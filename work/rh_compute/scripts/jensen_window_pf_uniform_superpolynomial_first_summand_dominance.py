#!/usr/bin/env python3
"""Upgrade compact-heat first-summand dominance to superpolynomial decay."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_superpolynomial_first_summand_dominance.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md"
)
M100_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json"
)
TRANSFER_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_lambda0_first_summand_dominance_transfer.json"
)


@dataclass(frozen=True)
class DominanceRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_diagnostics() -> dict:
    m100 = json.loads(M100_SOURCE.read_text(encoding="utf-8"))
    transfer = json.loads(TRANSFER_SOURCE.read_text(encoding="utf-8"))
    if m100.get("summary", {}).get("full_tail_power") != 6:
        raise RuntimeError("lambda=-100 first-summand source is not closed")
    if transfer.get("summary", {}).get("uniform_heat_intervals") != 1:
        raise RuntimeError("compact-heat monotonicity transfer is not closed")
    return {
        "m100_status": m100.get("status"),
        "m100_tail": m100.get("diagnostics", {}).get("full_tail_relative_bound"),
        "transfer_status": transfer.get("status"),
        "transfer_monotonicity": transfer.get("exact", {}).get("heat_monotonicity"),
        "transfer_uniform_bound": transfer.get("exact", {}).get("uniform_dominance"),
    }


def asymptotic_diagnostics() -> dict:
    k = sp.symbols("k", positive=True)
    p = sp.symbols("p", positive=True)
    L = sp.log(k)
    low_exponent = (
        -sp.Rational(8, 5) * k / (L + sp.Rational(4, 5))
        + sp.Rational(5, 2) * L
        + 1
        + sp.pi * (sp.exp(sp.Rational(2, 5)) - 1) * sp.sqrt(k)
    )
    high_exponent = sp.log(17) - 3 * sp.pi * sp.sqrt(k)
    low_scaled_limit = sp.limit((low_exponent + p * L) * L / k, k, sp.oo)
    high_scaled_limit = sp.limit(
        (high_exponent + p * L) / sp.sqrt(k), k, sp.oo
    )
    if low_scaled_limit != -sp.Rational(8, 5):
        raise RuntimeError("low-region superpolynomial limit changed")
    if high_scaled_limit != -3 * sp.pi:
        raise RuntimeError("high-region superpolynomial limit changed")
    return {
        "adaptive_split": "a(k)=log(k)/8",
        "low_region_bound": (
            "epsilon(0)*P_(k,100)(u<a(k))<=exp(B_low(k))"
        ),
        "low_exponent": (
            "B_low(k)=-8*k/(5*(log(k)+4/5))+(5/2)*log(k)+1+"
            "pi*(exp(2/5)-1)*sqrt(k)"
        ),
        "low_limit": (
            "lim_(k->infinity) (B_low(k)+p*log(k))*log(k)/k=-8/5"
        ),
        "high_region_bound": (
            "epsilon(a(k))<=17*exp(-3*pi*sqrt(k))"
        ),
        "high_limit": (
            "lim_(k->infinity) (log(17)-3*pi*sqrt(k)+p*log(k))/sqrt(k)=-3*pi"
        ),
        "m100_superpolynomial": (
            "for every p>0 there exists K_p such that "
            "0<=delta_k(100)<=2*k^(-p), k>=K_p"
        ),
        "compact_heat_superpolynomial": (
            "for every p>0 there exists K_p such that sup_(0<=T<=100) "
            "delta_k(T)<=2*k^(-p), k>=K_p"
        ),
        "log_tail": "e_k(T)=log(1+delta_k(T)), 0<=e_k(T)<=delta_k(T)",
        "difference_bound": (
            "|nabla^m e_k(T)|<=2^m*max_(0<=j<=m)e_(k+j)(T)"
        ),
        "local_difference_consequence": (
            "for every p>0 and fixed m, sup_(0<=T<=100)|nabla^m e_k(T)|="
            "O_p,m(k^(-p))"
        ),
        "determinant_consequence": (
            "the complete-kernel logarithmic ratio correction is o(k^-N) "
            "uniformly for every fixed N and every fixed local difference order"
        ),
        "symbolic_low_scaled_limit": str(low_scaled_limit),
        "symbolic_high_scaled_limit": str(high_scaled_limit),
    }


def build_artifact() -> dict:
    sources = source_diagnostics()
    exact = asymptotic_diagnostics()
    rows = [
        DominanceRow(
            id="usfsd_01_low_region_asymptotic",
            role="exact_asymptotic_lemma",
            readiness="ready_to_apply",
            claim="The inherited low-region probability exponent dominates every fixed logarithmic power target.",
            formula=exact["low_limit"],
            proof_boundary="Non-effective asymptotic consequence of the proved lambda=-100 low-region bound.",
            diagnostics={"bound": exact["low_region_bound"], "exponent": exact["low_exponent"]},
        ),
        DominanceRow(
            id="usfsd_02_high_region_asymptotic",
            role="exact_asymptotic_lemma",
            readiness="ready_to_apply",
            claim="The inherited pointwise high-region theta-tail bound dominates every fixed inverse power.",
            formula=exact["high_limit"],
            proof_boundary="Non-effective asymptotic consequence of the proved high-region bound.",
            diagnostics={"bound": exact["high_region_bound"]},
        ),
        DominanceRow(
            id="usfsd_03_m100_superpolynomial",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="At lambda=-100, the complete-to-first-summand relative moment tail is superpolynomially small.",
            formula=exact["m100_superpolynomial"],
            proof_boundary="All sufficiently large indices for each fixed power; thresholds are not effectivized.",
        ),
        DominanceRow(
            id="usfsd_04_uniform_heat_transfer",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Heat-parameter monotonicity transfers superpolynomial first-summand dominance uniformly to the full compact interval.",
            formula=exact["compact_heat_superpolynomial"],
            proof_boundary="Uniform coefficient-tail theorem only.",
            diagnostics={"monotonicity": sources["transfer_monotonicity"]},
        ),
        DominanceRow(
            id="usfsd_05_local_log_differences",
            role="exact_perturbation_theorem",
            readiness="ready_to_apply",
            claim="Every fixed local finite difference of the logarithmic higher-summand correction is uniformly superpolynomial.",
            formula=exact["local_difference_consequence"],
            proof_boundary="Fixed finite difference order only; sufficient for every fixed determinant size.",
            diagnostics={"difference_bound": exact["difference_bound"]},
        ),
        DominanceRow(
            id="usfsd_06_uniform_ratio_handoff",
            role="exact_theorem_handoff",
            readiness="ready_to_apply",
            claim="Higher theta summands are negligible to every fixed order in the compact-heat order-four ratio expansion.",
            formula=exact["determinant_consequence"],
            proof_boundary="Closes only the higher-summand remainder; the first-summand heat-tilt asymptotic remains separate.",
        ),
    ]
    source_paths = (M100_SOURCE, TRANSFER_SOURCE)
    return {
        "kind": "jensen_window_pf_uniform_superpolynomial_first_summand_dominance",
        "date": "2026-07-13",
        "status": "uniform compact-heat superpolynomial first-summand dominance theorem",
        "proof_boundary": (
            "This artifact proves non-effective superpolynomial higher-theta "
            "suppression and all fixed local log-difference bounds uniformly for "
            "-100<=lambda<=0. It does not prove the first-summand heat-tilt ratio "
            "expansion, uniform order-four positivity, forward order-four invariance, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            str(path.relative_to(REPO_ROOT)).replace("\\", "/") for path in source_paths
        ],
        "source_hashes": {
            str(path.relative_to(REPO_ROOT)).replace("\\", "/"): sha256(path)
            for path in source_paths
        },
        "source_diagnostics": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_rows": len(rows),
            "ready_to_apply_rows": len(rows),
            "open_analytic_rows": 0,
            "superpolynomial_tail_theorems": 1,
            "uniform_heat_transfers": 1,
            "local_log_difference_theorems": 1,
            "higher_theta_handoffs_closed": 1,
        },
        "remaining_target": (
            "prove the first-summand heat-tilt finite-difference estimate uniformly "
            "for 0<=T<=100"
        ),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_uniform_superpolynomial_first_summand_dominance.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_uniform_superpolynomial_first_summand_dominance.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Uniform Superpolynomial First-Summand Dominance",
        "",
        "Date: 2026-07-13",
        "",
        "Status: uniform compact-heat superpolynomial first-summand dominance",
        "theorem. This is not a proof of uniform order-four positivity, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.json",
        "python work/rh_compute/scripts/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_uniform_superpolynomial_first_summand_dominance.py",
        "```",
        "",
        "## Low And High Regions",
        "",
        "The proved lambda=-100 split uses `a(k)=log(k)/8`. Its low-region",
        "exponent satisfies",
        "",
        "```text",
        exact["low_exponent"],
        exact["low_limit"] + ".",
        "```",
        "",
        "The high-region estimate satisfies",
        "",
        "```text",
        exact["high_region_bound"],
        exact["high_limit"] + ".",
        "```",
        "",
        "Both limits are negative for every fixed `p>0`, so the two regions give",
        "",
        "```text",
        exact["m100_superpolynomial"] + ".",
        "```",
        "",
        "## Uniform Heat Transfer",
        "",
        "The covariance monotonicity theorem gives `delta_k(T)<=delta_k(100)`.",
        "Consequently",
        "",
        "```text",
        exact["compact_heat_superpolynomial"] + ".",
        "```",
        "",
        "For `e_k(T)=log(1+delta_k(T))`, the elementary stencil bound",
        "",
        "```text",
        exact["difference_bound"],
        "```",
        "",
        "proves every fixed local log difference is uniformly superpolynomial.",
        "Thus higher theta summands are negligible to every fixed order in the",
        "uniform order-four ratio expansion. The remaining asymptotic problem is",
        "only the first-summand heat tilt.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "outputs/jensen_window_pf_lambda0_first_summand_dominance_transfer.md",
        "outputs/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "proved uniform superpolynomial first-summand dominance: "
        f"{summary['rows']} rows, {summary['exact_rows']} exact rows, "
        f"{summary['superpolynomial_tail_theorems']} superpolynomial theorem, "
        f"{summary['local_log_difference_theorems']} local-difference theorem, "
        f"{summary['open_analytic_rows']} open rows"
    )


if __name__ == "__main__":
    main()
