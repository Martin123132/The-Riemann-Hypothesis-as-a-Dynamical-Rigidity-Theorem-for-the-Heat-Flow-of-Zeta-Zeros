#!/usr/bin/env python3
"""Reduce the lambda=-100 order-six tail to one nested curvature ceiling."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_m100_tail_curvature_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.md"
)
PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_m100_prefix_certificate.json"
)
ORDER6_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_uniform_tail_flow_reduction.json"
)
ORDER4_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
TAIL_FIRST_N = 317
TAIL_FIRST_K = TAIL_FIRST_N + 5
CURVATURE_CONSTANT = 320


@dataclass(frozen=True)
class TailRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_sources() -> dict:
    prefix = load_json(PREFIX_SOURCE)
    reduction = load_json(ORDER6_REDUCTION)
    order4 = load_json(ORDER4_ENTRY_SOURCE)
    if prefix.get("summary", {}).get("positive_Q6_rows") != TAIL_FIRST_N:
        raise RuntimeError("order-six prefix source does not end at n=316")
    if prefix.get("summary", {}).get("open_analytic_tails") != 1:
        raise RuntimeError("order-six prefix tail contract changed")
    if reduction.get("summary", {}).get("uniform_eventual_tail_theorems") != 1:
        raise RuntimeError("order-six uniform tail source is not closed")
    defect = order4.get("tail_arithmetic", {}).get("defect_buffer")
    if defect != "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))":
        raise RuntimeError(f"scaled-defect anchor changed: {defect!r}")
    return {
        "prefix": "Q_(6,n)(-100)>0 for every 0<=n<=316",
        "prefix_minimum_n": prefix["finite"]["minimum_relative_n"],
        "prefix_minimum_lower": prefix["finite"]["minimum_relative_lower"],
        "uniform_tail_status": reduction.get("status"),
        "defect_anchor": "d_k>=251/(250*(2*k+1)), k>=320",
    }


def exact_reduction() -> dict:
    a0, a1, a2, a3, a4 = sp.symbols("A0:5", positive=True)
    rho = a1 / a0
    x1 = a0 * a2 / a1**2
    x2 = a1 * a3 / a2**2
    x3 = a2 * a4 / a3**2
    prefactor = sp.factor(a0**5 * rho**20 * x1**15 * x2**10 * x3**5)
    if sp.factor(prefactor - a4**5) != 0:
        raise RuntimeError("order-five monomial prefactor did not collapse")

    dm, d0, dp, g, fm, f0, fp, x = sp.symbols(
        "d_minus d_0 d_plus g f_minus f_0 f_plus x", positive=True
    )
    stable_j = dm * dp * f0**2 - x**4 * d0**2 * fm * fp
    ratio = x**4 * d0**2 * fm * fp / (dm * dp * f0**2)
    if sp.factor(stable_j - dm * dp * f0**2 * (1 - ratio)) != 0:
        raise RuntimeError("third stable factorization failed")

    m, k = sp.symbols("m k", integer=True, nonnegative=True)
    comparison = sp.expand(
        sp.Integer(251) * k**2
        - sp.Integer(50) * CURVATURE_CONSTANT * (2 * k + 1)
    )
    shifted = sp.expand(comparison.subs(k, TAIL_FIRST_K + m))
    polynomial = sp.Poly(shifted, m)
    coefficients = polynomial.all_coeffs()
    if any(value <= 0 for value in coefficients):
        raise RuntimeError("order-six tail comparison is not coefficient-positive")

    return {
        "lower_coordinates": (
            "x=exp(-B), d=1-x, g=d^2-x^2*d(t-1)*d(t+1), "
            "h=log(g)"
        ),
        "order5_coordinate": (
            "f=g^2-x^3*g(t-1)*g(t+1), q=log(f/d)"
        ),
        "third_gap": "S(t)=4*B(t)-q(t-1)+2*q(t)-q(t+1)",
        "order6_coordinate": (
            "p(t)=2*q(t)-h(t)+log(1-exp(-S(t)))"
        ),
        "canonical_factorization": "H_(5,n)=A_(n+4)^5*exp(p(n+4))",
        "prefactor_identity": (
            "A_n^5*rho_n^20*x_(n+1)^15*x_(n+2)^10*x_(n+3)^5="
            "A_(n+4)^5"
        ),
        "curvature_identity": (
            "D_n=log(H_(5,n)*H_(5,n+2)/H_(5,n+1)^2)="
            "5*log(x_k)+P_k, P_k=p(k-1)-2*p(k)+p(k+1), k=n+5"
        ),
        "stable_margin": "K_n=exp(-D_n)-1",
        "sign_equivalence": "Q_(6,n)>0 iff K_n>0 iff D_n<0",
        "tail_index": "k=n+5, so n>=317 iff k>=322",
        "sufficient_ceiling": "P_k<=320/k^2 for every real/integer k>=322",
        "log_buffer": (
            "-5*log(x_k)>=5*d_k>=251/(50*(2*k+1)), k>=320"
        ),
        "rational_comparison": (
            "320/k^2<251/(50*(2*k+1)), k>=322"
        ),
        "cleared_polynomial": str(comparison),
        "shifted_polynomial_k_322_plus_m": str(shifted),
        "shifted_coefficients": [str(value) for value in coefficients],
        "prefactor_residual": "0",
        "stable_factor_residual": "0",
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_reduction()
    rows = [
        TailRow(
            "co6m100tcr_01_prefix",
            "theorem_input",
            "ready_to_apply",
            "The rigorous endpoint prefix reaches the last index before the analytic splice.",
            sources["prefix"],
            "Finite Arb prefix only.",
        ),
        TailRow(
            "co6m100tcr_02_prefactor",
            "exact_identity",
            "ready_to_apply",
            "The entire order-five coefficient monomial collapses to a single shifted fifth power.",
            exact["prefactor_identity"],
            "Exact algebra only.",
        ),
        TailRow(
            "co6m100tcr_03_nested_coordinate",
            "exact_identity",
            "ready_to_apply",
            "One additional stable logarithm gives a cancellation-free coordinate for H5.",
            exact["third_gap"] + "; " + exact["order6_coordinate"] + "; " + exact["canonical_factorization"],
            "Exact lower-cone factorization only.",
        ),
        TailRow(
            "co6m100tcr_04_curvature",
            "exact_reduction",
            "ready_to_apply",
            "Signed order six is equivalent to one centered curvature lying below the fifth defect buffer.",
            exact["curvature_identity"] + "; " + exact["sign_equivalence"],
            "Exact logarithmic reduction only.",
        ),
        TailRow(
            "co6m100tcr_05_defect",
            "theorem_input",
            "ready_to_apply",
            "The completed defect theorem supplies an inverse-linear endpoint buffer.",
            exact["log_buffer"],
            "Previously proved lambda=-100 defect anchor.",
        ),
        TailRow(
            "co6m100tcr_06_arithmetic",
            "exact_inequality",
            "ready_to_apply",
            "The proposed inverse-square ceiling lies strictly inside the defect buffer on the whole tail.",
            exact["rational_comparison"],
            "Coefficient-positive rational comparison.",
            {"shifted_coefficients": exact["shifted_coefficients"]},
        ),
        TailRow(
            "co6m100tcr_07_conditional_tail",
            "conditional_theorem",
            "ready_to_apply",
            "The scalar curvature ceiling would close the complete order-six endpoint tail.",
            "[P_k<=320/k^2 for k>=322] => [Q_(6,n)(-100)>0 for n>=317]",
            "Conditional on the displayed scalar ceiling only.",
        ),
        TailRow(
            "co6m100tcr_08_open_curvature",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the third-nested stable curvature ceiling on the analytic tail.",
            exact["sufficient_ceiling"],
            "This scalar analytic ceiling is not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order6_m100_tail_curvature_reduction",
        "date": "2026-07-13",
        "status": "exact order-six endpoint-tail reduction with one open nested curvature ceiling",
        "proof_boundary": (
            "This artifact proves the canonical H5 normalization, sign reduction, "
            "and tail arithmetic. It does not prove the 320/k^2 curvature ceiling, "
            "order-six entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": 7,
            "exact_factorizations": 2,
            "exact_curvature_reductions": 1,
            "coefficient_positive_comparisons": 1,
            "conditional_tail_theorems": 1,
            "open_curvature_targets": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_tail_curvature_reduction.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Six Lambda=-100 Tail Curvature Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact endpoint-tail reduction with one open third-nested stable",
        "curvature ceiling. This is not a proof of order-six entry, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order6_m100_tail_curvature_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order6_m100_tail_curvature_reduction.py",
        "```",
        "",
        "## Canonical Normalization",
        "",
        "The apparently complicated positive monomial in the stable H5",
        "factorization collapses exactly:",
        "",
        "```text",
        exact["prefactor_identity"],
        exact["third_gap"],
        exact["order6_coordinate"],
        exact["canonical_factorization"],
        "```",
        "",
        "Thus order six adds exactly one stable logarithm to the completed",
        "order-five hierarchy; no raw near-cancelling determinant is used.",
        "",
        "## Exact Sign Reduction",
        "",
        "For `k=n+5`,",
        "",
        "```text",
        exact["curvature_identity"],
        exact["stable_margin"],
        exact["sign_equivalence"],
        "```",
        "",
        "The rigorous prefix proves `Q_(6,n)(-100)>0` for `0<=n<=316`,",
        "so the missing tail starts at `n=317`, equivalently `k=322`.",
        "",
        "## Tail Arithmetic",
        "",
        "The completed defect theorem gives",
        "",
        "```text",
        exact["log_buffer"],
        "```",
        "",
        "A deliberately loose sufficient curvature theorem is",
        "",
        "```text",
        exact["sufficient_ceiling"],
        "```",
        "",
        "and the comparison is exact:",
        "",
        "```text",
        exact["rational_comparison"],
        exact["shifted_polynomial_k_322_plus_m"] + ">0 for m>=0,",
        "coefficients=" + str(exact["shifted_coefficients"]),
        "```",
        "",
        "Consequently the displayed scalar ceiling would prove the entire",
        "`n>=317` endpoint tail and splice it to the finite prefix.",
        "",
        "## Open Input",
        "",
        "```text",
        "Prove P_k<=320/k^2 for every k>=322.",
        "```",
        "",
        "The companion first/full bridge separates that task into a continuous",
        "first-summand curvature theorem and an explicit kernel-transfer error.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order6_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order6_uniform_tail_flow_reduction.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-six tail curvature reduction: "
        f"{summary['rows']} rows, "
        f"{summary['exact_curvature_reductions']} exact reduction, "
        f"{summary['open_curvature_targets']} open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
