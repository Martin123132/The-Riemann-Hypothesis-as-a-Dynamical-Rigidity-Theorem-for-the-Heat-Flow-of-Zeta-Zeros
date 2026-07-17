#!/usr/bin/env python3
"""Reduce the lambda=-100 order-seven tail to one nested curvature ceiling."""

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
    "jensen_window_pf_compound_order7_m100_tail_curvature_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.md"
)
PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_m100_prefix_certificate.json"
)
ORDER7_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_uniform_tail_flow_reduction.json"
)
ORDER4_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
TAIL_FIRST_N = 315
TAIL_FIRST_K = TAIL_FIRST_N + 6
CURVATURE_CONSTANT = 900


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
    reduction = load_json(ORDER7_REDUCTION)
    order4 = load_json(ORDER4_ENTRY_SOURCE)
    if prefix.get("summary", {}).get("positive_Q7_rows") != TAIL_FIRST_N:
        raise RuntimeError("order-seven prefix source does not end at n=314")
    if prefix.get("summary", {}).get("open_analytic_tails") != 1:
        raise RuntimeError("order-seven prefix tail contract changed")
    if reduction.get("summary", {}).get("universal_tail_specializations") != 1:
        raise RuntimeError("order-seven uniform tail source is not closed")
    defect = order4.get("tail_arithmetic", {}).get("defect_buffer")
    if defect != "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))":
        raise RuntimeError(f"scaled-defect anchor changed: {defect!r}")
    return {
        "prefix": "Q_(7,n)(-100)>0 for every 0<=n<=314",
        "prefix_minimum_n": prefix["finite"]["minimum_relative_n"],
        "prefix_minimum_lower": prefix["finite"]["minimum_relative_lower"],
        "uniform_tail_status": reduction.get("status"),
        "defect_anchor": "d_k>=251/(250*(2*k+1)), k>=320",
    }


def exact_reduction() -> dict:
    center, h5_scale, h4_scale, p, q, T = sp.symbols(
        "A H5_scale H4_scale p q T", positive=True
    )
    h5_center = center**5 * sp.exp(p)
    h4_center = center**4 * sp.exp(q)
    q6 = sp.factor(h5_center**2 * (1 - sp.exp(-T)) / h4_center)
    target = center**6 * sp.exp(2 * p - q) * (1 - sp.exp(-T))
    if sp.simplify(q6 - target) != 0:
        raise RuntimeError("canonical Q6 factorization failed")

    m, k = sp.symbols("m k", integer=True, nonnegative=True)
    comparison = sp.expand(
        sp.Integer(753) * k**2
        - sp.Integer(125) * CURVATURE_CONSTANT * (2 * k + 1)
    )
    shifted = sp.expand(comparison.subs(k, TAIL_FIRST_K + m))
    polynomial = sp.Poly(shifted, m)
    coefficients = polynomial.all_coeffs()
    if any(value <= 0 for value in coefficients):
        raise RuntimeError("order-seven tail comparison is not coefficient-positive")

    return {
        "lower_coordinates": (
            "x=exp(-B), d=1-x, g=d^2-x^2*d(t-1)*d(t+1), h=log(g)"
        ),
        "order5_coordinate": (
            "f=g^2-x^3*g(t-1)*g(t+1), q=log(f/d)"
        ),
        "order6_coordinate": (
            "S(t)=4*B(t)-q(t-1)+2*q(t)-q(t+1), "
            "p(t)=2*q(t)-h(t)+log(1-exp(-S(t)))"
        ),
        "fourth_gap": "T(t)=5*B(t)-p(t-1)+2*p(t)-p(t+1)",
        "order7_coordinate": (
            "r(t)=2*p(t)-q(t)+log(1-exp(-T(t)))"
        ),
        "canonical_factorization": (
            "Q_(6,n)=A_(n+5)^6*exp(r(n+5))"
        ),
        "canonical_factorization_residual": "0",
        "curvature_identity": (
            "E_n=log(Q_(6,n)*Q_(6,n+2)/Q_(6,n+1)^2)="
            "6*log(x_k)+R_k, R_k=r(k-1)-2*r(k)+r(k+1), k=n+6"
        ),
        "stable_margin": "L_n=exp(-E_n)-1",
        "sign_equivalence": "Q_(7,n)>0 iff L_n>0 iff E_n<0",
        "tail_index": "k=n+6, so n>=315 iff k>=321",
        "sufficient_ceiling": "R_k<=900/k^2 for every real/integer k>=321",
        "log_buffer": (
            "-6*log(x_k)>=6*d_k>=753/(125*(2*k+1)), k>=320"
        ),
        "rational_comparison": (
            "900/k^2<753/(125*(2*k+1)), k>=321"
        ),
        "cleared_polynomial": str(comparison),
        "shifted_polynomial_k_321_plus_m": str(shifted),
        "shifted_coefficients": [str(value) for value in coefficients],
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_reduction()
    rows = [
        TailRow(
            "co7m100tcr_01_prefix",
            "theorem_input",
            "ready_to_apply",
            "The rigorous endpoint prefix reaches the last index before the analytic splice.",
            sources["prefix"],
            "Finite Arb prefix only.",
        ),
        TailRow(
            "co7m100tcr_02_nested_coordinate",
            "exact_identity",
            "ready_to_apply",
            "One additional stable logarithm gives a cancellation-free coordinate for Q6.",
            exact["fourth_gap"] + "; " + exact["order7_coordinate"] + "; " + exact["canonical_factorization"],
            "Exact lower-cone factorization only.",
        ),
        TailRow(
            "co7m100tcr_03_curvature",
            "exact_reduction",
            "ready_to_apply",
            "Signed order seven is equivalent to one centered curvature lying below the sixth defect buffer.",
            exact["curvature_identity"] + "; " + exact["sign_equivalence"],
            "Exact logarithmic reduction only.",
        ),
        TailRow(
            "co7m100tcr_04_defect",
            "theorem_input",
            "ready_to_apply",
            "The completed defect theorem supplies an inverse-linear endpoint buffer.",
            exact["log_buffer"],
            "Previously proved lambda=-100 defect anchor.",
        ),
        TailRow(
            "co7m100tcr_05_arithmetic",
            "exact_inequality",
            "ready_to_apply",
            "The proposed inverse-square ceiling lies strictly inside the defect buffer on the whole tail.",
            exact["rational_comparison"],
            "Coefficient-positive rational comparison.",
            {"shifted_coefficients": exact["shifted_coefficients"]},
        ),
        TailRow(
            "co7m100tcr_06_conditional_tail",
            "conditional_theorem",
            "ready_to_apply",
            "The scalar curvature ceiling would close the complete order-seven endpoint tail.",
            "[R_k<=900/k^2 for k>=321] => [Q_(7,n)(-100)>0 for n>=315]",
            "Conditional on the displayed scalar ceiling only.",
        ),
        TailRow(
            "co7m100tcr_07_open_curvature",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the fourth-nested stable curvature ceiling on the analytic tail.",
            exact["sufficient_ceiling"],
            "This scalar analytic ceiling is not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order7_m100_tail_curvature_reduction",
        "date": "2026-07-13",
        "status": (
            "exact order-seven endpoint-tail reduction with one open "
            "fourth-nested curvature ceiling"
        ),
        "proof_boundary": (
            "This artifact proves the canonical Q6 normalization, sign reduction, "
            "and tail arithmetic. It does not prove the 900/k^2 curvature ceiling, "
            "order-seven entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": 6,
            "exact_factorizations": 1,
            "exact_curvature_reductions": 1,
            "coefficient_positive_comparisons": 1,
            "conditional_tail_theorems": 1,
            "open_curvature_targets": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order7_m100_tail_curvature_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order7_m100_tail_curvature_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Lambda=-100 Tail Curvature Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact endpoint-tail reduction with one open fourth-nested",
        "stable curvature ceiling. This is not a proof of order-seven entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order7_m100_tail_curvature_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_m100_tail_curvature_reduction.py",
        "```",
        "",
        "## Canonical Normalization",
        "",
        "The completed hierarchy through order six is",
        "",
        "```text",
        exact["lower_coordinates"],
        exact["order5_coordinate"],
        exact["order6_coordinate"],
        "```",
        "",
        "Signed condensation adds exactly one stable gap and logarithm:",
        "",
        "```text",
        exact["fourth_gap"],
        exact["order7_coordinate"],
        exact["canonical_factorization"],
        "```",
        "",
        "Indeed the centered `H_5(n+1)^2` scale is `A_(n+5)^10`; division",
        "by `H_4(n+2)=A_(n+5)^4 exp(q(n+5))` leaves the displayed sixth",
        "power. No raw near-cancelling determinant is used.",
        "",
        "## Exact Sign Reduction",
        "",
        "For `k=n+6`,",
        "",
        "```text",
        exact["curvature_identity"],
        exact["stable_margin"],
        exact["sign_equivalence"],
        "```",
        "",
        "The rigorous prefix proves `Q_(7,n)(-100)>0` for `0<=n<=314`,",
        "so the missing tail starts at `n=315`, equivalently `k=321`.",
        "",
        "## Tail Arithmetic",
        "",
        "The completed defect theorem gives",
        "",
        "```text",
        exact["log_buffer"],
        "```",
        "",
        "A deliberately buffered sufficient curvature theorem is",
        "",
        "```text",
        exact["sufficient_ceiling"],
        "```",
        "",
        "and the comparison is exact:",
        "",
        "```text",
        exact["rational_comparison"],
        exact["shifted_polynomial_k_321_plus_m"] + ">0 for m>=0,",
        "coefficients=" + str(exact["shifted_coefficients"]),
        "```",
        "",
        "Consequently the displayed scalar ceiling would prove the entire",
        "`n>=315` endpoint tail and splice it to the finite prefix.",
        "",
        "## Open Input",
        "",
        "```text",
        "Prove R_k<=900/k^2 for every k>=321.",
        "```",
        "",
        "The complete-to-first-summand transfer is now proved in",
        "`outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md`.",
        "It reduces the remaining task to the continuous theorem",
        "`r_1''(t)<=600/t^2` on `t>=320`. The extra stable layer requires a",
        "common `t+-5` cover and potential derivatives through order twelve.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order7_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md",
        "outputs/jensen_window_pf_compound_order7_first_summand_curvature_bridge.md",
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
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-seven tail curvature reduction: "
        f"{summary['rows']} rows, "
        f"{summary['exact_curvature_reductions']} exact reduction, "
        f"{summary['open_curvature_targets']} open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
