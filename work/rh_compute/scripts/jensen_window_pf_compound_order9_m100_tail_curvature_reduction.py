#!/usr/bin/env python3
"""Reduce the lambda=-100 order-nine tail to one nested curvature ceiling."""

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
    "jensen_window_pf_compound_order9_m100_tail_curvature_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md"
)
PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_prefix_certificate.json"
)
ORDER9_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_uniform_tail_flow_reduction.json"
)
ORDER4_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
TAIL_FIRST_N = 1241
TAIL_FIRST_K = TAIL_FIRST_N + 8
CURVATURE_CONSTANT = 4900


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
    reduction = load_json(ORDER9_REDUCTION)
    order4 = load_json(ORDER4_ENTRY_SOURCE)
    if prefix.get("summary", {}).get("positive_Q9_rows") != TAIL_FIRST_N:
        raise RuntimeError("order-nine prefix source does not end at n=1240")
    if prefix.get("summary", {}).get("open_analytic_tails") != 1:
        raise RuntimeError("order-nine prefix tail contract changed")
    if reduction.get("summary", {}).get("universal_tail_specializations") != 1:
        raise RuntimeError("order-nine uniform tail source is not closed")
    defect = order4.get("tail_arithmetic", {}).get("defect_buffer")
    if defect != "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))":
        raise RuntimeError(f"scaled-defect anchor changed: {defect!r}")
    return {
        "prefix": "Q_(9,n)(-100)>0 for every 0<=n<=1240",
        "prefix_minimum_n": prefix["finite"]["minimum_relative_n"],
        "prefix_minimum_lower": prefix["finite"]["minimum_relative_lower"],
        "uniform_tail_status": reduction.get("status"),
        "defect_anchor": "d_k>=251/(250*(2*k+1)), k>=320",
    }


def exact_reduction() -> dict:
    center, s, r, v = sp.symbols("A s r V", positive=True)
    q7_center = center**7 * sp.exp(s)
    q6_center = center**6 * sp.exp(r)
    q8 = sp.factor(q7_center**2 * (1 - sp.exp(-v)) / q6_center)
    target = center**8 * sp.exp(2 * s - r) * (1 - sp.exp(-v))
    if sp.simplify(q8 - target) != 0:
        raise RuntimeError("canonical Q8 factorization failed")

    m, k = sp.symbols("m k", integer=True, nonnegative=True)
    comparison = sp.expand(
        sp.Integer(1004) * k**2
        - sp.Integer(125) * CURVATURE_CONSTANT * (2 * k + 1)
    )
    shifted = sp.expand(comparison.subs(k, TAIL_FIRST_K + m))
    coefficients = sp.Poly(shifted, m).all_coeffs()
    if any(value <= 0 for value in coefficients):
        raise RuntimeError("order-nine tail comparison is not coefficient-positive")

    return {
        "completed_coordinates": (
            "U(t)=6*B(t)-r(t-1)+2*r(t)-r(t+1); "
            "s(t)=2*r(t)-p(t)+log(1-exp(-U(t)))"
        ),
        "sixth_gap": "V(t)=7*B(t)-s(t-1)+2*s(t)-s(t+1)",
        "order9_coordinate": "w(t)=2*s(t)-r(t)+log(1-exp(-V(t)))",
        "canonical_factorization": "Q_(8,n)=A_(n+7)^8*exp(w(n+7))",
        "canonical_factorization_residual": "0",
        "curvature_identity": (
            "E_n=log(Q_(8,n)*Q_(8,n+2)/Q_(8,n+1)^2)="
            "8*log(x_k)+Y_k, Y_k=w(k-1)-2*w(k)+w(k+1), k=n+8"
        ),
        "stable_margin": "M_n=exp(-E_n)-1",
        "sign_equivalence": "Q_(9,n)>0 iff M_n>0 iff E_n<0",
        "tail_index": "k=n+8, so n>=1241 iff k>=1249",
        "sufficient_ceiling": "Y_k<=4900/k^2 for every real/integer k>=1249",
        "log_buffer": (
            "-8*log(x_k)>=8*d_k>=1004/(125*(2*k+1)), k>=320"
        ),
        "rational_comparison": (
            "4900/k^2<1004/(125*(2*k+1)), k>=1249"
        ),
        "cleared_polynomial": str(comparison),
        "shifted_polynomial_k_1249_plus_m": str(shifted),
        "shifted_coefficients": [str(value) for value in coefficients],
        "continuous_derivative_requirement": "H derivatives through order sixteen",
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_reduction()
    rows = [
        TailRow(
            "co9m100tcr_01_prefix",
            "theorem_input",
            "ready_to_apply",
            "The rigorous endpoint prefix reaches the last index before the analytic splice.",
            sources["prefix"],
            "Finite Arb prefix only.",
        ),
        TailRow(
            "co9m100tcr_02_nested_coordinate",
            "exact_identity",
            "ready_to_apply",
            "One additional stable logarithm gives a cancellation-free coordinate for Q8.",
            exact["sixth_gap"] + "; " + exact["order9_coordinate"] + "; " + exact["canonical_factorization"],
            "Exact lower-cone factorization only.",
        ),
        TailRow(
            "co9m100tcr_03_curvature",
            "exact_reduction",
            "ready_to_apply",
            "Signed order nine is equivalent to one centered curvature lying below the eighth defect buffer.",
            exact["curvature_identity"] + "; " + exact["sign_equivalence"],
            "Exact logarithmic reduction only.",
        ),
        TailRow(
            "co9m100tcr_04_defect",
            "theorem_input",
            "ready_to_apply",
            "The completed defect theorem supplies an inverse-linear endpoint buffer.",
            exact["log_buffer"],
            "Previously proved lambda=-100 defect anchor.",
        ),
        TailRow(
            "co9m100tcr_05_arithmetic",
            "exact_inequality",
            "ready_to_apply",
            "The inverse-square ceiling lies strictly inside the defect buffer on the whole tail.",
            exact["rational_comparison"],
            "Coefficient-positive rational comparison.",
            {"shifted_coefficients": exact["shifted_coefficients"]},
        ),
        TailRow(
            "co9m100tcr_06_conditional_tail",
            "conditional_theorem",
            "ready_to_apply",
            "The scalar curvature ceiling would close the complete order-nine endpoint tail.",
            "[Y_k<=4900/k^2 for k>=1249] => [Q_(9,n)(-100)>0 for n>=1241]",
            "Conditional on the displayed scalar ceiling only.",
        ),
        TailRow(
            "co9m100tcr_07_open_curvature",
            "open_handoff",
            "not_ready_to_apply",
            "Prove the sixth-nested stable curvature ceiling on the analytic tail.",
            exact["sufficient_ceiling"],
            "This scalar analytic ceiling is not proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_m100_tail_curvature_reduction",
        "date": "2026-07-13",
        "status": (
            "exact order-nine endpoint-tail reduction with one open "
            "sixth-nested curvature ceiling"
        ),
        "proof_boundary": (
            "This artifact proves the canonical Q8 normalization, sign reduction, "
            "and tail arithmetic. It does not prove the 4900/k^2 curvature ceiling, "
            "order-nine entry, PF-infinity, RH, or Lambda<=0."
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
            "required_top_derivative_order": 16,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_m100_tail_curvature_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_m100_tail_curvature_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    lines = [
        "# Jensen-Window PF Compound Order-Nine Lambda=-100 Tail Curvature Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact endpoint-tail reduction with one open sixth-nested",
        "stable curvature ceiling. This is not a proof of order-nine entry,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_tail_curvature_reduction.py",
        "```",
        "",
        "## Canonical Normalization",
        "",
        "```text",
        exact["completed_coordinates"],
        exact["sixth_gap"],
        exact["order9_coordinate"],
        exact["canonical_factorization"],
        "```",
        "",
        "## Exact Sign Reduction",
        "",
        "```text",
        exact["curvature_identity"],
        exact["stable_margin"],
        exact["sign_equivalence"],
        exact["tail_index"],
        "```",
        "",
        "## Tail Arithmetic",
        "",
        "```text",
        exact["log_buffer"],
        exact["sufficient_ceiling"],
        exact["rational_comparison"],
        exact["shifted_polynomial_k_1249_plus_m"] + ">0 for m>=0,",
        "coefficients=" + str(exact["shifted_coefficients"]),
        "```",
        "",
        "The scalar ceiling would prove the complete `n>=1241` endpoint tail",
        "and splice it to the finite prefix.",
        "",
        "## Open Input",
        "",
        "```text",
        "Prove Y_k<=4900/k^2 for every k>=1249.",
        "```",
        "",
        "A continuous proof requires a common `t+-7` cover and `H` derivatives",
        "through order sixteen before tent integration.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order9_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-nine tail curvature reduction: "
        f"{summary['rows']} rows, "
        f"{summary['coefficient_positive_comparisons']} positive comparison, "
        f"{summary['open_curvature_targets']} open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
