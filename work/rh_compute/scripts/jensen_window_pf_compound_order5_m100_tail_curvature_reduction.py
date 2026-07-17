#!/usr/bin/env python3
"""Reduce the lambda=-100 order-five tail to one scalar curvature ceiling."""

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
    "jensen_window_pf_compound_order5_m100_tail_curvature_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.md"
)
PREFIX_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_m100_prefix_certificate.json"
)
ORDER5_REDUCTION = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order5_uniform_tail_flow_reduction.json"
)
ORDER4_ENTRY_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_m100_entry_certificate.json"
)
TAIL_FIRST_N = 317
TAIL_FIRST_K = TAIL_FIRST_N + 4
CURVATURE_CONSTANT = 100


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
    reduction = load_json(ORDER5_REDUCTION)
    order4_entry = load_json(ORDER4_ENTRY_SOURCE)
    if prefix.get("summary", {}).get("positive_H5_rows") != TAIL_FIRST_N:
        raise RuntimeError("order-five prefix source does not end at n=316")
    if prefix.get("summary", {}).get("open_analytic_tails") != 1:
        raise RuntimeError("order-five prefix tail contract changed")
    if reduction.get("summary", {}).get("uniform_eventual_tail_theorems") != 1:
        raise RuntimeError("order-five uniform tail source is not closed")
    if reduction.get("summary", {}).get("cooperative_flow_theorems") != 1:
        raise RuntimeError("order-five flow source is not closed")
    defect_buffer = order4_entry.get("tail_arithmetic", {}).get("defect_buffer")
    if defect_buffer != "-3*log(x_k)>3*d_k>=753/(250*(2*k+1))":
        raise RuntimeError(f"scaled-defect anchor changed: {defect_buffer!r}")
    finite = prefix.get("finite", {})
    return {
        "prefix_status": prefix.get("status"),
        "prefix_theorem": "H_(5,n)(-100)>0 for every 0<=n<=316",
        "uniform_tail_status": reduction.get("status"),
        "defect_anchor": "d_k>=251/(250*(2*k+1)), k>=320",
        "source_defect_buffer": defect_buffer,
        "boundary_relative_ball": finite.get("minimum_relative_ball"),
        "boundary_relative_lower": finite.get("minimum_relative_lower"),
    }


def exact_reduction() -> dict:
    m, k = sp.symbols("m k", integer=True, nonnegative=True)
    x = sp.symbols("x", positive=True)
    d_minus, d_mid, d_plus = sp.symbols(
        "d_minus d_mid d_plus", positive=True
    )
    f0, f1, f2 = sp.symbols("F0 F1 F2", positive=True)
    stable = d_minus * d_plus * f1**2 - x**4 * d_mid**2 * f0 * f2
    ratio = sp.factor(
        d_minus * d_plus * f1**2 / (x**4 * d_mid**2 * f0 * f2)
    )
    if sp.factor(stable - x**4 * d_mid**2 * f0 * f2 * (ratio - 1)) != 0:
        raise RuntimeError("stable ratio identity failed")

    comparison = sp.expand(
        sp.Integer(502) * k**2
        - sp.Integer(125) * CURVATURE_CONSTANT * (2 * k + 1)
    )
    shifted = sp.expand(comparison.subs(k, TAIL_FIRST_K + m))
    polynomial = sp.Poly(shifted, m)
    if any(coefficient <= 0 for coefficient in polynomial.all_coeffs()):
        raise RuntimeError("tail comparison polynomial is not coefficient-positive")
    return {
        "stable_margin": (
            "J_n=d_(n+3)*d_(n+5)*F_(n+1)^2-"
            "x_(n+4)^4*d_(n+4)^2*F_n*F_(n+2)"
        ),
        "tail_index": "k=n+4, so n>=317 iff k>=321",
        "curvature_quantity": (
            "C_n=log(F_n*F_(n+2)/F_(n+1)^2)+"
            "log(d_k^2/(d_(k-1)*d_(k+1))), k=n+4"
        ),
        "curvature_identity": (
            "C_n=Delta^2 log(F_n)-Delta^2 log(d_(n+3))"
        ),
        "sign_equivalence": "J_n>0 iff C_n<-4*log(x_k)",
        "sufficient_ceiling": "C_n<=100/k^2 for every k=n+4>=321",
        "log_buffer": (
            "-4*log(x_k)>=4*d_k>=502/(125*(2*k+1)), k>=320"
        ),
        "rational_comparison": (
            "100/k^2<502/(125*(2*k+1)), k>=321"
        ),
        "cleared_polynomial": str(comparison),
        "shifted_polynomial_k_321_plus_n": str(shifted),
        "shifted_coefficients": [
            str(coefficient) for coefficient in polynomial.all_coeffs()
        ],
        "stable_ratio_residual": "0",
    }


def conclusions(exact: dict, sources: dict) -> dict:
    if exact["shifted_coefficients"] != ["502", "297284", "43689082"]:
        raise RuntimeError("tail comparison coefficients changed")
    return {
        "conditional_tail": (
            "[C_n<=100/(n+4)^2 for every n>=317] => "
            "[J_n(-100)>0 and H_(5,n)(-100)>0 for every n>=317]"
        ),
        "conditional_entry": (
            "prefix 0<=n<=316 plus the curvature tail => "
            "H_(5,n)(-100)>0 for every n>=0"
        ),
        "conditional_forward": (
            "the completed uniform tail and cooperative flow then imply "
            "H_(5,n)(lambda)>0 for every n>=0 and -100<=lambda<=0"
        ),
        "boundary_scout": (
            "at n=316, (n+4)^2*C_n=3.5869277550969014082... "
            "while the sufficient cap is 100"
        ),
        "prefix_relative": sources["boundary_relative_ball"],
        "open_target": "prove C_n<=100/(n+4)^2 for every n>=317 at lambda=-100",
    }


def build_artifact() -> dict:
    sources = validate_sources()
    exact = exact_reduction()
    theorem = conclusions(exact, sources)
    rows = [
        TailRow(
            id="co5m100tcr_01_prefix_input",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="The rigorous endpoint prefix already proves order-five entry through n=316.",
            formula=sources["prefix_theorem"],
            proof_boundary="Finite Arb prefix only.",
        ),
        TailRow(
            id="co5m100tcr_02_log_equivalence",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The stable order-five sign is exactly one scalar logarithmic inequality.",
            formula=exact["sign_equivalence"] + "; " + exact["curvature_quantity"],
            proof_boundary="Exact algebra inside the completed positive lower cone.",
        ),
        TailRow(
            id="co5m100tcr_03_curvature_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The scalar target is a difference of adjacent log-curvatures of F and d.",
            formula=exact["curvature_identity"],
            proof_boundary="Exact finite-difference identity.",
        ),
        TailRow(
            id="co5m100tcr_04_defect_buffer",
            role="theorem_input",
            readiness="ready_to_apply",
            claim="The existing scaled-defect anchor supplies an order-one-over-k logarithmic buffer.",
            formula=exact["log_buffer"],
            proof_boundary="Previously proved lambda=-100 all-index anchor.",
        ),
        TailRow(
            id="co5m100tcr_05_rational_comparison",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="The coarse order-one-over-k-squared curvature ceiling lies strictly below the defect buffer on the complete tail.",
            formula=exact["rational_comparison"],
            proof_boundary="Coefficient-positive rational arithmetic for k>=321.",
            diagnostics={
                "shifted_polynomial": exact["shifted_polynomial_k_321_plus_n"],
                "coefficients": exact["shifted_coefficients"],
            },
        ),
        TailRow(
            id="co5m100tcr_06_conditional_tail",
            role="conditional_theorem",
            readiness="conditional_on_open_input",
            claim="The scalar curvature ceiling would prove the complete missing lambda=-100 order-five tail.",
            formula=theorem["conditional_tail"],
            proof_boundary="Conditional on the displayed scalar ceiling only.",
        ),
        TailRow(
            id="co5m100tcr_07_conditional_completion",
            role="conditional_theorem",
            readiness="conditional_on_open_input",
            claim="Composing the prefix, tail, uniform eventual theorem, and cooperative flow would complete order five through lambda zero.",
            formula=theorem["conditional_entry"] + "; " + theorem["conditional_forward"],
            proof_boundary="Conditional; not PF-infinity, RH, or Lambda<=0.",
        ),
        TailRow(
            id="co5m100tcr_08_open_curvature",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="Prove one zeta-specific stable log-curvature ceiling on the analytic tail.",
            formula=theorem["open_target"],
            proof_boundary="This scalar analytic ceiling is not proved here.",
            diagnostics={"boundary_scout": theorem["boundary_scout"]},
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order5_m100_tail_curvature_reduction",
        "date": "2026-07-13",
        "status": (
            "exact order-five lambda=-100 tail reduction to one open stable "
            "log-curvature ceiling"
        ),
        "proof_boundary": (
            "This artifact proves the exact scalar reduction and sufficient "
            "constant comparison. It does not prove the curvature ceiling, "
            "the order-five analytic tail, all-shift order five, PF-infinity, "
            "RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md",
            "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
            "outputs/jensen_window_pf_compound_order4_m100_entry_certificate.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/formal_core.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order5_m100_tail_curvature_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order5_m100_tail_curvature_reduction.py"
        ),
        "source_contract": sources,
        "exact": exact,
        "conclusions": theorem,
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": sum(
                row.readiness == "ready_to_apply" for row in rows
            ),
            "conditional_rows": sum(
                row.readiness == "conditional_on_open_input" for row in rows
            ),
            "open_rows": sum(
                row.readiness == "not_ready_to_apply" for row in rows
            ),
            "exact_identity_rows": 2,
            "rational_comparison_rows": 1,
            "conditional_tail_theorems": 1,
            "open_curvature_targets": 1,
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    theorem = artifact["conclusions"]
    lines = [
        "# Jensen-Window PF Compound Order-Five Lambda=-100 Tail Curvature Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact order-five `lambda=-100` tail reduction to one open",
        "stable log-curvature ceiling. This is not a proof of the analytic tail,",
        "all-shift order five, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order5_m100_tail_curvature_reduction`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order5_m100_tail_curvature_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order5_m100_tail_curvature_reduction.py",
        "```",
        "",
        "## Exact Scalar Target",
        "",
        "Put `k=n+4` and define",
        "",
        "```text",
        exact["curvature_quantity"],
        "```",
        "",
        "The stable factorization gives the exact equivalence",
        "",
        "```text",
        exact["sign_equivalence"],
        "```",
        "",
        "and the finite-difference form is",
        "",
        "```text",
        exact["curvature_identity"],
        "```",
        "",
        "Thus the raw nested determinant has been reduced to one difference of",
        "stable log-curvatures.",
        "",
        "## Coarse Sufficient Ceiling",
        "",
        "The proved lambda=-100 defect anchor gives",
        "",
        "```text",
        exact["log_buffer"],
        "```",
        "",
        "A very loose sufficient curvature theorem is",
        "",
        "```text",
        exact["sufficient_ceiling"],
        "```",
        "",
        "because exact rational arithmetic gives",
        "",
        "```text",
        exact["rational_comparison"],
        f"after k=321+m: {exact['shifted_polynomial_k_321_plus_n']}>0.",
        "```",
        "",
        "Every coefficient of the shifted polynomial is positive. Therefore the",
        "curvature ceiling implies `C_n<-4log(x_k)`, hence `J_n>0`, on the",
        "whole tail.",
        "",
        "## Boundary Calibration",
        "",
        "The Arb prefix gives at its final row",
        "",
        "```text",
        theorem["boundary_scout"],
        f"relative_316={theorem['prefix_relative']}.",
        "```",
        "",
        "The proposed constant `100` has a factor-above-27 reserve over the",
        "observed scaled curvature at the splice; it is chosen for analytic",
        "robustness, not numerical sharpness.",
        "",
        "## Conditional Completion",
        "",
        "```text",
        theorem["conditional_tail"],
        theorem["conditional_entry"],
        theorem["conditional_forward"],
        "```",
        "",
        "The sole new analytic target in this reduction is",
        "",
        "```text",
        theorem["open_target"],
        "```",
        "",
        "No finite data or eventual non-effective asymptotic is promoted into",
        "that ceiling.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order5_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order5_uniform_tail_flow_reduction.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
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
    with args.out.open("w", encoding="utf-8") as handle:
        json.dump(artifact, handle, indent=2, sort_keys=True)
        handle.write("\n")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-five tail curvature reduction: "
        f"{summary['rows']} rows, "
        f"{summary['conditional_tail_theorems']} conditional tail theorem, "
        f"{summary['open_curvature_targets']} open curvature target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
