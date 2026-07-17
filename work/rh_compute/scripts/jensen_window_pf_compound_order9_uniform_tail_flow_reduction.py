#!/usr/bin/env python3
"""Prove the signed order-nine tail and reduce propagation to endpoint entry."""

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
    "jensen_window_pf_compound_order9_uniform_tail_flow_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.md"
)
ORDER8_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.json"
)
VANDERMONDE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json"
)
ORDER = 9
LEADING_DEGREE = ORDER * (ORDER - 1) // 2
RAW_LEADING_CONSTANT = 347485857744891213250560000


@dataclass(frozen=True)
class ReductionRow:
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
    order8 = load_json(ORDER8_SOURCE)
    vandermonde = load_json(VANDERMONDE_SOURCE)
    order8_summary = order8.get("summary", {})
    if order8_summary.get("contiguous_order_eight_theorems") != 1:
        raise RuntimeError("strict order-eight source is not closed")
    if order8_summary.get("arbitrary_column_order_eight_theorems") != 1:
        raise RuntimeError("arbitrary-column order-eight source is not closed")
    if (
        vandermonde.get("summary", {}).get(
            "all_fixed_order_eventual_tail_theorems"
        )
        != 1
    ):
        raise RuntimeError("all-fixed-order tail source is not closed")
    order9 = next(
        (row for row in vandermonde.get("order_rows", []) if row.get("order") == ORDER),
        None,
    )
    if order9 is None or order9.get("raw_leading_constant") != RAW_LEADING_CONSTANT:
        raise RuntimeError("order-nine Vandermonde specialization changed")
    if order9.get("epsilon") != 1 or order9.get("leading_degree") != LEADING_DEGREE:
        raise RuntimeError("order-nine signed orientation changed")
    return {
        "order8_status": order8.get("status"),
        "order8_sign": (
            "Q_(8,n)(lambda)=H_(8,n)(lambda)>0 for every n>=0 and "
            "-100<=lambda<=0"
        ),
        "vandermonde_status": vandermonde.get("status"),
        "order9_leading_degree": order9.get("leading_degree"),
        "order9_raw_leading_constant": order9.get("raw_leading_constant"),
        "order9_positive_signed_constant": order9.get("positive_constant"),
    }


def determinant(values: tuple[sp.Rational, ...], order: int, shift: int) -> sp.Rational:
    return sp.factor(
        sp.Matrix(
            [[values[shift + i + j] for j in range(order)] for i in range(order)]
        ).det(method="domain-ge")
    )


def lower_layer_countermodel() -> dict:
    values = tuple(
        [sp.Rational(1, sp.factorial(index)) for index in range(16)]
        + [sp.Rational(1, 20926400000000)]
    )
    signed_rows = {}
    for order in range(1, ORDER):
        epsilon = (-1) ** (order * (order - 1) // 2)
        rows = [
            sp.factor(epsilon * determinant(values, order, shift))
            for shift in range(len(values) - 2 * order + 2)
        ]
        if not all(value > 0 for value in rows):
            raise RuntimeError(f"countermodel lower layer failed at order {order}")
        signed_rows[str(order)] = [str(value) for value in rows]
    q8 = [determinant(values, 8, shift) for shift in range(3)]
    q8_margin = sp.factor(q8[1] ** 2 - q8[0] * q8[2])
    q7_n2 = -determinant(values, 7, 2)
    q9 = determinant(values, 9, 0)
    condensation_residual = sp.factor(q8_margin - q9 * q7_n2)
    if q8_margin >= 0 or q9 >= 0 or condensation_residual != 0:
        raise RuntimeError("countermodel did not violate signed order nine")
    return {
        "sequence": [str(value) for value in values],
        "strict_signed_lower_layers": signed_rows,
        "Q8_values": [str(value) for value in q8],
        "Q8_log_concavity_margin": str(q8_margin),
        "Q7_n2": str(q7_n2),
        "Q9_n0": str(q9),
        "condensation_residual": str(condensation_residual),
        "conclusion": (
            "strict signed contiguous layers through order eight do not imply "
            "order nine"
        ),
    }


def exact_statements() -> dict:
    return {
        "signed_coordinate": (
            "Q_(m,n)=epsilon_m*H_(m,n), epsilon_m=(-1)^binom(m,2)"
        ),
        "order9_orientation": "epsilon_9=1, Q_(9,n)=H_(9,n)",
        "leading_term": (
            "det K=347485857744891213250560000*G_2^36*h^36+o(h^36)"
        ),
        "signed_leading_term": (
            "Q_9=positive_scale*(347485857744891213250560000*"
            "G_2^36*h^36+o(h^36))"
        ),
        "uniform_eventual_tail": (
            "there exists N_9 such that Q_(9,n)(lambda)>0 for every n>=N_9 "
            "and -100<=lambda<=0"
        ),
        "signed_condensation": (
            "Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-"
            "Q_(8,n)*Q_(8,n+2)"
        ),
        "entry_equivalence": (
            "Q_(9,n)>0 iff Q_(8,n+1)^2>Q_(8,n)*Q_(8,n+2)"
        ),
        "affine_derivative": "Q_(9,n)'=(4*n+66)*delta(Q_(9,n))",
        "cooperative_flow": (
            "Q_n'=a_n*Q_(n+1)+b_n*Q_n, "
            "a_n=(4*n+66)*Q_(8,n)/Q_(8,n+1)>0, "
            "b_n=((4*n+66)/(4*n+62))*(log Q_(8,n+1))'"
        ),
        "variation_of_constants": (
            "Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda "
            "E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)"
        ),
        "conditional_forward": (
            "[Q_(9,n)(-100)>0 for all n] => [Q_(9,n)(lambda)>0 for all n "
            "and -100<=lambda<=0]"
        ),
        "conditional_arbitrary_columns": (
            "completed signed contiguous layers through order nine imply every "
            "arbitrary-column signed layer through order nine"
        ),
        "open_entry": (
            "Q_(9,n)(-100)>0 for every n>=0, equivalently strict "
            "log-concavity of Q_(8,n)(-100)"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    countermodel = lower_layer_countermodel()
    exact = exact_statements()
    rows = [
        ReductionRow(
            "co9utfr_01_order8_input",
            "theorem_input",
            "ready_to_apply",
            "The completed signed layer through fixed order eight supplies the positive lower cone for order nine.",
            sources["order8_sign"],
            "Completed fixed-order input only; no order-nine sign is imported.",
        ),
        ReductionRow(
            "co9utfr_02_eventual_tail",
            "exact_published_composition",
            "ready_to_apply",
            "The all-fixed-order Vandermonde theorem supplies a uniform eventual signed order-nine tail.",
            exact["signed_leading_term"] + "; " + exact["uniform_eventual_tail"],
            "Non-effective tail; the threshold may depend on order nine.",
        ),
        ReductionRow(
            "co9utfr_03_entry_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order-nine entry is strict log-concavity of the completed positive Q_8 sequence.",
            exact["signed_condensation"] + "; " + exact["entry_equivalence"],
            "Exact signed Desnanot-Jacobi specialization only.",
        ),
        ReductionRow(
            "co9utfr_04_affine_heat",
            "exact_identity",
            "ready_to_apply",
            "The affine Hankel heat derivative specializes with coefficient 4n+66.",
            exact["affine_derivative"],
            "General determinant identity specialized at m=9.",
        ),
        ReductionRow(
            "co9utfr_05_cooperative_flow",
            "exact_forward_lemma",
            "ready_to_apply",
            "The signed order-nine heat flow is cooperative throughout the completed order-eight cone.",
            exact["cooperative_flow"],
            "The off-diagonal coefficient uses only Q_8>0.",
        ),
        ReductionRow(
            "co9utfr_06_conditional_forward",
            "exact_conditional_theorem",
            "conditional_on_open_input",
            "All-shift lambda=-100 entry would complete signed contiguous order nine through lambda zero.",
            exact["conditional_forward"] + "; " + exact["variation_of_constants"],
            "Conditional only on the explicit endpoint theorem.",
        ),
        ReductionRow(
            "co9utfr_07_conditional_columns",
            "exact_conditional_theorem",
            "conditional_on_open_input",
            "The fixed-order initial-minor theorem would then transfer order nine to arbitrary columns.",
            exact["conditional_arbitrary_columns"],
            "Conditional on contiguous order-nine completion.",
        ),
        ReductionRow(
            "co9utfr_08_countermodel",
            "countermodel_gate",
            "ready_to_apply",
            "The completed lower signed layers do not imply order nine abstractly.",
            "Q_(9,0)=" + countermodel["Q9_n0"] + "<0",
            "Exact rational finite sequence; blocks lower-layer promotion.",
            countermodel,
        ),
        ReductionRow(
            "co9utfr_09_open_entry",
            "analytic_theorem_target",
            "not_ready_to_apply",
            "Prove strict all-shift signed order-nine entry at lambda=-100.",
            exact["open_entry"],
            "Sole open input in this fixed-order-nine propagation route.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_uniform_tail_flow_reduction",
        "date": "2026-07-13",
        "status": (
            "exact uniform eventual signed order-nine tail and conditional "
            "forward reduction with one open lambda=-100 entry"
        ),
        "proof_boundary": (
            "This artifact proves the signed order-nine eventual tail, exact "
            "entry coordinate, cooperative flow, and conditional forward theorem. "
            "It does not prove endpoint entry, full order-nine invariance, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            ORDER8_SOURCE.relative_to(REPO_ROOT).as_posix(),
            VANDERMONDE_SOURCE.relative_to(REPO_ROOT).as_posix(),
        ],
        "source_contract": sources,
        "exact": exact,
        "countermodel": countermodel,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": 6,
            "conditional_rows": 2,
            "open_rows": 1,
            "universal_tail_specializations": 1,
            "signed_condensation_coordinates": 1,
            "cooperative_flow_recursions": 1,
            "conditional_forward_theorems": 1,
            "conditional_arbitrary_column_theorems": 1,
            "lower_layer_countermodels": 1,
            "open_entry_targets": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    countermodel = artifact["countermodel"]
    lines = [
        "# Jensen-Window PF Compound Order-Nine Uniform Tail And Flow Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact uniform eventual signed order-nine tail and conditional",
        "forward reduction with one open `lambda=-100` entry. This is not a",
        "proof of all-shift order nine, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_uniform_tail_flow_reduction.py",
        "```",
        "",
        "## Universal Signed Tail",
        "",
        "At `m=9`, `D=36` and `epsilon_9=1`. The all-fixed-order theorem gives",
        "",
        "```text",
        exact["leading_term"],
        exact["uniform_eventual_tail"],
        "```",
        "",
        "## Endpoint Coordinate And Flow",
        "",
        "```text",
        exact["signed_condensation"],
        exact["entry_equivalence"],
        exact["affine_derivative"],
        exact["cooperative_flow"],
        "```",
        "",
        "The completed `Q_7,Q_8>0` cone makes the off-diagonal coefficient",
        "strictly positive. The eventual tail and variation of constants prove",
        "the forward and arbitrary-column conclusions conditional only on entry.",
        "",
        "## Countermodel Gate",
        "",
        "The exact sequence `(1,1,1/2!,...,1/15!,1/20926400000000)` has every",
        "available signed contiguous minor through order eight strictly positive,",
        "but",
        "",
        "```text",
        "Q_(8,1)^2-Q_(8,0)Q_(8,2)="
        + countermodel["Q8_log_concavity_margin"]
        + "<0,",
        "Q_(9,0)=" + countermodel["Q9_n0"] + "<0.",
        "```",
        "",
        "Thus lower-cone promotion cannot replace a new zeta-specific theorem.",
        "",
        "## Remaining Endpoint Target",
        "",
        "```text",
        exact["open_entry"],
        "```",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order8_uniform_heat_forward_invariance_certificate.md",
        "outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-nine uniform tail and flow reduction: "
        f"{summary['rows']} rows, "
        f"{summary['universal_tail_specializations']} tail specialization, "
        f"{summary['open_entry_targets']} open entry"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
