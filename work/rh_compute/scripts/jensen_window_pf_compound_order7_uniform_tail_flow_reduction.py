#!/usr/bin/env python3
"""Prove the signed order-seven tail and reduce propagation to endpoint entry."""

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
    "jensen_window_pf_compound_order7_uniform_tail_flow_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.md"
)
ORDER6_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.json"
)
VANDERMONDE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json"
)
ORDER = 7
LEADING_DEGREE = ORDER * (ORDER - 1) // 2
RAW_LEADING_CONSTANT = -52183852646400


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
    order6 = load_json(ORDER6_SOURCE)
    vandermonde = load_json(VANDERMONDE_SOURCE)
    order6_summary = order6.get("summary", {})
    if order6_summary.get("contiguous_order_six_theorems") != 1:
        raise RuntimeError("strict order-six source is not closed")
    if order6_summary.get("arbitrary_column_order_six_theorems") != 1:
        raise RuntimeError("arbitrary-column order-six source is not closed")
    vandermonde_summary = vandermonde.get("summary", {})
    if vandermonde_summary.get("all_fixed_order_eventual_tail_theorems") != 1:
        raise RuntimeError("all-fixed-order tail source is not closed")
    order_rows = vandermonde.get("order_rows", [])
    order7 = next((row for row in order_rows if row.get("order") == ORDER), None)
    if order7 is None or order7.get("raw_leading_constant") != RAW_LEADING_CONSTANT:
        raise RuntimeError("order-seven Vandermonde specialization changed")
    return {
        "order6_status": order6.get("status"),
        "order6_sign": (
            "Q_(6,n)(lambda)=-H_(6,n)(lambda)>0 for every n>=0 and "
            "-100<=lambda<=0"
        ),
        "vandermonde_status": vandermonde.get("status"),
        "order7_leading_degree": order7.get("leading_degree"),
        "order7_raw_leading_constant": order7.get("raw_leading_constant"),
        "order7_positive_signed_constant": order7.get("positive_constant"),
    }


def determinant(values: tuple[sp.Rational, ...], order: int, shift: int) -> sp.Rational:
    return sp.factor(
        sp.Matrix(
            [[values[shift + i + j] for j in range(order)] for i in range(order)]
        ).det(method="domain-ge")
    )


def lower_layer_countermodel() -> dict:
    values = tuple(
        [sp.Rational(1, sp.factorial(index)) for index in range(12)]
        + [sp.Rational(1, 480100000)]
    )
    signed_rows = {}
    for order in range(1, ORDER):
        epsilon = (-1) ** (order * (order - 1) // 2)
        count = len(values) - 2 * order + 2
        rows = [
            sp.factor(epsilon * determinant(values, order, shift))
            for shift in range(count)
        ]
        if not all(value > 0 for value in rows):
            raise RuntimeError(f"countermodel lower layer failed at order {order}")
        signed_rows[str(order)] = [str(value) for value in rows]
    q6 = [-determinant(values, 6, shift) for shift in range(3)]
    q6_margin = sp.factor(q6[1] ** 2 - q6[0] * q6[2])
    h7 = determinant(values, 7, 0)
    q7 = -h7
    if q6_margin >= 0 or q7 >= 0 or q6_margin == 0 or q7 == 0:
        raise RuntimeError("countermodel did not violate signed order seven")
    return {
        "sequence": [str(value) for value in values],
        "strict_signed_lower_layers": signed_rows,
        "Q6_values": [str(value) for value in q6],
        "Q6_log_concavity_margin": str(q6_margin),
        "H7_n0": str(h7),
        "Q7_n0": str(q7),
        "conclusion": (
            "strict signed contiguous layers through order six do not imply order seven"
        ),
    }


def exact_statements() -> dict:
    return {
        "signed_coordinate": (
            "Q_(m,n)=epsilon_m*H_(m,n), epsilon_m=(-1)^binom(m,2)"
        ),
        "order7_orientation": "epsilon_7=-1, Q_(7,n)=-H_(7,n)",
        "leading_term": (
            "det K=-52183852646400*G_2^21*h^21+o(h^21)"
        ),
        "signed_leading_term": (
            "Q_7=positive_scale*(52183852646400*G_2^21*h^21+o(h^21))"
        ),
        "uniform_eventual_tail": (
            "there exists N_7 such that Q_(7,n)(lambda)>0 for every n>=N_7 "
            "and -100<=lambda<=0"
        ),
        "signed_condensation": (
            "Q_(7,n)*Q_(5,n+2)=Q_(6,n+1)^2-"
            "Q_(6,n)*Q_(6,n+2)"
        ),
        "order7_condensation": (
            "Q_(7,n)*H_(5,n+2)=Q_(6,n+1)^2-"
            "Q_(6,n)*Q_(6,n+2)"
        ),
        "entry_equivalence": (
            "Q_(7,n)>0 iff Q_(6,n+1)^2>Q_(6,n)*Q_(6,n+2)"
        ),
        "affine_derivative": "Q_(7,n)'=(4*n+50)*delta(Q_(7,n))",
        "cooperative_flow": (
            "Q_n'=a_n*Q_(n+1)+b_n*Q_n, "
            "a_n=(4*n+50)*Q_(6,n)/Q_(6,n+1)>0, "
            "b_n=((4*n+50)/(4*n+46))*(log Q_(6,n+1))'"
        ),
        "variation_of_constants": (
            "Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda "
            "E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)"
        ),
        "conditional_forward": (
            "[Q_(7,n)(-100)>0 for all n] => [Q_(7,n)(lambda)>0 for all n "
            "and -100<=lambda<=0]"
        ),
        "conditional_arbitrary_columns": (
            "completed signed contiguous layers through order seven imply every "
            "arbitrary-column signed layer through order seven"
        ),
        "open_entry": (
            "Q_(7,n)(-100)>0 for every n>=0, equivalently strict "
            "log-concavity of Q_(6,n)(-100)"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    countermodel = lower_layer_countermodel()
    exact = exact_statements()
    rows = [
        ReductionRow(
            "co7utfr_01_order6_input",
            "theorem_input",
            "ready_to_apply",
            "The completed signed layer through fixed order six supplies the positive lower cone for order seven.",
            sources["order6_sign"],
            "Completed fixed-order input only; no order-seven sign is imported.",
        ),
        ReductionRow(
            "co7utfr_02_eventual_tail",
            "exact_published_composition",
            "ready_to_apply",
            "The all-fixed-order Vandermonde theorem supplies a uniform eventual signed order-seven tail.",
            exact["signed_leading_term"] + "; " + exact["uniform_eventual_tail"],
            "Non-effective tail; the threshold may depend on order seven.",
        ),
        ReductionRow(
            "co7utfr_03_entry_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order-seven entry is strict log-concavity of the completed positive Q_6 sequence.",
            exact["order7_condensation"] + "; " + exact["entry_equivalence"],
            "Exact Desnanot-Jacobi specialization only.",
        ),
        ReductionRow(
            "co7utfr_04_affine_heat",
            "exact_identity",
            "ready_to_apply",
            "The affine Hankel heat derivative specializes with coefficient 4n+50.",
            exact["affine_derivative"],
            "General determinant identity specialized at m=7.",
        ),
        ReductionRow(
            "co7utfr_05_cooperative_flow",
            "exact_forward_lemma",
            "ready_to_apply",
            "The signed order-seven heat flow is cooperative throughout the completed order-six cone.",
            exact["cooperative_flow"],
            "The off-diagonal coefficient uses only Q_6>0.",
        ),
        ReductionRow(
            "co7utfr_06_conditional_forward",
            "exact_conditional_theorem",
            "conditional_on_open_input",
            "All-shift lambda=-100 entry would complete signed contiguous order seven through lambda zero.",
            exact["conditional_forward"] + "; " + exact["variation_of_constants"],
            "Conditional only on the explicit endpoint theorem.",
        ),
        ReductionRow(
            "co7utfr_07_conditional_columns",
            "exact_conditional_theorem",
            "conditional_on_open_input",
            "The fixed-order initial-minor theorem would then transfer order seven to arbitrary columns.",
            exact["conditional_arbitrary_columns"],
            "Conditional on contiguous order-seven completion.",
        ),
        ReductionRow(
            "co7utfr_08_countermodel",
            "countermodel_gate",
            "ready_to_apply",
            "The completed lower signed layers do not imply order seven abstractly.",
            "Q_(7,0)=" + countermodel["Q7_n0"] + "<0",
            "Exact rational finite sequence; blocks lower-layer promotion.",
            countermodel,
        ),
        ReductionRow(
            "co7utfr_09_open_entry",
            "analytic_theorem_target",
            "not_ready_to_apply",
            "Prove strict all-shift signed order-seven entry at lambda=-100.",
            exact["open_entry"],
            "Sole open input in this fixed-order-seven propagation route.",
        ),
    ]
    source_paths = (ORDER6_SOURCE, VANDERMONDE_SOURCE)
    return {
        "kind": "jensen_window_pf_compound_order7_uniform_tail_flow_reduction",
        "date": "2026-07-13",
        "status": (
            "exact uniform eventual signed order-seven tail and conditional "
            "forward reduction with one open lambda=-100 entry"
        ),
        "proof_boundary": (
            "This artifact proves the signed order-seven eventual tail, exact "
            "entry coordinate, cooperative flow, and conditional forward theorem. "
            "It does not prove endpoint entry, full order-seven invariance, "
            "PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [path.relative_to(REPO_ROOT).as_posix() for path in source_paths],
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
            "jensen_window_pf_compound_order7_uniform_tail_flow_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order7_uniform_tail_flow_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    countermodel = artifact["countermodel"]
    lines = [
        "# Jensen-Window PF Compound Order-Seven Uniform Tail And Flow Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact uniform eventual signed order-seven tail and conditional",
        "forward reduction with one open `lambda=-100` entry. This is not a",
        "proof of all-shift order seven, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order7_uniform_tail_flow_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order7_uniform_tail_flow_reduction.py",
        "```",
        "",
        "## Universal Signed Tail",
        "",
        "The graded-kernel Vandermonde lemma applies at `m=7`, where",
        "`D=binom(7,2)=21` and `epsilon_7=-1`. It gives",
        "",
        "```text",
        exact["leading_term"],
        exact["signed_leading_term"],
        exact["uniform_eventual_tail"],
        "```",
        "",
        "The threshold is finite but non-effective. This is an eventual tail,",
        "not an all-shift order-seven theorem.",
        "",
        "## Endpoint Coordinate",
        "",
        "Signed Desnanot-Jacobi gives",
        "",
        "```text",
        exact["signed_condensation"],
        exact["order7_condensation"],
        exact["entry_equivalence"],
        "```",
        "",
        "The completed positive `Q_5=H_5` and `Q_6=-H_6` layers make this an",
        "exact strict log-concavity target for the endpoint `Q_6` sequence.",
        "",
        "## Cooperative Heat Flow",
        "",
        "The general affine-Hankel and adjacent Plucker identities specialize to",
        "",
        "```text",
        exact["affine_derivative"],
        exact["cooperative_flow"],
        "```",
        "",
        "Since `Q_(6,n)>0` on the complete heat interval, `a_n>0`. The uniform",
        "eventual tail and variation of constants therefore prove",
        "",
        "```text",
        exact["conditional_forward"],
        exact["conditional_arbitrary_columns"],
        "```",
        "",
        "conditional only on all-shift endpoint entry.",
        "",
        "## Countermodel Gate",
        "",
        "The exact rational sequence",
        "",
        "```text",
        "(" + ",".join(countermodel["sequence"]) + ")",
        "```",
        "",
        "has every available signed contiguous minor through order six strict",
        "and positive, but",
        "",
        "```text",
        "Q_(6,1)^2-Q_(6,0)Q_(6,2)="
        + countermodel["Q6_log_concavity_margin"]
        + "<0,",
        "Q_(7,0)=" + countermodel["Q7_n0"] + "<0.",
        "```",
        "",
        "Thus no lower-cone promotion can replace the new endpoint theorem.",
        "",
        "## Remaining Endpoint Target",
        "",
        "The sole fixed-order-seven input still missing is",
        "",
        "```text",
        exact["open_entry"],
        "```",
        "",
        "A rigorous coefficient prefix and a cancellation-preserving analytic",
        "tail coordinate are the next two sub-obligations.",
        "",
        "```text",
        "outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md",
        "outputs/jensen_window_pf_compound_order6_uniform_heat_forward_invariance_certificate.md",
        "outputs/signed_hankel_jensen_bridge_target.md",
        "outputs/formal_core.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-seven uniform tail and flow reduction: "
        f"{summary['rows']} rows, "
        f"{summary['universal_tail_specializations']} tail specialization, "
        f"{summary['open_entry_targets']} open entry"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
