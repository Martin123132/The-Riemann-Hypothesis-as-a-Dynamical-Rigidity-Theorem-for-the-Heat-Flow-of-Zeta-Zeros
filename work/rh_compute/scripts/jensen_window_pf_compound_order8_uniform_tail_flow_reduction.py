#!/usr/bin/env python3
"""Prove the signed order-eight tail and reduce propagation to endpoint entry."""

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
    "jensen_window_pf_compound_order8_uniform_tail_flow_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.md"
)
ORDER7_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.json"
)
VANDERMONDE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.json"
)
ORDER = 8
LEADING_DEGREE = ORDER * (ORDER - 1) // 2
RAW_LEADING_CONSTANT = 33664847019245568000


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
    order7 = load_json(ORDER7_SOURCE)
    vandermonde = load_json(VANDERMONDE_SOURCE)
    order7_summary = order7.get("summary", {})
    if order7_summary.get("contiguous_order_seven_theorems") != 1:
        raise RuntimeError("strict order-seven source is not closed")
    if order7_summary.get("arbitrary_column_order_seven_theorems") != 1:
        raise RuntimeError("arbitrary-column order-seven source is not closed")
    vandermonde_summary = vandermonde.get("summary", {})
    if vandermonde_summary.get("all_fixed_order_eventual_tail_theorems") != 1:
        raise RuntimeError("all-fixed-order tail source is not closed")
    order_rows = vandermonde.get("order_rows", [])
    order8 = next((row for row in order_rows if row.get("order") == ORDER), None)
    if order8 is None or order8.get("raw_leading_constant") != RAW_LEADING_CONSTANT:
        raise RuntimeError("order-eight Vandermonde specialization changed")
    if order8.get("epsilon") != 1 or order8.get("leading_degree") != LEADING_DEGREE:
        raise RuntimeError("order-eight signed orientation changed")
    return {
        "order7_status": order7.get("status"),
        "order7_sign": (
            "Q_(7,n)(lambda)=-H_(7,n)(lambda)>0 for every n>=0 and "
            "-100<=lambda<=0"
        ),
        "vandermonde_status": vandermonde.get("status"),
        "order8_leading_degree": order8.get("leading_degree"),
        "order8_raw_leading_constant": order8.get("raw_leading_constant"),
        "order8_positive_signed_constant": order8.get("positive_constant"),
    }


def determinant(values: tuple[sp.Rational, ...], order: int, shift: int) -> sp.Rational:
    return sp.factor(
        sp.Matrix(
            [[values[shift + i + j] for j in range(order)] for i in range(order)]
        ).det(method="domain-ge")
    )


def lower_layer_countermodel() -> dict:
    values = tuple(
        [sp.Rational(1, sp.factorial(index)) for index in range(14)]
        + [sp.Rational(1, 87120000000)]
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
    q7 = [-determinant(values, 7, shift) for shift in range(3)]
    q7_margin = sp.factor(q7[1] ** 2 - q7[0] * q7[2])
    q6_n2 = -determinant(values, 6, 2)
    q8 = determinant(values, 8, 0)
    condensation_residual = sp.factor(q7_margin - q8 * q6_n2)
    if q7_margin >= 0 or q8 >= 0 or condensation_residual != 0:
        raise RuntimeError("countermodel did not violate signed order eight")
    return {
        "sequence": [str(value) for value in values],
        "strict_signed_lower_layers": signed_rows,
        "Q7_values": [str(value) for value in q7],
        "Q7_log_concavity_margin": str(q7_margin),
        "Q6_n2": str(q6_n2),
        "Q8_n0": str(q8),
        "condensation_residual": str(condensation_residual),
        "conclusion": (
            "strict signed contiguous layers through order seven do not imply order eight"
        ),
    }


def exact_statements() -> dict:
    return {
        "signed_coordinate": (
            "Q_(m,n)=epsilon_m*H_(m,n), epsilon_m=(-1)^binom(m,2)"
        ),
        "order8_orientation": "epsilon_8=1, Q_(8,n)=H_(8,n)",
        "leading_term": (
            "det K=33664847019245568000*G_2^28*h^28+o(h^28)"
        ),
        "signed_leading_term": (
            "Q_8=positive_scale*(33664847019245568000*G_2^28*h^28+o(h^28))"
        ),
        "uniform_eventual_tail": (
            "there exists N_8 such that Q_(8,n)(lambda)>0 for every n>=N_8 "
            "and -100<=lambda<=0"
        ),
        "signed_condensation": (
            "Q_(8,n)*Q_(6,n+2)=Q_(7,n+1)^2-"
            "Q_(7,n)*Q_(7,n+2)"
        ),
        "entry_equivalence": (
            "Q_(8,n)>0 iff Q_(7,n+1)^2>Q_(7,n)*Q_(7,n+2)"
        ),
        "affine_derivative": "Q_(8,n)'=(4*n+58)*delta(Q_(8,n))",
        "cooperative_flow": (
            "Q_n'=a_n*Q_(n+1)+b_n*Q_n, "
            "a_n=(4*n+58)*Q_(7,n)/Q_(7,n+1)>0, "
            "b_n=((4*n+58)/(4*n+54))*(log Q_(7,n+1))'"
        ),
        "variation_of_constants": (
            "Q_n(lambda)=E_n(lambda)*(Q_n(-100)+integral_(-100)^lambda "
            "E_n(s)^(-1)*a_n(s)*Q_(n+1)(s)ds)"
        ),
        "conditional_forward": (
            "[Q_(8,n)(-100)>0 for all n] => [Q_(8,n)(lambda)>0 for all n "
            "and -100<=lambda<=0]"
        ),
        "conditional_arbitrary_columns": (
            "completed signed contiguous layers through order eight imply every "
            "arbitrary-column signed layer through order eight"
        ),
        "open_entry": (
            "Q_(8,n)(-100)>0 for every n>=0, equivalently strict "
            "log-concavity of Q_(7,n)(-100)"
        ),
    }


def build_artifact() -> dict:
    sources = validate_sources()
    countermodel = lower_layer_countermodel()
    exact = exact_statements()
    rows = [
        ReductionRow(
            "co8utfr_01_order7_input",
            "theorem_input",
            "ready_to_apply",
            "The completed signed layer through fixed order seven supplies the positive lower cone for order eight.",
            sources["order7_sign"],
            "Completed fixed-order input only; no order-eight sign is imported.",
        ),
        ReductionRow(
            "co8utfr_02_eventual_tail",
            "exact_published_composition",
            "ready_to_apply",
            "The all-fixed-order Vandermonde theorem supplies a uniform eventual signed order-eight tail.",
            exact["signed_leading_term"] + "; " + exact["uniform_eventual_tail"],
            "Non-effective tail; the threshold may depend on order eight.",
        ),
        ReductionRow(
            "co8utfr_03_entry_coordinate",
            "exact_identity",
            "ready_to_apply",
            "Signed order-eight entry is strict log-concavity of the completed positive Q_7 sequence.",
            exact["signed_condensation"] + "; " + exact["entry_equivalence"],
            "Exact signed Desnanot-Jacobi specialization only.",
        ),
        ReductionRow(
            "co8utfr_04_affine_heat",
            "exact_identity",
            "ready_to_apply",
            "The affine Hankel heat derivative specializes with coefficient 4n+58.",
            exact["affine_derivative"],
            "General determinant identity specialized at m=8.",
        ),
        ReductionRow(
            "co8utfr_05_cooperative_flow",
            "exact_forward_lemma",
            "ready_to_apply",
            "The signed order-eight heat flow is cooperative throughout the completed order-seven cone.",
            exact["cooperative_flow"],
            "The off-diagonal coefficient uses only Q_7>0.",
        ),
        ReductionRow(
            "co8utfr_06_conditional_forward",
            "exact_conditional_theorem",
            "conditional_on_open_input",
            "All-shift lambda=-100 entry would complete signed contiguous order eight through lambda zero.",
            exact["conditional_forward"] + "; " + exact["variation_of_constants"],
            "Conditional only on the explicit endpoint theorem.",
        ),
        ReductionRow(
            "co8utfr_07_conditional_columns",
            "exact_conditional_theorem",
            "conditional_on_open_input",
            "The fixed-order initial-minor theorem would then transfer order eight to arbitrary columns.",
            exact["conditional_arbitrary_columns"],
            "Conditional on contiguous order-eight completion.",
        ),
        ReductionRow(
            "co8utfr_08_countermodel",
            "countermodel_gate",
            "ready_to_apply",
            "The completed lower signed layers do not imply order eight abstractly.",
            "Q_(8,0)=" + countermodel["Q8_n0"] + "<0",
            "Exact rational finite sequence; blocks lower-layer promotion.",
            countermodel,
        ),
        ReductionRow(
            "co8utfr_09_open_entry",
            "analytic_theorem_target",
            "not_ready_to_apply",
            "Prove strict all-shift signed order-eight entry at lambda=-100.",
            exact["open_entry"],
            "Sole open input in this fixed-order-eight propagation route.",
        ),
    ]
    source_paths = (ORDER7_SOURCE, VANDERMONDE_SOURCE)
    return {
        "kind": "jensen_window_pf_compound_order8_uniform_tail_flow_reduction",
        "date": "2026-07-13",
        "status": (
            "exact uniform eventual signed order-eight tail and conditional "
            "forward reduction with one open lambda=-100 entry"
        ),
        "proof_boundary": (
            "This artifact proves the signed order-eight eventual tail, exact "
            "entry coordinate, cooperative flow, and conditional forward theorem. "
            "It does not prove endpoint entry, full order-eight invariance, "
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
            "jensen_window_pf_compound_order8_uniform_tail_flow_reduction.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order8_uniform_tail_flow_reduction.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    countermodel = artifact["countermodel"]
    lines = [
        "# Jensen-Window PF Compound Order-Eight Uniform Tail And Flow Reduction",
        "",
        "Date: 2026-07-13",
        "",
        "Status: exact uniform eventual signed order-eight tail and conditional",
        "forward reduction with one open `lambda=-100` entry. This is not a",
        "proof of all-shift order eight, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order8_uniform_tail_flow_reduction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order8_uniform_tail_flow_reduction.py",
        "```",
        "",
        "## Universal Signed Tail",
        "",
        "The graded-kernel Vandermonde lemma applies at `m=8`, where",
        "`D=binom(8,2)=28` and `epsilon_8=1`. It gives",
        "",
        "```text",
        exact["leading_term"],
        exact["signed_leading_term"],
        exact["uniform_eventual_tail"],
        "```",
        "",
        "The threshold is finite but non-effective. This is an eventual tail,",
        "not an all-shift order-eight theorem.",
        "",
        "## Endpoint Coordinate",
        "",
        "Signed Desnanot-Jacobi gives",
        "",
        "```text",
        exact["signed_condensation"],
        exact["entry_equivalence"],
        "```",
        "",
        "The completed positive `Q_6=-H_6` and `Q_7=-H_7` layers make this",
        "an exact strict log-concavity target for the endpoint `Q_7` sequence.",
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
        "Since `Q_(7,n)>0` on the complete heat interval, `a_n>0`. The uniform",
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
        "has every available signed contiguous minor through order seven strict",
        "and positive, but",
        "",
        "```text",
        "Q_(7,1)^2-Q_(7,0)Q_(7,2)="
        + countermodel["Q7_log_concavity_margin"]
        + "<0,",
        "Q_(8,0)=" + countermodel["Q8_n0"] + "<0.",
        "```",
        "",
        "Thus no lower-cone promotion can replace the new endpoint theorem.",
        "",
        "## Remaining Endpoint Target",
        "",
        "The sole fixed-order-eight input still missing is",
        "",
        "```text",
        exact["open_entry"],
        "```",
        "",
        "The next decision is whether to build another zeta-specific stable",
        "curvature ladder or replace the layer-by-layer endpoint analysis with",
        "a uniform-in-order theorem.",
        "",
        "```text",
        "outputs/jensen_window_pf_graded_kernel_vandermonde_all_order_lemma.md",
        "outputs/jensen_window_pf_compound_order7_uniform_heat_forward_invariance_certificate.md",
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
        "wrote order-eight uniform tail and flow reduction: "
        f"{summary['rows']} rows, "
        f"{summary['universal_tail_specializations']} tail specialization, "
        f"{summary['open_entry_targets']} open entry"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
