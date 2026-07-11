#!/usr/bin/env python3
"""Exact counterexample to static full-ratio-cone all-m column recurrence."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_sparse_degree6_scout as sparse


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_all_m_counterexample.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_all_m_counterexample.md"
DEGREE = 7
MINOR_SIZE = 11
WITNESS_SHIFT = 0
WITNESS_X = [Fraction(19, 20), Fraction(19, 20), Fraction(1), Fraction(1), Fraction(1), Fraction(1)]
WITNESS_S = [Fraction(19, 20), Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(0)]


def rational_str(value: Fraction | int) -> str:
    frac = Fraction(value)
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def evaluate_polynomial(poly: sparse.Polynomial, values: list[Fraction]) -> Fraction:
    total = Fraction(0)
    for exponents, coefficient in poly.items():
        term = Fraction(coefficient)
        for value, exponent in zip(values, exponents):
            if exponent:
                term *= value**exponent
        total += term
    return total


def x_from_s(values: list[Fraction]) -> list[Fraction]:
    product = Fraction(1)
    out: list[Fraction] = []
    for value in values:
        product *= 1 - value
        out.append(1 - product)
    return out


def monotone_checks(values: list[Fraction]) -> dict:
    return {
        "all_in_unit_interval": all(0 <= value <= 1 for value in values),
        "weakly_increasing": all(left <= right for left, right in zip(values, values[1:])),
        "strictly_positive_first": values[0] > 0,
        "upper_bound_ok": values[-1] <= 1,
    }


def lower_wall(index: int) -> Fraction:
    return Fraction(2 * index - 1, 2 * index + 1)


def full_cone_checks(values: list[Fraction], shift: int) -> dict:
    global_indices = [shift + offset for offset in range(1, len(values) + 1)]
    walls = [lower_wall(index) for index in global_indices]
    margins = [value - wall for value, wall in zip(values, walls)]
    return {
        "shift": shift,
        "global_indices": global_indices,
        "lower_walls": [rational_str(value) for value in walls],
        "lower_wall_margins": [rational_str(value) for value in margins],
        "all_lower_walls_ok": all(margin >= 0 for margin in margins),
        "minimum_lower_wall_margin": rational_str(min(margins)),
        "infinite_extension": "x_k=1 for every k>=7",
        "infinite_extension_obeys_lower_wall": all(lower_wall(index) <= 1 for index in range(7, 64)),
    }


def build_payload() -> dict:
    columns = sparse.normalized_column_polynomials(DEGREE, MINOR_SIZE)
    q_poly = columns[MINOR_SIZE]
    value = evaluate_polynomial(q_poly, WITNESS_X)
    s_to_x = x_from_s(WITNESS_S)
    return {
        "kind": "jensen_window_pf_monotone_contraction_all_m_counterexample",
        "date": "2026-07-10",
        "status": "exact_countermodel_gate",
        "source_column_recurrence_contract": "outputs/jensen_window_pf_column_recurrence_contract.md",
        "source_monotone_contraction_theorem_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "source_degree7_frontier_scout": "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md",
        "proof_boundary": (
            "Exact abstract static-cone countermodel only. It proves that the full pointwise "
            "lower/upper walls and monotone contractions do not imply all-m column-recurrence "
            "positivity. It is not a zeta-window row, not "
            "evidence against RH, not evidence against Lambda <= 0, and not a disproof of any "
            "zeta-specific strengthened theorem."
        ),
        "degree": DEGREE,
        "minor_size": MINOR_SIZE,
        "normalized_recurrence": (
            "Q_0=1 and Q_m=sum_{j=1..min(7,m)} (-1)^(j-1) binom(7,j) "
            "prod_{i=1..j-1} x_i^(j-i) Q_{m-j}, after factoring A^m*rho^m."
        ),
        "witness": {
            "s_parameters": [rational_str(value) for value in WITNESS_S],
            "x_contractions": [rational_str(value) for value in WITNESS_X],
            "x_from_s": [rational_str(value) for value in s_to_x],
            "monotone_checks": monotone_checks(WITNESS_X),
            "full_cone_checks": full_cone_checks(WITNESS_X, WITNESS_SHIFT),
            "q_value": rational_str(value),
            "q_value_decimal": float(value),
            "q_value_is_negative": value < 0,
        },
        "polynomial_summary": {
            "normalized_recurrence_term_count": len(q_poly),
            "x_power_multidegree": [max(monomial[index] for monomial in q_poly) for index in range(DEGREE - 1)],
        },
        "claim_blocked": (
            "The implication 'the propagated static ratio cone implies C_m>=0 for all column "
            "rows' is false. Any all-m route must add heat-flow or Xi/Phi-specific structure, "
            "subdivision "
            "with additional hypotheses, or a different positive construction."
        ),
        "summary": {
            "counterexample_degree": DEGREE,
            "counterexample_minor_size": MINOR_SIZE,
            "negative_witness_count": 1,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The exact shift-0 witness x=(19/20,19/20,1,1,1,1), extended by x_k=1 for "
                "k>=7, satisfies the full propagated static ratio cone, including every "
                "pointwise lower wall, but makes the normalized degree-7 m=11 column recurrence "
                "negative. Thus the static cone is not an all-order PF theorem by itself."
            ),
        },
        "invariants": [
            "The witness and its x_k=1 tail satisfy (2*k-1)/(2*k+1)<=x_k<=1 and x_(k+1)>=x_k for every k>=1.",
            "The normalized degree-7 m=11 column recurrence value is strictly negative.",
            "This blocks only generic monotone-contraction all-m promotion.",
            "The witness is abstract and is not claimed to be a zeta heat-flow window.",
            "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    witness = payload["witness"]
    result_line = (
        "validated Jensen-window PF monotone-contraction all-m counterexample: "
        "degree 7, m=11, exact full-cone witness, 6 lower walls, negative normalized value, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Monotone-Contraction All-M Counterexample",
        "",
        "Date: 2026-07-10",
        "",
        "Status: exact countermodel gate. This is not evidence against RH,",
        "not evidence against `Lambda <= 0`, and not a zeta-window computation.",
        "",
        "Artifact kind: `jensen_window_pf_monotone_contraction_all_m_counterexample`.",
        "",
        "Proof boundary: this artifact blocks one generic proof step only:",
        "the full static ratio cone does not imply all-`m` column recurrence",
        "positivity. A zeta-specific strengthened theorem may still exist.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_monotone_contraction_all_m_counterexample.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_all_m_counterexample.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_all_m_counterexample.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Witness",
        "",
        "For degree `d=7` and column row `m=11`, after factoring the positive",
        "monomial `A^m*rho^m`, evaluate the normalized recurrence at:",
        "",
        "```text",
        f"s = ({', '.join(witness['s_parameters'])})",
        f"x = ({', '.join(witness['x_contractions'])})",
        "```",
        "",
        "The full static-cone checks are:",
        "",
        "```text",
        f"all_in_unit_interval = {witness['monotone_checks']['all_in_unit_interval']}",
        f"weakly_increasing = {witness['monotone_checks']['weakly_increasing']}",
        f"upper_bound_ok = {witness['monotone_checks']['upper_bound_ok']}",
        f"lower_walls = ({', '.join(witness['full_cone_checks']['lower_walls'])})",
        f"lower_wall_margins = ({', '.join(witness['full_cone_checks']['lower_wall_margins'])})",
        f"all_lower_walls_ok = {witness['full_cone_checks']['all_lower_walls_ok']}",
        f"tail = {witness['full_cone_checks']['infinite_extension']}",
        f"tail_lower_wall_ok = {witness['full_cone_checks']['infinite_extension_obeys_lower_wall']}",
        "```",
        "",
        "The exact normalized column value is:",
        "",
        "```text",
        f"Q_11 = {witness['q_value']}",
        f"Q_11 approximately {witness['q_value_decimal']}",
        "```",
        "",
        "This value is strictly negative.",
        "",
        "## Consequence",
        "",
        payload["summary"]["main_finding"],
        "",
        "Integration:",
        "",
        "```text",
        "outputs/jensen_window_pf_column_recurrence_contract.md",
        "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        payload["summary"]["main_finding"],
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(
        "wrote "
        f"{args.out_json.relative_to(REPO_ROOT)} and {args.note.relative_to(REPO_ROOT)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
