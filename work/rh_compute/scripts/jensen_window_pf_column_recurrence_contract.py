#!/usr/bin/env python3
"""Build the column-shape recurrence contract for Jensen-window PF.

Column-shape Jacobi-Trudi determinants s_(1^m) are elementary-symmetric
coefficients after normalizing h_j by h_0.  For Jensen windows,

    h_j = binom(d, j) A_{n+j},

the unnormalized column determinant C_m satisfies a finite recurrence.  This
script records that recurrence and checks the hard frontier column shapes
against the exact rational countermodel.  It is a theorem-search contract, not
a positivity proof.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_SHAPE_CONTRACT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json"


def expr_str(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr))


def rational_str(expr: sp.Expr) -> str:
    value = sp.Rational(expr)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def column_recurrence_polynomials(degree: int, max_m: int) -> tuple[tuple[sp.Symbol, ...], list[sp.Expr]]:
    h = sp.symbols(f"h0:{degree + 1}")
    columns = [sp.Integer(1)]
    for size in range(1, max_m + 1):
        value = sp.Integer(0)
        for offset in range(1, min(degree, size) + 1):
            value += (-1) ** (offset - 1) * h[0] ** (offset - 1) * h[offset] * columns[size - offset]
        columns.append(sp.factor(value))
    return h, columns


def direct_column_determinant(degree: int, size: int) -> sp.Expr:
    h = sp.symbols(f"h0:{degree + 1}")
    matrix = [
        [h[col - row] if 0 <= col - row <= degree else sp.Integer(0) for col in range(1, size + 1)]
        for row in range(size)
    ]
    return sp.factor(sp.Matrix(matrix).det())


def recurrence_terms(degree: int) -> list[dict]:
    h = sp.symbols(f"h0:{degree + 1}")
    rows = []
    for offset in range(1, degree + 1):
        sign = 1 if offset % 2 == 1 else -1
        rows.append(
            {
                "offset": offset,
                "sign": "+" if sign > 0 else "-",
                "coefficient": expr_str(h[0] ** (offset - 1) * h[offset]),
                "term": f"{'+' if sign > 0 else '-'} {expr_str(h[0] ** (offset - 1) * h[offset])} * C[m-{offset}]",
            }
        )
    return rows


def hard_row(degree: int, size: int, algebra: dict, shape_contract: dict) -> dict:
    h, columns = column_recurrence_polynomials(degree, size)
    recurrence_poly = columns[size]
    determinant_poly = direct_column_determinant(degree, size)
    a_values = [sp.Rational(item) for item in algebra["finite_countermodel"]["sequence_A0_to_A4"]]
    substitutions = {h[index]: sp.binomial(degree, index) * a_values[index] for index in range(degree + 1)}
    counter_value = sp.factor(recurrence_poly.subs(substitutions))
    expected = algebra["finite_countermodel"][f"d{degree}_first_negative_contiguous_toeplitz_minor"]["determinant"]
    shape_rows = {
        row["id"]: row for row in shape_contract.get("hard_frontier_column_shapes", []) if isinstance(row, dict)
    }
    shape_row = shape_rows[f"d{degree}_column_shape_m{size}"]
    return {
        "id": f"d{degree}_column_recurrence_m{size}",
        "source_shape_row": shape_row["id"],
        "degree": degree,
        "minor_size": size,
        "column_shape": shape_row["schur_label"],
        "recurrence_terms": recurrence_terms(degree),
        "recurrence_polynomial_h": expr_str(recurrence_poly),
        "direct_determinant_polynomial_h": expr_str(determinant_poly),
        "recurrence_matches_determinant": bool(sp.factor(recurrence_poly - determinant_poly) == 0),
        "countermodel_value": rational_str(counter_value),
        "expected_countermodel_value": expected,
        "countermodel_negative": bool(counter_value < 0),
        "matches_obligation_algebra_countermodel": rational_str(counter_value) == expected,
        "structural_obligation": (
            "A positive Schur, network, production-matrix, or determinant-integral route "
            "must prove C_m >= 0 for this column recurrence from Xi/Phi-specific input, "
            "not from generic log-concavity."
        ),
    }


def degree_contract_row(degree: int) -> dict:
    return {
        "degree": degree,
        "normalized_generating_identity": "E(t) = 1 / H(-t), H(t)=1+g1*t+...+gd*t^d, g_j=h_j/h0",
        "unnormalized_column_recurrence": (
            "C[m] = sum_{j=1..min(d,m)} (-1)^(j-1) h0^(j-1) h_j C[m-j], "
            "C[0]=1, C[m<0]=0"
        ),
        "recurrence_terms": recurrence_terms(degree),
        "necessary_condition": "C[m] >= 0 for every m >= 0 and every shift n",
    }


def build_payload(algebra: dict, shape_contract: dict) -> dict:
    hard_rows = [hard_row(3, 8, algebra, shape_contract), hard_row(4, 6, algebra, shape_contract)]
    return {
        "kind": "jensen_window_pf_column_recurrence_contract",
        "date": "2026-07-06",
        "target_obligation": "jwpf_06_sign_regular_to_jensen_pf_conversion",
        "source_algebra": "work/rh_compute/results/jensen_window_pf_obligation_algebra.json",
        "source_schur_shape_contract": "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json",
        "proof_boundary": (
            "Column-shape recurrence diagnostic only; not a proof of Schur positivity, "
            "Jensen-window PF-infinity, Jensen hyperbolicity, RH, or Lambda <= 0."
        ),
        "normalization": {
            "h_j": "binom(d,j) * A_{n+j}",
            "g_j": "h_j / h_0",
            "C_m": "det(h_{1+col-row}) for rows=0..m-1 and cols=1..m",
            "C_m_relation_to_elementary": "C_m = h_0^m * e_m(g_1,...,g_d)",
            "generating_identity": "sum_m e_m t^m = 1 / (1 - g_1*t + g_2*t^2 - ... + (-1)^d*g_d*t^d)",
        },
        "degree_contract_rows": [degree_contract_row(degree) for degree in (2, 3, 4, 5)],
        "hard_frontier_recurrence_rows": hard_rows,
        "summary": {
            "degree_rows": 4,
            "hard_frontier_rows": len(hard_rows),
            "negative_countermodel_rows": sum(1 for row in hard_rows if row["countermodel_negative"]),
            "target_closing": False,
            "main_finding": (
                "The hard column-shape frontier minors are elementary-symmetric "
                "recurrence coefficients for the finite-support Jensen-window "
                "specialization.  Generic log-concavity does not force the "
                "recurrence sequence C_m nonnegative: the exact rational "
                "countermodel makes both hard recurrence rows negative."
            ),
        },
        "invariants": [
            "Every hard recurrence polynomial must match the direct column determinant.",
            "Every hard recurrence countermodel value must match the obligation algebra.",
            "The recurrence is a necessary Schur-positivity condition only, not a proof of all shapes.",
            "No row may set target_closing=true.",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--shape-contract", type=Path, default=DEFAULT_SHAPE_CONTRACT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload(load_json(args.algebra), load_json(args.shape_contract))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        summary = payload["summary"]
        print(
            "wrote Jensen-window PF column recurrence contract: "
            f"{summary['degree_rows']} degree rows, "
            f"{summary['hard_frontier_rows']} hard frontier rows, "
            f"{summary['negative_countermodel_rows']} negative countermodel rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
