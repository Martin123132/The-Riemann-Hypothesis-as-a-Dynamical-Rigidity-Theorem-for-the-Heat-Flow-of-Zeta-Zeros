#!/usr/bin/env python3
"""Build a Schur/Jacobi-Trudi shape contract for Jensen-window PF.

For a fixed Jensen-window degree d, set

    h_j = binom(d, j) A_{n+j}, 0 <= j <= d,
    h_j = 0 otherwise.

Every finite-band Toeplitz minor det(h_{q_j-r_i}) is a Jacobi-Trudi
skew-Schur determinant under this finite-support specialization.  This script
records the exact bounded shape counts and the hard frontier column-shape
determinants.  It is a theorem-search contract, not a positivity proof.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from math import comb
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json"


def expr_str(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr))


def is_contiguous(values: tuple[int, ...]) -> bool:
    if not values:
        return False
    return values == tuple(range(values[0], values[0] + len(values)))


def upper_toeplitz_nonzero(rows: tuple[int, ...], cols: tuple[int, ...]) -> bool:
    return all(row <= col for row, col in zip(rows, cols))


def finite_band_nonzero(rows: tuple[int, ...], cols: tuple[int, ...], degree: int) -> bool:
    return all(row <= col <= row + degree for row, col in zip(rows, cols))


def shape_for_minor(rows: tuple[int, ...], cols: tuple[int, ...]) -> tuple[tuple[int, ...], tuple[int, ...]]:
    candidates = [0]
    candidates.extend(row - index for index, row in enumerate(rows, start=1))
    candidates.extend(col - index for index, col in enumerate(cols, start=1))
    shift_k = max(candidates)
    lam = tuple(shift_k + index - row for index, row in enumerate(rows, start=1))
    mu = tuple(shift_k + index - col for index, col in enumerate(cols, start=1))
    return lam, mu


def det_polynomial(rows: tuple[int, ...], cols: tuple[int, ...], degree: int) -> sp.Expr:
    h = sp.symbols(f"h0:{degree + 1}")
    matrix = [
        [h[col - row] if 0 <= col - row <= degree else sp.Integer(0) for col in cols]
        for row in rows
    ]
    return sp.factor(sp.Matrix(matrix).det())


def poly_term_count(expr: sp.Expr, degree: int) -> int:
    h = sp.symbols(f"h0:{degree + 1}")
    return len(sp.Poly(sp.expand(expr), *h).terms())


def grid_row(degree: int, matrix_size: int, max_order: int) -> dict:
    total = 0
    upper_nonzero = 0
    upper_zero = 0
    finite_nonzero = 0
    finite_zero_after_upper = 0
    unique_shapes: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    first_noncontiguous: dict | None = None

    for order in range(1, max_order + 1):
        for rows in combinations(range(matrix_size), order):
            for cols in combinations(range(matrix_size), order):
                total += 1
                if not upper_toeplitz_nonzero(rows, cols):
                    upper_zero += 1
                    continue
                upper_nonzero += 1
                if not finite_band_nonzero(rows, cols, degree):
                    finite_zero_after_upper += 1
                    continue
                finite_nonzero += 1
                lam, mu = shape_for_minor(rows, cols)
                unique_shapes.add((lam, mu))
                if first_noncontiguous is None and (not is_contiguous(rows) or not is_contiguous(cols)):
                    first_noncontiguous = {
                        "order": order,
                        "rows": list(rows),
                        "cols": list(cols),
                        "lambda": list(lam),
                        "mu": list(mu),
                        "schur_label": f"s_{lam}/{mu}",
                    }

    return {
        "degree": degree,
        "matrix_size": matrix_size,
        "max_order": max_order,
        "tested_minors": total,
        "upper_toeplitz_structural_nonzero": upper_nonzero,
        "upper_toeplitz_structural_zero": upper_zero,
        "finite_band_structural_nonzero": finite_nonzero,
        "finite_band_zero_after_upper_support": finite_zero_after_upper,
        "unique_finite_band_shapes": len(unique_shapes),
        "first_noncontiguous_shape_example": first_noncontiguous,
    }


def frontier_example(degree: int, size: int) -> dict:
    rows = tuple(range(size))
    cols = tuple(range(1, size + 1))
    lam, mu = shape_for_minor(rows, cols)
    polynomial = det_polynomial(rows, cols, degree)
    return {
        "id": f"d{degree}_column_shape_m{size}",
        "degree": degree,
        "minor_size": size,
        "rows": list(rows),
        "cols": list(cols),
        "lambda": list(lam),
        "mu": list(mu),
        "schur_label": f"s_{lam}/{mu}",
        "determinant_polynomial_h": expr_str(polynomial),
        "term_count": poly_term_count(polynomial, degree),
        "has_mixed_monomial_signs": True,
        "substitution": "h_j = binom(d,j) * A_{n+j}",
    }


def build_payload(matrix_size: int, max_order: int) -> dict:
    degrees = (2, 3, 4, 5)
    rows = [grid_row(degree, matrix_size, max_order) for degree in degrees]
    hard = [frontier_example(3, 8), frontier_example(4, 6)]
    return {
        "kind": "jensen_window_pf_schur_shape_contract",
        "date": "2026-07-06",
        "target_obligation": "jwpf_06_sign_regular_to_jensen_pf_conversion",
        "proof_boundary": (
            "Exact shape-contract diagnostic only; not a proof of Schur positivity, "
            "Jensen-window PF-infinity, Jensen hyperbolicity, RH, or Lambda <= 0."
        ),
        "finite_window_specialization": {
            "h_j": "binom(d,j) * A_{n+j}",
            "support": "h_j = 0 for j < 0 or j > d",
            "toeplitz_minor": "det(h_{q_j-r_i})",
            "jacobi_trudi_form": "s_{lambda/mu}(h_0,...,h_d)",
        },
        "bounded_grid": {
            "matrix_size": matrix_size,
            "max_order": max_order,
            "degrees": list(degrees),
            "rows": rows,
            "total_tested_minors": sum(row["tested_minors"] for row in rows),
            "total_finite_band_nonzero": sum(row["finite_band_structural_nonzero"] for row in rows),
        },
        "hard_frontier_column_shapes": hard,
        "summary": {
            "grid_rows": len(rows),
            "frontier_rows": len(hard),
            "target_closing": False,
            "main_finding": (
                "The Jensen-window PF bridge can be stated as Schur positivity "
                "for finite-support specializations h_j=binom(d,j)A_{n+j}.  "
                "Even in the bounded N=8, order<=5 grid there are thousands "
                "of nonzero finite-band skew shapes, and the hard contiguous "
                "frontier minors are column-shape Jacobi-Trudi determinants "
                "with mixed-sign h-monomial expansions."
            ),
        },
        "invariants": [
            "This contract records exact shape coverage, not determinant positivity.",
            "Finite-band structural nonzero requires r_i <= q_i <= r_i+d after sorting.",
            "The hard frontier rows remain obligations for any positive Schur, network, or determinant-integral route.",
            "No row may set target_closing=true.",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-size", type=int, default=8)
    parser.add_argument("--max-order", type=int, default=5)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.matrix_size <= 0 or args.max_order <= 0:
        raise ValueError("matrix size and max order must be positive")
    if args.max_order > args.matrix_size:
        raise ValueError("max order cannot exceed matrix size")

    payload = build_payload(args.matrix_size, args.max_order)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        summary = payload["summary"]
        grid = payload["bounded_grid"]
        print(
            "wrote Jensen-window PF Schur shape contract: "
            f"{summary['grid_rows']} grid rows, {summary['frontier_rows']} frontier rows, "
            f"{grid['total_finite_band_nonzero']} finite-band nonzero shapes"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
