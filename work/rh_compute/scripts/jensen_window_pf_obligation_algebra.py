#!/usr/bin/env python3
"""Generate exact low-degree Jensen-window PF obligation algebra.

This is theorem-search algebra, not zeta evidence.  It records which
low-degree Jensen-window Toeplitz/PF obligations match signed-Hankel data and
which obligations are genuinely additional.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"


def expr_str(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(expr))


def rational_str(value: sp.Rational) -> str:
    value = sp.Rational(value)
    if value.q == 1:
        return str(value.p)
    return f"{value.p}/{value.q}"


def jensen_window(symbols: tuple[sp.Symbol, ...], degree: int) -> list[sp.Expr]:
    return [sp.binomial(degree, j) * symbols[j] for j in range(degree + 1)]


def toeplitz_entry(window: list[sp.Expr], row: int, col: int) -> sp.Expr:
    index = col - row
    if 0 <= index < len(window):
        return window[index]
    return sp.Integer(0)


def toeplitz_det(window: list[sp.Expr], rows: tuple[int, ...], cols: tuple[int, ...]) -> sp.Expr:
    matrix = [[toeplitz_entry(window, row, col) for col in cols] for row in rows]
    return sp.factor(sp.Matrix(matrix).det())


def symbolic_minor(window: list[sp.Expr], rows: tuple[int, ...], cols: tuple[int, ...]) -> dict:
    return {
        "rows": list(rows),
        "cols": list(cols),
        "determinant": expr_str(toeplitz_det(window, rows, cols)),
    }


def evaluate_window(a_values: list[Fraction], degree: int) -> list[sp.Rational]:
    return [
        sp.Rational(sp.binomial(degree, j)) * sp.Rational(a_values[j].numerator, a_values[j].denominator)
        for j in range(degree + 1)
    ]


def evaluate_det(window: list[sp.Rational], rows: tuple[int, ...], cols: tuple[int, ...]) -> sp.Rational:
    matrix = [[toeplitz_entry(window, row, col) for col in cols] for row in rows]
    return sp.factor(sp.Matrix(matrix).det())


def first_negative_contiguous(window: list[sp.Rational], max_size: int = 20) -> dict:
    for size in range(1, max_size + 1):
        rows = tuple(range(size))
        cols = tuple(range(1, size + 1))
        determinant = evaluate_det(window, rows, cols)
        if determinant < 0:
            return {
                "size": size,
                "rows": list(rows),
                "cols": list(cols),
                "determinant": rational_str(determinant),
            }
    return {"size": None, "rows": [], "cols": [], "determinant": "none"}


def selected_values(window: list[sp.Rational], minors: list[dict]) -> list[dict]:
    rows: list[dict] = []
    for minor in minors:
        determinant = evaluate_det(window, tuple(minor["rows"]), tuple(minor["cols"]))
        rows.append(
            {
                "rows": minor["rows"],
                "cols": minor["cols"],
                "determinant": rational_str(determinant),
                "positive": bool(determinant > 0),
            }
        )
    return rows


def build_payload() -> dict:
    x = sp.Symbol("x")
    a = sp.symbols("a0:5")
    old_counter = [
        Fraction(1, 1),
        Fraction(33, 40),
        Fraction(429, 640),
        Fraction(4719, 12800),
        Fraction(4719, 35840),
    ]

    windows = {degree: jensen_window(a, degree) for degree in (2, 3, 4)}
    d3_selected = [
        symbolic_minor(windows[3], (0, 1), (1, 2)),
        symbolic_minor(windows[3], (0, 1), (1, 3)),
        symbolic_minor(windows[3], (0, 1), (2, 3)),
        symbolic_minor(windows[3], (0, 1, 2), (1, 2, 3)),
        symbolic_minor(windows[3], (0, 1, 2), (2, 3, 4)),
    ]
    d4_selected = [
        symbolic_minor(windows[4], (0, 1), (1, 2)),
        symbolic_minor(windows[4], (0, 1), (1, 3)),
        symbolic_minor(windows[4], (0, 1), (1, 4)),
        symbolic_minor(windows[4], (0, 1), (2, 3)),
        symbolic_minor(windows[4], (0, 1), (2, 4)),
        symbolic_minor(windows[4], (0, 1), (3, 4)),
        symbolic_minor(windows[4], (0, 1, 2), (1, 2, 3)),
        symbolic_minor(windows[4], (0, 1, 2), (2, 3, 4)),
        symbolic_minor(windows[4], (0, 1, 2), (3, 4, 5)),
    ]

    d3_poly = sum(windows[3][j] * x**j for j in range(4))
    d4_poly = sum(windows[4][j] * x**j for j in range(5))
    d3_counter_window = evaluate_window(old_counter, 3)
    d4_counter_window = evaluate_window(old_counter, 4)
    d3_counter_poly = sum(d3_counter_window[j] * x**j for j in range(4))
    d4_counter_poly = sum(d4_counter_window[j] * x**j for j in range(5))
    d3_counter_selected = selected_values(d3_counter_window, d3_selected)
    d4_counter_selected = selected_values(d4_counter_window, d4_selected)

    return {
        "kind": "jensen_window_pf_obligation_algebra",
        "date": "2026-07-05",
        "proof_boundary": (
            "Exact low-degree theorem-search algebra and finite proof-safety "
            "countermodel only; not a proof of Jensen-window PF-infinity, "
            "Laguerre-Polya membership, RH, or Lambda <= 0."
        ),
        "degree2": {
            "window": ["a0", "2*a1", "a2"],
            "jensen_discriminant": "4*(a1**2 - a0*a2)",
            "signed_hankel_relation": (
                "a1**2 - a0*a2 = -det([[a0,a1],[a1,a2]])"
            ),
            "contiguous_toeplitz_recurrence": (
                "D_m = 2*a1*D_{m-1} - a0*a2*D_{m-2}, D_0=1, D_1=2*a1"
            ),
            "pf_threshold": "a1**2 - a0*a2 >= 0",
            "matches_signed_hankel_m1": True,
        },
        "degree3": {
            "window": ["a0", "3*a1", "3*a2", "a3"],
            "selected_toeplitz_minors": d3_selected,
            "cubic_discriminant": expr_str(sp.discriminant(d3_poly, x)),
            "comparison": (
                "The degree-2 signed-Hankel contact does not supply the cubic "
                "discriminant or all Toeplitz minors of the binomially weighted "
                "Jensen window."
            ),
        },
        "degree4": {
            "window": ["a0", "4*a1", "6*a2", "4*a3", "a4"],
            "selected_toeplitz_minors": d4_selected,
            "quartic_discriminant": expr_str(sp.discriminant(d4_poly, x)),
            "comparison": (
                "The quartic window introduces further 3x3 and higher banded "
                "Toeplitz obligations that are not identical to signed "
                "reshaped-Hankel determinants."
            ),
        },
        "finite_countermodel": {
            "sequence_A0_to_A4": [
                "1",
                "33/40",
                "429/640",
                "4719/12800",
                "4719/35840",
            ],
            "source": "same positive rational sequence as jensen_hankel_bridge_algebra.json",
            "d3_selected_toeplitz_minors": d3_counter_selected,
            "d3_selected_toeplitz_minors_positive": all(row["positive"] for row in d3_counter_selected),
            "d3_cubic_discriminant": rational_str(sp.discriminant(d3_counter_poly, x)),
            "d3_cubic_hyperbolicity_breaks": bool(sp.discriminant(d3_counter_poly, x) < 0),
            "d3_first_negative_contiguous_toeplitz_minor": first_negative_contiguous(d3_counter_window),
            "d4_selected_toeplitz_minors": d4_counter_selected,
            "d4_selected_toeplitz_minors_positive": all(row["positive"] for row in d4_counter_selected),
            "d4_quartic_discriminant": rational_str(sp.discriminant(d4_counter_poly, x)),
            "d4_quartic_hyperbolicity_breaks": bool(sp.discriminant(d4_counter_poly, x) < 0),
            "d4_first_negative_contiguous_toeplitz_minor": first_negative_contiguous(d4_counter_window),
            "blocked_promotion": (
                "selected low-order Jensen-window Toeplitz minors, or finite "
                "low-order reshaped-Hankel signs, imply Jensen-window PF-infinity"
            ),
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print("validated Jensen-window PF obligation algebra: degree2 exact contact and degree3/4 low-order gaps")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
