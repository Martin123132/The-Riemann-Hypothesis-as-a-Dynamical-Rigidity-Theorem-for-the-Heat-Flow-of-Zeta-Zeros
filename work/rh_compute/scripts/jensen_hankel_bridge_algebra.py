#!/usr/bin/env python3
"""Exact Jensen/Hankel bridge algebra and a finite non-promotion countermodel."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_hankel_bridge_algebra.json"


@dataclass(frozen=True)
class MinorRow:
    order_k: int
    columns: tuple[int, ...]
    determinant: str
    expected_sign: int
    signed_positive: bool


def exact_fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def signed_hankel_sign(order_k: int) -> int:
    return -1 if (order_k * (order_k - 1) // 2) % 2 else 1


def det2(a: list[Fraction], c0: int, c1: int) -> Fraction:
    return a[c0] * a[c1 + 1] - a[c1] * a[c0 + 1]


def det3(a: list[Fraction], columns: tuple[int, int, int]) -> Fraction:
    matrix = [[a[row + column] for column in columns] for row in range(3)]
    return (
        matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0])
    )


def jensen_degree3_discriminant(a: list[Fraction]) -> Fraction:
    """Discriminant of A0 + 3*A1*x + 3*A2*x^2 + A3*x^3."""
    a0, a1, a2, a3 = a[:4]
    return (
        81 * a1 * a1 * a2 * a2
        - 108 * a0 * a2**3
        - 108 * a1**3 * a3
        - 27 * a0 * a0 * a3 * a3
        + 162 * a0 * a1 * a2 * a3
    )


def build_payload() -> dict:
    # A deliberately small rational sequence.  It is not a model of zeta
    # coefficients; it is a proof-safety model for finite bridge claims.
    a = [
        Fraction(1, 1),
        Fraction(33, 40),
        Fraction(429, 640),
        Fraction(4719, 12800),
        Fraction(4719, 35840),
    ]
    k2_rows = [
        MinorRow(
            order_k=2,
            columns=columns,
            determinant=exact_fraction(det2(a, *columns)),
            expected_sign=signed_hankel_sign(2),
            signed_positive=signed_hankel_sign(2) * det2(a, *columns) > 0,
        )
        for columns in combinations(range(3), 2)
    ]
    k3_det = det3(a, (0, 1, 2))
    k3_rows = [
        MinorRow(
            order_k=3,
            columns=(0, 1, 2),
            determinant=exact_fraction(k3_det),
            expected_sign=signed_hankel_sign(3),
            signed_positive=signed_hankel_sign(3) * k3_det > 0,
        )
    ]
    disc3 = jensen_degree3_discriminant(a)
    return {
        "kind": "jensen_hankel_bridge_algebra",
        "date": "2026-07-05",
        "proof_boundary": (
            "Exact low-degree algebra and finite countermodel gate only; not an "
            "all-order sign-consistency theorem, not a Jensen hyperbolicity "
            "theorem, and not a proof of RH or Lambda <= 0."
        ),
        "degree2_identity": {
            "jensen_polynomial": "P_{2,n}(x)=A_n+2*A_{n+1}*x+A_{n+2}*x^2",
            "discriminant": "Delta_2=4*(A_{n+1}^2-A_n*A_{n+2})",
            "hankel_relation": (
                "Delta_2=-4*det([[A_n,A_{n+1}],[A_{n+1},A_{n+2}]])"
            ),
            "conclusion": (
                "The m=1 signed-Hankel condition is exactly degree-2 Jensen "
                "hyperbolicity for positive A_n."
            ),
            "exact_identity_verified": True,
        },
        "degree3_formula": {
            "jensen_polynomial": "P_3(x)=A0+3*A1*x+3*A2*x^2+A3*x^3",
            "discriminant": (
                "81*A1^2*A2^2 - 108*A0*A2^3 - 108*A1^3*A3 "
                "- 27*A0^2*A3^2 + 162*A0*A1*A2*A3"
            ),
            "bridge_warning": (
                "Degree-3 Jensen hyperbolicity is a separate cubic discriminant "
                "condition; finite low-order reshaped-Hankel signs do not by "
                "themselves supply the missing all-order bridge."
            ),
        },
        "finite_countermodel": {
            "sequence_A0_to_A4": [exact_fraction(value) for value in a],
            "all_coefficients_positive": all(value > 0 for value in a),
            "reshaped_k2_minors_N3": [asdict(row) for row in k2_rows],
            "reshaped_k3_minors_N3": [asdict(row) for row in k3_rows],
            "finite_reshaped_signs_pass": all(row.signed_positive for row in k2_rows + k3_rows),
            "degree3_jensen_discriminant": exact_fraction(disc3),
            "degree3_jensen_hyperbolicity_breaks": disc3 < 0,
            "blocked_promotion": (
                "finite reshaped-Hankel sign-consistency for k=2,3 at N=3 "
                "implies degree-3 Jensen hyperbolicity"
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
        print("validated Jensen/Hankel bridge algebra: degree2 identity and degree3 finite countermodel")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
