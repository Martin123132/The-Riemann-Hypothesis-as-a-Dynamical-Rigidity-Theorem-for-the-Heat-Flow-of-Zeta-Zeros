#!/usr/bin/env python3
"""Check Taylor-model algebra and stable-log composition for order eleven."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_localized_final_gap_interval_core as order10  # noqa: E402
from jensen_window_pf_compound_order11_shifted_taylor_model_core import (  # noqa: E402
    TaylorModel,
    _stable_derivative_model,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DOMAIN_LEFT = Fraction(-1, 4)
DOMAIN_RIGHT = Fraction(1, 4)
SAMPLE_POINTS = tuple(Fraction(index, 16) for index in range(-4, 5))


def evaluate(model: TaylorModel, point: Fraction) -> flint.arb:
    value = model.coefficients[model.degree]
    point_arb = arb_rational(point)
    for index in range(model.degree - 1, -1, -1):
        value = value * point_arb + model.coefficients[index]
    return value + order10._symmetric(model.remainder)


def polynomial_value(coefficients: list[Fraction], point: Fraction) -> flint.arb:
    value = flint.arb(0)
    point_arb = arb_rational(point)
    for coefficient in reversed(coefficients):
        value = value * point_arb + arb_rational(coefficient)
    return value


def polynomial_fraction(coefficients: list[Fraction], point: Fraction) -> Fraction:
    value = Fraction(0)
    for coefficient in reversed(coefficients):
        value = value * point + coefficient
    return value


def stable_derivative(value: flint.arb, derivative: int) -> flint.arb:
    exponential = value.exp()
    if derivative == 0:
        return (flint.arb(1) - (-value).exp()).log()
    if derivative == 1:
        return 1 / (exponential - 1)
    if derivative == 2:
        return -exponential / (exponential - 1) ** 2
    raise ValueError("fixture supports stable derivatives only through two")


def multiplication_checks() -> list[str]:
    left_coefficients = [Fraction(3, 2), Fraction(-5, 8), Fraction(1, 16)]
    right_coefficients = [Fraction(-3, 4), Fraction(1, 8), Fraction(3, 32)]
    left = TaylorModel(
        [arb_rational(value) for value in left_coefficients],
        2,
        DOMAIN_LEFT,
        DOMAIN_RIGHT,
    )
    right = TaylorModel(
        [arb_rational(value) for value in right_coefficients],
        2,
        DOMAIN_LEFT,
        DOMAIN_RIGHT,
    )
    product = left * right
    issues = []
    for point in SAMPLE_POINTS:
        expected = polynomial_fraction(
            left_coefficients, point
        ) * polynomial_fraction(
            right_coefficients, point
        )
        expected_fmpq = flint.fmpq(expected.numerator, expected.denominator)
        if not bool(evaluate(product, point).contains(expected_fmpq)):
            issues.append(f"truncated product misses at x={point}")
    return issues


def stable_composition_checks() -> list[str]:
    coefficients = [Fraction(2), Fraction(1, 10), Fraction(1, 50)]
    perturbation_bound = Fraction(1, 10**12)
    issues = []
    for with_remainder in (False, True):
        remainder = perturbation_bound if with_remainder else Fraction(0)
        coordinate = TaylorModel(
            [arb_rational(value) for value in coefficients],
            6,
            DOMAIN_LEFT,
            DOMAIN_RIGHT,
            arb_rational(remainder),
        )
        models = [
            _stable_derivative_model(
                coordinate,
                derivative,
                arb_rational(Fraction(19, 10)),
                6,
            )
            for derivative in range(3)
        ]
        perturbations = (
            (-perturbation_bound / 2, perturbation_bound / 2)
            if with_remainder
            else (Fraction(0),)
        )
        for point in SAMPLE_POINTS:
            polynomial = polynomial_value(coefficients, point)
            for perturbation in perturbations:
                actual_coordinate = polynomial + arb_rational(perturbation)
                for derivative, model in enumerate(models):
                    expected = stable_derivative(actual_coordinate, derivative)
                    if not bool(evaluate(model, point).contains(expected)):
                        label = "perturbed" if with_remainder else "exact"
                        issues.append(
                            f"{label} F{derivative} composition misses at x={point}"
                        )
    return issues


def main() -> int:
    flint.ctx.prec = 256
    issues = multiplication_checks() + stable_composition_checks()
    if issues:
        print(f"order-eleven shifted Taylor-model algebra: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-eleven shifted Taylor-model algebra: "
        "9 products + 81 stable-log enclosures, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
