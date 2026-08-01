#!/usr/bin/env python3
"""Check the radius-ten Taylor algebra and ninth stable recurrence."""

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
import jensen_window_pf_compound_order12_shifted_taylor_model_core as source  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


TaylorModel = source.base.TaylorModel
DOMAIN_LEFT = Fraction(-1, 4)
DOMAIN_RIGHT = Fraction(1, 4)
EXPANSION_ANCHOR = Fraction(10000)
SAMPLE_POINTS = (Fraction(-1, 4), Fraction(0), Fraction(1, 4))


def evaluate(model: TaylorModel, point: Fraction) -> flint.arb:
    value = model.coefficients[model.degree]
    point_arb = arb_rational(point)
    for index in range(model.degree - 1, -1, -1):
        value = value * point_arb + model.coefficients[index]
    return value + order10._symmetric(model.remainder)


def polynomial_fraction(coefficients: list[Fraction], point: Fraction) -> Fraction:
    value = Fraction(0)
    for coefficient in reversed(coefficients):
        value = value * point + coefficient
    return value


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
        expected = (
            polynomial_fraction(left_coefficients, point)
            * polynomial_fraction(right_coefficients, point)
        )
        exact = flint.fmpq(expected.numerator, expected.denominator)
        if not bool(evaluate(product, point).contains(exact)):
            issues.append(f"truncated product misses at x={point}")
    return issues


def stable_triplet(coordinate: list[flint.arb]) -> list[flint.arb]:
    value, first, second = coordinate
    exponential = value.exp()
    stable_value = (flint.arb(1) - (-value).exp()).log()
    stable_first_derivative = 1 / (exponential - 1)
    stable_second_derivative = -exponential / (exponential - 1) ** 2
    return [
        stable_value,
        stable_first_derivative * first,
        stable_first_derivative * second
        + stable_second_derivative * first**2,
    ]


def b_model(shift: int) -> list[TaylorModel]:
    constant = Fraction(10) + Fraction(shift * shift, 100)
    linear = Fraction(1, 1000)
    quadratic = Fraction(1, 2000)
    value_degree, first_degree, second_degree = source.DERIVATIVE_MODEL_DEGREES
    return [
        TaylorModel(
            [
                arb_rational(constant),
                arb_rational(linear),
                arb_rational(quadratic),
            ],
            value_degree,
            DOMAIN_LEFT,
            DOMAIN_RIGHT,
        ),
        TaylorModel(
            [arb_rational(linear), arb_rational(2 * quadratic)],
            first_degree,
            DOMAIN_LEFT,
            DOMAIN_RIGHT,
        ),
        TaylorModel(
            [arb_rational(2 * quadratic)],
            second_degree,
            DOMAIN_LEFT,
            DOMAIN_RIGHT,
        ),
    ]


def b_triplet(shift: int, point: Fraction) -> list[flint.arb]:
    constant = Fraction(10) + Fraction(shift * shift, 100)
    linear = Fraction(1, 1000)
    quadratic = Fraction(1, 2000)
    return [
        arb_rational(constant + linear * point + quadratic * point**2),
        arb_rational(linear + 2 * quadratic * point),
        arb_rational(2 * quadratic),
    ]


def model_recurrence() -> dict[int, dict[int, list[TaylorModel]]]:
    b_models = {shift: b_model(shift) for shift in range(-9, 10)}
    current = {
        shift: source._stable_function_derivatives(
            b_models[shift],
            1,
            EXPANSION_ANCHOR + shift,
        )
        for shift in range(-9, 10)
    }
    previous = None
    stages = {1: current}
    for stage in range(2, 11):
        shifts = range(-(source.STENCIL_RADIUS - stage), source.STENCIL_RADIUS + 1 - stage)
        coordinates = {
            shift: [
                stage * b_models[shift][derivative]
                - (
                    current[shift - 1][derivative]
                    - 2 * current[shift][derivative]
                    + current[shift + 1][derivative]
                )
                for derivative in range(3)
            ]
            for shift in shifts
        }
        following = {}
        for shift in shifts:
            stable = source._stable_function_derivatives(
                coordinates[shift],
                stage,
                EXPANSION_ANCHOR + shift,
            )
            following[shift] = [
                2 * current[shift][derivative]
                - (previous[shift][derivative] if previous is not None else 0)
                + stable[derivative]
                for derivative in range(3)
            ]
        stages[stage] = following
        previous, current = current, following
    return stages


def direct_recurrence(point: Fraction) -> dict[int, dict[int, list[flint.arb]]]:
    b_values = {shift: b_triplet(shift, point) for shift in range(-9, 10)}
    current = {shift: stable_triplet(b_values[shift]) for shift in range(-9, 10)}
    previous = None
    stages = {1: current}
    for stage in range(2, 11):
        shifts = range(-(source.STENCIL_RADIUS - stage), source.STENCIL_RADIUS + 1 - stage)
        coordinates = {
            shift: [
                stage * b_values[shift][derivative]
                - (
                    current[shift - 1][derivative]
                    - 2 * current[shift][derivative]
                    + current[shift + 1][derivative]
                )
                for derivative in range(3)
            ]
            for shift in shifts
        }
        following = {}
        for shift in shifts:
            stable = stable_triplet(coordinates[shift])
            following[shift] = [
                2 * current[shift][derivative]
                - (previous[shift][derivative] if previous is not None else 0)
                + stable[derivative]
                for derivative in range(3)
            ]
        stages[stage] = following
        previous, current = current, following
    return stages


def ninth_recurrence_checks() -> tuple[list[str], int]:
    issues = []
    comparisons = 0
    models = model_recurrence()
    for point in SAMPLE_POINTS:
        direct = direct_recurrence(point)
        for stage in range(1, 11):
            if set(models[stage]) != set(direct[stage]):
                issues.append(f"stage-{stage} shift geometry changed")
                continue
            for shift in models[stage]:
                for derivative in range(3):
                    comparisons += 1
                    if not bool(
                        evaluate(models[stage][shift][derivative], point).contains(
                            direct[stage][shift][derivative]
                        )
                    ):
                        issues.append(
                            f"D{stage - 1} derivative {derivative} misses at "
                            f"shift={shift}, x={point}"
                        )
                        if len(issues) >= 20:
                            return issues, comparisons
    return issues, comparisons


def main() -> int:
    flint.ctx.prec = 256
    recurrence_issues, comparisons = ninth_recurrence_checks()
    issues = multiplication_checks() + recurrence_issues
    if issues:
        print(f"order-twelve shifted Taylor-model algebra: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-twelve shifted Taylor-model algebra: "
        f"3 products + {comparisons} D0-D9 value/derivative enclosures, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
