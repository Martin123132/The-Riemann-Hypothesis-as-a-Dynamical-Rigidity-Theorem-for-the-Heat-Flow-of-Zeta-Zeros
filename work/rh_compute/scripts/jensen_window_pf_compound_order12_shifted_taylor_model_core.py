#!/usr/bin/env python3
"""Cancellation-preserving quarter-block Taylor core for order twelve."""

from __future__ import annotations

from fractions import Fraction
import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order11_shifted_taylor_model_core as base  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
)


PRECISION_BITS = base.PRECISION_BITS
DERIVATIVE_MODEL_DEGREES = (18, 17, 16)
STABLE_TAYLOR_SURPLUS = base.STABLE_TAYLOR_SURPLUS
POINT_H_MAXIMUM_ORDER = base.POINT_H_MAXIMUM_ORDER
MAXIMUM_H_ORDER = max(
    derivative + degree + 1
    for derivative, degree in enumerate(DERIVATIVE_MODEL_DEGREES)
)
CURVATURE_CONSTANT = 8000
STENCIL_RADIUS = 10

INHERITED_CURVATURE_CONSTANTS = {
    **base.INHERITED_CURVATURE_CONSTANTS,
    10: Fraction(6000),
}


def _coordinate_floor(
    stage: int,
    target: Fraction,
    domain_left: Fraction,
    domain_right: Fraction,
) -> flint.arb:
    left = target + domain_left
    right = target + domain_right
    if not 0 < left < right:
        raise ValueError("invalid coordinate-floor interval")
    if stage == 1:
        floor = 1 / (2 * arb_rational(right) + 3)
    elif stage == 2:
        floor = 1 / (7 * arb_rational(right))
    else:
        try:
            inherited = INHERITED_CURVATURE_CONSTANTS[stage]
        except KeyError as exc:
            raise ValueError(f"unsupported stable stage {stage}") from exc
        floor = (
            flint.arb(stage) / (2 * arb_rational(right) + 3)
            - arb_rational(inherited) / (arb_rational(left) ** 2 - 1)
        )
    if not bool(floor > 0):
        raise RuntimeError(
            f"nonpositive analytic stage-{stage} floor on {left}..{right}"
        )
    return floor


def _stable_function_derivatives(
    coordinate: list[base.TaylorModel],
    stage: int,
    target: Fraction,
) -> list[base.TaylorModel]:
    domain_left = coordinate[0].domain_left
    domain_right = coordinate[0].domain_right
    floor = _coordinate_floor(stage, target, domain_left, domain_right)
    raw_constant = coordinate[0].coefficients[0]
    raw_lower = base.order10._lower_point(raw_constant)
    raw_upper = base.order10._upper_point(raw_constant)
    improved_lower = floor if bool(floor > raw_lower) else raw_lower
    if not bool(raw_upper >= improved_lower):
        raise RuntimeError(
            f"analytic stage-{stage} floor is disjoint from the point model "
            f"at t={target}"
        )
    if bool(improved_lower > raw_lower):
        improved_coefficients = list(coordinate[0].coefficients)
        improved_coefficients[0] = base.order10._interval_from_endpoints(
            improved_lower,
            raw_upper,
        )
        coordinate = [
            base.TaylorModel(
                improved_coefficients,
                coordinate[0].degree,
                domain_left,
                domain_right,
                coordinate[0].remainder,
            ),
            coordinate[1],
            coordinate[2],
        ]
    value_degree, first_degree, second_degree = DERIVATIVE_MODEL_DEGREES
    stable_value = base._stable_derivative_model(
        coordinate[0],
        0,
        floor,
        value_degree,
    )
    stable_first = base._stable_derivative_model(
        coordinate[0],
        1,
        floor,
        first_degree,
    ) * coordinate[1]
    coordinate_value_second = coordinate[0].truncate(second_degree)
    coordinate_first_second = coordinate[1].truncate(second_degree)
    stable_second = (
        base._stable_derivative_model(
            coordinate_value_second,
            1,
            floor,
            second_degree,
        )
        * coordinate[2]
        + base._stable_derivative_model(
            coordinate_value_second,
            2,
            floor,
            second_degree,
        )
        * coordinate_first_second.square()
    )
    return [stable_value, stable_first, stable_second]


def shifted_taylor_model_curvature_row(
    expansion_anchor: Fraction,
    block_left: Fraction,
    block_right: Fraction,
    h_rows: list[dict],
    *,
    point_h_source: dict[Fraction, tuple[list[flint.arb], dict]],
    require_pass: bool = True,
) -> dict:
    """Certify ``v_1''`` on one quarter block through the shifted recurrence."""
    flint.ctx.prec = PRECISION_BITS
    if (
        not block_left <= expansion_anchor <= block_right
        or not block_left < block_right
    ):
        raise ValueError("invalid order-twelve Taylor-model block")
    domain_left = block_left - expansion_anchor
    domain_right = block_right - expansion_anchor
    if max(abs(domain_left), abs(domain_right)) > Fraction(1, 4):
        raise ValueError("Taylor-model block radius exceeds one quarter")

    h_models: dict[int, list[base.TaylorModel]] = {}
    point_diagnostics = []
    for shift in range(-STENCIL_RADIUS, STENCIL_RADIUS + 1):
        target = expansion_anchor + shift
        try:
            point_series, diagnostics = point_h_source[target]
        except KeyError as exc:
            raise RuntimeError(f"exact point H source misses t={target}") from exc
        point_diagnostics.append(diagnostics)
        h_models[shift] = [
            base._h_derivative_model(
                target,
                derivative,
                degree,
                domain_left,
                domain_right,
                h_rows,
                point_series,
            )
            for derivative, degree in enumerate(DERIVATIVE_MODEL_DEGREES)
        ]

    b_models = {
        shift: [
            h_models[shift - 1][derivative]
            - 2 * h_models[shift][derivative]
            + h_models[shift + 1][derivative]
            for derivative in range(3)
        ]
        for shift in range(-(STENCIL_RADIUS - 1), STENCIL_RADIUS)
    }
    current = {
        shift: _stable_function_derivatives(
            b_models[shift],
            1,
            expansion_anchor + shift,
        )
        for shift in range(-(STENCIL_RADIUS - 1), STENCIL_RADIUS)
    }
    previous = None
    stage_diagnostics = []
    coordinate_names = {
        2: "J",
        3: "R",
        4: "S",
        5: "T",
        6: "U",
        7: "V",
        8: "W",
        9: "X",
        10: "D9",
    }

    for stage in range(2, 11):
        shifts = range(-(STENCIL_RADIUS - stage), STENCIL_RADIUS + 1 - stage)
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
        raw_lowers = []
        analytic_floors = []
        point_lowers = []
        output_remainders = [[], [], []]
        stable_value_diagnostics = {
            "derivative_floor": [],
            "polynomial_floor": [],
            "displacement_absolute_range": [],
            "coordinate_remainder": [],
            "algebraic_remainder": [],
            "analytic_tail": [],
            "input_remainder": [],
        }
        for shift in shifts:
            target = expansion_anchor + shift
            stable = _stable_function_derivatives(
                coordinates[shift],
                stage,
                target,
            )
            following[shift] = [
                2 * current[shift][derivative]
                - (
                    previous[shift][derivative]
                    if previous is not None
                    else 0
                )
                + stable[derivative]
                for derivative in range(3)
            ]
            raw_lowers.append(
                base.order10._lower_point(coordinates[shift][0].range())
            )
            analytic_floors.append(
                _coordinate_floor(stage, target, domain_left, domain_right)
            )
            point_lowers.append(
                base.order10._lower_point(coordinates[shift][0].coefficients[0])
            )
            for derivative in range(3):
                output_remainders[derivative].append(
                    following[shift][derivative].remainder
                )
            for name in stable_value_diagnostics:
                stable_value_diagnostics[name].append(
                    stable[0].composition_diagnostics[name]
                )
        stage_diagnostics.append(
            {
                "stage": stage,
                "coordinate": coordinate_names[stage],
                "shift_count": len(list(shifts)),
                "raw_coordinate_range_lower": arb_lower_text(
                    base._minimum(raw_lowers)
                ),
                "analytic_coordinate_floor_lower": arb_lower_text(
                    base._minimum(analytic_floors)
                ),
                "point_coordinate_lower": arb_lower_text(
                    base._minimum(point_lowers)
                ),
                "maximum_output_remainder": {
                    str(derivative): arb_upper_text(
                        base._maximum(output_remainders[derivative])
                    )
                    for derivative in range(3)
                },
                "stable_value_composition": {
                    "minimum_derivative_floor": arb_lower_text(
                        base._minimum(stable_value_diagnostics["derivative_floor"])
                    ),
                    "minimum_polynomial_floor": arb_lower_text(
                        base._minimum(stable_value_diagnostics["polynomial_floor"])
                    ),
                    "maximum_displacement_absolute_range": arb_upper_text(
                        base._maximum(
                            stable_value_diagnostics[
                                "displacement_absolute_range"
                            ]
                        )
                    ),
                    "maximum_coordinate_remainder": arb_upper_text(
                        base._maximum(
                            stable_value_diagnostics["coordinate_remainder"]
                        )
                    ),
                    "maximum_algebraic_remainder": arb_upper_text(
                        base._maximum(
                            stable_value_diagnostics["algebraic_remainder"]
                        )
                    ),
                    "maximum_analytic_tail": arb_upper_text(
                        base._maximum(stable_value_diagnostics["analytic_tail"])
                    ),
                    "maximum_input_remainder": arb_upper_text(
                        base._maximum(stable_value_diagnostics["input_remainder"])
                    ),
                },
            }
        )
        previous, current = current, following

    v_second = current[0][2]
    v_second_range = v_second.range()
    curvature_upper = base.order10._upper_point(v_second_range)
    if bool(curvature_upper > 0):
        scaled_upper = base.order10._upper_point(
            arb_rational(block_right) ** 2 * curvature_upper
        )
    else:
        scaled_upper = flint.arb(0)
    margin = flint.arb(CURVATURE_CONSTANT) - scaled_upper
    passed = bool(margin > 0)
    if require_pass and not passed:
        raise RuntimeError(
            f"order-twelve shifted Taylor model failed on "
            f"{block_left}..{block_right}: scaled={scaled_upper}"
        )

    return {
        "anchor": str(block_left),
        "expansion_anchor": str(expansion_anchor),
        "right": str(block_right),
        "width": str(block_right - block_left),
        "local_domain": [str(domain_left), str(domain_right)],
        "model_degrees": list(DERIVATIVE_MODEL_DEGREES),
        "stable_taylor_surplus": STABLE_TAYLOR_SURPLUS,
        "maximum_h_derivative_order": MAXIMUM_H_ORDER,
        "point_scaled_curvature": arb_upper_text(
            arb_rational(expansion_anchor) ** 2
            * base.order10._upper_point(v_second.coefficients[0])
        ),
        "curvature_range": [
            arb_lower_text(v_second_range),
            arb_upper_text(v_second_range),
        ],
        "curvature_upper": arb_upper_text(curvature_upper),
        "scaled_curvature_upper": arb_upper_text(scaled_upper),
        "curvature_margin_lower": arb_lower_text(margin),
        "final_polynomial_coefficients": [
            value.str(50).replace("e", "E")
            for value in v_second.coefficients
        ],
        "final_uniform_remainder_upper": arb_upper_text(v_second.remainder),
        "stage_diagnostics": stage_diagnostics,
        "quadrature": {
            "shift_count": len(point_diagnostics),
            "maximum_panel_error_upper": max(
                (row["maximum_panel_error_upper"] for row in point_diagnostics),
                key=float,
            ),
            "maximum_tail_moment_upper": max(
                (row["maximum_tail_moment_upper"] for row in point_diagnostics),
                key=float,
            ),
            "mode_brackets": [
                row["mode_bracket"] for row in point_diagnostics
            ],
        },
        "proof_formula": (
            "Delta^2 f(t+x)=f(t-1+x)-2f(t+x)+f(t+1+x); "
            "F'=1/(exp(C)-1); F''=-exp(C)/(exp(C)-1)^2; "
            "g''=2f''-p''+F'(C)C''+F''(C)(C')^2"
        ),
        "passed": passed,
    }
