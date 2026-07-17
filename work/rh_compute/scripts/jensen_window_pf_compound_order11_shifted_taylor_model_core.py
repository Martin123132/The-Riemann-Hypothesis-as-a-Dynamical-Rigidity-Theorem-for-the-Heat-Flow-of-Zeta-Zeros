#!/usr/bin/env python3
"""Cancellation-preserving Taylor-model core for order eleven.

The older localized hierarchy encloses every stable layer independently.  At
the eighth stable layer that destroys the shifted cancellations and inflates
the sixth normalized coefficient of ``X`` by roughly 26 decimal orders near
``t=1252``.  This module keeps one common local variable ``x`` through all
shifted recurrences instead.

For each exact half-grid anchor it builds Taylor models for
``H^(d)(anchor+j+x)``, ``d=0,1,2``.  Coefficients through H8 come from the
exact point cache; higher coefficients and the final Taylor remainder come
from the immutable localized H cache.  Addition and centered differences are
then coefficient-wise, so their cancellations survive.  Products and the
stable logarithm carry explicit uniform remainder bounds on the exact block.
"""

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

import jensen_window_pf_compound_order10_localized_final_gap_interval_core as order10  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    stable_log_series,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    upper_absolute,
)


PRECISION_BITS = order10.PRECISION_BITS
DERIVATIVE_MODEL_DEGREES = (16, 15, 14)
STABLE_TAYLOR_SURPLUS = 4
POINT_H_MAXIMUM_ORDER = 8
MAXIMUM_H_ORDER = max(
    derivative + degree + 1
    for derivative, degree in enumerate(DERIVATIVE_MODEL_DEGREES)
)
CURVATURE_CONSTANT = 6000

INHERITED_CURVATURE_CONSTANTS = {
    3: Fraction(7, 2),
    4: Fraction(60),
    5: Fraction(200),
    6: Fraction(600),
    7: Fraction(2500),
    8: Fraction(4200),
    9: Fraction(4200),
}


def _maximum(values: list[flint.arb]) -> flint.arb:
    if not values:
        raise ValueError("empty Arb maximum")
    result = values[0]
    for value in values[1:]:
        if bool(value > result):
            result = value
    return result


def _minimum(values: list[flint.arb]) -> flint.arb:
    if not values:
        raise ValueError("empty Arb minimum")
    result = values[0]
    for value in values[1:]:
        if bool(value < result):
            result = value
    return result


def _hull(values: list[flint.arb]) -> flint.arb:
    return order10.hull(values)


def _intersection(values: list[flint.arb]) -> flint.arb:
    if not values:
        raise ValueError("empty Arb intersection")
    lower = _maximum([order10._lower_point(value) for value in values])
    upper = _minimum([order10._upper_point(value) for value in values])
    if not bool(upper >= lower):
        raise RuntimeError("localized endpoint derivative intervals are disjoint")
    return order10._interval_from_endpoints(lower, upper)


class TaylorModel:
    """A polynomial plus a uniform absolute error on one rational interval."""

    def __init__(
        self,
        coefficients: list,
        degree: int,
        domain_left: Fraction,
        domain_right: Fraction,
        remainder=0,
    ) -> None:
        if degree < 0 or not domain_left < domain_right:
            raise ValueError("invalid Taylor-model degree or domain")
        if len(coefficients) > degree + 1:
            raise ValueError("too many Taylor-model coefficients")
        self.degree = degree
        self.domain_left = domain_left
        self.domain_right = domain_right
        self.radius = max(abs(domain_left), abs(domain_right))
        self.x_interval = (
            arb_rational((domain_left + domain_right) / 2)
            + order10._symmetric(arb_rational((domain_right - domain_left) / 2))
        )
        self.coefficients = [flint.arb(value) for value in coefficients]
        self.coefficients.extend(
            flint.arb(0) for _ in range(degree + 1 - len(self.coefficients))
        )
        self.remainder = upper_absolute(
            remainder if isinstance(remainder, flint.arb) else flint.arb(remainder)
        )

    @classmethod
    def zero(
        cls, degree: int, domain_left: Fraction, domain_right: Fraction
    ) -> "TaylorModel":
        return cls([0], degree, domain_left, domain_right)

    @classmethod
    def one(
        cls, degree: int, domain_left: Fraction, domain_right: Fraction
    ) -> "TaylorModel":
        return cls([1], degree, domain_left, domain_right)

    def _same_domain(self, other: "TaylorModel") -> None:
        if (
            self.degree != other.degree
            or self.domain_left != other.domain_left
            or self.domain_right != other.domain_right
        ):
            raise ValueError("Taylor-model degree/domain mismatch")

    def truncate(self, degree: int) -> "TaylorModel":
        if not 0 <= degree <= self.degree:
            raise ValueError("invalid Taylor-model truncation")
        radius = arb_rational(self.radius)
        discarded = flint.arb(0)
        for index in range(degree + 1, self.degree + 1):
            discarded += upper_absolute(self.coefficients[index]) * radius**index
        return TaylorModel(
            self.coefficients[: degree + 1],
            degree,
            self.domain_left,
            self.domain_right,
            self.remainder + discarded,
        )

    def polynomial_range(self) -> flint.arb:
        value = self.coefficients[self.degree]
        for index in range(self.degree - 1, -1, -1):
            value = value * self.x_interval + self.coefficients[index]
        return value

    def range(self) -> flint.arb:
        return self.polynomial_range() + order10._symmetric(self.remainder)

    def polynomial_absolute_range(self) -> flint.arb:
        return upper_absolute(self.polynomial_range())

    def absolute_range(self) -> flint.arb:
        return upper_absolute(self.range())

    def __add__(self, other) -> "TaylorModel":
        if not isinstance(other, TaylorModel):
            other = TaylorModel(
                [other],
                self.degree,
                self.domain_left,
                self.domain_right,
            )
        self._same_domain(other)
        return TaylorModel(
            [
                self.coefficients[index] + other.coefficients[index]
                for index in range(self.degree + 1)
            ],
            self.degree,
            self.domain_left,
            self.domain_right,
            self.remainder + other.remainder,
        )

    __radd__ = __add__

    def __neg__(self) -> "TaylorModel":
        return TaylorModel(
            [-value for value in self.coefficients],
            self.degree,
            self.domain_left,
            self.domain_right,
            self.remainder,
        )

    def __sub__(self, other) -> "TaylorModel":
        return self + (-other)

    def __rsub__(self, other) -> "TaylorModel":
        return TaylorModel(
            [other],
            self.degree,
            self.domain_left,
            self.domain_right,
        ) - self

    def __mul__(self, other) -> "TaylorModel":
        if not isinstance(other, TaylorModel):
            scalar = other if isinstance(other, flint.arb) else flint.arb(other)
            return TaylorModel(
                [scalar * value for value in self.coefficients],
                self.degree,
                self.domain_left,
                self.domain_right,
                upper_absolute(scalar) * self.remainder,
            )
        self._same_domain(other)
        degree = self.degree
        full = [flint.arb(0) for _ in range(2 * degree + 1)]
        for left_index in range(degree + 1):
            for right_index in range(degree + 1):
                full[left_index + right_index] += (
                    self.coefficients[left_index]
                    * other.coefficients[right_index]
                )
        radius = arb_rational(self.radius)
        polynomial_tail = flint.arb(0)
        for index in range(degree + 1, 2 * degree + 1):
            polynomial_tail += upper_absolute(full[index]) * radius**index
        remainder = (
            polynomial_tail
            + self.polynomial_absolute_range() * other.remainder
            + other.polynomial_absolute_range() * self.remainder
            + self.remainder * other.remainder
        )
        return TaylorModel(
            full[: degree + 1],
            degree,
            self.domain_left,
            self.domain_right,
            remainder,
        )

    __rmul__ = __mul__

    def square(self) -> "TaylorModel":
        return self * self


def _stable_derivative_model(
    coordinate: TaylorModel,
    derivative: int,
    analytic_floor: flint.arb,
    output_degree: int,
) -> TaylorModel:
    """Taylor-model ``F^(derivative)(coordinate)`` for F=log(1-exp(-x))."""
    if derivative < 0 or output_degree > coordinate.degree:
        raise ValueError("invalid stable-log derivative model")
    coordinate = coordinate.truncate(output_degree)
    surplus = STABLE_TAYLOR_SURPLUS
    stable_order = output_degree + surplus
    center = coordinate.coefficients[0]
    polynomial_coordinate = TaylorModel(
        coordinate.coefficients,
        output_degree,
        coordinate.domain_left,
        coordinate.domain_right,
    )
    displacement = TaylorModel(
        [0, *coordinate.coefficients[1:]],
        output_degree,
        coordinate.domain_left,
        coordinate.domain_right,
    )

    center_series = stable_log_series(
        [center, flint.arb(1)]
        + [flint.arb(0) for _ in range(derivative + stable_order)],
        derivative + stable_order + 1,
    )
    result = TaylorModel.zero(
        output_degree, coordinate.domain_left, coordinate.domain_right
    )
    power = TaylorModel.one(
        output_degree, coordinate.domain_left, coordinate.domain_right
    )
    for index in range(stable_order + 1):
        coefficient = (
            math.factorial(derivative + index)
            // math.factorial(index)
            * center_series[derivative + index]
        )
        result = result + power * coefficient
        if index != stable_order:
            power = power * displacement

    polynomial_lower = order10._lower_point(polynomial_coordinate.range())
    if not bool(polynomial_lower > 0):
        raise RuntimeError(
            "stable-log polynomial path has no positive lower bound: "
            f"polynomial={polynomial_coordinate.range()}, "
            f"remainder={coordinate.remainder}, floor={analytic_floor}"
        )
    input_path_lower = polynomial_lower - coordinate.remainder
    input_path_floor = (
        analytic_floor
        if bool(analytic_floor > input_path_lower)
        else input_path_lower
    )
    if not bool(input_path_floor > 0):
        raise RuntimeError("stable-log input path has no positive lower bound")

    tail_order = derivative + stable_order + 1
    floor_series = stable_log_series(
        [polynomial_lower, flint.arb(1)]
        + [flint.arb(0) for _ in range(tail_order)],
        tail_order,
    )
    derivative_upper = upper_absolute(
        math.factorial(tail_order) * floor_series[tail_order]
    )
    algebraic_remainder = result.remainder
    analytic_tail = (
        derivative_upper
        * displacement.absolute_range() ** (stable_order + 1)
        / math.factorial(stable_order + 1)
    )
    input_derivative_order = derivative + 1
    input_floor_series = stable_log_series(
        [order10._lower_point(input_path_floor), flint.arb(1)]
        + [flint.arb(0) for _ in range(input_derivative_order)],
        input_derivative_order,
    )
    input_lipschitz = upper_absolute(
        math.factorial(input_derivative_order)
        * input_floor_series[input_derivative_order]
    )
    input_remainder = input_lipschitz * coordinate.remainder
    result.remainder = upper_absolute(
        result.remainder + analytic_tail + input_remainder
    )
    result.composition_diagnostics = {
        "derivative_floor": input_path_floor,
        "polynomial_floor": polynomial_lower,
        "displacement_absolute_range": displacement.absolute_range(),
        "coordinate_remainder": coordinate.remainder,
        "algebraic_remainder": algebraic_remainder,
        "analytic_tail": analytic_tail,
        "input_remainder": input_remainder,
    }
    return result


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
    coordinate: list[TaylorModel],
    stage: int,
    target: Fraction,
) -> list[TaylorModel]:
    domain_left = coordinate[0].domain_left
    domain_right = coordinate[0].domain_right
    floor = _coordinate_floor(stage, target, domain_left, domain_right)
    raw_constant = coordinate[0].coefficients[0]
    raw_lower = order10._lower_point(raw_constant)
    raw_upper = order10._upper_point(raw_constant)
    improved_lower = floor if bool(floor > raw_lower) else raw_lower
    if not bool(raw_upper >= improved_lower):
        raise RuntimeError(
            f"analytic stage-{stage} floor is disjoint from the point model "
            f"at t={target}"
        )
    if bool(improved_lower > raw_lower):
        improved_coefficients = list(coordinate[0].coefficients)
        improved_coefficients[0] = order10._interval_from_endpoints(
            improved_lower, raw_upper
        )
        coordinate = [
            TaylorModel(
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
    stable_value = _stable_derivative_model(
        coordinate[0], 0, floor, value_degree
    )
    stable_first = _stable_derivative_model(
        coordinate[0], 1, floor, first_degree
    ) * coordinate[1]
    coordinate_value_second = coordinate[0].truncate(second_degree)
    coordinate_first_second = coordinate[1].truncate(second_degree)
    stable_second = (
        _stable_derivative_model(
            coordinate_value_second, 1, floor, second_degree
        )
        * coordinate[2]
        + _stable_derivative_model(
            coordinate_value_second, 2, floor, second_degree
        )
        * coordinate_first_second.square()
    )
    return [stable_value, stable_first, stable_second]


def _overlap_rows(
    h_rows: list[dict], left: Fraction, right: Fraction
) -> list[dict]:
    rows = [
        row
        for row in h_rows
        if row["target_t_right"] > left and row["target_t_left"] < right
    ]
    if not rows:
        raise RuntimeError(f"localized H source misses {left}..{right}")
    rows.sort(key=lambda row: row["target_t_left"])
    if rows[0]["target_t_left"] > left or rows[-1]["target_t_right"] < right:
        raise RuntimeError(f"localized H source incompletely covers {left}..{right}")
    for previous, current in zip(rows, rows[1:]):
        if previous["target_t_right"] != current["target_t_left"]:
            raise RuntimeError(f"localized H source has a gap on {left}..{right}")
    return rows


def _point_rows(h_rows: list[dict], target: Fraction) -> list[dict]:
    rows = [
        row
        for row in h_rows
        if row["target_t_left"] <= target <= row["target_t_right"]
    ]
    if not rows:
        raise RuntimeError(f"localized H source misses point {target}")
    return rows


def _h_derivative_model(
    target: Fraction,
    derivative: int,
    degree: int,
    domain_left: Fraction,
    domain_right: Fraction,
    h_rows: list[dict],
    point_series: list[flint.arb],
) -> TaylorModel:
    if len(point_series) <= POINT_H_MAXIMUM_ORDER:
        raise RuntimeError("exact point H source stops before H8")
    block_rows = _overlap_rows(
        h_rows, target + domain_left, target + domain_right
    )
    point_rows = _point_rows(h_rows, target)
    coefficients = []
    for local_order in range(degree + 1):
        h_order = derivative + local_order
        if h_order < len(point_series):
            coefficient = (
                point_series[h_order]
                * math.factorial(h_order)
                / math.factorial(local_order)
            )
        else:
            if any(h_order not in row["H"] for row in point_rows):
                raise RuntimeError(f"localized point source misses H{h_order}")
            coefficient = _intersection(
                [row["H"][h_order] for row in point_rows]
            ) / math.factorial(local_order)
        coefficients.append(coefficient)

    remainder_order = derivative + degree + 1
    if any(remainder_order not in row["H"] for row in block_rows):
        raise RuntimeError(f"localized block source misses H{remainder_order}")
    derivative_upper = _maximum(
        [upper_absolute(row["H"][remainder_order]) for row in block_rows]
    )
    radius = max(abs(domain_left), abs(domain_right))
    remainder = (
        derivative_upper
        * arb_rational(radius) ** (degree + 1)
        / math.factorial(degree + 1)
    )
    return TaylorModel(
        coefficients,
        degree,
        domain_left,
        domain_right,
        remainder,
    )


def shifted_taylor_model_curvature_row(
    expansion_anchor: Fraction,
    block_left: Fraction,
    block_right: Fraction,
    h_rows: list[dict],
    *,
    point_h_source: dict[Fraction, tuple[list[flint.arb], dict]],
    require_pass: bool = True,
) -> dict:
    """Certify ``y_1''`` on one block through the exact shifted recurrence."""
    flint.ctx.prec = PRECISION_BITS
    if not block_left <= expansion_anchor <= block_right or not block_left < block_right:
        raise ValueError("invalid order-eleven Taylor-model block")
    domain_left = block_left - expansion_anchor
    domain_right = block_right - expansion_anchor
    if max(abs(domain_left), abs(domain_right)) > Fraction(1, 4):
        raise ValueError("Taylor-model block radius exceeds one quarter")

    h_models: dict[int, list[TaylorModel]] = {}
    point_diagnostics = []
    for shift in range(-9, 10):
        target = expansion_anchor + shift
        try:
            point_series, diagnostics = point_h_source[target]
        except KeyError as exc:
            raise RuntimeError(f"exact point H source misses t={target}") from exc
        point_diagnostics.append(diagnostics)
        h_models[shift] = [
            _h_derivative_model(
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
        for shift in range(-8, 9)
    }
    current = {
        shift: _stable_function_derivatives(
            b_models[shift], 1, expansion_anchor + shift
        )
        for shift in range(-8, 9)
    }
    previous = None
    stage_diagnostics = []
    coordinate_names = {2: "J", 3: "R", 4: "S", 5: "T", 6: "U", 7: "V", 8: "W", 9: "X"}

    for stage in range(2, 10):
        shifts = range(-(9 - stage), 10 - stage)
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
                coordinates[shift], stage, target
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
            raw_lowers.append(order10._lower_point(coordinates[shift][0].range()))
            analytic_floors.append(
                _coordinate_floor(stage, target, domain_left, domain_right)
            )
            point_lowers.append(
                order10._lower_point(coordinates[shift][0].coefficients[0])
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
                "raw_coordinate_range_lower": arb_lower_text(_minimum(raw_lowers)),
                "analytic_coordinate_floor_lower": arb_lower_text(
                    _minimum(analytic_floors)
                ),
                "point_coordinate_lower": arb_lower_text(_minimum(point_lowers)),
                "maximum_output_remainder": {
                    str(derivative): arb_upper_text(
                        _maximum(output_remainders[derivative])
                    )
                    for derivative in range(3)
                },
                "stable_value_composition": {
                    "minimum_derivative_floor": arb_lower_text(
                        _minimum(stable_value_diagnostics["derivative_floor"])
                    ),
                    "minimum_polynomial_floor": arb_lower_text(
                        _minimum(stable_value_diagnostics["polynomial_floor"])
                    ),
                    "maximum_displacement_absolute_range": arb_upper_text(
                        _maximum(
                            stable_value_diagnostics[
                                "displacement_absolute_range"
                            ]
                        )
                    ),
                    "maximum_coordinate_remainder": arb_upper_text(
                        _maximum(stable_value_diagnostics["coordinate_remainder"])
                    ),
                    "maximum_algebraic_remainder": arb_upper_text(
                        _maximum(stable_value_diagnostics["algebraic_remainder"])
                    ),
                    "maximum_analytic_tail": arb_upper_text(
                        _maximum(stable_value_diagnostics["analytic_tail"])
                    ),
                    "maximum_input_remainder": arb_upper_text(
                        _maximum(stable_value_diagnostics["input_remainder"])
                    ),
                },
            }
        )
        previous, current = current, following

    y_second = current[0][2]
    y_second_range = y_second.range()
    curvature_upper = order10._upper_point(y_second_range)
    if bool(curvature_upper > 0):
        scaled_upper = order10._upper_point(
            arb_rational(block_right) ** 2 * curvature_upper
        )
    else:
        scaled_upper = flint.arb(0)
    margin = flint.arb(CURVATURE_CONSTANT) - scaled_upper
    passed = bool(margin > 0)
    if require_pass and not passed:
        raise RuntimeError(
            f"order-eleven shifted Taylor model failed on "
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
            * order10._upper_point(y_second.coefficients[0])
        ),
        "curvature_range": [
            arb_lower_text(y_second_range),
            arb_upper_text(y_second_range),
        ],
        "curvature_upper": arb_upper_text(curvature_upper),
        "scaled_curvature_upper": arb_upper_text(scaled_upper),
        "curvature_margin_lower": arb_lower_text(margin),
        "final_polynomial_coefficients": [
            value.str(50).replace("e", "E")
            for value in y_second.coefficients
        ],
        "final_uniform_remainder_upper": arb_upper_text(y_second.remainder),
        "stage_diagnostics": stage_diagnostics,
        "quadrature": {
            "shift_count": len(point_diagnostics),
            "maximum_panel_error_upper": max(
                (
                    row["maximum_panel_error_upper"]
                    for row in point_diagnostics
                ),
                key=float,
            ),
            "maximum_tail_moment_upper": max(
                (
                    row["maximum_tail_moment_upper"]
                    for row in point_diagnostics
                ),
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
