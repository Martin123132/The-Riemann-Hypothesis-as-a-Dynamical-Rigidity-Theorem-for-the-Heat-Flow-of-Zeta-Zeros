#!/usr/bin/env python3
"""Interval core for the localized first-summand order-four curvature bound."""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_interval,
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    exact_lower,
    exact_upper,
    positive_interval_from_upper,
    tail_polynomial,
    upper_absolute,
)


DEFAULT_PRECISION_BITS = 256
DEFAULT_MAX_MOMENT = 8


def derivative_power(curvature: flint.arb, order: int) -> flint.arb:
    return curvature ** (flint.arb(order) / 2)


def signed_hurwitz_gamma_derivative(order: int, argument: flint.arb) -> flint.arb:
    sign = 1 if order % 2 == 0 else -1
    return sign * math.factorial(order - 1) * flint.arb(order).zeta(argument)


def integrate_h_derivatives(
    left: Fraction,
    right: Fraction,
    panels: int,
    *,
    window_y: int,
    eighth_envelope_bound: Fraction,
    max_moment: int = DEFAULT_MAX_MOMENT,
) -> dict:
    """Enclose H^(2)..H^(max_moment) throughout a mode interval."""
    if panels <= 0 or panels % 2:
        raise ValueError("Simpson panels must be a positive even integer")
    if max_moment < 2:
        raise ValueError("max_moment must be at least two")

    mode = arb_interval(left, right)
    jet = potential_jet_arb(mode, 7)
    t, curvature = jet[1], jet[2]
    if not bool(curvature > 0):
        return {"passed": False, "failure": "nonpositive-curvature"}

    y_cut = flint.arb(window_y)
    panel_width = y_cut / panels
    eighth_envelope = arb_rational(eighth_envelope_bound)
    normalized = [flint.arb(0) for _ in range(8)]
    normalized[2] = flint.arb(1)
    for order in range(3, 8):
        normalized[order] = jet[order] / derivative_power(curvature, order)

    branch_moments = [
        [flint.arb(0) for _ in range(max_moment + 1)] for _ in range(2)
    ]
    for node in range(panels + 1):
        y = panel_width * node
        coefficient = 1 if node in (0, panels) else (4 if node % 2 else 2)
        for branch_index, sign in enumerate((1, -1)):
            potential = y**2 / 2
            for order in range(3, 8):
                potential += (
                    normalized[order]
                    * sign**order
                    * y**order
                    / math.factorial(order)
                )
            value_radius = eighth_envelope * y**8 / math.factorial(8)
            weight = (-potential + flint.arb(0, value_radius)).exp()
            for power in range(max_moment + 1):
                branch_moments[branch_index][power] += (
                    coefficient * y**power * weight
                )
    for branch_index in range(2):
        for power in range(max_moment + 1):
            branch_moments[branch_index][power] *= panel_width / 3

    simpson_errors = [
        [flint.arb(0) for _ in range(max_moment + 1)] for _ in range(2)
    ]
    pair_error_factor = (2 * panel_width) ** 5 / 2880
    for pair_start in range(0, panels, 2):
        y_left = panel_width * pair_start
        y_right = panel_width * (pair_start + 2)
        y_interval = (y_left + y_right) / 2 + flint.arb(
            0, (y_right - y_left) / 2
        )
        for branch_index, sign in enumerate((1, -1)):
            derivatives = [flint.arb(0) for _ in range(5)]
            derivatives[0] = y_interval**2 / 2
            derivatives[1] = y_interval
            derivatives[2] = flint.arb(1)
            for derivative_order in range(5):
                for order in range(max(3, derivative_order), 8):
                    derivatives[derivative_order] += (
                        normalized[order]
                        * sign**order
                        * y_interval ** (order - derivative_order)
                        / math.factorial(order - derivative_order)
                    )
                derivatives[derivative_order] += flint.arb(
                    0,
                    eighth_envelope
                    * y_right ** (8 - derivative_order)
                    / math.factorial(8 - derivative_order),
                )

            w0, w1, w2, w3, w4 = derivatives
            weight = (-w0).exp()
            exponential_derivatives = [
                weight,
                -w1 * weight,
                (w1**2 - w2) * weight,
                (-w1**3 + 3 * w1 * w2 - w3) * weight,
                (
                    w1**4
                    - 6 * w1**2 * w2
                    + 3 * w2**2
                    + 4 * w1 * w3
                    - w4
                )
                * weight,
            ]
            for power in range(max_moment + 1):
                fourth_derivative = flint.arb(0)
                for polynomial_order in range(min(power, 4) + 1):
                    polynomial_derivative = (
                        math.factorial(power)
                        // math.factorial(power - polynomial_order)
                    ) * y_interval ** (power - polynomial_order)
                    fourth_derivative += (
                        math.comb(4, polynomial_order)
                        * polynomial_derivative
                        * exponential_derivatives[4 - polynomial_order]
                    )
                simpson_errors[branch_index][power] += (
                    upper_absolute(fourth_derivative) * pair_error_factor
                )

    for branch_index in range(2):
        for power in range(max_moment + 1):
            branch_moments[branch_index][power] += flint.arb(
                0, simpson_errors[branch_index][power]
            )

    minimum_tail_slope: flint.arb | None = None
    maximum_tail_moment: flint.arb | None = None
    for branch_index, sign in enumerate((1, -1)):
        endpoint = y_cut**2 / 2
        outward_slope = y_cut
        for order in range(3, 8):
            endpoint += (
                normalized[order]
                * sign**order
                * y_cut**order
                / math.factorial(order)
            )
            outward_slope += (
                normalized[order]
                * sign**order
                * y_cut ** (order - 1)
                / math.factorial(order - 1)
            )
        endpoint -= eighth_envelope * y_cut**8 / math.factorial(8)
        outward_slope -= eighth_envelope * y_cut**7 / math.factorial(7)
        if not bool(endpoint > 0 and outward_slope > 0):
            return {"passed": False, "failure": "tail-endpoint-or-slope"}
        if minimum_tail_slope is None or outward_slope.lower() < minimum_tail_slope.lower():
            minimum_tail_slope = outward_slope
        for power in range(max_moment + 1):
            tail_upper = (-exact_lower(endpoint)).exp() * tail_polynomial(
                power, exact_lower(outward_slope), y_cut
            )
            branch_moments[branch_index][power] += positive_interval_from_upper(
                tail_upper
            )
            if maximum_tail_moment is None or tail_upper.upper() > maximum_tail_moment.upper():
                maximum_tail_moment = tail_upper

    assert minimum_tail_slope is not None and maximum_tail_moment is not None
    moments = [flint.arb(0) for _ in range(max_moment + 1)]
    for power in range(max_moment + 1):
        moments[power] = branch_moments[0][power] + (-1) ** power * branch_moments[1][power]
    if not bool(moments[0] > 0):
        return {"passed": False, "failure": "nonpositive-normalizer"}

    raw = [value / moments[0] for value in moments]
    cumulants = [flint.arb(0) for _ in range(max_moment + 1)]
    cumulants[1] = raw[1]
    for order in range(2, max_moment + 1):
        cumulants[order] = raw[order] - sum(
            math.comb(order - 1, previous - 1)
            * cumulants[previous]
            * raw[order - previous]
            for previous in range(1, order)
        )

    argument = t + flint.arb(1) / 2
    h_derivatives = {}
    for order in range(2, max_moment + 1):
        h_derivatives[order] = signed_hurwitz_gamma_derivative(
            order, argument
        ) - cumulants[order] / derivative_power(curvature, order)

    return {
        "passed": True,
        "failure": None,
        "mode_left": str(left),
        "mode_right": str(right),
        "panels": panels,
        "window_y": window_y,
        "eighth_envelope": str(eighth_envelope_bound),
        "t_ball": t.str(40).replace("e", "E"),
        "curvature_ball": curvature.str(40).replace("e", "E"),
        "h_derivatives": h_derivatives,
        "normalizer_value": moments[0],
        "raw_moments": {order: raw[order] for order in range(max_moment + 1)},
        "cumulants": {
            order: cumulants[order] for order in range(1, max_moment + 1)
        },
        "normalizer_lower": arb_lower_text(moments[0]),
        "minimum_tail_slope_lower": arb_lower_text(minimum_tail_slope),
        "maximum_tail_moment_upper": arb_upper_text(maximum_tail_moment),
        "maximum_simpson_error_upper": arb_upper_text(
            max(
                (error for branch in simpson_errors for error in branch),
                key=lambda value: value.upper(),
            )
        ),
    }


def integrate_h_jet_exact_potential(
    left: Fraction,
    right: Fraction,
    panels: int,
    *,
    window_y: int,
    max_moment: int = 2,
) -> dict:
    """Enclose H^(0)..H^(max_moment) using the exact saddle potential."""
    if panels <= 0 or panels % 2:
        raise ValueError("Simpson panels must be a positive even integer")
    if max_moment < 2:
        raise ValueError("max_moment must be at least two")

    mode = arb_interval(left, right)
    mode_jet = potential_jet_arb(mode, 4)
    potential_zero, t, curvature = mode_jet[:3]
    if not bool(curvature > 0):
        return {"passed": False, "failure": "nonpositive-curvature"}
    x_zero = 2 * mode.log()
    sqrt_curvature = curvature.sqrt()
    y_cut = flint.arb(window_y)
    panel_width = y_cut / panels

    def potential_data(y: flint.arb, sign: int) -> list[flint.arb]:
        if bool(y == 0):
            return [
                flint.arb(0),
                flint.arb(0),
                flint.arb(1),
                sign * mode_jet[3] / curvature ** flint.arb("1.5"),
                mode_jet[4] / curvature**2,
            ]
        x = x_zero + sign * y / sqrt_curvature
        shifted_mode = (x / 2).exp()
        jet = potential_jet_arb(shifted_mode, 4)
        return [
            jet[0] - potential_zero - t * (x - x_zero),
            sign * (jet[1] - t) / sqrt_curvature,
            jet[2] / curvature,
            sign * jet[3] / curvature ** flint.arb("1.5"),
            jet[4] / curvature**2,
        ]

    branch_moments = [
        [flint.arb(0) for _ in range(max_moment + 1)] for _ in range(2)
    ]
    for node in range(panels + 1):
        y = panel_width * node
        coefficient = 1 if node in (0, panels) else (4 if node % 2 else 2)
        for branch_index, sign in enumerate((1, -1)):
            potential = potential_data(y, sign)[0]
            weight = (-potential).exp()
            for power in range(max_moment + 1):
                branch_moments[branch_index][power] += (
                    coefficient * y**power * weight
                )
    for branch_index in range(2):
        for power in range(max_moment + 1):
            branch_moments[branch_index][power] *= panel_width / 3

    simpson_errors = [
        [flint.arb(0) for _ in range(max_moment + 1)] for _ in range(2)
    ]
    pair_error_factor = (2 * panel_width) ** 5 / 2880
    for pair_start in range(0, panels, 2):
        y_left = panel_width * pair_start
        y_right = panel_width * (pair_start + 2)
        y_interval = (y_left + y_right) / 2 + flint.arb(
            0, (y_right - y_left) / 2
        )
        for branch_index, sign in enumerate((1, -1)):
            w0, w1, w2, w3, w4 = potential_data(y_interval, sign)
            weight = (-w0).exp()
            exponential_derivatives = [
                weight,
                -w1 * weight,
                (w1**2 - w2) * weight,
                (-w1**3 + 3 * w1 * w2 - w3) * weight,
                (
                    w1**4
                    - 6 * w1**2 * w2
                    + 3 * w2**2
                    + 4 * w1 * w3
                    - w4
                )
                * weight,
            ]
            for power in range(max_moment + 1):
                fourth_derivative = flint.arb(0)
                for polynomial_order in range(min(power, 4) + 1):
                    polynomial_derivative = (
                        math.factorial(power)
                        // math.factorial(power - polynomial_order)
                    ) * y_interval ** (power - polynomial_order)
                    fourth_derivative += (
                        math.comb(4, polynomial_order)
                        * polynomial_derivative
                        * exponential_derivatives[4 - polynomial_order]
                    )
                simpson_errors[branch_index][power] += (
                    upper_absolute(fourth_derivative) * pair_error_factor
                )

    for branch_index in range(2):
        for power in range(max_moment + 1):
            branch_moments[branch_index][power] += flint.arb(
                0, simpson_errors[branch_index][power]
            )

    minimum_tail_slope: flint.arb | None = None
    maximum_tail_moment: flint.arb | None = None
    for branch_index, sign in enumerate((1, -1)):
        endpoint, outward_slope = potential_data(y_cut, sign)[:2]
        if not bool(endpoint > 0 and outward_slope > 0):
            return {"passed": False, "failure": "tail-endpoint-or-slope"}
        if minimum_tail_slope is None or outward_slope.lower() < minimum_tail_slope.lower():
            minimum_tail_slope = outward_slope
        for power in range(max_moment + 1):
            tail_upper = (-exact_lower(endpoint)).exp() * tail_polynomial(
                power,
                exact_lower(outward_slope),
                y_cut,
            )
            branch_moments[branch_index][power] += positive_interval_from_upper(
                tail_upper
            )
            if maximum_tail_moment is None or tail_upper.upper() > maximum_tail_moment.upper():
                maximum_tail_moment = tail_upper

    assert minimum_tail_slope is not None and maximum_tail_moment is not None
    moments = [
        branch_moments[0][power] + (-1) ** power * branch_moments[1][power]
        for power in range(max_moment + 1)
    ]
    if not bool(moments[0] > 0):
        return {"passed": False, "failure": "nonpositive-normalizer"}
    raw = [value / moments[0] for value in moments]
    cumulants = [flint.arb(0) for _ in range(max_moment + 1)]
    cumulants[1] = raw[1]
    for order in range(2, max_moment + 1):
        cumulants[order] = raw[order] - sum(
            math.comb(order - 1, previous - 1)
            * cumulants[previous]
            * raw[order - previous]
            for previous in range(1, order)
        )

    argument = t + flint.arb(1) / 2
    log_moment = (
        flint.arb.pi().log()
        + t * x_zero
        - potential_zero
        - curvature.log() / 2
        + moments[0].log()
    )
    derivatives = {
        0: argument.lgamma() - log_moment,
        1: argument.digamma() - x_zero - raw[1] / sqrt_curvature,
    }
    for order in range(2, max_moment + 1):
        derivatives[order] = signed_hurwitz_gamma_derivative(
            order, argument
        ) - cumulants[order] / derivative_power(curvature, order)
    return {
        "passed": True,
        "failure": None,
        "mode_left": str(left),
        "mode_right": str(right),
        "t_ball": t.str(40).replace("e", "E"),
        "h_derivatives": derivatives,
        "raw_moments": {order: raw[order] for order in range(max_moment + 1)},
        "normalizer_value": moments[0],
        "normalizer_lower": arb_lower_text(moments[0]),
        "minimum_tail_slope_lower": arb_lower_text(minimum_tail_slope),
        "maximum_tail_moment_upper": arb_upper_text(maximum_tail_moment),
        "maximum_simpson_error_upper": arb_upper_text(
            max(
                (error for branch in simpson_errors for error in branch),
                key=lambda value: value.upper(),
            )
        ),
        "tail_convexity_contract": (
            "V''>0 for u>=1/100 and V'<0 for 0<u<=1/100"
        ),
    }


def integrate_h_jet_taylor_quadrature(
    left: Fraction,
    right: Fraction,
    panels: int,
    *,
    window_y: int,
    taylor_order: int = 20,
    max_moment: int = 2,
) -> dict:
    """High-order rigorous H jet from exact-potential panel Taylor series."""
    if panels <= 0:
        raise ValueError("panel count must be positive")
    if taylor_order < 4 or taylor_order % 2:
        raise ValueError("Taylor order must be an even integer at least four")
    if max_moment < 2:
        raise ValueError("max_moment must be at least two")

    mode = arb_interval(left, right)
    mode_jet = potential_jet_arb(mode, taylor_order + 1)
    potential_zero, t, curvature = mode_jet[:3]
    if not bool(curvature > 0):
        return {"passed": False, "failure": "nonpositive-curvature"}
    x_zero = 2 * mode.log()
    sqrt_curvature = curvature.sqrt()
    panel_width = flint.arb(window_y) / panels
    half_width = panel_width / 2

    def multiply_series(
        left_series: list[flint.arb],
        right_series: list[flint.arb],
        order: int,
    ) -> list[flint.arb]:
        return [
            sum(
                (
                    left_series[index] * right_series[degree - index]
                    for index in range(degree + 1)
                ),
                flint.arb(0),
            )
            for degree in range(order + 1)
        ]

    def exponential_series(
        values: list[flint.arb], order: int
    ) -> list[flint.arb]:
        result = [values[0].exp()]
        for degree in range(1, order + 1):
            result.append(
                sum(
                    (
                        index * values[index] * result[degree - index]
                        for index in range(1, degree + 1)
                    ),
                    flint.arb(0),
                )
                / degree
            )
        return result

    def potential_series(
        y: flint.arb, sign: int, order: int
    ) -> list[flint.arb]:
        if bool(y == 0):
            values = [flint.arb(0), flint.arb(0)]
            values.extend(
                sign**degree
                * mode_jet[degree]
                / derivative_power(curvature, degree)
                / math.factorial(degree)
                for degree in range(2, order + 1)
            )
            return values
        x = x_zero + sign * y / sqrt_curvature
        shifted_mode = (x / 2).exp()
        jet = potential_jet_arb(shifted_mode, order)
        values = [
            jet[0] - potential_zero - t * (x - x_zero),
            sign * (jet[1] - t) / sqrt_curvature,
        ]
        values.extend(
            sign**degree
            * jet[degree]
            / derivative_power(curvature, degree)
            / math.factorial(degree)
            for degree in range(2, order + 1)
        )
        return values

    def integrand_series(
        y: flint.arb,
        sign: int,
        power: int,
        order: int,
    ) -> list[flint.arb]:
        potential = potential_series(y, sign, order)
        weight = exponential_series([-value for value in potential], order)
        polynomial = [flint.arb(0) for _ in range(order + 1)]
        for degree in range(min(power, order) + 1):
            polynomial[degree] = (
                math.comb(power, degree) * y ** (power - degree)
            )
        return multiply_series(polynomial, weight, order)

    branch_moments = [
        [flint.arb(0) for _ in range(max_moment + 1)] for _ in range(2)
    ]
    maximum_panel_error = flint.arb(0)
    for panel in range(panels):
        midpoint = panel_width * (flint.arb(panel) + flint.arb(1) / 2)
        panel_interval = midpoint + flint.arb(0, half_width)
        for branch_index, sign in enumerate((1, -1)):
            for power in range(max_moment + 1):
                point_series = integrand_series(
                    midpoint,
                    sign,
                    power,
                    taylor_order,
                )
                integral = sum(
                    (
                        2
                        * point_series[degree]
                        * half_width ** (degree + 1)
                        / (degree + 1)
                        for degree in range(0, taylor_order + 1, 2)
                    ),
                    flint.arb(0),
                )
                remainder_series = integrand_series(
                    panel_interval,
                    sign,
                    power,
                    taylor_order + 1,
                )
                error = (
                    2
                    * upper_absolute(remainder_series[taylor_order + 1])
                    * half_width ** (taylor_order + 2)
                    / (taylor_order + 2)
                )
                branch_moments[branch_index][power] += integral + flint.arb(
                    0, error
                )
                if error.upper() > maximum_panel_error.upper():
                    maximum_panel_error = error

    y_cut = flint.arb(window_y)
    minimum_tail_slope: flint.arb | None = None
    maximum_tail_moment: flint.arb | None = None
    for branch_index, sign in enumerate((1, -1)):
        endpoint_series = potential_series(y_cut, sign, 1)
        endpoint, outward_slope = endpoint_series
        if not bool(endpoint > 0 and outward_slope > 0):
            return {"passed": False, "failure": "tail-endpoint-or-slope"}
        if minimum_tail_slope is None or outward_slope.lower() < minimum_tail_slope.lower():
            minimum_tail_slope = outward_slope
        for power in range(max_moment + 1):
            tail_upper = (-exact_lower(endpoint)).exp() * tail_polynomial(
                power,
                exact_lower(outward_slope),
                y_cut,
            )
            branch_moments[branch_index][power] += positive_interval_from_upper(
                tail_upper
            )
            if maximum_tail_moment is None or tail_upper.upper() > maximum_tail_moment.upper():
                maximum_tail_moment = tail_upper

    assert minimum_tail_slope is not None and maximum_tail_moment is not None
    moments = [
        branch_moments[0][power] + (-1) ** power * branch_moments[1][power]
        for power in range(max_moment + 1)
    ]
    if not bool(moments[0] > 0):
        return {"passed": False, "failure": "nonpositive-normalizer"}
    raw = [value / moments[0] for value in moments]
    cumulants = [flint.arb(0) for _ in range(max_moment + 1)]
    cumulants[1] = raw[1]
    for order in range(2, max_moment + 1):
        cumulants[order] = raw[order] - sum(
            math.comb(order - 1, previous - 1)
            * cumulants[previous]
            * raw[order - previous]
            for previous in range(1, order)
        )

    argument = t + flint.arb(1) / 2
    log_moment = (
        flint.arb.pi().log()
        + t * x_zero
        - potential_zero
        - curvature.log() / 2
        + moments[0].log()
    )
    derivatives = {
        0: argument.lgamma() - log_moment,
        1: argument.digamma() - x_zero - raw[1] / sqrt_curvature,
    }
    for order in range(2, max_moment + 1):
        derivatives[order] = signed_hurwitz_gamma_derivative(
            order, argument
        ) - cumulants[order] / derivative_power(curvature, order)
    return {
        "passed": True,
        "failure": None,
        "mode_left": str(left),
        "mode_right": str(right),
        "t_ball": t.str(40).replace("e", "E"),
        "h_derivatives": derivatives,
        "raw_moments": {order: raw[order] for order in range(max_moment + 1)},
        "normalizer_value": moments[0],
        "normalizer_lower": arb_lower_text(moments[0]),
        "minimum_tail_slope_lower": arb_lower_text(minimum_tail_slope),
        "maximum_tail_moment_upper": arb_upper_text(maximum_tail_moment),
        "maximum_panel_error_upper": arb_upper_text(maximum_panel_error),
        "parameters": {
            "panels": panels,
            "window_y": window_y,
            "taylor_order": taylor_order,
            "max_moment": max_moment,
        },
        "tail_convexity_contract": (
            "V''>0 for u>=1/100 and V'<0 for 0<u<=1/100"
        ),
    }


def hull(values: list[flint.arb]) -> flint.arb:
    lower = flint.arb(min(value.lower() for value in values))
    upper = flint.arb(max(value.upper() for value in values))
    return (lower + upper) / 2 + flint.arb(0, (upper - lower) / 2)


def integrate_h_derivative_cover(
    left: Fraction,
    right: Fraction,
    width: Fraction,
    panels: int,
    *,
    window_y: int,
    eighth_envelope_bound: Fraction,
    max_moment: int = DEFAULT_MAX_MOMENT,
) -> dict:
    """Partition a mode interval, certify each piece, and hull H derivatives."""
    if width <= 0:
        raise ValueError("cover width must be positive")
    blocks = []
    cursor = left
    while cursor < right:
        endpoint = min(cursor + width, right)
        result = integrate_h_derivatives(
            cursor,
            endpoint,
            panels,
            window_y=window_y,
            eighth_envelope_bound=eighth_envelope_bound,
            max_moment=max_moment,
        )
        if not result.get("passed"):
            return {
                "passed": False,
                "failure": "h-cover-subblock",
                "failed_block": [str(cursor), str(endpoint)],
                "subblock_result": result,
            }
        blocks.append(result)
        cursor = endpoint
    h_derivatives = {
        order: hull([block["h_derivatives"][order] for block in blocks])
        for order in range(2, max_moment + 1)
    }
    return {
        "passed": True,
        "failure": None,
        "block_count": len(blocks),
        "width": str(width),
        "panels_per_block": panels,
        "mode_interval": [str(left), str(right)],
        "h_derivatives": h_derivatives,
        "minimum_normalizer_lower": min(
            (block["normalizer_lower"] for block in blocks), key=Decimal
        ),
        "minimum_tail_slope_lower": min(
            (block["minimum_tail_slope_lower"] for block in blocks), key=Decimal
        ),
        "maximum_tail_moment_upper": max(
            (block["maximum_tail_moment_upper"] for block in blocks), key=Decimal
        ),
        "maximum_simpson_error_upper": max(
            (block["maximum_simpson_error_upper"] for block in blocks), key=Decimal
        ),
        "selected_blocks": [
            {
                key: block[key]
                for key in (
                    "mode_left",
                    "mode_right",
                    "t_ball",
                    "curvature_ball",
                    "normalizer_lower",
                    "maximum_simpson_error_upper",
                )
            }
            for block in (blocks[0], blocks[len(blocks) // 2], blocks[-1])
        ],
    }


def ell_derivative_jet(b_derivatives: list[flint.arb], order: int = 6) -> list[flint.arb]:
    """Differentiate ell=log(1-exp(-B)) with truncated power series."""
    if len(b_derivatives) < order + 1:
        raise ValueError("insufficient B derivatives")
    b_series = [
        b_derivatives[index] / math.factorial(index) for index in range(order + 1)
    ]
    negative = [-value for value in b_series]
    exponential = [flint.arb(0) for _ in range(order + 1)]
    exponential[0] = negative[0].exp()
    for degree in range(1, order + 1):
        exponential[degree] = sum(
            index * negative[index] * exponential[degree - index]
            for index in range(1, degree + 1)
        ) / degree
    defect = [flint.arb(0) for _ in range(order + 1)]
    defect[0] = 1 - exponential[0]
    for degree in range(1, order + 1):
        defect[degree] = -exponential[degree]
    if not bool(defect[0] > 0):
        raise RuntimeError("defect interval is not positive")
    logarithm = [flint.arb(0) for _ in range(order + 1)]
    logarithm[0] = defect[0].log()
    for degree in range(1, order + 1):
        logarithm[degree] = (
            degree * defect[degree]
            - sum(
                index * logarithm[index] * defect[degree - index]
                for index in range(1, degree)
            )
        ) / (degree * defect[0])
    return [
        logarithm[index] * math.factorial(index) for index in range(order + 1)
    ]


def lower_absolute(value: flint.arb) -> flint.arb:
    if value.contains(0):
        return flint.arb(0)
    return flint.arb(min(abs(value.lower()), abs(value.upper())))


def positive_upper(value: flint.arb) -> flint.arb:
    return flint.arb(max(0, value.upper()))


def evaluate_localized_curvature_from_h_cover(
    left: Fraction,
    right: Fraction,
    h_derivatives: dict[int, flint.arb],
    *,
    cover_diagnostics: dict | None = None,
) -> dict:
    """Evaluate U(t)<=7/(2t^2) from a supplied H-derivative cover."""
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    central_mode = arb_interval(left, right)
    central_t = potential_jet_arb(central_mode, 1)[1]
    missing_orders = [order for order in range(2, 9) if order not in h_derivatives]
    if missing_orders:
        return {
            "passed": False,
            "failure": "missing-h-derivatives",
            "missing_orders": missing_orders,
        }
    b_derivatives = [h_derivatives[index + 2] for index in range(7)]
    ell = ell_derivative_jet(b_derivatives, 6)
    envelopes = [upper_absolute(ell[index + 4]) / 12 for index in range(3)]
    local_jets = [2 * b_derivatives[index] - ell[index + 2] for index in range(3)]

    j_lower = exact_lower(local_jets[0] - envelopes[0])
    j_upper = exact_upper(local_jets[0] + envelopes[0])
    if not bool(j_lower > 0):
        return {
            "passed": False,
            "failure": "nonpositive-localized-j",
            "j_lower": j_lower.str(40),
        }
    q_upper = positive_upper(local_jets[2] + envelopes[2])
    p_lower = lower_absolute(local_jets[1]) - exact_upper(envelopes[1])
    if p_lower < 0:
        p_lower = flint.arb(0)
    phi_upper = 1 / (j_lower.exp() - 1)
    chi_lower = j_upper.exp() / (j_upper.exp() - 1) ** 2
    localized_upper = (
        2 * exact_upper(ell[2])
        + phi_upper * q_upper
        - chi_lower * p_lower**2
    )
    target_lower = 7 / (2 * exact_upper(central_t) ** 2)
    margin = target_lower - localized_upper
    passed = bool(margin > 0)
    return {
        "passed": passed,
        "failure": None if passed else "nonpositive-localized-curvature-margin",
        "central_mode": [str(left), str(right)],
        "central_t_ball": central_t.str(40).replace("e", "E"),
        "j_lower": arb_lower_text(j_lower),
        "localized_upper": arb_upper_text(localized_upper),
        "target_lower": arb_lower_text(target_lower),
        "margin_lower": arb_lower_text(margin),
        "scaled_localized_upper": arb_upper_text(
            localized_upper * exact_upper(central_t) ** 2
        ),
        "ell_second_upper": arb_upper_text(ell[2]),
        "local_jets": [value.str(40).replace("e", "E") for value in local_jets],
        "envelopes": [value.str(40).replace("e", "E") for value in envelopes],
        "h_derivatives": {
            str(order): value.str(40).replace("e", "E")
            for order, value in h_derivatives.items()
        },
        "cover_diagnostics": cover_diagnostics or {},
    }


def certify_localized_curvature_block(
    left: Fraction,
    right: Fraction,
    outer_left: Fraction,
    outer_right: Fraction,
    panels: int,
    *,
    outer_width: Fraction = Fraction(1, 10_000),
    window_y: int = 6,
    eighth_envelope_bound: Fraction = Fraction(1, 50_000),
) -> dict:
    """Certify U(t)<=7/(2t^2) on one central mode block."""
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    lower_t = potential_jet_arb(arb_rational(left), 1)[1]
    upper_t = potential_jet_arb(arb_rational(right), 1)[1]
    outer_lower_t = potential_jet_arb(arb_rational(outer_left), 1)[1]
    outer_upper_t = potential_jet_arb(arb_rational(outer_right), 1)[1]
    if not bool(outer_lower_t < lower_t - 2 and outer_upper_t > upper_t + 2):
        return {"passed": False, "failure": "outer-mode-cover"}

    integrated = integrate_h_derivative_cover(
        outer_left,
        outer_right,
        outer_width,
        panels,
        window_y=window_y,
        eighth_envelope_bound=eighth_envelope_bound,
        max_moment=8,
    )
    if not integrated.get("passed"):
        return integrated
    result = evaluate_localized_curvature_from_h_cover(
        left,
        right,
        integrated["h_derivatives"],
        cover_diagnostics={
            key: value
            for key, value in integrated.items()
            if key not in {"h_derivatives", "passed", "failure"}
        },
    )
    result.update(
        {
            "outer_mode": [str(outer_left), str(outer_right)],
            "panels": panels,
            "outer_width": str(outer_width),
            "window_y": window_y,
            "eighth_envelope": str(eighth_envelope_bound),
            "quadrature": result.pop("cover_diagnostics"),
        }
    )
    return result


def parse_fraction(value: str) -> Fraction:
    return Fraction(value)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--left", type=parse_fraction, required=True)
    parser.add_argument("--right", type=parse_fraction, required=True)
    parser.add_argument("--outer-left", type=parse_fraction, required=True)
    parser.add_argument("--outer-right", type=parse_fraction, required=True)
    parser.add_argument("--panels", type=int, default=1000)
    parser.add_argument("--outer-width", type=parse_fraction, default=Fraction(1, 10_000))
    parser.add_argument("--window-y", type=int, default=6)
    parser.add_argument("--eighth-envelope", type=parse_fraction, default=Fraction(1, 50_000))
    args = parser.parse_args()
    result = certify_localized_curvature_block(
        args.left,
        args.right,
        args.outer_left,
        args.outer_right,
        args.panels,
        outer_width=args.outer_width,
        window_y=args.window_y,
        eighth_envelope_bound=args.eighth_envelope,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result.get("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
