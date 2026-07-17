#!/usr/bin/env python3
"""Certify the paired seventh-order remainder on a large compact mode range."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md"
)
DEFAULT_PRECISION_BITS = 192
DEFAULT_MODE_START = Fraction(579, 625)  # 0.9264
DEFAULT_MODE_END = Fraction(5, 1)
DEFAULT_BASE_WIDTH = Fraction(1, 1000)
DEFAULT_MAX_REFINEMENT_DEPTH = 4
DEFAULT_PRIMARY_PANELS = 300
DEFAULT_MID_PANELS = 600
DEFAULT_RETRY_PANELS = 900
DEFAULT_HIGH_PANEL_START = Fraction(19, 5)  # 3.8
DEFAULT_WINDOW_Y = 6
DEFAULT_EIGHTH_ENVELOPE = Fraction(1, 50_000)
DEFAULT_REMAINDER_FLOOR = Fraction(-79, 1000)
DEFAULT_AUXILIARY_STEP = Fraction(1, 10_000)


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str | None
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None
    gap: str | None = None


def arb_rational(value: Fraction) -> flint.arb:
    return flint.arb(value.numerator) / value.denominator


def arb_interval(left: Fraction, right: Fraction) -> flint.arb:
    left_arb = arb_rational(left)
    right_arb = arb_rational(right)
    return (left_arb + right_arb) / 2 + flint.arb(0, (right_arb - left_arb) / 2)


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def arb_upper_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.upper().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_plus(), "E")


def exact_upper(value: flint.arb) -> flint.arb:
    return flint.arb(value.upper())


def exact_lower(value: flint.arb) -> flint.arb:
    return flint.arb(value.lower())


def upper_absolute(value: flint.arb) -> flint.arb:
    return exact_upper(abs(value))


def scaled_curvature_quantities(
    t: flint.arb,
    curvature: flint.arb,
    raw1: flint.arb,
    raw2: flint.arb,
    raw3: flint.arb,
) -> dict[str, flint.arb]:
    variance = raw2 - raw1**2
    cumulant3 = raw3 - 3 * raw1 * raw2 + 2 * raw1**3
    x_variance = variance / curvature
    x_cumulant3 = cumulant3 / curvature ** flint.arb("1.5")
    gamma_argument = t + flint.arb(1) / 2
    gamma_second = flint.arb(2).zeta(gamma_argument)
    gamma_third = -2 * flint.arb(3).zeta(gamma_argument)
    h_second = gamma_second - x_variance
    h_third = gamma_third - x_cumulant3
    pointwise_derivative = 2 * h_second + (2 * t + 1) * h_third
    shifted_buffer = pointwise_derivative - 2 * upper_absolute(h_third)
    return {
        "variance": variance,
        "cumulant3": cumulant3,
        "h_second": h_second,
        "h_third": h_third,
        "pointwise_derivative": pointwise_derivative,
        "shifted_buffer": shifted_buffer,
    }


def fraction_grid(start: Fraction, end: Fraction, width: Fraction) -> list[tuple[Fraction, Fraction]]:
    rows = []
    left = start
    while left < end:
        right = min(left + width, end)
        rows.append((left, right))
        left = right
    return rows


def auxiliary_geometry_certificate() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    curvature_start = Fraction(1, 100)
    curvature_end = DEFAULT_MODE_START
    curvature_rows = fraction_grid(curvature_start, curvature_end, DEFAULT_AUXILIARY_STEP)
    positive_curvature = 0
    minimum_curvature: tuple[flint.arb, int] | None = None
    for index, (left, right) in enumerate(curvature_rows):
        value = potential_jet_arb(arb_interval(left, right), 2)[2]
        positive_curvature += int(bool(value > 0))
        if minimum_curvature is None or value.lower() < minimum_curvature[0].lower():
            minimum_curvature = (value, index)
    if positive_curvature != len(curvature_rows):
        raise RuntimeError("left curvature certificate failed")
    assert minimum_curvature is not None

    tiny = arb_rational(curvature_start)
    tiny_jet = potential_jet_arb(tiny, 1)
    tiny_slope = tiny_jet[1]
    tiny_q = flint.arb.pi() * (4 * tiny).exp()
    tiny_crude_upper = 100 * tiny**2 + 2 * tiny * tiny_q - flint.arb(1) / 2
    left_reference = arb_rational(Fraction(41, 50))
    left_reference_slope = potential_jet_arb(left_reference, 1)[1]
    if not bool(
        tiny_q < flint.arb(33) / 10
        and tiny_crude_upper < 0
        and tiny_slope < 0
        and left_reference_slope > 0
    ):
        raise RuntimeError("tiny-u potential slope is not negative")

    v9_start = Fraction(41, 50)  # 0.82
    v9_end = Fraction(51, 10)  # 5.1
    v9_rows = fraction_grid(v9_start, v9_end, DEFAULT_AUXILIARY_STEP)
    positive_v9 = 0
    minimum_v9: tuple[flint.arb, int] | None = None
    for index, (left, right) in enumerate(v9_rows):
        value = potential_jet_arb(arb_interval(left, right), 9)[9]
        positive_v9 += int(bool(value > 0))
        if minimum_v9 is None or value.lower() < minimum_v9[0].lower():
            minimum_v9 = (value, index)
    if positive_v9 != len(v9_rows):
        raise RuntimeError("V9 monotonicity certificate failed")
    assert minimum_v9 is not None

    ray_q = flint.arb.pi() * flint.arb(20).exp()
    ray_v9_margin = flint.arb(512) * 1_000_000_000 - 97_312_467_060
    if not bool(ray_q > 1_000_000_000 and ray_v9_margin > 0):
        raise RuntimeError("V9 ray gate failed")
    return {
        "positive_curvature_intervals": positive_curvature,
        "curvature_interval_count": len(curvature_rows),
        "minimum_curvature_lower": arb_lower_text(minimum_curvature[0]),
        "minimum_curvature_interval": minimum_curvature[1],
        "tiny_u_slope_upper": arb_upper_text(tiny_slope),
        "tiny_u_q_upper": arb_upper_text(tiny_q),
        "tiny_u_crude_slope_upper": arb_upper_text(tiny_crude_upper),
        "left_reference_u": "41/50",
        "left_reference_slope_lower": arb_lower_text(left_reference_slope),
        "tiny_u_argument": (
            "For 0<u<=1/100, q<3.3 gives V'(x)<0; for u>=1/100 the "
            "certified V''>0 makes V' increasing."
        ),
        "positive_v9_intervals": positive_v9,
        "v9_interval_count": len(v9_rows),
        "minimum_v9_lower": arb_lower_text(minimum_v9[0]),
        "minimum_v9_interval": minimum_v9[1],
        "v9_ray_q_lower": arb_lower_text(ray_q),
        "v9_ray_margin_lower": arb_lower_text(ray_v9_margin),
        "v9_ray_argument": (
            "For u>=5, q_9>=512*u^9*q and the absolute log-composition "
            "correction is at most 97312467060*u^9, so V^(9)>0."
        ),
    }


def eighth_derivative_envelope() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    rows = fraction_grid(DEFAULT_MODE_START, DEFAULT_MODE_END, DEFAULT_AUXILIARY_STEP)
    target = arb_rational(DEFAULT_EIGHTH_ENVELOPE)
    positive_rows = 0
    minimum_margin: tuple[flint.arb, int] | None = None
    selected = []
    minimum_window_u: flint.arb | None = None
    maximum_window_u: flint.arb | None = None
    selected_indices = {0, len(rows) // 4, len(rows) // 2, 3 * len(rows) // 4, len(rows) - 1}
    for index, (left, right) in enumerate(rows):
        mode = arb_interval(left, right)
        curvature = potential_jet_arb(mode, 2)[2]
        right_window = mode * (flint.arb(3) / curvature.sqrt()).exp()
        left_window = mode * (-flint.arb(3) / curvature.sqrt()).exp()
        if minimum_window_u is None or left_window.lower() < minimum_window_u.lower():
            minimum_window_u = left_window
        if maximum_window_u is None or right_window.upper() > maximum_window_u.upper():
            maximum_window_u = right_window
        right_point = exact_upper(right_window)
        eighth = potential_jet_arb(right_point, 8)[8]
        ratio = eighth / curvature**4
        margin = target - ratio
        positive_rows += int(bool(margin > 0))
        if minimum_margin is None or margin.lower() < minimum_margin[0].lower():
            minimum_margin = (margin, index)
        if index in selected_indices:
            selected.append(
                {
                    "index": index,
                    "mode_left": str(left),
                    "mode_right": str(right),
                    "right_window_upper": arb_upper_text(right_window, 35),
                    "normalized_eighth_upper": arb_upper_text(ratio, 35),
                    "margin_lower": arb_lower_text(margin, 35),
                }
            )
    if positive_rows != len(rows):
        raise RuntimeError("normalized eighth-derivative envelope failed")
    assert minimum_margin is not None and minimum_window_u is not None and maximum_window_u is not None
    if not bool(minimum_window_u > flint.arb(41) / 50):
        raise RuntimeError("central window leaves the V9-positive range")

    ray_margin = Fraction(1, 50_000) - Fraction(1_000 * 2**9, 1) / (
        Fraction(39, 10) ** 4 * 1_000_000_000**3
    )
    if ray_margin <= 0:
        raise RuntimeError("normalized eighth-derivative ray failed")
    return {
        "mode_interval": "0.9264<=u<=5",
        "interval_count": len(rows),
        "positive_interval_count": positive_rows,
        "minimum_margin_lower": arb_lower_text(minimum_margin[0]),
        "minimum_margin_interval": minimum_margin[1],
        "selected_rows": selected,
        "minimum_window_u_lower": arb_lower_text(minimum_window_u),
        "maximum_window_u_upper": arb_upper_text(maximum_window_u),
        "ray_margin_fraction": str(ray_margin),
        "ray_argument": (
            "V^(9)>0 puts the central-window maximum of V^(8) at y=+6. "
            "For u>=5, a>=(39/10)u^2q, q>=10^9, the window has u_+/u<1.00002, "
            "and V^(8)(u_+)<1000*u^8*q, which is far below a^4/50000."
        ),
        "certified_bound": "sup_|y|<=6 V^(8)(x_t+y/sqrt(a_t))/a_t^4<=1/50000",
    }


def tail_polynomial(power: int, slope: flint.arb, y_cut: flint.arb) -> flint.arb:
    return sum(
        flint.arb(math.comb(power, degree) * math.factorial(degree))
        * y_cut ** (power - degree)
        / slope ** (degree + 1)
        for degree in range(power + 1)
    )


def positive_interval_from_upper(upper: flint.arb) -> flint.arb:
    bound = exact_upper(upper)
    return bound / 2 + flint.arb(0, bound / 2)


def certify_mode_block(
    left: Fraction,
    right: Fraction,
    panels: int,
    *,
    window_y: int = DEFAULT_WINDOW_Y,
    eighth_envelope_bound: Fraction = DEFAULT_EIGHTH_ENVELOPE,
) -> dict:
    mode = arb_interval(left, right)
    jet = potential_jet_arb(mode, 7)
    t, curvature = jet[1], jet[2]
    if not bool(curvature > 0):
        return {"passed": False, "failure": "nonpositive-curvature"}
    y_cut = flint.arb(window_y)
    panel_width = y_cut / panels
    eighth_envelope = arb_rational(eighth_envelope_bound)
    normalized = [flint.arb(0) for _ in range(8)]
    # Build the powers explicitly: python-flint does not accept Fraction exponents.
    normalized[2] = flint.arb(1)
    normalized[3] = jet[3] / curvature ** flint.arb("1.5")
    normalized[4] = jet[4] / curvature**2
    normalized[5] = jet[5] / curvature ** flint.arb("2.5")
    normalized[6] = jet[6] / curvature**3
    normalized[7] = jet[7] / curvature ** flint.arb("3.5")

    moments = [flint.arb(0) for _ in range(4)]
    for panel in range(panels):
        y_left = panel_width * panel
        y_right = panel_width * (panel + 1)
        y_mid = (y_left + y_right) / 2
        y_interval = y_mid + flint.arb(0, (y_right - y_left) / 2)
        midpoint_weights: list[flint.arb] = []
        second_errors = [[flint.arb(0) for _ in range(4)] for _ in range(2)]
        for branch_index, sign in enumerate((1, -1)):
            midpoint_potential = y_mid**2 / 2
            interval_potential = y_interval**2 / 2
            first_derivative = y_interval
            second_derivative = flint.arb(1)
            for order in range(3, 8):
                parity = sign**order
                midpoint_potential += (
                    normalized[order] * parity * y_mid**order / math.factorial(order)
                )
                interval_potential += (
                    normalized[order] * parity * y_interval**order / math.factorial(order)
                )
                first_derivative += (
                    normalized[order]
                    * parity
                    * y_interval ** (order - 1)
                    / math.factorial(order - 1)
                )
                second_derivative += (
                    normalized[order]
                    * parity
                    * y_interval ** (order - 2)
                    / math.factorial(order - 2)
                )
            value_radius = (
                eighth_envelope * y_right**8 / math.factorial(8)
            )
            first_radius = (
                eighth_envelope * y_right**7 / math.factorial(7)
            )
            second_radius = (
                eighth_envelope * y_right**6 / math.factorial(6)
            )
            midpoint_weight = (
                -midpoint_potential + flint.arb(0, value_radius)
            ).exp()
            interval_weight = (
                -interval_potential + flint.arb(0, value_radius)
            ).exp()
            first_derivative += flint.arb(0, first_radius)
            second_derivative += flint.arb(0, second_radius)
            midpoint_weights.append(midpoint_weight)
            for power in range(4):
                second = (
                    y_interval**power
                    * (first_derivative**2 - second_derivative)
                    * interval_weight
                )
                if power >= 1:
                    second -= (
                        2
                        * power
                        * y_interval ** (power - 1)
                        * first_derivative
                        * interval_weight
                    )
                if power >= 2:
                    second += (
                        power
                        * (power - 1)
                        * y_interval ** (power - 2)
                        * interval_weight
                    )
                second_errors[branch_index][power] = (
                    upper_absolute(second) * panel_width**3 / 24
                )
        plus, minus = midpoint_weights
        for power in range(4):
            pair = plus + minus if power % 2 == 0 else plus - minus
            quadrature_error = second_errors[0][power] + second_errors[1][power]
            moments[power] += (
                y_mid**power * pair * panel_width + flint.arb(0, quadrature_error)
            )

    tail_uppers = [[flint.arb(0) for _ in range(4)] for _ in range(2)]
    minimum_tail_slope: flint.arb | None = None
    maximum_tail_third: flint.arb | None = None
    for branch_index, sign in enumerate((1, -1)):
        endpoint = y_cut**2 / 2
        outward_slope = y_cut
        for order in range(3, 8):
            parity = sign**order
            endpoint += (
                normalized[order] * parity * y_cut**order / math.factorial(order)
            )
            outward_slope += (
                normalized[order]
                * parity
                * y_cut ** (order - 1)
                / math.factorial(order - 1)
            )
        endpoint -= eighth_envelope * y_cut**8 / math.factorial(8)
        outward_slope -= eighth_envelope * y_cut**7 / math.factorial(7)
        if not bool(endpoint > 0 and outward_slope > 0):
            return {"passed": False, "failure": "tail-endpoint-or-slope"}
        if minimum_tail_slope is None or outward_slope.lower() < minimum_tail_slope.lower():
            minimum_tail_slope = outward_slope
        for power in range(4):
            upper = (-exact_lower(endpoint)).exp() * tail_polynomial(
                power, exact_lower(outward_slope), y_cut
            )
            tail_uppers[branch_index][power] = upper
            if power == 3 and (
                maximum_tail_third is None or upper.upper() > maximum_tail_third.upper()
            ):
                maximum_tail_third = upper
    assert minimum_tail_slope is not None and maximum_tail_third is not None
    for power in range(4):
        total_upper = tail_uppers[0][power] + tail_uppers[1][power]
        if power % 2 == 0:
            moments[power] += positive_interval_from_upper(total_upper)
        else:
            moments[power] += flint.arb(0, exact_upper(total_upper))

    if not bool(moments[0] > 0):
        return {"passed": False, "failure": "nonpositive-normalizer"}
    raw1 = moments[1] / moments[0]
    raw2 = moments[2] / moments[0]
    raw3 = moments[3] / moments[0]
    curvature_quantities = scaled_curvature_quantities(
        t, curvature, raw1, raw2, raw3
    )
    variance = curvature_quantities["variance"]
    cumulant3 = curvature_quantities["cumulant3"]
    h_second = curvature_quantities["h_second"]
    h_third = curvature_quantities["h_third"]
    pointwise_curvature_derivative = curvature_quantities["pointwise_derivative"]
    shifted_curvature_buffer = curvature_quantities["shifted_buffer"]

    alpha = normalized[3]
    beta = normalized[4]
    gamma = normalized[5]
    delta = normalized[6]
    epsilon = normalized[7]
    cubic = (8 * alpha**3 - 7 * alpha * beta + gamma) / 2
    fifth = (
        525 * alpha**5
        - 954 * alpha**3 * beta
        + 234 * alpha**2 * gamma
        + 298 * alpha * beta**2
        - 37 * alpha * delta
        - 57 * beta * gamma
        + 3 * epsilon
    ) / 24
    prefactor = t**2 / curvature ** flint.arb("1.5")
    remainder = prefactor * (cumulant3 + alpha + cubic + fifth)
    margin = remainder - arb_rational(DEFAULT_REMAINDER_FLOOR)
    passed_margin = bool(margin > 0)
    negative_cumulant = bool(cumulant3 < 0)
    failure = None
    if not passed_margin:
        failure = "nonpositive-remainder-margin"
    elif not negative_cumulant:
        failure = "nonnegative-or-inconclusive-third-cumulant"
    return {
        "passed": passed_margin and negative_cumulant,
        "failure": failure,
        "panels": panels,
        "cumulant3_ball": cumulant3.str(40).replace("e", "E"),
        "cumulant3_upper": arb_upper_text(cumulant3),
        "cumulant3_negative": negative_cumulant,
        "variance_positive": bool(variance > 0),
        "h_second_positive": bool(h_second > 0),
        "h_third_negative": bool(h_third < 0),
        "pointwise_curvature_derivative_lower": arb_lower_text(
            pointwise_curvature_derivative
        ),
        "shifted_curvature_buffer_positive": bool(shifted_curvature_buffer > 0),
        "shifted_curvature_buffer_lower": arb_lower_text(shifted_curvature_buffer),
        "margin_lower": arb_lower_text(margin),
        "remainder_ball": remainder.str(40).replace("e", "E"),
        "normalizer_lower": arb_lower_text(moments[0]),
        "minimum_tail_slope_lower": arb_lower_text(minimum_tail_slope),
        "maximum_tail_third_upper": arb_upper_text(maximum_tail_third),
    }


def certify_scaled_curvature_mode_block(
    left: Fraction,
    right: Fraction,
    panels: int,
    *,
    window_y: int = 8,
    eighth_envelope_bound: Fraction = Fraction(1, 10_000_000_000),
) -> dict:
    """Certify the buffered pointwise curvature derivative with Simpson error."""
    if panels <= 0 or panels % 2:
        raise ValueError("Simpson panels must be a positive even integer")

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
    normalized[3] = jet[3] / curvature ** flint.arb("1.5")
    normalized[4] = jet[4] / curvature**2
    normalized[5] = jet[5] / curvature ** flint.arb("2.5")
    normalized[6] = jet[6] / curvature**3
    normalized[7] = jet[7] / curvature ** flint.arb("3.5")

    branch_moments = [[flint.arb(0) for _ in range(4)] for _ in range(2)]
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
            for power in range(4):
                branch_moments[branch_index][power] += (
                    coefficient * y**power * weight
                )
    for branch_index in range(2):
        for power in range(4):
            branch_moments[branch_index][power] *= panel_width / 3

    simpson_errors = [[flint.arb(0) for _ in range(4)] for _ in range(2)]
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
                radius = (
                    eighth_envelope
                    * y_right ** (8 - derivative_order)
                    / math.factorial(8 - derivative_order)
                )
                derivatives[derivative_order] += flint.arb(0, radius)

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
            for power in range(4):
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
        for power in range(4):
            branch_moments[branch_index][power] += flint.arb(
                0, simpson_errors[branch_index][power]
            )

    minimum_tail_slope: flint.arb | None = None
    maximum_tail_third: flint.arb | None = None
    for branch_index, sign in enumerate((1, -1)):
        endpoint = y_cut**2 / 2
        outward_slope = y_cut
        for order in range(3, 8):
            endpoint += (
                normalized[order] * sign**order * y_cut**order / math.factorial(order)
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
        for power in range(4):
            tail_upper = (-exact_lower(endpoint)).exp() * tail_polynomial(
                power, exact_lower(outward_slope), y_cut
            )
            branch_moments[branch_index][power] += positive_interval_from_upper(
                tail_upper
            )
            if power == 3 and (
                maximum_tail_third is None
                or tail_upper.upper() > maximum_tail_third.upper()
            ):
                maximum_tail_third = tail_upper

    assert minimum_tail_slope is not None and maximum_tail_third is not None
    moments = [flint.arb(0) for _ in range(4)]
    for power in range(4):
        if power % 2 == 0:
            moments[power] = branch_moments[0][power] + branch_moments[1][power]
        else:
            moments[power] = branch_moments[0][power] - branch_moments[1][power]
    if not bool(moments[0] > 0):
        return {"passed": False, "failure": "nonpositive-normalizer"}

    raw1 = moments[1] / moments[0]
    raw2 = moments[2] / moments[0]
    raw3 = moments[3] / moments[0]
    quantities = scaled_curvature_quantities(t, curvature, raw1, raw2, raw3)
    transfer_margin = quantities["shifted_buffer"] - 64 / (t - 3) ** 5
    passed = bool(
        quantities["variance"] > 0
        and quantities["h_second"] > 0
        and quantities["h_third"] < 0
        and quantities["shifted_buffer"] > 0
        and transfer_margin > 0
    )
    return {
        "passed": passed,
        "failure": None if passed else "nonpositive-scaled-curvature-buffer",
        "panels": panels,
        "window_y": window_y,
        "eighth_envelope_bound": str(eighth_envelope_bound),
        "variance_positive": bool(quantities["variance"] > 0),
        "h_second_positive": bool(quantities["h_second"] > 0),
        "h_third_negative": bool(quantities["h_third"] < 0),
        "pointwise_curvature_derivative_lower": arb_lower_text(
            quantities["pointwise_derivative"]
        ),
        "shifted_curvature_buffer_lower": arb_lower_text(
            quantities["shifted_buffer"]
        ),
        "full_kernel_transfer_margin_lower": arb_lower_text(transfer_margin),
        "normalizer_lower": arb_lower_text(moments[0]),
        "maximum_simpson_error_upper": arb_upper_text(
            max(
                (error for branch in simpson_errors for error in branch),
                key=lambda value: value.upper(),
            )
        ),
        "minimum_tail_slope_lower": arb_lower_text(minimum_tail_slope),
        "maximum_tail_third_upper": arb_upper_text(maximum_tail_third),
    }


def compact_remainder_certificate() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    base_rows = fraction_grid(DEFAULT_MODE_START, DEFAULT_MODE_END, DEFAULT_BASE_WIDTH)
    accepted: list[dict] = []
    attempts = 0
    maximum_depth = 0
    start_time = perf_counter()

    def cover(left: Fraction, right: Fraction, depth: int) -> None:
        nonlocal attempts, maximum_depth
        maximum_depth = max(maximum_depth, depth)
        primary_panels = (
            DEFAULT_MID_PANELS if left >= DEFAULT_HIGH_PANEL_START else DEFAULT_PRIMARY_PANELS
        )
        retry_panels = (
            DEFAULT_RETRY_PANELS if left >= DEFAULT_HIGH_PANEL_START else DEFAULT_MID_PANELS
        )
        result = certify_mode_block(left, right, primary_panels)
        attempts += 1
        if not result["passed"]:
            result = certify_mode_block(left, right, retry_panels)
            attempts += 1
        if result["passed"]:
            accepted.append(
                {
                    "left": str(left),
                    "right": str(right),
                    "depth": depth,
                    **result,
                }
            )
            return
        if depth >= DEFAULT_MAX_REFINEMENT_DEPTH:
            raise RuntimeError(
                f"paired remainder block failed at {left}..{right}: {result['failure']}"
            )
        middle = (left + right) / 2
        cover(left, middle, depth + 1)
        cover(middle, right, depth + 1)

    for index, (left, right) in enumerate(base_rows):
        cover(left, right, 0)
        if (index + 1) % 500 == 0:
            print(
                f"paired remainder progress: {index + 1}/{len(base_rows)} base blocks, "
                f"{len(accepted)} accepted"
            )

    worst = min(accepted, key=lambda row: Decimal(row["margin_lower"]))
    closest_cumulant = max(accepted, key=lambda row: Decimal(row["cumulant3_upper"]))
    selected_indices = {0, len(accepted) // 4, len(accepted) // 2, 3 * len(accepted) // 4, len(accepted) - 1}
    selected = [accepted[index] for index in sorted(selected_indices)]
    covered_left = Fraction(accepted[0]["left"])
    covered_right = Fraction(accepted[-1]["right"])
    if covered_left != DEFAULT_MODE_START or covered_right != DEFAULT_MODE_END:
        raise RuntimeError("adaptive compact cover endpoints do not match")
    for previous, current in zip(accepted, accepted[1:]):
        if Fraction(previous["right"]) != Fraction(current["left"]):
            raise RuntimeError("adaptive compact cover has a gap")
    end_t = potential_jet_arb(arb_rational(DEFAULT_MODE_END), 1)[1]
    return {
        "mode_interval": "0.9264<=u<=5",
        "base_block_count": len(base_rows),
        "accepted_block_count": len(accepted),
        "attempt_count": attempts,
        "maximum_refinement_depth": maximum_depth,
        "primary_panels": DEFAULT_PRIMARY_PANELS,
        "mid_panels": DEFAULT_MID_PANELS,
        "retry_panels": DEFAULT_RETRY_PANELS,
        "high_panel_start": str(DEFAULT_HIGH_PANEL_START),
        "window": "|y|<=6 plus two analytic tails",
        "minimum_margin_lower": worst["margin_lower"],
        "minimum_margin_block": {key: worst[key] for key in ("left", "right", "depth", "panels")},
        "negative_cumulant_blocks": sum(row["cumulant3_negative"] for row in accepted),
        "maximum_cumulant3_upper": closest_cumulant["cumulant3_upper"],
        "maximum_cumulant3_block": {
            key: closest_cumulant[key] for key in ("left", "right", "depth", "panels")
        },
        "selected_blocks": selected,
        "elapsed_seconds": perf_counter() - start_time,
        "end_t_ball": end_t.str(50).replace("e", "E"),
        "end_t_lower": arb_lower_text(end_t),
        "open_ray": "u>=5, equivalently t>=V'(x(5)), remains to be proved analytically",
    }


def build_artifact() -> dict:
    geometry = auxiliary_geometry_certificate()
    eighth = eighth_derivative_envelope()
    compact = compact_remainder_certificate()
    rows = [
        CertificateRow(
            id="fsprc_01_standardized_pairing",
            role="exact_reduction",
            claim="Pair y and -y in standardized mode coordinates and retain normalization and centralization exactly.",
            formula=(
                "p_t(y) proportional exp(-W_t(y)); H_t=(t^2/a_t^(3/2))*"
                "(kappa_3(Y)+alpha+C3+C5)"
            ),
            readiness="available_exact",
            proof_boundary="Exact coordinate and cumulant reduction only.",
        ),
        CertificateRow(
            id="fsprc_02_left_tail_monotonicity",
            role="interval_analytic_lemma",
            claim="The left outward slope at y=-6 controls the entire farther-left tail.",
            formula="V''>0 for u>=1/100 and V'<0 below that threshold",
            readiness="interval_validated",
            proof_boundary="Potential-geometry theorem only.",
            diagnostics=geometry,
        ),
        CertificateRow(
            id="fsprc_03_eighth_derivative_envelope",
            role="interval_theorem",
            claim="The seventh-order Taylor model has a global central-window eighth-derivative envelope.",
            formula="sup_|y|<=6 V^(8)(x_t+y/sqrt(a_t))/a_t^4<=1/50000",
            readiness="interval_validated",
            proof_boundary="Central Taylor-remainder theorem only.",
            diagnostics=eighth,
        ),
        CertificateRow(
            id="fsprc_04_paired_midpoint_rule",
            role="interval_certificate",
            claim="A second-derivative midpoint rule encloses all four paired raw moments on |y|<=6.",
            formula="integral_panel f in h*f(mid)+[-sup|f''|*h^3/24,+sup|f''|*h^3/24]",
            readiness="interval_validated",
            proof_boundary="Finite adaptive compact quadrature only.",
        ),
        CertificateRow(
            id="fsprc_05_two_tail_budget",
            role="exact_interval_bound",
            claim="Endpoint values and outward slopes give explicit exponential moment tails beyond y=+/-6.",
            formula="W(6+z)>=W(6)+s*z; integrate (6+z)^j exp(-s*z) exactly",
            readiness="interval_validated",
            proof_boundary="Tail enclosure only.",
        ),
        CertificateRow(
            id="fsprc_06_compact_remainder_theorem",
            role="interval_theorem",
            claim="The seventh-order normalized remainder floor and strict negative third-cumulant sign hold throughout the compact mode range.",
            formula="H_t>=-79/1000 and kappa_3,t(2*log(U))<0 for 0.9264<=u_t<=5",
            readiness="interval_validated",
            proof_boundary="Large compact all-real-parameter theorem; not the u>=5 ray.",
            diagnostics=compact,
        ),
        CertificateRow(
            id="fsprc_07_far_ray_target",
            role="open_requirement",
            claim="Prove the same paired remainder floor on the asymptotic mode ray u>=5.",
            formula="H_t>=-79/1000 for u_t>=5",
            readiness="open_target",
            proof_boundary="Open analytic ray; finite compact coverage does not prove it.",
            gap="Needs a small-normalized-derivative paired perturbation theorem on q>=pi*exp(20).",
        ),
        CertificateRow(
            id="fsprc_08_wall_handoff",
            role="conditional_handoff",
            claim="The compact theorem plus the open ray would close the first-summand and full-kernel tails.",
            formula="compact+ray => cumulant floor => L_k^(1)>=1/(4k^2) => L_k>0",
            readiness="blocked_by_open_requirement",
            proof_boundary="Conditional handoff only; not cone entry, RH, or Lambda <= 0.",
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "positive_curvature_intervals": geometry["positive_curvature_intervals"],
        "positive_v9_intervals": geometry["positive_v9_intervals"],
        "positive_eighth_envelope_intervals": eighth["positive_interval_count"],
        "compact_base_blocks": compact["base_block_count"],
        "compact_accepted_blocks": compact["accepted_block_count"],
        "compact_negative_cumulant_blocks": compact["negative_cumulant_blocks"],
        "compact_maximum_cumulant3_upper": compact["maximum_cumulant3_upper"],
        "compact_minimum_margin_lower": compact["minimum_margin_lower"],
        "open_ray_rows": 1,
        "ready_to_apply_rows": 0,
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate",
        "date": "2026-07-10",
        "status": "interval compact paired-remainder and negative-skewness theorem with open asymptotic ray",
        "proof_boundary": (
            "This artifact proves the seventh-order paired remainder floor and negative third-cumulant sign on a large "
            "compact all-real mode range. It does not prove the u>=5 ray, the global "
            "cumulant bound, first-summand wall, cone entry, RH, or Lambda <= 0."
        ),
        "source_saddle_expansion": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md"
        ),
        "source_cumulant_bridge": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md"
        ),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py"
        ),
        "diagnostics": {
            "parameters": {
                "precision_bits": DEFAULT_PRECISION_BITS,
                "mode_start": str(DEFAULT_MODE_START),
                "mode_end": str(DEFAULT_MODE_END),
                "base_width": str(DEFAULT_BASE_WIDTH),
                "maximum_refinement_depth": DEFAULT_MAX_REFINEMENT_DEPTH,
                "primary_panels": DEFAULT_PRIMARY_PANELS,
                "mid_panels": DEFAULT_MID_PANELS,
                "retry_panels": DEFAULT_RETRY_PANELS,
                "high_panel_start": str(DEFAULT_HIGH_PANEL_START),
                "window_y": DEFAULT_WINDOW_Y,
                "eighth_envelope": str(DEFAULT_EIGHTH_ENVELOPE),
                "remainder_floor": str(DEFAULT_REMAINDER_FLOOR),
            },
            "geometry": geometry,
            "eighth_derivative_envelope": eighth,
            "compact_remainder": compact,
        },
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    compact = diagnostics["compact_remainder"]
    lines = [
        "# Jensen-Window PF Negative-Lambda First-Summand Paired-Remainder Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: interval compact paired-remainder and negative-skewness theorem with open asymptotic ray.",
        "This is not a proof of the global cumulant wall, cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.py",
        "```",
        "",
        "## Exact Paired Form",
        "",
        "At the log-variable mode put `y=sqrt(a_t)*(x-x_t)` and pair `y` with",
        "`-y`. The four raw moments are integrated together, then normalized and",
        "centralized before evaluating",
        "",
        "```text",
        "H_t=(t^2/a_t^(3/2))*(kappa_3(Y_t)+alpha_t+C3_t+C5_t).",
        "```",
        "",
        "The target is `H_t>=-79/1000`.",
        "",
        "## Taylor And Tails",
        "",
        f"Arb certifies `{summary['positive_v9_intervals']}` positive `V^(9)` intervals and",
        f"`{summary['positive_eighth_envelope_intervals']}` normalized eighth-derivative",
        "envelope intervals, proving",
        "",
        "```text",
        "sup_|y|<=6 V^(8)(x_t+y/sqrt(a_t))/a_t^4<=1/50000.",
        "```",
        "",
        "On `|y|<=6`, a paired midpoint rule uses interval second derivatives.",
        "Beyond each endpoint, the certified outward slope supplies an exact",
        "exponential moment tail. The left-tail use is legitimate because `V'` is",
        "increasing for `u>=1/100` and negative below that threshold.",
        "",
        "## Compact Theorem",
        "",
        "The adaptive cover proves",
        "",
        "```text",
        "H_t>=-79/1000 for every real mode parameter 0.9264<=u_t<=5.",
        "kappa_3,t(2*log(U))<0 for every real mode parameter 0.9264<=u_t<=5.",
        "```",
        "",
        f"- Base blocks: `{compact['base_block_count']}`.",
        f"- Accepted adaptive blocks: `{compact['accepted_block_count']}`.",
        f"- Maximum refinement depth: `{compact['maximum_refinement_depth']}`.",
        f"- Minimum outward-rounded margin: `{compact['minimum_margin_lower']}`.",
        f"- Negative third-cumulant blocks: `{compact['negative_cumulant_blocks']}`.",
        f"- Maximum outward-rounded third-cumulant upper bound: `{compact['maximum_cumulant3_upper']}`.",
        f"- Upper endpoint `t=V'(x(5))`: `{compact['end_t_ball']}`.",
        "",
        "## Remaining Ray",
        "",
        "The only unproved row in this artifact is",
        "",
        "```text",
        "H_t>=-79/1000 for u_t>=5.",
        "```",
        "",
        "Here `q>=pi*exp(20)>10^9` and every normalized derivative is very small.",
        "The next step is an analytic paired perturbation bound on this ray; finite",
        "compact quadrature is not promoted across the unbounded parameter range.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        (
            "Paired interval quadrature, a global eighth-derivative Taylor envelope, "
            "and explicit two-sided tail budgets prove the seventh-order normalized "
            "remainder floor -79/1000 and strict negative third-cumulant sign for every real mode 0.9264<=u<=5. The "
            "asymptotic mode ray u>=5 remains open."
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "validated Jensen-window PF negative-lambda first-summand paired-remainder certificate: "
        f"8 rows, 0 issues, {summary['positive_eighth_envelope_intervals']} eighth-envelope "
        f"intervals, {summary['compact_accepted_blocks']} compact remainder blocks, "
        "1 open asymptotic ray, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
