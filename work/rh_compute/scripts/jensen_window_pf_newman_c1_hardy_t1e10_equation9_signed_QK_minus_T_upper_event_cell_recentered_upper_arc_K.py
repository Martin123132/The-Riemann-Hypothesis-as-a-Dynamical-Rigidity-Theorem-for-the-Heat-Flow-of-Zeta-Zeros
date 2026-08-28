#!/usr/bin/env python3
"""Recenter the upper-arc Hardy operator with a stable action envelope.

The first-subcell donor is left untouched.  This companion evaluates the same
integrals and guards, but uses the envelope-theorem action enclosure already
computed at the local midpoint in the two exponential remainder estimates.
"""

from __future__ import annotations

import math
import os
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_quarter_disk_vertical_arc_tail_arb as donor
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined


EVENT_ATLAS = ROOT / "work/rh_compute/results" / (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_"
    "local_height_event_atlas_gate.json"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(65, more=True),
        "imag_ball": value.imag.str(65, more=True),
        "absolute_ball": abs(value).str(55, more=True),
    }


def real_record(value: arb) -> dict[str, str]:
    return {
        "ball": value.str(65, more=True),
        "lower": value.lower().str(50, more=True),
        "upper": value.upper().str(50, more=True),
    }


def add_complex_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def real_interval(left: arb, right: arb) -> arb:
    require(left.upper() <= right.lower(), "invalid real interval")
    midpoint = (left + right) / 2
    radius = (right - left) / 2
    return arb(midpoint, radius)


def build(
    center_text: str,
    radius_text: str,
    precision_bits: int,
    variant: str = "production",
) -> dict[str, Any]:
    started = time.perf_counter()
    resource_mode = joined.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    require(variant in ("production", "independent"), "unknown arc variant")
    ctx.prec = precision_bits
    ctx.threads = 1
    cfg = donor.configuration(variant)
    radius_box = arb(radius_text)
    require(radius_box > 0, "radius must be positive")
    t_center = arb(center_text)
    t = arb(t_center, radius_box)

    import json

    event = json.loads(EVENT_ATLAS.read_text(encoding="utf-8"))["certificate"]["primary_open_cell"]
    event_lower = arb(event["lower_ball"]["ball"])
    event_upper = arb(event["upper_ball"]["ball"])
    require(event_lower.upper() < t.lower(), "upper-arc box crosses lower event wall")
    require(t.upper() < event_upper.lower(), "upper-arc box crosses upper event wall")

    pi = arb.pi()
    imaginary = acb(0, 1)
    y = arb(donor.Y_TEXT)
    q_minus = arb(donor.Q_MINUS_TEXT)
    radial_scale = 2 * pi * y
    count = donor.COUNT
    sqrt_y = y.sqrt()
    degree = int(cfg["degree"])
    delta_min = arb(str(cfg["delta_min"]))
    delta_cut = arb(str(cfg["delta_cut"]))
    tolerance = arb(str(cfg["tolerance"]))
    initial_panels = int(cfg["initial_panels"])

    theta, theta_prime = joined.theta_and_derivative(t)
    H_log_prime = joined.direct_log_H_derivative(t)
    theta_center, _ = joined.theta_and_derivative(t_center)
    endpoint_phase_center = pi - t_center * y.log() - pi * y * y
    combined_phase_center = theta_center + endpoint_phase_center
    combined_phase_derivative = theta_prime - y.log()
    transport_distance = abs(t - t_center).upper()
    combined_phase = combined_phase_center + arb(
        0,
        transport_distance * abs(combined_phase_derivative).upper(),
    )
    log_H_center = joined.direct_log_H(t_center)
    H = (log_H_center + arb(0, transport_distance * abs(H_log_prime).upper())).exp()
    require(H.lower() > 0, "H lost positivity")

    slope = t + 2 * pi * y * y - q_minus * radial_scale
    require(slope > 0, "upper endpoint action must initially increase")
    exponent_coefficients = [acb(0) for _ in range(degree + 1)]
    w_coefficients = [acb(0) for _ in range(degree + 1)]
    exponent_coefficients[1] = slope + imaginary / 2
    for order in range(2, degree + 1):
        factorial = arb(math.factorial(order))
        exponent_coefficients[order] = (
            -imaginary * pi * y * y * (2 * imaginary) ** order
            + imaginary * q_minus * radial_scale * imaginary**order
        ) / factorial
    for order in range(1, degree + 1):
        w_coefficients[order] = (
            imaginary * radial_scale * imaginary**order / arb(math.factorial(order))
        )

    def components(z: acb) -> tuple[acb, acb]:
        return donor.polynomial(exponent_coefficients, z), donor.polynomial(w_coefficients, z)

    def full_integrand(z: acb, analytic: bool = True) -> acb:
        exponent, w_exponent = components(z)
        w = w_exponent.exp()
        return sqrt_y * exponent.exp() * (1 + (count * w_exponent).exp()) / (1 + w)

    def reduced_integrand(z: acb, analytic: bool = True) -> acb:
        exponent, w_exponent = components(z)
        w = w_exponent.exp()
        return sqrt_y * exponent.exp() / (1 + w)

    hardy_constant = imaginary * (theta_prime - y.log()) - H_log_prime

    def full_K_integrand(z: acb, analytic: bool = True) -> acb:
        return (z + hardy_constant) * full_integrand(z, analytic)

    def reduced_K_integrand(z: acb, analytic: bool = True) -> acb:
        return (z + hardy_constant) * reduced_integrand(z, analytic)

    def full_prime_integrand(z: acb, analytic: bool = True) -> acb:
        return z * full_integrand(z, analytic)

    def reduced_prime_integrand(z: acb, analytic: bool = True) -> acb:
        return z * reduced_integrand(z, analytic)

    initial = acb(0)
    initial_prime = acb(0)
    initial_K = acb(0)
    for index in range(initial_panels):
        left = delta_min * index / initial_panels
        right = delta_min * (index + 1) / initial_panels
        kwargs = {
            "abs_tol": tolerance,
            "rel_tol": tolerance,
            "eval_limit": 150_000,
            "depth_limit": 40,
        }
        initial += acb.integral(full_integrand, left, right, **kwargs)
        initial_prime += acb.integral(full_prime_integrand, left, right, **kwargs)
        initial_K += acb.integral(full_K_integrand, left, right, **kwargs)

    main = acb(0)
    main_prime = acb(0)
    main_K = acb(0)
    points = [arb(text) for text in cfg["main_points"]]
    require(
        (points[0] - delta_min).contains(0) and (points[-1] - delta_cut).contains(0),
        "panel endpoints drifted",
    )
    for left, right in zip(points, points[1:]):
        kwargs = {
            "abs_tol": tolerance,
            "rel_tol": tolerance,
            "eval_limit": 300_000,
            "depth_limit": 45,
        }
        main += acb.integral(reduced_integrand, left, right, **kwargs)
        main_prime += acb.integral(reduced_prime_integrand, left, right, **kwargs)
        main_K += acb.integral(reduced_K_integrand, left, right, **kwargs)

    def stationary_delta(t_value: arb) -> arb:
        constant = t_value - 2 * pi * y * y
        discriminant = (q_minus * radial_scale) ** 2 - 16 * pi * y * y * constant
        root_low = (q_minus * radial_scale - discriminant.sqrt()) / (8 * pi * y * y)
        root_high = (q_minus * radial_scale + discriminant.sqrt()) / (8 * pi * y * y)
        require(root_low < 0, "lower action-derivative root is not negative")
        require(0 < delta_cut.cos() < root_high < 1, "action monotonic tail guard failed")
        return root_high.acos()

    delta_star = stationary_delta(t)

    def action(t_value: arb, delta: arb) -> arb:
        return (
            t_value * delta
            + pi * y * y * (2 * delta).sin()
            - q_minus * radial_scale * delta.sin()
        )

    raw_action_maximum = action(t, delta_star)
    action_at_cut = action(t, delta_cut)
    require(action_at_cut < arb("-80"), "compact arc cutoff is not in the decaying tail")
    phase_deviation = radial_scale * (1 - delta_cut.cos())
    require(phase_deviation < pi / 2, "core geometric denominator left the positive-real sector")

    delta_anchor = stationary_delta(t_center)
    action_anchor = action(t_center, delta_anchor)
    stable_action = action_anchor + arb(0, transport_distance * delta_star.upper())
    require(stable_action.is_finite(), "stable local action is not finite")

    next_order = degree + 1
    factorial = arb(math.factorial(next_order))
    exponent_remainder = (
        pi * y * y * (2 * delta_cut).exp() * (2 * delta_cut) ** next_order
        + q_minus * radial_scale * delta_cut.exp() * delta_cut**next_order
    ) / factorial
    w_remainder = radial_scale * delta_cut.exp() * delta_cut**next_order / factorial
    polynomial_relative_error = (
        2 * exponent_remainder * exponent_remainder.exp()
        + 8
        * (count + 2)
        * w_remainder
        * ((count + 2) * w_remainder).exp()
    )
    polynomial_integral_error = (
        2 * delta_cut * sqrt_y * stable_action.exp() * polynomial_relative_error
    )
    omitted_numerator_error = (
        (delta_cut - delta_min)
        * sqrt_y
        * stable_action.exp()
        * (-count * radial_scale * delta_min.sin()).exp()
    )
    radial_decay = (-radial_scale * delta_cut.sin()).exp()
    quotient_tail = (1 + radial_decay**count) / (1 - radial_decay)
    compact_arc_tail = (
        sqrt_y
        * (pi / 2 - delta_cut)
        * action_at_cut.exp()
        * quotient_tail
    )
    total_added_error = polynomial_integral_error + omitted_numerator_error + compact_arc_tail
    require(total_added_error < arb("1e-20"), "analytic arc error budget is too large")

    core_delta = real_interval(arb(0), delta_cut)
    tail_delta = real_interval(delta_cut, pi / 2)
    core_weight = abs(acb(core_delta) + hardy_constant).upper()
    tail_weight = abs(acb(tail_delta) + hardy_constant).upper()
    core_error = polynomial_integral_error + omitted_numerator_error
    K_error = core_weight * core_error + tail_weight * compact_arc_tail
    prime_error = delta_cut * core_error + (pi / 2) * compact_arc_tail
    require(K_error < arb("1e-19"), "Hardy-operator arc error budget is too large")

    rotated_arc = add_complex_error(initial + main, total_added_error)
    rotated_prime = add_complex_error(initial_prime + main_prime, prime_error)
    rotated_K_direct = add_complex_error(initial_K + main_K, K_error)
    endpoint_phase = pi - t * y.log() - pi * y * y
    endpoint_factor = (imaginary * endpoint_phase).exp()
    actual_arc = endpoint_factor * rotated_arc
    actual_arc_prime = endpoint_factor * (
        rotated_prime - imaginary * y.log() * rotated_arc
    )
    actual_K_assembled = (
        actual_arc_prime + (imaginary * theta_prime - H_log_prime) * actual_arc
    )
    actual_K_direct = endpoint_factor * rotated_K_direct
    require(
        actual_K_direct.overlaps(actual_K_assembled),
        "direct and assembled upper-arc K miss",
    )

    hardy_contribution = (
        2 * ((imaginary * combined_phase).exp() * rotated_K_direct).real / H
    )
    return {
        "kind": "rh_diagnostic_recentered_upper_event_cell_upper_arc_K",
        "date": "2026-08-28",
        "status": "recentered_upper_arc_Hardy_operator_diagnostic",
        "passed": True,
        "resource_mode": resource_mode,
        "workers": 1,
        "precision_bits": precision_bits,
        "variant": variant,
        "center_height": center_text,
        "radius": radius_text,
        "height_ball": t.str(60, more=True),
        "configuration": {
            "degree": degree,
            "delta_min": str(delta_min),
            "delta_cut": str(delta_cut),
            "initial_panels": initial_panels,
            "main_panel_count": len(points) - 1,
            "tolerance": str(tolerance),
        },
        "action_guard": {
            "raw_dependent_action_ball": raw_action_maximum.str(60, more=True),
            "stable_local_action_ball": stable_action.str(60, more=True),
            "envelope_identity": "d A_max(t)/dt=delta_*(t)",
            "stable_action_used_in_both_exponential_remainder_terms": True,
        },
        "error_budgets": {
            "polynomial_integral_error": real_record(polynomial_integral_error),
            "omitted_numerator_error": real_record(omitted_numerator_error),
            "compact_arc_tail": real_record(compact_arc_tail),
            "total_arc_error": real_record(total_added_error),
            "rotated_derivative_error": real_record(prime_error),
            "rotated_K_error": real_record(K_error),
            "core_operator_weight_upper": core_weight.str(40, more=True),
            "tail_operator_weight_upper": tail_weight.str(40, more=True),
        },
        "phase_transport": {
            "combined_phase_center": combined_phase_center.str(60, more=True),
            "combined_phase_derivative": combined_phase_derivative.str(60, more=True),
            "combined_phase_ball": combined_phase.str(60, more=True),
            "maximum_transport_distance": transport_distance.str(40, more=True),
            "H_ball": H.str(60, more=True),
        },
        "rotated_arc": complex_record(rotated_arc),
        "rotated_arc_prime": complex_record(rotated_prime),
        "rotated_K_direct": complex_record(rotated_K_direct),
        "actual_arc": complex_record(actual_arc),
        "actual_K_direct": complex_record(actual_K_direct),
        "actual_K_assembled": complex_record(actual_K_assembled),
        "direct_and_assembled_K_overlap": True,
        "upper_arc_Q_derivative_contribution": real_record(hardy_contribution),
        "elapsed_seconds": time.perf_counter() - started,
        "diagnostic_boundary": (
            "One recentered upper-arc Hardy-operator box only. Other K_T components, connected "
            "height transport, event-wall handoffs, and RH are outside this certificate."
        ),
    }
