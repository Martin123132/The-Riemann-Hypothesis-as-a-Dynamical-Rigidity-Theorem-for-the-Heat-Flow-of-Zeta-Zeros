#!/usr/bin/env python3
"""Certify a generic rational-phase terminal and one adjacent recursive tile."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_rational_phase_third_level_adjacent_tile_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
TWO_SIDED_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_two_sided_wall_period_three_recursive_box_gate.py"
)
TWO_SIDED_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_two_sided_wall_period_three_recursive_box_gate.json"
)

A = 159_577
L = 2_481_422
K = L - 1
FIRST_P = Fraction(31_916)
FIRST_N = 992_568
OPEN_CHILD_N = 496_283
X0 = Fraction(2, 5)
S0 = Fraction(1, 3)
RIGHT_BASE_P = Fraction(95_747, 6)
RIGHT_BASE_A = Fraction(-1, 3)
LEFT_BASE_P = Fraction(47_873, 3)
LEFT_BASE_A = Fraction(1, 3)
ADJACENT_WIDTH = Fraction(1, 10**23)
SCALE_POWERS = (20, 21, 22, 23, 24, 26)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def fractional_part(value: Fraction) -> Fraction:
    return value - floor_fraction(value)


def minimal_phase_period(a: Fraction, tau: Fraction, limit: int = 131_072) -> int:
    """Return the exact period of e(a*k+tau*k^2)."""

    for period in range(1, limit + 1):
        if (2 * tau * period).denominator == 1 and (
            a * period + tau * period * period
        ).denominator == 1:
            return period
    raise RuntimeError("rational quadratic phase period not found")


def add_coefficients(
    left: dict[Fraction, Fraction],
    right: dict[Fraction, Fraction],
) -> dict[Fraction, Fraction]:
    result = dict(left)
    for phase, coefficient in right.items():
        result[phase] = result.get(phase, Fraction(0)) + coefficient
        if result[phase] == 0:
            del result[phase]
    return result


def scale_coefficients(
    values: dict[Fraction, Fraction], scalar: Fraction
) -> dict[Fraction, Fraction]:
    return {phase: scalar * coefficient for phase, coefficient in values.items() if scalar * coefficient}


def rational_phase_power_coefficients(
    n: int,
    degree: int,
    a0: Fraction,
    tau0: Fraction,
    small: Any,
) -> tuple[int, dict[Fraction, Fraction]]:
    """Group an exact power moment by root-of-unity phase."""

    require(n >= 0 and 0 <= degree <= 4, "unsupported rational-phase moment")
    period = minimal_phase_period(a0, tau0)
    coefficients: dict[Fraction, Fraction] = {}
    for residue in range(period):
        coefficient = Fraction(residue_power_sum(n, residue, period, degree, small))
        phase = fractional_part(a0 * residue + tau0 * residue * residue)
        coefficients[phase] = coefficients.get(phase, Fraction(0)) + coefficient
    return period, {phase: value for phase, value in coefficients.items() if value}


def coefficient_ball(values: dict[Fraction, Fraction], prior: Any) -> acb:
    result = acb(0)
    for phase, coefficient in sorted(values.items()):
        result += acb(prior.arb_rational(coefficient)) * prior.acb_cis_pi(2 * phase)
    return result


def coefficient_record(values: dict[Fraction, Fraction]) -> list[dict[str, str]]:
    return [
        {"phase_turns": str(phase), "coefficient": str(coefficient)}
        for phase, coefficient in sorted(values.items())
    ]


def residue_power_sum(
    n: int,
    residue: int,
    modulus: int,
    degree: int,
    small: Any,
) -> int:
    """Return sum k^degree for one residue class in 0..n."""

    require(modulus >= 1 and 0 <= residue < modulus, "invalid residue class")
    require(0 <= degree <= 5, "unsupported residue degree")
    if n < residue:
        return 0
    count = (n - residue) // modulus + 1
    return sum(
        math.comb(degree, j)
        * residue ** (degree - j)
        * modulus**j
        * small.power_sum(count, j)
        for j in range(degree + 1)
    )


def max_distance(value: Any, center: Fraction) -> Fraction:
    return max(abs(value.lo - center), abs(value.hi - center))


def rational_phase_taylor_box(
    p_box: Any,
    a_box: Any,
    tau_box: Any,
    n: int,
    base_p: Fraction,
    base_a: Fraction,
    base_tau: Fraction,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    """Enclose an affine theta current around any exact rational phase."""

    require(p_box.lo > 0, "rational-phase weight crossed zero")
    f_coefficients: list[dict[Fraction, Fraction]] = []
    periods: list[int] = []
    for degree in range(4):
        period, coefficients = rational_phase_power_coefficients(
            n, degree, base_a, base_tau, small
        )
        periods.append(period)
        f_coefficients.append(coefficients)
    require(len(set(periods)) == 1, "phase period drift between moments")
    m_coefficients = [
        add_coefficients(
            scale_coefficients(f_coefficients[degree], base_p),
            f_coefficients[degree + 1],
        )
        for degree in range(3)
    ]
    f_balls = [coefficient_ball(values, prior) for values in f_coefficients]
    m_balls = [coefficient_ball(values, prior) for values in m_coefficients]

    rho = box.arb_box(p_box) - prior.arb_rational(base_p)
    alpha = box.arb_box(a_box) - prior.arb_rational(base_a)
    beta = box.arb_box(tau_box) - prior.arb_rational(base_tau)
    b0 = m_balls[0] + acb(rho) * f_balls[0]
    b1 = m_balls[1] + acb(rho) * f_balls[1]
    b2 = m_balls[2] + acb(rho) * f_balls[2]
    first_order = b0 + 2 * acb.pi() * acb(0, 1) * (
        acb(alpha) * b1 + acb(beta) * b2
    )

    delta_p = max_distance(p_box, base_p)
    delta_a = max_distance(a_box, base_a)
    delta_tau = max_distance(tau_box, base_tau)
    masses = {degree: small.weighted_mass(p_box.hi, n, degree) for degree in range(2, 5)}
    remainder_mass = (
        delta_a * delta_a * masses[2]
        + 2 * delta_a * delta_tau * masses[3]
        + delta_tau * delta_tau * masses[4]
    )
    remainder_radius = 2 * arb.pi() ** 2 * prior.arb_rational(remainder_mass)
    enclosure = first_order + box.symmetric_complex_error(remainder_radius)

    dp = prior.arb_rational(delta_p)
    da = prior.arb_rational(delta_a)
    dt = prior.arb_rational(delta_tau)
    linear_radius = dp * f_balls[0].abs_upper() + 2 * arb.pi() * (
        da * (m_balls[1].abs_upper() + dp * f_balls[1].abs_upper())
        + dt * (m_balls[2].abs_upper() + dp * f_balls[2].abs_upper())
    )
    base_current, direct_period = affine.weighted_period_ball(
        base_p, base_a, base_tau, n, prior
    )
    require(direct_period == periods[0], "generic and direct phase periods disagree")
    require(m_balls[0].overlaps(base_current), "generic residue moment misses direct base current")
    anchor_in_box = (
        p_box.lo <= base_p <= p_box.hi
        and a_box.lo <= base_a <= a_box.hi
        and tau_box.lo <= base_tau <= tau_box.hi
    )
    if anchor_in_box:
        require(enclosure.overlaps(base_current), "Taylor box misses contained rational anchor")
    return enclosure, {
        "n": n,
        "base_p": str(base_p),
        "base_a": str(base_a),
        "base_tau": str(base_tau),
        "base_period": periods[0],
        "p_box": box.box_record(p_box),
        "a_box": box.box_record(a_box),
        "tau_box": box.box_record(tau_box),
        "anchor_is_inside_parameter_box": anchor_in_box,
        "delta_p_max": str(delta_p),
        "delta_a_max": str(delta_a),
        "delta_tau_max": str(delta_tau),
        "F_phase_coefficients": {
            f"F{degree}": coefficient_record(values)
            for degree, values in enumerate(f_coefficients)
        },
        "M_phase_coefficients": {
            f"M{degree}": coefficient_record(values)
            for degree, values in enumerate(m_coefficients)
        },
        "F_balls": {f"F{degree}": box.acb_record(value) for degree, value in enumerate(f_balls)},
        "M_balls": {f"M{degree}": box.acb_record(value) for degree, value in enumerate(m_balls)},
        "S2": str(masses[2]),
        "S3": str(masses[3]),
        "S4": str(masses[4]),
        "first_order_variation_radius_upper": box.arb_upper_text(linear_radius),
        "second_order_remainder_radius_upper": box.arb_upper_text(remainder_radius),
        "total_deviation_from_base_upper": box.arb_upper_text(linear_radius + remainder_radius),
        "first_order_enclosure": box.acb_record(first_order),
        "rational_phase_enclosure": box.acb_record(enclosure),
        "exact_base_current": box.acb_record(base_current),
        "residue_base_overlaps_direct_base": True,
        "passed": True,
    }


def exact_right_parent_boxes(x_box: Any, s_box: Any, box: Any) -> tuple[Any, Any]:
    a_box = box.RationalBox(
        FIRST_P / x_box.hi - Fraction(A, 2) - s_box.hi - 1,
        FIRST_P / x_box.lo - Fraction(A, 2) - s_box.lo - 1,
    )
    tau_box = box.RationalBox(
        Fraction(1, 2) / x_box.hi - 1,
        Fraction(1, 2) / x_box.lo - 1,
    )
    require(Fraction(-1, 2) < a_box.lo <= a_box.hi < Fraction(1, 2), "right parent crossed an a wall")
    require(Fraction(0) < tau_box.lo <= tau_box.hi < Fraction(1, 4), "adjacent tile left right branch")
    return a_box, tau_box


def exact_right_normalized_child_boxes(
    x_box: Any, s_box: Any, box: Any
) -> tuple[Any, Any, Any, dict[str, Any]]:
    """Use symbolic cancellation before interval evaluation."""

    c_lo = Fraction(A, 2) + s_box.lo + 1
    c_hi = Fraction(A, 2) + s_box.hi + 1
    p_box = box.RationalBox(c_lo - 2 * FIRST_P, c_hi - 2 * FIRST_P)
    raw_a = box.RationalBox(
        (FIRST_P - c_hi * x_box.hi) / (1 - 2 * x_box.hi),
        (FIRST_P - c_lo * x_box.lo) / (1 - 2 * x_box.lo),
    )
    raw_tau = box.RationalBox(
        -x_box.hi / (2 * (1 - 2 * x_box.hi)),
        -x_box.lo / (2 * (1 - 2 * x_box.lo)),
    )
    a_box = box.RationalBox(-raw_a.hi, -raw_a.lo)
    tau_box = box.RationalBox(-raw_tau.hi - 1, -raw_tau.lo - 1)
    require(a_box.lo > Fraction(-1, 2) and a_box.hi < Fraction(1, 2), "child nearest-integer wall crossed")
    require(tau_box.lo > 0, "adjacent child did not enter the conjugated small-tau branch")
    old_a, old_tau = exact_right_parent_boxes(x_box, s_box, box)
    old_p_box, _, _ = box_dependency_child_boxes(old_a, old_tau, box)
    require(old_p_box.lo <= p_box.lo <= p_box.hi <= old_p_box.hi, "exact p cancellation escaped legacy hull")
    return p_box, a_box, tau_box, {
        "raw_child_a_box": box.box_record(raw_a),
        "raw_child_tau_box": box.box_record(raw_tau),
        "normalization": [
            "tau -> tau+1",
            "conjugate (a,tau) -> (-a,-tau)",
            "nearest integer a shift remains zero",
        ],
        "exact_weight_identity": "p_2=A/2+s+1-2*31916",
        "legacy_dependency_p_box": box.box_record(old_p_box),
        "exact_cancellation_p_box": box.box_record(p_box),
        "p_radius_reduction_factor": str(old_p_box.radius / p_box.radius),
    }


def box_dependency_child_boxes(a_box: Any, tau_box: Any, box: Any) -> tuple[Any, Any, Any]:
    p_box = box.RationalBox(
        2 * tau_box.lo * FIRST_P - a_box.hi,
        2 * tau_box.hi * FIRST_P - a_box.lo,
    )
    a_child = box.RationalBox(
        a_box.lo / (2 * tau_box.hi),
        a_box.hi / (2 * tau_box.lo),
    )
    tau_child = box.RationalBox(
        -Fraction(1, 4) / tau_box.lo,
        -Fraction(1, 4) / tau_box.hi,
    )
    return p_box, a_child, tau_child


def right_parent_complete_box(
    x_box: Any,
    s_box: Any,
    settings: Any,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    a_parent, tau_parent = exact_right_parent_boxes(x_box, s_box, box)
    floor_box = box.RationalBox(2 * FIRST_N * tau_parent.lo, 2 * FIRST_N * tau_parent.hi)
    require(OPEN_CHILD_N < floor_box.lo <= floor_box.hi < OPEN_CHILD_N + 1, "second child floor wall crossed")
    p_child, a_child, tau_child, cancellation = exact_right_normalized_child_boxes(x_box, s_box, box)
    child_normalized, child_record = rational_phase_taylor_box(
        p_child,
        a_child,
        tau_child,
        OPEN_CHILD_N,
        RIGHT_BASE_P,
        RIGHT_BASE_A,
        Fraction(0),
        box,
        small,
        affine,
        prior,
    )
    child_original = child_normalized.conjugate()

    phase_box = box.RationalBox(
        Fraction(1, 4) - a_parent.hi * a_parent.hi / (2 * tau_parent.lo),
        Fraction(1, 4) - a_parent.lo * a_parent.lo / (2 * tau_parent.hi),
    )
    tt = box.arb_box(tau_parent)
    coefficient = box.cis_pi_interval(box.arb_box(phase_box)) / (2 * acb(tt) * acb(2 * tt).sqrt())
    main = coefficient * child_original

    c_lo = Fraction(A, 2) + s_box.lo + 1
    c_hi = Fraction(A, 2) + s_box.hi + 1
    lower_argument = box.RationalBox(
        (FIRST_P - Fraction(1, 2)) / x_box.hi - Fraction(A, 2) - s_box.hi + Fraction(1, 2),
        (FIRST_P - Fraction(1, 2)) / x_box.lo - Fraction(A, 2) - s_box.lo + Fraction(1, 2),
    )
    d = FIRST_P + FIRST_N + Fraction(1, 2)
    upper_argument = box.RationalBox(
        d / x_box.hi - c_hi - 2 * FIRST_N - OPEN_CHILD_N - Fraction(3, 2),
        d / x_box.lo - c_lo - 2 * FIRST_N - OPEN_CHILD_N - Fraction(3, 2),
    )
    mordell_tau = box.RationalBox(2 * tau_parent.lo, 2 * tau_parent.hi)
    lower_h, lower_hz, lower_record = box.mordell_target_parameter_box(
        lower_argument, mordell_tau, -1, settings, prior
    )
    upper_h, upper_hz, upper_record = box.mordell_target_parameter_box(
        upper_argument, mordell_tau, -1, settings, prior
    )
    lower_phase = box.RationalBox(
        c_lo - Fraction(1, 2) - (FIRST_P - Fraction(1, 4)) / x_box.lo,
        c_hi - Fraction(1, 2) - (FIRST_P - Fraction(1, 4)) / x_box.hi,
    )
    half_n = Fraction(FIRST_N) + Fraction(1, 2)
    upper_inner = box.RationalBox(
        (FIRST_P + half_n / 2) / x_box.hi - c_hi - half_n,
        (FIRST_P + half_n / 2) / x_box.lo - c_lo - half_n,
    )
    upper_phase = box.RationalBox(2 * half_n * upper_inner.lo, 2 * half_n * upper_inner.hi)
    pi_i = acb.pi() * acb(0, 1)
    endpoint = -acb(0, 1) * (
        box.cis_pi_interval(box.arb_box(lower_phase))
        * (acb(prior.arb_rational(2 * FIRST_P - 1)) * lower_h + lower_hz / pi_i)
        + ((-1) ** OPEN_CHILD_N)
        * box.cis_pi_interval(box.arb_box(upper_phase))
        * (
            acb(prior.arb_rational(2 * FIRST_P + 2 * FIRST_N + 1)) * upper_h
            + upper_hz / pi_i
        )
    ) / 4
    complete = main + endpoint

    third_floor = box.RationalBox(
        2 * OPEN_CHILD_N * tau_child.lo,
        2 * OPEN_CHILD_N * tau_child.hi,
    )
    require(0 < third_floor.lo <= third_floor.hi < 1, "third child floor wall crossed")
    third_p = box.RationalBox(
        2 * tau_child.lo * p_child.lo - a_child.hi,
        2 * tau_child.hi * p_child.hi - a_child.lo,
    )
    small_tt = box.arb_box(tau_child)
    third_coefficient_upper = (
        1 / (2 * small_tt * (2 * small_tt).sqrt())
    ).abs_upper()
    return complete, {
        "parent_p": str(FIRST_P),
        "parent_n": FIRST_N,
        "parent_a_box": box.box_record(a_parent),
        "parent_tau_box": box.box_record(tau_parent),
        "second_child_floor_argument_box": box.box_record(floor_box),
        "second_child_n": OPEN_CHILD_N,
        "normalized_child_p_box": box.box_record(p_child),
        "normalized_child_a_box": box.box_record(a_child),
        "normalized_child_tau_box": box.box_record(tau_child),
        "cancellation_audit": cancellation,
        "rational_phase_terminal": child_record,
        "original_child_current": box.acb_record(child_original),
        "second_main_phase_box": box.box_record(phase_box),
        "second_main_coefficient": box.acb_record(coefficient),
        "second_recursive_main": box.acb_record(main),
        "second_lower_argument_box": box.box_record(lower_argument),
        "second_upper_argument_box": box.box_record(upper_argument),
        "second_lower_Mordell": lower_record,
        "second_upper_Mordell": upper_record,
        "second_lower_phase_box": box.box_record(lower_phase),
        "second_upper_phase_box": box.box_record(upper_phase),
        "second_joined_endpoint": box.acb_record(endpoint),
        "complete_normalized_parent": box.acb_record(complete),
        "third_level_branch": {
            "normalized_nearest_integer_shift": 0,
            "third_floor_argument_box": box.box_record(third_floor),
            "third_child_n": 0,
            "third_child_p_box": box.box_record(third_p),
            "formal_main_coefficient_absolute_upper": box.arb_upper_text(third_coefficient_upper),
            "decision": "terminate_by_rational_phase_Taylor_before_small_tau_Mordell_cancellation",
        },
        "passed": True,
    }


def adjacent_source_tile(
    width: Fraction,
    settings: Any,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
    previous: dict[str, Any],
) -> dict[str, Any]:
    flint.ctx.prec = settings.precision_bits
    x_box = box.RationalBox(X0 + width, X0 + 2 * width)
    s_box = box.RationalBox(S0 - width, S0 + width)
    normalized_parent, parent_record = right_parent_complete_box(
        x_box, s_box, settings, box, small, affine, prior
    )
    original_transformed = normalized_parent.conjugate()

    first_phase = box.RationalBox(
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.lo) - FIRST_P * FIRST_P / x_box.lo,
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.hi) - FIRST_P * FIRST_P / x_box.hi,
    )
    xx = box.arb_box(x_box)
    first_multiplier = acb(2 / (xx * xx.sqrt())) * box.cis_pi_interval(box.arb_box(first_phase))
    first_main = first_multiplier * original_transformed
    first_endpoint, endpoint_record = box.physical_endpoint_parameter_box(
        "adjacent_right_third_level_first_joined_endpoint", x_box, s_box, settings, prior
    )
    complete_source = first_main + first_endpoint
    previous_source = parse_acb(previous["certified_two_sided_wall_box"]["complete_source_box"])
    require(complete_source.overlaps(previous_source), "adjacent source tile does not join prior wall box")
    return {
        "half_width_unit": str(width),
        "x_box": box.box_record(x_box),
        "s_box": box.box_record(s_box),
        "shares_boundary_with_previous_tile_at_x": str(X0 + width),
        "right_normalization": "conjugate after the first physical transform",
        "complete_normalized_parent": parent_record,
        "original_transformed_parent": box.acb_record(original_transformed),
        "first_phase_box": box.box_record(first_phase),
        "first_multiplier": box.acb_record(first_multiplier),
        "first_recursive_main": box.acb_record(first_main),
        "first_joined_endpoint": endpoint_record,
        "complete_source_box": box.acb_record(complete_source),
        "complete_radius_absolute_upper": box.arb_upper_text(complete_source.rad()),
        "overlaps_previous_two_sided_wall_source_box": True,
        "passed": True,
    }


def branch_wall_inventory(width: Fraction, box: Any) -> dict[str, Any]:
    n = OPEN_CHILD_N
    right_first = Fraction(2 * n + 1, 5 * n + 2)
    left_first = Fraction(2 * n + 1, 5 * n + 3)
    tile_right = box.RationalBox(X0 + width, X0 + 2 * width)
    right_t = box.RationalBox(
        (5 * tile_right.lo - 2) / (2 * (1 - 2 * tile_right.lo)),
        (5 * tile_right.hi - 2) / (2 * (1 - 2 * tile_right.hi)),
    )
    require(tile_right.hi < right_first, "adjacent tile reached the first third-floor wall")
    return {
        "open_second_child_n": n,
        "normalization_wall": str(X0),
        "right_small_tau_map": "t_R=(5*x-2)/(2*(1-2*x))",
        "left_small_tau_map": "t_L=(2-5*x)/(2*(3*x-1))",
        "right_qth_floor_wall": "x_R(q)=(2*n+q)/(5*n+2*q)",
        "left_qth_floor_wall": "x_L(q)=(2*n+q)/(5*n+3*q)",
        "right_first_positive_floor_wall": str(right_first),
        "left_first_positive_floor_wall": str(left_first),
        "right_first_wall_distance_from_2_over_5": str(right_first - X0),
        "left_first_wall_distance_from_2_over_5": str(X0 - left_first),
        "adjacent_width_units_to_first_right_floor_wall": str((right_first - X0) / width),
        "adjacent_right_tile": box.box_record(tile_right),
        "adjacent_right_small_tau_box": box.box_record(right_t),
        "adjacent_right_third_floor_argument_box": box.box_record(
            box.RationalBox(2 * n * right_t.lo, 2 * n * right_t.hi)
        ),
        "adjacent_tile_is_strictly_inside_m3_zero_branch": True,
    }


def generic_terminal_audit(
    box: Any, small: Any, affine: Any, prior: Any
) -> list[dict[str, Any]]:
    cases = (
        ("alternating_period_two", Fraction(39_894), Fraction(-1, 2), Fraction(0), 1_240_710, 2),
        ("right_period_three", Fraction(95_747, 6), Fraction(1, 3), Fraction(-1), 496_283, 3),
        ("left_period_three", Fraction(47_873, 3), Fraction(2, 3), Fraction(-1), 496_283, 3),
    )
    rows = []
    for name, p0, a0, tau0, n, expected_period in cases:
        point_p = box.RationalBox(p0, p0)
        point_a = box.RationalBox(a0, a0)
        point_tau = box.RationalBox(tau0, tau0)
        enclosure, record = rational_phase_taylor_box(
            point_p, point_a, point_tau, n, p0, a0, tau0, box, small, affine, prior
        )
        require(record["base_period"] == expected_period, f"{name} period drift")
        require(arb(record["second_order_remainder_radius_upper"]) == 0, f"{name} point remainder")
        rows.append(
            {
                "name": name,
                "expected_period": expected_period,
                "derived_period": record["base_period"],
                "base_current": record["exact_base_current"],
                "generic_enclosure": box.acb_record(enclosure),
                "generic_contains_anchor": True,
            }
        )
    return rows


def scale_audit(
    box: Any, small: Any, affine: Any, prior: Any
) -> list[dict[str, Any]]:
    rows = []
    for power in SCALE_POWERS:
        width = Fraction(1, 10**power)
        x_box = box.RationalBox(X0 + width, X0 + 2 * width)
        s_box = box.RationalBox(S0 - width, S0 + width)
        p_child, a_child, tau_child, cancellation = exact_right_normalized_child_boxes(
            x_box, s_box, box
        )
        _, record = rational_phase_taylor_box(
            p_child,
            a_child,
            tau_child,
            OPEN_CHILD_N,
            RIGHT_BASE_P,
            RIGHT_BASE_A,
            Fraction(0),
            box,
            small,
            affine,
            prior,
        )
        rows.append(
            {
                "width_decimal": f"1e-{power}",
                "width": str(width),
                "delta_p_max": record["delta_p_max"],
                "delta_a_max": record["delta_a_max"],
                "delta_tau_max": record["delta_tau_max"],
                "child_total_deviation_upper": record["total_deviation_from_base_upper"],
                "child_second_order_remainder_upper": record["second_order_remainder_radius_upper"],
                "legacy_to_exact_p_radius_factor": cancellation["p_radius_reduction_factor"],
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    tile = artifact["certified_adjacent_third_level_tile"]
    parent = tile["complete_normalized_parent"]
    child = parent["rational_phase_terminal"]
    walls = artifact["recursive_branch_wall_inventory"]
    rows = [
        "| width unit | delta p | delta a | delta tau | child deviation | quadratic remainder | p-radius reduction |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for item in artifact["scale_audit"]:
        rows.append(
            f"| `{item['width_decimal']}` | `{item['delta_p_max']}` | `{item['delta_a_max']}` | "
            f"`{item['delta_tau_max']}` | `{item['child_total_deviation_upper']}` | "
            f"`{item['child_second_order_remainder_upper']}` | `{item['legacy_to_exact_p_radius_factor']}` |"
        )
    return f"""# Rational-phase terminal and adjacent third-level tile

Date: 2026-08-25

Status: certified local recursive-tile lemma; finite-domain tiling remains open

## Generic rational-phase terminal

For rational `(a_0,tau_0)`, the phase `e(a_0 k+tau_0 k^2)` has exact period
`P` when

```text
(RP1) 2*tau_0*P is an integer,
      a_0*P+tau_0*P^2 is an integer.
```

Each power moment is grouped exactly by residues modulo `P`.  If
`rho=p-p_0`, `alpha=a-a_0`, and `beta=tau-tau_0`, the common evaluator is

```text
(RP2) B_j=M_j+rho F_j,

(RP3) W=B_0+2*pi*i(alpha B_1+beta B_2)+R_2,

(RP4) |R_2| <= 2*pi^2(delta_a^2 S_2
                       +2 delta_a delta_tau S_3
                       +delta_tau^2 S_4).
```

This reproduces the old alternating period-two terminal and both conjugate
period-three terminals.  The root-of-unity coefficients are stored exactly;
`pi` enters only through `e(u)=exp(2*pi*i*u)` and Taylor's exponential
remainder.

## Dependency cancellation

On the right branch the second affine weight simplifies before interval
evaluation:

```text
(RP5) p_2=2*tau_1*31916-a_1=A/2+s+1-2*31916.
```

Thus `p_2` is independent of `x`.  On the certified adjacent tile this reduces
the inherited interval radius by a factor
`{parent['cancellation_audit']['p_radius_reduction_factor']}`.

## Next branch walls

After integer-curvature normalization and conjugation, the two open
period-three children have

```text
(RP6) t_R=(5x-2)/(2(1-2x)),
      t_L=(2-5x)/(2(3x-1)).
```

For child length `n={walls['open_second_child_n']}`, their floor walls are

```text
(RP7) x_R(q)=(2n+q)/(5n+2q),
      x_L(q)=(2n+q)/(5n+3q).
```

The first positive right wall is `{walls['right_first_positive_floor_wall']}`,
at distance `{walls['right_first_wall_distance_from_2_over_5']}` from `2/5`.
The adjacent `10^-23` tile is therefore strictly inside the `m_3=0` branch.
A uniform march at this width would require more than
`{walls['adjacent_width_units_to_first_right_floor_wall']}` tile units merely
to reach that first wall, so fixed-width enumeration is not a viable cover
strategy.

## Adjacent source tile

The certified box is

```text
x in [2/5+10^-23,2/5+2*10^-23],
s in [1/3-10^-23,1/3+10^-23].
```

Its rational-phase terminal has period `{child['base_period']}` and deviation
below `{child['total_deviation_from_base_upper']}`.  The second recursive main,
its two joined Mordell endpoint currents, the first physical main, and both
physical endpoint currents are all retained.  The complete source radius is
below `{tile['complete_radius_absolute_upper']}` and its source box overlaps
the previously certified wall tile at their shared boundary.

Formally taking one more Mordell step would have child index zero, but its main
coefficient is already bounded only by
`{parent['third_level_branch']['formal_main_coefficient_absolute_upper']}`.
That `t^(-3/2)` inflation would require cancellation against the small-`t`
endpoint terms.  The rational-phase terminal avoids discarding that
cancellation; it does not claim a useful third Mordell decomposition.

## Scale audit

{chr(10).join(rows)}

## Proof boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(TWO_SIDED_BUILDER.is_file() and TWO_SIDED_RESULT.is_file() and CHECKER.is_file(), "missing dependency or checker")
    previous = json.loads(TWO_SIDED_RESULT.read_text(encoding="utf-8"))
    require(previous.get("passed") is True, "two-sided dependency not passed")
    two_sided = load_module("mordell_rational_phase_two_sided_dependency", TWO_SIDED_BUILDER)
    right = load_module("mordell_rational_phase_right_dependency", two_sided.RIGHT_BUILDER)
    small = load_module("mordell_rational_phase_small_dependency", right.SMALL_BUILDER)
    affine = load_module("mordell_rational_phase_affine_dependency", right.AFFINE_BUILDER)
    box = load_module("mordell_rational_phase_box_dependency", right.BOX_BUILDER)
    prior = load_module("mordell_rational_phase_point_dependency", right.PRIOR_BUILDER)
    settings = box.PRODUCTION

    generic_audit = generic_terminal_audit(box, small, affine, prior)
    walls = branch_wall_inventory(ADJACENT_WIDTH, box)
    tile = adjacent_source_tile(
        ADJACENT_WIDTH, settings, box, small, affine, prior, previous
    )
    audits = scale_audit(box, small, affine, prior)
    selected = next(row for row in audits if row["width_decimal"] == "1e-23")
    require(arb(selected["child_total_deviation_upper"]) < arb("0.001"), "adjacent terminal target failed")
    require(arb(tile["complete_radius_absolute_upper"]) < arb("0.02"), "adjacent source tile target failed")

    proof_boundary = (
        "One exact rational-phase Taylor terminal theorem, exact next-wall formulas for the two local "
        "period-three children, and one endpoint-complete adjacent right tile of width unit 1e-23 only. "
        "No finite recursive cover, useful third Mordell decomposition at small tau, physical quadrature, "
        "non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or "
        "prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "generic_rational_phase_terminal_and_adjacent_third_level_tile_certified_finite_recursive_cover_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "L": L,
            "K": K,
            "workers": 1,
            "adjacent_width_unit": str(ADJACENT_WIDTH),
        },
        "production_settings": settings.__dict__,
        "exact_theorem": {
            "phase_period_conditions": ["2*tau*P integer", "a*P+tau*P^2 integer"],
            "terminal": "B0+2*pi*i*(alpha*B1+beta*B2)+R2",
            "remainder": "2*pi^2*(delta_a^2*S2+2*delta_a*delta_tau*S3+delta_tau^2*S4)",
            "right_exact_child_weight": "A/2+s+1-2*31916",
            "left_exact_child_weight": "3*31916-3/2-A/2-s",
        },
        "generic_terminal_reproduction_audit": generic_audit,
        "recursive_branch_wall_inventory": walls,
        "certified_adjacent_third_level_tile": tile,
        "scale_audit": audits,
        "summary": {
            "generic_periods_reproduced": [2, 3],
            "adjacent_tile_third_child_index": 0,
            "selected_child_deviation_upper": selected["child_total_deviation_upper"],
            "selected_exact_p_radius_reduction_factor": selected["legacy_to_exact_p_radius_factor"],
            "complete_source_radius_upper": tile["complete_radius_absolute_upper"],
            "joins_previous_wall_tile": True,
        },
        "decision": {
            "generic_rational_phase_terminal_built": True,
            "period_two_and_three_special_cases_reproduced": True,
            "next_left_and_right_floor_walls_derived": True,
            "symbolic_child_weight_cancellation_used": True,
            "adjacent_endpoint_complete_source_tile_certified": True,
            "small_tau_third_Mordell_step_useful": False,
            "finite_recursive_cover_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Exploit symbolic cancellation in every child coordinate and endpoint argument, then measure "
            "the complete-source Lipschitz radius on geometrically widened tiles up to the first m3 wall. "
            "Develop a cancellation-preserving small-tau endpoint expansion only if direct rational-phase "
            "termination ceases to scale."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": previous["primary_source"],
        "dependencies": {
            "two_sided_result": {"path": relative(TWO_SIDED_RESULT), "sha256": file_hash(TWO_SIDED_RESULT)},
            "two_sided_builder": {"path": relative(TWO_SIDED_BUILDER), "sha256": file_hash(TWO_SIDED_BUILDER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified the rational-phase terminal and adjacent third-level source tile", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
