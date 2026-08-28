#!/usr/bin/env python3
"""Certify Mordell parameter boxes and audit direct current-box scaling."""

from __future__ import annotations

from dataclasses import dataclass
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
from mpmath import mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PRIOR_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_interval_one_step_current_gate.py"
)
PRIOR_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_interval_one_step_current_gate.json"
)

A = 159_577
L = 2_481_422
K = L - 1


@dataclass(frozen=True)
class Settings:
    precision_bits: int
    cutoff: int
    abs_tol: str
    deg_limit: int = 64
    eval_limit: int = 500_000
    depth_limit: int = 64


@dataclass(frozen=True)
class RationalBox:
    lo: Fraction
    hi: Fraction

    def __post_init__(self) -> None:
        if self.lo > self.hi:
            raise ValueError("reversed rational box")

    @property
    def mid(self) -> Fraction:
        return (self.lo + self.hi) / 2

    @property
    def radius(self) -> Fraction:
        return (self.hi - self.lo) / 2


PRODUCTION = Settings(precision_bits=320, cutoff=10, abs_tol="1e-55")
TARGET_BOXES = (
    {
        "name": "interior_lower_endpoint_neighbourhood",
        "z": RationalBox(Fraction(-503, 3000), Fraction(-497, 3000)),
        "tau": RationalBox(Fraction(399, 1000), Fraction(401, 1000)),
        "signed_tau": 1,
    },
    {
        "name": "interior_upper_endpoint_neighbourhood",
        "z": RationalBox(Fraction(-11, 30) - Fraction(1, 1000), Fraction(-11, 30) + Fraction(1, 1000)),
        "tau": RationalBox(Fraction(399, 1000), Fraction(401, 1000)),
        "signed_tau": -1,
    },
    {
        "name": "half_integer_recurrence_wall_neighbourhood",
        "z": RationalBox(Fraction(499, 1000), Fraction(501, 1000)),
        "tau": RationalBox(Fraction(499, 1000), Fraction(501, 1000)),
        "signed_tau": -1,
    },
)

PHYSICAL_ENDPOINT_BOXES = (
    {
        "name": "interior_physical_endpoint_box",
        "x": RationalBox(Fraction(2, 5) - Fraction(1, 10**8), Fraction(2, 5) + Fraction(1, 10**8)),
        "s": RationalBox(Fraction(1, 3) - Fraction(1, 10**4), Fraction(1, 3) + Fraction(1, 10**4)),
    },
    {
        "name": "physical_corner_recurrence_wall_box",
        "x": RationalBox(Fraction(1, 2) - Fraction(1, 10**9), Fraction(1, 2)),
        "s": RationalBox(Fraction(0), Fraction(1, 10**4)),
    },
)

LOCAL_WIDTH = Fraction(1, 10**28)
COMPLETE_BOXES = (
    {
        "name": "interior_complete_current_microbox",
        "center_x": Fraction(2, 5),
        "center_s": Fraction(1, 3),
        "x": RationalBox(Fraction(2, 5) - LOCAL_WIDTH, Fraction(2, 5) + LOCAL_WIDTH),
        "s": RationalBox(Fraction(1, 3) - LOCAL_WIDTH, Fraction(1, 3) + LOCAL_WIDTH),
    },
    {
        "name": "physical_corner_complete_current_microbox",
        "center_x": Fraction(1, 2),
        "center_s": Fraction(0),
        "x": RationalBox(Fraction(1, 2) - LOCAL_WIDTH, Fraction(1, 2)),
        "s": RationalBox(Fraction(0), LOCAL_WIDTH),
    },
)


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


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def nearest_shift(value: Fraction) -> tuple[int, Fraction]:
    shift = floor_fraction(value + Fraction(1, 2))
    return shift, value - shift


def arb_rational(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def arb_box(value: RationalBox) -> arb:
    return arb(arb_rational(value.mid), arb_rational(value.radius))


def box_record(value: RationalBox) -> dict[str, str]:
    return {"lo": str(value.lo), "hi": str(value.hi), "mid": str(value.mid), "radius": str(value.radius)}


def arb_upper_text(value: arb, digits: int = 40) -> str:
    midpoint, radius, exponent = value.mid_rad_10exp()
    upper = mp.mpf(abs(int(midpoint)) + int(radius)) * mp.power(10, int(exponent))
    return mp.nstr(upper, digits, min_fixed=-10, max_fixed=10)


def acb_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": str(value.real),
        "imag_ball": str(value.imag),
        "radius_upper": arb_upper_text(value.rad()),
        "absolute_upper": arb_upper_text(value.abs_upper()),
    }


def symmetric_complex_error(radius: arb) -> acb:
    upper = radius.abs_upper()
    return acb(arb(0, upper), arb(0, upper))


def cis_pi_interval(phase: arb) -> acb:
    return (acb(0, 1) * acb.pi() * acb(phase)).exp()


def central_mordell_parameter_box(
    z_box: RationalBox,
    tau_box: RationalBox,
    settings: Settings,
    prior: Any,
) -> tuple[acb, acb, dict[str, Any]]:
    """Enclose h and h_z simultaneously over a central (z,tau) box."""

    require(Fraction(-1, 2) <= z_box.lo <= z_box.hi <= Fraction(1, 2), "central z box escaped")
    require(tau_box.lo > 0, "tau box must be positive")
    flint.ctx.prec = settings.precision_bits
    pi = acb.pi()
    theta = (acb(0, 1) * pi / 4).exp()
    zz = acb(arb_box(z_box))
    tt = acb(arb_box(tau_box))
    tolerance = arb(settings.abs_tol)

    def common(y: acb) -> acb:
        return (-pi * tt * y * y).exp() / (pi * theta * y).cosh()

    def h_integrand(y: acb, analytic: bool) -> acb:
        return common(y) * (2 * pi * zz * theta * y).cosh()

    def hz_integrand(y: acb, analytic: bool) -> acb:
        return y * common(y) * (2 * pi * zz * theta * y).sinh()

    started = time.time()
    h_finite = 2 * theta * acb.integral(
        h_integrand,
        acb(0),
        acb(settings.cutoff),
        abs_tol=tolerance,
        rel_tol=tolerance,
        deg_limit=settings.deg_limit,
        eval_limit=settings.eval_limit,
        depth_limit=settings.depth_limit,
    )
    hz_finite = 4 * pi * acb(0, 1) * acb.integral(
        hz_integrand,
        acb(0),
        acb(settings.cutoff),
        abs_tol=tolerance,
        rel_tol=tolerance,
        deg_limit=settings.deg_limit,
        eval_limit=settings.eval_limit,
        depth_limit=settings.depth_limit,
    )
    require(h_finite.is_finite() and hz_finite.is_finite(), "parameter-box integration failed")
    h_tail, hz_tail, tail = prior.central_tail_bounds(tau_box.lo, settings.cutoff)
    h = h_finite + symmetric_complex_error(h_tail)
    hz = hz_finite + symmetric_complex_error(hz_tail)
    return h, hz, {
        "z_box": box_record(z_box),
        "tau_box": box_record(tau_box),
        "settings": settings.__dict__,
        "finite_h": acb_record(h_finite),
        "finite_h_z": acb_record(hz_finite),
        "tail_at_tau_lower": tail,
        "h": acb_record(h),
        "h_z": acb_record(hz),
        "elapsed_seconds": round(time.time() - started, 6),
    }


def split_at_half_integers(target: RationalBox) -> list[tuple[RationalBox, int]]:
    if target.lo == target.hi:
        shift, _ = nearest_shift(target.lo)
        return [(target, shift)]
    points = [target.lo]
    index = floor_fraction(target.lo - Fraction(1, 2)) + 1
    wall = Fraction(index) + Fraction(1, 2)
    while wall < target.hi:
        points.append(wall)
        index += 1
        wall = Fraction(index) + Fraction(1, 2)
    points.append(target.hi)
    pieces: list[tuple[RationalBox, int]] = []
    for lo, hi in zip(points, points[1:]):
        midpoint = (lo + hi) / 2
        shift, _ = nearest_shift(midpoint)
        pieces.append((RationalBox(lo, hi), shift))
    return pieces


def shifted_mordell_branch(
    target: RationalBox,
    shift: int,
    tau_box: RationalBox,
    signed_tau: int,
    settings: Settings,
    prior: Any,
) -> tuple[acb, acb, dict[str, Any]]:
    central = RationalBox(target.lo - shift, target.hi - shift)
    h, hz, central_record = central_mordell_parameter_box(central, tau_box, settings, prior)
    flint.ctx.prec = settings.precision_bits
    tt = acb(arb_box(tau_box))
    y = arb_box(central)
    rows = []
    if shift > 0:
        for _ in range(shift):
            phase = arb_rational(Fraction(1, 4)) + (y + arb_rational(Fraction(1, 2))) ** 2 / tt.real
            g = 2 * cis_pi_interval(phase) / tt.sqrt()
            gz = g * 2 * acb.pi() * acb(0, 1) * acb(y + arb_rational(Fraction(1, 2))) / tt
            h, hz = g - h, gz - hz
            rows.append({"direction": 1, "from": str(y), "phase": str(phase)})
            y += 1
    elif shift < 0:
        for _ in range(-shift):
            previous = y - 1
            phase = arb_rational(Fraction(1, 4)) + (previous + arb_rational(Fraction(1, 2))) ** 2 / tt.real
            g = 2 * cis_pi_interval(phase) / tt.sqrt()
            gz = g * 2 * acb.pi() * acb(0, 1) * acb(previous + arb_rational(Fraction(1, 2))) / tt
            h, hz = g - h, gz - hz
            rows.append({"direction": -1, "from": str(y), "phase": str(phase)})
            y = previous
    if signed_tau < 0:
        h = h.conjugate()
        hz = hz.conjugate()
    return h, hz, {
        "target_box": box_record(target),
        "central_box": box_record(central),
        "shift": shift,
        "signed_tau": signed_tau,
        "central_enclosure": central_record,
        "unit_recurrence": rows,
        "h": acb_record(h),
        "h_z": acb_record(hz),
    }


def mordell_target_parameter_box(
    target: RationalBox,
    tau_box: RationalBox,
    signed_tau: int,
    settings: Settings,
    prior: Any,
) -> tuple[acb, acb, dict[str, Any]]:
    require(signed_tau in (-1, 1), "signed_tau must be a sign")
    branches = []
    h_hull: acb | None = None
    hz_hull: acb | None = None
    for piece, shift in split_at_half_integers(target):
        h, hz, record = shifted_mordell_branch(piece, shift, tau_box, signed_tau, settings, prior)
        h_hull = h if h_hull is None else h_hull.union(h)
        hz_hull = hz if hz_hull is None else hz_hull.union(hz)
        branches.append(record)
    require(h_hull is not None and hz_hull is not None, "empty Mordell branch partition")
    return h_hull, hz_hull, {
        "target_box": box_record(target),
        "tau_box": box_record(tau_box),
        "signed_tau": signed_tau,
        "branch_count": len(branches),
        "branches": branches,
        "h": acb_record(h_hull),
        "h_z": acb_record(hz_hull),
    }


def physical_geometry(x_box: RationalBox, s_box: RationalBox) -> dict[str, Any]:
    center_x = x_box.mid
    center_s = s_box.mid
    center_a = center_x * (A + 2 * center_s) / 2
    shift, _ = nearest_shift(center_a)
    a_box = RationalBox(
        x_box.lo * (A + 2 * s_box.lo) / 2,
        x_box.hi * (A + 2 * s_box.hi) / 2,
    )
    require(a_box.lo >= shift - Fraction(1, 2), "nearest-shift lower wall crossed")
    require(a_box.hi < shift + Fraction(1, 2), "nearest-shift upper wall crossed")
    m = floor_fraction(K * center_x)
    require(K * x_box.lo >= m and K * x_box.hi < m + 1, "transformed-index wall crossed")
    z_box = RationalBox(a_box.lo - shift, a_box.hi - shift)
    lower = RationalBox(
        x_box.lo * (A + 2 * s_box.lo - 1) / 2 - shift + Fraction(1, 2),
        x_box.hi * (A + 2 * s_box.hi - 1) / 2 - shift + Fraction(1, 2),
    )
    upper = RationalBox(
        x_box.lo * (A + 2 * s_box.lo + 2 * K + 1) / 2 - shift - m - Fraction(1, 2),
        x_box.hi * (A + 2 * s_box.hi + 2 * K + 1) / 2 - shift - m - Fraction(1, 2),
    )
    return {
        "integer_shift_r": shift,
        "transformed_m": m,
        "a_box": a_box,
        "z_box": z_box,
        "lower_argument_box": lower,
        "upper_argument_box": upper,
        "nearest_shift_lower_margin": a_box.lo - (shift - Fraction(1, 2)),
        "nearest_shift_upper_margin": shift + Fraction(1, 2) - a_box.hi,
        "m_lower_margin": K * x_box.lo - m,
        "m_upper_margin": m + 1 - K * x_box.hi,
    }


def physical_endpoint_parameter_box(
    name: str,
    x_box: RationalBox,
    s_box: RationalBox,
    settings: Settings,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    geometry = physical_geometry(x_box, s_box)
    r = geometry["integer_shift_r"]
    m = geometry["transformed_m"]
    tau_box = RationalBox(x_box.lo / 2, x_box.hi / 2)
    mordell_tau_box = x_box
    lower_h, lower_hz, lower_record = mordell_target_parameter_box(
        geometry["lower_argument_box"], mordell_tau_box, -1, settings, prior
    )
    upper_h, upper_hz, upper_record = mordell_target_parameter_box(
        geometry["upper_argument_box"], mordell_tau_box, -1, settings, prior
    )

    flint.ctx.prec = settings.precision_bits
    xx = arb_box(x_box)
    ss = arb_box(s_box)
    cc = arb(A) + 2 * ss
    tau = xx / 2
    zz = xx * cc / 2 - r
    lower_phase = -(zz - tau / 2)
    upper_phase = 2 * (arb(K) + arb_rational(Fraction(1, 2))) * (
        zz + tau * (arb(K) + arb_rational(Fraction(1, 2)))
    )
    outer_phase = xx * cc * cc / 4
    pi_i = acb.pi() * acb(0, 1)
    lower_joined = (
        -acb(0, 1)
        * cis_pi_interval(lower_phase)
        * (acb(cc - 1) * lower_h + lower_hz / pi_i)
        / 2
    )
    upper_joined = (
        -acb(0, 1)
        * ((-1) ** m)
        * cis_pi_interval(upper_phase)
        * (acb(cc + 2 * K + 1) * upper_h + upper_hz / pi_i)
        / 2
    )
    outer = cis_pi_interval(outer_phase)
    lower_current = outer * lower_joined
    upper_current = outer * upper_joined
    endpoint = lower_current + upper_current
    return endpoint, {
        "name": name,
        "x_box": box_record(x_box),
        "s_box": box_record(s_box),
        "integer_shift_r": r,
        "transformed_m": m,
        "a_box": box_record(geometry["a_box"]),
        "z_box": box_record(geometry["z_box"]),
        "lower_argument_box": box_record(geometry["lower_argument_box"]),
        "upper_argument_box": box_record(geometry["upper_argument_box"]),
        "nearest_shift_lower_margin": str(geometry["nearest_shift_lower_margin"]),
        "nearest_shift_upper_margin": str(geometry["nearest_shift_upper_margin"]),
        "m_lower_margin": str(geometry["m_lower_margin"]),
        "m_upper_margin": str(geometry["m_upper_margin"]),
        "lower_Mordell": lower_record,
        "upper_Mordell": upper_record,
        "lower_endpoint_current": acb_record(lower_current),
        "upper_endpoint_current": acb_record(upper_current),
        "joined_endpoint_current": acb_record(endpoint),
    }


def polynomial_masses(r: int, m: int) -> tuple[int, int, int]:
    s0 = (m + 1) * r + m * (m + 1) // 2
    s1 = r * m * (m + 1) // 2 + m * (m + 1) * (2 * m + 1) // 6
    s2 = r * m * (m + 1) * (2 * m + 1) // 6 + (m * (m + 1) // 2) ** 2
    return s0, s1, s2


def main_variation_data(
    center_x: Fraction,
    center_s: Fraction,
    x_box: RationalBox,
    s_box: RationalBox,
    prior: Any,
) -> tuple[acb, arb, dict[str, Any]]:
    geometry = physical_geometry(x_box, s_box)
    r = geometry["integer_shift_r"]
    m = geometry["transformed_m"]
    c0 = Fraction(A) + 2 * center_s
    _, z0 = nearest_shift(center_x * c0 / 2)
    w0 = c0 / 2 - Fraction(r) / center_x
    sigma0 = -Fraction(1, 2) / center_x
    transformed0, period = prior.transformed_weighted_ball(w0, sigma0, m, r)
    phase0 = Fraction(1, 4) + r * c0 - Fraction(r * r) / center_x

    c_lo = Fraction(A) + 2 * s_box.lo
    c_hi = Fraction(A) + 2 * s_box.hi
    w_box = RationalBox(c_lo / 2 - Fraction(r) / x_box.lo, c_hi / 2 - Fraction(r) / x_box.hi)
    sigma_box = RationalBox(-Fraction(1, 2) / x_box.lo, -Fraction(1, 2) / x_box.hi)
    phase_box = RationalBox(
        Fraction(1, 4) + r * c_lo - Fraction(r * r) / x_box.lo,
        Fraction(1, 4) + r * c_hi - Fraction(r * r) / x_box.hi,
    )
    delta_w = max(abs(w_box.lo - w0), abs(w_box.hi - w0))
    delta_sigma = max(abs(sigma_box.lo - sigma0), abs(sigma_box.hi - sigma0))
    delta_phase = max(abs(phase_box.lo - phase0), abs(phase_box.hi - phase0))
    s0, s1, s2 = polynomial_masses(r, m)

    flint.ctx.prec = max(flint.ctx.prec, 320)
    raw_transform_error = 2 * arb.pi() * arb_rational(delta_w * s1 + delta_sigma * s2)
    cap = arb(2 * s0)
    transform_error = raw_transform_error if bool(raw_transform_error < cap) else cap
    q0 = 2 / (arb_rational(center_x) * arb_rational(center_x).sqrt())
    q_box = 2 / (arb_box(x_box) * arb_box(x_box).sqrt())
    prefactor_error = (q_box - q0).abs_upper() + q0.abs_upper() * arb.pi() * arb_rational(delta_phase)
    prefactor_upper = q_box.abs_upper()
    main_error = prefactor_upper * transform_error + prefactor_error * transformed0.abs_upper()
    main0 = acb(q0) * prior.acb_cis_pi(phase0) * transformed0
    main_box = main0 + symmetric_complex_error(main_error)
    return main_box, main_error, {
        "center_x": str(center_x),
        "center_s": str(center_s),
        "center_z": str(z0),
        "center_w": str(w0),
        "center_sigma": str(sigma0),
        "center_phase": str(phase0),
        "transformed_period_at_center": period,
        "w_box": box_record(w_box),
        "sigma_box": box_record(sigma_box),
        "phase_box": box_record(phase_box),
        "delta_w": str(delta_w),
        "delta_sigma": str(delta_sigma),
        "delta_phase": str(delta_phase),
        "S0": str(s0),
        "S1": str(s1),
        "S2": str(s2),
        "transform_first_difference_bound": arb_upper_text(transform_error),
        "prefactor_difference_bound": arb_upper_text(prefactor_error),
        "main_current_variation_bound": arb_upper_text(main_error),
        "center_transformed_current": acb_record(transformed0),
        "center_main_current": acb_record(main0),
        "main_current_box": acb_record(main_box),
    }


def complete_current_parameter_box(case: dict[str, Any], settings: Settings, prior: Any) -> dict[str, Any]:
    endpoint, endpoint_record = physical_endpoint_parameter_box(
        case["name"], case["x"], case["s"], settings, prior
    )
    main, main_error, main_record = main_variation_data(
        case["center_x"], case["center_s"], case["x"], case["s"], prior
    )
    complete = main + endpoint
    center_source, source_period = prior.source_current_ball(case["center_x"], case["center_s"])
    require(complete.contains(center_source), f"complete box misses its exact center source: {case['name']}")
    source_mass = arb(L) * acb(arb_rational(Fraction(A) + 2 * case["center_s"] + L - 1)).real
    normalized_radius = complete.rad() / source_mass
    return {
        "name": case["name"],
        "x_box": box_record(case["x"]),
        "s_box": box_record(case["s"]),
        "center_x": str(case["center_x"]),
        "center_s": str(case["center_s"]),
        "endpoint": endpoint_record,
        "main_variation": main_record,
        "complete_current_box": acb_record(complete),
        "center_source_current": acb_record(center_source),
        "center_source_period": source_period,
        "contains_exact_center_source": True,
        "complete_radius_absolute_upper": arb_upper_text(complete.rad()),
        "complete_radius_over_center_source_mass_upper": arb_upper_text(normalized_radius),
        "main_variation_bound": arb_upper_text(main_error),
        "passed": True,
    }


def target_box_probe_record(spec: dict[str, Any], settings: Settings, prior: Any) -> dict[str, Any]:
    h, hz, enclosure = mordell_target_parameter_box(
        spec["z"], spec["tau"], spec["signed_tau"], settings, prior
    )
    z_values = (spec["z"].lo, spec["z"].mid, spec["z"].hi)
    tau_values = (spec["tau"].lo, spec["tau"].mid, spec["tau"].hi)
    probes = []
    for z in z_values:
        for tau in tau_values:
            point_h, point_hz, _ = prior.mordell_h_and_derivative_balls(
                z, spec["signed_tau"] * tau, settings
            )
            require(h.contains(point_h), f"h parameter box misses point probe: {spec['name']}")
            require(hz.contains(point_hz), f"h_z parameter box misses point probe: {spec['name']}")
            probes.append({"z": str(z), "tau": str(tau), "h_contained": True, "h_z_contained": True})
    return {
        "name": spec["name"],
        "enclosure": enclosure,
        "probe_count": len(probes),
        "point_probes": probes,
        "passed": True,
    }


def width_box(center_x: Fraction, center_s: Fraction, width: Fraction) -> tuple[RationalBox, RationalBox]:
    if center_x == Fraction(1, 2):
        x_box = RationalBox(center_x - width, center_x)
    else:
        x_box = RationalBox(center_x - width, center_x + width)
    if center_s == 0:
        s_box = RationalBox(Fraction(0), width)
    else:
        s_box = RationalBox(center_s - width, center_s + width)
    return x_box, s_box


def scale_audit(prior: Any) -> list[dict[str, Any]]:
    rows = []
    for center_x, center_s in ((Fraction(1, 2), Fraction(0)), (Fraction(2, 5), Fraction(1, 3))):
        for exponent in (12, 16, 20, 24, 28):
            width = Fraction(1, 10**exponent)
            x_box, s_box = width_box(center_x, center_s, width)
            _, error, record = main_variation_data(center_x, center_s, x_box, s_box, prior)
            rows.append({
                "center_x": str(center_x),
                "center_s": str(center_s),
                "half_width": str(width),
                "half_width_decimal": f"1e-{exponent}",
                "main_variation_bound": arb_upper_text(error),
                "delta_w": record["delta_w"],
                "delta_sigma": record["delta_sigma"],
                "delta_phase": record["delta_phase"],
            })
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    summary = artifact["summary"]
    return f"""# Mordell parameter boxes and direct phase-scale audit

Date: 2026-08-24

Status: interval-certificate and method-barrier note; central Mordell boxes,
recurrence-wall splitting, and two microscopic complete-current boxes are
certified, but scalable recursion and physical quadrature remain open

## Parameter boxes

Let `Z=[z_-,z_+]` lie in `[-1/2,1/2]` and let
`T=[tau_-,tau_+]` have `tau_->0`.  Substituting the Arb balls `Z,T` directly
in Kuznetsov's rotated integral encloses `h(Z,T)` and `h_z(Z,T)` on `[0,Y]`.
The Section 11.464 tails are evaluated at `tau_-`, where they are largest:

```text
T_h <=coth(pi Y/sqrt(2))*erfc(sqrt(pi tau_-)Y)/sqrt(tau_-),
T_hz<=2*coth(pi Y/sqrt(2))*exp(-pi tau_- Y^2)/tau_-. (PB1)
```

Target boxes crossing a half integer are split there.  On each piece one
fixed exact unit recurrence transports the central enclosure; the branch
balls are then joined by interval hull.  Thus a recurrence wall is covered,
not sampled around or silently assigned to one side.                    (PB2)

Production at {artifact['production_settings']['precision_bits']} bits,
`Y={artifact['production_settings']['cutoff']}`, and tolerance
`{artifact['production_settings']['abs_tol']}` certifies
{summary['target_box_count']} target boxes, including one crossing `z=1/2`.
All {summary['dense_point_probe_count']} rigorous pointwise probes are
contained in their parent boxes.  The largest target-box radii are

```text
h:   {summary['maximum_h_box_radius_upper']},
h_z: {summary['maximum_h_z_box_radius_upper']}.       (PB3)
```

The physical map also certifies two joined endpoint-current boxes while
checking that `r=nearest(x(A+2s)/2)` and `m=floor(Kx)` remain fixed.  The
corner box crosses the Mordell recurrence wall and is split automatically.

## Complete-current microboxes

For fixed `r,m`, write

```text
T(w,sigma)=sum_(k=0)^m (r+k)exp(2pi i(wk+sigma k^2)).
```

Around an exact rational centre `(w_0,sigma_0)`, the elementary chord bound
`|exp(iu)-exp(iv)|<=|u-v|` gives

```text
|T-T_0|<=2pi(delta_w S_1+delta_sigma S_2),
S_1=sum_(k=0)^m(r+k)k,
S_2=sum_(k=0)^m(r+k)k^2.                             (PB4)
```

The centre `T_0` is enclosed by exact rational-period compression.  The
prefactor is treated as `q exp(i pi phi)`, where

```text
q=2/x^(3/2),  w=(A+2s)/2-r/x,
sigma=-1/(2x),  phi=1/4+r(A+2s)-r^2/x.              (PB5)
```

Combining (PB4) with interval endpoint currents preserves the exact centre
cancellation and bounds only its variation.  At half-width `1e-28`, both
full-roster boxes contain their independently enclosed exact centre source
currents.  Their largest absolute complete radius is
`{summary['maximum_complete_microbox_radius_upper']}`.

## Scaling decision

The ten-row width audit applies exactly the same proved first-difference
bound at half-widths `1e-12,1e-16,1e-20,1e-24,1e-28`.  Its largest and
smallest main-variation bounds are respectively
`{summary['maximum_scale_audit_main_variation_upper']}` and
`{summary['minimum_scale_audit_main_variation_upper']}`.  The `1e-24` rows
already have main-variation bounds above `{summary['minimum_1e_minus_24_main_variation_upper']}`.

This is a barrier for this direct first-difference box strategy, not an
impossibility theorem for all interval methods.  It shows that tiling the
physical domain with direct source-length boxes is not a credible next step.
The next implementation must recurse the truncated theta current, reduce its
active length, and propagate boxes through the transformed main and joined
endpoint currents before attempting quadrature.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Primary-source boundary: the rotated integral and exact unit recurrence are
from Alexey Kuznetsov, *Computing the truncated theta function via Mordell
integral*, arXiv:1306.4081v2.  The box substitution, wall partition, and
first-difference phase audit are derived here.  Practical Gauss--Laguerre
quadrature is not used.

Pi provenance: every `pi` in this gate comes from the inherited Fourier phase
or rotated Mordell representation.  No fitted circle or polygon constant is
introduced.

This gate proves three central/transported Mordell parameter boxes, two
physical joined endpoint boxes, and two microscopic complete-current boxes.
It does not prove a scalable parameter-box current evaluator, the small-`tau`
branch, recursive interval control, physical quadrature, the non-A bound,
joined `R_after_A`, `R_Dir`, `Q_K-T`, an all-height theorem, `Lambda<=0`,
PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(PRIOR_RESULT.is_file() and PRIOR_BUILDER.is_file(), "prior Mordell certificate missing")
    prior_artifact = json.loads(PRIOR_RESULT.read_text(encoding="utf-8"))
    require(prior_artifact.get("passed") is True, "prior Mordell certificate not passed")
    prior = load_module("mordell_point_interval_dependency", PRIOR_BUILDER)

    target_rows = [target_box_probe_record(spec, PRODUCTION, prior) for spec in TARGET_BOXES]
    endpoint_rows = []
    for spec in PHYSICAL_ENDPOINT_BOXES:
        endpoint, record = physical_endpoint_parameter_box(
            spec["name"], spec["x"], spec["s"], PRODUCTION, prior
        )
        require(endpoint.is_finite(), f"nonfinite physical endpoint box: {spec['name']}")
        endpoint_rows.append(record)
    complete_rows = [complete_current_parameter_box(case, PRODUCTION, prior) for case in COMPLETE_BOXES]
    audit_rows = scale_audit(prior)

    h_radii = [arb(row["enclosure"]["h"]["radius_upper"]) for row in target_rows]
    hz_radii = [arb(row["enclosure"]["h_z"]["radius_upper"]) for row in target_rows]
    complete_radii = [arb(row["complete_radius_absolute_upper"]) for row in complete_rows]
    audit_errors = [arb(row["main_variation_bound"]) for row in audit_rows]
    audit_24 = [arb(row["main_variation_bound"]) for row in audit_rows if row["half_width_decimal"] == "1e-24"]
    require(all(bool(radius < arb("0.02")) for radius in complete_radii), "microbox radius target failed")
    require(any(row["enclosure"]["branch_count"] == 2 for row in target_rows), "recurrence wall was not split")
    require(any(row["lower_Mordell"]["branch_count"] > 1 or row["upper_Mordell"]["branch_count"] > 1 for row in endpoint_rows), "physical wall box was not split")

    artifact = {
        "kind": STEM,
        "status": "Mordell_parameter_boxes_recurrence_wall_split_and_two_complete_microboxes_certified_direct_phase_scale_barrier_recorded_recursion_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "L": L,
            "K": K,
            "workers": 1,
            "complete_microbox_half_width": str(LOCAL_WIDTH),
        },
        "production_settings": PRODUCTION.__dict__,
        "target_parameter_boxes": target_rows,
        "physical_endpoint_boxes": endpoint_rows,
        "complete_current_microboxes": complete_rows,
        "direct_phase_scale_audit": audit_rows,
        "summary": {
            "target_box_count": len(target_rows),
            "dense_point_probe_count": sum(row["probe_count"] for row in target_rows),
            "physical_endpoint_box_count": len(endpoint_rows),
            "complete_microbox_count": len(complete_rows),
            "maximum_h_box_radius_upper": arb_upper_text(max(h_radii)),
            "maximum_h_z_box_radius_upper": arb_upper_text(max(hz_radii)),
            "maximum_complete_microbox_radius_upper": arb_upper_text(max(complete_radii)),
            "maximum_scale_audit_main_variation_upper": arb_upper_text(max(audit_errors)),
            "minimum_scale_audit_main_variation_upper": arb_upper_text(min(audit_errors)),
            "minimum_1e_minus_24_main_variation_upper": arb_upper_text(min(audit_24)),
        },
        "decision": {
            "uniform_central_Mordell_parameter_boxes_built": True,
            "half_integer_recurrence_wall_split_exactly": True,
            "physical_r_and_m_branch_guards_enforced": True,
            "joined_endpoint_parameter_boxes_built": True,
            "two_complete_current_microboxes_certified": True,
            "direct_source_length_first_difference_boxes_scale_to_physical_quadrature": False,
            "small_tau_branch_built": False,
            "recursive_interval_theta_evaluator_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Implement one recursive interval Mordell-current step that contracts the active theta length, "
            "including branch-stable parameter boxes and the joined endpoint current; then derive the tied "
            "small-tau branch before any physical tiling."
        ),
        "proof_boundary": (
            "Three Mordell target parameter boxes, exact half-integer wall splitting, two physical joined "
            "endpoint boxes, and two width-1e-28 complete-current boxes only. The direct phase audit is a "
            "barrier for the recorded first-difference strategy, not every possible interval method. No "
            "scalable current boxes, small-tau branch, recursive evaluator, physical quadrature, non-A bound, "
            "joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level "
            "conclusion is proved."
        ),
        "primary_source": prior_artifact["primary_source"],
        "dependencies": {
            "prior_result": {"path": relative(PRIOR_RESULT), "sha256": file_hash(PRIOR_RESULT)},
            "prior_builder": {"path": relative(PRIOR_BUILDER), "sha256": file_hash(PRIOR_BUILDER)},
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
    print(
        f"certified {len(target_rows)} Mordell parameter boxes and "
        f"{len(complete_rows)} complete-current microboxes; direct phase scaling rejects tiling",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
