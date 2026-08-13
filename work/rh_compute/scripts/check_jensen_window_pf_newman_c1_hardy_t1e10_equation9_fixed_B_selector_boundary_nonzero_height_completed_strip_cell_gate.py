#!/usr/bin/env python3
"""Independently validate the nonzero completed-strip selector cell."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
SCRIPT_ROOT = REPO_ROOT / "work/rh_compute/scripts"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate as point_gate


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.json"
PRECISION = 100
CHECK_RADIUS = 10
CHECK_PANELS = 40
MAJORANT_PANELS = 1200


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def interval_ball(left: arb, right: arb) -> arb:
    return arb((left + right) / 2, (right - left) / 2)


def derivative(beta: arb, pi: arb) -> arb:
    i = acb(0, 1)
    rotation = (i * pi / 6).exp()

    def integrand(radius: acb, _: bool) -> acb:
        z = rotation * radius
        scaled = z / beta
        exact = scaled.cosh() ** (-arb(3) / 2) * (i * beta**3 * point_gate.argument_minus_tanh_polynomial(scaled)).exp()
        return rotation * (-i * z) * (exact - (-radius**3 / 3).exp())

    ray = acb(0)
    radius = arb(CHECK_RADIUS)
    for index in range(CHECK_PANELS):
        ray += acb.integral(
            integrand,
            radius * index / CHECK_PANELS,
            radius * (index + 1) / CHECK_PANELS,
            abs_tol=arb("1e-29"),
            rel_tol=arb("1e-29"),
            eval_limit=400_000,
            depth_limit=50,
        )
    rho = radius / beta
    phase_tail = beta**3 * 3 * rho**18 / (1 - rho)
    replacement = radius**2 * phase_tail * phase_tail.exp()
    exact_tail = 2 * (8 * (-radius**3 / 4).exp() / (3 * radius) + (-beta**3 / 5).exp() * (5 / beta + 25 / beta**4))
    canonical_tail = 2 * (-radius**3 / 3).exp() / radius
    return 2 * ray.real + arb(0, replacement + exact_tail + canonical_tail)


def second_derivative_bound(beta: arb, lambda_radius: arb) -> arb:
    radius = arb(CHECK_RADIUS)
    compact = arb(0)
    for index in range(MAJORANT_PANELS):
        left = radius * index / MAJORANT_PANELS
        right = radius * (index + 1) / MAJORANT_PANELS
        r = interval_ball(left, right)
        rho = r / beta
        w = rho**2 * rho.cosh() / 2
        amplitude_error = arb(3) / 2 * w / (1 - w) ** (arb(5) / 2)
        phase_error = arb(0)
        for degree, coefficient in point_gate.TANH_COEFFICIENTS.items():
            if degree >= 5:
                phase_error += abs(arb(coefficient.numerator) / coefficient.denominator) * rho**degree
        phase_error = beta**3 * (phase_error + 3 * rho**18 / (1 - rho))
        majorant = (-r**3 / 3 + lambda_radius * r / 2).exp() * (amplitude_error + phase_error) * phase_error.exp()
        compact += (right - left) * (r**2 * majorant).upper()
    re = (10 * lambda_radius / 3).sqrt()
    ce = lambda_radius * re / 2 - re**3 / 20
    rc = (2 * lambda_radius).sqrt()
    cc = lambda_radius * rc / 2 - rc**3 / 12
    exact_tail = ce.exp() * arb(5) / 3 * (-radius**3 / 5).exp()
    canonical_tail = cc.exp() * arb(4) / 3 * (-radius**3 / 4).exp()
    rate = beta**2 / 5 - lambda_radius / 2
    far = (-rate * beta).exp() * (beta**2 / rate + 2 * beta / rate**2 + 2 / rate**3)
    return 2 * (compact + exact_tail + canonical_tail + far)


def main() -> None:
    priority = set_low_priority()
    require(RESULT.is_file(), "missing nonzero selector-cell result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "nonzero selector-cell artifact is not passed")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file() and file_hash(path) == record["sha256"], f"hash drift: {path}")

    ctx.dps = PRECISION
    p = point_gate.parameters()
    lambda_radius = p["pi"] / (16 * p["beta"])
    fresh_derivative = derivative(p["beta"], p["pi"])
    saved_derivative = arb(artifact["center_derivative"]["rigorous_derivative_ball"])
    require(fresh_derivative.overlaps(saved_derivative), "independent center derivative misses saved enclosure")
    formal = -arb(13) / 60 * 2 * p["pi"] * arb(0).airy_ai() / p["beta"] ** 2
    require(arb("0.999999999") < fresh_derivative / formal < arb("1.000000001"), "independent formal derivative comparison failed")

    fresh_second = second_derivative_bound(p["beta"], lambda_radius)
    saved_second = arb(artifact["second_derivative_majorant"]["two_ray_second_derivative_bound"])
    require(fresh_second <= saved_second, "independent finer majorant exceeds saved bound")
    point = json.loads(point_gate.RESULT.read_text(encoding="utf-8"))
    delta0 = arb(point["numerical_certificate"]["reduced_exact_minus_canonical_fold_ball"])
    fresh_uniform = abs(delta0) + lambda_radius * abs(fresh_derivative) + lambda_radius**2 * fresh_second / 2
    saved_uniform = arb(artifact["uniform_error"]["uniform_reduced_completed_strip_error_bound_ball"])
    require(fresh_uniform <= saved_uniform, "independent Taylor bound exceeds saved bound")
    require(fresh_uniform < arb("9.6e-12"), "independent reduced cell bound exceeds target")

    strip_scale = arb(2).sqrt() / p["pi"] * p["Y"]
    max_normalizer = (p["pi"] / (32 * (p["height"] - p["pi"] / 16))) ** (arb(1) / 4)
    require(max_normalizer * strip_scale * fresh_uniform < arb("8.9e-13"), "independent normalized cell bound exceeds target")
    require(artifact["event_geometry"]["adjacent_selector_first_event_mode"] == 40094, "adjacent event mode drift")
    require(artifact["decision"]["ordinary_Morse_join_proved"] is False, "ordinary-Morse join was overpromoted")
    print(f"validated nonzero completed-strip selector cell independently; priority={priority}")


if __name__ == "__main__":
    main()
