#!/usr/bin/env python3
"""Independently check the exact common-profile initial-slope gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
RADIUS = 11
PANELS = 22
PRECISION = 90


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval_ball(lower: arb, upper: arb) -> arb:
    return arb((lower + upper) / 2, (upper - lower) / 2)


def tanh_coefficients() -> dict[int, Fraction]:
    x = sp.symbols("x")
    series = sp.Poly(sp.series(sp.tanh(x), x, 0, 19).removeO(), x)
    result: dict[int, Fraction] = {}
    for (degree,), coefficient in series.terms():
        result[degree] = Fraction(int(sp.numer(coefficient)), int(sp.denom(coefficient)))
    return result


def polynomial(value: acb, coefficients: dict[int, Fraction]) -> acb:
    total = acb(0)
    for degree, coefficient in coefficients.items():
        total += arb(coefficient.numerator) / coefficient.denominator * value**degree
    return total


def argument_minus_tanh(value: acb, coefficients: dict[int, Fraction]) -> acb:
    total = acb(0)
    for degree, coefficient in coefficients.items():
        if degree >= 3:
            total -= arb(coefficient.numerator) / coefficient.denominator * value**degree
    return total


def beta4_multiplier(z: acb, epsilon: arb) -> acb:
    i = acb(0, 1)
    return (
        -i * z
        + epsilon * (arb(1) / 2 + arb(13) * i * z**3 / 12 - arb(2) * z**6 / 15)
        + epsilon**2 * (
            -arb(3) * z**2 / 8 - arb(137) * i * z**5 / 160
            + arb(25) * z**8 / 126 + arb(2) * i * z**11 / 225
        )
    )


def cubic_tail(radius: arb, power: int, rate: arb) -> arb:
    exponential = (-rate * radius**3).exp()
    if power <= 1:
        return radius ** (power - 2) * exponential / (3 * rate)
    return radius ** (power - 2) * exponential / (3 * rate) + arb(power - 2) / (3 * rate) * cubic_tail(radius, power - 3, rate)


def recompute() -> arb:
    ctx.dps = PRECISION
    pi = arb.pi()
    c = arb(159577)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    epsilon = beta**-2
    lambda_max = pi / (16 * beta)
    lam = interval_ball(arb(0), lambda_max)
    i = acb(0, 1)
    rotation = (i * pi / 6).exp()
    coefficients = tanh_coefficients()

    def integrand(radius: acb, _: bool) -> acb:
        z = rotation * radius
        scaled = z / beta
        amplitude = scaled.cosh() ** (-arb(3) / 2)
        ratio = amplitude * (i * (beta**3 * argument_minus_tanh(scaled, coefficients) - z**3 / 3)).exp()
        multiplier = arb(1) / (2 * beta**2) - i * beta * polynomial(scaled, coefficients)
        carrier = (i * (z**3 / 3 - lam * z)).exp()
        return rotation * carrier * (ratio * multiplier - beta4_multiplier(z, epsilon))

    radius = arb(RADIUS)
    ray = acb(0)
    for panel in range(PANELS):
        left = radius * panel / PANELS
        right = radius * (panel + 1) / PANELS
        ray += acb.integral(
            integrand, left, right, abs_tol=arb("1e-27"), rel_tol=arb("1e-27"),
            eval_limit=500_000, depth_limit=55,
        )

    rho = radius / beta
    tanh_tail = 3 * rho**18 / (1 - rho)
    phase_tail = beta**3 * tanh_tail
    tanh_bound = rho.sinh() / rho.cos()
    multiplier_bound = arb(1) / (2 * beta**2) + beta * tanh_bound
    replacement = 2 * radius * (lambda_max * radius / 2).exp() * phase_tail.exp() * (
        multiplier_bound * phase_tail * phase_tail.exp() + beta * tanh_tail
    )
    near_rate = radius**2 / 5
    near_exp = (-radius**3 / 5).exp()
    exact_near = 2 * (
        near_exp / near_rate / (2 * beta**2)
        + arb("2.2") * near_exp * (radius / near_rate + 1 / near_rate**2)
    )
    far_rate = beta**2 / 5 - lambda_max / 2
    exact_far = 2 * (arb(1) / (2 * beta**2) + 2 * beta) * (-far_rate * beta).exp() / far_rate
    coefficients4 = {
        0: epsilon / 2, 1: arb(1), 2: arb(3) * epsilon**2 / 8,
        3: arb(13) * epsilon / 12, 5: arb(137) * epsilon**2 / 160,
        6: arb(2) * epsilon / 15, 8: arb(25) * epsilon**2 / 126,
        11: arb(2) * epsilon**2 / 225,
    }
    canonical = 2 * sum(value * cubic_tail(radius, degree, arb(1) / 4) for degree, value in coefficients4.items())
    return 2 * ray.real + arb(0, replacement + exact_near + exact_far + canonical)


def main() -> None:
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    require(artifact["decision"]["uniform_initial_slope_below_5_56e_minus_22"] is True, "saved slope decision failed")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")

    saved = arb(artifact["numerical_certificate"]["uniform_exact_minus_beta4_initial_slope_ball"])
    independent = recompute()
    require(saved.overlaps(independent), "independently generated degree-17/radius-11 recomputation misses saved slope")
    require(abs(independent) < arb("5.56e-22"), "independent slope exceeds 5.56e-22")
    print("PASS: independently generated degree-17/radius-11 initial-slope recomputation overlaps the saved gate")


if __name__ == "__main__":
    main()
