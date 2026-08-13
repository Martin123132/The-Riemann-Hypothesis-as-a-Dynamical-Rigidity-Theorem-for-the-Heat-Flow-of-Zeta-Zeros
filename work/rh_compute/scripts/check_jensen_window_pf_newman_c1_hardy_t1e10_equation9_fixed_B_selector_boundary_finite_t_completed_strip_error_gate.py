#!/usr/bin/env python3
"""Independently validate the completed finite-t selector-strip certificate."""

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


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.json"
PRECISION = 100
LOWER_ALPHA = 159_577
CHECK_RADIUS = 10
CHECK_PANELS = 40

COEFFICIENTS = {
    1: Fraction(1), 3: Fraction(-1, 3), 5: Fraction(2, 15),
    7: Fraction(-17, 315), 9: Fraction(62, 2835),
    11: Fraction(-1382, 155925), 13: Fraction(21844, 6081075),
    15: Fraction(-929569, 638512875), 17: Fraction(6404582, 10854718875),
}


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


def polynomial(value: acb) -> acb:
    total = acb(0)
    for degree, coefficient in COEFFICIENTS.items():
        total += arb(coefficient.numerator) / coefficient.denominator * value**degree
    return total


def argument_minus_polynomial(value: acb) -> acb:
    total = acb(0)
    for degree, coefficient in COEFFICIENTS.items():
        if degree >= 3:
            total -= arb(coefficient.numerator) / coefficient.denominator * value**degree
    return total


def independent_difference(beta: arb, pi: arb) -> arb:
    i = acb(0, 1)
    rotation = (i * pi / 6).exp()

    def integrand(radius: acb, _: bool) -> acb:
        z = rotation * radius
        scaled = z / beta
        amplitude = scaled.cosh() ** (-arb(3) / 2)
        surrogate = amplitude * (i * beta**3 * argument_minus_polynomial(scaled)).exp()
        return rotation * (surrogate - (-radius**3 / 3).exp())

    total = acb(0)
    for index in range(CHECK_PANELS):
        left = arb(CHECK_RADIUS * index) / CHECK_PANELS
        right = arb(CHECK_RADIUS * (index + 1)) / CHECK_PANELS
        total += acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb("1e-30"),
            rel_tol=arb("1e-30"),
            eval_limit=400_000,
            depth_limit=50,
        )

    radius = arb(CHECK_RADIUS)
    rho = radius / beta
    polynomial_tail = 3 * rho**18 / (1 - rho)
    phase_tail = beta**3 * polynomial_tail
    replacement = 2 * radius * phase_tail * phase_tail.exp()
    exact_tail = 2 * (4 * (-radius**3 / 4).exp() / (3 * radius**2) + 5 * (-beta**3 / 5).exp() / beta**2)
    canonical_tail = 2 * (-radius**3 / 3).exp() / radius**2
    return 2 * total.real + arb(0, replacement + exact_tail + canonical_tail)


def check_decay_independently() -> None:
    root_three = arb(3).sqrt()

    def decay(rho: arb) -> arb:
        return rho / 2 - rho.sin() / ((root_three * rho).cosh() + rho.cos())

    for index in range(512):
        left = arb(1) / 2 + arb(index) / 1024
        right = arb(1) / 2 + arb(index + 1) / 1024
        rho = arb((left + right) / 2, (right - left) / 2)
        require((decay(rho) - rho**3 / 4).lower() > 0, f"independent cubic ray panel failed: {index}")
    for index in range(256):
        left = arb(1) + arb(index) / 1280
        right = arb(1) + arb(index + 1) / 1280
        rho = arb((left + right) / 2, (right - left) / 2)
        require((decay(rho) - rho / 5).lower() > 0, f"independent linear ray panel failed: {index}")
    endpoint = arb("1.2")
    require((3 * endpoint / 10 - 1 / ((root_three * endpoint).cosh() - 1)).lower() > 0, "independent analytic tail failed")


def main() -> None:
    priority = set_low_priority()
    require(RESULT.is_file(), "missing completed-strip result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "completed-strip artifact is not passed")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"source hash drift: {path}")

    ctx.dps = PRECISION
    pi = arb.pi()
    lower = arb(LOWER_ALPHA)
    height = pi * lower**2 / 8
    beta = height ** (arb(1) / 3)
    sigma = 4 * beta / (pi * lower)
    y_width = 2 / sigma
    h = 4 * beta / lower
    require((h * y_width - 2 * pi).contains(0), "independent hY identity failed")
    require(LOWER_ALPHA % 8 == 1 and ((LOWER_ALPHA**2 - 1) // 8) % 2 == 0, "independent source parity failed")
    check_decay_independently()

    fresh = independent_difference(beta, pi)
    saved = arb(artifact["numerical_certificate"]["reduced_exact_minus_canonical_fold_ball"])
    require(fresh.overlaps(saved), "independent radius-10 integral misses saved radius-12 enclosure")
    require(arb("1.2e-15") < fresh < arb("1.22e-15"), "independent completed correction leaves certified range")

    strip_scale = arb(2).sqrt() / pi * y_width
    normalizer = (pi / (32 * height)) ** (arb(1) / 4)
    normalized = normalizer * strip_scale * fresh
    saved_normalized = arb(artifact["numerical_certificate"]["equation9_normalized_finite_t_minus_canonical_ball"])
    require(normalized.overlaps(saved_normalized), "independent normalized correction misses saved enclosure")
    require(abs(normalized) < arb("1.2e-16"), "independent normalized correction exceeds target")

    leading = -arb(9) / 560 * 2 * pi * arb(0).airy_ai(derivative=1) / beta**4
    require(arb("0.999999") < fresh / leading < arb("1.000001"), "independent beta^-4 comparison failed")
    print(f"validated completed selector-strip finite-t error independently; priority={priority}")


if __name__ == "__main__":
    main()
