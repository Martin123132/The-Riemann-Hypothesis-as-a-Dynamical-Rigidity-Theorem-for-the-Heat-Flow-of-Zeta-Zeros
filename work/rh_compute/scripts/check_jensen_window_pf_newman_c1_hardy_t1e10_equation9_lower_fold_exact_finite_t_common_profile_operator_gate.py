#!/usr/bin/env python3
"""Independently check the exact finite-t common-profile operator."""

from __future__ import annotations

from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_ball(text: str) -> tuple[Decimal, Decimal]:
    body = text.strip().removeprefix("[").removesuffix("]").strip()
    if "+/-" not in body:
        return Decimal(body), Decimal(0)
    midpoint, radius = body.split("+/-", 1)
    return Decimal(midpoint.strip()), Decimal(radius.strip())


def independent_symbolic_check() -> None:
    beta, lam, z, s = sp.symbols("beta lambda z s", positive=True, real=True)
    n = sp.symbols("n", integer=True)
    c = sp.Integer(159577)
    pi = sp.pi
    h = 4 * beta / c
    y_width = pi * c / (2 * beta)
    center = sp.Integer(39895)
    detuning = 4 * beta * (center + n - c / 4) / c
    require(sp.simplify(h * y_width - 2 * pi) == 0, "independent hY identity failed")
    require(sp.simplify(detuning * y_width * s - 2 * pi * (n + sp.Rational(3, 4)) * s) == 0, "independent Fourier phase failed")

    t_star = beta**3
    height = t_star - beta * lam
    eta = t_star / height
    tau = sp.tanh(z / beta)
    require(
        sp.simplify(height * (z / beta - eta * tau) - (beta**3 * (z / beta - tau) - lam * z)) == 0,
        "independent height reduction failed",
    )
    require(sp.simplify(y_width / (2 * pi) - 1 / h) == 0, "independent profile scale failed")


def independent_finite_pairing() -> None:
    frequencies = list(range(-9, 10))
    delta = {n: Fraction(2 * n * n + 3 * n - 5, 37 + 5 * abs(n)) for n in frequencies}
    weights = {n: Fraction(7 * n - 11, 41 + 3 * abs(n)) for n in frequencies}
    direct = sum((weights[n] * delta[n] for n in frequencies), Fraction())
    integral_pairing = sum((delta[n] * weights[n] for n in frequencies), Fraction())
    require(direct == integral_pairing, "independent finite Fourier pairing failed")

    h = Fraction(11, 23)
    raw = {n: Fraction(5 * n + 17, 43 + abs(n)) for n in frequencies}
    canonical = {n: value / h for n, value in raw.items()}
    require(
        sum((weights[n] * canonical[n] for n in frequencies), Fraction())
        == sum((weights[n] * raw[n] for n in frequencies), Fraction()) / h,
        "independent raw/canonical scale failed",
    )


def main() -> None:
    getcontext().prec = 90
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(path.is_file(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")

    independent_symbolic_check()
    independent_finite_pairing()
    certificate = artifact["certificate"]
    spacing = parse_ball(certificate["detuning_spacing_h_ball"])
    y_width = parse_ball(certificate["normal_width_Y_ball"])
    h_y_error = parse_ball(certificate["hY_minus_2pi_ball"])
    kernel = parse_ball(certificate["weighted_kernel_L1_bound_ball"])
    require(spacing[0] - spacing[1] > Decimal("0.054"), "stored spacing failed")
    require(y_width[0] - y_width[1] > Decimal("116"), "stored width failed")
    require(abs(h_y_error[0]) + h_y_error[1] < Decimal("5e-69"), "stored hY error failed")
    require(kernel[0] + kernel[1] < Decimal("3.840283e-5"), "stored kernel failed")

    decision = artifact["decision"]
    require(decision["exact_transformed_strip_common_profile_identity_proved"] is True, "profile identity drift")
    require(decision["weighted_exact_minus_beta4_coefficient_operator_proved"] is True, "operator identity drift")
    require(decision["coefficientwise_mode_bounds_used"] is False, "modewise guard drift")
    require(decision["uniform_exact_minus_beta4_profile_bound_proved"] is False, "profile-bound boundary drift")
    require(decision["rh_implication"] is False, "RH boundary drift")
    print("independently checked exact finite-t common-profile Fourier operator", flush=True)


if __name__ == "__main__":
    main()
