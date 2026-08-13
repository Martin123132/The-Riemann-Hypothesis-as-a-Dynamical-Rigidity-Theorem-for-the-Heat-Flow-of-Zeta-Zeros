#!/usr/bin/env python3
"""Independent replay of the B-face truncation dictionary correction."""

from __future__ import annotations

import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_fresnel_boundary_truncation_dictionary_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball_upper(text: str) -> Decimal:
    match = re.fullmatch(r"\[([^\]]+) \+/- ([^\]]+)\]", text)
    require(match is not None, f"unparsed Arb ball: {text}")
    midpoint, radius = (Decimal(value) for value in match.groups())
    return midpoint + radius


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Re-derive the dictionary at one exact rational witness without importing
    # the production symbolic expression.
    x, endpoint, mode = sp.Rational(7, 29), sp.Integer(43), sp.Integer(5)
    pi, imaginary_unit = sp.pi, sp.I
    c = endpoint * x / 2
    positive = endpoint / (imaginary_unit * pi * (endpoint * x - 2 * mode)) - mode / (2 * pi**2 * (c - mode) ** 3) + 3 * imaginary_unit * mode * x / (4 * pi**3 * (c - mode) ** 5)
    negative = endpoint / (imaginary_unit * pi * (endpoint * x + 2 * mode)) + mode / (2 * pi**2 * (c + mode) ** 3) - 3 * imaginary_unit * mode * x / (4 * pi**3 * (c + mode) ** 5)
    d = c**2 - mode**2
    currents = -2 * (2 + imaginary_unit * pi * endpoint**2 * x) / (pi**2 * (endpoint**2 * x**2 - 4 * mode**2))
    currents += c / (pi**3 * endpoint) * ((-4 * pi * endpoint * c**3 + 4 * imaginary_unit * c**2) / d**3 + (5 * pi * endpoint * c - 3 * imaginary_unit) / d**2)
    currents += -imaginary_unit * c**2 / (pi**4 * endpoint**2) * ((-48 * pi * endpoint * c**5 + 48 * imaginary_unit * c**4) / d**5 + (84 * pi * endpoint * c**3 - 60 * imaginary_unit * c**2) / d**4 + (-35 * pi * endpoint * c + 15 * imaginary_unit) / d**3)
    correction = 48 * x**2 * (endpoint**4 * x**4 + 40 * endpoint**2 * mode**2 * x**2 + 80 * mode**4) / (pi**4 * (endpoint * x - 2 * mode) ** 5 * (endpoint * x + 2 * mode) ** 5)
    require(sp.simplify(currents - positive - negative - correction) == 0, "independent rational witness failed")

    mp.mp.dps = 80
    t, endpoint_f = mp.mpf(10_000_000_000), mp.mpf(5_122_421)
    x0 = (1 - mp.sqrt(1 - 8 * t / (mp.pi * endpoint_f**2))) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    x_low, x_high = x0 - 70 / mp.sqrt(hessian), x0 + 70 / mp.sqrt(hessian)
    total = mp.mpf(0)
    for m in range(1, 200_001):
        if m in (621, 622):
            continue
        numerator = endpoint_f**4 * x_high**4 + 40 * endpoint_f**2 * m**2 * x_high**2 + 80 * m**4
        if m <= 620:
            left, right = endpoint_f * x_low - 2 * m, endpoint_f * x_low + 2 * m
        else:
            left, right = 2 * m - endpoint_f * x_high, 2 * m + endpoint_f * x_low
        total += 48 * x_high**2 * numerator / (mp.pi**4 * left**5 * right**5)
    physical = 2 * (mp.pi / (32 * t)) ** mp.mpf("0.25") * (x_high - x_low) * (x_low * (1 - x_low)) ** mp.mpf("-0.25") * total
    require(physical < mp.mpf("2.1e-20"), "independent 200000-mode correction bound failed")

    certificate = artifact["interval_certificate"]
    require(ball_upper(certificate["physical_dictionary_correction_ball"]) < Decimal("2.1e-20"), "stored correction too large")
    decision = artifact["decision"]
    require(decision["direct_Fresnel_and_boundary_current_truncations_identical"] is False, "dictionary distinction lost")
    require(decision["complete_B_face_estimate_proved"] is False, "complete-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked Fresnel/boundary truncation dictionary correction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
