#!/usr/bin/env python3
"""Independently check the Airy-lattice cubic resonance gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_airy_lattice_cubic_resonance_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = (C - 1) // 4
COEFFICIENTS = (sp.Rational(3, 2), sp.Rational(3, 8), sp.Rational(-1, 16), sp.Rational(3, 128), sp.Rational(-3, 256))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def normal_form(j: int, u: mp.mpf) -> mp.mpf:
    xj = (8 * j - 2) / mp.mpf(C) + u
    x0 = -2 / mp.mpf(C) + u
    total = mp.mpf("0")
    for degree, coefficient in enumerate(COEFFICIENTS, start=1):
        total += mp.mpf(int(coefficient.p)) / int(coefficient.q) * (xj**degree - x0**degree)
    return mp.pi * C**2 / 12 * total - mp.pi * j * (C - 1 + 2 * j)


def exact_leading(j: int, u: mp.mpf) -> mp.mpf:
    xj = (8 * j - 2) / mp.mpf(C) + u
    x0 = -2 / mp.mpf(C) + u
    return mp.pi * C**2 / 12 * ((1 + xj) ** mp.mpf("1.5") - (1 + x0) ** mp.mpf("1.5"))


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    c, j = sp.symbols("C j", positive=True, real=True)
    center = -sp.pi * j * (16 * j**2 - 12 * j + 3) / (6 * c)
    center += sp.pi * ((8 * j - 2) ** 4 - 16) / (512 * c**2)
    center -= sp.pi * ((8 * j - 2) ** 5 + 32) / (1024 * c**3)
    xj = (8 * j - 2) / c
    x0 = -2 / c
    reconstructed = sp.pi * c**2 / 12 * sum(
        COEFFICIENTS[k - 1] * (xj**k - x0**k) for k in range(1, 6)
    ) - sp.pi * j * (c - 1 + 2 * j)
    require(sp.expand(center - reconstructed) == 0, "independent center polynomial failed")

    mp.mp.dps = 90
    u_max = mp.mpf(1) / (2 * C**2)
    bound_text = artifact["interval_certificate"]["uniform_binomial_degree5_remainder_ball"]
    bound = mp.mpf(bound_text.split("+/-")[0].strip("[ ")) + mp.mpf(bound_text.split("+/-")[1].split("]")[0].strip())
    for index in (-198, -150, -100, -17, -2, -1, 1, 2, 17, 100, 150, 200):
        for height in (mp.mpf("0"), u_max / 2, u_max):
            exact = exact_leading(index, height) - mp.pi * index * (C - 1 + 2 * index)
            error = abs(exact - normal_form(index, height))
            require(error < bound, f"quintic remainder failed at j={index}")
    exact_one_low = exact_leading(1, mp.mpf("0")) - 2 * mp.pi * (2 * Q + 1)
    exact_one_high = exact_leading(1, u_max) - 2 * mp.pi * (2 * Q + 1)
    require(exact_one_low < 0 and exact_one_high < 0, "central adjacent resonance changed sign")
    require(abs(exact_one_low) < mp.mpf("2.4e-5") and abs(exact_one_high) > mp.mpf("1.7e-5"), "central adjacent resonance escaped interval")
    require((C - 1) % 4 == 0, "selector parity drift")
    print("independently checked Airy-lattice cubic resonance normal form", flush=True)


if __name__ == "__main__":
    main()
