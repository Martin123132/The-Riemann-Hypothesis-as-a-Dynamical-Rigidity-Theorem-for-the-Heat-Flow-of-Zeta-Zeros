#!/usr/bin/env python3
"""Independently check the beta^-6 common-profile coefficient."""

from __future__ import annotations

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
PANELS = 8192


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_symbolic_check(stored: dict[str, str]) -> None:
    q, z, y, x, lam = sp.symbols("q z y X lambda", real=True)
    ratio = sp.cosh(q * z) ** (-sp.Rational(3, 2)) * (1 + q**2 * y / 2)
    remainder = (
        (q * z - sp.tanh(q * z)) / q**3 - z**3 / 3
        + y * (z - sp.tanh(q * z) / q) - q * y**2 * sp.tanh(q * z) / 4
    )
    c3 = sp.expand(sp.series(ratio * sp.exp(sp.I * remainder), q, 0, 8).removeO().coeff(q, 6))
    locals_map = {"z": z, "y": y, "X": x, "lam": lam, "I": sp.I}
    require(sp.expand(c3 - sp.sympify(stored["C3"], locals=locals_map)) == 0, "independent C3 drift")

    p, r = sp.Integer(1), sp.Integer(0)
    u3, v3 = sp.Integer(0), sp.Integer(0)
    for degree in range(16):
        coefficient = c3.coeff(z, degree)
        u3 += coefficient * sp.I**degree * p
        v3 += coefficient * sp.I**degree * r
        p, r = sp.expand(sp.diff(p, x) - x * r), sp.expand(p + sp.diff(r, x))
    stored_u3 = sp.sympify(stored["U3"].replace("lambda", "lam"), locals=locals_map)
    stored_v3 = sp.sympify(stored["V3"].replace("lambda", "lam"), locals=locals_map)
    require(sp.expand(u3.subs(x, lam + y) - stored_u3) == 0, "independent U3 drift")
    require(sp.expand(v3.subs(x, lam + y) - stored_v3) == 0, "independent V3 drift")


def uv3(lam: arb, y: arb) -> tuple[arb, arb]:
    u3 = (
        -49856 * lam**6 + 1344 * lam**5 * y + 3360 * lam**4 * y**2
        - 12920 * lam**3 * y**3 - 516535 * lam**3 + 7785 * lam**2 * y**4
        + 5445 * lam**2 * y - 1836 * lam * y**5 + 21195 * lam * y**2
        + 1674 * y**6 - 21285 * y**3 - 265860
    ) / 9072000
    v3 = -(
        3584 * lam**7 - 1792 * lam**6 * y + 1344 * lam**5 * y**2
        + 2240 * lam**4 * y**3 - 43880 * lam**4 - 1960 * lam**3 * y**4
        + 18580 * lam**3 * y + 1764 * lam**2 * y**5 + 5445 * lam**2 * y**2
        - 567 * lam * y**6 + 14130 * lam * y**3 - 127530 * lam
        + 189 * y**7 - 7605 * y**4 + 42570 * y
    ) / 9072000
    return u3, v3


def independent_interval_bound() -> tuple[arb, arb, arb]:
    ctx.dps = 90
    pi = arb.pi()
    c = arb(C)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    epsilon = beta**-2
    y_width = pi * c / (2 * beta)
    lambda_max = pi / (16 * beta)
    lam = arb(lambda_max / 2, lambda_max / 2)
    maximum = arb(0)
    for panel in range(PANELS):
        y = arb(y_width * (2 * panel + 1) / (2 * PANELS), y_width / (2 * PANELS))
        u, v = uv3(lam, y)
        x = lam + y
        ai = (-x).airy_ai()
        aip = (-x).airy_ai(derivative=1)
        bound = 2 * pi * epsilon**3 * (abs(u).upper() * abs(ai).upper() + abs(v).upper() * abs(aip).upper())
        if bound.upper() > maximum.upper():
            maximum = bound
    h = 4 * beta / c
    return maximum, maximum / h, h


def main() -> None:
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

    independent_symbolic_check(artifact["polynomial_certificate"])
    raw, profile, _ = independent_interval_bound()
    require(raw < arb("6.951e-10"), "independent reduced coefficient bound failed")
    require(profile < arb("1.288e-8"), "independent profile coefficient bound failed")
    stored = artifact["interval_certificate"]
    require(arb(stored["uniform_reduced_beta6_coefficient_bound_ball"]) < arb("6.951e-10"), "stored raw bound failed")
    require(arb(stored["uniform_canonical_beta6_profile_bound_ball"]) < arb("1.288e-8"), "stored profile bound failed")
    require(arb(stored["uniform_weighted_beta6_coefficient_bound_ball"]) < arb("4.95e-13"), "stored weighted bound failed")
    decision = artifact["decision"]
    require(decision["beta6_Airy_two_channel_contraction_exact"] is True, "contraction drift")
    require(decision["exact_post_beta6_Taylor_remainder_bound_proved"] is False, "tail boundary drift")
    require(decision["rh_implication"] is False, "RH boundary drift")
    print("independently checked beta^-6 common-profile leading coefficient", flush=True)


if __name__ == "__main__":
    main()
