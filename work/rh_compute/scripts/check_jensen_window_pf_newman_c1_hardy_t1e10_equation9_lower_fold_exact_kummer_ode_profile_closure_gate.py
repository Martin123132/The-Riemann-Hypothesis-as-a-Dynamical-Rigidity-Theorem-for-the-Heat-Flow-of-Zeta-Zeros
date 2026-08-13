#!/usr/bin/env python3
"""Independently check the exact Kummer-ODE profile closure."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_kummer_ode_profile_closure_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
PANELS = 16_384
PRECISION = 65


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_algebra_check() -> None:
    beta, lam, y = sp.symbols("beta lambda y", positive=True, real=True)
    rho = 1 + y / (2 * beta**2)
    t = beta**3 - beta * lam
    potential = sp.simplify((beta**3 * rho**2 - t) / beta)
    require(sp.expand(potential - (lam + y + y**2 / (4 * beta**2))) == 0, "independent ODE potential check failed")

    # Rebuild the residual directly from a formal Airy jet A''=-X A.
    e = sp.symbols("e", positive=True, real=True)
    X = lam + y
    U = 1 + e * (-(13 * lam + 3 * y) / 60) + e**2 * (
        -(448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2) / 50400
    )
    V = e * ((8 * lam**2 - 4 * lam * y + 3 * y**2) / 60) + e**2 * (
        (40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27) / 1680
    )
    d1u, d1v = sp.diff(U, y) - X * V, U + sp.diff(V, y)
    d2u, d2v = sp.diff(d1u, y) - X * d1v, d1u + sp.diff(d1v, y)
    ru = sp.expand(d2u + (X + e * y**2 / 4) * U)
    rv = sp.expand(d2v + (X + e * y**2 / 4) * V)
    require(all(ru.coeff(e, k) == 0 and rv.coeff(e, k) == 0 for k in range(3)), "lower residual order survived")
    require(sp.degree(ru, e) == 3 and sp.degree(rv, e) == 3, "residual is not exactly cubic in beta^-2")


def numerical_recheck(artifact: dict) -> None:
    ctx.dps = PRECISION
    pi = arb.pi()
    c = arb(159577)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    epsilon = beta**-2
    y_width = pi * c / (2 * beta)
    h = 4 * beta / c
    lambda_max = pi / (16 * beta)
    lam = arb(lambda_max / 2, lambda_max / 2)
    width = y_width / PANELS
    residual_l1 = arb(0)
    modulus_upper = arb(0)
    for panel in range(PANELS):
        y = arb(y_width * (2 * panel + 1) / (2 * PANELS), width / 2)
        x = lam + y
        airy = (-x).airy_ai()
        airy_x = -(-x).airy_ai(derivative=1)
        bi = (-x).airy_bi()
        p3 = -y**2 * (448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2 - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2) / 201600
        q3 = -y**2 * (-40 * lam**3 + 20 * lam**2 * y - lam * y**2 + 9 * y**3 - 27) / 6720
        residual_l1 += width * (
            2 * pi * epsilon**3 * (abs(p3).upper() * abs(airy).upper() + abs(q3).upper() * abs(airy_x).upper())
        ).upper()
        modulus = (airy * airy + bi * bi).sqrt().upper()
        if modulus > modulus_upper:
            modulus_upper = modulus

    saved = artifact["interval_certificate"]
    require(residual_l1 < arb(saved["residual_L1_upper_bound"]), "16384-panel residual bound is not below saved 8192-panel bound")
    require(modulus_upper < arb(saved["Airy_modulus_upper_bound"]), "16384-panel Airy modulus is not below saved bound")

    derivative_modulus = ((-lam).airy_ai(derivative=1)**2 + (-lam).airy_bi(derivative=1)**2).sqrt().upper()
    kernel = (pi * modulus_upper**2).upper()
    dkernel = (pi * modulus_upper * derivative_modulus).upper()
    qint = (epsilon * y_width**3 / 12).upper()
    denominator = (1 - kernel * qint).lower()
    e0 = arb(saved["initial_value_uniform_bound"])
    e1 = arb(saved["initial_slope_uniform_bound"])
    raw = ((dkernel * e0 + kernel * e1 + kernel * residual_l1) / denominator).upper()
    kernel_l1 = arb(saved["weighted_kernel_L1_upper_bound"])
    weighted = (raw / h * kernel_l1).upper()
    require(weighted < arb("1.21e-11"), "independent weighted correction exceeds 1.21e-11")


def main() -> None:
    require(RESULT.is_file(), "missing result artifact")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "saved gate did not pass")
    require(artifact["decision"]["full_exact_minus_beta4_profile_bound_proved"] is True, "saved profile decision failed")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(path.is_file() and file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    exact_algebra_check()
    numerical_recheck(artifact)
    print("PASS: independent exact algebra and 16384-panel Kummer-ODE profile closure check")


if __name__ == "__main__":
    main()
