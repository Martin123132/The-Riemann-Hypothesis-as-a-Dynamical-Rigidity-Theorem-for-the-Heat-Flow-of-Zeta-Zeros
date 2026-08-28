#!/usr/bin/env python3
"""Independently check the affine-wedge remainder identity certificate."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_remainder_identity_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
HEIGHT = 10_000_000_000
A = 159_577
B = 5_122_421
KEYS = ((B, 621), (B, 622), (A, 39_852), (A, 39_853), (A, 39_894), (A, 39_895))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def check_symbolic() -> None:
    u, v, sigma = sp.symbols("u v sigma", positive=True, real=True)
    m, endpoint, t = sp.symbols("m D t", positive=True, real=True)
    p = sp.sqrt(u) - 1 / sp.sqrt(u)
    z = 2 * m * u / endpoint
    P = sp.sqrt(sp.pi * m * endpoint) * p
    r = t / (2 * sp.pi * m**2)
    x = 1 / (1 + r * v)
    dz_dP = sp.diff(z, u) / sp.diff(P, u)
    minus_dx_dS = -sp.diff(x, v) * v * sigma / (v - 1) * sp.sqrt(2 / t)
    jacobian = sp.simplify(
        z ** (-sp.Rational(1, 2))
        * (2 + sp.I * sp.pi * endpoint**2 * z)
        * x ** (-sp.Rational(7, 4))
        * (1 - x) ** (-sp.Rational(1, 4))
        * dz_dP
        * minus_dx_dS
    )
    K = sp.pi * endpoint * m
    W0 = 4 * (K - sp.I) * r ** sp.Rational(3, 4) / (endpoint * sp.sqrt(sp.pi * t))
    E = 2 * u * (K * u - sp.I) / ((u + 1) * (K - sp.I))
    O = v ** sp.Rational(3, 4) * sigma / (v - 1)
    require(sp.simplify(jacobian / sp.I - W0 * E * O) == 0, "fresh amplitude factorization failed")

    chi_e, chi_t, amp, aff = sp.symbols("chi_e chi_t amp aff")
    require(sp.expand(chi_e * (amp - aff) + (chi_e - chi_t) * aff - (chi_e * amp - chi_t * aff)) == 0, "fresh signed split failed")
    C, C_a, C_h, lam, mu, rho = sp.symbols("C C_a C_h lambda mu rho")
    require(sp.simplify(C + lam * (-sp.I * C_a) + mu * sp.I * (C_h - rho * C_a) - (C + sp.I * mu * C_h - sp.I * (lam + mu * rho) * C_a)) == 0, "fresh affine moment failed")

    v_series = 1 + sigma + sigma**2 / 3 + sigma**3 / 36 - sigma**4 / 270 + sigma**5 / 4320
    require(sp.series(2 * (v_series - 1 - sp.log(v_series)) - sigma**2, sigma, 0, 6).removeO().expand() == 0, "fresh Morse series failed")
    require(
        sp.simplify(sp.series(v_series ** sp.Rational(3, 4) * sigma / (v_series - 1), sigma, 0, 3).removeO() - (1 + 5 * sigma / 12 - sigma**2 / 96)) == 0,
        "fresh outer slope failed",
    )


def fresh_row(endpoint_int: int, mode_int: int) -> dict[str, arb]:
    pi, t = arb.pi(), arb(HEIGHT)
    endpoint, mode = arb(endpoint_int), arb(mode_int)
    denominator = 2 * pi * mode**2 + t
    alpha = 2 * mode + t / (pi * mode)
    a = (pi * mode / alpha).sqrt() * (endpoint - alpha)
    rho = -(2 * t).sqrt() * (pi * endpoint * mode + denominator) / (2 * denominator ** arb("1.5"))
    curvature = -(
        8 * pi**2 * endpoint * mode**3
        - 5 * pi * endpoint * mode * t
        + 16 * pi**2 * mode**4
        + 10 * pi * mode**2 * t
        + t**2
    ) / (6 * denominator ** arb("2.5"))
    K = pi * endpoint * mode
    scale = (pi * mode * endpoint).sqrt()
    lambda_re = (3 * K**2 + 1) / (2 * (K**2 + 1) * scale)
    lambda_im = K / ((K**2 + 1) * scale)
    r = t / (2 * pi * mode**2)
    w0_scale = 4 * r ** arb("0.75") / (endpoint * (pi * t).sqrt())
    return {
        "a_face_detuning_ball": a,
        "rho_signed_ball": rho,
        "c_characteristic_defect_ball": rho**2 - 1,
        "face_second_derivative_ball": curvature,
        "lambda_P_real_ball": lambda_re,
        "lambda_P_imag_ball": lambda_im,
        "mu_S_ball": arb(5) / 12 * (arb(2) / t).sqrt(),
        "W0_real_ball": w0_scale * K,
        "W0_imag_ball": -w0_scale,
    }


def main() -> int:
    ctx.dps = 120
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True, "artifact is not passed")
    require(artifact.get("kind") == STEM, "kind drift")
    decisions = artifact["decision"]
    for key in (
        "exact_transformed_amplitude_factorized",
        "mandatory_affine_P_current_retained",
        "outer_affine_S_current_retained",
        "affine_wedge_reduced_to_C_and_boundary_derivatives",
        "exact_defect_split_into_signed_amplitude_and_face_remainders",
        "phase_remainder_identically_zero",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in ("global_amplitude_remainder_bound_proved", "global_curved_face_remainder_bound_proved", "R_Dir_bound_proved", "rh_implication"):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    check_symbolic()
    rows = {(row["endpoint"], row["mode"]): row for row in artifact["geometry_certificate"]["rows"]}
    require(set(rows) == set(KEYS), "geometry key drift")
    for key in KEYS:
        fresh = fresh_row(*key)
        saved = rows[key]
        for field, value in fresh.items():
            require(value.overlaps(arb(saved[field])), f"saved interval misses fresh {field} at {key}")
        require(abs(fresh["face_second_derivative_ball"]).upper() < arb("7e-6"), f"curvature cap failed at {key}")
        require(fresh["lambda_P_real_ball"].upper() < arb("2e-5"), f"P-slope cap failed at {key}")
        require(abs(fresh["lambda_P_imag_ball"]).upper() < arb("1e-14"), f"imaginary P-slope cap failed at {key}")
        require(fresh["mu_S_ball"].upper() < arb("6e-6"), f"S-slope cap failed at {key}")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("R_exact-aff=R_amp+R_face" in note, "note identity missing")
    require("No global `R_amp` or `R_face` bound" in note, "note proof boundary missing")
    print("independently checked exact affine-wedge remainder identity", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
