#!/usr/bin/env python3
"""Independently check the selected detuning Airy Green-kernel gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

from flint import arb, ctx
import mpmath as mp
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = (C - 1) // 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    beta, lam, d = sp.symbols("beta lambda d", positive=True, real=True)
    i = sp.I
    chi = 2 * beta**2 * d + beta * d**2
    r = beta**2 + lam + 2 * beta * d
    hfun = sp.Function("h")
    g = sp.exp(-i * chi) * hfun(r)
    lhs = sp.diff(g, d, 2) / (4 * beta**2)
    lhs += i * (1 + d / beta) * sp.diff(g, d)
    lhs += (lam - d**2 + i / (2 * beta)) * g
    rhs = sp.exp(-i * chi) * (sp.diff(hfun(r), d, 2) / (4 * beta**2) + r * hfun(r))
    require(sp.simplify(lhs - rhs) == 0, "independent gauge conjugation failed")

    ctx.dps = 100
    pi = arb.pi()
    beta_ball = (pi * arb(C) ** 2 / 8) ** (arb(1) / 3)
    lambda_max = pi / (16 * beta_ball)
    r_min = beta_ball**2 * (arb(8) * 39_696 / C - 1)
    r_max = beta_ball**2 * (arb(8) * 40_094 / C - 1) + lambda_max
    certificate = artifact["interval_certificate"]
    require(r_min.overlaps(arb(certificate["r_min_ball"])), "r minimum drift")
    require(r_max.overlaps(arb(certificate["r_max_ball"])), "r maximum drift")
    require((8 * beta_ball**2 / C).overlaps(arb(certificate["exact_lattice_step_ball"])), "lattice step drift")
    require((1 / r_min.sqrt()).overlaps(arb(certificate["green_kernel_absolute_envelope_ball"])), "kernel envelope drift")

    mp.mp.dps = 70
    sample_s = mp.mpf("4620000.25")
    sample_r = mp.mpf("4650000.75")

    def kernel(rv: mp.mpf, sv: mp.mpf) -> mp.mpf:
        return mp.pi * (mp.airyai(-rv) * mp.airybi(-sv) - mp.airybi(-rv) * mp.airyai(-sv))

    diagonal = kernel(sample_s, sample_s)
    diagonal_derivative = mp.diff(lambda rv: kernel(rv, sample_s), sample_s)
    equation_residual = mp.diff(lambda rv: kernel(rv, sample_s), sample_r, 2) + sample_r * kernel(sample_r, sample_s)
    require(abs(diagonal) < mp.mpf("1e-60"), "Green diagonal is nonzero")
    require(abs(diagonal_derivative - 1) < mp.mpf("1e-55"), "Green diagonal derivative drift")
    require(abs(equation_residual) < mp.mpf("1e-50"), "Green Airy equation residual drift")

    step = 8 * beta_ball**2 / C
    center = beta_ball**2 * (1 - arb(2) / C) + arb(lambda_max / 2, lambda_max / 2)
    for j in (1, 17, 198):
        r_minus = beta_ball**2 * (arb(8) * (Q - j) / C - 1) + arb(lambda_max / 2, lambda_max / 2)
        r_plus = beta_ball**2 * (arb(8) * (Q + j) / C - 1) + arb(lambda_max / 2, lambda_max / 2)
        require((r_minus - (center - j * step)).contains(0), f"minus lattice symmetry failed at {j}")
        require((r_plus - (center + j * step)).contains(0), f"plus lattice symmetry failed at {j}")
    print("independently checked selected detuning Airy Green kernel", flush=True)


if __name__ == "__main__":
    main()
