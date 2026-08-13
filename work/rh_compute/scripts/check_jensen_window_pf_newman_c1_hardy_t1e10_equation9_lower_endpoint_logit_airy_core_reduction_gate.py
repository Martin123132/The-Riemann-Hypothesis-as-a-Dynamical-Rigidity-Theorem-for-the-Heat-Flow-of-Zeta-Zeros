#!/usr/bin/env python3
"""Independently check the exact lower-endpoint logit Airy core."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.md"
T = 10_000_000_000
A = 159577
CORE_Z = 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    artifact = load(RESULT)
    decision = artifact["decision"]
    require(artifact["passed"] is True, "artifact is not passed")
    require(decision["exact_logit_fold_phase_proved"] is True, "fold phase decision drift")
    require(decision["canonical_Airy_scaling_proved"] is True, "Airy scaling decision drift")
    require(decision["Airy_core_phase_remainder_below_0_001"] is True, "core phase decision drift")
    require(decision["portcullis_roots_recovered"] is True, "root decision drift")
    require(decision["separated_endpoint_mode_series_summable"] is False, "grouping guard drift")
    require(decision["grouped_amplitude_remainder_proved"] is False, "remainder overpromotion")

    u, t, C, m, z = sp.symbols("u t C m z", real=True, positive=True)
    pi = sp.pi
    x = 1 / (1 + sp.exp(u))
    eta = pi * C**2 / (8 * t)
    mu = eta - 1
    phase = pi * C**2 * x / 4 - pi * m * C + t * u / 2
    constant = pi * C**2 / 8 - pi * m * C
    fold = t * (u / 2 - eta * sp.tanh(u / 2))
    require(sp.simplify(sp.expand_trig(phase - constant - fold).rewrite(sp.exp)) == 0, "independent fold phase failed")
    derivative = (1 - eta * sp.sech(u / 2) ** 2) / 2
    root = 2 * sp.acosh(sp.sqrt(eta))
    require(sp.simplify(derivative.subs(u, root)) == 0, "independent stationary root failed")
    scale = 2 / (t * eta) ** sp.Rational(1, 3)
    airy_lambda = mu * t ** sp.Rational(2, 3) * eta ** (-sp.Rational(1, 3))
    cubic = -t * mu * u / 2 + t * eta * u**3 / 24
    require(sp.simplify(sp.powsimp(cubic.subs(u, scale * z) - z**3 / 3 + airy_lambda * z, force=True)) == 0, "independent Airy scaling failed")
    measure = sp.simplify((-sp.diff(x, u)) * (x * (1 - x)) ** (-sp.Rational(1, 4)))
    expected_measure = 2 ** (-sp.Rational(3, 2)) * sp.cosh(u / 2) ** (-sp.Rational(3, 2))
    require(sp.simplify(sp.expand_trig(measure - expected_measure).rewrite(sp.exp)) == 0, "independent measure failed")

    ctx.dps = 130
    tt = arb(T)
    cc = arb(A)
    api = arb.pi()
    eeta = api * cc**2 / (8 * tt)
    mmu = eeta - 1
    cscale = 2 / (tt * eeta) ** (arb(1) / 3)
    lam = mmu * tt ** (arb(2) / 3) * eeta ** (-arb(1) / 3)
    ustar = 2 * eeta.sqrt().acosh()
    xlower = 1 / (1 + ustar.exp())
    xupper = 1 / (1 + (-ustar).exp())
    nlower = cc * xlower / 2
    nupper = cc * xupper / 2
    bound = arb(64) / 15 * tt ** (-arb(2) / 3) * eeta ** (-arb(2) / 3) * arb(CORE_Z) ** 5
    saved = artifact["certified_core"]
    for field, value in (
        ("eta_ball", eeta),
        ("mu_ball", mmu),
        ("Airy_scale_ball", cscale),
        ("Airy_lambda_ball", lam),
        ("positive_stationary_u_ball", ustar),
        ("lower_stationary_x_ball", xlower),
        ("upper_stationary_x_ball", xupper),
        ("lower_continuous_mode_ball", nlower),
        ("upper_continuous_mode_ball", nupper),
        ("conservative_phase_remainder_bound_ball", bound),
    ):
        require(arb(saved[field]).overlaps(value), f"saved ball drift: {field}")
    require(arb(5) < lam < arb(6), "independent lambda classification failed")
    require(arb(39852) < nlower < arb(39853), "independent lower root failed")
    require(arb(39936) < nupper < arb(39937), "independent upper root failed")
    require(bound < arb("0.001"), "independent phase bound failed")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {source['path']}")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "phase is exactly",
        "canonical Airy polynomial exactly",
        "reproduce the prior portcullis roots",
        "cannot be summed separately",
        "canonical fold phase is no longer the unknown",
    ):
        require(token in note, f"note token missing: {token}")

    print("validated lower-endpoint logit Airy core independently: lambda=5.1099, roots=2, phase<0.001")


if __name__ == "__main__":
    main()
