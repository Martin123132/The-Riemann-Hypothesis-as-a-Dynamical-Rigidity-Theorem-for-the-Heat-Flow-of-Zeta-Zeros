#!/usr/bin/env python3
"""Independently check symmetric finite-Poisson Kummer interchange."""

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


RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.md"
A = 159577
B = 5122421
L = (B - A) // 2


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
    require(decision["symmetric_poisson_interchange_proved"] is True, "interchange decision drift")
    require(decision["auxiliary_endpoint_regulator_required"] is False, "regulator decision drift")
    require(decision["paired_one_over_m_endpoint_current_cancels"] is True, "endpoint cancellation drift")
    require(decision["tail_bound_uniform_in_height"] is True, "height-uniform decision drift")
    require(decision["raw_c2_tail_bound_quantitatively_viable"] is False, "scale guard drift")
    require(decision["stationary_phase_remainder_proved"] is False, "stationary theorem overpromotion")

    u, x, m = sp.symbols("u x m", real=True)
    alpha_u = sp.Integer(A) + 2 * u
    phase = sp.exp(sp.I * sp.pi * x * alpha_u**2 / 4)
    f = alpha_u * phase
    f1 = sp.simplify(sp.diff(f, u) / phase)
    f2 = sp.simplify(sp.diff(f, u, 2) / phase)
    require(sp.simplify(f1 - (2 + sp.I * sp.pi * x * alpha_u**2)) == 0, "independent first derivative failed")
    require(sp.simplify(f2 - (6 * sp.I * sp.pi * x * alpha_u - sp.pi**2 * x**2 * alpha_u**3)) == 0, "independent second derivative failed")

    D = 2 * sp.pi * sp.I * m
    df, df1, rp, rn = sp.symbols("df df1 rp rn")
    ip = -df / D - df1 / D**2 + rp / D**2
    im = ip.subs({m: -m, rp: rn})
    pair = sp.simplify(ip + im)
    require(sp.simplify(sp.diff(pair, df)) == 0, "independent 1/m cancellation failed")
    require(sp.simplify(pair - (-2 * df1 + rp + rn) / D**2) == 0, "independent paired identity failed")

    beta0 = sp.beta(sp.Rational(3, 4), sp.Rational(3, 4))
    require(sp.simplify(sp.expand_func(sp.beta(sp.Rational(7, 4), sp.Rational(3, 4)) - beta0 / 2)) == 0, "independent beta-1 recurrence failed")
    require(sp.simplify(sp.expand_func(sp.beta(sp.Rational(11, 4), sp.Rational(3, 4)) - 7 * beta0 / 20)) == 0, "independent beta-2 recurrence failed")

    ctx.dps = 130
    q = arb(3) / 4
    beta_ball = q.gamma() ** 2 / (arb(3) / 2).gamma()
    a = arb(A)
    b = arb(B)
    pi = arb.pi()
    c_ab = beta_ball * (4 + pi * (5 * b**2 - a**2) / 4 + 7 * pi**2 * (b**4 - a**4) / 160)
    coefficient = c_ab / (2 * pi**2)
    cutoff = coefficient / arb("0.005")
    saved = artifact["certified_constants"]
    require(saved["alpha_count"] == L + 1 == 2481423, "roster count drift")
    require(arb(saved["beta_3_4_3_4_ball"]).overlaps(beta_ball), "beta ball drift")
    require(arb(saved["integrated_C_AB_ball"]).overlaps(c_ab), "C_AB ball drift")
    require(arb(saved["tail_coefficient_C_AB_over_2pi2_ball"]).overlaps(coefficient), "tail coefficient drift")
    require(arb(saved["raw_cutoff_required_for_bound_below_0_005_ball"]).overlaps(cutoff), "cutoff ball drift")
    require(coefficient > arb("1e25") and cutoff > arb("1e27"), "absolute-tail scale guard failed")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {source['path']}")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "The apparent `Delta f_x/(2*pi*i*m)` endpoint current cancels exactly",
        "every real `t` and every integer",
        "no auxiliary endpoint regulator is required",
        "not a viable quantitative Hardy-error estimate",
        "stationary-phase remainder",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "validated symmetric finite-Poisson interchange independently: "
        "alpha=2481423, uniform_tail=True, quantitative=False"
    )


if __name__ == "__main__":
    main()
