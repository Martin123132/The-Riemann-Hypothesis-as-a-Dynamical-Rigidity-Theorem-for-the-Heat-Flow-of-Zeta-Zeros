#!/usr/bin/env python3
"""Independent replay of the completed-B outer Fresnel remainder derivatives."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_fresnel_remainder_derivative_tail_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tail(exponent: int) -> arb:
    start = arb(B)
    return start ** (-exponent) + start ** (1 - exponent) / (exponent - 1)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    q = sp.symbols("q", nonzero=True, real=True)
    z = sp.Function("z")(q)
    pi = sp.pi
    p3 = sp.I / (pi * q) + 1 / (pi**2 * q**3) - 3 * sp.I / (pi**3 * q**5)
    z_prime = -1 - sp.I * pi * q * z
    r_prime = sp.simplify(z_prime - sp.diff(p3, q)).subs(z, sp.Symbol("r") + p3)
    require(sp.simplify(r_prime + sp.I * pi * q * sp.Symbol("r") + 15 * sp.I / (pi**3 * q**6)) == 0, "independent remainder ODE failed")

    ctx.dps = 110
    ctx.threads = 1
    pi_ball = arb.pi()
    alpha, beta, gamma = arb(3) / 4, arb(5) / 4, arb(13) / 4
    delta, xmax = arb("1e-4"), arb(1) / 2
    a0 = arb(15) / (4 * pi_ball**4 * alpha**7)
    a11 = arb(45) * beta / (8 * pi_ball**3 * alpha**6)
    one0 = a0 * xmax**2 * tail(6)
    one1 = (3 * a0 / 2) * xmax * tail(6) + a11 * tail(4)
    one2 = (
        (15 * a0 / 4) * tail(6)
        + 3 * a11 * delta**-1 * tail(4)
        + 45 * beta**2 / (8 * pi_ball**2 * alpha**5) * delta**-2 * tail(2)
        + 15 * beta**2 / (2 * pi_ball**3 * alpha**7) * delta**-1 * tail(4)
        + 45 * gamma / (16 * pi_ball**3 * alpha**6) * delta**-1 * tail(4)
    )
    raw = {0: 2 * one0, 1: 2 * one1, 2: 2 * one2}

    endpoint, t = arb(B), arb(T)
    x0 = (1 - (1 - 8 * t / (pi_ball * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    weight = (delta * (1 - delta)) ** (-arb(1) / 4)
    l1 = (1 - 2 * delta) / (4 * delta * (1 - delta))
    w2 = (12 * delta**2 - 12 * delta + 5) / (16 * delta**2 * (1 - delta) ** 2)
    normalized = {
        0: weight * raw[0],
        1: weight * (raw[1] + l1 * raw[0]) / hessian.sqrt(),
        2: weight * (raw[2] + 2 * l1 * raw[1] + w2 * raw[0]) / hessian,
    }

    stored_raw = artifact["interval_certificate"]["two_sign_remainder_x_derivative_tail_balls"]
    stored_norm = artifact["interval_certificate"]["weighted_normalized_derivative_tail_balls"]
    for order in range(3):
        require(arb(stored_raw[str(order)]).overlaps(raw[order]), f"raw order-{order} replay mismatch")
        require(arb(stored_norm[str(order)]).overlaps(normalized[order]), f"normalized order-{order} replay mismatch")

    decision = artifact["decision"]
    require(decision["Fresnel_to_boundary_dictionary_derivative_tail_bounded"] is False, "dictionary overclaim")
    require(decision["complete_outer_mode_derivative_block_bounded"] is False, "outer-block overclaim")
    require(decision["finite_mode_block_below_outer_start_bounded"] is False, "finite-block overclaim")
    require(decision["completed_exterior_current_bound_proved"] is False, "completed-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(
        "independently checked outer Fresnel remainder derivatives: "
        + ", ".join(f"order{order}={normalized[order]}" for order in range(3)),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
