#!/usr/bin/env python3
"""Independent replay of the completed-B outer dictionary derivative tail."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_dictionary_derivative_tail_gate"
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


def rows(order: int, p: int, k: int) -> tuple[tuple[int, int, int], ...]:
    candidates = {
        0: ((1, p, k),),
        1: ((p, p - 1, k), (2 * k, p + 1, k + 1)),
        2: ((p * (p - 1), p - 2, k), (2 * k * (2 * p + 1), p, k + 1), (4 * k * (k + 1), p + 2, k + 2)),
    }[order]
    return tuple(row for row in candidates if row[0])


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

    x, endpoint, mode, c = sp.symbols("x B m c", positive=True, real=True)
    original = 48 * x**2 * (endpoint**4 * x**4 + 40 * endpoint**2 * mode**2 * x**2 + 80 * mode**4) / (
        sp.pi**4 * (endpoint * x - 2 * mode) ** 5 * (endpoint * x + 2 * mode) ** 5
    )
    reduced = 3 * c**2 * (c**4 + 10 * mode**2 * c**2 + 5 * mode**4) / (
        sp.pi**4 * endpoint**2 * (c**2 - mode**2) ** 5
    )
    require(sp.simplify(original.subs(x, 2 * c / endpoint) - reduced) == 0, "independent dictionary reduction failed")

    ctx.dps = 110
    ctx.threads = 1
    endpoint_ball, pi = arb(B), arb.pi()
    cmax, rho = endpoint_ball / 4, arb(15) / 16
    raw = {order: arb(0) for order in range(3)}
    for coefficient, p, mode_power, k in ((1, 6, 0, 5), (10, 4, 2, 5), (5, 2, 4, 5)):
        for order in range(3):
            for factor, cp, kp in rows(order, p, k):
                raw[order] += (
                    3 * coefficient / (pi**4 * endpoint_ball**2)
                    * factor * cmax**cp * rho**(-kp)
                    * tail(2 * kp - mode_power) * (endpoint_ball / 2) ** order
                )

    delta, t = arb("1e-4"), arb(T)
    x0 = (1 - (1 - 8 * t / (pi * endpoint_ball**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    weight = (delta * (1 - delta)) ** (-arb(1) / 4)
    l1 = (1 - 2 * delta) / (4 * delta * (1 - delta))
    w2 = (12 * delta**2 - 12 * delta + 5) / (16 * delta**2 * (1 - delta) ** 2)
    norm = {
        0: weight * raw[0],
        1: weight * (raw[1] + l1 * raw[0]) / hessian.sqrt(),
        2: weight * (raw[2] + 2 * l1 * raw[1] + w2 * raw[0]) / hessian,
    }
    stored_raw = artifact["interval_certificate"]["dictionary_x_derivative_tail_balls"]
    stored_norm = artifact["interval_certificate"]["weighted_normalized_derivative_tail_balls"]
    for order in range(3):
        require(arb(stored_raw[str(order)]).overlaps(raw[order]), f"raw order-{order} replay mismatch")
        require(arb(stored_norm[str(order)]).overlaps(norm[order]), f"normalized order-{order} replay mismatch")

    decision = artifact["decision"]
    require(decision["complete_outer_mode_derivative_block_bounded"] is True, "outer completion decision drift")
    require(decision["finite_mode_block_below_outer_start_bounded"] is False, "finite-block overclaim")
    require(decision["completed_exterior_current_bound_proved"] is False, "completed-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(
        "independently checked outer dictionary derivatives: "
        + ", ".join(f"order{order}={norm[order]}" for order in range(3)),
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
