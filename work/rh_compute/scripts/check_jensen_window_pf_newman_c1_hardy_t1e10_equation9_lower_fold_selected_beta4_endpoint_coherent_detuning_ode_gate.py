#!/usr/bin/env python3
"""Independently check the selected beta^-4 endpoint-coherent detuning ODE."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
C = 159_577
Q = (C - 1) // 4


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def weights(lam: arb, y: arb, beta: arb) -> tuple[arb, arb, arb, arb]:
    e = 1 / beta**2
    u = 1 - e * (13 * lam + 3 * y) / 60
    u -= e**2 * (
        448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2
        - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    ) / 50400
    v = e * (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    v += e**2 * (40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27) / 1680
    uy = -e / 20 - e**2 * (
        840 * lam**2 * y**2 - 420 * lam * y**3 + 30 * lam + 315 * y**4 - 810 * y
    ) / 50400
    vy = e * (-4 * lam + 6 * y) / 60 + e**2 * (-20 * lam**2 + 2 * lam * y - 27 * y**2) / 1680
    return u, v, uy, vy


def branch(lam: arb, y: arb, beta: arb, sigma: int) -> tuple[acb, acb]:
    x = lam + y
    ai, aip, bi, bip = (-x).airy()
    w = acb(ai, -sigma * bi)
    wx = acb(-aip, sigma * bip)
    u, v, uy, vy = weights(lam, y, beta)
    return (u * w + v * wx) / 2, ((uy - x * v) * w + (u + vy) * wx) / 2


def direct_endpoint_sums(lam: arb, y: arb, beta: arb, upper: bool) -> acb:
    h = 4 * beta / C
    total = acb(0)
    for mode in range(39_696, 40_095):
        if mode == Q:
            continue
        d = h * (arb(mode) - arb(C) / 4)
        sigma = -1 if mode < Q else 1
        value, derivative = branch(lam, y, beta, sigma)
        if upper:
            phase = acb(0, y**2 / (4 * beta) - d * y).exp()
            total += phase * (acb(0, y / (2 * beta) - d) * value - derivative)
        else:
            total += derivative + acb(0, d) * value
    return total


def stored_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> None:
    require(RESULT.is_file(), "missing result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "stored gate is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path}")

    lam, y, e = sp.symbols("lambda y epsilon", real=True)
    x = lam + y
    u = 1 - e * (13 * lam + 3 * y) / 60
    u -= e**2 * (
        448 * lam**5 + 280 * lam**2 * y**3 + 4565 * lam**2
        - 105 * lam * y**4 + 30 * lam * y + 63 * y**5 - 405 * y**2
    ) / 50400
    v = e * (8 * lam**2 - 4 * lam * y + 3 * y**2) / 60
    v += e**2 * (40 * lam**3 - 20 * lam**2 * y + lam * y**2 - 9 * y**3 + 27) / 1680
    q0 = sp.expand(sp.diff(u, y, 2) - v - 2 * x * sp.diff(v, y))
    q1 = sp.expand(2 * sp.diff(u, y) + sp.diff(v, y, 2))
    require(sp.factor(q0 + e * y**2 / 4 - e**2 * (13 * lam * y**2 / 240 + y**3 / 80)) == 0, "independent W forcing failed")
    require(sp.factor(q1 + e**2 * y**2 * (8 * lam**2 - 4 * lam * y + 3 * y**2) / 240) == 0, "independent WX forcing failed")

    ctx.dps = 100
    pi = arb.pi()
    beta = (pi * arb(C) ** 2 / 8) ** (arb(1) / 3)
    y_max = pi * C / (2 * beta)
    lambda_max = pi / (16 * beta)
    lambda_ball = arb(lambda_max / 2, lambda_max / 2)
    lower_direct = direct_endpoint_sums(lambda_ball, arb(0), beta, False)
    upper_direct = direct_endpoint_sums(lambda_ball, y_max, beta, True)
    certificate = artifact["interval_certificate"]
    lower_stored = stored_complex(certificate["remaining_398_lower_endpoint_source_ball"])
    upper_stored = stored_complex(certificate["remaining_398_upper_endpoint_source_ball"])
    require(lower_direct.real.overlaps(lower_stored.real) and lower_direct.imag.overlaps(lower_stored.imag), "lower source reconstruction missed stored ball")
    require(upper_direct.real.overlaps(upper_stored.real) and upper_direct.imag.overlaps(upper_stored.imag), "upper source reconstruction missed stored ball")
    require(abs(lower_direct).overlaps(arb(certificate["remaining_398_lower_endpoint_source_absolute_ball"])), "lower absolute enclosure drift")
    require(abs(upper_direct).overlaps(arb(certificate["remaining_398_upper_endpoint_source_absolute_ball"])), "upper absolute enclosure drift")
    require((4 * beta / C * y_max - 2 * pi).contains(0), "independent hY identity failed")
    print("independently checked endpoint-coherent beta^-4 detuning ODE", flush=True)


if __name__ == "__main__":
    main()
