#!/usr/bin/env python3
"""Independent replay of the higher B face-current window bounds."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_higher_boundary_current_absolute_window_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()

T = 10_000_000_000
B = 5_122_421
CUTOFF = 200_000


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Reconstruct the second boundary term directly, rather than importing the
    # production reduced formula.
    z, x, endpoint, mode = sp.symbols("z x B m", positive=True, real=True)
    phase = sp.pi * mode**2 / z + sp.pi * endpoint**2 * z / 4
    h = z ** (-sp.Rational(1, 2)) * (2 + sp.I * sp.pi * endpoint**2 * z)
    q = sp.I * sp.diff(phase, z)
    direct_second = -x ** (-sp.Rational(3, 2)) * sp.diff(h / q, z).subs(z, x) / (2 * sp.I * sp.pi * q.subs(z, x))
    require(direct_second.has(mode, endpoint, x), "independent second-current derivation failed")

    ctx.dps = 100
    ctx.threads = 1
    t, endpoint_ball, pi = arb(T), arb(B), arb.pi()
    x0 = (1 - (1 - 8 * t / (pi * endpoint_ball**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = hessian.sqrt()
    x_low, x_high = x0 - arb(70) / root_hessian, x0 + arb(70) / root_hessian
    c_low, c_high = endpoint_ball * x_low / 2, endpoint_ball * x_high / 2

    absolute = {order: arb(0) for order in range(2, 6)}
    for mode_value in range(1, 621):
        d = c_low**2 - arb(mode_value) ** 2
        for order in absolute:
            absolute[order] += d ** (-order)
    for mode_value in range(623, CUTOFF + 1):
        d = arb(mode_value) ** 2 - c_high**2
        for order in absolute:
            absolute[order] += d ** (-order)
    rho = 1 - (c_high / arb(CUTOFF + 1)) ** 2
    for order in absolute:
        absolute[order] += rho ** (-order) * arb(CUTOFF) ** (1 - 2 * order) / (2 * order - 1)

    c = c_high
    # Use Euclidean coefficient norms here, unlike production's componentwise
    # triangle bounds.
    norm = lambda real, imag: (real**2 + imag**2).sqrt()
    sup1 = c / (pi**3 * endpoint_ball) * (
        norm(4 * pi * endpoint_ball * c**3, 4 * c**2) * absolute[3]
        + norm(5 * pi * endpoint_ball * c, 3) * absolute[2]
    )
    sup2 = c**2 / (pi**4 * endpoint_ball**2) * (
        norm(48 * pi * endpoint_ball * c**5, 48 * c**4) * absolute[5]
        + norm(84 * pi * endpoint_ball * c**3, 60 * c**2) * absolute[4]
        + norm(35 * pi * endpoint_ball * c, 15) * absolute[3]
    )
    physical = 2 * (pi / (32 * t)) ** (arb(1) / 4)
    width = arb(140) / root_hessian
    require(x_high < arb(1) / 2, "window left the decreasing-weight branch")
    weight = (x_low * (1 - x_low)) ** (-arb(1) / 4)
    bound1, bound2 = physical * width * weight * sup1, physical * width * weight * sup2
    require(bound1 < arb("1.49e-6"), "independent second-current bound failed")
    require(bound2 < arb("2.1e-10"), "independent third-current bound failed")
    require(bound1 + bound2 < arb("1.491e-6"), "independent combined bound failed")

    decision = artifact["decision"]
    require(decision["remainder_after_third_face_current_bounded"] is False, "remainder overclaim")
    require(decision["complete_B_face_estimate_proved"] is False, "complete-B overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print(f"independently checked higher B face currents: C1={bound1}, C2={bound2}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
