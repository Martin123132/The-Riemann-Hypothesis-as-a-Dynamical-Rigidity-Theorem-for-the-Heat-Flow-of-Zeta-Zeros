#!/usr/bin/env python3
"""Independent replay of the outer Volterra-kernel bound."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_volterra_incomplete_gamma_outer_kernel_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")

    # Recompute the Jacobian cancellation independently.
    x = sp.symbols("x", positive=True, real=True)
    w = (1 - x) / x
    dx_dw = -1 / (1 + sp.symbols("w", positive=True)) ** 2
    w_symbol = sp.symbols("w", positive=True)
    x_of_w = 1 / (1 + w_symbol)
    transformed = sp.simplify(
        x_of_w ** (-sp.Rational(3, 2))
        * (x_of_w * (1 - x_of_w)) ** (-sp.Rational(1, 4))
        * (-dx_dw)
    )
    require(sp.simplify(transformed - w_symbol ** (-sp.Rational(1, 4))) == 0, "Jacobian replay failed")
    require(sp.simplify(w - (1 / x - 1)) == 0, "logistic identity failed")

    mp.mp.dps = 100
    t = mp.mpf(10) ** 10
    pi = mp.pi
    start = mp.mpf(39_895)
    nu = mp.sqrt(t / (2 * pi))
    first = 2 / (pi * (start**2 - nu**2))
    tail = mp.log((start + nu) / (start - nu)) / (pi * nu)
    bound = first + tail
    require(first < mp.mpf("0.00001034"), "first-mode bound replay failed")
    require(bound < mp.mpf("0.000103"), "outer l1 bound replay failed")

    decision = artifact["decision"]
    require(decision["Volterra_kernel_reduced_exactly_to_incomplete_Gamma_path"] is True, "kernel decision drift")
    require(decision["outer_kernel_family_l1_bound_proved"] is True, "l1 decision drift")
    require(decision["endpoint_driver_may_be_bounded_by_raw_absolute_value"] is False, "driver overclaim")
    require(decision["outer_pair_residual_quantitatively_closed"] is False, "residual overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked outer Volterra kernel by alternate Jacobian and 100-digit tail replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
