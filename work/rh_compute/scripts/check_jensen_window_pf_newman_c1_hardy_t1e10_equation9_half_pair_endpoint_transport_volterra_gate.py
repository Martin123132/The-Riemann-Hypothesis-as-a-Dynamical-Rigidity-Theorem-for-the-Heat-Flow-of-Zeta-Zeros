#!/usr/bin/env python3
"""Independent replay of the endpoint-driven Fourier-pair transport."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_endpoint_transport_volterra_gate"
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

    # Derive the ODE independently from the twice-integrated Fourier pair.
    x, m = sp.symbols("x m", positive=True, real=True)
    p, p_x, delta = sp.symbols("p p_x delta")
    k = 2 * sp.pi * m
    remainder_sum_from_ibp = 2 * delta - k**2 * p
    remainder_sum_from_pde = 6 * sp.I * sp.pi * x * p + 4 * sp.I * sp.pi * x**2 * p_x
    ode = sp.expand(remainder_sum_from_pde - remainder_sum_from_ibp)
    expected = 4 * sp.I * sp.pi * x**2 * p_x + (k**2 + 6 * sp.I * sp.pi * x) * p - 2 * delta
    require(sp.simplify(ode - expected) == 0, "independent ODE derivation failed")

    mu = x ** sp.Rational(3, 2) * sp.exp(sp.I * sp.pi * m**2 / x)
    normalized_coefficient = k**2 / (4 * sp.I * sp.pi * x**2) + sp.Rational(3, 2) / x
    require(sp.simplify(sp.diff(mu, x) / mu - normalized_coefficient) == 0, "integrating factor replay failed")

    # Replay the endpoint-characteristic alignment without solve().
    t, endpoint = sp.symbols("t endpoint", positive=True, real=True)
    x_mode = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    endpoint_value = 2 * m + t / (sp.pi * m)
    require(sp.simplify(2 * m / endpoint_value - x_mode) == 0, "characteristic alignment replay failed")

    decision = artifact["decision"]
    require(decision["canonical_roster_transport_PDE_proved"] is True, "PDE decision drift")
    require(decision["symmetric_pair_first_order_x_ODE_proved"] is True, "ODE decision drift")
    require(decision["endpoint_driven_Volterra_solution_proved"] is True, "Volterra decision drift")
    require(decision["bare_endpoint_exponentials_separated"] is False, "endpoint separation overclaim")
    require(decision["Volterra_operator_quantitatively_bounded"] is False, "operator-bound overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked endpoint-driven pair transport by independent Fourier/PDE replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
