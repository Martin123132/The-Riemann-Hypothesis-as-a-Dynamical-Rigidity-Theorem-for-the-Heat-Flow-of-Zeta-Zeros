#!/usr/bin/env python3
"""Independent symbolic replay of the finite endpoint-tail decomposition."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_finite_endpoint_tail_decomposition_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "endpoint-tail artifact is not passed")
    require(artifact["sources"]["builder"]["sha256"] == file_hash(BUILDER), "builder hash drift")
    require(artifact["sources"]["checker"]["sha256"] == file_hash(CHECKER), "checker hash drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(dependency["sha256"] == file_hash(path), f"dependency hash drift: {path}")

    f_a, f_b, e_a, e_b, coefficient, p = sp.symbols("f_a f_b e_a e_b coefficient p")
    full = 1 + sp.I
    finite = coefficient * (f_b - f_a) + (e_b - e_a) / (sp.I * p)
    pieces = (
        coefficient * full
        - e_a / (sp.I * p)
        - coefficient * (f_a + full / 2)
        + e_b / (sp.I * p)
        - coefficient * (full / 2 - f_b)
    )
    require(sp.expand(finite - pieces) == 0, "independent current decomposition failed")

    x, d, m = sp.symbols("x d m", positive=True)
    q_squared_over_two = x * (d - 2 * m / x) ** 2 / 4
    require(sp.expand(-m**2 / x + q_squared_over_two - (x * d**2 / 4 - m * d)) == 0, "phase replay failed")
    t = sp.symbols("t", positive=True)
    endpoint_phase = sp.pi * d**2 * x / 4 + t * sp.log((1 - x) / x) / 2
    target_derivative = sp.pi * d**2 / 4 - t / (2 * x * (1 - x))
    require(sp.simplify(sp.diff(endpoint_phase, x) - target_derivative) == 0, "endpoint derivative replay failed")

    require(artifact["decision"]["bare_endpoint_current_may_be_summed_independently"] is False, "summation guard overclaim")
    require(artifact["decision"]["complete_symmetric_endpoint_tail_bound_proved"] is False, "endpoint-bound overclaim")
    require(artifact["decision"]["complete_T_upper_proved"] is False, "T_upper overclaim")
    print("checked exact endpoint-tail decomposition and symmetric summation guard", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
