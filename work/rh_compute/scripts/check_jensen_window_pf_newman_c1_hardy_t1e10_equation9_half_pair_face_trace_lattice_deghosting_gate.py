#!/usr/bin/env python3
"""Independent replay of the exact face-trace lattice geometry."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_face_trace_lattice_deghosting_gate"
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

    # Replay the inverse-root identity directly from the quadratic formula.
    endpoint, t = sp.symbols("endpoint t", positive=True, real=True)
    root = (endpoint - sp.sqrt(endpoint**2 - 8 * t / sp.pi)) / 4
    alpha = 2 * root + t / (sp.pi * root)
    require(sp.simplify(alpha - endpoint) == 0, "inverse-root replay failed")

    mp.mp.dps = 100
    T = mp.mpf(10) ** 10
    pi = mp.pi

    def center(D: int) -> mp.mpf:
        endpoint_value = mp.mpf(D)
        return (endpoint_value - mp.sqrt(endpoint_value**2 - 8 * T / pi)) / 4

    b = center(5_122_421)
    a = center(159_577)
    require(mp.mpf("621.55") < b < mp.mpf("621.56"), "B center replay failed")
    require(mp.mpf("39852.39") < a < mp.mpf("39852.40"), "A center replay failed")
    normal_257 = pi * mp.mpf(5_122_421) ** 2 / 4 * (1 - (mp.mpf(257) / b) ** 2)
    require(abs(normal_257) > mp.mpf("1.70e13"), "mode-257 normal replay failed")

    decision = artifact["decision"]
    require(decision["triangle_face_trace_is_mode_independent"] is True, "trace decision drift")
    require(decision["B_face_local_integer_roster_is_621_622"] is True, "B roster decision drift")
    require(decision["mode_257_is_true_face_critical_point"] is False, "mode-257 overclaim")
    require(decision["B_nonlocal_normal_IBP_bound_proved"] is False, "IBP overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked face-trace lattice geometry and mode-257 deghosting at 100 digits", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
