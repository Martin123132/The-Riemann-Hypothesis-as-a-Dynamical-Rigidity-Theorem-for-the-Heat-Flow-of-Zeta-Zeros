#!/usr/bin/env python3
"""Independent replay of the signed-pair triangle geometry."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_separable_triangle_geometry_gate"
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

    m, t, endpoint = sp.symbols("m t endpoint", positive=True, real=True)
    z_d = 2 * m / endpoint
    x_m = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    alpha_m = 2 * m + t / (sp.pi * m)
    gap = 2 * sp.pi * m**2 * (endpoint - alpha_m) / (endpoint * (t + 2 * sp.pi * m**2))
    require(sp.simplify(x_m - z_d - gap) == 0, "independent gap identity failed")

    rows = {row["mode"]: row for row in artifact["saved_height_ledger"]["rows"]}
    require(float(rows[621]["B_gap_x_minus_z"]) < 0 < float(rows[622]["B_gap_x_minus_z"]), "B face signs failed")
    require(float(rows[39_852]["A_gap_x_minus_z"]) < 0 < float(rows[39_853]["A_gap_x_minus_z"]), "A face signs failed")
    require(float(rows[39_894]["x_m"]) < 0.5 < float(rows[39_895]["x_m"]), "half face signs failed")

    decision = artifact["decision"]
    require(decision["pair_Volterra_integral_reduced_to_separable_triangle"] is True, "triangle decision drift")
    require(decision["endpoint_characteristic_is_triangle_face_crossing"] is True, "face decision drift")
    require(decision["tangent_B_profile_is_complete_pair_triangle_theorem"] is False, "tangent overclaim")
    require(decision["uniform_triangle_face_corner_estimate_proved"] is False, "uniform-bound overclaim")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")
    print("checked pair-triangle geometry by independent gap and face-sign replay", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
