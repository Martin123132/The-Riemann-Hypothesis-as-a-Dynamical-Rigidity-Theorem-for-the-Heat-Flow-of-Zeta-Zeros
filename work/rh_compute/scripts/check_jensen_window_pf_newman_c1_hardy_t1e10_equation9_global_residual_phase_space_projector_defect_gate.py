#!/usr/bin/env python3
"""Independent check of the phase-space projector-defect gate."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_phase_space_projector_defect_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sign(value: arb) -> int:
    if value.lower() > 0:
        return 1
    if value.upper() < 0:
        return -1
    return 0


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "production artifact is not passed")
    require(NOTE.is_file(), "missing theorem note")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file() and file_hash(path) == dependency["sha256"], f"dependency hash drift: {path.name}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    # Independent aggregate algebra.
    h, outside, target_interval, target_full = sp.symbols("h outside target_interval target_full")
    all_interval = outside + target_interval
    target_exterior = target_full - target_interval
    require(sp.expand((h + all_interval - target_full) - (h + outside - target_exterior)) == 0, "independent projector identity failed")

    m, t, endpoint = sp.symbols("m t endpoint", positive=True)
    alpha = 2 * m + t / (sp.pi * m)
    x_m = 2 * sp.pi * m**2 / (t + 2 * sp.pi * m**2)
    z_d = 2 * m / endpoint
    face_formula = 2 * sp.pi * m**2 * (endpoint - alpha) / (endpoint * (t + 2 * sp.pi * m**2))
    require(sp.simplify(x_m - z_d - face_formula) == 0, "independent face formula failed")
    require(sp.simplify(x_m - sp.Rational(1, 2) - (2 * sp.pi * m**2 - t) / (2 * (t + 2 * sp.pi * m**2))) == 0, "independent half-boundary formula failed")

    for row in artifact["mode_partition_certificate"]["rows"]:
        cutoff = row["cutoff"]
        target_count = 39_894 - 622 + 1
        low = cutoff + 622
        high = max(0, cutoff - 39_894)
        require(row["target_count"] == target_count == 39_273, "target count drift")
        require(row["low_outside_count"] == low and row["high_outside_count"] == high, "outside block count drift")
        require(row["outside_count"] == low + high == 2 * cutoff + 1 - target_count, "outside total drift")

    # Recompute all transition signs from the defining equations.
    ctx.dps = 90
    pi = arb.pi()
    height = arb(10_000_000_000)

    def alpha_value(mode: int) -> arb:
        mode_ball = arb(mode)
        return 2 * mode_ball + height / (pi * mode_ball)

    def face_gap(endpoint_value: int, mode: int) -> arb:
        mode_ball = arb(mode)
        return 2 * pi * mode_ball**2 * (arb(endpoint_value) - alpha_value(mode)) / (
            arb(endpoint_value) * (height + 2 * pi * mode_ball**2)
        )

    def half_gap(mode: int) -> arb:
        mode_ball = arb(mode)
        return (2 * pi * mode_ball**2 - height) / (2 * (height + 2 * pi * mode_ball**2))

    recomputed = {
        "B_face_m621": face_gap(5_122_421, 621),
        "B_face_m622": face_gap(5_122_421, 622),
        "A_face_m39852": face_gap(159_577, 39_852),
        "A_face_m39853": face_gap(159_577, 39_853),
        "half_boundary_m39894": half_gap(39_894),
        "half_boundary_m39895": half_gap(39_895),
        "ownership_A_face_m39694": face_gap(159_577, 39_694),
        "ownership_A_face_m39695": face_gap(159_577, 39_695),
        "ownership_half_m39694": half_gap(39_694),
        "ownership_half_m39695": half_gap(39_695),
    }
    expected_signs = {
        "B_face_m621": -1,
        "B_face_m622": 1,
        "A_face_m39852": -1,
        "A_face_m39853": 1,
        "half_boundary_m39894": -1,
        "half_boundary_m39895": 1,
        "ownership_A_face_m39694": -1,
        "ownership_A_face_m39695": -1,
        "ownership_half_m39694": -1,
        "ownership_half_m39695": -1,
    }
    for key, value in recomputed.items():
        require(sign(value) == expected_signs[key] == artifact["transition_certificate"]["signs"][key], f"transition sign drift: {key}")
        require(value.overlaps(arb(artifact["transition_certificate"]["values"][key])), f"transition interval drift: {key}")
    require(4 * 39_894 == 159_577 - 1 and 4 * 39_895 == 159_577 + 3, "null-face lattice bracket drift")

    symmetric_path = REPO_ROOT / artifact["dependencies"]["symmetric_Poisson"]["path"]
    symmetric = json.loads(symmetric_path.read_text(encoding="utf-8"))
    coefficient = arb(symmetric["certified_constants"]["tail_coefficient_C_AB_over_2pi2_ball"])
    at_target = coefficient / arb(39_894)
    at_outer = coefficient / arb(5_122_421)
    required = coefficient / arb("0.00014057919999999995")
    obstruction = artifact["raw_norm_obstruction"]
    require(at_target.overlaps(arb(obstruction["bound_at_M_39894_ball"])), "target-edge obstruction interval drift")
    require(at_outer.overlaps(arb(obstruction["bound_at_M_5122421_ball"])), "outer-threshold obstruction interval drift")
    require(required.overlaps(arb(obstruction["cutoff_ratio_needed_below_working_target_ball"])), "required-cutoff interval drift")
    require(at_outer.lower() > arb("0.00014057919999999995") and required.lower() > arb("1e29"), "raw-majorant obstruction failed")

    decision = artifact["decision"]
    require(decision["Gamma_normalized_residual_is_one_endpoint_completed_oriented_projector_defect"] is True, "projector decision drift")
    require(decision["three_saved_height_stationary_occupancy_interfaces_identified"] is True, "occupancy-interface decision drift")
    require(decision["ordinary_fold_ownership_interface_is_not_stationary_transition"] is True, "ownership-interface decision drift")
    require(decision["raw_C2_interchange_majorant_closes_working_target"] is False, "raw bound overclaim")
    require(decision["phase_adapted_signed_projector_bound_proved"] is False, "signed bound overclaim")
    require(decision["compressed_R_Dir_target_proved"] is False, "R_Dir overclaim")
    require(decision["rh_implication"] is False, "RH overclaim")
    print("independently checked endpoint-completed projector defect and transition audit", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
