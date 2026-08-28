#!/usr/bin/env python3
"""Independently replay the carrier-oriented affine A-transition roster."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_transition_roster_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER_CUTOFF = 22


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_builder_module():
    spec = importlib.util.spec_from_file_location("A_roster_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_complex(module, record: dict[str, str]):
    return module.acb(module.arb(record["real_ball"]), module.arb(record["imag_ball"]))


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "complete_84_mode_A_transition_roster_evaluated",
        "target_and_outer_sides_each_have_42_modes",
        "delta_negative_and_positive_branches_joined_without_c_denominator_instability",
        "all_steepest_rays_point_away_from_their_quadratic_saddles",
        "common_phase_and_paired_endpoint_orientation_retained",
        "signed_affine_tangent_wedge_sum_certified_before_norms",
        "termwise_absolute_value_shown_to_lose_cancellation",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "curved_face_remainder_bound_proved",
        "transformed_amplitude_remainder_bound_proved",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    fresh = module.evaluate_roster(CHECKER_CUTOFF, precision=115)
    saved = artifact["certificate"]
    require(fresh["mode_count"] == saved["mode_count"] == 84, "roster count drift")
    require(fresh["target_side_count"] == saved["target_side_count"] == 42, "target-side count drift")
    require(fresh["outer_side_count"] == saved["outer_side_count"] == 42, "outer-side count drift")
    fresh_rows = {row["mode"]: row for row in fresh["rows"]}
    saved_rows = {row["mode"]: row for row in saved["rows"]}
    require(set(fresh_rows) == set(saved_rows) == set(range(39_853, 39_937)), "mode roster mismatch")
    corner = load_json(
        REPO_ROOT
        / artifact["dependencies"]["orientation"]["path"]
    )
    corner_source = load_json(
        REPO_ROOT
        / corner["dependencies"]["corner"]["path"]
    )
    certified_corner_rows = {row["mode"]: row for row in corner_source["rows"]}
    for mode in fresh_rows:
        for field in ("canonical_wedge_C", "affine_wedge_C_aff", "normalized_affine_half_mode"):
            require(parse_complex(module, fresh_rows[mode][field]).overlaps(parse_complex(module, saved_rows[mode][field])), f"longer-ray replay misses {field} at {mode}")
        for field in (
            "physical_scalar_A_tangent_wedge_ball",
            "physical_affine_A_tangent_wedge_ball",
            "physical_affine_correction_ball",
        ):
            require(module.arb(fresh_rows[mode][field]).overlaps(module.arb(saved_rows[mode][field])), f"longer-ray replay misses {field} at {mode}")
        require(saved_rows[mode]["paired_endpoint_projector_coefficient"] == 1, f"endpoint orientation drift at {mode}")
    for mode in (39_894, 39_895):
        for field in ("canonical_wedge_C", "affine_wedge_C_aff"):
            require(
                parse_complex(module, saved_rows[mode][field]).overlaps(
                    parse_complex(module, certified_corner_rows[mode][field])
                ),
                f"full-roster branch misses prior corner certificate for {field} at {mode}",
            )

    for field in (
        "target_side_physical_affine_A_tangent_wedge_ball",
        "outer_side_physical_affine_A_tangent_wedge_ball",
        "complete_physical_scalar_A_tangent_wedge_ball",
        "complete_physical_affine_A_tangent_wedge_ball",
        "complete_physical_affine_correction_ball",
        "termwise_absolute_physical_affine_sum_ball",
        "signed_to_termwise_absolute_ratio_ball",
    ):
        require(module.arb(fresh["signed_sums"][field]).overlaps(module.arb(saved["signed_sums"][field])), f"longer-ray replay misses aggregate {field}")
    require(module.arb(fresh["maximum_ray_tail_absolute_bound"]).upper() < module.arb("1e-90"), "longer-ray tail target failed")
    require(module.arb(saved["signed_sums"]["signed_to_termwise_absolute_ratio_ball"]).upper() < 1, "cancellation ratio drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("already contains the Abel" in note, "stable positive-delta half-jump convention missing")
    require("No exact curved-face" in note, "proof boundary missing")
    print("independently checked carrier-oriented affine A-transition roster", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
