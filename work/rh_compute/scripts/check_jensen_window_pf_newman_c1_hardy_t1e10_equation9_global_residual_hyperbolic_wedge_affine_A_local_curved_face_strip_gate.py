#!/usr/bin/env python3
"""Independently replay the finite-box affine A curved-face strip."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_curved_face_strip_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER_PRECISION = 70
CHECKER_SERIES_TERMS = 28
CHECKER_TOL = "3e-15"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_builder_module():
    spec = importlib.util.spec_from_file_location("A_face_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "logistic_null_coordinate_regularized_by_certified_power_series",
        "all_84_half_boundaries_contained_in_finite_Airy_box",
        "exact_curved_face_and_tangent_used_without_face_Taylor_truncation",
        "affine_inner_strip_reduced_to_exact_Fresnel_primitives",
        "signed_local_curved_face_strip_summed_before_norms",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    require(decisions.get("local_curved_face_strip_is_small_at_R_Dir_target_scale") is False, "local-strip scale guard drift")
    for key in (
        "exterior_curved_face_strip_bound_proved",
        "transformed_amplitude_remainder_bound_proved",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    priority = module.set_low_priority()
    require(priority == "below_normal", f"checker priority drift: {priority}")
    module.ctx.dps = CHECKER_PRECISION
    module.ctx.threads = 1
    saved_rows = {row["mode"]: row for row in artifact["certificate"]["rows"]}
    require(set(saved_rows) == set(range(39_853, 39_937)), "saved mode roster drift")
    fresh_rows: list[dict[str, Any]] = []
    for index, mode in enumerate(range(39_853, 39_937), start=1):
        fresh = module.evaluate_mode(
            mode,
            CHECKER_PRECISION,
            CHECKER_SERIES_TERMS,
            CHECKER_TOL,
        )
        fresh_rows.append(fresh)
        require(
            module.arb(fresh["physical_local_curved_face_strip_ball"]).overlaps(
                module.arb(saved_rows[mode]["physical_local_curved_face_strip_ball"])
            ),
            f"higher-precision face strip misses at mode {mode}",
        )
        require(
            module.acb(
                module.arb(fresh["local_strip_canonical_ball"]["real_ball"]),
                module.arb(fresh["local_strip_canonical_ball"]["imag_ball"]),
            ).overlaps(
                module.acb(
                    module.arb(saved_rows[mode]["local_strip_canonical_ball"]["real_ball"]),
                    module.arb(saved_rows[mode]["local_strip_canonical_ball"]["imag_ball"]),
                )
            ),
            f"higher-precision canonical strip misses at mode {mode}",
        )
        if index % 8 == 0:
            print(f"checked local curved-face rows {index}/84", flush=True)

    fresh_sums = module.aggregate_rows(fresh_rows)
    saved_sums = artifact["certificate"]["signed_sums"]
    for field in fresh_sums:
        require(
            module.arb(fresh_sums[field]).overlaps(module.arb(saved_sums[field])),
            f"higher-precision aggregate misses {field}",
        )
    require(module.arb(fresh_sums["complete_physical_local_curved_face_strip_ball"]) < module.arb("-0.009"), "signed strip scale drift")

    cache_path = REPO_ROOT / artifact["cache"]["path"]
    require(cache_path.is_file(), "missing resumable cache")
    require(file_hash(cache_path) == artifact["cache"]["sha256"], "cache hash drift")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("not a perturbation" in note, "local-strip scale interpretation missing")
    require("No exterior face" in note, "proof boundary missing")
    print("independently checked exact local affine A curved-face strip", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
