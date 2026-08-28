#!/usr/bin/env python3
"""Independently replay the affine A-corner wedge evaluation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_null_coordinate_affine_corner_evaluation_gate"
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
    spec = importlib.util.spec_from_file_location("affine_corner_builder", BUILDER)
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
        "epsilon_extrapolation_avoided",
        "Abel_half_jump_retained_by_lower_indentation",
        "canonical_scalar_wedge_values_certified",
        "canonical_affine_wedge_values_certified",
        "mandatory_affine_current_quantitatively_nonzero",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in ("curved_face_remainder_bound_proved", "transformed_amplitude_remainder_bound_proved", "R_Dir_bound_proved", "rh_implication"):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    module.ctx.dps = 115
    affine_record = load_json(REPO_ROOT / artifact["dependencies"]["affine_remainder"]["path"])
    fresh_rows = module.evaluate_rows(affine_record, CHECKER_CUTOFF)
    fresh = {row["mode"]: row for row in fresh_rows}
    saved = {row["mode"]: row for row in artifact["rows"]}
    require(set(fresh) == set(saved) == {39_894, 39_895}, "mode roster drift")
    for mode in fresh:
        for field in (
            "canonical_wedge_C",
            "partial_a_C",
            "partial_h_C",
            "affine_wedge_C_aff",
            "affine_minus_scalar_wedge",
            "W0_times_affine_wedge",
        ):
            require(
                parse_complex(module, fresh[mode][field]).overlaps(parse_complex(module, saved[mode][field])),
                f"independent longer-ray replay misses {field} at mode {mode}",
            )
        require(module.arb(fresh[mode]["first_ray_tail_absolute_bound"]).upper() < module.arb("1e-80"), f"longer first-ray tail failed at {mode}")
    require(abs(module.arb(saved[39_894]["affine_minus_scalar_wedge"]["real_ball"])).lower() > module.arb("6e-4"), "mandatory affine correction drift")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("avoids epsilon extrapolation" in note, "note method missing")
    require("No exact curved-face" in note, "note proof boundary missing")
    print("independently checked null-coordinate affine A-corner wedge values", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
