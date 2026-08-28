#!/usr/bin/env python3
"""Independently check the complete translated A-face assembly."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from flint import arb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_complete_absolute_assembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> int:
    ctx.dps = 150
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decision = artifact["decision"]
    for key in (
        "common_x16_handoff_verified",
        "common_physical_normalization_verified",
        "complete_translated_42_mode_A_face_absolute_bound_below_0_point_00364_proved",
        "remaining_R_after_A_allowance_above_0_point_03317_bookkept",
    ):
        require(decision.get(key) is True, f"missing assembly decision: {key}")
    for key in (
        "oscillatory_cancellation_used",
        "other_post_A_channels_bounded",
        "R_after_A_bound_proved",
        "R_Dir_bound_proved",
        "QK_minus_T_bound_proved",
        "rh_implication",
    ):
        require(decision.get(key) is False, f"proof-boundary drift: {key}")

    dependencies = {}
    for name, record in artifact["dependencies"].items():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {name}")
        dependencies[name] = load_json(path)
    core = dependencies["pre_endpoint_core"]
    endpoint = dependencies["endpoint_layer"]
    require(core["scope"]["height"] == endpoint["scope"]["height"] == 10_000_000_000, "height mismatch")
    require(core["scope"]["A"] == endpoint["scope"]["A"] == 159_577, "A mismatch")
    x_core = arb(core["interval_certificate"]["exact_x_16_ball"])
    x_endpoint = arb(endpoint["interval_certificate"]["exact_x_16_ball"])
    require(x_core.overlaps(x_endpoint), "independent x_16 handoff failed")

    core_bound = arb(core["interval_certificate"]["physical_complete_pre_endpoint_core_absolute_bound"])
    endpoint_bound = arb(endpoint["interval_certificate"]["physical_complete_endpoint_bound"])
    independent_total = core_bound + endpoint_bound
    independent_remaining = arb("0.0368147039947") - independent_total
    certificate = artifact["assembly_certificate"]
    recorded_total = arb(certificate["complete_translated_A_face_absolute_bound"])
    require(recorded_total.overlaps(independent_total), "assembly total drift")
    require(independent_total.upper() < arb("0.00364"), "independent assembly threshold failed")
    require(independent_remaining.lower() > arb("0.03317"), "independent remaining allowance failed")

    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(BUILDER.is_file() and NOTE.is_file(), "builder or note missing")
    note = NOTE.read_text(encoding="utf-8")
    require("|I_(A,42)|" in note and "< 0.00364" in note, "assembly theorem missing")
    require("budget bookkeeping only" in note, "remaining-budget boundary missing")
    require("No" in note and "joined `R_after_A`" in note and "RH" in note, "proof boundary missing")
    print(
        "independently checked complete translated 42-mode A-face assembly; "
        f"bound {independent_total.str(14, more=True)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
