#!/usr/bin/env python3
"""Validate the all-call exceptional-mode recombination certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_exceptional_mode_recombination_full_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["status"] == "rigorous_all_374_exceptional_mode_recombination_validated",
        "exceptional-mode full status drift",
    )
    require(len(artifact["rows"]) == 374, "exceptional-mode full row count drift")
    aggregate = artifact["aggregate"]
    require(aggregate["call_count"] == 374, "exceptional-mode full call count drift")
    require(aggregate["w3_call_count"] > 0 and aggregate["w4_call_count"] > 0, "selector coverage drift")
    require(aggregate["w3_call_count"] + aggregate["w4_call_count"] == 374, "selector count drift")
    require(aggregate["model_identity_count"] == 374, "model identity count drift")
    require(aggregate["recombination_identity_count"] == 374, "recombination identity count drift")
    require(float(aggregate["maximum_model_identity_gap_upper"]) < 1e-40, "model identity gap drift")
    require(float(aggregate["maximum_total_recombination_gap_upper"]) < 1e-10, "recombination gap drift")
    require(float(aggregate["maximum_exceptional_mode_tail_upper"]) < 1e-20, "mode tail drift")
    require(float(aggregate["maximum_component_triangle_to_correction_ratio_upper"]) > 1, "triangle diagnostic drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing exceptional-mode full dependency: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == dependency["sha256"], f"hash drift: {path}")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("all 374", "global `Phi2`", "no `1/delta`", "not a fitted theorem constant"):
        require(token in note, f"exceptional-mode full note token missing: {token}")
    print("validated exceptional-mode recombination full gate: 374/374 model and recombination identities")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
