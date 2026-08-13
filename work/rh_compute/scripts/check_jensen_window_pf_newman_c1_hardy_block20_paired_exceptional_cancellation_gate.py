#!/usr/bin/env python3
"""Validate the paired exceptional-cancellation audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_paired_exceptional_cancellation_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact["status"] == "rigorous_all_374_paired_exceptional_cancellation_audited", "paired status drift")
    require(len(artifact["rows"]) == 374, "paired row count drift")
    aggregate = artifact["aggregate"]
    require(aggregate["call_count"] == 374, "paired call count drift")
    require(aggregate["paired_identity_count"] == 374, "paired identity count drift")
    require(float(aggregate["maximum_paired_identity_gap_upper"]) < 1e-10, "paired identity gap drift")
    require(float(aggregate["maximum_raw_component_ratio_upper"]) > 1e5, "raw cancellation diagnostic drift")
    ratio = float(aggregate["maximum_paired_triangle_to_correction_ratio_upper"])
    require(200 < ratio < 211, "paired cancellation ratio drift")
    require(float(aggregate["raw_to_paired_worst_ratio_reduction_lower"]) > 600, "paired reduction drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing paired dependency: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == dependency["sha256"], f"hash drift: {path}")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("all 374", "above 210", "not fitted constants", "not a height-uniform"):
        require(token in note, f"paired note token missing: {token}")
    print("validated paired exceptional-cancellation gate: 374/374 identities, worst paired ratio below 211")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
