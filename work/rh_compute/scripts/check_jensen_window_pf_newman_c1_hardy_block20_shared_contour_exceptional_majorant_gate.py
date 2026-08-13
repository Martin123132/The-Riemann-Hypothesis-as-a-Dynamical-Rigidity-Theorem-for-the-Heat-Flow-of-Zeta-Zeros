#!/usr/bin/env python3
"""Validate the shared-contour exceptional-mode majorant certificate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_shared_contour_exceptional_majorant_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["status"] == "rigorous_gap_uniform_W2_W3_shared_contour_majorants_validated_all_374",
        "shared-contour status drift",
    )
    require(len(artifact["rows"]) == 374, "shared-contour row count drift")
    aggregate = artifact["aggregate"]
    require(aggregate["w2_call_count"] == aggregate["w2_majorant_count"] == 374, "W2 majorant count drift")
    require(aggregate["w3_call_count"] == aggregate["w3_majorant_count"] == 165, "W3 majorant count drift")
    require(float(aggregate["minimum_w2_majorant_slack_lower"]) > 0, "W2 majorant slack drift")
    require(float(aggregate["minimum_w3_majorant_slack_lower"]) > 0, "W3 majorant slack drift")
    require(float(aggregate["maximum_w2_majorant_to_actual_ratio_upper"]) > 1e4, "W2 conservatism diagnostic drift")
    require(float(aggregate["maximum_w3_majorant_to_actual_ratio_upper"]) > 1e5, "W3 conservatism diagnostic drift")
    require("delta" not in artifact["theorem"]["b"], "W2 theorem reintroduced delta")
    require("eta" not in artifact["theorem"]["c_W3"], "W3 theorem reintroduced eta")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing shared-contour dependency: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == dependency["sha256"], f"hash drift: {path}")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("Common-contour lemma", "7*pi/12", "5*pi/12", "Fourier normalization", "fitted approximations"):
        require(token in note, f"shared-contour note token missing: {token}")
    print("validated shared-contour exceptional majorants: 374/374 W2 and 165/165 W3 bounds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
