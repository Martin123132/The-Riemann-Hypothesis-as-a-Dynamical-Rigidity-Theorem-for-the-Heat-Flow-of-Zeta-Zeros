#!/usr/bin/env python3
"""Validate the W4 contour and complete endpoint Hurwitz majorant gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_w4_contour_endpoint_hurwitz_majorant_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["status"] == "rigorous_W4_contour_and_complete_selector_gap_uniform_endpoint_majorants_validated_all_374",
        "W4/Hurwitz status drift",
    )
    require(len(artifact["rows"]) == 374, "W4/Hurwitz row count drift")
    aggregate = artifact["aggregate"]
    require(aggregate["call_count"] == 374, "W4/Hurwitz call count drift")
    require(aggregate["w4_contour_identity_count"] == 209, "W4 contour count drift")
    require(aggregate["upper_majorant_count"] == 374, "upper majorant count drift")
    require(aggregate["lower_majorant_count"] == 374, "lower majorant count drift")
    require(aggregate["complete_majorant_count"] == 374, "complete majorant count drift")
    require(float(aggregate["maximum_w4_contour_identity_gap_upper"]) < 1e-40, "W4 contour identity gap drift")
    require(float(aggregate["maximum_w4_half_ray_tail_upper"]) < 1e-30, "W4 half-ray tail drift")
    require(float(aggregate["minimum_upper_majorant_slack_lower"]) > 0, "upper majorant slack drift")
    require(float(aggregate["minimum_lower_majorant_slack_lower"]) > 0, "lower majorant slack drift")
    require(float(aggregate["minimum_complete_majorant_slack_lower"]) > 0, "complete majorant slack drift")
    ratio = float(aggregate["maximum_complete_majorant_to_correction_ratio_upper"])
    require(2700 < ratio < 2900, "complete majorant ratio drift")
    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing W4/Hurwitz dependency: {path}")
        require(hashlib.sha256(path.read_bytes()).hexdigest() == dependency["sha256"], f"hash drift: {path}")
    note = NOTE.read_text(encoding="utf-8")
    for token in ("W4 is an endpoint half-ray difference", "Hurwitz", "Fourier phase", "not fitted"):
        require(token in note, f"W4/Hurwitz note token missing: {token}")
    print("validated W4/Hurwitz endpoint majorants: 209/209 W4 contours and 374/374 complete bounds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
