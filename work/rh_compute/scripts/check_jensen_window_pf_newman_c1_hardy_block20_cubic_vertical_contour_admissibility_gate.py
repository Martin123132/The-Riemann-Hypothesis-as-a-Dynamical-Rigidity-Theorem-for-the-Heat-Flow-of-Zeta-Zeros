#!/usr/bin/env python3
"""Validate the cubic vertical-contour admissibility gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing cubic contour artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_cubic_vertical_contour_admissibility_gate",
        "cubic contour kind drift",
    )
    require(
        artifact["status"] == "exact_cubic_far_field_obstruction_to_stated_vertical_nonsaddle_contours",
        "cubic contour status drift",
    )
    scope = artifact["scope"]
    aggregate = artifact["aggregate"]
    rows = artifact["rows"]
    require(scope["recursive_call_count"] == len(rows) == 374, "cubic contour row count drift")
    require(scope["precision_ladder_decimal_digits"] == [90, 150], "cubic contour precision drift")
    require(scope["precision_overlap_count"] == 374, "cubic contour overlap count drift")
    require(aggregate["negative_phi3_count"] == aggregate["failing_67b_count"] == 187, "cubic contour 67b count drift")
    require(aggregate["positive_phi3_count"] == aggregate["failing_67c_count"] == 187, "cubic contour 67c count drift")
    require(aggregate["zero_phi3_count"] == 0, "cubic contour zero Phi3 drift")
    require(Decimal(aggregate["minimum_crossover_radius_lower"]) > Decimal("100"), "cubic contour minimum radius drift")
    require(Decimal(aggregate["maximum_crossover_radius_upper"]) < Decimal("1e7"), "cubic contour maximum radius drift")
    require(Decimal(aggregate["maximum_witness_identity_gap_upper"]) < Decimal("1e-130"), "cubic contour identity drift")
    require(Decimal(aggregate["maximum_witness_imaginary_phase_upper"]) < 0, "cubic contour growth witness drift")

    chain_roster = [row["chain"] for row in rows]
    require(chain_roster == sorted(chain_roster) and len(set(chain_roster)) == 374, "cubic contour row order drift")
    require(sum(row["failing_family"] == "67b" for row in rows) == 187, "cubic contour row 67b drift")
    require(sum(row["failing_family"] == "67c" for row in rows) == 187, "cubic contour row 67c drift")
    for row in rows:
        require(row["precision_overlap"], "cubic contour row precision drift")
        require(row["phi3_sign"] in ("negative", "positive"), "cubic contour row sign drift")
        require(
            (row["phi3_sign"] == "negative" and row["failing_family"] == "67b")
            or (row["phi3_sign"] == "positive" and row["failing_family"] == "67c"),
            "cubic contour sign/family mismatch",
        )
        require(Decimal(row["witness_identity_gap_upper"]) < Decimal("1e-130"), "cubic contour row identity drift")

    require(artifact["paper"]["pages"] == [27, 28, 29, 30], "cubic contour paper page drift")
    require(artifact["paper"]["equations"] == [82, 83, 84, 85, 87, 90, 91], "cubic contour equation drift")
    require(artifact["repair_target"]["candidate_ray_angle"] == "pi/6", "cubic contour repair ray drift")
    require(artifact["repair_target"]["leading_decay"] == "sin(3*pi/6)=1", "cubic contour decay identity drift")
    for group in (artifact["paper"], artifact["source"], *artifact["sources"].values()):
        path = REPO_ROOT / group["path"]
        require(path.is_file(), f"cubic contour dependency missing: {group['path']}")
        require(file_hash(path) == group["sha256"], f"cubic contour dependency hash drift: {group['path']}")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "cubic contour builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "cubic contour checker hash drift")

    for token in (
        "Status: exact cubic far-field obstruction",
        "Im g_n(x+iR)",
        "Im h_n(x+iR)",
        "does not even tend to zero",
        "187 calls; stated (67b)",
        "187 calls; stated (67c)",
        "pi/6",
        "This gate identifies a missing",
    ):
        require(token in note, f"cubic contour note token missing: {token}")
    require("does not supply the repaired contour" in artifact["proof_boundary"], "cubic contour boundary drift")
    print("validated cubic contour admissibility gate: 187 failing 67b, 187 failing 67c")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
