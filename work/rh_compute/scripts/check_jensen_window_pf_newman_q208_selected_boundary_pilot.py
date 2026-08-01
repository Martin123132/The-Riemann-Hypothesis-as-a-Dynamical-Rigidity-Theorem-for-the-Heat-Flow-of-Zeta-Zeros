#!/usr/bin/env python3
"""Validate the stored Q208 selected-boundary pilot."""

from __future__ import annotations

import json
from pathlib import Path

import jensen_window_pf_newman_q208_selected_boundary_pilot as pilot


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_q208_selected_boundary_pilot.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_q208_selected_boundary_pilot.md"
)


def validate() -> list[str]:
    issues: list[str] = []
    pilot.bridge.compact.flint.ctx.prec = 352
    if not RESULT.exists():
        return ["missing stored pilot result"]
    if not NOTE.exists():
        issues.append("missing rendered pilot note")
    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    if stored.get("kind") != (
        "jensen_window_pf_newman_q208_selected_boundary_pilot"
    ):
        issues.append("kind drifted")
    if stored.get("status") != (
        "rigorous complete Q208 right-strip theorem and "
        "selected-bottom boundary diagnostic"
    ):
        issues.append("status drifted")

    config = stored.get("config", {})
    if config.get("q207") != "[1/1035,1/4]x[0,245]":
        issues.append("Q207 source rectangle drifted")
    if config.get("q208") != "[1/1040,1/4]x[0,246]":
        issues.append("Q208 pilot rectangle drifted")
    expected_bottom = [
        [str(lower), str(upper)]
        for lower, upper in pilot.BOTTOM_PANELS
    ]
    expected_right = [
        [str(lower), str(upper)]
        for lower, upper in pilot.RIGHT_PANELS
    ]
    if config.get("bottom_panels") != expected_bottom:
        issues.append("bottom selected panels drifted")
    if config.get("right_panels") != expected_right:
        issues.append("right selected panels drifted")
    if config.get("precision_bits") != 352:
        issues.append("precision drifted")
    if config.get("retained_terms") != 6:
        issues.append("retained count drifted")

    source_paths = {
        "q207_builder": Path(pilot.bridge.__file__).resolve(),
        "q207_result": pilot.bridge.DEFAULT_OUT.resolve(),
        "six_term_tail": pilot.bridge.TAIL_SOURCE.resolve(),
        "outer_tail": pilot.bridge.OUTER_SOURCE.resolve(),
    }
    expected_hashes = {
        key: pilot.file_hash(path)
        for key, path in source_paths.items()
    }
    if stored.get("source_sha256") != expected_hashes:
        issues.append("source hashes drifted")

    expected_tasks = pilot.tasks()
    records = stored.get("records", [])
    if len(records) != 8:
        issues.append(f"expected 8 panel records, got {len(records)}")
    if [record.get("task") for record in records] != expected_tasks:
        issues.append("panel tasks or order drifted")
    expected_sequence = list(range(1, len(records) + 1))
    if [record.get("sequence") for record in records] != expected_sequence:
        issues.append("panel sequence drifted")

    leaf_count = 0
    value_count = 0
    derivative_count = 0
    for record in records:
        result = record.get("result", {})
        if result.get("status") != "certified":
            issues.append(
                f"panel {record.get('sequence')} is not certified"
            )
        if result.get("unresolved_boxes") != 0:
            issues.append(
                f"panel {record.get('sequence')} has unresolved boxes"
            )
        leaves = result.get("certified_leaves", [])
        if len(leaves) != result.get("certified_leaf_boxes"):
            issues.append(
                f"panel {record.get('sequence')} leaf count mismatch"
            )
        leaf_count += len(leaves)
        for leaf in leaves:
            ratio = pilot.bridge.compact.arb(
                leaf["certified_ratio_lower"]
            )
            if not ratio.is_finite() or ratio.lower() <= 1:
                issues.append(
                    f"panel {record.get('sequence')} has a nonpassing ratio"
                )
            branch = leaf.get("branch")
            if branch == "value":
                value_count += 1
            elif branch == "derivative":
                derivative_count += 1
            else:
                issues.append(
                    f"panel {record.get('sequence')} has invalid branch"
                )

    summary = stored.get("summary", {})
    if leaf_count != 46:
        issues.append(f"expected 46 certified leaves, got {leaf_count}")
    if summary.get("panels") != len(records):
        issues.append("summary panel count mismatch")
    if summary.get("certified_leaves") != leaf_count:
        issues.append("summary leaf count mismatch")
    if summary.get("unresolved_boxes") != 0:
        issues.append("summary has unresolved boxes")
    if summary.get("branch_counts") != {
        "value": value_count,
        "derivative": derivative_count,
    }:
        issues.append("summary branch counts mismatch")
    if not summary.get("all_selected_panels_certified"):
        issues.append("selected-panel completion flag is false")
    right = summary.get("right_strip", {})
    right_records = [
        record
        for record in records
        if record.get("task", {}).get("region")
        == "full_time_right_strip"
    ]
    right_leaves = [
        leaf
        for record in right_records
        for leaf in record.get("result", {}).get("certified_leaves", [])
    ]
    right_lowers = [
        pilot.full_derivative_lower(leaf)
        for leaf in right_leaves
    ]
    if len(right_leaves) != 40:
        issues.append(
            f"expected 40 complete right-strip cells, got {len(right_leaves)}"
        )
    if not right_lowers or not all(value > 0 for value in right_lowers):
        issues.append("full derivative is not positive on the right strip")
    if right.get("domain") != "[1/1040,1/5]x[245,246]":
        issues.append("right-strip theorem domain drifted")
    if right.get("certified_cells") != len(right_leaves):
        issues.append("right-strip certified-cell count mismatch")
    if not right.get("full_derivative_positive"):
        issues.append("right-strip derivative theorem flag is false")
    if right_lowers:
        recorded_lower = pilot.bridge.compact.arb(
            right.get("minimum_full_derivative_lower")
        )
        actual_lower = min(right_lowers)
        if not recorded_lower.contains(actual_lower):
            issues.append("right-strip minimum derivative lower drifted")

    resource = stored.get("resource", {})
    if resource.get("worker_count") != 1:
        issues.append("pilot did not record one worker")
    baseline = resource.get("baseline_samples", [])
    if len(baseline) != pilot.CPU_BASELINE_SECONDS:
        issues.append("baseline sample count drifted")
    if resource.get("baseline_mean", 100) > pilot.CPU_BASELINE_LIMIT:
        issues.append("pilot launched above its baseline CPU limit")

    boundary = stored.get("proof_boundary", "")
    for phrase in [
        "complete new right strip",
        "do not cover the complete bottom edge",
        "compute the closed boundary winding",
        "prove Lambda<=0",
        "prove RH",
    ]:
        if phrase not in boundary:
            issues.append(f"proof boundary missing phrase: {phrase}")

    if NOTE.exists():
        note = NOTE.read_text(encoding="utf-8")
        required = [
            "not a proof",
            "does not certify Q208",
            "complete cover",
            "open upper half-plane",
            "compute a closed-boundary winding",
            "boundary phase-unwrapping",
        ]
        for phrase in required:
            if phrase not in note:
                issues.append(f"note missing required phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    summary = stored["summary"]
    print(
        "validated Q208 selected-boundary pilot: "
        f"{summary['panels']} panels, "
        f"{summary['certified_leaves']} certified leaves, "
        f"{summary['unresolved_boxes']} unresolved, 0 global promotions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
