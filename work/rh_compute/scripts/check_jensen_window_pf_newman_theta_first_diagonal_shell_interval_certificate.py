#!/usr/bin/env python3
"""Validate the first diagonal-shell interval certificate."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_first_diagonal_shell_interval_certificate as target  # noqa: E402


EXPECTED_IDS = [
    "ntfdic_01_interval_method",
    "ntfdic_02_edge_partition",
    "ntfdic_03_edge_contact_exclusion",
    "ntfdic_04_compact_composition",
    "ntfdic_05_above_boundary_simplicity",
    "ntfdic_06_first_stage_theorem",
    "ntfdic_07_winding_composition",
    "ntfdic_08_second_stage_handoff",
]


def validate(payload: dict, note: str) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != target.STEM:
        issues.append("artifact kind mismatch")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if any(
        row.get("readiness") != "ready_to_apply" for row in rows[:-1]
    ):
        issues.append("certified row was demoted")
    if not rows or rows[-1].get("readiness") != "not_ready_to_apply":
        issues.append("second-stage handoff was promoted")
    if payload.get("exact") != target.build_exact():
        issues.append("exact payload drifted")
    if payload.get("source_audit") != target.source_audit():
        issues.append("source audit drifted")

    certificate = payload.get("certificate", {})
    expected_summary = {
        "precision_bits": 160,
        "integration_cutoff": "2",
        "tail_radius": "1e-800",
        "edge_time": "1/5",
        "x_lower": "38",
        "x_upper": "39",
        "initial_x_step": "1/5",
        "initial_boxes": 5,
        "evaluated_boxes": 5,
        "certified_boxes": 5,
        "subdivisions": 0,
        "unresolved_boxes": 0,
        "maximum_depth": 0,
        "branch_counts": {"value": 0, "derivative": 5},
    }
    for key, expected in expected_summary.items():
        if certificate.get(key) != expected:
            issues.append(
                f"certificate {key} drifted: "
                f"{certificate.get(key)!r} != {expected!r}"
            )
    if certificate.get("below_normal_priority_applied") is not True:
        issues.append("below-normal priority was not applied")
    try:
        minimum = target.compact.arb(
            certificate["minimum_certified_ratio_lower"]
        ).lower()
        if not minimum > target.compact.arb("1.2"):
            issues.append("minimum ratio is not above 6/5")
    except (KeyError, ValueError) as exc:
        issues.append(f"invalid minimum ratio: {exc}")

    records = certificate.get("records", [])
    if len(records) != 5:
        issues.append("edge record count mismatch")
    cursor = target.X_LOWER
    for index, row in enumerate(records):
        if row.get("certified") is not True:
            issues.append(f"record {index} is not certified")
        if row.get("branch") != "derivative":
            issues.append(f"record {index} changed branch")
        try:
            t_low = Fraction(row["t_low"])
            t_high = Fraction(row["t_high"])
            x_low = Fraction(row["x_low"])
            x_high = Fraction(row["x_high"])
        except (KeyError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"record {index} has invalid bounds: {exc}")
            continue
        if t_low != target.EDGE_TIME or t_high != target.EDGE_TIME:
            issues.append(f"record {index} leaves the exact time edge")
        if x_low != cursor or x_high - x_low != target.X_STEP:
            issues.append(f"record {index} breaks the edge partition")
        cursor = x_high
        try:
            ratio = target.compact.arb(
                row["certified_ratio_lower"]
            ).lower()
            if not ratio > 1:
                issues.append(f"record {index} lost strict separation")
        except (KeyError, ValueError) as exc:
            issues.append(f"record {index} has invalid ratio: {exc}")
    if cursor != target.X_UPPER:
        issues.append("edge partition does not end at x=39")

    status = payload.get("status", "")
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "closing j=1",
        "j>=2",
        "Lambda<=0",
        "RH remain open",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "Q_1 only",
        "does not certify Q_2",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    required_note = (
        "Certified Edge",
        "First Stage",
        "Proof Boundary",
        "38<=x<=39",
        "certified boxes=5",
        "unresolved boxes=0",
        "minimum normalized disjunction ratio",
        "first member of the independent diagonal exhaustion",
        "Stages `j>=2`,",
    )
    for marker in required_note:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    try:
        payload = json.loads(target.DEFAULT_OUT.read_text(encoding="utf-8"))
        note = target.DEFAULT_NOTE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as exc:
        print(f"first diagonal-shell certificate: source failure: {exc}")
        return 1
    issues = validate(payload, note)
    rebuilt = target.build_payload()
    if payload != rebuilt:
        issues.append("stored certificate differs from full Arb replay")
    success = target.success_line(payload)
    if success not in note:
        issues.append("success line missing from note")
    if issues:
        print(f"first diagonal-shell certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
