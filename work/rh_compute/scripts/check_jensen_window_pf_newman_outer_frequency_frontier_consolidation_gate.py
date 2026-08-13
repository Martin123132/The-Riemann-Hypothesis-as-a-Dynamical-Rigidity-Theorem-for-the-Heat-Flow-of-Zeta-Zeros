#!/usr/bin/env python3
"""Validate the Newman outer-frequency frontier consolidation gate."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_outer_frequency_frontier_consolidation_gate as target  # noqa: E402


def validate(payload: dict, note: str) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != target.STEM:
        issues.append("artifact kind drifted")
    if payload.get("schema_version") != 1:
        issues.append("schema version drifted")
    status = payload.get("status", "")
    for phrase in (
        "compact-base and outer-frequency frontier consolidation",
        "outer Xi joint-avoidance theorem remains open",
    ):
        if phrase not in status:
            issues.append(f"status boundary missing: {phrase}")

    exact = payload.get("exact", {})
    compact = exact.get("compact_core", {})
    if compact != {
        "domain": "0<=t<=1/5 and |x|<=38",
        "certified_boxes": 1900,
        "unresolved_boxes": 0,
        "contact_free": True,
        "new_boundary_scout_required": False,
    }:
        issues.append("compact-core summary drifted")
    finite = exact.get("finite_cofinal_base", {})
    if finite.get("proved_through") != 208:
        issues.append("proved cofinal base drifted")
    if finite.get("first_unresolved_linear_stage") != 209:
        issues.append("first unresolved stage drifted")
    if finite.get("q209_open_shell_regions") != 2:
        issues.append("Q209 open-shell count drifted")
    if exact.get("direct_remaining_domain") != "0<t<=1/5 and x>38":
        issues.append("direct outer domain drifted")

    rows = payload.get("rows", [])
    expected_ids = [f"ofc_{index:02d}_{suffix}" for index, suffix in (
        (1, "compact_core_closed"),
        (2, "q208_finite_base"),
        (3, "q209_exact_frontier"),
        (4, "direct_outer_reduction"),
        (5, "local_band_placement"),
        (6, "compact_rescout_redundant"),
        (7, "cofinal_successor_target"),
        (8, "global_outer_target"),
    )]
    if [row.get("id") for row in rows] != expected_ids:
        issues.append("row ids or ordering drifted")
    if [row.get("readiness") for row in rows] != [
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "open",
        "open",
    ]:
        issues.append("row readiness boundary drifted")

    summary = payload.get("summary", {})
    expected_summary = {
        "rows": 8,
        "compact_certified_boxes": 1900,
        "compact_rescout_required": False,
        "proved_cofinal_base": 208,
        "first_unresolved_linear_stage": 209,
        "q209_open_shell_regions": 2,
        "ready_rows": 6,
        "open_rows": 2,
    }
    if summary != expected_summary:
        issues.append("summary drifted")

    for key, record in payload.get("source_audit", {}).items():
        path = REPO_ROOT / record.get("path", "")
        if not path.exists():
            issues.append(f"source missing: {key}")
            continue
        if target.sha256_path(path) != record.get("sha256"):
            issues.append(f"source hash drifted: {key}")
        try:
            source = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            issues.append(f"source parse failure {key}: {exc}")
            continue
        if source.get("kind") != record.get("kind"):
            issues.append(f"source kind drifted: {key}")
        if source.get("status") != record.get("status"):
            issues.append(f"source status drifted: {key}")
    if set(payload.get("source_audit", {})) != set(target.SOURCES):
        issues.append("source audit roster drifted")

    for phrase in (
        "# Newman Outer-Frequency Frontier Consolidation Gate",
        "## Closed Compact Base",
        "1900 certified Arb/Taylor boxes, 0 unresolved boxes",
        "A new boundary scout on this rectangle would be weaker and redundant",
        "## Finite Cofinal Frontier",
        "first unresolved linear stage: Q_209",
        "## Live Outer Target",
        "0<t<=1/5 and x>38",
        "## Proof Boundary",
        "not a proof of RH",
    ):
        if phrase.lower() not in note.lower():
            issues.append(f"note phrase missing: {phrase}")
    if target.success_line(payload) not in note:
        issues.append("success line missing from note")
    return issues


def main() -> int:
    try:
        payload = json.loads(target.RESULT.read_text(encoding="utf-8"))
        note = target.NOTE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as exc:
        print(f"outer-frequency frontier gate: source failure: {exc}")
        return 1

    issues = validate(payload, note)
    rebuilt = target.build_payload()
    if payload != rebuilt:
        issues.append("stored result differs from deterministic rebuild")
    if issues:
        print(f"outer-frequency frontier gate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(target.success_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
