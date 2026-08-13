#!/usr/bin/env python3
"""Validate the Newman C1 remainder supersession frontier gate."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_c1_remainder_supersession_frontier_gate as target  # noqa: E402


def validate(payload: dict, note: str) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != target.STEM:
        issues.append("artifact kind drifted")
    if payload.get("schema_version") != 1:
        issues.append("schema version drifted")
    status = payload.get("status", "")
    for phrase in (
        "C1 remainder supersession proved",
        "signed Xi theorem remains open",
    ):
        if phrase not in status:
            issues.append(f"status boundary missing: {phrase}")

    exact = payload.get("exact", {})
    if exact.get("closed_domain") != (
        "L>=50 and 0<tL<=25, including adjacent cutoffs"
    ):
        issues.append("closed C1 domain drifted")
    if exact.get("global_remainder") != {
        "eta_0": 100000,
        "eta_1": 200000,
        "value": "|r_[1]|<100000*exp(-5L/4)",
        "derivative": "|r_[1],x|<200000*L*exp(-5L/4)",
    }:
        issues.append("global C1 remainder budget drifted")
    if exact.get("boundary_transfer", {}).get("required_c2_residual_bounds") != 0:
        issues.append("C2 boundary-transfer count drifted")
    if exact.get("direct_remaining_domain") != "0<t<=1/5 and x>38":
        issues.append("outer domain drifted")

    rows = payload.get("rows", [])
    expected_ids = [
        "c1rsf_01_old_target_identified",
        "c1rsf_02_remainder_superseded",
        "c1rsf_03_global_c1_budget",
        "c1rsf_04_contact_box",
        "c1rsf_05_c1_boundary_transfer",
        "c1rsf_06_live_c1_frontier",
        "c1rsf_07_remote_pair_boundary",
        "c1rsf_08_ordinary_alternative",
        "c1rsf_09_outer_domain",
        "c1rsf_10_route_guard",
    ]
    if [row.get("id") for row in rows] != expected_ids:
        issues.append("row ids or ordering drifted")
    if [row.get("readiness") for row in rows] != [
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "ready_to_apply",
        "open",
        "open",
        "open",
        "ready_to_apply",
        "ready_to_apply",
    ]:
        issues.append("row readiness boundary drifted")

    expected_summary = {
        "rows": 10,
        "superseded_remainder_targets": 1,
        "cutoff_uniform_c1_remainders": 1,
        "eta_0": 100000,
        "eta_1": 200000,
        "required_c2_boundary_bounds": 0,
        "ready_rows": 7,
        "open_rows": 3,
    }
    if payload.get("summary") != expected_summary:
        issues.append("summary drifted")

    audit = payload.get("source_audit", {})
    if set(audit) != set(target.SOURCES):
        issues.append("source audit roster drifted")
    for key, record in audit.items():
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

    for phrase in (
        "# Newman C1 Remainder Supersession and Live Arithmetic Frontier",
        "## Superseded Target",
        "|r_[1]|<100000*exp(-5L/4)",
        "including adjacent cutoffs",
        "## Contact Box",
        "requires no C2 residual estimate",
        "## Live C1 Target",
        "carrier-plus-near",
        "## Independent Alternative",
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
        print(f"C1 remainder supersession frontier gate: source failure: {exc}")
        return 1

    issues = validate(payload, note)
    rebuilt = target.build_payload()
    if payload != rebuilt:
        issues.append("stored result differs from deterministic rebuild")
    if issues:
        print(f"C1 remainder supersession frontier gate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(target.success_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
