#!/usr/bin/env python3
"""Validate the positive-boundary delta-localization gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_newman_positive_boundary_delta_localization_gate as gate


EXPECTED_IDS = [
    "npbdlg_01_boundary_attainment",
    "npbdlg_02_positive_time_equivalence",
    "npbdlg_03_delta_contradiction",
    "npbdlg_04_dominant_localization",
    "npbdlg_05_cofinal_sequence",
    "npbdlg_06_compact_minimum",
    "npbdlg_07_delta_target",
    "npbdlg_08_quadratic_calibration",
    "npbdlg_09_uniform_floor_guard",
    "npbdlg_10_proof_handoff",
]


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_positive_boundary_"
        "delta_localization_gate"
    ):
        issues.append("artifact kind mismatch")
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if artifact.get("exact") != gate.build_exact():
        issues.append("exact payload drifted")
    if artifact.get("source_audit") != gate.source_audit():
        issues.append("source audit drifted")
    if artifact.get("cofinal_examples") != gate.cofinal_examples():
        issues.append("cofinal examples drifted")

    expected_roles = {
        "published_input": 1,
        "exact_equivalence": 2,
        "exact_reduction": 2,
        "exact_composition": 1,
        "exact_countermodel": 1,
        "nonpromotion_gate": 1,
        "open_theorem_target": 1,
        "proof_guard": 1,
    }
    role_counts = {
        role: sum(row.get("role") == role for row in rows)
        for role in expected_roles
    }
    if role_counts != expected_roles:
        issues.append(f"role counts drifted: {role_counts}")

    open_row = next(
        (row for row in rows if row.get("role") == "open_theorem_target"),
        {},
    )
    if open_row.get("readiness") != "not_ready_to_apply":
        issues.append("delta-dependent target was promoted")
    guard_rows = [
        row
        for row in rows
        if row.get("role") in {"exact_countermodel", "nonpromotion_gate"}
    ]
    if len(guard_rows) != 2 or any(
        row.get("readiness") != "guard_validated" for row in guard_rows
    ):
        issues.append("endpoint-uniformity guard is not active")

    exact_text = json.dumps(artifact.get("exact", {}))
    for marker in (
        "R_delta=4*pi*exp(25/delta)",
        "delta_j=1/(5j)",
        "x^2-2t",
        "4t^2",
        "8t/A^2",
        "not logically necessary",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    status = artifact.get("status", "")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "delta-localized",
        "remains open",
        "not a proof",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "does not certify",
        "positive-time simplicity",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = gate.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Exact Localization",
        "Endpoint-Uniformity Guard",
        "Live Handoff",
        "4*pi*exp(125)",
        "not a proof of `Lambda <= 0` or RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=gate.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman positive-boundary delta-localization gate: "
        f"{len(artifact['rows'])} rows, 0 issues, 5 exact reductions, "
        "1 cofinal-strip criterion, 1 quadratic endpoint-uniformity "
        "countermodel, 1 delta-dependent open target"
    )


if __name__ == "__main__":
    main()
