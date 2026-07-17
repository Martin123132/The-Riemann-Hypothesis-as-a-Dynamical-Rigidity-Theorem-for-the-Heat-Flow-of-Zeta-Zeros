#!/usr/bin/env python3
"""Validate the adjacent C0 recurrence and transition scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp

import jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate.py"
)
EXPECTED_IDS = {
    "np15c0rtg_01_recurrence",
    "np15c0rtg_02_signed_jump",
    "np15c0rtg_03_corrected_jump",
    "np15c0rtg_04_boundary_ratio",
    "np15c0rtg_05_transition_scout",
    "np15c0rtg_06_scale_guard",
    "np15c0rtg_07_live_bound",
    "np15c0rtg_08_absorption",
    "np15c0rtg_09_no_large_subtraction",
    "np15c0rtg_10_nonpromotion",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_endpoint_C0_recurrence_transition_gate":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 10:
        issues.append("expected ten rows")
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("row ids drifted")
    if artifact.get("exact") != gate.build_exact():
        issues.append("exact recurrence payload drifted")
    if sum(row.get("role") == "open_theorem_target" for row in rows) != 1:
        issues.append("expected one open complex-collar target")
    live = next(
        (row for row in rows if row.get("id") == "np15c0rtg_07_live_bound"),
        {},
    )
    if live.get("readiness") != "not_ready_to_apply":
        issues.append("complex-collar target was promoted")

    diagnostics = artifact.get("diagnostics", {})
    point_rows = diagnostics.get("rows", [])
    if diagnostics.get("row_count") != len(point_rows) or len(point_rows) < 10:
        issues.append("transition diagnostic count drifted")
    try:
        precision = mp.mpf(diagnostics["max_relative_delta"])
        max_value = mp.mpf(diagnostics["max_scaled_value"])
        max_first = mp.mpf(diagnostics["max_scaled_first"])
    except (KeyError, ValueError):
        issues.append("malformed transition diagnostics")
    else:
        if precision >= mp.mpf("1e-25"):
            issues.append("transition precision guard failed")
        if max_value >= 1 or max_first >= 1:
            issues.append("transition next-saddle guard failed")
    for row in point_rows:
        if mp.mpf(row["t"]) > mp.mpf("0.5"):
            issues.append(f"out-of-scope time at n={row.get('n')} c={row.get('c')}")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "Exact Recurrence",
        "C0(p+2)+C0(p)",
        "Boundary Scout",
        "Live Estimate",
        "remains finite point evidence",
        "This is not a proof of",
        "`Lambda <= 0` or RH",
    )
    for marker in required:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    status = artifact.get("status", "")
    if "not a complex-collar transition theorem" not in status or "RH" not in status:
        issues.append("status nonpromotion boundary drifted")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    artifact = load_artifact(args.rebuild)
    issues = validate(artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        raise SystemExit(1)
    print(
        "validated Newman Polymath-15 endpoint C0 recurrence transition gate: "
        f"10 rows, 0 issues, 1 exact recurrence, "
        f"{artifact['diagnostics']['row_count']} transition diagnostics, "
        "1 open complex-collar bound, 1 nonpromotion gate"
    )


if __name__ == "__main__":
    main()
