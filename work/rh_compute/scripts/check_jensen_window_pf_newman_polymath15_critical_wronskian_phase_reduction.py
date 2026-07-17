#!/usr/bin/env python3
"""Validate the corrected critical Wronskian phase reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp

import jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction as reduction


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.py"
)
EXPECTED_IDS = {
    "np15cwpr_01_complex_split",
    "np15cwpr_02_cartesian_jet",
    "np15cwpr_03_polar_jet",
    "np15cwpr_04_crossing_wronskian",
    "np15cwpr_05_collision_cap",
    "np15cwpr_06_thin_collar",
    "np15cwpr_07_explicit_target",
    "np15cwpr_08_phase_critical_values",
    "np15cwpr_09_sign_diagnostic",
    "np15cwpr_10_nonpromotion",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 10:
        issues.append("expected ten rows")
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("row ids drifted")
    if artifact.get("exact") != reduction.build_exact():
        issues.append("exact Wronskian payload drifted")
    target_domains = reduction.build_exact().get("target_domains", "")
    if "c_*=4911678521/1933561194" not in target_domains:
        issues.append("live asymptotic domain drifted")
    if "bounded L" not in target_domains:
        issues.append("bounded-L compact obligation missing")
    fine = reduction.diagnostics(reduction.FINE_DPS)
    if artifact.get("diagnostics") != fine:
        issues.append("finite Wronskian diagnostics drifted")
    root_rows = fine.get("root_rows", [])
    if len(root_rows) != 2:
        issues.append("expected two corrected crossing diagnostics")
    elif not (
        mp.mpf(root_rows[0]["phase_velocity"]) < 0
        and mp.mpf(root_rows[1]["phase_velocity"]) > 0
    ):
        issues.append("crossing phase sign diagnostic drifted")
    if any(
        mp.mpf(row.get("abs_corrected_value_residual", "inf")) >= mp.mpf("1e-45")
        for row in root_rows
    ):
        issues.append("corrected crossing residual failed")
    lehmer_pair = fine.get("lehmer_pair", {})
    if len(lehmer_pair.get("roots", [])) != 2:
        issues.append("expected two corrected Lehmer crossings")
    if not mp.mpf(lehmer_pair.get("gap", "0")) > 0:
        issues.append("corrected Lehmer gap failed")
    open_rows = [row for row in rows if row.get("role") == "open_theorem_target"]
    if len(open_rows) != 1 or open_rows[0].get("readiness") != "not_ready_to_apply":
        issues.append("thin-collar target was promoted")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "Complex Main",
        "Robust Collision Gate",
        "Sign Shortcut Rejected",
        "Dynamical Bridge",
        "c_*=4911678521/1933561194",
        "denominator of the zero-flow generator",
        "phase-critical-value avoidance",
        "do not prove or disprove",
        "RH-level and is not established here",
    )
    for marker in required:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    status = artifact.get("status", "")
    if "open arithmetic separation" not in status or "not Lambda<=0 or RH" not in status:
        issues.append("status theorem boundary drifted")
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
        "validated Newman Polymath-15 critical Wronskian phase reduction: "
        "10 rows, 0 issues, 6 exact identities, 1 robust collision gate, "
        "2 precision-stable crossing diagnostics, 1 open thin-collar target"
    )


if __name__ == "__main__":
    main()
