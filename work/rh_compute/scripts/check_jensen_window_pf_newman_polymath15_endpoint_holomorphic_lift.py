#!/usr/bin/env python3
"""Validate the Polymath-15 endpoint holomorphic-lift artifact."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp

import jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift as lift


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.py"
)
EXPECTED_IDS = {
    "np15ehl_01_phase_cancellation",
    "np15ehl_02_slow_factor",
    "np15ehl_03_holomorphic_lift",
    "np15ehl_04_real_axis_match",
    "np15ehl_05_slow_derivative",
    "np15ehl_06_defect_ode",
    "np15ehl_07_saddle_gain",
    "np15ehl_08_C0_derivative_scout",
    "np15ehl_09_live_bound",
    "np15ehl_10_nonpromotion",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 10:
        issues.append("expected ten rows")
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("row ids drifted")
    if artifact.get("exact") != lift.build_exact():
        issues.append("exact lift payload drifted")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("expected one open strip handoff")
    open_row = next(
        (row for row in rows if row.get("id") == "np15ehl_09_live_bound"),
        {},
    )
    if open_row.get("readiness") != "not_ready_to_apply":
        issues.append("strip bound was promoted")

    diagnostics = artifact.get("diagnostics", {})
    try:
        phase_delta = mp.mpf(
            diagnostics["max_phase_cancellation_relative_delta"]
        )
        derivative_max = mp.mpf(diagnostics["C0_prime_grid_max"])
    except (KeyError, ValueError):
        issues.append("malformed diagnostics")
    else:
        if phase_delta >= mp.mpf("1e-60"):
            issues.append("phase cancellation diagnostic failed")
        if derivative_max >= mp.mpf("0.7"):
            issues.append("C0 derivative grid guard failed")
    if len(diagnostics.get("phase_rows", [])) != 3:
        issues.append("expected three phase diagnostics")
    if "diagnostic only" not in diagnostics.get("C0_prime_status", ""):
        issues.append("C0 derivative nonpromotion drifted")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "Exact Phase Cancellation",
        "Holomorphic Lift",
        "Defect Equation",
        "only a diagnostic",
        "not a `C1` remainder certificate",
        "`Lambda <= 0` or RH",
    )
    for marker in required:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    if "not a C1 remainder certificate" not in artifact.get("status", ""):
        issues.append("status boundary drifted")
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
        "validated Newman Polymath-15 endpoint holomorphic lift: "
        "10 rows, 0 issues, 4 exact identities, 1 defect equation, "
        "3 phase diagnostics, 1 open strip bound"
    )


if __name__ == "__main__":
    main()
