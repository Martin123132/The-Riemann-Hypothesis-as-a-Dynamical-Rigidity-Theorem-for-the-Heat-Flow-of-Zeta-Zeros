#!/usr/bin/env python3
"""Validate the corrected critical C1 fixed-cell remainder certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate as cert


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.py"
)
EXPECTED_IDS = {
    "np15c1crc_01_region",
    "np15c1crc_02_coefficient_sum",
    "np15c1crc_03_published_eAB",
    "np15c1crc_04_published_eC",
    "np15c1crc_05_holomorphic_lift",
    "np15c1crc_06_raw_collar",
    "np15c1crc_07_C1_transfer",
    "np15c1crc_08_transversality_budget",
    "np15c1crc_09_transition_handoff",
    "np15c1crc_10_nonpromotion",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 10:
        issues.append("expected ten rows")
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("row ids drifted")
    if artifact.get("exact") != cert.build_exact():
        issues.append("exact payload drifted")
    interval = cert.interval_budget()
    if artifact.get("interval") != interval:
        issues.append("interval budget drifted")
    if interval.get("combined_raw_constant_lt") != "2500":
        issues.append("raw collar constant drifted")
    if interval.get("first_constant") != "5000":
        issues.append("first-jet constant drifted")
    if interval.get("norm_squared_constant_lt") != "32000000":
        issues.append("error norm constant drifted")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("expected one transition handoff")
    transition = next(
        (row for row in rows if row.get("id") == "np15c1crc_09_transition_handoff"),
        {},
    )
    if transition.get("readiness") != "not_ready_to_apply":
        issues.append("transition gap was promoted")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "Critical C1 Cell Remainder Certificate",
        "e_A+e_B<1000*exp(-3L/4)",
        "e_C<100*exp(-3L/4)",
        "r^2+(r'/L)^2<32000000*exp(-3L/2)",
        "Neither is supplied by this certificate",
        "This is not a proof of",
        "`Lambda <= 0` or RH",
    )
    for marker in required:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    status = artifact.get("status", "")
    if "not corrected transversality" not in status or "RH" not in status:
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
        "validated Newman Polymath-15 critical C1 cell remainder certificate: "
        "10 rows, 0 issues, 5 interval constants, 1 explicit C1 budget, "
        "1 transition handoff, 1 nonpromotion gate"
    )


if __name__ == "__main__":
    main()
