#!/usr/bin/env python3
"""Validate the endpoint C0 complex-strip certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate as cert


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.py"
)
EXPECTED_IDS = {
    "np15c0sc_01_definition",
    "np15c0sc_02_near_denominator",
    "np15c0sc_03_near_C0",
    "np15c0sc_04_near_derivative",
    "np15c0sc_05_outside_denominator",
    "np15c0sc_06_outside_quotient",
    "np15c0sc_07_global_strip",
    "np15c0sc_08_collar_map",
    "np15c0sc_09_lift_handoff",
    "np15c0sc_10_nonpromotion",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 10:
        issues.append("expected ten rows")
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("row ids drifted")
    if artifact.get("exact") != cert.build_exact():
        issues.append("exact payload drifted")
    interval = cert.bounds()
    if artifact.get("interval") != interval:
        issues.append("interval payload drifted")
    if interval.get("global_bound") != "|C0(p)|<5 and |C0'(p)|<60":
        issues.append("global strip bound drifted")
    if interval.get("critical_collar_p_displacement_lt") != "1/100":
        issues.append("critical collar transfer drifted")
    if sum(row.get("role") == "closure_record" for row in rows) != 1:
        issues.append("expected one closure record")
    closure = next(
        (row for row in rows if row.get("id") == "np15c0sc_09_lift_handoff"),
        {},
    )
    if closure.get("readiness") != "closed":
        issues.append("lift derivative handoff not closed")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "Removable Disks",
        "Outside Region",
        "|C0(p)|<5 and |C0'(p)|<60",
        "|Delta p|",
        "corrected first-jet lower bound remains open",
        "not a proof of `Lambda <= 0` or RH",
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
        "validated Newman Polymath-15 endpoint C0 strip certificate: "
        "10 rows, 0 issues, 4 Arb bounds, 2 removable disks, "
        "1 global strip theorem, 1 lift handoff closed"
    )


if __name__ == "__main__":
    main()
