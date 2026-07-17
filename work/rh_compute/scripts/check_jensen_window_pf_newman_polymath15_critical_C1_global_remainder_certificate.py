#!/usr/bin/env python3
"""Validate the global corrected critical C1 remainder certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import mpmath as mp

import jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate as cert


REPO_ROOT = Path(__file__).resolve().parents[3]
ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.md"
GENERATOR = (
    REPO_ROOT
    / "work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.py"
)
EXPECTED_IDS = {
    "np15c1grc_01_transition_geometry",
    "np15c1grc_02_cubic_saddle",
    "np15c1grc_03_log_ratio",
    "np15c1grc_04_ratio_bound",
    "np15c1grc_05_endpoint_size",
    "np15c1grc_06_adjacent_bound",
    "np15c1grc_07_absorption",
    "np15c1grc_08_global_C1",
    "np15c1grc_09_live_arithmetic",
    "np15c1grc_10_nonpromotion",
}


def load_artifact(rebuild: bool) -> dict:
    if rebuild:
        subprocess.run([sys.executable, str(GENERATOR)], cwd=REPO_ROOT, check=True)
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def validate(artifact: dict) -> list[str]:
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate":
        issues.append("kind drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 10:
        issues.append("expected ten rows")
    if {row.get("id") for row in rows} != EXPECTED_IDS:
        issues.append("row ids drifted")
    if artifact.get("exact") != cert.build_exact():
        issues.append("exact transition payload drifted")
    interval = cert.interval_budget()
    if artifact.get("interval") != interval:
        issues.append("interval budget drifted")
    ratio_audit = artifact.get("numerical_log_ratio_audit", {})
    if ratio_audit != cert.log_ratio_audit():
        issues.append("numerical log-ratio audit drifted")
    if ratio_audit.get("row_count") != 33:
        issues.append("expected 33 complex log-ratio audit points")
    if mp.mpf(ratio_audit.get("max_relative_identity_delta", "inf")) >= mp.mpf("1e-80"):
        issues.append("direct log-ratio identity precision failed")
    if mp.mpf(ratio_audit.get("max_abs_log_ratio_times_abs_T", "inf")) >= 1:
        issues.append("sampled scaled log-ratio exceeded one")
    if interval.get("global_raw_constant_lt") != "2500":
        issues.append("global raw constant drifted")
    if interval.get("global_first_constant") != "5000":
        issues.append("global first constant drifted")
    if interval.get("global_norm_constant_lt") != "32000000":
        issues.append("global norm constant drifted")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("expected one open arithmetic handoff")
    live = next(
        (row for row in rows if row.get("id") == "np15c1grc_09_live_arithmetic"),
        {},
    )
    if live.get("readiness") != "not_ready_to_apply":
        issues.append("arithmetic target was promoted")

    note = NOTE.read_text(encoding="utf-8")
    required = (
        "Adjacent Ratio",
        "Transition Bound",
        "Global C1 Budget",
        "|log(main_n/endpoint_n)|<3/|T|",
        "including every cutoff transition",
        "direct complex-point audit",
        "only guards signs and branches",
        "That arithmetic lower bound is not proved here",
        "This is not a proof of",
        "`Lambda <= 0` or RH",
    )
    for marker in required:
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    status = artifact.get("status", "")
    if "including every cutoff" not in status or "not corrected transversality" not in status:
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
        "validated Newman Polymath-15 critical C1 global remainder certificate: "
        "10 rows, 0 issues, 2 exact transition identities, 1 global cutoff theorem, "
        "1 explicit C1 budget, 1 open arithmetic target"
    )


if __name__ == "__main__":
    main()
