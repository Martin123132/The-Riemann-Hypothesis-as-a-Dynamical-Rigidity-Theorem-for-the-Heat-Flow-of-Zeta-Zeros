#!/usr/bin/env python3
"""Validate the published/candidate Newman upper-bound provenance audit."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_lambda01965_provenance_audit import (
    ARCHIVE_MD5,
    CANDIDATE_UPPER,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MANIFEST_BYTES,
    MANIFEST_ENTRIES,
    PORTABLE_RUNS,
    PT_HEIGHT,
    PUBLISHED_UPPER,
    build_exact,
)


EXPECTED_IDS = [
    f"np15lpa_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "published_lower"),
        (2, "published_upper"),
        (3, "candidate_source"),
        (4, "archive_integrity"),
        (5, "exact_row"),
        (6, "portable_rerun"),
        (7, "compiled_boundary"),
        (8, "review_boundary"),
        (9, "interval_effect"),
        (10, "handoff"),
    )
]


def validate(path: Path, note_path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_lambda01965_provenance_audit"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")

    exact = build_exact()
    if artifact.get("exact") != exact:
        issues.append("exact payload drifted")
    published = exact["published_state"]
    candidate = exact["candidate_record"]
    integrity = exact["archive_integrity_audit"]
    rerun = exact["portable_rerun_audit"]
    boundary = exact["unreproduced_boundary"]

    if published.get("platt_trudgian_height") != PT_HEIGHT:
        issues.append("published height drifted")
    if Fraction(published["platt_trudgian_upper"]["exact"]) != PUBLISHED_UPPER:
        issues.append("published upper bound drifted")
    if Fraction(candidate["claimed_upper"]["exact"]) != CANDIDATE_UPPER:
        issues.append("candidate upper bound drifted")
    if candidate.get("pt_height_margin", {}).get("exact") != "350479773/2":
        issues.append("candidate height margin drifted")
    if candidate.get("improvement_over_published_0p2", {}).get("exact") != "7/2000":
        issues.append("candidate improvement drifted")

    if integrity.get("archive_md5_observed") != ARCHIVE_MD5:
        issues.append("archive MD5 drifted")
    if integrity.get("manifest_entries") != MANIFEST_ENTRIES:
        issues.append("manifest entry count drifted")
    if integrity.get("manifest_bytes_hashed") != MANIFEST_BYTES:
        issues.append("manifest byte count drifted")
    if integrity.get("manifest_bad") != 0 or integrity.get("manifest_missing") != 0:
        issues.append("manifest integrity is not clean")

    if rerun.get("invocation_count") != len(PORTABLE_RUNS):
        issues.append("portable invocation count drifted")
    if rerun.get("pass_count") != len(PORTABLE_RUNS) or rerun.get("fail_count") != 0:
        issues.append("portable rerun is not clean")
    if any(row.get("exit_code") != 0 for row in rerun.get("invocations", [])):
        issues.append("nonzero portable exit code recorded")
    if len(boundary.get("compiled_packages_not_live_rerun", [])) != 3:
        issues.append("compiled-package boundary drifted")

    by_id = {row.get("id"): row for row in artifact.get("rows", [])}
    if by_id.get("np15lpa_07_compiled_boundary", {}).get("readiness") != (
        "guard_validated"
    ):
        issues.append("compiled nonpromotion gate missing")
    if by_id.get("np15lpa_08_review_boundary", {}).get("readiness") != (
        "guard_validated"
    ):
        issues.append("review nonpromotion gate missing")
    if by_id.get("np15lpa_10_handoff", {}).get("readiness") != (
        "not_ready_to_apply"
    ):
        issues.append("positive-boundary handoff was promoted")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "peer-reviewed interval",
        "June 2026 candidate",
        "does not claim a full FLINT/Arb rerun",
        "theorem-to-code correspondence",
        "positive Newman boundary",
        "Lambda<=0",
        "RH",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "Lambda 0.1965 Provenance Audit",
        "0 <= Lambda <= 1/5 = 0.2",
        "3,000,175,332,800".replace(",", ""),
        "393/2000=0.1965",
        "DRAFT v0.3 - not for distribution",
        "PRE-FREEZE DRAFT",
        ARCHIVE_MD5,
        "portable invocations=22",
        "portable pass=22",
        "portable fail=0",
        "reproduced unrefereed candidate",
        "[0,1/5] -> [0,393/2000]",
        "No `Lambda <= 0` or RH conclusion is supplied",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    issues = validate(args.artifact, args.note)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman Polymath-15 Lambda 0.1965 provenance audit: "
        f"{len(artifact['rows'])} rows, 0 issues, 1 published interval, "
        f"{len(PORTABLE_RUNS)} portable passes, 3 compiled-package boundaries, "
        "1 candidate nonpromotion gate"
    )


if __name__ == "__main__":
    main()
