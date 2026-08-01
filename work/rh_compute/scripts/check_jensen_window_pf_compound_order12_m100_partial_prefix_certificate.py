#!/usr/bin/env python3
"""Validate the order-twelve partial endpoint-prefix certificate."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_m100_partial_prefix_certificate as generator  # noqa: E402


def main() -> int:
    artifact = json.loads(generator.DEFAULT_OUT.read_text(encoding="utf-8"))
    issues = []
    if artifact.get("kind") != "jensen_window_pf_compound_order12_m100_partial_prefix_certificate":
        issues.append("bad artifact kind")
    if artifact.get("source", {}).get("sha256") != generator.sha256(generator.SOURCE):
        issues.append("source hash changed")
    rebuilt = generator.build_artifact()
    if artifact != rebuilt:
        issues.append("rebuilt artifact differs")
    finite = artifact.get("finite", {})
    if finite.get("negative_indices") != [0, 1, 2, 3]:
        issues.append("negative set changed")
    if finite.get("positive_range") != [4, 1240]:
        issues.append("positive range changed")
    if finite.get("inconclusive_indices") != []:
        issues.append("inconclusive rows present")
    rows = finite.get("rows", [])
    if [row.get("n") for row in rows] != list(range(1241)):
        issues.append("rows are not contiguous")
    expected_summary = {
        "rows": 4,
        "ready_rows": 3,
        "open_rows": 1,
        "checked_shifts": 1241,
        "negative_Q12_rows": 4,
        "positive_Q12_rows": 1237,
        "inconclusive_Q12_rows": 0,
        "remaining_positive_shifts": 252,
        "all_shift_endpoint_theorems": 0,
        "rh_claims": 0,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("summary changed")
    note = generator.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "4 negative, 1237 positive",
        "Q_(12,n)(-100)<0 for n=0,1,2,3",
        "Q_(12,n)(-100)>0 for every 4<=n<=1240",
        "remaining 252 positive shifts",
        "must not be promoted",
    ):
        if marker not in note:
            issues.append(f"note misses {marker}")
    if issues:
        for issue in issues:
            print(f"issue: {issue}")
        return 1
    print("validated order-twelve partial endpoint prefix: 1241 shifts, 4 negative, 1237 positive, 0 inconclusive, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
