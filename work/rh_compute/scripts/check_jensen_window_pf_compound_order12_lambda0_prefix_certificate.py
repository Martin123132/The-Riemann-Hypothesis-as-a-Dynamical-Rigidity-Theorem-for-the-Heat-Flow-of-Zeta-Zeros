#!/usr/bin/env python3
"""Validate the lambda-zero order-twelve prefix certificate."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_lambda0_prefix_certificate as generator  # noqa: E402


def main() -> int:
    artifact = json.loads(generator.DEFAULT_OUT.read_text(encoding="utf-8"))
    issues = []
    if artifact.get("kind") != "jensen_window_pf_compound_order12_lambda0_prefix_certificate":
        issues.append("bad artifact kind")
    rebuilt = generator.build_artifact()
    if artifact != rebuilt:
        issues.append("rebuilt artifact differs")
    finite = artifact.get("finite", {})
    if finite.get("theorem") != "Q_(12,n)(0)>0 for every integer 0<=n<=3":
        issues.append("theorem changed")
    rows = finite.get("rows", [])
    if [row.get("n") for row in rows] != [0, 1, 2, 3]:
        issues.append("prefix rows changed")
    if any(row.get("classification") != "positive" for row in rows):
        issues.append("nonpositive row present")
    expected_summary = {
        "coefficient_rows": 26,
        "prefix_rows": 4,
        "positive_Q12_rows": 4,
        "inconclusive_rows": 0,
        "orders_above_12": 0,
        "rh_claims": 0,
    }
    if artifact.get("summary") != expected_summary:
        issues.append("summary changed")
    note = generator.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Fresh 520-digit outward-rounded determinant rebuilds",
        "Q_(12,n)(0)>0 for every integer 0<=n<=3",
        "not all-shift order twelve",
        "do not prove that ray",
    ):
        if marker not in note:
            issues.append(f"note misses {marker}")
    if issues:
        for issue in issues:
            print(f"issue: {issue}")
        return 1
    print("validated lambda-zero order-twelve prefix: 26 coefficients, 4 positive rows, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
