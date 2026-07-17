#!/usr/bin/env python3
"""Validate the high-precision multiplier complete-monotonicity frontier."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_multiplier_complete_monotonicity_frontier_scout as scout  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
EXPECTED = (
    "validated Jensen-window PF multiplier complete-monotonicity frontier scout: "
    "7980 positive intervals, 0 inconclusive, orders 0..55, 5 lambdas, 0 issues"
)


def validate(stored: dict, note_path: Path) -> list[str]:
    issues: list[str] = []
    rebuilt = scout.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from independent 250-dps Arb reconstruction")
    if stored.get("kind") != "jensen_window_pf_multiplier_complete_monotonicity_frontier_scout":
        issues.append(f"bad kind: {stored.get('kind')!r}")
    summary = stored.get("summary", {})
    expected_summary = {
        "lambda_count": 5,
        "coefficient_rows": 290,
        "contraction_rows": 280,
        "difference_order_rows": 280,
        "checked_intervals": 7980,
        "strictly_positive_intervals": 7980,
        "inconclusive_intervals": 0,
        "strictly_negative_intervals": 0,
        "fully_certified_orders": list(range(56)),
        "max_certified_order": 55,
        "target_closing": False,
        "ready_to_apply_rows": 0,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted: {summary.get(key)!r}")
    weak = summary.get("weakest_certified_row", {})
    if (weak.get("lam"), weak.get("order"), weak.get("k")) != ("0.0", 29, 27):
        issues.append(f"weakest row drifted: {weak!r}")
    rows = stored.get("rows", [])
    if len(rows) != 280:
        issues.append(f"bad row count: {len(rows)}")
    for row in rows:
        if row.get("all_strictly_positive") is not True:
            issues.append(f"nonpositive row: {row.get('lam')} order {row.get('order')}")
        if row.get("contains_zero_count") != 0 or row.get("strictly_negative_count") != 0:
            issues.append(f"bad sign counts: {row!r}")
    for ref in stored.get("source_enclosures", []):
        if not (REPO_ROOT / ref).exists():
            issues.append(f"missing source: {ref}")
    if not note_path.exists():
        issues.append(f"missing note: {note_path}")
    else:
        text = note_path.read_text(encoding="utf-8")
        required = (
            EXPECTED,
            "Status: finite high-precision interval frontier diagnostic.",
            "(-1)^m*Delta^m y_k>0,",
            "0<=m<=55, 1<=k<=56-m,",
            "lambda=0.0, order=29, k=27",
            "insufficient working precision",
            "does not construct the required unit-atomic",
        )
        issues.extend(
            f"missing note text: {item}" for item in required if item not in text
        )
    boundary = str(stored.get("proof_boundary", "")).lower()
    for marker in (
        "finite alternating difference",
        "does not prove all-k",
        "counting measure",
        "pf-infinity",
        "rh",
        "lambda <= 0",
    ):
        if marker not in boundary:
            issues.append(f"weak proof boundary: {marker}")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=scout.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=scout.DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    stored = json.loads(args.result.read_text(encoding="utf-8"))
    issues = validate(stored, args.note)
    for item in issues:
        print(f"ISSUE {item}")
    print(EXPECTED.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
