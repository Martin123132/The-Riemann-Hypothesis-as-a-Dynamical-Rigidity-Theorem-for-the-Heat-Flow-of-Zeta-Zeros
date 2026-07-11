#!/usr/bin/env python3
"""Validate the defect complete-monotonicity scout and exact guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_defect_complete_monotonicity_scout as scout


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_defect_complete_monotonicity_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md"
EXPECTED = (
    "validated Jensen-window PF defect complete-monotonicity scout: "
    "3284 defect positives, 3288 log positives, 838 inconclusive, both certified through order 8, 5 lambdas, 1 exact all-shape countermodel, 0 issues"
)


def validate_payload(stored: dict, rebuilt: dict) -> list[str]:
    issues: list[str] = []
    if stored != rebuilt:
        return ["stored payload differs from independent Arb/exact reconstruction"]
    summary = stored.get("summary", {})
    expected_counts = {
        "lambda_count": 5,
        "order_rows": 65,
        "checked_intervals": 3705,
        "strictly_positive_intervals": 3284,
        "inconclusive_intervals": 421,
        "strictly_negative_intervals": 0,
        "fully_certified_orders": list(range(9)),
        "negative_log_strictly_positive_intervals": 3288,
        "negative_log_inconclusive_intervals": 417,
        "negative_log_strictly_negative_intervals": 0,
        "negative_log_fully_certified_orders": list(range(9)),
        "exact_countermodels": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_counts.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted: {summary.get(key)!r}")
    for row in stored.get("finite_rows", []):
        if row.get("strictly_negative_count") != 0:
            issues.append(f"finite row has a negative interval: {row.get('lam')} order {row.get('order')}")
        if row.get("order", 99) <= 8 and (
            row.get("all_strictly_positive") is not True or row.get("contains_zero_count") != 0
        ):
            issues.append(f"certified-order row failed: {row.get('lam')} order {row.get('order')}")
    for row in stored.get("negative_log_finite_rows", []):
        if row.get("strictly_negative_count") != 0:
            issues.append(f"log row has a negative interval: {row.get('lam')} order {row.get('order')}")
        if row.get("order", 99) <= 8 and (
            row.get("all_strictly_positive") is not True or row.get("contains_zero_count") != 0
        ):
            issues.append(f"certified-order log row failed: {row.get('lam')} order {row.get('order')}")
    guard = stored.get("exact_all_shape_countermodel", {})
    if guard.get("discriminant") != "-27/16" or guard.get("discriminant_is_negative") is not True:
        issues.append("exact cubic discriminant guard drifted")
    checks = guard.get("full_static_cone_checks", {})
    for key in ("tail_lower_walls_ok", "upper_walls_ok", "monotone_contractions_ok"):
        if checks.get(key) is not True:
            issues.append(f"countermodel cone check failed: {key}")
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = (
        EXPECTED,
        "Status: finite interval diagnostic with exact countermodel gate.",
        "(-1)^m*Delta^m d_k > 0, 0<=m<=8.",
        "421` intervals containing zero",
        "(-1)^m*Delta^m y_k > 0, 0<=m<=8.",
        "3288` log intervals",
        "417` contain zero",
        "(1/2)*delta_0",
        "1+3*z+(3/2)*z^2+(1/4)*z^3",
        "discriminant = -27/16 < 0.",
        "cannot be promoted directly to the all-shape bridge",
    )
    return [f"missing note text: {item}" for item in required if item not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    stored = json.loads(args.result.read_text(encoding="utf-8"))
    rebuilt = scout.build_payload()
    issues = validate_payload(stored, rebuilt)
    issues.extend(validate_note(args.note))
    for item in issues:
        print(f"ISSUE {item}")
    print(EXPECTED.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
