#!/usr/bin/env python3
"""Validate the sparse degree-6 monotone-contraction column scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_sparse_degree6_scout as scout


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree6_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_sparse_degree6_scout.md"

EXPECTED_RESULT = (
    "validated Jensen-window PF monotone-contraction sparse degree-6 scout: "
    "10 degree-6 rows, 63347 Bernstein coefficients, m<=10, "
    "0 negative Bernstein rows, 0 zero Bernstein rows, 0 issues"
)
EXPECTED_COUNTS = {
    1: (1, [], "6"),
    2: (2, [1], "21"),
    3: (8, [3, 1], "56"),
    4: (56, [6, 3, 1], "126"),
    5: (616, [10, 6, 3, 1], "252"),
    6: (9856, [15, 10, 6, 3, 1], "462"),
    7: (9856, [15, 10, 6, 3, 1], "792"),
    8: (10472, [16, 10, 6, 3, 1], "1287"),
    9: (12768, [18, 11, 6, 3, 1], "2002"),
    10: (19712, [21, 13, 7, 3, 1], "14794/7"),
}


def load_payload(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_payload(payload: dict, expected_payload: dict) -> list[str]:
    issues: list[str] = []
    if payload != expected_payload:
        issues.append("stored payload differs from exact regenerated sparse certificate")
        return issues

    summary = payload.get("summary", {})
    if summary.get("degree6_rows") != 10:
        issues.append("degree6 row count is not 10")
    if summary.get("total_bernstein_coefficients") != 63347:
        issues.append("total Bernstein coefficient count is not 63347")
    if summary.get("negative_bernstein_rows") != 0:
        issues.append("negative Bernstein rows are present")
    if summary.get("zero_bernstein_rows") != 0:
        issues.append("zero Bernstein rows are present")
    if summary.get("ready_to_apply_rows") != 0 or summary.get("target_closing") is not False:
        issues.append("proof-boundary readiness flags drifted")

    rows = payload.get("rows", [])
    if len(rows) != 10:
        issues.append("row list length is not 10")
    for row in rows:
        size = row.get("minor_size")
        expected = EXPECTED_COUNTS.get(size)
        if expected is None:
            issues.append(f"unexpected row size {size!r}")
            continue
        count, multidegree, min_coeff = expected
        if row.get("degree") != 6:
            issues.append(f"{row.get('id')} has non-degree-6 row")
        if row.get("bernstein_coefficient_count") != count:
            issues.append(f"{row.get('id')} has wrong Bernstein count")
        if row.get("bernstein_multidegree") != multidegree:
            issues.append(f"{row.get('id')} has wrong Bernstein multidegree")
        if row.get("bernstein_min_coefficient") != min_coeff:
            issues.append(f"{row.get('id')} has wrong minimum Bernstein coefficient")
        if row.get("bernstein_negative_coefficient_count") != 0:
            issues.append(f"{row.get('id')} has negative Bernstein coefficients")
        if row.get("bernstein_zero_coefficient_count") != 0:
            issues.append(f"{row.get('id')} has zero Bernstein coefficients")
        if row.get("bernstein_coefficients_strictly_positive") is not True:
            issues.append(f"{row.get('id')} is not marked strictly positive")
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = [
        EXPECTED_RESULT,
        "Status: exact sparse degree-6 diagnostic. This is not a proof",
        "degree 6: 0 <= x1 <= x2 <= x3 <= x4 <= x5 <= 1",
        "mcs6_d6_m10: m=10",
        "min=14794/7",
        "not prove all column rows or all",
        "outputs/signed_hankel_jensen_dependency_graph.md",
    ]
    return [f"missing note text: {needle}" for needle in required if needle not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-path", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = load_payload(args.json_path)
    expected_payload = scout.build_payload()
    issues = validate_payload(payload, expected_payload)
    issues.extend(validate_note(args.note))
    for issue in issues:
        print(f"ISSUE {issue}")
    if not issues:
        print(EXPECTED_RESULT)
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
