#!/usr/bin/env python3
"""Validate the sparse degree-7 monotone-contraction frontier scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout as scout


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_sparse_degree7_frontier_scout.md"

EXPECTED_RESULT = (
    "validated Jensen-window PF monotone-contraction sparse degree-7 frontier scout: "
    "9 positive rows, 1 certificate-obstruction row, 932691 Bernstein coefficients, "
    "first obstruction m=10, 126 negative Bernstein coefficients, 0 zero Bernstein coefficients, 0 issues"
)
EXPECTED_COUNTS = {
    1: (1, [], "7", 0, "strict_positive_global_bernstein_certificate"),
    2: (2, [1], "28", 0, "strict_positive_global_bernstein_certificate"),
    3: (8, [3, 1], "84", 0, "strict_positive_global_bernstein_certificate"),
    4: (56, [6, 3, 1], "210", 0, "strict_positive_global_bernstein_certificate"),
    5: (616, [10, 6, 3, 1], "462", 0, "strict_positive_global_bernstein_certificate"),
    6: (9856, [15, 10, 6, 3, 1], "924", 0, "strict_positive_global_bernstein_certificate"),
    7: (216832, [21, 15, 10, 6, 3, 1], "1716", 0, "strict_positive_global_bernstein_certificate"),
    8: (216832, [21, 15, 10, 6, 3, 1], "3003", 0, "strict_positive_global_bernstein_certificate"),
    9: (226688, [22, 15, 10, 6, 3, 1], "72835/22", 0, "strict_positive_global_bernstein_certificate"),
    10: (261800, [24, 16, 10, 6, 3, 1], "-4928", 126, "global_bernstein_certificate_obstruction"),
}


def load_payload(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_payload(payload: dict, expected_payload: dict) -> list[str]:
    issues: list[str] = []
    if payload != expected_payload:
        issues.append("stored payload differs from exact regenerated sparse degree-7 frontier certificate")
        return issues

    summary = payload.get("summary", {})
    expected_summary = {
        "degree7_rows": 10,
        "positive_certificate_rows": 9,
        "certificate_obstruction_rows": 1,
        "first_certificate_obstruction_m": 10,
        "first_certificate_obstruction_min_coefficient": "-4928",
        "total_bernstein_coefficients": 932691,
        "positive_row_bernstein_coefficients": 670891,
        "obstruction_row_bernstein_coefficients": 261800,
        "negative_bernstein_coefficients": 126,
        "zero_bernstein_coefficients": 0,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(f"summary {key} expected {expected!r}, got {summary.get(key)!r}")
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
        count, multidegree, min_coeff, negative_count, classification = expected
        if row.get("degree") != 7:
            issues.append(f"{row.get('id')} has non-degree-7 row")
        if row.get("bernstein_coefficient_count") != count:
            issues.append(f"{row.get('id')} has wrong Bernstein count")
        if row.get("bernstein_multidegree") != multidegree:
            issues.append(f"{row.get('id')} has wrong Bernstein multidegree")
        if row.get("bernstein_min_coefficient") != min_coeff:
            issues.append(f"{row.get('id')} has wrong minimum Bernstein coefficient")
        if row.get("bernstein_negative_coefficient_count") != negative_count:
            issues.append(f"{row.get('id')} has wrong negative Bernstein count")
        if row.get("bernstein_zero_coefficient_count") != 0:
            issues.append(f"{row.get('id')} has zero Bernstein coefficients")
        if row.get("frontier_classification") != classification:
            issues.append(f"{row.get('id')} has wrong frontier classification")
    obstruction = [row for row in rows if row.get("minor_size") == 10]
    if obstruction:
        if obstruction[0].get("bernstein_min_index") != [22, 0, 10, 0, 0, 0]:
            issues.append("m=10 minimum Bernstein index drifted")
        examples = obstruction[0].get("bernstein_negative_examples", [])
        if not examples or examples[0] != {"index": [22, 0, 9, 0, 0, 0], "coefficient": "-36267/230"}:
            issues.append("m=10 first negative Bernstein example drifted")
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = [
        EXPECTED_RESULT,
        "Status: exact sparse degree-7 frontier diagnostic. This is not a proof",
        "degree 7: 0 <= x1 <= x2 <= x3 <= x4 <= x5 <= x6 <= 1",
        "mcs7_d7_m10: m=10",
        "min=-4928",
        "index=[22, 0, 9, 0, 0, 0], coefficient=-36267/230",
        "That failure is not a proof that the",
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
