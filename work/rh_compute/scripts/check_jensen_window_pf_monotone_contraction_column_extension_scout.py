#!/usr/bin/env python3
"""Validate the monotone-contraction column extension scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_column_extension_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_column_extension_scout.md"


EXPECTED_D3_IDS = tuple(f"mccx_d3_m{size}" for size in range(1, 11))
EXPECTED_D4_IDS = tuple(f"mccx_d4_m{size}" for size in range(1, 8))
EXPECTED_D5_IDS = tuple(f"mccx_d5_m{size}" for size in range(1, 9))
EXPECTED_BEYOND_IDS = {"mccx_d3_m9", "mccx_d3_m10", "mccx_d4_m7"}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Monotone-Contraction Column Extension Scout",
    "Status: exact bounded column-extension diagnostic",
    "Artifact kind: `jensen_window_pf_monotone_contraction_column_extension_scout`.",
    "work/rh_compute/results/jensen_window_pf_monotone_contraction_column_extension_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_column_extension_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_column_extension_scout.py",
    "validated Jensen-window PF monotone-contraction column extension scout: 25 column rows, 3329 Bernstein coefficients, 3 beyond-frontier rows, 0 negative Bernstein rows, 0 issues",
    "degree 3: 0 <= x1 <= x2 <= 1",
    "degree 4: 0 <= x1 <= x2 <= x3 <= 1",
    "degree 5: 0 <= x1 <= x2 <= x3 <= x4 <= 1",
    "degree 3: m=9, m=10",
    "degree 4: m=7",
    "degree 5: m=1..8",
    "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md",
    "outputs/jensen_window_pf_column_recurrence_contract.md",
)

FORBIDDEN_TEXT = (
    "therefore rh",
    "we have proved lambda <= 0",
    "proves lambda <= 0",
    "jensen-window pf-infinity is proved",
    "all schur shapes are proved",
    "monotone contractions are proved for the zeta coefficients",
)


@dataclass(frozen=True)
class ScoutIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def require(condition: bool, issues: list[ScoutIssue], row_id: str, issue: str, detail: str) -> None:
    if not condition:
        issues.append(ScoutIssue(row_id, issue, detail))


def validate_row(row: dict, expected_degree: int) -> list[ScoutIssue]:
    issues: list[ScoutIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    require(row.get("degree") == expected_degree, issues, row_id, "bad-degree", repr(row.get("degree")))
    require(isinstance(row.get("minor_size"), int), issues, row_id, "bad-minor-size", repr(row.get("minor_size")))
    require(row.get("bernstein_coefficient_count", 0) > 0, issues, row_id, "bad-bernstein-count", repr(row.get("bernstein_coefficient_count")))
    require(row.get("bernstein_negative_coefficient_count") == 0, issues, row_id, "negative-bernstein", repr(row.get("bernstein_negative_coefficient_count")))
    require(row.get("bernstein_coefficients_positive") is True, issues, row_id, "not-positive-bernstein", repr(row.get("bernstein_coefficients_positive")))
    require(bool(str(row.get("monotone_region_polynomial", "")).strip()), issues, row_id, "missing-polynomial", row_id)
    require("Bounded exact column-recurrence certificate" in str(row.get("proof_boundary", "")), issues, row_id, "weak-proof-boundary", str(row.get("proof_boundary", "")))
    combined = " ".join(str(row.get(key, "")) for key in ("monotone_region_polynomial", "proof_boundary")).lower()
    for forbidden in FORBIDDEN_TEXT:
        if forbidden in combined:
            issues.append(ScoutIssue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_note(path: Path) -> list[ScoutIssue]:
    if not path.exists():
        return [ScoutIssue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ScoutIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        require(required in text, issues, "note", "missing-text", required)
    lowered = text.lower()
    for forbidden in FORBIDDEN_TEXT:
        if forbidden in lowered:
            issues.append(ScoutIssue("note", "forbidden-text", forbidden))
    return issues


def validate(result_path: Path, note_path: Path) -> tuple[list[ScoutIssue], dict]:
    result = load_json(result_path)
    issues: list[ScoutIssue] = []
    require(
        result.get("kind") == "jensen_window_pf_monotone_contraction_column_extension_scout",
        issues,
        "<result>",
        "bad-kind",
        repr(result.get("kind")),
    )
    require(
        result.get("status") == "exact_bounded_column_extension_diagnostic",
        issues,
        "<result>",
        "bad-status",
        repr(result.get("status")),
    )
    boundary = str(result.get("proof_boundary", "")).lower()
    for required in ("bounded", "does not prove", "lambda <= 0"):
        require(required in boundary, issues, "<result>", "weak-proof-boundary", result.get("proof_boundary", ""))

    rows_by_degree = result.get("rows_by_degree", {})
    d3_rows = rows_by_degree.get("3", [])
    d4_rows = rows_by_degree.get("4", [])
    d5_rows = rows_by_degree.get("5", [])
    require(tuple(row.get("id") for row in d3_rows) == EXPECTED_D3_IDS, issues, "rows", "bad-d3-ids", repr([row.get("id") for row in d3_rows]))
    require(tuple(row.get("id") for row in d4_rows) == EXPECTED_D4_IDS, issues, "rows", "bad-d4-ids", repr([row.get("id") for row in d4_rows]))
    require(tuple(row.get("id") for row in d5_rows) == EXPECTED_D5_IDS, issues, "rows", "bad-d5-ids", repr([row.get("id") for row in d5_rows]))
    all_rows = d3_rows + d4_rows + d5_rows
    for row in d3_rows:
        issues.extend(validate_row(row, 3))
    for row in d4_rows:
        issues.extend(validate_row(row, 4))
    for row in d5_rows:
        issues.extend(validate_row(row, 5))

    beyond_ids = {row.get("id") for row in all_rows if row.get("beyond_first_hard_frontier") is True}
    require(beyond_ids == EXPECTED_BEYOND_IDS, issues, "rows", "bad-beyond-frontier-ids", repr(sorted(beyond_ids)))

    expected_counts = {
        "degree3_rows": 10,
        "degree4_rows": 7,
        "degree5_rows": 8,
        "total_column_rows": 25,
        "total_bernstein_coefficients": 3329,
        "beyond_first_hard_frontier_rows": 3,
        "higher_degree_extension_rows": 8,
        "negative_bernstein_rows": 0,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    summary = result.get("summary", {})
    for key, value in expected_counts.items():
        require(summary.get(key) == value, issues, "summary", f"bad-{key}", repr(summary.get(key)))
    require(summary.get("total_column_rows") == len(all_rows), issues, "summary", "row-count-mismatch", repr(len(all_rows)))
    require(
        summary.get("total_bernstein_coefficients")
        == sum(row.get("bernstein_coefficient_count", 0) for row in all_rows),
        issues,
        "summary",
        "bernstein-count-mismatch",
        repr(summary.get("total_bernstein_coefficients")),
    )

    invariants = "\n".join(str(item) for item in result.get("invariants", []))
    for required in (
        "Every row has strictly positive Bernstein coefficients.",
        "The scout covers only degree 3 m<=10, degree 4 m<=7, and degree 5 m<=8.",
        "Rows beyond the first hard frontier are diagnostic extensions, not an all-m theorem.",
        "Degree-5 rows are higher-degree diagnostic extensions, not an all-degree theorem.",
        "No row is ready_to_apply.",
        "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
    ):
        require(required in invariants, issues, "invariants", "missing-invariant", required)

    for path in (
        result.get("source_frontier_scout"),
        result.get("source_column_recurrence_contract"),
        result.get("source_monotone_contraction_target"),
    ):
        require(isinstance(path, str) and (REPO_ROOT / path).exists(), issues, "paths", "missing-source", repr(path))
    issues.extend(validate_note(note_path))
    return issues, result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, result = validate(args.result, args.note)
    summary = result.get("summary", {})
    if args.json:
        print(json.dumps({"ok": not issues, "issues": [asdict(issue) for issue in issues]}, indent=2, sort_keys=True))
    else:
        for issue in issues:
            print(f"JENSEN-WINDOW-PF-MONOTONE-COLUMN-EXTENSION {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated Jensen-window PF monotone-contraction column extension scout: "
            f"{summary.get('total_column_rows', 0)} column rows, "
            f"{summary.get('total_bernstein_coefficients', 0)} Bernstein coefficients, "
            f"{summary.get('beyond_first_hard_frontier_rows', 0)} beyond-frontier rows, "
            f"{summary.get('negative_bernstein_rows', 0)} negative Bernstein rows, "
            f"{len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
