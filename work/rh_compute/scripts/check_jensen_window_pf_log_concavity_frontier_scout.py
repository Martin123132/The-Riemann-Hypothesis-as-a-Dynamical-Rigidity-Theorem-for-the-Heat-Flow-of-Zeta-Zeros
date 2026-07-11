#!/usr/bin/env python3
"""Validate the Jensen-window PF log-concavity frontier scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json"
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_log_concavity_frontier_scout.md"


REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Log-Concavity Frontier Scout",
    "Status: symbolic frontier diagnostic",
    "This is not a proof of a positive kernel",
    "work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_log_concavity_frontier_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py",
    "validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues",
    "degree 3 size 6",
    "degree 3 size 8",
    "degree 4 size 5",
    "degree 4 size 6",
    "kernel_identity_found=false",
    "Ratio-Condition Audit",
    "outputs/jensen_window_pf_ratio_condition_scout.md",
    "work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py",
    "validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction",
)


@dataclass(frozen=True)
class FrontierIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def issue(row_id: str, name: str, detail: str) -> FrontierIssue:
    return FrontierIssue(row_id=row_id, issue=name, detail=detail)


def validate_top_level(scout: dict) -> list[FrontierIssue]:
    issues: list[FrontierIssue] = []
    if scout.get("kind") != "jensen_window_pf_log_concavity_frontier_scout":
        issues.append(issue("<scout>", "bad-kind", repr(scout.get("kind"))))
    if scout.get("target_ansatz") != "ansatz_02_positive_cauchy_binet_kernel":
        issues.append(issue("<scout>", "bad-target-ansatz", repr(scout.get("target_ansatz"))))
    boundary = str(scout.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<scout>", "weak-proof-boundary", str(scout.get("proof_boundary", ""))))
    summary = scout.get("summary", {})
    if summary.get("degree3_rows") != 8:
        issues.append(issue("<summary>", "bad-degree3-row-count", repr(summary.get("degree3_rows"))))
    if summary.get("degree4_rows") != 6:
        issues.append(issue("<summary>", "bad-degree4-row-count", repr(summary.get("degree4_rows"))))
    if summary.get("total_rows") != 14:
        issues.append(issue("<summary>", "bad-total-row-count", repr(summary.get("total_rows"))))
    if summary.get("kernel_identity_found") is not False:
        issues.append(issue("<summary>", "kernel-overclaim", repr(summary.get("kernel_identity_found"))))
    if summary.get("target_closing") is not False:
        issues.append(issue("<summary>", "target-closing-overclaim", repr(summary.get("target_closing"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("degree 3 size 6", "degree 3 size 8", "degree 4 size 5", "degree 4 size 6"):
        if required not in finding:
            issues.append(issue("<summary>", "missing-frontier-finding", required))
    return issues


def validate_frontiers(scout: dict, algebra: dict) -> list[FrontierIssue]:
    issues: list[FrontierIssue] = []
    frontiers = scout.get("frontiers", {})
    expected = {
        "3": {
            "first_negative_bernstein_size": 6,
            "first_countermodel_negative_size": 8,
            "first_countermodel_negative_value": algebra["finite_countermodel"]["d3_first_negative_contiguous_toeplitz_minor"]["determinant"],
        },
        "4": {
            "first_negative_bernstein_size": 5,
            "first_countermodel_negative_size": 6,
            "first_countermodel_negative_value": algebra["finite_countermodel"]["d4_first_negative_contiguous_toeplitz_minor"]["determinant"],
        },
    }
    for degree, checks in expected.items():
        row = frontiers.get(degree, {})
        for key, value in checks.items():
            if row.get(key) != value:
                issues.append(issue(f"degree{degree}", f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
        if row.get("first_negative_bernstein_min_coefficient") in (None, "0"):
            issues.append(issue(f"degree{degree}", "missing-negative-bernstein-min", repr(row.get("first_negative_bernstein_min_coefficient"))))
    return issues


def validate_rows(scout: dict) -> list[FrontierIssue]:
    issues: list[FrontierIssue] = []
    rows_by_degree = scout.get("rows_by_degree", {})
    expected_sizes = {"3": 8, "4": 6}
    expected_bernstein_first = {"3": 6, "4": 5}
    expected_counter_first = {"3": 8, "4": 6}
    for degree, count in expected_sizes.items():
        rows = rows_by_degree.get(degree, [])
        if len(rows) != count:
            issues.append(issue(f"degree{degree}", "bad-row-count", str(len(rows))))
            continue
        for row in rows:
            row_id = str(row.get("id", "<missing-id>"))
            size = row.get("minor_size")
            for key in (
                "determinant",
                "normalized_ratio_polynomial",
                "bernstein_min_coefficient",
                "countermodel_value",
                "proof_boundary",
            ):
                if key not in row:
                    issues.append(issue(row_id, "missing-field", key))
            if size < expected_bernstein_first[degree] and row.get("bernstein_coefficients_nonnegative") is not True:
                issues.append(issue(row_id, "early-bernstein-failure", repr(row.get("bernstein_min_coefficient"))))
            if size == expected_bernstein_first[degree] and row.get("bernstein_coefficients_nonnegative") is not False:
                issues.append(issue(row_id, "missing-first-bernstein-failure", repr(row.get("bernstein_min_coefficient"))))
            if size < expected_counter_first[degree] and row.get("countermodel_positive") is not True:
                issues.append(issue(row_id, "early-countermodel-negative", repr(row.get("countermodel_value"))))
            if size == expected_counter_first[degree]:
                if row.get("countermodel_positive") is not False or row.get("countermodel_sign") != -1:
                    issues.append(issue(row_id, "missing-first-countermodel-negative", repr(row.get("countermodel_value"))))
            boundary = str(row.get("proof_boundary", "")).lower()
            if "not a positive kernel identity" not in boundary:
                issues.append(issue(row_id, "weak-row-boundary", str(row.get("proof_boundary", ""))))
    return issues


def validate_algebra_check(scout: dict, algebra: dict) -> list[FrontierIssue]:
    issues: list[FrontierIssue] = []
    row = scout.get("algebra_countermodel_check", {})
    expected = {
        "d3_expected_first_negative_size": algebra["finite_countermodel"]["d3_first_negative_contiguous_toeplitz_minor"]["size"],
        "d3_expected_first_negative_value": algebra["finite_countermodel"]["d3_first_negative_contiguous_toeplitz_minor"]["determinant"],
        "d4_expected_first_negative_size": algebra["finite_countermodel"]["d4_first_negative_contiguous_toeplitz_minor"]["size"],
        "d4_expected_first_negative_value": algebra["finite_countermodel"]["d4_first_negative_contiguous_toeplitz_minor"]["determinant"],
    }
    for key, value in expected.items():
        if row.get(key) != value:
            issues.append(issue("algebra_countermodel_check", f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
    interpretation = str(scout.get("interpretation", "")).lower()
    if "strictly weaker than jensen-window pf-infinity" not in interpretation:
        issues.append(issue("<scout>", "missing-interpretation", str(scout.get("interpretation", ""))))
    return issues


def validate_refs(scout: dict) -> list[FrontierIssue]:
    source = scout.get("source_algebra")
    if not isinstance(source, str) or not (REPO_ROOT / source).exists():
        return [issue("<scout>", "missing-source-algebra", repr(source))]
    return []


def validate_note(path: Path) -> list[FrontierIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[FrontierIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "kernel identity is proved", "cauchy-binet identity is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(scout_path: Path, algebra_path: Path, note_path: Path) -> tuple[list[FrontierIssue], int]:
    scout = load_json(scout_path)
    algebra = load_json(algebra_path)
    issues: list[FrontierIssue] = []
    issues.extend(validate_top_level(scout))
    issues.extend(validate_refs(scout))
    issues.extend(validate_frontiers(scout, algebra))
    issues.extend(validate_rows(scout))
    issues.extend(validate_algebra_check(scout, algebra))
    issues.extend(validate_note(note_path))
    return issues, scout.get("summary", {}).get("total_rows", 0)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scout", type=Path, default=DEFAULT_SCOUT)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, row_count = validate(args.scout, args.algebra, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "contiguous_rows": row_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-LC-FRONTIER {item.row_id} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF log-concavity frontier scout: "
            f"{row_count} contiguous rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
