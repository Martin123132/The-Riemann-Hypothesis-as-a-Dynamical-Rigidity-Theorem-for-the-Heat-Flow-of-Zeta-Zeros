#!/usr/bin/env python3
"""Validate the Jensen-window PF ratio-condition scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SCOUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_ratio_condition_scout.md"

ALLOWED_STATUSES = {
    "rejected_by_countermodel",
    "rejected_by_constructed_extension",
    "tautological_window_condition",
    "open_candidate_not_validated",
}

REQUIRED_IDS = {
    "rc_01_adjacent_log_concavity",
    "rc_02_decreasing_ratio_contractions",
    "rc_03_second_order_log_concavity",
    "rc_04_selected_low_degree_bernstein_positivity",
    "rc_05_degree3_discriminant",
    "rc_06_degree4_discriminant",
    "rc_07_contraction_log_concavity",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Ratio-Condition Scout",
    "Status: ratio-condition theorem-search diagnostic",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_ratio_condition_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py",
    "validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction",
    "x1 = 65/66",
    "x2 = 44/65",
    "x3 = 50/77",
    "x3 = 1/3",
    "second-order log-concavity",
    "degree 3 size 8",
    "degree 4 size 6",
    "x2^2 >= x1*x3",
    "rejected_by_constructed_extension",
)


@dataclass(frozen=True)
class RatioIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def issue(row_id: str, name: str, detail: str) -> RatioIssue:
    return RatioIssue(row_id=row_id, issue=name, detail=detail)


def validate_refs(scout: dict) -> list[RatioIssue]:
    issues: list[RatioIssue] = []
    for key in ("source_algebra", "source_low_degree_scout", "source_frontier_scout", "source_contraction_log_concavity_scout"):
        ref = scout.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("<scout>", f"missing-{key}", repr(ref)))
    return issues


def validate_top_level(scout: dict) -> list[RatioIssue]:
    issues: list[RatioIssue] = []
    if scout.get("kind") != "jensen_window_pf_ratio_condition_scout":
        issues.append(issue("<scout>", "bad-kind", repr(scout.get("kind"))))
    if scout.get("target_ansatz") != "ansatz_02_positive_cauchy_binet_kernel":
        issues.append(issue("<scout>", "bad-target-ansatz", repr(scout.get("target_ansatz"))))
    boundary = str(scout.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<scout>", "weak-proof-boundary", str(scout.get("proof_boundary", ""))))

    summary = scout.get("summary", {})
    expected = {
        "candidate_rows": 7,
        "rejected_by_countermodel_rows": 4,
        "rejected_by_constructed_extension_rows": 1,
        "tautological_window_condition_rows": 2,
        "open_candidate_not_validated_rows": 0,
        "target_closing_rows": 0,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("<summary>", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    main_finding = str(summary.get("main_finding", "")).lower()
    if "insufficient bridge conditions" not in main_finding or "contraction log-concavity" not in main_finding:
        issues.append(issue("<summary>", "missing-main-finding", str(summary.get("main_finding", ""))))
    return issues


def validate_countermodel_point(scout: dict) -> list[RatioIssue]:
    issues: list[RatioIssue] = []
    point = scout.get("countermodel_ratio_point", {})
    expected = {"x1": "65/66", "x2": "44/65", "x3": "50/77"}
    for key, value in expected.items():
        if point.get(key) != value:
            issues.append(issue("countermodel_ratio_point", f"bad-{key}", f"{point.get(key)!r} != {value!r}"))
    values = point.get("condition_values", {})
    positive_keys = [
        "1 - x1",
        "1 - x2",
        "1 - x3",
        "x1 - x2",
        "x2 - x3",
        "(1 - x2)^2 - x2^2*(1 - x1)*(1 - x3)",
    ]
    for key in positive_keys:
        if key not in values or values[key].startswith("-"):
            issues.append(issue("countermodel_ratio_point", "expected-positive-condition-value", f"{key}={values.get(key)!r}"))
    if values.get("x2^2 - x1*x3") != "-1946249/10735725":
        issues.append(issue("countermodel_ratio_point", "bad-contraction-log-concavity-value", repr(values.get("x2^2 - x1*x3"))))
    return issues


def validate_rows(scout: dict) -> list[RatioIssue]:
    rows = scout.get("candidate_rows", [])
    issues: list[RatioIssue] = []
    if not isinstance(rows, list):
        return [issue("candidate_rows", "bad-rows", repr(type(rows)))]
    seen: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("candidate_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        seen.add(row_id)
        status = row.get("status")
        if status not in ALLOWED_STATUSES:
            issues.append(issue(row_id, "bad-status", repr(status)))
        if row.get("target_closing") is not False:
            issues.append(issue(row_id, "target-closing-overclaim", repr(row.get("target_closing"))))
        for key in ("condition", "conclusion", "next_action", "proof_boundary"):
            if not str(row.get(key, "")).strip():
                issues.append(issue(row_id, f"missing-{key}", row_id))
        boundary = str(row.get("proof_boundary", "")).lower()
        if "not a proof" not in boundary:
            issues.append(issue(row_id, "weak-row-boundary", str(row.get("proof_boundary", ""))))
        if status == "rejected_by_countermodel":
            if row.get("countermodel_satisfies") is not True:
                issues.append(issue(row_id, "rejected-row-not-satisfied", repr(row.get("countermodel_satisfies"))))
            failures = row.get("target_failures_at_countermodel", [])
            failure_ids = {failure.get("id") for failure in failures if isinstance(failure, dict)}
            if {"d3_contiguous_m8", "d4_contiguous_m6"} - failure_ids:
                issues.append(issue(row_id, "missing-target-failures", repr(sorted(failure_ids))))
        if status == "rejected_by_constructed_extension":
            if row.get("countermodel_satisfies") is not False:
                issues.append(issue(row_id, "constructed-row-should-fail-original-countermodel", repr(row.get("countermodel_satisfies"))))
            if row.get("constructed_extension_satisfies") is not True:
                issues.append(issue(row_id, "constructed-extension-not-satisfied", repr(row.get("constructed_extension_satisfies"))))
            constructed_values = row.get("condition_values_at_constructed_extension", {})
            expected_values = {
                "x3": "1/3",
                "x2^2 - x1*x3": "108703/836550",
            }
            for key, value in expected_values.items():
                if constructed_values.get(key) != value:
                    issues.append(issue(row_id, f"bad-constructed-{key}", repr(constructed_values.get(key))))
            failures = row.get("target_failures_at_constructed_extension", [])
            failure_ids = {failure.get("id") for failure in failures if isinstance(failure, dict)}
            if {"d3_contiguous_m8", "d4_contiguous_m6"} - failure_ids:
                issues.append(issue(row_id, "missing-constructed-target-failures", repr(sorted(failure_ids))))
            for failure in failures:
                if isinstance(failure, dict) and failure.get("sign") != -1:
                    issues.append(issue(row_id, "constructed-target-failure-not-negative", repr(failure)))
            if "not a standalone bridge" not in str(row.get("conclusion", "")).lower():
                issues.append(issue(row_id, "missing-constructed-rejection-warning", str(row.get("conclusion", ""))))
        if status == "tautological_window_condition":
            if row.get("countermodel_satisfies") is not False:
                issues.append(issue(row_id, "tautological-row-should-fail-countermodel", repr(row.get("countermodel_satisfies"))))
            if "target itself" not in str(row.get("conclusion", "")).lower() and "window hyperbolicity" not in str(row.get("conclusion", "")).lower():
                issues.append(issue(row_id, "missing-tautology-warning", str(row.get("conclusion", ""))))
        if status == "open_candidate_not_validated":
            if row.get("countermodel_satisfies") is not False:
                issues.append(issue(row_id, "open-row-should-not-be-countermodel-satisfied", repr(row.get("countermodel_satisfies"))))
            if "not reject" not in str(row.get("conclusion", "")).lower():
                issues.append(issue(row_id, "missing-open-warning", str(row.get("conclusion", ""))))
    missing = REQUIRED_IDS - seen
    for row_id in sorted(missing):
        issues.append(issue(row_id, "missing-required-row", row_id))
    return issues


def validate_note(path: Path) -> list[RatioIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[RatioIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "kernel identity is proved", "cauchy-binet identity is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(scout_path: Path, note_path: Path) -> tuple[list[RatioIssue], int, int, int]:
    scout = load_json(scout_path)
    issues: list[RatioIssue] = []
    issues.extend(validate_top_level(scout))
    issues.extend(validate_refs(scout))
    issues.extend(validate_countermodel_point(scout))
    issues.extend(validate_rows(scout))
    issues.extend(validate_note(note_path))
    summary = scout.get("summary", {})
    return (
        issues,
        int(summary.get("candidate_rows", 0)),
        int(summary.get("rejected_by_countermodel_rows", 0)),
        int(summary.get("rejected_by_constructed_extension_rows", 0)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scout", type=Path, default=DEFAULT_SCOUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, row_count, rejected_count, construction_count = validate(args.scout, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "candidate_rows": row_count,
                    "rejected_by_countermodel_rows": rejected_count,
                    "rejected_by_constructed_extension_rows": construction_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-RATIO-CONDITION {item.row_id} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF ratio-condition scout: "
            f"{row_count} candidate rows, {len(issues)} issues, "
            f"{rejected_count} rejected by countermodel, "
            f"{construction_count} rejected by construction"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
