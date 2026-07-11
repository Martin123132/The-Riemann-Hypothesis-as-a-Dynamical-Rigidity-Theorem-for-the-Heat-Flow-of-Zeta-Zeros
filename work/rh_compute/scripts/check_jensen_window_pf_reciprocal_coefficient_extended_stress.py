#!/usr/bin/env python3
"""Validate the extended reciprocal coefficient stress artifact."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_stress.json"
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_lamgrid_n0_n32_d2_d12_m1_m40_dps520.jsonl"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md"

EXPECTED_LAMBDAS = ["0", "1e-6", "1e-4", "1e-2", "1e-1"]
EXPECTED_SHIFTS = list(range(33))
EXPECTED_DEGREES = list(range(2, 13))
EXPECTED_SIZES = list(range(1, 41))
EXPECTED_SUMMARY_ROWS = 5 * 33 * 11
EXPECTED_ROWS = EXPECTED_SUMMARY_ROWS * 40
EXPECTED_PER_DEGREE = 5 * 33 * 40

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal Coefficient Extended Stress",
    "Status: extended finite Arb reciprocal-coefficient diagnostic",
    "This is not a proof",
    "E(t) = 1 / H(-t)",
    "work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_stress.json",
    "work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_lamgrid_n0_n32_d2_d12_m1_m40_dps520.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_reciprocal_coefficient_extended_stress.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_coefficient_extended_stress.py",
    "validated Jensen-window PF reciprocal coefficient extended stress: 72600 rows, 0 issues",
    "degrees d=2..12",
    "sizes m=1..40",
    "shifts n=0..32",
    "lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}",
    "72,600 / 72,600",
    "finite evidence only",
    "not an all-order recurrence theorem",
)


@dataclass(frozen=True)
class ExtendedStressIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ExtendedStressIssue:
    return ExtendedStressIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_summary(summary: dict) -> list[ExtendedStressIssue]:
    issues: list[ExtendedStressIssue] = []
    if summary.get("kind") != "jensen_window_pf_reciprocal_coefficient_extended_stress":
        issues.append(issue("<summary>", "bad-kind", repr(summary.get("kind"))))
    if summary.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<summary>", "bad-target-obligation", repr(summary.get("target_obligation"))))
    if summary.get("lambdas") != EXPECTED_LAMBDAS:
        issues.append(issue("<summary>", "bad-lambdas", repr(summary.get("lambdas"))))
    if summary.get("shifts") != EXPECTED_SHIFTS:
        issues.append(issue("<summary>", "bad-shifts", repr(summary.get("shifts"))))
    if summary.get("degrees") != EXPECTED_DEGREES:
        issues.append(issue("<summary>", "bad-degrees", repr(summary.get("degrees"))))
    if summary.get("coefficient_sizes") != EXPECTED_SIZES:
        issues.append(issue("<summary>", "bad-sizes", repr(summary.get("coefficient_sizes"))))
    if summary.get("dps") != 520:
        issues.append(issue("<summary>", "bad-dps", repr(summary.get("dps"))))
    if summary.get("needed_max_k") != 44:
        issues.append(issue("<summary>", "bad-needed-max-k", repr(summary.get("needed_max_k"))))
    if summary.get("summary_rows") != EXPECTED_SUMMARY_ROWS:
        issues.append(issue("<summary>", "bad-summary-rows", repr(summary.get("summary_rows"))))
    if summary.get("rows") != EXPECTED_ROWS:
        issues.append(issue("<summary>", "bad-rows", repr(summary.get("rows"))))
    if summary.get("ok") != EXPECTED_ROWS or summary.get("failed_or_inconclusive") != 0:
        issues.append(issue("<summary>", "not-all-ok", f"ok={summary.get('ok')} failed={summary.get('failed_or_inconclusive')}"))
    if summary.get("all_ok") is not True:
        issues.append(issue("<summary>", "all-ok-false", repr(summary.get("all_ok"))))
    boundary = str(summary.get("proof_boundary", ""))
    if "not all degrees" not in boundary or "Lambda <= 0" not in boundary:
        issues.append(issue("<summary>", "weak-proof-boundary", boundary))
    row_log = summary.get("row_log")
    if not isinstance(row_log, str) or not (REPO_ROOT / row_log).exists():
        issues.append(issue("<summary>", "missing-row-log", repr(row_log)))
    for ref in summary.get("enclosure_jsonl", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("<summary>", "missing-enclosure-ref", repr(ref)))

    ok_by_degree = summary.get("ok_by_degree", {})
    for degree in EXPECTED_DEGREES:
        if ok_by_degree.get(str(degree)) != EXPECTED_PER_DEGREE:
            issues.append(issue("<summary>", f"bad-ok-by-degree-{degree}", repr(ok_by_degree.get(str(degree)))))

    summaries = summary.get("summaries", [])
    if len(summaries) != EXPECTED_SUMMARY_ROWS:
        issues.append(issue("<summary>", "bad-subsummary-count", repr(len(summaries))))
    for sub in summaries:
        row_id = f"lambda={sub.get('lam')} n={sub.get('shift_n')} d={sub.get('degree_d')}"
        if sub.get("all_ok") is not True:
            issues.append(issue(row_id, "subsummary-not-ok", repr(sub)))
        if sub.get("tests") != len(EXPECTED_SIZES) or sub.get("positive") != len(EXPECTED_SIZES):
            issues.append(issue(row_id, "bad-subsummary-counts", repr(sub)))
        if sub.get("negative") or sub.get("zero") or sub.get("inconclusive"):
            issues.append(issue(row_id, "subsummary-has-bad-classification", repr(sub)))

    block = summary.get("summary", {})
    if block.get("ready_to_apply_rows") != 0 or block.get("target_closing") is not False:
        issues.append(issue("summary", "bad-readiness", repr(block)))
    main_finding = str(block.get("main_finding", "")).lower()
    for required in ("72,600", "positive", "finite", "all-order theorem"):
        if required.lower() not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_rows(rows_path: Path) -> list[ExtendedStressIssue]:
    issues: list[ExtendedStressIssue] = []
    row_count = 0
    ok_count = 0
    per_degree = {degree: 0 for degree in EXPECTED_DEGREES}
    seen: set[tuple] = set()
    with rows_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row_count += 1
            row = json.loads(raw)
            row_id = f"lambda={row.get('lam')} n={row.get('shift_n')} d={row.get('degree_d')} m={row.get('coefficient_m')}"
            key = (str(row.get("lam")), row.get("shift_n"), row.get("degree_d"), row.get("coefficient_m"))
            if key in seen:
                issues.append(issue(row_id, "duplicate-row", repr(key)))
            seen.add(key)
            if row.get("kind") != "jensen_window_pf_reciprocal_coefficient_extended_stress_row":
                issues.append(issue(row_id, "bad-kind", repr(row.get("kind"))))
            degree = row.get("degree_d")
            if degree not in EXPECTED_DEGREES:
                issues.append(issue(row_id, "bad-degree", repr(degree)))
            elif row.get("ok") is True:
                per_degree[degree] += 1
            if row.get("shift_n") not in EXPECTED_SHIFTS:
                issues.append(issue(row_id, "bad-shift", repr(row.get("shift_n"))))
            if row.get("coefficient_m") not in EXPECTED_SIZES:
                issues.append(issue(row_id, "bad-size", repr(row.get("coefficient_m"))))
            expected_max_index = int(row.get("shift_n", -1000)) + int(row.get("degree_d", -1000))
            if row.get("max_coefficient_index") != expected_max_index:
                issues.append(issue(row_id, "bad-max-coefficient-index", repr(row)))
            if row.get("classification") != "positive" or row.get("contains_zero") is not False or row.get("ok") is not True:
                issues.append(issue(row_id, "row-not-positive", repr(row)))
            else:
                ok_count += 1
    if row_count != EXPECTED_ROWS:
        issues.append(issue("<rows>", "bad-row-count", repr(row_count)))
    if ok_count != EXPECTED_ROWS:
        issues.append(issue("<rows>", "bad-ok-count", repr(ok_count)))
    for degree, count in per_degree.items():
        if count != EXPECTED_PER_DEGREE:
            issues.append(issue("<rows>", f"bad-degree-count-{degree}", repr(count)))
    return issues


def validate_note(path: Path) -> list[ExtendedStressIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ExtendedStressIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "all-order recurrence theorem is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, rows_path: Path, note_path: Path) -> list[ExtendedStressIssue]:
    issues: list[ExtendedStressIssue] = []
    issues.extend(validate_summary(load_json(summary_path)))
    issues.extend(validate_rows(rows_path))
    issues.extend(validate_note(note_path))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues = validate(args.summary, args.rows, args.note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"RECIPROCAL-EXTENDED-STRESS {item.section} [{item.issue}] {item.detail}")
        print(f"validated Jensen-window PF reciprocal coefficient extended stress: {EXPECTED_ROWS} rows, {len(issues)} issues")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
