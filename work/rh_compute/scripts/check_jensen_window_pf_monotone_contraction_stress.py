#!/usr/bin/env python3
"""Validate the Jensen-window PF monotone-contraction stress artifact."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_stress.md"

EXPECTED_DEGREES = list(range(3, 13))
EXPECTED_ROWS = sum(5 * (65 - degree) for degree in EXPECTED_DEGREES)
EXPECTED_OK_BY_DEGREE = {str(degree): 5 * (65 - degree) for degree in EXPECTED_DEGREES}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Monotone-Contraction Stress",
    "Status: finite Arb ratio-curvature stress",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_monotone_contraction_stress_summary`",
    "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json",
    "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_monotone_contraction_stress.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_stress.py",
    "validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues",
    "A_{k+2}*A_k**3 >= A_{k+1}**3*A_{k-1}",
    "degrees d=3..12",
    "k<=64",
    "2875 finite Arb-classified rows",
    "global min monotone gap sample",
    "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md",
)


@dataclass(frozen=True)
class StressIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> StressIssue:
    return StressIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[StressIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_summary(summary: dict) -> list[StressIssue]:
    issues: list[StressIssue] = []
    if summary.get("kind") != "jensen_window_pf_monotone_contraction_stress_summary":
        issues.append(issue("summary", "bad-kind", repr(summary.get("kind"))))
    if summary.get("status") != "finite_arb_ratio_curvature_stress":
        issues.append(issue("summary", "bad-status", repr(summary.get("status"))))
    if summary.get("degrees") != EXPECTED_DEGREES:
        issues.append(issue("summary", "bad-degrees", repr(summary.get("degrees"))))
    if summary.get("max_coefficient_index") != 64:
        issues.append(issue("summary", "bad-max-coefficient-index", repr(summary.get("max_coefficient_index"))))
    if summary.get("rows") != EXPECTED_ROWS:
        issues.append(issue("summary", "bad-row-count", repr(summary.get("rows"))))
    if summary.get("ok") != EXPECTED_ROWS:
        issues.append(issue("summary", "bad-ok-count", repr(summary.get("ok"))))
    if summary.get("failed_or_inconclusive") != 0:
        issues.append(issue("summary", "nonzero-failures", repr(summary.get("failed_or_inconclusive"))))
    if summary.get("all_ok") is not True:
        issues.append(issue("summary", "all-ok-not-true", repr(summary.get("all_ok"))))
    if summary.get("ok_by_degree") != EXPECTED_OK_BY_DEGREE:
        issues.append(issue("summary", "bad-ok-by-degree", repr(summary.get("ok_by_degree"))))
    condition = summary.get("condition", {})
    if condition.get("equivalent_polynomial_inequality") != "A_{k+2}*A_k**3 >= A_{k+1}**3*A_{k-1}":
        issues.append(issue("condition", "bad-equivalent-inequality", repr(condition.get("equivalent_polynomial_inequality"))))
    for key in ("row_jsonl", "source_frontier_scout"):
        issues.extend(validate_ref("summary", summary.get(key)))
    for ref in summary.get("enclosure_jsonl", []):
        issues.extend(validate_ref("summary", ref))
    boundary = str(summary.get("proof_boundary", "")).lower()
    for required in ("finite arb stress", "not all shifts", "not all lambda", "not an analytic", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("summary", "weak-proof-boundary", required))
    inner = summary.get("summary", {})
    expected_inner = {
        "stress_rows": EXPECTED_ROWS,
        "positive_rows": EXPECTED_ROWS,
        "failed_or_inconclusive_rows": 0,
        "degree_count": 10,
        "target_closing": False,
    }
    for key, value in expected_inner.items():
        if inner.get(key) != value:
            issues.append(issue("inner-summary", f"bad-{key}", f"{inner.get(key)!r} != {value!r}"))
    return issues


def validate_rows(summary: dict) -> tuple[list[StressIssue], int]:
    row_ref = summary.get("row_jsonl")
    if not isinstance(row_ref, str):
        return [issue("rows", "missing-row-jsonl", repr(row_ref))], 0
    row_path = REPO_ROOT / row_ref
    issues: list[StressIssue] = []
    count = 0
    ok_by_degree: Counter[str] = Counter()
    with row_path.open("r", encoding="utf-8") as handle:
        for line_no, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            count += 1
            row = json.loads(raw)
            section = f"row:{line_no}"
            if row.get("kind") != "jensen_window_pf_monotone_contraction_stress_row":
                issues.append(issue(section, "bad-kind", repr(row.get("kind"))))
            degree = row.get("degree_d")
            if degree not in EXPECTED_DEGREES:
                issues.append(issue(section, "bad-degree", repr(degree)))
            if row.get("contraction_count") != int(degree) - 1:
                issues.append(issue(section, "bad-contraction-count", repr(row.get("contraction_count"))))
            if row.get("max_coefficient_index") != int(row.get("shift_n")) + int(degree):
                issues.append(issue(section, "bad-max-index", repr(row)))
            for key in ("adjacent_log_concavity_ok", "monotone_contractions_ok", "ok"):
                if row.get(key) is not True:
                    issues.append(issue(section, f"{key}-not-true", repr(row.get(key))))
            if row.get("contains_zero") is not False:
                issues.append(issue(section, "contains-zero", repr(row.get("contains_zero"))))
            if degree in EXPECTED_DEGREES and row.get("ok") is True:
                ok_by_degree[str(degree)] += 1
    if count != EXPECTED_ROWS:
        issues.append(issue("rows", "bad-row-count", f"{count} != {EXPECTED_ROWS}"))
    if dict(sorted(ok_by_degree.items(), key=lambda item: int(item[0]))) != EXPECTED_OK_BY_DEGREE:
        issues.append(issue("rows", "bad-row-ok-by-degree", repr(ok_by_degree)))
    return issues, count


def validate_note(path: Path) -> list[StressIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[StressIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "monotone contractions are proved for all zeta windows",
        "analytic monotone-contraction theorem is proved",
        "cauchy-binet identity is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, note_path: Path) -> tuple[list[StressIssue], int, int]:
    summary = load_json(summary_path)
    issues: list[StressIssue] = []
    issues.extend(validate_summary(summary))
    row_issues, row_count = validate_rows(summary)
    issues.extend(row_issues)
    issues.extend(validate_note(note_path))
    return issues, row_count, int(summary.get("ok", 0))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, rows, positive_rows = validate(args.summary, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "rows": rows,
                    "positive_rows": positive_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-MONOTONE-STRESS {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF monotone contraction stress: "
            f"{rows} rows, {positive_rows} positive rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
