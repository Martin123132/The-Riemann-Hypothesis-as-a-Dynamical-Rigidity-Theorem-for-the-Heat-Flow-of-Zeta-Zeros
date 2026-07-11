#!/usr/bin/env python3
"""Validate the Jensen-window PF column recurrence finite coverage map."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_COVERAGE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_column_recurrence_finite_coverage.md"

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Column Recurrence Finite Coverage",
    "Status: finite coverage map",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json",
    "python work/rh_compute/scripts/jensen_window_pf_column_recurrence_finite_coverage.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_finite_coverage.py",
    "validated Jensen-window PF column recurrence finite coverage: 1470 direct positive rows, 210 hard recurrence rows, 315 Sturm/PF windows, 0 issues",
    "d=3, m=1..8",
    "d=4, m=1..6",
    "d=3,m=8",
    "d=4,m=6",
    "315 checked Sturm/PF windows",
    "outputs/arb_jensen_window_column_recurrence_stress.md",
    "work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520_summary.json",
    "work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520.jsonl",
    "python work/rh_compute/scripts/check_arb_jensen_window_column_recurrence_stress.py",
    "validated 12600 Arb Jensen-window column recurrence stress rows with 0 issues",
    "finite evidence only",
    "not an all-order recurrence theorem",
)


@dataclass(frozen=True)
class CoverageIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> CoverageIssue:
    return CoverageIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_refs(coverage: dict) -> list[CoverageIssue]:
    issues: list[CoverageIssue] = []
    refs = [coverage.get("source_column_recurrence_contract")]
    direct = coverage.get("direct_arb_determinant_coverage", {})
    refs.extend([direct.get("source_summary"), direct.get("source_rows")])
    sturm = coverage.get("sturm_pf_window_coverage", {})
    refs.extend(sturm.get("source_summaries", []))
    stress = coverage.get("arb_column_recurrence_stress_coverage", {})
    refs.extend([stress.get("source_summary"), stress.get("source_rows")])
    for ref in refs:
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("refs", "missing-ref", repr(ref)))
    return issues


def validate_top_level(coverage: dict) -> list[CoverageIssue]:
    issues: list[CoverageIssue] = []
    if coverage.get("kind") != "jensen_window_pf_column_recurrence_finite_coverage":
        issues.append(issue("<coverage>", "bad-kind", repr(coverage.get("kind"))))
    if coverage.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<coverage>", "bad-target-obligation", repr(coverage.get("target_obligation"))))
    boundary = str(coverage.get("proof_boundary", "")).lower()
    if "finite coverage map only" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<coverage>", "weak-proof-boundary", str(coverage.get("proof_boundary", ""))))

    summary = coverage.get("summary", {})
    expected = {
        "direct_checked_rows": 1470,
        "direct_positive_rows": 1470,
        "direct_hard_checked_rows": 210,
        "direct_hard_ok_rows": 210,
        "sturm_pf_checked_windows": 315,
        "sturm_pf_ok_windows": 315,
        "stress_checked_rows": 12600,
        "stress_ok_rows": 12600,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("<summary>", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("1470 direct", "210 hard", "315 checked", "12600", "finite evidence only"):
        if required not in finding:
            issues.append(issue("<summary>", "missing-main-finding-text", required))
    return issues


def validate_direct(coverage: dict) -> list[CoverageIssue]:
    direct = coverage.get("direct_arb_determinant_coverage", {})
    issues: list[CoverageIssue] = []
    expected = {
        "kind": "arb_jensen_window_pf_obligation_summary",
        "checked_rows": 1470,
        "positive_rows": 1470,
        "failed_or_inconclusive_rows": 0,
        "all_ok": True,
        "hard_recurrence_rows_checked": 210,
        "hard_recurrence_rows_ok": 210,
    }
    for key, value in expected.items():
        if direct.get(key) != value:
            issues.append(issue("direct_arb_determinant_coverage", f"bad-{key}", f"{direct.get(key)!r} != {value!r}"))
    if direct.get("degree_sizes") != {"3": list(range(1, 9)), "4": list(range(1, 7))}:
        issues.append(issue("direct_arb_determinant_coverage", "bad-degree-sizes", repr(direct.get("degree_sizes"))))
    if direct.get("lambdas") != ["0", "1e-6", "1e-4", "1e-2", "1e-1"]:
        issues.append(issue("direct_arb_determinant_coverage", "bad-lambdas", repr(direct.get("lambdas"))))
    if direct.get("shifts") != list(range(21)):
        issues.append(issue("direct_arb_determinant_coverage", "bad-shifts", repr(direct.get("shifts"))))

    hard = direct.get("hard_recurrence_rows", {})
    for row_id in ("d3_column_recurrence_m8", "d4_column_recurrence_m6"):
        row = hard.get(row_id, {})
        if row.get("checked_rows") != 105 or row.get("positive_rows") != 105 or row.get("ok_rows") != 105:
            issues.append(issue(row_id, "bad-hard-row-counts", repr(row)))
        if row.get("failed_or_inconclusive_rows") != 0:
            issues.append(issue(row_id, "hard-row-not-all-ok", repr(row)))
    return issues


def validate_sturm(coverage: dict) -> list[CoverageIssue]:
    sturm = coverage.get("sturm_pf_window_coverage", {})
    issues: list[CoverageIssue] = []
    expected = {
        "degrees": [3, 4, 5],
        "checked_windows": 315,
        "ok_windows": 315,
        "failed_or_inconclusive_windows": 0,
        "rows_by_degree": {"3": 105, "4": 105, "5": 105},
    }
    for key, value in expected.items():
        if sturm.get(key) != value:
            issues.append(issue("sturm_pf_window_coverage", f"bad-{key}", f"{sturm.get(key)!r} != {value!r}"))
    if "including all column recurrence minors" not in str(sturm.get("finite_consequence", "")):
        issues.append(issue("sturm_pf_window_coverage", "missing-finite-consequence", str(sturm.get("finite_consequence", ""))))
    if "not an all-degree" not in str(sturm.get("proof_boundary", "")):
        issues.append(issue("sturm_pf_window_coverage", "weak-proof-boundary", str(sturm.get("proof_boundary", ""))))
    return issues


def validate_stress(coverage: dict) -> list[CoverageIssue]:
    stress = coverage.get("arb_column_recurrence_stress_coverage", {})
    issues: list[CoverageIssue] = []
    expected = {
        "degrees": list(range(3, 9)),
        "sizes": list(range(1, 21)),
        "checked_rows": 12600,
        "ok_rows": 12600,
        "failed_or_inconclusive_rows": 0,
        "ok_by_degree": {str(degree): 2100 for degree in range(3, 9)},
    }
    for key, value in expected.items():
        if stress.get(key) != value:
            issues.append(issue("arb_column_recurrence_stress_coverage", f"bad-{key}", f"{stress.get(key)!r} != {value!r}"))
    if stress.get("lambdas") != ["0", "1e-6", "1e-4", "1e-2", "1e-1"]:
        issues.append(issue("arb_column_recurrence_stress_coverage", "bad-lambdas", repr(stress.get("lambdas"))))
    if stress.get("shifts") != list(range(21)):
        issues.append(issue("arb_column_recurrence_stress_coverage", "bad-shifts", repr(stress.get("shifts"))))
    if "not an all-order recurrence theorem" not in str(stress.get("proof_boundary", "")):
        issues.append(issue("arb_column_recurrence_stress_coverage", "weak-proof-boundary", str(stress.get("proof_boundary", ""))))
    return issues


def validate_gap(coverage: dict) -> list[CoverageIssue]:
    gap = coverage.get("gap_analysis", {})
    issues: list[CoverageIssue] = []
    not_covered = " ".join(gap.get("not_covered", []))
    for required in ("all degrees d", "all shifts n", "all lambda values", "analytic proof"):
        if required not in not_covered:
            issues.append(issue("gap_analysis", "missing-not-covered", required))
    return issues


def validate_note(path: Path) -> list[CoverageIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[CoverageIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "recurrence theorem is proved", "all-order recurrence theorem is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(coverage_path: Path, note_path: Path) -> tuple[list[CoverageIssue], int, int, int]:
    coverage = load_json(coverage_path)
    issues: list[CoverageIssue] = []
    issues.extend(validate_top_level(coverage))
    issues.extend(validate_refs(coverage))
    issues.extend(validate_direct(coverage))
    issues.extend(validate_sturm(coverage))
    issues.extend(validate_stress(coverage))
    issues.extend(validate_gap(coverage))
    issues.extend(validate_note(note_path))
    summary = coverage.get("summary", {})
    return (
        issues,
        int(summary.get("direct_positive_rows", 0)),
        int(summary.get("direct_hard_ok_rows", 0)),
        int(summary.get("sturm_pf_ok_windows", 0)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coverage", type=Path, default=DEFAULT_COVERAGE)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, direct_rows, hard_rows, sturm_windows = validate(args.coverage, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "direct_positive_rows": direct_rows,
                    "direct_hard_ok_rows": hard_rows,
                    "sturm_pf_ok_windows": sturm_windows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-COLUMN-COVERAGE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF column recurrence finite coverage: "
            f"{direct_rows} direct positive rows, {hard_rows} hard recurrence rows, "
            f"{sturm_windows} Sturm/PF windows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
