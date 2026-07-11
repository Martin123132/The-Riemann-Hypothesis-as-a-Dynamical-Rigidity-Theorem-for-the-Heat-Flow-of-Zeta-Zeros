#!/usr/bin/env python3
"""Validate the Jensen-window PF reciprocal fraction scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_reciprocal_fraction_scout.json"
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_fraction_sign_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_reciprocal_fraction_scout.md"

EXPECTED_FINITE_ROWS = 5 * 21 * 7

REQUIRED_SYMBOLIC_IDS = {
    "symbolic_standard_s_fraction_first_obstruction",
    "symbolic_standard_j_fraction_first_obstruction",
    "signed_or_modified_fraction_remaining_candidate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal Fraction Scout",
    "Status: reciprocal continued-fraction sign diagnostic",
    "This is not a proof of the all-order column recurrence",
    "a_2 = -g_2/g_1",
    "lambda_1 = -g_2",
    "work/rh_compute/results/jensen_window_pf_reciprocal_fraction_scout.json",
    "work/rh_compute/results/jensen_window_pf_reciprocal_fraction_sign_lamgrid_n0_n20_d2_d8_dps520.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_reciprocal_fraction_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_fraction_scout.py",
    "validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues",
    "rp_03_positive_stieltjes_or_j_fraction",
    "rp_04_companion_or_production_matrix_total_positivity",
    "standard positive Stieltjes/Jacobi",
    "signed or modified continued fraction",
)


@dataclass(frozen=True)
class FractionScoutIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> FractionScoutIssue:
    return FractionScoutIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def load_jsonl(path: Path) -> list[dict]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def validate_refs(summary: dict) -> list[FractionScoutIssue]:
    issues: list[FractionScoutIssue] = []
    grid = summary.get("finite_grid", {})
    for ref in grid.get("source_enclosures", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("finite_grid", "missing-source-enclosure", repr(ref)))
    row_log = grid.get("row_log")
    if not isinstance(row_log, str) or not (REPO_ROOT / row_log).exists():
        issues.append(issue("finite_grid", "missing-row-log", repr(row_log)))
    return issues


def validate_summary(summary: dict) -> list[FractionScoutIssue]:
    issues: list[FractionScoutIssue] = []
    if summary.get("kind") != "jensen_window_pf_reciprocal_fraction_scout":
        issues.append(issue("<summary>", "bad-kind", repr(summary.get("kind"))))
    boundary = str(summary.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<summary>", "weak-proof-boundary", str(summary.get("proof_boundary", ""))))
    targets = set(summary.get("target_route_rows", []))
    expected_targets = {
        "rp_03_positive_stieltjes_or_j_fraction",
        "rp_04_companion_or_production_matrix_total_positivity",
    }
    if targets != expected_targets:
        issues.append(issue("<summary>", "bad-target-route-rows", repr(sorted(targets))))

    symbolic = summary.get("symbolic_rows", [])
    if not isinstance(symbolic, list):
        issues.append(issue("symbolic_rows", "bad-symbolic-rows", repr(type(symbolic))))
        symbolic = []
    by_id = {row.get("id"): row for row in symbolic if isinstance(row, dict)}
    for row_id in sorted(REQUIRED_SYMBOLIC_IDS - set(by_id)):
        issues.append(issue(row_id, "missing-symbolic-row", row_id))
    s_row = by_id.get("symbolic_standard_s_fraction_first_obstruction", {})
    if s_row.get("first_parameter") != "g1":
        issues.append(issue("symbolic_s_fraction", "bad-first-parameter", repr(s_row.get("first_parameter"))))
    if s_row.get("second_parameter") != "-g2/g1":
        issues.append(issue("symbolic_s_fraction", "bad-second-parameter", repr(s_row.get("second_parameter"))))
    j_row = by_id.get("symbolic_standard_j_fraction_first_obstruction", {})
    if j_row.get("first_beta") != "g1":
        issues.append(issue("symbolic_j_fraction", "bad-first-beta", repr(j_row.get("first_beta"))))
    if j_row.get("first_lambda") != "-g2":
        issues.append(issue("symbolic_j_fraction", "bad-first-lambda", repr(j_row.get("first_lambda"))))
    signed_row = by_id.get("signed_or_modified_fraction_remaining_candidate", {})
    if "separate" not in str(signed_row.get("obstruction", "")).lower():
        issues.append(issue("symbolic_signed_fraction", "missing-theorem-warning", str(signed_row.get("obstruction", ""))))

    grid = summary.get("finite_grid", {})
    expected_grid = {
        "dps": 520,
        "rows": EXPECTED_FINITE_ROWS,
        "standard_positive_fraction_obstructed_rows": EXPECTED_FINITE_ROWS,
        "all_rows_obstruct_standard_positive_fraction": True,
    }
    for key, value in expected_grid.items():
        if grid.get(key) != value:
            issues.append(issue("finite_grid", f"bad-{key}", f"{grid.get(key)!r} != {value!r}"))
    if grid.get("degrees") != [2, 8]:
        issues.append(issue("finite_grid", "bad-degrees", repr(grid.get("degrees"))))
    if grid.get("shifts") != [0, 20]:
        issues.append(issue("finite_grid", "bad-shifts", repr(grid.get("shifts"))))
    if grid.get("lambdas") != ["0", "1e-6", "1e-4", "1e-2", "1e-1"]:
        issues.append(issue("finite_grid", "bad-lambdas", repr(grid.get("lambdas"))))

    summary_block = summary.get("summary", {})
    expected_summary = {
        "symbolic_rows": 3,
        "finite_rows": EXPECTED_FINITE_ROWS,
        "standard_positive_fraction_obstructed_rows": EXPECTED_FINITE_ROWS,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary_block.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary_block.get(key)!r} != {value!r}"))
    main_finding = str(summary_block.get("main_finding", "")).lower()
    for required in ("wrong", "a2=-g2/g1", "lambda1=-g2", "signed or modified"):
        if required not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_rows(rows: list[dict]) -> list[FractionScoutIssue]:
    issues: list[FractionScoutIssue] = []
    if len(rows) != EXPECTED_FINITE_ROWS:
        issues.append(issue("rows", "bad-row-count", str(len(rows))))
    seen_keys: set[tuple[str, int, int]] = set()
    for row in rows:
        row_id = f"{row.get('lam')}:{row.get('shift_n')}:{row.get('degree_d')}"
        key = (str(row.get("lam")), int(row.get("shift_n", -1)), int(row.get("degree_d", -1)))
        if key in seen_keys:
            issues.append(issue(row_id, "duplicate-row", repr(key)))
        seen_keys.add(key)
        if row.get("kind") != "jensen_window_pf_reciprocal_fraction_sign_row":
            issues.append(issue(row_id, "bad-kind", repr(row.get("kind"))))
        expected_classifications = {
            "h0_classification": "positive",
            "h1_classification": "positive",
            "h2_classification": "positive",
            "standard_s_a1_classification": "positive",
            "standard_s_a2_classification": "negative",
            "standard_j_beta0_classification": "positive",
            "standard_j_lambda1_classification": "negative",
        }
        for field, value in expected_classifications.items():
            if row.get(field) != value:
                issues.append(issue(row_id, f"bad-{field}", repr(row.get(field))))
        if row.get("standard_positive_fraction_obstructed") is not True:
            issues.append(issue(row_id, "not-obstructed", repr(row)))
    return issues


def validate_note(path: Path) -> list[FractionScoutIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[FractionScoutIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "continued fraction proof is complete", "jwpf_06 is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, rows_path: Path, note_path: Path) -> tuple[list[FractionScoutIssue], int, int]:
    summary = load_json(summary_path)
    rows = load_jsonl(rows_path)
    issues: list[FractionScoutIssue] = []
    issues.extend(validate_summary(summary))
    issues.extend(validate_refs(summary))
    issues.extend(validate_rows(rows))
    issues.extend(validate_note(note_path))
    symbolic_count = len(summary.get("symbolic_rows", [])) if isinstance(summary.get("symbolic_rows"), list) else 0
    return issues, symbolic_count, len(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, symbolic_count, finite_count = validate(args.summary, args.rows, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "symbolic_rows": symbolic_count,
                    "finite_rows": finite_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-RECIPROCAL-FRACTION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF reciprocal fraction scout: "
            f"{symbolic_count} symbolic rows, {finite_count} finite rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
