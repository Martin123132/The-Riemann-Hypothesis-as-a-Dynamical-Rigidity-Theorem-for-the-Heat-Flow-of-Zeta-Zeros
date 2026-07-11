#!/usr/bin/env python3
"""Validate the Jensen-window PF reciprocal Motzkin path obstruction scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json"
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md"

EXPECTED_MU2_ROWS = 5 * 21 * 7
EXPECTED_BETA1_ROWS = 5 * 21 * 6
EXPECTED_TOTAL_ROWS = EXPECTED_MU2_ROWS + EXPECTED_BETA1_ROWS

REQUIRED_SYMBOLIC_IDS = {
    "ordinary_j_fraction_motzkin_model",
    "diagonal_conjugation_obstruction",
    "first_mu2_cancellation",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal Motzkin Path Obstruction Scout",
    "Status: ordinary Motzkin-path obstruction diagnostic",
    "This is not a proof against every modified signed continued-fraction",
    "mu_2 = beta_0^2 + lambda_1",
    "beta_0^2 > 0",
    "lambda_1 < 0",
    "mu_2 > 0",
    "beta_1 < 0",
    "diagonal sign conjugation",
    "raw ordinary Motzkin",
    "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json",
    "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.py",
    "validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",
    "rp_04_companion_or_production_matrix_total_positivity",
    "rp_09_signed_or_modified_continued_fraction",
    "finite evidence only",
)


@dataclass(frozen=True)
class MotzkinIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> MotzkinIssue:
    return MotzkinIssue(section=section, issue=name, detail=detail)


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


def validate_refs(summary: dict) -> list[MotzkinIssue]:
    issues: list[MotzkinIssue] = []
    grid = summary.get("finite_grid", {})
    for ref in grid.get("source_enclosures", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("finite_grid", "missing-source-enclosure", repr(ref)))
    row_log = grid.get("row_log")
    if not isinstance(row_log, str) or not (REPO_ROOT / row_log).exists():
        issues.append(issue("finite_grid", "missing-row-log", repr(row_log)))
    return issues


def validate_summary(summary: dict) -> list[MotzkinIssue]:
    issues: list[MotzkinIssue] = []
    if summary.get("kind") != "jensen_window_pf_reciprocal_motzkin_path_obstruction_scout":
        issues.append(issue("<summary>", "bad-kind", repr(summary.get("kind"))))
    targets = set(summary.get("target_route_rows", []))
    for target in ("rp_04_companion_or_production_matrix_total_positivity", "rp_09_signed_or_modified_continued_fraction"):
        if target not in targets:
            issues.append(issue("<summary>", "missing-target-route-row", target))
    boundary = str(summary.get("proof_boundary", "")).lower()
    for required in ("not a proof against every", "modified", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<summary>", "weak-proof-boundary", str(summary.get("proof_boundary", ""))))
            break

    symbolic = summary.get("symbolic_rows", [])
    if not isinstance(symbolic, list):
        issues.append(issue("symbolic_rows", "bad-symbolic-rows", repr(type(symbolic))))
        symbolic = []
    by_id = {row.get("id"): row for row in symbolic if isinstance(row, dict)}
    for row_id in sorted(REQUIRED_SYMBOLIC_IDS - set(by_id)):
        issues.append(issue(row_id, "missing-symbolic-row", row_id))
    diagonal_text = " ".join(str(value) for value in by_id.get("diagonal_conjugation_obstruction", {}).values()).lower()
    for required in ("diagonal", "beta", "lambda", "nonnegative tridiagonal"):
        if required not in diagonal_text:
            issues.append(issue("diagonal_conjugation_obstruction", "missing-symbolic-text", required))
    mu2_text = " ".join(str(value) for value in by_id.get("first_mu2_cancellation", {}).values())
    for required in ("mu_2=beta_0^2+lambda_1", "lambda_1<0", "mu_2>0"):
        if required not in mu2_text:
            issues.append(issue("first_mu2_cancellation", "missing-symbolic-text", required))

    grid = summary.get("finite_grid", {})
    expected_grid = {
        "dps": 520,
        "motzkin_mu2_cancellation_rows": EXPECTED_MU2_ROWS,
        "motzkin_mu2_cancellation_ok_rows": EXPECTED_MU2_ROWS,
        "beta1_diagonal_obstruction_rows": EXPECTED_BETA1_ROWS,
        "beta1_diagonal_obstruction_ok_rows": EXPECTED_BETA1_ROWS,
        "all_mu2_rows_show_negative_path_cancellation": True,
        "all_beta1_rows_block_nonnegative_diagonal": True,
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
        "motzkin_mu2_cancellation_rows": EXPECTED_MU2_ROWS,
        "beta1_diagonal_obstruction_rows": EXPECTED_BETA1_ROWS,
        "total_finite_rows": EXPECTED_TOTAL_ROWS,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary_block.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary_block.get(key)!r} != {value!r}"))
    main_finding = str(summary_block.get("main_finding", "")).lower()
    for required in ("motzkin", "lambda_1<0", "beta_1", "modified path"):
        if required not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_rows(rows: list[dict]) -> list[MotzkinIssue]:
    issues: list[MotzkinIssue] = []
    if len(rows) != EXPECTED_TOTAL_ROWS:
        issues.append(issue("rows", "bad-row-count", str(len(rows))))
    mu2_rows = 0
    beta1_rows = 0
    seen: set[tuple] = set()
    for row in rows:
        kind = row.get("kind")
        if kind == "jensen_window_pf_reciprocal_motzkin_mu2_cancellation_row":
            key = (
                kind,
                str(row.get("lam")),
                int(row.get("shift_n", -1)),
                int(row.get("degree_d", -1)),
            )
            mu2_rows += 1
            expected = {
                "beta0_classification": "positive",
                "beta0_contains_zero": False,
                "beta0_square_classification": "positive",
                "beta0_square_contains_zero": False,
                "lambda1_classification": "negative",
                "lambda1_contains_zero": False,
                "mu2_classification": "positive",
                "mu2_contains_zero": False,
                "has_negative_path_weight": True,
                "ok": True,
            }
            for field, value in expected.items():
                if row.get(field) != value:
                    issues.append(issue(str(key), f"bad-{field}", repr(row)))
        elif kind == "jensen_window_pf_reciprocal_motzkin_beta1_diagonal_obstruction_row":
            key = (
                kind,
                str(row.get("lam")),
                int(row.get("shift_n", -1)),
                int(row.get("degree_d", -1)),
            )
            beta1_rows += 1
            degree = int(row.get("degree_d", -1))
            if degree < 3:
                issues.append(issue(str(key), "unexpected-degree", repr(row)))
            expected = {
                "beta1_classification": "negative",
                "beta1_contains_zero": False,
                "diagonal_conjugation_changes_diagonal": False,
                "ok": True,
            }
            for field, value in expected.items():
                if row.get(field) != value:
                    issues.append(issue(str(key), f"bad-{field}", repr(row)))
        else:
            issues.append(issue("rows", "bad-kind", repr(row)))
            continue
        if key in seen:
            issues.append(issue(str(key), "duplicate-row", repr(key)))
        seen.add(key)
    if mu2_rows != EXPECTED_MU2_ROWS:
        issues.append(issue("rows", "bad-mu2-row-count", str(mu2_rows)))
    if beta1_rows != EXPECTED_BETA1_ROWS:
        issues.append(issue("rows", "bad-beta1-row-count", str(beta1_rows)))
    return issues


def validate_note(path: Path) -> list[MotzkinIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[MotzkinIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all production matrices are impossible",
        "the signed route is dead",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, rows_path: Path, note_path: Path) -> tuple[list[MotzkinIssue], int, int, int]:
    summary = load_json(summary_path)
    rows = load_jsonl(rows_path)
    issues: list[MotzkinIssue] = []
    issues.extend(validate_summary(summary))
    issues.extend(validate_refs(summary))
    issues.extend(validate_rows(rows))
    issues.extend(validate_note(note_path))
    summary_block = summary.get("summary", {})
    return (
        issues,
        int(summary_block.get("symbolic_rows", 0)),
        int(summary_block.get("motzkin_mu2_cancellation_rows", 0)),
        int(summary_block.get("beta1_diagonal_obstruction_rows", 0)),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, symbolic_count, mu2_count, beta1_count = validate(args.summary, args.rows, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "symbolic_rows": symbolic_count,
                    "motzkin_mu2_cancellation_rows": mu2_count,
                    "beta1_diagonal_obstruction_rows": beta1_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"MOTZKIN {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF reciprocal Motzkin path obstruction scout: "
            f"{symbolic_count} symbolic rows, {mu2_count} mu2 cancellation rows, "
            f"{beta1_count} beta1 diagonal obstruction rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
