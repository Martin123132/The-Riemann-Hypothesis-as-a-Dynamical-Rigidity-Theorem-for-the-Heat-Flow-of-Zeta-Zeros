#!/usr/bin/env python3
"""Validate the Jensen-window PF reciprocal signed Jacobi beta scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json"
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md"

EXPECTED_BETA_ROWS = 5 * 21 * sum(range(2, 9))
EXPECTED_NEGATIVE_ROWS = 5 * 21 * 6
EXPECTED_TERMINAL_ZERO_ROWS = 5 * 21
EXPECTED_POSITIVE_ROWS = EXPECTED_BETA_ROWS - EXPECTED_NEGATIVE_ROWS - EXPECTED_TERMINAL_ZERO_ROWS

REQUIRED_SYMBOLIC_IDS = {
    "last_column_shifted_beta_formula",
    "degree_two_formula_sanity_check",
    "finite_beta_signature",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal Signed Jacobi Beta Scout",
    "Status: signed Jacobi beta-parameter diagnostic",
    "This is not a proof of a signed continued-fraction theorem",
    "Delta_r^*",
    "Q_r = Delta_r^*/Delta_r",
    "beta_n = Q_{n+1} - Q_n",
    "not the fully shifted determinant",
    "beta_0 > 0",
    "beta_1 < 0",
    "beta_n > 0 for n >= 2",
    "terminal degree-2 beta_1 contains zero",
    "work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json",
    "work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_lamgrid_n0_n20_d2_d8_dps520.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_jacobi_beta_scout.py",
    "validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues",
    "rp_09_signed_or_modified_continued_fraction",
    "finite evidence only",
)


@dataclass(frozen=True)
class BetaIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> BetaIssue:
    return BetaIssue(section=section, issue=name, detail=detail)


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


def validate_refs(summary: dict) -> list[BetaIssue]:
    issues: list[BetaIssue] = []
    grid = summary.get("finite_grid", {})
    for ref in grid.get("source_enclosures", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("finite_grid", "missing-source-enclosure", repr(ref)))
    row_log = grid.get("row_log")
    if not isinstance(row_log, str) or not (REPO_ROOT / row_log).exists():
        issues.append(issue("finite_grid", "missing-row-log", repr(row_log)))
    return issues


def validate_summary(summary: dict) -> list[BetaIssue]:
    issues: list[BetaIssue] = []
    if summary.get("kind") != "jensen_window_pf_reciprocal_signed_jacobi_beta_scout":
        issues.append(issue("<summary>", "bad-kind", repr(summary.get("kind"))))
    if summary.get("target_route_row") != "rp_09_signed_or_modified_continued_fraction":
        issues.append(issue("<summary>", "bad-target-route-row", repr(summary.get("target_route_row"))))
    boundary = str(summary.get("proof_boundary", "")).lower()
    for required in ("not a proof", "production-matrix", "lambda <= 0"):
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
    formula_text = " ".join(str(value) for value in by_id.get("last_column_shifted_beta_formula", {}).values())
    for required in ("Delta_r^*", "Q_r=Delta_r^*/Delta_r", "beta_n=Q_{n+1}-Q_n", "fully shifted"):
        if required not in formula_text:
            issues.append(issue("last_column_shifted_beta_formula", "missing-formula-text", required))
    sanity_text = " ".join(str(value) for value in by_id.get("degree_two_formula_sanity_check", {}).values())
    for required in ("beta_0=g_1", "lambda_1=-g_2", "beta_1=0"):
        if required not in sanity_text:
            issues.append(issue("degree_two_formula_sanity_check", "missing-sanity-text", required))

    grid = summary.get("finite_grid", {})
    expected_grid = {
        "dps": 520,
        "beta_rows": EXPECTED_BETA_ROWS,
        "ok_beta_rows": EXPECTED_BETA_ROWS,
        "positive_beta_rows": EXPECTED_POSITIVE_ROWS,
        "negative_beta_rows": EXPECTED_NEGATIVE_ROWS,
        "terminal_zero_beta_rows": EXPECTED_TERMINAL_ZERO_ROWS,
        "expected_positive_beta_rows": EXPECTED_POSITIVE_ROWS,
        "expected_negative_beta_rows": EXPECTED_NEGATIVE_ROWS,
        "expected_terminal_zero_beta_rows": EXPECTED_TERMINAL_ZERO_ROWS,
        "all_beta_rows_match_signature": True,
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
        "beta_rows": EXPECTED_BETA_ROWS,
        "positive_beta_rows": EXPECTED_POSITIVE_ROWS,
        "negative_beta_rows": EXPECTED_NEGATIVE_ROWS,
        "terminal_zero_beta_rows": EXPECTED_TERMINAL_ZERO_ROWS,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary_block.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary_block.get(key)!r} != {value!r}"))
    main_finding = str(summary_block.get("main_finding", "")).lower()
    for required in ("beta_0", "beta_1", "finite", "missing all-order theorem"):
        if required not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def expected_signature(degree_d: int, beta_index: int) -> str:
    if degree_d == 2 and beta_index == 1:
        return "terminal_degree2_zero"
    if beta_index == 1:
        return "negative_second"
    return "positive"


def validate_rows(rows: list[dict]) -> list[BetaIssue]:
    issues: list[BetaIssue] = []
    if len(rows) != EXPECTED_BETA_ROWS:
        issues.append(issue("rows", "bad-row-count", str(len(rows))))
    counts = {"positive": 0, "negative": 0, "terminal": 0}
    seen: set[tuple] = set()
    for row in rows:
        if row.get("kind") != "jensen_window_pf_reciprocal_signed_jacobi_beta_row":
            issues.append(issue("rows", "bad-kind", repr(row)))
            continue
        degree_d = int(row.get("degree_d", -1))
        beta_index = int(row.get("beta_index", -1))
        key = (
            str(row.get("kind")),
            str(row.get("lam")),
            int(row.get("shift_n", -1)),
            degree_d,
            beta_index,
        )
        if key in seen:
            issues.append(issue(str(key), "duplicate-row", repr(key)))
        seen.add(key)
        expected = expected_signature(degree_d, beta_index)
        if row.get("expected_signature") != expected:
            issues.append(issue(str(key), "bad-expected-signature", repr(row)))
        classification = row.get("beta_classification")
        contains_zero = row.get("beta_contains_zero")
        if expected == "terminal_degree2_zero":
            counts["terminal"] += 1
            if classification not in {"zero", "inconclusive_contains_zero"}:
                issues.append(issue(str(key), "terminal-beta-not-zero-containing-class", repr(row)))
            if contains_zero is not True:
                issues.append(issue(str(key), "terminal-beta-does-not-contain-zero", repr(row)))
        elif expected == "negative_second":
            counts["negative"] += 1
            if classification != "negative":
                issues.append(issue(str(key), "second-beta-not-negative", repr(row)))
            if contains_zero is not False:
                issues.append(issue(str(key), "second-beta-contains-zero", repr(row)))
        else:
            counts["positive"] += 1
            if classification != "positive":
                issues.append(issue(str(key), "beta-not-positive", repr(row)))
            if contains_zero is not False:
                issues.append(issue(str(key), "positive-beta-contains-zero", repr(row)))
        if row.get("ok") is not True:
            issues.append(issue(str(key), "row-not-ok", repr(row)))
    if counts["positive"] != EXPECTED_POSITIVE_ROWS:
        issues.append(issue("rows", "bad-positive-count", str(counts["positive"])))
    if counts["negative"] != EXPECTED_NEGATIVE_ROWS:
        issues.append(issue("rows", "bad-negative-count", str(counts["negative"])))
    if counts["terminal"] != EXPECTED_TERMINAL_ZERO_ROWS:
        issues.append(issue("rows", "bad-terminal-count", str(counts["terminal"])))
    return issues


def validate_note(path: Path) -> list[BetaIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[BetaIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "signed jacobi proof is complete", "jwpf_06 is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, rows_path: Path, note_path: Path) -> tuple[list[BetaIssue], int, int, int, int, int]:
    summary = load_json(summary_path)
    rows = load_jsonl(rows_path)
    issues: list[BetaIssue] = []
    issues.extend(validate_summary(summary))
    issues.extend(validate_refs(summary))
    issues.extend(validate_rows(rows))
    issues.extend(validate_note(note_path))
    symbolic_count = len(summary.get("symbolic_rows", [])) if isinstance(summary.get("symbolic_rows"), list) else 0
    summary_block = summary.get("summary", {})
    return (
        issues,
        symbolic_count,
        int(summary_block.get("beta_rows", 0)),
        int(summary_block.get("positive_beta_rows", 0)),
        int(summary_block.get("negative_beta_rows", 0)),
        int(summary_block.get("terminal_zero_beta_rows", 0)),
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
    issues, symbolic_count, beta_count, positive_count, negative_count, terminal_count = validate(
        args.summary, args.rows, args.note
    )
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "symbolic_rows": symbolic_count,
                    "beta_rows": beta_count,
                    "positive_beta_rows": positive_count,
                    "negative_beta_rows": negative_count,
                    "terminal_zero_beta_rows": terminal_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"BETA {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF reciprocal signed Jacobi beta scout: "
            f"{symbolic_count} symbolic rows, {beta_count} beta rows, "
            f"{positive_count} positive rows, {negative_count} negative rows, "
            f"{terminal_count} terminal-zero rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
