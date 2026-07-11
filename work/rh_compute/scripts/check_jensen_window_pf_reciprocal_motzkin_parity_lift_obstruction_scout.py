#!/usr/bin/env python3
"""Validate the Jensen-window PF reciprocal Motzkin parity-lift obstruction scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json"
)
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_lamgrid_n0_n20_d2_d8_m2_m8_dps520.jsonl"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md"

EXPECTED_ROWS = 5 * 21 * 7 * 7

REQUIRED_SYMBOLIC_IDS = {
    "same_length_mixed_sign_paths",
    "global_length_sign_obstruction",
    "closed_path_conjugation_invariance",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal Motzkin Parity-Lift Obstruction Scout",
    "Status: global parity/sign-lift obstruction diagnostic",
    "This is not a proof against state-space doubled",
    "H0^m",
    "U D H0^(m-2)",
    "same length",
    "global length-parity sign",
    "diagonal sign conjugation",
    "closed excursion",
    "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json",
    "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_lamgrid_n0_n20_d2_d8_m2_m8_dps520.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.py",
    "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",
    "rp_04_companion_or_production_matrix_total_positivity",
    "rp_09_signed_or_modified_continued_fraction",
    "finite evidence only",
)


@dataclass(frozen=True)
class ParityIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ParityIssue:
    return ParityIssue(section=section, issue=name, detail=detail)


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


def validate_refs(summary: dict) -> list[ParityIssue]:
    issues: list[ParityIssue] = []
    grid = summary.get("finite_grid", {})
    for ref in grid.get("source_enclosures", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("finite_grid", "missing-source-enclosure", repr(ref)))
    row_log = grid.get("row_log")
    if not isinstance(row_log, str) or not (REPO_ROOT / row_log).exists():
        issues.append(issue("finite_grid", "missing-row-log", repr(row_log)))
    return issues


def validate_summary(summary: dict) -> list[ParityIssue]:
    issues: list[ParityIssue] = []
    if summary.get("kind") != "jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout":
        issues.append(issue("<summary>", "bad-kind", repr(summary.get("kind"))))
    targets = set(summary.get("target_route_rows", []))
    for target in ("rp_04_companion_or_production_matrix_total_positivity", "rp_09_signed_or_modified_continued_fraction"):
        if target not in targets:
            issues.append(issue("<summary>", "missing-target-route-row", target))
    boundary = str(summary.get("proof_boundary", "")).lower()
    for required in ("not a proof", "state-space", "lambda <= 0"):
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
    mixed_text = " ".join(str(value) for value in by_id.get("same_length_mixed_sign_paths", {}).values())
    for required in ("H0^m", "U D H0^(m-2)", "positive", "negative"):
        if required not in mixed_text:
            issues.append(issue("same_length_mixed_sign_paths", "missing-symbolic-text", required))
    global_text = " ".join(str(value) for value in by_id.get("global_length_sign_obstruction", {}).values()).lower()
    for required in ("same length", "sign depending only on m", "negative"):
        if required not in global_text:
            issues.append(issue("global_length_sign_obstruction", "missing-symbolic-text", required))
    conjugation_text = " ".join(str(value) for value in by_id.get("closed_path_conjugation_invariance", {}).values()).lower()
    for required in ("diagonal", "closed", "invariant"):
        if required not in conjugation_text:
            issues.append(issue("closed_path_conjugation_invariance", "missing-symbolic-text", required))

    grid = summary.get("finite_grid", {})
    expected_grid = {
        "dps": 520,
        "parity_lift_obstruction_rows": EXPECTED_ROWS,
        "parity_lift_obstruction_ok_rows": EXPECTED_ROWS,
        "all_rows_have_same_length_mixed_sign_witnesses": True,
    }
    for key, value in expected_grid.items():
        if grid.get(key) != value:
            issues.append(issue("finite_grid", f"bad-{key}", f"{grid.get(key)!r} != {value!r}"))
    if grid.get("degrees") != [2, 8]:
        issues.append(issue("finite_grid", "bad-degrees", repr(grid.get("degrees"))))
    if grid.get("lengths") != [2, 8]:
        issues.append(issue("finite_grid", "bad-lengths", repr(grid.get("lengths"))))
    if grid.get("shifts") != [0, 20]:
        issues.append(issue("finite_grid", "bad-shifts", repr(grid.get("shifts"))))
    if grid.get("lambdas") != ["0", "1e-6", "1e-4", "1e-2", "1e-1"]:
        issues.append(issue("finite_grid", "bad-lambdas", repr(grid.get("lambdas"))))

    summary_block = summary.get("summary", {})
    expected_summary = {
        "symbolic_rows": 3,
        "parity_lift_obstruction_rows": EXPECTED_ROWS,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary_block.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary_block.get(key)!r} != {value!r}"))
    main_finding = str(summary_block.get("main_finding", "")).lower()
    for required in ("same-length", "global length-parity", "diagonal sign", "state space"):
        if required not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_rows(rows: list[dict]) -> list[ParityIssue]:
    issues: list[ParityIssue] = []
    if len(rows) != EXPECTED_ROWS:
        issues.append(issue("rows", "bad-row-count", str(len(rows))))
    seen: set[tuple] = set()
    for row in rows:
        if row.get("kind") != "jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_row":
            issues.append(issue("rows", "bad-kind", repr(row)))
            continue
        key = (
            str(row.get("kind")),
            str(row.get("lam")),
            int(row.get("shift_n", -1)),
            int(row.get("degree_d", -1)),
            int(row.get("path_length_m", -1)),
        )
        if key in seen:
            issues.append(issue(str(key), "duplicate-row", repr(key)))
        seen.add(key)
        expected = {
            "positive_path": "H0^m",
            "positive_path_classification": "positive",
            "positive_path_contains_zero": False,
            "negative_path": "U D H0^(m-2)",
            "negative_path_classification": "negative",
            "negative_path_contains_zero": False,
            "same_length_and_endpoints": True,
            "global_length_sign_can_fix": False,
            "diagonal_conjugation_changes_closed_path_sign": False,
            "ok": True,
        }
        for field, value in expected.items():
            if row.get(field) != value:
                issues.append(issue(str(key), f"bad-{field}", repr(row)))
        if int(row.get("path_length_m", -1)) < 2:
            issues.append(issue(str(key), "bad-path-length", repr(row)))
    return issues


def validate_note(path: Path) -> list[ParityIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ParityIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all parity lifts are impossible",
        "all modified models are impossible",
        "the signed route is dead",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, rows_path: Path, note_path: Path) -> tuple[list[ParityIssue], int, int]:
    summary = load_json(summary_path)
    rows = load_jsonl(rows_path)
    issues: list[ParityIssue] = []
    issues.extend(validate_summary(summary))
    issues.extend(validate_refs(summary))
    issues.extend(validate_rows(rows))
    issues.extend(validate_note(note_path))
    summary_block = summary.get("summary", {})
    return issues, int(summary_block.get("symbolic_rows", 0)), int(summary_block.get("parity_lift_obstruction_rows", 0))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, symbolic_count, row_count = validate(args.summary, args.rows, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "symbolic_rows": symbolic_count,
                    "parity_lift_obstruction_rows": row_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"PARITY {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: "
            f"{symbolic_count} symbolic rows, {row_count} mixed-sign witness rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
