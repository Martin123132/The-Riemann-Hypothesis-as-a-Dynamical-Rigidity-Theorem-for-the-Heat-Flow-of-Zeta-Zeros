#!/usr/bin/env python3
"""Validate the Jensen-window PF state-space sign-lift obstruction scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_scout.json"
)
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md"

EXPECTED_ROWS = 5 * 21 * 7
REQUIRED_SYMBOLIC_IDS = {
    "absolute_value_sign_state_cover",
    "mu2_absolute_lift_gap",
    "surviving_state_space_requirement",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF State-Space Sign-Lift Obstruction Scout",
    "Status: state-space sign-lift obstruction diagnostic",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_scout.json",
    "work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_state_space_sign_lift_obstruction_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_state_space_sign_lift_obstruction_scout.py",
    "validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues",
    "mu_2 = beta_0^2 + lambda_1",
    "kappa_1 = -lambda_1 > 0",
    "absolute_lift_mu2 - mu_2 = 2*kappa_1 = -2*lambda_1 > 0",
    "735 / 735",
    "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
    "msm_04_state_space_doubled_model",
    "rp_04_companion_or_production_matrix_total_positivity",
    "rp_09_signed_or_modified_continued_fraction",
    "does not rule out a genuinely modified state-space doubled model",
)


@dataclass(frozen=True)
class SignLiftIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> SignLiftIssue:
    return SignLiftIssue(section=section, issue=name, detail=detail)


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


def validate_ref(section: str, ref: str) -> list[SignLiftIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_summary(summary: dict) -> list[SignLiftIssue]:
    issues: list[SignLiftIssue] = []
    if summary.get("kind") != "jensen_window_pf_state_space_sign_lift_obstruction_scout":
        issues.append(issue("<summary>", "bad-kind", repr(summary.get("kind"))))
    if summary.get("target_model_row") != "msm_04_state_space_doubled_model":
        issues.append(issue("<summary>", "bad-target-model-row", repr(summary.get("target_model_row"))))
    targets = set(summary.get("target_route_rows", []))
    for target in ("rp_04_companion_or_production_matrix_total_positivity", "rp_09_signed_or_modified_continued_fraction"):
        if target not in targets:
            issues.append(issue("<summary>", "missing-target-route-row", target))
    boundary = str(summary.get("proof_boundary", "")).lower()
    for required in ("rejects only", "absolute-value", "not every state-space", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<summary>", "weak-proof-boundary", required))
    for ref in summary.get("source_artifacts", []):
        if isinstance(ref, str):
            issues.extend(validate_ref("source_artifacts", ref))
        else:
            issues.append(issue("source_artifacts", "bad-ref", repr(ref)))

    symbolic = summary.get("symbolic_rows", [])
    if not isinstance(symbolic, list):
        issues.append(issue("symbolic_rows", "bad-symbolic-rows", repr(type(symbolic))))
        symbolic = []
    by_id = {row.get("id"): row for row in symbolic if isinstance(row, dict)}
    for row_id in sorted(REQUIRED_SYMBOLIC_IDS - set(by_id)):
        issues.append(issue(row_id, "missing-symbolic-row", row_id))
    gap_text = " ".join(str(value) for value in by_id.get("mu2_absolute_lift_gap", {}).values()).lower()
    for required in ("beta_0^2-kappa_1", "beta_0^2+kappa_1", "2*kappa_1", "lambda_1"):
        if required not in gap_text:
            issues.append(issue("mu2_absolute_lift_gap", "missing-symbolic-text", required))
    requirement_text = " ".join(str(value) for value in by_id.get("surviving_state_space_requirement", {}).values()).lower()
    for required in ("exactly e(t)=1/h(-t)", "manifestly nonnegative", "not a construction"):
        if required not in requirement_text:
            issues.append(issue("surviving_state_space_requirement", "missing-symbolic-text", required))

    grid = summary.get("finite_grid", {})
    expected_grid = {
        "mu2_sign_lift_obstruction_rows": EXPECTED_ROWS,
        "mu2_sign_lift_obstruction_ok_rows": EXPECTED_ROWS,
        "all_absolute_lift_gaps_positive": True,
        "all_absolute_lift_rows_fail_exact_mu2": True,
        "lambdas": ["0", "1e-6", "1e-4", "1e-2", "1e-1"],
        "shifts": [0, 20],
        "degrees": [2, 8],
    }
    for key, value in expected_grid.items():
        if grid.get(key) != value:
            issues.append(issue("finite_grid", f"bad-{key}", f"{grid.get(key)!r} != {value!r}"))
    for key in ("derived_from", "row_log"):
        ref = grid.get(key)
        if isinstance(ref, str):
            issues.extend(validate_ref("finite_grid", ref))
        else:
            issues.append(issue("finite_grid", f"missing-{key}", repr(ref)))

    summary_block = summary.get("summary", {})
    expected_summary = {
        "symbolic_rows": 3,
        "mu2_sign_lift_obstruction_rows": EXPECTED_ROWS,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary_block.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary_block.get(key)!r} != {value!r}"))
    main_finding = str(summary_block.get("main_finding", "")).lower()
    for required in ("absolute", "mu_2", "2*kappa_1", "state-space doubled", "not just split"):
        if required not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_rows(rows: list[dict]) -> list[SignLiftIssue]:
    issues: list[SignLiftIssue] = []
    if len(rows) != EXPECTED_ROWS:
        issues.append(issue("rows", "bad-row-count", str(len(rows))))
    seen: set[tuple[str, str, int, int]] = set()
    for row in rows:
        key = (
            str(row.get("kind")),
            str(row.get("lam")),
            int(row.get("shift_n", -1)),
            int(row.get("degree_d", -1)),
        )
        if key in seen:
            issues.append(issue(str(key), "duplicate-row", repr(key)))
        seen.add(key)
        expected = {
            "kind": "jensen_window_pf_state_space_sign_lift_mu2_obstruction_row",
            "source_kind": "jensen_window_pf_reciprocal_motzkin_mu2_cancellation_row",
            "raw_mu2_formula": "mu_2=beta_0^2+lambda_1=beta_0^2-kappa_1",
            "absolute_sign_lift_formula": "beta_0^2+kappa_1",
            "gap_formula": "absolute_lift_mu2-mu_2=2*kappa_1=-2*lambda_1",
            "lambda1_classification": "negative",
            "kappa1_classification": "positive",
            "mu2_classification": "positive",
            "absolute_lift_gap_classification": "positive",
            "absolute_lift_matches_mu2": False,
            "ok": True,
        }
        for field, value in expected.items():
            if row.get(field) != value:
                issues.append(issue(str(key), f"bad-{field}", repr(row)))
    return issues


def validate_note(path: Path) -> list[SignLiftIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[SignLiftIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "all state-space models are impossible",
        "the state-space route is dead",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, rows_path: Path, note_path: Path) -> tuple[list[SignLiftIssue], int, int]:
    summary = load_json(summary_path)
    rows = load_jsonl(rows_path)
    issues: list[SignLiftIssue] = []
    issues.extend(validate_summary(summary))
    issues.extend(validate_rows(rows))
    issues.extend(validate_note(note_path))
    symbolic_rows = int(summary.get("summary", {}).get("symbolic_rows", 0))
    finite_rows = int(summary.get("summary", {}).get("mu2_sign_lift_obstruction_rows", 0))
    return issues, symbolic_rows, finite_rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, symbolic_rows, finite_rows = validate(args.summary, args.rows, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "symbolic_rows": symbolic_rows,
                    "mu2_sign_lift_obstruction_rows": finite_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"SIGN-LIFT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF state-space sign-lift obstruction scout: "
            f"{symbolic_rows} symbolic rows, {finite_rows} mu2 sign-lift obstruction rows, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
