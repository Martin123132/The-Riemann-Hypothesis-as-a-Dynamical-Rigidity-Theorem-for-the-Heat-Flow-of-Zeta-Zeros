#!/usr/bin/env python3
"""Validate the Jensen-window PF reciprocal signed J-fraction scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json"
DEFAULT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md"

EXPECTED_HANKEL_ROWS = 5 * 21 * sum(range(2, 9))
EXPECTED_LAMBDA_ROWS = 5 * 21 * sum(degree - 1 for degree in range(2, 9))
EXPECTED_TOTAL_ROWS = EXPECTED_HANKEL_ROWS + EXPECTED_LAMBDA_ROWS

REQUIRED_SYMBOLIC_IDS = {
    "signed_hankel_signature_to_signed_j_lambda",
    "endpoint_real_rooted_model_signature",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal Signed J-Fraction Scout",
    "Status: signed J-fraction Hankel-signature diagnostic",
    "This is not a proof of a signed continued-fraction theorem",
    "(-1)^(r(r-1)/2) Delta_r > 0",
    "lambda_n = Delta_{n+1} Delta_{n-1} / Delta_n^2",
    "kappa_n = -lambda_n > 0",
    "work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json",
    "work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_lamgrid_n0_n20_d2_d8_dps520.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_reciprocal_signed_j_fraction_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_signed_j_fraction_scout.py",
    "validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",
    "rp_09_signed_or_modified_continued_fraction",
    "finite evidence only",
)


@dataclass(frozen=True)
class SignedJIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> SignedJIssue:
    return SignedJIssue(section=section, issue=name, detail=detail)


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


def validate_refs(summary: dict) -> list[SignedJIssue]:
    issues: list[SignedJIssue] = []
    grid = summary.get("finite_grid", {})
    for ref in grid.get("source_enclosures", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(issue("finite_grid", "missing-source-enclosure", repr(ref)))
    row_log = grid.get("row_log")
    if not isinstance(row_log, str) or not (REPO_ROOT / row_log).exists():
        issues.append(issue("finite_grid", "missing-row-log", repr(row_log)))
    return issues


def validate_summary(summary: dict) -> list[SignedJIssue]:
    issues: list[SignedJIssue] = []
    if summary.get("kind") != "jensen_window_pf_reciprocal_signed_j_fraction_scout":
        issues.append(issue("<summary>", "bad-kind", repr(summary.get("kind"))))
    if summary.get("target_route_row") != "rp_09_signed_or_modified_continued_fraction":
        issues.append(issue("<summary>", "bad-target-route-row", repr(summary.get("target_route_row"))))
    boundary = str(summary.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<summary>", "weak-proof-boundary", str(summary.get("proof_boundary", ""))))

    symbolic = summary.get("symbolic_rows", [])
    if not isinstance(symbolic, list):
        issues.append(issue("symbolic_rows", "bad-symbolic-rows", repr(type(symbolic))))
        symbolic = []
    by_id = {row.get("id"): row for row in symbolic if isinstance(row, dict)}
    for row_id in sorted(REQUIRED_SYMBOLIC_IDS - set(by_id)):
        issues.append(issue(row_id, "missing-symbolic-row", row_id))
    sign_row = by_id.get("signed_hankel_signature_to_signed_j_lambda", {})
    for required in ("Delta_{n+1}Delta_{n-1}/Delta_n^2", "kappa_n = -lambda_n > 0"):
        if required not in " ".join(str(value) for value in sign_row.values()):
            issues.append(issue("signed_hankel_signature_to_signed_j_lambda", "missing-symbolic-text", required))
    endpoint_row = by_id.get("endpoint_real_rooted_model_signature", {})
    if "circular" not in str(endpoint_row.get("proof_boundary", "")).lower():
        issues.append(issue("endpoint_real_rooted_model_signature", "missing-circular-warning", str(endpoint_row)))

    grid = summary.get("finite_grid", {})
    expected_grid = {
        "dps": 520,
        "signed_hankel_rows": EXPECTED_HANKEL_ROWS,
        "signed_hankel_positive_rows": EXPECTED_HANKEL_ROWS,
        "ordinary_lambda_rows": EXPECTED_LAMBDA_ROWS,
        "ordinary_lambda_negative_rows": EXPECTED_LAMBDA_ROWS,
        "signed_lambda_positive_rows": EXPECTED_LAMBDA_ROWS,
        "all_signed_hankel_rows_positive": True,
        "all_ordinary_lambda_rows_negative": True,
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
        "symbolic_rows": 2,
        "signed_hankel_rows": EXPECTED_HANKEL_ROWS,
        "signed_lambda_rows": EXPECTED_LAMBDA_ROWS,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary_block.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary_block.get(key)!r} != {value!r}"))
    main_finding = str(summary_block.get("main_finding", "")).lower()
    for required in ("signed", "hankel", "negative", "finite"):
        if required not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_rows(rows: list[dict]) -> list[SignedJIssue]:
    issues: list[SignedJIssue] = []
    if len(rows) != EXPECTED_TOTAL_ROWS:
        issues.append(issue("rows", "bad-row-count", str(len(rows))))
    hankel_rows = 0
    lambda_rows = 0
    seen: set[tuple] = set()
    for row in rows:
        kind = row.get("kind")
        if kind == "jensen_window_pf_reciprocal_signed_hankel_row":
            key = (
                kind,
                str(row.get("lam")),
                int(row.get("shift_n", -1)),
                int(row.get("degree_d", -1)),
                int(row.get("determinant_order_r", -1)),
            )
            hankel_rows += 1
            if row.get("signed_delta_classification") != "positive":
                issues.append(issue(str(key), "signed-delta-not-positive", repr(row)))
            if row.get("signed_delta_contains_zero") is not False:
                issues.append(issue(str(key), "signed-delta-contains-zero", repr(row)))
            order = int(row.get("determinant_order_r", -1))
            expected_sign = 1 if (order * (order - 1) // 2) % 2 == 0 else -1
            if row.get("expected_delta_sign") != expected_sign:
                issues.append(issue(str(key), "bad-expected-delta-sign", repr(row)))
        elif kind == "jensen_window_pf_reciprocal_signed_j_lambda_row":
            key = (
                kind,
                str(row.get("lam")),
                int(row.get("shift_n", -1)),
                int(row.get("degree_d", -1)),
                int(row.get("lambda_index", -1)),
            )
            lambda_rows += 1
            if row.get("ordinary_lambda_classification") != "negative":
                issues.append(issue(str(key), "ordinary-lambda-not-negative", repr(row)))
            if row.get("ordinary_lambda_contains_zero") is not False:
                issues.append(issue(str(key), "ordinary-lambda-contains-zero", repr(row)))
            if row.get("signed_lambda_classification") != "positive":
                issues.append(issue(str(key), "signed-lambda-not-positive", repr(row)))
            if row.get("signed_lambda_contains_zero") is not False:
                issues.append(issue(str(key), "signed-lambda-contains-zero", repr(row)))
        else:
            issues.append(issue("rows", "bad-kind", repr(row)))
            continue
        if key in seen:
            issues.append(issue(str(key), "duplicate-row", repr(key)))
        seen.add(key)
        if row.get("ok") is not True:
            issues.append(issue(str(key), "row-not-ok", repr(row)))
    if hankel_rows != EXPECTED_HANKEL_ROWS:
        issues.append(issue("rows", "bad-hankel-row-count", str(hankel_rows)))
    if lambda_rows != EXPECTED_LAMBDA_ROWS:
        issues.append(issue("rows", "bad-lambda-row-count", str(lambda_rows)))
    return issues


def validate_note(path: Path) -> list[SignedJIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[SignedJIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "signed j-fraction proof is complete", "jwpf_06 is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(summary_path: Path, rows_path: Path, note_path: Path) -> tuple[list[SignedJIssue], int, int, int]:
    summary = load_json(summary_path)
    rows = load_jsonl(rows_path)
    issues: list[SignedJIssue] = []
    issues.extend(validate_summary(summary))
    issues.extend(validate_refs(summary))
    issues.extend(validate_rows(rows))
    issues.extend(validate_note(note_path))
    symbolic_count = len(summary.get("symbolic_rows", [])) if isinstance(summary.get("symbolic_rows"), list) else 0
    return issues, symbolic_count, int(summary.get("summary", {}).get("signed_hankel_rows", 0)), int(summary.get("summary", {}).get("signed_lambda_rows", 0))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, symbolic_count, hankel_count, lambda_count = validate(args.summary, args.rows, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "symbolic_rows": symbolic_count,
                    "signed_hankel_rows": hankel_count,
                    "signed_lambda_rows": lambda_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-RECIPROCAL-SIGNED-J {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF reciprocal signed J-fraction scout: "
            f"{symbolic_count} symbolic rows, {hankel_count} signed Hankel rows, "
            f"{lambda_count} signed-lambda rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
