#!/usr/bin/env python3
"""Validate shifted Arb reshaped-Hankel sign-consistency certificate artifacts."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import gzip
import json
from math import comb
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n12_k2_k5_N18_dps520_summary.json"
)
DEFAULT_JSONL = (
    REPO_ROOT
    / "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n12_k2_k5_N18_dps520.jsonl.gz"
)
EXPECTED_LAMBDAS = ("0", "1e-6", "1e-4", "1e-2", "1e-1")
EXPECTED_SHIFTS = tuple(range(13))
EXPECTED_ORDERS = (2, 3, 4, 5)
EXPECTED_N_COLS = 18
EXPECTED_DPS = 520
EXPECTED_NEEDED_MAX_K = 33


@dataclass(frozen=True)
class ShiftedManifestIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def open_row_input(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def expected_rows() -> int:
    return len(EXPECTED_LAMBDAS) * len(EXPECTED_SHIFTS) * sum(
        comb(EXPECTED_N_COLS, order) for order in EXPECTED_ORDERS
    )


def validate_summary(summary: dict) -> list[ShiftedManifestIssue]:
    issues: list[ShiftedManifestIssue] = []
    if summary.get("kind") != "arb_shifted_hankel_sign_consistency_summary":
        issues.append(ShiftedManifestIssue("<summary>", "bad-kind", repr(summary.get("kind"))))
    if "not all-order sign consistency" not in str(summary.get("proof_boundary", "")):
        issues.append(ShiftedManifestIssue("<summary>", "weak-proof-boundary", str(summary.get("proof_boundary"))))
    checks = {
        "lambdas": list(EXPECTED_LAMBDAS),
        "shifts": list(EXPECTED_SHIFTS),
        "orders": list(EXPECTED_ORDERS),
        "n_cols": EXPECTED_N_COLS,
        "dps": EXPECTED_DPS,
        "needed_max_k": EXPECTED_NEEDED_MAX_K,
        "summary_rows": len(EXPECTED_LAMBDAS) * len(EXPECTED_SHIFTS) * len(EXPECTED_ORDERS),
        "rows": expected_rows(),
        "ok": expected_rows(),
        "failed_or_inconclusive": 0,
        "all_ok": True,
    }
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            issues.append(ShiftedManifestIssue("<summary>", f"bad-{key}", f"{actual!r} != {expected!r}"))

    subrows = summary.get("summaries", [])
    if not isinstance(subrows, list) or len(subrows) != checks["summary_rows"]:
        issues.append(ShiftedManifestIssue("<summary>", "bad-subsummary-count", repr(len(subrows))))
        return issues

    seen: set[tuple[str, int, int]] = set()
    for row in subrows:
        lam = str(row.get("lam"))
        shift_n = int(row.get("shift_n", -1))
        order = int(row.get("order_k", -1))
        row_id = f"lambda={lam},n={shift_n},k={order}"
        seen.add((lam, shift_n, order))
        tests = comb(EXPECTED_N_COLS, order) if order in EXPECTED_ORDERS else None
        if lam not in EXPECTED_LAMBDAS:
            issues.append(ShiftedManifestIssue(row_id, "unexpected-lambda", lam))
        if shift_n not in EXPECTED_SHIFTS:
            issues.append(ShiftedManifestIssue(row_id, "unexpected-shift", str(shift_n)))
        if order not in EXPECTED_ORDERS:
            issues.append(ShiftedManifestIssue(row_id, "unexpected-order", str(order)))
        if row.get("n_cols") != EXPECTED_N_COLS:
            issues.append(ShiftedManifestIssue(row_id, "bad-n-cols", repr(row.get("n_cols"))))
        if row.get("tests") != tests:
            issues.append(ShiftedManifestIssue(row_id, "bad-tests", f"{row.get('tests')} != {tests}"))
        if row.get("positive") != tests:
            issues.append(ShiftedManifestIssue(row_id, "not-all-positive", repr(row.get("positive"))))
        for key in ("negative", "zero", "inconclusive"):
            if row.get(key) != 0:
                issues.append(ShiftedManifestIssue(row_id, f"nonzero-{key}", repr(row.get(key))))
        if row.get("first_bad_columns") is not None:
            issues.append(ShiftedManifestIssue(row_id, "first-bad-present", repr(row.get("first_bad_columns"))))
        if row.get("all_ok") is not True:
            issues.append(ShiftedManifestIssue(row_id, "not-all-ok", repr(row.get("all_ok"))))
        expected_max_index = shift_n + EXPECTED_N_COLS + order - 2
        if row.get("max_coefficient_index") != expected_max_index:
            issues.append(
                ShiftedManifestIssue(
                    row_id,
                    "bad-max-index",
                    f"{row.get('max_coefficient_index')} != {expected_max_index}",
                )
            )

    expected_seen = {
        (lam, shift_n, order)
        for lam in EXPECTED_LAMBDAS
        for shift_n in EXPECTED_SHIFTS
        for order in EXPECTED_ORDERS
    }
    for lam, shift_n, order in sorted(expected_seen - seen):
        issues.append(ShiftedManifestIssue(f"lambda={lam},n={shift_n},k={order}", "missing-subsummary", "expected row missing"))
    return issues


def validate_jsonl(path: Path) -> list[ShiftedManifestIssue]:
    issues: list[ShiftedManifestIssue] = []
    row_count = 0
    with open_row_input(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            row_count += 1
            row = json.loads(line)
            row_id = f"line={row_count},lambda={row.get('lam')},n={row.get('shift_n')},k={row.get('order_k')}"
            if row.get("kind") != "arb_shifted_hankel_sign_consistency_row":
                issues.append(ShiftedManifestIssue(row_id, "bad-kind", repr(row.get("kind"))))
            if row.get("classification") != "positive":
                issues.append(ShiftedManifestIssue(row_id, "bad-classification", repr(row.get("classification"))))
            if row.get("contains_zero") is not False:
                issues.append(ShiftedManifestIssue(row_id, "contains-zero", repr(row.get("contains_zero"))))
            if row.get("ok") is not True:
                issues.append(ShiftedManifestIssue(row_id, "not-ok", repr(row.get("ok"))))
            shift_n = int(row.get("shift_n", -1))
            order = int(row.get("order_k", -1))
            columns = row.get("columns", [])
            if shift_n not in EXPECTED_SHIFTS:
                issues.append(ShiftedManifestIssue(row_id, "unexpected-shift", str(shift_n)))
            if len(columns) != order:
                issues.append(ShiftedManifestIssue(row_id, "bad-column-count", repr(columns)))
            if row_count > expected_rows() + 1 and len(issues) > 20:
                break
    if row_count != expected_rows():
        issues.append(ShiftedManifestIssue("<jsonl>", "bad-row-count", f"{row_count} != {expected_rows()}"))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    summary = load_json(args.summary)
    issues = validate_summary(summary)
    issues.extend(validate_jsonl(args.jsonl))
    if args.json:
        print(
            json.dumps(
                {
                    "ok": not issues,
                    "rows": summary.get("rows"),
                    "summary_rows": summary.get("summary_rows"),
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for issue in issues:
            print(f"SHIFTED-ARB-HANKEL {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated "
            f"{summary.get('rows')} shifted Arb reshaped-Hankel sign-consistency finite certificates "
            f"with {len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
