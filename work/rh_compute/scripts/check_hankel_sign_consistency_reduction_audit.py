#!/usr/bin/env python3
"""Validate the finite reshaped-Hankel sign-consistency point audit."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_AUDIT = (
    REPO_ROOT
    / "work/rh_compute/results/hankel_sign_consistency_reduction_lamgrid_k2_k5_n18.json"
)
EXPECTED_LAMBDAS = ("0.0", "0.000001", "0.0001", "0.01", "0.1")
EXPECTED_ORDERS = (2, 3, 4, 5)
EXPECTED_N_COLS = 18


@dataclass(frozen=True)
class AuditIssue:
    row_id: str
    issue: str
    detail: str


def load_audit(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_audit(audit: dict) -> list[AuditIssue]:
    issues: list[AuditIssue] = []
    if audit.get("kind") != "hankel_sign_consistency_reduction_point_audit":
        issues.append(AuditIssue("<audit>", "bad-kind", repr(audit.get("kind"))))
    if "not an interval certificate" not in str(audit.get("proof_boundary", "")):
        issues.append(AuditIssue("<audit>", "weak-proof-boundary", str(audit.get("proof_boundary"))))
    if audit.get("orders") != list(EXPECTED_ORDERS):
        issues.append(AuditIssue("<audit>", "bad-orders", repr(audit.get("orders"))))
    if audit.get("n_cols") != EXPECTED_N_COLS:
        issues.append(AuditIssue("<audit>", "bad-n-cols", repr(audit.get("n_cols"))))
    if audit.get("ok") is not True:
        issues.append(AuditIssue("<audit>", "not-ok", repr(audit.get("ok"))))

    rows = audit.get("rows", [])
    expected_count = len(EXPECTED_LAMBDAS) * len(EXPECTED_ORDERS)
    if not isinstance(rows, list) or len(rows) != expected_count:
        issues.append(AuditIssue("<audit>", "bad-row-count", f"{len(rows)} != {expected_count}"))
        return issues

    seen: set[tuple[str, int]] = set()
    for row in rows:
        lam = str(row.get("lam"))
        order_k = int(row.get("order_k", -1))
        row_id = f"lambda={lam},k={order_k}"
        seen.add((lam, order_k))

        if lam not in EXPECTED_LAMBDAS:
            issues.append(AuditIssue(row_id, "unexpected-lambda", lam))
        if order_k not in EXPECTED_ORDERS:
            issues.append(AuditIssue(row_id, "unexpected-order", str(order_k)))
        tests = math.comb(EXPECTED_N_COLS, order_k)
        if row.get("tests") != tests:
            issues.append(AuditIssue(row_id, "bad-test-count", f"{row.get('tests')} != {tests}"))
        if row.get("positive") != tests:
            issues.append(AuditIssue(row_id, "not-all-positive", repr(row.get("positive"))))
        if row.get("negative") != 0:
            issues.append(AuditIssue(row_id, "negative-minors", repr(row.get("negative"))))
        if row.get("zero") != 0:
            issues.append(AuditIssue(row_id, "zero-minors", repr(row.get("zero"))))
        if row.get("first_bad_columns") is not None:
            issues.append(AuditIssue(row_id, "first-bad-present", repr(row.get("first_bad_columns"))))
        if row.get("minimum_signed_value") in (None, ""):
            issues.append(AuditIssue(row_id, "missing-minimum", repr(row.get("minimum_signed_value"))))
        if int(row.get("coefficient_common_denominator_digits", 0)) <= 0:
            issues.append(
                AuditIssue(
                    row_id,
                    "missing-common-denominator-size",
                    repr(row.get("coefficient_common_denominator_digits")),
                )
            )
        expected_max_index = EXPECTED_N_COLS + order_k - 2
        if row.get("max_coefficient_index") != expected_max_index:
            issues.append(
                AuditIssue(
                    row_id,
                    "bad-max-index",
                    f"{row.get('max_coefficient_index')} != {expected_max_index}",
                )
            )

    expected_seen = {(lam, order) for lam in EXPECTED_LAMBDAS for order in EXPECTED_ORDERS}
    missing = sorted(expected_seen - seen)
    for lam, order in missing:
        issues.append(AuditIssue(f"lambda={lam},k={order}", "missing-row", "expected row missing"))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    audit = load_audit(args.audit)
    issues = validate_audit(audit)
    rows = audit.get("rows", [])
    ok_rows = sum(
        1
        for row in rows
        if row.get("positive") == row.get("tests")
        and row.get("negative") == 0
        and row.get("zero") == 0
    )

    if args.json:
        print(
            json.dumps(
                {
                    "ok": not issues,
                    "rows": len(rows),
                    "ok_rows": ok_rows,
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for issue in issues:
            print(f"AUDIT {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated "
            f"{ok_rows} reshaped Hankel sign-consistency point audits "
            f"with {len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
