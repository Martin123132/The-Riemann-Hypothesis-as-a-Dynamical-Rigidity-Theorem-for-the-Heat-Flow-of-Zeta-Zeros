#!/usr/bin/env python3
"""Validate shifted Arb reshaped-Hankel order-7 certificate artifacts."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import check_arb_shifted_hankel_k6_manifest as shifted_checker


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k7_N15_dps520_summary.json"
)
DEFAULT_JSONL = (
    REPO_ROOT
    / "work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k7_N15_dps520.jsonl.gz"
)


def configure_order_7_spec() -> None:
    shifted_checker.DEFAULT_SUMMARY = DEFAULT_SUMMARY
    shifted_checker.DEFAULT_JSONL = DEFAULT_JSONL
    shifted_checker.EXPECTED_ORDERS = (7,)
    shifted_checker.EXPECTED_N_COLS = 15
    shifted_checker.EXPECTED_NEEDED_MAX_K = 40


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--jsonl", type=Path, default=DEFAULT_JSONL)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    configure_order_7_spec()
    args = build_parser().parse_args()
    summary = shifted_checker.load_json(args.summary)
    issues = shifted_checker.validate_summary(summary)
    issues.extend(shifted_checker.validate_jsonl(args.jsonl))
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
            print(f"SHIFTED-ARB-HANKEL-K7 {issue.row_id} [{issue.issue}] {issue.detail}")
        print(
            "validated "
            f"{summary.get('rows')} shifted Arb reshaped-Hankel order-7 finite certificates "
            f"with {len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
