#!/usr/bin/env python3
"""Validate the promoted shifted Arb reshaped-Hankel staircase certificates."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]


@dataclass(frozen=True)
class StaircaseSpec:
    name: str
    checker: str
    summary: str
    orders: tuple[int, ...]
    n_cols: int
    rows: int
    needed_max_k: int
    expected_line: str


@dataclass(frozen=True)
class StaircaseIssue:
    spec: str
    issue: str
    detail: str


SPECS: tuple[StaircaseSpec, ...] = (
    StaircaseSpec(
        name="k2-k5_N18",
        checker="work/rh_compute/scripts/check_arb_shifted_hankel_sign_consistency_manifest.py",
        summary="work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k2_k5_N18_dps520_summary.json",
        orders=(2, 3, 4, 5),
        n_cols=18,
        rows=1_322_685,
        needed_max_k=41,
        expected_line="validated 1322685 shifted Arb reshaped-Hankel sign-consistency finite certificates with 0 issues",
    ),
    StaircaseSpec(
        name="k6_N16",
        checker="work/rh_compute/scripts/check_arb_shifted_hankel_k6_manifest.py",
        summary="work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k6_N16_dps520_summary.json",
        orders=(6,),
        n_cols=16,
        rows=840_840,
        needed_max_k=40,
        expected_line="validated 840840 shifted Arb reshaped-Hankel order-6 finite certificates with 0 issues",
    ),
    StaircaseSpec(
        name="k7_N15",
        checker="work/rh_compute/scripts/check_arb_shifted_hankel_k7_manifest.py",
        summary="work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k7_N15_dps520_summary.json",
        orders=(7,),
        n_cols=15,
        rows=675_675,
        needed_max_k=40,
        expected_line="validated 675675 shifted Arb reshaped-Hankel order-7 finite certificates with 0 issues",
    ),
    StaircaseSpec(
        name="k8_N14",
        checker="work/rh_compute/scripts/check_arb_shifted_hankel_k8_manifest.py",
        summary="work/rh_compute/results/arb_shifted_hankel_sign_consistency_lamgrid_n0_n20_k8_N14_dps520_summary.json",
        orders=(8,),
        n_cols=14,
        rows=315_315,
        needed_max_k=40,
        expected_line="validated 315315 shifted Arb reshaped-Hankel order-8 finite certificates with 0 issues",
    ),
)

EXPECTED_TOTAL_ROWS = 3_154_515
EXPECTED_LAMBDAS = ["0", "1e-6", "1e-4", "1e-2", "1e-1"]
EXPECTED_SHIFTS = list(range(21))
EXPECTED_DPS = 520


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def run_checker(spec: StaircaseSpec) -> list[StaircaseIssue]:
    completed = subprocess.run(
        [sys.executable, spec.checker],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    combined = completed.stdout + "\n" + completed.stderr
    issues: list[StaircaseIssue] = []
    if completed.returncode != 0:
        issues.append(StaircaseIssue(spec.name, "checker-failed", combined.strip()))
    if spec.expected_line not in combined:
        issues.append(StaircaseIssue(spec.name, "missing-expected-line", spec.expected_line))
    return issues


def validate_summary(spec: StaircaseSpec) -> tuple[int, list[StaircaseIssue]]:
    summary_path = REPO_ROOT / spec.summary
    if not summary_path.exists():
        return 0, [StaircaseIssue(spec.name, "missing-summary", spec.summary)]
    summary = load_json(summary_path)
    checks = {
        "lambdas": EXPECTED_LAMBDAS,
        "shifts": EXPECTED_SHIFTS,
        "orders": list(spec.orders),
        "n_cols": spec.n_cols,
        "dps": EXPECTED_DPS,
        "needed_max_k": spec.needed_max_k,
        "rows": spec.rows,
        "ok": spec.rows,
        "failed_or_inconclusive": 0,
        "all_ok": True,
    }
    issues: list[StaircaseIssue] = []
    if "not all-order sign consistency" not in str(summary.get("proof_boundary", "")):
        issues.append(StaircaseIssue(spec.name, "weak-proof-boundary", str(summary.get("proof_boundary"))))
    for key, expected in checks.items():
        actual = summary.get(key)
        if actual != expected:
            issues.append(StaircaseIssue(spec.name, f"bad-{key}", f"{actual!r} != {expected!r}"))
    return int(summary.get("rows", 0)), issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues: list[StaircaseIssue] = []
    total_rows = 0
    for spec in SPECS:
        issues.extend(run_checker(spec))
        rows, summary_issues = validate_summary(spec)
        total_rows += rows
        issues.extend(summary_issues)
    if total_rows != EXPECTED_TOTAL_ROWS:
        issues.append(StaircaseIssue("<total>", "bad-total-rows", f"{total_rows} != {EXPECTED_TOTAL_ROWS}"))

    if args.json:
        print(
            json.dumps(
                {
                    "ok": not issues,
                    "total_rows": total_rows,
                    "issues": [asdict(issue) for issue in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for issue in issues:
            print(f"SHIFTED-STAIRCASE {issue.spec} [{issue.issue}] {issue.detail}")
        print(
            "validated "
            f"{total_rows} shifted Arb reshaped-Hankel staircase finite certificates "
            f"with {len(issues)} issues"
        )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
