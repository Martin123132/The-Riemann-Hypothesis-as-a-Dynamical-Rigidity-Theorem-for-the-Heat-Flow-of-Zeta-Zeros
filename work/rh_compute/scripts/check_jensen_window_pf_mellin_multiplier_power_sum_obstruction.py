#!/usr/bin/env python3
"""Validate the Mellin multiplier power-sum obstruction certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_mellin_multiplier_power_sum_obstruction as obstruction


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_mellin_multiplier_power_sum_obstruction.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md"


def exact_issues() -> list[str]:
    issues: list[str] = []
    s, alpha = sp.symbols("s alpha", positive=True)
    atom_log = (s - 1) * sp.log(alpha / (alpha + 1)) + sp.log(alpha + s) - sp.log(alpha + 1)
    for m in range(2, 9):
        derivative = sp.simplify(sp.diff(atom_log, s, m).subs(s, 0))
        expected = (-1) ** (m - 1) * sp.factorial(m - 1) / alpha**m
        if sp.simplify(derivative - expected) != 0:
            issues.append(f"atom log derivative failed at m={m}")
    return issues


def validate_payload(stored: dict) -> list[str]:
    issues: list[str] = []
    rebuilt = obstruction.build_payload()
    if stored != rebuilt:
        return ["stored Mellin obstruction differs from rigorous reconstruction"]
    summary = stored.get("summary", {})
    expected = {
        "rows": 6,
        "raw_log_moments": 11,
        "power_sum_candidates": 9,
        "hankel_determinants": 6,
        "positive_hankel_determinants": 3,
        "negative_hankel_determinants": 3,
        "inconclusive_hankel_determinants": 0,
        "continuous_product_ruled_out": True,
        "integer_only_product_ruled_out": False,
        "ready_to_apply_rows": 1,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted")
    if not stored.get("tail_diagnostics", {}).get("all_tail_errors_below_inflation"):
        issues.append("analytic tail budget is not closed")
    for row in stored.get("hankel_rows", []):
        if row.get("classification") in {"positive", "negative"} and row.get("contains_zero"):
            issues.append(f"promoted determinant contains zero: {row}")
    if sum(row.get("classification") == "negative" for row in stored.get("hankel_rows", [])) != 3:
        issues.append("negative Hankel count drifted")
    for ref in stored.get("source_artifacts", []):
        if not (REPO_ROOT / ref).exists():
            issues.append(f"missing source artifact: {ref}")
    issues.extend(exact_issues())
    return issues


def validate_note(path: Path, expected_line: str) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = (
        expected_line,
        "Status: interval-certified continuous-interpolation obstruction.",
        "A(k)=A_k",
        "sum_j alpha_j^(-m)",
        "positive semidefinite",
        "shift-2 size-4 determinant is strictly negative",
        "does not by itself force equality with this Mellin",
        "discrete counting-measure route therefore remains open",
        "not a",
    )
    return [f"missing note text: {value}" for value in required if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    stored = json.loads(args.result.read_text(encoding="utf-8"))
    issues = validate_payload(stored)
    expected_line = obstruction.result_line(stored)
    issues.extend(validate_note(args.note, expected_line))
    for value in issues:
        print(f"ISSUE {value}")
    print(expected_line.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
