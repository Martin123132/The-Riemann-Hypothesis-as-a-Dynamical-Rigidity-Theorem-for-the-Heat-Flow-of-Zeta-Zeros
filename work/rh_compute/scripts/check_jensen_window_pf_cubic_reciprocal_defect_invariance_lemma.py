#!/usr/bin/env python3
"""Validate the cubic reciprocal-defect invariance lemma."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_cubic_reciprocal_defect_invariance_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.md"


def symbolic_issues() -> list[str]:
    issues: list[str] = []
    p, q = sp.symbols("p q", positive=True)
    x = 1 - p**-2
    y = 1 - q**-2
    frontier = x**2 * y**2 - 6 * x * y + 4 * x + 4 * y - 3
    expected = (p - q - 1) * (p - q + 1) * (p + q - 1) * (p + q + 1)
    if sp.factor(frontier * p**4 * q**4 - expected) != 0:
        issues.append("reciprocal-defect frontier factorization failed")

    t, z, n = sp.symbols("t z n", positive=True)
    xb = 1 - t**2
    yb = (1 + 2 * t) / (1 + t) ** 2
    if sp.factor((xb**2 * yb**2 - 6 * xb * yb + 4 * xb + 4 * yb - 3)) != 0:
        issues.append("frontier parameterization failed")
    X, Y = sp.symbols("X Y")
    F = X**2 * Y**2 - 6 * X * Y + 4 * X + 4 * Y - 3
    x_dot_over_r = 2 * X * ((2 * n + 5) * X * Y - 2 * (2 * n + 3) * X + (2 * n + 1))
    y_dot_over_r = 2 * X * Y * ((2 * n + 7) * Y * z - 2 * (2 * n + 5) * Y + (2 * n + 3))
    derivative = sp.factor(
        (sp.diff(F, X) * x_dot_over_r + sp.diff(F, Y) * y_dot_over_r).subs({X: xb, Y: yb})
    )
    threshold = (1 + t) * (1 + 3 * t) / (1 + 2 * t) ** 2
    coefficient = 8 * t**3 * (2 * n + 7) * (1 - t) * (1 + 2 * t) ** 2 / (1 + t) ** 3
    if sp.factor(derivative - coefficient * (z - threshold)) != 0:
        issues.append("heat frontier factorization failed")
    if sp.factor(1 - threshold - t**2 / (1 + 2 * t) ** 2) != 0:
        issues.append("next-defect threshold failed")
    if sp.factor(yb - xb) != t**3 * (t + 2) / (t + 1) ** 2:
        issues.append("frontier monotone-contraction identity failed")
    return issues


def validate_payload(stored: dict) -> list[str]:
    issues: list[str] = []
    rebuilt = lemma.build_payload()
    if stored != rebuilt:
        return ["stored cubic reciprocal-defect payload differs from reconstruction"]
    summary = stored.get("summary", {})
    expected = {
        "rows": 10,
        "exact_coordinate_rows": 6,
        "conditional_invariance_rows": 1,
        "finite_certificate_rows": 2,
        "open_handoff_rows": 1,
        "lambda_minus_100_positive_margins": 318,
        "nonnegative_grid_positive_margins": 310,
        "failed_or_inconclusive_margins": 0,
        "ready_to_apply_rows": 1,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted")
    for key in ("source_hierarchy", "source_flow_cone", "source_m100_entry", "source_rank_two_family"):
        ref = stored.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(f"missing source reference {key}")
    issues.extend(symbolic_issues())
    return issues


def validate_note(path: Path, expected_line: str) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = (
        expected_line,
        "Status: exact cubic reciprocal-defect invariance reduction",
        "0 <= q_(k+1)-q_k <= 1",
        "partial_lambda F/r_(n+1)=C_n(t)*(z-Z(t))",
        "threshold is independent of the",
        "finite-active-set first-crossing argument",
        "all `318` repaired-prefix margins",
        "all `310`",
        "uniformly on compact lambda intervals",
        "Higher degrees still need additional minor coordinates",
        "not a proof",
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
    expected_line = lemma.result_line(stored)
    issues.extend(validate_note(args.note, expected_line))
    for value in issues:
        print(f"ISSUE {value}")
    print(expected_line.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
