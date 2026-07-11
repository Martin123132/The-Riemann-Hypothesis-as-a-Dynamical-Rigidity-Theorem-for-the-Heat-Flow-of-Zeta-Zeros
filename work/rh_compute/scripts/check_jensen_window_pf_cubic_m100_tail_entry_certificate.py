#!/usr/bin/env python3
"""Validate the lambda=-100 all-k cubic tail-entry certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_cubic_m100_tail_entry_certificate as certificate


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cubic_m100_tail_entry_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cubic_m100_tail_entry_certificate.md"


def symbolic_issues() -> list[str]:
    issues: list[str] = []
    k, m, n = sp.symbols("k m n", integer=True, positive=True)
    defect = sp.expand(
        (sp.Rational(1, 4) / m - 16 / (m - 1) ** 6 - sp.Rational(16, 5) / (m - 1) ** 5 - sp.Rational(1, 5) / m)
        * 20
        * m
        * (m - 1) ** 6
    )
    expected_defect = [1053765790738561, 19820046304378, 155329516751, 649235180, 1526415, 1914, 1]
    actual_defect = certificate.positive_shift_coefficients(defect, m, 320)
    if actual_defect != expected_defect:
        issues.append("defect-floor shifted polynomial drifted")
    epsilon = sp.expand((k - 1) ** 6 - 16 * k**2)
    expected_epsilon = [1034100431206448, 19511328911200, 153390950624, 643148640, 1516860, 1908, 1]
    if certificate.positive_shift_coefficients(epsilon, k, 319) != expected_epsilon:
        issues.append("epsilon shifted polynomial drifted")
    final = sp.expand(k**4 * 99_999**4 - (5 * k + 6) ** 3 * 100_000**4)
    expected_final = [
        625120211416829925149906901121,
        9139382719646770221315447036,
        49046657772633715774210566,
        115094896076559489601276,
        99996000059999600001,
    ]
    if certificate.positive_shift_coefficients(final, k, 319) != expected_final:
        issues.append("final shifted quartic drifted")
    return issues


def validate_payload(stored: dict) -> list[str]:
    issues: list[str] = []
    rebuilt = certificate.build_payload()
    if stored != rebuilt:
        return ["stored cubic tail-entry payload differs from reconstruction"]
    summary = stored.get("summary", {})
    expected = {
        "rows": 10,
        "compact_negative_cumulant_blocks": 4074,
        "analytic_ray_negative_skewness": True,
        "lambda_minus_100_prefix_margins": 318,
        "lambda_minus_100_analytic_tail_start": 319,
        "full_cubic_entry": True,
        "entry_spatial_tail": True,
        "forward_uniform_tail": False,
        "open_handoff_rows": 1,
        "ready_to_apply_rows": 9,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted")
    for key in (
        "source_compact_skewness",
        "source_ray_remainder",
        "source_cumulant_bridge",
        "source_dominance",
        "source_cubic_prefix",
    ):
        ref = stored.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(f"missing source reference {key}")
    issues.extend(symbolic_issues())
    return issues


def validate_note(path: Path, expected_line: str) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = (
        expected_line,
        "Status: all-k lambda=-100 cubic Jensen entry theorem",
        "all `4074`",
        "V''<=5*u^2*q",
        "V'''>=7*u^3*q",
        "49*q>1800000",
        "d_m=1-x_m>=1/(5*m+1)",
        "At `k=319`",
        "final quartic has positive",
        "every shifted",
        "degree-3 Jensen polynomial is hyperbolic",
        "Uniformity of the",
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
    expected_line = certificate.result_line(stored)
    issues.extend(validate_note(args.note, expected_line))
    for value in issues:
        print(f"ISSUE {value}")
    print(expected_line.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
