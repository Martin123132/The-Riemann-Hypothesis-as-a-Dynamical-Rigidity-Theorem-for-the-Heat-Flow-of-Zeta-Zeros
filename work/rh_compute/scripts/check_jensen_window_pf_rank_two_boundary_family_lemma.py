#!/usr/bin/env python3
"""Validate the rank-two all-degree boundary-family lemma."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import comb
from pathlib import Path

import sympy as sp

import jensen_window_pf_rank_two_boundary_family_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_rank_two_boundary_family_lemma.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_rank_two_boundary_family_lemma.md"
EXPECTED = (
    "validated Jensen-window PF rank-two boundary-family lemma: "
    "11 rows, 0 issues, 4 exact identities, 1 all-degree factorization, "
    "1 integer-product closure, 2 exact countermodels, 1 open structural handoff"
)


def family_value(u: Fraction, k: int) -> Fraction:
    a = (2 * u - 1) / u
    b = 1 - u
    c = 2 * u - 1
    return a ** (k - 1) * (c + k * b) / u


def cubic_discriminant(coefficients: list[Fraction]) -> Fraction:
    d, c, b, a = coefficients
    return b * b * c * c - 4 * a * c**3 - 4 * b**3 * d - 27 * a * a * d * d + 18 * a * b * c * d


def exact_issues() -> list[str]:
    issues: list[str] = []
    d, j, a, b, cn = sp.symbols("d j a b cn", positive=True)
    q = a * (cn + d * b) / cn
    factor_coefficient = (d - j) * a**j / d + j * q * a ** (j - 1) / d
    window_coefficient = a**j * (cn + j * b) / cn
    if sp.factor(factor_coefficient - window_coefficient) != 0:
        issues.append("generic all-degree coefficient factorization failed")
    k, c = sp.symbols("k c", positive=True)
    contraction = sp.factor(((c + (k + 1) * b) * (c + (k - 1) * b)) / (c + k * b) ** 2)
    if sp.factor(contraction - (1 - b**2 / (c + k * b) ** 2)) != 0:
        issues.append("ratio contraction identity failed")
    alpha = sp.symbols("alpha", positive=True)
    q_alpha = alpha / (alpha + 1)
    atom_alpha = q_alpha ** (k - 1) * (alpha + k) / (alpha + 1)
    egf_coefficient = q_alpha**k + k * q_alpha ** (k - 1) / (alpha + 1)
    if sp.factor(atom_alpha - egf_coefficient) != 0:
        issues.append("elementary multiplier EGF coefficient identity failed")
    for u in (Fraction(3, 5), Fraction(2, 3)):
        if family_value(u, 0) != 1 or family_value(u, 1) != 1:
            issues.append(f"family normalization failed at u={u}")
    u1, u2 = Fraction(3, 5), Fraction(2, 3)
    coefficients = [
        Fraction(comb(3, index), 2) * (family_value(u1, index) + family_value(u2, index))
        for index in range(4)
    ]
    if coefficients != [Fraction(1), Fraction(3), Fraction(47, 24), Fraction(41, 108)]:
        issues.append(f"mixture coefficients drifted: {coefficients}")
    if cubic_discriminant(coefficients) != Fraction(-937, 3456):
        issues.append("mixture discriminant drifted")
    radical_discriminant = -432 * sp.sqrt(2) - 324 * sp.sqrt(3) + 378 + 324 * sp.sqrt(6)
    rational_upper_bound = (
        -432 * sp.Rational(707, 500)
        - 324 * sp.Rational(433, 250)
        + 378
        + 324 * sp.Rational(49, 20)
    )
    if rational_upper_bound != sp.Rational(-27, 125):
        issues.append("fractional-power rational upper bound drifted")
    if not bool(sp.N(radical_discriminant, 50) < rational_upper_bound < 0):
        issues.append("fractional-power cubic countermodel failed")
    return issues


def validate_payload(stored: dict) -> list[str]:
    issues: list[str] = []
    if stored != lemma.build_payload():
        return ["stored boundary-family payload differs from reconstruction"]
    expected_summary = {
        "rows": 11,
        "exact_identity_rows": 4,
        "all_degree_factorization_rows": 1,
        "integer_product_closure_rows": 1,
        "exact_countermodels": 2,
        "open_structural_handoffs": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    summary = stored.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(f"summary {key} drifted")
    for key in ("source_hierarchy", "source_bridge_target"):
        ref = stored.get(key)
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            issues.append(f"missing source reference {key}")
    issues.extend(exact_issues())
    return issues


def validate_note(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    required = (
        EXPECTED,
        "Status: exact all-degree boundary-family lemma with mixture countermodel.",
        "J_(d,n)(z)=A_n*(1+a*z)^(d-1)",
        "d_k=1/(k+c/b)^2",
        "Disc(J_3)=-937/3456<0.",
        "(1+z/(alpha+1))*exp(alpha*z/(alpha+1))",
        "-log x_k=sum_j -log(1-1/(k+alpha_j)^2)",
        "<-27/125<0.",
        "discrete counting-measure factorization",
        "not, by",
        "itself, the missing zeta kernel representation",
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
    issues.extend(validate_note(args.note))
    for value in issues:
        print(f"ISSUE {value}")
    print(EXPECTED.replace("0 issues", f"{len(issues)} issues"))
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
