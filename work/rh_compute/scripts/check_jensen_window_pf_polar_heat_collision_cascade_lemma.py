#!/usr/bin/env python3
"""Independently validate the polar heat-collision cascade lemma."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_polar_heat_collision_cascade_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_polar_heat_collision_cascade_lemma.md"


EXPECTED_IDS = [f"phcc_{index:02d}_{suffix}" for index, suffix in enumerate(
    [
        "adjacent_polar_heat_hierarchy",
        "multiple_root_taylor_recurrence",
        "heat_jet_identity",
        "double_root_splitting",
        "higher_multiplicity_viability",
        "polar_multiplicity_lift",
        "infinite_multiplicity_cascade",
        "exponential_polynomial_classification",
        "unbounded_degree_escape",
        "scaled_tail_handoff",
    ],
    start=1,
)]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = lemma.build_payload()
    if stored != rebuilt:
        issues.append("stored polar cascade payload differs from reconstruction")

    rows = stored.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    if len(rows) != 10:
        issues.append("expected 10 theorem rows")

    exact = stored.get("exact", {})
    checks = exact.get("checks", {})
    if checks.get("sample_degree_range") != [2, 10]:
        issues.append("sample degree range drifted")
    if checks.get("taylor_recurrence_checks") != 210:
        issues.append("Taylor recurrence count drifted")
    if checks.get("heat_jet_checks") != 165:
        issues.append("heat-jet count drifted")
    factor_checks = checks.get("exponential_polynomial_factor_checks", [])
    if len(factor_checks) != 8:
        issues.append("exponential-polynomial factor count drifted")
    for row in factor_checks:
        if row.get("forced_multiplicity") != row.get("degree") - 2:
            issues.append("sample cascade multiplicity failed")

    d, ell, xi, q0 = sp.symbols("d ell xi q0", integer=True, positive=True)
    recurrence = sp.simplify(
        (d + 1 - ell) * sp.binomial(d + 1, ell) * q0 / xi**ell
        - xi * (ell + 1) * sp.binomial(d + 1, ell + 1) * q0 / xi ** (ell + 1)
    )
    if recurrence != 0:
        issues.append("independent Taylor recurrence failed")

    n = sp.symbols("n", integer=True, nonnegative=True)
    for degree in (2, 3, 5, 8, 13):
        for order in range(degree - 1):
            q_r1 = (
                sp.factorial(order + 1)
                * sp.binomial(degree + 1, order + 1)
                * q0
                / xi ** (order + 1)
            )
            q_r2 = (
                sp.factorial(order + 2)
                * sp.binomial(degree + 1, order + 2)
                * q0
                / xi ** (order + 2)
            )
            observed = sp.factor(
                ((4 * n + 2 + 4 * order) * q_r1 + 4 * xi * q_r2)
                / (degree + 1)
            )
            expected = sp.factor(
                (4 * n + 4 * degree + 2)
                * sp.factorial(degree)
                / sp.factorial(degree - order)
                * q0
                / xi ** (order + 1)
            )
            if sp.simplify(observed - expected) != 0:
                issues.append(f"independent heat jet failed at d={degree}, r={order}")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "Status: exact polar heat-collision cascade",
        "not a proof of PF-infinity, RH, or `Lambda <= 0`",
        "least failing degree must tend to infinity",
        "quantitative estimate in the rescaled `D->infinity`",
    ]
    for phrase in required:
        if phrase not in note:
            issues.append(f"note missing required phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Jensen-window PF polar heat-collision cascade lemma: "
        "10 rows, 0 issues, 3 exact local identities, 1 double-root criterion, "
        "1 higher-multiplicity gate, 1 infinite polar cascade, "
        "1 exponential-polynomial classification, 1 unbounded-degree escape theorem, "
        "1 open scaled-tail handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
