#!/usr/bin/env python3
"""Independently validate the Newman local odd-count reduction lemma."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_local_odd_count_reduction_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_local_odd_count_reduction_lemma.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_local_odd_count_reduction_lemma.md"


EXPECTED_IDS = [
    "nloc_01_positive_zero_stieltjes_field",
    "nloc_02_outer_counting_bound",
    "nloc_03_published_uniform_counting",
    "nloc_04_log_squared_localization",
    "nloc_05_odd_count_identity",
    "nloc_06_paired_gap_norm",
    "nloc_07_counting_only_guard",
    "nloc_08_classical_field_birth_guard",
    "nloc_09_fixed_time_scope",
    "nloc_10_xi_collision_handoff",
]


def independent_odd_integral(
    left: list[sp.Rational], right: list[sp.Rational], cutoff: sp.Rational
) -> sp.Expr:
    events: list[tuple[sp.Rational, int]] = [
        *((distance, 1) for distance in left),
        *((distance, -1) for distance in right),
    ]
    grouped: dict[sp.Rational, int] = {}
    for distance, jump in events:
        grouped[distance] = grouped.get(distance, 0) + jump
    discrepancy = 0
    previous: sp.Rational | None = None
    integral = sp.Rational(0)
    for distance in sorted(grouped):
        if previous is not None:
            integral += discrepancy * (1 / previous - 1 / distance)
        discrepancy += grouped[distance]
        previous = distance
    if previous is not None:
        integral += discrepancy * (1 / previous - 1 / cutoff)
    return sp.factor(discrepancy / cutoff + integral)


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = lemma.build_payload()
    if stored != rebuilt:
        issues.append("stored local odd-count payload differs from reconstruction")
    if stored.get("status") != "exact local odd-count reduction with counting and birth guards":
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 theorem rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("expected exactly one open collision handoff")
    if sum(row.get("role") == "exact_countermodel" for row in rows) != 1:
        issues.append("expected exactly one exact heat-birth countermodel")

    c, y, width = sp.symbols("c y width", positive=True)
    kernel = 2 * c / (c**2 - y**2)
    if sp.simplify(kernel - 1 / (c - y) - 1 / (c + y)) != 0:
        issues.append("independent kernel split failed")
    derivative = sp.factor(sp.diff(kernel, y))
    if sp.simplify(derivative - 4 * c * y / (c**2 - y**2) ** 2) != 0:
        issues.append("independent kernel derivative failed")

    left_scaled = sp.factor(width * kernel.subs(y, c - width))
    right_scaled = sp.factor(-width * kernel.subs(y, c + width))
    if sp.simplify(left_scaled.subs(width, c / 2) - sp.Rational(4, 3)) != 0:
        issues.append("left boundary extremal constant failed")
    if sp.limit(left_scaled, width, 0, dir="+") != 1:
        issues.append("left boundary small-window limit failed")
    if sp.simplify(right_scaled.subs(width, c / 2) - sp.Rational(4, 5)) != 0:
        issues.append("right boundary endpoint constant failed")
    if sp.limit(right_scaled, width, 0, dir="+") != 1:
        issues.append("right boundary small-window limit failed")
    if 2 * sp.Rational(4, 3) + 2 != sp.Rational(14, 3):
        issues.append("outer left-plus-mid coefficient failed")

    radius = sp.symbols("radius", positive=True)
    tail = sp.log(2 * radius) / (2 * radius**2) + 1 / (4 * radius**2)
    if sp.simplify(sp.diff(tail, radius) + sp.log(2 * radius) / radius**3) != 0:
        issues.append("independent tail primitive failed")
    if sp.limit(tail, radius, sp.oo) != 0:
        issues.append("independent tail boundary failed")
    tail_at_2c = sp.simplify(tail.subs(radius, 2 * c))
    expected_tail = sp.log(4 * c) / (8 * c**2) + 1 / (16 * c**2)
    if sp.simplify(tail_at_2c - expected_tail) != 0:
        issues.append("tail evaluation at 2c failed")
    tail_coefficient = sp.expand(sp.Rational(64, 9) * c * tail_at_2c)
    expected_coefficient = 8 * sp.log(4 * c) / (9 * c) + sp.Rational(4, 9) / c
    if sp.simplify(tail_coefficient - expected_coefficient) != 0:
        issues.append("tail derivative-bound coefficient failed")

    # Validate Stieltjes integration by parts on exact signed atomic measures.
    c_value = sp.Rational(10)
    for atoms in [
        [(sp.Rational(2), 1), (sp.Rational(5), -1), (sp.Rational(7), 1)],
        [(sp.Rational(4), -1), (sp.Rational(6), 2), (sp.Rational(8), -1)],
        [(sp.Rational(12), 1), (sp.Rational(15), -2), (sp.Rational(18), 1)],
    ]:
        direct = sum(
            weight * 2 * c_value / (c_value**2 - location**2)
            for location, weight in atoms
        )
        ordered = sorted(atoms)
        cumulative = 0
        lower = ordered[0][0]
        upper = ordered[-1][0]
        boundary = 0
        integral = 0
        previous = lower
        for location, weight in ordered:
            if location > previous:
                primitive_change = (
                    2 * c_value / (c_value**2 - location**2)
                    - 2 * c_value / (c_value**2 - previous**2)
                )
                integral += cumulative * primitive_change
            cumulative += weight
            previous = location
        boundary = cumulative * 2 * c_value / (c_value**2 - upper**2)
        first_boundary = 0
        by_parts = boundary - first_boundary - integral
        if sp.simplify(direct - by_parts) != 0:
            issues.append(f"atomic Stieltjes integration failed for {atoms}")

    odd_examples = [
        ([sp.Rational(1, 5), sp.Rational(7, 10)], [sp.Rational(1, 4)], sp.Rational(1)),
        ([sp.Rational(1, 3)], [sp.Rational(2, 5), sp.Rational(4, 5)], sp.Rational(1)),
        (
            [sp.Rational(1, 6), sp.Rational(5, 6), sp.Rational(3, 2)],
            [sp.Rational(1, 5), sp.Rational(4, 5), sp.Rational(7, 4)],
            sp.Rational(2),
        ),
    ]
    for left, right, cutoff in odd_examples:
        direct = sp.factor(sum(1 / item for item in left) - sum(1 / item for item in right))
        by_count = independent_odd_integral(left, right, cutoff)
        if sp.simplify(direct - by_count) != 0:
            issues.append(f"independent odd-count identity failed for {left}, {right}")

    # Verify the exact even heat-flow birth with classical field.
    z, tau = sp.symbols("z tau", real=True)
    a_squared = 1 + 16 / (8 + sp.pi)
    initial = sp.expand((z**2 - 1) ** 2 * (z**2 - a_squared))
    heat = sp.expand(
        sum(
            (-tau) ** order
            / sp.factorial(order)
            * sp.diff(initial, z, 2 * order)
            for order in range(4)
        )
    )
    if sp.simplify(sp.diff(heat, tau) + sp.diff(heat, z, 2)) != 0:
        issues.append("independent even heat-birth PDE failed")
    cofactor = sp.cancel(initial / (z - 1) ** 2)
    field = sp.simplify(sp.diff(cofactor, z).subs(z, 1) / cofactor.subs(z, 1))
    if field != -sp.pi / 8:
        issues.append("independent classical countermodel field failed")
    eps = sp.symbols("eps", positive=True)
    for sign in (-1, 1):
        branch = 1 + sign * sp.sqrt(2) * eps - sp.pi * eps**2 / 4
        residual = sp.series(
            heat.subs({tau: eps**2, z: branch}), eps, 0, 4
        ).removeO()
        if sp.simplify(residual) != 0:
            issues.append(f"independent heat-birth branch failed for sign={sign}")

    exact = stored.get("exact", {})
    published = exact.get("published_counting_input", {})
    if published.get("source") != "https://arxiv.org/abs/1904.12438":
        issues.append("published counting source drifted")
    if "absolute implied constant" not in published.get("statement", ""):
        issues.append("uniform counting constant warning missing")
    if exact.get("local_window") != "H=log(4c)^2":
        issues.append("local window drifted")
    if "field balance alone is compatible" not in exact.get("open_handoff", ""):
        issues.append("field-only non-exclusion warning missing")
    if len(exact.get("checks", {}).get("odd_count_checks", [])) != 3:
        issues.append("stored odd-count check count drifted")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "All mesoscopic and far-field uncertainty is now `O(1/log c)`",
        "inverse-square weight at small `u`",
        "exact even solution of `F_tau=-F_zz`",
        "classical field `-pi/8` and drift `-pi/4`",
        "would not, by itself, exclude a positive Newman",
        "additional Xi-specific global",
        "https://arxiv.org/abs/1904.12438",
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
        "validated Jensen-window PF Newman local odd-count reduction lemma: "
        "10 rows, 0 issues, 3 exact Stieltjes identities, 1 explicit outer bound, "
        "1 log-squared localization theorem, 1 odd-count formula, "
        "3 finite reciprocal-gap checks, 1 published uniform counting input, "
        "1 exact classical-field birth countermodel, 1 open Xi collision-exclusion handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
