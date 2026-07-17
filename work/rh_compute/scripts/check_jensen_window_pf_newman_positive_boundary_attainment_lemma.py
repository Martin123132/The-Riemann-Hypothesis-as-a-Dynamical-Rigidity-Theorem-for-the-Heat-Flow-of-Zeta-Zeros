#!/usr/bin/env python3
"""Independently validate the Newman positive-boundary attainment lemma."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_positive_boundary_attainment_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_positive_boundary_attainment_lemma.md"


EXPECTED_IDS = [
    "npba_01_debruijn_strip",
    "npba_02_uniform_high_zero_reality",
    "npba_03_positive_boundary_compactness",
    "npba_04_boundary_multiplicity_limit",
    "npba_05_positive_boundary_attainment",
    "npba_06_positive_time_simplicity_equivalence",
    "npba_07_hermite_cluster_split",
    "npba_08_hermite_reciprocal_stiffness",
    "npba_09_cluster_energy_blowup",
    "npba_10_positive_simplicity_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = lemma.build_payload()
    if stored != rebuilt:
        issues.append("stored positive-boundary payload differs from reconstruction")
    if stored.get("status") != "exact positive-boundary attainment and arbitrary-multiplicity energy lemma":
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 theorem rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    for role, expected in [
        ("published_theorem", 3),
        ("exact_theorem", 1),
        ("exact_equivalence", 1),
        ("open_handoff", 1),
    ]:
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    exact = stored.get("exact", {})
    strip = exact.get("published_strip", {})
    high = exact.get("published_high_zero", {})
    expected_source = "https://arxiv.org/abs/1904.12438"
    if strip.get("source") != expected_source or high.get("source") != expected_source:
        issues.append("positive-time primary source drifted")
    if "Theorem 3.2" not in strip.get("theorem", ""):
        issues.append("strip theorem number missing")
    if "Theorem 1.5" not in high.get("theorem", ""):
        issues.append("high-zero theorem number missing")
    if high.get("absolute_constants") is not True:
        issues.append("absolute-constant uniformity flag missing")
    if "0<t<=1/2" not in high.get("statement", ""):
        issues.append("positive-time range missing")

    compactness = exact.get("positive_boundary_compactness", "")
    if "Lambda/2<=t_n<Lambda" not in compactness:
        issues.append("boundary time sequence missing")
    if "exp(2*C/Lambda)" not in compactness:
        issues.append("uniform compact threshold missing")
    if "|Im z_n|<=1" not in compactness:
        issues.append("vertical compactness missing")
    if "locally uniformly" not in exact.get("local_uniform_convergence", ""):
        issues.append("local-uniform convergence step missing")
    multiplicity = exact.get("multiplicity_argument", "")
    for phrase in ("conjugate(z_n)", "Rouche", ">=2"):
        if phrase not in multiplicity:
            issues.append(f"multiplicity step missing: {phrase}")
    if "cannot be realized solely" not in exact.get("attainment_theorem", ""):
        issues.append("height-escape exclusion missing")
    if "if and only if" not in exact.get("positive_time_simplicity_equivalence", ""):
        issues.append("positive-time simplicity equivalence missing")
    if "every 0<t<=1/5" not in exact.get("positive_time_simplicity_equivalence", ""):
        issues.append("published reduced simplicity window missing")

    # Re-derive the Hermite ODE and reciprocal-square identities algebraically.
    x = sp.symbols("x", real=True)
    numeric_rows: list[dict] = []
    for m in range(2, 11):
        polynomial = sp.hermite_prob(m, x)
        if sp.expand(sp.diff(polynomial, x, 2) - x * sp.diff(polynomial, x) + m * polynomial) != 0:
            issues.append(f"independent Hermite ODE failed for m={m}")
            continue
        derivative = sp.diff(polynomial, x)
        third = sp.diff(polynomial, x, 3)
        first_field = sp.diff(polynomial, x, 2) / (2 * derivative)
        stiffness = first_field**2 - third / (3 * derivative)
        target_stiffness = (4 * (m - 1) - x**2) / 12
        numerator = sp.cancel(stiffness - target_stiffness).as_numer_denom()[0]
        remainder = sp.rem(sp.Poly(numerator, x), sp.Poly(polynomial, x)).as_expr()
        if sp.expand(remainder) != 0:
            issues.append(f"rootwise reciprocal-square identity failed for m={m}")

        coefficient = sp.Poly(polynomial, x).nth(m - 2)
        square_sum = sp.factor(-2 * coefficient)
        if square_sum != m * (m - 1):
            issues.append(f"root-square sum failed for m={m}")

        roots = [complex(root) for root in sp.nroots(polynomial, n=40, maxsteps=200)]
        ordered = sum(
            1 / (roots[a] - roots[b]) ** 2
            for a in range(m)
            for b in range(m)
            if a != b
        )
        expected_ordered = m * (m - 1) / 4
        if abs(ordered.imag) > 1e-10 or abs(ordered.real - expected_ordered) > 1e-9:
            issues.append(f"numeric ordered Hermite sum failed for m={m}")
        numeric_rows.append(
            {
                "m": m,
                "ordered": ordered.real,
                "expected": expected_ordered,
            }
        )

    stored_rows = exact.get("checks", {}).get("hermite_rows", [])
    if len(stored_rows) != 9:
        issues.append("expected nine stored Hermite checks")
    else:
        for row in stored_rows:
            m = int(row.get("multiplicity", 0))
            expected_ordered = sp.Rational(m * (m - 1), 4)
            if sp.sympify(row.get("ordered_reciprocal_gap_square_sum", "nan")) != expected_ordered:
                issues.append(f"stored ordered Hermite sum failed for m={m}")
            expected_energy = sp.Rational(m * (m - 1), 8)
            if sp.sympify(row.get("ordered_cluster_energy_coefficient", "nan")) != expected_energy:
                issues.append(f"stored cluster energy coefficient failed for m={m}")

    tau = sp.symbols("tau", positive=True)
    for m in range(2, 11):
        reciprocal_sum = sp.Rational(m * (m - 1), 4)
        leading = reciprocal_sum / (2 * tau)
        if sp.simplify(leading - sp.Rational(m * (m - 1), 8) / tau) != 0:
            issues.append(f"cluster energy leading term failed for m={m}")
        cutoff = sp.symbols(f"cutoff_{m}", positive=True)
        integral = sp.integrate(
            sp.Rational(m * (m - 1), 8) / tau,
            (tau, cutoff, 1),
        )
        if sp.limit(integral, cutoff, 0, dir="+") != sp.oo:
            issues.append(f"cluster logarithmic divergence failed for m={m}")

    cluster = exact.get("cluster_energy", {})
    if "m*(m-1)/(8*tau)" not in cluster.get("ordered_inverse_square_sum", ""):
        issues.append("cluster coefficient text missing")
    if "every m>=2" not in cluster.get("integral", ""):
        issues.append("arbitrary multiplicity warning missing")
    if "finite truncation containing it" not in exact.get("conditional_global_exclusion", ""):
        issues.append("finite-cluster exclusion criterion missing")
    if "positive-time simplicity" not in exact.get("open_handoff", ""):
        issues.append("positive-time simplicity handoff missing")
    if "every 0<t<=1/5" not in exact.get("open_handoff", ""):
        issues.append("reduced positive-time handoff missing")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "force every such choice into",
        "Rouche zero counting",
        "closes the high-index escape loophole",
        "positive-time simplicity equivalence",
        "probabilists' Hermite",
        "No separate",
        "endpoint-integrability",
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
        "validated Jensen-window PF Newman positive-boundary attainment lemma: "
        "10 rows, 0 issues, 2 published compactness inputs, "
        "1 finite-boundary attainment theorem, 1 positive-time simplicity equivalence, "
        "1 arbitrary-multiplicity Hermite split, 9 exact Hermite checks, "
        "1 cluster-energy blow-up, 1 open Xi endpoint handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
