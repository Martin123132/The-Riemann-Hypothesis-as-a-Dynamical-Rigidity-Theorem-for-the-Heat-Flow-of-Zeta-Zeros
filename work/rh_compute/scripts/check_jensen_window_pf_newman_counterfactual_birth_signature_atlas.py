#!/usr/bin/env python3
"""Independently validate the Newman counterfactual birth-signature atlas."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_counterfactual_birth_signature_atlas as atlas


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_counterfactual_birth_signature_atlas.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_counterfactual_birth_signature_atlas.md"
)
STRICT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_strict_laguerre_correlation_target.json"
)
BOUNDARY_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_positive_boundary_attainment_lemma.json"
)
SUZUKI_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_fixed_omega_phase_diagram.json"
)

EXPECTED_IDS = [
    "ncbsa_01_quartet_coordinate_map",
    "ncbsa_02_symmetric_quartet_heat_model",
    "ncbsa_03_universal_laguerre_birth_jet",
    "ncbsa_04_correlation_zero_mode",
    "ncbsa_05_quadratic_jensen_threshold",
    "ncbsa_06_li_quartet_signature",
    "ncbsa_07_suzuki_shift_signature",
    "ncbsa_08_diagonal_escape_theorem",
    "ncbsa_09_finite_sensor_nonpromotion",
    "ncbsa_10_primary_route_selection",
    "ncbsa_11_phase_critical_handoff",
    "ncbsa_12_new_math_alternative",
]


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.exists():
        issues.append(f"missing {label}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "stored atlas", issues)
    if not NOTE.exists():
        issues.append("missing rendered note")
    if not stored:
        return issues

    rebuilt = atlas.build_payload()
    if stored != rebuilt:
        issues.append("stored atlas payload differs from reconstruction")
    if stored.get("kind") != (
        "jensen_window_pf_newman_counterfactual_birth_signature_atlas"
    ):
        issues.append("kind drifted")
    if stored.get("status") != (
        "exact counterfactual birth-signature atlas and "
        "finite-sensor nonpromotion gate"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 12:
        issues.append("expected 12 atlas rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "exact_identity": 1,
        "exact_countermodel": 2,
        "exact_lemma": 3,
        "exact_composition": 2,
        "nonpromotion_gate": 1,
        "route_decision": 1,
        "open_handoff": 2,
    }
    for role, expected in expected_roles.items():
        observed = sum(row.get("role") == role for row in rows)
        if observed != expected:
            issues.append(
                f"expected {expected} rows with role {role}, got {observed}"
            )

    # Re-derive the zeta-to-H and squared-zero coordinate map.
    epsilon, gamma = sp.symbols("epsilon gamma", positive=True)
    h_zero = 2 * gamma - 2 * sp.I * epsilon
    squared_zero = sp.expand(-(h_zero**2))
    expected_squared = (
        -4 * (gamma**2 - epsilon**2)
        + 8 * sp.I * gamma * epsilon
    )
    if sp.simplify(squared_zero - expected_squared) != 0:
        issues.append("independent squared-zero coordinate map failed")
    sector_sine_squared = sp.factor(
        4 * gamma**2 * epsilon**2 / (gamma**2 + epsilon**2) ** 2
    )

    # Reconstruct the exact four-zero heat model without using builder output.
    z, t = sp.symbols("z t", real=True)
    x, y = sp.symbols("x y", positive=True)
    initial = sp.expand(
        ((z - x) ** 2 + y**2) * ((z + x) ** 2 + y**2)
    )
    heat = sp.expand(
        initial - t * sp.diff(initial, z, 2)
        + t**2 * sp.diff(initial, z, 4) / 2
    )
    if sp.simplify(sp.diff(heat, t) + sp.diff(heat, z, 2)) != 0:
        issues.append("independent quartet heat equation failed")

    w = sp.symbols("w", real=True)
    heat_w = (
        w**2
        + 2 * (y**2 - x**2 - 6 * t) * w
        + (x**2 + y**2) ** 2
        + 4 * t * (x**2 - y**2)
        + 12 * t**2
    )
    if sp.simplify(heat - heat_w.subs(w, z**2)) != 0:
        issues.append("independent squared-variable heat polynomial failed")
    discriminant = sp.factor(sp.discriminant(heat_w, w))
    expected_discriminant = 16 * (
        6 * t**2 + 2 * t * (x**2 - y**2) - x**2 * y**2
    )
    if sp.simplify(discriminant - expected_discriminant) != 0:
        issues.append("independent quartet discriminant failed")

    root = sp.sqrt(x**4 + 4 * x**2 * y**2 + y**4)
    collision_time = (y**2 - x**2 + root) / 6
    if sp.simplify(
        heat.subs(t, collision_time) - (z**2 - root) ** 2
    ) != 0:
        issues.append("independent collision factorization failed")
    if not bool(
        sp.N(collision_time.subs({x: sp.Rational(7, 3), y: sp.Rational(2, 5)}))
        > 0
    ):
        issues.append("collision time is not positive in a rational test case")

    # Check the universal double-zero birth jet directly.
    q, h, amplitude = sp.symbols("q h amplitude", real=True)
    local = amplitude * (q**2 / 2 - h)
    laguerre = sp.expand(
        sp.diff(local, q) ** 2 - local * sp.diff(local, q, 2)
    )
    if sp.simplify(
        laguerre - amplitude**2 * (q**2 / 2 + h)
    ) != 0:
        issues.append("independent universal Laguerre birth jet failed")

    # Re-derive the quadratic Jensen threshold from exponential coefficients.
    delta, radius = sp.symbols("delta radius", positive=True)
    degree = sp.symbols("degree", integer=True, positive=True)
    variable = sp.symbols("variable", real=True)
    jensen = (
        1
        + 2 * degree * sp.cos(delta) * variable / radius
        + degree * (degree - 1) * variable**2 / radius**2
    )
    jensen_discriminant = sp.factor(sp.discriminant(jensen, variable))
    expected_jensen = (
        4 * degree * (1 - degree * sp.sin(delta) ** 2) / radius**2
    )
    if sp.simplify(sp.trigsimp(jensen_discriminant - expected_jensen)) != 0:
        issues.append("independent Jensen discriminant failed")
    if sp.signsimp(
        expected_jensen.subs(
            {degree: 4, delta: sp.pi / 6, radius: 3}
        )
    ) != 0:
        issues.append("Jensen threshold equality test failed")
    if not bool(
        expected_jensen.subs(
            {degree: 5, delta: sp.pi / 6, radius: 3}
        )
        < 0
    ):
        issues.append("Jensen first-failure side test failed")
    continuous_threshold = sp.factor(1 / sector_sine_squared)

    # Derive the reciprocal Li quartet law and its amplification rate.
    rho = sp.Rational(1, 2) + epsilon + sp.I * gamma
    multiplier = sp.simplify((rho - 1) / rho)
    radius_squared = sp.simplify(multiplier * sp.conjugate(multiplier))
    expected_radius_squared = (
        gamma**2 + (sp.Rational(1, 2) - epsilon) ** 2
    ) / (
        gamma**2 + (sp.Rational(1, 2) + epsilon) ** 2
    )
    if sp.simplify(radius_squared - expected_radius_squared) != 0:
        issues.append("independent Li multiplier radius failed")
    if not bool(
        expected_radius_squared.subs(
            {epsilon: sp.Rational(1, 10), gamma: 7}
        )
        < 1
    ):
        issues.append("Li multiplier radius is not below one")
    r, theta = sp.symbols("r theta", positive=True, real=True)
    n = sp.symbols("n", integer=True, positive=True)
    quartet_sum = sp.expand(
        (1 - r**n * sp.exp(sp.I * n * theta))
        + (1 - r**n * sp.exp(-sp.I * n * theta))
        + (1 - r ** (-n) * sp.exp(sp.I * n * theta))
        + (1 - r ** (-n) * sp.exp(-sp.I * n * theta))
    )
    expected_quartet_sum = (
        4 - 2 * (r**n + r ** (-n)) * sp.cos(n * theta)
    )
    if sp.simplify(sp.expand_complex(quartet_sum) - expected_quartet_sum) != 0:
        issues.append("independent Li quartet contribution failed")

    # The diagonal family must evade all five fixed sensor cutoffs at once.
    m = sp.symbols("m", positive=True)
    diagonal_collision = sp.factor(
        collision_time.subs({x: 2 * m, y: 2 / m})
    )
    diagonal_jensen = sp.factor(
        continuous_threshold.subs({gamma: m, epsilon: 1 / m})
    )
    growth = sp.log(1 / expected_radius_squared) / 2
    diagonal_growth = sp.factor(
        growth.subs({gamma: m, epsilon: 1 / m})
    )
    diagonal_multiplier = sp.simplify(
        multiplier.subs({gamma: m, epsilon: 1 / m})
    )
    diagonal_li_contribution = sp.simplify(
        4
        - diagonal_multiplier**n
        - sp.conjugate(diagonal_multiplier) ** n
        - diagonal_multiplier ** (-n)
        - sp.conjugate(diagonal_multiplier) ** (-n)
    )
    limits = {
        "m2_collision_time": sp.limit(
            m**2 * diagonal_collision, m, sp.oo
        ),
        "jensen_threshold_over_m4": sp.limit(
            diagonal_jensen / m**4, m, sp.oo
        ),
        "m3_li_growth_rate": sp.limit(
            m**3 * diagonal_growth, m, sp.oo
        ),
        "li_multiplier": sp.limit(
            diagonal_multiplier, m, sp.oo
        ),
        "fixed_index_li_contribution": sp.limit(
            m**2 * diagonal_li_contribution, m, sp.oo
        ),
    }
    expected_limits = {
        "m2_collision_time": sp.Integer(2),
        "jensen_threshold_over_m4": sp.Rational(1, 4),
        "m3_li_growth_rate": sp.Integer(1),
        "li_multiplier": sp.Integer(1),
        "fixed_index_li_contribution": 2 * n**2,
    }
    if limits != expected_limits:
        issues.append(f"independent diagonal escape limits failed: {limits}")
    for cutoff in (sp.Integer(2), sp.Integer(5), sp.Integer(11)):
        if not bool(
            diagonal_jensen.subs(m, cutoff)
            > cutoff
        ):
            issues.append(
                f"diagonal Jensen growth test failed at m={cutoff}"
            )

    # Confirm every Xi-specific composition against its source artifact.
    strict = load_json(STRICT_RESULT, "strict-Laguerre source result", issues)
    boundary = load_json(BOUNDARY_RESULT, "boundary-attainment source result", issues)
    suzuki = load_json(SUZUKI_RESULT, "Suzuki source result", issues)
    strict_identity = (
        strict.get("exact", {})
        .get("correlation_identity", "")
    )
    if "L_t" not in strict_identity or "xi/2" not in strict_identity:
        issues.append("strict-Laguerre correlation identity source missing")
    attainment = (
        boundary.get("exact", {})
        .get("attainment_theorem", "")
    )
    if "multiple zero" not in attainment or "finite" not in attainment:
        issues.append("positive-boundary finite-attainment source missing")
    suzuki_exact = suzuki.get("exact", {})
    phase_diagram = suzuki_exact.get("phase_diagram", "")
    cancellation_set = suzuki_exact.get("cancellation_set", "")
    local_exclusion = suzuki_exact.get("local_exclusion", "")
    if (
        "delta_xi" not in phase_diagram
        or "countable" not in cancellation_set
        or "0<omega<epsilon_rho" not in local_exclusion
    ):
        issues.append("Suzuki fixed-shift phase diagram source missing")

    proof_boundary = stored.get("proof_boundary", "")
    forbidden_promotions = [
        "Therefore RH",
        "This proves RH",
        "This proves Lambda<=0",
        "Xi has an off-axis zero.",
    ]
    for phrase in forbidden_promotions:
        if phrase in proof_boundary:
            issues.append(f"proof boundary contains forbidden promotion: {phrase}")
    required_boundary = [
        "does not assert that Xi has an off-axis zero",
        "does not",
        "prove RH",
    ]
    for phrase in required_boundary:
        if phrase not in proof_boundary:
            issues.append(f"proof boundary missing phrase: {phrase}")

    if NOTE.exists():
        note = NOTE.read_text(encoding="utf-8")
        required_note = [
            "not a proof of `Lambda<=0` or RH",
            "Assume one off-axis quartet provisionally",
            "It is not the Xi heat flow.",
            "does not identify the",
            "first negative complete Xi Li coefficient",
            "no fixed real band, finite degree, finite Li index",
            "phase-critical-value theorem as the primary next attack",
            "Neither is proved here.",
        ]
        for phrase in required_note:
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
        "validated Newman counterfactual birth-signature atlas: "
        "12 rows, 0 issues, 1 symmetric quartet heat collision, "
        "1 exact Jensen threshold, 1 exact Li quartet law, "
        "1 Suzuki shift law, 1 diagonal finite-sensor escape theorem, "
        "2 open global handoffs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
