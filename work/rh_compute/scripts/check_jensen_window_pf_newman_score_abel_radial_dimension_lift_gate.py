#!/usr/bin/env python3
"""Independently validate the score-Abel radial dimension-lift gate."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_score_abel_radial_dimension_lift_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_score_abel_radial_dimension_lift_gate.md"
)

EXPECTED_IDS = [
    "sardl_01_imported_abel_coordinate",
    "sardl_02_radial_probability",
    "sardl_03_planar_abel_marginal",
    "sardl_04_bessel_derivative",
    "sardl_05_radial_fourier",
    "sardl_06_dimension_walk",
    "sardl_07_jensen_identification",
    "sardl_08_schoenberg_scope",
    "sardl_09_gaussian_flow_shape",
    "sardl_10_gaussian_radial_lift",
    "sardl_11_gaussian_jensen_failure",
    "sardl_12_xi_quadratic_closure",
    "sardl_13_xi_cubic_closure",
    "sardl_14_imported_degree361_closure",
    "sardl_15_nonpromotion",
]

EXPECTED_ROLES = {
    "imported_exact_bridge",
    "exact_probability_lift",
    "exact_radon_identity",
    "exact_bessel_identity",
    "exact_positive_definite_theorem",
    "exact_dimension_walk",
    "exact_window_identification",
    "exact_scope_guard",
    "exact_newman_flow_countermodel",
    "exact_radial_countermodel",
    "exact_jensen_countermodel",
    "imported_exact_quadratic_closure",
    "imported_exact_cubic_closure",
    "imported_exact_bounded_degree_closure",
    "nonpromotion_gate",
}


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = gate.build_payload()
    if stored != rebuilt:
        issues.append("stored radial dimension-lift payload differs from reconstruction")
    if stored.get("kind") != (
        "jensen_window_pf_newman_score_abel_"
        "radial_dimension_lift_gate"
    ):
        issues.append("artifact kind drifted")
    if "https://doi.org/10.2307/1968466" not in stored.get(
        "sources", []
    ):
        issues.append("Schoenberg source DOI drifted")
    if "https://arxiv.org/abs/2408.11612" not in stored.get(
        "sources", []
    ):
        issues.append("radial dimension-walk source drifted")
    for source in (
        "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
        "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md",
        "outputs/jensen_window_pf_strong_logconcave_local_quartic_countermodel.md",
        "outputs/jensen_window_pf_quartic_double_root_threshold_lemma.md",
        "outputs/jensen_window_pf_quartic_signed_hankel_branch_exclusion_lemma.md",
        "outputs/jensen_window_pf_quartic_outer_threshold_order4_nonpromotion_gate.md",
        "outputs/jensen_window_pf_quartic_outer_branch_length13_obstruction.md",
        "outputs/jensen_window_pf_quartic_outer_branch_alternate_length13_survivor_gate.md",
        "outputs/jensen_window_pf_quartic_outer_contact_length14_survivor_gate.md",
        "outputs/jensen_window_pf_newman_zero_slab_degree71_sector_certificate.md",
        "outputs/jensen_window_pf_newman_real_zero_band_degree361_sector_certificate.md",
    ):
        if source not in stored.get("sources", []):
            issues.append(f"low-degree reconciliation source missing: {source}")

    rows = stored.get("rows", [])
    if len(rows) != 15:
        issues.append("expected 15 gate rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    roles = {row.get("role") for row in rows}
    if roles != EXPECTED_ROLES:
        issues.append(f"role set drifted: {sorted(roles)}")
    if any(
        sum(row.get("role") == role for row in rows) != 1
        for role in EXPECTED_ROLES
    ):
        issues.append("each expected role must occur exactly once")

    bounded_degree_row = next(
        (
            row
            for row in rows
            if row.get("role")
            == "imported_exact_bounded_degree_closure"
        ),
        {},
    )
    if bounded_degree_row.get("readiness") != "ready_to_apply":
        issues.append("Xi degree-361 closure is not ready to apply")
    if any(
        row.get("readiness") != "guard_validated"
        for row in rows
        if "countermodel" in row.get("role", "")
    ):
        issues.append("one or more Gaussian countermodel guards are inactive")

    # Recompute the derivative/0F1 coefficients without using gate helpers.
    for shift in range(10):
        for power in range(11):
            derivative = 1 / (
                sp.factorial(shift + power)
                * sp.factorial(power)
            )
            hyper = 1 / (
                sp.factorial(shift)
                * sp.rf(shift + 1, power)
                * sp.factorial(power)
            )
            if sp.simplify(derivative - hyper) != 0:
                issues.append(
                    "independent derivative/0F1 coefficient failed "
                    f"at n={shift}, k={power}"
                )

    # Verify polar normalization over a wider dimension range.
    for shift in range(13):
        sphere = 2 * sp.pi ** (shift + 1) / sp.factorial(shift)
        jacobian = 2 ** (2 * shift + 1)
        coefficient = sp.simplify(
            sphere
            * jacobian
            / (4 * sp.pi) ** (shift + 1)
        )
        if coefficient != 1 / sp.factorial(shift):
            issues.append(
                f"independent radial normalization failed at n={shift}"
            )

    inner_abel = sp.expand_func(
        sp.beta(sp.Rational(1, 2), sp.Rational(1, 2))
    ) / 2
    if sp.simplify(inner_abel - sp.pi / 2) != 0:
        issues.append("independent planar Abel marginal integral failed")

    # Match the radial characteristic series to F^(n)/A_n.
    for shift in range(9):
        for power in range(9):
            a_n, a_nk = sp.symbols(
                f"B_{shift} B_{shift}_{power}", positive=True
            )
            tilted_moment = (
                sp.factorial(shift + power)
                * a_nk
                / (sp.factorial(shift) * a_n)
            )
            radial = sp.simplify(
                tilted_moment
                / (
                    sp.rf(shift + 1, power)
                    * sp.factorial(power)
                )
            )
            expected = a_nk / (a_n * sp.factorial(power))
            if sp.simplify(radial - expected) != 0:
                issues.append(
                    "independent radial characteristic failed "
                    f"at n={shift}, k={power}"
                )

    x = sp.symbols("x", nonzero=True, real=True)
    a_n, a_np1, f_np1 = sp.symbols(
        "C_n C_np1 G_np1", positive=True
    )
    phi_prime = -2 * x * f_np1 / a_n
    walked = -a_n * phi_prime / (2 * x * a_np1)
    if sp.simplify(walked - f_np1 / a_np1) != 0:
        issues.append("independent dimension walk failed")

    # Check the Gaussian heat-flow identity symbolically.
    u, t, c = sp.symbols("u t c", real=True)
    gaussian = (
        sp.exp(-(c**2 - t) * u**2)
        + sp.exp(-(2 * c**2 - t) * u**2)
    )
    if sp.simplify(sp.diff(gaussian, t) - u**2 * gaussian) != 0:
        issues.append("two-Gaussian Newman flow identity failed")

    # Check the shape, Abel moments, zero lattice, and strict D=2 failure
    # on several exact parameter pairs.
    parameter_pairs = [
        (sp.Integer(1), sp.Integer(0)),
        (sp.Integer(1), sp.Rational(1, 5)),
        (sp.Rational(3, 2), sp.Rational(1, 10)),
        (sp.Integer(2), sp.Rational(1, 5)),
    ]
    for c_value, t_value in parameter_pairs:
        a = c_value**2 - t_value
        b = 2 * c_value**2 - t_value
        if not (a > 0 and b > a):
            issues.append(
                f"invalid Gaussian parameters c={c_value}, t={t_value}"
            )
            continue

        curvature = (
            (2 - 4 / sp.E) * c_value**2 - sp.Rational(2, 5)
        )
        if not bool(sp.N(curvature, 40) > 0):
            issues.append(
                "Gaussian strong-log-concavity bound failed "
                f"at c={c_value}, t={t_value}"
            )

        alpha = 1 / (4 * a)
        beta = 1 / (4 * b)
        p = sp.sqrt(sp.pi / a)
        q = sp.sqrt(sp.pi / b)
        delta = sp.simplify(alpha - beta)
        if sp.simplify(delta - c_value**2 / (4 * a * b)) != 0:
            issues.append(
                f"Gaussian exponential gap failed at c={c_value}, t={t_value}"
            )

        for lattice_index in (-2, -1, 0, 1, 2):
            zero = (
                2 * a * b * sp.log(a / b) / c_value**2
                + 4
                * sp.pi
                * sp.I
                * a
                * b
                * (2 * lattice_index + 1)
                / c_value**2
            )
            target_log = (
                sp.log(a / b) / 2
                + sp.I * sp.pi * (2 * lattice_index + 1)
            )
            if sp.simplify(delta * zero - target_log) != 0:
                issues.append(
                    "Gaussian zero-lattice equation failed "
                    f"at c={c_value}, t={t_value}, k={lattice_index}"
                )
            if sp.im(zero) == 0:
                issues.append(
                    "Gaussian zero unexpectedly real "
                    f"at c={c_value}, t={t_value}, k={lattice_index}"
                )

        for power in range(11):
            abel_moment = sp.simplify(
                4
                * sp.sqrt(sp.pi * a)
                * sp.factorial(power)
                / (4 * a) ** (power + 1)
                + 4
                * sp.sqrt(sp.pi * b)
                * sp.factorial(power)
                / (4 * b) ** (power + 1)
            )
            coefficient_moment = sp.factorial(power) * (
                p * alpha**power + q * beta**power
            )
            if sp.simplify(abel_moment - coefficient_moment) != 0:
                issues.append(
                    "Gaussian Abel moment failed "
                    f"at c={c_value}, t={t_value}, k={power}"
                )

        for shift in range(11):
            a0 = p * alpha**shift + q * beta**shift
            a1 = (
                p * alpha ** (shift + 1)
                + q * beta ** (shift + 1)
            )
            a2 = (
                p * alpha ** (shift + 2)
                + q * beta ** (shift + 2)
            )
            discriminant = sp.factor(4 * (a1**2 - a0 * a2))
            expected = sp.factor(
                -4
                * p
                * q
                * (alpha * beta) ** shift
                * (alpha - beta) ** 2
            )
            if sp.simplify(discriminant - expected) != 0:
                issues.append(
                    "Gaussian quadratic Jensen identity failed "
                    f"at c={c_value}, t={t_value}, n={shift}"
                )
            if not bool(sp.N(discriminant, 40) < 0):
                issues.append(
                    "Gaussian quadratic discriminant is not negative "
                    f"at c={c_value}, t={t_value}, n={shift}"
                )

    exact = stored.get("exact", {})
    checks = exact.get("checks", {})
    expected_check_lengths = {
        "derivative_0f1_coefficients": 56,
        "radial_normalizations": 9,
        "radial_characteristic_coefficients": 49,
        "gaussian_abel_moments": 9,
        "gaussian_quadratic_discriminants": 9,
        "xi_quadratic_abel_scalings": 12,
    }
    for key, expected_length in expected_check_lengths.items():
        if len(checks.get(key, [])) != expected_length:
            issues.append(
                f"stored check count drifted for {key}: "
                f"{len(checks.get(key, []))}"
            )

    boundary = stored.get("proof_boundary", "")
    for phrase in (
        "do not imply a quadratic Jensen inequality",
        "closes every shifted window through degree 361",
        "does not prove degree 362",
        "unbounded cofinal degree sequence",
        "does not prove",
        "PF-infinity",
        "Lambda<=0",
        "RH",
    ):
        if phrase not in boundary:
            issues.append(f"proof boundary missing phrase: {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    required_note_phrases = [
        "Radial Probability Ladder",
        "planar density whose one-coordinate marginal",
        "Bessel-Fourier Identity",
        "dimension `2n+2`",
        "Schoenberg Scope Guard",
        "fixed-profile all-dimension hypothesis is absent",
        "Two-Gaussian Newman Countermodel",
        "arbitrarily large certified strong-concavity",
        "fails the degree-two Jensen test for every shift",
        "Quadratic And Cubic Reconciliation",
        "strictly log-convex, not log-concave",
        "every shifted degree-two Jensen polynomial is hyperbolic",
        "every shifted degree-three Jensen polynomial is hyperbolic",
        "Degree-361 Reconciliation",
        "strong-log-concave local Mellin countermodel",
        "outer-root quartic threshold u<=U(a,p)",
        "all 120, 126, and 56 supported signed-Hankel conditions",
        "all 364, 715, and 792 finite signed minors",
        "4,043 supported signed minors",
        "x_14 compatibility is negative, but only for that selected tail",
        "nonhyperbolic adjacent quintic",
        "through degree 361",
        "finite-degree frontier begins at degree 362",
        "https://doi.org/10.2307/1968466",
        "https://arxiv.org/abs/2408.11612",
        "not a proof of PF-infinity, RH, or `Lambda <= 0`",
    ]
    for phrase in required_note_phrases:
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
        "validated score-Abel radial dimension-lift gate: "
        "15 rows, 0 issues, 1 imported Abel coordinate, "
        "1 radial probability ladder, 1 planar marginal, "
        "1 Bessel derivative identity, "
        "1 own-dimension positive-definiteness theorem, "
        "1 dimension walk, 1 Jensen identification, "
        "1 Schoenberg scope guard, 3 exact Gaussian countermodel rows, "
        "2 imported Xi low-degree closures, "
        "1 imported degree-361 closure, 1 nonpromotion gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
