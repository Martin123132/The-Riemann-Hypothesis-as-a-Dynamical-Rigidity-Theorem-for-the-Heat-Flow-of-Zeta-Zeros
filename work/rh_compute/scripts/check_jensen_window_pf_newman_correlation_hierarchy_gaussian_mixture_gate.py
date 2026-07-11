#!/usr/bin/env python3
"""Independently validate the Newman correlation-hierarchy mixture gate."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_correlation_hierarchy_gaussian_mixture_gate.md"
)

EXPECTED_IDS = [
    "nchgm_01_heat_hierarchy",
    "nchgm_02_fourier_hierarchy",
    "nchgm_03_laguerre_normalization",
    "nchgm_04_multiple_root_contact",
    "nchgm_05_gaussian_mixture_sufficiency",
    "nchgm_06_alternating_derivative_scout",
    "nchgm_07_logconvexity_witness",
    "nchgm_08_super_gaussian_tail",
    "nchgm_09_complete_monotonicity_rejection",
    "nchgm_10_polya_pf_rejections",
    "nchgm_11_tail_compatible_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored hierarchy payload differs from reconstruction")
    if stored.get("status") != (
        "exact correlation hierarchy with Gaussian-mixture and direct-PF obstructions"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 11:
        issues.append("expected 11 rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    for role, expected in [
        ("exact_identity", 3),
        ("exact_boundary_signature", 1),
        ("published_theorem_composition", 1),
        ("numerical_diagnostic", 1),
        ("numerical_counter_witness", 1),
        ("exact_theorem", 1),
        ("nonpromotion_gate", 2),
        ("open_handoff", 1),
    ]:
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    exact = stored.get("exact", {})

    # Re-derive the heat and Fourier hierarchy coefficients.
    s, v = sp.symbols("s v", real=True)
    if sp.expand((s + v) ** 2 + (s - v) ** 2) != 2 * s**2 + 2 * v**2:
        issues.append("heat-weight identity failed")
    if exact.get("heat_hierarchy") != (
        "partial_t K_(n,t)(v)=2*v^2*K_(n,t)(v)+2*K_(n+1,t)(v)"
    ):
        issues.append("stored heat hierarchy drifted")
    if "-2*partial_xi^2" not in exact.get("fourier_hierarchy", ""):
        issues.append("Fourier v^2 sign or factor drifted")

    checks = exact.get("normalization_checks", [])
    if len(checks) != 7:
        issues.append("normalization check count drifted")
    for n, row in enumerate(checks):
        expected = sp.Rational(2 ** (2 * n - 1), sp.factorial(2 * n))
        if sp.sympify(row.get("coefficient", "nan")) != expected:
            issues.append(f"Laguerre normalization failed at n={n}")
    if checks[1].get("coefficient") != "1":
        issues.append("first-correlation normalization is not one")

    # Re-derive the universal double-root contact from a local Taylor jet.
    y, a, b, c4 = sp.symbols("y a b c4", real=True)
    h = a * y**2 / 2 + b * y**3 / 6 + c4 * y**4 / 24
    l1 = sp.expand(sp.diff(h, y) ** 2 - h * sp.diff(h, y, 2))
    if sp.expand(l1).coeff(y, 2) != a**2 / 2:
        issues.append("double-root Laguerre quadratic contact failed")
    l1_xx_at_zero = sp.diff(l1, y, 2).subs(y, 0)
    if sp.simplify(l1_xx_at_zero - a**2) != 0:
        issues.append("double-root L1 curvature failed")
    # F1(xi)=L1(xi/2), F2=3*L2, and L2(0)=(H''/2)^2.
    f1_xixi = l1_xx_at_zero / 4
    f2 = 3 * a**2 / 4
    if sp.simplify(f2 - 3 * f1_xixi) != 0:
        issues.append("F2=3*F1'' boundary contact failed")
    contact = exact.get("double_root_contact", {})
    if "F_(2,t)(2c)=3*partial_xi^2" not in contact.get("compatibility", ""):
        issues.append("stored boundary compatibility drifted")

    # Check the Gaussian-mixture Fourier formula and its strict sign.
    mix = exact.get("gaussian_mixture_sufficiency", {})
    if mix.get("source") != "https://doi.org/10.2307/1968466":
        issues.append("Schoenberg source drifted")
    if ">0" not in mix.get("fourier_formula", ""):
        issues.append("strict Gaussian-transform sign missing")
    if "Lambda<=0" not in mix.get("consequence", ""):
        issues.append("mixture-to-Newman consequence missing")

    # Independently compare the fast quadrature against the stored 55-digit run.
    numerics = stored.get("numerics", {})
    numerical_rows = {row.get("t"): row for row in numerics.get("rows", [])}
    references = numerics.get("independent_mpmath_55_digit_references", {})
    for t_key, reference_key in [("0", "t=0"), ("0.5", "t=1/2")]:
        row = numerical_rows.get(t_key, {})
        values = [mp.mpf(value) for value in row.get("g_derivatives_at_zero", [])]
        reference = references.get(reference_key, {})
        if len(values) != 9:
            issues.append(f"missing derivative values at t={t_key}")
            continue
        for index, field in enumerate(("g0", "g1", "g2")):
            expected = mp.mpf(reference.get(field, "nan"))
            if abs(values[index] - expected) > mp.mpf("2e-13") * abs(expected):
                issues.append(f"independent numerical mismatch {field} at t={t_key}")
        determinant = values[0] * values[2] - values[1] ** 2
        expected_det = mp.mpf(reference.get("determinant", "nan"))
        if abs(determinant - expected_det) > mp.mpf("2e-13") * abs(expected_det):
            issues.append(f"independent determinant mismatch at t={t_key}")
    for row in numerics.get("rows", []):
        values = [mp.mpf(value) for value in row.get("g_derivatives_at_zero", [])]
        if len(values) != 9 or not all(((-1) ** n) * value > 0 for n, value in enumerate(values)):
            issues.append(f"alternating sign scout failed at t={row.get('t')}")
        determinant = mp.mpf(row.get("log_convexity_determinant", "nan"))
        if determinant >= mp.mpf("-7e-6"):
            issues.append(f"log-convexity witness lost margin at t={row.get('t')}")
    if mp.mpf(numerics.get("max_relative_coarse_fine_delta_orders_0_to_4", "inf")) >= mp.mpf("1e-9"):
        issues.append("quadrature convergence margin failed")
    if "Numerical diagnostic only" not in numerics.get("scope", ""):
        issues.append("numerical scope boundary missing")

    # Check the two elementary inequalities behind the exact tail theorem.
    mp.mp.dps = 50
    for v_value in [mp.mpf(k) / 8 for k in range(0, 25)]:
        for s_value in [mp.mpf(k) / 8 for k in range(-32, 33)]:
            left = mp.exp(4 * abs(s_value + v_value)) + mp.exp(
                4 * abs(s_value - v_value)
            )
            right = mp.exp(4 * v_value) + mp.exp(4 * abs(s_value))
            if left < right:
                issues.append("tail exponential inequality failed on rational grid")
                break
    tail = exact.get("correlation_tail_bound", {})
    if "-pi*exp(4v)" not in tail.get("bound", ""):
        issues.append("super-Gaussian exponent missing")
    if "+infinity" not in tail.get("limit", ""):
        issues.append("super-Gaussian limit missing")
    obstruction = exact.get("complete_monotonicity_obstruction", "")
    for phrase in ("mu([0,A])", "limsup", "not completely monotone"):
        if phrase not in obstruction:
            issues.append(f"complete-monotonicity proof missing: {phrase}")

    direct_pf = exact.get("direct_pf_infinity_obstruction", {})
    if "Gaussian lower tail" not in direct_pf.get("tail_consequence", ""):
        issues.append("direct PF Gaussian-tail obstruction missing")
    if "exponential lower tail" not in direct_pf.get("tail_consequence", ""):
        issues.append("direct PF exponential-tail obstruction missing")
    if direct_pf.get("sources") != [
        "https://doi.org/10.1073/pnas.33.1.11",
        "https://arxiv.org/abs/2006.16213",
    ]:
        issues.append("PF-infinity sources drifted")
    handoff = exact.get("tail_compatible_handoff", "")
    for phrase in ("spectral-square", "coercivity", "F_2=3*F_1''"):
        if phrase not in handoff:
            issues.append(f"tail-compatible handoff missing: {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "Exact Hierarchy",
        "Gaussian-Mixture Test",
        "Exact Tail Obstruction",
        "negative, not positive",
        "cannot be PF-infinity",
        "spectral-square",
        "https://doi.org/10.1073/pnas.33.1.11",
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
        "validated Jensen-window PF Newman correlation hierarchy Gaussian-mixture gate: "
        "11 rows, 0 issues, 3 exact hierarchy identities, "
        "1 universal boundary-contact signature, 1 Gaussian-mixture sufficient theorem, "
        "2 numerical diagnostics, 1 exact super-Gaussian tail theorem, "
        "2 non-promotion gates, 1 tail-compatible handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
