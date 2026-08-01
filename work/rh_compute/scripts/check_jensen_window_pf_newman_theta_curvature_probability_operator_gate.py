#!/usr/bin/env python3
"""Check the theta-curvature probability/operator gate independently."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARTIFACT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_curvature_probability_operator_gate.json"
)

EXPECTED_IDS = [
    "ntcpo_01_theta_primitive",
    "ntcpo_02_summand_monotonicity",
    "ntcpo_03_summand_convexity",
    "ntcpo_04_probability_law",
    "ntcpo_04a_fixed_theta_weights",
    "ntcpo_04b_dominant_component_budget",
    "ntcpo_05_triangular_mixture",
    "ntcpo_06_primitive_fourier_positivity",
    "ntcpo_07_oscillator_factorization",
    "ntcpo_08_positive_normalizer",
    "ntcpo_08a_doob_evolution",
    "ntcpo_08b_sturm_nodal_guard",
    "ntcpo_08c_drift_curvature_guard",
    "ntcpo_09_characteristic_transfer",
    "ntcpo_09a_componentwise_contact",
    "ntcpo_10_contact_reduction",
    "ntcpo_10a_dominant_block_disjunction",
    "ntcpo_11_endpoint_laplace_reference",
    "ntcpo_12_generic_double_contact_guard",
    "ntcpo_13_open_xi_transversality",
]

SUCCESS = (
    "validated Newman theta-curvature probability/operator gate: "
    "20 rows, 0 issues, 3 theta-primitive identities, "
    "2 monotone-convex inequalities, 1 probability law, "
    "1 fixed-weight theta mixture, 1 uniform C2 dominant-block budget, "
    "1 dominant-mass endpoint-jet guard, "
    "2 transform identities, 1 oscillator factorization, "
    "1 Doob diffusion identity, 1 Sturm nodal-loss guard, "
    "1 asymptotic drift-curvature guard, "
    "1 characteristic contact reduction, "
    "1 componentwise contact decomposition, 1 explicit C1 tail disjunction, "
    "1 endpoint Laplace comparison, "
    "1 generic double-contact guard, 1 open Xi transversality handoff"
)


def symbolic_audit() -> list[str]:
    issues: list[str] = []
    x, t = sp.symbols("x t", nonzero=True)
    a_fun = sp.Function("A")(x)
    c_fun = (1 - a_fun) / (2 * x**2)
    d_c = (
        -4 * t**2 * sp.diff(c_fun, x, 2)
        + 4 * t * x * sp.diff(c_fun, x)
        + (2 * t - 1 - x**2) * c_fun
    )

    def q(f):
        return 2 * t * sp.diff(f, x) - x * f

    if sp.simplify(d_c + q(q(c_fun)) + c_fun) != 0:
        issues.append("oscillator factorization residual is nonzero")

    transfer = sp.expand(
        sp.simplify(16 * x**4 * (sp.Rational(1, 16) + d_c / 8))
    )
    expected = (
        (x**4 + (6 * t + 1) * x**2 + 24 * t**2) * a_fun
        - 4 * t * x * (x**2 + 4 * t) * sp.diff(a_fun, x)
        + 4 * t**2 * x**2 * sp.diff(a_fun, x, 2)
        - ((6 * t + 1) * x**2 + 24 * t**2)
    )
    if sp.simplify(transfer - expected) != 0:
        issues.append("characteristic transfer residual is nonzero")

    w = sp.symbols("w", real=True)
    component_c = (w - a_fun) / (2 * x**2)
    component_d = (
        -4 * t**2 * sp.diff(component_c, x, 2)
        + 4 * t * x * sp.diff(component_c, x)
        + (2 * t - 1 - x**2) * component_c
    )
    component_transfer = sp.expand(
        sp.simplify(16 * x**4 * (w / 16 + component_d / 8))
    )
    p_poly = x**4 + (6 * t + 1) * x**2 + 24 * t**2
    r_poly = (6 * t + 1) * x**2 + 24 * t**2
    component_expected = (
        4 * t**2 * x**2 * sp.diff(a_fun, x, 2)
        - 4 * t * x * (x**2 + 4 * t) * sp.diff(a_fun, x)
        + p_poly * a_fun
        - r_poly * w
    )
    if sp.simplify(component_transfer - component_expected) != 0:
        issues.append("componentwise characteristic transfer residual is nonzero")
    component_derivative = (
        4 * t**2 * x**2 * sp.diff(a_fun, x, 3)
        - 4 * t * x * (x**2 + 2 * t) * sp.diff(a_fun, x, 2)
        + (p_poly - 4 * t * (3 * x**2 + 4 * t))
        * sp.diff(a_fun, x)
        + (4 * x**3 + 2 * (6 * t + 1) * x) * a_fun
        - 2 * (6 * t + 1) * x * w
    )
    if sp.simplify(sp.diff(component_expected, x) - component_derivative) != 0:
        issues.append("componentwise transfer derivative residual is nonzero")
    if not (
        sp.Rational(6, 49**3 * 2800)
        < sp.Rational(1, 54_900_000)
    ):
        issues.append("third-moment theta-tail certificate failed")

    endpoint = x**2 * ((1 + x**2) * a_fun - 1)
    if sp.simplify(expected.subs(t, 0) - endpoint) != 0:
        issues.append("endpoint Laplace comparison residual is nonzero")

    c0, c1, c2, g0, g1, g2 = sp.symbols("c0 c1 c2 g0 g1 g2")
    product = sp.expand(
        (c1 * g0 + c0 * g1) ** 2
        - c0 * g0 * (c2 * g0 + 2 * c1 * g1 + c0 * g2)
    )
    target = c0**2 * (g1**2 - g0 * g2) + g0**2 * (
        c1**2 - c0 * c2
    )
    if sp.simplify(product - target) != 0:
        issues.append("product Laguerre residual is nonzero")

    c_ratio = sp.Function("C")(x)
    h_ratio = sp.Function("H")(x)
    g_ratio = 8 * h_ratio / c_ratio
    ratio_time = 8 * (
        -sp.diff(h_ratio, x, 2) * c_ratio
        + h_ratio * sp.diff(c_ratio, x, 2)
    ) / c_ratio**2
    ratio_residual = sp.simplify(
        ratio_time
        + sp.diff(g_ratio, x, 2)
        + 2
        * sp.diff(c_ratio, x)
        / c_ratio
        * sp.diff(g_ratio, x)
    )
    if ratio_residual != 0:
        issues.append("positive-normalizer Doob residual is nonzero")
    divergence_residual = sp.simplify(
        sp.diff(g_ratio, x, 2)
        + 2
        * sp.diff(c_ratio, x)
        / c_ratio
        * sp.diff(g_ratio, x)
        - sp.diff(c_ratio**2 * sp.diff(g_ratio, x), x) / c_ratio**2
    )
    if divergence_residual != 0:
        issues.append("positive-normalizer divergence residual is nonzero")

    a_guard = sp.symbols("a_guard", real=True)
    h_guard = (x**2 + a_guard - 2 * t) / 8
    g_guard = 8 * h_guard
    if sp.simplify(sp.diff(h_guard, t) + sp.diff(h_guard, x, 2)) != 0:
        issues.append("Sturm guard heat residual is nonzero")
    if sp.simplify(sp.diff(g_guard, t) + sp.diff(g_guard, x, 2)) != 0:
        issues.append("Sturm guard ratio residual is nonzero")

    r1, r3 = sp.symbols("r1 r3", real=True)
    third_jet = (r3 + 6 * t * r1).subs(
        {r1: -sp.Rational(1, 2), r3: -sp.Rational(1, 2)}
    )
    if sp.simplify(third_jet + sp.Rational(1, 2) + 3 * t) != 0:
        issues.append("primitive third-jet transfer failed")
    free = sp.symbols("free", real=True)
    c_asymptotic = (
        sp.Rational(1, 2) / x**2
        - (sp.Rational(1, 2) + 3 * t) / x**4
        + free / x**6
    )
    log_curvature_lead = sp.simplify(
        sp.limit(
            x**4
            * (
                sp.diff(sp.log(c_asymptotic), x, 2)
                - 2 / x**2
            ),
            x,
            sp.oo,
        )
    )
    if log_curvature_lead != -6 * (1 + 6 * t):
        issues.append("primitive-transform log-curvature asymptotic failed")

    y = sp.symbols("y", positive=True)
    margin = sp.Rational(121, 144) * y**2 - 4 * y
    if sp.simplify(margin.subs(y, 12)) <= 0:
        issues.append("convexity floor at y=12 is not positive")

    exp15_lower = sum(
        sp.Rational(15) ** k / sp.factorial(k) for k in range(31)
    )
    if exp15_lower <= 2_250_000:
        issues.append("theta-tail ratio Taylor certificate failed")
    exp4pi_lower = sum(
        sp.Rational(666, 53) ** k / sp.factorial(k) for k in range(41)
    )
    if exp4pi_lower <= sp.Rational(32 * 22 * 2801, 7):
        issues.append("first omitted theta-weight certificate failed")
    if not (
        sp.Rational(1_000_000, 2801 * 999_999)
        < sp.Rational(1, 2800)
    ):
        issues.append("theta-weight geometric sum failed")

    z, b = sp.symbols("z b", real=True)
    fifth = (
        -1024 * z**5
        + 11520 * z**4
        - 33920 * z**3
        + 26400 * z**2
        - 3124 * z
        + 1
    )
    monotonicity = sp.diff(fifth, z) - fifth
    lo = sp.Rational(333, 106)
    hi = sp.Rational(22, 7)
    unit_poly = sp.Poly(
        sp.expand(monotonicity.subs(z, lo + (hi - lo) * b)), b
    )
    degree = unit_poly.degree()
    powers = [unit_poly.nth(i) for i in range(degree + 1)]
    bernstein = [
        sp.factor(
            sum(
                powers[i]
                * sp.binomial(k, i)
                / sp.binomial(degree, i)
                for i in range(k + 1)
            )
        )
        for k in range(degree + 1)
    ]
    if any(value <= 0 for value in bernstein):
        issues.append("fifth-jet monotonicity certificate failed")
    if fifth.subs(z, lo) <= 7200:
        issues.append("fifth-jet polynomial floor failed")
    exp_upper_partial = sum(hi**k / sp.factorial(k) for k in range(19))
    exp_upper_term = hi**19 / sp.factorial(19)
    exp_upper = sp.simplify(
        exp_upper_partial + exp_upper_term / (1 - hi / 20)
    )
    if exp_upper >= 24:
        issues.append("exp(pi) upper certificate failed")

    xi = sp.symbols("xi", real=True)
    transform = (
        8
        * sp.sin(xi / 2) ** 2
        / xi**2
        * sp.sqrt(2 * sp.pi)
        * sp.exp(-2 * xi**2)
    )
    point = 2 * sp.pi
    if sp.simplify(transform.subs(xi, point)) != 0:
        issues.append("guard transform does not vanish at 2*pi")
    if sp.simplify(sp.diff(transform, xi).subs(xi, point)) != 0:
        issues.append("guard transform slope does not vanish at 2*pi")
    if sp.simplify(sp.diff(transform, xi, 2).subs(xi, point)) == 0:
        issues.append("guard transform lacks a nonzero second derivative")

    u, v = sp.symbols("u v", positive=True)
    green_left = sp.exp(-(v - u)) + sp.exp(-(u + v))
    green_right = sp.exp(-(u - v)) + sp.exp(-(u + v))
    if sp.simplify(sp.diff(green_left, u).subs(u, 0)) != 0:
        issues.append("Neumann Green kernel has nonzero endpoint slope")
    jump = sp.simplify(
        sp.diff(green_right, u).subs(u, v)
        - sp.diff(green_left, u).subs(u, v)
    )
    if jump != -2:
        issues.append("Neumann Green kernel derivative jump is not -2")
    return issues


def content_audit(payload: dict) -> list[str]:
    issues: list[str] = []
    if payload.get("kind") != (
        "jensen_window_pf_newman_theta_curvature_probability_operator_gate"
    ):
        issues.append("unexpected artifact kind")
    rows = payload.get("rows")
    if not isinstance(rows, list) or len(rows) != 20:
        issues.append("expected exactly 20 rows")
        return issues
    ids = [row.get("id") for row in rows]
    if ids != EXPECTED_IDS:
        issues.append("row ids or deterministic order changed")
    if len(set(ids)) != len(ids):
        issues.append("duplicate row id")

    exact = payload.get("exact", {})
    required = [
        "theta_primitive",
        "monotone_convex_proof",
        "probability_law",
        "fixed_theta_mixture",
        "dominant_component_c2_budget",
        "positive_transform",
        "operator_factorization",
        "normalizer_evolution",
        "sturm_nodal_guard",
        "normalizer_weight_asymptotic",
        "characteristic_transfer",
        "componentwise_contact",
        "dominant_block_disjunction",
        "endpoint_laplace_comparison",
        "laguerre_product",
        "generic_guard",
        "open_handoff",
        "checks",
    ]
    for key in required:
        if key not in exact:
            issues.append(f"missing exact section: {key}")

    primitive = exact.get("theta_primitive", {})
    for phrase in ("R''(u)-R(u)=8*Phi(u)", "R'(0)=-1/2"):
        if phrase not in primitive.values():
            issues.append(f"missing primitive identity: {phrase}")

    law = exact.get("probability_law", {})
    if law.get("definition") != "dmu_t(v)=2*S_t''(v)dv on [0,infinity)":
        issues.append("probability density formula changed")
    if not str(law.get("mass", "")).endswith("=1"):
        issues.append("probability mass identity missing")

    mixture = exact.get("fixed_theta_mixture", {})
    if "2*(4*pi*n^2-1)*exp(-pi*n^2)" not in mixture.get(
        "component_mass", ""
    ):
        issues.append("fixed theta-component weight missing")
    if "w_1>2799/2800" not in mixture.get("dominant_weight", ""):
        issues.append("dominant first-component weight missing")

    budget = exact.get("dominant_component_c2_budget", {})
    for phrase in (
        "1/2800",
        "1/137200",
        "1/3361400",
        "1/54900000",
    ):
        if phrase not in " ".join(str(value) for value in budget.values()):
            issues.append(f"dominant-component budget missing: {phrase}")
    for phrase in (
        "R^((2k+1))(0)=R'(0)=-1/2",
        ">300",
        "<-601/2",
    ):
        if phrase not in " ".join(str(value) for value in budget.values()):
            issues.append(f"endpoint-jet guard missing: {phrase}")
    if "modularly coupled oscillatory decay" not in budget.get("scope", ""):
        issues.append("global-frequency scope guard missing")

    transform = exact.get("positive_transform", {})
    if "C_t(x)>0 for every real x" not in transform.get(
        "strict_positivity", ""
    ):
        issues.append("strict primitive-transform positivity missing")
    if "not for the first correlation K_(1,t)" not in transform.get(
        "scope", ""
    ):
        issues.append("first-correlation scope guard missing")

    operator = exact.get("operator_factorization", {})
    if operator.get("factorization") != "D_t=-(Q_t^2+1)":
        issues.append("operator factorization changed")
    if "C_t(x)>0" not in operator.get("positive_normalizer", ""):
        issues.append("positive normalizer statement missing")

    evolution = exact.get("normalizer_evolution", {})
    evolution_text = " ".join(str(value) for value in evolution.values())
    for phrase in (
        "partial_t G_t=-partial_x^2 G_t",
        "-C_t^(-2)*partial_x(C_t^2*partial_x G_t)",
        "partial_tau G=",
        "<f,A_t f>_(C_t^2)=integral C_t^2*(f')^2 dx",
        "Sturm zero number may drop",
    ):
        if phrase not in evolution_text:
            issues.append(f"missing normalizer-evolution phrase: {phrase}")

    sturm = exact.get("sturm_nodal_guard", {})
    sturm_text = " ".join(str(value) for value in sturm.values())
    for phrase in (
        "C_t(x)=1",
        "G_t(x)=x^2+a-2t",
        "two simple real zeros",
        "positive heat semigroup loses two real nodes",
        "do not exclude contact",
        "Xi-specific conserved nodal count",
    ):
        if phrase not in sturm_text:
            issues.append(f"missing Sturm nodal guard phrase: {phrase}")

    weight_asymptotic = exact.get("normalizer_weight_asymptotic", {})
    weight_text = " ".join(
        str(value) for value in weight_asymptotic.values()
    )
    for phrase in (
        "S_t'''(0)=R'''(0)+6*t*R'(0)=-1/2-3*t",
        "C_t(x)=1/(2*x^2)-(1/2+3*t)/x^4+O_t(x^-6)",
        "(log C_t)''=2/x^2-6*(1+6*t)/x^4+O(x^-6)>0",
        "uniformly on 0<=t<=1/5",
        "eventually log-convex",
        "Bakry-Emery positive-curvature",
    ):
        if phrase not in weight_text:
            issues.append(f"missing drift-curvature guard phrase: {phrase}")

    contact = exact.get("characteristic_transfer", {})
    for phrase in ("J_t=4*t^2*x^2*A_t''", "H_t(c)=H_t'(c)=0"):
        if phrase not in " ".join(str(value) for value in contact.values()):
            issues.append(f"missing contact phrase: {phrase}")

    component = exact.get("componentwise_contact", {})
    component_text = " ".join(str(value) for value in component.values())
    for phrase in (
        "J_t=sum_(n>=1)J_(n,t)",
        "J_(n,t)=4*t^2*x^2*A_(n,t)''",
        "J_(n,t)'=4*t^2*x^2*A_(n,t)'''",
        "Gamma((1+i*x)/4,pi)",
    ):
        if phrase not in component_text:
            issues.append(f"missing componentwise contact phrase: {phrase}")

    disjunction = exact.get("dominant_block_disjunction", {})
    disjunction_text = " ".join(
        str(value) for value in disjunction.values()
    )
    for phrase in (
        "delta_3=1/54900000",
        "B_0(t,x)=",
        "B_1(t,x)=",
        "|J_t-J_(1,t)|<B_0(t,x)",
        "|J_(1,t)(x)|>B_0(t,x) or",
        "(J_t(x),J_t'(x))!=(0,0)",
        "boundary-contact criterion gives Lambda<=0",
        "proved pointwise exclusion test",
        "modular endpoint cancellation",
    ):
        if phrase not in disjunction_text:
            issues.append(f"missing dominant-block disjunction phrase: {phrase}")
    disjunction_rows = [
        row
        for row in rows
        if row.get("id") == "ntcpo_10a_dominant_block_disjunction"
    ]
    if (
        len(disjunction_rows) != 1
        or disjunction_rows[0].get("readiness") != "ready_to_apply"
    ):
        issues.append("dominant-block pointwise exclusion is not ready to apply")

    endpoint = exact.get("endpoint_laplace_comparison", {})
    if endpoint.get("formula") != "16*x^2*H_0(x)=(1+x^2)*A_0(x)-1":
        issues.append("endpoint transfer formula changed")
    if "unit exponential law" not in endpoint.get("reference", ""):
        issues.append("Laplace reference missing")

    guard = exact.get("generic_guard", {})
    for phrase in (
        "g''-g=8*f",
        "multiple real zero at 2*pi",
        "do not alone exclude contact",
    ):
        if phrase not in " ".join(str(value) for value in guard.values()):
            issues.append(f"missing generic guard phrase: {phrase}")

    if rows[-1].get("readiness") != "open":
        issues.append("Xi transversality handoff is not marked open")
    boundary = payload.get("proof_boundary", "")
    for phrase in ("does not prove", "Lambda<=0", "RH"):
        if phrase not in boundary:
            issues.append(f"proof boundary missing: {phrase}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_ARTIFACT)
    args = parser.parse_args()
    payload = json.loads(args.artifact.read_text(encoding="utf-8"))
    issues = symbolic_audit() + content_audit(payload)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        print(f"theta-curvature probability/operator gate failed: {len(issues)} issues")
        return 1
    print(SUCCESS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
