#!/usr/bin/env python3
"""Build the theta-curvature probability/operator gate for the Newman target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_curvature_probability_operator_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_curvature_probability_operator_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    x, t = sp.symbols("x t", nonzero=True)
    a_fun = sp.Function("A")(x)
    c_fun = (1 - a_fun) / (2 * x**2)

    d_c = (
        -4 * t**2 * sp.diff(c_fun, x, 2)
        + 4 * t * x * sp.diff(c_fun, x)
        + (2 * t - 1 - x**2) * c_fun
    )
    q = lambda f: 2 * t * sp.diff(f, x) - x * f
    oscillator_residual = sp.simplify(d_c + q(q(c_fun)) + c_fun)
    if oscillator_residual != 0:
        raise RuntimeError("harmonic-oscillator factorization failed")

    transfer = sp.expand(
        sp.simplify(16 * x**4 * (sp.Rational(1, 16) + d_c / 8))
    )
    transfer_expected = (
        (x**4 + (6 * t + 1) * x**2 + 24 * t**2) * a_fun
        - 4 * t * x * (x**2 + 4 * t) * sp.diff(a_fun, x)
        + 4 * t**2 * x**2 * sp.diff(a_fun, x, 2)
        - ((6 * t + 1) * x**2 + 24 * t**2)
    )
    if sp.simplify(transfer - transfer_expected) != 0:
        raise RuntimeError("characteristic-function transfer failed")

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
    component_transfer_expected = (
        4 * t**2 * x**2 * sp.diff(a_fun, x, 2)
        - 4 * t * x * (x**2 + 4 * t) * sp.diff(a_fun, x)
        + p_poly * a_fun
        - r_poly * w
    )
    if sp.simplify(component_transfer - component_transfer_expected) != 0:
        raise RuntimeError("componentwise characteristic transfer failed")

    component_derivative_expected = (
        4 * t**2 * x**2 * sp.diff(a_fun, x, 3)
        - 4 * t * x * (x**2 + 2 * t) * sp.diff(a_fun, x, 2)
        + (p_poly - 4 * t * (3 * x**2 + 4 * t))
        * sp.diff(a_fun, x)
        + (4 * x**3 + 2 * (6 * t + 1) * x) * a_fun
        - 2 * (6 * t + 1) * x * w
    )
    if (
        sp.simplify(
            sp.diff(component_transfer_expected, x)
            - component_derivative_expected
        )
        != 0
    ):
        raise RuntimeError("componentwise transfer derivative failed")

    tail_c3_check = (
        sp.Rational(6, 49**3 * 2800)
        < sp.Rational(1, 54_900_000)
    )
    if tail_c3_check is not sp.true:
        raise RuntimeError("third-moment theta-tail certificate failed")

    endpoint_transfer = sp.factor(transfer_expected.subs(t, 0))
    endpoint_expected = x**2 * ((1 + x**2) * a_fun - 1)
    if sp.simplify(endpoint_transfer - endpoint_expected) != 0:
        raise RuntimeError("endpoint Laplace comparison failed")

    c0, c1, c2, g0, g1, g2 = sp.symbols("c0 c1 c2 g0 g1 g2")
    product_laguerre = sp.expand(
        (c1 * g0 + c0 * g1) ** 2
        - c0 * g0 * (c2 * g0 + 2 * c1 * g1 + c0 * g2)
    )
    product_expected = c0**2 * (g1**2 - g0 * g2) + g0**2 * (
        c1**2 - c0 * c2
    )
    if sp.simplify(product_laguerre - product_expected) != 0:
        raise RuntimeError("product Laguerre identity failed")

    c_ratio = sp.Function("C")(x)
    h_ratio = sp.Function("H")(x)
    g_ratio = 8 * h_ratio / c_ratio
    ratio_time = 8 * (
        -sp.diff(h_ratio, x, 2) * c_ratio
        + h_ratio * sp.diff(c_ratio, x, 2)
    ) / c_ratio**2
    ratio_pde_residual = sp.simplify(
        ratio_time
        + sp.diff(g_ratio, x, 2)
        + 2
        * sp.diff(c_ratio, x)
        / c_ratio
        * sp.diff(g_ratio, x)
    )
    if ratio_pde_residual != 0:
        raise RuntimeError("positive-normalizer Doob evolution failed")
    ratio_divergence_residual = sp.simplify(
        sp.diff(g_ratio, x, 2)
        + 2
        * sp.diff(c_ratio, x)
        / c_ratio
        * sp.diff(g_ratio, x)
        - sp.diff(c_ratio**2 * sp.diff(g_ratio, x), x) / c_ratio**2
    )
    if ratio_divergence_residual != 0:
        raise RuntimeError("positive-normalizer divergence form failed")

    a_guard = sp.symbols("a_guard", real=True)
    h_guard = (x**2 + a_guard - 2 * t) / 8
    c_guard = sp.Integer(1)
    g_guard = sp.simplify(8 * h_guard / c_guard)
    guard_heat_residual = sp.simplify(
        sp.diff(h_guard, t) + sp.diff(h_guard, x, 2)
    )
    guard_ratio_residual = sp.simplify(
        sp.diff(g_guard, t) + sp.diff(g_guard, x, 2)
    )
    if guard_heat_residual != 0 or guard_ratio_residual != 0:
        raise RuntimeError("Sturm nodal-loss guard failed")

    r1, r3 = sp.symbols("r1 r3", real=True)
    primitive_third_jet = sp.expand(r3 + 6 * t * r1)
    primitive_third_modular = sp.simplify(
        primitive_third_jet.subs(
            {r1: -sp.Rational(1, 2), r3: -sp.Rational(1, 2)}
        )
    )
    if primitive_third_modular != -sp.Rational(1, 2) - 3 * t:
        raise RuntimeError("primitive third-jet transfer failed")
    asymptotic_free = sp.symbols("asymptotic_free", real=True)
    c_asymptotic = (
        sp.Rational(1, 2) / x**2
        - (sp.Rational(1, 2) + 3 * t) / x**4
        + asymptotic_free / x**6
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
        raise RuntimeError("primitive-transform log-curvature asymptotic failed")

    y = sp.symbols("y", positive=True)
    convex_margin = sp.factor(sp.Rational(121, 144) * y**2 - 4 * y)
    convex_floor = sp.simplify(convex_margin.subs(y, 12))
    if convex_floor <= 0:
        raise RuntimeError("convexity floor failed")

    exp15_lower = sum(
        sp.Rational(15) ** k / sp.factorial(k) for k in range(31)
    )
    if exp15_lower <= 2_250_000:
        raise RuntimeError("geometric theta-tail ratio certificate failed")
    exp4pi_lower = sum(
        sp.Rational(666, 53) ** k / sp.factorial(k) for k in range(41)
    )
    exp4pi_target = sp.Rational(32 * 22 * 2801, 7)
    if exp4pi_lower <= exp4pi_target:
        raise RuntimeError("first omitted theta-weight certificate failed")
    geometric_tail_check = sp.simplify(
        sp.Rational(1_000_000, 2801 * 999_999)
        < sp.Rational(1, 2800)
    )
    if geometric_tail_check is not sp.true:
        raise RuntimeError("theta-weight geometric sum failed")

    z, b = sp.symbols("z b", real=True)
    fifth_jet_polynomial = (
        -1024 * z**5
        + 11520 * z**4
        - 33920 * z**3
        + 26400 * z**2
        - 3124 * z
        + 1
    )
    jet_monotonicity = sp.diff(fifth_jet_polynomial, z) - fifth_jet_polynomial
    pi_lo = sp.Rational(333, 106)
    pi_hi = sp.Rational(22, 7)
    jet_poly_on_unit = sp.Poly(
        sp.expand(jet_monotonicity.subs(z, pi_lo + (pi_hi - pi_lo) * b)),
        b,
    )
    degree = jet_poly_on_unit.degree()
    power_coefficients = [
        jet_poly_on_unit.nth(i) for i in range(degree + 1)
    ]
    jet_bernstein = [
        sp.factor(
            sum(
                power_coefficients[i]
                * sp.binomial(k, i)
                / sp.binomial(degree, i)
                for i in range(k + 1)
            )
        )
        for k in range(degree + 1)
    ]
    if any(value <= 0 for value in jet_bernstein):
        raise RuntimeError("fifth-jet monotonicity certificate failed")
    if fifth_jet_polynomial.subs(z, pi_lo) <= 7200:
        raise RuntimeError("fifth-jet polynomial floor failed")
    exp_pi_upper_partial = sum(
        pi_hi**k / sp.factorial(k) for k in range(19)
    )
    exp_pi_upper_term = pi_hi**19 / sp.factorial(19)
    exp_pi_upper = sp.simplify(
        exp_pi_upper_partial
        + exp_pi_upper_term / (1 - pi_hi / 20)
    )
    if exp_pi_upper >= 24:
        raise RuntimeError("exp(pi) rational upper certificate failed")

    xi = sp.symbols("xi", real=True)
    guard_transform = (
        8
        * sp.sin(xi / 2) ** 2
        / xi**2
        * sp.sqrt(2 * sp.pi)
        * sp.exp(-2 * xi**2)
    )
    guard_point = 2 * sp.pi
    guard_zero = sp.simplify(guard_transform.subs(xi, guard_point))
    guard_slope = sp.simplify(sp.diff(guard_transform, xi).subs(xi, guard_point))
    guard_curvature = sp.simplify(
        sp.diff(guard_transform, xi, 2).subs(xi, guard_point)
    )
    if guard_zero != 0 or guard_slope != 0 or guard_curvature == 0:
        raise RuntimeError("generic double-contact guard failed")

    u, v = sp.symbols("u v", positive=True)
    green_left = sp.exp(-(v - u)) + sp.exp(-(u + v))
    green_right = sp.exp(-(u - v)) + sp.exp(-(u + v))
    green_endpoint_slope = sp.simplify(
        sp.diff(green_left, u).subs(u, 0)
    )
    green_derivative_jump = sp.simplify(
        sp.diff(green_right, u).subs(u, v)
        - sp.diff(green_left, u).subs(u, v)
    )
    if green_endpoint_slope != 0 or green_derivative_jump != -2:
        raise RuntimeError("Neumann Green-kernel audit failed")

    return {
        "theta_primitive": {
            "definition": (
                "R(u)=sum_(n>=1) exp(u-pi*n^2*exp(4u)), u>=0"
            ),
            "differential_identity": "R''(u)-R(u)=8*Phi(u)",
            "modular_endpoint": "R'(0)=-1/2",
            "newman_primitive": "S_t(u)=exp(t*u^2)*R(u), 0<=t<=1/5",
            "summand": (
                "s_(n,t)(u)=exp(t*u^2+u-pi*n^2*exp(4u))"
            ),
        },
        "monotone_convex_proof": {
            "variables": "y=4*pi*n^2*exp(4u), y>=4*pi>12",
            "log_derivative": "s_(n,t)'/s_(n,t)=1+2*t*u-y",
            "second_derivative": (
                "s_(n,t)''/s_(n,t)=(y-1-2*t*u)^2+2*t-4*y"
            ),
            "elementary_bound": (
                "1+2*t*u<=1+(2/5)u<=y/12, hence "
                "y-1-2*t*u>=11*y/12"
            ),
            "strict_margin": (
                "s_(n,t)''/s_(n,t)>=(121/144)*y^2-4*y"
                "=y*((121/144)*y-4)>0"
            ),
            "conclusion": (
                "Every s_(n,t), and therefore S_t, is positive, strictly "
                "decreasing, and strictly convex on u>=0."
            ),
        },
        "probability_law": {
            "definition": "dmu_t(v)=2*S_t''(v)dv on [0,infinity)",
            "mass": (
                "mu_t([0,infinity))=2*(S_t'(infinity)-S_t'(0))=1"
            ),
            "triangular_mixture": (
                "S_t(u)=(1/2)*integral_[0,infinity)(v-u)_+ dmu_t(v)"
            ),
            "density_at_endpoint": (
                "S_t'(0)=R'(0)=-1/2 and S_t'(infinity)=0"
            ),
        },
        "fixed_theta_mixture": {
            "component_measure": (
                "dmu_(n,t)(v)=2*s_(n,t)''(v)dv"
            ),
            "component_mass": (
                "w_n=mu_(n,t)([0,infinity))="
                "2*(4*pi*n^2-1)*exp(-pi*n^2)"
            ),
            "time_independence": (
                "The weights w_n are independent of t and sum_n w_n=1."
            ),
            "normalized_components": (
                "dnu_(n,t)=dmu_(n,t)/w_n, "
                "A_t(x)=sum_(n>=1)w_n*a_(n,t)(x)"
            ),
            "dominant_weight": (
                "sum_(n>=2)w_n<1/2800, hence w_1>2799/2800"
            ),
            "exact_first_weight": (
                "w_1=2*(4*pi-1)*exp(-pi)="
                + str(sp.N(2 * (4 * sp.pi - 1) * sp.exp(-sp.pi), 32))
            ),
            "proof_certificate": (
                "Using 333/106<pi<22/7, the 40-term exponential Taylor "
                "lower sum proves 32*pi*exp(-4*pi)<1/2801. For "
                "a_n=n^2*exp(-pi*n^2), a_(n+1)/a_n<10^-6 on n>=2; "
                "the resulting geometric sum is below 1/2800."
            ),
        },
        "dominant_component_c2_budget": {
            "component_first_moment": (
                "integral v*dmu_(n,t)=2*exp(-pi*n^2)"
            ),
            "component_second_moment": (
                "integral v^2*dmu_(n,t)=4*integral_0^infinity "
                "s_(n,t)(v)dv"
            ),
            "component_third_moment": (
                "integral v^3*dmu_(n,t)=12*integral_0^infinity "
                "v*s_(n,t)(v)dv"
            ),
            "uniform_summand_bound": (
                "s_(n,t)(v)<=exp(-pi*n^2)"
                "*exp(-(4*pi*n^2-1)*v)"
            ),
            "tail_c0": (
                "|A_t-w_1*a_(1,t)|<1/2800"
            ),
            "tail_c1": (
                "|A_t'-w_1*a_(1,t)'|<1/137200"
            ),
            "tail_c2": (
                "|A_t''-w_1*a_(1,t)''|<1/3361400"
            ),
            "tail_c3": (
                "|A_t'''-w_1*a_(1,t)'''|<1/54900000"
            ),
            "modular_odd_jets": (
                "R^((2k+1))(0)=R'(0)=-1/2 for every k>=0"
            ),
            "first_component_fifth_jet": (
                "s_(1,0)^(5)(0)=exp(-pi)*"
                "(-1024*pi^5+11520*pi^4-33920*pi^3"
                "+26400*pi^2-3124*pi+1)>300"
            ),
            "tiny_mass_large_jet": (
                "sum_(n>=2)s_(n,0)^(5)(0)="
                "-1/2-s_(1,0)^(5)(0)<-601/2"
            ),
            "scope": (
                "The C0-C2 bounds form the dominant-block budget; the C3 "
                "bound supplies the derivative error below. Polynomial "
                "coefficients in J_t grow with x, and the tiny-mass tail cancels "
                "large first-block endpoint jets. A global proof therefore "
                "requires modularly coupled oscillatory decay rather than this "
                "compact-frequency budget alone."
            ),
        },
        "positive_transform": {
            "definition": (
                "C_t(x)=integral_0^infinity S_t(u)*cos(x*u)du"
            ),
            "characteristic_function": (
                "A_t(x)=integral_[0,infinity)cos(x*v)dmu_t(v)"
            ),
            "triangular_transform": (
                "C_t(x)=(1/(2*x^2))*integral(1-cos(x*v))dmu_t(v)"
                "=(1-A_t(x))/(2*x^2), x!=0"
            ),
            "strict_positivity": (
                "C_t(x)>0 for every real x; for x!=0, 0<x^2*C_t(x)<1"
            ),
            "scope": (
                "This is a strict Fourier-positivity theorem for the theta "
                "primitive S_t, not for the first correlation K_(1,t)."
            ),
        },
        "operator_factorization": {
            "operator": (
                "D_t=-4*t^2*d_x^2+4*t*x*d_x+(2*t-1-x^2)"
            ),
            "first_order": "Q_t=2*t*d_x-x",
            "factorization": "D_t=-(Q_t^2+1)",
            "heat_transform": (
                "H_t(x)=1/16+D_t[C_t](x)/8"
                "=(1-2*(Q_t^2+1)C_t(x))/16"
            ),
            "positive_normalizer": (
                "G_t(x)=8*H_t(x)/C_t(x) is globally real analytic because "
                "C_t(x)>0 on R."
            ),
        },
        "normalizer_evolution": {
            "heat_equations": (
                "partial_t H_t=-partial_x^2 H_t and "
                "partial_t C_t=-partial_x^2 C_t"
            ),
            "ratio": "G_t=8*H_t/C_t",
            "backward_pde": (
                "partial_t G_t=-partial_x^2 G_t"
                "-2*(partial_x log C_t)*partial_x G_t"
                "=-C_t^(-2)*partial_x(C_t^2*partial_x G_t)"
            ),
            "reverse_time": (
                "For tau=T-t, partial_tau G="
                "C^(-2)*partial_x(C^2*partial_x G)."
            ),
            "instantaneous_operator": (
                "A_t f=-C_t^(-2)*partial_x(C_t^2*partial_x f) "
                "is nonnegative in L2(C_t^2 dx): "
                "<f,A_t f>_(C_t^2)=integral C_t^2*(f')^2 dx."
            ),
            "nodal_interpretation": (
                "A common zero G_t(c)=partial_x G_t(c)=0 is a multiple "
                "node at which the reverse-time Sturm zero number may drop."
            ),
        },
        "sturm_nodal_guard": {
            "model": (
                "C_t(x)=1, H_t(x)=(x^2+a-2t)/8, "
                "G_t(x)=x^2+a-2t"
            ),
            "heat_identity": (
                "partial_t H_t=-partial_x^2 H_t and "
                "partial_t G_t=-partial_x^2 G_t"
            ),
            "nodal_transition": (
                "At t_*=a/2, G has a double zero at 0; for t>t_* "
                "it has two simple real zeros, while for t<t_* it has none."
            ),
            "reverse_time_consequence": (
                "Under tau=T-t, the positive heat semigroup loses two real "
                "nodes at the multiple zero. This is allowed by the Sturm "
                "zero-number theorem."
            ),
            "scope": (
                "A positive denominator, divergence form, nonnegative "
                "instantaneous operator, and zero-number monotonicity do not "
                "exclude contact. A completion needs an Xi-specific conserved "
                "nodal count, boundary flux, or quantitative transversality."
            ),
        },
        "normalizer_weight_asymptotic": {
            "endpoint_jets": (
                "S_t'(0)=-1/2 and "
                "S_t'''(0)=R'''(0)+6*t*R'(0)=-1/2-3*t"
            ),
            "cosine_expansion": (
                "C_t(x)=1/(2*x^2)-(1/2+3*t)/x^4+O_t(x^-6)"
            ),
            "uniformity": (
                "The O(x^-6) remainder is uniform for 0<=t<=1/5 "
                "because the required S_t derivatives have uniform L1 tails."
            ),
            "log_curvature": (
                "(log C_t)''=2/x^2-6*(1+6*t)/x^4+O(x^-6)>0 "
                "for all sufficiently large x, uniformly on 0<=t<=1/5."
            ),
            "doob_drift": (
                "2*(log C_t)'=-4/x+4*(1+6*t)/x^3+O(x^-5)"
            ),
            "scope": (
                "The positive normalizer weight is eventually log-convex, "
                "not log-concave. Therefore a global contracting-drift, "
                "Bakry-Emery positive-curvature, or weight-log-concavity "
                "argument cannot close the contact problem."
            ),
        },
        "characteristic_transfer": {
            "definition": "J_t(x)=16*x^4*H_t(x)",
            "formula": (
                "J_t=4*t^2*x^2*A_t''-4*t*x*(x^2+4*t)*A_t'"
                "+(x^4+(6*t+1)*x^2+24*t^2)*A_t"
                "-((6*t+1)*x^2+24*t^2)"
            ),
            "expectation_form": (
                "J_t=E_mu[((x^4+(6t+1)x^2+24t^2)-4t^2x^2V^2)"
                "*cos(xV)+4tx(x^2+4t)V*sin(xV)]"
                "-((6t+1)x^2+24t^2)"
            ),
            "contact_reduction": (
                "Because H_t(0)>0, a boundary multiple zero c is nonzero; "
                "H_t(c)=H_t'(c)=0 iff J_t(c)=J_t'(c)=0."
            ),
        },
        "componentwise_contact": {
            "definitions": (
                "A_(n,t)(x)=integral cos(xv)dmu_(n,t)(v), "
                "C_(n,t)=(w_n-A_(n,t))/(2*x^2), "
                "H_(n,t)=w_n/16+D_t[C_(n,t)]/8"
            ),
            "additivity": (
                "J_t=sum_(n>=1)J_(n,t), "
                "J_(n,t)=16*x^4*H_(n,t)"
            ),
            "polynomials": (
                "P_t=x^4+(6*t+1)*x^2+24*t^2, "
                "R_t=(6*t+1)*x^2+24*t^2"
            ),
            "component_formula": (
                "J_(n,t)=4*t^2*x^2*A_(n,t)''"
                "-4*t*x*(x^2+4*t)*A_(n,t)'"
                "+P_t*A_(n,t)-R_t*w_n"
            ),
            "component_derivative": (
                "J_(n,t)'=4*t^2*x^2*A_(n,t)'''"
                "-4*t*x*(x^2+2*t)*A_(n,t)''"
                "+(P_t-4*t*(3*x^2+4*t))*A_(n,t)'"
                "+(4*x^3+2*(6*t+1)*x)*A_(n,t)"
                "-2*(6*t+1)*x*w_n"
            ),
            "endpoint_first_primitive": (
                "C_(1,0)(x)=Re[(1/4)*pi^(-(1+i*x)/4)"
                "*Gamma((1+i*x)/4,pi)]"
            ),
        },
        "dominant_block_disjunction": {
            "tail_constants": (
                "delta_0=1/2800, delta_1=1/137200, "
                "delta_2=1/3361400, delta_3=1/54900000"
            ),
            "tail_value_bound": (
                "B_0(t,x)=4*t^2*x^2*delta_2"
                "+4*t*x*(x^2+4*t)*delta_1"
                "+(x^4+2*((6*t+1)*x^2+24*t^2))*delta_0"
            ),
            "tail_derivative_bound": (
                "B_1(t,x)=4*t^2*x^2*delta_3"
                "+4*t*x*(x^2+2*t)*delta_2"
                "+|P_t-4*t*(3*x^2+4*t)|*delta_1"
                "+(4*x^3+4*(6*t+1)*x)*delta_0"
            ),
            "bounds": (
                "|J_t-J_(1,t)|<B_0(t,x) and "
                "|J_t'-J_(1,t)'|<B_1(t,x), x>0"
            ),
            "pointwise_exclusion": (
                "At any 0<t<=1/5 and x>0, if "
                "|J_(1,t)(x)|>B_0(t,x) or "
                "|J_(1,t)'(x)|>B_1(t,x), then "
                "(J_t(x),J_t'(x))!=(0,0)."
            ),
            "global_corollary": (
                "If the pointwise disjunction holds for every "
                "0<t<=1/5 and x>0, then the boundary-contact criterion gives "
                "Lambda<=0."
            ),
            "scope": (
                "The displayed disjunction is a proved pointwise exclusion "
                "test, not an established global inequality. Since its raw "
                "moment bars contain delta_0*x^4 and 4*delta_0*x^3 terms, a "
                "completion needs a compact-frequency certificate plus a "
                "separate high-frequency estimate that retains modular "
                "endpoint cancellation."
            ),
        },
        "endpoint_laplace_comparison": {
            "formula": "16*x^2*H_0(x)=(1+x^2)*A_0(x)-1",
            "reference": (
                "1/(1+x^2) is the cosine characteristic function of the "
                "unit exponential law on [0,infinity)."
            ),
            "double_contact": (
                "H_0(c)=H_0'(c)=0 iff A_0(c)=1/(1+c^2) and "
                "A_0'(c)=-2*c/(1+c^2)^2."
            ),
        },
        "laguerre_product": {
            "normalization": "H_t=(C_t/8)*G_t",
            "identity": (
                "L[H_t]=(C_t^2/64)*L[G_t]+(G_t^2/64)*L[C_t]"
            ),
            "zero_value": (
                "At G_t(c)=0, L[H_t](c)=C_t(c)^2*G_t'(c)^2/64."
            ),
        },
        "generic_guard": {
            "kernel": (
                "f=(1-|.|)_+ * exp(-(.^2)/8), a smooth positive even "
                "Gaussian-tail kernel"
            ),
            "fourier_transform": str(guard_transform),
            "double_zero": "Fourier[f](2*pi)=Fourier[f]'(2*pi)=0",
            "neumann_lift": (
                "g(u)=-4*integral_0^infinity"
                " (exp(-|u-v|)+exp(-(u+v)))*f(v)dv"
            ),
            "lift_identities": "g''-g=8*f, g'(0)=0, g<0",
            "perturbation": (
                "R_epsilon(u)=exp(-u)/2+epsilon*g(u). For sufficiently "
                "small epsilon>0, R_epsilon is positive and decreasing, "
                "R_epsilon''=R_epsilon+8*epsilon*f>0, and "
                "R_epsilon'(0)=-1/2."
            ),
            "consequence": (
                "Its curvature law 2*R_epsilon''du is a probability and its "
                "primitive transform is strictly positive, while the associated "
                "H_epsilon=epsilon*integral_0^infinity f(u)cos(xu)du has a "
                "multiple real zero at 2*pi."
            ),
            "scope": (
                "Theta-curvature probability, primitive Fourier positivity, "
                "and the positive normalizer do not alone exclude contact. "
                "The remaining estimate must use the arithmetic density mu_t "
                "or the full characteristic differential expression J_t."
            ),
        },
        "open_handoff": (
            "Prove uniformly for 0<t<=1/5 that the Xi theta-curvature "
            "characteristic expression J_t has no common real zero with J_t'. "
            "Equivalently, use the explicit positive probability density "
            "dmu_t=2*d_u^2[exp(tu^2)R(u)]du to obtain a quantitative C1 "
            "separation for the displayed oscillatory expectation. The generic "
            "Neumann-lift guard forbids promotion from probability positivity, "
            "convexity, or C_t>0 alone. Any successful estimate must preserve "
            "the theta arithmetic and remain uniform as t tends to zero. The "
            "separate compact-transversality certificate now proves the displayed "
            "first-block disjunction on 1/4<=x<=38 and combines it with an exact "
            "origin collar. The live obligation is therefore x>38, where the raw "
            "B_0/B_1 bars must be abandoned in favor of the corrected "
            "Riemann-Siegel phase-aware partition."
        ),
        "checks": {
            "oscillator_residual": str(oscillator_residual),
            "characteristic_transfer_residual": str(
                sp.simplify(transfer - transfer_expected)
            ),
            "component_transfer_residual": str(
                sp.simplify(
                    component_transfer - component_transfer_expected
                )
            ),
            "component_derivative_residual": str(
                sp.simplify(
                    sp.diff(component_transfer_expected, x)
                    - component_derivative_expected
                )
            ),
            "ratio_pde_residual": str(ratio_pde_residual),
            "ratio_divergence_residual": str(ratio_divergence_residual),
            "sturm_guard_heat_residual": str(guard_heat_residual),
            "sturm_guard_ratio_residual": str(guard_ratio_residual),
            "primitive_third_modular": str(primitive_third_modular),
            "log_curvature_asymptotic_coefficient": str(
                log_curvature_lead
            ),
            "third_moment_tail_check": str(tail_c3_check),
            "endpoint_transfer": str(endpoint_transfer),
            "product_laguerre_residual": str(
                sp.simplify(product_laguerre - product_expected)
            ),
            "convex_floor_at_y12": str(convex_floor),
            "exp15_taylor_margin": str(exp15_lower - 2_250_000),
            "exp4pi_taylor_margin": str(exp4pi_lower - exp4pi_target),
            "geometric_tail_check": str(geometric_tail_check),
            "fifth_jet_polynomial": str(fifth_jet_polynomial),
            "fifth_jet_bernstein_min": str(min(jet_bernstein)),
            "fifth_jet_polynomial_floor_at_pi_lo": str(
                fifth_jet_polynomial.subs(z, pi_lo) - 7200
            ),
            "exp_pi_upper_margin": str(24 - exp_pi_upper),
            "guard_zero": str(guard_zero),
            "guard_slope": str(guard_slope),
            "guard_curvature": str(guard_curvature),
            "green_endpoint_slope": str(green_endpoint_slope),
            "green_derivative_jump": str(green_derivative_jump),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            "ntcpo_01_theta_primitive",
            "exact_identity",
            "ready_to_apply",
            "The positive theta primitive differentiates to the Xi kernel.",
            "R''-R=8*Phi and R'(0)=-1/2",
            "Exact theta-series differentiation and modular endpoint identity.",
        ),
        GateRow(
            "ntcpo_02_summand_monotonicity",
            "exact_inequality",
            "ready_to_apply",
            "Every Newman-weighted primitive summand decreases strictly.",
            "s_(n,t)'/s_(n,t)=1+2tu-y<0",
            "Uniform only on 0<=t<=1/5 and u>=0.",
        ),
        GateRow(
            "ntcpo_03_summand_convexity",
            "exact_inequality",
            "ready_to_apply",
            "Every Newman-weighted primitive summand is strictly convex.",
            exact["monotone_convex_proof"]["strict_margin"],
            "Uses y>=4*pi>12 and 1+2tu<=y/12.",
        ),
        GateRow(
            "ntcpo_04_probability_law",
            "exact_theorem",
            "ready_to_apply",
            "Theta curvature defines a probability measure.",
            "dmu_t=2*S_t''du, integral dmu_t=1",
            "Uses strict convexity and the exact endpoint slopes.",
        ),
        GateRow(
            "ntcpo_04a_fixed_theta_weights",
            "exact_theorem",
            "ready_to_apply",
            "Theta curvature is a fixed-weight positive component mixture.",
            exact["fixed_theta_mixture"]["component_mass"],
            "The weights are exact and independent of Newman time.",
        ),
        GateRow(
            "ntcpo_04b_dominant_component_budget",
            "exact_bound",
            "ready_to_apply",
            "The first theta component has a uniform rational C2 tail budget.",
            (
                exact["dominant_component_c2_budget"]["tail_c0"]
                + ", "
                + exact["dominant_component_c2_budget"]["tail_c1"]
                + ", "
                + exact["dominant_component_c2_budget"]["tail_c2"]
            ),
            "Absolute bounds alone do not control the polynomially amplified global-frequency tail.",
        ),
        GateRow(
            "ntcpo_05_triangular_mixture",
            "exact_identity",
            "ready_to_apply",
            "The primitive is a positive mixture of triangular kernels.",
            exact["probability_law"]["triangular_mixture"],
            "Tonelli applies because the curvature density is positive.",
        ),
        GateRow(
            "ntcpo_06_primitive_fourier_positivity",
            "exact_theorem",
            "ready_to_apply",
            "The primitive cosine transform is strictly positive.",
            exact["positive_transform"]["triangular_transform"],
            "This does not assert Fourier positivity of K_(1,t).",
        ),
        GateRow(
            "ntcpo_07_oscillator_factorization",
            "exact_identity",
            "ready_to_apply",
            "The endpoint-subtracted transform operator factors exactly.",
            "D_t=-(Q_t^2+1), Q_t=2t*d_x-x",
            "Pure differential algebra.",
        ),
        GateRow(
            "ntcpo_08_positive_normalizer",
            "exact_reduction",
            "ready_to_apply",
            "The strict-Laguerre target has a global positive normalizer.",
            "G_t=8*H_t/C_t, C_t>0 on R",
            "The normalized contact problem remains RH-strength.",
        ),
        GateRow(
            "ntcpo_08a_doob_evolution",
            "exact_identity",
            "ready_to_apply",
            "The positive-normalizer ratio obeys a weighted Doob diffusion.",
            exact["normalizer_evolution"]["backward_pde"],
            "The instantaneous reverse-time generator is nonnegative in L2(C_t^2 dx).",
        ),
        GateRow(
            "ntcpo_08b_sturm_nodal_guard",
            "countermodel_gate",
            "guard_validated",
            "Positive operator evolution and Sturm zero-number monotonicity permit nodal loss at contact.",
            exact["sturm_nodal_guard"]["nodal_transition"],
            exact["sturm_nodal_guard"]["scope"],
        ),
        GateRow(
            "ntcpo_08c_drift_curvature_guard",
            "asymptotic_guard",
            "guard_validated",
            "The positive Doob weight is eventually log-convex.",
            exact["normalizer_weight_asymptotic"]["log_curvature"],
            exact["normalizer_weight_asymptotic"]["scope"],
        ),
        GateRow(
            "ntcpo_09_characteristic_transfer",
            "exact_identity",
            "ready_to_apply",
            "H_t is an explicit second-order expression in a probability characteristic function.",
            exact["characteristic_transfer"]["formula"],
            "The apparent singularity at x=0 is removable; contacts occur at c!=0.",
        ),
        GateRow(
            "ntcpo_09a_componentwise_contact",
            "exact_identity",
            "ready_to_apply",
            "The contact functional splits exactly over the fixed theta components.",
            exact["componentwise_contact"]["component_formula"],
            "Both the value and derivative decompositions are audited symbolically.",
        ),
        GateRow(
            "ntcpo_10_contact_reduction",
            "exact_reduction",
            "ready_to_apply",
            "A positive Newman boundary is a common zero of J_t and J_t'.",
            exact["characteristic_transfer"]["contact_reduction"],
            "No zero simplicity or RH assumption enters.",
        ),
        GateRow(
            "ntcpo_10a_dominant_block_disjunction",
            "exact_reduction",
            "ready_to_apply",
            "Explicit first-block C1 alternatives exclude contact pointwise.",
            exact["dominant_block_disjunction"]["pointwise_exclusion"],
            "The raw moment bars are compact-frequency tools; no global inequality is claimed.",
        ),
        GateRow(
            "ntcpo_11_endpoint_laplace_reference",
            "exact_identity",
            "ready_to_apply",
            "At t=0, Xi is a transverse-crossing problem against the exponential characteristic function.",
            exact["endpoint_laplace_comparison"]["double_contact"],
            "This is a reformulation, not a proof of transversality.",
        ),
        GateRow(
            "ntcpo_12_generic_double_contact_guard",
            "countermodel_gate",
            "guard_validated",
            "Positive curvature probability and primitive Fourier positivity do not exclude a multiple spectral zero.",
            exact["generic_guard"]["consequence"],
            "The guard is generic and does not reproduce theta arithmetic.",
        ),
        GateRow(
            "ntcpo_13_open_xi_transversality",
            "open_theorem_target",
            "open",
            "The remaining theorem is arithmetic C1 separation for J_t.",
            exact["open_handoff"],
            "No strict-Laguerre, Lambda<=0, RH, or Clay-prize claim is made.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_theta_curvature_probability_operator_gate",
        "date": "2026-07-24",
        "status": (
            "exact theta-curvature probability/operator and dominant-component "
            "reduction with an explicit C1 tail-disjunction target and a "
            "generic double-contact guard and one open Xi transversality handoff"
        ),
        "proof_boundary": (
            "The artifact proves strict positivity of the theta-primitive "
            "transform C_t and an exact probability-characteristic reduction. "
            "It does not prove strict positivity of the first correlation "
            "transform, positive-time simplicity, Lambda<=0, or RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    return "\n".join(
        [
            "# Newman Theta-Curvature Probability/Operator Gate",
            "",
            "Date: 2026-07-24",
            "",
            "Status: exact theta-curvature probability/operator reduction with",
            "a generic double-contact guard. This is not a proof of",
            "`Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_curvature_probability_operator_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_curvature_probability_operator_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_curvature_probability_operator_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Newman theta-curvature probability/operator gate: 20 rows, 0 issues, 3 theta-primitive identities, 2 monotone-convex inequalities, 1 probability law, 1 fixed-weight theta mixture, 1 uniform C2 dominant-block budget, 1 dominant-mass endpoint-jet guard, 2 transform identities, 1 oscillator factorization, 1 Doob diffusion identity, 1 Sturm nodal-loss guard, 1 asymptotic drift-curvature guard, 1 characteristic contact reduction, 1 componentwise contact decomposition, 1 explicit C1 tail disjunction, 1 endpoint Laplace comparison, 1 generic double-contact guard, 1 open Xi transversality handoff",
            "```",
            "",
            "## Curvature Probability",
            "",
            "For the positive theta primitive, set",
            "",
            "```text",
            exact["theta_primitive"]["definition"],
            exact["theta_primitive"]["differential_identity"],
            exact["theta_primitive"]["modular_endpoint"],
            exact["theta_primitive"]["newman_primitive"],
            "```",
            "",
            "Writing `y=4*pi*n^2*exp(4u)`, every summand obeys",
            "",
            "```text",
            exact["monotone_convex_proof"]["log_derivative"],
            exact["monotone_convex_proof"]["second_derivative"],
            exact["monotone_convex_proof"]["elementary_bound"],
            exact["monotone_convex_proof"]["strict_margin"],
            "```",
            "",
            "Thus `S_t` is positive, strictly decreasing, and strictly convex",
            "uniformly for `0<=t<=1/5`. Its curvature is a probability law:",
            "",
            "```text",
            exact["probability_law"]["definition"],
            exact["probability_law"]["mass"],
            exact["probability_law"]["triangular_mixture"],
            "```",
            "",
            "The probability law has a fixed arithmetic component split:",
            "",
            "```text",
            exact["fixed_theta_mixture"]["component_measure"],
            exact["fixed_theta_mixture"]["component_mass"],
            exact["fixed_theta_mixture"]["time_independence"],
            exact["fixed_theta_mixture"]["normalized_components"],
            exact["fixed_theta_mixture"]["dominant_weight"],
            exact["fixed_theta_mixture"]["exact_first_weight"],
            "```",
            "",
            "Integration by parts gives a uniform dominant-block budget:",
            "",
            "```text",
            exact["dominant_component_c2_budget"]["component_first_moment"],
            exact["dominant_component_c2_budget"]["component_second_moment"],
            exact["dominant_component_c2_budget"]["component_third_moment"],
            exact["dominant_component_c2_budget"]["uniform_summand_bound"],
            exact["dominant_component_c2_budget"]["tail_c0"],
            exact["dominant_component_c2_budget"]["tail_c1"],
            exact["dominant_component_c2_budget"]["tail_c2"],
            exact["dominant_component_c2_budget"]["tail_c3"],
            "```",
            "",
            "The mass hierarchy cannot be promoted to a global one-block",
            "approximation. Modular evenness gives the exact endpoint identities",
            "",
            "```text",
            exact["dominant_component_c2_budget"]["modular_odd_jets"],
            exact["dominant_component_c2_budget"]["first_component_fifth_jet"],
            exact["dominant_component_c2_budget"]["tiny_mass_large_jet"],
            "```",
            "",
            "Thus a tail carrying less than `1/2800` of the probability mass",
            "cancels a fifth endpoint jet of magnitude above `300`.",
            "",
            "The polynomial coefficients in `J_t` still amplify fixed absolute",
            "errors at high frequency, so this is a rigorous compact-frequency",
            "handoff rather than a global positivity proof.",
            "",
            "## Positive Primitive Transform",
            "",
            "Let",
            "",
            "```text",
            exact["positive_transform"]["definition"],
            exact["positive_transform"]["characteristic_function"],
            "```",
            "",
            "Twice integrating the triangular mixture gives",
            "",
            "```text",
            exact["positive_transform"]["triangular_transform"],
            exact["positive_transform"]["strict_positivity"],
            "```",
            "",
            "This is a genuine strict Fourier-positivity theorem, but for the",
            "theta primitive rather than for `K_(1,t)`.",
            "",
            "## Operator And Contact",
            "",
            "The endpoint-subtracted operator has the exact factorization",
            "",
            "```text",
            exact["operator_factorization"]["operator"],
            exact["operator_factorization"]["first_order"],
            exact["operator_factorization"]["factorization"],
            exact["operator_factorization"]["heat_transform"],
            "```",
            "",
            "Because both numerator and denominator solve the same backward",
            "heat equation, the positive-normalizer ratio also has an exact",
            "evolution law:",
            "",
            "```text",
            exact["normalizer_evolution"]["heat_equations"],
            exact["normalizer_evolution"]["ratio"],
            exact["normalizer_evolution"]["backward_pde"],
            exact["normalizer_evolution"]["reverse_time"],
            exact["normalizer_evolution"]["instantaneous_operator"],
            "```",
            "",
            "This weighted Sturm form does not itself exclude nodal loss. The",
            "exact calibration model",
            "",
            "```text",
            exact["sturm_nodal_guard"]["model"],
            exact["sturm_nodal_guard"]["nodal_transition"],
            exact["sturm_nodal_guard"]["reverse_time_consequence"],
            "```",
            "",
            "shows that a nonnegative generator and the zero-number theorem",
            "permit two nodes to disappear at a multiple zero. An Xi-specific",
            "conserved nodal flux or quantitative transversality is still",
            "required.",
            "",
            "The modular endpoint jets also determine the Doob weight at high",
            "frequency:",
            "",
            "```text",
            exact["normalizer_weight_asymptotic"]["endpoint_jets"],
            exact["normalizer_weight_asymptotic"]["cosine_expansion"],
            exact["normalizer_weight_asymptotic"]["log_curvature"],
            exact["normalizer_weight_asymptotic"]["doob_drift"],
            "```",
            "",
            "Thus `C_t` is eventually log-convex, uniformly on the target",
            "window. A global contracting-drift or positive Bakry-Emery",
            "curvature argument is therefore unavailable.",
            "",
            "Since `C_t>0`, it supplies a global positive normalizer. More",
            "concretely, with `J_t=16*x^4*H_t`,",
            "",
            "```text",
            exact["characteristic_transfer"]["formula"],
            exact["characteristic_transfer"]["contact_reduction"],
            "```",
            "",
            "The same identity can be read as the explicit oscillatory",
            "expectation",
            "",
            "```text",
            exact["characteristic_transfer"]["expectation_form"],
            "```",
            "",
            "The fixed theta mixture makes the contact equation additive:",
            "",
            "```text",
            exact["componentwise_contact"]["definitions"],
            exact["componentwise_contact"]["additivity"],
            exact["componentwise_contact"]["polynomials"],
            exact["componentwise_contact"]["component_formula"],
            exact["componentwise_contact"]["component_derivative"],
            exact["componentwise_contact"]["endpoint_first_primitive"],
            "```",
            "",
            "Writing `J_(1,t)` for the unnormalized first component, the",
            "remaining theta tail obeys the explicit error bars",
            "",
            "```text",
            exact["dominant_block_disjunction"]["tail_constants"],
            exact["dominant_block_disjunction"]["tail_value_bound"],
            exact["dominant_block_disjunction"]["tail_derivative_bound"],
            exact["dominant_block_disjunction"]["bounds"],
            "```",
            "",
            "Consequently the following is a rigorous pointwise exclusion test:",
            "",
            "```text",
            exact["dominant_block_disjunction"]["pointwise_exclusion"],
            exact["dominant_block_disjunction"]["global_corollary"],
            "```",
            "",
            "The implication is proved, but the raw moment bars grow too",
            "quickly to serve as an all-frequency estimate. The open theorem",
            "requires a compact certificate joined to a modularly coupled",
            "high-frequency argument.",
            "",
            "At the undeformed endpoint this collapses to a Laplace reference:",
            "",
            "```text",
            exact["endpoint_laplace_comparison"]["formula"],
            exact["endpoint_laplace_comparison"]["reference"],
            exact["endpoint_laplace_comparison"]["double_contact"],
            "```",
            "",
            "Thus a multiple endpoint zero is exactly tangential contact between",
            "the theta-curvature characteristic function and that of the unit",
            "exponential law.",
            "",
            "## Nonpromotion Guard",
            "",
            "The positive normalization is not itself the missing theorem.",
            "Let",
            "",
            "```text",
            exact["generic_guard"]["kernel"],
            exact["generic_guard"]["neumann_lift"],
            exact["generic_guard"]["lift_identities"],
            exact["generic_guard"]["perturbation"],
            "```",
            "",
            "The kernel has a double Fourier zero at `2*pi`. For small positive",
            "`epsilon`, its Neumann lift has all the curvature-probability and",
            "primitive-transform positivity properties above, while the",
            "associated transform still has that multiple zero. Hence generic",
            "probability positivity, convexity, and `C_t>0` cannot be promoted",
            "to the Newman conclusion.",
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Newman theta-curvature probability/operator gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
