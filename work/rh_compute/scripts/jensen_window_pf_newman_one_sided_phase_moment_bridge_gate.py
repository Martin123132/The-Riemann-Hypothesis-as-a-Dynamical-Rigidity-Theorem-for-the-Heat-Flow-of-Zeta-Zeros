#!/usr/bin/env python3
"""Build the exact one-sided phase/probability bridge for the Newman kernel."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.md"
)
COARSE_DPS = 35
FINE_DPS = 45
THETA_TERMS = 8
FIRST_XI_ZERO_X = (
    "28.2694502834693875809145039671249405415685142313984864"
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
    f0, x, sine, cosine, sine_prime, ell = sp.symbols(
        "f0 x B C B_prime ell", positive=True, real=True
    )
    chi = cosine + sp.I * sine
    lift = sp.expand(sp.I * f0 * (1 - chi) / x)
    if sp.simplify(sp.re(lift) - f0 * sine / x) != 0:
        raise RuntimeError("one-sided lift real-part identity failed")
    if sp.simplify(sp.im(lift) - f0 * (1 - cosine) / x) != 0:
        raise RuntimeError("one-sided lift imaginary-part identity failed")

    h = f0 * sine / x
    h_prime = sp.diff(h, x) + f0 * sine_prime / x
    expected_h_prime = f0 * (sine_prime - sine / x) / x
    if sp.simplify(h_prime - expected_h_prime) != 0:
        raise RuntimeError("characteristic first-jet identity failed")
    scaled_jet = sp.expand(h**2 + (expected_h_prime / ell) ** 2)
    expected_scaled_jet = (
        f0**2
        / x**2
        * (sine**2 + (sine_prime - sine / x) ** 2 / ell**2)
    )
    if sp.simplify(scaled_jet - expected_scaled_jet) != 0:
        raise RuntimeError("scaled characteristic first-jet identity failed")

    k, raw_moment, odd_moment = sp.symbols(
        "k raw_moment odd_moment", positive=True
    )
    odd_relation = (2 * k + 1) * raw_moment / (2 * f0)
    coefficient_from_odd = (
        2
        * f0
        * odd_moment
        / sp.gamma(2 * k + 2)
    )
    coefficient_from_raw = raw_moment / sp.gamma(2 * k + 1)
    if (
        sp.simplify(
            coefficient_from_odd.subs(odd_moment, odd_relation)
            - coefficient_from_raw
        )
        != 0
    ):
        raise RuntimeError("odd-moment coefficient identity failed")
    z = sp.symbols("z")
    for degree in range(1, 5):
        for shift in range(4):
            finite_kernel = sp.Add(
                *[
                    sp.binomial(degree, index)
                    * z**index
                    / sp.factorial(shift + index)
                    for index in range(degree + 1)
                ]
            )
            laguerre_kernel = (
                sp.factorial(degree)
                / sp.factorial(degree + shift)
                * sp.assoc_laguerre(degree, shift, -z)
            )
            if sp.simplify(finite_kernel - laguerre_kernel) != 0:
                raise RuntimeError("shifted Laguerre kernel identity failed")
    for index in range(5):
        beta_moment = (
            sp.factorial(index) ** 2 / sp.factorial(2 * index + 1)
        )
        score_odd = (
            (2 * index + 1) * raw_moment / (2 * f0)
        )
        rho_moment = sp.simplify(
            2 * f0 * score_odd * beta_moment
        )
        expected_rho_moment = (
            sp.factorial(index) ** 2
            * raw_moment
            / sp.factorial(2 * index)
        )
        if sp.simplify(rho_moment - expected_rho_moment) != 0:
            raise RuntimeError("score-Beta measure moment identity failed")
        next_moment = sp.symbols(f"M_{index + 1}", positive=True)
        flow_moment = (
            4 * next_moment
            - 2 * next_moment / (index + 1)
        )
        expected_flow_moment = (
            2
            * (2 * index + 1)
            * next_moment
            / (index + 1)
        )
        if sp.simplify(flow_moment - expected_flow_moment) != 0:
            raise RuntimeError("score-Beta measure flow identity failed")
    chi_function = sp.Function("chi")
    lift_function = sp.I * f0 * (1 - chi_function(x)) / x
    solved_chi_time = sp.simplify(
        x * sp.diff(lift_function, x, 2) / (sp.I * f0)
    )
    expected_chi_time = (
        -sp.diff(chi_function(x), x, 2)
        + 2 * sp.diff(chi_function(x), x) / x
        + 2 * (1 - chi_function(x)) / x**2
    )
    if sp.simplify(solved_chi_time - expected_chi_time) != 0:
        raise RuntimeError("score characteristic flow identity failed")
    gram = sp.Matrix([[1, 0], [sp.Rational(1, 2), 1]])
    gram = gram.T * gram
    transfer_min = (9 - sp.sqrt(17)) / 8
    transfer_max = (9 + sp.sqrt(17)) / 8
    if sp.simplify(sp.trace(gram) - transfer_min - transfer_max) != 0:
        raise RuntimeError("normalizer-transfer trace identity failed")
    if sp.simplify(gram.det() - transfer_min * transfer_max) != 0:
        raise RuntimeError("normalizer-transfer determinant identity failed")

    return {
        "imported_shape": (
            "For 0<=t<=1/5, f_t(u)=exp(t*u^2)*Phi(u) is positive, even, "
            "C-infinity, strictly decreasing for u>0, and super-exponentially "
            "decaying; (log f_t)''<=-(kappa-2t), "
            "kappa=-Phi''(0)/Phi(0)>74.9076"
        ),
        "score_probability": (
            "dnu_t(u)=-f_t'(u)du/f_t(0), u>0, is a nondegenerate probability "
            "law with a continuous positive density"
        ),
        "complex_lift": (
            "F_t(x)=integral_0^infinity f_t(u)exp(i*x*u)du=H_t(x)+i*Y_t(x); "
            "chi_t(x)=integral_0^infinity exp(i*x*u)dnu_t(u); "
            "F_t(x)=i*f_t(0)*(1-chi_t(x))/x for x>0"
        ),
        "cartesian_lift": (
            "Writing chi_t=C_t+i*B_t gives "
            "H_t=f_t(0)*B_t/x and Y_t=f_t(0)*(1-C_t)/x>0 for x>0"
        ),
        "phase_contact": (
            "For Theta_t=arg(F_t) in (0,pi), at H_t(x)=0 one has "
            "Theta_t=pi/2, Theta_t'=-H_t'/Y_t=-B_t'/(1-C_t); "
            "B_t=E_nu[sin(xU)] and B_t'=E_nu[U*cos(xU)]; "
            "therefore H_t=H_t'=0 iff B_t=B_t'=0"
        ),
        "scaled_first_jet": (
            "T_L[H_t]=H_t^2+(H_t'/L)^2="
            "f_t(0)^2/x^2*(B_t^2+(B_t'-B_t/x)^2/L^2)"
        ),
        "coefficient_bridge": (
            "If mu_(2k)(t)=integral_R u^(2k)f_t(u)du and "
            "c_k(t)=mu_(2k)(t)/(2k)!, then "
            "E_nu[U^(2k+1)]=(2k+1)*mu_(2k)/(2*f_t(0)), "
            "c_k=2*f_t(0)*E_nu[U^(2k+1)]/(2k+1)!, "
            "A_k=k!*c_k, and 2*H_t(sqrt(-z))=sum_(k>=0)c_k*z^k"
        ),
        "shifted_laguerre_bridge": (
            "For q=s*(1-s) and "
            "P_(D,n)(w)=sum_(k=0)^D binom(D,k)A_(n+k)w^k, "
            "P_(D,n)(w)=2*f_t(0)*D!/(D+n)!*integral_0^1 "
            "E_nu[U^(2n+1)*q^n*L_D^(n)(-q*U^2*w)]ds. "
            "For every U>0 and 0<s<1, the fixed kernel has D simple "
            "negative w-roots; "
            "the remaining issue is preservation under this specific positive "
            "score/Beta scale mixture"
        ),
        "score_beta_abel_measure": (
            "Let S be uniform on (0,1), independent of U~nu_t, "
            "Q=S*(1-S), and define the finite measure "
            "rho_t(E)=2*f_t(0)*E_nu[U*1_(Q*U^2 in E)]. Then "
            "rho_t has density r_t(v)=4*integral_(2*sqrt(v))^infinity "
            "(-f_t'(u))/sqrt(u^2-4*v)du, v>0, and "
            "M_k(t)=integral_0^infinity v^k*r_t(v)dv=k!*A_k(t)"
        ),
        "abel_bessel_laguerre_bridge": (
            "2*H_t(sqrt(-z))=integral_0^infinity "
            "I_0(2*sqrt(v*z))*r_t(v)dv and "
            "P_(D,n)(w)=D!/(D+n)!*integral_0^infinity "
            "v^n*L_D^(n)(-v*w)*r_t(v)dv"
        ),
        "abel_flow_concentration_guard": (
            "Writing R_t(v)=integral_v^infinity r_t(s)ds, "
            "partial_t r_t(v)=4*v*r_t(v)-2*R_t(v), equivalently "
            "M_k'=2*(2*k+1)*M_(k+1)/(k+1). Moreover "
            "disc P_(2,n)=4*(A_(n+1)^2-A_n*A_(n+2))>=0 iff "
            "M_(n+1)^2/(M_n*M_(n+2))>=(n+1)/(n+2); positivity of rho_t "
            "alone gives only the opposite-side Cauchy-Schwarz upper bound "
            "M_(n+1)^2/(M_n*M_(n+2))<=1"
        ),
        "full_scale_interlacing_guard": (
            "Because -f_t'(u)>0 for u>0, r_t(v)>0 for every v>0. If "
            "0<ell_1<...<ell_D are the roots of L_D^(n), then the roots of "
            "L_D^(n)(-v*w) are -ell_D/v<...<-ell_1/v. For D>=2 the family "
            "{L_D^(n)(-v*w):v>0} has no common interlacer: as v tends to 0 "
            "every root tends to -infinity, while as v tends to infinity "
            "every root tends to 0 from below. Hence a direct global "
            "common-interlacing mixture theorem cannot close the Xi target; "
            "a weighted Xi-specific total-positive or variation-diminishing "
            "connection would still be sufficient"
        ),
        "abel_kernel_sign_guard": (
            "For K(x,y)=(y-x)_+^(-1/2), the 2x2 minor at "
            "x=(0,1), y=(2,3) is 1/2-1/sqrt(3)<0, whereas the minor at "
            "x=(0,2), y=(1,3) is 1>0. Thus the bare Abel kernel is neither "
            "TP_2 nor sign-regular of order 2; any variation-diminishing "
            "closure must use the Xi weight or a larger composed kernel"
        ),
        "score_flow": (
            "Writing m_t(u)=-f_t'(u)/f_t(0) and "
            "bar_nu_t(u)=integral_u^infinity m_t(v)dv=f_t(u)/f_t(0), "
            "partial_t m_t=u^2*m_t-2*u*bar_nu_t. Equivalently, "
            "partial_t chi_t=-chi_t''+2*chi_t'/x+2*(1-chi_t)/x^2 and "
            "partial_t B_t=-B_t''+2*B_t'/x-2*B_t/x^2"
        ),
        "normalizer_transfer": (
            "Let mathcal_A_t=|M_t((1-i*x)/2)|, Z_t=H_t/mathcal_A_t, "
            "g_t=x*mathcal_A_t/f_t(0), and a_t=(log mathcal_A_t)'. Then "
            "B_t=g_t*Z_t and B_t'-B_t/x=g_t*(Z_t'+a_t*Z_t). "
            "On L>=50, 0<tL<=25, |a_t|/L<1/2, so with "
            "m_-=(9-sqrt(17))/8 and m_+=(9+sqrt(17))/8, "
            "m_-*g_t^2*T_L[Z_t] <= "
            "B_t^2+(B_t'-B_t/x)^2/L^2 <= "
            "m_+*g_t^2*T_L[Z_t]. Equivalently, "
            "m_-*mathcal_A_t^2*T_L[Z_t] <= T_L[H_t] <= "
            "m_+*mathcal_A_t^2*T_L[Z_t]"
        ),
        "phase_conditioning": (
            "Uniformly for 0<=t<=1/5, "
            "Y_t(x)=f_t(0)/x-f_t''(0)/x^3+O(x^-5)="
            "f_t(0)/x*(1+(kappa-2t)/x^2+O(x^-4)), while H_t and H_t' "
            "are rapidly decreasing; hence pi/2-Theta_t and Theta_t' are "
            "O_A(x^-A) for every A"
        ),
        "conditioning_guard": (
            "The exact one-sided phase is globally nonvanishing but "
            "asymptotically flat. No fixed polynomial absolute lower bound for "
            "|Theta_t'| can be transferred from the corrected Riemann-Siegel "
            "phase; quantitative margins must retain their own lift amplitude."
        ),
        "countermodel": (
            "For f=K_(1,2)=(1-|.|)_+*exp(-(.^2)/8), f is positive, even, "
            "strictly decreasing and strongly log-concave, while "
            "Fourier[f](xi)=8*sqrt(2*pi)*exp(-2*xi^2)*sin(xi/2)^2/xi^2 "
            "has double zeros at xi=2*pi*n; its score probability therefore "
            "has B=B'=0 at those points"
        ),
        "countermodel_quadratic_coefficient": (
            "lim_(xi->2*pi) Fourier[f](xi)/(xi-2*pi)^2="
            "sqrt(2*pi)*exp(-8*pi^2)/(2*pi^2)>0"
        ),
        "live_target": (
            "Prove for the Xi score laws nu_t that "
            "(B_t(x),B_t'(x))!=(0,0) for every x>38 and 0<t<=1/5, "
            "using theta arithmetic or an all-order PF/sign-regular theorem; "
            "probability, strong log-concavity, and the zero-free complex lift "
            "alone are insufficient"
        ),
    }


def phi_value(u: mp.mpf) -> mp.mpf:
    return mp.fsum(
        (
            2 * mp.pi**2 * n**4 * mp.exp(9 * u)
            - 3 * mp.pi * n**2 * mp.exp(5 * u)
        )
        * mp.exp(-mp.pi * n**2 * mp.exp(4 * u))
        for n in range(1, THETA_TERMS + 1)
    )


def phi_first(u: mp.mpf) -> mp.mpf:
    total = mp.mpf(0)
    for n in range(1, THETA_TERMS + 1):
        saddle = mp.pi * n**2
        e4 = mp.exp(4 * u)
        e5 = mp.exp(5 * u)
        e9 = mp.exp(9 * u)
        main = 2 * saddle**2 * e9 - 3 * saddle * e5
        main_first = 18 * saddle**2 * e9 - 15 * saddle * e5
        total += (
            main_first - 4 * saddle * e4 * main
        ) * mp.exp(-saddle * e4)
    return total


def integrate_half_line(function) -> mp.mpc:
    points = [
        mp.mpf("0"),
        mp.mpf("0.2"),
        mp.mpf("0.4"),
        mp.mpf("0.6"),
        mp.mpf("0.8"),
        mp.mpf("1"),
        mp.mpf("1.25"),
        mp.mpf("1.5"),
        mp.mpf("2"),
    ]
    return mp.fsum(
        mp.quad(function, [left, right])
        for left, right in zip(points[:-1], points[1:], strict=True)
    )


def completed_xi(s: mp.mpc) -> mp.mpc:
    return (
        mp.mpf("0.5")
        * s
        * (s - 1)
        * mp.pi ** (-s / 2)
        * mp.gamma(s / 2)
        * mp.zeta(s)
    )


def diagnostics(dps: int) -> dict:
    mp.mp.dps = dps
    f0 = phi_value(mp.mpf(0))
    cases = [
        ("origin_scale", mp.mpf(0), mp.mpf(1)),
        ("first_xi_crossing", mp.mpf(0), mp.mpf(FIRST_XI_ZERO_X)),
        ("positive_time_edge", mp.mpf("0.2"), mp.mpf(38)),
    ]
    rows = []
    tolerance = mp.power(10, -dps + 8)
    for label, time, x in cases:
        def deformed(u: mp.mpf) -> mp.mpf:
            return mp.exp(time * u**2) * phi_value(u)

        def score_density(u: mp.mpf) -> mp.mpf:
            return (
                -mp.exp(time * u**2)
                * (2 * time * u * phi_value(u) + phi_first(u))
                / f0
            )

        lift = integrate_half_line(
            lambda u: deformed(u) * mp.exp(mp.j * x * u)
        )
        lift_first = integrate_half_line(
            lambda u: mp.j * u * deformed(u) * mp.exp(mp.j * x * u)
        )
        chi = integrate_half_line(
            lambda u: score_density(u) * mp.exp(mp.j * x * u)
        )
        chi_first = integrate_half_line(
            lambda u: mp.j * u * score_density(u) * mp.exp(mp.j * x * u)
        )
        lift_from_chi = mp.j * f0 * (1 - chi) / x
        a_value = mp.im(chi)
        a_first = mp.im(chi_first)
        h_from_chi = f0 * a_value / x
        h_first_from_chi = f0 * (a_first - a_value / x) / x
        phase_direct = mp.im(lift_first / lift)
        phase_from_chi = mp.im(-chi_first / (1 - chi))
        row = {
            "label": label,
            "t": mp.nstr(time, 20, strip_zeros=False),
            "x": mp.nstr(x, 35, strip_zeros=False),
            "F_real": mp.nstr(mp.re(lift), 30, strip_zeros=False),
            "F_imag": mp.nstr(mp.im(lift), 30, strip_zeros=False),
            "B": mp.nstr(a_value, 30, strip_zeros=False),
            "B_prime": mp.nstr(a_first, 30, strip_zeros=False),
            "phase_velocity": mp.nstr(
                phase_direct, 30, strip_zeros=False
            ),
            "complex_lift_abs_residual": mp.nstr(
                abs(lift - lift_from_chi), 18, strip_zeros=False
            ),
            "H_characteristic_abs_residual": mp.nstr(
                abs(mp.re(lift) - h_from_chi), 18, strip_zeros=False
            ),
            "H_first_characteristic_abs_residual": mp.nstr(
                abs(mp.re(lift_first) - h_first_from_chi),
                18,
                strip_zeros=False,
            ),
            "phase_formula_abs_residual": mp.nstr(
                abs(phase_direct - phase_from_chi),
                18,
                strip_zeros=False,
            ),
        }
        if time == 0:
            s = (1 + mp.j * x) / 2
            exact_h = mp.re(completed_xi(s)) / 8
            row["exact_Xi_H"] = mp.nstr(
                exact_h, 30, strip_zeros=False
            )
            row["Xi_integral_abs_residual"] = mp.nstr(
                abs(mp.re(lift) - exact_h), 18, strip_zeros=False
            )
        if mp.im(lift) <= 0:
            raise RuntimeError(f"one-sided imaginary lift failed at {label}")
        for key in (
            "complex_lift_abs_residual",
            "H_characteristic_abs_residual",
            "H_first_characteristic_abs_residual",
            "phase_formula_abs_residual",
        ):
            if mp.mpf(row[key]) >= tolerance:
                raise RuntimeError(f"{label} diagnostic failed: {key}")
        rows.append(row)
    return {
        "role": "normalization_and_conditioning_diagnostics_only",
        "proof_boundary": (
            "The rows numerically cross-check exact identities and the first "
            "Xi crossing. They do not establish the cofinal no-contact target."
        ),
        "dps": dps,
        "theta_terms": THETA_TERMS,
        "rows": rows,
    }


def compare_diagnostics(coarse: dict, fine: dict) -> dict:
    maximum = mp.mpf(0)
    fields = ("F_real", "F_imag", "B", "B_prime", "phase_velocity")
    for left, right in zip(coarse["rows"], fine["rows"], strict=True):
        if left["label"] != right["label"]:
            raise RuntimeError("diagnostic case ordering drifted")
        for field in fields:
            maximum = max(
                maximum,
                abs(mp.mpf(left[field]) - mp.mpf(right[field])),
            )
    if maximum >= mp.mpf("1e-24"):
        raise RuntimeError("one-sided phase diagnostics are precision-unstable")
    return {
        "max_abs_selected_field_delta": mp.nstr(
            maximum, 18, strip_zeros=False
        )
    }


def build_artifact() -> dict:
    exact = build_exact()
    coarse = diagnostics(COARSE_DPS)
    fine = diagnostics(FINE_DPS)
    convergence = compare_diagnostics(coarse, fine)
    rows = [
        GateRow(
            id="nospmb_01_imported_shape",
            role="imported_exact_theorem",
            readiness="ready_to_apply",
            claim="The positive-time Xi half-kernel has the shape needed for a score probability.",
            formula=exact["imported_shape"],
            proof_boundary="Imports the already checked strong-log-concavity gate.",
        ),
        GateRow(
            id="nospmb_02_score_probability",
            role="exact_probability_identity",
            readiness="ready_to_apply",
            claim="The negative logarithmic slope defines a probability law on the positive half-line.",
            formula=exact["score_probability"],
            proof_boundary="Uses strict decrease and the endpoint values of f_t.",
        ),
        GateRow(
            id="nospmb_03_zero_free_lift",
            role="exact_complex_identity",
            readiness="ready_to_apply",
            claim="The one-sided Fourier lift is an exact characteristic-function defect and never vanishes on the positive real axis.",
            formula=f"{exact['complex_lift']}; {exact['cartesian_lift']}",
            proof_boundary="Real x>0; strict positivity uses the continuous nondegenerate score law.",
        ),
        GateRow(
            id="nospmb_04_phase_contact",
            role="exact_contact_equivalence",
            readiness="ready_to_apply",
            claim="A Newman double contact is exactly a stationary pi/2 crossing of the globally defined one-sided phase.",
            formula=exact["phase_contact"],
            proof_boundary="No exceptional F_t=0 case remains.",
        ),
        GateRow(
            id="nospmb_05_scaled_first_jet",
            role="exact_target_reduction",
            readiness="conditional_ready",
            claim="The full first-jet target is an exact two-observable small-ball problem for the score characteristic function.",
            formula=exact["scaled_first_jet"],
            proof_boundary="The displayed identity is exact; positivity on the open outer region remains unproved.",
        ),
        GateRow(
            id="nospmb_06_signed_hankel_bridge",
            role="exact_coefficient_bridge",
            readiness="ready_to_apply",
            claim="The signed-Hankel coefficients and every shifted Jensen window are odd-moment/Laguerre mixtures of the same score probability.",
            formula=(
                f"{exact['coefficient_bridge']}; "
                f"{exact['shifted_laguerre_bridge']}"
            ),
            proof_boundary=(
                "This identification does not prove mixture preservation, "
                "PF-infinity, or sign regularity."
            ),
        ),
        GateRow(
            id="nospmb_06a_score_beta_abel_measure",
            role="exact_abel_measure_bridge",
            readiness="ready_to_apply",
            claim="The score/Beta mixture is a single positive Abel measure whose moments, Bessel transform, shifted Laguerre averages, and Newman-time flow are all exact.",
            formula=(
                f"{exact['score_beta_abel_measure']}; "
                f"{exact['abel_bessel_laguerre_bridge']}; "
                f"{exact['abel_flow_concentration_guard']}"
            ),
            proof_boundary=(
                "Measure positivity and the closed flow do not prove the "
                "degree-two concentration threshold at all times, the "
                "higher-degree Jensen inequalities, or PF-infinity."
            ),
        ),
        GateRow(
            id="nospmb_06b_full_scale_interlacing_guard",
            role="exact_interlacing_countergate",
            readiness="guard_validated",
            claim="Full support of the Abel measure rules out a common interlacer for the entire fixed-scale Laguerre family in every degree at least two.",
            formula=exact["full_scale_interlacing_guard"],
            proof_boundary=(
                "This rejects the direct global common-interlacing shortcut, "
                "not a weighted Xi-specific total-positive connection."
            ),
        ),
        GateRow(
            id="nospmb_06c_abel_kernel_sign_guard",
            role="exact_abel_kernel_countergate",
            readiness="guard_validated",
            claim="The fractional Abel kernel underlying the score pushforward has both signs of two-by-two minors.",
            formula=exact["abel_kernel_sign_guard"],
            proof_boundary=(
                "This rejects total positivity or sign regularity of the bare "
                "Abel operator, not a weighted Xi-specific composition."
            ),
        ),
        GateRow(
            id="nospmb_07_score_flow",
            role="exact_dynamical_identity",
            readiness="ready_to_apply",
            claim="The score probability and its sine characteristic obey a closed tail-source/radial backward-heat evolution.",
            formula=exact["score_flow"],
            proof_boundary=(
                "The closed evolution preserves the exact coordinate but does "
                "not itself prevent a radial-heat collision."
            ),
        ),
        GateRow(
            id="nospmb_08_normalizer_transfer",
            role="exact_quantitative_transfer",
            readiness="ready_to_apply",
            claim="After retaining the lift amplitude, the score first jet and the normalized Newman first jet are uniformly equivalent on the corrected overlap.",
            formula=exact["normalizer_transfer"],
            proof_boundary=(
                "Imports the proved normalizer derivative bound; it does not "
                "supply a lower bound for either first jet."
            ),
        ),
        GateRow(
            id="nospmb_09_phase_conditioning",
            role="exact_asymptotic_guard",
            readiness="guard_validated",
            claim="The exact zero-free lift is asymptotically flat in phase and has a different quantitative conditioning from the corrected Riemann-Siegel lift.",
            formula=f"{exact['phase_conditioning']}; {exact['conditioning_guard']}",
            proof_boundary="Uniform on the compact time interval; it rejects only lift-independent phase-speed floors.",
        ),
        GateRow(
            id="nospmb_10_shape_countermodel",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Probability, monotonicity, and strong log-concavity do not exclude a score-characteristic double zero.",
            formula=(
                f"{exact['countermodel']}; "
                f"{exact['countermodel_quadratic_coefficient']}"
            ),
            proof_boundary="Generic shape countermodel, not the Xi kernel.",
        ),
        GateRow(
            id="nospmb_11_xi_diagnostics",
            role="finite_diagnostics",
            readiness="diagnostic_only",
            claim="Independent quadrature checks the lift, characteristic, first-jet, and phase formulas, including at the first Xi crossing.",
            formula="F_t=i*f_t(0)*(1-chi_t)/x and Theta_t'=Im(-chi_t'/(1-chi_t))",
            proof_boundary="Three finite rows and no cofinal conclusion.",
            diagnostics={
                "convergence": convergence,
                "rows": fine["rows"],
            },
        ),
        GateRow(
            id="nospmb_12_live_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="The remaining theorem is Xi-specific joint avoidance for the sine characteristic observable and its derivative.",
            formula=exact["live_target"],
            proof_boundary="This is still equivalent to the missing positive-time transversality endgame.",
        ),
        GateRow(
            id="nospmb_13_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The one-sided phase bridge and its coefficient identification do not prove the required joint avoidance.",
            formula="exact lift plus score probability plus finite diagnostics != Lambda<=0 or RH",
            proof_boundary="No Clay-prize conclusion is asserted.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_one_sided_phase_moment_bridge_gate",
        "date": "2026-07-25",
        "status": (
            "exact one-sided zero-free phase lift, score-probability and "
            "signed-Hankel/Abel moment bridge, with conditioning and shape "
            "guards plus full-scale common-interlacing and bare-Abel sign "
            "rejections; the Xi joint-avoidance theorem remains open"
        ),
        "proof_boundary": (
            "This artifact proves exact probability, complex-lift, phase, "
            "first-jet, coefficient, Abel-measure, flow, concentration-"
            "threshold, and asymptotic identities. It does not prove the Xi "
            "score-characteristic joint-avoidance target, corrected Riemann-"
            "Siegel transversality, PF-infinity, Lambda<=0, or RH."
        ),
        "exact": exact,
        "diagnostics": fine,
        "convergence": convergence,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/jensen_window_pf_newman_positive_time_strong_logconcavity_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_wronskian_phase_reduction.md",
            "outputs/jensen_window_pf_newman_theta_forward_sqrt_to_corrected_rs_C1_transfer_gate.md",
            "outputs/signed_hankel_jensen_audit.md",
            "outputs/jensen_window_pf_newman_theta_curvature_probability_operator_gate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    diagnostics_payload = artifact["diagnostics"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman One-Sided Phase/Moment Bridge Gate",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact zero-free one-sided phase lift and score-probability",
            "bridge to the signed-Hankel coefficients and their positive Abel",
            "measure. The Xi joint-avoidance theorem remains open;",
            "this is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_one_sided_phase_moment_bridge_gate.py",
            "```",
            "",
            "## Score Probability",
            "",
            "Import the proved positive-time kernel shape:",
            "",
            "```text",
            exact["imported_shape"],
            exact["score_probability"],
            "```",
            "",
            "The probability is not an auxiliary model. It is obtained directly",
            "from the logarithmic slope of the positive half-kernel.",
            "",
            "## Zero-Free Complex Lift",
            "",
            "Integration by parts gives",
            "",
            "```text",
            exact["complex_lift"],
            exact["cartesian_lift"],
            "```",
            "",
            "Because the score law has a continuous positive density, `C_t(x)<1`",
            "for every `x>0`. Hence `Y_t(x)>0`, so this exact complex lift never",
            "vanishes on the positive real axis.",
            "",
            "## Exact Contact System",
            "",
            "```text",
            exact["phase_contact"],
            exact["scaled_first_jet"],
            "```",
            "",
            "This removes the exceptional `E=0` branch from the phase reduction.",
            "The surviving obligation is the joint small-ball problem for the",
            "sine characteristic component and its derivative.",
            "",
            "## Signed-Hankel Coordinate",
            "",
            "The same probability law encodes the coefficient workstream:",
            "",
            "```text",
            exact["coefficient_bridge"],
            exact["shifted_laguerre_bridge"],
            "```",
            "",
            "Thus the phase-critical and signed-Hankel programmes are two",
            "coordinates on the same score law. The identity itself does not",
            "promote the positive Laguerre scale mixture or the finite signed",
            "minors to PF-infinity.",
            "",
            "## Score-Beta Abel Measure",
            "",
            "The two-variable score/Beta mixture can be pushed to one",
            "positive measure on the Laguerre scale:",
            "",
            "```text",
            exact["score_beta_abel_measure"],
            exact["abel_bessel_laguerre_bridge"],
            exact["abel_flow_concentration_guard"],
            "```",
            "",
            "This makes the first Jensen obstruction a sharp concentration",
            "threshold for the tilted Abel measure. Ordinary moment positivity",
            "supplies only Cauchy-Schwarz log-convexity, so it cannot be used",
            "as the missing all-degree preservation theorem.",
            "",
            "## Interlacing Guard",
            "",
            "```text",
            exact["full_scale_interlacing_guard"],
            "```",
            "",
            "The support geometry therefore retires global common interlacing",
            "of the component family. It does not rule out a theorem using the",
            "specific Abel weights and a stronger total-positive connection.",
            "",
            "The bare Abel operator cannot supply that connection by itself:",
            "",
            "```text",
            exact["abel_kernel_sign_guard"],
            "```",
            "",
            "## Closed Score Flow",
            "",
            "The probability coordinate remains closed under Newman time:",
            "",
            "```text",
            exact["score_flow"],
            "```",
            "",
            "This is an exact radial backward-heat equation for the sine",
            "observable. It does not by itself prevent collision; the missing",
            "input remains Xi-specific sign regularity or arithmetic phase",
            "separation.",
            "",
            "## Normalizer-Compatible Jet",
            "",
            "The score coordinate is quantitatively compatible with the exact",
            "normalization used by the corrected Riemann-Siegel branch:",
            "",
            "```text",
            exact["normalizer_transfer"],
            "```",
            "",
            "Thus the asymptotically flat raw phase is a conditioning effect,",
            "not a new obstruction after the lift amplitude is restored. The",
            "existing corrected lower-bound target remains unchanged.",
            "",
            "## Conditioning Guard",
            "",
            "Repeated integration by parts gives",
            "",
            "```text",
            exact["phase_conditioning"],
            "```",
            "",
            exact["conditioning_guard"],
            "",
            "The corrected Riemann-Siegel phase and this exact phase therefore",
            "cannot share an absolute velocity threshold without an explicit",
            "amplitude conversion.",
            "",
            "## Shape Countermodel",
            "",
            "```text",
            exact["countermodel"],
            exact["countermodel_quadratic_coefficient"],
            "```",
            "",
            "So probability positivity, monotonicity, and strong log-concavity",
            "still do not prove transversality.",
            "",
            "## Numerical Crosscheck",
            "",
            "| label | t | x | Re F | Im F | B | phase velocity |",
            "|---|---:|---:|---:|---:|---:|---:|",
            *[
                "| {label} | {t} | {x} | {F_real} | {F_imag} | {B} | {phase_velocity} |".format(
                    **row
                )
                for row in diagnostics_payload["rows"]
            ],
            "",
            "The quadrature checks the exact lift and first-jet identities at",
            "two endpoint rows, including the first Xi crossing, and one",
            "positive-time row. These are diagnostics only.",
            "",
            "## Live Target",
            "",
            "```text",
            exact["live_target"],
            "```",
            "",
            "A successful next theorem must use the theta arithmetic of this",
            "specific score law or close the all-order PF/sign-regular route.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman one-sided phase/moment bridge gate: "
        f"{len(artifact['rows'])} rows, 10 exact identities/bridges, "
        "1 conditioning guard, 1 full-scale interlacing guard, "
        "1 Abel-kernel sign guard, "
        "1 exact shape countermodel, "
        "3 precision-stable diagnostics, 1 open Xi joint-avoidance target"
    )


if __name__ == "__main__":
    main()
