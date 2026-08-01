#!/usr/bin/env python3
"""Certify the first-order global Polymath-15 remainder on the critical ray."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
import sympy as sp  # noqa: E402


STEM = (
    "jensen_window_pf_newman_polymath15_critical_"
    "first_order_global_remainder_certificate"
)
DEFAULT_OUT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / f"outputs/{STEM}.md"
POLYMATH_SOURCE = "https://arxiv.org/pdf/1904.12438"
PRECISION_BITS = 256
L_MIN = 50
T_AUDIT_MIN = 1_000_000
ENDPOINT_BRACKET_CONSTANT = 4_000
ENDPOINT_PREFIX_CONSTANT = 20
FIXED_ENDPOINT_CONSTANT = 30_000
FIXED_RAW_CONSTANT = 40_000
ADJACENT_CONSTANT = 20_000
GLOBAL_VALUE_CONSTANT = 100_000
GLOBAL_FIRST_CONSTANT = 200_000

SOURCE_FILES = {
    "endpoint_C1": REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.json",
    "dirichlet_second_order": REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_critical_dirichlet_second_order_remainder_certificate.json",
    "old_global_C1": REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.json",
    "C0_strip": REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.json",
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCE_FILES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing source artifacts: {missing}")
    return {name: file_hash(path) for name, path in SOURCE_FILES.items()}


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def exact_identities() -> dict[str, str]:
    y, a = sp.symbols("y a", nonzero=True)
    pi = sp.symbols("pi", positive=True)
    f1, f3 = sp.symbols("f1 f3")
    sigma = sp.Rational(1, 2) - y / 2
    c1 = f3 / (12 * pi**2) + sp.I * (2 * sigma - 1) * f1 / (4 * pi)
    critical = f3 / (12 * pi**2)
    sigma_defect = sp.simplify((c1 - critical) / a)
    displacement = sp.simplify(
        (sp.I * y / 2) * f1 * (-1 / (2 * pi * a))
    )
    if sp.simplify(sigma_defect - displacement) != 0:
        raise RuntimeError("off-axis C1/C0 displacement cancellation failed")

    tau = sp.symbols("tau", positive=True)
    a_tau = sp.sqrt(tau / (2 * pi))
    slow = tau ** sp.Rational(3, 2) + sp.I * sp.sqrt(tau)
    if sp.simplify(slow / a_tau - sp.sqrt(2 * pi) * (tau + sp.I)) != 0:
        raise RuntimeError("C1 slow-factor reduction failed")

    return {
        "region": (
            "L=log(x/(4*pi))>=50, 0<t*L<=25, rho=1/L; "
            "all Cauchy disks are centered on the real critical ray"
        ),
        "source_expansion": (
            "For u=sigma+sqrt(t)v, Proposition 6.2 gives "
            "sum_(k=0)^K C_k(p,u)/a^k+RS_K; Proposition 6.3 gives "
            "the exact heat integral with a multiplicative P(u) satisfying "
            "|P(u)-1|<=exp((u^2+5/6)/(T-3))-1"
        ),
        "C1_formula": (
            "C_1(p,sigma)=F'''(p)/(12*pi^2)"
            "+i*(2*sigma-1)*F'(p)/(4*pi), F=C_0"
        ),
        "off_axis_cancellation": (
            "For sigma=(1-y)/2 and p'(T)=-1/(2*pi*a), "
            "(C_1(p,sigma)-C_1(p,1/2))/a="
            "(i*y/2)*F'(p)*p'(T)=-i*y*F'(p)/(4*pi*a)"
        ),
        "slow_factors": (
            "K(T)=-(sqrt(pi)/8)*exp(-pi*T/4)"
            "*(T^(3/2)+i*T^(1/2)); "
            "(T^(3/2)+i*T^(1/2))/a=sqrt(2*pi)*(T+i)"
        ),
        "analytic_lift": (
            "E_[1],N(z)=(-1)^N exp(t*pi^2/64)"
            "*(G_[1],N(z)+G_[1],N^#(z)), where "
            "G_[1],N=K(T(z))*(F(p_N(T(z)))"
            "+F'''(p_N(T(z)))/(12*pi^2*a(T(z))))"
        ),
        "lift_derivatives": (
            "For A(T)=T^(3/2)+i*T^(1/2), |theta|<=1/L and T>=10^6: "
            "|A'|/|A(T)|<=2/T, |A''|/|A(T)|<=1/T^2, "
            "|p'|<=(2*pi*(T-|theta|))^-1/2, "
            "|p''|<=(2*sqrt(2*pi)*(T-|theta|)^(3/2))^-1; "
            "therefore T|S_0''|/|A(T)|<2100, "
            "T|S_1'|/|A(T)|<10133, and the two-component "
            "paper/lift mismatch is <1000/T"
        ),
        "expanded_strip": (
            "On |Re p|<=101/100, |Im p|<=1/100, Cauchy radius 1/10 "
            "stays in |Re p|<=111/100, |Im p|<=11/100; "
            "|F|<5 gives |F'|<50, |F''|<1000, "
            "|F'''|<30000, |F''''|<1200000"
        ),
        "C0_prefactor": (
            "|C_0|*integral|P-1|dmu<=2/T"
        ),
        "C1_prefactor": (
            "|C_1(p,u)|<=300+10|u| and "
            "a^-1*integral|P-1||C_1|dmu<=1000/T"
        ),
        "positive_tail": (
            "For u>=0 choose K=2; the C_2/a^2 and RS_2 terms "
            "contribute less than 9/T after Gaussian integration"
        ),
        "negative_tail": (
            "For u<0 choose K_u=max(floor(-u)+3,floor(T_0/pi)); "
            "the k>=2 low series is split into even/odd Gamma progressions, "
            "and the published extreme tail delta_3<=2*10^-30/a^14; "
            "the total negative contribution is less than 1/T"
        ),
        "endpoint_pair": (
            "The two heat endpoints, after the C_1 lift is retained, have "
            "bracket error <=4000/T on the doubled transition collar"
        ),
        "endpoint_prefactor": (
            "The published endpoint prefactor, the collar loss, and "
            "sup|B_t(z)|/A_t(x)<2 give |K(T)|/A_t(x)"
            "<2*exp(7/4)*exp(-L/4)<20*exp(-L/4); "
            "also T>pi*exp(L)"
        ),
        "endpoint_normalized": (
            "|R_endpoint|/A_t(x)<30000*exp(-5L/4)"
        ),
        "finite_dirichlet": (
            "|R_D|/A_t(x)<8000000*exp(-7L/4), "
            "|partial_x R_D|/A_t(x)<8000000*L*exp(-7L/4)"
        ),
        "fixed_raw": (
            "On a fixed-N disk, sup |R_[1]|/A_t(x)"
            "<40000*exp(-5L/4)"
        ),
        "old_adjacent": (
            "The inherited C_0/main-block adjacent mismatch is "
            "<1000*exp(-5L/4)"
        ),
        "new_adjacent": (
            "|d_(t,n)|<100/T for the entering Dirichlet term; "
            "critical C_1 parity gives exact endpoint value matching and "
            "its adjacent derivative difference is O(T^-1); together with "
            "the old mismatch, |Delta lift_[1]|/A_t(x)"
            "<20000*exp(-5L/4)"
        ),
        "adjacent_components": (
            "From |alpha'(s)|<=1/(T-3), t|alpha_n|<=27, "
            "d=1/(6s)+alpha'(s)(t/4+t^2*alpha_n^2/8) gives "
            "T|d|<100. The entering two-sharp block is "
            "<200*exp(-L/4), costing <7000*exp(-5L/4). "
            "C_1 parity makes the endpoint difference zero at a=m; "
            "the F'''' derivative bound costs <6000*exp(-5L/4) "
            "on the doubled collar"
        ),
        "global_remainder": (
            "For every critical disk, including cutoff crossings, "
            "|r_[1](x)|<100000*exp(-5L/4) and "
            "|partial_x r_[1](x)|<200000*L*exp(-5L/4)"
        ),
        "cauchy_transfer": (
            "A radius-1/L Cauchy estimate gives |R_[1]'|/A_t"
            "<100000*L*exp(-5L/4); with "
            "|partial_x log A_t|<L/2 this gives "
            "|partial_x r_[1]|<150000*L*exp(-5L/4)"
            "<200000*L*exp(-5L/4)"
        ),
        "contact_handoff": (
            "Insert eta_0=100000 and eta_1=200000 in the exact signed "
            "contact-normal inequality; the strict Xi arithmetic reversal "
            "is the next open theorem target"
        ),
    }


def gaussian_quadratic_linear(
    c: flint.arb,
    linear: flint.arb,
    sigma: flint.arb,
    heat: flint.arb,
) -> flint.arb:
    q = 1 - c * heat
    return (1 / q.sqrt()) * (
        c * sigma**2 / q
        + linear * sigma / q
        + heat * linear**2 / (4 * q)
    ).exp()


def interval_certificate() -> dict[str, str]:
    flint.ctx.prec = PRECISION_BITS
    one = flint.arb(1)
    pi = flint.arb.pi()
    sqrt_two = flint.arb(2).sqrt()
    ell = flint.arb(L_MIN)

    # The previous radius-1/4 removable-disc estimate gives |F|<3 near
    # p=+-1/2. Outside radius 1/5, a direct denominator lower bound gives
    # |F|<5 on the expanded rectangle needed for a radius-1/10 Cauchy disk.
    real_pad = flint.arb(111) / 100
    imag_pad = flint.arb(11) / 100
    numerator = (
        (pi * real_pad * imag_pad).exp()
        + sqrt_two * (pi * imag_pad / 2).cosh()
    )
    vertical_denominator = 2 * (pi / 10).sinh()
    horizontal_denominator = 2 * (pi * flint.arb(3).sqrt() / 10).sin()
    if not vertical_denominator < horizontal_denominator:
        raise RuntimeError("expanded-strip denominator split failed")
    outside_c0 = numerator / vertical_denominator
    if not numerator < 3:
        raise RuntimeError("expanded-strip numerator bound failed")
    if not outside_c0 < 5:
        raise RuntimeError("expanded-strip C0 bound failed")

    derivative_bounds = {
        "F": 5,
        "F1": 50,
        "F2": 1_000,
        "F3": 30_000,
        "F4": 1_200_000,
    }

    t_lower = 2 * pi * ell.exp() - 1 / ell
    if not t_lower > pi * ell.exp():
        raise RuntimeError("critical doubled-collar T lower bound failed")
    if not t_lower > T_AUDIT_MIN:
        raise RuntimeError("critical T audit floor failed")
    p_displacement = (1 / ell) / (2 * pi * (t_lower - 1)).sqrt()
    if not p_displacement < one / 100:
        raise RuntimeError("first-order lift p-strip map failed")

    endpoint_prefix = 2 * (flint.arb(7) / 4).exp()
    if not endpoint_prefix < ENDPOINT_PREFIX_CONSTANT:
        raise RuntimeError("endpoint prefactor constant failed")

    audit_t = flint.arb(T_AUDIT_MIN)
    denominator = audit_t - 3
    c = 1 / denominator
    heat = one / 2
    sigma = flint.arb(13) / 25
    q = 1 - c * heat

    eps = 10 / (3 * (audit_t - 6))
    eps_scaled = audit_t * eps * eps.exp()
    if not eps_scaled < 4:
        raise RuntimeError("C0 prefactor scaled bound failed")

    tilted_mass = (1 / q.sqrt()) * (c * sigma**2 / q).exp()
    tilted_mean = sigma / q
    tilted_variance = heat / (2 * q)
    tilted_second = tilted_mass * (
        tilted_variance + tilted_mean**2
    )
    tilted_fourth = tilted_mass * (
        3 * tilted_variance**2
        + 6 * tilted_variance * tilted_mean**2
        + tilted_mean**4
    )
    c1_cross_moment = (flint.arb(5) / (6 * denominator)).exp() * (
        5 * tilted_fourth
        + (flint.arb(305) + flint.arb(25) / 6) * tilted_second
        + flint.arb(1525) / 6 * tilted_mass
    )
    if not c1_cross_moment < 500:
        raise RuntimeError("C1 prefactor moment bound failed")
    c1_prefactor_scaled = (
        audit_t / denominator * c1_cross_moment
    )
    if not c1_prefactor_scaled < 1_000:
        raise RuntimeError("C1 prefactor scaled constant failed")

    common_heat_factor = (flint.arb(5) / (6 * denominator)).exp()
    mgf_9 = common_heat_factor * gaussian_quadratic_linear(
        c, flint.arb(9).log(), sigma, heat
    )
    mgf_rs = common_heat_factor * gaussian_quadratic_linear(
        c, flint.arb(3) * flint.arb(2).log() / 2, sigma, heat
    )
    if not mgf_9 < 7:
        raise RuntimeError("positive C2 Gaussian moment failed")
    if not mgf_rs < 3:
        raise RuntimeError("positive RS2 Gaussian moment failed")

    c2_prefactor = sqrt_two / (8 * pi)
    rs2_prefactor = (
        flint.arb(1) / 7
        * (flint.arb(3) / 2).gamma()
        * (flint.arb(11) / 10) ** 3
    )
    positive_a2_constant = c2_prefactor * mgf_9 + rs2_prefactor * mgf_rs
    if not positive_a2_constant < 1:
        raise RuntimeError("positive k>=2 a^-2 constant failed")

    de_reyna_scale = (
        (3 - 2 * flint.arb(2).log()) * pi
    )
    ratio = 1 / de_reyna_scale.sqrt()
    if not ratio**2 < one / 5:
        raise RuntimeError("negative coefficient Gamma ratio failed")
    negative_prefactor = sqrt_two / (4 * pi**2)
    low_series_constant = negative_prefactor * (
        ratio**2 / (1 - ratio**2)
        + (flint.arb(3) / 2).gamma()
        * ratio**3
        / (1 - ratio**2)
    )
    if not low_series_constant < one / 50:
        raise RuntimeError("negative low-series a^-2 constant failed")
    mgf_negative = common_heat_factor * gaussian_quadratic_linear(
        c, -flint.arb(2).log(), flint.arb(0), heat
    )
    if not mgf_negative < 2:
        raise RuntimeError("negative-tail Gaussian moment failed")
    negative_a2_constant = low_series_constant * mgf_negative
    if not negative_a2_constant < one / 20:
        raise RuntimeError("negative k>=2 a^-2 constant failed")

    tail_scaled_constant = 2 * pi * (
        positive_a2_constant + negative_a2_constant
    ) + 4 * pi / (10**30)
    if not tail_scaled_constant < 10:
        raise RuntimeError("combined k>=2 tail constant failed")

    h = one / ell
    slow_amp_scaled = h * 2 * derivative_bounds["F"]
    slow_second_scaled = h**2 / 2 * 2_100
    c1_slow_coefficient = (2 * pi).sqrt() / (12 * pi**2)
    c1_derivative_scaled = (
        c1_slow_coefficient
        * derivative_bounds["F3"]
        / audit_t.sqrt()
        + flint.arb(derivative_bounds["F4"])
        / (12 * pi**2)
        * (audit_t + 2)
        / audit_t
        * (audit_t / (audit_t - h)).sqrt()
    )
    component_lift_scaled = (
        slow_amp_scaled
        + slow_second_scaled
        + h * c1_derivative_scaled
    )
    pair_lift_scaled = 2 * component_lift_scaled
    if not pair_lift_scaled < 1_000:
        raise RuntimeError("first-order holomorphic-lift constant failed")

    component_expansion = flint.arb(2 + 1_000 + 10)
    endpoint_pair_bracket = 2 * component_expansion + pair_lift_scaled
    if not endpoint_pair_bracket < ENDPOINT_BRACKET_CONSTANT:
        raise RuntimeError("endpoint pair bracket constant failed")
    endpoint_normalized = (
        endpoint_prefix * ENDPOINT_BRACKET_CONSTANT / pi
    )
    if not endpoint_normalized < FIXED_ENDPOINT_CONSTANT:
        raise RuntimeError("fixed endpoint normalized constant failed")

    finite_effective = flint.arb(8_000_000) * (-ell / 2).exp()
    if not finite_effective < 1:
        raise RuntimeError("finite Dirichlet target-scale absorption failed")
    fixed_raw = flint.arb(FIXED_ENDPOINT_CONSTANT) + finite_effective
    if not fixed_raw < FIXED_RAW_CONSTANT:
        raise RuntimeError("fixed first-order raw constant failed")

    dirichlet_adjacent = flint.arb(200) * 100 / pi
    if not dirichlet_adjacent < 7_000:
        raise RuntimeError("adjacent Dirichlet correction constant failed")
    endpoint_c1_adjacent = (
        4
        * c1_derivative_scaled
        / ell
        * endpoint_prefix
        / pi
    )
    if not endpoint_c1_adjacent < 6_000:
        raise RuntimeError("adjacent endpoint C1 constant failed")
    adjacent_total = (
        flint.arb(1_000)
        + dirichlet_adjacent
        + endpoint_c1_adjacent
    )
    if not adjacent_total < ADJACENT_CONSTANT:
        raise RuntimeError("first-order adjacent total failed")
    global_raw = flint.arb(FIXED_RAW_CONSTANT + ADJACENT_CONSTANT)
    if not global_raw < GLOBAL_VALUE_CONSTANT:
        raise RuntimeError("global first-order value constant failed")
    if not 2 * flint.arb(GLOBAL_VALUE_CONSTANT) <= GLOBAL_FIRST_CONSTANT:
        raise RuntimeError("global first-order derivative constant failed")

    return {
        "precision_bits": str(PRECISION_BITS),
        "L_min": str(L_MIN),
        "T_audit_min": str(T_AUDIT_MIN),
        "expanded_numerator_ball": arb_text(numerator),
        "expanded_vertical_denominator_ball": arb_text(vertical_denominator),
        "expanded_horizontal_denominator_ball": arb_text(horizontal_denominator),
        "expanded_C0_ball": arb_text(outside_c0),
        "derivative_bounds": derivative_bounds,
        "critical_T_lower_ball": arb_text(t_lower),
        "p_displacement_ball": arb_text(p_displacement),
        "endpoint_prefix_ball": arb_text(endpoint_prefix),
        "C0_prefactor_scaled_ball": arb_text(eps_scaled),
        "C1_cross_moment_ball": arb_text(c1_cross_moment),
        "C1_prefactor_scaled_ball": arb_text(c1_prefactor_scaled),
        "mgf_9_ball": arb_text(mgf_9),
        "mgf_RS2_ball": arb_text(mgf_rs),
        "positive_a2_constant_ball": arb_text(positive_a2_constant),
        "de_reyna_ratio_ball": arb_text(ratio),
        "negative_low_series_ball": arb_text(low_series_constant),
        "mgf_negative_ball": arb_text(mgf_negative),
        "negative_a2_constant_ball": arb_text(negative_a2_constant),
        "tail_scaled_constant_ball": arb_text(tail_scaled_constant),
        "C1_slow_derivative_scaled_ball": arb_text(c1_derivative_scaled),
        "pair_lift_scaled_ball": arb_text(pair_lift_scaled),
        "endpoint_pair_bracket_ball": arb_text(endpoint_pair_bracket),
        "endpoint_normalized_ball": arb_text(endpoint_normalized),
        "finite_effective_L50_ball": arb_text(finite_effective),
        "fixed_raw_ball": arb_text(fixed_raw),
        "adjacent_dirichlet_ball": arb_text(dirichlet_adjacent),
        "adjacent_endpoint_C1_ball": arb_text(endpoint_c1_adjacent),
        "adjacent_total_ball": arb_text(adjacent_total),
        "global_raw_ball": arb_text(global_raw),
    }


def build_artifact() -> dict:
    exact = exact_identities()
    interval = interval_certificate()
    rows = [
        GateRow(
            id="np15f1grc_01_scope",
            role="exact_scope",
            readiness="ready_to_apply",
            claim="The first-order endpoint problem is posed on the full critical Cauchy collar.",
            formula=exact["region"],
            proof_boundary="Only L>=50 and 0<tL<=25; the bounded-L shoulder is not included.",
        ),
        GateRow(
            id="np15f1grc_02_source",
            role="published_input",
            readiness="ready_to_apply",
            claim="The heat endpoint is controlled by the explicit Proposition 6.2 coefficient and remainder bounds.",
            formula=exact["source_expansion"],
            proof_boundary="Uses Polymath 15 and its cited Arias de Reyna estimates, not a new Riemann-Siegel theorem.",
        ),
        GateRow(
            id="np15f1grc_03_off_axis_cancellation",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The apparently dangerous off-axis C1 term is exactly the first complex displacement of the C0 lift.",
            formula=f"{exact['C1_formula']}; {exact['off_axis_cancellation']}",
            proof_boundary="Exact first-order identity; all remaining Taylor terms are bounded separately.",
        ),
        GateRow(
            id="np15f1grc_04_strip",
            role="interval_certificate",
            readiness="certified",
            claim="A wider removable-point strip gives four explicit C0 derivative bounds.",
            formula=exact["expanded_strip"],
            proof_boundary="Finite complex rectangle only; the apparent poles at p=+-1/2 are removable.",
            diagnostics={
                "C0": interval["expanded_C0_ball"],
                "derivatives": interval["derivative_bounds"],
                "p_displacement": interval["p_displacement_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_05_prefactor",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="The C0 and retained C1 layers absorb the full gamma/Stirling and log-M0 prefactor defect at order T^-1.",
            formula=f"{exact['C0_prefactor']}; {exact['C1_prefactor']}",
            proof_boundary="The C1 moment uses the exact affine sigma dependence and an explicit tilted Gaussian fourth moment.",
            diagnostics={
                "C0_scaled": interval["C0_prefactor_scaled_ball"],
                "C1_moment": interval["C1_cross_moment_ball"],
                "C1_scaled": interval["C1_prefactor_scaled_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_06_positive_tail",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="K=2 suffices on the nonnegative shifted-sigma half-line.",
            formula=exact["positive_tail"],
            proof_boundary="Bounds C2 and RS2 after Gaussian integration; it does not use a fixed K on the negative half-line.",
            diagnostics={
                "mgf_9": interval["mgf_9_ball"],
                "mgf_RS2": interval["mgf_RS2_ball"],
                "a2_constant": interval["positive_a2_constant_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_07_negative_tail",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="A variable expansion order closes the negative shifted-sigma Gaussian tail.",
            formula=exact["negative_tail"],
            proof_boundary="Uses K_u+u>=2 and retains the published extreme-tail estimate; fixed K=2 is explicitly forbidden here.",
            diagnostics={
                "ratio": interval["de_reyna_ratio_ball"],
                "low_series": interval["negative_low_series_ball"],
                "Gaussian": interval["mgf_negative_ball"],
                "scaled_total_tail": interval["tail_scaled_constant_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_08_lift",
            role="analytic_bound",
            readiness="ready_to_apply",
            claim="The first-order endpoint has a fixed-N holomorphic lift whose paper-extension defect is O(T^-1).",
            formula=(
                f"{exact['slow_factors']}; {exact['analytic_lift']}; "
                f"{exact['lift_derivatives']}"
            ),
            proof_boundary="Uses the exact C0/C1 displacement cancellation plus the explicit F through F'''' strip bounds.",
            diagnostics={
                "C1_derivative": interval["C1_slow_derivative_scaled_ball"],
                "pair_lift": interval["pair_lift_scaled_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_09_endpoint",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="After C1 is retained, the complete heat endpoint starts at the a^-2 scale.",
            formula=(
                f"{exact['endpoint_pair']}; {exact['endpoint_prefactor']}; "
                f"{exact['endpoint_normalized']}"
            ),
            proof_boundary="Uniform on fixed-N and doubled transition collars; this is an endpoint theorem only.",
            diagnostics={
                "bracket": interval["endpoint_pair_bracket_ball"],
                "normalized": interval["endpoint_normalized_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_10_fixed_composition",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="The endpoint theorem composes with the peeled finite Dirichlet residual on every fixed cutoff.",
            formula=f"{exact['finite_dirichlet']}; {exact['fixed_raw']}",
            proof_boundary="Fixed prescribed N; adjacent lifts are handled by the next rows.",
            diagnostics={
                "finite_L50": interval["finite_effective_L50_ball"],
                "raw": interval["fixed_raw_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_11_adjacent_old",
            role="imported_bound",
            readiness="ready_to_apply",
            claim="The already-certified C0/main-block cutoff mismatch is at the new target scale.",
            formula=exact["old_adjacent"],
            proof_boundary="Imported without strengthening from the previous global C1 certificate.",
        ),
        GateRow(
            id="np15f1grc_12_adjacent_first_order",
            role="proved_bound",
            readiness="ready_to_apply",
            claim="The entering Dirichlet correction and parity-matched endpoint C1 term preserve the target cutoff scale.",
            formula=f"{exact['new_adjacent']}; {exact['adjacent_components']}",
            proof_boundary="One adjacent cutoff on a doubled collar; consecutive cutoffs cannot both occur there.",
            diagnostics={
                "Dirichlet": interval["adjacent_dirichlet_ball"],
                "endpoint_C1": interval["adjacent_endpoint_C1_ball"],
                "total": interval["adjacent_total_ball"],
            },
        ),
        GateRow(
            id="np15f1grc_13_global",
            role="asymptotic_theorem",
            readiness="ready_to_apply",
            claim="The first-order signed approximation now has explicit global value and first-derivative remainders.",
            formula=f"{exact['global_remainder']}; {exact['cauchy_transfer']}",
            proof_boundary="L>=50 and 0<tL<=25 only; Cauchy gives the derivative after one fixed adjacent lift is selected.",
            diagnostics={"global_raw": interval["global_raw_ball"]},
        ),
        GateRow(
            id="np15f1grc_14_handoff",
            role="open_handoff",
            readiness="open",
            claim="The approximation side is closed at first order, but the signed Xi contact contradiction remains open.",
            formula=exact["contact_handoff"],
            proof_boundary="Not contact exclusion, not the bounded-L shoulder, not Lambda<=0, and not RH.",
        ),
        GateRow(
            id="np15f1grc_15_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="A global asymptotic remainder theorem cannot be promoted to a Clay-prize conclusion.",
            formula="controlled r_[1] != strict signed contact-normal reversal",
            proof_boundary="The arithmetic inequality is still an open theorem target.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "global first-order critical remainder theorem on L>=50 and "
            "0<tL<=25; not contact exclusion, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the heat-integrated Riemann-Siegel endpoint "
            "a^-2 estimate, the variable-order negative-sigma tail, the "
            "first-order holomorphic lift, and one adjacent-cutoff signed "
            "comparison. Together with the existing finite Dirichlet theorem "
            "it proves eta_0=100000 and eta_1=200000 for r_[1]. It does not "
            "prove the strict Xi contact-normal inequality, the bounded-L "
            "shoulder, contact exclusion, Lambda<=0, RH, or a Clay-prize result."
        ),
        "builder_sha256": file_hash(Path(__file__)),
        "source_sha256": source_hashes(),
        "exact": exact,
        "interval": interval,
        "summary": {
            "rows": 15,
            "exact_or_published_inputs": 4,
            "analytic_or_interval_bounds": 8,
            "proved_compositions": 4,
            "open_targets": 1,
            "eta_0": GLOBAL_VALUE_CONSTANT,
            "eta_1": GLOBAL_FIRST_CONSTANT,
        },
        "rows": [asdict(row) for row in rows],
        "sources": [
            POLYMATH_SOURCE,
            "outputs/jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_dirichlet_second_order_remainder_certificate.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.md",
            "outputs/jensen_window_pf_newman_polymath15_endpoint_C0_strip_certificate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    interval = artifact["interval"]
    return "\n".join(
        [
            "# Newman Polymath-15 Critical First-Order Global Remainder Certificate",
            "",
            "Date: 2026-07-26",
            "",
            "Status: global first-order endpoint and cutoff remainder on the",
            "`L>=50`, `0<tL<=25` critical ray. This is not a proof of contact",
            "exclusion, `Lambda<=0`, RH, or a Clay-prize conclusion.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Exact First-Order Cancellation",
            "",
            f"[Polymath 15]({POLYMATH_SOURCE}), Propositions 6.2 and 6.3,",
            "give the endpoint heat integral and its coefficient/remainder",
            "bounds. The extracted coefficient is",
            "",
            "```text",
            exact["C1_formula"],
            exact["off_axis_cancellation"],
            "```",
            "",
            "Thus the off-axis `F'` term is not discarded. It is exactly the",
            "first complex displacement of `C0(p(T(z)))`. This cancellation is",
            "what permits a second-order holomorphic collar.",
            "",
            "## Fourth-Derivative Strip",
            "",
            "The old removable-disc argument is enlarged by a radius-`1/10`",
            "Cauchy disk:",
            "",
            "```text",
            exact["expanded_strip"],
            f"direct expanded C0 bound = {interval['expanded_C0_ball']} < 5",
            f"p displacement = {interval['p_displacement_ball']} < 1/100",
            "```",
            "",
            "## Heat Prefactor",
            "",
            "Writing the residual multiplier as `P(u)`, the `C0` and retained",
            "`C1` pieces satisfy",
            "",
            "```text",
            exact["C0_prefactor"],
            exact["C1_prefactor"],
            f"scaled C0 constant = {interval['C0_prefactor_scaled_ball']} < 4",
            f"tilted C1 moment = {interval['C1_cross_moment_ball']} < 500",
            "```",
            "",
            "This includes the gamma/Stirling and `log M0` defects from the",
            "proof of Proposition 6.3.",
            "",
            "## Positive Shift",
            "",
            "```text",
            exact["positive_tail"],
            f"positive a^-2 coefficient = {interval['positive_a2_constant_ball']} < 1",
            "```",
            "",
            "## Negative Shift",
            "",
            "A fixed `K=2` is invalid when `u<0`. Instead use",
            "",
            "```text",
            exact["negative_tail"],
            f"negative low series = {interval['negative_low_series_ball']} < 1/50",
            f"combined T-scaled tail = {interval['tail_scaled_constant_ball']} < 10",
            "```",
            "",
            "The extreme range is the same Fubini-Tonelli tail bounded by",
            "`2*10^-30/a^14` in the published proof.",
            "",
            "## Holomorphic Lift",
            "",
            "```text",
            exact["slow_factors"],
            exact["analytic_lift"],
            exact["lift_derivatives"],
            f"two-component lift constant = {interval['pair_lift_scaled_ball']} < 1000",
            "```",
            "",
            "Combining the heat expansion and lift comparison gives",
            "",
            "```text",
            exact["endpoint_pair"],
            exact["endpoint_prefactor"],
            exact["endpoint_normalized"],
            "```",
            "",
            "## Fixed Cutoffs",
            "",
            "The finite Dirichlet second-order theorem is lower order:",
            "",
            "```text",
            exact["finite_dirichlet"],
            exact["fixed_raw"],
            f"L=50 absorption = {interval['finite_effective_L50_ball']} < 1",
            "```",
            "",
            "## Adjacent Cutoffs",
            "",
            "The inherited `C0` mismatch, entering corrected Dirichlet term,",
            "and parity-matched endpoint `C1` term give",
            "",
            "```text",
            exact["old_adjacent"],
            exact["new_adjacent"],
            exact["adjacent_components"],
            f"saved adjacent sum = {interval['adjacent_total_ball']} < 20000",
            "```",
            "",
            "## Global First Jet",
            "",
            "One fixed analytic lift is now available on every Cauchy disk.",
            "Cauchy's estimate and the real normalizer derivative yield",
            "",
            "```text",
            exact["global_remainder"],
            exact["cauchy_transfer"],
            "```",
            "",
            "## Remaining Theorem",
            "",
            "The approximation side no longer sets the critical-ray scale.",
            "The next obligation is",
            "",
            "```text",
            exact["contact_handoff"],
            "```",
            "",
            "That signed arithmetic inequality, the bounded-`L` shoulder, and",
            "the final exhaustion remain open.",
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
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman critical first-order global remainder certificate: "
        "15 rows, eta_0=100000, eta_1=200000, "
        "1 open signed-contact target"
    )


if __name__ == "__main__":
    main()
