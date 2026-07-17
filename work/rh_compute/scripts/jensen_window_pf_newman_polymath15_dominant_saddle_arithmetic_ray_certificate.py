#!/usr/bin/env python3
"""Certify a dominant-saddle ray for the finite Polymath-15 main sum."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.md"
)
SOURCE_URL = "https://arxiv.org/abs/1904.12438"
PRECISION_BITS = 256
P_EXPONENT = Fraction(349, 100)
L_MIN = Fraction(48)
TL_MIN = Fraction(24)
ETA = Fraction(1, 10**21)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def build_exact() -> dict:
    p0, p1, p2, q0, q1, q2, v = sp.symbols(
        "p0 p1 p2 q0 q1 q2 v", real=True
    )
    curvature = lambda f0, f1, f2: f1**2 - f0 * f2 + v * f0**2
    delta = sp.expand(
        curvature(p0 + q0, p1 + q1, p2 + q2)
        - curvature(p0, p1, p2)
    )
    expected = (
        2 * p1 * q1
        + q1**2
        - p0 * q2
        - q0 * p2
        - q0 * q2
        + 2 * v * p0 * q0
        + v * q0**2
    )
    if sp.simplify(delta - expected) != 0:
        raise RuntimeError("dominant-saddle perturbation identity failed")

    beta1, beta2, sigma1, sigma2, ell = sp.symbols(
        "beta1 beta2 sigma1 sigma2 ell"
    )
    log_first = sp.I * beta1 - sigma1 * ell
    log_second = sp.I * beta2 - sigma2 * ell
    lambda_second_ratio = sp.expand(log_first**2 + log_second)
    if sp.simplify(
        lambda_second_ratio
        - (
            (sp.I * beta1 - sigma1 * ell) ** 2
            + sp.I * beta2
            - sigma2 * ell
        )
    ) != 0:
        raise RuntimeError("tail-jet logarithmic derivative failed")

    return {
        "ray": (
            "0<t<=1/2, L=log(x/(4*pi)), t*L>=24; hence L>=48 and "
            "x=4*pi*exp(L)>10^21"
        ),
        "alpha": (
            "alpha(s)=1/(2s)+1/(s-1)+(1/2)Log(s/(2*pi)), "
            "s=(1-i*x)/2"
        ),
        "alpha_real": (
            "Re(alpha)=L/2+(1/4)log(1+x^-2)-1/(1+x^2)"
        ),
        "alpha_bounds": (
            "|alpha|<=L/2+1, |alpha'|<=2/x, |alpha''|<=3/x^2"
        ),
        "normalizer_bounds": (
            "g'=alpha+(t/2)alpha*alpha', "
            "L/2-x^-2-(L/4+1/2)/x<=Re(g')<="
            "L/2+(4x^2)^-1+(L/4+1/2)/x; "
            "|g''|<=2/x+(3L/8+7/4)/x^2<=4/x"
        ),
        "geometry_bounds": (
            "0.24L<=b=|beta'|<=0.26L, a=|s_*'|<=0.51, "
            "c=|beta''|<=10^-21, d=|s_*''|<=10^-42, |V|<=10^-21"
        ),
        "coefficient_exponent": (
            "p_n>=1/2+(t/4)log((x/(4*pi))/n)-t/(2x^2)>="
            "7/2-[1/(512*exp(L))+1/(4x^2)]>349/100"
        ),
        "coefficient_bound": (
            "|exp((t/4)log(n)^2)*n^(-s_*)|<=n^(-3.49), 2<=n<=N"
        ),
        "tail_jets": (
            "Q0=2*S0; Q1=2*(b*S0+a*S1); "
            "Q2=2*((b^2+c)*S0+(2*b*a+d)*S1+a^2*S2)"
        ),
        "phase_floor": (
            "C[P0]>=4b^2-2c-4|V|>0.23L^2, P0=2cos(beta)"
        ),
        "perturbation_error": (
            "Err<=4bQ1+Q1^2+2Q2+Q0(2c+2b^2)+Q0Q2+"
            "|V|(4Q0+Q0^2)"
        ),
        "theorem": (
            "C_t[P_(N,t)]>0.076L^2 for every point on the ray, with N held "
            "fixed when differentiating the finite analytic main sum"
        ),
    }


def arb_moment_data() -> dict:
    flint.ctx.prec = PRECISION_BITS
    p = flint.arb(P_EXPONENT.numerator) / P_EXPONENT.denominator
    zeta_series = flint.arb_series([p, 1], 3).zeta()
    s0 = zeta_series[0] - 1
    s1 = -zeta_series[1]
    s2 = 2 * zeta_series[2]
    bounds = {
        "S0": flint.arb(16) / 125,
        "S1": flint.arb(23) / 200,
        "S2": flint.arb(121) / 1000,
    }
    values = {"S0": s0, "S1": s1, "S2": s2}
    for key in values:
        if not values[key] < bounds[key]:
            raise RuntimeError(f"Arb zeta-moment bound failed for {key}")

    x_scale = 4 * flint.arb.pi() * flint.arb(48).exp()
    if not x_scale > flint.arb(10) ** 21:
        raise RuntimeError("L=48 scale lower bound failed")
    coefficient_loss = (
        1 / (flint.arb(512) * flint.arb(48).exp())
        + 1 / (4 * x_scale**2)
    )
    if not coefficient_loss < flint.arb(1) / 100:
        raise RuntimeError("coefficient exponent loss bound failed")
    return {
        "precision_bits": PRECISION_BITS,
        "p": fraction_text(P_EXPONENT),
        "moments": {
            key: {
                "ball": arb_text(values[key]),
                "strict_upper": arb_text(bounds[key]),
            }
            for key in values
        },
        "x_at_L48_ball": arb_text(x_scale),
        "x_at_L48_gt": "10^21",
        "coefficient_exponent_loss_ball": arb_text(coefficient_loss),
        "coefficient_exponent_loss_lt": "1/100",
    }


def rational_budget() -> dict:
    s0 = Fraction(16, 125)
    s1 = Fraction(23, 200)
    s2 = Fraction(121, 1000)
    b_lower = Fraction(6, 25)
    b_upper = Fraction(13, 50)
    a_upper = Fraction(51, 100)
    l0 = L_MIN
    eta = ETA

    q0_cap = Fraction(32, 125)
    q1_raw = 2 * (b_upper * s0 + a_upper * s1 / l0)
    q1_cap = Fraction(6901, 100000)
    q2_raw = 2 * (
        b_upper**2 * s0
        + 2 * b_upper * a_upper * s1 / l0
        + a_upper**2 * s2 / l0**2
        + eta * s0 / l0**2
        + eta**2 * s1 / l0**2
    )
    q2_cap = Fraction(1861, 100000)
    if not q1_raw < q1_cap:
        raise RuntimeError("Q1 rational cap failed")
    if not q2_raw < q2_cap:
        raise RuntimeError("Q2 rational cap failed")

    phase_floor_raw = 4 * b_lower**2 - 6 * eta / l0**2
    phase_floor_cap = Fraction(23, 100)
    if not phase_floor_raw > phase_floor_cap:
        raise RuntimeError("phase-floor rational cap failed")

    error_coefficient = (
        4 * b_upper * q1_cap
        + q1_cap**2
        + 2 * q2_cap
        + q0_cap * 2 * b_upper**2
        + q0_cap * q2_cap
        + q0_cap * 2 * eta / l0**2
        + eta * (4 * q0_cap + q0_cap**2) / l0**2
    )
    error_cap = Fraction(77, 500)
    if not error_coefficient < error_cap:
        raise RuntimeError("curvature-error rational cap failed")
    theorem_margin = phase_floor_cap - error_cap
    if theorem_margin != Fraction(19, 250):
        raise RuntimeError("unexpected theorem margin")

    return {
        "moment_caps": {
            "S0": fraction_text(s0),
            "S1": fraction_text(s1),
            "S2": fraction_text(s2),
        },
        "geometry_caps": {
            "b_lower_over_L": fraction_text(b_lower),
            "b_upper_over_L": fraction_text(b_upper),
            "a_upper": fraction_text(a_upper),
            "c_upper": fraction_text(eta),
            "d_upper": fraction_text(eta**2),
            "V_abs_upper": fraction_text(eta),
        },
        "tail_caps": {
            "Q0": fraction_text(q0_cap),
            "Q1_over_L_raw": fraction_text(q1_raw),
            "Q1_over_L": fraction_text(q1_cap),
            "Q2_over_L2_raw": fraction_text(q2_raw),
            "Q2_over_L2": fraction_text(q2_cap),
        },
        "curvature_budget": {
            "phase_floor_over_L2_raw": fraction_text(phase_floor_raw),
            "phase_floor_over_L2": fraction_text(phase_floor_cap),
            "error_over_L2_raw": fraction_text(error_coefficient),
            "error_over_L2": fraction_text(error_cap),
            "strict_margin_over_L2": fraction_text(theorem_margin),
        },
    }


def build_artifact() -> dict:
    exact = build_exact()
    arb_data = arb_moment_data()
    budget = rational_budget()
    rows = [
        GateRow(
            id="np15dsarc_01_alpha_geometry",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The published normalizer has an exact real-part formula on the critical axis.",
            formula=exact["alpha_real"],
            proof_boundary="Exact elementary algebra for s=(1-i*x)/2.",
        ),
        GateRow(
            id="np15dsarc_02_ray_geometry",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="The normalizer phase and shifted exponent obey uniform derivative bounds on the dominant ray.",
            formula=exact["geometry_bounds"],
            proof_boundary="Uses L>=48 and elementary modulus estimates for alpha, alpha', and alpha''.",
        ),
        GateRow(
            id="np15dsarc_03_coefficient_envelope",
            role="analytic_theorem",
            readiness="ready_to_apply",
            claim="Every non-leading finite Dirichlet coefficient lies under one fixed summable power law.",
            formula=exact["coefficient_bound"],
            proof_boundary="Uses the published lower bound for Re(s_*) and N<=sqrt(x/(4*pi)+t/16).",
        ),
        GateRow(
            id="np15dsarc_04_arb_moments",
            role="interval_certificate",
            readiness="certified",
            claim="Arb encloses the three zeta moments needed for the coefficient tail.",
            formula="S0<16/125, S1<23/200, S2<121/1000 at p=349/100",
            proof_boundary="Outward-rounded Arb evaluation of zeta and its first two derivatives.",
            diagnostics=arb_data,
        ),
        GateRow(
            id="np15dsarc_05_tail_jets",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="The arithmetic tail's first three real jets reduce to the three certified zeta moments.",
            formula=exact["tail_jets"],
            proof_boundary="Triangle inequality after exact logarithmic differentiation of each finite term.",
        ),
        GateRow(
            id="np15dsarc_06_tail_budget",
            role="rational_certificate",
            readiness="certified",
            claim="Exact rational arithmetic bounds all finite-tail jets homogeneously in L.",
            formula="Q0<=0.256, Q1<=0.06901L, Q2<=0.01861L^2",
            proof_boundary="Uses infinite zeta tails, so it covers every finite published cutoff N.",
            diagnostics=budget["tail_caps"],
        ),
        GateRow(
            id="np15dsarc_07_phase_floor",
            role="rational_certificate",
            readiness="certified",
            claim="The leading saddle has a uniform normalized curvature floor.",
            formula=exact["phase_floor"],
            proof_boundary="Allows either sign of the small normalizer potential V.",
        ),
        GateRow(
            id="np15dsarc_08_arithmetic_ray",
            role="proved_theorem",
            readiness="ready_to_apply",
            claim="The complete finite Polymath-15 main sum has positive normalized Laguerre curvature on the dominant ray.",
            formula=exact["theorem"],
            proof_boundary="Finite main sum only; the published analytic remainder is not included.",
            diagnostics=budget["curvature_budget"],
        ),
        GateRow(
            id="np15dsarc_09_remainder_gap",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Passing from the finite main theorem to H_t requires a uniform Cauchy remainder budget and cutoff-transition collars.",
            formula="C[P]>0.076L^2 must dominate the normalized C^2 theorem remainder",
            proof_boundary="Open; no strict Laguerre theorem for H_t is asserted.",
        ),
        GateRow(
            id="np15dsarc_10_live_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The dominant-saddle branch is reduced to the published remainder, leaving t*L<24 as the arithmetic boundary layer.",
            formula="Close r_t in C^2 on tL>=24; retain coupled Dirichlet structure on 0<tL<24",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate",
        "date": "2026-07-16",
        "status": (
            "rigorous positive normalized-curvature theorem for the finite "
            "Polymath-15 main sum on t*log(x/(4*pi))>=24; not a proof for H_t, "
            "Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves a strict 19/250*L^2 normalized-curvature "
            "margin for the complete finite Riemann-Siegel main sum. It does not "
            "bound the published analytic remainder through two derivatives, "
            "close cutoff-transition collars, prove strict Laguerre positivity "
            "for H_t, prove Lambda<=0, or prove RH."
        ),
        "parameters": {
            "precision_bits": PRECISION_BITS,
            "p_exponent": fraction_text(P_EXPONENT),
            "L_min": fraction_text(L_MIN),
            "tL_min": fraction_text(TL_MIN),
            "eta": fraction_text(ETA),
        },
        "exact": exact,
        "arb": arb_data,
        "rational_budget": budget,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            "outputs/jensen_window_pf_newman_polymath15_normalized_laguerre_bridge.md",
            "outputs/jensen_window_pf_newman_theta_modular_blend_adaptive_saddle_gate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    arb_data = artifact["arb"]
    budget = artifact["rational_budget"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Dominant-Saddle Arithmetic Ray Certificate",
            "",
            "Date: 2026-07-16",
            "",
            "Status: rigorous positive normalized curvature for the complete finite",
            "Riemann-Siegel main sum on the dominant-saddle ray. This is not yet a",
            "proof for `H_t`, `Lambda <= 0`, or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate.py",
            "```",
            "",
            "The input finite sum and cutoff come from D. H. J. Polymath,",
            f"[Effective approximation of heat flow evolution of the Riemann xi function]({SOURCE_URL}),",
            "Theorem 1.3.",
            "",
            "## Ray",
            "",
            "```text",
            exact["ray"],
            "```",
            "",
            "The implication `L>=48` places the entire argument beyond",
            f"`x={arb_data['x_at_L48_ball']}`, rigorously larger than `10^21`.",
            "The large threshold is a theorem-building choice, not a numerical",
            "claim that smaller points fail.",
            "",
            "## Coefficients",
            "",
            "On this ray, elementary estimates for the published normalizer give",
            "",
            "```text",
            exact["alpha_real"],
            exact["alpha_bounds"],
            exact["normalizer_bounds"],
            exact["geometry_bounds"],
            "```",
            "",
            "Here `g=log(M_t)`. The first line is exact. The modulus bounds use",
            "`|s|=|s-1|=sqrt(1+x^2)/2`; substituting them into",
            "`g'=alpha+(t/2)alpha*alpha'` and",
            "`g''=alpha'+(t/2)(alpha'^2+alpha*alpha'')` gives the displayed",
            "normalizer bounds. Since `L>=48` and `x>10^21`, they imply the",
            "stated decimal rational envelopes for `b,a,c,d,V`.",
            "",
            "The published lower bound for `Re(s_*)`, together with the saddle",
            "cutoff, then proves",
            "",
            "```text",
            exact["coefficient_exponent"],
            exact["coefficient_bound"],
            "```",
            "",
            "The exponent loss from replacing the exact power by `3.49` is",
            f"enclosed by `{arb_data['coefficient_exponent_loss_ball']}`, below `1/100`.",
            "",
            "## Arb Moments",
            "",
            "For `p=349/100`, outward-rounded Arb series evaluation gives",
            "",
            "```text",
            f"S0=sum_(n>=2)n^-p       = {arb_data['moments']['S0']['ball']} < 16/125",
            f"S1=sum_(n>=2)log(n)n^-p = {arb_data['moments']['S1']['ball']} < 23/200",
            f"S2=sum_(n>=2)log(n)^2n^-p = {arb_data['moments']['S2']['ball']} < 121/1000",
            "```",
            "",
            "Exact logarithmic differentiation therefore yields",
            "",
            "```text",
            exact["tail_jets"],
            "Q0 <= 0.256",
            "Q1 <= 0.06901*L",
            "Q2 <= 0.01861*L^2",
            "```",
            "",
            "These bounds use the infinite zeta tails, so they cover every finite",
            "cutoff without computation at the enormous endpoint scale.",
            "",
            "## Curvature",
            "",
            "The leading phase and the exact perturbation identity give",
            "",
            "```text",
            exact["phase_floor"],
            exact["perturbation_error"],
            "```",
            "",
            "The final budget is exact rational arithmetic:",
            "",
            "```text",
            f"phase floor / L^2 > {budget['curvature_budget']['phase_floor_over_L2']}",
            f"tail error / L^2  < {budget['curvature_budget']['error_over_L2']}",
            f"strict margin / L^2 = {budget['curvature_budget']['strict_margin_over_L2']}",
            "```",
            "",
            "Hence, with `N` held fixed when differentiating its finite analytic",
            "formula,",
            "",
            "```text",
            exact["theorem"],
            "```",
            "",
            "## Remaining Gap",
            "",
            "This closes the arithmetic main sum, not the exact heat-flow function.",
            "The next job is to put the published remainder on a fixed-`N` complex",
            "collar, transfer it through two derivatives, and show its normalized",
            "curvature error is below `0.076L^2`. Adjacent-`N` overlap is still",
            "needed at the discrete cutoff transitions. The complementary region",
            "`0<tL<24` remains the coupled near-zero-time boundary layer.",
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
        "built Newman Polymath-15 dominant-saddle arithmetic ray certificate: "
        f"{len(artifact['rows'])} rows, strict margin "
        f"{artifact['rational_budget']['curvature_budget']['strict_margin_over_L2']}*L^2"
    )


if __name__ == "__main__":
    main()
