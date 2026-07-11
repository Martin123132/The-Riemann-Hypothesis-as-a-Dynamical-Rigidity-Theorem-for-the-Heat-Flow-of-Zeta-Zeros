#!/usr/bin/env python3
"""Certify all-k first-summand dominance for the lambda=-100 Newman moments."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md"
)
DEFAULT_PRECISION_BITS = 1024
DEFAULT_T = 100
DEFAULT_K = 300
DEFAULT_EPSILON_ZERO_CAP = "0.0022"
DEFAULT_SERIES_CUTOFF = 20


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    formula: str | None = None
    diagnostics: dict | None = None
    gap: str | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def arb_upper_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.upper().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_plus(), "E")


def ratio_term(n_int: int, q: flint.arb) -> flint.arb:
    n2 = flint.arb(n_int * n_int)
    return n2 * (2 * n2 * q - 3) / (2 * q - 3) * (-(n2 - 1) * q).exp()


def ratio_sum_with_tail(q: flint.arb, cutoff: int = DEFAULT_SERIES_CUTOFF) -> dict:
    finite = flint.arb(0)
    for n_int in range(2, cutoff + 1):
        finite += ratio_term(n_int, q)
    n = flint.arb(cutoff + 1)
    comparison_constant = 2 * q / (2 * q - 3)
    first_power_term = n**4 * (-((n**2 - 1) * q)).exp()
    geometric_ratio = ((n + 1) / n) ** 4 * (-(2 * n + 1) * q).exp()
    if not bool(geometric_ratio < 1 and geometric_ratio > 0):
        raise RuntimeError("Gaussian n-tail ratio did not certify in (0,1)")
    tail = comparison_constant * first_power_term / (1 - geometric_ratio)
    return {
        "finite": finite,
        "tail": tail,
        "total": finite + tail,
        "comparison_constant": comparison_constant,
        "geometric_ratio": geometric_ratio,
    }


def log_phi1(u: flint.arb) -> flint.arb:
    q = flint.arb.pi() * (4 * u).exp()
    return flint.arb.pi().log() + 5 * u + (2 * q - 3).log() - q


def log_integrand_derivative(k: flint.arb, u: flint.arb) -> flint.arb:
    q = flint.arb.pi() * (4 * u).exp()
    return 2 * k / u - 2 * DEFAULT_T * u + 5 + 8 * q / (2 * q - 3) - 4 * q


def build_diagnostics() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    k = flint.arb(DEFAULT_K)
    L = k.log()
    sqrt_k = k.sqrt()
    pi = flint.arb.pi()
    epsilon_zero = ratio_sum_with_tail(pi)
    epsilon_zero_cap = flint.arb(DEFAULT_EPSILON_ZERO_CAP)

    q_min = pi * sqrt_k
    high_geometric = (flint.arb(3) / 2) ** 4 * (-5 * q_min).exp()
    high_prefactor = 16 * (2 * q_min / (2 * q_min - 3)) / (1 - high_geometric)
    high_prefactor_margin = 17 - high_prefactor
    high_power_endpoint = flint.arb(17).log() - 3 * pi * sqrt_k + 6 * L
    high_power_derivative_margin = flint.arb(3) * pi * sqrt_k / 2 - 6

    ga = 16 * k / L - 25 * L + 5 - 4 * pi * sqrt_k
    ga_margin = ga - 100
    ga_sqrt_derivative_margin = 8 * (L - 1) / L**2 - 2 * pi / sqrt_k
    ga_k_derivative_margin = 8 * (L - 1) / L**2 - 25 / k

    offset_c = flint.arb(22) / 25
    exp_c = (flint.arb(11) / 25).exp()
    gc = 16 * k / (L + offset_c) - 25 * (L + offset_c) + 5 - 4 * pi * exp_c * sqrt_k
    gc_sqrt_derivative_margin = (
        8 * (L - flint.arb(3) / 25) / (L + offset_c) ** 2 - 2 * pi * exp_c / sqrt_k
    )
    gc_k_derivative_margin = 8 * (L - flint.arb(3) / 25) / (L + offset_c) ** 2 - 25 / k

    exp_h = (flint.arb(2) / 5).exp()
    B = (
        -flint.arb(8) * k / (5 * (L + flint.arb(4) / 5))
        + flint.arb(5) * L / 2
        + 1
        + pi * (exp_h - 1) * sqrt_k
    )
    low_power_endpoint = B + epsilon_zero_cap.log() + 6 * L
    half_negative = (
        flint.arb(4)
        * (L - flint.arb(1) / 5)
        / (5 * (L + flint.arb(4) / 5) ** 2)
    )
    low_sqrt_derivative_margin = half_negative - pi * (exp_h - 1) / (2 * sqrt_k)
    low_k_derivative_margin = half_negative - flint.arb(17) / (2 * k)

    a = L / 8
    b = a + flint.arb(1) / 10
    c = b + flint.arb(1) / 100
    D_a = log_integrand_derivative(k, a)
    D_c = log_integrand_derivative(k, c)
    log_low_exact = (
        2 * k * (a / b).log()
        - DEFAULT_T * (a**2 - b**2)
        + log_phi1(a)
        - log_phi1(b)
        - ((c - b) * D_a).log()
    )
    low_probability_exact = log_low_exact.exp()

    positive_gates = {
        "epsilon_zero_below_cap": bool(epsilon_zero["total"] < epsilon_zero_cap),
        "high_prefactor_below_17": bool(high_prefactor_margin > 0),
        "high_endpoint_below_k_minus_6": bool(high_power_endpoint < 0),
        "high_comparison_decreasing": bool(high_power_derivative_margin > 0),
        "D_a_lower_above_100": bool(ga_margin > 0),
        "D_a_lower_increasing_sqrt_side": bool(ga_sqrt_derivative_margin > 0),
        "D_a_lower_increasing_k_side": bool(ga_k_derivative_margin > 0),
        "D_c_lower_positive": bool(gc > 0),
        "D_c_lower_increasing_sqrt_side": bool(gc_sqrt_derivative_margin > 0),
        "D_c_lower_increasing_k_side": bool(gc_k_derivative_margin > 0),
        "low_endpoint_below_k_minus_6": bool(low_power_endpoint < 0),
        "low_comparison_decreasing_sqrt_side": bool(low_sqrt_derivative_margin > 0),
        "low_comparison_decreasing_k_side": bool(low_k_derivative_margin > 0),
        "exact_D_a_positive": bool(D_a > 0),
        "exact_D_c_positive": bool(D_c > 0),
    }
    if not all(positive_gates.values()):
        failed = [key for key, value in positive_gates.items() if not value]
        raise RuntimeError(f"first-summand dominance gates failed: {failed}")

    return {
        "parameters": {
            "T": DEFAULT_T,
            "K": DEFAULT_K,
            "adaptive_a": "a(k)=log(k)/8",
            "adaptive_b": "b(k)=a(k)+1/10",
            "adaptive_c": "c(k)=b(k)+1/100",
            "series_cutoff": DEFAULT_SERIES_CUTOFF,
            "precision_bits": DEFAULT_PRECISION_BITS,
        },
        "epsilon_zero": {
            "ball": arb_text(epsilon_zero["total"]),
            "upper": arb_upper_text(epsilon_zero["total"]),
            "cap": DEFAULT_EPSILON_ZERO_CAP,
            "finite_ball": arb_text(epsilon_zero["finite"]),
            "tail_upper": arb_upper_text(epsilon_zero["tail"]),
        },
        "high_region": {
            "q_min_ball": arb_text(q_min),
            "prefactor_ball": arb_text(high_prefactor),
            "prefactor_margin_lower": arb_lower_text(high_prefactor_margin),
            "endpoint_log_margin_ball": arb_text(high_power_endpoint),
            "endpoint_log_margin_upper": arb_upper_text(high_power_endpoint),
            "derivative_margin_lower": arb_lower_text(high_power_derivative_margin),
        },
        "saddle_geometry": {
            "a_at_K_ball": arb_text(a),
            "b_at_K_ball": arb_text(b),
            "c_at_K_ball": arb_text(c),
            "D_a_ball": arb_text(D_a),
            "D_c_ball": arb_text(D_c),
            "G_a_lower_ball": arb_text(ga),
            "G_a_minus_100_lower": arb_lower_text(ga_margin),
            "G_c_lower_ball": arb_text(gc),
            "G_a_sqrt_derivative_margin_lower": arb_lower_text(ga_sqrt_derivative_margin),
            "G_a_k_derivative_margin_lower": arb_lower_text(ga_k_derivative_margin),
            "G_c_sqrt_derivative_margin_lower": arb_lower_text(gc_sqrt_derivative_margin),
            "G_c_k_derivative_margin_lower": arb_lower_text(gc_k_derivative_margin),
        },
        "low_region": {
            "coarse_B_ball": arb_text(B),
            "endpoint_log_margin_ball": arb_text(low_power_endpoint),
            "endpoint_log_margin_upper": arb_upper_text(low_power_endpoint),
            "sqrt_derivative_margin_lower": arb_lower_text(low_sqrt_derivative_margin),
            "k_derivative_margin_lower": arb_lower_text(low_k_derivative_margin),
            "exact_log_probability_bound_ball": arb_text(log_low_exact),
            "exact_probability_bound_upper": arb_upper_text(low_probability_exact),
        },
        "half_line_propagation": {
            "G_a_sqrt_ratio_log_derivative": "1/2+1/(L-1)-2/L>0 for L>=log(300)",
            "G_a_k_ratio_log_derivative": "1+1/(L-1)-2/L>0 for L>=log(300)",
            "G_c_sqrt_ratio_log_derivative": "1/2+1/(L-3/25)-2/(L+22/25)>0 for L>=log(300)",
            "G_c_k_ratio_log_derivative": "1+1/(L-3/25)-2/(L+22/25)>0 for L>=log(300)",
            "low_sqrt_ratio_log_derivative": "1/2+1/(L-1/5)-2/(L+4/5)>0 for L>=log(300)",
            "low_k_ratio_log_derivative": "1+1/(L-1/5)-2/(L+4/5)>0 for L>=log(300)",
            "high_log_comparison_derivative": "d_k(log(17)-3*pi*sqrt(k)+6*log(k))=(6-(3*pi/2)*sqrt(k))/k<0",
            "low_log_comparison_derivative": (
                "Q'(k)=-(8/5)*(L-1/5)/(L+4/5)^2+17/(2*k)"
                "+pi*(exp(2/5)-1)/(2*sqrt(k))<0"
            ),
            "comparison_rule": (
                "Each positive term in -Q' and in the G_a/G_c derivative lower bounds is split "
                "into two endpoint comparisons; the displayed ratio log-derivatives prove those "
                "comparisons strengthen for k>=300."
            ),
        },
        "positive_gates": positive_gates,
        "all_positive_gates": all(positive_gates.values()),
        "full_tail_relative_bound": "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<=2/k^6 for every integer k>=300",
        "adjacent_log_wall_error_bound": (
            "for k>=301, |L_k-L_k^(1)|<=16/(k-1)^6, "
            "where L_k=log(x_(k+1)/x_k)"
        ),
    }


def build_artifact() -> dict:
    diagnostics = build_diagnostics()
    rows = [
        CertificateRow(
            id="fsdc_01_exact_original_ratio",
            role="exact_reduction",
            readiness="available_exact",
            claim="The nth-to-first Newman summand ratio has an explicit positive q-form.",
            formula="r_n(q)=n^2*(2*n^2*q-3)/(2*q-3)*exp(-(n^2-1)*q), q=pi*exp(4u)",
            proof_boundary="Exact pointwise kernel ratio only.",
        ),
        CertificateRow(
            id="fsdc_02_ratio_monotonicity",
            role="exact_lemma",
            readiness="available_exact",
            claim="For every n>=2, r_n is strictly decreasing in q and hence in u>=0.",
            formula="d_q log(r_n)=-(n^2-1)*(1+6/((2*n^2*q-3)*(2*q-3)))<0",
            proof_boundary="Exact all-n pointwise monotonicity; no moment comparison yet.",
        ),
        CertificateRow(
            id="fsdc_03_zero_endpoint_tail_sum",
            role="interval_analytic_bound",
            readiness="interval_validated",
            claim="A finite n=2..20 Arb sum plus a geometric n-tail proves epsilon(0)<0.0022.",
            formula="epsilon(u)=sum_(n>=2) r_n(pi*exp(4u))",
            diagnostics=diagnostics["epsilon_zero"],
            proof_boundary="Full pointwise n-tail bound at u=0 only.",
        ),
        CertificateRow(
            id="fsdc_04_first_integrand_strict_concavity",
            role="exact_lemma",
            readiness="available_exact",
            claim="The logarithm S_k(u) of the n=1 moment integrand is strictly concave for every k>0 and T>0.",
            formula=(
                "S_k''(u)=-2*k/u^2-2*T-16*q-96*q/(2*q-3)^2<0, "
                "q=pi*exp(4u)"
            ),
            proof_boundary="Exact dominant-integrand geometry only; not the adjacent moment inequality.",
        ),
        CertificateRow(
            id="fsdc_05_adaptive_saddle_bracket",
            role="exact_analytic_bound",
            readiness="interval_validated",
            claim="For every k>=300, S_k is increasing through c(k)=log(k)/8+11/100, so the low-region integral admits a one-sided concavity bound.",
            formula=(
                "S_k'(u)=2*k/u-200*u+5+8*q/(2*q-3)-4*q; "
                "S_k'(a(k))>100 and S_k'(c(k))>0; the endpoint derivative "
                "splits propagate because sqrt(k)*(L-alpha)/(L+beta)^2 and "
                "k*(L-alpha)/(L+beta)^2 increase"
            ),
            diagnostics=diagnostics["saddle_geometry"],
            proof_boundary="Analytic half-line endpoint propagation for the chosen adaptive bracket only.",
        ),
        CertificateRow(
            id="fsdc_06_low_region_probability_bound",
            role="exact_analytic_bound",
            readiness="interval_validated",
            claim="The n=1 tilted probability of 0<u<a(k) times epsilon(0) is below k^(-6) for every k>=300.",
            formula=(
                "P_k(u<a)<=exp(S_k(a)-S_k(b))/(0.01*S_k'(a)); "
                "epsilon(0)*P_k(u<a)<k^(-6); Q'(k)<0 by two increasing-ratio comparisons"
            ),
            diagnostics=diagnostics["low_region"],
            proof_boundary="Uniform low-region first-summand mass comparison only.",
        ),
        CertificateRow(
            id="fsdc_07_high_region_kernel_tail_bound",
            role="exact_analytic_bound",
            readiness="interval_validated",
            claim="On u>=a(k), the complete pointwise n>=2 kernel ratio is below k^(-6) for every k>=300.",
            formula="epsilon(a(k))<=17*exp(-3*pi*sqrt(k))<k^(-6)",
            diagnostics=diagnostics["high_region"],
            proof_boundary="Uniform high-region pointwise kernel bound only.",
        ),
        CertificateRow(
            id="fsdc_08_full_moment_dominance",
            role="interval_theorem",
            readiness="interval_validated",
            claim="For lambda=-100 and every integer k>=300, the complete n>=2 moment tail is at most 2/k^6 of the n=1 moment.",
            formula="M_k=M_k^(1)*(1+delta_k), 0<=delta_k<=2/k^6",
            diagnostics={"full_tail_relative_bound": diagnostics["full_tail_relative_bound"]},
            proof_boundary="All-k coefficient dominance theorem only; not adjacent-k monotonicity.",
        ),
        CertificateRow(
            id="fsdc_09_adjacent_wall_stability",
            role="exact_reduction",
            readiness="available_exact",
            claim="For k>=301, the full-kernel adjacent log wall differs from the n=1 wall by at most 16/(k-1)^6.",
            formula=(
                "L_k-L_k^(1)=e_(k+2)-3*e_(k+1)+3*e_k-e_(k-1), "
                "e_j=log(1+delta_j)"
            ),
            diagnostics={"error_bound": diagnostics["adjacent_log_wall_error_bound"]},
            proof_boundary="Exact perturbation transfer conditional on a positive n=1 wall margin.",
        ),
        CertificateRow(
            id="fsdc_10_dominant_wall_handoff",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim="It now suffices to prove L_k^(1)>=1/(4*k^2) for every k>=301; this would dominate the certified full-kernel error and close the lambda=-100 adjacent wall after the finite collar.",
            formula="L_k^(1)=log(x_(k+1)^(1)/x_k^(1))",
            gap="A referee-grade dominant n=1 saddle inequality remains open.",
            proof_boundary="Open dominant-summand theorem only; not cone entry, RH, or Lambda <= 0.",
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "available_exact_rows": sum(row.readiness == "available_exact" for row in rows),
        "interval_validated_rows": sum(row.readiness == "interval_validated" for row in rows),
        "open_requirement_rows": sum(row.readiness == "not_ready_to_apply" for row in rows),
        "positive_analytic_gates": sum(diagnostics["positive_gates"].values()),
        "full_tail_power": 6,
        "full_tail_constant": 2,
        "tail_start_k": DEFAULT_K,
        "wall_transfer_start_k": DEFAULT_K + 1,
        "ready_to_apply_rows": 0,
        "main_finding": (
            "Exact monotonicity of every original-variable tail ratio, strict concavity of the "
            "first-summand tilted integrand, and 15 positive Arb endpoint gates prove that at "
            "lambda=-100 the complete n>=2 moment tail satisfies 0<=delta_k<=2/k^6 for every "
            "integer k>=300. Consequently the full adjacent log wall differs from the n=1 wall "
            "by at most 16/(k-1)^6 for k>=301. The shifted far-tail source is discharged; the "
            "only remaining tail inequality is a quantitative dominant n=1 saddle wall."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_dominance_certificate",
        "status": "all-k analytic first-summand dominance certificate",
        "date": "2026-07-10",
        "proof_boundary": (
            "This analytic/interval certificate proves all-k n=1 dominance for the lambda=-100 "
            "moments and an exact adjacent-wall perturbation bound. It does not prove the "
            "dominant n=1 adjacent-wall lower bound, does not by itself prove cone entry, and "
            "does not prove RH or Lambda <= 0."
        ),
        "source_shift_lemma": "outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md",
        "source_k300_audit": "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py",
        "diagnostics": diagnostics,
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda First-Summand Dominance Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: all-k analytic first-summand dominance certificate. This is not",
        "a proof of cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_dominance_certificate`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_dominance_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF negative-lambda first-summand dominance certificate: 10 rows, 0 issues, 4 exact rows, 5 interval rows, 15 positive analytic gates, 1 open dominant-wall row, 0 ready-to-apply rows",
        "```",
        "",
        "## Exact Kernel Ratio",
        "",
        "For `q=pi*exp(4u)`, divide the nth positive Newman summand by the",
        "first:",
        "",
        "```text",
        "r_n(q)=n^2*(2*n^2*q-3)/(2*q-3)*exp(-(n^2-1)*q).",
        "d_q log(r_n)=-(n^2-1)*(1+6/((2*n^2*q-3)*(2*q-3)))<0.",
        "```",
        "",
        "Thus every tail ratio and their sum `epsilon(u)` decrease with `u`.",
        "A finite Arb sum through `n=20` plus a geometric tail gives",
        "",
        "```text",
        f"epsilon(0) = {diagnostics['epsilon_zero']['ball']}",
        f"epsilon(0) < {diagnostics['epsilon_zero']['cap']}.",
        "```",
        "",
        "## Adaptive Split",
        "",
        "For every integer `k>=300`, set",
        "",
        "```text",
        "a(k)=log(k)/8",
        "b(k)=a(k)+1/10",
        "c(k)=b(k)+1/100.",
        "```",
        "",
        "If `S_k` is the logarithm of the `n=1`, `lambda=-100` moment",
        "integrand, then",
        "",
        "```text",
        "S_k'(u)=2*k/u-200*u+5+8*q/(2*q-3)-4*q",
        "S_k''(u)=-2*k/u^2-200-16*q-96*q/(2*q-3)^2<0.",
        "```",
        "",
        "Endpoint bounds and monotone comparison functions prove",
        "`S_k'(a(k))>100` and `S_k'(c(k))>0` for the complete half-line",
        "`k>=300`. Concavity therefore gives",
        "",
        "```text",
        "P_k(u<a(k)) <= exp(S_k(a)-S_k(b))/(0.01*S_k'(a)).",
        "epsilon(0)*P_k(u<a(k)) < k^(-6).",
        "```",
        "",
        "On the complementary region, the decreasing kernel ratio and",
        "`q(a(k))=pi*sqrt(k)` give",
        "",
        "```text",
        "epsilon(a(k)) <= 17*exp(-3*pi*sqrt(k)) < k^(-6).",
        "```",
        "",
        "All 15 endpoint and monotonic-propagation gates are Arb-positive.",
        "The half-line propagation is exact: with `L=log(k)`, the comparison",
        "ratios have logarithmic derivatives",
        "",
        "```text",
        "1/2+1/(L-alpha)-2/(L+beta)>0",
        "1  +1/(L-alpha)-2/(L+beta)>0",
        "```",
        "",
        "for the three pairs `(alpha,beta)=(1,0),(3/25,22/25),`",
        "`(1/5,4/5)` and every `L>=log(300)`. Hence each split endpoint",
        "comparison only strengthens with `k`. Also",
        "",
        "```text",
        "d_k[log(17)-3*pi*sqrt(k)+6*log(k)]",
        "  =(6-(3*pi/2)*sqrt(k))/k<0.",
        "```",
        "",
        "## All-k Consequence",
        "",
        "Let `M_k^(1)` be the first-summand moment and `M_k` the full kernel",
        "moment. The two regions compose to prove",
        "",
        "```text",
        "M_k=M_k^(1)*(1+delta_k), 0<=delta_k<=2/k^6, k>=300.",
        "```",
        "",
        "For the adjacent log wall",
        "",
        "```text",
        "L_k=log(x_(k+1)/x_k),",
        "|L_k-L_k^(1)|<=16/(k-1)^6, k>=301.",
        "```",
        "",
        "The shifted `v>=3/2` far-tail source from the earlier shift lemma is",
        "therefore discharged. It remains to prove a quantitative first-summand",
        "wall, for example",
        "",
        "```text",
        "L_k^(1)>=1/(4*k^2), k>=301.",
        "```",
        "",
        "That bound would dominate the certified perturbation by many orders of",
        "magnitude and splice to the repaired finite collar.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_kernel_summand_shift_lemma.md",
        "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF negative-lambda first-summand dominance certificate: "
        "10 rows, 0 issues, 4 exact rows, 5 interval rows, 15 positive analytic gates, "
        "1 open dominant-wall row, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
