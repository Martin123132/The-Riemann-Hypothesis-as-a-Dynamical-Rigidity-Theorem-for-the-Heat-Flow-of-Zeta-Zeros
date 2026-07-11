#!/usr/bin/env python3
"""Certify kernel sqrt-log-concavity and the Mellin upper cone wall.

The interval part proves that y -> Phi(sqrt(y)) is strictly log-concave on
[0, infinity).  The final coefficient inequality uses the classical
Berwald-Borell theorem for Gamma-normalized Mellin transforms of log-concave
measures.  It does not prove the adjacent-k monotonicity wall.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import REPO_ROOT  # noqa: E402
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)
from jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate import (  # noqa: E402
    full_disk_majorant,
)


DEFAULT_EVENNESS_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.json"
)
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_kernel_mellin_upper_wall_certificate.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md"

DEFAULT_PRECISION_BITS = 1024
DEFAULT_COEFFICIENT_MAX_DEGREE = 40
DEFAULT_COEFFICIENT_CUTOFF_N = 80
DEFAULT_COMPACT_Y_MAX = "0.04"
DEFAULT_COMPACT_SUBINTERVALS = 200
DEFAULT_RAY_U_START = "0.2"
DEFAULT_RAY_TAIL_CUTOFF_N = 40
DEFAULT_C0_LOWER = "0.44"
DEFAULT_FIRST_OMITTED_Y_DEGREE = 21


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    source_artifacts: list[str]
    diagnostics: dict | None = None
    gap: str | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


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


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def polynomial_value(coefficients: list[flint.arb], x: flint.arb) -> flint.arb:
    total = flint.arb(0)
    for coefficient in reversed(coefficients):
        total = total * x + coefficient
    return total


def cauchy_y_tail_bounds(disk: dict[str, flint.arb | bool]) -> dict[str, flint.arb]:
    radius = disk["radius"]
    majorant = disk["majorant"]
    assert isinstance(radius, flint.arb)
    assert isinstance(majorant, flint.arb)
    normalized_majorant = majorant / flint.arb(DEFAULT_C0_LOWER)
    y_max = flint.arb(DEFAULT_COMPACT_Y_MAX)
    z = y_max / (radius * radius)
    m = DEFAULT_FIRST_OMITTED_Y_DEGREE
    if not arb_positive(1 - z):
        raise RuntimeError("compact y-range did not fit inside the Cauchy disk")
    value = normalized_majorant * z**m / (1 - z)
    first = normalized_majorant / (radius * radius) * (
        flint.arb(m) * z ** (m - 1) / (1 - z) + z**m / (1 - z) ** 2
    )
    second = normalized_majorant / radius**4 * (
        flint.arb(m * (m - 1)) * z ** (m - 2) / (1 - z)
        + flint.arb(2 * m) * z ** (m - 1) / (1 - z) ** 2
        + flint.arb(2) * z**m / (1 - z) ** 3
    )
    return {
        "normalized_majorant": normalized_majorant,
        "disk_ratio": z,
        "value": value,
        "first": first,
        "second": second,
    }


def compact_log_concavity_certificate() -> dict:
    ratio_rows, ratios = build_ratio_rows(
        DEFAULT_COEFFICIENT_MAX_DEGREE,
        DEFAULT_COEFFICIENT_CUTOFF_N,
    )
    disk = full_disk_majorant()
    tails = cauchy_y_tail_bounds(disk)
    first_coefficients = [flint.arb(j) * ratios[j] for j in range(1, len(ratios))]
    second_coefficients = [
        flint.arb(j * (j - 1)) * ratios[j] for j in range(2, len(ratios))
    ]
    y_max = flint.arb(DEFAULT_COMPACT_Y_MAX)
    count = DEFAULT_COMPACT_SUBINTERVALS
    interval_rows = []
    minimum_q: tuple[flint.arb, int] | None = None
    minimum_g: tuple[flint.arb, int] | None = None
    positive_rows = 0
    positive_g_rows = 0
    for index in range(count):
        left = flint.arb(index) * y_max / count
        right = flint.arb(index + 1) * y_max / count
        y = (left + right) / 2 + flint.arb(0, (right - left) / 2)
        g = polynomial_value(ratios, y) + flint.arb(0, tails["value"])
        g_first = polynomial_value(first_coefficients, y) + flint.arb(0, tails["first"])
        g_second = polynomial_value(second_coefficients, y) + flint.arb(0, tails["second"])
        q_value = g_first * g_first - g * g_second
        q_positive = arb_positive(q_value)
        g_positive = arb_positive(g)
        positive_rows += int(q_positive)
        positive_g_rows += int(g_positive)
        if minimum_q is None or q_value.lower() < minimum_q[0].lower():
            minimum_q = (q_value, index)
        if minimum_g is None or g.lower() < minimum_g[0].lower():
            minimum_g = (g, index)
        if index in {0, count // 4, count // 2, 3 * count // 4, count - 1}:
            interval_rows.append(
                {
                    "index": index,
                    "y_left": arb_text(left, 30),
                    "y_right": arb_text(right, 30),
                    "g_ball": arb_text(g, 40),
                    "Q_ball": arb_text(q_value, 40),
                    "g_positive": g_positive,
                    "Q_positive": q_positive,
                }
            )
    if minimum_q is None or minimum_g is None:
        raise RuntimeError("compact subdivision produced no rows")
    if positive_rows != count or positive_g_rows != count:
        raise RuntimeError("compact log-concavity subdivision failed")
    return {
        "coefficient_ball_rows": [asdict(row) for row in ratio_rows],
        "coefficient_ball_count": len(ratio_rows),
        "disk_radius": arb_text(disk["radius"], 40),
        "disk_majorant": arb_text(disk["majorant"], 50),
        "cauchy_disk_ratio_upper": arb_upper_text(tails["disk_ratio"], 60),
        "cauchy_value_tail_radius_upper": arb_upper_text(tails["value"], 60),
        "cauchy_first_tail_radius_upper": arb_upper_text(tails["first"], 60),
        "cauchy_second_tail_radius_upper": arb_upper_text(tails["second"], 60),
        "subinterval_count": count,
        "positive_g_subintervals": positive_g_rows,
        "positive_Q_subintervals": positive_rows,
        "minimum_g_lower": arb_lower_text(minimum_g[0], 60),
        "minimum_g_subinterval": minimum_g[1],
        "minimum_Q_lower": arb_lower_text(minimum_q[0], 60),
        "minimum_Q_subinterval": minimum_q[1],
        "selected_subinterval_rows": interval_rows,
        "compact_log_concavity_certified": positive_rows == count,
    }


def gaussian_n_power_tail(power: int, q: flint.arb) -> tuple[flint.arb, flint.arb]:
    cutoff = DEFAULT_RAY_TAIL_CUTOFF_N
    total = flint.arb(0)
    for n_int in range(2, cutoff + 1):
        n = flint.arb(n_int)
        total += n**power * (-(flint.arb(n_int * n_int - 1) * q)).exp()
    n = flint.arb(cutoff + 1)
    first = n**power * (-((n * n - 1) * q)).exp()
    ratio = ((n + 1) / n) ** power * (-(q * (2 * n + 1))).exp()
    if not arb_positive(1 - ratio):
        raise RuntimeError("ray n-tail ratio did not certify below one")
    return total + first / (1 - ratio), ratio


def ray_log_concavity_certificate() -> dict:
    u0 = flint.arb(DEFAULT_RAY_U_START)
    q0 = flint.arb.pi() * (4 * u0).exp()
    denominator = 2 * q0 - 3
    h0 = 8 * q0 / denominator
    g0 = 96 * q0 / denominator**2
    ratio_constant = 2 * q0 / denominator

    log_first = 5 - 4 * q0 + h0
    log_second = -16 * q0 - g0
    n1_D = log_first - u0 * log_second
    third_derivative_margin = denominator**3 - 6 * (2 * q0 + 3)
    if not arb_positive(third_derivative_margin):
        raise RuntimeError("n=1 ray monotonicity margin failed")

    sum4, ratio4 = gaussian_n_power_tail(4, q0)
    sum6, ratio6 = gaussian_n_power_tail(6, q0)
    sum8, ratio8 = gaussian_n_power_tail(8, q0)
    s1 = ratio_constant * (4 * q0 * sum6 + h0 * sum4)
    s2 = ratio_constant * (
        16 * q0 * q0 * sum8
        + (8 * q0 * h0 + 16 * q0) * sum6
        + (h0 * h0 + g0) * sum4
    )
    perturbation = s1 + u0 * (s2 + s1 * s1)
    full_D = n1_D - perturbation
    endpoint_decay_margin = 4 * flint.arb(3) * q0 - flint.arb(13)
    if not all(
        arb_positive(value)
        for value in (n1_D, endpoint_decay_margin, full_D)
    ):
        raise RuntimeError("analytic ray log-concavity margin failed")
    return {
        "u_start": arb_text(u0, 30),
        "y_start": arb_text(u0 * u0, 30),
        "q_start": arb_text(q0, 60),
        "n1_D_lower": arb_lower_text(n1_D, 60),
        "n1_third_derivative_margin_lower": arb_lower_text(
            third_derivative_margin, 60
        ),
        "endpoint_decay_margin_lower": arb_lower_text(endpoint_decay_margin, 60),
        "S1_upper": arb_upper_text(s1, 60),
        "S2_upper": arb_upper_text(s2, 60),
        "perturbation_upper": arb_upper_text(perturbation, 60),
        "full_D_lower": arb_lower_text(full_D, 60),
        "tail_ratio_power4_upper": arb_upper_text(ratio4, 40),
        "tail_ratio_power6_upper": arb_upper_text(ratio6, 40),
        "tail_ratio_power8_upper": arb_upper_text(ratio8, 40),
        "ray_log_concavity_certified": arb_positive(full_D),
        "proof_formula": (
            "D(u)=L'(u)-u*L''(u)>0 implies d_y^2 log(Phi(sqrt(y)))=-D(u)/(4*u^3)<0"
        ),
    }


def build_diagnostics(evenness: dict) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    if evenness["summary"].get("full_kernel_evenness_certified") is not True:
        raise RuntimeError("source full-kernel evenness was not certified")
    if evenness["summary"].get("residual_order_42_certified") is not True:
        raise RuntimeError("source order-42 residual zero was not certified")
    compact = compact_log_concavity_certificate()
    ray = ray_log_concavity_certificate()
    global_certificate = bool(
        compact["compact_log_concavity_certified"]
        and ray["ray_log_concavity_certified"]
    )
    return {
        "parameters": {
            "precision_bits": DEFAULT_PRECISION_BITS,
            "coefficient_max_degree": DEFAULT_COEFFICIENT_MAX_DEGREE,
            "coefficient_cutoff_n": DEFAULT_COEFFICIENT_CUTOFF_N,
            "compact_y_interval": "0<=y<=1/25",
            "compact_subintervals": DEFAULT_COMPACT_SUBINTERVALS,
            "cauchy_disk_radius": "0.38",
            "ray_u_interval": "u>=1/5",
            "ray_y_interval": "y>=1/25",
            "ray_tail_cutoff_n": DEFAULT_RAY_TAIL_CUTOFF_N,
        },
        "source_full_kernel_evenness_certified": True,
        "source_order42_residual_zero_certified": True,
        "compact": compact,
        "ray": ray,
        "global_kernel_sqrt_log_concavity_certified": global_certificate,
        "mellin_reindexing": (
            "M_k(lambda)=int_0^infty y^(k-1/2)*exp(lambda*y)*Phi(sqrt(y))dy; "
            "A_k(lambda)=sqrt(pi)*4^(-k)*M_k(lambda)/Gamma(k+1/2)"
        ),
        "berwald_borell_input": (
            "For a log-concave measure on [0,infinity), its Mellin transform divided by Gamma is log-concave."
        ),
        "all_real_lambda_integrability": (
            "The exp(-pi*exp(4*sqrt(y))) kernel decay dominates exp(lambda*y) for every fixed real lambda."
        ),
        "all_k_upper_wall_certified": global_certificate,
        "remaining_open_cone_clause": "x_(k+1)(lambda)>=x_k(lambda)",
    }


def build_rows(diagnostics: dict) -> list[CertificateRow]:
    note = "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md"
    evenness_note = (
        "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md"
    )
    cone_target = "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md"
    defect_target = "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md"
    return [
        CertificateRow(
            id="kmuw_01_exact_y_mellin_reindexing",
            role="exact_reduction",
            readiness="proved",
            claim=(
                "The substitution y=u^2 and Legendre duplication give "
                "A_k=sqrt(pi)*4^(-k)*M_k/Gamma(k+1/2)."
            ),
            proof_boundary="Exact coefficient reindexing only.",
            source_artifacts=[cone_target],
        ),
        CertificateRow(
            id="kmuw_02_full_kernel_evenness_input",
            role="exact_input",
            readiness="available_exact",
            claim=(
                "Exact evenness makes Phi(sqrt(y)) analytic at y=0 and gives an even Taylor series."
            ),
            proof_boundary="Imports the exact full-kernel lemma; no finite truncation parity is assumed.",
            source_artifacts=[evenness_note],
        ),
        CertificateRow(
            id="kmuw_03_compact_interval_log_concavity",
            role="interval_certificate",
            readiness="interval_validated",
            claim=(
                "On 0<=y<=1/25, Arb Taylor balls plus a full-kernel Cauchy tail prove "
                "g'(y)^2-g(y)g''(y)>0 on every subdivision."
            ),
            proof_boundary="Compact full-kernel interval theorem only.",
            source_artifacts=[note, evenness_note],
            diagnostics={
                "subintervals": diagnostics["compact"]["subinterval_count"],
                "minimum_Q_lower": diagnostics["compact"]["minimum_Q_lower"],
            },
        ),
        CertificateRow(
            id="kmuw_04_n1_analytic_ray",
            role="exact_analytic_bound",
            readiness="proved",
            claim=(
                "For u>=1/5 the n=1 summand has D_1=L_1'-uL_1'' increasing and positive."
            ),
            proof_boundary="Dominant-summand ray bound; the full tail is handled separately.",
            source_artifacts=[note],
            diagnostics={"n1_D_lower": diagnostics["ray"]["n1_D_lower"]},
        ),
        CertificateRow(
            id="kmuw_05_full_n_tail_ray",
            role="interval_analytic_bound",
            readiness="interval_validated",
            claim=(
                "Explicit geometric n-tail sums bound the log-derivative perturbation below the n=1 ray margin."
            ),
            proof_boundary="Full infinite-kernel ray composition only.",
            source_artifacts=[note],
            diagnostics={
                "perturbation_upper": diagnostics["ray"]["perturbation_upper"],
                "full_D_lower": diagnostics["ray"]["full_D_lower"],
            },
        ),
        CertificateRow(
            id="kmuw_06_global_kernel_sqrt_log_concavity",
            role="interval_theorem",
            readiness="interval_validated",
            claim="The function y -> Phi(sqrt(y)) is strictly log-concave on [0,infinity).",
            proof_boundary="Kernel theorem only; it is not adjacent-k monotonicity.",
            source_artifacts=[note, evenness_note],
        ),
        CertificateRow(
            id="kmuw_07_berwald_borell_upper_wall",
            role="literature_theorem_application",
            readiness="proved",
            claim=(
                "Berwald-Borell applied to exp(lambda*y)*Phi(sqrt(y)) proves x_k(lambda)<=1 "
                "for every real lambda and every k>=1."
            ),
            proof_boundary=(
                "This proves the upper cone wall only; the lower wall is the separate Cauchy-Schwarz lemma."
            ),
            source_artifacts=[note, "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md"],
        ),
        CertificateRow(
            id="kmuw_08_remaining_monotone_wall",
            role="non_promotion_gate",
            readiness="open",
            claim=(
                "Full cone entry still requires x_(k+1)>=x_k; neither log-concavity nor the Mellin theorem supplies it."
            ),
            proof_boundary="No claim of cone entry, Jensen-window PF-infinity, RH, or Lambda<=0.",
            source_artifacts=[cone_target, defect_target],
            gap="Prove the zeta-specific adjacent-k monotone-contraction inequality for all k.",
        ),
    ]


def build_artifact(evenness_path: Path) -> dict:
    evenness = load_json(evenness_path)
    diagnostics = build_diagnostics(evenness)
    rows = build_rows(diagnostics)
    summary = {
        "certificate_rows": len(rows),
        "ready_to_apply_rows": 0,
        "compact_subintervals": diagnostics["compact"]["subinterval_count"],
        "positive_compact_subintervals": diagnostics["compact"]["positive_Q_subintervals"],
        "compact_minimum_Q_lower": diagnostics["compact"]["minimum_Q_lower"],
        "ray_full_D_lower": diagnostics["ray"]["full_D_lower"],
        "global_kernel_sqrt_log_concavity_certified": diagnostics[
            "global_kernel_sqrt_log_concavity_certified"
        ],
        "all_k_upper_wall_certified": diagnostics["all_k_upper_wall_certified"],
        "remaining_open_cone_clauses": 1,
        "main_finding": (
            "A compact Arb/Cauchy certificate and an analytic n=1-dominant ray prove that "
            "y->Phi(sqrt(y)) is strictly log-concave. Berwald-Borell then proves the upper "
            "ratio-cone wall x_k(lambda)<=1 for every real lambda and k>=1. The adjacent-k "
            "monotonicity x_(k+1)>=x_k remains open, so this is not cone entry or Lambda<=0."
        ),
    }
    return {
        "kind": "jensen_window_pf_kernel_mellin_upper_wall_certificate",
        "status": "interval theorem certificate",
        "date": "2026-07-10",
        "proof_boundary": (
            "This interval-backed theorem certificate proves global log-concavity of the full "
            "kernel after y=u^2 and, via the classical Berwald-Borell Mellin theorem, the "
            "all-k upper cone wall. It does not prove the adjacent-k monotone wall, full cone "
            "entry, Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
        "source_evenness_json": evenness_path.relative_to(REPO_ROOT).as_posix(),
        "source_evenness_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md"
        ),
        "source_boundary_threshold_note": "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "literature_sources": [
            "https://doi.org/10.1007/BF01362702",
            "https://www.math.tau.ac.il/~klartagb/papers/log_concave_bernstein.pdf",
        ],
        "generator": "work/rh_compute/scripts/jensen_window_pf_kernel_mellin_upper_wall_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py",
        "summary": summary,
        "diagnostics": diagnostics,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    compact = artifact["diagnostics"]["compact"]
    ray = artifact["diagnostics"]["ray"]
    lines = [
        "# Jensen-Window PF Kernel Mellin Upper-Wall Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: interval theorem certificate. This is not a proof of full cone entry,",
        "Jensen-window PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_kernel_mellin_upper_wall_certificate`.",
        "",
        "Proof boundary: this certificate proves the all-`k` upper wall `x_k<=1`.",
        "It does not prove the remaining adjacent-`k` wall `x_(k+1)>=x_k`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_kernel_mellin_upper_wall_certificate.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_kernel_mellin_upper_wall_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_kernel_mellin_upper_wall_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF kernel Mellin upper-wall certificate: "
            f"{summary['certificate_rows']} rows, 0 issues, "
            f"{summary['positive_compact_subintervals']} positive compact intervals, "
            "1 positive analytic ray, 1 remaining open cone clause, 0 ready-to-apply rows"
        ),
        "```",
        "",
        "## Exact Mellin Reduction",
        "",
        "With `y=u^2` and `g(y)=Phi(sqrt(y))`,",
        "",
        "```text",
        "M_k(lambda) = integral_0^infinity y^(k-1/2)*exp(lambda*y)*g(y) dy",
        "A_k(lambda) = sqrt(pi)*4^(-k)*M_k(lambda)/Gamma(k+1/2)",
        "x_k = (A_(k+1)/A_k)/(A_k/A_(k-1))",
        "```",
        "",
        "The geometric factor `4^(-k)` cancels from `x_k`. The Berwald-Borell",
        "inequality says that the Mellin transform of a log-concave measure on",
        "`[0,infinity)` divided by `Gamma` is log-concave. The primary references",
        "used for that theorem are:",
        "",
        "- Christer Borell, *Complements of Lyapunov's inequality*: https://doi.org/10.1007/BF01362702",
        "- Bo'az Klartag and Joseph Lehec, *Poisson processes and a log-concave Bernstein theorem*: https://www.math.tau.ac.il/~klartagb/papers/log_concave_bernstein.pdf",
        "",
        "Thus it remains to prove that `g(y)` is log-concave. Multiplication by",
        "`exp(lambda*y)` preserves log-concavity for every real `lambda`.",
        "",
        "## Compact Certificate",
        "",
        "Exact full-kernel evenness supplies an analytic series",
        "`g(y)/g(0)=sum r_j*y^j`. On `0<=y<=1/25`, coefficient balls through",
        "`j=20` are combined with a degree-21 Cauchy tail and its first two",
        "`y` derivatives. Arb proves",
        "",
        "```text",
        "Q(y)=g'(y)^2-g(y)*g''(y) > 0",
        "```",
        "",
        f"on `{compact['positive_Q_subintervals']}/{compact['subinterval_count']}` rational subintervals.",
        "",
        f"- Minimum outward-rounded `Q` lower bound: `{compact['minimum_Q_lower']}`.",
        f"- Minimum outward-rounded `g` lower bound: `{compact['minimum_g_lower']}`.",
        f"- Cauchy value-tail radius upper: `{compact['cauchy_value_tail_radius_upper']}`.",
        f"- Cauchy first-derivative radius upper: `{compact['cauchy_first_tail_radius_upper']}`.",
        f"- Cauchy second-derivative radius upper: `{compact['cauchy_second_tail_radius_upper']}`.",
        "",
        "## Analytic Ray",
        "",
        "For `u>=1/5`, write `Phi(u)=t_1(u)*(1+R(u))` and",
        "`L=log(Phi)`. Log-concavity in `y=u^2` is equivalent to",
        "",
        "```text",
        "D(u)=L'(u)-u*L''(u) > 0.",
        "```",
        "",
        "The `n=1` contribution is increasing on the ray. Explicit geometric",
        "bounds for the `n>=2` log-derivative perturbation give:",
        "",
        f"- `n=1` margin lower at `u=1/5`: `{ray['n1_D_lower']}`.",
        f"- Full perturbation upper: `{ray['perturbation_upper']}`.",
        f"- Full ray `D` lower: `{ray['full_D_lower']}`.",
        "",
        "The compact interval and ray meet at `y=1/25`, so",
        "`y -> Phi(sqrt(y))` is strictly log-concave on the full half-line.",
        "The kernel's double-exponential tail makes all fixed-real-`lambda`",
        "Mellin moments finite.",
        "",
        "## Cone Consequence",
        "",
        "Berwald-Borell now gives, for every real `lambda` and `k>=1`,",
        "",
        "```text",
        "A_k(lambda)^2 >= A_(k-1)(lambda)*A_(k+1)(lambda)",
        "x_k(lambda) <= 1.",
        "```",
        "",
        "The separate Cauchy-Schwarz moment lemma already gives the lower wall",
        "`(2*k-1)/(2*k+1)<=x_k`. The only remaining ratio-cone clause is",
        "`x_(k+1)>=x_k`. That clause is still open and is not implied by ordinary",
        "log-concavity or by Berwald-Borell.",
        "",
        "## Integration",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md",
        "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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
    parser.add_argument("--evenness-json", type=Path, default=DEFAULT_EVENNESS_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact(args.evenness_json)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "validated Jensen-window PF kernel Mellin upper-wall certificate: "
        f"{summary['certificate_rows']} rows, 0 issues, "
        f"{summary['positive_compact_subintervals']} positive compact intervals, "
        "1 positive analytic ray, 1 remaining open cone clause, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
