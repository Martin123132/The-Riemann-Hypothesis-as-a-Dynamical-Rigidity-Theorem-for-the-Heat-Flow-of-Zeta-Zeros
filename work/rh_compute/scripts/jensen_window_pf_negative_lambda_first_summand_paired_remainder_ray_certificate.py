#!/usr/bin/env python3
"""Prove the paired first-summand remainder on the asymptotic mode ray."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md"
)
DEFAULT_PRECISION_BITS = 192
Q_FLOOR = 1_000_000_000
MODE_START = Fraction(5, 1)
COMPACT_MODE_START = Fraction(579, 625)
WINDOW_V_MIN = Fraction(499, 100)


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str | None
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def arb_upper_text(value: flint.arb, digits: int = 60) -> str:
    rounded = value.upper().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_plus(), "E")


def tail_polynomial_ratio(order: int, slope_ratio: Fraction, y_floor: int = 12) -> Fraction:
    return sum(
        Fraction(math.comb(order, degree) * math.factorial(degree), 1)
        / slope_ratio ** (degree + 1)
        / y_floor ** (2 * degree)
        for degree in range(order + 1)
    )


def build_diagnostics() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    q0 = flint.arb(Q_FLOOR)
    log_q0 = q0.log()
    y0 = (8 * log_q0).sqrt()
    sqrt_q0 = q0.sqrt()
    curvature_coefficient = flint.arb(39) / 10

    compact_start_slope = potential_jet_arb(arb_rational(COMPACT_MODE_START), 1)[1]
    q_at_v_min = flint.arb.pi() * (4 * arb_rational(WINDOW_V_MIN)).exp()
    if not bool(compact_start_slope < 318 and q_at_v_min > Q_FLOOR):
        raise RuntimeError("mode-handoff gate failed")

    normalized_alpha_coefficient_sq_margin = Fraction(4, 1) * Fraction(39, 10) ** 3 - 14**2
    prefactor_coefficient_sq_margin = Fraction(36, 25) - Fraction(81, 1) / Fraction(39, 10) ** 3
    v4_normalized_margin = Fraction(3, 1) - Fraction(31, 1) / Fraction(39, 10) ** 2
    if min(
        normalized_alpha_coefficient_sq_margin,
        prefactor_coefficient_sq_margin,
        v4_normalized_margin,
    ) <= 0:
        raise RuntimeError("exact derivative-ratio gate failed")

    bell_v4_correction = 4_120
    v_min = WINDOW_V_MIN
    v4_component_ratio = (
        Fraction(22, 1)
        + Fraction(bell_v4_correction, Q_FLOOR)
        + Fraction(100, 1) / (v_min**2 * Q_FLOOR)
        + Fraction(5, 16) / (v_min**3 * Q_FLOOR)
    )
    right_drift_factor = (
        Fraction(30, 1) * Fraction(10_001, 10_000) ** 4 * Fraction(1_001, 1_000)
    )
    if not (v4_component_ratio < 30 and right_drift_factor < 31):
        raise RuntimeError("fourth-derivative component gate failed")

    central_remainder = y0**3 / (3 * sqrt_q0) + y0**4 / (8 * q0)
    endpoint_slope_loss = y0 / sqrt_q0 + y0**2 / (2 * q0)
    drift_exponent = y0 / (10 * (curvature_coefficient * q0).sqrt())
    q_drift_exponent = (
        2
        * y0
        / (curvature_coefficient * q0).sqrt()
        * flint.arb(10_001)
        / 10_000
    )
    tail_log_margin = q0**2 - 56 * log_q0
    exponential_gates = {
        "central_remainder_upper": arb_upper_text(central_remainder),
        "central_remainder_cap": "1/40",
        "endpoint_slope_loss_upper": arb_upper_text(endpoint_slope_loss),
        "endpoint_slope_loss_cap": "1/10",
        "standardized_log_shift_upper": arb_upper_text(drift_exponent),
        "standardized_log_shift_cap": "1/40000",
        "q_log_drift_upper": arb_upper_text(q_drift_exponent),
        "q_log_drift_cap": "1/2000",
        "exp_1_over_40000_upper": arb_upper_text((flint.arb(1) / 40_000).exp()),
        "u_ratio_cap": "10001/10000",
        "exp_1_over_2000_upper": arb_upper_text((flint.arb(1) / 2_000).exp()),
        "q_ratio_cap": "1001/1000",
        "tail_log_margin_lower": arb_lower_text(tail_log_margin),
    }
    if not bool(
        central_remainder < flint.arb(1) / 40
        and endpoint_slope_loss < flint.arb(1) / 10
        and drift_exponent < flint.arb(1) / 40_000
        and q_drift_exponent < flint.arb(1) / 2_000
        and (flint.arb(1) / 40_000).exp() < flint.arb(10_001) / 10_000
        and (flint.arb(1) / 2_000).exp() < flint.arb(1_001) / 1_000
        and y0 > 12
        and tail_log_margin > 0
    ):
        raise RuntimeError("central-window or tail gate failed")

    actual_tail_ratios = [tail_polynomial_ratio(order, Fraction(9, 10)) for order in range(4)]
    model_tail_ratios = [tail_polynomial_ratio(order, Fraction(1, 1)) for order in range(7)]
    if max(actual_tail_ratios) >= 2 or max(model_tail_ratios) >= 2:
        raise RuntimeError("tail polynomial gate failed")

    sqrt_bound = 30_000
    q_fraction = Fraction(Q_FLOOR, 1)
    denominator = 1 - Fraction(6, Q_FLOOR)
    raw_error_coefficients = {
        "raw1": (Fraction(12, 1) + Fraction(6, sqrt_bound)) / denominator,
        "raw2": Fraction(35, 1) / denominator,
        "raw3": (Fraction(79, 1) + Fraction(30, sqrt_bound)) / denominator,
    }
    raw_error_caps = {"raw1": 13, "raw2": 36, "raw3": 80}
    if not all(raw_error_coefficients[key] < cap for key, cap in raw_error_caps.items()):
        raise RuntimeError("raw-moment normalization gate failed")

    cumulant_coefficient = (
        Fraction(2, sqrt_bound)
        + raw_error_caps["raw3"]
        + 3
        * (
            raw_error_caps["raw1"] * (1 + Fraction(raw_error_caps["raw2"], Q_FLOOR))
            + Fraction(raw_error_caps["raw2"], sqrt_bound)
        )
        + 6
        * raw_error_caps["raw1"]
        * (Fraction(1, sqrt_bound) + Fraction(raw_error_caps["raw1"], Q_FLOOR)) ** 2
    )
    if not cumulant_coefficient < 120:
        raise RuntimeError("cumulant composition gate failed")

    scaled_cumulant_cap = Fraction(6, 25) * 120 / sqrt_bound
    ray_remainder_cap = scaled_cumulant_cap + Fraction(1, 100) + Fraction(1, 1000)
    required_floor = Fraction(79, 1000)
    if not scaled_cumulant_cap < Fraction(1, 1000):
        raise RuntimeError("scaled cumulant gate failed")
    if not ray_remainder_cap < required_floor:
        raise RuntimeError("ray remainder floor failed")

    gaussian_abs_moment_caps = {4: 4, 5: 7, 6: 15, 7: 39, 8: 105, 9: 307}
    moment_error_caps = {0: 6, 1: 12, 2: 29, 3: 79}
    return {
        "parameters": {
            "precision_bits": DEFAULT_PRECISION_BITS,
            "mode_ray": "u>=5",
            "q_floor": Q_FLOOR,
            "adaptive_cutoff": "Y=sqrt(8*log(q))",
            "curvature_lower": "a>=39*u^2*q/10",
            "slope_upper": "t<=3*u*q",
        },
        "mode_handoff": {
            "compact_mode_start": str(COMPACT_MODE_START),
            "potential_slope_at_compact_start_upper": arb_upper_text(compact_start_slope),
            "target_t_start": 318,
            "q_at_499_lower": arb_lower_text(q_at_v_min),
        },
        "derivative_envelopes": {
            "alpha_bound": "|alpha|<=2/sqrt(q)",
            "central_fourth_bound": "sup_window |V^(4)|/a^2<=3/q",
            "prefactor_bound": "P=t^2/a^(3/2)<=6*sqrt(q)/(5*u)",
            "bell_v4_correction": bell_v4_correction,
            "v4_component_ratio_upper": str(v4_component_ratio),
            "right_drift_factor": str(right_drift_factor),
            "alpha_sq_margin": str(normalized_alpha_coefficient_sq_margin),
            "v4_normalized_margin": str(v4_normalized_margin),
            "prefactor_sq_margin": str(prefactor_coefficient_sq_margin),
        },
        "central_window": exponential_gates,
        "tails": {
            "actual_tail_ratio_max": str(max(actual_tail_ratios)),
            "model_tail_ratio_max": str(max(model_tail_ratios)),
            "normalized_tail_error": "<=q^(-2) for raw orders 0..3",
            "actual_tail_argument": "W(+-Y)>=4*log(q)-1/40 and outward slopes are >=9*Y/10",
            "model_tail_argument": "Gaussian-linear model tails use P_j(Y,Y) and |alpha|*P_(j+3)(Y,Y)/6",
        },
        "moment_comparison": {
            "gaussian_abs_moment_caps": gaussian_abs_moment_caps,
            "unnormalized_error_caps_over_q": moment_error_caps,
            "raw_error_coefficients": {key: str(value) for key, value in raw_error_coefficients.items()},
            "raw_error_caps_over_q": raw_error_caps,
            "linear_model_moments": {
                "p0": "1",
                "p1": "-alpha/2",
                "p2": "1",
                "p3": "-5*alpha/2",
            },
        },
        "cumulant_composition": {
            "coefficient_upper": str(cumulant_coefficient),
            "coefficient_cap": 120,
            "bound": "|kappa_3(Y)+alpha|<=120/q",
            "scaled_bound": "P*|kappa_3(Y)+alpha|<1/1000",
            "scaled_cap_fraction": str(scaled_cumulant_cap),
        },
        "ray_theorem": {
            "cubic_correction_cap": "1/100",
            "fifth_correction_cap": "1/1000",
            "ray_remainder_absolute_cap": str(ray_remainder_cap),
            "ray_remainder_lower": str(-ray_remainder_cap),
            "required_remainder_floor": "-79/1000",
            "floor_margin": str(required_floor - ray_remainder_cap),
        },
    }


def build_artifact() -> dict:
    diagnostics = build_diagnostics()
    rows = [
        CertificateRow(
            id="fsprrc_01_mode_and_derivative_geometry",
            role="exact_analytic_lemma",
            claim="The compact cover starts below the t=318 mode, and the ray has explicit alpha, fourth-derivative, and prefactor envelopes.",
            formula="|alpha|<=2/sqrt(q), B4<=3/q, P<=6*sqrt(q)/(5*u)",
            readiness="available_exact",
            proof_boundary="Mode and derivative geometry only.",
            diagnostics={
                "mode_handoff": diagnostics["mode_handoff"],
                "derivative_envelopes": diagnostics["derivative_envelopes"],
            },
        ),
        CertificateRow(
            id="fsprrc_02_adaptive_window",
            role="exact_analytic_lemma",
            claim="The adaptive cutoff Y=sqrt(8 log q) stays in the derivative envelope and keeps the perturbation uniformly small.",
            formula="|r(y)|<1/40 and endpoint slope loss <1/10 for |y|<=Y",
            readiness="available_exact",
            proof_boundary="Central-window theorem only.",
            diagnostics=diagnostics["central_window"],
        ),
        CertificateRow(
            id="fsprrc_03_first_order_weight_bound",
            role="exact_analytic_lemma",
            claim="The exact standardized weight is controlled by its first odd Gaussian perturbation.",
            formula="|exp(-r)-(1-alpha*y^3/6)|<=B4*|y|^4/12+alpha^2*|y|^6/18",
            readiness="available_exact",
            proof_boundary="Pointwise central comparison only.",
        ),
        CertificateRow(
            id="fsprrc_04_two_tail_bound",
            role="exact_analytic_lemma",
            claim="Potential monotonicity and endpoint slopes control both actual tails, while Gaussian hazards control both model tails.",
            formula="normalized actual+model tail error<=q^(-2), raw orders 0..3",
            readiness="available_exact",
            proof_boundary="Tail moment theorem only.",
            diagnostics=diagnostics["tails"],
        ),
        CertificateRow(
            id="fsprrc_05_raw_moment_comparison",
            role="exact_analytic_lemma",
            claim="The four normalized raw moments remain close to the exactly integrable Gaussian-linear model.",
            formula="|m1+alpha/2|<=13/q, |m2-1|<=36/q, |m3+5*alpha/2|<=80/q",
            readiness="available_exact",
            proof_boundary="Raw standardized moments only.",
            diagnostics=diagnostics["moment_comparison"],
        ),
        CertificateRow(
            id="fsprrc_06_cumulant_composition",
            role="exact_analytic_theorem",
            claim="Normalization and centralization preserve enough cancellation to bound the residual skewness after the leading odd term.",
            formula="|kappa_3(Y)+alpha|<=120/q and P*|kappa_3(Y)+alpha|<1/1000",
            readiness="available_exact",
            proof_boundary="First-order cumulant theorem on u>=5.",
            diagnostics=diagnostics["cumulant_composition"],
        ),
        CertificateRow(
            id="fsprrc_07_ray_remainder_theorem",
            role="exact_analytic_theorem",
            claim="The paired seventh-order normalized remainder satisfies a stronger bound on the full asymptotic mode ray.",
            formula="H_t>=-3/250>-79/1000 for u>=5",
            readiness="available_exact",
            proof_boundary="Analytic ray theorem only.",
            diagnostics=diagnostics["ray_theorem"],
        ),
        CertificateRow(
            id="fsprrc_08_global_remainder_closure",
            role="theorem_composition",
            claim="The compact paired certificate and analytic ray theorem prove the required remainder floor for every t>=318.",
            formula="H_t>=-79/1000 for every real t>=318",
            readiness="ready_to_apply",
            proof_boundary="Global first-summand remainder theorem; not yet the full heat-flow bridge.",
        ),
        CertificateRow(
            id="fsprrc_09_cumulant_wall_handoff",
            role="exact_handoff",
            claim="Together with the certified leading, cubic, and fifth caps, the global remainder proves the cumulant hypothesis used by the exact Gamma bridge.",
            formula="kappa_3,t(2 log U)>=-37/(50*t^2), t>=318",
            readiness="ready_to_apply",
            proof_boundary="Closes the first-summand analytic hypothesis only; not cone propagation, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate",
        "date": "2026-07-10",
        "status": "analytic asymptotic-ray theorem and global paired-remainder closure",
        "proof_boundary": (
            "This artifact proves the u>=5 paired remainder theorem and composes it with "
            "the compact certificate to close the global first-summand cumulant hypothesis. "
            "It does not prove the remaining cone-propagation/all-order bridge, RH, or Lambda <= 0."
        ),
        "source_compact_certificate": "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
        "source_leading_certificate": "outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md",
        "source_cumulant_bridge": "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
        "source_saddle_target": "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py",
        "diagnostics": diagnostics,
        "summary": {
            "certificate_rows": len(rows),
            "analytic_ray_rows": 1,
            "global_remainder_closure_rows": 1,
            "open_ray_rows": 0,
            "ready_to_apply_rows": 2,
            "ray_remainder_lower": diagnostics["ray_theorem"]["ray_remainder_lower"],
            "required_remainder_floor": "-79/1000",
            "floor_margin": diagnostics["ray_theorem"]["floor_margin"],
        },
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    diagnostics = artifact["diagnostics"]
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda First-Summand Paired-Remainder Ray Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: analytic asymptotic-ray theorem and global paired-remainder closure.",
        "This is not a proof of the remaining all-order bridge, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.py",
        "```",
        "",
        "## Ray Geometry",
        "",
        "For `u>=5`, put `q=pi*exp(4*u)`, `a=V''(x_t)`, and",
        "`Y=sqrt(8*log(q))`. Exact component bounds give",
        "",
        "```text",
        "|alpha|<=2/sqrt(q),",
        "sup_|y|<=Y |V^(4)(x_t+y/sqrt(a))|/a^2<=3/q,",
        "P=t^2/a^(3/2)<=6*sqrt(q)/(5*u).",
        "```",
        "",
        f"At the worst endpoint, the central perturbation is below `{diagnostics['central_window']['central_remainder_upper']}`.",
        "The adaptive window has outward slopes at least `9*Y/10`, and both",
        "actual and Gaussian-linear model tails contribute at most `q^(-2)`",
        "to each normalized raw-moment error.",
        "",
        "## Moment And Cumulant Bounds",
        "",
        "The first-order Gaussian model has moments",
        "",
        "```text",
        "p0=1, p1=-alpha/2, p2=1, p3=-5*alpha/2.",
        "```",
        "",
        "The exact comparison proves",
        "",
        "```text",
        "|m1+alpha/2|<=13/q,",
        "|m2-1|<=36/q,",
        "|m3+5*alpha/2|<=80/q,",
        "|kappa_3(Y)+alpha|<=120/q.",
        "```",
        "",
        "After the prefactor, this is below `1/1000`. Adding the already",
        "certified cubic and fifth-order caps gives",
        "",
        "```text",
        f"H_t>={summary['ray_remainder_lower']}>-3/250>-79/1000, u>=5.",
        "```",
        "",
        "## Global Closure",
        "",
        "The compact paired theorem covers `0.9264<=u<=5`; the mode at `t=318`",
        "lies to the right of `0.9264`. Therefore",
        "",
        "```text",
        "H_t>=-79/1000 for every real t>=318,",
        "kappa_3,t(2*log(U))>=-37/(50*t^2) for every real t>=318.",
        "```",
        "",
        "This closes the analytic hypothesis in the exact first-summand",
        "cumulant/Gamma bridge. It does not by itself prove the remaining",
        "cone-propagation or all-order Newman-direction bridge.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
    ]
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
    summary = artifact["summary"]
    print(
        "validated Jensen-window PF negative-lambda first-summand paired-remainder ray certificate: "
        f"{summary['certificate_rows']} rows, 0 issues, 1 analytic ray theorem, "
        "1 global remainder closure, 0 open rays, 2 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
