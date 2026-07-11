#!/usr/bin/env python3
"""Build a worst-row far-tail split certificate for interval integration."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FIRST_OMITTED_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json"
)
DEFAULT_QUADRATURE_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_SPLIT_Y = 200
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PRECISION_BITS = 512
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_C0_LOWER = "0.44"

getcontext().prec = 120


@dataclass(frozen=True)
class MatrixRow:
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


def arb_upper_text(value: flint.arb, digits: int = 40) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 40) -> str:
    return value.lower().str(digits, radius=False).replace("e", "E")


def arb_mid_text(value: flint.arb, digits: int = 40) -> str:
    return value.str(digits, radius=False).replace("e", "E")


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def abs_upper(value: flint.arb) -> flint.arb:
    lower_abs = abs(flint.arb(value.lower()))
    upper_abs = abs(flint.arb(value.upper()))
    return lower_abs if lower_abs > upper_abs else upper_abs


def source_paths(first_omitted_path: Path, quadrature_route_path: Path) -> dict[str, str]:
    return {
        "first_omitted_json": first_omitted_path.relative_to(REPO_ROOT).as_posix(),
        "first_omitted_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
        ),
        "quadrature_route_json": quadrature_route_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"
        ),
        "node_c0_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md",
        "phi_tail_grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
        "cancellation_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md"
        ),
        "intervalization_target": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def build_phi_bounds(x0: flint.arb, term_count: int, c0_lower: flint.arb) -> dict:
    pi = flint.arb.pi()
    exp4 = (flint.arb(4) * x0).exp()
    value_sum = flint.arb(0)
    x_phip_sum = flint.arb(0)
    for n in range(1, term_count + 1):
        n2 = flint.arb(n * n)
        n4 = flint.arb(n * n * n * n)
        q = pi * n2
        e5 = (flint.arb(5) * x0).exp()
        e9 = (flint.arb(9) * x0).exp()
        decay = (-q * exp4).exp()
        value_sum += (flint.arb(2) * pi * pi * n4 * e9 + flint.arb(3) * pi * n2 * e5) * decay
        prefactor_abs = flint.arb(2) * pi * pi * n4 * e9 + flint.arb(3) * pi * n2 * e5
        prefactor_derivative_abs = flint.arb(18) * pi * pi * n4 * e9 + flint.arb(15) * pi * n2 * e5
        x_phip_sum += x0 * (prefactor_derivative_abs + flint.arb(4) * q * exp4 * prefactor_abs) * decay
    return {
        "phi_abs_bound_at_split": value_sum,
        "x_phip_abs_bound_at_split": x_phip_sum,
        "phi_over_c0_abs_bound_at_split": value_sum / c0_lower,
        "x_phip_over_2c0_abs_bound_at_split": x_phip_sum / (flint.arb(2) * c0_lower),
    }


def build_tail_moment_bounds(
    ratios: list[flint.arb],
    split_y: flint.arb,
    alpha: flint.arb,
    u: flint.arb,
    polynomial_m: int,
) -> dict:
    gamma_alpha = (alpha + 1).gamma()
    tail_mass = split_y.gamma_upper(alpha + 1) / gamma_alpha
    value_poly = flint.arb(0)
    derivative_poly = flint.arb(0)
    moment_rows: list[dict] = []
    for j in range(polynomial_m + 1):
        tail_moment = split_y.gamma_upper(alpha + flint.arb(j) + 1) / gamma_alpha
        ratio_abs = abs_upper(ratios[j])
        value_contribution = ratio_abs * (u**j) * tail_moment
        value_poly += value_contribution
        derivative_contribution = flint.arb(0)
        if j >= 1:
            derivative_contribution = flint.arb(j) * ratio_abs * (u**j) * tail_moment
            derivative_poly += derivative_contribution
        if j in {0, 1, 2, 10, 20}:
            moment_rows.append(
                {
                    "j": j,
                    "tail_moment_upper": arb_upper_text(tail_moment),
                    "ratio_abs_upper": arb_upper_text(ratio_abs),
                    "value_contribution_upper": arb_upper_text(value_contribution),
                    "derivative_contribution_upper": arb_upper_text(derivative_contribution),
                    "proof_boundary": "Sampled polynomial tail-moment contribution only.",
                }
            )
    return {
        "tail_mass": tail_mass,
        "value_polynomial_tail_bound": value_poly,
        "derivative_polynomial_tail_bound": derivative_poly,
        "moment_sample_rows": moment_rows,
    }


def build_diagnostics(first_omitted: dict, quadrature_route: dict) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    _ratio_rows, ratios = build_ratio_rows(2 * (DEFAULT_POLYNOMIAL_M + 1), DEFAULT_RATIO_CUTOFF_N)
    split_y = flint.arb(DEFAULT_SPLIT_Y)
    T = flint.arb(DEFAULT_T)
    u = flint.arb(1) / T
    alpha = flint.arb(DEFAULT_ALPHA)
    x0 = (split_y / T).sqrt()
    c0_lower = flint.arb(DEFAULT_C0_LOWER)
    pi = flint.arb.pi()
    exp4 = (flint.arb(4) * x0).exp()
    value_monotone_margin = flint.arb(4) * pi * exp4 - flint.arb(9)
    derivative_monotone_margin = flint.arb(4) * pi * exp4 - flint.arb(13) - (flint.arb(1) / x0)
    if not bool(value_monotone_margin > 0):
        raise RuntimeError("value Phi majorant monotonicity margin did not certify positive")
    if not bool(derivative_monotone_margin > 0):
        raise RuntimeError("derivative Phi majorant monotonicity margin did not certify positive")

    phi = build_phi_bounds(x0, DEFAULT_PHI_TERM_COUNT, c0_lower)
    moments = build_tail_moment_bounds(ratios, split_y, alpha, u, DEFAULT_POLYNOMIAL_M)
    value_phi_tail = phi["phi_over_c0_abs_bound_at_split"] * moments["tail_mass"]
    derivative_phi_tail = phi["x_phip_over_2c0_abs_bound_at_split"] * moments["tail_mass"]
    value_total = value_phi_tail + moments["value_polynomial_tail_bound"]
    derivative_total = derivative_phi_tail + moments["derivative_polynomial_tail_bound"]
    value_scaled = value_total / (u**3)
    derivative_scaled = derivative_total / (u**2)
    value_denom_lower = flint.arb(first_omitted["summary"]["minimum_value_denominator_lower"])
    derivative_denom_lower = flint.arb(first_omitted["summary"]["minimum_derivative_denominator_lower"])
    value_ratio = value_scaled / value_denom_lower
    derivative_ratio = derivative_scaled / derivative_denom_lower
    quadrature_cap = flint.arb(quadrature_route["summary"]["quadrature_ratio_radius_cap"])
    value_unscaled_cap = flint.arb(quadrature_route["summary"]["value_unscaled_expectation_error_cap"])
    derivative_unscaled_cap = flint.arb(quadrature_route["summary"]["derivative_unscaled_expectation_error_cap"])
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "u": "1/10000",
        "alpha": DEFAULT_ALPHA,
        "split_y": DEFAULT_SPLIT_Y,
        "split_x": arb_upper_text(x0),
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "polynomial_M": DEFAULT_POLYNOMIAL_M,
        "precision_bits": DEFAULT_PRECISION_BITS,
        "c0_lower": DEFAULT_C0_LOWER,
        "value_majorant_monotonicity_margin_lower": arb_lower_text(value_monotone_margin),
        "derivative_majorant_monotonicity_margin_lower": arb_lower_text(derivative_monotone_margin),
        "phi_abs_bound_at_split": arb_upper_text(phi["phi_abs_bound_at_split"]),
        "x_phip_abs_bound_at_split": arb_upper_text(phi["x_phip_abs_bound_at_split"]),
        "tail_mass_upper": arb_upper_text(moments["tail_mass"]),
        "value_phi_tail_bound": arb_upper_text(value_phi_tail),
        "derivative_phi_tail_bound": arb_upper_text(derivative_phi_tail),
        "value_polynomial_tail_bound": arb_upper_text(moments["value_polynomial_tail_bound"]),
        "derivative_polynomial_tail_bound": arb_upper_text(moments["derivative_polynomial_tail_bound"]),
        "value_total_tail_bound": arb_upper_text(value_total),
        "derivative_total_tail_bound": arb_upper_text(derivative_total),
        "value_scaled_tail_bound": arb_upper_text(value_scaled),
        "derivative_scaled_tail_bound": arb_upper_text(derivative_scaled),
        "value_ratio_to_first_omitted_upper": arb_upper_text(value_ratio),
        "derivative_ratio_to_first_omitted_upper": arb_upper_text(derivative_ratio),
        "quadrature_ratio_radius_cap": quadrature_route["summary"]["quadrature_ratio_radius_cap"],
        "value_tail_ratio_below_quadrature_cap": bool(value_ratio < quadrature_cap),
        "derivative_tail_ratio_below_quadrature_cap": bool(derivative_ratio < quadrature_cap),
        "value_unscaled_tail_below_quadrature_cap": bool(value_total < value_unscaled_cap),
        "derivative_unscaled_tail_below_quadrature_cap": bool(derivative_total < derivative_unscaled_cap),
        "value_tail_to_unscaled_cap_ratio_upper": arb_upper_text(value_total / value_unscaled_cap),
        "derivative_tail_to_unscaled_cap_ratio_upper": arb_upper_text(derivative_total / derivative_unscaled_cap),
        "moment_sample_rows": moments["moment_sample_rows"],
        "remaining_compact_interval": "0<=y<=200",
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrftsc_01_split_and_monotonicity",
            role="arb_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "For the worst row T=10000, F_21, split the normalized Gamma expectation at y=200. "
                "Arb margins certify the finite n<=30 Phi and x Phi' majorants are decreasing beyond the split."
            ),
            diagnostics={
                "split_y": diagnostics["split_y"],
                "split_x": diagnostics["split_x"],
                "value_majorant_monotonicity_margin_lower": diagnostics[
                    "value_majorant_monotonicity_margin_lower"
                ],
                "derivative_majorant_monotonicity_margin_lower": diagnostics[
                    "derivative_majorant_monotonicity_margin_lower"
                ],
                "phi_term_count": diagnostics["phi_term_count"],
            },
            source_artifacts=[paths["node_c0_note"], paths["cancellation_grid_note"]],
            proof_boundary="Far-tail split monotonicity certificate only; not compact-interval integration.",
        ),
        MatrixRow(
            id="nlrgwrftsc_02_finite_phi_tail_bound",
            role="arb_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "The finite n<=30 Phi and derivative-core contributions on y>=200 are bounded by the "
                "monotone majorant at the split times the Arb upper Gamma tail mass."
            ),
            diagnostics={
                "tail_mass_upper": diagnostics["tail_mass_upper"],
                "phi_abs_bound_at_split": diagnostics["phi_abs_bound_at_split"],
                "x_phip_abs_bound_at_split": diagnostics["x_phip_abs_bound_at_split"],
                "value_phi_tail_bound": diagnostics["value_phi_tail_bound"],
                "derivative_phi_tail_bound": diagnostics["derivative_phi_tail_bound"],
                "c0_lower": diagnostics["c0_lower"],
            },
            source_artifacts=[paths["node_c0_note"], paths["phi_tail_grid_note"]],
            proof_boundary="Finite n<=30 far-tail Phi bound only; not the omitted n>30 source or compact interval.",
        ),
        MatrixRow(
            id="nlrgwrftsc_03_polynomial_tail_moments",
            role="arb_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "The cancellation polynomial tail on y>=200 is bounded by Arb upper incomplete-Gamma "
                "moments and Arb coefficient-ratio enclosures."
            ),
            diagnostics={
                "value_polynomial_tail_bound": diagnostics["value_polynomial_tail_bound"],
                "derivative_polynomial_tail_bound": diagnostics["derivative_polynomial_tail_bound"],
                "moment_sample_rows": diagnostics["moment_sample_rows"],
            },
            source_artifacts=[paths["first_omitted_json"], paths["first_omitted_note"]],
            proof_boundary="Polynomial far-tail moment certificate only; not compact-interval integration.",
        ),
        MatrixRow(
            id="nlrgwrftsc_04_tail_vs_quadrature_cap",
            role="arb_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "The combined finite n<=30 value and derivative far tails on y>=200 are far below the "
                "unscaled and first-omitted ratio caps required by the quadrature route matrix."
            ),
            diagnostics={
                "value_total_tail_bound": diagnostics["value_total_tail_bound"],
                "derivative_total_tail_bound": diagnostics["derivative_total_tail_bound"],
                "value_scaled_tail_bound": diagnostics["value_scaled_tail_bound"],
                "derivative_scaled_tail_bound": diagnostics["derivative_scaled_tail_bound"],
                "value_ratio_to_first_omitted_upper": diagnostics["value_ratio_to_first_omitted_upper"],
                "derivative_ratio_to_first_omitted_upper": diagnostics[
                    "derivative_ratio_to_first_omitted_upper"
                ],
                "quadrature_ratio_radius_cap": diagnostics["quadrature_ratio_radius_cap"],
                "value_tail_ratio_below_quadrature_cap": diagnostics[
                    "value_tail_ratio_below_quadrature_cap"
                ],
                "derivative_tail_ratio_below_quadrature_cap": diagnostics[
                    "derivative_tail_ratio_below_quadrature_cap"
                ],
                "value_tail_to_unscaled_cap_ratio_upper": diagnostics[
                    "value_tail_to_unscaled_cap_ratio_upper"
                ],
                "derivative_tail_to_unscaled_cap_ratio_upper": diagnostics[
                    "derivative_tail_to_unscaled_cap_ratio_upper"
                ],
            },
            source_artifacts=[paths["quadrature_route_json"], paths["quadrature_route_note"]],
            proof_boundary="Far-tail comparison only; it does not certify the compact interval 0<=y<=200.",
        ),
        MatrixRow(
            id="nlrgwrftsc_05_interval_integration_handoff",
            role="intervalization_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "After this far-tail certificate, the independent interval-integration route is reduced to "
                "certifying the compact worst-row integral on 0<=y<=200 plus already separated tail sources."
            ),
            diagnostics={
                "remaining_compact_interval": diagnostics["remaining_compact_interval"],
                "retired_component": "finite n<=30 cancellation-reduced far tail y>=200 for T=10000, F_21",
                "remaining_sources": [
                    "Arb or interval-adaptive integration on 0<=y<=200",
                    "rounding and aggregation with the other intervalization sources",
                    "all recorded rows/orders, not just the worst row",
                    "finite-grid to full-collar coverage",
                ],
            },
            source_artifacts=[paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Handoff only; not a finite-grid interval certificate or uniform collar theorem.",
        ),
        MatrixRow(
            id="nlrgwrftsc_06_far_tail_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The y>=200 far-tail certificate proves the quadrature-remainder source for the worst row.",
            gap=(
                "It certifies only the continuum far tail after the split. The compact interval 0<=y<=200 "
                "still needs interval integration or a derivative/remainder theorem."
            ),
            source_artifacts=[paths["quadrature_route_note"], paths["dependency_graph"]],
            proof_boundary="Rejected promotion only; not quadrature remainder, RH, or Lambda <= 0.",
        ),
        MatrixRow(
            id="nlrgwrftsc_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted worst-row interval-integration certificate may use this far-tail split only "
                "together with a compact-interval enclosure and aggregation ledger."
            ),
            source_artifacts=[paths["dependency_graph"]],
            proof_boundary="Proof-hygiene gate only; not a finite-grid interval certificate, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(first_omitted_path: Path, quadrature_route_path: Path) -> dict:
    first_omitted = load_json(first_omitted_path)
    quadrature_route = load_json(quadrature_route_path)
    paths = source_paths(first_omitted_path, quadrature_route_path)
    diagnostics = build_diagnostics(first_omitted, quadrature_route)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "split_y": diagnostics["split_y"],
        "split_x": diagnostics["split_x"],
        "phi_term_count": diagnostics["phi_term_count"],
        "polynomial_M": diagnostics["polynomial_M"],
        "precision_bits": diagnostics["precision_bits"],
        "tail_mass_upper": diagnostics["tail_mass_upper"],
        "value_total_tail_bound": diagnostics["value_total_tail_bound"],
        "derivative_total_tail_bound": diagnostics["derivative_total_tail_bound"],
        "value_scaled_tail_bound": diagnostics["value_scaled_tail_bound"],
        "derivative_scaled_tail_bound": diagnostics["derivative_scaled_tail_bound"],
        "value_ratio_to_first_omitted_upper": diagnostics["value_ratio_to_first_omitted_upper"],
        "derivative_ratio_to_first_omitted_upper": diagnostics[
            "derivative_ratio_to_first_omitted_upper"
        ],
        "quadrature_ratio_radius_cap": diagnostics["quadrature_ratio_radius_cap"],
        "value_tail_ratio_below_quadrature_cap": diagnostics["value_tail_ratio_below_quadrature_cap"],
        "derivative_tail_ratio_below_quadrature_cap": diagnostics[
            "derivative_tail_ratio_below_quadrature_cap"
        ],
        "value_unscaled_tail_below_quadrature_cap": diagnostics[
            "value_unscaled_tail_below_quadrature_cap"
        ],
        "derivative_unscaled_tail_below_quadrature_cap": diagnostics[
            "derivative_unscaled_tail_below_quadrature_cap"
        ],
        "value_tail_to_unscaled_cap_ratio_upper": diagnostics["value_tail_to_unscaled_cap_ratio_upper"],
        "derivative_tail_to_unscaled_cap_ratio_upper": diagnostics[
            "derivative_tail_to_unscaled_cap_ratio_upper"
        ],
        "remaining_compact_interval": diagnostics["remaining_compact_interval"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst-row interval-integration route now has an Arb-certified far-tail split: on y>=200, "
            "finite n<=30 Phi majorants are monotone decreasing, the polynomial part is bounded by upper "
            "incomplete-Gamma moments, and the combined value and derivative tails are far below the "
            "quadrature route matrix caps. This retires only the far-tail part; compact integration on "
            "0<=y<=200, aggregation, all-row coverage, and grid-to-collar coverage remain open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate",
        "date": "2026-07-07",
        "status": "worst-row far-tail split certificate",
        "source_first_omitted_denominator_certificate": paths["first_omitted_note"],
        "source_first_omitted_denominator_json": paths["first_omitted_json"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_quadrature_remainder_route_json": paths["quadrature_route_json"],
        "source_node_c0_range_certificate": paths["node_c0_note"],
        "source_phi_tail_grid_certificate": paths["phi_tail_grid_note"],
        "source_cancellation_reduced_grid_scout": paths["cancellation_grid_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py"
        ),
        "proof_boundary": (
            "Worst-row far-tail split certificate only. It certifies the finite n<=30 cancellation-reduced "
            "continuum tail y>=200 for T=10000, F_21 against the quadrature route caps, but it does not "
            "integrate the compact interval 0<=y<=200, does not aggregate rounding, does not cover all rows "
            "or quadrature orders, does not prove a finite-grid interval certificate, does not prove a "
            "uniform collar theorem, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "Only the worst-row finite n<=30 far tail y>=200 is certified here.",
            "The compact interval 0<=y<=200 remains open.",
            "The omitted n>30 source, aggregation, all-row coverage, and grid-to-collar coverage remain separate.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, split y={summary['split_y']}, "
        "2 tail ratios below quadrature cap, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Far-Tail Split Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row far-tail split certificate. This is not a proof",
        "of a compact interval-integration certificate, finite-grid interval",
        "certificate, uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate`.",
        "",
        "Proof boundary: this artifact certifies only the finite `n<=30`",
        "cancellation-reduced continuum far tail `y>=200` for `T=10000`,",
        "`F_21`. It does not integrate the compact interval `0<=y<=200`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Tail Bounds",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"split y: {summary['split_y']}",
        f"split x: {summary['split_x']}",
        f"tail mass upper: {summary['tail_mass_upper']}",
        f"value total tail bound: {summary['value_total_tail_bound']}",
        f"derivative total tail bound: {summary['derivative_total_tail_bound']}",
        f"value scaled tail bound: {summary['value_scaled_tail_bound']}",
        f"derivative scaled tail bound: {summary['derivative_scaled_tail_bound']}",
        f"value ratio to first omitted upper: {summary['value_ratio_to_first_omitted_upper']}",
        f"derivative ratio to first omitted upper: {summary['derivative_ratio_to_first_omitted_upper']}",
        f"quadrature ratio radius cap: {summary['quadrature_ratio_radius_cap']}",
        f"value tail ratio below quadrature cap: {summary['value_tail_ratio_below_quadrature_cap']}",
        f"derivative tail ratio below quadrature cap: {summary['derivative_tail_ratio_below_quadrature_cap']}",
        f"remaining compact interval: {summary['remaining_compact_interval']}",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_first_omitted_denominator_certificate"],
        artifact["source_first_omitted_denominator_json"],
        artifact["source_quadrature_remainder_route_matrix"],
        artifact["source_quadrature_remainder_route_json"],
        artifact["source_node_c0_range_certificate"],
        artifact["source_phi_tail_grid_certificate"],
        artifact["source_cancellation_reduced_grid_scout"],
        artifact["source_intervalization_target"],
        artifact["source_uniform_remainder_target"],
        artifact["source_dependency_graph"],
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-omitted-json", type=Path, default=DEFAULT_FIRST_OMITTED_JSON)
    parser.add_argument("--quadrature-route-json", type=Path, default=DEFAULT_QUADRATURE_ROUTE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    first_omitted_path = (
        args.first_omitted_json if args.first_omitted_json.is_absolute() else REPO_ROOT / args.first_omitted_json
    )
    quadrature_route_path = (
        args.quadrature_route_json
        if args.quadrature_route_json.is_absolute()
        else REPO_ROOT / args.quadrature_route_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(first_omitted_path, quadrature_route_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
