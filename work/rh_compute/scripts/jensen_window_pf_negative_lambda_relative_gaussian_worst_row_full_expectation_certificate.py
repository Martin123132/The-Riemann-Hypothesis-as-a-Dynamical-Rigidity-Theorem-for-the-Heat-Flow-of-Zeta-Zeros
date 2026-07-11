#!/usr/bin/env python3
"""Build a complete worst-row relative-Gaussian expectation certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_COMPACT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.json"
)
DEFAULT_FAR_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json"
)
DEFAULT_FIRST_OMITTED_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_TAIL_START_N = 31
DEFAULT_TAIL_SPLIT_X = 1
DEFAULT_C0_LOWER_TARGET = "0.44"
DEFAULT_PRECISION_BITS = 8192


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


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 70) -> str:
    return value.lower().str(digits, radius=False).replace("e", "E")


def arb_upper_text(value: flint.arb, digits: int = 70) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def abs_upper(value: flint.arb) -> flint.arb:
    lower_abs = abs(flint.arb(value.lower()))
    upper_abs = abs(flint.arb(value.upper()))
    return lower_abs if lower_abs > upper_abs else upper_abs


def parse_arb(text: str) -> flint.arb:
    return flint.arb(text.replace("E", "e"))


def source_paths(compact_path: Path, far_tail_path: Path, first_omitted_path: Path) -> dict[str, str]:
    return {
        "compact_json": compact_path.relative_to(REPO_ROOT).as_posix(),
        "compact_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.md"
        ),
        "endpoint_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.md"
        ),
        "far_tail_json": far_tail_path.relative_to(REPO_ROOT).as_posix(),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "first_omitted_json": first_omitted_path.relative_to(REPO_ROOT).as_posix(),
        "first_omitted_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
        ),
        "coefficient_core_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md"
        ),
        "phi_tail_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md"
        ),
        "node_c0_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md"
        ),
        "finite_plus_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md"
        ),
        "quadrature_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"
        ),
        "intervalization_target": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "formal_tail_obstruction": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"
        ),
    }


def geometric_tail(first_n: int, degree: int, first_term: flint.arb) -> tuple[flint.arb, flint.arb]:
    pi = flint.arb.pi()
    ratio = (
        (flint.arb(first_n + 1) / flint.arb(first_n)) ** degree
        * (-pi * flint.arb(2 * first_n + 1)).exp()
    )
    if ratio.contains(1) or not bool(ratio < 1):
        raise RuntimeError("n-tail geometric ratio did not certify below one")
    return first_term / (1 - ratio), ratio


def global_n_tail_bounds() -> dict[str, flint.arb | bool]:
    """Arb n>=31 bounds, extended from x<=1 to x>=0 by monotonicity."""
    pi = flint.arb.pi()
    n = flint.arb(DEFAULT_TAIL_START_N)
    n2 = n * n
    exp4 = flint.arb(4).exp()
    exp5 = flint.arb(5).exp()
    exp9 = flint.arb(9).exp()
    value_first = (
        2 * pi * pi * n2 * n2 * exp9 + 3 * pi * n2 * exp5
    ) * (-pi * n2).exp()
    derivative_first = (
        18 * pi * pi * n2 * n2 * exp9
        + 15 * pi * n2 * exp5
        + 4
        * pi
        * n2
        * exp4
        * (2 * pi * pi * n2 * n2 * exp9 + 3 * pi * n2 * exp5)
    ) * (-pi * n2).exp()
    c0_first = (2 * pi * pi * n2 * n2 + 3 * pi * n2) * (-pi * n2).exp()
    value_bound, value_ratio = geometric_tail(DEFAULT_TAIL_START_N, 4, value_first)
    derivative_bound, derivative_ratio = geometric_tail(
        DEFAULT_TAIL_START_N, 6, derivative_first
    )
    c0_bound, c0_ratio = geometric_tail(DEFAULT_TAIL_START_N, 4, c0_first)

    value_monotonicity_margin = 4 * pi * n2 * exp4 - 9
    derivative_x_monotonicity_margin = 4 * pi * n2 * exp4 - 14
    global_extension_certified = bool(
        value_monotonicity_margin > 0 and derivative_x_monotonicity_margin > 0
    )
    if not global_extension_certified:
        raise RuntimeError("global x>=1 n-tail monotonicity did not certify")
    return {
        "value_first": value_first,
        "derivative_first": derivative_first,
        "c0_first": c0_first,
        "value_bound": value_bound,
        "derivative_bound": derivative_bound,
        "c0_bound": c0_bound,
        "value_ratio": value_ratio,
        "derivative_ratio": derivative_ratio,
        "c0_ratio": c0_ratio,
        "value_monotonicity_margin": value_monotonicity_margin,
        "derivative_x_monotonicity_margin": derivative_x_monotonicity_margin,
        "global_extension_certified": global_extension_certified,
    }


def build_diagnostics(compact: dict, far_tail: dict, first_omitted: dict) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    compact_summary = compact["summary"]
    compact_diagnostics = compact["diagnostics"]
    far_summary = far_tail["summary"]
    far_phi = far_tail["matrix_rows"][1]["diagnostics"]
    n_tail = global_n_tail_bounds()

    c0_finite = parse_arb(compact_diagnostics["finite_c0_ball"])
    c0_finite_lower = flint.arb(c0_finite.lower())
    c0_finite_upper = flint.arb(c0_finite.upper())
    c0_tail = n_tail["c0_bound"]
    assert isinstance(c0_tail, flint.arb)
    c0_full_lower = c0_finite_lower - c0_tail
    c0_target = flint.arb(DEFAULT_C0_LOWER_TARGET)
    c0_lowers_certified = bool(c0_finite_lower > c0_target and c0_full_lower > c0_target)
    if not c0_lowers_certified:
        raise RuntimeError("finite/full Phi(0) lower bounds did not certify above 0.44")

    compact_mass = parse_arb(compact_diagnostics["compact_gamma_mass_ball"])
    disk_majorant = flint.arb(compact_summary["finite_phi_complex_disk_majorant_upper"])
    q = flint.arb(compact_summary["x_to_disk_radius_ratio_upper"])
    compact_phi_abs_integral_bound = disk_majorant * compact_mass
    compact_x_phip_over_2_abs_integral_bound = (
        disk_majorant * q * compact_mass / (2 * ((1 - q) ** 2))
    )
    tail_mass = flint.arb(far_phi["tail_mass_upper"])
    far_phi_abs_integral_bound = flint.arb(far_phi["phi_abs_bound_at_split"]) * tail_mass
    far_x_phip_over_2_abs_integral_bound = (
        flint.arb(far_phi["x_phip_abs_bound_at_split"]) * tail_mass / 2
    )
    finite_phi_abs_integral_bound = (
        compact_phi_abs_integral_bound + far_phi_abs_integral_bound
    )
    finite_x_phip_over_2_abs_integral_bound = (
        compact_x_phip_over_2_abs_integral_bound
        + far_x_phip_over_2_abs_integral_bound
    )

    value_numerator_tail = n_tail["value_bound"]
    derivative_numerator_tail = n_tail["derivative_bound"]
    assert isinstance(value_numerator_tail, flint.arb)
    assert isinstance(derivative_numerator_tail, flint.arb)
    value_normalized_correction = (
        value_numerator_tail / c0_full_lower
        + finite_phi_abs_integral_bound
        * c0_tail
        / (c0_finite_lower * c0_full_lower)
    )
    derivative_normalized_correction = (
        derivative_numerator_tail / (2 * c0_full_lower)
        + finite_x_phip_over_2_abs_integral_bound
        * c0_tail
        / (c0_finite_lower * c0_full_lower)
    )

    compact_value = parse_arb(compact_summary["value_compact_total_ball"])
    compact_derivative = parse_arb(compact_summary["derivative_compact_total_ball"])
    far_value_bound = flint.arb(far_summary["value_total_tail_bound"])
    far_derivative_bound = flint.arb(far_summary["derivative_total_tail_bound"])
    finite_entire_value = compact_value + flint.arb(0, far_value_bound)
    finite_entire_derivative = compact_derivative + flint.arb(0, far_derivative_bound)
    full_value = finite_entire_value + flint.arb(0, value_normalized_correction)
    full_derivative = finite_entire_derivative + flint.arb(0, derivative_normalized_correction)

    quadrature_ratio_cap = flint.arb(far_summary["quadrature_ratio_radius_cap"])
    value_error_cap = flint.arb(compact_summary["value_unscaled_expectation_error_cap"])
    derivative_error_cap = flint.arb(
        compact_summary["derivative_unscaled_expectation_error_cap"]
    )
    value_first_omitted = value_error_cap / quadrature_ratio_cap
    derivative_first_omitted = derivative_error_cap / quadrature_ratio_cap
    value_ratio = abs_upper(full_value) / value_first_omitted
    derivative_ratio = abs_upper(full_derivative) / derivative_first_omitted
    value_margin = 1 - value_ratio
    derivative_margin = 1 - derivative_ratio
    both_below_one = bool(
        full_value < 0
        and full_derivative < 0
        and value_ratio < 1
        and derivative_ratio < 1
    )
    if not both_below_one:
        raise RuntimeError("full worst-row expectation ratios did not certify below one")

    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "precision_bits": DEFAULT_PRECISION_BITS,
        "tail_start_n": DEFAULT_TAIL_START_N,
        "tail_split_x": DEFAULT_TAIL_SPLIT_X,
        "value_n_tail_first_term_upper": arb_upper_text(n_tail["value_first"]),
        "derivative_n_tail_first_term_upper": arb_upper_text(n_tail["derivative_first"]),
        "c0_n_tail_first_term_upper": arb_upper_text(n_tail["c0_first"]),
        "value_n_tail_geometric_ratio_upper": arb_upper_text(n_tail["value_ratio"]),
        "derivative_n_tail_geometric_ratio_upper": arb_upper_text(
            n_tail["derivative_ratio"]
        ),
        "c0_n_tail_geometric_ratio_upper": arb_upper_text(n_tail["c0_ratio"]),
        "value_n_tail_global_bound_upper": arb_upper_text(value_numerator_tail),
        "derivative_n_tail_global_bound_upper": arb_upper_text(derivative_numerator_tail),
        "c0_n_tail_bound_upper": arb_upper_text(c0_tail),
        "value_x_ge_1_monotonicity_margin_lower": arb_lower_text(
            n_tail["value_monotonicity_margin"]
        ),
        "derivative_x_ge_1_monotonicity_margin_lower": arb_lower_text(
            n_tail["derivative_x_monotonicity_margin"]
        ),
        "global_n_tail_extension_certified": n_tail["global_extension_certified"],
        "global_tail_argument": (
            "For 0<=x<=1 use exp(a*x)<=exp(a) and exp(-pi*n^2*exp(4x))<=exp(-pi*n^2). "
            "For x>=1, every value monomial exp(a*x-pi*n^2*exp(4x)) and derivative-core "
            "monomial x*exp(a*x-pi*n^2*exp(4x)) is decreasing for n>=31 by the recorded margins."
        ),
        "finite_c0_ball": arb_text(c0_finite, 80),
        "finite_c0_lower": arb_lower_text(c0_finite_lower),
        "finite_c0_upper": arb_upper_text(c0_finite_upper),
        "full_c0_lower": arb_lower_text(c0_full_lower),
        "c0_lower_target": DEFAULT_C0_LOWER_TARGET,
        "finite_and_full_c0_lowers_certified": c0_lowers_certified,
        "compact_phi_abs_integral_bound_upper": arb_upper_text(
            compact_phi_abs_integral_bound
        ),
        "far_phi_abs_integral_bound_upper": arb_upper_text(far_phi_abs_integral_bound),
        "finite_phi_abs_integral_bound_upper": arb_upper_text(
            finite_phi_abs_integral_bound
        ),
        "compact_x_phip_over_2_abs_integral_bound_upper": arb_upper_text(
            compact_x_phip_over_2_abs_integral_bound
        ),
        "far_x_phip_over_2_abs_integral_bound_upper": arb_upper_text(
            far_x_phip_over_2_abs_integral_bound
        ),
        "finite_x_phip_over_2_abs_integral_bound_upper": arb_upper_text(
            finite_x_phip_over_2_abs_integral_bound
        ),
        "value_normalized_n_tail_correction_upper": arb_upper_text(
            value_normalized_correction
        ),
        "derivative_normalized_n_tail_correction_upper": arb_upper_text(
            derivative_normalized_correction
        ),
        "normalization_correction_formula": (
            "|Phi/Phi0-Phi30/Phi30_0| <= |Phi-Phi30|/Phi0_lower "
            "+ |Phi30|*|Phi0-Phi30_0|/(Phi30_0_lower*Phi0_lower), integrated channelwise."
        ),
        "source_compact_value_ball": compact_summary["value_compact_total_ball"],
        "source_compact_derivative_ball": compact_summary["derivative_compact_total_ball"],
        "source_far_value_bound": far_summary["value_total_tail_bound"],
        "source_far_derivative_bound": far_summary["derivative_total_tail_bound"],
        "finite_entire_value_ball": arb_text(finite_entire_value, 80),
        "finite_entire_derivative_ball": arb_text(finite_entire_derivative, 80),
        "full_value_expectation_ball": arb_text(full_value, 80),
        "full_derivative_expectation_ball": arb_text(full_derivative, 80),
        "full_value_certified_negative": bool(full_value < 0),
        "full_derivative_certified_negative": bool(full_derivative < 0),
        "quadrature_ratio_radius_cap": far_summary["quadrature_ratio_radius_cap"],
        "value_first_omitted_expectation_ball": arb_text(value_first_omitted, 70),
        "derivative_first_omitted_expectation_ball": arb_text(
            derivative_first_omitted, 70
        ),
        "value_full_ratio_to_first_omitted_upper": arb_upper_text(value_ratio),
        "derivative_full_ratio_to_first_omitted_upper": arb_upper_text(derivative_ratio),
        "value_remaining_margin_below_one_lower": arb_lower_text(value_margin),
        "derivative_remaining_margin_below_one_lower": arb_lower_text(derivative_margin),
        "both_full_expectation_ratios_below_one": both_below_one,
        "complete_worst_row_expectation_certified": True,
        "quadrature_needed_for_worst_row_expectation": False,
        "remaining_worst_row_integral_sources": [],
        "source_minimum_value_denominator_lower": first_omitted["summary"][
            "minimum_value_denominator_lower"
        ],
        "source_minimum_derivative_denominator_lower": first_omitted["summary"][
            "minimum_derivative_denominator_lower"
        ],
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrfec_01_compact_source_import",
            role="interval_certificate_import",
            readiness="available_interval_certificate",
            claim=(
                "Import the degree-128 exact-moment/Cauchy certificate for the finite n<=30 core on "
                "0<=y<=200 at T=10000, F_21."
            ),
            diagnostics={
                "source_compact_value_ball": diagnostics["source_compact_value_ball"],
                "source_compact_derivative_ball": diagnostics[
                    "source_compact_derivative_ball"
                ],
            },
            source_artifacts=[paths["compact_json"], paths["compact_note"]],
            proof_boundary="Finite n<=30 compact source import for one worst row only.",
        ),
        MatrixRow(
            id="nlrgwrfec_02_far_tail_source_import",
            role="tail_certificate_import",
            readiness="available_tail_certificate",
            claim=(
                "Import the Arb finite n<=30 cancellation-reduced far-tail bounds on y>=200 for the same row."
            ),
            diagnostics={
                "source_far_value_bound": diagnostics["source_far_value_bound"],
                "source_far_derivative_bound": diagnostics["source_far_derivative_bound"],
            },
            source_artifacts=[paths["far_tail_json"], paths["far_tail_note"]],
            proof_boundary="Finite n<=30 far-tail source import for one row only.",
        ),
        MatrixRow(
            id="nlrgwrfec_03_global_n_tail_certificate",
            role="analytic_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "Arb geometric majorants and x>=1 monotonicity extend the n>=31 Phi and x*Phi' tail "
                "bounds from 0<=x<=1 to every x>=0."
            ),
            diagnostics={
                "tail_start_n": diagnostics["tail_start_n"],
                "global_tail_argument": diagnostics["global_tail_argument"],
                "value_n_tail_global_bound_upper": diagnostics[
                    "value_n_tail_global_bound_upper"
                ],
                "derivative_n_tail_global_bound_upper": diagnostics[
                    "derivative_n_tail_global_bound_upper"
                ],
                "c0_n_tail_bound_upper": diagnostics["c0_n_tail_bound_upper"],
                "value_x_ge_1_monotonicity_margin_lower": diagnostics[
                    "value_x_ge_1_monotonicity_margin_lower"
                ],
                "derivative_x_ge_1_monotonicity_margin_lower": diagnostics[
                    "derivative_x_ge_1_monotonicity_margin_lower"
                ],
                "global_n_tail_extension_certified": diagnostics[
                    "global_n_tail_extension_certified"
                ],
            },
            source_artifacts=[paths["phi_tail_grid_note"], paths["node_c0_note"]],
            proof_boundary="Global n>=31 analytic tail certificate only; not all-row coverage.",
        ),
        MatrixRow(
            id="nlrgwrfec_04_normalization_correction",
            role="normalization_certificate",
            readiness="available_normalization_certificate",
            claim=(
                "The global n>=31 numerator bounds and Phi(0) tail are propagated through the exact "
                "finite/full normalization identity in both value and derivative channels."
            ),
            diagnostics={
                "normalization_correction_formula": diagnostics[
                    "normalization_correction_formula"
                ],
                "finite_c0_lower": diagnostics["finite_c0_lower"],
                "full_c0_lower": diagnostics["full_c0_lower"],
                "finite_and_full_c0_lowers_certified": diagnostics[
                    "finite_and_full_c0_lowers_certified"
                ],
                "value_normalized_n_tail_correction_upper": diagnostics[
                    "value_normalized_n_tail_correction_upper"
                ],
                "derivative_normalized_n_tail_correction_upper": diagnostics[
                    "derivative_normalized_n_tail_correction_upper"
                ],
            },
            source_artifacts=[paths["compact_note"], paths["phi_tail_grid_note"], paths["node_c0_note"]],
            proof_boundary="Worst-row normalization correction only; not continuum-in-T or all-index coverage.",
        ),
        MatrixRow(
            id="nlrgwrfec_05_complete_expectation_composition",
            role="full_expectation_certificate",
            readiness="available_worst_row_expectation_certificate",
            claim=(
                "Composing compact finite core, finite-core far tail, and the normalized global n>=31 "
                "correction certifies the complete true relative-Gaussian expectation for the worst row."
            ),
            diagnostics={
                "finite_entire_value_ball": diagnostics["finite_entire_value_ball"],
                "finite_entire_derivative_ball": diagnostics[
                    "finite_entire_derivative_ball"
                ],
                "full_value_expectation_ball": diagnostics["full_value_expectation_ball"],
                "full_derivative_expectation_ball": diagnostics[
                    "full_derivative_expectation_ball"
                ],
                "full_value_certified_negative": diagnostics["full_value_certified_negative"],
                "full_derivative_certified_negative": diagnostics[
                    "full_derivative_certified_negative"
                ],
                "complete_worst_row_expectation_certified": diagnostics[
                    "complete_worst_row_expectation_certified"
                ],
                "quadrature_needed_for_worst_row_expectation": diagnostics[
                    "quadrature_needed_for_worst_row_expectation"
                ],
            },
            source_artifacts=[
                paths["compact_note"],
                paths["far_tail_note"],
                paths["coefficient_core_note"],
            ],
            proof_boundary="Complete true expectation for T=10000, F_21 only; not the other grid rows or collar.",
        ),
        MatrixRow(
            id="nlrgwrfec_06_below_first_omitted_certificate",
            role="ratio_certificate",
            readiness="available_worst_row_ratio_certificate",
            claim=(
                "Both complete true worst-row expectation magnitudes are strictly below their corresponding "
                "first omitted terms, with certified positive margins."
            ),
            diagnostics={
                "value_first_omitted_expectation_ball": diagnostics[
                    "value_first_omitted_expectation_ball"
                ],
                "derivative_first_omitted_expectation_ball": diagnostics[
                    "derivative_first_omitted_expectation_ball"
                ],
                "value_full_ratio_to_first_omitted_upper": diagnostics[
                    "value_full_ratio_to_first_omitted_upper"
                ],
                "derivative_full_ratio_to_first_omitted_upper": diagnostics[
                    "derivative_full_ratio_to_first_omitted_upper"
                ],
                "value_remaining_margin_below_one_lower": diagnostics[
                    "value_remaining_margin_below_one_lower"
                ],
                "derivative_remaining_margin_below_one_lower": diagnostics[
                    "derivative_remaining_margin_below_one_lower"
                ],
                "both_full_expectation_ratios_below_one": diagnostics[
                    "both_full_expectation_ratios_below_one"
                ],
            },
            source_artifacts=[paths["first_omitted_json"], paths["first_omitted_note"]],
            proof_boundary="One-row first-omitted ratio certificate only; not a uniform residual theorem.",
        ),
        MatrixRow(
            id="nlrgwrfec_07_all_row_handoff",
            role="intervalization_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "The complete worst-row expectation source is retired; the next finite task is to apply the "
                "same direct certificate to every remaining recorded T/index row before a collar bridge."
            ),
            diagnostics={
                "retired_component": "complete true relative-Gaussian expectation at T=10000, F_21",
                "remaining_worst_row_integral_sources": diagnostics[
                    "remaining_worst_row_integral_sources"
                ],
                "remaining_requirements": [
                    "direct compact/tail certificates for the other recorded T/index rows",
                    "all-source aggregation and rounding accounting on the finite grid",
                    "finite-grid to full-collar coverage",
                    "uniform scaled-curvature and cone-entry handoff",
                ],
                "target_closing": diagnostics["target_closing"],
            },
            source_artifacts=[paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Handoff only; not an all-row finite-grid certificate or uniform collar theorem.",
        ),
        MatrixRow(
            id="nlrgwrfec_08_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "This one-row expectation certificate cannot be promoted to the finite grid, uniform collar, "
                "scaled-curvature theorem, cone entry, RH, or Lambda <= 0."
            ),
            diagnostics={
                "forbidden_promotions": [
                    "all-row finite-grid interval certificate",
                    "uniform collar theorem",
                    "scaled-curvature monotonicity",
                    "cone entry",
                    "RH",
                    "Lambda <= 0",
                ]
            },
            source_artifacts=[paths["dependency_graph"], paths["formal_tail_obstruction"]],
            proof_boundary="Proof-hygiene gate only; not RH and not Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(compact_path: Path, far_tail_path: Path, first_omitted_path: Path) -> dict:
    compact = load_json(compact_path)
    far_tail = load_json(far_tail_path)
    first_omitted = load_json(first_omitted_path)
    paths = source_paths(compact_path, far_tail_path, first_omitted_path)
    diagnostics = build_diagnostics(compact, far_tail, first_omitted)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "precision_bits": diagnostics["precision_bits"],
        "tail_start_n": diagnostics["tail_start_n"],
        "global_n_tail_extension_certified": diagnostics[
            "global_n_tail_extension_certified"
        ],
        "finite_and_full_c0_lowers_certified": diagnostics[
            "finite_and_full_c0_lowers_certified"
        ],
        "value_n_tail_global_bound_upper": diagnostics[
            "value_n_tail_global_bound_upper"
        ],
        "derivative_n_tail_global_bound_upper": diagnostics[
            "derivative_n_tail_global_bound_upper"
        ],
        "c0_n_tail_bound_upper": diagnostics["c0_n_tail_bound_upper"],
        "value_normalized_n_tail_correction_upper": diagnostics[
            "value_normalized_n_tail_correction_upper"
        ],
        "derivative_normalized_n_tail_correction_upper": diagnostics[
            "derivative_normalized_n_tail_correction_upper"
        ],
        "finite_entire_value_ball": diagnostics["finite_entire_value_ball"],
        "finite_entire_derivative_ball": diagnostics["finite_entire_derivative_ball"],
        "full_value_expectation_ball": diagnostics["full_value_expectation_ball"],
        "full_derivative_expectation_ball": diagnostics[
            "full_derivative_expectation_ball"
        ],
        "full_value_certified_negative": diagnostics["full_value_certified_negative"],
        "full_derivative_certified_negative": diagnostics[
            "full_derivative_certified_negative"
        ],
        "value_full_ratio_to_first_omitted_upper": diagnostics[
            "value_full_ratio_to_first_omitted_upper"
        ],
        "derivative_full_ratio_to_first_omitted_upper": diagnostics[
            "derivative_full_ratio_to_first_omitted_upper"
        ],
        "value_remaining_margin_below_one_lower": diagnostics[
            "value_remaining_margin_below_one_lower"
        ],
        "derivative_remaining_margin_below_one_lower": diagnostics[
            "derivative_remaining_margin_below_one_lower"
        ],
        "both_full_expectation_ratios_below_one": diagnostics[
            "both_full_expectation_ratios_below_one"
        ],
        "composed_source_count": 3,
        "below_one_ratio_count": 2,
        "complete_worst_row_expectation_certified": True,
        "quadrature_needed_for_worst_row_expectation": False,
        "remaining_worst_row_integral_source_count": 0,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst-row full-expectation certificate composes three rigorous sources: the degree-128 "
            "finite-core compact x-moment integral, the finite-core y>=200 tail, and a new Arb global n>=31 "
            "normalization correction. For T=10000, F_21, the complete true value and derivative expectations "
            "are certified negative, with first-omitted ratio uppers below 0.971 and positive margins above "
            "0.029. This removes generalized Gauss-Laguerre quadrature from the worst-row expectation proof. "
            "It remains a one-row certificate: the other finite-grid rows, all-source aggregation, and the "
            "finite-grid-to-collar theorem remain open."
        ),
    }
    return {
        "kind": (
            "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate"
        ),
        "date": "2026-07-09",
        "status": "worst-row full-expectation certificate",
        "source_compact_x_moment_taylor_certificate": paths["compact_note"],
        "source_compact_x_moment_taylor_json": paths["compact_json"],
        "source_endpoint_x_moment_taylor_certificate": paths["endpoint_note"],
        "source_far_tail_split_certificate": paths["far_tail_note"],
        "source_far_tail_split_json": paths["far_tail_json"],
        "source_first_omitted_denominator_certificate": paths["first_omitted_note"],
        "source_first_omitted_denominator_json": paths["first_omitted_json"],
        "source_coefficient_core_certificate": paths["coefficient_core_note"],
        "source_phi_tail_grid_certificate": paths["phi_tail_grid_note"],
        "source_node_c0_range_certificate": paths["node_c0_note"],
        "source_finite_plus_tail_budget_certificate": paths["finite_plus_tail_note"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py"
        ),
        "proof_boundary": (
            "Complete true relative-Gaussian expectation certificate for the single T=10000, F_21 worst row. "
            "It composes the finite n<=30 compact and far-tail integrals with a global n>=31 Phi/Phi'/Phi(0) "
            "normalization correction and certifies both first-omitted ratios below one. It does not certify "
            "the other finite-grid rows, does not aggregate every intervalization source over the grid, does "
            "not prove finite-grid-to-collar coverage, does not prove scaled-curvature monotonicity or cone "
            "entry, and does not prove RH or Lambda <= 0."
        ),
        "diagnostics": diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "The complete expectation certificate applies only to T=10000, F_21.",
            "The finite compact core, finite far tail, and global n>=31 normalization correction are all included.",
            "The global n>=31 extension uses explicit x>=1 monotonicity margins, not a finite-node assumption.",
            "Both full expectation balls are negative and strictly below one first omitted term in magnitude.",
            "Generalized Gauss-Laguerre quadrature is not used to certify the final worst-row expectation.",
            "No claim about the other finite-grid rows or continuum collar is made.",
            "No row is ready_to_apply for the full intervalization or uniform-collar target.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation "
        f"certificate: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['composed_source_count']} composed sources, "
        f"{summary['below_one_ratio_count']} below-one full ratios, "
        f"{summary['remaining_worst_row_integral_source_count']} open worst-row integral sources, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Full-Expectation Certificate",
        "",
        "Date: 2026-07-09",
        "",
        "Status: worst-row full-expectation certificate. This is not a proof",
        "of an all-row finite-grid certificate, a uniform collar theorem, RH,",
        "or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate`.",
        "",
        "Proof boundary: this artifact certifies the complete true relative-Gaussian",
        "value and derivative expectations only for `T=10000`, `F_21`. It does not",
        "certify the other finite-grid rows or bridge the finite grid to the collar.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Global N-Tail",
        "",
        "```text",
        diagnostics["global_tail_argument"],
        f"tail start n: {summary['tail_start_n']}",
        f"global extension certified: {summary['global_n_tail_extension_certified']}",
        f"value x>=1 monotonicity margin lower: {diagnostics['value_x_ge_1_monotonicity_margin_lower']}",
        f"derivative x>=1 monotonicity margin lower: {diagnostics['derivative_x_ge_1_monotonicity_margin_lower']}",
        f"value n-tail bound upper: {summary['value_n_tail_global_bound_upper']}",
        f"derivative n-tail bound upper: {summary['derivative_n_tail_global_bound_upper']}",
        f"Phi(0) n-tail bound upper: {summary['c0_n_tail_bound_upper']}",
        "```",
        "",
        "## Normalization",
        "",
        "```text",
        diagnostics["normalization_correction_formula"],
        f"finite Phi_30(0) lower: {diagnostics['finite_c0_lower']}",
        f"full Phi(0) lower: {diagnostics['full_c0_lower']}",
        f"both c0 lower bounds certified: {summary['finite_and_full_c0_lowers_certified']}",
        f"value normalized n-tail correction upper: {summary['value_normalized_n_tail_correction_upper']}",
        f"derivative normalized n-tail correction upper: {summary['derivative_normalized_n_tail_correction_upper']}",
        "```",
        "",
        "## Complete Expectations",
        "",
        "```text",
        f"finite entire value ball: {summary['finite_entire_value_ball']}",
        f"finite entire derivative ball: {summary['finite_entire_derivative_ball']}",
        f"full value expectation ball: {summary['full_value_expectation_ball']}",
        f"full derivative expectation ball: {summary['full_derivative_expectation_ball']}",
        f"full value certified negative: {summary['full_value_certified_negative']}",
        f"full derivative certified negative: {summary['full_derivative_certified_negative']}",
        f"value ratio / first omitted upper: {summary['value_full_ratio_to_first_omitted_upper']}",
        f"derivative ratio / first omitted upper: {summary['derivative_full_ratio_to_first_omitted_upper']}",
        f"value margin below one lower: {summary['value_remaining_margin_below_one_lower']}",
        f"derivative margin below one lower: {summary['derivative_remaining_margin_below_one_lower']}",
        f"both full ratios below one: {summary['both_full_expectation_ratios_below_one']}",
        f"complete worst-row expectation certified: {summary['complete_worst_row_expectation_certified']}",
        f"quadrature needed for worst-row expectation: {summary['quadrature_needed_for_worst_row_expectation']}",
        "```",
        "",
        "## Remaining Work",
        "",
        "```text",
        "Apply the same direct compact/global-tail certificate to the other recorded T/index rows.",
        "Aggregate every finite-grid source and its rounding allowance.",
        "Prove the finite-grid-to-full-collar and scaled-curvature handoff.",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_compact_x_moment_taylor_certificate"],
        artifact["source_compact_x_moment_taylor_json"],
        artifact["source_far_tail_split_certificate"],
        artifact["source_far_tail_split_json"],
        artifact["source_first_omitted_denominator_certificate"],
        artifact["source_first_omitted_denominator_json"],
        artifact["source_coefficient_core_certificate"],
        artifact["source_phi_tail_grid_certificate"],
        artifact["source_node_c0_range_certificate"],
        artifact["source_finite_plus_tail_budget_certificate"],
        artifact["source_quadrature_remainder_route_matrix"],
        artifact["source_intervalization_target"],
        artifact["source_uniform_remainder_target"],
        artifact["source_dependency_graph"],
        artifact["source_formal_tail_obstruction"],
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
    parser.add_argument("--compact-json", type=Path, default=DEFAULT_COMPACT_JSON)
    parser.add_argument("--far-tail-json", type=Path, default=DEFAULT_FAR_TAIL_JSON)
    parser.add_argument("--first-omitted-json", type=Path, default=DEFAULT_FIRST_OMITTED_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    compact_path = args.compact_json if args.compact_json.is_absolute() else REPO_ROOT / args.compact_json
    far_tail_path = args.far_tail_json if args.far_tail_json.is_absolute() else REPO_ROOT / args.far_tail_json
    first_omitted_path = (
        args.first_omitted_json
        if args.first_omitted_json.is_absolute()
        else REPO_ROOT / args.first_omitted_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(compact_path, far_tail_path, first_omitted_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row full-expectation "
        f"certificate: {out_json.relative_to(REPO_ROOT).as_posix()} and "
        f"{note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
