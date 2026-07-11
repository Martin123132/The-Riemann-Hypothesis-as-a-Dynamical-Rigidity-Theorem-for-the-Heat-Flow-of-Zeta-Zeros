#!/usr/bin/env python3
"""Build a rigorous compact x-moment Taylor certificate for the worst row."""

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

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)
from jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix import (  # noqa: E402
    finite_phi_series,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ENDPOINT_CERTIFICATE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.json"
)
DEFAULT_PHI_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json"
)
DEFAULT_FAR_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_Y_RIGHT = 200
DEFAULT_DISK_RADIUS = "0.38"
DEFAULT_TAYLOR_DEGREE = 128
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PRECISION_BITS = 8192
SELECTED_DEGREES = {0, 1, 2, 20, 40, 41, 42, 64, 80, 96, 112, 128}


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


def arb_text(value: flint.arb, digits: int = 60) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    return value.lower().str(digits, radius=False).replace("e", "E")


def arb_upper_text(value: flint.arb, digits: int = 60) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def abs_upper(value: flint.arb) -> flint.arb:
    lower_abs = abs(flint.arb(value.lower()))
    upper_abs = abs(flint.arb(value.upper()))
    return lower_abs if lower_abs > upper_abs else upper_abs


def source_paths(
    endpoint_certificate_path: Path,
    phi_tail_path: Path,
    far_tail_path: Path,
) -> dict[str, str]:
    return {
        "endpoint_certificate_json": endpoint_certificate_path.relative_to(REPO_ROOT).as_posix(),
        "endpoint_certificate_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.md"
        ),
        "endpoint_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.md"
        ),
        "interpolation_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md"
        ),
        "arb_interpolant_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md"
        ),
        "compact_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md"
        ),
        "phi_tail_json": phi_tail_path.relative_to(REPO_ROOT).as_posix(),
        "phi_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md"
        ),
        "far_tail_json": far_tail_path.relative_to(REPO_ROOT).as_posix(),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
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


def configure_precision() -> None:
    flint.ctx.prec = DEFAULT_PRECISION_BITS


def compact_gamma_moment(degree: int) -> flint.arb:
    """Integral of x^degree over 0<=y<=200 after y=T*x^2."""
    alpha = flint.arb(DEFAULT_ALPHA)
    shape = alpha + 1 + flint.arb(degree) / 2
    scale = flint.arb(DEFAULT_T) ** (flint.arb(degree) / 2)
    return flint.arb(DEFAULT_Y_RIGHT).gamma_lower(shape) / ((alpha + 1).gamma() * scale)


def complex_disk_majorant() -> dict[str, flint.arb | bool]:
    radius = flint.arb(DEFAULT_DISK_RADIUS)
    pi = flint.arb.pi()
    radius_admissible = bool(4 * radius < pi / 2)
    real_exp_lower = (-4 * radius).exp() * (4 * radius).cos()
    if not radius_admissible or real_exp_lower.contains(0) or not bool(real_exp_lower > 0):
        raise RuntimeError("compact complex-disk damping lower bound was not certified positive")
    majorant = flint.arb(0)
    for n in range(1, DEFAULT_PHI_TERM_COUNT + 1):
        n2 = flint.arb(n * n)
        prefactor = (
            2 * pi * pi * n2 * n2 * (9 * radius).exp()
            + 3 * pi * n2 * (5 * radius).exp()
        )
        majorant += prefactor * (-pi * n2 * real_exp_lower).exp()
    return {
        "radius": radius,
        "radius_admissible": radius_admissible,
        "real_exp_lower": real_exp_lower,
        "majorant": majorant,
    }


def build_moment_rows(
    coefficients: list[flint.arb],
    ratios: list[flint.arb],
    c0: flint.arb,
) -> tuple[list[dict], flint.arb, flint.arb]:
    rows: list[dict] = []
    value_integral = flint.arb(0)
    derivative_integral = flint.arb(0)
    for degree, phi_coefficient in enumerate(coefficients):
        value_coefficient = phi_coefficient / c0
        derivative_coefficient = flint.arb(degree) * phi_coefficient / (2 * c0)
        subtracted_ratio = None
        if degree % 2 == 0 and degree // 2 <= DEFAULT_POLYNOMIAL_M:
            ratio_index = degree // 2
            subtracted_ratio = ratios[ratio_index]
            value_coefficient -= subtracted_ratio
            if degree:
                derivative_coefficient -= flint.arb(ratio_index) * subtracted_ratio
        moment = compact_gamma_moment(degree)
        value_contribution = value_coefficient * moment
        derivative_contribution = derivative_coefficient * moment
        value_integral += value_contribution
        derivative_integral += derivative_contribution
        rows.append(
            {
                "degree": degree,
                "phi_taylor_coefficient_ball": arb_text(phi_coefficient),
                "subtracted_even_ratio_ball": (
                    None if subtracted_ratio is None else arb_text(subtracted_ratio)
                ),
                "value_core_coefficient_ball": arb_text(value_coefficient),
                "derivative_core_coefficient_ball": arb_text(derivative_coefficient),
                "compact_transformed_gamma_moment_ball": arb_text(moment),
                "value_moment_contribution_ball": arb_text(value_contribution),
                "derivative_moment_contribution_ball": arb_text(derivative_contribution),
                "proof_boundary": (
                    "Arb coefficient and exact 0<=y<=200 transformed-moment row for the finite n<=30 core only."
                ),
            }
        )
    return rows, value_integral, derivative_integral


def build_diagnostics(endpoint_certificate: dict, phi_tail: dict, far_tail: dict) -> dict:
    configure_precision()
    coefficients = finite_phi_series(DEFAULT_PHI_TERM_COUNT, DEFAULT_TAYLOR_DEGREE)
    _ratio_rows, ratios = build_ratio_rows(2 * DEFAULT_POLYNOMIAL_M, DEFAULT_RATIO_CUTOFF_N)
    c0 = coefficients[0]
    if c0.contains(0) or not bool(c0 > 0):
        raise RuntimeError("finite Phi(0) did not certify positive")
    moment_rows, value_polynomial_integral, derivative_polynomial_integral = build_moment_rows(
        coefficients, ratios, c0
    )
    disk = complex_disk_majorant()
    radius = disk["radius"]
    majorant = disk["majorant"]
    assert isinstance(radius, flint.arb)
    assert isinstance(majorant, flint.arb)
    x_right = (flint.arb(DEFAULT_Y_RIGHT) / flint.arb(DEFAULT_T)).sqrt()
    q = x_right / radius
    if not bool(q < 1):
        raise RuntimeError("compact x interval is not inside the Cauchy disk")
    c0_lower = flint.arb(c0.lower())
    first_tail_degree = DEFAULT_TAYLOR_DEGREE + 1
    value_tail_sup = majorant * (q**first_tail_degree) / (c0_lower * (1 - q))
    derivative_tail_sup = (
        majorant
        * (q**first_tail_degree)
        * (flint.arb(first_tail_degree) - flint.arb(DEFAULT_TAYLOR_DEGREE) * q)
        / (2 * c0_lower * ((1 - q) ** 2))
    )
    compact_mass = compact_gamma_moment(0)
    value_tail_integral_radius = abs_upper(value_tail_sup * compact_mass)
    derivative_tail_integral_radius = abs_upper(derivative_tail_sup * compact_mass)
    value_total = value_polynomial_integral + flint.arb(0, value_tail_integral_radius)
    derivative_total = derivative_polynomial_integral + flint.arb(0, derivative_tail_integral_radius)

    endpoint_summary = endpoint_certificate["summary"]
    value_cap = flint.arb(endpoint_summary["value_unscaled_expectation_error_cap"])
    derivative_cap = flint.arb(endpoint_summary["derivative_unscaled_expectation_error_cap"])
    value_tail_ratio = value_tail_integral_radius / value_cap
    derivative_tail_ratio = derivative_tail_integral_radius / derivative_cap
    both_tail_radii_below_caps = bool(value_tail_ratio < 1 and derivative_tail_ratio < 1)
    if not both_tail_radii_below_caps:
        raise RuntimeError("compact Taylor remainder radii did not fit the source caps")
    if value_total.contains(0) or derivative_total.contains(0):
        raise RuntimeError("compact finite-core total balls did not separate zero")

    selected_rows = [row for row in moment_rows if row["degree"] in SELECTED_DEGREES]
    far_summary = far_tail["summary"]
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "alpha": DEFAULT_ALPHA,
        "y_interval": "0<=y<=200",
        "x_interval": "0<=x<=sqrt(0.02)",
        "x_right_endpoint": arb_upper_text(x_right, 70),
        "precision_bits": DEFAULT_PRECISION_BITS,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
        "polynomial_M": DEFAULT_POLYNOMIAL_M,
        "taylor_degree": DEFAULT_TAYLOR_DEGREE,
        "first_tail_degree": first_tail_degree,
        "finite_c0_ball": arb_text(c0, 80),
        "finite_c0_lower": arb_lower_text(c0, 70),
        "disk_radius": DEFAULT_DISK_RADIUS,
        "four_disk_radius_upper": arb_upper_text(4 * radius, 40),
        "pi_over_two_lower": arb_lower_text(flint.arb.pi() / 2, 40),
        "disk_radius_admissible": disk["radius_admissible"],
        "real_exp4_lower_on_disk": arb_lower_text(disk["real_exp_lower"], 70),
        "finite_phi_complex_disk_majorant_upper": arb_upper_text(majorant, 70),
        "x_to_disk_radius_ratio_upper": arb_upper_text(q, 70),
        "compact_gamma_mass_ball": arb_text(compact_mass, 70),
        "compact_gamma_mass_upper": arb_upper_text(compact_mass, 70),
        "compact_moment_formula": (
            "mu_k(200)=T^(-k/2)*lower_gamma(alpha+1+k/2,200)/Gamma(alpha+1)"
        ),
        "complex_disk_bound_formula": (
            "For |z|<=R<pi/8, Re(exp(4z))>=exp(-4R)*cos(4R)>0; hence "
            "|Phi_30(z)| is bounded by the recorded positive termwise majorant."
        ),
        "value_cauchy_tail_formula": "M*q^(D+1)/(c0_lower*(1-q))",
        "derivative_cauchy_tail_formula": (
            "M*q^(D+1)*((D+1)-D*q)/(2*c0_lower*(1-q)^2)"
        ),
        "value_cauchy_tail_sup_upper": arb_upper_text(value_tail_sup, 70),
        "derivative_cauchy_tail_sup_upper": arb_upper_text(derivative_tail_sup, 70),
        "value_cauchy_tail_integral_radius_upper": arb_upper_text(value_tail_integral_radius, 70),
        "derivative_cauchy_tail_integral_radius_upper": arb_upper_text(
            derivative_tail_integral_radius, 70
        ),
        "compact_taylor_moment_rows": moment_rows,
        "compact_taylor_moment_row_count": len(moment_rows),
        "selected_compact_taylor_moment_rows": selected_rows,
        "value_taylor_moment_integral_ball": arb_text(value_polynomial_integral, 80),
        "derivative_taylor_moment_integral_ball": arb_text(
            derivative_polynomial_integral, 80
        ),
        "value_compact_total_ball": arb_text(value_total, 80),
        "derivative_compact_total_ball": arb_text(derivative_total, 80),
        "value_compact_certified_negative": bool(value_total < 0),
        "derivative_compact_certified_negative": bool(derivative_total < 0),
        "value_unscaled_expectation_error_cap": endpoint_summary[
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": endpoint_summary[
            "derivative_unscaled_expectation_error_cap"
        ],
        "value_tail_radius_to_cap_ratio_upper": arb_upper_text(value_tail_ratio, 70),
        "derivative_tail_radius_to_cap_ratio_upper": arb_upper_text(derivative_tail_ratio, 70),
        "both_tail_radii_below_caps": both_tail_radii_below_caps,
        "finite_core_compact_interval_certified": True,
        "panel_interpolation_needed_for_finite_core": False,
        "source_endpoint_first_panel_certified": endpoint_summary[
            "finite_core_first_panel_certified"
        ],
        "source_phi_tail_certified": phi_tail["summary"]["finite_grid_tail_source_certified"],
        "source_far_tail_value_ratio": far_summary[
            "value_tail_to_unscaled_cap_ratio_upper"
        ],
        "source_far_tail_derivative_ratio": far_summary[
            "derivative_tail_to_unscaled_cap_ratio_upper"
        ],
        "full_n_tail_composed_here": False,
        "far_y_tail_composed_here": False,
        "remaining_compact_panels": [],
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrcxmtc_01_compact_x_reduction",
            role="scope_reduction",
            readiness="available_scope_reduction",
            claim=(
                "For T=10000, y=T*x^2 maps the complete compact interval 0<=y<=200 into "
                "0<=x<=sqrt(0.02), entirely inside the certified complex disk |z|<=0.38."
            ),
            diagnostics={
                "T": diagnostics["T"],
                "index": diagnostics["index"],
                "y_interval": diagnostics["y_interval"],
                "x_interval": diagnostics["x_interval"],
                "x_right_endpoint": diagnostics["x_right_endpoint"],
                "disk_radius": diagnostics["disk_radius"],
            },
            source_artifacts=[paths["endpoint_certificate_json"], paths["endpoint_certificate_note"]],
            proof_boundary="Worst-row compact finite-core change of variables only; not all rows or the full collar.",
        ),
        MatrixRow(
            id="nlrgwrcxmtc_02_exact_compact_moments",
            role="exact_moment_reduction",
            readiness="available_exact_reduction",
            claim=(
                "Each x^k term is integrated over 0<=y<=200 by the exact lower-incomplete-Gamma moment "
                "identity, evaluated with Arb for all degrees 0 through 128."
            ),
            diagnostics={
                "compact_moment_formula": diagnostics["compact_moment_formula"],
                "compact_taylor_moment_rows": diagnostics["compact_taylor_moment_row_count"],
                "compact_gamma_mass_ball": diagnostics["compact_gamma_mass_ball"],
                "selected_compact_taylor_moment_rows": diagnostics[
                    "selected_compact_taylor_moment_rows"
                ],
            },
            source_artifacts=[paths["endpoint_certificate_note"], paths["finite_part_note"]],
            proof_boundary="Exact 0<=y<=200 transformed moments for this worst-row finite core only.",
        ),
        MatrixRow(
            id="nlrgwrcxmtc_03_compact_complex_disk_majorant",
            role="analytic_remainder_certificate",
            readiness="available_analytic_certificate",
            claim=(
                "The finite Phi_30 kernel has the recorded rigorous majorant on |z|<=0.38, with "
                "4R<pi/2 and a positive lower bound for Re(exp(4z))."
            ),
            diagnostics={
                "complex_disk_bound_formula": diagnostics["complex_disk_bound_formula"],
                "disk_radius_admissible": diagnostics["disk_radius_admissible"],
                "four_disk_radius_upper": diagnostics["four_disk_radius_upper"],
                "pi_over_two_lower": diagnostics["pi_over_two_lower"],
                "real_exp4_lower_on_disk": diagnostics["real_exp4_lower_on_disk"],
                "finite_phi_complex_disk_majorant_upper": diagnostics[
                    "finite_phi_complex_disk_majorant_upper"
                ],
                "x_to_disk_radius_ratio_upper": diagnostics[
                    "x_to_disk_radius_ratio_upper"
                ],
            },
            source_artifacts=[paths["endpoint_certificate_note"], paths["interpolation_route_note"]],
            proof_boundary="Finite n<=30 disk majorant only; the infinite n-tail remains a separate source.",
        ),
        MatrixRow(
            id="nlrgwrcxmtc_04_compact_taylor_moment_interval",
            role="arb_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The degree-128 finite-core Taylor polynomial, after subtracting the degree-40 ratio "
                "polynomial, is integrated over the full compact interval by exact transformed moments."
            ),
            diagnostics={
                "precision_bits": diagnostics["precision_bits"],
                "taylor_degree": diagnostics["taylor_degree"],
                "value_taylor_moment_integral_ball": diagnostics[
                    "value_taylor_moment_integral_ball"
                ],
                "derivative_taylor_moment_integral_ball": diagnostics[
                    "derivative_taylor_moment_integral_ball"
                ],
            },
            source_artifacts=[paths["arb_interpolant_note"], paths["compact_scout_note"]],
            proof_boundary="Finite-core compact Taylor-polynomial integral only; Cauchy tails are added separately.",
        ),
        MatrixRow(
            id="nlrgwrcxmtc_05_compact_remainder_certificate",
            role="compact_interval_certificate",
            readiness="available_compact_finite_core_certificate",
            claim=(
                "Adding the post-degree-128 Cauchy tails certifies the true finite n<=30 value and "
                "derivative integrals on all of 0<=y<=200, with both remainder radii below the source caps."
            ),
            diagnostics={
                "value_cauchy_tail_integral_radius_upper": diagnostics[
                    "value_cauchy_tail_integral_radius_upper"
                ],
                "derivative_cauchy_tail_integral_radius_upper": diagnostics[
                    "derivative_cauchy_tail_integral_radius_upper"
                ],
                "value_compact_total_ball": diagnostics["value_compact_total_ball"],
                "derivative_compact_total_ball": diagnostics["derivative_compact_total_ball"],
                "value_tail_radius_to_cap_ratio_upper": diagnostics[
                    "value_tail_radius_to_cap_ratio_upper"
                ],
                "derivative_tail_radius_to_cap_ratio_upper": diagnostics[
                    "derivative_tail_radius_to_cap_ratio_upper"
                ],
                "both_tail_radii_below_caps": diagnostics["both_tail_radii_below_caps"],
                "panel_interpolation_needed_for_finite_core": diagnostics[
                    "panel_interpolation_needed_for_finite_core"
                ],
            },
            source_artifacts=[
                paths["endpoint_certificate_note"],
                paths["interpolation_route_note"],
                paths["quadrature_route_note"],
            ],
            proof_boundary=(
                "Complete 0<=y<=200 certificate for the finite n<=30 worst-row core only; not the n>30 "
                "tail, y>200 tail, all-row grid, or uniform collar."
            ),
        ),
        MatrixRow(
            id="nlrgwrcxmtc_06_tail_splice_handoff",
            role="intervalization_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "The finite-core compact integration source is retired for the worst row; the next splice "
                "must compose the separately certified n>30 Phi tail and y>200 far tail without changing normalization."
            ),
            diagnostics={
                "finite_core_compact_interval_certified": diagnostics[
                    "finite_core_compact_interval_certified"
                ],
                "source_phi_tail_certified": diagnostics["source_phi_tail_certified"],
                "source_far_tail_value_ratio": diagnostics["source_far_tail_value_ratio"],
                "source_far_tail_derivative_ratio": diagnostics[
                    "source_far_tail_derivative_ratio"
                ],
                "full_n_tail_composed_here": diagnostics["full_n_tail_composed_here"],
                "far_y_tail_composed_here": diagnostics["far_y_tail_composed_here"],
                "target_closing": diagnostics["target_closing"],
            },
            source_artifacts=[
                paths["phi_tail_json"],
                paths["phi_tail_note"],
                paths["far_tail_json"],
                paths["far_tail_note"],
                paths["finite_plus_tail_note"],
                paths["intervalization_target"],
            ],
            proof_boundary="Tail-splice handoff only; not the complete worst-row expectation or uniform theorem.",
        ),
        MatrixRow(
            id="nlrgwrcxmtc_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "The compact finite-core certificate cannot be promoted past its one-row scope until tail "
                "composition, all-row coverage, rounding aggregation, and grid-to-collar coverage are proved."
            ),
            diagnostics={
                "forbidden_promotions": [
                    "full n>30 tail composition",
                    "full y>200 far-tail composition",
                    "complete worst-row expectation",
                    "all-row finite-grid certificate",
                    "uniform collar theorem",
                    "scaled-curvature monotonicity",
                    "cone entry",
                    "RH",
                    "Lambda <= 0",
                ]
            },
            source_artifacts=[
                paths["dependency_graph"],
                paths["formal_tail_obstruction"],
                paths["uniform_remainder_target"],
            ],
            proof_boundary="Proof-hygiene gate only; not RH and not Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(
    endpoint_certificate_path: Path,
    phi_tail_path: Path,
    far_tail_path: Path,
) -> dict:
    endpoint_certificate = load_json(endpoint_certificate_path)
    phi_tail = load_json(phi_tail_path)
    far_tail = load_json(far_tail_path)
    paths = source_paths(endpoint_certificate_path, phi_tail_path, far_tail_path)
    diagnostics = build_diagnostics(endpoint_certificate, phi_tail, far_tail)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "alpha": diagnostics["alpha"],
        "y_interval": diagnostics["y_interval"],
        "x_interval": diagnostics["x_interval"],
        "x_right_endpoint": diagnostics["x_right_endpoint"],
        "precision_bits": diagnostics["precision_bits"],
        "phi_term_count": diagnostics["phi_term_count"],
        "polynomial_M": diagnostics["polynomial_M"],
        "taylor_degree": diagnostics["taylor_degree"],
        "first_tail_degree": diagnostics["first_tail_degree"],
        "compact_taylor_moment_row_count": diagnostics[
            "compact_taylor_moment_row_count"
        ],
        "disk_radius": diagnostics["disk_radius"],
        "disk_radius_admissible": diagnostics["disk_radius_admissible"],
        "real_exp4_lower_on_disk": diagnostics["real_exp4_lower_on_disk"],
        "finite_phi_complex_disk_majorant_upper": diagnostics[
            "finite_phi_complex_disk_majorant_upper"
        ],
        "x_to_disk_radius_ratio_upper": diagnostics["x_to_disk_radius_ratio_upper"],
        "compact_gamma_mass_upper": diagnostics["compact_gamma_mass_upper"],
        "value_cauchy_tail_integral_radius_upper": diagnostics[
            "value_cauchy_tail_integral_radius_upper"
        ],
        "derivative_cauchy_tail_integral_radius_upper": diagnostics[
            "derivative_cauchy_tail_integral_radius_upper"
        ],
        "value_taylor_moment_integral_ball": diagnostics[
            "value_taylor_moment_integral_ball"
        ],
        "derivative_taylor_moment_integral_ball": diagnostics[
            "derivative_taylor_moment_integral_ball"
        ],
        "value_compact_total_ball": diagnostics["value_compact_total_ball"],
        "derivative_compact_total_ball": diagnostics["derivative_compact_total_ball"],
        "value_compact_certified_negative": diagnostics[
            "value_compact_certified_negative"
        ],
        "derivative_compact_certified_negative": diagnostics[
            "derivative_compact_certified_negative"
        ],
        "value_unscaled_expectation_error_cap": diagnostics[
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "value_tail_radius_to_cap_ratio_upper": diagnostics[
            "value_tail_radius_to_cap_ratio_upper"
        ],
        "derivative_tail_radius_to_cap_ratio_upper": diagnostics[
            "derivative_tail_radius_to_cap_ratio_upper"
        ],
        "both_tail_radii_below_caps": diagnostics["both_tail_radii_below_caps"],
        "finite_core_compact_interval_certified": True,
        "panel_interpolation_needed_for_finite_core": False,
        "compact_certificate_rows": 1,
        "remaining_compact_panel_count": 0,
        "full_n_tail_composed_here": False,
        "far_y_tail_composed_here": False,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst-row compact x-moment Taylor certificate replaces the six-panel interpolation "
            "obligation for the finite n<=30 core. For T=10000 and F_21, the full range 0<=y<=200 maps "
            "inside |x|<=sqrt(0.02)<0.38; a degree-128 Arb Taylor model is integrated by exact transformed "
            "Gamma moments, and the true-function Cauchy remainder radii consume less than 1.3e-9 of their "
            "source caps. Both compact integrals are certified negative. This closes only the finite-core "
            "compact source for one row; n>30 tail composition, y>200 far-tail composition, all-row coverage, "
            "and the uniform collar theorem remain open."
        ),
    }
    return {
        "kind": (
            "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate"
        ),
        "date": "2026-07-09",
        "status": "worst-row compact x-moment Taylor certificate",
        "source_endpoint_x_moment_taylor_certificate": paths["endpoint_certificate_note"],
        "source_endpoint_x_moment_taylor_json": paths["endpoint_certificate_json"],
        "source_endpoint_x_panel_route_matrix": paths["endpoint_route_note"],
        "source_interpolation_remainder_route_matrix": paths["interpolation_route_note"],
        "source_arb_chebyshev_interpolant_moment_scout": paths["arb_interpolant_note"],
        "source_compact_interval_integration_scout": paths["compact_scout_note"],
        "source_phi_tail_grid_certificate": paths["phi_tail_note"],
        "source_phi_tail_grid_json": paths["phi_tail_json"],
        "source_far_tail_split_certificate": paths["far_tail_note"],
        "source_far_tail_split_json": paths["far_tail_json"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_finite_plus_tail_budget_certificate": paths["finite_plus_tail_note"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.py"
        ),
        "proof_boundary": (
            "Worst-row compact x-moment Taylor certificate for the finite n<=30 core only. It certifies the "
            "T=10000, F_21 value and derivative integrals on 0<=y<=200 by exact transformed Gamma moments "
            "and a complex-disk Cauchy remainder. It does not compose the n>30 Phi tail or y>200 far tail "
            "in this artifact, does not certify all finite-grid rows, does not aggregate all error sources, "
            "does not prove the uniform collar theorem, does not prove scaled-curvature monotonicity or cone "
            "entry, and does not prove RH or Lambda <= 0."
        ),
        "diagnostics": diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "The certificate applies only to the finite n<=30 cancellation-reduced core at T=10000, F_21.",
            "All 129 Taylor monomials are integrated over 0<=y<=200 by exact transformed Gamma moments.",
            "The Cauchy disk satisfies 4R<pi/2 and contains the full compact x interval.",
            "Both compact total balls include the full post-degree-128 true-function Cauchy tail.",
            "Panel interpolation is not required for this finite-core compact certificate.",
            "The separate n>30 Phi tail and y>200 far tail are not silently composed here.",
            "No row is ready_to_apply for the full intervalization or uniform-collar target.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact x-moment Taylor "
        f"certificate: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['compact_taylor_moment_row_count']} exact-moment rows, "
        f"{summary['compact_certificate_rows']} certified compact interval, "
        f"{summary['remaining_compact_panel_count']} open compact panels, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Compact X-Moment Taylor Certificate",
        "",
        "Date: 2026-07-09",
        "",
        "Status: worst-row compact x-moment Taylor certificate. This is not a proof",
        "of a complete worst-row expectation, an all-row finite-grid certificate,",
        "a uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate`.",
        "",
        "Proof boundary: this artifact certifies the finite `n<=30` cancellation-reduced",
        "value and derivative cores on all of `0<=y<=200` for `T=10000`, `F_21`.",
        "The separate `n>30` Phi tail and `y>200` far tail are not composed here.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Certified Setup",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"alpha: {summary['alpha']}",
        f"y interval: {summary['y_interval']}",
        f"x interval: {summary['x_interval']}",
        f"x right endpoint: {summary['x_right_endpoint']}",
        f"disk radius: {summary['disk_radius']}",
        f"Taylor degree: {summary['taylor_degree']}",
        f"exact transformed-moment rows: {summary['compact_taylor_moment_row_count']}",
        f"Arb precision bits: {summary['precision_bits']}",
        "```",
        "",
        "Exact compact moment:",
        "",
        "```text",
        diagnostics["compact_moment_formula"],
        "```",
        "",
        "## Complex-Disk Remainder",
        "",
        "```text",
        diagnostics["complex_disk_bound_formula"],
        f"4R upper: {diagnostics['four_disk_radius_upper']}",
        f"pi/2 lower: {diagnostics['pi_over_two_lower']}",
        f"4R<pi/2 certified: {summary['disk_radius_admissible']}",
        f"Re(exp(4z)) lower: {summary['real_exp4_lower_on_disk']}",
        f"finite Phi_30 disk majorant upper: {summary['finite_phi_complex_disk_majorant_upper']}",
        f"x_max/R upper: {summary['x_to_disk_radius_ratio_upper']}",
        f"value integrated Cauchy radius upper: {summary['value_cauchy_tail_integral_radius_upper']}",
        f"derivative integrated Cauchy radius upper: {summary['derivative_cauchy_tail_integral_radius_upper']}",
        f"value tail-radius/cap ratio upper: {summary['value_tail_radius_to_cap_ratio_upper']}",
        f"derivative tail-radius/cap ratio upper: {summary['derivative_tail_radius_to_cap_ratio_upper']}",
        f"both tail radii below caps: {summary['both_tail_radii_below_caps']}",
        "```",
        "",
        "## Compact Enclosures",
        "",
        "```text",
        f"value Taylor-moment integral: {summary['value_taylor_moment_integral_ball']}",
        f"derivative Taylor-moment integral: {summary['derivative_taylor_moment_integral_ball']}",
        f"value compact total ball: {summary['value_compact_total_ball']}",
        f"derivative compact total ball: {summary['derivative_compact_total_ball']}",
        f"value certified negative: {summary['value_compact_certified_negative']}",
        f"derivative certified negative: {summary['derivative_compact_certified_negative']}",
        f"finite-core compact interval certified: {summary['finite_core_compact_interval_certified']}",
        f"panel interpolation needed for finite core: {summary['panel_interpolation_needed_for_finite_core']}",
        "```",
        "",
        "Selected coefficient/moment rows:",
        "",
        "```text",
    ]
    for row in diagnostics["selected_compact_taylor_moment_rows"]:
        lines.append(
            f"k={row['degree']}: moment={row['compact_transformed_gamma_moment_ball']}; "
            f"value contribution={row['value_moment_contribution_ball']}; "
            f"derivative contribution={row['derivative_moment_contribution_ball']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Remaining Splice",
            "",
            "```text",
            "Compose the separately certified n>30 Phi/Phi'/Phi(0) tail with matching normalization.",
            "Compose the separately certified y>200 far-tail bound.",
            "Then extend the compact/tail certificate from this worst row to every required grid row and the collar.",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_endpoint_x_moment_taylor_certificate"],
            artifact["source_endpoint_x_moment_taylor_json"],
            artifact["source_interpolation_remainder_route_matrix"],
            artifact["source_arb_chebyshev_interpolant_moment_scout"],
            artifact["source_compact_interval_integration_scout"],
            artifact["source_phi_tail_grid_certificate"],
            artifact["source_phi_tail_grid_json"],
            artifact["source_far_tail_split_certificate"],
            artifact["source_far_tail_split_json"],
            artifact["source_finite_plus_tail_budget_certificate"],
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
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--endpoint-certificate-json", type=Path, default=DEFAULT_ENDPOINT_CERTIFICATE_JSON
    )
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--far-tail-json", type=Path, default=DEFAULT_FAR_TAIL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    endpoint_certificate_path = (
        args.endpoint_certificate_json
        if args.endpoint_certificate_json.is_absolute()
        else REPO_ROOT / args.endpoint_certificate_json
    )
    phi_tail_path = (
        args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    )
    far_tail_path = (
        args.far_tail_json if args.far_tail_json.is_absolute() else REPO_ROOT / args.far_tail_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(endpoint_certificate_path, phi_tail_path, far_tail_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row compact x-moment Taylor "
        f"certificate: {out_json.relative_to(REPO_ROOT).as_posix()} and "
        f"{note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
