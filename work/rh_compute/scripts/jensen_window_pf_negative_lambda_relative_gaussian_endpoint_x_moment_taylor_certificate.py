#!/usr/bin/env python3
"""Build a rigorous endpoint x-moment Taylor certificate for the worst row."""

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
DEFAULT_ENDPOINT_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.json"
)
DEFAULT_PHI_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_X_RIGHT = "0.01"
DEFAULT_DISK_RADIUS = "0.2"
DEFAULT_TAYLOR_DEGREE = 64
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PRECISION_BITS = 8192
SELECTED_DEGREES = {0, 1, 2, 20, 40, 41, 42, 48, 56, 64}


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


def arb_mid_text(value: flint.arb, digits: int = 60) -> str:
    return value.str(digits, radius=False).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    return value.lower().str(digits, radius=False).replace("e", "E")


def arb_upper_text(value: flint.arb, digits: int = 60) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def abs_upper(value: flint.arb) -> flint.arb:
    lower_abs = abs(flint.arb(value.lower()))
    upper_abs = abs(flint.arb(value.upper()))
    return lower_abs if lower_abs > upper_abs else upper_abs


def source_paths(endpoint_route_path: Path, phi_tail_path: Path) -> dict[str, str]:
    return {
        "endpoint_route_json": endpoint_route_path.relative_to(REPO_ROOT).as_posix(),
        "endpoint_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.md"
        ),
        "phi_tail_json": phi_tail_path.relative_to(REPO_ROOT).as_posix(),
        "phi_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md"
        ),
        "endpoint_parity_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md"
        ),
        "interpolation_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
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


def transformed_gamma_moment(degree: int) -> flint.arb:
    """Integral of x^degree against the transformed first-panel Gamma density."""
    alpha = flint.arb(DEFAULT_ALPHA)
    shape = alpha + flint.arb(1) + flint.arb(degree) / flint.arb(2)
    scale = flint.arb(DEFAULT_T) ** (flint.arb(degree) / flint.arb(2))
    return flint.arb(1).gamma_lower(shape) / ((alpha + 1).gamma() * scale)


def complex_disk_majorant() -> dict[str, flint.arb | bool]:
    """Bound the finite Phi sum on |z|<=R using positive real damping."""
    radius = flint.arb(DEFAULT_DISK_RADIUS)
    pi = flint.arb.pi()
    radius_admissible = bool(flint.arb(4) * radius < pi / 2)
    real_exp_lower = (-flint.arb(4) * radius).exp() * (flint.arb(4) * radius).cos()
    if not radius_admissible or real_exp_lower.contains(0) or not bool(real_exp_lower > 0):
        raise RuntimeError("complex-disk damping lower bound was not certified positive")
    majorant = flint.arb(0)
    for n in range(1, DEFAULT_PHI_TERM_COUNT + 1):
        n2 = flint.arb(n * n)
        prefactor = (
            flint.arb(2) * pi * pi * n2 * n2 * (flint.arb(9) * radius).exp()
            + flint.arb(3) * pi * n2 * (flint.arb(5) * radius).exp()
        )
        majorant += prefactor * (-pi * n2 * real_exp_lower).exp()
    return {
        "radius": radius,
        "radius_admissible": radius_admissible,
        "real_exp_lower": real_exp_lower,
        "finite_phi_disk_majorant": majorant,
    }


def taylor_moment_rows(
    coefficients: list[flint.arb],
    ratios: list[flint.arb],
    c0: flint.arb,
) -> tuple[list[dict], flint.arb, flint.arb]:
    rows: list[dict] = []
    value_integral = flint.arb(0)
    derivative_integral = flint.arb(0)
    for degree, phi_coefficient in enumerate(coefficients):
        value_coefficient = phi_coefficient / c0
        derivative_coefficient = flint.arb(degree) * phi_coefficient / (flint.arb(2) * c0)
        subtracted_ratio = None
        if degree % 2 == 0 and degree // 2 <= DEFAULT_POLYNOMIAL_M:
            ratio_index = degree // 2
            subtracted_ratio = ratios[ratio_index]
            value_coefficient -= subtracted_ratio
            if degree:
                derivative_coefficient -= flint.arb(ratio_index) * subtracted_ratio
        moment = transformed_gamma_moment(degree)
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
                "transformed_gamma_moment_ball": arb_text(moment),
                "value_moment_contribution_ball": arb_text(value_contribution),
                "derivative_moment_contribution_ball": arb_text(derivative_contribution),
                "proof_boundary": (
                    "Arb Taylor-coefficient and exact transformed-moment row for the finite n<=30 core only."
                ),
            }
        )
    return rows, value_integral, derivative_integral


def build_diagnostics(endpoint_route: dict, phi_tail: dict) -> dict:
    configure_precision()
    coefficients = finite_phi_series(DEFAULT_PHI_TERM_COUNT, DEFAULT_TAYLOR_DEGREE)
    _ratio_rows, ratios = build_ratio_rows(2 * DEFAULT_POLYNOMIAL_M, DEFAULT_RATIO_CUTOFF_N)
    c0 = coefficients[0]
    if c0.contains(0) or not bool(c0 > 0):
        raise RuntimeError("finite Phi(0) did not certify positive")

    moment_rows, value_polynomial_integral, derivative_polynomial_integral = taylor_moment_rows(
        coefficients, ratios, c0
    )
    disk = complex_disk_majorant()
    x_right = flint.arb(DEFAULT_X_RIGHT)
    radius = disk["radius"]
    assert isinstance(radius, flint.arb)
    q = x_right / radius
    if not bool(q < 1):
        raise RuntimeError("Taylor evaluation radius is not inside the Cauchy disk")
    majorant = disk["finite_phi_disk_majorant"]
    assert isinstance(majorant, flint.arb)
    c0_lower = flint.arb(c0.lower())
    first_tail_degree = DEFAULT_TAYLOR_DEGREE + 1
    value_tail_sup = (
        majorant
        * (q**first_tail_degree)
        / (c0_lower * (flint.arb(1) - q))
    )
    derivative_tail_sup = (
        majorant
        * (q**first_tail_degree)
        * (
            flint.arb(first_tail_degree)
            - flint.arb(DEFAULT_TAYLOR_DEGREE) * q
        )
        / (
            flint.arb(2)
            * c0_lower
            * ((flint.arb(1) - q) ** 2)
        )
    )
    first_panel_mass = transformed_gamma_moment(0)
    value_tail_integral_radius = abs_upper(value_tail_sup * first_panel_mass)
    derivative_tail_integral_radius = abs_upper(derivative_tail_sup * first_panel_mass)
    value_total = value_polynomial_integral + flint.arb(0, value_tail_integral_radius)
    derivative_total = derivative_polynomial_integral + flint.arb(0, derivative_tail_integral_radius)
    value_abs = abs_upper(value_total)
    derivative_abs = abs_upper(derivative_total)

    route_summary = endpoint_route["summary"]
    value_cap = flint.arb(route_summary["value_unscaled_expectation_error_cap"])
    derivative_cap = flint.arb(route_summary["derivative_unscaled_expectation_error_cap"])
    value_ratio = value_abs / value_cap
    derivative_ratio = derivative_abs / derivative_cap
    value_below_cap = bool(value_ratio < 1)
    derivative_below_cap = bool(derivative_ratio < 1)
    both_below_caps = value_below_cap and derivative_below_cap
    if not both_below_caps:
        raise RuntimeError("first-panel Taylor certificate did not fit the allocated caps")

    selected_rows = [row for row in moment_rows if row["degree"] in SELECTED_DEGREES]
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "alpha": DEFAULT_ALPHA,
        "y_interval": "0<=y<=1",
        "x_interval": "0<=x<=0.01",
        "x_right_endpoint": DEFAULT_X_RIGHT,
        "precision_bits": DEFAULT_PRECISION_BITS,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
        "polynomial_M": DEFAULT_POLYNOMIAL_M,
        "taylor_degree": DEFAULT_TAYLOR_DEGREE,
        "first_tail_degree": first_tail_degree,
        "finite_c0_ball": arb_text(c0, 80),
        "finite_c0_lower": arb_lower_text(c0, 70),
        "disk_radius": DEFAULT_DISK_RADIUS,
        "four_disk_radius": arb_mid_text(flint.arb(4) * radius, 30),
        "pi_over_two_lower": arb_lower_text(flint.arb.pi() / 2, 30),
        "disk_radius_admissible": disk["radius_admissible"],
        "real_exp4_lower_on_disk": arb_lower_text(disk["real_exp_lower"], 60),
        "finite_phi_complex_disk_majorant_upper": arb_upper_text(majorant, 60),
        "x_to_disk_radius_ratio": arb_upper_text(q, 40),
        "first_panel_mass_ball": arb_text(first_panel_mass, 70),
        "first_panel_mass_upper": arb_upper_text(first_panel_mass, 70),
        "x_moment_formula": (
            "mu_k=T^(-k/2)*lower_gamma(alpha+1+k/2,1)/Gamma(alpha+1)"
        ),
        "complex_disk_bound_formula": (
            "For |z|<=R<pi/8, Re(exp(4z))>=exp(-4R)*cos(4R)>0; hence "
            "|Phi_30(z)| is bounded termwise by the recorded positive majorant."
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
        "taylor_moment_rows": moment_rows,
        "taylor_moment_row_count": len(moment_rows),
        "selected_taylor_moment_rows": selected_rows,
        "value_taylor_moment_integral_ball": arb_text(value_polynomial_integral, 80),
        "derivative_taylor_moment_integral_ball": arb_text(
            derivative_polynomial_integral, 80
        ),
        "value_first_panel_total_ball": arb_text(value_total, 80),
        "derivative_first_panel_total_ball": arb_text(derivative_total, 80),
        "value_first_panel_abs_upper": arb_upper_text(value_abs, 70),
        "derivative_first_panel_abs_upper": arb_upper_text(derivative_abs, 70),
        "value_first_panel_certified_negative": bool(value_total < 0),
        "derivative_first_panel_certified_negative": bool(derivative_total < 0),
        "value_unscaled_expectation_error_cap": route_summary[
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": route_summary[
            "derivative_unscaled_expectation_error_cap"
        ],
        "value_abs_to_cap_ratio_upper": arb_upper_text(value_ratio, 60),
        "derivative_abs_to_cap_ratio_upper": arb_upper_text(derivative_ratio, 60),
        "value_below_cap": value_below_cap,
        "derivative_below_cap": derivative_below_cap,
        "both_channels_below_caps": both_below_caps,
        "source_first_panel_mass_upper": route_summary["first_panel_mass_upper"],
        "source_tail_certificate_rows": phi_tail["summary"]["tail_source_rows"],
        "source_tail_certified": phi_tail["summary"]["finite_grid_tail_source_certified"],
        "finite_core_first_panel_certified": True,
        "full_n_tail_composed_here": False,
        "remaining_compact_panels": ["1<=y<=5", "5<=y<=20", "20<=y<=50", "50<=y<=100", "100<=y<=200"],
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgexmtc_01_route_scope_import",
            role="scope_reduction",
            readiness="available_scope_reduction",
            claim=(
                "The endpoint route maps the first Gamma panel 0<=y<=1 to 0<=x<=0.01 for the "
                "T=10000, F_21 worst row and isolates its finite n<=30 cancellation-reduced core."
            ),
            diagnostics={
                "T": diagnostics["T"],
                "index": diagnostics["index"],
                "y_interval": diagnostics["y_interval"],
                "x_interval": diagnostics["x_interval"],
                "phi_term_count": diagnostics["phi_term_count"],
                "polynomial_M": diagnostics["polynomial_M"],
            },
            source_artifacts=[paths["endpoint_route_json"], paths["endpoint_route_note"]],
            proof_boundary="Worst-row first-panel finite-core reduction only; not all compact panels or all rows.",
        ),
        MatrixRow(
            id="nlrgexmtc_02_exact_transformed_moments",
            role="exact_moment_reduction",
            readiness="available_exact_reduction",
            claim=(
                "Every x^k Taylor monomial is integrated against the transformed first-panel Gamma density "
                "by an exact incomplete-gamma identity, implemented with Arb ball arithmetic."
            ),
            diagnostics={
                "x_moment_formula": diagnostics["x_moment_formula"],
                "taylor_moment_rows": diagnostics["taylor_moment_row_count"],
                "first_panel_mass_ball": diagnostics["first_panel_mass_ball"],
                "selected_taylor_moment_rows": diagnostics["selected_taylor_moment_rows"],
            },
            source_artifacts=[paths["endpoint_route_note"], paths["finite_part_note"]],
            proof_boundary="Exact transformed-moment reduction for this first panel only.",
        ),
        MatrixRow(
            id="nlrgexmtc_03_complex_disk_majorant",
            role="analytic_remainder_certificate",
            readiness="available_analytic_certificate",
            claim=(
                "On |z|<=0.2, the finite n<=30 Phi kernel has a rigorous positive termwise complex-disk "
                "majorant because Re(exp(4z)) has the recorded positive lower bound."
            ),
            diagnostics={
                "complex_disk_bound_formula": diagnostics["complex_disk_bound_formula"],
                "disk_radius": diagnostics["disk_radius"],
                "four_disk_radius": diagnostics["four_disk_radius"],
                "pi_over_two_lower": diagnostics["pi_over_two_lower"],
                "disk_radius_admissible": diagnostics["disk_radius_admissible"],
                "real_exp4_lower_on_disk": diagnostics["real_exp4_lower_on_disk"],
                "finite_phi_complex_disk_majorant_upper": diagnostics[
                    "finite_phi_complex_disk_majorant_upper"
                ],
            },
            source_artifacts=[paths["endpoint_parity_note"], paths["interpolation_route_note"]],
            proof_boundary="Finite n<=30 complex-disk bound only; the infinite n-tail is a separate source.",
        ),
        MatrixRow(
            id="nlrgexmtc_04_taylor_moment_interval",
            role="arb_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The degree-64 finite-core Taylor polynomial, including subtraction of the degree-40 even "
                "ratio polynomial, is integrated by the exact transformed moments in both channels."
            ),
            diagnostics={
                "precision_bits": diagnostics["precision_bits"],
                "taylor_degree": diagnostics["taylor_degree"],
                "ratio_cutoff_n": diagnostics["ratio_cutoff_n"],
                "value_taylor_moment_integral_ball": diagnostics[
                    "value_taylor_moment_integral_ball"
                ],
                "derivative_taylor_moment_integral_ball": diagnostics[
                    "derivative_taylor_moment_integral_ball"
                ],
            },
            source_artifacts=[paths["finite_part_note"], paths["quadrature_route_note"]],
            proof_boundary="Finite-core first-panel Taylor-polynomial integral only; Cauchy tails are added separately.",
        ),
        MatrixRow(
            id="nlrgexmtc_05_first_panel_remainder_certificate",
            role="first_panel_certificate",
            readiness="available_first_panel_certificate",
            claim=(
                "Adding the rigorous degree-64 Cauchy tails certifies the complete finite n<=30 first-panel "
                "value and derivative integrals, with both absolute enclosures far below their global caps."
            ),
            diagnostics={
                "value_cauchy_tail_formula": diagnostics["value_cauchy_tail_formula"],
                "derivative_cauchy_tail_formula": diagnostics["derivative_cauchy_tail_formula"],
                "value_cauchy_tail_integral_radius_upper": diagnostics[
                    "value_cauchy_tail_integral_radius_upper"
                ],
                "derivative_cauchy_tail_integral_radius_upper": diagnostics[
                    "derivative_cauchy_tail_integral_radius_upper"
                ],
                "value_first_panel_total_ball": diagnostics["value_first_panel_total_ball"],
                "derivative_first_panel_total_ball": diagnostics[
                    "derivative_first_panel_total_ball"
                ],
                "value_abs_to_cap_ratio_upper": diagnostics["value_abs_to_cap_ratio_upper"],
                "derivative_abs_to_cap_ratio_upper": diagnostics[
                    "derivative_abs_to_cap_ratio_upper"
                ],
                "both_channels_below_caps": diagnostics["both_channels_below_caps"],
            },
            source_artifacts=[paths["endpoint_route_note"], paths["quadrature_route_note"]],
            proof_boundary=(
                "Complete finite n<=30 first-panel integral certificate only; not the n>30 tail composition, "
                "remaining compact panels, or a full quadrature-remainder theorem."
            ),
        ),
        MatrixRow(
            id="nlrgexmtc_06_tail_and_splice_handoff",
            role="compact_interval_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "The endpoint branch obstruction is retired for the finite n<=30 first panel; the compact "
                "route must now compose the separate n>30 tail source and certify the five panels with y>=1."
            ),
            diagnostics={
                "finite_core_first_panel_certified": diagnostics[
                    "finite_core_first_panel_certified"
                ],
                "source_tail_certified": diagnostics["source_tail_certified"],
                "full_n_tail_composed_here": diagnostics["full_n_tail_composed_here"],
                "remaining_compact_panels": diagnostics["remaining_compact_panels"],
                "target_closing": diagnostics["target_closing"],
            },
            source_artifacts=[
                paths["phi_tail_json"],
                paths["phi_tail_note"],
                paths["interpolation_route_note"],
                paths["intervalization_target"],
            ],
            proof_boundary="Handoff only; the full compact interval and uniform collar remain open.",
        ),
        MatrixRow(
            id="nlrgexmtc_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "This first-panel certificate may be used only as one local component of a compact-interval "
                "proof; it does not promote the remaining panels or collar-level targets."
            ),
            diagnostics={
                "forbidden_promotions": [
                    "full n>30 tail composition",
                    "remaining five compact panels",
                    "full compact interval certificate",
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


def build_artifact(endpoint_route_path: Path, phi_tail_path: Path) -> dict:
    endpoint_route = load_json(endpoint_route_path)
    phi_tail = load_json(phi_tail_path)
    paths = source_paths(endpoint_route_path, phi_tail_path)
    diagnostics = build_diagnostics(endpoint_route, phi_tail)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "alpha": diagnostics["alpha"],
        "y_interval": diagnostics["y_interval"],
        "x_interval": diagnostics["x_interval"],
        "precision_bits": diagnostics["precision_bits"],
        "phi_term_count": diagnostics["phi_term_count"],
        "polynomial_M": diagnostics["polynomial_M"],
        "taylor_degree": diagnostics["taylor_degree"],
        "first_tail_degree": diagnostics["first_tail_degree"],
        "taylor_moment_row_count": diagnostics["taylor_moment_row_count"],
        "disk_radius": diagnostics["disk_radius"],
        "disk_radius_admissible": diagnostics["disk_radius_admissible"],
        "real_exp4_lower_on_disk": diagnostics["real_exp4_lower_on_disk"],
        "finite_phi_complex_disk_majorant_upper": diagnostics[
            "finite_phi_complex_disk_majorant_upper"
        ],
        "x_to_disk_radius_ratio": diagnostics["x_to_disk_radius_ratio"],
        "first_panel_mass_upper": diagnostics["first_panel_mass_upper"],
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
        "value_first_panel_total_ball": diagnostics["value_first_panel_total_ball"],
        "derivative_first_panel_total_ball": diagnostics[
            "derivative_first_panel_total_ball"
        ],
        "value_first_panel_abs_upper": diagnostics["value_first_panel_abs_upper"],
        "derivative_first_panel_abs_upper": diagnostics[
            "derivative_first_panel_abs_upper"
        ],
        "value_first_panel_certified_negative": diagnostics[
            "value_first_panel_certified_negative"
        ],
        "derivative_first_panel_certified_negative": diagnostics[
            "derivative_first_panel_certified_negative"
        ],
        "value_unscaled_expectation_error_cap": diagnostics[
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "value_abs_to_cap_ratio_upper": diagnostics["value_abs_to_cap_ratio_upper"],
        "derivative_abs_to_cap_ratio_upper": diagnostics[
            "derivative_abs_to_cap_ratio_upper"
        ],
        "both_channels_below_caps": diagnostics["both_channels_below_caps"],
        "finite_core_first_panel_certified": diagnostics[
            "finite_core_first_panel_certified"
        ],
        "first_panel_certificate_rows": 1,
        "full_n_tail_composed_here": diagnostics["full_n_tail_composed_here"],
        "remaining_compact_panel_count": len(diagnostics["remaining_compact_panels"]),
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The endpoint x-moment Taylor certificate closes the finite n<=30 first-panel integration "
            "problem for the T=10000, F_21 worst row. A degree-64 Arb Taylor model is integrated by exact "
            "transformed Gamma moments, and a rigorous |z|<=0.2 Cauchy majorant bounds the true-function "
            "tail. Both first-panel channels are certified negative and far below their global error caps. "
            "This retires the endpoint branch obstruction for the finite core only; the separate n>30 tail "
            "composition, five y>=1 panels, all-row coverage, and uniform collar theorem remain open."
        ),
    }
    return {
        "kind": (
            "jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate"
        ),
        "date": "2026-07-09",
        "status": "endpoint x-moment Taylor certificate",
        "source_endpoint_x_panel_route_matrix": paths["endpoint_route_note"],
        "source_endpoint_x_panel_route_json": paths["endpoint_route_json"],
        "source_phi_tail_grid_certificate": paths["phi_tail_note"],
        "source_phi_tail_grid_json": paths["phi_tail_json"],
        "source_endpoint_parity_repair_matrix": paths["endpoint_parity_note"],
        "source_interpolation_remainder_route_matrix": paths["interpolation_route_note"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py"
        ),
        "proof_boundary": (
            "Endpoint x-moment Taylor certificate for the finite n<=30 first-panel core only. It certifies "
            "the T=10000, F_21 value and derivative integrals on 0<=y<=1 by exact transformed Gamma moments "
            "and a complex-disk Cauchy remainder. It does not compose the n>30 tail in this artifact, does "
            "not certify the five panels with y>=1, does not prove the full compact-interval integral, does "
            "not certify all finite-grid rows, does not prove the uniform collar theorem, does not prove "
            "scaled-curvature monotonicity or cone entry, and does not prove RH or Lambda <= 0."
        ),
        "diagnostics": diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "The x-panel certificate applies only to the finite n<=30 cancellation-reduced core.",
            "The transformed moment identity is evaluated with Arb for every degree 0 through 64.",
            "The Cauchy disk satisfies 4R<pi/2, so the damping majorant has a positive real-part lower bound.",
            "Both first-panel total balls include the full post-degree-64 Cauchy tail.",
            "The separate n>30 Phi-tail source is not silently absorbed into this finite-core certificate.",
            "The five panels with y>=1 and the full compact interval remain open.",
            "No row is ready_to_apply for the full intervalization or uniform-collar target.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor "
        f"certificate: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['taylor_moment_row_count']} exact-moment rows, "
        f"{summary['first_panel_certificate_rows']} certified first panel, "
        f"{summary['remaining_compact_panel_count']} open compact panels, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint X-Moment Taylor Certificate",
        "",
        "Date: 2026-07-09",
        "",
        "Status: endpoint x-moment Taylor certificate. This is not a proof",
        "of the full compact interval, a finite-grid interval certificate, a uniform",
        "collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate`.",
        "",
        "Proof boundary: this artifact certifies the finite `n<=30` cancellation-reduced",
        "value and derivative cores on the first panel `0<=y<=1` for `T=10000`,",
        "`F_21`. It does not compose the separate `n>30` tail here and does not",
        "certify the five panels with `y>=1`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_moment_taylor_certificate.py",
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
        f"Phi finite terms: {summary['phi_term_count']}",
        f"subtracted polynomial M: {summary['polynomial_M']}",
        f"Taylor degree: {summary['taylor_degree']}",
        f"Arb precision bits: {summary['precision_bits']}",
        f"exact transformed-moment rows: {summary['taylor_moment_row_count']}",
        "```",
        "",
        "Exact transformed moment:",
        "",
        "```text",
        diagnostics["x_moment_formula"],
        "```",
        "",
        "## Complex-Disk Remainder",
        "",
        "```text",
        diagnostics["complex_disk_bound_formula"],
        f"disk radius R: {summary['disk_radius']}",
        f"4R<pi/2 certified: {summary['disk_radius_admissible']}",
        f"Re(exp(4z)) lower: {summary['real_exp4_lower_on_disk']}",
        f"finite Phi_30 disk majorant upper: {summary['finite_phi_complex_disk_majorant_upper']}",
        f"x/R: {summary['x_to_disk_radius_ratio']}",
        f"value integrated Cauchy radius upper: {summary['value_cauchy_tail_integral_radius_upper']}",
        f"derivative integrated Cauchy radius upper: {summary['derivative_cauchy_tail_integral_radius_upper']}",
        "```",
        "",
        "## First-Panel Enclosures",
        "",
        "```text",
        f"value Taylor-moment integral: {summary['value_taylor_moment_integral_ball']}",
        f"derivative Taylor-moment integral: {summary['derivative_taylor_moment_integral_ball']}",
        f"value total ball: {summary['value_first_panel_total_ball']}",
        f"derivative total ball: {summary['derivative_first_panel_total_ball']}",
        f"value certified negative: {summary['value_first_panel_certified_negative']}",
        f"derivative certified negative: {summary['derivative_first_panel_certified_negative']}",
        f"value absolute/cap ratio upper: {summary['value_abs_to_cap_ratio_upper']}",
        f"derivative absolute/cap ratio upper: {summary['derivative_abs_to_cap_ratio_upper']}",
        f"both channels below caps: {summary['both_channels_below_caps']}",
        "```",
        "",
        "Selected coefficient/moment rows:",
        "",
        "```text",
    ]
    for row in diagnostics["selected_taylor_moment_rows"]:
        lines.append(
            f"k={row['degree']}: moment={row['transformed_gamma_moment_ball']}; "
            f"value contribution={row['value_moment_contribution_ball']}; "
            f"derivative contribution={row['derivative_moment_contribution_ball']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Remaining Work",
            "",
            "```text",
            "Compose the already separate n>30 Phi-tail source without changing normalization.",
            "Certify the five compact panels 1<=y<=200 by analytic-domain/Taylor or Bernstein bounds.",
            "Lift the worst-row compact certificate to all required rows and then to the full collar.",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_endpoint_x_panel_route_matrix"],
            artifact["source_endpoint_x_panel_route_json"],
            artifact["source_phi_tail_grid_certificate"],
            artifact["source_phi_tail_grid_json"],
            artifact["source_endpoint_parity_repair_matrix"],
            artifact["source_interpolation_remainder_route_matrix"],
            artifact["source_finite_part_weighted_sum_interval_certificate"],
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
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint-route-json", type=Path, default=DEFAULT_ENDPOINT_ROUTE_JSON)
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    endpoint_route_path = (
        args.endpoint_route_json
        if args.endpoint_route_json.is_absolute()
        else REPO_ROOT / args.endpoint_route_json
    )
    phi_tail_path = (
        args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(endpoint_route_path, phi_tail_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian endpoint x-moment Taylor "
        f"certificate: {out_json.relative_to(REPO_ROOT).as_posix()} and "
        f"{note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
