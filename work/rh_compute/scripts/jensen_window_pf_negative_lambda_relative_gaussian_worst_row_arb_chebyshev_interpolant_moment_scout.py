#!/usr/bin/env python3
"""Build a worst-row Arb Chebyshev interpolant-moment scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
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
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate import (  # noqa: E402
    finite_phi,
    finite_phip,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FLOATING_SCOUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.json"
)
DEFAULT_QUADRATURE_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_PRECISION_BITS = 1024
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_PANELS = ((0, 1), (1, 5), (5, 20), (20, 50), (50, 100), (100, 200))
DEFAULT_DEGREES = (16, 20, 24, 32)
REFERENCE_DEGREE = DEFAULT_DEGREES[-1]


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


def arb_text(value: flint.arb, digits: int = 50) -> str:
    return value.str(digits).replace("e", "E")


def arb_mid_text(value: flint.arb, digits: int = 50) -> str:
    return value.str(digits, radius=False).replace("e", "E")


def arb_upper_text(value: flint.arb, digits: int = 50) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def abs_upper(value: flint.arb) -> flint.arb:
    lower_abs = abs(flint.arb(value.lower()))
    upper_abs = abs(flint.arb(value.upper()))
    return lower_abs if lower_abs > upper_abs else upper_abs


def source_paths(floating_scout_path: Path, quadrature_route_path: Path) -> dict[str, str]:
    return {
        "floating_scout_json": floating_scout_path.relative_to(REPO_ROOT).as_posix(),
        "floating_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md"
        ),
        "quadrature_route_json": quadrature_route_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"
        ),
        "compact_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md"
        ),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
        ),
        "intervalization_target": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


CHEB_POWER_CACHE: dict[int, list[flint.arb]] = {}


def chebyshev_power_coefficients(degree: int) -> list[flint.arb]:
    if degree in CHEB_POWER_CACHE:
        return CHEB_POWER_CACHE[degree]
    if degree == 0:
        coefficients = [flint.arb(1)]
    elif degree == 1:
        coefficients = [flint.arb(0), flint.arb(1)]
    else:
        previous = [flint.arb(1)]
        current = [flint.arb(0), flint.arb(1)]
        for _ in range(1, degree):
            twice_t_current = [flint.arb(0)] + [2 * item for item in current]
            if len(previous) < len(twice_t_current):
                previous = previous + [flint.arb(0)] * (len(twice_t_current) - len(previous))
            next_coefficients = [twice_t_current[i] - previous[i] for i in range(len(twice_t_current))]
            previous, current = current, next_coefficients
        coefficients = current
    CHEB_POWER_CACHE[degree] = coefficients
    return coefficients


def configure_precision() -> None:
    flint.ctx.prec = DEFAULT_PRECISION_BITS


def panel_chebyshev_moments(left: int, right: int, degree: int) -> list[flint.arb]:
    a = flint.arb(left)
    b = flint.arb(right)
    center = (a + b) / 2
    half_width = (b - a) / 2
    alpha = flint.arb(DEFAULT_ALPHA)
    gamma_alpha = (alpha + 1).gamma()
    y_moments: list[flint.arb] = []
    for p in range(degree + 1):
        shape = alpha + flint.arb(p) + 1
        y_moments.append((b.gamma_lower(shape) - a.gamma_lower(shape)) / gamma_alpha)
    t_moments: list[flint.arb] = []
    for m in range(degree + 1):
        total = flint.arb(0)
        for p in range(m + 1):
            total += flint.arb(math.comb(m, p)) * ((-center) ** (m - p)) * y_moments[p]
        t_moments.append(total / (half_width**m))
    chebyshev_moments: list[flint.arb] = []
    for k in range(degree + 1):
        coefficients = chebyshev_power_coefficients(k)
        chebyshev_moments.append(
            sum((coefficients[m] * t_moments[m] for m in range(len(coefficients))), flint.arb(0))
        )
    return chebyshev_moments


def make_core_functions(ratios: list[flint.arb], c0_finite: flint.arb):
    t = flint.arb(DEFAULT_T)

    def value_core(y_value: flint.arb) -> flint.arb:
        v = y_value / t
        x = v.sqrt()
        polynomial = flint.arb(0)
        for j in range(DEFAULT_POLYNOMIAL_M + 1):
            polynomial += ratios[j] * (v**j)
        return finite_phi(x, DEFAULT_PHI_TERM_COUNT) / c0_finite - polynomial

    def derivative_core(y_value: flint.arb) -> flint.arb:
        v = y_value / t
        x = v.sqrt()
        polynomial = flint.arb(0)
        for j in range(1, DEFAULT_POLYNOMIAL_M + 1):
            polynomial += flint.arb(j) * ratios[j] * (v**j)
        return x * finite_phip(x, DEFAULT_PHI_TERM_COUNT) / (flint.arb(2) * c0_finite) - polynomial

    return value_core, derivative_core


def chebyshev_coefficients(func, left: int, right: int, degree: int) -> list[flint.arb]:
    n_nodes = degree + 1
    a = flint.arb(left)
    b = flint.arb(right)
    pi = flint.arb.pi()
    theta_values: list[flint.arb] = []
    function_values: list[flint.arb] = []
    for j in range(n_nodes):
        theta = pi * (flint.arb(2 * j + 1) / flint.arb(2 * n_nodes))
        t = theta.cos()
        y = (a + b) / 2 + (b - a) * t / 2
        theta_values.append(theta)
        function_values.append(func(y))
    coefficients: list[flint.arb] = []
    for k in range(n_nodes):
        total = flint.arb(0)
        for j in range(n_nodes):
            total += function_values[j] * (flint.arb(k) * theta_values[j]).cos()
        coefficients.append((flint.arb(2) / flint.arb(n_nodes)) * total)
    return coefficients


def panel_integral_from_coefficients(coefficients: list[flint.arb], moments: list[flint.arb]) -> flint.arb:
    return coefficients[0] * moments[0] / 2 + sum(
        (coefficients[k] * moments[k] for k in range(1, len(coefficients))),
        flint.arb(0),
    )


def panel_integral(func, left: int, right: int, degree: int) -> tuple[flint.arb, dict]:
    coefficients = chebyshev_coefficients(func, left, right, degree)
    moments = panel_chebyshev_moments(left, right, degree)
    integral = panel_integral_from_coefficients(coefficients, moments)
    max_coefficient_abs = max(abs_upper(item) for item in coefficients)
    max_moment_abs = max(abs_upper(item) for item in moments)
    return integral, {
        "max_interpolant_coefficient_abs_upper": arb_upper_text(max_coefficient_abs, 30),
        "max_panel_chebyshev_moment_abs_upper": arb_upper_text(max_moment_abs, 30),
    }


def build_ladder_rows(value_core, derivative_core) -> tuple[list[dict], dict[int, tuple[flint.arb, flint.arb]]]:
    rows: list[dict] = []
    estimates: dict[int, tuple[flint.arb, flint.arb]] = {}
    for degree in DEFAULT_DEGREES:
        value_total = flint.arb(0)
        derivative_total = flint.arb(0)
        for left, right in DEFAULT_PANELS:
            value_panel, _value_meta = panel_integral(value_core, left, right, degree)
            derivative_panel, _derivative_meta = panel_integral(derivative_core, left, right, degree)
            value_total += value_panel
            derivative_total += derivative_panel
        estimates[degree] = (value_total, derivative_total)
        rows.append(
            {
                "degree": degree,
                "panel_count": len(DEFAULT_PANELS),
                "value_interpolant_integral_ball": arb_text(value_total, 60),
                "derivative_interpolant_integral_ball": arb_text(derivative_total, 60),
                "proof_boundary": (
                    "Arb-enclosed Chebyshev interpolant integral only; not an interpolation-remainder "
                    "or true-function integral certificate."
                ),
            }
        )
    return rows, estimates


def build_cauchy_rows(
    estimates: dict[int, tuple[flint.arb, flint.arb]],
    value_cap: flint.arb,
    derivative_cap: flint.arb,
) -> list[dict]:
    rows: list[dict] = []
    for lower, upper in zip(DEFAULT_DEGREES, DEFAULT_DEGREES[1:]):
        value_delta = abs_upper(estimates[upper][0] - estimates[lower][0])
        derivative_delta = abs_upper(estimates[upper][1] - estimates[lower][1])
        value_ratio = value_delta / value_cap
        derivative_ratio = derivative_delta / derivative_cap
        rows.append(
            {
                "pair": f"{lower}->{upper}",
                "lower_degree": lower,
                "upper_degree": upper,
                "value_delta_upper": arb_upper_text(value_delta),
                "derivative_delta_upper": arb_upper_text(derivative_delta),
                "value_delta_to_cap_ratio_upper": arb_upper_text(value_ratio, 40),
                "derivative_delta_to_cap_ratio_upper": arb_upper_text(derivative_ratio, 40),
                "both_deltas_below_caps": bool(value_ratio < 1 and derivative_ratio < 1),
                "proof_boundary": (
                    "Arb Cauchy-ladder comparison for interpolants only; not a bound for the "
                    "interpolation remainder."
                ),
            }
        )
    return rows


def build_panel_rows(value_core, derivative_core) -> list[dict]:
    rows: list[dict] = []
    for left, right in DEFAULT_PANELS:
        value_panel, value_meta = panel_integral(value_core, left, right, REFERENCE_DEGREE)
        derivative_panel, derivative_meta = panel_integral(derivative_core, left, right, REFERENCE_DEGREE)
        rows.append(
            {
                "panel": f"{left}<=y<={right}",
                "degree": REFERENCE_DEGREE,
                "value_interpolant_contribution_ball": arb_text(value_panel, 50),
                "derivative_interpolant_contribution_ball": arb_text(derivative_panel, 50),
                "value_max_interpolant_coefficient_abs_upper": value_meta[
                    "max_interpolant_coefficient_abs_upper"
                ],
                "derivative_max_interpolant_coefficient_abs_upper": derivative_meta[
                    "max_interpolant_coefficient_abs_upper"
                ],
                "max_panel_chebyshev_moment_abs_upper": value_meta[
                    "max_panel_chebyshev_moment_abs_upper"
                ],
                "proof_boundary": "Reference panel contribution for the interpolant only.",
            }
        )
    return rows


def build_diagnostics(floating_scout: dict, quadrature_route: dict) -> dict:
    configure_precision()
    _ratio_rows, ratios = build_ratio_rows(2 * (DEFAULT_POLYNOMIAL_M + 1), DEFAULT_RATIO_CUTOFF_N)
    c0_finite = finite_phi(flint.arb(0), DEFAULT_PHI_TERM_COUNT)
    if c0_finite.contains(0) or not bool(c0_finite > 0):
        raise RuntimeError("finite Phi(0) interval did not certify positive")
    value_core, derivative_core = make_core_functions(ratios, c0_finite)
    value_cap = flint.arb(quadrature_route["summary"]["value_unscaled_expectation_error_cap"])
    derivative_cap = flint.arb(quadrature_route["summary"]["derivative_unscaled_expectation_error_cap"])
    ladder_rows, estimates = build_ladder_rows(value_core, derivative_core)
    cauchy_rows = build_cauchy_rows(estimates, value_cap, derivative_cap)
    panel_rows = build_panel_rows(value_core, derivative_core)
    cap_safe_pairs = [row for row in cauchy_rows if row["both_deltas_below_caps"]]
    reference_value = estimates[REFERENCE_DEGREE][0]
    reference_derivative = estimates[REFERENCE_DEGREE][1]
    floating_summary = floating_scout["summary"]
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "u": "1/10000",
        "alpha": DEFAULT_ALPHA,
        "precision_bits": DEFAULT_PRECISION_BITS,
        "polynomial_M": DEFAULT_POLYNOMIAL_M,
        "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "panel_count": len(DEFAULT_PANELS),
        "degree_count": len(DEFAULT_DEGREES),
        "degrees": list(DEFAULT_DEGREES),
        "cauchy_pair_count": len(cauchy_rows),
        "cap_safe_pair_count": len(cap_safe_pairs),
        "first_cap_safe_pair": cap_safe_pairs[0]["pair"] if cap_safe_pairs else None,
        "reference_degree": REFERENCE_DEGREE,
        "compact_interval": floating_summary["compact_interval"],
        "floating_reference_degree": floating_summary["reference_degree"],
        "floating_reference_value_integral_estimate": floating_summary[
            "reference_value_integral_estimate"
        ],
        "floating_reference_derivative_integral_estimate": floating_summary[
            "reference_derivative_integral_estimate"
        ],
        "value_unscaled_expectation_error_cap": quadrature_route["summary"][
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": quadrature_route["summary"][
            "derivative_unscaled_expectation_error_cap"
        ],
        "reference_value_interpolant_integral_ball": arb_text(reference_value, 60),
        "reference_derivative_interpolant_integral_ball": arb_text(reference_derivative, 60),
        "reference_value_interpolant_mid": arb_mid_text(reference_value, 60),
        "reference_derivative_interpolant_mid": arb_mid_text(reference_derivative, 60),
        "ladder_rows": ladder_rows,
        "cauchy_rows": cauchy_rows,
        "panel_rows": panel_rows,
        "recommended_upgrade": (
            "Complete this route by bounding the difference between each Arb Chebyshev interpolant and "
            "the true cancellation-reduced core on its panel."
        ),
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwracims_01_floating_scout_import",
            role="scope_reduction",
            readiness="diagnostic_only",
            claim=(
                "The floating Chebyshev panel-moment scout identifies the six-panel compact route and "
                "degree ladder to promote into Arb arithmetic."
            ),
            diagnostics={
                "compact_interval": diagnostics["compact_interval"],
                "floating_reference_degree": diagnostics["floating_reference_degree"],
                "floating_reference_value_integral_estimate": diagnostics[
                    "floating_reference_value_integral_estimate"
                ],
                "floating_reference_derivative_integral_estimate": diagnostics[
                    "floating_reference_derivative_integral_estimate"
                ],
            },
            source_artifacts=[paths["floating_scout_json"], paths["floating_scout_note"]],
            proof_boundary="Scope import only; not a compact interval-integration proof.",
        ),
        MatrixRow(
            id="nlrgwracims_02_arb_interpolant_moment_arithmetic",
            role="arb_interpolant_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "For each panel and degree, Arb evaluates the cancellation-reduced core at Chebyshev "
                "nodes, forms the interpolant coefficients, and integrates that interpolant against "
                "Arb incomplete-Gamma panel moments."
            ),
            diagnostics={
                "precision_bits": diagnostics["precision_bits"],
                "panel_count": diagnostics["panel_count"],
                "degrees": diagnostics["degrees"],
                "moment_formula": (
                    "Chebyshev basis moments are expanded into monomials in y and evaluated by "
                    "lower incomplete-Gamma differences."
                ),
            },
            source_artifacts=[paths["quadrature_route_note"], paths["finite_part_note"]],
            proof_boundary="Arb interpolant arithmetic only; not a true-function panel enclosure.",
        ),
        MatrixRow(
            id="nlrgwracims_03_arb_degree_ladder",
            role="arb_interpolant_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "The Arb interpolant-integral ladder for degrees 16,20,24,32 agrees with the floating "
                "Chebyshev calibration while preserving ball arithmetic for coefficients and moments."
            ),
            diagnostics={"ladder_rows": diagnostics["ladder_rows"]},
            source_artifacts=[paths["floating_scout_note"], paths["quadrature_route_json"]],
            proof_boundary="Arb interpolant ladder only; interpolation remainders are not bounded.",
        ),
        MatrixRow(
            id="nlrgwracims_04_cauchy_delta_vs_caps",
            role="arb_interpolant_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "All consecutive Arb interpolant Cauchy deltas in the degree ladder are below the "
                "unscaled quadrature route caps in both value and derivative channels."
            ),
            diagnostics={
                "cauchy_rows": diagnostics["cauchy_rows"],
                "cap_safe_pair_count": diagnostics["cap_safe_pair_count"],
                "first_cap_safe_pair": diagnostics["first_cap_safe_pair"],
            },
            source_artifacts=[paths["quadrature_route_json"], paths["quadrature_route_note"]],
            proof_boundary="Cauchy-delta diagnostic only; not an interpolation-remainder theorem.",
        ),
        MatrixRow(
            id="nlrgwracims_05_reference_panel_contributions",
            role="arb_interpolant_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "At reference degree 32, Arb panel contribution balls identify the central panels as the "
                "dominant compact-interval contribution for the interpolant."
            ),
            diagnostics={
                "reference_degree": diagnostics["reference_degree"],
                "reference_value_interpolant_integral_ball": diagnostics[
                    "reference_value_interpolant_integral_ball"
                ],
                "reference_derivative_interpolant_integral_ball": diagnostics[
                    "reference_derivative_interpolant_integral_ball"
                ],
                "panel_rows": diagnostics["panel_rows"],
            },
            source_artifacts=[paths["compact_scout_note"], paths["far_tail_note"]],
            proof_boundary="Reference interpolant contribution only; not a certified true integral.",
        ),
        MatrixRow(
            id="nlrgwracims_06_remainder_gap",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The Arb Chebyshev interpolant-moment ladder proves the compact-interval integration source.",
            gap=(
                "The interpolant arithmetic is Arb-enclosed, but the difference between each interpolant "
                "and the true cancellation-reduced core on its panel is not yet bounded."
            ),
            source_artifacts=[paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Rejected promotion only; not a compact interval certificate, RH, or Lambda <= 0.",
        ),
        MatrixRow(
            id="nlrgwracims_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted compact-interval certificate may reuse this Arb interpolant arithmetic only "
                "after a panel remainder theorem supplies true-function enclosures."
            ),
            diagnostics={"recommended_upgrade": diagnostics["recommended_upgrade"]},
            source_artifacts=[paths["dependency_graph"], paths["intervalization_target"]],
            proof_boundary="Proof-hygiene gate only; not a finite-grid certificate, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(floating_scout_path: Path, quadrature_route_path: Path) -> dict:
    floating_scout = load_json(floating_scout_path)
    quadrature_route = load_json(quadrature_route_path)
    paths = source_paths(floating_scout_path, quadrature_route_path)
    diagnostics = build_diagnostics(floating_scout, quadrature_route)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "u": diagnostics["u"],
        "alpha": diagnostics["alpha"],
        "precision_bits": diagnostics["precision_bits"],
        "polynomial_M": diagnostics["polynomial_M"],
        "ratio_cutoff_n": diagnostics["ratio_cutoff_n"],
        "phi_term_count": diagnostics["phi_term_count"],
        "panel_count": diagnostics["panel_count"],
        "degree_count": diagnostics["degree_count"],
        "degrees": diagnostics["degrees"],
        "cauchy_pair_count": diagnostics["cauchy_pair_count"],
        "cap_safe_pair_count": diagnostics["cap_safe_pair_count"],
        "first_cap_safe_pair": diagnostics["first_cap_safe_pair"],
        "reference_degree": diagnostics["reference_degree"],
        "compact_interval": diagnostics["compact_interval"],
        "value_unscaled_expectation_error_cap": diagnostics["value_unscaled_expectation_error_cap"],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "reference_value_interpolant_integral_ball": diagnostics[
            "reference_value_interpolant_integral_ball"
        ],
        "reference_derivative_interpolant_integral_ball": diagnostics[
            "reference_derivative_interpolant_integral_ball"
        ],
        "reference_value_interpolant_mid": diagnostics["reference_value_interpolant_mid"],
        "reference_derivative_interpolant_mid": diagnostics["reference_derivative_interpolant_mid"],
        "ready_to_apply_rows": 0,
        "target_closing": diagnostics["target_closing"],
        "recommended_upgrade": diagnostics["recommended_upgrade"],
        "main_finding": (
            "The Arb Chebyshev interpolant-moment route upgrades the floating scout's arithmetic: "
            "degree 16,20,24,32 interpolant integrals are Arb-enclosed, all consecutive Cauchy deltas "
            "are below the unscaled quadrature caps, and the degree-32 interpolant balls are "
            f"{diagnostics['reference_value_interpolant_integral_ball']} and "
            f"{diagnostics['reference_derivative_interpolant_integral_ball']}. This still does not "
            "close the compact interval because panel interpolation remainders remain unbounded."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout",
        "date": "2026-07-08",
        "status": "worst-row Arb Chebyshev interpolant-moment scout",
        "source_floating_chebyshev_panel_moment_scout": paths["floating_scout_note"],
        "source_floating_chebyshev_panel_moment_scout_json": paths["floating_scout_json"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_quadrature_remainder_route_json": paths["quadrature_route_json"],
        "source_compact_interval_integration_scout": paths["compact_scout_note"],
        "source_far_tail_split_certificate": paths["far_tail_note"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py"
        ),
        "proof_boundary": (
            "Worst-row Arb Chebyshev interpolant-moment scout only. It certifies the arithmetic of "
            "Chebyshev interpolant coefficients and incomplete-Gamma panel moments on 0<=y<=200, but "
            "it does not bound the interpolation remainder, does not prove the true compact integral, "
            "does not prove a quadrature-remainder theorem, does not cover all rows or quadrature orders, "
            "does not prove a finite-grid interval certificate, does not prove a uniform collar theorem, "
            "and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "Arb interpolant arithmetic is not promoted to a true-function integral enclosure.",
            "Panel interpolation remainders remain open.",
            "The compact interval 0<=y<=200 remains open until those remainders are certified.",
            "All-row coverage, rounding aggregation, and grid-to-collar coverage remain separate.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev "
        f"interpolant-moment scout: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['degree_count']} degrees, {summary['cauchy_pair_count']} Cauchy pairs, "
        f"{summary['cap_safe_pair_count']} cap-safe pairs, {summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Arb Chebyshev Interpolant-Moment Scout",
        "",
        "Date: 2026-07-08",
        "",
        "Status: worst-row Arb Chebyshev interpolant-moment scout. This is not a proof",
        "of a compact interval-integration certificate, quadrature-remainder",
        "theorem, finite-grid interval certificate, uniform collar theorem, RH,",
        "or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout`.",
        "",
        "Proof boundary: this artifact certifies Arb arithmetic for Chebyshev",
        "interpolants and incomplete-Gamma panel moments. It does not bound the",
        "interpolation remainder between the interpolant and the true core.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Arb Interpolant Ladder",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"compact interval: {summary['compact_interval']}",
        f"panel count: {summary['panel_count']}",
        f"degrees: {summary['degrees']}",
        f"value unscaled cap: {summary['value_unscaled_expectation_error_cap']}",
        f"derivative unscaled cap: {summary['derivative_unscaled_expectation_error_cap']}",
        f"first cap-safe pair: {summary['first_cap_safe_pair']}",
        f"cap-safe pair count: {summary['cap_safe_pair_count']}",
        f"reference degree: {summary['reference_degree']}",
        f"reference value interpolant integral ball: {summary['reference_value_interpolant_integral_ball']}",
        f"reference derivative interpolant integral ball: {summary['reference_derivative_interpolant_integral_ball']}",
        f"target closing: {summary['target_closing']}",
        "```",
        "",
        "Cauchy rows:",
        "",
        "```text",
    ]
    for row in artifact["matrix_rows"][3]["diagnostics"]["cauchy_rows"]:
        lines.append(
            f"{row['pair']}: value delta/cap {row['value_delta_to_cap_ratio_upper']}; "
            f"derivative delta/cap {row['derivative_delta_to_cap_ratio_upper']}; "
            f"below caps {row['both_deltas_below_caps']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Required upgrade:",
            "",
            "```text",
            summary["recommended_upgrade"],
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_floating_chebyshev_panel_moment_scout"],
            artifact["source_floating_chebyshev_panel_moment_scout_json"],
            artifact["source_quadrature_remainder_route_matrix"],
            artifact["source_quadrature_remainder_route_json"],
            artifact["source_compact_interval_integration_scout"],
            artifact["source_far_tail_split_certificate"],
            artifact["source_finite_part_weighted_sum_interval_certificate"],
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
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--floating-scout-json", type=Path, default=DEFAULT_FLOATING_SCOUT_JSON)
    parser.add_argument("--quadrature-route-json", type=Path, default=DEFAULT_QUADRATURE_ROUTE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    floating_scout_path = (
        args.floating_scout_json if args.floating_scout_json.is_absolute() else REPO_ROOT / args.floating_scout_json
    )
    quadrature_route_path = (
        args.quadrature_route_json
        if args.quadrature_route_json.is_absolute()
        else REPO_ROOT / args.quadrature_route_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(floating_scout_path, quadrature_route_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row Arb Chebyshev "
        f"interpolant-moment scout: {out_json.relative_to(REPO_ROOT).as_posix()} and "
        f"{note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
