#!/usr/bin/env python3
"""Certify the fixed-k stencil system uniformly on 1156<=T<=10000."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate import (  # noqa: E402
    DEFAULT_PHI_TERM_COUNT,
    DEFAULT_RATIO_CUTOFF_N,
    arb_lower_text,
    arb_text,
    arb_upper_text,
)
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)
from jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix import (  # noqa: E402
    finite_phi_series,
)
from jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate import (  # noqa: E402
    perturbation_ledger,
    finite_collar_data,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate import (  # noqa: E402
    abs_upper,
    build_phi_bounds,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate import (  # noqa: E402
    global_n_tail_bounds,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RECORDED_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.md"
)

DEFAULT_T_LEFT = 1156
DEFAULT_T_TRANSITION_NUMERATOR = 15625
DEFAULT_T_TRANSITION_DENOMINATOR = 8
DEFAULT_T_RIGHT = 10000
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_TAYLOR_DEGREE = 384
DEFAULT_DISK_RADIUS = "0.38"
DEFAULT_X_CAP_NUMERATOR = 8
DEFAULT_X_CAP_DENOMINATOR = 25
DEFAULT_FIXED_SPLIT_Y = 200
DEFAULT_C0_LOWER = "0.44"
DEFAULT_PRECISION_BITS = 4096
VALUE_BUDGET = "0.5"
DERIVATIVE_BUDGET = "0.009"
SELECTED_DEGREES = {0, 1, 2, 20, 40, 41, 42, 64, 128, 192, 256, 320, 384}


@dataclass(frozen=True)
class MatrixRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    source_artifacts: list[str]
    diagnostics: dict | None = None
    formula: str | None = None
    gap: str | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(recorded_grid_path: Path) -> dict[str, str]:
    return {
        "recorded_grid_json": recorded_grid_path.relative_to(REPO_ROOT).as_posix(),
        "recorded_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.md"
        ),
        "all_row_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.md"
        ),
        "compact_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.md"
        ),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "degree40_budget_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md"
        ),
        "degree40_ladder_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "formal_tail_obstruction": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"
        ),
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def configure_precision() -> None:
    flint.ctx.prec = DEFAULT_PRECISION_BITS


def complex_disk_majorant() -> dict[str, flint.arb | bool]:
    radius = flint.arb(DEFAULT_DISK_RADIUS)
    pi = flint.arb.pi()
    admissible = bool(4 * radius < pi / 2)
    real_exp_lower = (-4 * radius).exp() * (4 * radius).cos()
    if not admissible or not bool(real_exp_lower > 0):
        raise RuntimeError("common disk failed")
    majorant = flint.arb(0)
    for n in range(1, DEFAULT_PHI_TERM_COUNT + 1):
        n2 = flint.arb(n * n)
        majorant += (
            2 * pi * pi * n2 * n2 * (9 * radius).exp()
            + 3 * pi * n2 * (5 * radius).exp()
        ) * (-pi * n2 * real_exp_lower).exp()
    return {
        "radius": radius,
        "admissible": admissible,
        "real_exp_lower": real_exp_lower,
        "majorant": majorant,
    }


def max_power_on_interval(exponent: flint.arb, left: flint.arb, right: flint.arb) -> flint.arb:
    return right**exponent if bool(exponent >= 0) else left**exponent


def compact_polynomial_majorant_rows(
    index: int,
    coefficients: list[flint.arb],
    ratios: list[flint.arb],
    c0: flint.arb,
    T_left: flint.arb,
    T_right: flint.arb,
) -> dict:
    alpha = flint.arb(index) - flint.arb("0.5")
    gamma_base = (alpha + 1).gamma()
    value_total = flint.arb(0)
    derivative_total = flint.arb(0)
    selected_rows: list[dict] = []
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
        moment_constant = (alpha + 1 + flint.arb(degree) / 2).gamma() / gamma_base
        value_exponent = flint.arb(3) - flint.arb(degree) / 2
        derivative_exponent = flint.arb(2) - flint.arb(degree) / 2
        value_T_factor = max_power_on_interval(value_exponent, T_left, T_right)
        derivative_T_factor = max_power_on_interval(
            derivative_exponent, T_left, T_right
        )
        value_contribution = abs_upper(value_coefficient) * moment_constant * value_T_factor
        derivative_contribution = (
            abs_upper(derivative_coefficient) * moment_constant * derivative_T_factor
        )
        value_total += value_contribution
        derivative_total += derivative_contribution
        if degree in SELECTED_DEGREES:
            selected_rows.append(
                {
                    "degree": degree,
                    "value_core_coefficient_ball": arb_text(value_coefficient, 50),
                    "derivative_core_coefficient_ball": arb_text(
                        derivative_coefficient, 50
                    ),
                    "subtracted_even_ratio_ball": (
                        None if subtracted_ratio is None else arb_text(subtracted_ratio, 50)
                    ),
                    "full_gamma_moment_constant_ball": arb_text(moment_constant, 50),
                    "value_scaled_contribution_upper": arb_upper_text(
                        value_contribution, 50
                    ),
                    "derivative_scaled_contribution_upper": arb_upper_text(
                        derivative_contribution, 50
                    ),
                    "proof_boundary": (
                        "Selected absolute full-Gamma moment majorant row; the generator sums every degree 0 through 384."
                    ),
                }
            )
    return {
        "value": value_total,
        "derivative": derivative_total,
        "selected_rows": selected_rows,
    }


def cauchy_scaled_majorants(
    disk: dict[str, flint.arb | bool], c0: flint.arb
) -> dict[str, flint.arb]:
    radius = disk["radius"]
    majorant = disk["majorant"]
    assert isinstance(radius, flint.arb)
    assert isinstance(majorant, flint.arb)
    x_cap = flint.arb(DEFAULT_X_CAP_NUMERATOR) / DEFAULT_X_CAP_DENOMINATOR
    q = x_cap / radius
    T_transition = flint.arb(DEFAULT_T_TRANSITION_NUMERATOR) / DEFAULT_T_TRANSITION_DENOMINATOR
    c0_lower = flint.arb(c0.lower())
    first_tail_degree = DEFAULT_TAYLOR_DEGREE + 1
    value_sup = majorant * q**first_tail_degree / (c0_lower * (1 - q))
    derivative_sup = (
        majorant
        * q**first_tail_degree
        * (flint.arb(first_tail_degree) - flint.arb(DEFAULT_TAYLOR_DEGREE) * q)
        / (2 * c0_lower * (1 - q) ** 2)
    )
    return {
        "q": q,
        "value_sup": value_sup,
        "derivative_sup": derivative_sup,
        "value_scaled": value_sup * T_transition**3,
        "derivative_scaled": derivative_sup * T_transition**2,
    }


def gamma_upper_ratio(shape: flint.arb, split_y: flint.arb, gamma_base: flint.arb) -> flint.arb:
    return split_y.gamma_upper(shape) / gamma_base


def far_regime_bound(
    index: int,
    ratios: list[flint.arb],
    x_majorant_point: flint.arb,
    split_y_lower: flint.arb,
    u_upper: flint.arb,
    value_scale_upper: flint.arb,
    derivative_scale_upper: flint.arb,
) -> dict:
    alpha = flint.arb(index) - flint.arb("0.5")
    gamma_base = (alpha + 1).gamma()
    phi = build_phi_bounds(
        x_majorant_point, DEFAULT_PHI_TERM_COUNT, flint.arb(DEFAULT_C0_LOWER)
    )
    tail_mass = gamma_upper_ratio(alpha + 1, split_y_lower, gamma_base)
    value_unscaled = phi["phi_over_c0_abs_bound_at_split"] * tail_mass
    derivative_unscaled = phi["x_phip_over_2c0_abs_bound_at_split"] * tail_mass
    value_polynomial = flint.arb(0)
    derivative_polynomial = flint.arb(0)
    for j in range(DEFAULT_POLYNOMIAL_M + 1):
        tail_moment = gamma_upper_ratio(alpha + j + 1, split_y_lower, gamma_base)
        ratio_abs = abs_upper(ratios[j])
        term = ratio_abs * u_upper**j * tail_moment
        value_polynomial += term
        if j:
            derivative_polynomial += flint.arb(j) * term
    value_unscaled += value_polynomial
    derivative_unscaled += derivative_polynomial
    return {
        "tail_mass": tail_mass,
        "value_phi_unscaled": phi["phi_over_c0_abs_bound_at_split"] * tail_mass,
        "derivative_phi_unscaled": (
            phi["x_phip_over_2c0_abs_bound_at_split"] * tail_mass
        ),
        "value_polynomial_unscaled": value_polynomial,
        "derivative_polynomial_unscaled": derivative_polynomial,
        "value_scaled": value_unscaled * value_scale_upper,
        "derivative_scaled": derivative_unscaled * derivative_scale_upper,
    }


def uniform_far_majorants(index: int, ratios: list[flint.arb]) -> dict:
    T_left = flint.arb(DEFAULT_T_LEFT)
    T_transition = flint.arb(DEFAULT_T_TRANSITION_NUMERATOR) / DEFAULT_T_TRANSITION_DENOMINATOR
    T_right = flint.arb(DEFAULT_T_RIGHT)
    x_cap = flint.arb(DEFAULT_X_CAP_NUMERATOR) / DEFAULT_X_CAP_DENOMINATOR
    split_one_lower = T_left * flint.arb(64) / 625
    regime_one = far_regime_bound(
        index,
        ratios,
        x_cap,
        split_one_lower,
        1 / T_left,
        T_transition**3,
        T_transition**2,
    )
    split_two = flint.arb(DEFAULT_FIXED_SPLIT_Y)
    x_min = (split_two / T_right).sqrt()
    regime_two = far_regime_bound(
        index,
        ratios,
        x_min,
        split_two,
        1 / T_transition,
        T_right**3,
        T_right**2,
    )
    pi = flint.arb.pi()
    exp4_x_min = (4 * x_min).exp()
    value_margin = 4 * pi * exp4_x_min - 9
    derivative_margin = 4 * pi * exp4_x_min - 13 - 1 / x_min
    if not bool(value_margin > 0 and derivative_margin > 0):
        raise RuntimeError("finite-Phi real-tail monotonicity failed at the smallest split x")
    return {
        "regime_one": regime_one,
        "regime_two": regime_two,
        "value": max(regime_one["value_scaled"], regime_two["value_scaled"]),
        "derivative": max(
            regime_one["derivative_scaled"], regime_two["derivative_scaled"]
        ),
        "x_min": x_min,
        "value_monotonicity_margin": value_margin,
        "derivative_monotonicity_margin": derivative_margin,
    }


def normalization_scaled_majorants(
    disk: dict[str, flint.arb | bool], c0: flint.arb, n_tail: dict[str, flint.arb | bool]
) -> dict[str, flint.arb]:
    majorant = disk["majorant"]
    radius = disk["radius"]
    assert isinstance(majorant, flint.arb)
    assert isinstance(radius, flint.arb)
    x_cap = flint.arb(DEFAULT_X_CAP_NUMERATOR) / DEFAULT_X_CAP_DENOMINATOR
    q = x_cap / radius
    c0_lower = flint.arb(c0.lower())
    c0_tail = n_tail["c0_bound"]
    value_tail = n_tail["value_bound"]
    derivative_tail = n_tail["derivative_bound"]
    assert isinstance(c0_tail, flint.arb)
    assert isinstance(value_tail, flint.arb)
    assert isinstance(derivative_tail, flint.arb)
    full_c0_lower = c0_lower - c0_tail
    finite_phi_abs_expectation = majorant
    finite_derivative_abs_expectation = majorant * q / (2 * (1 - q) ** 2)
    value = (
        value_tail / full_c0_lower
        + finite_phi_abs_expectation * c0_tail / (c0_lower * full_c0_lower)
    )
    derivative = (
        derivative_tail / (2 * full_c0_lower)
        + finite_derivative_abs_expectation
        * c0_tail
        / (c0_lower * full_c0_lower)
    )
    T_right = flint.arb(DEFAULT_T_RIGHT)
    return {
        "value_unscaled": value,
        "derivative_unscaled": derivative,
        "value_scaled": value * T_right**3,
        "derivative_scaled": derivative * T_right**2,
        "finite_phi_abs_expectation_bound": finite_phi_abs_expectation,
        "finite_derivative_abs_expectation_bound": finite_derivative_abs_expectation,
        "full_c0_lower": full_c0_lower,
    }


def build_uniform_rows(
    coefficients: list[flint.arb],
    ratios: list[flint.arb],
    c0: flint.arb,
    disk: dict[str, flint.arb | bool],
    n_tail: dict[str, flint.arb | bool],
) -> tuple[list[dict], dict]:
    T_left = flint.arb(DEFAULT_T_LEFT)
    T_right = flint.arb(DEFAULT_T_RIGHT)
    cauchy = cauchy_scaled_majorants(disk, c0)
    normalization = normalization_scaled_majorants(disk, c0, n_tail)
    A = flint.arb(VALUE_BUDGET)
    B = flint.arb(DERIVATIVE_BUDGET)
    rows = []
    for index in DEFAULT_INDICES:
        compact = compact_polynomial_majorant_rows(
            index, coefficients, ratios, c0, T_left, T_right
        )
        far = uniform_far_majorants(index, ratios)
        value_total = (
            compact["value"]
            + cauchy["value_scaled"]
            + far["value"]
            + normalization["value_scaled"]
        )
        derivative_total = (
            compact["derivative"]
            + cauchy["derivative_scaled"]
            + far["derivative"]
            + normalization["derivative_scaled"]
        )
        value_fraction = value_total / A
        derivative_fraction = derivative_total / B
        certified = bool(value_total < A and derivative_total < B)
        if not certified:
            raise RuntimeError(f"F_{index} did not fit the finite-collar residual budgets")
        rows.append(
            {
                "index": index,
                "row_label": f"F_{index}",
                "T_interval": "1156<=T<=10000",
                "compact_polynomial_value_scaled_upper": arb_upper_text(
                    compact["value"], 70
                ),
                "compact_polynomial_derivative_scaled_upper": arb_upper_text(
                    compact["derivative"], 70
                ),
                "selected_compact_majorant_rows": compact["selected_rows"],
                "selected_compact_majorant_row_count": len(compact["selected_rows"]),
                "cauchy_value_scaled_upper": arb_upper_text(
                    cauchy["value_scaled"], 70
                ),
                "cauchy_derivative_scaled_upper": arb_upper_text(
                    cauchy["derivative_scaled"], 70
                ),
                "regime_one_far_value_scaled_upper": arb_upper_text(
                    far["regime_one"]["value_scaled"], 70
                ),
                "regime_one_far_derivative_scaled_upper": arb_upper_text(
                    far["regime_one"]["derivative_scaled"], 70
                ),
                "regime_two_far_value_scaled_upper": arb_upper_text(
                    far["regime_two"]["value_scaled"], 70
                ),
                "regime_two_far_derivative_scaled_upper": arb_upper_text(
                    far["regime_two"]["derivative_scaled"], 70
                ),
                "uniform_far_value_scaled_upper": arb_upper_text(far["value"], 70),
                "uniform_far_derivative_scaled_upper": arb_upper_text(
                    far["derivative"], 70
                ),
                "normalization_value_scaled_upper": arb_upper_text(
                    normalization["value_scaled"], 70
                ),
                "normalization_derivative_scaled_upper": arb_upper_text(
                    normalization["derivative_scaled"], 70
                ),
                "value_residual_scaled_uniform_upper": arb_upper_text(value_total, 70),
                "derivative_residual_scaled_uniform_upper": arb_upper_text(
                    derivative_total, 70
                ),
                "value_fraction_of_budget_upper": arb_upper_text(value_fraction, 70),
                "derivative_fraction_of_budget_upper": arb_upper_text(
                    derivative_fraction, 70
                ),
                "value_budget_certified_uniformly": bool(value_total < A),
                "derivative_budget_certified_uniformly": bool(derivative_total < B),
                "uniform_residual_budget_certified": certified,
                "proof_boundary": (
                    "Uniform residual majorant on 1156<=T<=10000 for one F_i only; not T>10000."
                ),
            }
        )
    common = {
        "cauchy": cauchy,
        "normalization": normalization,
    }
    return rows, common


def build_diagnostics(recorded_grid: dict) -> dict:
    configure_precision()
    if recorded_grid["summary"]["all_5_recorded_T_stencil_systems_certified"] is not True:
        raise RuntimeError("recorded-grid source is not certified")
    coefficients = finite_phi_series(DEFAULT_PHI_TERM_COUNT, DEFAULT_TAYLOR_DEGREE)
    _ratio_rows, ratios = build_ratio_rows(
        2 * (DEFAULT_POLYNOMIAL_M + 1), DEFAULT_RATIO_CUTOFF_N
    )
    c0 = coefficients[0]
    disk = complex_disk_majorant()
    n_tail = global_n_tail_bounds()
    uniform_rows, common = build_uniform_rows(coefficients, ratios, c0, disk, n_tail)
    collar = finite_collar_data()
    perturbations = perturbation_ledger(collar)
    worst_value = max(
        uniform_rows, key=lambda row: flint.arb(row["value_fraction_of_budget_upper"])
    )
    worst_derivative = max(
        uniform_rows,
        key=lambda row: flint.arb(row["derivative_fraction_of_budget_upper"]),
    )
    x_cap = flint.arb(DEFAULT_X_CAP_NUMERATOR) / DEFAULT_X_CAP_DENOMINATOR
    T_transition = flint.arb(DEFAULT_T_TRANSITION_NUMERATOR) / DEFAULT_T_TRANSITION_DENOMINATOR
    return {
        "parameters": {
            "T_interval": "1156<=T<=10000",
            "T_left": DEFAULT_T_LEFT,
            "T_transition": "15625/8",
            "T_right": DEFAULT_T_RIGHT,
            "indices": list(DEFAULT_INDICES),
            "polynomial_M": DEFAULT_POLYNOMIAL_M,
            "taylor_degree": DEFAULT_TAYLOR_DEGREE,
            "phi_term_count": DEFAULT_PHI_TERM_COUNT,
            "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
            "disk_radius": DEFAULT_DISK_RADIUS,
            "x_cap": "8/25",
            "fixed_split_y": DEFAULT_FIXED_SPLIT_Y,
            "precision_bits": DEFAULT_PRECISION_BITS,
            "value_budget_A": VALUE_BUDGET,
            "derivative_budget_B": DERIVATIVE_BUDGET,
        },
        "split_rule": (
            "For 1156<=T<=15625/8 use x_*=8/25 and y_*=(64/625)T; "
            "for 15625/8<=T<=10000 use y_*=200 and x_*=sqrt(200/T)."
        ),
        "T_transition_ball": arb_text(T_transition, 50),
        "x_cap_ball": arb_text(x_cap, 50),
        "disk_radius_admissible": disk["admissible"],
        "finite_phi_disk_majorant_upper": arb_upper_text(disk["majorant"], 70),
        "cauchy_q_upper": arb_upper_text(common["cauchy"]["q"], 70),
        "cauchy_value_scaled_uniform_upper": arb_upper_text(
            common["cauchy"]["value_scaled"], 70
        ),
        "cauchy_derivative_scaled_uniform_upper": arb_upper_text(
            common["cauchy"]["derivative_scaled"], 70
        ),
        "normalization_value_scaled_uniform_upper": arb_upper_text(
            common["normalization"]["value_scaled"], 70
        ),
        "normalization_derivative_scaled_uniform_upper": arb_upper_text(
            common["normalization"]["derivative_scaled"], 70
        ),
        "global_n_tail_extension_certified": n_tail["global_extension_certified"],
        "uniform_residual_rows": uniform_rows,
        "uniform_residual_row_count": len(uniform_rows),
        "maximum_value_fraction_of_budget_upper": worst_value[
            "value_fraction_of_budget_upper"
        ],
        "maximum_value_fraction_location": worst_value["row_label"],
        "maximum_derivative_fraction_of_budget_upper": worst_derivative[
            "derivative_fraction_of_budget_upper"
        ],
        "maximum_derivative_fraction_location": worst_derivative["row_label"],
        "all_four_value_budgets_certified_uniformly": all(
            row["value_budget_certified_uniformly"] for row in uniform_rows
        ),
        "all_four_derivative_budgets_certified_uniformly": all(
            row["derivative_budget_certified_uniformly"] for row in uniform_rows
        ),
        "all_four_residual_budgets_certified_uniformly": all(
            row["uniform_residual_budget_certified"] for row in uniform_rows
        ),
        "perturbation_rows": perturbations["rows"],
        "all_retained_perturbation_margins_positive": perturbations[
            "all_retained_margins_positive"
        ],
        "finite_collar_segment_stencil_system_certified": bool(
            perturbations["all_retained_margins_positive"]
            and all(row["uniform_residual_budget_certified"] for row in uniform_rows)
        ),
        "remaining_finite_segment_stencil_sources": [],
        "remaining_obligations": [
            "extend the residual majorants from T=10000 to every T>10000",
            "extend fixed k=22 to the complete required k range",
            "complete cone invariance and the sign-regularity/Newman bridges",
        ],
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgfcssc_01_two_regime_split",
            role="exact_reduction",
            readiness="available_exact",
            claim=(
                "Split the finite T segment at T=15625/8: use fixed x_*=8/25 below the transition and fixed y_*=200 above it."
            ),
            formula=diagnostics["split_rule"],
            source_artifacts=[paths["compact_note"], paths["far_tail_note"]],
            proof_boundary="Exact split geometry on 1156<=T<=10000 only.",
        ),
        MatrixRow(
            id="nlrgfcssc_02_compact_polynomial_majorant",
            role="analytic_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "Absolute full-Gamma moments bound every degree-0..384 compact Taylor contribution uniformly in T after the value T^3 and derivative T^2 scalings."
            ),
            diagnostics={
                "uniform_residual_rows": diagnostics["uniform_residual_row_count"],
                "taylor_degree": diagnostics["parameters"]["taylor_degree"],
            },
            source_artifacts=[paths["all_row_note"], paths["compact_note"]],
            proof_boundary="Compact polynomial majorant for four indices on the finite T segment only.",
        ),
        MatrixRow(
            id="nlrgfcssc_03_uniform_cauchy_remainder",
            role="analytic_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The common |z|<=0.38 Cauchy disk gives scaled degree-384 remainders uniformly across both T regimes."
            ),
            diagnostics={
                "cauchy_q_upper": diagnostics["cauchy_q_upper"],
                "value_scaled_uniform_upper": diagnostics[
                    "cauchy_value_scaled_uniform_upper"
                ],
                "derivative_scaled_uniform_upper": diagnostics[
                    "cauchy_derivative_scaled_uniform_upper"
                ],
            },
            source_artifacts=[paths["compact_note"], paths["formal_tail_obstruction"]],
            proof_boundary="Convergent compact Cauchy remainder only; not formal full-moment tail summation.",
        ),
        MatrixRow(
            id="nlrgfcssc_04_two_regime_real_tail",
            role="analytic_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "Endpoint-monotone finite-Phi majorants and upper incomplete-Gamma moments bound the real tails uniformly in each T regime."
            ),
            diagnostics={
                "split_rule": diagnostics["split_rule"],
                "uniform_residual_rows": diagnostics["uniform_residual_row_count"],
            },
            source_artifacts=[paths["far_tail_note"], paths["all_row_note"]],
            proof_boundary="Finite n<=30 real-tail majorants on 1156<=T<=10000 only.",
        ),
        MatrixRow(
            id="nlrgfcssc_05_global_normalization_tail",
            role="normalization_certificate",
            readiness="available_normalization_certificate",
            claim=(
                "The global n>=31 numerator and Phi(0) correction remains negligible after worst-case T^3 and T^2 scaling on the finite segment."
            ),
            diagnostics={
                "global_n_tail_extension_certified": diagnostics[
                    "global_n_tail_extension_certified"
                ],
                "value_scaled_uniform_upper": diagnostics[
                    "normalization_value_scaled_uniform_upper"
                ],
                "derivative_scaled_uniform_upper": diagnostics[
                    "normalization_derivative_scaled_uniform_upper"
                ],
            },
            source_artifacts=[paths["all_row_note"], paths["far_tail_note"]],
            proof_boundary="Global omitted-n normalization correction only.",
        ),
        MatrixRow(
            id="nlrgfcssc_06_uniform_residual_budget_certificate",
            role="real_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "For every real 1156<=T<=10000 and i=21..24, |R_i(u)|<=(1/2)u^3 and |R_i'(u)|<=(9/1000)u."
            ),
            diagnostics={
                "uniform_residual_rows": diagnostics["uniform_residual_rows"],
                "maximum_value_fraction_of_budget_upper": diagnostics[
                    "maximum_value_fraction_of_budget_upper"
                ],
                "maximum_derivative_fraction_of_budget_upper": diagnostics[
                    "maximum_derivative_fraction_of_budget_upper"
                ],
                "all_four_residual_budgets_certified_uniformly": diagnostics[
                    "all_four_residual_budgets_certified_uniformly"
                ],
            },
            source_artifacts=[paths["degree40_budget_note"], paths["recorded_grid_note"]],
            proof_boundary="Four-index residual theorem on the bounded real-T segment only; not T>10000.",
        ),
        MatrixRow(
            id="nlrgfcssc_07_finite_collar_segment_stencil_certificate",
            role="real_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The exact rational residual budgets and Arb perturbation ledger certify the fixed-k=22 normalizers and three stencil inequalities for every 1156<=T<=10000."
            ),
            diagnostics={
                "perturbation_rows": diagnostics["perturbation_rows"],
                "finite_collar_segment_stencil_system_certified": diagnostics[
                    "finite_collar_segment_stencil_system_certified"
                ],
                "remaining_finite_segment_stencil_sources": diagnostics[
                    "remaining_finite_segment_stencil_sources"
                ],
            },
            source_artifacts=[paths["degree40_ladder_note"], paths["recorded_grid_note"]],
            proof_boundary="Fixed k=22 on 1156<=T<=10000 only; not the unbounded collar or all k.",
        ),
        MatrixRow(
            id="nlrgfcssc_08_T_gt_10000_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The bounded segment certificate proves the same residual and stencil bounds for every T>10000.",
            gap=(
                "The crude normalization and fixed-degree Cauchy bounds were scaled only through T=10000; an asymptotic or adaptive split is still required on the unbounded ray."
            ),
            source_artifacts=[paths["uniform_remainder_target"], paths["cone_entry_target"]],
            proof_boundary="Rejected bounded-to-unbounded T promotion only.",
        ),
        MatrixRow(
            id="nlrgfcssc_09_all_k_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The fixed-k=22 segment certificate proves cone entry for all k.",
            gap="Only F_21 through F_24 and one k-window are controlled.",
            source_artifacts=[paths["cone_entry_target"], paths["dependency_graph"]],
            proof_boundary="Rejected fixed-k-to-all-k promotion only.",
        ),
        MatrixRow(
            id="nlrgfcssc_10_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "Promotion requires a certified T>10000 tail and all required k windows; neither finite endpoints nor RH-assuming inputs may replace them."
            ),
            source_artifacts=[paths["uniform_remainder_target"], paths["dependency_graph"]],
            proof_boundary="Proof-hygiene gate only; not cone entry, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(recorded_grid_path: Path) -> dict:
    recorded_grid = load_json(recorded_grid_path)
    paths = source_paths(recorded_grid_path)
    diagnostics = build_diagnostics(recorded_grid)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "uniform_residual_rows": diagnostics["uniform_residual_row_count"],
        "T_regimes": 2,
        "positive_retained_perturbation_margins": sum(
            1
            for row in diagnostics["perturbation_rows"]
            if row["certified_positive_after_perturbation"]
        ),
        "maximum_value_fraction_of_budget_upper": diagnostics[
            "maximum_value_fraction_of_budget_upper"
        ],
        "maximum_derivative_fraction_of_budget_upper": diagnostics[
            "maximum_derivative_fraction_of_budget_upper"
        ],
        "finite_collar_segment_stencil_system_certified": diagnostics[
            "finite_collar_segment_stencil_system_certified"
        ],
        "remaining_finite_segment_stencil_source_count": len(
            diagnostics["remaining_finite_segment_stencil_sources"]
        ),
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "A two-regime exact-moment/Cauchy and incomplete-Gamma majorant proves the rational residual "
            "budgets uniformly for every real 1156<=T<=10000 and i=21..24. Composed with the Arb "
            "degree-40 perturbation ledger, this certifies the complete fixed-k=22 stencil system on the "
            "whole bounded collar segment. The unbounded ray T>10000, remaining k coverage, cone entry, "
            "RH, and Lambda <= 0 remain open."
        ),
    }
    return {
        "kind": (
            "jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate"
        ),
        "date": "2026-07-10",
        "status": "bounded real-T fixed-k stencil certificate",
        "source_recorded_grid_stencil_certificate": paths["recorded_grid_note"],
        "source_recorded_grid_stencil_json": paths["recorded_grid_json"],
        "source_all_row_expectation_certificate": paths["all_row_note"],
        "source_compact_certificate": paths["compact_note"],
        "source_far_tail_certificate": paths["far_tail_note"],
        "source_degree40_budget": paths["degree40_budget_note"],
        "source_degree40_ladder": paths["degree40_ladder_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_cone_entry_target": paths["cone_entry_target"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py"
        ),
        "proof_boundary": (
            "Uniform fixed-k=22 residual and stencil certificate on the bounded real interval "
            "1156<=T<=10000 only. It does not cover T>10000, does not cover all k, does not prove the "
            "full cone-entry or sign-regularity theorem, and does not prove RH or Lambda <= 0."
        ),
        "diagnostics": diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "The two T regimes cover every real T in [1156,10000] with a shared endpoint.",
            "No floating quadrature node, weight, or threshold is a proof input.",
            "The compact Cauchy tail is convergent and does not use the obstructed formal full-moment tail.",
            "The bounded T segment is not promoted to T>10000.",
            "Fixed k=22 is not promoted to all k.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian finite-collar-segment stencil certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, {summary['uniform_residual_rows']} uniform residual rows, "
        f"{summary['T_regimes']} T regimes, "
        f"{summary['positive_retained_perturbation_margins']} positive perturbation margins, "
        f"{summary['remaining_finite_segment_stencil_source_count']} open finite-segment stencil sources, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    diagnostics = artifact["diagnostics"]
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Finite-Collar-Segment Stencil Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: bounded real-T fixed-k stencil certificate. This is not a proof",
        "of the T>10000 ray, an all-k cone-entry theorem, RH, or `Lambda <= 0`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        artifact["generator"],
        artifact["checker"],
        "```",
        "",
        "## Result",
        "",
        summary["main_finding"],
        "",
        "The exact split is",
        "",
        "```text",
        diagnostics["split_rule"],
        "```",
        "",
        "Both regimes keep the compact core inside `x<=8/25<0.38`. The",
        "compact Taylor polynomial is bounded by absolute full-Gamma moments,",
        "the degree-384 remainder by a convergent Cauchy estimate, and the real",
        "tails by endpoint-monotone finite-Phi majorants and upper",
        "incomplete-Gamma moments. The global `n>=31` normalization correction",
        "is then scaled at the worst finite endpoint `T=10000`.",
        "",
        "## Uniform Residual Bounds",
        "",
        "| Index | value scaled upper | fraction of 1/2 | derivative scaled upper | fraction of 9/1000 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in diagnostics["uniform_residual_rows"]:
        lines.append(
            f"| `F_{row['index']}` | `{row['value_residual_scaled_uniform_upper']}` | "
            f"`{row['value_fraction_of_budget_upper']}` | "
            f"`{row['derivative_residual_scaled_uniform_upper']}` | "
            f"`{row['derivative_fraction_of_budget_upper']}` |"
        )
    lines.extend(
        [
            "",
            "All four rows therefore satisfy, for every real `1156<=T<=10000`,",
            "",
            "```text",
            "|R_i(1/T)|  <= (1/2) T^(-3)",
            "|R_i'(1/T)| <= (9/1000) T^(-1).",
            "```",
            "",
            "## Stencil Composition",
            "",
            "| Target | finite margin lower | perturbation upper | retained lower |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in diagnostics["perturbation_rows"]:
        lines.append(
            f"| `{row['name']}` | `{row['finite_margin_lower']}` | "
            f"`{row['perturbation_bound_upper']}` | `{row['retained_margin_lower']}` |"
        )
    lines.extend(
        [
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
            "The numerical-integration and between-grid-T obligations are now closed",
            "for this fixed k on the bounded segment. The next T obligation is the",
            "unbounded ray `T>10000`; the separate all-k obligation also remains.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            result_line(artifact),
            "```",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recorded-grid-json", type=Path, default=DEFAULT_RECORDED_GRID_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = build_artifact(args.recorded_grid_json)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, args.note)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
