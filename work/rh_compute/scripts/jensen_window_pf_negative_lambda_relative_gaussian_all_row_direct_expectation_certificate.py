#!/usr/bin/env python3
"""Certify the complete recorded relative-Gaussian expectation grid directly."""

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

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)
from jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix import (  # noqa: E402
    finite_phi_series,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate import (  # noqa: E402
    abs_upper,
    build_phi_bounds,
    build_tail_moment_bounds,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate import (  # noqa: E402
    global_n_tail_bounds,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FIRST_OMITTED_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json"
)
DEFAULT_FLOATING_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.md"
)

DEFAULT_T_GRID = (1156, 1500, 2000, 5000, 10000)
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_TAYLOR_DEGREE = 384
DEFAULT_DISK_RADIUS = "0.38"
DEFAULT_CORE_X_CAP_SQUARED_NUMERATOR = 64
DEFAULT_CORE_X_CAP_SQUARED_DENOMINATOR = 625
DEFAULT_MAX_SPLIT_Y = 200
DEFAULT_C0_LOWER = "0.44"
DEFAULT_PRECISION_BITS = 4096
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


def parse_arb(text: str) -> flint.arb:
    return flint.arb(text.replace("E", "e"))


def positive_abs_ball(value: flint.arb) -> flint.arb:
    if value.contains(0):
        raise RuntimeError("expected a sign-separated coefficient ball")
    return -value if bool(value < 0) else value


def source_paths(first_omitted_path: Path, floating_grid_path: Path) -> dict[str, str]:
    return {
        "first_omitted_json": first_omitted_path.relative_to(REPO_ROOT).as_posix(),
        "first_omitted_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
        ),
        "floating_grid_json": floating_grid_path.relative_to(REPO_ROOT).as_posix(),
        "floating_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md"
        ),
        "worst_compact_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_x_moment_taylor_certificate.md"
        ),
        "worst_full_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_full_expectation_certificate.md"
        ),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "phi_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md"
        ),
        "coefficient_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md"
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


def split_y(T: int) -> flint.arb:
    disk_limited = (
        flint.arb(T)
        * flint.arb(DEFAULT_CORE_X_CAP_SQUARED_NUMERATOR)
        / flint.arb(DEFAULT_CORE_X_CAP_SQUARED_DENOMINATOR)
    )
    maximum = flint.arb(DEFAULT_MAX_SPLIT_Y)
    return disk_limited if bool(disk_limited < maximum) else maximum


def complex_disk_majorant() -> dict[str, flint.arb | bool]:
    radius = flint.arb(DEFAULT_DISK_RADIUS)
    pi = flint.arb.pi()
    admissible = bool(4 * radius < pi / 2)
    real_exp_lower = (-4 * radius).exp() * (4 * radius).cos()
    if not admissible or real_exp_lower.contains(0) or not bool(real_exp_lower > 0):
        raise RuntimeError("the common complex disk did not certify Re(exp(4z))>0")
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
        "admissible": admissible,
        "real_exp_lower": real_exp_lower,
        "majorant": majorant,
    }


def compact_moment(T: int, index: int, right_y: flint.arb, degree: int) -> flint.arb:
    alpha = flint.arb(index) - flint.arb("0.5")
    shape = alpha + 1 + flint.arb(degree) / 2
    scale = flint.arb(T) ** (flint.arb(degree) / 2)
    return right_y.gamma_lower(shape) / ((alpha + 1).gamma() * scale)


def compact_integrals(
    T: int,
    index: int,
    right_y: flint.arb,
    coefficients: list[flint.arb],
    ratios: list[flint.arb],
    c0: flint.arb,
) -> tuple[flint.arb, flint.arb, list[dict]]:
    value_integral = flint.arb(0)
    derivative_integral = flint.arb(0)
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
        moment = compact_moment(T, index, right_y, degree)
        value_contribution = value_coefficient * moment
        derivative_contribution = derivative_coefficient * moment
        value_integral += value_contribution
        derivative_integral += derivative_contribution
        if degree in SELECTED_DEGREES:
            selected_rows.append(
                {
                    "degree": degree,
                    "phi_taylor_coefficient_ball": arb_text(phi_coefficient, 50),
                    "subtracted_even_ratio_ball": (
                        None if subtracted_ratio is None else arb_text(subtracted_ratio, 50)
                    ),
                    "compact_moment_ball": arb_text(moment, 50),
                    "value_contribution_ball": arb_text(value_contribution, 50),
                    "derivative_contribution_ball": arb_text(derivative_contribution, 50),
                    "proof_boundary": "Selected exact-moment audit row; the generator sums every degree 0 through 384.",
                }
            )
    return value_integral, derivative_integral, selected_rows


def source_denominator_rows(first_omitted: dict) -> dict[tuple[int, int], dict]:
    rows = first_omitted["matrix_rows"][2]["diagnostics"]["denominator_rows"]
    return {(int(row["T"]), int(row["index"])): row for row in rows}


def source_floating_rows(floating_grid: dict) -> dict[tuple[int, int], dict]:
    rows = floating_grid["matrix_rows"][1]["diagnostics"]["grid_rows"]
    return {(int(row["T"]), int(row["index"])): row for row in rows}


def build_grid_row(
    T: int,
    index: int,
    coefficients: list[flint.arb],
    ratios: list[flint.arb],
    c0: flint.arb,
    disk: dict[str, flint.arb | bool],
    n_tail: dict[str, flint.arb | bool],
    denominator_source: dict,
    floating_source: dict,
) -> dict:
    right_y = split_y(T)
    x_right = (right_y / flint.arb(T)).sqrt()
    radius = disk["radius"]
    majorant = disk["majorant"]
    assert isinstance(radius, flint.arb)
    assert isinstance(majorant, flint.arb)
    q = x_right / radius
    if not bool(q < 1):
        raise RuntimeError(f"row T={T}, F_{index} escaped the common disk")

    c0_lower = flint.arb(c0.lower())
    c0_tail = n_tail["c0_bound"]
    assert isinstance(c0_tail, flint.arb)
    full_c0_lower = c0_lower - c0_tail
    if not bool(c0_lower > flint.arb(DEFAULT_C0_LOWER)) or not bool(
        full_c0_lower > flint.arb(DEFAULT_C0_LOWER)
    ):
        raise RuntimeError("finite/full c0 lower bound failed")

    value_integral, derivative_integral, selected_rows = compact_integrals(
        T, index, right_y, coefficients, ratios, c0
    )
    compact_mass = compact_moment(T, index, right_y, 0)
    first_tail_degree = DEFAULT_TAYLOR_DEGREE + 1
    value_tail_sup = majorant * (q**first_tail_degree) / (c0_lower * (1 - q))
    derivative_tail_sup = (
        majorant
        * (q**first_tail_degree)
        * (flint.arb(first_tail_degree) - flint.arb(DEFAULT_TAYLOR_DEGREE) * q)
        / (2 * c0_lower * ((1 - q) ** 2))
    )
    value_taylor_radius = abs_upper(value_tail_sup * compact_mass)
    derivative_taylor_radius = abs_upper(derivative_tail_sup * compact_mass)

    alpha = flint.arb(index) - flint.arb("0.5")
    u = flint.arb(1) / flint.arb(T)
    pi = flint.arb.pi()
    exp4 = (4 * x_right).exp()
    value_monotonicity_margin = 4 * pi * exp4 - 9
    derivative_monotonicity_margin = 4 * pi * exp4 - 13 - 1 / x_right
    if not bool(value_monotonicity_margin > 0) or not bool(
        derivative_monotonicity_margin > 0
    ):
        raise RuntimeError(f"row T={T}, F_{index} real-tail monotonicity failed")
    phi_tail = build_phi_bounds(x_right, DEFAULT_PHI_TERM_COUNT, flint.arb(DEFAULT_C0_LOWER))
    moment_tail = build_tail_moment_bounds(
        ratios, right_y, alpha, u, DEFAULT_POLYNOMIAL_M
    )
    value_far_phi = phi_tail["phi_over_c0_abs_bound_at_split"] * moment_tail["tail_mass"]
    derivative_far_phi = (
        phi_tail["x_phip_over_2c0_abs_bound_at_split"] * moment_tail["tail_mass"]
    )
    value_far_bound = value_far_phi + moment_tail["value_polynomial_tail_bound"]
    derivative_far_bound = derivative_far_phi + moment_tail["derivative_polynomial_tail_bound"]

    compact_phi_abs_integral = majorant * compact_mass
    compact_x_phip_over_2_abs_integral = majorant * q * compact_mass / (2 * (1 - q) ** 2)
    far_phi_abs_integral = phi_tail["phi_abs_bound_at_split"] * moment_tail["tail_mass"]
    far_x_phip_over_2_abs_integral = (
        phi_tail["x_phip_abs_bound_at_split"] * moment_tail["tail_mass"] / 2
    )
    finite_phi_abs_integral = compact_phi_abs_integral + far_phi_abs_integral
    finite_x_phip_over_2_abs_integral = (
        compact_x_phip_over_2_abs_integral + far_x_phip_over_2_abs_integral
    )
    value_numerator_tail = n_tail["value_bound"]
    derivative_numerator_tail = n_tail["derivative_bound"]
    assert isinstance(value_numerator_tail, flint.arb)
    assert isinstance(derivative_numerator_tail, flint.arb)
    value_normalization_radius = (
        value_numerator_tail / full_c0_lower
        + finite_phi_abs_integral * c0_tail / (c0_lower * full_c0_lower)
    )
    derivative_normalization_radius = (
        derivative_numerator_tail / (2 * full_c0_lower)
        + finite_x_phip_over_2_abs_integral * c0_tail / (c0_lower * full_c0_lower)
    )

    value_total_radius = value_taylor_radius + value_far_bound + value_normalization_radius
    derivative_total_radius = (
        derivative_taylor_radius + derivative_far_bound + derivative_normalization_radius
    )
    full_value = value_integral + flint.arb(0, value_total_radius)
    full_derivative = derivative_integral + flint.arb(0, derivative_total_radius)

    first_ratio_abs = positive_abs_ball(ratios[DEFAULT_POLYNOMIAL_M + 1])
    rising = (alpha + 1).rising(DEFAULT_POLYNOMIAL_M + 1)
    first_value = first_ratio_abs * rising * (u ** (DEFAULT_POLYNOMIAL_M + 1))
    first_derivative = flint.arb(DEFAULT_POLYNOMIAL_M + 1) * first_value
    source_scaled_value_ball = parse_arb(denominator_source["value_denominator_ball"])
    source_scaled_derivative_ball = parse_arb(denominator_source["derivative_denominator_ball"])
    source_scaled_value_lower = flint.arb(source_scaled_value_ball.lower())
    source_scaled_derivative_lower = flint.arb(source_scaled_derivative_ball.lower())
    first_value_lower = source_scaled_value_lower * (u**3)
    first_derivative_lower = source_scaled_derivative_lower * (u**2)
    value_ratio = abs_upper(full_value) / first_value_lower
    derivative_ratio = abs_upper(full_derivative) / first_derivative_lower
    value_margin = 1 - value_ratio
    derivative_margin = 1 - derivative_ratio

    denominator_crosscheck = bool(
        first_value_lower > 0
        and first_derivative_lower > 0
        and flint.arb(first_value.lower()) <= flint.arb(source_scaled_value_ball.upper()) * (u**3)
        and first_value_lower <= flint.arb(first_value.upper())
        and flint.arb(first_derivative.lower())
        <= flint.arb(source_scaled_derivative_ball.upper()) * (u**2)
        and first_derivative_lower <= flint.arb(first_derivative.upper())
    )
    floating_value_ratio = parse_arb(floating_source["selected_value_ratio_to_first_omitted"])
    floating_derivative_ratio = parse_arb(
        floating_source["selected_derivative_ratio_to_first_omitted"]
    )
    floating_tolerance = flint.arb("0.0002")
    floating_value_gap = abs_upper(value_ratio - floating_value_ratio)
    floating_derivative_gap = abs_upper(derivative_ratio - floating_derivative_ratio)
    floating_consistent = bool(
        floating_value_gap < floating_tolerance
        and floating_derivative_gap < floating_tolerance
    )
    certified = bool(
        full_value < 0
        and full_derivative < 0
        and value_ratio < 1
        and derivative_ratio < 1
        and denominator_crosscheck
        and floating_consistent
    )
    if not certified:
        raise RuntimeError(
            f"row T={T}, F_{index} did not certify: "
            f"value_negative={bool(full_value < 0)}, "
            f"derivative_negative={bool(full_derivative < 0)}, "
            f"value_ratio_below_one={bool(value_ratio < 1)}, "
            f"derivative_ratio_below_one={bool(derivative_ratio < 1)}, "
            f"denominator_crosscheck={denominator_crosscheck}, "
            f"floating_consistent={floating_consistent}"
        )

    return {
        "T": T,
        "index": index,
        "row_label": f"T={T}, F_{index}",
        "alpha": f"{2 * index - 1}/2",
        "u": f"1/{T}",
        "split_y_ball": arb_text(right_y, 50),
        "split_y_upper": arb_upper_text(right_y, 50),
        "split_x_upper": arb_upper_text(x_right, 50),
        "x_to_disk_radius_ratio_upper": arb_upper_text(q, 50),
        "compact_gamma_mass_ball": arb_text(compact_mass, 50),
        "selected_compact_moment_rows": selected_rows,
        "selected_compact_moment_row_count": len(selected_rows),
        "value_compact_polynomial_integral_ball": arb_text(value_integral, 80),
        "derivative_compact_polynomial_integral_ball": arb_text(derivative_integral, 80),
        "value_taylor_remainder_radius_upper": arb_upper_text(value_taylor_radius, 70),
        "derivative_taylor_remainder_radius_upper": arb_upper_text(
            derivative_taylor_radius, 70
        ),
        "value_real_tail_monotonicity_margin_lower": arb_lower_text(
            value_monotonicity_margin, 50
        ),
        "derivative_real_tail_monotonicity_margin_lower": arb_lower_text(
            derivative_monotonicity_margin, 50
        ),
        "gamma_tail_mass_upper": arb_upper_text(moment_tail["tail_mass"], 70),
        "value_finite_phi_far_tail_bound_upper": arb_upper_text(value_far_phi, 70),
        "derivative_finite_phi_far_tail_bound_upper": arb_upper_text(
            derivative_far_phi, 70
        ),
        "value_polynomial_far_tail_bound_upper": arb_upper_text(
            moment_tail["value_polynomial_tail_bound"], 70
        ),
        "derivative_polynomial_far_tail_bound_upper": arb_upper_text(
            moment_tail["derivative_polynomial_tail_bound"], 70
        ),
        "value_total_far_tail_bound_upper": arb_upper_text(value_far_bound, 70),
        "derivative_total_far_tail_bound_upper": arb_upper_text(
            derivative_far_bound, 70
        ),
        "value_n_tail_normalization_radius_upper": arb_upper_text(
            value_normalization_radius, 70
        ),
        "derivative_n_tail_normalization_radius_upper": arb_upper_text(
            derivative_normalization_radius, 70
        ),
        "value_total_added_radius_upper": arb_upper_text(value_total_radius, 70),
        "derivative_total_added_radius_upper": arb_upper_text(
            derivative_total_radius, 70
        ),
        "full_value_expectation_ball": arb_text(full_value, 80),
        "full_derivative_expectation_ball": arb_text(full_derivative, 80),
        "full_value_certified_negative": bool(full_value < 0),
        "full_derivative_certified_negative": bool(full_derivative < 0),
        "value_first_omitted_expectation_ball": arb_text(first_value, 70),
        "derivative_first_omitted_expectation_ball": arb_text(first_derivative, 70),
        "value_first_omitted_expectation_lower": arb_lower_text(first_value_lower, 70),
        "derivative_first_omitted_expectation_lower": arb_lower_text(
            first_derivative_lower, 70
        ),
        "value_full_ratio_to_first_omitted_upper": arb_upper_text(value_ratio, 70),
        "derivative_full_ratio_to_first_omitted_upper": arb_upper_text(
            derivative_ratio, 70
        ),
        "value_remaining_margin_below_one_lower": arb_lower_text(value_margin, 70),
        "derivative_remaining_margin_below_one_lower": arb_lower_text(
            derivative_margin, 70
        ),
        "source_scaled_value_denominator_lower": denominator_source[
            "value_denominator_lower"
        ],
        "source_scaled_value_denominator_ball": denominator_source[
            "value_denominator_ball"
        ],
        "source_scaled_derivative_denominator_lower": denominator_source[
            "derivative_denominator_lower"
        ],
        "source_scaled_derivative_denominator_ball": denominator_source[
            "derivative_denominator_ball"
        ],
        "source_denominator_lowers_contained_by_computed_balls": denominator_crosscheck,
        "source_floating_value_ratio": floating_source[
            "selected_value_ratio_to_first_omitted"
        ],
        "source_floating_derivative_ratio": floating_source[
            "selected_derivative_ratio_to_first_omitted"
        ],
        "floating_value_ratio_gap_upper": arb_upper_text(floating_value_gap, 50),
        "floating_derivative_ratio_gap_upper": arb_upper_text(
            floating_derivative_gap, 50
        ),
        "floating_diagnostic_consistency_tolerance": "0.0002",
        "floating_diagnostic_agrees_within_tolerance": floating_consistent,
        "row_direct_expectation_certified": certified,
        "quadrature_used_as_proof_input": False,
        "proof_boundary": (
            "One recorded finite-grid expectation row only; no interpolation in T and no uniform collar claim."
        ),
    }


def build_diagnostics(first_omitted: dict, floating_grid: dict) -> dict:
    configure_precision()
    coefficients = finite_phi_series(DEFAULT_PHI_TERM_COUNT, DEFAULT_TAYLOR_DEGREE)
    _ratio_rows, ratios = build_ratio_rows(
        2 * (DEFAULT_POLYNOMIAL_M + 1), DEFAULT_RATIO_CUTOFF_N
    )
    c0 = coefficients[0]
    if c0.contains(0) or not bool(c0 > flint.arb(DEFAULT_C0_LOWER)):
        raise RuntimeError("finite Phi_30(0) did not certify above the c0 floor")
    disk = complex_disk_majorant()
    n_tail = global_n_tail_bounds()
    denominator_rows = source_denominator_rows(first_omitted)
    floating_rows = source_floating_rows(floating_grid)
    expected_keys = {(T, index) for T in DEFAULT_T_GRID for index in DEFAULT_INDICES}
    if set(denominator_rows) != expected_keys or set(floating_rows) != expected_keys:
        raise RuntimeError("source grid keys do not match the required 20-row grid")

    grid_rows = [
        build_grid_row(
            T,
            index,
            coefficients,
            ratios,
            c0,
            disk,
            n_tail,
            denominator_rows[(T, index)],
            floating_rows[(T, index)],
        )
        for T in DEFAULT_T_GRID
        for index in DEFAULT_INDICES
    ]
    value_worst = max(
        grid_rows, key=lambda row: parse_arb(row["value_full_ratio_to_first_omitted_upper"])
    )
    derivative_worst = max(
        grid_rows,
        key=lambda row: parse_arb(row["derivative_full_ratio_to_first_omitted_upper"]),
    )
    minimum_value_margin = min(
        parse_arb(row["value_remaining_margin_below_one_lower"]) for row in grid_rows
    )
    minimum_derivative_margin = min(
        parse_arb(row["derivative_remaining_margin_below_one_lower"]) for row in grid_rows
    )
    return {
        "parameters": {
            "T_grid": list(DEFAULT_T_GRID),
            "indices": list(DEFAULT_INDICES),
            "grid_row_count": len(grid_rows),
            "phi_term_count": DEFAULT_PHI_TERM_COUNT,
            "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
            "polynomial_M": DEFAULT_POLYNOMIAL_M,
            "taylor_degree": DEFAULT_TAYLOR_DEGREE,
            "complex_disk_radius": DEFAULT_DISK_RADIUS,
            "core_x_cap_squared": (
                f"{DEFAULT_CORE_X_CAP_SQUARED_NUMERATOR}/{DEFAULT_CORE_X_CAP_SQUARED_DENOMINATOR}"
            ),
            "core_x_cap": "8/25",
            "maximum_split_y": DEFAULT_MAX_SPLIT_Y,
            "precision_bits": DEFAULT_PRECISION_BITS,
            "compact_moment_formula": (
                "T^(-k/2)*lower_gamma(i+1/2+k/2,split_y)/Gamma(i+1/2)"
            ),
            "first_omitted_value_formula": (
                "abs(r_21)*(i+1/2)_21*T^(-21)"
            ),
            "first_omitted_derivative_formula": (
                "21*abs(r_21)*(i+1/2)_21*T^(-21)"
            ),
        },
        "finite_c0_ball": arb_text(c0, 80),
        "finite_c0_lower": arb_lower_text(flint.arb(c0.lower()), 70),
        "disk_radius_admissible": disk["admissible"],
        "real_exp4_lower_on_disk": arb_lower_text(disk["real_exp_lower"], 70),
        "finite_phi_complex_disk_majorant_upper": arb_upper_text(
            disk["majorant"], 70
        ),
        "global_n_tail_extension_certified": n_tail["global_extension_certified"],
        "value_n_tail_global_bound_upper": arb_upper_text(n_tail["value_bound"], 70),
        "derivative_n_tail_global_bound_upper": arb_upper_text(
            n_tail["derivative_bound"], 70
        ),
        "c0_n_tail_bound_upper": arb_upper_text(n_tail["c0_bound"], 70),
        "grid_rows": grid_rows,
        "value_worst_ratio_location": {
            "T": value_worst["T"],
            "index": value_worst["index"],
        },
        "derivative_worst_ratio_location": {
            "T": derivative_worst["T"],
            "index": derivative_worst["index"],
        },
        "maximum_value_ratio_to_first_omitted_upper": value_worst[
            "value_full_ratio_to_first_omitted_upper"
        ],
        "maximum_derivative_ratio_to_first_omitted_upper": derivative_worst[
            "derivative_full_ratio_to_first_omitted_upper"
        ],
        "minimum_value_margin_below_one_lower": arb_lower_text(minimum_value_margin, 70),
        "minimum_derivative_margin_below_one_lower": arb_lower_text(
            minimum_derivative_margin, 70
        ),
        "all_value_expectations_certified_negative": all(
            row["full_value_certified_negative"] for row in grid_rows
        ),
        "all_derivative_expectations_certified_negative": all(
            row["full_derivative_certified_negative"] for row in grid_rows
        ),
        "all_value_ratios_below_one": all(
            parse_arb(row["value_full_ratio_to_first_omitted_upper"]) < 1
            for row in grid_rows
        ),
        "all_derivative_ratios_below_one": all(
            parse_arb(row["derivative_full_ratio_to_first_omitted_upper"]) < 1
            for row in grid_rows
        ),
        "all_denominator_crosschecks_pass": all(
            row["source_denominator_lowers_contained_by_computed_balls"] for row in grid_rows
        ),
        "all_floating_diagnostics_agree_within_tolerance": all(
            row["floating_diagnostic_agrees_within_tolerance"]
            for row in grid_rows
        ),
        "complete_recorded_grid_expectations_certified": all(
            row["row_direct_expectation_certified"] for row in grid_rows
        ),
        "quadrature_needed_for_recorded_grid_expectations": False,
        "remaining_recorded_grid_integral_sources": [],
        "remaining_obligations": [
            "compose the certified expectation rows with the signed finite-degree stencil margins",
            "replace the five sampled T values by a certified real-T collar argument",
            "complete the downstream sign-regularity and Newman-direction theorem bridges",
        ],
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgardec_01_uniform_split_geometry",
            role="scope_reduction",
            readiness="available_scope_reduction",
            claim=(
                "Every recorded row is split at y=min(200,(64/625)T), so the exact-moment core lies in "
                "0<=x<=8/25 inside the common certified disk |z|<=0.38."
            ),
            diagnostics={
                "core_x_cap": diagnostics["parameters"]["core_x_cap"],
                "complex_disk_radius": diagnostics["parameters"]["complex_disk_radius"],
                "disk_radius_admissible": diagnostics["disk_radius_admissible"],
                "real_exp4_lower_on_disk": diagnostics["real_exp4_lower_on_disk"],
            },
            source_artifacts=[paths["worst_compact_note"], paths["far_tail_note"]],
            proof_boundary="Recorded 20-row split geometry only; not continuous in T.",
        ),
        MatrixRow(
            id="nlrgardec_02_degree384_exact_moment_core",
            role="analytic_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "Arb sums the degree-384 finite-Phi Taylor model against exact lower-incomplete-Gamma "
                "moments on each compact core and bounds the analytic remainder by one common Cauchy disk."
            ),
            diagnostics={
                "taylor_degree": diagnostics["parameters"]["taylor_degree"],
                "compact_moment_formula": diagnostics["parameters"]["compact_moment_formula"],
                "finite_phi_complex_disk_majorant_upper": diagnostics[
                    "finite_phi_complex_disk_majorant_upper"
                ],
            },
            source_artifacts=[paths["worst_compact_note"], paths["coefficient_note"]],
            proof_boundary="Finite n<=30 compact cores for the recorded grid only.",
        ),
        MatrixRow(
            id="nlrgardec_03_rowwise_real_tail_certificate",
            role="analytic_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "Beyond each split, monotone finite-Phi majorants and upper incomplete-Gamma moments "
                "bound both cancellation-reduced real tails."
            ),
            diagnostics={
                "grid_rows": len(diagnostics["grid_rows"]),
                "minimum_split_y": diagnostics["grid_rows"][0]["split_y_upper"],
                "maximum_split_y": DEFAULT_MAX_SPLIT_Y,
            },
            source_artifacts=[paths["far_tail_note"], paths["first_omitted_note"]],
            proof_boundary="Finite n<=30 real tails for 20 rows only.",
        ),
        MatrixRow(
            id="nlrgardec_04_global_n_tail_normalization",
            role="normalization_certificate",
            readiness="available_normalization_certificate",
            claim=(
                "The global n>=31 Phi, x*Phi', and Phi(0) tail bounds are propagated through the "
                "finite/full normalization identity for every row."
            ),
            diagnostics={
                "global_n_tail_extension_certified": diagnostics[
                    "global_n_tail_extension_certified"
                ],
                "value_n_tail_global_bound_upper": diagnostics[
                    "value_n_tail_global_bound_upper"
                ],
                "derivative_n_tail_global_bound_upper": diagnostics[
                    "derivative_n_tail_global_bound_upper"
                ],
                "c0_n_tail_bound_upper": diagnostics["c0_n_tail_bound_upper"],
            },
            source_artifacts=[paths["phi_tail_note"], paths["worst_full_note"]],
            proof_boundary="Global omitted-n correction only; it does not interpolate between T rows.",
        ),
        MatrixRow(
            id="nlrgardec_05_complete_recorded_grid_certificate",
            role="finite_grid_interval_certificate",
            readiness="available_finite_grid_certificate",
            claim=(
                "All 20 recorded full value and derivative expectation balls are strictly negative and "
                "their absolute values are rigorously below the corresponding first omitted terms."
            ),
            diagnostics={
                "grid_rows": len(diagnostics["grid_rows"]),
                "value_worst_ratio_location": diagnostics["value_worst_ratio_location"],
                "derivative_worst_ratio_location": diagnostics[
                    "derivative_worst_ratio_location"
                ],
                "maximum_value_ratio_to_first_omitted_upper": diagnostics[
                    "maximum_value_ratio_to_first_omitted_upper"
                ],
                "maximum_derivative_ratio_to_first_omitted_upper": diagnostics[
                    "maximum_derivative_ratio_to_first_omitted_upper"
                ],
                "minimum_value_margin_below_one_lower": diagnostics[
                    "minimum_value_margin_below_one_lower"
                ],
                "minimum_derivative_margin_below_one_lower": diagnostics[
                    "minimum_derivative_margin_below_one_lower"
                ],
                "complete_recorded_grid_expectations_certified": diagnostics[
                    "complete_recorded_grid_expectations_certified"
                ],
                "quadrature_needed_for_recorded_grid_expectations": diagnostics[
                    "quadrature_needed_for_recorded_grid_expectations"
                ],
            },
            source_artifacts=[paths["first_omitted_json"], paths["worst_full_note"]],
            proof_boundary=(
                "Complete expectation certificate for five T values and four F indices only; not a real-T collar."
            ),
        ),
        MatrixRow(
            id="nlrgardec_06_floating_grid_noninput_check",
            role="independent_consistency_check",
            readiness="available_consistency_check",
            claim=(
                "Every earlier floating quadrature ratio agrees with the corresponding rigorous upper "
                "bound within 2e-4, while quadrature is excluded from the proof inputs."
            ),
            diagnostics={
                "all_floating_diagnostics_agree_within_tolerance": diagnostics[
                    "all_floating_diagnostics_agree_within_tolerance"
                ],
                "quadrature_used_as_proof_input": False,
            },
            source_artifacts=[paths["floating_grid_json"], paths["floating_grid_note"]],
            proof_boundary="Consistency check only; floating quadrature proves nothing here.",
        ),
        MatrixRow(
            id="nlrgardec_07_real_T_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The 20 certified rows prove the expectation inequality for every real T>=1156.",
            gap=(
                "Only five isolated T values are certified. No interval-in-T enclosure or monotonicity theorem "
                "has yet filled the four intervening intervals or the ray beyond T=10000."
            ),
            source_artifacts=[paths["uniform_remainder_target"], paths["dependency_graph"]],
            proof_boundary="Rejected finite-grid-to-collar promotion only.",
        ),
        MatrixRow(
            id="nlrgardec_08_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "Downstream promotion requires a signed-stencil aggregation certificate and a genuine real-T "
                "collar argument; neither may be replaced by sampling or an RH-assuming input."
            ),
            source_artifacts=[
                paths["intervalization_target"],
                paths["uniform_remainder_target"],
                paths["formal_tail_obstruction"],
            ],
            proof_boundary="Proof-hygiene gate only; not RH or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(first_omitted_path: Path, floating_grid_path: Path) -> dict:
    first_omitted = load_json(first_omitted_path)
    floating_grid = load_json(floating_grid_path)
    paths = source_paths(first_omitted_path, floating_grid_path)
    diagnostics = build_diagnostics(first_omitted, floating_grid)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "grid_rows": len(diagnostics["grid_rows"]),
        "certified_value_rows": sum(
            1 for row in diagnostics["grid_rows"] if row["full_value_certified_negative"]
        ),
        "certified_derivative_rows": sum(
            1 for row in diagnostics["grid_rows"] if row["full_derivative_certified_negative"]
        ),
        "below_one_value_rows": sum(
            1
            for row in diagnostics["grid_rows"]
            if parse_arb(row["value_full_ratio_to_first_omitted_upper"]) < 1
        ),
        "below_one_derivative_rows": sum(
            1
            for row in diagnostics["grid_rows"]
            if parse_arb(row["derivative_full_ratio_to_first_omitted_upper"]) < 1
        ),
        "value_worst_ratio_location": diagnostics["value_worst_ratio_location"],
        "derivative_worst_ratio_location": diagnostics["derivative_worst_ratio_location"],
        "maximum_value_ratio_to_first_omitted_upper": diagnostics[
            "maximum_value_ratio_to_first_omitted_upper"
        ],
        "maximum_derivative_ratio_to_first_omitted_upper": diagnostics[
            "maximum_derivative_ratio_to_first_omitted_upper"
        ],
        "minimum_value_margin_below_one_lower": diagnostics[
            "minimum_value_margin_below_one_lower"
        ],
        "minimum_derivative_margin_below_one_lower": diagnostics[
            "minimum_derivative_margin_below_one_lower"
        ],
        "complete_recorded_grid_expectations_certified": diagnostics[
            "complete_recorded_grid_expectations_certified"
        ],
        "quadrature_needed_for_recorded_grid_expectations": diagnostics[
            "quadrature_needed_for_recorded_grid_expectations"
        ],
        "remaining_recorded_grid_integral_source_count": len(
            diagnostics["remaining_recorded_grid_integral_sources"]
        ),
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "A uniform degree-384 exact-moment/Cauchy core plus rowwise incomplete-Gamma real tail "
            "certifies all 20 recorded full relative-Gaussian expectations directly. Both channels are "
            "negative and below one first omitted term in every row; the worst row remains T=10000, F_21. "
            "This closes the recorded-grid integration source, not the real-T collar, signed-stencil "
            "aggregation, sign-regularity bridge, RH, or Lambda <= 0."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate",
        "date": "2026-07-10",
        "status": "complete recorded-grid direct expectation certificate",
        "source_first_omitted_denominator_certificate": paths["first_omitted_note"],
        "source_first_omitted_denominator_json": paths["first_omitted_json"],
        "source_floating_grid_scout": paths["floating_grid_note"],
        "source_floating_grid_json": paths["floating_grid_json"],
        "source_worst_row_compact_certificate": paths["worst_compact_note"],
        "source_worst_row_full_expectation_certificate": paths["worst_full_note"],
        "source_far_tail_certificate": paths["far_tail_note"],
        "source_phi_tail_certificate": paths["phi_tail_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.py"
        ),
        "proof_boundary": (
            "Complete interval certificate for the 20 recorded expectations at T=1156,1500,2000,5000,10000 "
            "and F_21 through F_24 only. It does not certify any T interval, does not yet compose the rows "
            "with signed stencil margins, does not prove the uniform collar or sign-regularity bridge, and "
            "does not prove RH or Lambda <= 0."
        ),
        "diagnostics": diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "All 20 recorded value and derivative expectations are Arb-certified negative.",
            "All 40 certified absolute ratios are below one first omitted term.",
            "Floating quadrature is an independent consistency check and is not a proof input.",
            "The artifact makes no interpolation or monotonicity claim between sampled T values.",
            "Signed-stencil aggregation and the real-T collar remain separate obligations.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian all-row direct expectation certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, {summary['grid_rows']} grid rows, "
        f"{summary['certified_value_rows']} negative values, "
        f"{summary['certified_derivative_rows']} negative derivatives, "
        f"{summary['remaining_recorded_grid_integral_source_count']} open recorded-grid integral sources, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    diagnostics = artifact["diagnostics"]
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian All-Row Direct Expectation Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: complete recorded-grid direct expectation certificate. This is not a proof",
        "of a real-`T` collar theorem, signed-stencil aggregation theorem, RH,",
        "or proof of `Lambda <= 0`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        artifact["generator"],
        artifact["checker"],
        "```",
        "",
        "## Main result",
        "",
        summary["main_finding"],
        "",
        f"- Certified rows: `{summary['grid_rows']}`.",
        f"- Negative value balls: `{summary['certified_value_rows']}/{summary['grid_rows']}`.",
        f"- Negative derivative balls: `{summary['certified_derivative_rows']}/{summary['grid_rows']}`.",
        f"- Value ratios below one: `{summary['below_one_value_rows']}/{summary['grid_rows']}`.",
        f"- Derivative ratios below one: `{summary['below_one_derivative_rows']}/{summary['grid_rows']}`.",
        f"- Maximum value ratio: `{summary['maximum_value_ratio_to_first_omitted_upper']}` at `{summary['value_worst_ratio_location']}`.",
        f"- Maximum derivative ratio: `{summary['maximum_derivative_ratio_to_first_omitted_upper']}` at `{summary['derivative_worst_ratio_location']}`.",
        f"- Minimum value margin below one: `{summary['minimum_value_margin_below_one_lower']}`.",
        f"- Minimum derivative margin below one: `{summary['minimum_derivative_margin_below_one_lower']}`.",
        "- Generalized Gauss-Laguerre quadrature used as a proof input: `False`.",
        "- Remaining recorded-grid integration sources: `0`.",
        "",
        "## Construction",
        "",
        "For each row, split at",
        "",
        "```text",
        "y_* = min(200, (64/625) T),    x_* = sqrt(y_*/T) <= 8/25.",
        "```",
        "",
        "The compact core lies inside `|z|<=0.38<pi/8`. Arb integrates the",
        "degree-384 finite-`Phi_30` Taylor model by exact transformed",
        "lower-incomplete-Gamma moments and bounds the analytic remainder by",
        "Cauchy's estimate. On `y>=y_*`, monotone real majorants and upper",
        "incomplete-Gamma moments bound the finite-`Phi_30` and cancellation",
        "polynomial tails. A global `n>=31` bound then propagates the numerator",
        "and `Phi(0)` normalization corrections.",
        "",
        "The earlier floating quadrature ratios agree with the rigorous row bounds",
        "within `2e-4`, but they are not used to derive any enclosure.",
        "",
        "## Row ledger",
        "",
        "| Row | split y | value ball | value / first | derivative ball | derivative / first |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in diagnostics["grid_rows"]:
        lines.append(
            f"| `{row['row_label']}` | `{row['split_y_upper']}` | "
            f"`{row['full_value_expectation_ball']}` | "
            f"`{row['value_full_ratio_to_first_omitted_upper']}` | "
            f"`{row['full_derivative_expectation_ball']}` | "
            f"`{row['derivative_full_ratio_to_first_omitted_upper']}` |"
        )
    lines.extend(
        [
            "",
            "## Proof boundary",
            "",
            artifact["proof_boundary"],
            "",
            "The next live obligation is no longer interval integration at the 20",
            "recorded rows. It is to compose these rowwise expectation enclosures",
            "with the signed finite-degree stencil margins and then replace the five",
            "sampled `T` values by a genuine interval or monotonicity argument on",
            "the complete real collar `T>=1156`.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.py",
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
    parser.add_argument("--first-omitted-json", type=Path, default=DEFAULT_FIRST_OMITTED_JSON)
    parser.add_argument("--floating-grid-json", type=Path, default=DEFAULT_FLOATING_GRID_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = build_artifact(args.first_omitted_json, args.floating_grid_json)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, args.note)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
