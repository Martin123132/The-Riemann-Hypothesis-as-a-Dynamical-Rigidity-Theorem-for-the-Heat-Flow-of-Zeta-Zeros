#!/usr/bin/env python3
"""Certify the fixed-k=22 stencil system for every real T>=1156."""

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
    arb_lower_text,
    arb_text,
    arb_upper_text,
)
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)
from jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate import (  # noqa: E402
    finite_collar_data,
    perturbation_ledger,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate import (  # noqa: E402
    abs_upper,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_EVENNESS_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.json"
)
DEFAULT_FINITE_SEGMENT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.md"
)

DEFAULT_RAY_START_T = 10000
DEFAULT_COLLAR_START_T = 1156
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_K = 22
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_FIRST_RESIDUAL_DEGREE = 42
DEFAULT_DISK_RADIUS = "0.38"
DEFAULT_SPLIT_X = "0.2"
DEFAULT_FULL_KERNEL_SUM_N = 80
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_C0_LOWER = "0.44"
DEFAULT_PRECISION_BITS = 4096
VALUE_BUDGET = "0.5"
DERIVATIVE_BUDGET = "0.009"


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


def source_paths(evenness_path: Path, finite_segment_path: Path) -> dict[str, str]:
    return {
        "evenness_json": evenness_path.relative_to(REPO_ROOT).as_posix(),
        "evenness_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_full_kernel_evenness_cauchy_lemma.md"
        ),
        "finite_segment_json": finite_segment_path.relative_to(REPO_ROOT).as_posix(),
        "finite_segment_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_finite_collar_segment_stencil_certificate.md"
        ),
        "degree40_ladder_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md"
        ),
        "degree40_budget_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md"
        ),
        "node_c0_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def configure_precision() -> None:
    flint.ctx.prec = DEFAULT_PRECISION_BITS


def geometric_tail(
    first_n: int, degree: int, decay_coefficient: flint.arb, first_term: flint.arb
) -> tuple[flint.arb, flint.arb]:
    n = flint.arb(first_n)
    ratio = ((n + 1) / n) ** degree * (
        -decay_coefficient * (2 * n + 1)
    ).exp()
    if not bool(ratio < 1):
        raise RuntimeError("n-tail geometric ratio did not certify below one")
    return first_term / (1 - ratio), ratio


def full_disk_majorant() -> dict[str, flint.arb | bool]:
    R = flint.arb(DEFAULT_DISK_RADIUS)
    pi = flint.arb.pi()
    real_exp_lower = (-4 * R).exp() * (4 * R).cos()
    admissible = bool(4 * R < pi / 2 and real_exp_lower > 0)
    if not admissible:
        raise RuntimeError("full-kernel disk is inadmissible")
    decay = pi * real_exp_lower
    total = flint.arb(0)
    for n in range(1, DEFAULT_FULL_KERNEL_SUM_N + 1):
        n2 = flint.arb(n * n)
        total += (
            2 * pi * pi * n2 * n2 * (9 * R).exp()
            + 3 * pi * n2 * (5 * R).exp()
        ) * (-decay * n2).exp()
    first_n = DEFAULT_FULL_KERNEL_SUM_N + 1
    n = flint.arb(first_n)
    n2 = n * n
    degree4_first = 2 * pi * pi * n2 * n2 * (9 * R).exp() * (-decay * n2).exp()
    degree2_first = 3 * pi * n2 * (5 * R).exp() * (-decay * n2).exp()
    degree4_tail, degree4_ratio = geometric_tail(
        first_n, 4, decay, degree4_first
    )
    degree2_tail, degree2_ratio = geometric_tail(
        first_n, 2, decay, degree2_first
    )
    total += degree4_tail + degree2_tail
    return {
        "radius": R,
        "real_exp_lower": real_exp_lower,
        "decay_coefficient": decay,
        "majorant": total,
        "degree4_tail": degree4_tail,
        "degree2_tail": degree2_tail,
        "degree4_ratio": degree4_ratio,
        "degree2_ratio": degree2_ratio,
        "admissible": admissible,
    }


def full_real_majorants() -> dict[str, flint.arb | bool]:
    x = flint.arb(DEFAULT_SPLIT_X)
    pi = flint.arb.pi()
    exp4 = (4 * x).exp()
    exp5 = (5 * x).exp()
    exp9 = (9 * x).exp()
    decay = pi * exp4
    value = flint.arb(0)
    derivative = flint.arb(0)
    for n in range(1, DEFAULT_FULL_KERNEL_SUM_N + 1):
        n2 = flint.arb(n * n)
        n4 = n2 * n2
        exponential = (-decay * n2).exp()
        prefactor = 2 * pi * pi * n4 * exp9 + 3 * pi * n2 * exp5
        value += prefactor * exponential
        derivative += x * (
            18 * pi * pi * n4 * exp9
            + 15 * pi * n2 * exp5
            + 4 * pi * n2 * exp4 * prefactor
        ) * exponential
    first_n = DEFAULT_FULL_KERNEL_SUM_N + 1
    n = flint.arb(first_n)
    n2 = n * n
    n4 = n2 * n2
    exponential = (-decay * n2).exp()
    prefactor = 2 * pi * pi * n4 * exp9 + 3 * pi * n2 * exp5
    value_first = prefactor * exponential
    derivative_first = x * (
        18 * pi * pi * n4 * exp9
        + 15 * pi * n2 * exp5
        + 4 * pi * n2 * exp4 * prefactor
    ) * exponential
    value_tail, value_ratio = geometric_tail(first_n, 4, decay, value_first)
    derivative_tail, derivative_ratio = geometric_tail(
        first_n, 6, decay, derivative_first
    )
    value += value_tail
    derivative += derivative_tail
    value_monotonicity_margin = 4 * pi * exp4 - 9
    derivative_monotonicity_margin = 4 * pi * exp4 - 13 - 1 / x
    monotonicity = bool(
        value_monotonicity_margin > 0 and derivative_monotonicity_margin > 0
    )
    if not monotonicity:
        raise RuntimeError("full real-kernel majorants are not monotone after x=0.2")
    return {
        "x": x,
        "decay_coefficient": decay,
        "value": value,
        "derivative": derivative,
        "value_tail": value_tail,
        "derivative_tail": derivative_tail,
        "value_ratio": value_ratio,
        "derivative_ratio": derivative_ratio,
        "value_monotonicity_margin": value_monotonicity_margin,
        "derivative_monotonicity_margin": derivative_monotonicity_margin,
        "monotonicity_certified": monotonicity,
    }


def compact_factored_bounds(index: int, disk: dict[str, flint.arb | bool]) -> dict:
    R = disk["radius"]
    M = disk["majorant"]
    assert isinstance(R, flint.arb)
    assert isinstance(M, flint.arb)
    x = flint.arb(DEFAULT_SPLIT_X)
    q = x / R
    c0_lower = flint.arb(DEFAULT_C0_LOWER)
    first_degree = flint.arb(DEFAULT_FIRST_RESIDUAL_DEGREE)
    value_factor = M / (c0_lower * R**first_degree * (1 - q))
    derivative_factor = M / (2 * c0_lower * R**first_degree) * (
        first_degree / (1 - q) + q / (1 - q) ** 2
    )
    alpha = flint.arb(index) - flint.arb("0.5")
    moment_constant = (alpha + 1).rising(DEFAULT_FIRST_RESIDUAL_DEGREE // 2)
    T0 = flint.arb(DEFAULT_RAY_START_T)
    value_scaled = value_factor * moment_constant * T0 ** (
        3 - DEFAULT_FIRST_RESIDUAL_DEGREE // 2
    )
    derivative_scaled = derivative_factor * moment_constant * T0 ** (
        2 - DEFAULT_FIRST_RESIDUAL_DEGREE // 2
    )
    return {
        "q": q,
        "value_factor": value_factor,
        "derivative_factor": derivative_factor,
        "moment_constant": moment_constant,
        "value_scaled": value_scaled,
        "derivative_scaled": derivative_scaled,
        "value_T_exponent": -18,
        "derivative_T_exponent": -19,
    }


def gamma_upper_ratio(shape: flint.arb, y: flint.arb, gamma_base: flint.arb) -> flint.arb:
    return y.gamma_upper(shape) / gamma_base


def scaled_tail_monotonicity_margin(
    y0: flint.arb, shape: flint.arb, scale_power: int
) -> flint.arb:
    positive_power = max(scale_power, 0)
    return y0 - (shape - 1) - positive_power


def real_tail_bounds(
    index: int, ratios: list[flint.arb], real: dict[str, flint.arb | bool]
) -> dict:
    T0 = flint.arb(DEFAULT_RAY_START_T)
    x = real["x"]
    assert isinstance(x, flint.arb)
    y0 = T0 * x * x
    alpha = flint.arb(index) - flint.arb("0.5")
    gamma_base = (alpha + 1).gamma()
    c0_lower = flint.arb(DEFAULT_C0_LOWER)
    tail_mass = gamma_upper_ratio(alpha + 1, y0, gamma_base)
    value_phi_unscaled = real["value"] / c0_lower * tail_mass
    derivative_phi_unscaled = real["derivative"] / (2 * c0_lower) * tail_mass
    value_polynomial_unscaled = flint.arb(0)
    derivative_polynomial_unscaled = flint.arb(0)
    monotonicity_margins = [
        scaled_tail_monotonicity_margin(y0, alpha + 1, 3),
        scaled_tail_monotonicity_margin(y0, alpha + 1, 2),
    ]
    moment_rows = []
    for j in range(DEFAULT_POLYNOMIAL_M + 1):
        shape = alpha + j + 1
        tail_moment = gamma_upper_ratio(shape, y0, gamma_base)
        ratio_abs = abs_upper(ratios[j])
        value_term = ratio_abs * T0 ** (-j) * tail_moment
        value_polynomial_unscaled += value_term
        value_margin = scaled_tail_monotonicity_margin(y0, shape, 3 - j)
        monotonicity_margins.append(value_margin)
        derivative_term = flint.arb(0)
        derivative_margin = None
        if j:
            derivative_term = flint.arb(j) * value_term
            derivative_polynomial_unscaled += derivative_term
            derivative_margin = scaled_tail_monotonicity_margin(y0, shape, 2 - j)
            monotonicity_margins.append(derivative_margin)
        if j in {0, 1, 2, 10, 20}:
            moment_rows.append(
                {
                    "j": j,
                    "shape": arb_text(shape, 30),
                    "tail_moment_ball": arb_text(tail_moment, 50),
                    "value_scaled_term_upper_at_T0": arb_upper_text(
                        value_term * T0**3, 50
                    ),
                    "derivative_scaled_term_upper_at_T0": arb_upper_text(
                        derivative_term * T0**2, 50
                    ),
                    "value_monotonicity_margin_lower": arb_lower_text(
                        value_margin, 40
                    ),
                    "derivative_monotonicity_margin_lower": (
                        None
                        if derivative_margin is None
                        else arb_lower_text(derivative_margin, 40)
                    ),
                    "proof_boundary": "Selected upper-Gamma ray-tail moment row.",
                }
            )
    minimum_margin = min(monotonicity_margins)
    monotonicity_certified = bool(minimum_margin > 0)
    if not monotonicity_certified:
        raise RuntimeError(f"F_{index} scaled upper-Gamma tail did not decrease on the ray")
    value_scaled = (
        value_phi_unscaled + value_polynomial_unscaled
    ) * T0**3
    derivative_scaled = (
        derivative_phi_unscaled + derivative_polynomial_unscaled
    ) * T0**2
    return {
        "y0": y0,
        "tail_mass": tail_mass,
        "value_phi_scaled": value_phi_unscaled * T0**3,
        "derivative_phi_scaled": derivative_phi_unscaled * T0**2,
        "value_polynomial_scaled": value_polynomial_unscaled * T0**3,
        "derivative_polynomial_scaled": derivative_polynomial_unscaled * T0**2,
        "value_scaled": value_scaled,
        "derivative_scaled": derivative_scaled,
        "minimum_monotonicity_margin": minimum_margin,
        "monotonicity_certified": monotonicity_certified,
        "selected_moment_rows": moment_rows,
        "upper_gamma_hazard_bound": (
            "Gamma(s,y)<=y^s*exp(-y)/(y-s+1), so y^s*exp(-y)/Gamma(s,y)>=y-s+1."
        ),
    }


def build_ray_rows(
    ratios: list[flint.arb],
    disk: dict[str, flint.arb | bool],
    real: dict[str, flint.arb | bool],
) -> list[dict]:
    A = flint.arb(VALUE_BUDGET)
    B = flint.arb(DERIVATIVE_BUDGET)
    rows = []
    for index in DEFAULT_INDICES:
        compact = compact_factored_bounds(index, disk)
        tail = real_tail_bounds(index, ratios, real)
        value_total = compact["value_scaled"] + tail["value_scaled"]
        derivative_total = compact["derivative_scaled"] + tail["derivative_scaled"]
        value_fraction = value_total / A
        derivative_fraction = derivative_total / B
        certified = bool(value_total < A and derivative_total < B)
        if not certified:
            raise RuntimeError(f"F_{index} ray residual budget failed")
        rows.append(
            {
                "index": index,
                "row_label": f"F_{index}",
                "T_interval": "T>=10000",
                "compact_value_scaled_upper_at_T0": arb_upper_text(
                    compact["value_scaled"], 70
                ),
                "compact_derivative_scaled_upper_at_T0": arb_upper_text(
                    compact["derivative_scaled"], 70
                ),
                "compact_value_T_exponent": compact["value_T_exponent"],
                "compact_derivative_T_exponent": compact[
                    "derivative_T_exponent"
                ],
                "real_tail_value_scaled_upper_at_T0": arb_upper_text(
                    tail["value_scaled"], 70
                ),
                "real_tail_derivative_scaled_upper_at_T0": arb_upper_text(
                    tail["derivative_scaled"], 70
                ),
                "upper_gamma_minimum_monotonicity_margin_lower": arb_lower_text(
                    tail["minimum_monotonicity_margin"], 70
                ),
                "selected_real_tail_moment_rows": tail["selected_moment_rows"],
                "selected_real_tail_moment_row_count": len(
                    tail["selected_moment_rows"]
                ),
                "value_residual_scaled_uniform_upper": arb_upper_text(
                    value_total, 70
                ),
                "derivative_residual_scaled_uniform_upper": arb_upper_text(
                    derivative_total, 70
                ),
                "value_fraction_of_budget_upper": arb_upper_text(value_fraction, 70),
                "derivative_fraction_of_budget_upper": arb_upper_text(
                    derivative_fraction, 70
                ),
                "value_budget_certified_on_ray": bool(value_total < A),
                "derivative_budget_certified_on_ray": bool(derivative_total < B),
                "upper_gamma_scaled_tails_decrease_on_ray": tail[
                    "monotonicity_certified"
                ],
                "ray_residual_budget_certified": certified,
                "proof_boundary": "Uniform residual row on T>=10000 for one F_i only.",
            }
        )
    return rows


def build_diagnostics(evenness: dict, finite_segment: dict) -> dict:
    configure_precision()
    if evenness["summary"]["residual_order_42_certified"] is not True:
        raise RuntimeError("source order-42 evenness lemma is not certified")
    if finite_segment["summary"]["finite_collar_segment_stencil_system_certified"] is not True:
        raise RuntimeError("source finite segment is not certified")
    disk = full_disk_majorant()
    real = full_real_majorants()
    _ratio_rows, ratios = build_ratio_rows(
        2 * DEFAULT_POLYNOMIAL_M, DEFAULT_RATIO_CUTOFF_N
    )
    ray_rows = build_ray_rows(ratios, disk, real)
    collar = finite_collar_data()
    perturbations = perturbation_ledger(collar)
    worst_value = max(
        ray_rows, key=lambda row: flint.arb(row["value_fraction_of_budget_upper"])
    )
    worst_derivative = max(
        ray_rows,
        key=lambda row: flint.arb(row["derivative_fraction_of_budget_upper"]),
    )
    ray_stencil = bool(
        perturbations["all_retained_margins_positive"]
        and all(row["ray_residual_budget_certified"] for row in ray_rows)
    )
    full_collar = bool(
        ray_stencil
        and finite_segment["summary"][
            "finite_collar_segment_stencil_system_certified"
        ]
    )
    return {
        "parameters": {
            "full_T_interval": "T>=1156",
            "finite_segment": "1156<=T<=10000",
            "unbounded_ray": "T>=10000",
            "ray_start_T": DEFAULT_RAY_START_T,
            "indices": list(DEFAULT_INDICES),
            "k": DEFAULT_K,
            "polynomial_M": DEFAULT_POLYNOMIAL_M,
            "first_residual_degree": DEFAULT_FIRST_RESIDUAL_DEGREE,
            "disk_radius": DEFAULT_DISK_RADIUS,
            "split_x": DEFAULT_SPLIT_X,
            "split_y_at_ray_start": "400",
            "full_kernel_sum_n": DEFAULT_FULL_KERNEL_SUM_N,
            "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
            "precision_bits": DEFAULT_PRECISION_BITS,
            "value_budget_A": VALUE_BUDGET,
            "derivative_budget_B": DERIVATIVE_BUDGET,
        },
        "source_full_kernel_evenness_certified": evenness["summary"][
            "full_kernel_evenness_certified"
        ],
        "source_order42_residual_zero_certified": evenness["summary"][
            "residual_order_42_certified"
        ],
        "source_finite_segment_certified": finite_segment["summary"][
            "finite_collar_segment_stencil_system_certified"
        ],
        "full_disk_admissible": disk["admissible"],
        "full_disk_majorant_upper": arb_upper_text(disk["majorant"], 70),
        "full_disk_degree4_tail_upper": arb_upper_text(disk["degree4_tail"], 70),
        "full_disk_degree2_tail_upper": arb_upper_text(disk["degree2_tail"], 70),
        "full_disk_degree4_ratio_upper": arb_upper_text(disk["degree4_ratio"], 70),
        "full_disk_degree2_ratio_upper": arb_upper_text(disk["degree2_ratio"], 70),
        "full_real_value_majorant_upper_at_split": arb_upper_text(
            real["value"], 70
        ),
        "full_real_derivative_majorant_upper_at_split": arb_upper_text(
            real["derivative"], 70
        ),
        "full_real_value_n_tail_upper": arb_upper_text(real["value_tail"], 70),
        "full_real_derivative_n_tail_upper": arb_upper_text(
            real["derivative_tail"], 70
        ),
        "full_real_value_monotonicity_margin_lower": arb_lower_text(
            real["value_monotonicity_margin"], 70
        ),
        "full_real_derivative_monotonicity_margin_lower": arb_lower_text(
            real["derivative_monotonicity_margin"], 70
        ),
        "full_real_majorants_decrease_after_split": real[
            "monotonicity_certified"
        ],
        "upper_gamma_hazard_bound": (
            "Gamma(s,y)<=y^s*exp(-y)/(y-s+1) for y>s-1; this makes every scaled tail term decrease once y-s+1 exceeds its positive T-power."
        ),
        "ray_residual_rows": ray_rows,
        "ray_residual_row_count": len(ray_rows),
        "maximum_ray_value_fraction_of_budget_upper": worst_value[
            "value_fraction_of_budget_upper"
        ],
        "maximum_ray_value_fraction_location": worst_value["row_label"],
        "maximum_ray_derivative_fraction_of_budget_upper": worst_derivative[
            "derivative_fraction_of_budget_upper"
        ],
        "maximum_ray_derivative_fraction_location": worst_derivative["row_label"],
        "all_four_ray_residual_budgets_certified": all(
            row["ray_residual_budget_certified"] for row in ray_rows
        ),
        "all_upper_gamma_scaled_tails_decrease_on_ray": all(
            row["upper_gamma_scaled_tails_decrease_on_ray"] for row in ray_rows
        ),
        "perturbation_rows": perturbations["rows"],
        "ray_fixed_k_stencil_system_certified": ray_stencil,
        "full_real_T_fixed_k_stencil_system_certified": full_collar,
        "remaining_full_T_fixed_k_stencil_sources": [],
        "remaining_obligations": [
            "extend fixed k=22 to the complete required k range",
            "activate the cone-invariance theorem with an all-k entry statement",
            "complete the sign-regularity and Newman-direction bridges",
        ],
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgfrtfksc_01_order42_source",
            role="exact_analytic_import",
            readiness="available_exact",
            claim=(
                "Import exact full-kernel evenness, disk analyticity, and the order-42 zero of the degree-40 normalized residual."
            ),
            diagnostics={
                "full_kernel_evenness_certified": diagnostics[
                    "source_full_kernel_evenness_certified"
                ],
                "order42_residual_zero_certified": diagnostics[
                    "source_order42_residual_zero_certified"
                ],
            },
            source_artifacts=[paths["evenness_json"], paths["evenness_note"]],
            proof_boundary="Exact analytic source import only.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_02_full_disk_majorant",
            role="analytic_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "An Arb n<=80 sum plus geometric n-tail bounds gives a full infinite-kernel majorant on |z|<=0.38."
            ),
            diagnostics={
                "full_disk_majorant_upper": diagnostics["full_disk_majorant_upper"],
                "degree4_tail_upper": diagnostics["full_disk_degree4_tail_upper"],
                "degree2_tail_upper": diagnostics["full_disk_degree2_tail_upper"],
            },
            source_artifacts=[paths["evenness_note"], paths["node_c0_note"]],
            proof_boundary="Full-kernel local disk majorant only.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_03_factored_compact_ray_bound",
            role="analytic_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The order-42 Cauchy factors integrate against the full Gamma moment, yielding value T^-18 and derivative T^-19 scaled decay on T>=10000."
            ),
            diagnostics={
                "first_residual_degree": diagnostics["parameters"][
                    "first_residual_degree"
                ],
                "ray_residual_rows": diagnostics["ray_residual_row_count"],
            },
            source_artifacts=[paths["evenness_note"], paths["uniform_remainder_target"]],
            proof_boundary="Compact full-kernel residual bound on T>=10000 only.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_04_full_real_tail_majorant",
            role="analytic_tail_certificate",
            readiness="available_tail_certificate",
            claim=(
                "Full-kernel value and x*Phi' majorants, including geometric n>80 tails, decrease for x>=1/5."
            ),
            diagnostics={
                "value_majorant_upper": diagnostics[
                    "full_real_value_majorant_upper_at_split"
                ],
                "derivative_majorant_upper": diagnostics[
                    "full_real_derivative_majorant_upper_at_split"
                ],
                "value_monotonicity_margin_lower": diagnostics[
                    "full_real_value_monotonicity_margin_lower"
                ],
                "derivative_monotonicity_margin_lower": diagnostics[
                    "full_real_derivative_monotonicity_margin_lower"
                ],
            },
            source_artifacts=[paths["evenness_note"], paths["degree40_budget_note"]],
            proof_boundary="Full-kernel real-tail majorant on x>=1/5 only.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_05_scaled_upper_gamma_monotonicity",
            role="exact_analytic_lemma",
            readiness="available_exact",
            claim=(
                "The upper-Gamma hazard inequality proves every T-scaled real-tail term decreases from T=10000 onward."
            ),
            formula=diagnostics["upper_gamma_hazard_bound"],
            diagnostics={
                "all_scaled_tails_decrease": diagnostics[
                    "all_upper_gamma_scaled_tails_decrease_on_ray"
                ]
            },
            source_artifacts=[paths["uniform_remainder_target"]],
            proof_boundary="Exact monotonicity lemma for the recorded shapes and powers on T>=10000.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_06_ray_residual_budget_certificate",
            role="real_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "For every real T>=10000 and i=21..24, the full residuals satisfy |R_i|<=(1/2)T^-3 and |R_i'|<=(9/1000)T^-1."
            ),
            diagnostics={
                "ray_residual_rows": diagnostics["ray_residual_rows"],
                "maximum_value_fraction_of_budget_upper": diagnostics[
                    "maximum_ray_value_fraction_of_budget_upper"
                ],
                "maximum_derivative_fraction_of_budget_upper": diagnostics[
                    "maximum_ray_derivative_fraction_of_budget_upper"
                ],
                "all_four_ray_residual_budgets_certified": diagnostics[
                    "all_four_ray_residual_budgets_certified"
                ],
            },
            source_artifacts=[paths["evenness_note"], paths["degree40_budget_note"]],
            proof_boundary="Four-index residual theorem on T>=10000 only.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_07_ray_stencil_certificate",
            role="real_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The rational residual budgets and Arb perturbation ledger certify the fixed-k=22 stencil system for every T>=10000."
            ),
            diagnostics={
                "perturbation_rows": diagnostics["perturbation_rows"],
                "ray_fixed_k_stencil_system_certified": diagnostics[
                    "ray_fixed_k_stencil_system_certified"
                ],
            },
            source_artifacts=[paths["degree40_ladder_note"], paths["degree40_budget_note"]],
            proof_boundary="Fixed k=22 unbounded-ray stencil theorem only.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_08_full_real_T_fixed_k_composition",
            role="real_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The bounded segment 1156<=T<=10000 and unbounded ray T>=10000 compose to certify the fixed-k=22 stencil system for every real T>=1156."
            ),
            diagnostics={
                "finite_segment_certified": diagnostics[
                    "source_finite_segment_certified"
                ],
                "ray_certified": diagnostics["ray_fixed_k_stencil_system_certified"],
                "full_real_T_fixed_k_stencil_system_certified": diagnostics[
                    "full_real_T_fixed_k_stencil_system_certified"
                ],
                "remaining_sources": diagnostics[
                    "remaining_full_T_fixed_k_stencil_sources"
                ],
            },
            source_artifacts=[paths["finite_segment_json"], paths["finite_segment_note"]],
            proof_boundary="Full real-T theorem for one fixed k-window only; not all k or cone entry.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_09_all_k_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The full real-T fixed-k=22 theorem proves the required cone inequalities for all k.",
            gap="Only F_21 through F_24 are controlled; no all-k propagation theorem is supplied.",
            source_artifacts=[paths["cone_entry_target"], paths["dependency_graph"]],
            proof_boundary="Rejected fixed-k-to-all-k promotion only.",
        ),
        MatrixRow(
            id="nlrgfrtfksc_10_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "Any cone-entry promotion must supply the remaining k coverage and then pass the separate heat-flow invariance and sign-regularity bridges."
            ),
            source_artifacts=[paths["cone_entry_target"], paths["dependency_graph"]],
            proof_boundary="Proof-hygiene gate only; not cone entry, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(evenness_path: Path, finite_segment_path: Path) -> dict:
    evenness = load_json(evenness_path)
    finite_segment = load_json(finite_segment_path)
    paths = source_paths(evenness_path, finite_segment_path)
    diagnostics = build_diagnostics(evenness, finite_segment)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "ray_residual_rows": diagnostics["ray_residual_row_count"],
        "full_kernel_n_tail_channels": 4,
        "positive_retained_perturbation_margins": sum(
            1
            for row in diagnostics["perturbation_rows"]
            if row["certified_positive_after_perturbation"]
        ),
        "maximum_ray_value_fraction_of_budget_upper": diagnostics[
            "maximum_ray_value_fraction_of_budget_upper"
        ],
        "maximum_ray_derivative_fraction_of_budget_upper": diagnostics[
            "maximum_ray_derivative_fraction_of_budget_upper"
        ],
        "ray_fixed_k_stencil_system_certified": diagnostics[
            "ray_fixed_k_stencil_system_certified"
        ],
        "full_real_T_fixed_k_stencil_system_certified": diagnostics[
            "full_real_T_fixed_k_stencil_system_certified"
        ],
        "remaining_full_T_fixed_k_stencil_source_count": len(
            diagnostics["remaining_full_T_fixed_k_stencil_sources"]
        ),
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "Exact full-kernel evenness gives an order-42 residual zero. A full-kernel disk majorant, "
            "factored Cauchy bound, full real-tail majorant, and upper-Gamma hazard monotonicity certify "
            "the rational residual budgets on every T>=10000. Combined with the bounded segment theorem, "
            "the complete fixed-k=22 stencil system is now certified for every real T>=1156. Remaining "
            "k coverage, cone entry, sign regularity, RH, and Lambda <= 0 remain open."
        ),
    }
    return {
        "kind": (
            "jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate"
        ),
        "date": "2026-07-10",
        "status": "full real-T fixed-k stencil certificate",
        "source_full_kernel_evenness_lemma": paths["evenness_note"],
        "source_full_kernel_evenness_json": paths["evenness_json"],
        "source_finite_segment_certificate": paths["finite_segment_note"],
        "source_finite_segment_json": paths["finite_segment_json"],
        "source_degree40_ladder": paths["degree40_ladder_note"],
        "source_degree40_budget": paths["degree40_budget_note"],
        "source_node_c0_certificate": paths["node_c0_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_cone_entry_target": paths["cone_entry_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py"
        ),
        "proof_boundary": (
            "Uniform full real-T fixed-k=22 residual and stencil theorem for every T>=1156. It does "
            "not cover all k, does not by itself prove cone entry or sign regularity, and does not prove "
            "RH or Lambda <= 0."
        ),
        "diagnostics": diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "The bounded segment and unbounded ray share T=10000 and cover every real T>=1156.",
            "Full-kernel evenness, not finite-truncation near-evenness, supplies the order-42 zero.",
            "All full-kernel n tails are included by geometric Arb majorants.",
            "Every scaled upper-Gamma tail term is proved decreasing on T>=10000.",
            "Fixed k=22 is not promoted to all k or to cone entry.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian full-real-T fixed-k stencil certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, {summary['ray_residual_rows']} ray residual rows, "
        f"{summary['full_kernel_n_tail_channels']} full-kernel n-tail channels, "
        f"{summary['positive_retained_perturbation_margins']} positive perturbation margins, "
        f"{summary['remaining_full_T_fixed_k_stencil_source_count']} open full-T fixed-k stencil sources, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    diagnostics = artifact["diagnostics"]
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Full-Real-T Fixed-k Stencil Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: full real-T fixed-k stencil certificate. This is not a proof",
        "of all-k cone entry, sign regularity, RH, or `Lambda <= 0`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.json",
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
        "The bounded segment certificate covers `1156<=T<=10000`. For the",
        "unbounded ray, fix `x_*=1/5`, so the real-tail split starts at",
        "`y_*=T/25>=400`.",
        "",
        "## Full-Kernel Ray Bounds",
        "",
        f"- Full disk majorant on `|z|<=0.38`: `{diagnostics['full_disk_majorant_upper']}`.",
        f"- Full real value majorant at `x=1/5`: `{diagnostics['full_real_value_majorant_upper_at_split']}`.",
        f"- Full real derivative majorant at `x=1/5`: `{diagnostics['full_real_derivative_majorant_upper_at_split']}`.",
        f"- Value monotonicity margin: `{diagnostics['full_real_value_monotonicity_margin_lower']}`.",
        f"- Derivative monotonicity margin: `{diagnostics['full_real_derivative_monotonicity_margin_lower']}`.",
        "",
        "Exact evenness and subtraction through degree 40 factor the compact",
        "residual by `x^42`. After Gamma integration, the scaled compact value",
        "and derivative bounds decay as `T^-18` and `T^-19`.",
        "",
        "For the real tail,",
        "",
        "```text",
        diagnostics["upper_gamma_hazard_bound"],
        "```",
        "",
        "proves every scaled tail term decreases from `T=10000` onward.",
        "",
        "## Ray Residual Rows",
        "",
        "| Index | value scaled upper | fraction of 1/2 | derivative scaled upper | fraction of 9/1000 | hazard margin |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in diagnostics["ray_residual_rows"]:
        lines.append(
            f"| `F_{row['index']}` | `{row['value_residual_scaled_uniform_upper']}` | "
            f"`{row['value_fraction_of_budget_upper']}` | "
            f"`{row['derivative_residual_scaled_uniform_upper']}` | "
            f"`{row['derivative_fraction_of_budget_upper']}` | "
            f"`{row['upper_gamma_minimum_monotonicity_margin_lower']}` |"
        )
    lines.extend(
        [
            "",
            "## Full T Composition",
            "",
            "```text",
            "[1156,10000] union [10000,infinity) = [1156,infinity).",
            "```",
            "",
            "The exact rational residual budgets therefore preserve all four",
            "normalizers, the B product, companion product, and weighted-gap",
            "derivative for fixed `k=22` at every real `T>=1156`.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
            "The missing cone-entry work has shifted from real-T control to k",
            "coverage: this theorem controls only the window `F_21,...,F_24`.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_full_real_T_fixed_k_stencil_certificate.py",
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
    parser.add_argument("--evenness-json", type=Path, default=DEFAULT_EVENNESS_JSON)
    parser.add_argument("--finite-segment-json", type=Path, default=DEFAULT_FINITE_SEGMENT_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = build_artifact(args.evenness_json, args.finite_segment_json)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, args.note)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
