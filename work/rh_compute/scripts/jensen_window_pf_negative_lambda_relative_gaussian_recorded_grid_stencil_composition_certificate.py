#!/usr/bin/env python3
"""Compose the certified expectation grid with rigorous degree-40 stencil budgets."""

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
    DEFAULT_INDICES,
    DEFAULT_T_GRID,
    arb_lower_text,
    arb_text,
    arb_upper_text,
    parse_arb,
)
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_diagnostics as build_arb_collar_diagnostics,
    build_ratio_rows,
    pderivative,
    truncated_multiplier_polynomial,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate import (  # noqa: E402
    abs_upper,
)
from jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget import (  # noqa: E402
    max_abs_bernstein_upper,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_EXPECTATION_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.json"
)
DEFAULT_RESIDUAL_BUDGET_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.md"
)

DEFAULT_FINITE_DEGREE = 40
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_K = 22
DEFAULT_COLLAR_START_T = 1156
DEFAULT_PRECISION_BITS = 4096
VALUE_BUDGET_NUMERATOR = 1
VALUE_BUDGET_DENOMINATOR = 2
DERIVATIVE_BUDGET_NUMERATOR = 9
DERIVATIVE_BUDGET_DENOMINATOR = 1000


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


def source_paths(expectation_path: Path, residual_budget_path: Path) -> dict[str, str]:
    return {
        "expectation_json": expectation_path.relative_to(REPO_ROOT).as_posix(),
        "expectation_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_all_row_direct_expectation_certificate.md"
        ),
        "residual_budget_json": residual_budget_path.relative_to(REPO_ROOT).as_posix(),
        "residual_budget_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md"
        ),
        "degree40_ladder_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md"
        ),
        "stencil_obligations_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md"
        ),
        "intervalization_target": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def configure_precision() -> None:
    flint.ctx.prec = DEFAULT_PRECISION_BITS


def minimum_ball_lower(values: list[str]) -> flint.arb:
    balls = [parse_arb(value) for value in values]
    lowers = [flint.arb(value.lower()) for value in balls]
    return min(lowers)


def finite_collar_data() -> dict:
    collar = build_arb_collar_diagnostics(
        DEFAULT_FINITE_DEGREE,
        DEFAULT_RATIO_CUTOFF_N,
        DEFAULT_PRECISION_BITS,
        DEFAULT_K,
        DEFAULT_POLYNOMIAL_M,
        DEFAULT_COLLAR_START_T,
    )
    normalizer_margin = minimum_ball_lower(
        [row["min_bernstein_lower"] for row in collar["normalizer_rows"]]
    )
    stencil_margins = {
        row["name"]: flint.arb(parse_arb(row["min_bernstein_lower"]).lower())
        for row in collar["stencil_rows"]
    }
    _ratio_rows, ratios = build_ratio_rows(
        DEFAULT_FINITE_DEGREE, DEFAULT_RATIO_CUTOFF_N
    )
    finite_polys = {
        index: truncated_multiplier_polynomial(index, ratios, DEFAULT_POLYNOMIAL_M)
        for index in DEFAULT_INDICES
    }
    pmax = flint.arb(0)
    dmax = flint.arb(0)
    for poly in finite_polys.values():
        pmax = max(pmax, flint.arb(max_abs_bernstein_upper(poly, DEFAULT_COLLAR_START_T)))
        dmax = max(
            dmax,
            flint.arb(max_abs_bernstein_upper(pderivative(poly), DEFAULT_COLLAR_START_T)),
        )
    return {
        "normalizer_margin": normalizer_margin,
        "B_product_margin": stencil_margins["B_product"],
        "companion_product_margin": stencil_margins["companion_product"],
        "weighted_gap_derivative_margin": stencil_margins["weighted_gap_derivative"],
        "finite_normalizer_abs_upper": pmax,
        "finite_derivative_abs_upper": dmax,
    }


def perturbation_ledger(collar: dict) -> dict:
    U = flint.arb(1) / DEFAULT_COLLAR_START_T
    A = flint.arb(VALUE_BUDGET_NUMERATOR) / VALUE_BUDGET_DENOMINATOR
    B = flint.arb(DERIVATIVE_BUDGET_NUMERATOR) / DERIVATIVE_BUDGET_DENOMINATOR
    P = collar["finite_normalizer_abs_upper"]
    D = collar["finite_derivative_abs_upper"]
    weight_sum = flint.arb(368)

    normalizer_bound = A * (U**3)
    b_product_bound = 4 * P * A * U + 2 * A * A * (U**4)
    companion_bound = 2 * (
        4 * (P**3) * A
        + 6 * (P**2) * (A**2) * (U**3)
        + 4 * P * (A**3) * (U**6)
        + (A**4) * (U**9)
    )
    weighted_value_bound = weight_sum * D * (
        3 * (P**2) * A * (U**2)
        + 3 * P * (A**2) * (U**5)
        + (A**3) * (U**8)
    )
    weighted_derivative_bound = weight_sum * B * ((P + A * (U**3)) ** 3)
    weighted_total_bound = weighted_value_bound + weighted_derivative_bound

    entries = [
        (
            "normalizer",
            collar["normalizer_margin"],
            normalizer_bound,
            "|R_i| <= A*u^3",
        ),
        (
            "B_product",
            collar["B_product_margin"],
            b_product_bound,
            "4*P*A*U + 2*A^2*U^4",
        ),
        (
            "companion_product",
            collar["companion_product_margin"],
            companion_bound,
            "2*(4*P^3*A + 6*P^2*A^2*U^3 + 4*P*A^3*U^6 + A^4*U^9)",
        ),
        (
            "weighted_gap_derivative",
            collar["weighted_gap_derivative_margin"],
            weighted_total_bound,
            "368*(D*(3*P^2*A*U^2+3*P*A^2*U^5+A^3*U^8)+B*(P+A*U^3)^3)",
        ),
    ]
    rows = []
    for name, finite_margin, bound, formula in entries:
        retained = flint.arb(finite_margin.lower()) - flint.arb(bound.upper())
        certified = bool(retained > 0)
        if not certified:
            raise RuntimeError(f"{name} perturbation budget did not retain a positive margin")
        rows.append(
            {
                "name": name,
                "finite_margin_lower": arb_lower_text(finite_margin, 70),
                "perturbation_bound_upper": arb_upper_text(bound, 70),
                "retained_margin_lower": arb_lower_text(retained, 70),
                "perturbation_formula": formula,
                "certified_positive_after_perturbation": certified,
                "proof_boundary": "Arb perturbation inequality for the rational A,B budgets on 0<=u<=1/1156.",
            }
        )
    return {
        "U": U,
        "A": A,
        "B": B,
        "P": P,
        "D": D,
        "normalizer_bound": normalizer_bound,
        "B_product_bound": b_product_bound,
        "companion_bound": companion_bound,
        "weighted_value_bound": weighted_value_bound,
        "weighted_derivative_bound": weighted_derivative_bound,
        "weighted_total_bound": weighted_total_bound,
        "rows": rows,
        "all_retained_margins_positive": all(
            row["certified_positive_after_perturbation"] for row in rows
        ),
    }


def expectation_rows_by_key(expectation: dict) -> dict[tuple[int, int], dict]:
    rows = expectation["diagnostics"]["grid_rows"]
    return {(int(row["T"]), int(row["index"])): row for row in rows}


def build_grid_rows(expectation: dict, perturbations: dict) -> tuple[list[dict], list[dict]]:
    by_key = expectation_rows_by_key(expectation)
    expected_keys = {(T, index) for T in DEFAULT_T_GRID for index in DEFAULT_INDICES}
    if set(by_key) != expected_keys:
        raise RuntimeError("expectation source does not contain the required 20 rows")
    A = perturbations["A"]
    B = perturbations["B"]
    detail_rows: list[dict] = []
    t_rows: list[dict] = []
    for T in DEFAULT_T_GRID:
        value_bounds: list[flint.arb] = []
        derivative_bounds: list[flint.arb] = []
        all_rows_certified = True
        for index in DEFAULT_INDICES:
            source = by_key[(T, index)]
            u = flint.arb(1) / T
            value_ball = parse_arb(source["full_value_expectation_ball"])
            derivative_ball = parse_arb(source["full_derivative_expectation_ball"])
            value_scaled = abs_upper(value_ball) / (u**3)
            derivative_scaled = abs_upper(derivative_ball) / (u**2)
            value_fraction = value_scaled / A
            derivative_fraction = derivative_scaled / B
            certified = bool(
                source["row_direct_expectation_certified"]
                and value_scaled < A
                and derivative_scaled < B
            )
            if not certified:
                raise RuntimeError(f"T={T}, F_{index} did not fit the rational residual budgets")
            value_bounds.append(value_scaled)
            derivative_bounds.append(derivative_scaled)
            all_rows_certified = all_rows_certified and certified
            detail_rows.append(
                {
                    "T": T,
                    "index": index,
                    "row_label": f"T={T}, F_{index}",
                    "u": f"1/{T}",
                    "source_value_expectation_ball": source["full_value_expectation_ball"],
                    "source_derivative_expectation_ball": source[
                        "full_derivative_expectation_ball"
                    ],
                    "value_residual_scaled_upper": arb_upper_text(value_scaled, 70),
                    "derivative_residual_scaled_upper": arb_upper_text(
                        derivative_scaled, 70
                    ),
                    "value_fraction_of_rational_budget_upper": arb_upper_text(
                        value_fraction, 70
                    ),
                    "derivative_fraction_of_rational_budget_upper": arb_upper_text(
                        derivative_fraction, 70
                    ),
                    "value_budget_certified": bool(value_scaled < A),
                    "derivative_budget_certified": bool(derivative_scaled < B),
                    "row_residual_budget_certified": certified,
                    "identity_used": (
                        "value expectation=R_i(u); derivative expectation=u*R_i'(u)"
                    ),
                    "proof_boundary": "One finite T,F_i residual-budget row only.",
                }
            )
        max_value = max(value_bounds)
        max_derivative = max(derivative_bounds)
        t_certified = bool(
            all_rows_certified and perturbations["all_retained_margins_positive"]
        )
        t_rows.append(
            {
                "T": T,
                "u": f"1/{T}",
                "certified_indices": list(DEFAULT_INDICES),
                "maximum_value_residual_scaled_upper": arb_upper_text(max_value, 70),
                "maximum_derivative_residual_scaled_upper": arb_upper_text(
                    max_derivative, 70
                ),
                "all_four_normalizers_certified_positive": t_certified,
                "B_product_certified_positive": t_certified,
                "companion_product_certified_positive": t_certified,
                "weighted_gap_derivative_certified_positive": t_certified,
                "recorded_T_stencil_system_certified": t_certified,
                "proof_boundary": (
                    "Fixed k=22 stencil system at one recorded T only; not an interval in T or an all-k cone theorem."
                ),
            }
        )
    return detail_rows, t_rows


def build_diagnostics(expectation: dict, residual_budget: dict) -> dict:
    configure_precision()
    if expectation["summary"]["complete_recorded_grid_expectations_certified"] is not True:
        raise RuntimeError("source expectation grid is not complete")
    collar = finite_collar_data()
    perturbations = perturbation_ledger(collar)
    detail_rows, t_rows = build_grid_rows(expectation, perturbations)
    worst_value = max(
        detail_rows, key=lambda row: parse_arb(row["value_fraction_of_rational_budget_upper"])
    )
    worst_derivative = max(
        detail_rows,
        key=lambda row: parse_arb(row["derivative_fraction_of_rational_budget_upper"]),
    )
    return {
        "parameters": {
            "T_grid": list(DEFAULT_T_GRID),
            "indices": list(DEFAULT_INDICES),
            "k": DEFAULT_K,
            "finite_degree": DEFAULT_FINITE_DEGREE,
            "polynomial_M": DEFAULT_POLYNOMIAL_M,
            "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
            "collar_start_T": DEFAULT_COLLAR_START_T,
            "precision_bits": DEFAULT_PRECISION_BITS,
            "value_budget_A": "1/2",
            "derivative_budget_B": "9/1000",
            "value_residual_identity": "expectation_value=R_i(u)",
            "derivative_residual_identity": "expectation_derivative=u*R_i'(u)",
        },
        "source_float_budget_A": residual_budget["summary"][
            "value_residual_half_safety_budget_A"
        ],
        "source_float_budget_B": residual_budget["summary"][
            "derivative_residual_half_safety_budget_B"
        ],
        "source_float_budgets_used_as_proof_inputs": False,
        "rational_value_budget_ball": arb_text(perturbations["A"], 50),
        "rational_derivative_budget_ball": arb_text(perturbations["B"], 50),
        "finite_normalizer_abs_upper": arb_upper_text(
            collar["finite_normalizer_abs_upper"], 70
        ),
        "finite_derivative_abs_upper": arb_upper_text(
            collar["finite_derivative_abs_upper"], 70
        ),
        "finite_margin_balls": {
            "normalizer": arb_text(collar["normalizer_margin"], 70),
            "B_product": arb_text(collar["B_product_margin"], 70),
            "companion_product": arb_text(collar["companion_product_margin"], 70),
            "weighted_gap_derivative": arb_text(
                collar["weighted_gap_derivative_margin"], 70
            ),
        },
        "perturbation_rows": perturbations["rows"],
        "perturbation_row_count": len(perturbations["rows"]),
        "weighted_gap_value_perturbation_upper": arb_upper_text(
            perturbations["weighted_value_bound"], 70
        ),
        "weighted_gap_derivative_residual_perturbation_upper": arb_upper_text(
            perturbations["weighted_derivative_bound"], 70
        ),
        "all_retained_margins_positive": perturbations[
            "all_retained_margins_positive"
        ],
        "grid_residual_rows": detail_rows,
        "grid_residual_row_count": len(detail_rows),
        "recorded_T_stencil_rows": t_rows,
        "recorded_T_stencil_row_count": len(t_rows),
        "worst_value_budget_fraction_location": {
            "T": worst_value["T"],
            "index": worst_value["index"],
        },
        "worst_derivative_budget_fraction_location": {
            "T": worst_derivative["T"],
            "index": worst_derivative["index"],
        },
        "maximum_value_fraction_of_rational_budget_upper": worst_value[
            "value_fraction_of_rational_budget_upper"
        ],
        "maximum_derivative_fraction_of_rational_budget_upper": worst_derivative[
            "derivative_fraction_of_rational_budget_upper"
        ],
        "all_20_value_residual_budgets_certified": all(
            row["value_budget_certified"] for row in detail_rows
        ),
        "all_20_derivative_residual_budgets_certified": all(
            row["derivative_budget_certified"] for row in detail_rows
        ),
        "all_5_recorded_T_stencil_systems_certified": all(
            row["recorded_T_stencil_system_certified"] for row in t_rows
        ),
        "remaining_recorded_grid_stencil_sources": [],
        "remaining_obligations": [
            "certify the stencil system for every real T>=1156 rather than five points",
            "extend the fixed k=22 result to the complete required k range",
            "complete the downstream cone-invariance and sign-regularity bridges",
        ],
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgrgscc_01_residual_expectation_identity",
            role="exact_reduction",
            readiness="available_exact",
            claim=(
                "For each i, the certified value expectation is R_i(u), while the certified derivative "
                "expectation is u*R_i'(u) for F_i=P_i^(40)+R_i."
            ),
            formula="E[value_core]=R_i(u); E[derivative_core]=u*R_i'(u)",
            source_artifacts=[paths["expectation_note"], paths["residual_budget_note"]],
            proof_boundary="Exact coordinate identity only.",
        ),
        MatrixRow(
            id="nlrgrgscc_02_rational_residual_budgets",
            role="exact_sufficient_condition",
            readiness="available_exact",
            claim=(
                "Use the exact rational budgets |R_i(u)|<=(1/2)u^3 and "
                "|R_i'(u)|<=(9/1000)u, avoiding the source artifact's floating threshold decimals."
            ),
            diagnostics={
                "value_budget_A": diagnostics["parameters"]["value_budget_A"],
                "derivative_budget_B": diagnostics["parameters"]["derivative_budget_B"],
                "source_float_budgets_used_as_proof_inputs": diagnostics[
                    "source_float_budgets_used_as_proof_inputs"
                ],
            },
            source_artifacts=[paths["residual_budget_json"], paths["residual_budget_note"]],
            proof_boundary="Exact rational sufficient budgets only; satisfaction is certified separately.",
        ),
        MatrixRow(
            id="nlrgrgscc_03_arb_perturbation_ledger",
            role="analytic_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "Arb recomputation proves the rational A,B budgets preserve positive normalizer, "
                "B-product, companion-product, and weighted-gap derivative margins."
            ),
            diagnostics={
                "perturbation_rows": diagnostics["perturbation_rows"],
                "all_retained_margins_positive": diagnostics[
                    "all_retained_margins_positive"
                ],
            },
            source_artifacts=[paths["degree40_ladder_note"], paths["stencil_obligations_note"]],
            proof_boundary="Fixed k=22 degree-40 perturbation ledger on 0<=u<=1/1156.",
        ),
        MatrixRow(
            id="nlrgrgscc_04_all_row_budget_composition",
            role="finite_grid_interval_certificate",
            readiness="available_finite_grid_certificate",
            claim=(
                "All 20 certified expectation rows satisfy both exact rational residual budgets after "
                "the derivative expectation is divided by u^2."
            ),
            diagnostics={
                "grid_residual_rows": diagnostics["grid_residual_row_count"],
                "maximum_value_fraction_of_rational_budget_upper": diagnostics[
                    "maximum_value_fraction_of_rational_budget_upper"
                ],
                "maximum_derivative_fraction_of_rational_budget_upper": diagnostics[
                    "maximum_derivative_fraction_of_rational_budget_upper"
                ],
                "all_20_value_residual_budgets_certified": diagnostics[
                    "all_20_value_residual_budgets_certified"
                ],
                "all_20_derivative_residual_budgets_certified": diagnostics[
                    "all_20_derivative_residual_budgets_certified"
                ],
            },
            source_artifacts=[paths["expectation_json"], paths["expectation_note"]],
            proof_boundary="Twenty recorded expectation rows only; not continuous in T.",
        ),
        MatrixRow(
            id="nlrgrgscc_05_recorded_T_stencil_certificate",
            role="finite_grid_interval_certificate",
            readiness="available_finite_grid_certificate",
            claim=(
                "At every recorded T, all four full normalizers and the fixed-k B-product, companion "
                "product, and weighted-gap derivative are rigorously positive."
            ),
            diagnostics={
                "recorded_T_stencil_rows": diagnostics["recorded_T_stencil_rows"],
                "all_5_recorded_T_stencil_systems_certified": diagnostics[
                    "all_5_recorded_T_stencil_systems_certified"
                ],
                "remaining_recorded_grid_stencil_sources": diagnostics[
                    "remaining_recorded_grid_stencil_sources"
                ],
            },
            source_artifacts=[paths["expectation_note"], paths["degree40_ladder_note"]],
            proof_boundary="Fixed k=22 at five T values only; not a real-T collar or all-k theorem.",
        ),
        MatrixRow(
            id="nlrgrgscc_06_real_T_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The five certified stencil systems prove the same inequalities for every T>=1156.",
            gap="No interval subdivision or monotonicity theorem fills the real-T gaps.",
            source_artifacts=[paths["uniform_remainder_target"], paths["cone_entry_target"]],
            proof_boundary="Rejected finite-grid-to-collar promotion only.",
        ),
        MatrixRow(
            id="nlrgrgscc_07_all_k_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The fixed k=22 stencil certificate proves the required inequalities for all k.",
            gap="The composition uses only F_21 through F_24 and therefore certifies one k-window.",
            source_artifacts=[paths["cone_entry_target"], paths["dependency_graph"]],
            proof_boundary="Rejected fixed-k-to-all-k promotion only.",
        ),
        MatrixRow(
            id="nlrgrgscc_08_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "Promotion requires both a real-T collar proof and the remaining k coverage; finite "
                "sampling, floating thresholds, or RH-assuming inputs cannot replace either bridge."
            ),
            source_artifacts=[paths["uniform_remainder_target"], paths["dependency_graph"]],
            proof_boundary="Proof-hygiene gate only; not cone entry, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(expectation_path: Path, residual_budget_path: Path) -> dict:
    expectation = load_json(expectation_path)
    residual_budget = load_json(residual_budget_path)
    paths = source_paths(expectation_path, residual_budget_path)
    diagnostics = build_diagnostics(expectation, residual_budget)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "grid_residual_rows": diagnostics["grid_residual_row_count"],
        "recorded_T_stencil_rows": diagnostics["recorded_T_stencil_row_count"],
        "positive_retained_perturbation_margins": sum(
            1
            for row in diagnostics["perturbation_rows"]
            if row["certified_positive_after_perturbation"]
        ),
        "certified_recorded_T_stencil_systems": sum(
            1
            for row in diagnostics["recorded_T_stencil_rows"]
            if row["recorded_T_stencil_system_certified"]
        ),
        "maximum_value_fraction_of_rational_budget_upper": diagnostics[
            "maximum_value_fraction_of_rational_budget_upper"
        ],
        "maximum_derivative_fraction_of_rational_budget_upper": diagnostics[
            "maximum_derivative_fraction_of_rational_budget_upper"
        ],
        "all_5_recorded_T_stencil_systems_certified": diagnostics[
            "all_5_recorded_T_stencil_systems_certified"
        ],
        "remaining_recorded_grid_stencil_source_count": len(
            diagnostics["remaining_recorded_grid_stencil_sources"]
        ),
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The 20 direct expectation balls satisfy the exact rational residual budgets A=1/2 and "
            "B=9/1000. An independent Arb perturbation ledger preserves all degree-40 fixed-k=22 "
            "normalizer, B-product, companion-product, and weighted-gap derivative margins. Therefore "
            "the complete fixed-k stencil system is certified at all five recorded T values. The real-T "
            "collar, remaining k coverage, cone entry, RH, and Lambda <= 0 remain open."
        ),
    }
    return {
        "kind": (
            "jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate"
        ),
        "date": "2026-07-10",
        "status": "recorded-grid fixed-k stencil composition certificate",
        "source_all_row_expectation_certificate": paths["expectation_note"],
        "source_all_row_expectation_json": paths["expectation_json"],
        "source_degree40_residual_budget": paths["residual_budget_note"],
        "source_degree40_residual_budget_json": paths["residual_budget_json"],
        "source_degree40_ladder": paths["degree40_ladder_note"],
        "source_stencil_obligations": paths["stencil_obligations_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_cone_entry_target": paths["cone_entry_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py"
        ),
        "proof_boundary": (
            "Complete fixed-k=22 stencil composition certificate at T=1156,1500,2000,5000,10000 "
            "only. It does not certify an interval in T, does not cover all k, does not prove the full "
            "cone-entry theorem or sign-regularity bridge, and does not prove RH or Lambda <= 0."
        ),
        "diagnostics": diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "The source floating A,B thresholds are not proof inputs.",
            "The proof budgets are the exact rationals A=1/2 and B=9/1000.",
            "All perturbation comparisons are recomputed with Arb ball arithmetic.",
            "Five isolated T values are not promoted to the real-T collar.",
            "Fixed k=22 is not promoted to all k.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian recorded-grid stencil composition certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, {summary['grid_residual_rows']} residual rows, "
        f"{summary['positive_retained_perturbation_margins']} positive perturbation margins, "
        f"{summary['certified_recorded_T_stencil_systems']} certified T systems, "
        f"{summary['remaining_recorded_grid_stencil_source_count']} open recorded-grid stencil sources, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    diagnostics = artifact["diagnostics"]
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Recorded-Grid Stencil Composition Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: recorded-grid fixed-k stencil composition certificate. This is not a proof",
        "of a real-T collar, an all-k cone-entry theorem, RH, or `Lambda <= 0`.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.json",
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
        "The exact residual identities are",
        "",
        "```text",
        "E[value core]      = R_i(u)",
        "E[derivative core] = u R_i'(u).",
        "```",
        "",
        "The composition uses the exact rational sufficient budgets",
        "",
        "```text",
        "|R_i(u)|  <= (1/2) u^3",
        "|R_i'(u)| <= (9/1000) u.",
        "```",
        "",
        "The older floating threshold decimals are recorded for comparison but are",
        "not proof inputs.",
        "",
        "## Perturbation Ledger",
        "",
        "| Target | finite margin lower | perturbation upper | retained lower |",
        "|---|---:|---:|---:|",
    ]
    for row in diagnostics["perturbation_rows"]:
        lines.append(
            f"| `{row['name']}` | `{row['finite_margin_lower']}` | "
            f"`{row['perturbation_bound_upper']}` | `{row['retained_margin_lower']}` |"
        )
    lines.extend(
        [
            "",
            "## Recorded T Systems",
            "",
            "| T | max |R|/u^3 | max |R'|/u | normalizers | B | companion | weighted gap |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in diagnostics["recorded_T_stencil_rows"]:
        lines.append(
            f"| `{row['T']}` | `{row['maximum_value_residual_scaled_upper']}` | "
            f"`{row['maximum_derivative_residual_scaled_upper']}` | "
            f"`{row['all_four_normalizers_certified_positive']}` | "
            f"`{row['B_product_certified_positive']}` | "
            f"`{row['companion_product_certified_positive']}` | "
            f"`{row['weighted_gap_derivative_certified_positive']}` |"
        )
    lines.extend(
        [
            "",
            f"Maximum fraction of the value budget: `{summary['maximum_value_fraction_of_rational_budget_upper']}`.",
            f"Maximum fraction of the derivative budget: `{summary['maximum_derivative_fraction_of_rational_budget_upper']}`.",
            "",
            "## Proof Boundary",
            "",
            artifact["proof_boundary"],
            "",
            "This closes the fixed-k stencil composition source at the five recorded",
            "T values. It does not interpolate between them and does not extend the",
            "single k-window to all k.",
            "",
            "## Reproduction",
            "",
            "```powershell",
            "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_recorded_grid_stencil_composition_certificate.py",
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
    parser.add_argument("--expectation-json", type=Path, default=DEFAULT_EXPECTATION_JSON)
    parser.add_argument("--residual-budget-json", type=Path, default=DEFAULT_RESIDUAL_BUDGET_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    artifact = build_artifact(args.expectation_json, args.residual_budget_json)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, args.note)
    print(result_line(artifact))


if __name__ == "__main__":
    main()
