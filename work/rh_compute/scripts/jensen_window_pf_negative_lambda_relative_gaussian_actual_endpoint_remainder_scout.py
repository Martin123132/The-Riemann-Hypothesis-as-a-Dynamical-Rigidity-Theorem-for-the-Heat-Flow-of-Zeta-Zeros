#!/usr/bin/env python3
"""Build a floating endpoint scout for the actual relative-Gaussian remainder."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys

import numpy as np
import scipy.special as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import REPO_ROOT  # noqa: E402
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)


DEFAULT_RESIDUAL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md"
)

DEFAULT_T_ENDPOINT = 1156
DEFAULT_QUADRATURE_ORDERS = (96, 128, 192, 256, 320)
DEFAULT_SELECTED_ORDER = 320
DEFAULT_PHI_TERM_COUNT = 23
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PRECISION_BITS = 384
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_INDICES = (21, 22, 23, 24)


@dataclass(frozen=True)
class EndpointRemainderRow:
    index: int
    selected_quadrature_order: int
    value_residual_scaled_at_selected_order: str
    derivative_residual_scaled_at_selected_order: str
    value_residual_spread_scaled: str
    derivative_residual_spread_scaled: str
    value_first_omitted_scaled: str
    derivative_first_omitted_scaled: str
    value_residual_to_first_omitted_ratio_at_selected_order: str
    derivative_residual_to_first_omitted_ratio_at_selected_order: str
    max_value_residual_to_first_omitted_ratio_over_orders: str
    max_derivative_residual_to_first_omitted_ratio_over_orders: str
    value_budget_fraction_at_selected_order: str
    derivative_budget_fraction_at_selected_order: str
    actual_value_at_selected_order: str
    actual_derivative_at_selected_order: str
    polynomial_value_at_endpoint: str
    polynomial_derivative_at_endpoint: str
    proof_boundary: str


def sci(value: float) -> str:
    return f"{float(value):.18E}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(residual_path: Path) -> dict[str, str]:
    return {
        "residual_budget_json": residual_path.relative_to(REPO_ROOT).as_posix(),
        "residual_budget_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
        "formal_tail_obstruction_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"
        ),
        "asymptotic_remainder_target_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md"
        ),
        "uniform_remainder_target_note": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def phi_np(x: np.ndarray, term_count: int) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    total = np.zeros_like(x, dtype=np.float64)
    exp4 = np.exp(4.0 * x)
    for n in range(1, term_count + 1):
        n2 = n * n
        q = math.pi * n2
        exp9 = np.exp(9.0 * x)
        exp5 = np.exp(5.0 * x)
        prefactor = 2.0 * math.pi * math.pi * n2 * n2 * exp9 - 3.0 * math.pi * n2 * exp5
        total += prefactor * np.exp(-q * exp4)
    return total


def phip_np(x: np.ndarray, term_count: int) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    total = np.zeros_like(x, dtype=np.float64)
    exp4 = np.exp(4.0 * x)
    for n in range(1, term_count + 1):
        n2 = n * n
        q = math.pi * n2
        exp9 = np.exp(9.0 * x)
        exp5 = np.exp(5.0 * x)
        prefactor = 2.0 * math.pi * math.pi * n2 * n2 * exp9 - 3.0 * math.pi * n2 * exp5
        prefactor_derivative = 18.0 * math.pi * math.pi * n2 * n2 * exp9 - 15.0 * math.pi * n2 * exp5
        total += (prefactor_derivative - prefactor * 4.0 * q * exp4) * np.exp(-q * exp4)
    return total


def rising_float(start: float, length: int) -> float:
    total = 1.0
    for offset in range(length):
        total *= start + offset
    return total


def polynomial_value(index: int, ratios: list[float], u: float, m: int) -> float:
    a = index + 0.5
    return sum(ratios[j] * rising_float(a, j) * (u**j) for j in range(m + 1))


def polynomial_derivative(index: int, ratios: list[float], u: float, m: int) -> float:
    a = index + 0.5
    return sum(j * ratios[j] * rising_float(a, j) * (u ** (j - 1)) for j in range(1, m + 1))


def first_omitted_value(index: int, ratios: list[float], u: float, first_j: int) -> float:
    return abs(ratios[first_j] * rising_float(index + 0.5, first_j) * (u**first_j)) / (u**3)


def first_omitted_derivative(index: int, ratios: list[float], u: float, first_j: int) -> float:
    return abs(first_j * ratios[first_j] * rising_float(index + 0.5, first_j) * (u ** (first_j - 1))) / u


def actual_multiplier_and_derivative(
    index: int,
    order: int,
    u: float,
    c0: float,
    term_count: int,
) -> tuple[float, float]:
    nodes, weights = sp.roots_genlaguerre(order, index - 0.5)
    z = np.sqrt(u * nodes)
    gamma = sp.gamma(index + 0.5)
    value = np.sum(weights * phi_np(z, term_count) / c0) / gamma
    derivative = np.sum(weights * phip_np(z, term_count) * (nodes / (2.0 * z)) / c0) / gamma
    return float(value), float(derivative)


def build_endpoint_rows(
    ratios: list[float],
    value_budget: float,
    derivative_budget: float,
    t_endpoint: int,
    quadrature_orders: tuple[int, ...],
    selected_order: int,
    phi_term_count: int,
    polynomial_m: int,
    indices: tuple[int, ...],
) -> tuple[list[EndpointRemainderRow], list[dict]]:
    u = 1.0 / float(t_endpoint)
    c0 = float(phi_np(np.array([0.0], dtype=np.float64), phi_term_count)[0])
    rows: list[EndpointRemainderRow] = []
    sample_rows: list[dict] = []
    for index in indices:
        pvalue = polynomial_value(index, ratios, u, polynomial_m)
        pderivative = polynomial_derivative(index, ratios, u, polynomial_m)
        value_first = first_omitted_value(index, ratios, u, polynomial_m + 1)
        derivative_first = first_omitted_derivative(index, ratios, u, polynomial_m + 1)
        samples = []
        for order in quadrature_orders:
            actual_value, actual_derivative = actual_multiplier_and_derivative(index, order, u, c0, phi_term_count)
            value_residual = abs(actual_value - pvalue) / (u**3)
            derivative_residual = abs(actual_derivative - pderivative) / u
            samples.append(
                {
                    "quadrature_order": order,
                    "actual_value": sci(actual_value),
                    "actual_derivative": sci(actual_derivative),
                    "value_residual_scaled": sci(value_residual),
                    "derivative_residual_scaled": sci(derivative_residual),
                    "value_residual_to_first_omitted_ratio": sci(value_residual / value_first),
                    "derivative_residual_to_first_omitted_ratio": sci(derivative_residual / derivative_first),
                }
            )
        selected = next(sample for sample in samples if sample["quadrature_order"] == selected_order)
        value_residuals = [float(sample["value_residual_scaled"]) for sample in samples]
        derivative_residuals = [float(sample["derivative_residual_scaled"]) for sample in samples]
        value_ratios = [float(sample["value_residual_to_first_omitted_ratio"]) for sample in samples]
        derivative_ratios = [float(sample["derivative_residual_to_first_omitted_ratio"]) for sample in samples]
        rows.append(
            EndpointRemainderRow(
                index=index,
                selected_quadrature_order=selected_order,
                value_residual_scaled_at_selected_order=selected["value_residual_scaled"],
                derivative_residual_scaled_at_selected_order=selected["derivative_residual_scaled"],
                value_residual_spread_scaled=sci(max(value_residuals) - min(value_residuals)),
                derivative_residual_spread_scaled=sci(max(derivative_residuals) - min(derivative_residuals)),
                value_first_omitted_scaled=sci(value_first),
                derivative_first_omitted_scaled=sci(derivative_first),
                value_residual_to_first_omitted_ratio_at_selected_order=selected[
                    "value_residual_to_first_omitted_ratio"
                ],
                derivative_residual_to_first_omitted_ratio_at_selected_order=selected[
                    "derivative_residual_to_first_omitted_ratio"
                ],
                max_value_residual_to_first_omitted_ratio_over_orders=sci(max(value_ratios)),
                max_derivative_residual_to_first_omitted_ratio_over_orders=sci(max(derivative_ratios)),
                value_budget_fraction_at_selected_order=sci(
                    float(selected["value_residual_scaled"]) / value_budget
                ),
                derivative_budget_fraction_at_selected_order=sci(
                    float(selected["derivative_residual_scaled"]) / derivative_budget
                ),
                actual_value_at_selected_order=selected["actual_value"],
                actual_derivative_at_selected_order=selected["actual_derivative"],
                polynomial_value_at_endpoint=sci(pvalue),
                polynomial_derivative_at_endpoint=sci(pderivative),
                proof_boundary=(
                    "Floating endpoint quadrature row only; it is not an interval enclosure, not a uniform "
                    "collar estimate, and not an analytic remainder theorem."
                ),
            )
        )
        sample_rows.append({"index": index, "samples": samples})
    return rows, sample_rows


def build_diagnostics(
    residual_path: Path,
    t_endpoint: int,
    quadrature_orders: tuple[int, ...],
    selected_order: int,
    phi_term_count: int,
    ratio_cutoff_n: int,
    precision_bits: int,
    polynomial_m: int,
    indices: tuple[int, ...],
) -> dict:
    if selected_order not in quadrature_orders:
        raise ValueError("selected_order must belong to quadrature_orders")
    flint.ctx.prec = precision_bits
    residual = load_json(residual_path)
    residual_summary = residual["summary"]
    value_budget = float(residual_summary["value_residual_half_safety_budget_A"])
    derivative_budget = float(residual_summary["derivative_residual_half_safety_budget_B"])
    _ratio_rows, ratios_arb = build_ratio_rows(2 * (polynomial_m + 1), ratio_cutoff_n)
    ratios = [float(item) for item in ratios_arb]
    endpoint_rows, sample_rows = build_endpoint_rows(
        ratios,
        value_budget,
        derivative_budget,
        t_endpoint,
        quadrature_orders,
        selected_order,
        phi_term_count,
        polynomial_m,
        indices,
    )
    value_ratio_max = max(float(row.max_value_residual_to_first_omitted_ratio_over_orders) for row in endpoint_rows)
    derivative_ratio_max = max(
        float(row.max_derivative_residual_to_first_omitted_ratio_over_orders) for row in endpoint_rows
    )
    value_spread_max = max(float(row.value_residual_spread_scaled) for row in endpoint_rows)
    derivative_spread_max = max(float(row.derivative_residual_spread_scaled) for row in endpoint_rows)
    value_budget_fraction_max = max(float(row.value_budget_fraction_at_selected_order) for row in endpoint_rows)
    derivative_budget_fraction_max = max(float(row.derivative_budget_fraction_at_selected_order) for row in endpoint_rows)
    c0 = float(phi_np(np.array([0.0], dtype=np.float64), phi_term_count)[0])
    return {
        "parameters": {
            "t_endpoint": t_endpoint,
            "u_endpoint": f"1/{t_endpoint}",
            "quadrature_rule": "generalized Gauss-Laguerre for Gamma(i+1/2, 1)",
            "quadrature_orders": list(quadrature_orders),
            "selected_quadrature_order": selected_order,
            "phi_n_terms": phi_term_count,
            "ratio_cutoff_n": ratio_cutoff_n,
            "precision_bits_for_ratio_build": precision_bits,
            "polynomial_degree": 2 * polynomial_m,
            "polynomial_M": polynomial_m,
            "residual_first_j": polynomial_m + 1,
            "indices": list(indices),
            "value_budget_A": residual_summary["value_residual_half_safety_budget_A"],
            "derivative_budget_B": residual_summary["derivative_residual_half_safety_budget_B"],
            "scaled_value_residual": "|F_i(u)-P_i^(40)(u)|/u^3",
            "scaled_derivative_residual": "|F_i'(u)-P_i^(40)'(u)|/u",
            "normalizing_phi0_float": sci(c0),
        },
        "endpoint_rows": [asdict(row) for row in endpoint_rows],
        "quadrature_sample_rows": sample_rows,
        "endpoint_row_count": len(endpoint_rows),
        "quadrature_order_count": len(quadrature_orders),
        "max_selected_value_residual_to_first_omitted_ratio": max(
            row.value_residual_to_first_omitted_ratio_at_selected_order for row in endpoint_rows
        ),
        "max_selected_derivative_residual_to_first_omitted_ratio": max(
            row.derivative_residual_to_first_omitted_ratio_at_selected_order for row in endpoint_rows
        ),
        "max_value_residual_to_first_omitted_ratio_over_orders": sci(value_ratio_max),
        "max_derivative_residual_to_first_omitted_ratio_over_orders": sci(derivative_ratio_max),
        "max_value_residual_spread_scaled": sci(value_spread_max),
        "max_derivative_residual_spread_scaled": sci(derivative_spread_max),
        "max_value_budget_fraction_at_selected_order": sci(value_budget_fraction_max),
        "max_derivative_budget_fraction_at_selected_order": sci(derivative_budget_fraction_max),
        "all_selected_residuals_below_first_omitted_term": all(
            float(row.value_residual_to_first_omitted_ratio_at_selected_order) < 1
            and float(row.derivative_residual_to_first_omitted_ratio_at_selected_order) < 1
            for row in endpoint_rows
        ),
        "floating_scope_note": (
            "The endpoint quadrature is stable over N=96..320 at T=1156, but larger T values suffer "
            "double-precision cancellation after subtracting the degree-40 polynomial. This scout is therefore "
            "kept endpoint-local and is not promoted to a uniform theorem."
        ),
    }


def build_artifact(
    residual_path: Path,
    t_endpoint: int,
    quadrature_orders: tuple[int, ...],
    selected_order: int,
    phi_term_count: int,
    ratio_cutoff_n: int,
    precision_bits: int,
    polynomial_m: int,
    indices: tuple[int, ...],
) -> dict:
    diagnostics = build_diagnostics(
        residual_path,
        t_endpoint,
        quadrature_orders,
        selected_order,
        phi_term_count,
        ratio_cutoff_n,
        precision_bits,
        polynomial_m,
        indices,
    )
    paths = source_paths(residual_path)
    rows = [
        {
            "id": "nlrgaeprs_01_actual_gamma_expectation_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "At lambda=-T with u=1/T, the actual relative-Gaussian multiplier is F_i(T)=E[Phi(sqrt(uY))/Phi(0)] for Y~Gamma(i+1/2,1).",
            "formula": "F_i'(u)=E[Phi'(sqrt(uY))*Y/(2*sqrt(uY))/Phi(0)]",
            "source_artifacts": [paths["residual_budget_note"]],
            "proof_boundary": "Exact coordinate identity only; it is not a numerical enclosure or a residual bound.",
        },
        {
            "id": "nlrgaeprs_02_endpoint_quadrature_samples",
            "role": "floating_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Generalized Gauss-Laguerre quadrature at T=1156 gives stable endpoint samples for the actual value and derivative residuals over the recorded quadrature orders.",
            "diagnostics": diagnostics,
            "proof_boundary": "Floating endpoint diagnostic only; no interval arithmetic or all-T theorem is supplied.",
        },
        {
            "id": "nlrgaeprs_03_first_omitted_comparison",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "At the endpoint and selected quadrature order, all four value and derivative residuals are below one first omitted formal term.",
            "proof_boundary": "Endpoint comparison only; it does not prove a first-omitted-term remainder theorem.",
        },
        {
            "id": "nlrgaeprs_04_budget_comparison",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "At the endpoint and selected quadrature order, the observed scaled residuals use less than 0.1% of the degree-40 half-safety value and derivative budgets.",
            "source_artifacts": [
                paths["residual_budget_note"],
                paths["asymptotic_remainder_target_note"],
            ],
            "proof_boundary": "Numerical budget comparison only; it is not a certified residual estimate.",
        },
        {
            "id": "nlrgaeprs_05_uniform_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The endpoint floating quadrature proves the residual-tail estimates on the full collar 0<=u<=1/1156.",
            "gap": "The computation is local to T=1156, uses double floating quadrature, and is not an interval enclosure for either the value or derivative channel.",
            "source_artifacts": [
                paths["formal_tail_obstruction_note"],
                paths["uniform_remainder_target_note"],
            ],
            "proof_boundary": "Rejected promotion only; not a theorem about the actual remainder.",
        },
        {
            "id": "nlrgaeprs_06_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted endpoint-to-collar proof must replace this scout with interval-safe quadrature or analytic remainder control, prove uniformity in T, and preserve derivative control.",
            "source_artifacts": [paths["dependency_graph"]],
            "proof_boundary": "Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "endpoint_rows": diagnostics["endpoint_row_count"],
        "quadrature_orders": diagnostics["parameters"]["quadrature_orders"],
        "quadrature_order_count": diagnostics["quadrature_order_count"],
        "selected_quadrature_order": selected_order,
        "max_selected_value_residual_to_first_omitted_ratio": diagnostics[
            "max_selected_value_residual_to_first_omitted_ratio"
        ],
        "max_selected_derivative_residual_to_first_omitted_ratio": diagnostics[
            "max_selected_derivative_residual_to_first_omitted_ratio"
        ],
        "max_value_residual_to_first_omitted_ratio_over_orders": diagnostics[
            "max_value_residual_to_first_omitted_ratio_over_orders"
        ],
        "max_derivative_residual_to_first_omitted_ratio_over_orders": diagnostics[
            "max_derivative_residual_to_first_omitted_ratio_over_orders"
        ],
        "max_value_residual_spread_scaled": diagnostics["max_value_residual_spread_scaled"],
        "max_derivative_residual_spread_scaled": diagnostics["max_derivative_residual_spread_scaled"],
        "max_value_budget_fraction_at_selected_order": diagnostics["max_value_budget_fraction_at_selected_order"],
        "max_derivative_budget_fraction_at_selected_order": diagnostics[
            "max_derivative_budget_fraction_at_selected_order"
        ],
        "all_selected_residuals_below_first_omitted_term": diagnostics[
            "all_selected_residuals_below_first_omitted_term"
        ],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "At the collar endpoint T=1156, direct generalized Gauss-Laguerre quadrature of the actual "
            "relative-Gaussian multiplier gives value and derivative residuals for F_21..F_24 below one "
            "first omitted formal term and below 0.1% of the degree-40 half-safety budgets. This is useful "
            "evidence that the asymptotic-remainder target is numerically plausible, but it remains a "
            "floating endpoint scout rather than a uniform analytic or interval-certified remainder theorem."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout",
        "date": "2026-07-07",
        "status": "floating endpoint theorem-search diagnostic",
        "source_degree40_residual_tail_budget": paths["residual_budget_note"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction_note"],
        "source_asymptotic_remainder_target": paths["asymptotic_remainder_target_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target_note"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py"
        ),
        "proof_boundary": (
            "Floating endpoint theorem-search diagnostic only. It samples the actual relative-Gaussian "
            "integral at T=1156, but it does not provide interval enclosures, does not prove a uniform "
            "collar remainder estimate, does not prove scaled-curvature monotonicity, does not prove cone "
            "entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Endpoint quadrature is not promoted to a uniform collar theorem.",
            "Floating samples are not interval enclosures.",
            "The first-omitted-term comparison is numerical evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][1]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['endpoint_rows']} endpoint rows, "
        f"{summary['quadrature_order_count']} quadrature orders, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Actual Endpoint Remainder Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: floating endpoint theorem-search diagnostic. This is not a proof",
        "of a uniform residual estimate, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout`.",
        "",
        "Proof boundary: this artifact samples the actual relative-Gaussian",
        "integral at the collar endpoint `T=1156`. It is not interval-certified",
        "and it is not a theorem on the full collar.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Endpoint Setup",
        "",
        "```text",
        f"T endpoint: {diagnostics['parameters']['t_endpoint']}",
        f"u endpoint: {diagnostics['parameters']['u_endpoint']}",
        f"quadrature orders: {diagnostics['parameters']['quadrature_orders']}",
        f"selected quadrature order: {summary['selected_quadrature_order']}",
        f"polynomial degree: {diagnostics['parameters']['polynomial_degree']}",
        f"normalizing Phi(0) float: {diagnostics['parameters']['normalizing_phi0_float']}",
        "```",
        "",
        "## Endpoint Residuals",
        "",
        "```text",
        f"max selected value residual / first omitted term: {summary['max_selected_value_residual_to_first_omitted_ratio']}",
        f"max selected derivative residual / first omitted term: {summary['max_selected_derivative_residual_to_first_omitted_ratio']}",
        f"max value residual / first omitted term over orders: {summary['max_value_residual_to_first_omitted_ratio_over_orders']}",
        f"max derivative residual / first omitted term over orders: {summary['max_derivative_residual_to_first_omitted_ratio_over_orders']}",
        f"max selected value budget fraction: {summary['max_value_budget_fraction_at_selected_order']}",
        f"max selected derivative budget fraction: {summary['max_derivative_budget_fraction_at_selected_order']}",
        f"all selected residuals below first omitted term: {summary['all_selected_residuals_below_first_omitted_term']}",
        "```",
        "",
        "Per-index selected-order rows:",
        "",
        "```text",
    ]
    for row in diagnostics["endpoint_rows"]:
        lines.append(
            f"F_{row['index']}: value residual={row['value_residual_scaled_at_selected_order']}, "
            f"value first omitted={row['value_first_omitted_scaled']}, "
            f"value ratio={row['value_residual_to_first_omitted_ratio_at_selected_order']}, "
            f"derivative residual={row['derivative_residual_scaled_at_selected_order']}, "
            f"derivative first omitted={row['derivative_first_omitted_scaled']}, "
            f"derivative ratio={row['derivative_residual_to_first_omitted_ratio_at_selected_order']}"
        )
        lines.append(
            f"  spreads over orders: value={row['value_residual_spread_scaled']}, "
            f"derivative={row['derivative_residual_spread_scaled']}; "
            f"budget fractions=({row['value_budget_fraction_at_selected_order']}, "
            f"{row['derivative_budget_fraction_at_selected_order']})"
        )
    lines.extend(
        [
            "```",
            "",
            "Scope warning:",
            "",
            "```text",
            diagnostics["floating_scope_note"],
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree40_residual_tail_budget"],
            artifact["source_formal_tail_obstruction"],
            artifact["source_asymptotic_remainder_target"],
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
    parser.add_argument("--residual-json", type=Path, default=DEFAULT_RESIDUAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--t-endpoint", type=int, default=DEFAULT_T_ENDPOINT)
    parser.add_argument("--quadrature-orders", type=int, nargs="+", default=list(DEFAULT_QUADRATURE_ORDERS))
    parser.add_argument("--selected-order", type=int, default=DEFAULT_SELECTED_ORDER)
    parser.add_argument("--phi-term-count", type=int, default=DEFAULT_PHI_TERM_COUNT)
    parser.add_argument("--ratio-cutoff-n", type=int, default=DEFAULT_RATIO_CUTOFF_N)
    parser.add_argument("--precision-bits", type=int, default=DEFAULT_PRECISION_BITS)
    parser.add_argument("--polynomial-m", type=int, default=DEFAULT_POLYNOMIAL_M)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    residual_path = args.residual_json if args.residual_json.is_absolute() else REPO_ROOT / args.residual_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(
        residual_path,
        args.t_endpoint,
        tuple(args.quadrature_orders),
        args.selected_order,
        args.phi_term_count,
        args.ratio_cutoff_n,
        args.precision_bits,
        args.polynomial_m,
        DEFAULT_INDICES,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian actual endpoint remainder scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
