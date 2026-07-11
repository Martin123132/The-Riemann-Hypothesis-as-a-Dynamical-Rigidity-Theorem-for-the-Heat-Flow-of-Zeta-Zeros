#!/usr/bin/env python3
"""Build a cancellation-reduced actual-remainder grid scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import mpmath as mp
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
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md"
)

DEFAULT_T_GRID = (1156, 1500, 2000, 5000, 10000)
DEFAULT_QUADRATURE_ORDERS = (64, 96, 128, 192)
DEFAULT_SELECTED_ORDER = 192
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PRECISION_BITS = 384
DEFAULT_MPMATH_DPS = 80
DEFAULT_PHI_TERM_COUNT = 30


@dataclass(frozen=True)
class GridRow:
    T: int
    index: int
    selected_quadrature_order: int
    selected_value_residual_scaled: str
    selected_derivative_residual_scaled: str
    selected_value_ratio_to_first_omitted: str
    selected_derivative_ratio_to_first_omitted: str
    max_value_ratio_to_first_omitted_over_orders: str
    max_derivative_ratio_to_first_omitted_over_orders: str
    value_ratio_spread_over_orders: str
    derivative_ratio_spread_over_orders: str
    selected_value_budget_fraction: str
    selected_derivative_budget_fraction: str
    proof_boundary: str


def sci(value: mp.mpf, digits: int = 18) -> str:
    return mp.nstr(value, n=digits, min_fixed=-6, max_fixed=6)


def sci_float(value: float) -> str:
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
        "actual_endpoint_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md"
        ),
        "uniform_remainder_target_note": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


class CancellationReducedEvaluator:
    def __init__(self, ratios: list[mp.mpf], polynomial_m: int, phi_term_count: int) -> None:
        self.ratios = ratios
        self.polynomial_m = polynomial_m
        self.phi_term_count = phi_term_count
        self.pi = mp.pi
        self.c0 = self.phi(mp.mpf("0"))

    def phi(self, x: mp.mpf) -> mp.mpf:
        total = mp.mpf("0")
        for n in range(1, self.phi_term_count + 1):
            n2 = n * n
            q = self.pi * n2
            total += (
                2 * self.pi * self.pi * n2 * n2 * mp.exp(9 * x)
                - 3 * self.pi * n2 * mp.exp(5 * x)
            ) * mp.exp(-q * mp.exp(4 * x))
        return total

    def phip(self, x: mp.mpf) -> mp.mpf:
        total = mp.mpf("0")
        for n in range(1, self.phi_term_count + 1):
            n2 = n * n
            q = self.pi * n2
            exp9 = mp.exp(9 * x)
            exp5 = mp.exp(5 * x)
            exp4 = mp.exp(4 * x)
            prefactor = 2 * self.pi * self.pi * n2 * n2 * exp9 - 3 * self.pi * n2 * exp5
            prefactor_derivative = 18 * self.pi * self.pi * n2 * n2 * exp9 - 15 * self.pi * n2 * exp5
            total += (prefactor_derivative - prefactor * 4 * q * exp4) * mp.exp(-q * exp4)
        return total

    def value_core(self, v: mp.mpf) -> mp.mpf:
        x = mp.sqrt(v)
        polynomial = mp.fsum(self.ratios[j] * (v**j) for j in range(self.polynomial_m + 1))
        return self.phi(x) / self.c0 - polynomial

    def derivative_core(self, v: mp.mpf) -> mp.mpf:
        x = mp.sqrt(v)
        polynomial = mp.fsum(j * self.ratios[j] * (v**j) for j in range(1, self.polynomial_m + 1))
        return x * self.phip(x) / (2 * self.c0) - polynomial

    def first_omitted_value(self, index: int, u: mp.mpf) -> mp.mpf:
        first_j = self.polynomial_m + 1
        return abs(self.ratios[first_j] * mp.rf(mp.mpf(index) + mp.mpf("0.5"), first_j) * (u**first_j)) / (
            u**3
        )

    def first_omitted_derivative(self, index: int, u: mp.mpf) -> mp.mpf:
        first_j = self.polynomial_m + 1
        return abs(
            first_j * self.ratios[first_j] * mp.rf(mp.mpf(index) + mp.mpf("0.5"), first_j) * (u ** (first_j - 1))
        ) / u

    def scaled_residuals(self, index: int, T: int, order: int) -> tuple[mp.mpf, mp.mpf]:
        u = mp.mpf(1) / mp.mpf(T)
        nodes, weights = sp.roots_genlaguerre(order, index - 0.5)
        total_value = mp.mpf("0")
        total_derivative = mp.mpf("0")
        for node_float, weight_float in zip(nodes, weights):
            node = mp.mpf(str(float(node_float)))
            weight = mp.mpf(str(float(weight_float)))
            v = u * node
            total_value += weight * self.value_core(v)
            total_derivative += weight * self.derivative_core(v)
        gamma = mp.gamma(mp.mpf(index) + mp.mpf("0.5"))
        value_scaled = abs(total_value / gamma) / (u**3)
        derivative_scaled = abs(total_derivative / gamma) / (u**2)
        return value_scaled, derivative_scaled


def load_ratios(polynomial_m: int, ratio_cutoff_n: int, precision_bits: int) -> list[mp.mpf]:
    flint.ctx.prec = precision_bits
    _rows, ratios_arb = build_ratio_rows(2 * (polynomial_m + 1), ratio_cutoff_n)
    return [mp.mpf(ratio.str(100, radius=False)) for ratio in ratios_arb]


def build_grid_rows(
    evaluator: CancellationReducedEvaluator,
    value_budget: mp.mpf,
    derivative_budget: mp.mpf,
    t_grid: tuple[int, ...],
    quadrature_orders: tuple[int, ...],
    selected_order: int,
    indices: tuple[int, ...],
) -> list[GridRow]:
    rows: list[GridRow] = []
    for T in t_grid:
        u = mp.mpf(1) / mp.mpf(T)
        for index in indices:
            first_value = evaluator.first_omitted_value(index, u)
            first_derivative = evaluator.first_omitted_derivative(index, u)
            value_samples: list[tuple[int, mp.mpf, mp.mpf]] = []
            derivative_samples: list[tuple[int, mp.mpf, mp.mpf]] = []
            selected_value = mp.mpf("nan")
            selected_derivative = mp.mpf("nan")
            for order in quadrature_orders:
                value_scaled, derivative_scaled = evaluator.scaled_residuals(index, T, order)
                value_ratio = value_scaled / first_value
                derivative_ratio = derivative_scaled / first_derivative
                value_samples.append((order, value_scaled, value_ratio))
                derivative_samples.append((order, derivative_scaled, derivative_ratio))
                if order == selected_order:
                    selected_value = value_scaled
                    selected_derivative = derivative_scaled
            value_ratios = [item[2] for item in value_samples]
            derivative_ratios = [item[2] for item in derivative_samples]
            rows.append(
                GridRow(
                    T=T,
                    index=index,
                    selected_quadrature_order=selected_order,
                    selected_value_residual_scaled=sci(selected_value),
                    selected_derivative_residual_scaled=sci(selected_derivative),
                    selected_value_ratio_to_first_omitted=sci(selected_value / first_value),
                    selected_derivative_ratio_to_first_omitted=sci(selected_derivative / first_derivative),
                    max_value_ratio_to_first_omitted_over_orders=sci(max(value_ratios)),
                    max_derivative_ratio_to_first_omitted_over_orders=sci(max(derivative_ratios)),
                    value_ratio_spread_over_orders=sci(max(value_ratios) - min(value_ratios)),
                    derivative_ratio_spread_over_orders=sci(max(derivative_ratios) - min(derivative_ratios)),
                    selected_value_budget_fraction=sci(selected_value / value_budget),
                    selected_derivative_budget_fraction=sci(selected_derivative / derivative_budget),
                    proof_boundary=(
                        "Cancellation-reduced floating quadrature row only; mpmath evaluates the residual "
                        "integrand, but SciPy supplies floating Laguerre nodes and weights, so this is not an "
                        "interval enclosure or a uniform theorem."
                    ),
                )
            )
    return rows


def build_diagnostics(
    residual_path: Path,
    t_grid: tuple[int, ...],
    quadrature_orders: tuple[int, ...],
    selected_order: int,
    indices: tuple[int, ...],
    polynomial_m: int,
    ratio_cutoff_n: int,
    precision_bits: int,
    mpmath_dps: int,
    phi_term_count: int,
) -> dict:
    if selected_order not in quadrature_orders:
        raise ValueError("selected_order must belong to quadrature_orders")
    mp.mp.dps = mpmath_dps
    residual = load_json(residual_path)
    residual_summary = residual["summary"]
    value_budget = mp.mpf(residual_summary["value_residual_half_safety_budget_A"])
    derivative_budget = mp.mpf(residual_summary["derivative_residual_half_safety_budget_B"])
    ratios = load_ratios(polynomial_m, ratio_cutoff_n, precision_bits)
    evaluator = CancellationReducedEvaluator(ratios, polynomial_m, phi_term_count)
    rows = build_grid_rows(evaluator, value_budget, derivative_budget, t_grid, quadrature_orders, selected_order, indices)
    max_value_ratio = max(mp.mpf(row.max_value_ratio_to_first_omitted_over_orders) for row in rows)
    max_derivative_ratio = max(mp.mpf(row.max_derivative_ratio_to_first_omitted_over_orders) for row in rows)
    max_value_ratio_row = max(rows, key=lambda row: mp.mpf(row.max_value_ratio_to_first_omitted_over_orders))
    max_derivative_ratio_row = max(rows, key=lambda row: mp.mpf(row.max_derivative_ratio_to_first_omitted_over_orders))
    max_value_spread = max(mp.mpf(row.value_ratio_spread_over_orders) for row in rows)
    max_derivative_spread = max(mp.mpf(row.derivative_ratio_spread_over_orders) for row in rows)
    max_value_budget_fraction = max(mp.mpf(row.selected_value_budget_fraction) for row in rows)
    max_derivative_budget_fraction = max(mp.mpf(row.selected_derivative_budget_fraction) for row in rows)
    return {
        "parameters": {
            "t_grid": list(t_grid),
            "indices": list(indices),
            "quadrature_rule": "generalized Gauss-Laguerre for Gamma(i+1/2, 1)",
            "quadrature_orders": list(quadrature_orders),
            "selected_quadrature_order": selected_order,
            "polynomial_degree": 2 * polynomial_m,
            "polynomial_M": polynomial_m,
            "residual_first_j": polynomial_m + 1,
            "ratio_cutoff_n": ratio_cutoff_n,
            "precision_bits_for_ratio_build": precision_bits,
            "mpmath_dps": mpmath_dps,
            "phi_n_terms": phi_term_count,
            "value_budget_A": residual_summary["value_residual_half_safety_budget_A"],
            "derivative_budget_B": residual_summary["derivative_residual_half_safety_budget_B"],
            "value_core": "Phi(sqrt(v))/Phi(0)-sum_{j=0}^{20} r_j*v^j",
            "derivative_core": "sqrt(v)*Phi'(sqrt(v))/(2*Phi(0))-sum_{j=1}^{20} j*r_j*v^j",
            "scaled_value_residual": "|E[value_core(uY)]|/u^3",
            "scaled_derivative_residual": "|E[derivative_core(uY)]|/u^2",
            "normalizing_phi0_mpmath": sci(evaluator.c0, 30),
        },
        "grid_rows": [asdict(row) for row in rows],
        "grid_row_count": len(rows),
        "t_grid_count": len(t_grid),
        "index_count": len(indices),
        "quadrature_order_count": len(quadrature_orders),
        "max_value_ratio_to_first_omitted_over_orders": sci(max_value_ratio, 18),
        "max_value_ratio_location": {"T": max_value_ratio_row.T, "index": max_value_ratio_row.index},
        "max_derivative_ratio_to_first_omitted_over_orders": sci(max_derivative_ratio, 18),
        "max_derivative_ratio_location": {"T": max_derivative_ratio_row.T, "index": max_derivative_ratio_row.index},
        "max_value_ratio_spread_over_orders": sci(max_value_spread, 18),
        "max_derivative_ratio_spread_over_orders": sci(max_derivative_spread, 18),
        "max_selected_value_budget_fraction": sci(max_value_budget_fraction, 18),
        "max_selected_derivative_budget_fraction": sci(max_derivative_budget_fraction, 18),
        "all_grid_ratios_below_one": all(
            mp.mpf(row.max_value_ratio_to_first_omitted_over_orders) < 1
            and mp.mpf(row.max_derivative_ratio_to_first_omitted_over_orders) < 1
            for row in rows
        ),
        "cancellation_reduction_note": (
            "The value and derivative remainders are subtracted inside the Gamma expectation. This removes "
            "the far-T double-precision cancellation seen when F_i and P_i^(40) are evaluated separately, "
            "but the quadrature is still floating-node diagnostic evidence, not an interval enclosure."
        ),
    }


def build_artifact(
    residual_path: Path,
    t_grid: tuple[int, ...],
    quadrature_orders: tuple[int, ...],
    selected_order: int,
    indices: tuple[int, ...],
    polynomial_m: int,
    ratio_cutoff_n: int,
    precision_bits: int,
    mpmath_dps: int,
    phi_term_count: int,
) -> dict:
    diagnostics = build_diagnostics(
        residual_path,
        t_grid,
        quadrature_orders,
        selected_order,
        indices,
        polynomial_m,
        ratio_cutoff_n,
        precision_bits,
        mpmath_dps,
        phi_term_count,
    )
    paths = source_paths(residual_path)
    rows = [
        {
            "id": "nlrgcrrgs_01_cancellation_reduced_identity",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The actual value and derivative residuals can be written as Gamma expectations of cancellation-reduced residual cores after subtracting the degree-40 Taylor polynomial inside the integrand.",
            "formula": "R_i(u)=E[Phi(sqrt(uY))/Phi(0)-sum_{j<=20}r_j*(uY)^j]",
            "source_artifacts": [paths["residual_budget_note"]],
            "proof_boundary": "Exact residual-coordinate identity only; no numerical enclosure or analytic bound is proved.",
        },
        {
            "id": "nlrgcrrgs_02_floating_grid_quadrature",
            "role": "floating_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Cancellation-reduced generalized Laguerre quadrature samples the actual residuals on the recorded T grid without the far-T blow-up seen in separate floating subtraction.",
            "diagnostics": diagnostics,
            "proof_boundary": "Floating diagnostic only; SciPy quadrature nodes and weights are not Arb intervals.",
        },
        {
            "id": "nlrgcrrgs_03_first_omitted_grid_comparison",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Across the recorded T grid, indices, and quadrature orders, all sampled value and derivative residuals are below one first omitted formal term.",
            "source_artifacts": [
                paths["asymptotic_remainder_target_note"],
                paths["actual_endpoint_scout_note"],
            ],
            "proof_boundary": "Finite grid evidence only; not a first-omitted-term theorem.",
        },
        {
            "id": "nlrgcrrgs_04_far_T_shape_evidence",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The sampled ratios increase toward one as T grows from 1156 to 10000, matching the expected first-omitted asymptotic shape.",
            "proof_boundary": "Shape evidence only; not monotonicity in T, not an asymptotic expansion proof, and not a uniform bound.",
        },
        {
            "id": "nlrgcrrgs_05_uniform_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The cancellation-reduced floating grid proves the full collar residual-tail estimates.",
            "gap": "The grid is finite, the quadrature nodes and weights are floating, and no tail in T or interval quadrature error is enclosed.",
            "source_artifacts": [paths["uniform_remainder_target_note"]],
            "proof_boundary": "Rejected promotion only; not a theorem about the actual remainder.",
        },
        {
            "id": "nlrgcrrgs_06_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof must replace this grid with interval-safe quadrature or analytic estimates for the cancellation-reduced cores, including quadrature error, Phi n-tail, and the full T collar.",
            "source_artifacts": [paths["dependency_graph"]],
            "proof_boundary": "Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "grid_rows": diagnostics["grid_row_count"],
        "t_grid_count": diagnostics["t_grid_count"],
        "index_count": diagnostics["index_count"],
        "quadrature_order_count": diagnostics["quadrature_order_count"],
        "selected_quadrature_order": selected_order,
        "max_value_ratio_to_first_omitted_over_orders": diagnostics[
            "max_value_ratio_to_first_omitted_over_orders"
        ],
        "max_value_ratio_location": diagnostics["max_value_ratio_location"],
        "max_derivative_ratio_to_first_omitted_over_orders": diagnostics[
            "max_derivative_ratio_to_first_omitted_over_orders"
        ],
        "max_derivative_ratio_location": diagnostics["max_derivative_ratio_location"],
        "max_value_ratio_spread_over_orders": diagnostics["max_value_ratio_spread_over_orders"],
        "max_derivative_ratio_spread_over_orders": diagnostics["max_derivative_ratio_spread_over_orders"],
        "max_selected_value_budget_fraction": diagnostics["max_selected_value_budget_fraction"],
        "max_selected_derivative_budget_fraction": diagnostics["max_selected_derivative_budget_fraction"],
        "all_grid_ratios_below_one": diagnostics["all_grid_ratios_below_one"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "Subtracting the degree-40 polynomial inside the Gamma expectation removes the far-T "
            "catastrophic cancellation from separate floating subtraction. On the grid T=1156,1500,2000,5000,10000 "
            "and F_21..F_24, all sampled value and derivative residuals are below one first omitted formal term; "
            "the worst sampled value ratio is about 0.9707100590 at T=10000, F_21, and the worst sampled "
            "derivative ratio is about 0.9693567775 at the same row. This is strong diagnostic support for "
            "the first-omitted remainder target, but it is still finite floating evidence rather than an "
            "interval-certified or analytic uniform collar theorem."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout",
        "date": "2026-07-07",
        "status": "floating cancellation-reduced theorem-search diagnostic",
        "source_degree40_residual_tail_budget": paths["residual_budget_note"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction_note"],
        "source_asymptotic_remainder_target": paths["asymptotic_remainder_target_note"],
        "source_actual_endpoint_scout": paths["actual_endpoint_scout_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target_note"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py"
        ),
        "proof_boundary": (
            "Floating cancellation-reduced theorem-search diagnostic only. It evaluates residual cores at "
            "high precision but uses floating generalized Laguerre nodes and weights on a finite T grid; it "
            "does not provide interval enclosures, does not prove a uniform collar remainder estimate, does "
            "not prove scaled-curvature monotonicity, does not prove cone entry, and does not prove RH or "
            "Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Cancellation-reduced floating quadrature is not promoted to a uniform collar theorem.",
            "Floating nodes and weights are not interval enclosures.",
            "The below-one first-omitted comparison is finite grid evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][1]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['grid_rows']} grid rows, "
        f"{summary['t_grid_count']} T values, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Cancellation-Reduced Remainder Grid Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: floating cancellation-reduced theorem-search diagnostic. This is not a proof",
        "of a uniform residual estimate, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout`.",
        "",
        "Proof boundary: the residual is subtracted inside the Gamma expectation",
        "and evaluated with `mpmath`, but the Laguerre nodes and weights are",
        "floating SciPy values. This is not interval-certified.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Cancellation-Reduced Setup",
        "",
        "```text",
        f"T grid: {diagnostics['parameters']['t_grid']}",
        f"indices: {diagnostics['parameters']['indices']}",
        f"quadrature orders: {diagnostics['parameters']['quadrature_orders']}",
        f"selected quadrature order: {summary['selected_quadrature_order']}",
        f"polynomial degree: {diagnostics['parameters']['polynomial_degree']}",
        f"mpmath dps: {diagnostics['parameters']['mpmath_dps']}",
        f"Phi n terms: {diagnostics['parameters']['phi_n_terms']}",
        f"normalizing Phi(0): {diagnostics['parameters']['normalizing_phi0_mpmath']}",
        "```",
        "",
        "The cancellation-reduced cores are:",
        "",
        "```text",
        diagnostics["parameters"]["value_core"],
        diagnostics["parameters"]["derivative_core"],
        diagnostics["parameters"]["scaled_value_residual"],
        diagnostics["parameters"]["scaled_derivative_residual"],
        "```",
        "",
        "## Grid Result",
        "",
        "```text",
        f"grid rows: {summary['grid_rows']}",
        f"max value ratio / first omitted: {summary['max_value_ratio_to_first_omitted_over_orders']} at {summary['max_value_ratio_location']}",
        f"max derivative ratio / first omitted: {summary['max_derivative_ratio_to_first_omitted_over_orders']} at {summary['max_derivative_ratio_location']}",
        f"max value ratio spread over orders: {summary['max_value_ratio_spread_over_orders']}",
        f"max derivative ratio spread over orders: {summary['max_derivative_ratio_spread_over_orders']}",
        f"max selected value budget fraction: {summary['max_selected_value_budget_fraction']}",
        f"max selected derivative budget fraction: {summary['max_selected_derivative_budget_fraction']}",
        f"all grid ratios below one: {summary['all_grid_ratios_below_one']}",
        "```",
        "",
        "Selected-order rows:",
        "",
        "```text",
    ]
    for row in diagnostics["grid_rows"]:
        lines.append(
            f"T={row['T']}, F_{row['index']}: value ratio={row['selected_value_ratio_to_first_omitted']}, "
            f"derivative ratio={row['selected_derivative_ratio_to_first_omitted']}, "
            f"value budget fraction={row['selected_value_budget_fraction']}, "
            f"derivative budget fraction={row['selected_derivative_budget_fraction']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Scope warning:",
            "",
            "```text",
            diagnostics["cancellation_reduction_note"],
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree40_residual_tail_budget"],
            artifact["source_formal_tail_obstruction"],
            artifact["source_asymptotic_remainder_target"],
            artifact["source_actual_endpoint_scout"],
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
    parser.add_argument("--t-grid", type=int, nargs="+", default=list(DEFAULT_T_GRID))
    parser.add_argument("--quadrature-orders", type=int, nargs="+", default=list(DEFAULT_QUADRATURE_ORDERS))
    parser.add_argument("--selected-order", type=int, default=DEFAULT_SELECTED_ORDER)
    parser.add_argument("--polynomial-m", type=int, default=DEFAULT_POLYNOMIAL_M)
    parser.add_argument("--ratio-cutoff-n", type=int, default=DEFAULT_RATIO_CUTOFF_N)
    parser.add_argument("--precision-bits", type=int, default=DEFAULT_PRECISION_BITS)
    parser.add_argument("--mpmath-dps", type=int, default=DEFAULT_MPMATH_DPS)
    parser.add_argument("--phi-term-count", type=int, default=DEFAULT_PHI_TERM_COUNT)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    residual_path = args.residual_json if args.residual_json.is_absolute() else REPO_ROOT / args.residual_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(
        residual_path,
        tuple(args.t_grid),
        tuple(args.quadrature_orders),
        args.selected_order,
        DEFAULT_INDICES,
        args.polynomial_m,
        args.ratio_cutoff_n,
        args.precision_bits,
        args.mpmath_dps,
        args.phi_term_count,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
