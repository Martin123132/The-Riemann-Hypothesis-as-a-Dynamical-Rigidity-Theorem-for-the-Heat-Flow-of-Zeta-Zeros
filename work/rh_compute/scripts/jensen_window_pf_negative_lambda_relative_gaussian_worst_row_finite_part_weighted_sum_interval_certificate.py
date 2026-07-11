#!/usr/bin/env python3
"""Build a worst-row finite-part weighted-sum interval certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
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
    arb_rising,
    build_ratio_rows,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate import (  # noqa: E402
    christoffel_constant,
    laguerre_taylor_interval,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate import (  # noqa: E402
    DEFAULT_INDEX,
    DEFAULT_ORDER,
    DEFAULT_T,
    REPO_ROOT,
    build_root_brackets,
)


DEFAULT_WEIGHT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.json"
)
DEFAULT_PHI_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json"
)
DEFAULT_COEFFICIENT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.json"
)
DEFAULT_FIRST_OMITTED_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
)

DEFAULT_REFINEMENT_STEPS = 120
DEFAULT_DECIMAL_PRECISION = 220
DEFAULT_ARB_PRECISION_BITS = 4096
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_FIRST_J = DEFAULT_POLYNOMIAL_M + 1
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PHI_TERM_COUNT = 30


@dataclass(frozen=True)
class FinitePartNodeRow:
    root_index: int
    node_left_endpoint: str
    node_right_endpoint: str
    node_width: str
    weight_interval: str
    relative_weight_width: str
    phi_over_c0_interval: str
    derivative_phi_core_interval: str
    proof_boundary: str


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


def dec_text(value: Decimal) -> str:
    return format(value, "f")


def arb_text(value: flint.arb, digits: int = 40) -> str:
    return value.str(digits)


def arb_mid_text(value: flint.arb, digits: int = 40) -> str:
    return value.str(digits, radius=False)


def arb_lower_text(value: flint.arb, digits: int = 40) -> str:
    return value.lower().str(digits, radius=False)


def arb_upper_text(value: flint.arb, digits: int = 40) -> str:
    return value.upper().str(digits, radius=False)


def abs_upper(value: flint.arb) -> flint.arb:
    lower_abs = abs(flint.arb(value.lower()))
    upper_abs = abs(flint.arb(value.upper()))
    return lower_abs if lower_abs > upper_abs else upper_abs


def relative_width(value: flint.arb) -> flint.arb:
    return (value.upper() - value.lower()) / abs(value.mid())


def interval_from_mid_radius(midpoint: Decimal, radius: Decimal) -> flint.arb:
    return flint.arb(f"{dec_text(midpoint)} +/- {dec_text(radius)}")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(
    weight_path: Path,
    phi_tail_path: Path,
    coefficient_path: Path,
    first_omitted_path: Path,
    interval_path: Path,
) -> dict[str, str]:
    return {
        "weight_json": weight_path.relative_to(REPO_ROOT).as_posix(),
        "weight_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.md"
        ),
        "root_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md"
        ),
        "phi_tail_json": phi_tail_path.relative_to(REPO_ROOT).as_posix(),
        "phi_tail_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
        "coefficient_json": coefficient_path.relative_to(REPO_ROOT).as_posix(),
        "coefficient_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md"
        ),
        "first_omitted_json": first_omitted_path.relative_to(REPO_ROOT).as_posix(),
        "first_omitted_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
        ),
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "quadrature_ladder_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def finite_phi(x_value: flint.arb, term_count: int = DEFAULT_PHI_TERM_COUNT) -> flint.arb:
    pi = flint.arb.pi()
    total = flint.arb(0)
    for n in range(1, term_count + 1):
        n2 = n * n
        q = pi * flint.arb(n2)
        exp9 = (flint.arb(9) * x_value).exp()
        exp5 = (flint.arb(5) * x_value).exp()
        exp4 = (flint.arb(4) * x_value).exp()
        total += (
            flint.arb(2) * pi * pi * flint.arb(n2 * n2) * exp9
            - flint.arb(3) * pi * flint.arb(n2) * exp5
        ) * (-q * exp4).exp()
    return total


def finite_phip(x_value: flint.arb, term_count: int = DEFAULT_PHI_TERM_COUNT) -> flint.arb:
    pi = flint.arb.pi()
    total = flint.arb(0)
    for n in range(1, term_count + 1):
        n2 = n * n
        q = pi * flint.arb(n2)
        exp9 = (flint.arb(9) * x_value).exp()
        exp5 = (flint.arb(5) * x_value).exp()
        exp4 = (flint.arb(4) * x_value).exp()
        prefactor = (
            flint.arb(2) * pi * pi * flint.arb(n2 * n2) * exp9
            - flint.arb(3) * pi * flint.arb(n2) * exp5
        )
        prefactor_derivative = (
            flint.arb(18) * pi * pi * flint.arb(n2 * n2) * exp9
            - flint.arb(15) * pi * flint.arb(n2) * exp5
        )
        total += (prefactor_derivative - prefactor * flint.arb(4) * q * exp4) * (-q * exp4).exp()
    return total


def polynomial_moments(
    ratios: list[flint.arb],
    index: int,
    u_value: flint.arb,
    polynomial_m: int,
) -> tuple[flint.arb, flint.arb]:
    value_polynomial = flint.arb(0)
    derivative_polynomial = flint.arb(0)
    for j in range(polynomial_m + 1):
        value_polynomial += ratios[j] * arb_rising(index, j) * (u_value**j)
    for j in range(1, polynomial_m + 1):
        derivative_polynomial += flint.arb(j) * ratios[j] * arb_rising(index, j) * (u_value**j)
    return value_polynomial, derivative_polynomial


def first_omitted_denominators(
    ratios: list[flint.arb],
    index: int,
    u_value: flint.arb,
    first_j: int,
) -> tuple[flint.arb, flint.arb]:
    first_ratio_abs = abs(ratios[first_j])
    value_denominator = first_ratio_abs * arb_rising(index, first_j) * (u_value**first_j) / (u_value**3)
    derivative_denominator = (
        flint.arb(first_j) * first_ratio_abs * arb_rising(index, first_j) * (u_value ** (first_j - 2))
    )
    return value_denominator, derivative_denominator


def build_finite_part_certificate(refinement_steps: int) -> tuple[list[FinitePartNodeRow], dict]:
    root_rows, root_summary = build_root_brackets(DEFAULT_ORDER, DEFAULT_INDEX, DEFAULT_T, refinement_steps)
    alpha = flint.arb(2 * DEFAULT_INDEX - 1) / flint.arb(2)
    constant = christoffel_constant(DEFAULT_ORDER, alpha)
    u_value = flint.arb(1) / flint.arb(DEFAULT_T)
    gamma_mass = (flint.arb(DEFAULT_INDEX) + flint.arb("0.5")).gamma()
    ratio_rows, ratios = build_ratio_rows(2 * DEFAULT_FIRST_J, DEFAULT_RATIO_CUTOFF_N)
    c0_finite = finite_phi(flint.arb(0))
    if c0_finite.contains(0) or not bool(c0_finite > 0):
        raise RuntimeError("finite Phi(0) interval did not certify positive")

    weighted_phi_sum = flint.arb(0)
    weighted_phip_sum = flint.arb(0)
    node_rows: list[FinitePartNodeRow] = []
    max_relative_weight_width: tuple[flint.arb, int] | None = None
    max_node_width: tuple[Decimal, int] | None = None

    for root_row in root_rows:
        left = Decimal(root_row.left_endpoint)
        right = Decimal(root_row.right_endpoint)
        midpoint = (left + right) / Decimal(2)
        radius = (right - left) / Decimal(2)
        node_interval = interval_from_mid_radius(midpoint, radius)
        denominator = laguerre_taylor_interval(DEFAULT_ORDER + 1, alpha, midpoint, radius)
        if denominator.contains(0):
            raise RuntimeError(f"Taylor denominator contains zero at root {root_row.root_index}")
        weight = constant * node_interval / (denominator * denominator)
        if weight.contains(0) or not bool(weight > 0):
            raise RuntimeError(f"weight did not certify positive at root {root_row.root_index}")
        v_interval = u_value * node_interval
        x_interval = v_interval.sqrt()
        phi_over_c0 = finite_phi(x_interval) / c0_finite
        derivative_phi_core = x_interval * finite_phip(x_interval) / (flint.arb(2) * c0_finite)
        weighted_phi_sum += weight * phi_over_c0
        weighted_phip_sum += weight * derivative_phi_core

        weight_rel_width = relative_width(weight)
        if max_relative_weight_width is None or weight_rel_width > max_relative_weight_width[0]:
            max_relative_weight_width = (weight_rel_width, root_row.root_index)
        node_width = right - left
        if max_node_width is None or node_width > max_node_width[0]:
            max_node_width = (node_width, root_row.root_index)

        node_rows.append(
            FinitePartNodeRow(
                root_index=root_row.root_index,
                node_left_endpoint=root_row.left_endpoint,
                node_right_endpoint=root_row.right_endpoint,
                node_width=dec_text(node_width),
                weight_interval=arb_text(weight),
                relative_weight_width=arb_mid_text(weight_rel_width),
                phi_over_c0_interval=arb_text(phi_over_c0),
                derivative_phi_core_interval=arb_text(derivative_phi_core),
                proof_boundary=(
                    "Refined worst-row node/weight/Phi finite-part row only; not quadrature remainder, "
                    "not all-row coverage, and not Lambda <= 0."
                ),
            )
        )

    assert max_relative_weight_width is not None
    assert max_node_width is not None
    phi_expectation = weighted_phi_sum / gamma_mass
    phip_expectation = weighted_phip_sum / gamma_mass
    value_polynomial, derivative_polynomial = polynomial_moments(
        ratios,
        DEFAULT_INDEX,
        u_value,
        DEFAULT_POLYNOMIAL_M,
    )
    value_residual = phi_expectation - value_polynomial
    derivative_residual = phip_expectation - derivative_polynomial
    value_scaled_abs = abs_upper(value_residual) / (u_value**3)
    derivative_scaled_abs = abs_upper(derivative_residual) / (u_value**2)
    value_first_omitted, derivative_first_omitted = first_omitted_denominators(
        ratios,
        DEFAULT_INDEX,
        u_value,
        DEFAULT_FIRST_J,
    )
    value_ratio = value_scaled_abs / value_first_omitted
    derivative_ratio = derivative_scaled_abs / derivative_first_omitted
    below_one = bool(value_ratio < 1 and derivative_ratio < 1)
    if not below_one:
        raise RuntimeError("finite-part ratios did not certify below one")

    diagnostics = {
        "quadrature_order": DEFAULT_ORDER,
        "index": DEFAULT_INDEX,
        "T": DEFAULT_T,
        "u": "1/10000",
        "alpha": f"{2 * DEFAULT_INDEX - 1}/2",
        "refinement_steps": refinement_steps,
        "decimal_precision": DEFAULT_DECIMAL_PRECISION,
        "arb_precision_bits": DEFAULT_ARB_PRECISION_BITS,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "polynomial_M": DEFAULT_POLYNOMIAL_M,
        "first_j": DEFAULT_FIRST_J,
        "node_rows": len(node_rows),
        "root_summary": root_summary,
        "widest_refined_node_width": f"{max_node_width[0]:.18E}",
        "widest_refined_node_root_index": max_node_width[1],
        "finite_c0_interval": arb_text(c0_finite),
        "finite_c0_positive": True,
        "weighted_phi_sum_interval": arb_text(weighted_phi_sum),
        "weighted_phip_sum_interval": arb_text(weighted_phip_sum),
        "phi_expectation_interval": arb_text(phi_expectation),
        "phip_expectation_interval": arb_text(phip_expectation),
        "value_polynomial_exact_moment_interval": arb_text(value_polynomial),
        "derivative_polynomial_exact_moment_interval": arb_text(derivative_polynomial),
        "value_residual_interval": arb_text(value_residual),
        "derivative_residual_interval": arb_text(derivative_residual),
        "value_residual_separated_from_zero": not value_residual.contains(0),
        "derivative_residual_separated_from_zero": not derivative_residual.contains(0),
        "value_scaled_abs_upper": arb_upper_text(value_scaled_abs),
        "derivative_scaled_abs_upper": arb_upper_text(derivative_scaled_abs),
        "value_first_omitted_denominator_interval": arb_text(value_first_omitted),
        "derivative_first_omitted_denominator_interval": arb_text(derivative_first_omitted),
        "value_ratio_to_first_omitted_interval": arb_text(value_ratio),
        "derivative_ratio_to_first_omitted_interval": arb_text(derivative_ratio),
        "value_ratio_to_first_omitted_upper": arb_upper_text(value_ratio),
        "derivative_ratio_to_first_omitted_upper": arb_upper_text(derivative_ratio),
        "both_ratios_certified_below_one": below_one,
        "maximum_relative_weight_width": arb_mid_text(max_relative_weight_width[0]),
        "maximum_relative_weight_width_root_index": max_relative_weight_width[1],
        "exact_moment_rewrite": (
            "The polynomial part is not summed against interval weights; it is subtracted via exact Gamma "
            "moments sum_j r_j (i+1/2)_j u^j, avoiding artificial interval loss from cancellation."
        ),
        "ratio_rows_source_count": len(ratio_rows),
    }
    return node_rows, diagnostics


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrfpwsic_01_exact_moment_rewrite",
            role="exact_reduction",
            readiness="available_exact",
            claim=(
                "The finite n<=30 value and derivative residuals are evaluated as a certified weighted Phi "
                "expectation minus the polynomial part through exact Gamma moments, avoiding interval loss "
                "from summing already-cancelled node cores."
            ),
            diagnostics={
                "exact_moment_rewrite": diagnostics["sum_diagnostics"]["exact_moment_rewrite"],
                "polynomial_M": diagnostics["sum_diagnostics"]["polynomial_M"],
                "first_j": diagnostics["sum_diagnostics"]["first_j"],
            },
            source_artifacts=[paths["coefficient_note"], paths["quadrature_ladder_note"]],
            proof_boundary="Exact finite-part rewrite only; not a quadrature-remainder theorem.",
        ),
        MatrixRow(
            id="nlrgwrfpwsic_02_refined_nodes_and_weights",
            role="arb_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "Refining the worst-row L_320^(41/2) node brackets to 120 bisection steps gives node and "
                "Christoffel-weight intervals tight enough for finite-part weighted summation."
            ),
            diagnostics={
                "refinement_steps": diagnostics["sum_diagnostics"]["refinement_steps"],
                "node_rows": diagnostics["sum_diagnostics"]["node_rows"],
                "widest_refined_node_width": diagnostics["sum_diagnostics"]["widest_refined_node_width"],
                "widest_refined_node_root_index": diagnostics["sum_diagnostics"][
                    "widest_refined_node_root_index"
                ],
                "maximum_relative_weight_width": diagnostics["sum_diagnostics"]["maximum_relative_weight_width"],
                "maximum_relative_weight_width_root_index": diagnostics["sum_diagnostics"][
                    "maximum_relative_weight_width_root_index"
                ],
            },
            source_artifacts=[paths["root_note"], paths["weight_note"], paths["weight_json"]],
            proof_boundary="Worst-row refined node/weight intervals only; not all recorded rows or orders.",
        ),
        MatrixRow(
            id="nlrgwrfpwsic_03_finite_phi_weighted_sum",
            role="arb_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The finite n<=30 Phi and Phi' contributions are evaluated on the refined node intervals and "
                "summed with certified Christoffel-weight intervals for the T=10000, F_21, N=320 row."
            ),
            diagnostics={
                "node_evaluation_rows": diagnostics["node_rows"],
                "finite_c0_interval": diagnostics["sum_diagnostics"]["finite_c0_interval"],
                "weighted_phi_sum_interval": diagnostics["sum_diagnostics"]["weighted_phi_sum_interval"],
                "weighted_phip_sum_interval": diagnostics["sum_diagnostics"]["weighted_phip_sum_interval"],
                "phi_expectation_interval": diagnostics["sum_diagnostics"]["phi_expectation_interval"],
                "phip_expectation_interval": diagnostics["sum_diagnostics"]["phip_expectation_interval"],
            },
            source_artifacts=[paths["phi_tail_note"], paths["interval_note"]],
            proof_boundary=(
                "Finite n<=30 weighted-sum interval only; n-tail composition and quadrature remainder remain separate."
            ),
        ),
        MatrixRow(
            id="nlrgwrfpwsic_04_first_omitted_comparison",
            role="arb_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "The scaled worst-row finite-part value and derivative residual intervals are both certified "
                "below their first omitted denominators."
            ),
            diagnostics={
                "value_residual_interval": diagnostics["sum_diagnostics"]["value_residual_interval"],
                "derivative_residual_interval": diagnostics["sum_diagnostics"]["derivative_residual_interval"],
                "value_scaled_abs_upper": diagnostics["sum_diagnostics"]["value_scaled_abs_upper"],
                "derivative_scaled_abs_upper": diagnostics["sum_diagnostics"]["derivative_scaled_abs_upper"],
                "value_first_omitted_denominator_interval": diagnostics["sum_diagnostics"][
                    "value_first_omitted_denominator_interval"
                ],
                "derivative_first_omitted_denominator_interval": diagnostics["sum_diagnostics"][
                    "derivative_first_omitted_denominator_interval"
                ],
                "value_ratio_to_first_omitted_upper": diagnostics["sum_diagnostics"][
                    "value_ratio_to_first_omitted_upper"
                ],
                "derivative_ratio_to_first_omitted_upper": diagnostics["sum_diagnostics"][
                    "derivative_ratio_to_first_omitted_upper"
                ],
                "both_ratios_certified_below_one": diagnostics["sum_diagnostics"][
                    "both_ratios_certified_below_one"
                ],
            },
            source_artifacts=[paths["first_omitted_note"], paths["coefficient_note"]],
            proof_boundary=(
                "Worst-row finite-part first-omitted comparison only; not the full residual numerator or collar theorem."
            ),
        ),
        MatrixRow(
            id="nlrgwrfpwsic_05_intervalization_handoff",
            role="intervalization_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "This retires the worst-row finite n<=30 weighted-sum source for T=10000, F_21, N=320, "
                "but the intervalization target remains open until tails, quadrature remainder, rounding, "
                "all-row coverage, and grid-to-collar coverage are aggregated."
            ),
            diagnostics={
                "retired_component": "worst-row finite n<=30 weighted Phi/Phi' sum at T=10000, F_21, N=320",
                "remaining_sources": [
                    "compose n>30 Phi/Phi'/Phi(0) tail source with this finite-part interval",
                    "quadrature-remainder theorem or interval adaptive integration beyond the N=320 sum",
                    "rounding and cross-source aggregation",
                    "all recorded rows/orders, not just the worst row",
                    "finite-grid to full-collar coverage",
                ],
            },
            source_artifacts=[paths["interval_json"], paths["interval_note"], paths["uniform_remainder_target"]],
            proof_boundary="Handoff only; not a finite-grid interval certificate or uniform collar theorem.",
        ),
        MatrixRow(
            id="nlrgwrfpwsic_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted finite-grid certificate may use this artifact only as the worst-row finite-part "
                "weighted-sum component, with all remaining intervalization sources certified separately."
            ),
            source_artifacts=[paths["dependency_graph"]],
            proof_boundary=(
                "Proof-hygiene gate only; not a quadrature-remainder theorem, not RH, and not Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(
    weight_path: Path,
    phi_tail_path: Path,
    coefficient_path: Path,
    first_omitted_path: Path,
    interval_path: Path,
    refinement_steps: int = DEFAULT_REFINEMENT_STEPS,
) -> dict:
    getcontext().prec = DEFAULT_DECIMAL_PRECISION
    flint.ctx.prec = DEFAULT_ARB_PRECISION_BITS
    paths = source_paths(weight_path, phi_tail_path, coefficient_path, first_omitted_path, interval_path)
    weight_artifact = load_json(weight_path)
    phi_tail = load_json(phi_tail_path)
    coefficient = load_json(coefficient_path)
    first_omitted = load_json(first_omitted_path)
    interval = load_json(interval_path)
    node_rows, sum_diagnostics = build_finite_part_certificate(refinement_steps)
    diagnostics = {
        "source_weight_interval_rows": weight_artifact["summary"]["weight_interval_rows"],
        "source_phi_tail_certificate_rows": phi_tail["summary"]["matrix_rows"],
        "source_coefficient_rows": coefficient["summary"]["coefficient_rows"],
        "source_first_omitted_denominator_rows": first_omitted["summary"]["denominator_rows"],
        "source_interval_obligations": interval["summary"]["obligation_rows"],
        "sum_diagnostics": sum_diagnostics,
        "node_rows": [asdict(row) for row in node_rows],
    }
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "refined_node_rows": len(node_rows),
        "interval_weight_rows": len(node_rows),
        "quadrature_order": DEFAULT_ORDER,
        "index": DEFAULT_INDEX,
        "T": DEFAULT_T,
        "refinement_steps": refinement_steps,
        "precision_bits": DEFAULT_ARB_PRECISION_BITS,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "widest_refined_node_width": sum_diagnostics["widest_refined_node_width"],
        "maximum_relative_weight_width": sum_diagnostics["maximum_relative_weight_width"],
        "value_residual_separated_from_zero": sum_diagnostics["value_residual_separated_from_zero"],
        "derivative_residual_separated_from_zero": sum_diagnostics["derivative_residual_separated_from_zero"],
        "value_ratio_to_first_omitted_upper": sum_diagnostics["value_ratio_to_first_omitted_upper"],
        "derivative_ratio_to_first_omitted_upper": sum_diagnostics["derivative_ratio_to_first_omitted_upper"],
        "below_one_ratio_count": 2 if sum_diagnostics["both_ratios_certified_below_one"] else 0,
        "both_ratios_certified_below_one": sum_diagnostics["both_ratios_certified_below_one"],
        "available_interval_certificate_rows": 3,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst-row finite-part weighted-sum interval certificate refines the T=10000, F_21, N=320 "
            "Laguerre node brackets to 120 bisection steps, evaluates the finite n<=30 Phi/Phi' contributions "
            "with Arb on those node intervals, sums them with certified Christoffel weights, and subtracts the "
            "polynomial part through exact Gamma moments. The resulting scaled value and derivative residuals "
            "are both certified below one first omitted term. This is still only the worst-row finite-part "
            "quadrature sum: n-tail composition, quadrature remainder, rounding aggregation, all-row coverage, "
            "and grid-to-collar coverage remain open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate",
        "date": "2026-07-07",
        "status": "worst-row finite-part weighted-sum interval certificate",
        "source_weight_interval_certificate": paths["weight_note"],
        "source_weight_interval_json": paths["weight_json"],
        "source_root_bracket_certificate": paths["root_note"],
        "source_phi_tail_grid_certificate": paths["phi_tail_note"],
        "source_phi_tail_grid_json": paths["phi_tail_json"],
        "source_coefficient_core_certificate": paths["coefficient_note"],
        "source_coefficient_core_json": paths["coefficient_json"],
        "source_first_omitted_denominator_certificate": paths["first_omitted_note"],
        "source_first_omitted_denominator_json": paths["first_omitted_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_quadrature_ladder_scout": paths["quadrature_ladder_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py"
        ),
        "proof_boundary": (
            "Worst-row finite-part weighted-sum interval certificate only. It certifies the finite n<=30 "
            "weighted Phi/Phi' quadrature sum for T=10000, F_21, N=320 and compares it with the first omitted "
            "term using exact polynomial moments, but it does not compose the n>30 tail source, does not prove "
            "a quadrature-remainder theorem, does not cover all recorded rows/orders, does not aggregate "
            "rounding, does not bridge the finite grid to a uniform collar, does not prove scaled-curvature "
            "monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The certificate covers only the worst row T=10000, F_21 and quadrature order N=320.",
            "Only the finite n<=30 weighted-sum component is certified here.",
            "The n>30 tail source, quadrature remainder, rounding, all-row coverage, and grid-to-collar coverage remain open.",
            "The finite row is not promoted to a full finite-grid interval certificate or a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval "
        f"certificate: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['refined_node_rows']} refined nodes, "
        f"{summary['interval_weight_rows']} interval weights, "
        f"{summary['below_one_ratio_count']} below-one ratios, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    sum_diagnostics = artifact["matrix_rows"][3]["diagnostics"]
    node_rows = artifact["matrix_rows"][2]["diagnostics"]["node_evaluation_rows"]
    sample_indices = [0, 1, 42, 99, 199, 289, 319]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Finite-Part Weighted-Sum Interval Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row finite-part weighted-sum interval certificate. This is not a proof",
        "of a quadrature-remainder theorem, a finite-grid interval certificate,",
        "a uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate`.",
        "",
        "Proof boundary: this artifact certifies only the finite `n<=30`",
        "weighted Phi/Phi' quadrature sum for `T=10000`, `F_21`, `N=320`.",
        "It does not compose the `n>30` tail source, prove quadrature",
        "remainder, cover all rows, aggregate rounding, or bridge the grid to",
        "a full collar.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Refined Inputs",
        "",
        "```text",
        f"quadrature order: {summary['quadrature_order']}",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"refinement steps: {summary['refinement_steps']}",
        f"precision bits: {summary['precision_bits']}",
        f"Phi finite terms: {summary['phi_term_count']}",
        f"widest refined node width: {summary['widest_refined_node_width']}",
        f"maximum relative weight width: {summary['maximum_relative_weight_width']}",
        "```",
        "",
        "## First-Omitted Comparison",
        "",
        "```text",
        f"value residual interval: {sum_diagnostics['value_residual_interval']}",
        f"derivative residual interval: {sum_diagnostics['derivative_residual_interval']}",
        f"value scaled abs upper: {sum_diagnostics['value_scaled_abs_upper']}",
        f"derivative scaled abs upper: {sum_diagnostics['derivative_scaled_abs_upper']}",
        f"value ratio upper: {summary['value_ratio_to_first_omitted_upper']}",
        f"derivative ratio upper: {summary['derivative_ratio_to_first_omitted_upper']}",
        f"both ratios certified below one: {summary['both_ratios_certified_below_one']}",
        "```",
        "",
        "Sample node rows:",
        "",
        "```text",
    ]
    for sample_index in sample_indices:
        row = node_rows[sample_index]
        lines.append(
            f"root {row['root_index']}: width={row['node_width']}, "
            f"rel_weight_width={row['relative_weight_width']}, "
            f"Phi/c0={row['phi_over_c0_interval']}, "
            f"xPhi'/(2c0)={row['derivative_phi_core_interval']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Remaining sources:",
            "",
            "```text",
        ]
    )
    for item in artifact["matrix_rows"][4]["diagnostics"]["remaining_sources"]:
        lines.append(item)
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_weight_interval_certificate"],
            artifact["source_root_bracket_certificate"],
            artifact["source_phi_tail_grid_certificate"],
            artifact["source_coefficient_core_certificate"],
            artifact["source_first_omitted_denominator_certificate"],
            artifact["source_intervalization_target"],
            artifact["source_quadrature_ladder_scout"],
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
    parser.add_argument("--weight-json", type=Path, default=DEFAULT_WEIGHT_JSON)
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--coefficient-json", type=Path, default=DEFAULT_COEFFICIENT_JSON)
    parser.add_argument("--first-omitted-json", type=Path, default=DEFAULT_FIRST_OMITTED_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--refinement-steps", type=int, default=DEFAULT_REFINEMENT_STEPS)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    weight_path = args.weight_json if args.weight_json.is_absolute() else REPO_ROOT / args.weight_json
    phi_tail_path = args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    coefficient_path = (
        args.coefficient_json if args.coefficient_json.is_absolute() else REPO_ROOT / args.coefficient_json
    )
    first_omitted_path = (
        args.first_omitted_json
        if args.first_omitted_json.is_absolute()
        else REPO_ROOT / args.first_omitted_json
    )
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(
        weight_path,
        phi_tail_path,
        coefficient_path,
        first_omitted_path,
        interval_path,
        args.refinement_steps,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum "
        f"interval certificate: {out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
