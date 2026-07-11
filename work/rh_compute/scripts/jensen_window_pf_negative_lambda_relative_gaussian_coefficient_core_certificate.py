#!/usr/bin/env python3
"""Build a coefficient-core propagation certificate for the relative-Gaussian grid."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import REPO_ROOT  # noqa: E402
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    arb_rising,
    build_ratio_rows,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md"
)

DEFAULT_T_GRID = (1156, 1500, 2000, 5000, 10000)
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_FIRST_J = DEFAULT_POLYNOMIAL_M + 1
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PRECISION_BITS = 384
DEFAULT_RATIO_CAP = Decimal("1.0e-6")
DEFAULT_INTERVALIZATION_SOURCE_CAP = Decimal("2.0e-3")

getcontext().prec = 110


@dataclass(frozen=True)
class CorePropagationRow:
    T: int
    index: int
    polynomial_M: int
    first_j: int
    value_coefficient_scaled_radius_upper: str
    value_first_omitted_denominator_lower: str
    value_coefficient_ratio_to_first_omitted_upper: str
    derivative_coefficient_scaled_radius_upper: str
    derivative_first_omitted_denominator_lower: str
    derivative_coefficient_ratio_to_first_omitted_upper: str
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


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def dec_text(value: flint.arb, digits: int = 80, upper: bool = False) -> str:
    endpoint = value.upper() if upper else value.lower()
    return endpoint.str(digits, radius=False)


def ball_text(value: flint.arb, digits: int = 40) -> str:
    return value.str(digits)


def decimal_from_text(text: str) -> Decimal:
    return Decimal(text.replace("e", "E"))


def source_paths() -> dict[str, str]:
    return {
        "degree16_certificate_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md"
        ),
        "cancellation_reduced_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md"
        ),
        "intervalization_target_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md"
        ),
        "first_omitted_denominator_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
        ),
        "quadrature_ladder_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md"
        ),
        "node_c0_range_certificate_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md"
        ),
        "uniform_remainder_target_note": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "formal_tail_obstruction_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"
        ),
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def coefficient_row_table(max_degree: int, cutoff_n: int) -> tuple[list[dict], list[flint.arb]]:
    ratio_rows, ratios = build_ratio_rows(max_degree, cutoff_n)
    table = []
    for index, row in enumerate(ratio_rows):
        radius = ratios[index].rad()
        table.append(
            {
                **asdict(row),
                "coefficient_index": index,
                "ratio_to_c0_ball_80_digits": ratios[index].str(80),
                "radius_upper": dec_text(flint.arb(radius), digits=80, upper=True),
                "source": "build_ratio_rows finite n<=80 sum plus geometric Arb tail-radius bound",
            }
        )
    return table, ratios


def pointwise_radius_envelopes(ratios: list[flint.arb], polynomial_m: int) -> dict:
    vmax = flint.arb(809) / flint.arb(1156)
    padded = flint.arb(1)
    value_vmax = flint.arb(0)
    derivative_vmax = flint.arb(0)
    value_padded = flint.arb(0)
    derivative_padded = flint.arb(0)
    for j in range(polynomial_m + 1):
        radius = ratios[j].rad()
        value_vmax += radius * (vmax**j)
        value_padded += radius * (padded**j)
        if j:
            derivative_vmax += flint.arb(j) * radius * (vmax**j)
            derivative_padded += flint.arb(j) * radius * (padded**j)
    return {
        "vmax_source": "809/1156 from the node-c0 range certificate at order 192, F_24, T=1156",
        "vmax": "809/1156",
        "value_core_radius_on_v_le_809_1156_upper": dec_text(value_vmax, upper=True),
        "derivative_core_radius_on_v_le_809_1156_upper": dec_text(derivative_vmax, upper=True),
        "value_core_radius_on_v_le_1_upper": dec_text(value_padded, upper=True),
        "derivative_core_radius_on_v_le_1_upper": dec_text(derivative_padded, upper=True),
        "proof_boundary": (
            "Pointwise coefficient-radius envelope only; it is not a Phi-node enclosure, "
            "not a quadrature-remainder theorem, and not a residual numerator bound."
        ),
    }


def scaled_value_coefficient_radius(
    ratios: list[flint.arb],
    polynomial_m: int,
    index: int,
    u: flint.arb,
) -> flint.arb:
    total = flint.arb(0)
    for j in range(polynomial_m + 1):
        total += ratios[j].rad() * arb_rising(index, j) * (u**j) / (u**3)
    return total


def scaled_derivative_coefficient_radius(
    ratios: list[flint.arb],
    polynomial_m: int,
    index: int,
    u: flint.arb,
) -> flint.arb:
    total = flint.arb(0)
    for j in range(1, polynomial_m + 1):
        total += flint.arb(j) * ratios[j].rad() * arb_rising(index, j) * (u**j) / (u**2)
    return total


def first_omitted_denominators(
    first_ratio_abs: flint.arb,
    first_j: int,
    index: int,
    u: flint.arb,
) -> tuple[flint.arb, flint.arb]:
    rising = arb_rising(index, first_j)
    value_denominator = first_ratio_abs * rising * (u**first_j) / (u**3)
    derivative_denominator = flint.arb(first_j) * first_ratio_abs * rising * (u ** (first_j - 2))
    if not arb_positive(value_denominator) or not arb_positive(derivative_denominator):
        raise RuntimeError(f"first omitted denominator did not certify positive for F_{index}")
    return value_denominator, derivative_denominator


def build_propagation_rows(
    ratios: list[flint.arb],
    polynomial_m: int,
    first_j: int,
    t_grid: tuple[int, ...],
    indices: tuple[int, ...],
) -> tuple[list[CorePropagationRow], dict]:
    first_ratio = ratios[first_j]
    if not arb_negative(first_ratio):
        raise RuntimeError("first omitted ratio did not certify negative")
    first_ratio_abs = -first_ratio
    if not arb_positive(first_ratio_abs):
        raise RuntimeError("absolute first omitted ratio did not certify positive")

    rows: list[CorePropagationRow] = []
    max_value: tuple[Decimal, CorePropagationRow] | None = None
    max_derivative: tuple[Decimal, CorePropagationRow] | None = None
    for T in t_grid:
        u = flint.arb(1) / flint.arb(T)
        for index in indices:
            value_radius = scaled_value_coefficient_radius(ratios, polynomial_m, index, u)
            derivative_radius = scaled_derivative_coefficient_radius(ratios, polynomial_m, index, u)
            value_denominator, derivative_denominator = first_omitted_denominators(
                first_ratio_abs,
                first_j,
                index,
                u,
            )
            value_denominator_lower = flint.arb(value_denominator.lower())
            derivative_denominator_lower = flint.arb(derivative_denominator.lower())
            value_ratio = value_radius / value_denominator_lower
            derivative_ratio = derivative_radius / derivative_denominator_lower
            row = CorePropagationRow(
                T=T,
                index=index,
                polynomial_M=polynomial_m,
                first_j=first_j,
                value_coefficient_scaled_radius_upper=dec_text(value_radius, upper=True),
                value_first_omitted_denominator_lower=dec_text(value_denominator, upper=False),
                value_coefficient_ratio_to_first_omitted_upper=dec_text(value_ratio, upper=True),
                derivative_coefficient_scaled_radius_upper=dec_text(derivative_radius, upper=True),
                derivative_first_omitted_denominator_lower=dec_text(derivative_denominator, upper=False),
                derivative_coefficient_ratio_to_first_omitted_upper=dec_text(derivative_ratio, upper=True),
                proof_boundary=(
                    "Coefficient-radius propagation through exact Gamma moments only; this row does not "
                    "enclose Phi finite-node values, Laguerre nodes/weights, quadrature error, or rounding aggregation."
                ),
            )
            rows.append(row)
            value_decimal = decimal_from_text(row.value_coefficient_ratio_to_first_omitted_upper)
            derivative_decimal = decimal_from_text(row.derivative_coefficient_ratio_to_first_omitted_upper)
            if max_value is None or value_decimal > max_value[0]:
                max_value = (value_decimal, row)
            if max_derivative is None or derivative_decimal > max_derivative[0]:
                max_derivative = (derivative_decimal, row)
    assert max_value is not None
    assert max_derivative is not None
    summary = {
        "propagation_rows": len(rows),
        "t_grid": list(t_grid),
        "indices": list(indices),
        "polynomial_M": polynomial_m,
        "polynomial_degree": 2 * polynomial_m,
        "first_j": first_j,
        "value_core_formula": "Phi(sqrt(v))/Phi(0)-sum_{j=0}^{20} r_j*v^j",
        "derivative_core_formula": "sqrt(v)*Phi'(sqrt(v))/(2*Phi(0))-sum_{j=1}^{20} j*r_j*v^j",
        "scaled_value_coefficient_radius_formula": "sum_{j=0}^{20} rad(r_j)*(i+1/2)_j*u^(j-3)",
        "scaled_derivative_coefficient_radius_formula": "sum_{j=1}^{20} j*rad(r_j)*(i+1/2)_j*u^(j-2)",
        "maximum_value_coefficient_ratio_to_first_omitted_upper": f"{max_value[0]:.18E}",
        "maximum_value_coefficient_ratio_location": f"T={max_value[1].T}, F_{max_value[1].index}",
        "maximum_derivative_coefficient_ratio_to_first_omitted_upper": f"{max_derivative[0]:.18E}",
        "maximum_derivative_coefficient_ratio_location": f"T={max_derivative[1].T}, F_{max_derivative[1].index}",
        "ratio_cap": f"{DEFAULT_RATIO_CAP:.18E}",
        "intervalization_per_source_cap": f"{DEFAULT_INTERVALIZATION_SOURCE_CAP:.18E}",
        "all_coefficient_ratios_below_ratio_cap": max_value[0] < DEFAULT_RATIO_CAP
        and max_derivative[0] < DEFAULT_RATIO_CAP,
        "all_coefficient_ratios_below_intervalization_source_cap": max_value[0] < DEFAULT_INTERVALIZATION_SOURCE_CAP
        and max_derivative[0] < DEFAULT_INTERVALIZATION_SOURCE_CAP,
    }
    return rows, summary


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgccc_01_ratio_ball_source",
            role="arb_coefficient_certificate",
            readiness="available_arb",
            claim=(
                "The coefficient ratios r_0 through r_21 are rebuilt as Arb balls from the finite "
                "n<=80 coefficient sums plus the geometric Arb tail-radius bound."
            ),
            diagnostics={
                "coefficient_rows": diagnostics["coefficient_rows"],
                "coefficient_rows_count": diagnostics["coefficient_rows_count"],
                "maximum_ratio_ball_radius_upper": diagnostics["maximum_ratio_ball_radius_upper"],
                "maximum_ratio_ball_radius_index": diagnostics["maximum_ratio_ball_radius_index"],
                "pointwise_radius_envelopes": diagnostics["pointwise_radius_envelopes"],
            },
            source_artifacts=[
                paths["degree16_certificate_note"],
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            proof_boundary=(
                "Arb coefficient source only; it does not evaluate Phi at Laguerre nodes or prove a residual bound."
            ),
        ),
        MatrixRow(
            id="nlrgccc_02_value_core_gamma_propagation",
            role="coefficient_radius_propagation",
            readiness="available_for_intervalization",
            claim=(
                "The value-core coefficient-ball uncertainty propagates through the exact Gamma moments "
                "to a maximum first-omitted ratio far below the 1e-6 ratio-radius target on the recorded grid."
            ),
            diagnostics={
                "propagation_rows": diagnostics["propagation_rows"],
                "propagation_summary": diagnostics["propagation_summary"],
            },
            source_artifacts=[
                paths["cancellation_reduced_grid_note"],
                paths["first_omitted_denominator_note"],
                paths["intervalization_target_note"],
            ],
            proof_boundary=(
                "Value coefficient-radius handoff only; finite Phi node evaluation and quadrature enclosure remain open."
            ),
        ),
        MatrixRow(
            id="nlrgccc_03_derivative_core_gamma_propagation",
            role="coefficient_radius_propagation",
            readiness="available_for_intervalization",
            claim=(
                "The derivative-core coefficient-ball uncertainty propagates through the exact Gamma moments "
                "to a maximum first-omitted ratio far below the 1e-6 ratio-radius target on the recorded grid."
            ),
            diagnostics={
                "propagation_rows": diagnostics["propagation_rows"],
                "propagation_summary": diagnostics["propagation_summary"],
            },
            source_artifacts=[
                paths["cancellation_reduced_grid_note"],
                paths["first_omitted_denominator_note"],
                paths["intervalization_target_note"],
            ],
            proof_boundary=(
                "Derivative coefficient-radius handoff only; finite Phi node evaluation and quadrature enclosure remain open."
            ),
        ),
        MatrixRow(
            id="nlrgccc_04_intervalization_handoff",
            role="exact_budget_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "This certificate discharges the coefficient-ball component of nlrgit_05 for the recorded finite grid, "
                "but nlrgit_05 is still open until the signed residual numerator is enclosed with all other sources."
            ),
            diagnostics={
                "available_for_intervalization_rows": diagnostics["available_for_intervalization_rows"],
                "remaining_open_sources": [
                    "finite n<=30 Phi/Phi' node evaluation",
                    "Laguerre node and weight interval enclosures",
                    "quadrature-remainder error",
                    "rounding and cross-source aggregation",
                    "finite-grid to full-collar coverage",
                ],
            },
            source_artifacts=[
                paths["intervalization_target_note"],
                paths["phi_tail_grid_certificate_note"] if "phi_tail_grid_certificate_note" in paths else paths[
                    "node_c0_range_certificate_note"
                ],
                paths["quadrature_ladder_note"],
            ],
            proof_boundary=(
                "Budget handoff only; it is not a finite-grid interval certificate and not a uniform collar theorem."
            ),
        ),
        MatrixRow(
            id="nlrgccc_05_coefficient_only_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim=(
                "Small coefficient-ball radii alone prove the cancellation-reduced first-omitted residual comparison."
            ),
            gap=(
                "The propagated radii cover only uncertainty in the precomputed r_j balls. The actual residual "
                "numerator still needs Phi evaluation, quadrature, aggregation, and continuum coverage."
            ),
            source_artifacts=[
                paths["formal_tail_obstruction_note"],
                paths["intervalization_target_note"],
            ],
            proof_boundary=(
                "Rejected promotion only; not a first-omitted residual theorem, scaled-curvature monotonicity, "
                "cone entry, RH, or Lambda <= 0."
            ),
        ),
        MatrixRow(
            id="nlrgccc_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted finite-grid certificate must combine this coefficient-radius source with certified "
                "Phi node evaluation, Laguerre node/weight intervals, quadrature remainder, and rounding aggregation."
            ),
            source_artifacts=[
                paths["dependency_graph"],
                paths["uniform_remainder_target_note"],
                paths["intervalization_target_note"],
            ],
            proof_boundary=(
                "Proof-hygiene gate only; not a finite-grid interval certificate, not full-collar coverage, "
                "and not Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(
    polynomial_m: int = DEFAULT_POLYNOMIAL_M,
    first_j: int = DEFAULT_FIRST_J,
    ratio_cutoff_n: int = DEFAULT_RATIO_CUTOFF_N,
    precision_bits: int = DEFAULT_PRECISION_BITS,
    t_grid: tuple[int, ...] = DEFAULT_T_GRID,
    indices: tuple[int, ...] = DEFAULT_INDICES,
) -> dict:
    if first_j != polynomial_m + 1:
        raise ValueError("first_j must equal polynomial_m + 1")
    flint.ctx.prec = precision_bits
    paths = source_paths()
    paths["phi_tail_grid_certificate_note"] = (
        "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md"
    )
    coefficient_rows, ratios = coefficient_row_table(2 * first_j, ratio_cutoff_n)
    max_radius_index = max(range(len(ratios)), key=lambda index: decimal_from_text(coefficient_rows[index]["radius_upper"]))
    pointwise = pointwise_radius_envelopes(ratios, polynomial_m)
    propagation_rows, propagation_summary = build_propagation_rows(ratios, polynomial_m, first_j, t_grid, indices)
    diagnostics = {
        "parameters": {
            "polynomial_M": polynomial_m,
            "polynomial_degree": 2 * polynomial_m,
            "first_j": first_j,
            "ratio_cutoff_n": ratio_cutoff_n,
            "precision_bits": precision_bits,
            "t_grid": list(t_grid),
            "indices": list(indices),
            "gamma_moment_identity": "E[Y^j]=(i+1/2)_j for Y~Gamma(i+1/2,1)",
        },
        "coefficient_rows": coefficient_rows,
        "coefficient_rows_count": len(coefficient_rows),
        "maximum_ratio_ball_radius_upper": coefficient_rows[max_radius_index]["radius_upper"],
        "maximum_ratio_ball_radius_index": max_radius_index,
        "pointwise_radius_envelopes": pointwise,
        "propagation_rows": [asdict(row) for row in propagation_rows],
        "propagation_summary": propagation_summary,
        "available_for_intervalization_rows": 2,
        "ready_to_apply_rows": 0,
        "proof_boundary_note": (
            "This propagates coefficient-ball uncertainty only. It does not enclose the residual numerator, "
            "does not certify Laguerre nodes or weights, and does not bridge the finite grid to a real collar."
        ),
    }
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "coefficient_rows": diagnostics["coefficient_rows_count"],
        "propagation_rows": len(propagation_rows),
        "available_for_intervalization_rows": diagnostics["available_for_intervalization_rows"],
        "ready_to_apply_rows": diagnostics["ready_to_apply_rows"],
        "target_closing": False,
        "maximum_value_coefficient_ratio_to_first_omitted_upper": propagation_summary[
            "maximum_value_coefficient_ratio_to_first_omitted_upper"
        ],
        "maximum_value_coefficient_ratio_location": propagation_summary["maximum_value_coefficient_ratio_location"],
        "maximum_derivative_coefficient_ratio_to_first_omitted_upper": propagation_summary[
            "maximum_derivative_coefficient_ratio_to_first_omitted_upper"
        ],
        "maximum_derivative_coefficient_ratio_location": propagation_summary[
            "maximum_derivative_coefficient_ratio_location"
        ],
        "all_coefficient_ratios_below_ratio_cap": propagation_summary["all_coefficient_ratios_below_ratio_cap"],
        "all_coefficient_ratios_below_intervalization_source_cap": propagation_summary[
            "all_coefficient_ratios_below_intervalization_source_cap"
        ],
        "main_finding": (
            "The Arb coefficient-ratio balls r_0..r_21 propagate through the exact Gamma-moment scaling with "
            "negligible first-omitted impact on the recorded finite grid: the worst value coefficient-radius ratio "
            "is at T=10000, F_21 and the worst derivative coefficient-radius ratio is also at T=10000, F_21, "
            "both far below the 1e-6 ratio-radius target. This retires the coefficient-ball source for the "
            "finite-grid intervalization ledger, but not Phi node evaluation, quadrature, rounding, or grid-to-collar coverage."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate",
        "date": "2026-07-07",
        "status": "coefficient-core propagation certificate",
        "source_degree16_certificate": paths["degree16_certificate_note"],
        "source_cancellation_reduced_grid_scout": paths["cancellation_reduced_grid_note"],
        "source_intervalization_target": paths["intervalization_target_note"],
        "source_first_omitted_denominator_certificate": paths["first_omitted_denominator_note"],
        "source_quadrature_ladder_scout": paths["quadrature_ladder_note"],
        "source_node_c0_range_certificate": paths["node_c0_range_certificate_note"],
        "source_phi_tail_grid_certificate": paths["phi_tail_grid_certificate_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target_note"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction_note"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py"
        ),
        "proof_boundary": (
            "Coefficient-core propagation certificate only. It bounds the contribution of Arb coefficient-ratio "
            "ball radii to the scaled finite-grid first-omitted comparisons, but it does not enclose the Phi "
            "finite-node residual numerator, does not prove Laguerre node/weight intervals, does not prove a "
            "quadrature-remainder theorem, does not aggregate rounding sources, does not bridge the finite grid "
            "to a uniform collar, does not prove scaled-curvature monotonicity, does not prove cone entry, and "
            "does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Only coefficient-ratio ball uncertainty is propagated.",
            "The first-omitted denominator side is imported as a certified divisor, not as a numerator proof.",
            "The finite grid is not promoted to a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['coefficient_rows']} coefficient rows, "
        f"{summary['propagation_rows']} propagation rows, "
        f"{summary['available_for_intervalization_rows']} intervalization rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][0]["diagnostics"]
    propagation = artifact["matrix_rows"][1]["diagnostics"]["propagation_summary"]
    pointwise = diagnostics["pointwise_radius_envelopes"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Coefficient-Core Propagation Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: coefficient-core propagation certificate. This is not a proof",
        "of the first-omitted residual theorem, a finite-grid interval certificate,",
        "a uniform collar theorem, scaled-curvature monotonicity, cone entry, RH,",
        "or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate`.",
        "",
        "Proof boundary: this artifact propagates Arb coefficient-ratio ball",
        "uncertainty through exact Gamma moments for the recorded finite grid.",
        "It does not enclose finite Phi node values, Laguerre nodes or weights,",
        "quadrature error, rounding aggregation, or grid-to-collar coverage.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Coefficient Source",
        "",
        "```text",
        f"polynomial M: {propagation['polynomial_M']}",
        f"polynomial degree: {propagation['polynomial_degree']}",
        f"first omitted j: {propagation['first_j']}",
        f"coefficient rows: {summary['coefficient_rows']}",
        f"maximum ratio-ball radius upper: {diagnostics['maximum_ratio_ball_radius_upper']} at r_{diagnostics['maximum_ratio_ball_radius_index']}",
        "source: build_ratio_rows with finite n<=80 sums plus geometric Arb tail-radius bound",
        "```",
        "",
        "## Radius Envelopes",
        "",
        "```text",
        f"value core radius on v<=809/1156: {pointwise['value_core_radius_on_v_le_809_1156_upper']}",
        f"derivative core radius on v<=809/1156: {pointwise['derivative_core_radius_on_v_le_809_1156_upper']}",
        f"value core radius on v<=1: {pointwise['value_core_radius_on_v_le_1_upper']}",
        f"derivative core radius on v<=1: {pointwise['derivative_core_radius_on_v_le_1_upper']}",
        "```",
        "",
        "## Scaled Grid Handoff",
        "",
        "```text",
        propagation["scaled_value_coefficient_radius_formula"],
        propagation["scaled_derivative_coefficient_radius_formula"],
        f"maximum value coefficient ratio: {summary['maximum_value_coefficient_ratio_to_first_omitted_upper']} at {summary['maximum_value_coefficient_ratio_location']}",
        f"maximum derivative coefficient ratio: {summary['maximum_derivative_coefficient_ratio_to_first_omitted_upper']} at {summary['maximum_derivative_coefficient_ratio_location']}",
        f"ratio cap: {propagation['ratio_cap']}",
        f"intervalization per-source cap: {propagation['intervalization_per_source_cap']}",
        f"all coefficient ratios below ratio cap: {summary['all_coefficient_ratios_below_ratio_cap']}",
        "```",
        "",
        "Worst-row detail:",
        "",
        "```text",
    ]
    rows = artifact["matrix_rows"][1]["diagnostics"]["propagation_rows"]
    worst_value = max(rows, key=lambda row: Decimal(row["value_coefficient_ratio_to_first_omitted_upper"]))
    worst_derivative = max(rows, key=lambda row: Decimal(row["derivative_coefficient_ratio_to_first_omitted_upper"]))
    for row in (worst_value, worst_derivative):
        lines.append(
            f"T={row['T']}, F_{row['index']}: value ratio={row['value_coefficient_ratio_to_first_omitted_upper']}, "
            f"derivative ratio={row['derivative_coefficient_ratio_to_first_omitted_upper']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree16_certificate"],
            artifact["source_first_omitted_denominator_certificate"],
            artifact["source_intervalization_target"],
            artifact["source_quadrature_ladder_scout"],
            artifact["source_phi_tail_grid_certificate"],
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
    parser.add_argument("--polynomial-m", type=int, default=DEFAULT_POLYNOMIAL_M)
    parser.add_argument("--first-j", type=int, default=DEFAULT_FIRST_J)
    parser.add_argument("--ratio-cutoff-n", type=int, default=DEFAULT_RATIO_CUTOFF_N)
    parser.add_argument("--precision-bits", type=int, default=DEFAULT_PRECISION_BITS)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(
        polynomial_m=args.polynomial_m,
        first_j=args.first_j,
        ratio_cutoff_n=args.ratio_cutoff_n,
        precision_bits=args.precision_bits,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
