#!/usr/bin/env python3
"""Build a worst-row Christoffel-weight interval certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys

import scipy.special as sp


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate import (  # noqa: E402
    DEFAULT_INDEX,
    DEFAULT_ORDER,
    DEFAULT_T,
    REPO_ROOT,
    laguerre_value,
)


DEFAULT_ROOT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.json"
)
DEFAULT_MIDPOINT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.md"
)

DEFAULT_DECIMAL_PRECISION = 120
DEFAULT_ARB_PRECISION_BITS = 1024


@dataclass(frozen=True)
class WeightIntervalRow:
    root_index: int
    node_left_endpoint: str
    node_right_endpoint: str
    denominator_interval: str
    denominator_lower: str
    denominator_upper: str
    weight_interval: str
    weight_lower: str
    weight_upper: str
    relative_weight_width: str
    scipy_float_weight: str
    scipy_underflowed_to_zero: bool
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


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(root_path: Path, midpoint_path: Path, interval_path: Path) -> dict[str, str]:
    return {
        "root_json": root_path.relative_to(REPO_ROOT).as_posix(),
        "root_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md",
        "midpoint_json": midpoint_path.relative_to(REPO_ROOT).as_posix(),
        "midpoint_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "quadrature_ladder_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
        "node_c0_certificate": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md",
        "phi_tail_grid_certificate": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
        "coefficient_core_certificate": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md"
        ),
        "first_omitted_denominator_certificate": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def bracket_midpoint_and_radius(row: dict) -> tuple[Decimal, Decimal]:
    left = Decimal(row["left_endpoint"])
    right = Decimal(row["right_endpoint"])
    midpoint = (left + right) / Decimal(2)
    radius = (right - left) / Decimal(2)
    return midpoint, radius


def christoffel_constant(order: int, alpha: flint.arb) -> flint.arb:
    return (flint.arb(order) + alpha + 1).gamma() / (
        flint.arb(order + 1).gamma() * flint.arb((order + 1) ** 2)
    )


def laguerre_taylor_interval(order: int, alpha: flint.arb, midpoint: Decimal, radius: Decimal) -> flint.arb:
    """Evaluate L_order^(alpha)(midpoint + t), |t| <= radius, by exact Taylor identity.

    The derivative identity is d^m/dx^m L_n^(alpha)(x) =
    (-1)^m L_(n-m)^(alpha+m)(x), so the degree-order Taylor series is exact.
    """

    x_mid = flint.arb(dec_text(midpoint))
    t_interval = flint.arb(f"0 +/- {dec_text(radius)}")
    total = flint.arb(0)
    t_power = flint.arb(1)
    factorial = flint.arb(1)
    for derivative_order in range(order + 1):
        if derivative_order:
            t_power *= t_interval
            factorial *= flint.arb(derivative_order)
        coefficient = laguerre_value(
            order - derivative_order,
            alpha + flint.arb(derivative_order),
            x_mid,
        ) / factorial
        if derivative_order % 2:
            coefficient = -coefficient
        total += coefficient * t_power
    return total


def direct_interval_denominator_contains_zero(row: dict, order: int, alpha: flint.arb) -> bool:
    midpoint, radius = bracket_midpoint_and_radius(row)
    x_interval = flint.arb(f"{dec_text(midpoint)} +/- {dec_text(radius)}")
    denominator = laguerre_value(order + 1, alpha, x_interval)
    return bool(denominator.contains(0))


def relative_width(value: flint.arb) -> flint.arb:
    return (value.upper() - value.lower()) / abs(value.mid())


def build_weight_interval_rows(root_rows: list[dict], order: int, index: int) -> tuple[list[WeightIntervalRow], dict]:
    alpha = flint.arb(2 * index - 1) / flint.arb(2)
    constant = christoffel_constant(order, alpha)
    _nodes, scipy_weights = sp.roots_genlaguerre(order, index - 0.5)

    weights: list[flint.arb] = []
    rows: list[WeightIntervalRow] = []
    zero_underflows_repaired = 0
    direct_obstruction_count = 0
    taylor_denominator_contains_zero = 0
    min_weight: tuple[flint.arb, int] | None = None
    max_weight: tuple[flint.arb, int] | None = None
    max_relative_weight_width: tuple[flint.arb, int] | None = None
    max_relative_denominator_width: tuple[flint.arb, int] | None = None

    for root_row in root_rows:
        root_index = int(root_row["root_index"])
        midpoint, radius = bracket_midpoint_and_radius(root_row)
        node_interval = flint.arb(f"{dec_text(midpoint)} +/- {dec_text(radius)}")
        denominator = laguerre_taylor_interval(order + 1, alpha, midpoint, radius)
        if denominator.contains(0):
            taylor_denominator_contains_zero += 1
        if direct_interval_denominator_contains_zero(root_row, order, alpha):
            direct_obstruction_count += 1
        weight = constant * node_interval / (denominator * denominator)
        if not bool(weight > 0 and not weight.contains(0)):
            raise RuntimeError(f"weight interval did not certify positivity at root {root_index}")
        weights.append(weight)
        scipy_float_weight = float(scipy_weights[root_index - 1])
        if scipy_float_weight == 0.0:
            zero_underflows_repaired += 1

        rel_weight_width = relative_width(weight)
        rel_denominator_width = relative_width(denominator)
        if min_weight is None or weight.lower() < min_weight[0].lower():
            min_weight = (weight, root_index)
        if max_weight is None or weight.upper() > max_weight[0].upper():
            max_weight = (weight, root_index)
        if max_relative_weight_width is None or rel_weight_width > max_relative_weight_width[0]:
            max_relative_weight_width = (rel_weight_width, root_index)
        if max_relative_denominator_width is None or rel_denominator_width > max_relative_denominator_width[0]:
            max_relative_denominator_width = (rel_denominator_width, root_index)

        rows.append(
            WeightIntervalRow(
                root_index=root_index,
                node_left_endpoint=root_row["left_endpoint"],
                node_right_endpoint=root_row["right_endpoint"],
                denominator_interval=arb_text(denominator),
                denominator_lower=arb_lower_text(denominator),
                denominator_upper=arb_upper_text(denominator),
                weight_interval=arb_text(weight),
                weight_lower=arb_lower_text(weight),
                weight_upper=arb_upper_text(weight),
                relative_weight_width=arb_mid_text(rel_weight_width),
                scipy_float_weight=str(scipy_float_weight),
                scipy_underflowed_to_zero=bool(scipy_float_weight == 0.0),
                proof_boundary=(
                    "Worst-row Christoffel-weight interval only; this row certifies the weight for one "
                    "root bracket, not Phi node values, quadrature error, all-row coverage, or Lambda <= 0."
                ),
            )
        )

    assert min_weight is not None
    assert max_weight is not None
    assert max_relative_weight_width is not None
    assert max_relative_denominator_width is not None
    weight_sum = sum(weights, flint.arb(0))
    mass = (alpha + 1).gamma()
    relative_mass_error = abs(weight_sum - mass) / mass
    summary = {
        "weight_interval_rows": len(rows),
        "quadrature_order": order,
        "index": index,
        "T": DEFAULT_T,
        "alpha": f"{2 * index - 1}/2",
        "denominator_formula": "Taylor enclosure for L_(N+1)^(alpha) on each certified root bracket",
        "weight_formula": "w_j=Gamma(N+alpha+1)/(Gamma(N+1)*(N+1)^2)*x_j/[L_(N+1)^(alpha)(x_j)]^2",
        "taylor_identity": "d^m/dx^m L_n^(alpha)(x)=(-1)^m L_(n-m)^(alpha+m)(x)",
        "all_taylor_denominators_separated": taylor_denominator_contains_zero == 0,
        "taylor_denominator_contains_zero_rows": taylor_denominator_contains_zero,
        "direct_interval_denominator_contains_zero_rows": direct_obstruction_count,
        "all_weight_intervals_positive": True,
        "zero_scipy_float_weights_repaired": zero_underflows_repaired,
        "minimum_weight_root_index": min_weight[1],
        "minimum_weight_interval": arb_text(min_weight[0]),
        "maximum_weight_root_index": max_weight[1],
        "maximum_weight_interval": arb_text(max_weight[0]),
        "maximum_relative_weight_width_root_index": max_relative_weight_width[1],
        "maximum_relative_weight_width": arb_mid_text(max_relative_weight_width[0]),
        "maximum_relative_denominator_width_root_index": max_relative_denominator_width[1],
        "maximum_relative_denominator_width": arb_mid_text(max_relative_denominator_width[0]),
        "weight_sum_interval": arb_text(weight_sum),
        "exact_mass_gamma_alpha_plus_1": arb_text(mass),
        "mass_contained_by_weight_sum_interval": bool(weight_sum.contains(mass)),
        "relative_weight_sum_mass_error_interval": arb_text(relative_mass_error),
        "proof_boundary": (
            "This certifies worst-row Christoffel-weight intervals only. It does not certify Phi/Phi' "
            "node evaluation, quadrature remainder, rounding aggregation, all recorded rows/orders, "
            "grid-to-collar coverage, RH, or Lambda <= 0."
        ),
    }
    return rows, summary


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwcwic_01_taylor_denominator_formula",
            role="exact_formula",
            readiness="available_exact",
            claim=(
                "The exact Taylor identity d^m L_n^(alpha)/dx^m=(-1)^m L_(n-m)^(alpha+m) "
                "gives a centered interval enclosure for L_(N+1)^(alpha) on each certified root bracket."
            ),
            diagnostics={
                "taylor_identity": diagnostics["weight_summary"]["taylor_identity"],
                "denominator_formula": diagnostics["weight_summary"]["denominator_formula"],
                "quadrature_order": diagnostics["weight_summary"]["quadrature_order"],
                "alpha": diagnostics["weight_summary"]["alpha"],
            },
            source_artifacts=[paths["root_note"], paths["midpoint_note"]],
            proof_boundary="Exact polynomial enclosure formula only; not by itself a quadrature row certificate.",
        ),
        MatrixRow(
            id="nlrgwcwic_02_taylor_denominator_intervals",
            role="arb_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "Taylor-centered Arb evaluation separates zero from every L_321^(41/2) denominator interval "
                "over the 320 certified root brackets."
            ),
            diagnostics={
                "taylor_denominator_contains_zero_rows": diagnostics["weight_summary"][
                    "taylor_denominator_contains_zero_rows"
                ],
                "direct_interval_denominator_contains_zero_rows": diagnostics["weight_summary"][
                    "direct_interval_denominator_contains_zero_rows"
                ],
                "maximum_relative_denominator_width_root_index": diagnostics["weight_summary"][
                    "maximum_relative_denominator_width_root_index"
                ],
                "maximum_relative_denominator_width": diagnostics["weight_summary"][
                    "maximum_relative_denominator_width"
                ],
            },
            source_artifacts=[paths["root_json"], paths["root_note"]],
            proof_boundary=(
                "Worst-row denominator intervals only; not Phi node evaluation, quadrature error, or all-row coverage."
            ),
        ),
        MatrixRow(
            id="nlrgwcwic_03_christoffel_weight_intervals",
            role="arb_interval_certificate",
            readiness="available_interval_certificate",
            claim=(
                "Combining the certified node brackets with the separated Taylor denominator intervals gives "
                "320 positive Christoffel-weight intervals for the T=10000, F_21, N=320 row."
            ),
            diagnostics={
                "weight_summary": diagnostics["weight_summary"],
                "weight_interval_rows": diagnostics["weight_rows"],
            },
            source_artifacts=[paths["root_note"], paths["midpoint_note"], paths["quadrature_ladder_note"]],
            proof_boundary=(
                "Worst-row weight intervals only; not weighted Phi/Phi' summation, quadrature-remainder proof, "
                "all-row intervalization, or Lambda <= 0."
            ),
        ),
        MatrixRow(
            id="nlrgwcwic_04_underflow_repair_promoted_to_interval",
            role="finite_consistency_check",
            readiness="available_interval_certificate",
            claim=(
                "The interval weights repair the 30 SciPy double tail-weight underflows in the same row, now "
                "without relying on midpoint-only evidence."
            ),
            diagnostics={
                "zero_scipy_float_weights_repaired": diagnostics["weight_summary"][
                    "zero_scipy_float_weights_repaired"
                ],
                "minimum_weight_root_index": diagnostics["weight_summary"]["minimum_weight_root_index"],
                "minimum_weight_interval": diagnostics["weight_summary"]["minimum_weight_interval"],
            },
            source_artifacts=[paths["midpoint_note"], paths["quadrature_ladder_note"]],
            proof_boundary=(
                "Underflow repair for one worst row only; this does not certify the full quadrature ladder."
            ),
        ),
        MatrixRow(
            id="nlrgwcwic_05_mass_sum_cross_check",
            role="finite_consistency_check",
            readiness="not_ready_to_apply",
            claim=(
                "The independent weight intervals sum to an interval containing Gamma(43/2), providing a "
                "mass-sum consistency check for the worst-row weight table."
            ),
            diagnostics={
                "weight_sum_interval": diagnostics["weight_summary"]["weight_sum_interval"],
                "exact_mass_gamma_alpha_plus_1": diagnostics["weight_summary"][
                    "exact_mass_gamma_alpha_plus_1"
                ],
                "mass_contained_by_weight_sum_interval": diagnostics["weight_summary"][
                    "mass_contained_by_weight_sum_interval"
                ],
                "relative_weight_sum_mass_error_interval": diagnostics["weight_summary"][
                    "relative_weight_sum_mass_error_interval"
                ],
            },
            source_artifacts=[paths["node_c0_certificate"], paths["interval_note"]],
            proof_boundary=(
                "Mass-sum consistency check only; not a quadrature-remainder theorem or finite-grid closure."
            ),
        ),
        MatrixRow(
            id="nlrgwcwic_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "The next promoted row certificate must combine these interval weights with interval Phi/Phi' "
                "node evaluation, weighted summation radii, quadrature remainder, and rounding aggregation."
            ),
            source_artifacts=[
                paths["phi_tail_grid_certificate"],
                paths["coefficient_core_certificate"],
                paths["first_omitted_denominator_certificate"],
                paths["uniform_remainder_target"],
                paths["dependency_graph"],
            ],
            proof_boundary=(
                "Proof-hygiene gate only; not a finite-grid interval certificate, not a uniform collar theorem, "
                "not RH, and not Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(root_path: Path, midpoint_path: Path, interval_path: Path) -> dict:
    getcontext().prec = DEFAULT_DECIMAL_PRECISION
    flint.ctx.prec = DEFAULT_ARB_PRECISION_BITS
    root_artifact = load_json(root_path)
    midpoint_artifact = load_json(midpoint_path)
    interval = load_json(interval_path)
    paths = source_paths(root_path, midpoint_path, interval_path)
    root_rows = root_artifact["matrix_rows"][1]["diagnostics"]["root_bracket_rows"]
    weight_rows, weight_summary = build_weight_interval_rows(root_rows, DEFAULT_ORDER, DEFAULT_INDEX)
    diagnostics = {
        "source_root_bracket_rows": root_artifact["summary"]["root_bracket_rows"],
        "source_midpoint_weight_rows": midpoint_artifact["summary"]["midpoint_weight_rows"],
        "source_interval_open_requirement": "nlrgit_02_laguerre_node_weight_intervals",
        "intervalization_per_source_cap": interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
        "weight_rows": [asdict(row) for row in weight_rows],
        "weight_summary": weight_summary,
        "precision_bits": DEFAULT_ARB_PRECISION_BITS,
    }
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "weight_interval_rows": len(weight_rows),
        "quadrature_order": DEFAULT_ORDER,
        "index": DEFAULT_INDEX,
        "T": DEFAULT_T,
        "precision_bits": DEFAULT_ARB_PRECISION_BITS,
        "taylor_denominator_contains_zero_rows": weight_summary["taylor_denominator_contains_zero_rows"],
        "direct_interval_denominator_contains_zero_rows": weight_summary[
            "direct_interval_denominator_contains_zero_rows"
        ],
        "zero_scipy_float_weights_repaired": weight_summary["zero_scipy_float_weights_repaired"],
        "minimum_weight_root_index": weight_summary["minimum_weight_root_index"],
        "minimum_weight_interval": weight_summary["minimum_weight_interval"],
        "maximum_weight_root_index": weight_summary["maximum_weight_root_index"],
        "maximum_weight_interval": weight_summary["maximum_weight_interval"],
        "maximum_relative_weight_width_root_index": weight_summary["maximum_relative_weight_width_root_index"],
        "maximum_relative_weight_width": weight_summary["maximum_relative_weight_width"],
        "weight_sum_interval": weight_summary["weight_sum_interval"],
        "mass_contained_by_weight_sum_interval": weight_summary["mass_contained_by_weight_sum_interval"],
        "relative_weight_sum_mass_error_interval": weight_summary["relative_weight_sum_mass_error_interval"],
        "available_interval_certificate_rows": 3,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst-row Christoffel-weight interval certificate replaces the midpoint-only obstruction with "
            "a Taylor-centered Arb enclosure of L_321^(41/2) on each certified L_320^(41/2) root bracket. "
            "All 320 denominator intervals avoid zero, all 320 Christoffel-weight intervals are positive, "
            "and the interval weight sum contains Gamma(43/2). This certifies the worst-row weights only; "
            "Phi/Phi' node evaluation, quadrature remainder, all-row coverage, rounding aggregation, and the "
            "grid-to-collar bridge remain open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate",
        "date": "2026-07-07",
        "status": "worst-row Christoffel-weight interval certificate",
        "source_worst_row_laguerre_root_bracket_certificate": paths["root_note"],
        "source_worst_row_laguerre_root_bracket_json": paths["root_json"],
        "source_worst_row_christoffel_weight_midpoint_scout": paths["midpoint_note"],
        "source_worst_row_christoffel_weight_midpoint_json": paths["midpoint_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_quadrature_ladder_scout": paths["quadrature_ladder_note"],
        "source_node_c0_certificate": paths["node_c0_certificate"],
        "source_phi_tail_grid_certificate": paths["phi_tail_grid_certificate"],
        "source_coefficient_core_certificate": paths["coefficient_core_certificate"],
        "source_first_omitted_denominator_certificate": paths["first_omitted_denominator_certificate"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py"
        ),
        "proof_boundary": (
            "Worst-row Christoffel-weight interval certificate only. It certifies positive interval weights "
            "for the T=10000, F_21, N=320 generalized Laguerre row using Taylor-centered Arb denominator "
            "enclosures, but it does not evaluate Phi or Phi' on node intervals, does not prove a "
            "quadrature-remainder theorem, does not cover all recorded rows/orders, does not aggregate "
            "rounding, does not bridge the finite grid to a uniform collar, does not prove scaled-curvature "
            "monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The certificate covers only the worst row T=10000, F_21 and quadrature order N=320.",
            "Weights are certified only for the existing root brackets, not for all recorded rows or orders.",
            "Phi/Phi' node values, quadrature remainder, rounding, and grid-to-collar coverage remain open.",
            "The finite row is not promoted to a full finite-grid interval certificate or a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval "
        f"certificate: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['weight_interval_rows']} interval weights, "
        f"{summary['taylor_denominator_contains_zero_rows']} Taylor denominator obstructions, "
        f"{summary['zero_scipy_float_weights_repaired']} repaired floating underflows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    weight_summary = artifact["matrix_rows"][2]["diagnostics"]["weight_summary"]
    weight_rows = artifact["matrix_rows"][2]["diagnostics"]["weight_interval_rows"]
    sample_indices = [0, 1, 42, 99, 289, 290, 319]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Christoffel-Weight Interval Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row Christoffel-weight interval certificate. This is not a proof",
        "of Phi/Phi' node interval evaluation, a quadrature-remainder theorem, a",
        "finite-grid interval certificate, a uniform collar theorem, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate`.",
        "",
        "Proof boundary: this artifact certifies positive Christoffel-weight",
        "intervals for the certified `L_320^(41/2)` root brackets in the single",
        "worst row `T=10000`, `F_21`, `N=320`. It does not certify Phi/Phi'",
        "node values, quadrature error, all recorded rows/orders, or",
        "grid-to-collar coverage.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Interval Weights",
        "",
        "```text",
        f"formula: {weight_summary['weight_formula']}",
        f"Taylor identity: {weight_summary['taylor_identity']}",
        f"precision bits: {summary['precision_bits']}",
        f"quadrature order: {summary['quadrature_order']}",
        f"index: F_{summary['index']}",
        f"interval weights: {summary['weight_interval_rows']}",
        f"Taylor denominator obstructions: {summary['taylor_denominator_contains_zero_rows']}",
        f"direct recurrence denominator obstructions: {summary['direct_interval_denominator_contains_zero_rows']}",
        f"repaired floating underflows: {summary['zero_scipy_float_weights_repaired']}",
        f"minimum weight interval: {summary['minimum_weight_interval']} at root {summary['minimum_weight_root_index']}",
        f"maximum weight interval: {summary['maximum_weight_interval']} at root {summary['maximum_weight_root_index']}",
        f"maximum relative weight width: {summary['maximum_relative_weight_width']} at root {summary['maximum_relative_weight_width_root_index']}",
        "```",
        "",
        "Sample intervals:",
        "",
        "```text",
    ]
    for sample_index in sample_indices:
        row = weight_rows[sample_index]
        lines.append(
            f"root {row['root_index']}: weight={row['weight_interval']}, "
            f"denominator={row['denominator_interval']}, "
            f"SciPy float={row['scipy_float_weight']}, underflow={row['scipy_underflowed_to_zero']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Mass Check",
            "",
            "```text",
            f"weight sum interval: {summary['weight_sum_interval']}",
            f"mass Gamma(43/2) contained: {summary['mass_contained_by_weight_sum_interval']}",
            f"relative weight-sum mass error interval: {summary['relative_weight_sum_mass_error_interval']}",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_worst_row_laguerre_root_bracket_certificate"],
            artifact["source_worst_row_christoffel_weight_midpoint_scout"],
            artifact["source_intervalization_target"],
            artifact["source_quadrature_ladder_scout"],
            artifact["source_coefficient_core_certificate"],
            artifact["source_first_omitted_denominator_certificate"],
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
    parser.add_argument("--root-json", type=Path, default=DEFAULT_ROOT_JSON)
    parser.add_argument("--midpoint-json", type=Path, default=DEFAULT_MIDPOINT_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root_path = args.root_json if args.root_json.is_absolute() else REPO_ROOT / args.root_json
    midpoint_path = args.midpoint_json if args.midpoint_json.is_absolute() else REPO_ROOT / args.midpoint_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(root_path, midpoint_path, interval_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval "
        f"certificate: {out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
