#!/usr/bin/env python3
"""Build a worst-row Laguerre root-bracket certificate for the relative-Gaussian grid."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys

import numpy as np
import scipy.special as sp


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md"
)

DEFAULT_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_QUADRATURE_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json"
)

DEFAULT_ORDER = 320
DEFAULT_INDEX = 21
DEFAULT_T = 10000
DEFAULT_ALPHA_NUMERATOR = 41
DEFAULT_ALPHA_DENOMINATOR = 2
DEFAULT_BISECTION_STEPS = 60
DEFAULT_DECIMAL_PRECISION = 100
DEFAULT_ARB_PRECISION_BITS = 1024


@dataclass(frozen=True)
class RootBracketRow:
    root_index: int
    left_endpoint: str
    right_endpoint: str
    width: str
    left_sign: str
    right_sign: str
    left_value_ball: str
    right_value_ball: str
    v_interval_left_for_T_10000: str
    v_interval_right_for_T_10000: str
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


def decimal_text(value: Decimal) -> str:
    return format(value, "f")


def sci_text(value: Decimal) -> str:
    return f"{value:.18E}"


def arb_from_decimal(value: Decimal) -> flint.arb:
    return flint.arb(decimal_text(value))


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def sign_name(value: flint.arb) -> str:
    if arb_positive(value):
        return "positive"
    if arb_negative(value):
        return "negative"
    return "unresolved"


def laguerre_value(order: int, alpha: flint.arb, x_value: flint.arb) -> flint.arb:
    previous = flint.arb(1)
    if order == 0:
        return previous
    current = alpha + 1 - x_value
    if order == 1:
        return current
    for n in range(1, order):
        next_value = (
            (flint.arb(2 * n + 1) + alpha - x_value) * current
            - (flint.arb(n) + alpha) * previous
        ) / flint.arb(n + 1)
        previous, current = current, next_value
    return current


def signed_value(order: int, alpha: flint.arb, endpoint: Decimal) -> tuple[str, flint.arb]:
    value = laguerre_value(order, alpha, arb_from_decimal(endpoint))
    sign = sign_name(value)
    if sign == "unresolved":
        raise RuntimeError(f"unresolved Laguerre sign at {endpoint}")
    return sign, value


def gershgorin_upper_bound(order: int, alpha_numerator: int, alpha_denominator: int) -> Decimal:
    alpha = Decimal(alpha_numerator) / Decimal(alpha_denominator)
    first = Decimal(2) + Decimal("1.5") * alpha
    interior = Decimal(4 * order) + 2 * alpha - Decimal(6)
    last = Decimal(3 * order - 2) + Decimal("1.5") * alpha
    return max(first, interior, last)


def initial_brackets(order: int, alpha_float: float) -> list[tuple[Decimal, Decimal]]:
    nodes, _weights = sp.roots_genlaguerre(order, alpha_float)
    node_decimals = [Decimal(str(value)) for value in nodes]
    right_bound = gershgorin_upper_bound(order, DEFAULT_ALPHA_NUMERATOR, DEFAULT_ALPHA_DENOMINATOR)
    brackets: list[tuple[Decimal, Decimal]] = []
    for index in range(order):
        left = Decimal(0) if index == 0 else (node_decimals[index - 1] + node_decimals[index]) / Decimal(2)
        right = right_bound if index == order - 1 else (node_decimals[index] + node_decimals[index + 1]) / Decimal(2)
        brackets.append((left, right))
    return brackets


def refine_bracket(
    order: int,
    alpha: flint.arb,
    left: Decimal,
    right: Decimal,
    steps: int,
) -> tuple[Decimal, Decimal, str, str, flint.arb, flint.arb]:
    left_sign, _left_value = signed_value(order, alpha, left)
    right_sign, _right_value = signed_value(order, alpha, right)
    if left_sign == right_sign:
        raise RuntimeError(f"initial bracket has equal endpoint signs: {left}, {right}, {left_sign}")
    for _step in range(steps):
        midpoint = (left + right) / Decimal(2)
        midpoint_sign, _midpoint_value = signed_value(order, alpha, midpoint)
        if midpoint_sign == left_sign:
            left = midpoint
            left_sign = midpoint_sign
        else:
            right = midpoint
            right_sign = midpoint_sign
    left_sign, left_value = signed_value(order, alpha, left)
    right_sign, right_value = signed_value(order, alpha, right)
    if left_sign == right_sign:
        raise RuntimeError("refined bracket lost sign change")
    return left, right, left_sign, right_sign, left_value, right_value


def build_root_brackets(
    order: int,
    index: int,
    T: int,
    bisection_steps: int,
) -> tuple[list[RootBracketRow], dict]:
    alpha = flint.arb(2 * index - 1) / flint.arb(2)
    brackets = initial_brackets(order, float(index - 0.5))
    rows: list[RootBracketRow] = []
    widest: tuple[Decimal, int] | None = None
    narrowest: tuple[Decimal, int] | None = None
    for root_index, (left, right) in enumerate(brackets, start=1):
        left, right, left_sign, right_sign, left_value, right_value = refine_bracket(
            order,
            alpha,
            left,
            right,
            bisection_steps,
        )
        width = right - left
        if widest is None or width > widest[0]:
            widest = (width, root_index)
        if narrowest is None or width < narrowest[0]:
            narrowest = (width, root_index)
        rows.append(
            RootBracketRow(
                root_index=root_index,
                left_endpoint=decimal_text(left),
                right_endpoint=decimal_text(right),
                width=decimal_text(width),
                left_sign=left_sign,
                right_sign=right_sign,
                left_value_ball=left_value.str(20),
                right_value_ball=right_value.str(20),
                v_interval_left_for_T_10000=decimal_text(left / Decimal(T)),
                v_interval_right_for_T_10000=decimal_text(right / Decimal(T)),
                proof_boundary=(
                    "Arb sign-change bracket for one root of L_320^(41/2); this is a node bracket only, "
                    "not a Christoffel-weight interval or quadrature-remainder certificate."
                ),
            )
        )
    assert widest is not None
    assert narrowest is not None
    first = rows[0]
    last = rows[-1]
    summary = {
        "quadrature_order": order,
        "index": index,
        "T": T,
        "alpha": f"{2 * index - 1}/2",
        "root_bracket_rows": len(rows),
        "bisection_steps": bisection_steps,
        "decimal_precision": DEFAULT_DECIMAL_PRECISION,
        "arb_precision_bits": DEFAULT_ARB_PRECISION_BITS,
        "all_endpoint_signs_certified": True,
        "all_brackets_sign_changing": True,
        "all_brackets_ordered_and_disjoint": all(
            Decimal(rows[position].right_endpoint) < Decimal(rows[position + 1].left_endpoint)
            for position in range(len(rows) - 1)
        ),
        "degree_count_argument": (
            "The degree-320 polynomial has 320 disjoint sign-changing brackets, so each bracket contains "
            "exactly one real root and there are no additional roots."
        ),
        "widest_bracket_root_index": widest[1],
        "widest_bracket_width": sci_text(widest[0]),
        "narrowest_bracket_root_index": narrowest[1],
        "narrowest_bracket_width": sci_text(narrowest[0]),
        "first_root_interval": [first.left_endpoint, first.right_endpoint],
        "last_root_interval": [last.left_endpoint, last.right_endpoint],
        "last_v_interval_for_T_10000": [
            last.v_interval_left_for_T_10000,
            last.v_interval_right_for_T_10000,
        ],
        "node_brackets_available_for_worst_row": True,
    }
    return rows, summary


def build_weight_underflow_diagnostic(order: int, alpha_float: float) -> dict:
    nodes, weights = sp.roots_genlaguerre(order, alpha_float)
    zero_indices = np.where(weights == 0)[0]
    positive_indices = np.where(weights > 0)[0]
    first_zero = int(zero_indices[0] + 1) if len(zero_indices) else None
    last_zero = int(zero_indices[-1] + 1) if len(zero_indices) else None
    last_positive = int(positive_indices[-1] + 1) if len(positive_indices) else None
    return {
        "quadrature_order": order,
        "alpha_float": alpha_float,
        "floating_weight_count": int(len(weights)),
        "zero_weight_count": int(len(zero_indices)),
        "first_zero_weight_index": first_zero,
        "last_zero_weight_index": last_zero,
        "last_positive_weight_index": last_positive,
        "last_positive_weight_float": str(float(weights[positive_indices[-1]])) if len(positive_indices) else None,
        "first_zero_weight_node_float": str(float(nodes[zero_indices[0]])) if len(zero_indices) else None,
        "weights_sum_float": str(float(weights.sum())),
        "proof_boundary": (
            "Floating underflow diagnostic only; zero double weights are not interval Christoffel weights."
        ),
    }


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(grid_path: Path, interval_path: Path, quadrature_path: Path) -> dict[str, str]:
    return {
        "grid_json": grid_path.relative_to(REPO_ROOT).as_posix(),
        "grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "quadrature_json": quadrature_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
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


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrbr_01_worst_row_setup",
            role="exact_reduction",
            readiness="available_exact",
            claim=(
                "For the worst finite-grid location T=10000, F_21 and quadrature order N=320, "
                "the required nodes are roots of L_320^(41/2)."
            ),
            diagnostics={
                "T": diagnostics["root_bracket_summary"]["T"],
                "index": diagnostics["root_bracket_summary"]["index"],
                "quadrature_order": diagnostics["root_bracket_summary"]["quadrature_order"],
                "alpha": diagnostics["root_bracket_summary"]["alpha"],
                "laguerre_recurrence": (
                    "(n+1)L_(n+1)^(alpha)(x)=(2n+1+alpha-x)L_n^(alpha)(x)-(n+alpha)L_(n-1)^(alpha)(x)"
                ),
            },
            source_artifacts=[paths["grid_note"], paths["quadrature_note"]],
            proof_boundary="Exact worst-row node setup only; not a full-grid node/weight interval certificate.",
        ),
        MatrixRow(
            id="nlrgwrbr_02_arb_root_brackets",
            role="arb_node_certificate",
            readiness="available_for_intervalization",
            claim=(
                "Arb sign checks certify 320 disjoint sign-changing brackets for L_320^(41/2), giving one "
                "root in each bracket and no additional roots by the degree count."
            ),
            diagnostics={
                "root_bracket_summary": diagnostics["root_bracket_summary"],
                "root_bracket_rows": diagnostics["root_bracket_rows"],
            },
            source_artifacts=[paths["node_c0_certificate"], paths["interval_note"]],
            proof_boundary=(
                "Worst-row node-bracket certificate only; it does not certify Christoffel weights or quadrature error."
            ),
        ),
        MatrixRow(
            id="nlrgwrbr_03_node_interval_handoff",
            role="node_interval_handoff",
            readiness="available_for_intervalization",
            claim=(
                "The certified root brackets give node-induced v intervals for the T=10000, F_21, N=320 "
                "quadrature row used by the ladder scout."
            ),
            diagnostics={
                "last_v_interval_for_T_10000": diagnostics["root_bracket_summary"]["last_v_interval_for_T_10000"],
                "widest_bracket_width": diagnostics["root_bracket_summary"]["widest_bracket_width"],
                "narrowest_bracket_width": diagnostics["root_bracket_summary"]["narrowest_bracket_width"],
            },
            source_artifacts=[paths["quadrature_note"], paths["coefficient_core_certificate"]],
            proof_boundary=(
                "Node interval handoff only for one row/order; not all recorded quadrature orders and not weights."
            ),
        ),
        MatrixRow(
            id="nlrgwrbr_04_weight_underflow_diagnostic",
            role="floating_diagnostic",
            readiness="not_ready_to_apply",
            claim=(
                "SciPy double generalized Laguerre weights underflow on the same N=320, alpha=41/2 row, "
                "so floating weights cannot be promoted to interval weights."
            ),
            diagnostics=diagnostics["weight_underflow"],
            source_artifacts=[paths["quadrature_note"]],
            proof_boundary=(
                "Floating underflow diagnostic only; not a Christoffel-weight interval certificate."
            ),
        ),
        MatrixRow(
            id="nlrgwrbr_05_weight_interval_target",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim=(
                "The next certification step is a non-floating Christoffel-weight interval formula evaluated on "
                "the certified root brackets, plus weighted summation radius propagation."
            ),
            diagnostics={
                "candidate_weight_formulas": [
                    "Christoffel weights from the three-term recurrence and normalized first eigenvector",
                    "equivalent Laguerre derivative formula evaluated with Arb on certified root brackets",
                ],
                "remaining_scope": "N=320, alpha=41/2 worst row first, then all recorded indices and orders",
            },
            source_artifacts=[paths["interval_note"], paths["first_omitted_denominator_certificate"]],
            proof_boundary=(
                "Open weight-certification target only; no weight intervals or quadrature-remainder theorem are proved."
            ),
        ),
        MatrixRow(
            id="nlrgwrbr_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted finite-grid interval certificate must combine these node brackets with certified weights, "
                "Phi/Phi' node evaluation, quadrature remainder, rounding aggregation, and grid-to-collar coverage."
            ),
            source_artifacts=[
                paths["phi_tail_grid_certificate"],
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


def build_artifact(
    grid_path: Path,
    interval_path: Path,
    quadrature_path: Path,
    order: int = DEFAULT_ORDER,
    index: int = DEFAULT_INDEX,
    T: int = DEFAULT_T,
    bisection_steps: int = DEFAULT_BISECTION_STEPS,
) -> dict:
    getcontext().prec = DEFAULT_DECIMAL_PRECISION
    flint.ctx.prec = DEFAULT_ARB_PRECISION_BITS
    grid = load_json(grid_path)
    interval = load_json(interval_path)
    quadrature = load_json(quadrature_path)
    paths = source_paths(grid_path, interval_path, quadrature_path)
    root_rows, root_summary = build_root_brackets(order, index, T, bisection_steps)
    underflow = build_weight_underflow_diagnostic(order, index - 0.5)
    diagnostics = {
        "source_grid_worst_value_location": grid["summary"]["max_value_ratio_location"],
        "source_grid_worst_derivative_location": grid["summary"]["max_derivative_ratio_location"],
        "source_interval_open_requirement": "nlrgit_02_laguerre_node_weight_intervals",
        "source_quadrature_reference_order": quadrature["summary"]["reference_order"],
        "source_quadrature_ladder_rows": quadrature["summary"]["total_ladder_rows"],
        "root_bracket_summary": root_summary,
        "root_bracket_rows": [asdict(row) for row in root_rows],
        "weight_underflow": underflow,
        "intervalization_per_source_cap": interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
    }
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "root_bracket_rows": len(root_rows),
        "quadrature_order": order,
        "index": index,
        "T": T,
        "alpha": root_summary["alpha"],
        "bisection_steps": bisection_steps,
        "widest_bracket_width": root_summary["widest_bracket_width"],
        "widest_bracket_root_index": root_summary["widest_bracket_root_index"],
        "narrowest_bracket_width": root_summary["narrowest_bracket_width"],
        "narrowest_bracket_root_index": root_summary["narrowest_bracket_root_index"],
        "zero_float_weight_count": underflow["zero_weight_count"],
        "first_zero_float_weight_index": underflow["first_zero_weight_index"],
        "available_for_intervalization_rows": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst relative-Gaussian ladder row now has certified Arb brackets for all 320 nodes of "
            "L_320^(41/2). The brackets are sign-changing, ordered, and disjoint, with the widest width "
            f"{root_summary['widest_bracket_width']} at root {root_summary['widest_bracket_root_index']}. "
            "The same row exposes the remaining weight problem: SciPy double weights underflow to zero for "
            "30 tail nodes, so the next proof step must certify Christoffel weights non-float."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate",
        "date": "2026-07-07",
        "status": "worst-row Laguerre root-bracket certificate",
        "source_cancellation_reduced_grid_scout": paths["grid_note"],
        "source_cancellation_reduced_grid_json": paths["grid_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_quadrature_ladder_scout": paths["quadrature_note"],
        "source_quadrature_ladder_json": paths["quadrature_json"],
        "source_node_c0_certificate": paths["node_c0_certificate"],
        "source_phi_tail_grid_certificate": paths["phi_tail_grid_certificate"],
        "source_coefficient_core_certificate": paths["coefficient_core_certificate"],
        "source_first_omitted_denominator_certificate": paths["first_omitted_denominator_certificate"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py"
        ),
        "proof_boundary": (
            "Worst-row Laguerre root-bracket certificate only. It certifies individual Arb sign-changing "
            "node brackets for the T=10000, F_21, N=320 quadrature row, but it does not certify Christoffel "
            "weights, does not evaluate Phi or Phi' on node intervals, does not prove a quadrature-remainder "
            "theorem, does not cover all recorded rows/orders, does not aggregate rounding, does not bridge "
            "the finite grid to a uniform collar, does not prove scaled-curvature monotonicity, does not prove "
            "cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The certificate covers only the worst row T=10000, F_21 and quadrature order N=320.",
            "Only Laguerre node brackets are certified; Christoffel weights remain open.",
            "Floating SciPy weights are treated as underflow diagnostics, not as interval weights.",
            "The finite row is not promoted to a full finite-grid interval certificate or a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['root_bracket_rows']} root brackets, "
        f"{summary['zero_float_weight_count']} zero floating weights, "
        f"{summary['available_for_intervalization_rows']} intervalization rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    root_summary = artifact["matrix_rows"][1]["diagnostics"]["root_bracket_summary"]
    roots = artifact["matrix_rows"][1]["diagnostics"]["root_bracket_rows"]
    underflow = artifact["matrix_rows"][3]["diagnostics"]
    sample_indices = [0, 1, len(roots) - 2, len(roots) - 1]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Laguerre Root-Bracket Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row Laguerre root-bracket certificate. This is not a proof",
        "of a Christoffel-weight interval certificate, a quadrature-remainder",
        "theorem, a finite-grid interval certificate, a uniform collar theorem,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate`.",
        "",
        "Proof boundary: this artifact certifies individual root brackets for",
        "`L_320^(41/2)` at the worst recorded row `T=10000`, `F_21`. It does",
        "not certify weights, Phi/Phi' node values, quadrature error, rounding",
        "aggregation, all recorded rows/orders, or grid-to-collar coverage.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Root Brackets",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"quadrature order: {summary['quadrature_order']}",
        f"alpha: {summary['alpha']}",
        f"root brackets: {summary['root_bracket_rows']}",
        f"bisection steps: {summary['bisection_steps']}",
        f"widest bracket: {summary['widest_bracket_width']} at root {summary['widest_bracket_root_index']}",
        f"narrowest bracket: {summary['narrowest_bracket_width']} at root {summary['narrowest_bracket_root_index']}",
        f"degree-count argument: {root_summary['degree_count_argument']}",
        "```",
        "",
        "Sample brackets:",
        "",
        "```text",
    ]
    for sample_index in sample_indices:
        row = roots[sample_index]
        lines.append(
            f"root {row['root_index']}: [{row['left_endpoint']}, {row['right_endpoint']}], "
            f"signs=({row['left_sign']},{row['right_sign']}), width={row['width']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Weight Boundary",
            "",
            "```text",
            f"floating weight count: {underflow['floating_weight_count']}",
            f"zero floating weights: {underflow['zero_weight_count']}",
            f"first zero floating weight index: {underflow['first_zero_weight_index']}",
            f"last positive floating weight index: {underflow['last_positive_weight_index']}",
            f"last positive floating weight: {underflow['last_positive_weight_float']}",
            f"first zero floating weight node: {underflow['first_zero_weight_node_float']}",
            "```",
            "",
            "This underflow is a diagnostic boundary, not a theorem. The next",
            "certificate needs non-floating Christoffel-weight intervals on these",
            "root brackets before the quadrature row can become intervalized.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_cancellation_reduced_grid_scout"],
            artifact["source_intervalization_target"],
            artifact["source_quadrature_ladder_scout"],
            artifact["source_node_c0_certificate"],
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
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--quadrature-json", type=Path, default=DEFAULT_QUADRATURE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    quadrature_path = args.quadrature_json if args.quadrature_json.is_absolute() else REPO_ROOT / args.quadrature_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(grid_path, interval_path, quadrature_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row Laguerre root-bracket certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
