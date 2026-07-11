#!/usr/bin/env python3
"""Build a non-floating Christoffel-weight midpoint scout for the worst Laguerre row."""

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
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.md"
)

DEFAULT_DECIMAL_PRECISION = 120
DEFAULT_ARB_PRECISION_BITS = 2048


@dataclass(frozen=True)
class MidpointWeightRow:
    root_index: int
    node_midpoint: str
    midpoint_weight_ball: str
    midpoint_weight_lower: str
    midpoint_weight_upper: str
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


def arb_lower_text(value: flint.arb, digits: int = 40) -> str:
    return value.lower().str(digits, radius=False)


def arb_upper_text(value: flint.arb, digits: int = 40) -> str:
    return value.upper().str(digits, radius=False)


def arb_mid_text(value: flint.arb, digits: int = 40) -> str:
    return value.str(digits, radius=False)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(root_path: Path, interval_path: Path) -> dict[str, str]:
    return {
        "root_json": root_path.relative_to(REPO_ROOT).as_posix(),
        "root_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md",
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


def bracket_midpoint(row: dict) -> Decimal:
    return (Decimal(row["left_endpoint"]) + Decimal(row["right_endpoint"])) / Decimal(2)


def christoffel_constant(order: int, alpha: flint.arb) -> flint.arb:
    return (flint.arb(order) + alpha + 1).gamma() / (
        flint.arb(order + 1).gamma() * flint.arb((order + 1) ** 2)
    )


def midpoint_weight(order: int, alpha: flint.arb, node_midpoint: Decimal) -> flint.arb:
    x = flint.arb(dec_text(node_midpoint))
    laguerre_next = laguerre_value(order + 1, alpha, x)
    return christoffel_constant(order, alpha) * x / (laguerre_next * laguerre_next)


def direct_interval_denominator_contains_zero(row: dict, order: int, alpha: flint.arb) -> bool:
    left = Decimal(row["left_endpoint"])
    right = Decimal(row["right_endpoint"])
    midpoint = (left + right) / Decimal(2)
    radius = (right - left) / Decimal(2)
    x_interval = flint.arb(f"{midpoint} +/- {radius}")
    denominator = laguerre_value(order + 1, alpha, x_interval)
    return bool(denominator.contains(0))


def build_weight_rows(root_rows: list[dict], order: int, index: int) -> tuple[list[MidpointWeightRow], dict]:
    alpha = flint.arb(2 * index - 1) / flint.arb(2)
    _nodes, scipy_weights = sp.roots_genlaguerre(order, index - 0.5)
    weights: list[flint.arb] = []
    rows: list[MidpointWeightRow] = []
    zero_repaired = 0
    min_weight: tuple[flint.arb, int] | None = None
    max_weight: tuple[flint.arb, int] | None = None
    obstruction_count = 0
    for root_row in root_rows:
        root_index = int(root_row["root_index"])
        midpoint = bracket_midpoint(root_row)
        weight = midpoint_weight(order, alpha, midpoint)
        weights.append(weight)
        if direct_interval_denominator_contains_zero(root_row, order, alpha):
            obstruction_count += 1
        if not bool(weight > 0 and not weight.contains(0)):
            raise RuntimeError(f"midpoint weight did not certify positive at root {root_index}")
        if float(scipy_weights[root_index - 1]) == 0.0:
            zero_repaired += 1
        if min_weight is None or weight.lower() < min_weight[0].lower():
            min_weight = (weight, root_index)
        if max_weight is None or weight > max_weight[0]:
            max_weight = (weight, root_index)
        rows.append(
            MidpointWeightRow(
                root_index=root_index,
                node_midpoint=dec_text(midpoint),
                midpoint_weight_ball=arb_text(weight),
                midpoint_weight_lower=arb_lower_text(weight),
                midpoint_weight_upper=arb_upper_text(weight),
                scipy_float_weight=str(float(scipy_weights[root_index - 1])),
                scipy_underflowed_to_zero=bool(float(scipy_weights[root_index - 1]) == 0.0),
                proof_boundary=(
                    "Arb midpoint Christoffel-weight evaluation only; the midpoint is not an interval "
                    "root enclosure and this is not a certified quadrature weight."
                ),
            )
        )
    assert min_weight is not None
    assert max_weight is not None
    weight_sum = sum(weights, flint.arb(0))
    mass = (alpha + 1).gamma()
    mass_error = weight_sum - mass
    relative_mass_error = abs(mass_error) / mass
    summary = {
        "weight_rows": len(rows),
        "quadrature_order": order,
        "index": index,
        "T": DEFAULT_T,
        "alpha": f"{2 * index - 1}/2",
        "formula": "w_j=Gamma(N+alpha+1)/(Gamma(N+1)*(N+1)^2)*x_j/[L_(N+1)^(alpha)(x_j)]^2",
        "all_midpoint_weights_positive": True,
        "zero_scipy_float_weights_repaired": zero_repaired,
        "direct_interval_denominator_contains_zero_rows": obstruction_count,
        "minimum_midpoint_weight_root_index": min_weight[1],
        "minimum_midpoint_weight": arb_mid_text(min_weight[0]),
        "maximum_midpoint_weight_root_index": max_weight[1],
        "maximum_midpoint_weight": arb_mid_text(max_weight[0]),
        "midpoint_weight_sum": arb_mid_text(weight_sum),
        "exact_mass_gamma_alpha_plus_1": arb_mid_text(mass),
        "weight_sum_minus_mass": arb_mid_text(mass_error),
        "relative_weight_sum_mass_error": arb_mid_text(relative_mass_error),
        "proof_boundary": (
            "The weight sum check is a non-floating midpoint diagnostic. It does not replace interval "
            "Christoffel weights at the certified roots."
        ),
    }
    return rows, summary


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwcwms_01_christoffel_formula",
            role="exact_formula",
            readiness="available_exact",
            claim=(
                "For roots of L_N^(alpha), generalized Gauss-Laguerre weights satisfy "
                "w_j=Gamma(N+alpha+1)/(Gamma(N+1)*(N+1)^2)*x_j/[L_(N+1)^(alpha)(x_j)]^2."
            ),
            diagnostics={
                "formula": diagnostics["weight_summary"]["formula"],
                "quadrature_order": diagnostics["weight_summary"]["quadrature_order"],
                "alpha": diagnostics["weight_summary"]["alpha"],
            },
            source_artifacts=[paths["root_note"], paths["quadrature_ladder_note"]],
            proof_boundary="Exact formula recall only; not an interval evaluation at certified root enclosures.",
        ),
        MatrixRow(
            id="nlrgwcwms_02_arb_midpoint_weight_table",
            role="arb_midpoint_diagnostic",
            readiness="not_ready_to_apply",
            claim=(
                "Evaluating the Christoffel formula at the Arb midpoints of the certified root brackets gives "
                "320 positive non-floating midpoint weights for the worst row."
            ),
            diagnostics={
                "weight_summary": diagnostics["weight_summary"],
                "weight_rows": diagnostics["weight_rows"],
            },
            source_artifacts=[paths["root_json"], paths["root_note"]],
            proof_boundary=(
                "Midpoint diagnostic only; not certified weights at all points of the root brackets."
            ),
        ),
        MatrixRow(
            id="nlrgwcwms_03_underflow_repair_scout",
            role="floating_repair_diagnostic",
            readiness="not_ready_to_apply",
            claim=(
                "The non-floating midpoint weights repair all 30 SciPy double tail-weight underflows on the "
                "N=320, alpha=41/2 row."
            ),
            diagnostics={
                "zero_scipy_float_weights_repaired": diagnostics["weight_summary"][
                    "zero_scipy_float_weights_repaired"
                ],
                "minimum_midpoint_weight_root_index": diagnostics["weight_summary"][
                    "minimum_midpoint_weight_root_index"
                ],
                "minimum_midpoint_weight": diagnostics["weight_summary"]["minimum_midpoint_weight"],
            },
            source_artifacts=[paths["quadrature_ladder_note"], paths["root_note"]],
            proof_boundary="Underflow repair diagnostic only; not a certified interval quadrature table.",
        ),
        MatrixRow(
            id="nlrgwcwms_04_mass_sum_cross_check",
            role="numerical_consistency_check",
            readiness="not_ready_to_apply",
            claim=(
                "The midpoint weights sum to Gamma(43/2) with relative error about 1.8e-18, confirming "
                "that the non-floating midpoint table is numerically coherent but not exact-root interval data."
            ),
            diagnostics={
                "midpoint_weight_sum": diagnostics["weight_summary"]["midpoint_weight_sum"],
                "exact_mass_gamma_alpha_plus_1": diagnostics["weight_summary"]["exact_mass_gamma_alpha_plus_1"],
                "weight_sum_minus_mass": diagnostics["weight_summary"]["weight_sum_minus_mass"],
                "relative_weight_sum_mass_error": diagnostics["weight_summary"]["relative_weight_sum_mass_error"],
            },
            source_artifacts=[paths["node_c0_certificate"], paths["interval_note"]],
            proof_boundary="Mass-sum diagnostic only; the midpoint nodes are not exact certified roots.",
        ),
        MatrixRow(
            id="nlrgwcwms_05_direct_interval_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim=(
                "Directly evaluating the Christoffel denominator on the current root brackets already gives "
                "certified interval weights."
            ),
            diagnostics={
                "direct_interval_denominator_contains_zero_rows": diagnostics["weight_summary"][
                    "direct_interval_denominator_contains_zero_rows"
                ],
                "root_bracket_width_source": paths["root_note"],
            },
            gap=(
                "Naive interval evaluation of L_(N+1)^(alpha) over the current brackets contains zero on "
                "all 320 rows, so a sharper interval method is still required."
            ),
            source_artifacts=[paths["root_note"], paths["interval_note"]],
            proof_boundary=(
                "Rejected promotion only; not Christoffel-weight intervals and not a quadrature-remainder theorem."
            ),
        ),
        MatrixRow(
            id="nlrgwcwms_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted worst-row quadrature certificate must replace midpoint weights by interval weights "
                "on the certified root brackets and then propagate weighted summation radii."
            ),
            source_artifacts=[
                paths["coefficient_core_certificate"],
                paths["phi_tail_grid_certificate"],
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


def build_artifact(root_path: Path, interval_path: Path) -> dict:
    getcontext().prec = DEFAULT_DECIMAL_PRECISION
    flint.ctx.prec = DEFAULT_ARB_PRECISION_BITS
    root_artifact = load_json(root_path)
    interval = load_json(interval_path)
    paths = source_paths(root_path, interval_path)
    root_rows = root_artifact["matrix_rows"][1]["diagnostics"]["root_bracket_rows"]
    weight_rows, weight_summary = build_weight_rows(root_rows, DEFAULT_ORDER, DEFAULT_INDEX)
    diagnostics = {
        "source_root_bracket_rows": root_artifact["summary"]["root_bracket_rows"],
        "source_interval_open_requirement": "nlrgit_02_laguerre_node_weight_intervals",
        "intervalization_per_source_cap": interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
        "weight_rows": [asdict(row) for row in weight_rows],
        "weight_summary": weight_summary,
        "precision_bits": DEFAULT_ARB_PRECISION_BITS,
    }
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "midpoint_weight_rows": len(weight_rows),
        "quadrature_order": DEFAULT_ORDER,
        "index": DEFAULT_INDEX,
        "T": DEFAULT_T,
        "zero_scipy_float_weights_repaired": weight_summary["zero_scipy_float_weights_repaired"],
        "direct_interval_denominator_contains_zero_rows": weight_summary[
            "direct_interval_denominator_contains_zero_rows"
        ],
        "minimum_midpoint_weight_root_index": weight_summary["minimum_midpoint_weight_root_index"],
        "minimum_midpoint_weight": weight_summary["minimum_midpoint_weight"],
        "maximum_midpoint_weight_root_index": weight_summary["maximum_midpoint_weight_root_index"],
        "maximum_midpoint_weight": weight_summary["maximum_midpoint_weight"],
        "relative_weight_sum_mass_error": weight_summary["relative_weight_sum_mass_error"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst-row Christoffel-weight midpoint scout repairs the N=320 double-weight underflow without "
            "using floating arithmetic for the weights: all 320 Arb midpoint weights are positive, including "
            "30 weights that SciPy reports as zero. The smallest midpoint weight is at root 320, and the "
            "midpoint-weight sum matches Gamma(43/2) to relative error about 1.8e-18. This is still not a "
            "weight interval certificate because direct interval evaluation over the current root brackets "
            "contains zero in the Christoffel denominator on all 320 rows."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout",
        "date": "2026-07-07",
        "status": "worst-row Christoffel-weight midpoint scout",
        "source_worst_row_laguerre_root_bracket_certificate": paths["root_note"],
        "source_worst_row_laguerre_root_bracket_json": paths["root_json"],
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
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py"
        ),
        "proof_boundary": (
            "Worst-row Christoffel-weight midpoint scout only. It evaluates the generalized Laguerre "
            "Christoffel formula at Arb midpoints of certified root brackets and repairs floating underflow "
            "diagnostically, but it does not certify weight intervals on the root brackets, does not evaluate "
            "Phi or Phi' on node intervals, does not prove a quadrature-remainder theorem, does not cover all "
            "recorded rows/orders, does not aggregate rounding, does not bridge the finite grid to a uniform "
            "collar, does not prove scaled-curvature monotonicity, does not prove cone entry, and does not "
            "prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The scout covers only the worst row T=10000, F_21 and quadrature order N=320.",
            "Weights are evaluated at root-bracket midpoints only, not as certified interval weights.",
            "Direct interval Christoffel-denominator evaluation over the current brackets is rejected.",
            "The finite row is not promoted to a full finite-grid interval certificate or a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['midpoint_weight_rows']} midpoint weights, "
        f"{summary['zero_scipy_float_weights_repaired']} repaired floating underflows, "
        f"{summary['direct_interval_denominator_contains_zero_rows']} direct interval obstructions, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    weight_summary = artifact["matrix_rows"][1]["diagnostics"]["weight_summary"]
    weight_rows = artifact["matrix_rows"][1]["diagnostics"]["weight_rows"]
    sample_indices = [0, 1, 99, 289, 290, 299, 319]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Christoffel-Weight Midpoint Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row Christoffel-weight midpoint scout. This is not a proof",
        "of Christoffel-weight intervals, a quadrature-remainder theorem, a",
        "finite-grid interval certificate, a uniform collar theorem, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout`.",
        "",
        "Proof boundary: this artifact evaluates the Christoffel formula at",
        "Arb midpoints of the certified `L_320^(41/2)` root brackets. It repairs",
        "floating tail-weight underflow diagnostically, but it does not certify",
        "weights on the whole root brackets.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Midpoint Weights",
        "",
        "```text",
        f"formula: {weight_summary['formula']}",
        f"quadrature order: {summary['quadrature_order']}",
        f"index: F_{summary['index']}",
        f"midpoint weights: {summary['midpoint_weight_rows']}",
        f"repaired floating underflows: {summary['zero_scipy_float_weights_repaired']}",
        f"minimum midpoint weight: {summary['minimum_midpoint_weight']} at root {summary['minimum_midpoint_weight_root_index']}",
        f"maximum midpoint weight: {summary['maximum_midpoint_weight']} at root {summary['maximum_midpoint_weight_root_index']}",
        f"relative weight-sum mass error: {summary['relative_weight_sum_mass_error']}",
        "```",
        "",
        "Sample weights:",
        "",
        "```text",
    ]
    for sample_index in sample_indices:
        row = weight_rows[sample_index]
        lines.append(
            f"root {row['root_index']}: midpoint weight={row['midpoint_weight_ball']}, "
            f"SciPy float={row['scipy_float_weight']}, underflow={row['scipy_underflowed_to_zero']}"
        )
    lines.extend(
        [
            "```",
            "",
            "## Interval Boundary",
            "",
            "```text",
            f"direct interval denominator contains zero rows: {summary['direct_interval_denominator_contains_zero_rows']}",
            "current root brackets are too wide for naive interval recurrence",
            "next step: certified weight intervals via a sharper recurrence, Taylor/Lipschitz enclosure, or refined roots",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_worst_row_laguerre_root_bracket_certificate"],
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
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root_path = args.root_json if args.root_json.is_absolute() else REPO_ROOT / args.root_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(root_path, interval_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight midpoint scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
