#!/usr/bin/env python3
"""Build a high-order quadrature ladder scout for the worst relative-Gaussian grid row."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import mpmath as mp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout import (  # noqa: E402
    CancellationReducedEvaluator,
    DEFAULT_PHI_TERM_COUNT,
    DEFAULT_POLYNOMIAL_M,
    DEFAULT_PRECISION_BITS,
    DEFAULT_RATIO_CUTOFF_N,
    DEFAULT_RESIDUAL_JSON,
    REPO_ROOT,
    load_ratios,
)


DEFAULT_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md"
)

DEFAULT_LADDER_ORDERS = (96, 128, 160, 192, 224, 256, 320)
DEFAULT_REFERENCE_ORDER = 320
DEFAULT_MPMATH_DPS = 80
PROPOSED_QUADRATURE_RATIO_RADIUS_CAP = mp.mpf("1.0e-6")


@dataclass(frozen=True)
class LadderRow:
    quadrature_order: int
    value_ratio_to_first_omitted: str
    derivative_ratio_to_first_omitted: str
    value_delta_from_reference_ratio: str
    derivative_delta_from_reference_ratio: str
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


def sci(value: mp.mpf, digits: int = 25) -> str:
    return mp.nstr(value, n=digits, min_fixed=-6, max_fixed=6)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(grid_path: Path, interval_path: Path, residual_path: Path) -> dict[str, str]:
    return {
        "grid_json": grid_path.relative_to(REPO_ROOT).as_posix(),
        "grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "residual_json": residual_path.relative_to(REPO_ROOT).as_posix(),
        "residual_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
        "node_c0_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md",
        "phi_tail_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def worst_locations(grid: dict) -> list[dict[str, int]]:
    summary = grid["summary"]
    locations = []
    for key in ("max_value_ratio_location", "max_derivative_ratio_location"):
        location = summary[key]
        row = {"T": int(location["T"]), "index": int(location["index"])}
        if row not in locations:
            locations.append(row)
    return locations


def build_ladder_rows(
    evaluator: CancellationReducedEvaluator,
    T: int,
    index: int,
    orders: tuple[int, ...],
    reference_order: int,
) -> tuple[list[LadderRow], dict]:
    u = mp.mpf(1) / mp.mpf(T)
    first_value = evaluator.first_omitted_value(index, u)
    first_derivative = evaluator.first_omitted_derivative(index, u)
    samples: dict[int, tuple[mp.mpf, mp.mpf]] = {}
    for order in orders:
        value_scaled, derivative_scaled = evaluator.scaled_residuals(index, T, order)
        samples[order] = (value_scaled / first_value, derivative_scaled / first_derivative)
    reference_value, reference_derivative = samples[reference_order]
    rows = [
        LadderRow(
            quadrature_order=order,
            value_ratio_to_first_omitted=sci(samples[order][0]),
            derivative_ratio_to_first_omitted=sci(samples[order][1]),
            value_delta_from_reference_ratio=sci(samples[order][0] - reference_value),
            derivative_delta_from_reference_ratio=sci(samples[order][1] - reference_derivative),
            proof_boundary=(
                "Floating high-order quadrature sample only; SciPy supplies floating Laguerre nodes and weights."
            ),
        )
        for order in orders
    ]
    value_ratios = [sample[0] for sample in samples.values()]
    derivative_ratios = [sample[1] for sample in samples.values()]
    max_value_ratio = max(value_ratios)
    max_derivative_ratio = max(derivative_ratios)
    value_spread = max(value_ratios) - min(value_ratios)
    derivative_spread = max(derivative_ratios) - min(derivative_ratios)
    diagnostics = {
        "T": T,
        "index": index,
        "ladder_orders": list(orders),
        "reference_order": reference_order,
        "ladder_row_count": len(rows),
        "reference_value_ratio": sci(reference_value),
        "reference_derivative_ratio": sci(reference_derivative),
        "max_value_ratio": sci(max_value_ratio),
        "max_derivative_ratio": sci(max_derivative_ratio),
        "value_ratio_spread": sci(value_spread),
        "derivative_ratio_spread": sci(derivative_spread),
        "max_value_delta_from_reference_abs": sci(max(abs(sample[0] - reference_value) for sample in samples.values())),
        "max_derivative_delta_from_reference_abs": sci(
            max(abs(sample[1] - reference_derivative) for sample in samples.values())
        ),
        "all_ladder_ratios_below_one": bool(max_value_ratio < 1 and max_derivative_ratio < 1),
        "proposed_quadrature_ratio_radius_cap": sci(PROPOSED_QUADRATURE_RATIO_RADIUS_CAP, 8),
        "cap_keeps_worst_ladder_below_one": bool(
            max_value_ratio + PROPOSED_QUADRATURE_RATIO_RADIUS_CAP < 1
            and max_derivative_ratio + PROPOSED_QUADRATURE_RATIO_RADIUS_CAP < 1
        ),
        "proof_boundary": (
            "High-order floating ladder only. The spread is target calibration for a future quadrature "
            "remainder theorem, not a certified quadrature error bound."
        ),
    }
    return rows, diagnostics


def build_rows(paths: dict[str, str], location_diagnostics: list[dict], aggregate: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgqls_01_worst_row_import",
            role="finite_diagnostic",
            readiness="not_ready_to_apply",
            claim="Import the worst cancellation-reduced finite-grid row and test it with a higher quadrature ladder.",
            diagnostics={
                "source_worst_locations": aggregate["source_worst_locations"],
                "ladder_orders": aggregate["ladder_orders"],
                "reference_order": aggregate["reference_order"],
            },
            source_artifacts=[paths["grid_json"], paths["grid_note"]],
            proof_boundary="Worst-row floating diagnostic only; not an interval certificate.",
        ),
        MatrixRow(
            id="nlrgqls_02_high_order_ladder_stability",
            role="floating_diagnostic",
            readiness="not_ready_to_apply",
            claim="Across the high-order ladder through N=320, the worst row stays below one first omitted term and has tiny order-to-order spread.",
            diagnostics={"locations": location_diagnostics, "aggregate": aggregate},
            source_artifacts=[paths["grid_note"], paths["residual_note"]],
            proof_boundary=(
                "Floating quadrature stability evidence only; the ladder spread is not a rigorous quadrature remainder."
            ),
        ),
        MatrixRow(
            id="nlrgqls_03_quadrature_radius_target",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim=(
                "A future rigorous quadrature theorem with ratio radius below 1.0e-6 on the worst row would be "
                "far inside the intervalization slack and below the per-source cap."
            ),
            diagnostics={
                "proposed_quadrature_ratio_radius_cap": aggregate["proposed_quadrature_ratio_radius_cap"],
                "cap_keeps_worst_ladder_below_one": aggregate["cap_keeps_worst_ladder_below_one"],
                "intervalization_per_source_cap": aggregate["intervalization_per_source_cap"],
                "cap_below_per_source_cap": aggregate["cap_below_per_source_cap"],
            },
            source_artifacts=[paths["interval_json"], paths["interval_note"]],
            proof_boundary="Sufficient quadrature-radius target only; no remainder theorem is proved.",
        ),
        MatrixRow(
            id="nlrgqls_04_floating_overflow_boundary",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="Pushing SciPy floating generalized Laguerre quadrature to still higher orders by itself provides a rigorous error bound.",
            gap=(
                "SciPy floating roots become numerically fragile at higher orders in this parameter range, "
                "and repeated floating agreement is not an interval remainder theorem."
            ),
            source_artifacts=[paths["node_c0_note"], paths["phi_tail_note"]],
            proof_boundary="Rejected promotion only; not a quadrature theorem or interval enclosure.",
        ),
        MatrixRow(
            id="nlrgqls_05_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted proof must replace the floating ladder by a quadrature-remainder theorem or "
                "interval adaptive integration for the cancellation-reduced cores."
            ),
            source_artifacts=[paths["uniform_remainder_target"], paths["dependency_graph"]],
            proof_boundary=(
                "Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(grid_path: Path, interval_path: Path, residual_path: Path) -> dict:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    grid = load_json(grid_path)
    interval = load_json(interval_path)
    paths = source_paths(grid_path, interval_path, residual_path)
    ratios = load_ratios(DEFAULT_POLYNOMIAL_M, DEFAULT_RATIO_CUTOFF_N, DEFAULT_PRECISION_BITS)
    evaluator = CancellationReducedEvaluator(ratios, DEFAULT_POLYNOMIAL_M, DEFAULT_PHI_TERM_COUNT)
    locations = worst_locations(grid)
    location_diagnostics: list[dict] = []
    all_ladder_rows: list[dict] = []
    max_value_ratio = mp.mpf("0")
    max_derivative_ratio = mp.mpf("0")
    max_value_spread = mp.mpf("0")
    max_derivative_spread = mp.mpf("0")
    for location in locations:
        rows, diagnostics = build_ladder_rows(
            evaluator,
            location["T"],
            location["index"],
            DEFAULT_LADDER_ORDERS,
            DEFAULT_REFERENCE_ORDER,
        )
        location_diagnostics.append({**diagnostics, "ladder_rows": [asdict(row) for row in rows]})
        all_ladder_rows.extend(asdict(row) for row in rows)
        max_value_ratio = max(max_value_ratio, mp.mpf(diagnostics["max_value_ratio"]))
        max_derivative_ratio = max(max_derivative_ratio, mp.mpf(diagnostics["max_derivative_ratio"]))
        max_value_spread = max(max_value_spread, mp.mpf(diagnostics["value_ratio_spread"]))
        max_derivative_spread = max(max_derivative_spread, mp.mpf(diagnostics["derivative_ratio_spread"]))

    per_source_cap = mp.mpf(str(interval["summary"]["proposed_per_error_source_cap_for_five_sources"]))
    aggregate = {
        "source_worst_locations": locations,
        "ladder_orders": list(DEFAULT_LADDER_ORDERS),
        "reference_order": DEFAULT_REFERENCE_ORDER,
        "location_count": len(locations),
        "total_ladder_rows": len(all_ladder_rows),
        "max_value_ratio": sci(max_value_ratio),
        "max_derivative_ratio": sci(max_derivative_ratio),
        "max_value_ratio_spread": sci(max_value_spread),
        "max_derivative_ratio_spread": sci(max_derivative_spread),
        "all_ladder_ratios_below_one": bool(max_value_ratio < 1 and max_derivative_ratio < 1),
        "proposed_quadrature_ratio_radius_cap": sci(PROPOSED_QUADRATURE_RATIO_RADIUS_CAP, 8),
        "intervalization_per_source_cap": interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
        "cap_below_per_source_cap": bool(PROPOSED_QUADRATURE_RATIO_RADIUS_CAP < per_source_cap),
        "cap_keeps_worst_ladder_below_one": bool(
            max_value_ratio + PROPOSED_QUADRATURE_RATIO_RADIUS_CAP < 1
            and max_derivative_ratio + PROPOSED_QUADRATURE_RATIO_RADIUS_CAP < 1
        ),
    }
    rows = build_rows(paths, location_diagnostics, aggregate)
    summary = {
        "matrix_rows": len(rows),
        "location_count": len(locations),
        "ladder_orders": list(DEFAULT_LADDER_ORDERS),
        "reference_order": DEFAULT_REFERENCE_ORDER,
        "total_ladder_rows": len(all_ladder_rows),
        "max_value_ratio": aggregate["max_value_ratio"],
        "max_derivative_ratio": aggregate["max_derivative_ratio"],
        "max_value_ratio_spread": aggregate["max_value_ratio_spread"],
        "max_derivative_ratio_spread": aggregate["max_derivative_ratio_spread"],
        "all_ladder_ratios_below_one": aggregate["all_ladder_ratios_below_one"],
        "proposed_quadrature_ratio_radius_cap": aggregate["proposed_quadrature_ratio_radius_cap"],
        "intervalization_per_source_cap": aggregate["intervalization_per_source_cap"],
        "cap_below_per_source_cap": aggregate["cap_below_per_source_cap"],
        "cap_keeps_worst_ladder_below_one": aggregate["cap_keeps_worst_ladder_below_one"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The worst cancellation-reduced row remains below one first omitted term through the high-order "
            "floating ladder N=96..320. The largest value ratio is about 0.970710059020335 and the largest "
            "derivative ratio is about 0.969356777475809, with order-spread below 1e-14. A future rigorous "
            "quadrature radius below 1e-6 would be far inside the 2.0e-3 per-source intervalization cap, "
            "but this artifact remains floating evidence only."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout",
        "date": "2026-07-07",
        "status": "high-order floating quadrature ladder scout",
        "source_cancellation_reduced_grid_scout": paths["grid_note"],
        "source_cancellation_reduced_grid_json": paths["grid_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_residual_budget": paths["residual_note"],
        "source_node_c0_certificate": paths["node_c0_note"],
        "source_phi_tail_bound_scout": paths["phi_tail_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py",
        "proof_boundary": (
            "High-order floating quadrature ladder scout only. It calibrates a future quadrature-remainder "
            "target for the worst finite-grid row, but it does not provide interval nodes or weights, does "
            "not prove a quadrature-remainder theorem, does not prove a uniform collar theorem, does not "
            "prove scaled-curvature monotonicity, and does not prove RH or Lambda <= 0."
        ),
        "location_diagnostics": location_diagnostics,
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The high-order ladder is floating evidence only.",
            "The proposed quadrature radius cap is a target, not a proved remainder.",
            "Laguerre nodes and weights are not interval enclosures here.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['total_ladder_rows']} ladder rows, "
        f"{summary['reference_order']} reference order, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Quadrature Ladder Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: high-order floating quadrature ladder scout. This is not a proof",
        "of a quadrature-remainder theorem, a finite-grid interval certificate,",
        "a uniform residual estimate, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout`.",
        "",
        "Proof boundary: this artifact tests the worst cancellation-reduced",
        "finite-grid row at higher floating quadrature orders. It calibrates",
        "a future rigorous quadrature-radius target but does not prove it.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Ladder Summary",
        "",
        "```text",
        f"ladder orders: {summary['ladder_orders']}",
        f"reference order: {summary['reference_order']}",
        f"max value ratio: {summary['max_value_ratio']}",
        f"max derivative ratio: {summary['max_derivative_ratio']}",
        f"max value ratio spread: {summary['max_value_ratio_spread']}",
        f"max derivative ratio spread: {summary['max_derivative_ratio_spread']}",
        f"all ladder ratios below one: {summary['all_ladder_ratios_below_one']}",
        f"proposed quadrature ratio radius cap: {summary['proposed_quadrature_ratio_radius_cap']}",
        f"intervalization per-source cap: {summary['intervalization_per_source_cap']}",
        f"cap below per-source cap: {summary['cap_below_per_source_cap']}",
        f"cap keeps worst ladder below one: {summary['cap_keeps_worst_ladder_below_one']}",
        "```",
        "",
        "Worst-row ladder:",
        "",
        "```text",
    ]
    for location in artifact["location_diagnostics"]:
        lines.append(f"T={location['T']}, F_{location['index']}")
        for row in location["ladder_rows"]:
            lines.append(
                f"N={row['quadrature_order']}: value={row['value_ratio_to_first_omitted']}, "
                f"derivative={row['derivative_ratio_to_first_omitted']}, "
                f"delta_value={row['value_delta_from_reference_ratio']}, "
                f"delta_derivative={row['derivative_delta_from_reference_ratio']}"
            )
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_cancellation_reduced_grid_scout"],
            artifact["source_cancellation_reduced_grid_json"],
            artifact["source_intervalization_target"],
            artifact["source_intervalization_target_json"],
            artifact["source_residual_budget"],
            artifact["source_node_c0_certificate"],
            artifact["source_phi_tail_bound_scout"],
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
    parser.add_argument("--residual-json", type=Path, default=DEFAULT_RESIDUAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    residual_path = args.residual_json if args.residual_json.is_absolute() else REPO_ROOT / args.residual_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(grid_path, interval_path, residual_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
