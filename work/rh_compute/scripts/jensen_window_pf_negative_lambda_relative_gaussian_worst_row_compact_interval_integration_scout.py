#!/usr/bin/env python3
"""Build a worst-row compact-interval integration scout."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)
from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate import (  # noqa: E402
    finite_phi,
    finite_phip,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FAR_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json"
)
DEFAULT_QUADRATURE_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_PRECISION_BITS = 512
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_PANELS = ((0, 1), (1, 5), (5, 20), (20, 50), (50, 100), (100, 200))


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


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def arb_upper_text(value: flint.arb, digits: int = 40) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 40) -> str:
    return value.lower().str(digits, radius=False).replace("e", "E")


def arb_mid_text(value: flint.arb, digits: int = 40) -> str:
    return value.str(digits, radius=False).replace("e", "E")


def interval_width(value: flint.arb) -> flint.arb:
    return flint.arb(value.upper()) - flint.arb(value.lower())


def interval_hull(left: flint.arb, right: flint.arb) -> flint.arb:
    return left.union(right)


def source_paths(far_tail_path: Path, quadrature_route_path: Path) -> dict[str, str]:
    return {
        "far_tail_json": far_tail_path.relative_to(REPO_ROOT).as_posix(),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "quadrature_route_json": quadrature_route_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
        ),
        "cancellation_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md"
        ),
        "intervalization_target": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def panel_mass(alpha: flint.arb, left: int, right: int) -> flint.arb:
    gamma_alpha = (alpha + 1).gamma()
    return (flint.arb(left).gamma_upper(alpha + 1) - flint.arb(right).gamma_upper(alpha + 1)) / gamma_alpha


def panel_interval_core_widths(
    ratios: list[flint.arb],
    c0_finite: flint.arb,
    alpha: flint.arb,
    left: int,
    right: int,
) -> dict:
    t = flint.arb(DEFAULT_T)
    v_interval = interval_hull(flint.arb(left) / t, flint.arb(right) / t)
    x_left = flint.arb(0) if left == 0 else (flint.arb(left) / t).sqrt()
    x_right = (flint.arb(right) / t).sqrt()
    x_interval = interval_hull(x_left, x_right)
    phi_over_c0 = finite_phi(x_interval, DEFAULT_PHI_TERM_COUNT) / c0_finite
    derivative_phi_core = x_interval * finite_phip(x_interval, DEFAULT_PHI_TERM_COUNT) / (
        flint.arb(2) * c0_finite
    )
    polynomial = flint.arb(0)
    derivative_polynomial = flint.arb(0)
    for j in range(DEFAULT_POLYNOMIAL_M + 1):
        polynomial += ratios[j] * (v_interval**j)
        if j >= 1:
            derivative_polynomial += flint.arb(j) * ratios[j] * (v_interval**j)
    value_core = phi_over_c0 - polynomial
    derivative_core = derivative_phi_core - derivative_polynomial
    mass = panel_mass(alpha, left, right)
    mass_upper = flint.arb(mass.upper())
    value_core_width = interval_width(value_core)
    derivative_core_width = interval_width(derivative_core)
    value_width_contribution = value_core_width * mass_upper
    derivative_width_contribution = derivative_core_width * mass_upper
    return {
        "panel": f"{left}<=y<={right}",
        "left": left,
        "right": right,
        "gamma_mass_upper": arb_upper_text(mass),
        "x_interval_lower": arb_lower_text(x_interval, 18),
        "x_interval_upper": arb_upper_text(x_interval, 18),
        "value_core_ball": arb_mid_text(value_core, 18),
        "derivative_core_ball": arb_mid_text(derivative_core, 18),
        "value_core_width_upper": arb_upper_text(value_core_width),
        "derivative_core_width_upper": arb_upper_text(derivative_core_width),
        "value_width_contribution_upper": arb_upper_text(value_width_contribution),
        "derivative_width_contribution_upper": arb_upper_text(derivative_width_contribution),
        "endpoint_inflation_touches_negative_x": bool(flint.arb(x_interval.lower()) < flint.arb(0)),
        "proof_boundary": (
            "Raw interval-width diagnostic for this panel only; it is not an interval integration proof."
        ),
    }


def build_diagnostics(far_tail: dict, quadrature_route: dict) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    _ratio_rows, ratios = build_ratio_rows(2 * (DEFAULT_POLYNOMIAL_M + 1), DEFAULT_RATIO_CUTOFF_N)
    c0_finite = finite_phi(flint.arb(0), DEFAULT_PHI_TERM_COUNT)
    if c0_finite.contains(0) or not bool(c0_finite > 0):
        raise RuntimeError("finite Phi(0) interval did not certify positive")
    alpha = flint.arb(DEFAULT_ALPHA)
    panel_rows = [
        panel_interval_core_widths(ratios, c0_finite, alpha, left, right)
        for left, right in DEFAULT_PANELS
    ]
    value_raw_width = flint.arb(0)
    derivative_raw_width = flint.arb(0)
    for row in panel_rows:
        value_raw_width += flint.arb(row["value_width_contribution_upper"])
        derivative_raw_width += flint.arb(row["derivative_width_contribution_upper"])
    value_cap = flint.arb(quadrature_route["summary"]["value_unscaled_expectation_error_cap"])
    derivative_cap = flint.arb(quadrature_route["summary"]["derivative_unscaled_expectation_error_cap"])
    value_width_to_cap = value_raw_width / value_cap
    derivative_width_to_cap = derivative_raw_width / derivative_cap
    far_tail_summary = far_tail["summary"]
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "u": "1/10000",
        "alpha": DEFAULT_ALPHA,
        "precision_bits": DEFAULT_PRECISION_BITS,
        "polynomial_M": DEFAULT_POLYNOMIAL_M,
        "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "panel_count": len(panel_rows),
        "compact_interval": "0<=y<=200",
        "far_tail_split_y": far_tail_summary["split_y"],
        "far_tail_split_imported": far_tail_summary["remaining_compact_interval"] == "0<=y<=200",
        "far_tail_value_ratio_to_first_omitted_upper": far_tail_summary[
            "value_ratio_to_first_omitted_upper"
        ],
        "far_tail_derivative_ratio_to_first_omitted_upper": far_tail_summary[
            "derivative_ratio_to_first_omitted_upper"
        ],
        "quadrature_ratio_radius_cap": quadrature_route["summary"]["quadrature_ratio_radius_cap"],
        "value_unscaled_expectation_error_cap": quadrature_route["summary"][
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": quadrature_route["summary"][
            "derivative_unscaled_expectation_error_cap"
        ],
        "value_raw_interval_width_bound": arb_upper_text(value_raw_width),
        "derivative_raw_interval_width_bound": arb_upper_text(derivative_raw_width),
        "value_raw_width_to_cap_ratio_upper": arb_upper_text(value_width_to_cap),
        "derivative_raw_width_to_cap_ratio_upper": arb_upper_text(derivative_width_to_cap),
        "plain_interval_riemann_rejected": bool(
            value_width_to_cap > flint.arb(1) and derivative_width_to_cap > flint.arb(1)
        ),
        "endpoint_inflation_panel_count": sum(
            1 for row in panel_rows if row["endpoint_inflation_touches_negative_x"]
        ),
        "panel_rows": panel_rows,
        "recommended_upgrade": (
            "Use local Taylor or Chebyshev panel models for the cancellation-reduced core, with Arb "
            "coefficient balls and exact Gamma-weighted panel moments."
        ),
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrcis_01_compact_interval_import",
            role="scope_reduction",
            readiness="diagnostic_only",
            claim=(
                "The far-tail split reduces the independent interval-integration route for the worst row "
                "to the compact Gamma interval 0<=y<=200."
            ),
            diagnostics={
                "compact_interval": diagnostics["compact_interval"],
                "far_tail_split_y": diagnostics["far_tail_split_y"],
                "far_tail_split_imported": diagnostics["far_tail_split_imported"],
            },
            source_artifacts=[paths["far_tail_json"], paths["far_tail_note"]],
            proof_boundary="Scope-reduction diagnostic only; not a compact interval-integration proof.",
        ),
        MatrixRow(
            id="nlrgwrcis_02_raw_panel_width_scout",
            role="arb_panel_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "Six raw Arb panel hulls evaluate the finite n<=30 cancellation-reduced value and "
                "derivative cores on 0<=y<=200 and accumulate Gamma-mass-weighted interval widths."
            ),
            diagnostics={
                "panel_count": diagnostics["panel_count"],
                "panels": [row["panel"] for row in diagnostics["panel_rows"]],
                "value_raw_interval_width_bound": diagnostics["value_raw_interval_width_bound"],
                "derivative_raw_interval_width_bound": diagnostics[
                    "derivative_raw_interval_width_bound"
                ],
                "endpoint_inflation_panel_count": diagnostics["endpoint_inflation_panel_count"],
            },
            source_artifacts=[paths["finite_part_note"], paths["cancellation_grid_note"]],
            proof_boundary=(
                "Raw interval-width scout only; it measures dependency/cancellation loss and does not "
                "certify an integral value."
            ),
        ),
        MatrixRow(
            id="nlrgwrcis_03_width_vs_quadrature_caps",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="Plain interval-Riemann panel hulls close the compact-interval source below the quadrature caps.",
            diagnostics={
                "value_unscaled_expectation_error_cap": diagnostics[
                    "value_unscaled_expectation_error_cap"
                ],
                "derivative_unscaled_expectation_error_cap": diagnostics[
                    "derivative_unscaled_expectation_error_cap"
                ],
                "value_raw_width_to_cap_ratio_upper": diagnostics[
                    "value_raw_width_to_cap_ratio_upper"
                ],
                "derivative_raw_width_to_cap_ratio_upper": diagnostics[
                    "derivative_raw_width_to_cap_ratio_upper"
                ],
                "plain_interval_riemann_rejected": diagnostics["plain_interval_riemann_rejected"],
            },
            gap=(
                "The raw interval-width envelope is many orders of magnitude larger than the required "
                "unscaled caps, so ordinary panel hulls do not control the cancellation-reduced integral."
            ),
            source_artifacts=[paths["quadrature_route_json"], paths["quadrature_route_note"]],
            proof_boundary="Rejected diagnostic only; not a quadrature-remainder theorem or interval certificate.",
        ),
        MatrixRow(
            id="nlrgwrcis_04_panel_samples",
            role="arb_panel_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "The widest compact-panel contribution comes from the central Gamma-mass panels, not the "
                "already retired far tail."
            ),
            diagnostics={"panel_rows": diagnostics["panel_rows"]},
            source_artifacts=[paths["cancellation_grid_note"], paths["far_tail_note"]],
            proof_boundary="Panel-local diagnostics only; no finite grid, collar, or RH claim follows.",
        ),
        MatrixRow(
            id="nlrgwrcis_05_taylor_chebyshev_upgrade",
            role="live_route",
            readiness="not_ready_to_apply",
            claim=(
                "A viable compact-interval route should replace raw hulls by local Taylor or Chebyshev "
                "models, then integrate coefficient balls against exact Gamma-weighted panel moments."
            ),
            diagnostics={"recommended_upgrade": diagnostics["recommended_upgrade"]},
            source_artifacts=[paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Live route specification only; no local model certificate is produced here.",
        ),
        MatrixRow(
            id="nlrgwrcis_06_endpoint_hygiene",
            role="interval_hygiene",
            readiness="diagnostic_only",
            claim=(
                "The first panel touches y=0; the scout avoids taking square roots of midpoint-radius "
                "balls that slightly cross below zero, and records endpoint inflation explicitly."
            ),
            diagnostics={
                "endpoint_inflation_panel_count": diagnostics["endpoint_inflation_panel_count"],
                "inflation_explanation": (
                    "Hull endpoints in python-flint are outward rounded; for y=0 this can introduce a "
                    "tiny negative x lower endpoint. The scout evaluates Phi on the x hull rather than "
                    "taking sqrt of a sign-ambiguous v ball."
                ),
            },
            source_artifacts=[paths["finite_part_note"]],
            proof_boundary="Numerical hygiene diagnostic only; not an analytic endpoint theorem.",
        ),
        MatrixRow(
            id="nlrgwrcis_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted compact-interval certificate must beat the unscaled caps, compose with the "
                "far-tail split, and remain separate from all-row, rounding, and grid-to-collar coverage."
            ),
            source_artifacts=[paths["dependency_graph"], paths["intervalization_target"]],
            proof_boundary="Proof-hygiene gate only; not a finite-grid certificate, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(far_tail_path: Path, quadrature_route_path: Path) -> dict:
    far_tail = load_json(far_tail_path)
    quadrature_route = load_json(quadrature_route_path)
    paths = source_paths(far_tail_path, quadrature_route_path)
    diagnostics = build_diagnostics(far_tail, quadrature_route)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "u": diagnostics["u"],
        "alpha": diagnostics["alpha"],
        "precision_bits": diagnostics["precision_bits"],
        "polynomial_M": diagnostics["polynomial_M"],
        "ratio_cutoff_n": diagnostics["ratio_cutoff_n"],
        "phi_term_count": diagnostics["phi_term_count"],
        "panel_count": diagnostics["panel_count"],
        "compact_interval": diagnostics["compact_interval"],
        "far_tail_split_y": diagnostics["far_tail_split_y"],
        "far_tail_split_imported": diagnostics["far_tail_split_imported"],
        "quadrature_ratio_radius_cap": diagnostics["quadrature_ratio_radius_cap"],
        "value_unscaled_expectation_error_cap": diagnostics["value_unscaled_expectation_error_cap"],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "value_raw_interval_width_bound": diagnostics["value_raw_interval_width_bound"],
        "derivative_raw_interval_width_bound": diagnostics["derivative_raw_interval_width_bound"],
        "value_raw_width_to_cap_ratio_upper": diagnostics["value_raw_width_to_cap_ratio_upper"],
        "derivative_raw_width_to_cap_ratio_upper": diagnostics[
            "derivative_raw_width_to_cap_ratio_upper"
        ],
        "plain_interval_riemann_rejected": diagnostics["plain_interval_riemann_rejected"],
        "endpoint_inflation_panel_count": diagnostics["endpoint_inflation_panel_count"],
        "ready_to_apply_rows": 0,
        "target_closing": diagnostics["target_closing"],
        "recommended_upgrade": diagnostics["recommended_upgrade"],
        "main_finding": (
            "A six-panel raw Arb interval scout on the remaining compact interval 0<=y<=200 is far too "
            "wide: the value and derivative width-to-cap ratios are respectively "
            f"{diagnostics['value_raw_width_to_cap_ratio_upper']} and "
            f"{diagnostics['derivative_raw_width_to_cap_ratio_upper']}. Plain interval-Riemann hulls "
            "are therefore rejected; the live route is a local Taylor/Chebyshev panel model with exact "
            "Gamma-weighted moments. This does not close the compact interval."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout",
        "date": "2026-07-07",
        "status": "worst-row compact-interval integration scout",
        "source_far_tail_split_certificate": paths["far_tail_note"],
        "source_far_tail_split_json": paths["far_tail_json"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_quadrature_remainder_route_json": paths["quadrature_route_json"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_cancellation_reduced_grid_scout": paths["cancellation_grid_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.py"
        ),
        "proof_boundary": (
            "Worst-row compact-interval integration scout only. It imports the far-tail split and "
            "diagnoses raw Arb panel widths on 0<=y<=200, but it does not prove a compact "
            "interval-integration certificate, does not prove a quadrature-remainder theorem, does not "
            "aggregate rounding, does not cover all rows or quadrature orders, does not prove a "
            "finite-grid interval certificate, does not prove a uniform collar theorem, and does not "
            "prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "Plain interval-Riemann panel hulls are explicitly rejected here.",
            "The compact interval 0<=y<=200 remains open.",
            "The far-tail split is imported as support only, not as a quadrature-remainder theorem.",
            "All-row coverage, rounding aggregation, and grid-to-collar coverage remain separate.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval "
        f"integration scout: {summary['matrix_rows']} rows, 0 issues, {summary['panel_count']} panels, "
        "plain interval Riemann rejected, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Compact-Interval Integration Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row compact-interval integration scout. This is not a proof",
        "of a compact interval-integration certificate, quadrature-remainder",
        "theorem, finite-grid interval certificate, uniform collar theorem, RH,",
        "or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout`.",
        "",
        "Proof boundary: this artifact diagnoses only raw Arb panel hulls",
        "on the remaining compact interval `0<=y<=200` for `T=10000`,",
        "`F_21`. Plain interval-Riemann hulls are explicitly rejected.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Compact Scout",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"compact interval: {summary['compact_interval']}",
        f"far-tail split y: {summary['far_tail_split_y']}",
        f"panel count: {summary['panel_count']}",
        f"value unscaled cap: {summary['value_unscaled_expectation_error_cap']}",
        f"derivative unscaled cap: {summary['derivative_unscaled_expectation_error_cap']}",
        f"value raw interval width bound: {summary['value_raw_interval_width_bound']}",
        f"derivative raw interval width bound: {summary['derivative_raw_interval_width_bound']}",
        f"value raw width to cap ratio upper: {summary['value_raw_width_to_cap_ratio_upper']}",
        f"derivative raw width to cap ratio upper: {summary['derivative_raw_width_to_cap_ratio_upper']}",
        f"plain interval Riemann rejected: {summary['plain_interval_riemann_rejected']}",
        f"target closing: {summary['target_closing']}",
        "```",
        "",
        "Required upgrade:",
        "",
        "```text",
        summary["recommended_upgrade"],
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_far_tail_split_certificate"],
        artifact["source_far_tail_split_json"],
        artifact["source_quadrature_remainder_route_matrix"],
        artifact["source_quadrature_remainder_route_json"],
        artifact["source_finite_part_weighted_sum_interval_certificate"],
        artifact["source_cancellation_reduced_grid_scout"],
        artifact["source_intervalization_target"],
        artifact["source_uniform_remainder_target"],
        artifact["source_dependency_graph"],
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--far-tail-json", type=Path, default=DEFAULT_FAR_TAIL_JSON)
    parser.add_argument("--quadrature-route-json", type=Path, default=DEFAULT_QUADRATURE_ROUTE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    far_tail_path = args.far_tail_json if args.far_tail_json.is_absolute() else REPO_ROOT / args.far_tail_json
    quadrature_route_path = (
        args.quadrature_route_json
        if args.quadrature_route_json.is_absolute()
        else REPO_ROOT / args.quadrature_route_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(far_tail_path, quadrature_route_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row compact-interval integration scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
