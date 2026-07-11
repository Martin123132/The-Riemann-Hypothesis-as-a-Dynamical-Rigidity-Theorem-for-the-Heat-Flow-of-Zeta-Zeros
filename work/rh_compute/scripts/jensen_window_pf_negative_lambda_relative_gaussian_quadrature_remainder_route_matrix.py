#!/usr/bin/env python3
"""Build a quadrature-remainder route matrix for the relative-Gaussian worst row."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FIRST_OMITTED_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_QUADRATURE_LADDER_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json"
)
DEFAULT_FINITE_PLUS_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA_NUMERATOR = 41
DEFAULT_ALPHA_DENOMINATOR = 2
DEFAULT_QUADRATURE_ORDER = 320
DEFAULT_DERIVATIVE_ORDER = 2 * DEFAULT_QUADRATURE_ORDER
DEFAULT_MPMATH_DPS = 120

getcontext().prec = 120


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


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def fixed(value: Decimal) -> str:
    return format(value, "f")


def sci_decimal(value: Decimal, digits: int = 18) -> str:
    return f"{value:.{digits}E}"


def sci_mpf(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(value, n=digits, min_fixed=-6, max_fixed=6).replace("e", "E")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(
    first_omitted_path: Path,
    interval_path: Path,
    quadrature_ladder_path: Path,
    finite_plus_tail_path: Path,
) -> dict[str, str]:
    return {
        "first_omitted_json": first_omitted_path.relative_to(REPO_ROOT).as_posix(),
        "first_omitted_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
        ),
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "quadrature_ladder_json": quadrature_ladder_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_ladder_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md"
        ),
        "finite_plus_tail_json": finite_plus_tail_path.relative_to(REPO_ROOT).as_posix(),
        "finite_plus_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md"
        ),
        "cancellation_reduced_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md"
        ),
        "phi_tail_grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def laguerre_remainder_prefactor(quadrature_order: int, alpha: mp.mpf) -> mp.mpf:
    """Return N! Gamma(N+alpha+1)/(Gamma(alpha+1) (2N)!)."""

    return (
        mp.gamma(quadrature_order + 1)
        * mp.gamma(quadrature_order + alpha + 1)
        / (mp.gamma(alpha + 1) * mp.gamma(2 * quadrature_order + 1))
    )


def build_diagnostics(first_omitted: dict, interval: dict, ladder: dict, finite_plus_tail: dict) -> dict:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    alpha = mp.mpf(DEFAULT_ALPHA_NUMERATOR) / mp.mpf(DEFAULT_ALPHA_DENOMINATOR)
    prefactor_mpf = laguerre_remainder_prefactor(DEFAULT_QUADRATURE_ORDER, alpha)
    prefactor = dec(sci_mpf(prefactor_mpf, 80))
    u = Decimal(1) / Decimal(DEFAULT_T)
    ratio_cap = dec(ladder["summary"]["proposed_quadrature_ratio_radius_cap"])
    per_source_cap = dec(interval["summary"]["proposed_per_error_source_cap_for_five_sources"])

    value_scaled_cap = dec(first_omitted["summary"]["finest_value_scaled_absolute_radius_cap"])
    derivative_scaled_cap = dec(first_omitted["summary"]["finest_derivative_scaled_absolute_radius_cap"])
    value_unscaled_cap = value_scaled_cap * (u**3)
    derivative_unscaled_cap = derivative_scaled_cap * (u**2)
    value_derivative_sup_cap = value_unscaled_cap / prefactor
    derivative_derivative_sup_cap = derivative_unscaled_cap / prefactor

    value_finite_plus_tail = dec(finite_plus_tail["summary"]["value_composed_ratio_upper_using_full_tail_cap"])
    derivative_finite_plus_tail = dec(
        finite_plus_tail["summary"]["derivative_composed_ratio_upper_using_full_tail_cap"]
    )
    value_with_quadrature_cap = value_finite_plus_tail + ratio_cap
    derivative_with_quadrature_cap = derivative_finite_plus_tail + ratio_cap

    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "u": "1/10000",
        "alpha": f"{DEFAULT_ALPHA_NUMERATOR}/{DEFAULT_ALPHA_DENOMINATOR}",
        "quadrature_order": DEFAULT_QUADRATURE_ORDER,
        "derivative_order": DEFAULT_DERIVATIVE_ORDER,
        "laguerre_error_prefactor_formula": "N!*Gamma(N+alpha+1)/(Gamma(alpha+1)*(2N)!)",
        "laguerre_error_prefactor": sci_decimal(prefactor, 60),
        "laguerre_error_prefactor_log10": sci_mpf(mp.log10(prefactor_mpf), 40),
        "monic_norm_formula": "N!*Gamma(N+alpha+1)",
        "normalized_expectation_error_bound": "prefactor * sup_y |f^(2N)(y)|",
        "quadrature_ratio_radius_cap": fixed(ratio_cap),
        "intervalization_per_source_cap": interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
        "quadrature_cap_below_per_source_cap": bool(ratio_cap < per_source_cap),
        "value_scaled_absolute_radius_cap": sci_decimal(value_scaled_cap, 18),
        "derivative_scaled_absolute_radius_cap": sci_decimal(derivative_scaled_cap, 18),
        "value_unscaled_expectation_error_cap": sci_decimal(value_unscaled_cap, 18),
        "derivative_unscaled_expectation_error_cap": sci_decimal(derivative_unscaled_cap, 18),
        "value_required_640th_derivative_supremum_cap": sci_decimal(value_derivative_sup_cap, 40),
        "derivative_required_640th_derivative_supremum_cap": sci_decimal(
            derivative_derivative_sup_cap, 40
        ),
        "finite_plus_tail_value_ratio_upper": fixed(value_finite_plus_tail),
        "finite_plus_tail_derivative_ratio_upper": fixed(derivative_finite_plus_tail),
        "value_ratio_upper_after_quadrature_cap": fixed(value_with_quadrature_cap),
        "derivative_ratio_upper_after_quadrature_cap": fixed(derivative_with_quadrature_cap),
        "both_ratios_still_below_one_after_quadrature_cap": bool(
            value_with_quadrature_cap < Decimal(1) and derivative_with_quadrature_cap < Decimal(1)
        ),
        "ladder_order_spread_is_proof": False,
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgqrrm_01_gauss_laguerre_remainder_formula",
            role="exact_formula_target",
            readiness="available_exact_formula",
            claim=(
                "For normalized generalized Gauss-Laguerre expectation with alpha=41/2 and N=320, "
                "the classical Gaussian-quadrature remainder prefactor is "
                "N!*Gamma(N+alpha+1)/(Gamma(alpha+1)*(2N)!)."
            ),
            diagnostics={
                "alpha": diagnostics["alpha"],
                "quadrature_order": diagnostics["quadrature_order"],
                "derivative_order": diagnostics["derivative_order"],
                "laguerre_error_prefactor_formula": diagnostics["laguerre_error_prefactor_formula"],
                "laguerre_error_prefactor": diagnostics["laguerre_error_prefactor"],
                "monic_norm_formula": diagnostics["monic_norm_formula"],
                "normalized_expectation_error_bound": diagnostics["normalized_expectation_error_bound"],
            },
            source_artifacts=[paths["quadrature_ladder_note"], paths["cancellation_reduced_grid_note"]],
            proof_boundary=(
                "Exact quadrature-remainder formula target only; no 640th-derivative bound is proved here."
            ),
        ),
        MatrixRow(
            id="nlrgqrrm_02_ratio_cap_to_derivative_sup_caps",
            role="sufficient_condition",
            readiness="open_requirement",
            claim=(
                "A 1e-6 quadrature ratio-radius certificate at the worst row would follow from explicit "
                "640th-derivative supremum bounds for the value and derivative cancellation-reduced cores."
            ),
            diagnostics={
                "quadrature_ratio_radius_cap": diagnostics["quadrature_ratio_radius_cap"],
                "value_scaled_absolute_radius_cap": diagnostics["value_scaled_absolute_radius_cap"],
                "derivative_scaled_absolute_radius_cap": diagnostics["derivative_scaled_absolute_radius_cap"],
                "value_unscaled_expectation_error_cap": diagnostics["value_unscaled_expectation_error_cap"],
                "derivative_unscaled_expectation_error_cap": diagnostics[
                    "derivative_unscaled_expectation_error_cap"
                ],
                "value_required_640th_derivative_supremum_cap": diagnostics[
                    "value_required_640th_derivative_supremum_cap"
                ],
                "derivative_required_640th_derivative_supremum_cap": diagnostics[
                    "derivative_required_640th_derivative_supremum_cap"
                ],
            },
            source_artifacts=[paths["first_omitted_json"], paths["first_omitted_note"]],
            proof_boundary="Sufficient derivative-sup condition only; not a derivative estimate.",
        ),
        MatrixRow(
            id="nlrgqrrm_03_finite_plus_tail_plus_quadrature_budget",
            role="budget_composition",
            readiness="available_budget_arithmetic",
            claim=(
                "Adding the proposed 1e-6 quadrature ratio-radius cap to the certified worst-row "
                "finite-plus-tail ratio uppers still keeps both channels below one first omitted term."
            ),
            diagnostics={
                "finite_plus_tail_value_ratio_upper": diagnostics["finite_plus_tail_value_ratio_upper"],
                "finite_plus_tail_derivative_ratio_upper": diagnostics[
                    "finite_plus_tail_derivative_ratio_upper"
                ],
                "quadrature_ratio_radius_cap": diagnostics["quadrature_ratio_radius_cap"],
                "value_ratio_upper_after_quadrature_cap": diagnostics[
                    "value_ratio_upper_after_quadrature_cap"
                ],
                "derivative_ratio_upper_after_quadrature_cap": diagnostics[
                    "derivative_ratio_upper_after_quadrature_cap"
                ],
                "both_ratios_still_below_one_after_quadrature_cap": diagnostics[
                    "both_ratios_still_below_one_after_quadrature_cap"
                ],
            },
            source_artifacts=[paths["finite_plus_tail_json"], paths["finite_plus_tail_note"], paths["interval_note"]],
            proof_boundary=(
                "Budget arithmetic only; the quadrature radius remains an unproved target until a remainder "
                "or interval-integration theorem supplies it."
            ),
        ),
        MatrixRow(
            id="nlrgqrrm_04_cauchy_derivative_route",
            role="live_route",
            readiness="not_ready_to_apply",
            claim=(
                "One rigorous route is to prove complex-analytic Cauchy bounds for the 640th y-derivatives "
                "of the cancellation-reduced value and derivative cores, using the even-v expansion and "
                "separate n-tail control."
            ),
            diagnostics={
                "required_value_derivative_sup_cap": diagnostics[
                    "value_required_640th_derivative_supremum_cap"
                ],
                "required_derivative_derivative_sup_cap": diagnostics[
                    "derivative_required_640th_derivative_supremum_cap"
                ],
                "polynomial_part_note": "The degree-20 polynomial part has zero 640th y-derivative.",
                "tail_warning": (
                    "Node-only x<=1 tail control is not enough for a quadrature theorem; the continuum "
                    "integral tail or analytic derivative tail must be bounded separately."
                ),
            },
            source_artifacts=[paths["phi_tail_grid_note"], paths["uniform_remainder_target"]],
            proof_boundary="Live analytic route only; no Cauchy-radius or tail derivative theorem is proved.",
        ),
        MatrixRow(
            id="nlrgqrrm_05_interval_adaptive_integration_route",
            role="live_route",
            readiness="not_ready_to_apply",
            claim=(
                "A second rigorous route is independent Arb interval integration of the cancellation-reduced "
                "Gamma expectation, split into a finite interval and a separately bounded far Gamma tail."
            ),
            diagnostics={
                "required_value_unscaled_error": diagnostics["value_unscaled_expectation_error_cap"],
                "required_derivative_unscaled_error": diagnostics[
                    "derivative_unscaled_expectation_error_cap"
                ],
                "split_warning": (
                    "The finite-grid node range x<=1 does not cover the continuum y>T tail; a split proof "
                    "must include the far-tail integral."
                ),
            },
            source_artifacts=[paths["cancellation_reduced_grid_note"], paths["uniform_remainder_target"]],
            proof_boundary="Live interval-integration route only; no integration certificate is produced here.",
        ),
        MatrixRow(
            id="nlrgqrrm_06_order_spread_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim=(
                "The observed N=96..320 order-spread below 1e-14 proves the quadrature remainder is below "
                "the required 1e-6 ratio-radius cap."
            ),
            diagnostics={
                "ladder_order_spread_is_proof": diagnostics["ladder_order_spread_is_proof"],
                "quadrature_ladder_role": "floating stability scout",
            },
            gap=(
                "Order agreement is calibration evidence only. A referee-usable proof needs a Gaussian "
                "quadrature remainder theorem, a derivative-sup certificate, or independent interval integration."
            ),
            source_artifacts=[paths["quadrature_ladder_json"], paths["quadrature_ladder_note"]],
            proof_boundary="Rejected promotion only; not a quadrature-remainder theorem.",
        ),
        MatrixRow(
            id="nlrgqrrm_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted intervalized worst-row proof may use the quadrature source only after one of the "
                "live routes supplies a machine-checkable radius below the recorded cap."
            ),
            source_artifacts=[paths["interval_json"], paths["interval_note"], paths["dependency_graph"]],
            proof_boundary=(
                "Proof-hygiene gate only; not a finite-grid interval certificate, not RH, and not Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(
    first_omitted_path: Path,
    interval_path: Path,
    quadrature_ladder_path: Path,
    finite_plus_tail_path: Path,
) -> dict:
    first_omitted = load_json(first_omitted_path)
    interval = load_json(interval_path)
    ladder = load_json(quadrature_ladder_path)
    finite_plus_tail = load_json(finite_plus_tail_path)
    paths = source_paths(first_omitted_path, interval_path, quadrature_ladder_path, finite_plus_tail_path)
    diagnostics = build_diagnostics(first_omitted, interval, ladder, finite_plus_tail)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "quadrature_order": diagnostics["quadrature_order"],
        "derivative_order": diagnostics["derivative_order"],
        "laguerre_error_prefactor": diagnostics["laguerre_error_prefactor"],
        "quadrature_ratio_radius_cap": diagnostics["quadrature_ratio_radius_cap"],
        "intervalization_per_source_cap": diagnostics["intervalization_per_source_cap"],
        "quadrature_cap_below_per_source_cap": diagnostics["quadrature_cap_below_per_source_cap"],
        "value_unscaled_expectation_error_cap": diagnostics["value_unscaled_expectation_error_cap"],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "value_required_640th_derivative_supremum_cap": diagnostics[
            "value_required_640th_derivative_supremum_cap"
        ],
        "derivative_required_640th_derivative_supremum_cap": diagnostics[
            "derivative_required_640th_derivative_supremum_cap"
        ],
        "finite_plus_tail_value_ratio_upper": diagnostics["finite_plus_tail_value_ratio_upper"],
        "finite_plus_tail_derivative_ratio_upper": diagnostics[
            "finite_plus_tail_derivative_ratio_upper"
        ],
        "value_ratio_upper_after_quadrature_cap": diagnostics["value_ratio_upper_after_quadrature_cap"],
        "derivative_ratio_upper_after_quadrature_cap": diagnostics[
            "derivative_ratio_upper_after_quadrature_cap"
        ],
        "both_ratios_still_below_one_after_quadrature_cap": diagnostics[
            "both_ratios_still_below_one_after_quadrature_cap"
        ],
        "live_route_rows": 2,
        "rejected_route_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The classical N=320 generalized Gauss-Laguerre remainder formula reduces the worst-row "
            "quadrature source to 640th-derivative control with normalized prefactor about 2.79e-159. "
            "A 1e-6 ratio-radius certificate would follow from value and derivative 640th-derivative "
            "supremum bounds below the recorded caps, and adding that cap to the finite-plus-tail budget "
            "still leaves both channels below one. This is a route matrix only: no derivative or interval "
            "integration bound is proved."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix",
        "date": "2026-07-07",
        "status": "quadrature-remainder route matrix",
        "source_first_omitted_denominator_certificate": paths["first_omitted_note"],
        "source_first_omitted_denominator_json": paths["first_omitted_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_quadrature_ladder_scout": paths["quadrature_ladder_note"],
        "source_quadrature_ladder_json": paths["quadrature_ladder_json"],
        "source_finite_plus_tail_budget_certificate": paths["finite_plus_tail_note"],
        "source_finite_plus_tail_budget_json": paths["finite_plus_tail_json"],
        "source_cancellation_reduced_grid_scout": paths["cancellation_reduced_grid_note"],
        "source_phi_tail_grid_certificate": paths["phi_tail_grid_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py"
        ),
        "proof_boundary": (
            "Quadrature-remainder route matrix only. It records the exact classical Gauss-Laguerre remainder "
            "prefactor, converts the desired 1e-6 ratio-radius target into derivative-supremum and interval "
            "integration obligations, and checks budget compatibility, but it does not prove any derivative "
            "bound, does not produce interval adaptive integration, does not close the finite-grid interval "
            "certificate, does not prove a uniform collar theorem, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The route matrix does not promote floating order-spread to a quadrature theorem.",
            "The derivative-supremum caps are sufficient conditions, not proved estimates.",
            "Node-only x<=1 tail control is not a continuum quadrature-tail theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"derivative order {summary['derivative_order']}, "
        "2 derivative-sup caps, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Quadrature-Remainder Route Matrix",
        "",
        "Date: 2026-07-07",
        "",
        "Status: quadrature-remainder route matrix. This is not a proof",
        "of a quadrature-remainder theorem, a finite-grid interval certificate,",
        "a uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix`.",
        "",
        "Proof boundary: this artifact records the exact classical",
        "Gauss-Laguerre remainder prefactor and the derivative/interval",
        "obligations needed to turn the floating quadrature ladder into a",
        "certificate. It does not prove those obligations.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Remainder Constants",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"quadrature order: {summary['quadrature_order']}",
        f"derivative order: {summary['derivative_order']}",
        f"Laguerre error prefactor: {summary['laguerre_error_prefactor']}",
        f"quadrature ratio radius cap: {summary['quadrature_ratio_radius_cap']}",
        f"intervalization per-source cap: {summary['intervalization_per_source_cap']}",
        f"quadrature cap below per-source cap: {summary['quadrature_cap_below_per_source_cap']}",
        f"value unscaled expectation error cap: {summary['value_unscaled_expectation_error_cap']}",
        f"derivative unscaled expectation error cap: {summary['derivative_unscaled_expectation_error_cap']}",
        f"value required 640th derivative sup cap: {summary['value_required_640th_derivative_supremum_cap']}",
        f"derivative required 640th derivative sup cap: {summary['derivative_required_640th_derivative_supremum_cap']}",
        "```",
        "",
        "## Budget Compatibility",
        "",
        "```text",
        f"finite-plus-tail value ratio upper: {summary['finite_plus_tail_value_ratio_upper']}",
        f"finite-plus-tail derivative ratio upper: {summary['finite_plus_tail_derivative_ratio_upper']}",
        f"value ratio upper after quadrature cap: {summary['value_ratio_upper_after_quadrature_cap']}",
        f"derivative ratio upper after quadrature cap: {summary['derivative_ratio_upper_after_quadrature_cap']}",
        f"both ratios still below one after quadrature cap: {summary['both_ratios_still_below_one_after_quadrature_cap']}",
        "```",
        "",
        "Live routes:",
        "",
        "```text",
        "Cauchy derivative route: prove 640th y-derivative bounds for the cancellation-reduced cores.",
        "Interval adaptive route: integrate the cancellation-reduced Gamma expectation with a separate far-tail bound.",
        "Rejected route: do not promote floating N=96..320 order-spread into a remainder theorem.",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_first_omitted_denominator_certificate"],
        artifact["source_first_omitted_denominator_json"],
        artifact["source_intervalization_target"],
        artifact["source_intervalization_target_json"],
        artifact["source_quadrature_ladder_scout"],
        artifact["source_quadrature_ladder_json"],
        artifact["source_finite_plus_tail_budget_certificate"],
        artifact["source_finite_plus_tail_budget_json"],
        artifact["source_cancellation_reduced_grid_scout"],
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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-omitted-json", type=Path, default=DEFAULT_FIRST_OMITTED_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--quadrature-ladder-json", type=Path, default=DEFAULT_QUADRATURE_LADDER_JSON)
    parser.add_argument("--finite-plus-tail-json", type=Path, default=DEFAULT_FINITE_PLUS_TAIL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    first_omitted_path = (
        args.first_omitted_json if args.first_omitted_json.is_absolute() else REPO_ROOT / args.first_omitted_json
    )
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    quadrature_ladder_path = (
        args.quadrature_ladder_json
        if args.quadrature_ladder_json.is_absolute()
        else REPO_ROOT / args.quadrature_ladder_json
    )
    finite_plus_tail_path = (
        args.finite_plus_tail_json
        if args.finite_plus_tail_json.is_absolute()
        else REPO_ROOT / args.finite_plus_tail_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(first_omitted_path, interval_path, quadrature_ladder_path, finite_plus_tail_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
