#!/usr/bin/env python3
"""Build a padded-range Phi-tail bound scout for the relative-Gaussian grid."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path

import mpmath as mp
import scipy.special as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md"
)

DEFAULT_T_GRID = (1156, 1500, 2000, 5000, 10000)
DEFAULT_QUADRATURE_ORDERS = (64, 96, 128, 192)
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_PHI_TRUNCATION_N = 30
DEFAULT_X_RANGE_UPPER = "1"
DEFAULT_C0_LOWER_PROXY = "0.44"
DEFAULT_MPMATH_DPS = 100


@dataclass(frozen=True)
class MatrixRow:
    id: str
    role: str
    readiness: str
    claim: str
    source_artifacts: list[str]
    proof_boundary: str
    diagnostics: dict | None = None
    gap: str | None = None


def sci(value: mp.mpf | float, digits: int = 19) -> str:
    return mp.nstr(mp.mpf(value), n=digits, min_fixed=-6, max_fixed=6)


def sci_e(value: mp.mpf | float, digits: int = 18) -> str:
    if isinstance(value, float):
        return f"{value:.{digits}E}"
    text = f"{Decimal(str(value)):.{digits}E}"
    coeff, exponent = text.split("E")
    sign, magnitude = exponent[0], exponent[1:]
    return f"{coeff}E{sign}{magnitude.rjust(2, '0')}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(grid_path: Path, interval_path: Path) -> dict[str, str]:
    return {
        "grid_json": grid_path.relative_to(REPO_ROOT).as_posix(),
        "grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "actual_endpoint_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md",
        "asymptotic_target_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md"
        ),
        "uniform_remainder_target_note": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def observed_node_range(
    t_grid: tuple[int, ...],
    orders: tuple[int, ...],
    indices: tuple[int, ...],
    x_range_upper: mp.mpf,
) -> dict:
    max_row: dict | None = None
    row_count = 0
    for T in t_grid:
        for order in orders:
            for index in indices:
                nodes, _weights = sp.roots_genlaguerre(order, index - 0.5)
                node = float(nodes[-1])
                x_value = (node / float(T)) ** 0.5
                row_count += 1
                if max_row is None or x_value > max_row["x_float"]:
                    max_row = {
                        "x_float": x_value,
                        "node_float": node,
                        "T": T,
                        "quadrature_order": order,
                        "index": index,
                    }
    assert max_row is not None
    return {
        "node_range_row_count": row_count,
        "max_observed_x": sci_e(max_row["x_float"]),
        "max_observed_node": sci_e(max_row["node_float"]),
        "max_observed_location": {
            "T": max_row["T"],
            "quadrature_order": max_row["quadrature_order"],
            "index": max_row["index"],
        },
        "padded_x_range_upper": sci_e(x_range_upper),
        "observed_slack_to_padded_range": sci_e(float(x_range_upper) - max_row["x_float"]),
        "all_observed_nodes_inside_padded_range": bool(max_row["x_float"] < float(x_range_upper)),
        "proof_boundary": (
            "Floating SciPy node-range scout only. The tail majorant below is stated for the padded "
            "range 0<=x<=1, but a promoted proof still needs interval-certified node enclosures."
        ),
    }


def polynomial_gaussian_tail_bound(first_n: int, degree: int, first_term: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    rho = (mp.mpf(first_n + 1) / mp.mpf(first_n)) ** degree * mp.e ** (-mp.pi * (2 * first_n + 1))
    return first_term / (1 - rho), rho


def tail_bounds(first_n: int, x_upper: mp.mpf, c0_lower_proxy: mp.mpf) -> dict:
    pi = mp.pi
    exp4 = mp.e ** (4 * x_upper)
    exp5 = mp.e ** (5 * x_upper)
    exp9 = mp.e ** (9 * x_upper)
    n = mp.mpf(first_n)

    value_first_term = (2 * pi * pi * n**4 * exp9 + 3 * pi * n**2 * exp5) * mp.e ** (-pi * n**2)
    value_bound, value_rho = polynomial_gaussian_tail_bound(first_n, 4, value_first_term)

    derivative_first_term = (
        18 * pi * pi * n**4 * exp9
        + 15 * pi * n**2 * exp5
        + 4 * pi * n**2 * exp4 * (2 * pi * pi * n**4 * exp9 + 3 * pi * n**2 * exp5)
    ) * mp.e ** (-pi * n**2)
    derivative_bound, derivative_rho = polynomial_gaussian_tail_bound(first_n, 6, derivative_first_term)

    c0_first_term = (2 * pi * pi * n**4 + 3 * pi * n**2) * mp.e ** (-pi * n**2)
    c0_bound, c0_rho = polynomial_gaussian_tail_bound(first_n, 4, c0_first_term)

    c0_truncated = mp.fsum((2 * pi * pi * k**4 - 3 * pi * k**2) * mp.e ** (-pi * k**2) for k in range(1, first_n))

    normalized_value = value_bound / c0_lower_proxy
    normalized_derivative_core = x_upper * derivative_bound / (2 * c0_lower_proxy)
    denominator_relative = c0_bound / c0_lower_proxy
    return {
        "phi_truncation_n": first_n - 1,
        "tail_start_n": first_n,
        "x_range_upper": sci_e(x_upper),
        "value_tail_majorant_degree": 4,
        "derivative_tail_majorant_degree": 6,
        "value_tail_first_term": sci(value_first_term, 80),
        "value_tail_geometric_ratio_bound": sci(value_rho, 80),
        "value_phi_tail_bound": sci(value_bound, 80),
        "derivative_tail_first_term": sci(derivative_first_term, 80),
        "derivative_tail_geometric_ratio_bound": sci(derivative_rho, 80),
        "derivative_phi_prime_tail_bound": sci(derivative_bound, 80),
        "c0_tail_first_term": sci(c0_first_term, 80),
        "c0_tail_geometric_ratio_bound": sci(c0_rho, 80),
        "c0_tail_bound": sci(c0_bound, 80),
        "c0_truncated_through_n30": sci(c0_truncated, 80),
        "c0_lower_proxy": sci_e(c0_lower_proxy),
        "c0_truncated_margin_over_proxy": sci(c0_truncated - c0_lower_proxy, 80),
        "normalized_value_tail_bound_using_c0_proxy": sci(normalized_value, 80),
        "normalized_derivative_core_tail_bound_using_c0_proxy": sci(normalized_derivative_core, 80),
        "denominator_relative_tail_bound_using_c0_proxy": sci(denominator_relative, 80),
        "proof_boundary": (
            "Analytic positive-majorant calculation for 0<=x<=1 only. The c0 lower bound is a "
            "floating proxy here; interval promotion still needs a certified lower enclosure for Phi(0)."
        ),
    }


def build_rows(paths: dict[str, str], node_diagnostics: dict, tail_diagnostics: dict, budget: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgtb_01_padded_node_range_scout",
            role="finite_node_range_scout",
            readiness="not_ready_to_apply",
            claim=(
                "The finite cancellation-reduced grid has floating node-induced x=sqrt(node/T) values below "
                "the padded range endpoint x=1."
            ),
            source_artifacts=[paths["grid_json"], paths["grid_note"]],
            diagnostics=node_diagnostics,
            proof_boundary=(
                "Floating node-range scout only; not an interval enclosure of generalized Laguerre roots."
            ),
        ),
        MatrixRow(
            id="nlrgtb_02_value_phi_tail_majorant",
            role="analytic_tail_bound",
            readiness="not_ready_to_apply",
            claim=(
                "For 0<=x<=1, the n>30 value tail of Phi is bounded by a positive polynomial-Gaussian "
                "majorant starting at n=31."
            ),
            source_artifacts=[paths["interval_note"], paths["asymptotic_target_note"]],
            diagnostics={
                "value_phi_tail_bound": tail_diagnostics["value_phi_tail_bound"],
                "normalized_value_tail_bound_using_c0_proxy": tail_diagnostics[
                    "normalized_value_tail_bound_using_c0_proxy"
                ],
                "per_source_intervalization_cap": budget["per_source_cap"],
                "ratio_to_per_source_cap": budget["value_ratio_to_per_source_cap"],
            },
            proof_boundary=(
                "Padded-range tail majorant only; not a full interval certificate until the node range and "
                "Phi(0) lower bound are interval-certified."
            ),
        ),
        MatrixRow(
            id="nlrgtb_03_derivative_phi_prime_tail_majorant",
            role="analytic_tail_bound",
            readiness="not_ready_to_apply",
            claim=(
                "For 0<=x<=1, the n>30 Phi' tail gives a normalized derivative-core tail bound far below "
                "the per-source intervalization cap."
            ),
            source_artifacts=[paths["interval_note"], paths["asymptotic_target_note"]],
            diagnostics={
                "derivative_phi_prime_tail_bound": tail_diagnostics["derivative_phi_prime_tail_bound"],
                "normalized_derivative_core_tail_bound_using_c0_proxy": tail_diagnostics[
                    "normalized_derivative_core_tail_bound_using_c0_proxy"
                ],
                "per_source_intervalization_cap": budget["per_source_cap"],
                "ratio_to_per_source_cap": budget["derivative_ratio_to_per_source_cap"],
            },
            proof_boundary=(
                "Padded-range derivative-tail majorant only; not a quadrature-remainder theorem, "
                "not coefficient propagation, and not a uniform collar theorem."
            ),
        ),
        MatrixRow(
            id="nlrgtb_04_c0_tail_and_lower_proxy",
            role="denominator_tail_scout",
            readiness="not_ready_to_apply",
            claim=(
                "The n>30 absolute tail for Phi(0) is negligible relative to the floating lower proxy "
                "Phi_30(0)>0.44 used only for this diagnostic normalization."
            ),
            source_artifacts=[paths["interval_note"], paths["grid_note"]],
            diagnostics={
                "c0_tail_bound": tail_diagnostics["c0_tail_bound"],
                "c0_truncated_through_n30": tail_diagnostics["c0_truncated_through_n30"],
                "c0_lower_proxy": tail_diagnostics["c0_lower_proxy"],
                "denominator_relative_tail_bound_using_c0_proxy": tail_diagnostics[
                    "denominator_relative_tail_bound_using_c0_proxy"
                ],
            },
            proof_boundary=(
                "Floating denominator scout only; a promoted proof must replace the proxy by an interval "
                "lower bound for Phi(0)."
            ),
        ),
        MatrixRow(
            id="nlrgtb_05_intervalization_handoff",
            role="conditional_budget_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "The Phi/Phi'/Phi(0) n-tail component of nlrgit_03 is small enough to fit the "
                "2.0e-3 per-source ratio budget once node range and c0 lower enclosures are certified."
            ),
            source_artifacts=[paths["interval_json"], paths["interval_note"], paths["dependency_graph"]],
            diagnostics=budget,
            proof_boundary=(
                "Conditional handoff only; it does not close the Laguerre node/weight, quadrature, "
                "coefficient-propagation, rounding, or grid-to-collar obligations."
            ),
        ),
        MatrixRow(
            id="nlrgtb_06_full_obligation_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The padded Phi-tail scout by itself completes the finite-grid interval certificate.",
            source_artifacts=[paths["interval_note"], paths["uniform_remainder_target_note"]],
            gap=(
                "It only addresses the n>30 Phi-tail source conditionally. Laguerre node/weight intervals, "
                "quadrature error, coefficient-ball propagation, rounding aggregation, and the grid-to-collar "
                "bridge remain open."
            ),
            proof_boundary=(
                "Rejected promotion only; not a proof of the finite-grid certificate, the uniform collar "
                "remainder theorem, scaled-curvature monotonicity, RH, or Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(grid_path: Path, interval_path: Path) -> dict:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    grid = load_json(grid_path)
    interval = load_json(interval_path)
    paths = source_paths(grid_path, interval_path)
    first_n = DEFAULT_PHI_TRUNCATION_N + 1
    x_range_upper = mp.mpf(DEFAULT_X_RANGE_UPPER)
    c0_lower_proxy = mp.mpf(DEFAULT_C0_LOWER_PROXY)
    node_diagnostics = observed_node_range(DEFAULT_T_GRID, DEFAULT_QUADRATURE_ORDERS, DEFAULT_INDICES, x_range_upper)
    tail_diagnostics = tail_bounds(first_n, x_range_upper, c0_lower_proxy)

    per_source_cap = mp.mpf(str(interval["summary"]["proposed_per_error_source_cap_for_five_sources"]))
    value_ratio = mp.mpf(tail_diagnostics["normalized_value_tail_bound_using_c0_proxy"]) / per_source_cap
    derivative_ratio = mp.mpf(tail_diagnostics["normalized_derivative_core_tail_bound_using_c0_proxy"]) / per_source_cap
    denominator_ratio = mp.mpf(tail_diagnostics["denominator_relative_tail_bound_using_c0_proxy"]) / per_source_cap
    budget = {
        "source_obligation_id": "nlrgit_03_phi_and_c0_interval_tail",
        "per_source_cap": interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
        "total_ratio_cap": interval["summary"]["proposed_total_ratio_error_cap"],
        "value_ratio_to_per_source_cap": sci(value_ratio, 80),
        "derivative_ratio_to_per_source_cap": sci(derivative_ratio, 80),
        "denominator_ratio_to_per_source_cap": sci(denominator_ratio, 80),
        "tail_bounds_below_per_source_cap": bool(value_ratio < 1 and derivative_ratio < 1 and denominator_ratio < 1),
        "conditional_requirements": [
            "Replace the floating SciPy node range by interval enclosures proving x<=1 for all grid nodes.",
            "Replace the floating c0 lower proxy by an interval-certified lower bound Phi(0)>=0.44.",
        ],
        "remaining_non_tail_obligations": [
            "Laguerre node/weight intervals beyond the x-range check.",
            "Generalized Gauss-Laguerre quadrature-remainder bound.",
            "Coefficient-ratio and first-omitted denominator interval propagation.",
            "Rounding aggregation below the total ratio cap.",
            "Continuum-in-T or grid-to-collar bridge.",
        ],
    }
    rows = build_rows(paths, node_diagnostics, tail_diagnostics, budget)
    summary = {
        "matrix_rows": len(rows),
        "tail_bound_rows": 3,
        "conditional_requirement_rows": 2,
        "ready_to_apply_rows": 0,
        "source_grid_rows": grid["summary"]["grid_rows"],
        "source_t_grid_count": grid["summary"]["t_grid_count"],
        "source_index_count": grid["summary"]["index_count"],
        "source_selected_quadrature_order": grid["summary"]["selected_quadrature_order"],
        "node_range_row_count": node_diagnostics["node_range_row_count"],
        "max_observed_x": node_diagnostics["max_observed_x"],
        "padded_x_range_upper": node_diagnostics["padded_x_range_upper"],
        "observed_slack_to_padded_range": node_diagnostics["observed_slack_to_padded_range"],
        "all_observed_nodes_inside_padded_range": node_diagnostics["all_observed_nodes_inside_padded_range"],
        "phi_truncation_n": tail_diagnostics["phi_truncation_n"],
        "tail_start_n": tail_diagnostics["tail_start_n"],
        "value_phi_tail_bound": tail_diagnostics["value_phi_tail_bound"],
        "derivative_phi_prime_tail_bound": tail_diagnostics["derivative_phi_prime_tail_bound"],
        "c0_tail_bound": tail_diagnostics["c0_tail_bound"],
        "normalized_value_tail_bound_using_c0_proxy": tail_diagnostics[
            "normalized_value_tail_bound_using_c0_proxy"
        ],
        "normalized_derivative_core_tail_bound_using_c0_proxy": tail_diagnostics[
            "normalized_derivative_core_tail_bound_using_c0_proxy"
        ],
        "denominator_relative_tail_bound_using_c0_proxy": tail_diagnostics[
            "denominator_relative_tail_bound_using_c0_proxy"
        ],
        "per_source_intervalization_cap": budget["per_source_cap"],
        "tail_bounds_below_per_source_cap": budget["tail_bounds_below_per_source_cap"],
        "target_closing": False,
        "main_finding": (
            "On the padded range 0<=x<=1, the omitted n>30 tails in Phi, Phi', and Phi(0) are "
            "far below the 2.0e-3 per-source intervalization ratio cap after normalization by the "
            "diagnostic c0 lower proxy 0.44. This narrows the nlrgit_03 tail source to two concrete "
            "certification tasks: interval-prove the grid node range x<=1 and interval-prove Phi(0)>=0.44."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout",
        "date": "2026-07-07",
        "status": "analytic padded-range tail scout",
        "source_cancellation_reduced_grid_scout": paths["grid_note"],
        "source_cancellation_reduced_grid_json": paths["grid_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_actual_endpoint_scout": paths["actual_endpoint_note"],
        "source_asymptotic_remainder_target": paths["asymptotic_target_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target_note"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py",
        "proof_boundary": (
            "Analytic padded-range tail scout only. It bounds a Phi-tail source conditionally on x<=1 "
            "and Phi(0)>=0.44, but it does not interval-certify Laguerre nodes or weights, does not "
            "bound quadrature remainder, does not propagate coefficient balls, does not prove a uniform "
            "collar theorem, does not prove scaled-curvature monotonicity, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "tail_diagnostics": tail_diagnostics,
        "invariants": [
            "No row is ready_to_apply.",
            "The padded Phi-tail bound is conditional on future interval node and Phi(0) lower certificates.",
            "The cancellation-reduced floating grid is not promoted to a finite-grid interval certificate.",
            "The finite grid is not promoted to a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['tail_bound_rows']} tail bounds below 1e-1000, "
        f"{summary['conditional_requirement_rows']} conditional requirements, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Phi Tail Bound Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: analytic padded-range tail scout. This is not a proof",
        "of a finite-grid interval certificate, a uniform residual estimate,",
        "scaled-curvature monotonicity, cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout`.",
        "",
        "Proof boundary: this artifact bounds the omitted `n>30` Phi tails on",
        "the padded range `0<=x<=1`, conditionally on later interval proofs",
        "of the node range and `Phi(0)>=0.44`. It does not prove those",
        "interval facts or any quadrature-remainder theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Node Range",
        "",
        "```text",
        f"max observed node x: {summary['max_observed_x']}",
        f"padded x range: {summary['padded_x_range_upper']}",
        f"observed slack to padded range: {summary['observed_slack_to_padded_range']}",
        f"all observed nodes inside padded range: {summary['all_observed_nodes_inside_padded_range']}",
        "```",
        "",
        "## Tail Bounds",
        "",
        "```text",
        f"Phi truncation n: {summary['phi_truncation_n']}",
        f"tail start n: {summary['tail_start_n']}",
        f"value Phi tail bound: {summary['value_phi_tail_bound']}",
        f"derivative Phi-prime tail bound: {summary['derivative_phi_prime_tail_bound']}",
        f"c0 tail bound: {summary['c0_tail_bound']}",
        f"normalized value tail bound using c0 proxy: {summary['normalized_value_tail_bound_using_c0_proxy']}",
        (
            "normalized derivative-core tail bound using c0 proxy: "
            f"{summary['normalized_derivative_core_tail_bound_using_c0_proxy']}"
        ),
        f"denominator relative tail bound using c0 proxy: {summary['denominator_relative_tail_bound_using_c0_proxy']}",
        f"per-source intervalization cap: {summary['per_source_intervalization_cap']}",
        f"tail bounds below per-source cap: {summary['tail_bounds_below_per_source_cap']}",
        "```",
        "",
        "Conditional requirements:",
        "",
        "```text",
    ]
    for item in artifact["matrix_rows"][4]["diagnostics"]["conditional_requirements"]:
        lines.append(item)
    lines.extend(
        [
            "```",
            "",
            "Remaining non-tail obligations:",
            "",
            "```text",
        ]
    )
    for item in artifact["matrix_rows"][4]["diagnostics"]["remaining_non_tail_obligations"]:
        lines.append(item)
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
            artifact["source_actual_endpoint_scout"],
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
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(grid_path, interval_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
