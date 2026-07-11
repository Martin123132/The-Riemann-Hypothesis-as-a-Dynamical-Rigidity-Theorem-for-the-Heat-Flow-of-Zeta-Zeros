#!/usr/bin/env python3
"""Build an endpoint x-panel route matrix for the relative-Gaussian compact route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INTERPOLATION_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json"
)
DEFAULT_ENDPOINT_PARITY_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.md"

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_X_INTERVAL = "0<=x<=0.01"
DEFAULT_Y_INTERVAL = "0<=y<=1"
DEFAULT_DEGREES = (8, 12, 16, 20, 24, 32)
DEFAULT_RHOS = ("1.5", "2.0", "3.0")

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


def sci(value: Decimal, digits: int = 30) -> str:
    return f"{value:.{digits}E}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(interpolation_route_path: Path, endpoint_parity_path: Path) -> dict[str, str]:
    return {
        "interpolation_route_json": interpolation_route_path.relative_to(REPO_ROOT).as_posix(),
        "interpolation_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md"
        ),
        "endpoint_parity_json": endpoint_parity_path.relative_to(REPO_ROOT).as_posix(),
        "endpoint_parity_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md"
        ),
        "arb_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md"
        ),
        "compact_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md"
        ),
        "quadrature_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
        ),
        "intervalization_target": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "formal_tail_obstruction": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"
        ),
    }


def bernstein_constant(degree: int, rho: Decimal) -> Decimal:
    return Decimal(4) * (rho ** Decimal(-degree)) / (rho - Decimal(1))


def budget_rows(value_cap: Decimal, derivative_cap: Decimal, first_panel_mass: Decimal) -> list[dict]:
    panel_count = Decimal(6)
    rows: list[dict] = []
    for degree in DEFAULT_DEGREES:
        for rho_text in DEFAULT_RHOS:
            rho = Decimal(rho_text)
            constant = bernstein_constant(degree, rho)
            value_budget = value_cap / (panel_count * first_panel_mass * constant)
            derivative_budget = derivative_cap / (panel_count * first_panel_mass * constant)
            rows.append(
                {
                    "degree": degree,
                    "rho": rho_text,
                    "bernstein_error_constant": sci(constant),
                    "value_sup_norm_budget_on_first_x_panel": sci(value_budget),
                    "derivative_sup_norm_budget_on_first_x_panel": sci(derivative_budget),
                    "proof_boundary": (
                        "First x-panel sufficient-condition budget only; requires an independent "
                        "analytic-domain and sup-norm certificate in the x variable."
                    ),
                }
            )
    return rows


def build_diagnostics(interpolation_route: dict, endpoint_parity: dict) -> dict:
    interpolation_summary = interpolation_route["summary"]
    value_cap = dec(interpolation_summary["value_unscaled_expectation_error_cap"])
    derivative_cap = dec(interpolation_summary["derivative_unscaled_expectation_error_cap"])
    panel_rows = interpolation_route["matrix_rows"][1]["diagnostics"]["panel_masses"]
    first_panel = next(row for row in panel_rows if row["panel"] == DEFAULT_Y_INTERVAL)
    heaviest_panel = next(row for row in panel_rows if row["panel"] == interpolation_summary["heaviest_panel"])
    first_mass = dec(first_panel["gamma_mass_upper"])
    heaviest_mass = dec(heaviest_panel["gamma_mass_upper"])
    budgets = budget_rows(value_cap, derivative_cap, first_mass)
    rho_2_degree_32 = next(row for row in budgets if row["rho"] == "2.0" and row["degree"] == 32)
    rho_3_degree_32 = next(row for row in budgets if row["rho"] == "3.0" and row["degree"] == 32)
    x_weight_power = int(Decimal(2) * dec(DEFAULT_ALPHA) + Decimal(1))
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "alpha": DEFAULT_ALPHA,
        "x_interval": DEFAULT_X_INTERVAL,
        "y_interval": DEFAULT_Y_INTERVAL,
        "x_right_endpoint": "0.01",
        "compact_interval": interpolation_summary["compact_interval"],
        "panel_count": interpolation_summary["panel_count"],
        "first_panel_mass_upper": first_panel["gamma_mass_upper"],
        "heaviest_panel": heaviest_panel["panel"],
        "heaviest_panel_mass_upper": heaviest_panel["gamma_mass_upper"],
        "first_to_heaviest_mass_ratio": sci(first_mass / heaviest_mass),
        "value_unscaled_expectation_error_cap": interpolation_summary["value_unscaled_expectation_error_cap"],
        "derivative_unscaled_expectation_error_cap": interpolation_summary[
            "derivative_unscaled_expectation_error_cap"
        ],
        "first_panel_value_sup_budget_without_bernstein": sci(value_cap / (Decimal(6) * first_mass)),
        "first_panel_derivative_sup_budget_without_bernstein": sci(
            derivative_cap / (Decimal(6) * first_mass)
        ),
        "degrees": list(DEFAULT_DEGREES),
        "rhos": list(DEFAULT_RHOS),
        "bernstein_x_budget_rows": budgets,
        "rho_2_degree_32_value_sup_budget": rho_2_degree_32["value_sup_norm_budget_on_first_x_panel"],
        "rho_2_degree_32_derivative_sup_budget": rho_2_degree_32[
            "derivative_sup_norm_budget_on_first_x_panel"
        ],
        "rho_3_degree_32_value_sup_budget": rho_3_degree_32["value_sup_norm_budget_on_first_x_panel"],
        "rho_3_degree_32_derivative_sup_budget": rho_3_degree_32[
            "derivative_sup_norm_budget_on_first_x_panel"
        ],
        "x_transform_formula": (
            "y=T*x^2, dy=2*T*x dx, and the Gamma(alpha+1,1) density becomes "
            "2*T^(alpha+1)*x^(2*alpha+1)*exp(-T*x^2)/Gamma(alpha+1) dx."
        ),
        "x_weight_power": x_weight_power,
        "x_moment_formula": (
            "For monomial x^k on the first panel, the transformed Gamma moment equals "
            "T^(-k/2)*lower_gamma(alpha+1+k/2,1)/Gamma(alpha+1)."
        ),
        "endpoint_parity_ready_to_apply_rows": endpoint_parity["summary"]["ready_to_apply_rows"],
        "tiny_mass_is_proof": False,
        "y_plane_branch_resolved": False,
        "x_panel_certificate_produced": False,
        "target_closing": False,
        "required_missing_input": (
            "Certify the first-panel value and derivative cores in x, including an x-domain/sup-norm "
            "or Taylor-model remainder bound and exact transformed moments; then splice this with "
            "the y-panel route away from zero."
        ),
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgexrm_01_endpoint_route_import",
            role="scope_reduction",
            readiness="diagnostic_only",
            claim=(
                "The endpoint parity-repair matrix identifies the x-variable first-panel route as an "
                "admissible way to avoid the y=0 sqrt-branch obstruction."
            ),
            diagnostics={
                "source_ready_to_apply_rows": diagnostics["endpoint_parity_ready_to_apply_rows"],
                "x_interval": diagnostics["x_interval"],
                "y_interval": diagnostics["y_interval"],
            },
            source_artifacts=[paths["endpoint_parity_json"], paths["endpoint_parity_note"]],
            proof_boundary="Scope import only; not an x-panel certificate or endpoint theorem.",
        ),
        MatrixRow(
            id="nlrgexrm_02_exact_x_change_of_variables",
            role="change_of_variables",
            readiness="open_requirement",
            claim=(
                "On the first panel, y=T*x^2 maps 0<=y<=1 to 0<=x<=0.01 and removes the "
                "sqrt(y/T) branch from the finite core by working directly in x."
            ),
            diagnostics={
                "x_transform_formula": diagnostics["x_transform_formula"],
                "x_weight_power": diagnostics["x_weight_power"],
                "x_moment_formula": diagnostics["x_moment_formula"],
            },
            source_artifacts=[paths["interpolation_route_json"], paths["finite_part_note"]],
            proof_boundary="Change-of-variables route only; transformed x moments and remainders are not certified here.",
        ),
        MatrixRow(
            id="nlrgexrm_03_first_panel_mass_budget",
            role="budget_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "The first panel has extremely small Gamma mass compared with the heaviest compact panel, "
                "so its admissible first-panel sup-norm budgets can be much larger after a valid x-panel "
                "remainder theorem is supplied."
            ),
            diagnostics={
                "first_panel_mass_upper": diagnostics["first_panel_mass_upper"],
                "heaviest_panel": diagnostics["heaviest_panel"],
                "heaviest_panel_mass_upper": diagnostics["heaviest_panel_mass_upper"],
                "first_to_heaviest_mass_ratio": diagnostics["first_to_heaviest_mass_ratio"],
                "first_panel_value_sup_budget_without_bernstein": diagnostics[
                    "first_panel_value_sup_budget_without_bernstein"
                ],
                "first_panel_derivative_sup_budget_without_bernstein": diagnostics[
                    "first_panel_derivative_sup_budget_without_bernstein"
                ],
            },
            source_artifacts=[paths["interpolation_route_json"], paths["quadrature_route_note"]],
            proof_boundary="Budget diagnostic only; tiny mass is not a proof of the first-panel remainder.",
        ),
        MatrixRow(
            id="nlrgexrm_04_x_bernstein_budget_table",
            role="route_budget",
            readiness="not_ready_to_apply",
            claim=(
                "For selected x-panel degrees and Bernstein radii, the table records sufficient-condition "
                "sup-norm budgets for value and derivative channels on the first panel."
            ),
            diagnostics={
                "degree_count": len(diagnostics["degrees"]),
                "rho_count": len(diagnostics["rhos"]),
                "bernstein_x_budget_rows": diagnostics["bernstein_x_budget_rows"],
                "rho_2_degree_32_value_sup_budget": diagnostics["rho_2_degree_32_value_sup_budget"],
                "rho_3_degree_32_value_sup_budget": diagnostics["rho_3_degree_32_value_sup_budget"],
            },
            source_artifacts=[paths["interpolation_route_note"], paths["compact_scout_note"]],
            proof_boundary="Route sizing table only; no x-domain or x-panel sup-norm theorem is proved.",
        ),
        MatrixRow(
            id="nlrgexrm_05_x_moment_taylor_obligation",
            role="open_requirement",
            readiness="open_requirement",
            claim=(
                "A usable x-panel certificate must integrate the x-local model against exact transformed "
                "moments and certify a true-function x Taylor/Chebyshev remainder in both value and derivative channels."
            ),
            diagnostics={
                "x_moment_formula": diagnostics["x_moment_formula"],
                "required_missing_input": diagnostics["required_missing_input"],
            },
            source_artifacts=[paths["arb_scout_note"], paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Open requirement only; the x moment implementation and remainder proof are not produced.",
        ),
        MatrixRow(
            id="nlrgexrm_06_tiny_mass_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The tiny first-panel Gamma mass alone proves the endpoint interpolation remainder is harmless.",
            diagnostics={
                "tiny_mass_is_proof": diagnostics["tiny_mass_is_proof"],
                "y_plane_branch_resolved": diagnostics["y_plane_branch_resolved"],
            },
            gap=(
                "Mass allocation only multiplies a future pointwise or analytic remainder bound. It does not "
                "supply the missing endpoint analyticity, x-domain, sup-norm, or Taylor-model certificate."
            ),
            source_artifacts=[paths["endpoint_parity_note"], paths["formal_tail_obstruction"]],
            proof_boundary="Rejected promotion only; small measure does not prove the endpoint theorem.",
        ),
        MatrixRow(
            id="nlrgexrm_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "Before the compact-interval route can be promoted, the first panel must be certified in x "
                "or replaced by an exact even-core plus tail-charge theorem, and then spliced to the nonzero y-panels."
            ),
            diagnostics={
                "ready_to_apply_rows": 0,
                "target_closing": diagnostics["target_closing"],
                "x_panel_certificate_produced": diagnostics["x_panel_certificate_produced"],
                "required_missing_input": diagnostics["required_missing_input"],
            },
            source_artifacts=[paths["dependency_graph"], paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Proof-hygiene gate only; not a compact certificate, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(interpolation_route_path: Path, endpoint_parity_path: Path) -> dict:
    interpolation_route = load_json(interpolation_route_path)
    endpoint_parity = load_json(endpoint_parity_path)
    paths = source_paths(interpolation_route_path, endpoint_parity_path)
    diagnostics = build_diagnostics(interpolation_route, endpoint_parity)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "alpha": diagnostics["alpha"],
        "x_interval": diagnostics["x_interval"],
        "y_interval": diagnostics["y_interval"],
        "x_right_endpoint": diagnostics["x_right_endpoint"],
        "compact_interval": diagnostics["compact_interval"],
        "panel_count": diagnostics["panel_count"],
        "first_panel_mass_upper": diagnostics["first_panel_mass_upper"],
        "heaviest_panel": diagnostics["heaviest_panel"],
        "heaviest_panel_mass_upper": diagnostics["heaviest_panel_mass_upper"],
        "first_to_heaviest_mass_ratio": diagnostics["first_to_heaviest_mass_ratio"],
        "value_unscaled_expectation_error_cap": diagnostics["value_unscaled_expectation_error_cap"],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "first_panel_value_sup_budget_without_bernstein": diagnostics[
            "first_panel_value_sup_budget_without_bernstein"
        ],
        "first_panel_derivative_sup_budget_without_bernstein": diagnostics[
            "first_panel_derivative_sup_budget_without_bernstein"
        ],
        "x_weight_power": diagnostics["x_weight_power"],
        "degree_count": len(diagnostics["degrees"]),
        "rho_count": len(diagnostics["rhos"]),
        "bernstein_x_budget_row_count": len(diagnostics["bernstein_x_budget_rows"]),
        "rho_2_degree_32_value_sup_budget": diagnostics["rho_2_degree_32_value_sup_budget"],
        "rho_2_degree_32_derivative_sup_budget": diagnostics[
            "rho_2_degree_32_derivative_sup_budget"
        ],
        "rho_3_degree_32_value_sup_budget": diagnostics["rho_3_degree_32_value_sup_budget"],
        "rho_3_degree_32_derivative_sup_budget": diagnostics[
            "rho_3_degree_32_derivative_sup_budget"
        ],
        "ready_to_apply_rows": 0,
        "target_closing": diagnostics["target_closing"],
        "tiny_mass_is_proof": diagnostics["tiny_mass_is_proof"],
        "y_plane_branch_resolved": diagnostics["y_plane_branch_resolved"],
        "x_panel_certificate_produced": diagnostics["x_panel_certificate_produced"],
        "required_missing_input": diagnostics["required_missing_input"],
        "main_finding": (
            "The endpoint x-panel route matrix quantifies the first-panel repair suggested by the "
            "endpoint parity matrix. The change y=T*x^2 maps 0<=y<=1 to 0<=x<=0.01 and gives "
            "transformed density power x^42 for alpha=20.5. The first-panel Gamma mass upper is "
            f"{diagnostics['first_panel_mass_upper']}, only "
            f"{diagnostics['first_to_heaviest_mass_ratio']} of the heaviest panel mass. This makes "
            "the first panel numerically lightweight, but it remains a route matrix only: no x-domain, "
            "sup-norm, exact-moment implementation, or true x-panel interpolation-remainder certificate "
            "is proved."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix",
        "date": "2026-07-08",
        "status": "endpoint x-panel route matrix",
        "source_interpolation_remainder_route_matrix": paths["interpolation_route_note"],
        "source_interpolation_remainder_route_json": paths["interpolation_route_json"],
        "source_endpoint_parity_repair_matrix": paths["endpoint_parity_note"],
        "source_endpoint_parity_repair_json": paths["endpoint_parity_json"],
        "source_arb_chebyshev_interpolant_moment_scout": paths["arb_scout_note"],
        "source_compact_interval_integration_scout": paths["compact_scout_note"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py"
        ),
        "proof_boundary": (
            "Endpoint x-panel route matrix only. It records the y=T*x^2 change of variables, "
            "first-panel Gamma mass and Bernstein sufficient-condition budgets, exact-moment "
            "obligations, and rejected shortcuts. It does not prove an x-domain theorem, does not "
            "prove a sup-norm bound, does not implement exact x moments, does not prove a "
            "Taylor-model or Chebyshev interpolation remainder, does not prove a compact interval "
            "certificate, does not prove a finite-grid interval certificate, does not prove a uniform "
            "collar theorem, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The first-panel mass is not promoted into a proof.",
            "The y=0 branch is avoided only after a certified x-panel theorem is supplied.",
            "Bernstein x-panel budgets are sufficient-condition route sizing, not proved estimates.",
            "The x moment formula is an obligation until implemented and checked with Arb-safe enclosures.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"x interval {summary['x_interval']}, "
        f"{summary['bernstein_x_budget_row_count']} Bernstein budgets, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    budget_rows_text = artifact["matrix_rows"][3]["diagnostics"]["bernstein_x_budget_rows"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint X-Panel Route Matrix",
        "",
        "Date: 2026-07-08",
        "",
        "Status: endpoint x-panel route matrix. This is not a proof",
        "of an x-panel certificate, an interpolation-remainder theorem, a compact interval",
        "certificate, a finite-grid interval certificate, a uniform collar theorem, RH,",
        "or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix`.",
        "",
        "Proof boundary: this artifact records the endpoint x change of variables,",
        "first-panel mass and Bernstein budgets, exact-moment obligations, and rejected",
        "shortcuts. It does not prove the x-panel remainder.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Route Constants",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"alpha: {summary['alpha']}",
        f"y interval: {summary['y_interval']}",
        f"x interval: {summary['x_interval']}",
        f"x right endpoint: {summary['x_right_endpoint']}",
        f"x weight power: {summary['x_weight_power']}",
        f"first panel mass upper: {summary['first_panel_mass_upper']}",
        f"heaviest panel: {summary['heaviest_panel']}",
        f"heaviest panel mass upper: {summary['heaviest_panel_mass_upper']}",
        f"first/heaviest mass ratio: {summary['first_to_heaviest_mass_ratio']}",
        f"value unscaled cap: {summary['value_unscaled_expectation_error_cap']}",
        f"derivative unscaled cap: {summary['derivative_unscaled_expectation_error_cap']}",
        f"first-panel value sup budget without Bernstein: {summary['first_panel_value_sup_budget_without_bernstein']}",
        f"first-panel derivative sup budget without Bernstein: {summary['first_panel_derivative_sup_budget_without_bernstein']}",
        f"rho=2 degree=32 value sup budget: {summary['rho_2_degree_32_value_sup_budget']}",
        f"rho=3 degree=32 value sup budget: {summary['rho_3_degree_32_value_sup_budget']}",
        "```",
        "",
        "Endpoint x transform:",
        "",
        "```text",
        artifact["matrix_rows"][1]["diagnostics"]["x_transform_formula"],
        artifact["matrix_rows"][1]["diagnostics"]["x_moment_formula"],
        "```",
        "",
        "Selected first x-panel Bernstein budgets:",
        "",
        "```text",
    ]
    for row in budget_rows_text:
        if row["degree"] in {16, 24, 32}:
            lines.append(
                f"N={row['degree']}, rho={row['rho']}: "
                f"value M<={row['value_sup_norm_budget_on_first_x_panel']}; "
                f"derivative M<={row['derivative_sup_norm_budget_on_first_x_panel']}"
            )
    lines.extend(
        [
            "```",
            "",
            "Rejected shortcut:",
            "",
            "```text",
            "Tiny first-panel Gamma mass is not a proof of the endpoint remainder.",
            "```",
            "",
            "Required upgrade:",
            "",
            "```text",
            summary["required_missing_input"],
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_interpolation_remainder_route_matrix"],
            artifact["source_interpolation_remainder_route_json"],
            artifact["source_endpoint_parity_repair_matrix"],
            artifact["source_endpoint_parity_repair_json"],
            artifact["source_arb_chebyshev_interpolant_moment_scout"],
            artifact["source_compact_interval_integration_scout"],
            artifact["source_quadrature_remainder_route_matrix"],
            artifact["source_finite_part_weighted_sum_interval_certificate"],
            artifact["source_intervalization_target"],
            artifact["source_uniform_remainder_target"],
            artifact["source_dependency_graph"],
            artifact["source_formal_tail_obstruction"],
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
    parser.add_argument("--interpolation-route-json", type=Path, default=DEFAULT_INTERPOLATION_ROUTE_JSON)
    parser.add_argument("--endpoint-parity-json", type=Path, default=DEFAULT_ENDPOINT_PARITY_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    interpolation_route_path = (
        args.interpolation_route_json
        if args.interpolation_route_json.is_absolute()
        else REPO_ROOT / args.interpolation_route_json
    )
    endpoint_parity_path = (
        args.endpoint_parity_json
        if args.endpoint_parity_json.is_absolute()
        else REPO_ROOT / args.endpoint_parity_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(interpolation_route_path, endpoint_parity_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
