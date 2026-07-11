#!/usr/bin/env python3
"""Build a panel interpolation-remainder route matrix for the worst row."""

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


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ARB_SCOUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.json"
)
DEFAULT_QUADRATURE_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_PRECISION_BITS = 384
DEFAULT_PANELS = ((0, 1), (1, 5), (5, 20), (20, 50), (50, 100), (100, 200))
DEFAULT_DEGREES = (32, 64, 96, 128, 160)
DEFAULT_RHOS = ("1.25", "1.5", "2.0", "3.0")
DEFAULT_ASSUMED_SUP_NORMS = ("1e-12", "1e-6", "1", "1e6")
DEFAULT_MIN_DEGREE_START = 16
DEFAULT_MIN_DEGREE_LIMIT = 256

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


def sci(value: Decimal, digits: int = 18) -> str:
    return f"{value:.{digits}E}"


def fixed(value: Decimal) -> str:
    return format(value, "f")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(arb_scout_path: Path, quadrature_route_path: Path) -> dict[str, str]:
    return {
        "arb_scout_json": arb_scout_path.relative_to(REPO_ROOT).as_posix(),
        "arb_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md"
        ),
        "quadrature_route_json": quadrature_route_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"
        ),
        "floating_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md"
        ),
        "compact_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md"
        ),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
        ),
        "intervalization_target": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def bernstein_constant(degree: int, rho: Decimal) -> Decimal:
    return Decimal(4) * (rho ** Decimal(-degree)) / (rho - Decimal(1))


def panel_masses() -> list[dict]:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    alpha = flint.arb(DEFAULT_ALPHA)
    gamma_alpha = (alpha + 1).gamma()
    rows: list[dict] = []
    for left, right in DEFAULT_PANELS:
        mass = (flint.arb(left).gamma_upper(alpha + 1) - flint.arb(right).gamma_upper(alpha + 1)) / gamma_alpha
        rows.append(
            {
                "panel": f"{left}<=y<={right}",
                "left": left,
                "right": right,
                "gamma_mass_upper": mass.upper().str(50, radius=False).replace("e", "E"),
            }
        )
    return rows


def budget_rows(value_cap: Decimal, derivative_cap: Decimal, heaviest_mass: Decimal) -> list[dict]:
    rows: list[dict] = []
    panel_count = Decimal(len(DEFAULT_PANELS))
    for degree in DEFAULT_DEGREES:
        for rho_text in DEFAULT_RHOS:
            rho = Decimal(rho_text)
            constant = bernstein_constant(degree, rho)
            value_sup_budget = value_cap / (panel_count * heaviest_mass * constant)
            derivative_sup_budget = derivative_cap / (panel_count * heaviest_mass * constant)
            rows.append(
                {
                    "degree": degree,
                    "rho": rho_text,
                    "bernstein_error_constant": sci(constant, 30),
                    "value_sup_norm_budget_on_heaviest_panel": sci(value_sup_budget, 30),
                    "derivative_sup_norm_budget_on_heaviest_panel": sci(derivative_sup_budget, 30),
                    "proof_boundary": (
                        "Heaviest-panel sufficient-condition budget only; requires an independent "
                        "analytic-domain and sup-norm certificate."
                    ),
                }
            )
    return rows


def minimal_degree_rows(value_cap: Decimal, derivative_cap: Decimal, heaviest_mass: Decimal) -> list[dict]:
    rows: list[dict] = []
    panel_count = Decimal(len(DEFAULT_PANELS))
    for assumed_text in DEFAULT_ASSUMED_SUP_NORMS:
        assumed = Decimal(assumed_text)
        for rho_text in DEFAULT_RHOS:
            rho = Decimal(rho_text)
            found = None
            value_budget_at_found = None
            derivative_budget_at_found = None
            for degree in range(DEFAULT_MIN_DEGREE_START, DEFAULT_MIN_DEGREE_LIMIT + 1):
                constant = bernstein_constant(degree, rho)
                value_sup_budget = value_cap / (panel_count * heaviest_mass * constant)
                derivative_sup_budget = derivative_cap / (panel_count * heaviest_mass * constant)
                if value_sup_budget >= assumed and derivative_sup_budget >= assumed:
                    found = degree
                    value_budget_at_found = value_sup_budget
                    derivative_budget_at_found = derivative_sup_budget
                    break
            rows.append(
                {
                    "assumed_joint_sup_norm": assumed_text,
                    "rho": rho_text,
                    "minimal_degree_in_16_to_256": found,
                    "value_sup_norm_budget_at_minimal_degree": (
                        None if value_budget_at_found is None else sci(value_budget_at_found, 24)
                    ),
                    "derivative_sup_norm_budget_at_minimal_degree": (
                        None if derivative_budget_at_found is None else sci(derivative_budget_at_found, 24)
                    ),
                    "proof_boundary": (
                        "Conditional degree sizing only; it is not a proof that the transformed panel "
                        "core has this rho-domain or sup norm."
                    ),
                }
            )
    return rows


def build_diagnostics(arb_scout: dict, quadrature_route: dict) -> dict:
    value_cap = dec(quadrature_route["summary"]["value_unscaled_expectation_error_cap"])
    derivative_cap = dec(quadrature_route["summary"]["derivative_unscaled_expectation_error_cap"])
    masses = panel_masses()
    heaviest = max(masses, key=lambda row: dec(row["gamma_mass_upper"]))
    heaviest_mass = dec(heaviest["gamma_mass_upper"])
    budgets = budget_rows(value_cap, derivative_cap, heaviest_mass)
    minimal = minimal_degree_rows(value_cap, derivative_cap, heaviest_mass)
    rho_2_degree_128 = next(row for row in budgets if row["rho"] == "2.0" and row["degree"] == 128)
    rho_2_degree_160 = next(row for row in budgets if row["rho"] == "2.0" and row["degree"] == 160)
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "alpha": DEFAULT_ALPHA,
        "compact_interval": arb_scout["summary"]["compact_interval"],
        "arb_reference_degree": arb_scout["summary"]["reference_degree"],
        "arb_cap_safe_pair_count": arb_scout["summary"]["cap_safe_pair_count"],
        "arb_reference_value_interpolant_mid": arb_scout["summary"]["reference_value_interpolant_mid"],
        "arb_reference_derivative_interpolant_mid": arb_scout["summary"]["reference_derivative_interpolant_mid"],
        "value_unscaled_expectation_error_cap": quadrature_route["summary"][
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": quadrature_route["summary"][
            "derivative_unscaled_expectation_error_cap"
        ],
        "panel_count": len(DEFAULT_PANELS),
        "panel_masses": masses,
        "heaviest_panel": heaviest["panel"],
        "heaviest_panel_mass_upper": heaviest["gamma_mass_upper"],
        "degrees": list(DEFAULT_DEGREES),
        "rhos": list(DEFAULT_RHOS),
        "assumed_sup_norms": list(DEFAULT_ASSUMED_SUP_NORMS),
        "bernstein_budget_rows": budgets,
        "minimal_degree_rows": minimal,
        "rho_2_degree_128_value_sup_budget": rho_2_degree_128["value_sup_norm_budget_on_heaviest_panel"],
        "rho_2_degree_128_derivative_sup_budget": rho_2_degree_128[
            "derivative_sup_norm_budget_on_heaviest_panel"
        ],
        "rho_2_degree_160_value_sup_budget": rho_2_degree_160["value_sup_norm_budget_on_heaviest_panel"],
        "rho_2_degree_160_derivative_sup_budget": rho_2_degree_160[
            "derivative_sup_norm_budget_on_heaviest_panel"
        ],
        "route_formula": "panel_error <= gamma_panel_mass * 4*rho^(-N)/(rho-1) * M",
        "required_missing_input": (
            "For every panel, certify an analytic Bernstein-ellipse domain for the transformed "
            "cancellation-reduced value and derivative cores and an Arb-safe sup norm M on that domain."
        ),
        "endpoint_route_condition": (
            "The first panel touches y=0 while the finite core is evaluated through sqrt(y/T); a "
            "proof must either work in the x-variable near zero or certify the even analytic y-core "
            "after separating the tail/parity terms."
        ),
        "arb_cauchy_delta_is_proof": False,
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrirm_01_arb_interpolant_import",
            role="scope_reduction",
            readiness="diagnostic_only",
            claim=(
                "The Arb Chebyshev interpolant-moment scout supplies ball arithmetic for panel "
                "interpolants, but explicitly leaves the true-function interpolation remainder open."
            ),
            diagnostics={
                "compact_interval": diagnostics["compact_interval"],
                "arb_reference_degree": diagnostics["arb_reference_degree"],
                "arb_cap_safe_pair_count": diagnostics["arb_cap_safe_pair_count"],
                "reference_value_interpolant_mid": diagnostics["arb_reference_value_interpolant_mid"],
                "reference_derivative_interpolant_mid": diagnostics[
                    "arb_reference_derivative_interpolant_mid"
                ],
            },
            source_artifacts=[paths["arb_scout_json"], paths["arb_scout_note"], paths["floating_scout_note"]],
            proof_boundary="Scope import only; not an interpolation-remainder or compact-integral proof.",
        ),
        MatrixRow(
            id="nlrgwrirm_02_panel_mass_allocation",
            role="budget_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "The compact interval 0<=y<=200 is split into six Gamma-weighted panels; the heaviest "
                "panel for alpha=20.5 is 20<=y<=50."
            ),
            diagnostics={
                "panel_count": diagnostics["panel_count"],
                "panel_masses": diagnostics["panel_masses"],
                "heaviest_panel": diagnostics["heaviest_panel"],
                "heaviest_panel_mass_upper": diagnostics["heaviest_panel_mass_upper"],
                "equal_panel_allocation": "Each panel is assigned one sixth of the value/derivative caps.",
            },
            source_artifacts=[paths["compact_scout_note"], paths["far_tail_note"], paths["quadrature_route_json"]],
            proof_boundary="Budget allocation only; it does not bound any interpolation remainder.",
        ),
        MatrixRow(
            id="nlrgwrirm_03_bernstein_sufficient_condition",
            role="sufficient_condition",
            readiness="open_requirement",
            claim=(
                "A conservative Bernstein-ellipse route would bound degree-N panel interpolation error "
                "by gamma_panel_mass * 4*rho^(-N)/(rho-1) * M once an analytic rho-domain and sup norm "
                "M are independently certified."
            ),
            diagnostics={
                "route_formula": diagnostics["route_formula"],
                "required_missing_input": diagnostics["required_missing_input"],
                "value_cap": diagnostics["value_unscaled_expectation_error_cap"],
                "derivative_cap": diagnostics["derivative_unscaled_expectation_error_cap"],
            },
            source_artifacts=[paths["quadrature_route_note"], paths["uniform_remainder_target"]],
            proof_boundary="Sufficient-condition formula only; no analytic-domain or sup-norm certificate is proved.",
        ),
        MatrixRow(
            id="nlrgwrirm_04_heaviest_panel_budget_table",
            role="route_budget",
            readiness="not_ready_to_apply",
            claim=(
                "The heaviest-panel budget table quantifies how large a certified panel sup norm M could "
                "be for selected degrees and Bernstein radii while staying under the unscaled caps."
            ),
            diagnostics={
                "degree_count": len(diagnostics["degrees"]),
                "rho_count": len(diagnostics["rhos"]),
                "bernstein_budget_rows": diagnostics["bernstein_budget_rows"],
                "rho_2_degree_128_value_sup_budget": diagnostics["rho_2_degree_128_value_sup_budget"],
                "rho_2_degree_160_value_sup_budget": diagnostics["rho_2_degree_160_value_sup_budget"],
            },
            source_artifacts=[paths["arb_scout_note"], paths["quadrature_route_json"]],
            proof_boundary="Route sizing table only; not a panel sup-norm or interpolation-remainder theorem.",
        ),
        MatrixRow(
            id="nlrgwrirm_05_minimal_degree_table",
            role="route_budget",
            readiness="not_ready_to_apply",
            claim=(
                "For assumed joint sup norms 1e-12, 1e-6, 1, and 1e6, the table records the first "
                "degree in 16..256 that would satisfy both value and derivative heaviest-panel budgets."
            ),
            diagnostics={
                "assumed_sup_norms": diagnostics["assumed_sup_norms"],
                "rhos": diagnostics["rhos"],
                "minimal_degree_rows": diagnostics["minimal_degree_rows"],
            },
            source_artifacts=[paths["quadrature_route_note"], paths["intervalization_target"]],
            proof_boundary="Conditional degree sizing only; it is not a certified analytic estimate.",
        ),
        MatrixRow(
            id="nlrgwrirm_06_endpoint_and_branch_condition",
            role="route_constraint",
            readiness="open_requirement",
            claim=(
                "The endpoint y=0 and sqrt(y/T) composition make the first panel a branch-sensitive "
                "case unless the proof works in x or establishes an even analytic y-core plus tail split."
            ),
            diagnostics={"endpoint_route_condition": diagnostics["endpoint_route_condition"]},
            source_artifacts=[paths["finite_part_note"], paths["uniform_remainder_target"]],
            proof_boundary="Route constraint only; no endpoint analytic continuation theorem is proved.",
        ),
        MatrixRow(
            id="nlrgwrirm_07_arb_cauchy_delta_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim=(
                "The Arb Cauchy deltas between degree-16,20,24,32 interpolant integrals prove the "
                "true-function interpolation remainder is below the caps."
            ),
            diagnostics={
                "arb_cauchy_delta_is_proof": diagnostics["arb_cauchy_delta_is_proof"],
                "cap_safe_pair_count": diagnostics["arb_cap_safe_pair_count"],
            },
            gap=(
                "Arb Cauchy deltas compare interpolants to interpolants. They do not compare the "
                "interpolant to the true cancellation-reduced core, and therefore do not prove an "
                "interpolation remainder bound."
            ),
            source_artifacts=[paths["arb_scout_json"], paths["arb_scout_note"]],
            proof_boundary="Rejected promotion only; Arb Cauchy deltas do not prove the interpolation remainder.",
        ),
        MatrixRow(
            id="nlrgwrirm_08_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted compact-interval certificate may reuse the Arb interpolant arithmetic only "
                "after panel analytic-domain/sup-norm or Taylor-model remainder certificates close the "
                "interpolation-remainder gap."
            ),
            diagnostics={
                "ready_to_apply_rows": 0,
                "target_closing": diagnostics["target_closing"],
                "required_missing_input": diagnostics["required_missing_input"],
            },
            source_artifacts=[paths["dependency_graph"], paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Proof-hygiene gate only; not a compact certificate, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(arb_scout_path: Path, quadrature_route_path: Path) -> dict:
    arb_scout = load_json(arb_scout_path)
    quadrature_route = load_json(quadrature_route_path)
    paths = source_paths(arb_scout_path, quadrature_route_path)
    diagnostics = build_diagnostics(arb_scout, quadrature_route)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "alpha": diagnostics["alpha"],
        "compact_interval": diagnostics["compact_interval"],
        "panel_count": diagnostics["panel_count"],
        "heaviest_panel": diagnostics["heaviest_panel"],
        "heaviest_panel_mass_upper": diagnostics["heaviest_panel_mass_upper"],
        "degree_count": len(diagnostics["degrees"]),
        "rho_count": len(diagnostics["rhos"]),
        "bernstein_budget_row_count": len(diagnostics["bernstein_budget_rows"]),
        "minimal_degree_row_count": len(diagnostics["minimal_degree_rows"]),
        "value_unscaled_expectation_error_cap": diagnostics["value_unscaled_expectation_error_cap"],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "rho_2_degree_128_value_sup_budget": diagnostics["rho_2_degree_128_value_sup_budget"],
        "rho_2_degree_128_derivative_sup_budget": diagnostics[
            "rho_2_degree_128_derivative_sup_budget"
        ],
        "rho_2_degree_160_value_sup_budget": diagnostics["rho_2_degree_160_value_sup_budget"],
        "rho_2_degree_160_derivative_sup_budget": diagnostics[
            "rho_2_degree_160_derivative_sup_budget"
        ],
        "ready_to_apply_rows": 0,
        "target_closing": diagnostics["target_closing"],
        "arb_cauchy_delta_is_proof": diagnostics["arb_cauchy_delta_is_proof"],
        "required_missing_input": diagnostics["required_missing_input"],
        "endpoint_route_condition": diagnostics["endpoint_route_condition"],
        "main_finding": (
            "The interpolation-remainder route matrix isolates the remaining compact-interval theorem: "
            "the heaviest Gamma panel is 20<=y<=50 with mass upper "
            f"{diagnostics['heaviest_panel_mass_upper']}, and the Bernstein route formula "
            "panel_error <= mass * 4*rho^(-N)/(rho-1) * M gives concrete degree/rho budgets. "
            "At rho=2 and degree 128 the heaviest-panel value sup-norm budget is "
            f"{diagnostics['rho_2_degree_128_value_sup_budget']}, while degree 160 raises it to "
            f"{diagnostics['rho_2_degree_160_value_sup_budget']}. This is still a route matrix only: "
            "no analytic-domain, sup-norm, endpoint, or true interpolation-remainder certificate is proved."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix",
        "date": "2026-07-08",
        "status": "worst-row interpolation-remainder route matrix",
        "source_arb_chebyshev_interpolant_moment_scout": paths["arb_scout_note"],
        "source_arb_chebyshev_interpolant_moment_scout_json": paths["arb_scout_json"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_quadrature_remainder_route_json": paths["quadrature_route_json"],
        "source_floating_chebyshev_panel_moment_scout": paths["floating_scout_note"],
        "source_compact_interval_integration_scout": paths["compact_scout_note"],
        "source_far_tail_split_certificate": paths["far_tail_note"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py"
        ),
        "proof_boundary": (
            "Worst-row interpolation-remainder route matrix only. It records panel Gamma masses, "
            "Bernstein-ellipse sufficient-condition budgets, endpoint route constraints, and rejected "
            "shortcuts for the compact interval 0<=y<=200. It does not prove an analytic Bernstein "
            "ellipse domain, does not prove a sup-norm bound, does not prove a Taylor-model remainder, "
            "does not prove the true compact integral, does not prove a finite-grid interval certificate, "
            "does not prove a uniform collar theorem, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Bernstein budgets are sufficient-condition route sizing, not proved estimates.",
            "Arb Cauchy deltas do not prove the interpolation remainder.",
            "The endpoint y=0 branch/composition condition remains open.",
            "The compact interval 0<=y<=200 remains open until true-function panel remainders are certified.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row "
        f"interpolation-remainder route matrix: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['panel_count']} panel masses, "
        f"{summary['bernstein_budget_row_count']} Bernstein budgets, "
        f"{summary['minimal_degree_row_count']} minimal-degree rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    budget_rows_text = artifact["matrix_rows"][3]["diagnostics"]["bernstein_budget_rows"]
    minimal_rows_text = artifact["matrix_rows"][4]["diagnostics"]["minimal_degree_rows"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Interpolation-Remainder Route Matrix",
        "",
        "Date: 2026-07-08",
        "",
        "Status: worst-row interpolation-remainder route matrix. This is not a proof",
        "of an interpolation-remainder theorem, compact interval-integration certificate,",
        "finite-grid interval certificate, uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix`.",
        "",
        "Proof boundary: this artifact records panel Gamma masses, Bernstein",
        "sufficient-condition budgets, endpoint route constraints, and rejected",
        "shortcuts. It does not prove an analytic-domain bound or a true",
        "interpolation remainder.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py",
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
        f"compact interval: {summary['compact_interval']}",
        f"panel count: {summary['panel_count']}",
        f"heaviest panel: {summary['heaviest_panel']}",
        f"heaviest panel mass upper: {summary['heaviest_panel_mass_upper']}",
        f"value unscaled cap: {summary['value_unscaled_expectation_error_cap']}",
        f"derivative unscaled cap: {summary['derivative_unscaled_expectation_error_cap']}",
        f"rho=2 degree=128 value sup budget: {summary['rho_2_degree_128_value_sup_budget']}",
        f"rho=2 degree=160 value sup budget: {summary['rho_2_degree_160_value_sup_budget']}",
        f"target closing: {summary['target_closing']}",
        "```",
        "",
        "Bernstein route formula:",
        "",
        "```text",
        "panel_error <= gamma_panel_mass * 4*rho^(-N)/(rho-1) * M",
        "```",
        "",
        "Selected heaviest-panel budgets:",
        "",
        "```text",
    ]
    for row in budget_rows_text:
        if row["rho"] == "2.0":
            lines.append(
                f"N={row['degree']}, rho={row['rho']}: value M<={row['value_sup_norm_budget_on_heaviest_panel']}; "
                f"derivative M<={row['derivative_sup_norm_budget_on_heaviest_panel']}"
            )
    lines.extend(
        [
            "```",
            "",
            "Minimal degree rows:",
            "",
            "```text",
        ]
    )
    for row in minimal_rows_text:
        if row["rho"] in {"2.0", "3.0"}:
            lines.append(
                f"M={row['assumed_joint_sup_norm']}, rho={row['rho']}: "
                f"N={row['minimal_degree_in_16_to_256']}"
            )
    lines.extend(
        [
            "```",
            "",
            "Rejected shortcut:",
            "",
            "```text",
            "Arb Cauchy deltas do not prove the interpolation remainder.",
            "```",
            "",
            "Required upgrade:",
            "",
            "```text",
            summary["required_missing_input"],
            summary["endpoint_route_condition"],
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_arb_chebyshev_interpolant_moment_scout"],
            artifact["source_arb_chebyshev_interpolant_moment_scout_json"],
            artifact["source_quadrature_remainder_route_matrix"],
            artifact["source_quadrature_remainder_route_json"],
            artifact["source_floating_chebyshev_panel_moment_scout"],
            artifact["source_compact_interval_integration_scout"],
            artifact["source_far_tail_split_certificate"],
            artifact["source_finite_part_weighted_sum_interval_certificate"],
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
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arb-scout-json", type=Path, default=DEFAULT_ARB_SCOUT_JSON)
    parser.add_argument("--quadrature-route-json", type=Path, default=DEFAULT_QUADRATURE_ROUTE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    arb_scout_path = args.arb_scout_json if args.arb_scout_json.is_absolute() else REPO_ROOT / args.arb_scout_json
    quadrature_route_path = (
        args.quadrature_route_json
        if args.quadrature_route_json.is_absolute()
        else REPO_ROOT / args.quadrature_route_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(arb_scout_path, quadrature_route_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row "
        f"interpolation-remainder route matrix: {out_json.relative_to(REPO_ROOT).as_posix()} and "
        f"{note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
