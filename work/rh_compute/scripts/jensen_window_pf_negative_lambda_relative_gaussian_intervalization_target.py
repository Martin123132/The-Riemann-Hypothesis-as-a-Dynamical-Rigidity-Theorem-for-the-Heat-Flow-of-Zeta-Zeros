#!/usr/bin/env python3
"""Build the intervalization target for the cancellation-reduced remainder grid."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md"

getcontext().prec = 80


@dataclass(frozen=True)
class CertificationBudgetRow:
    channel: str
    observed_worst_ratio: str
    observed_worst_location: str
    slack_to_one: str
    half_slack: str
    proposed_total_ratio_error_cap: str
    closed_if_cap_met: bool
    proof_boundary: str


@dataclass(frozen=True)
class ObligationRow:
    id: str
    role: str
    readiness: str
    claim: str
    required_upgrade: str
    target_error_cap: str | None
    proof_boundary: str


def dec(value: object) -> Decimal:
    return Decimal(str(value))


def fmt(value: Decimal, places: int = 18) -> str:
    return format(value, f".{places}E")


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(grid_path: Path) -> dict[str, str]:
    return {
        "grid_json": grid_path.relative_to(REPO_ROOT).as_posix(),
        "grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
        "actual_endpoint_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md",
        "asymptotic_target_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md",
        "formal_tail_obstruction_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"
        ),
        "residual_budget_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def build_budget_rows(grid: dict) -> tuple[list[CertificationBudgetRow], dict]:
    summary = grid["summary"]
    value_ratio = dec(summary["max_value_ratio_to_first_omitted_over_orders"])
    derivative_ratio = dec(summary["max_derivative_ratio_to_first_omitted_over_orders"])
    cap = Decimal("1.0E-2")
    rows: list[CertificationBudgetRow] = []
    for channel, ratio, location in (
        ("value", value_ratio, summary["max_value_ratio_location"]),
        ("derivative", derivative_ratio, summary["max_derivative_ratio_location"]),
    ):
        slack = Decimal(1) - ratio
        half_slack = slack / Decimal(2)
        rows.append(
            CertificationBudgetRow(
                channel=channel,
                observed_worst_ratio=fmt(ratio),
                observed_worst_location=f"T={location['T']}, F_{location['index']}",
                slack_to_one=fmt(slack),
                half_slack=fmt(half_slack),
                proposed_total_ratio_error_cap=fmt(cap),
                closed_if_cap_met=bool(ratio + cap < Decimal(1)),
                proof_boundary=(
                    "Certification budget only; the intervalized quadrature and analytic error terms are not yet proved."
                ),
            )
        )
    common_slack = min(Decimal(1) - value_ratio, Decimal(1) - derivative_ratio)
    common_half_slack = common_slack / Decimal(2)
    channel_count = Decimal(5)
    return rows, {
        "value_slack_to_one": fmt(Decimal(1) - value_ratio),
        "derivative_slack_to_one": fmt(Decimal(1) - derivative_ratio),
        "common_slack_to_one": fmt(common_slack),
        "common_half_slack": fmt(common_half_slack),
        "proposed_total_ratio_error_cap": fmt(cap),
        "proposed_per_error_source_cap_for_five_sources": fmt(cap / channel_count),
        "closed_if_total_cap_met": bool(value_ratio + cap < Decimal(1) and derivative_ratio + cap < Decimal(1)),
    }


def build_obligations(diagnostics: dict) -> list[ObligationRow]:
    per_source_cap = diagnostics["proposed_per_error_source_cap_for_five_sources"]
    total_cap = diagnostics["proposed_total_ratio_error_cap"]
    return [
        ObligationRow(
            id="nlrgit_01_residual_core_identity",
            role="exact_reduction",
            readiness="available_exact",
            claim="Use the cancellation-reduced Gamma-expectation residual cores for value and derivative channels.",
            required_upgrade="None for the identity; downstream proof must preserve the same normalization and degree-40 polynomial.",
            target_error_cap=None,
            proof_boundary="Exact coordinate identity only; not a numerical enclosure.",
        ),
        ObligationRow(
            id="nlrgit_02_laguerre_node_weight_intervals",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim="Replace SciPy floating generalized Laguerre nodes and weights by interval enclosures for every recorded index and order.",
            required_upgrade="Certified roots and Christoffel weights for L_N^(i-1/2), including root separation and weight-radius propagation.",
            target_error_cap=per_source_cap,
            proof_boundary="Open numerical-certification requirement only.",
        ),
        ObligationRow(
            id="nlrgit_03_phi_and_c0_interval_tail",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim="Enclose Phi, Phi', Phi(0), and their n>30 tails on all node-induced v ranges used by the finite grid.",
            required_upgrade="Arb or analytic exponential-tail enclosure for value and derivative cores after the n=30 truncation.",
            target_error_cap=per_source_cap,
            proof_boundary="Open analytic/numerical tail requirement only.",
        ),
        ObligationRow(
            id="nlrgit_04_quadrature_remainder_error",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim="Bound the generalized Gauss-Laguerre quadrature remainder for the cancellation-reduced cores, not just the spread across orders.",
            required_upgrade="A rigorous quadrature remainder theorem or an interval adaptive integral cross-check for the Gamma expectation.",
            target_error_cap=per_source_cap,
            proof_boundary="Open quadrature-error requirement only.",
        ),
        ObligationRow(
            id="nlrgit_05_ratio_and_coefficient_ball_propagation",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim="Propagate Arb coefficient-ratio balls and first-omitted-term denominators through the ratio comparisons.",
            required_upgrade="Interval-safe denominator lower bounds and signed residual absolute-value enclosure in first-omitted units.",
            target_error_cap=per_source_cap,
            proof_boundary="Open interval propagation requirement only.",
        ),
        ObligationRow(
            id="nlrgit_06_rounding_and_aggregation_budget",
            role="open_requirement",
            readiness="not_ready_to_apply",
            claim="Account for arithmetic rounding, summation, and cross-source aggregation so the total ratio error remains below the proposed cap.",
            required_upgrade="A machine-checkable error ledger whose total ratio radius is below the recorded total cap.",
            target_error_cap=total_cap,
            proof_boundary="Open certification-ledger requirement only.",
        ),
        ObligationRow(
            id="nlrgit_07_finite_grid_not_uniform_collar",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="An intervalized finite T grid would by itself prove the full collar residual theorem on 0<=u<=1/1156.",
            required_upgrade="A separate T-continuity, monotonicity, interval subdivision, or far-tail theorem is still required.",
            target_error_cap=None,
            proof_boundary="Rejected promotion only; finite grid certification is not a uniform collar theorem.",
        ),
        ObligationRow(
            id="nlrgit_08_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim="A promoted proof must state all interval sources, cap their total ratio error, and separately bridge from finite grid certification to the full collar or all-k theorem.",
            required_upgrade="Intervalized grid certificate plus a continuum-in-T theorem or subdivision plan with certified coverage.",
            target_error_cap=total_cap,
            proof_boundary="Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        ),
    ]


def build_diagnostics(grid: dict) -> dict:
    budget_rows, budget_summary = build_budget_rows(grid)
    grid_summary = grid["summary"]
    obligations = build_obligations(budget_summary)
    return {
        "source_grid_summary": {
            "grid_rows": grid_summary["grid_rows"],
            "t_grid_count": grid_summary["t_grid_count"],
            "index_count": grid_summary["index_count"],
            "quadrature_order_count": grid_summary["quadrature_order_count"],
            "selected_quadrature_order": grid_summary["selected_quadrature_order"],
            "max_value_ratio_to_first_omitted_over_orders": grid_summary[
                "max_value_ratio_to_first_omitted_over_orders"
            ],
            "max_derivative_ratio_to_first_omitted_over_orders": grid_summary[
                "max_derivative_ratio_to_first_omitted_over_orders"
            ],
            "max_value_ratio_spread_over_orders": grid_summary["max_value_ratio_spread_over_orders"],
            "max_derivative_ratio_spread_over_orders": grid_summary["max_derivative_ratio_spread_over_orders"],
        },
        "certification_budget_rows": [asdict(row) for row in budget_rows],
        "certification_budget_row_count": len(budget_rows),
        "budget_summary": budget_summary,
        "obligation_rows": [asdict(row) for row in obligations],
        "obligation_row_count": len(obligations),
        "open_requirement_rows": sum(row.role == "open_requirement" for row in obligations),
        "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in obligations),
        "proof_boundary_note": (
            "This target is a certification roadmap for the finite cancellation-reduced grid. It does not "
            "provide node/weight intervals, quadrature-remainder bounds, Phi n-tail bounds, or a continuum-in-T theorem."
        ),
    }


def build_artifact(grid_path: Path) -> dict:
    grid = load_json(grid_path)
    paths = source_paths(grid_path)
    diagnostics = build_diagnostics(grid)
    rows = [
        {
            "id": "nlrgit_01_budget_import",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "Import the cancellation-reduced finite-grid ratios and compute remaining ratio slack to the first-omitted threshold one.",
            "diagnostics": diagnostics,
            "source_artifacts": [paths["grid_note"], paths["asymptotic_target_note"]],
            "proof_boundary": "Exact budget bookkeeping only; not an interval certificate.",
        },
        {
            "id": "nlrgit_02_common_error_cap",
            "role": "exact_sufficient_condition",
            "readiness": "not_ready_to_apply",
            "claim": "If all intervalization and quadrature-certification errors are bounded by the proposed total ratio cap in both channels, the finite grid remains below one first omitted term.",
            "proof_boundary": "Sufficient finite-grid certification condition only; the cap is not yet proved.",
        },
        {
            "id": "nlrgit_03_interval_sources",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "claim": "The required interval sources are Laguerre node/weight balls, Phi and Phi' n-tail balls, coefficient-ratio propagation, quadrature-remainder error, and aggregation/rounding control.",
            "proof_boundary": "Open numerical-certification requirement only.",
        },
        {
            "id": "nlrgit_04_grid_to_collar_gap",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "claim": "A certified finite T grid still requires a separate continuum-in-T or far-tail theorem before it can become the full real-collar residual estimate.",
            "source_artifacts": [paths["uniform_remainder_target"]],
            "proof_boundary": "Open collar-uniformity requirement only.",
        },
        {
            "id": "nlrgit_05_floating_grid_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The existing floating cancellation-reduced grid already proves the first-omitted-term theorem.",
            "gap": "The grid lacks interval node/weight enclosures, quadrature-remainder bounds, Phi n-tail bounds, and T-continuum coverage.",
            "source_artifacts": [paths["grid_note"], paths["formal_tail_obstruction_note"]],
            "proof_boundary": "Rejected promotion only; not a proof of the actual remainder.",
        },
        {
            "id": "nlrgit_06_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any promoted remainder proof must pass the interval-source ledger and keep total ratio radii below the recorded cap before using the grid as certified evidence.",
            "source_artifacts": [paths["dependency_graph"]],
            "proof_boundary": "Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "certification_budget_rows": diagnostics["certification_budget_row_count"],
        "obligation_rows": diagnostics["obligation_row_count"],
        "open_requirement_rows": diagnostics["open_requirement_rows"],
        "ready_to_apply_rows": diagnostics["ready_to_apply_rows"],
        "observed_worst_value_ratio": diagnostics["source_grid_summary"]["max_value_ratio_to_first_omitted_over_orders"],
        "observed_worst_derivative_ratio": diagnostics["source_grid_summary"][
            "max_derivative_ratio_to_first_omitted_over_orders"
        ],
        "common_slack_to_one": diagnostics["budget_summary"]["common_slack_to_one"],
        "proposed_total_ratio_error_cap": diagnostics["budget_summary"]["proposed_total_ratio_error_cap"],
        "proposed_per_error_source_cap_for_five_sources": diagnostics["budget_summary"][
            "proposed_per_error_source_cap_for_five_sources"
        ],
        "closed_if_total_cap_met": diagnostics["budget_summary"]["closed_if_total_cap_met"],
        "target_closing": False,
        "main_finding": (
            "The cancellation-reduced finite grid leaves a common ratio slack of about 2.928994097967e-2 "
            "before reaching one first omitted term. A future interval certificate with total ratio error "
            "below 1.0e-2 in both value and derivative channels would keep the finite grid below the "
            "first-omitted threshold. The remaining work is to intervalize Laguerre nodes/weights, Phi tails, "
            "coefficient propagation, quadrature error, and then separately bridge the finite grid to the full collar."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target",
        "date": "2026-07-07",
        "status": "open numerical-certification target",
        "source_cancellation_reduced_grid_scout": paths["grid_note"],
        "source_actual_endpoint_scout": paths["actual_endpoint_note"],
        "source_asymptotic_remainder_target": paths["asymptotic_target_note"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction_note"],
        "source_residual_budget": paths["residual_budget_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py",
        "proof_boundary": (
            "Open numerical-certification target only. It gives a sufficient error budget for certifying the "
            "finite cancellation-reduced grid, but it does not provide interval enclosures, does not prove a "
            "uniform collar remainder theorem, does not prove scaled-curvature monotonicity, does not prove "
            "cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The floating grid is not promoted to a certified theorem.",
            "The intervalization target is finite-grid only until a continuum-in-T bridge is proved.",
            "The proposed error cap is sufficient target bookkeeping, not a proved bound.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][0]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian intervalization target: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['obligation_rows']} obligations, "
        f"{summary['open_requirement_rows']} open requirements, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Intervalization Target",
        "",
        "Date: 2026-07-07",
        "",
        "Status: open numerical-certification target. This is not a proof",
        "of a uniform residual estimate, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target`.",
        "",
        "Proof boundary: this artifact translates the cancellation-reduced",
        "floating grid into explicit interval-certification obligations and",
        "ratio-error budgets. It does not prove those obligations.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Error Budget",
        "",
        "```text",
        f"observed worst value ratio: {summary['observed_worst_value_ratio']}",
        f"observed worst derivative ratio: {summary['observed_worst_derivative_ratio']}",
        f"common slack to one: {summary['common_slack_to_one']}",
        f"proposed total ratio error cap: {summary['proposed_total_ratio_error_cap']}",
        f"proposed per-source cap for five sources: {summary['proposed_per_error_source_cap_for_five_sources']}",
        f"closed if total cap met: {summary['closed_if_total_cap_met']}",
        "```",
        "",
        "Per-channel budget rows:",
        "",
        "```text",
    ]
    for row in diagnostics["certification_budget_rows"]:
        lines.append(
            f"{row['channel']}: observed={row['observed_worst_ratio']} at {row['observed_worst_location']}, "
            f"slack={row['slack_to_one']}, half_slack={row['half_slack']}, "
            f"cap={row['proposed_total_ratio_error_cap']}, closes={row['closed_if_cap_met']}"
        )
    lines.extend(["```", "", "## Certification Obligations", "", "```text"])
    for row in diagnostics["obligation_rows"]:
        lines.append(
            f"{row['id']}: {row['role']} / {row['readiness']} / cap={row['target_error_cap']} / "
            f"{row['claim']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Open gap:",
            "",
            diagnostics["proof_boundary_note"],
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_cancellation_reduced_grid_scout"],
            artifact["source_actual_endpoint_scout"],
            artifact["source_asymptotic_remainder_target"],
            artifact["source_formal_tail_obstruction"],
            artifact["source_residual_budget"],
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
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(grid_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian intervalization target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
