#!/usr/bin/env python3
"""Build a finite-grid Phi-tail source certificate for the relative-Gaussian intervalization route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PHI_TAIL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json"
)
DEFAULT_NODE_C0_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md"

getcontext().prec = 120


@dataclass(frozen=True)
class TailSourceRow:
    channel: str
    tail_bound: str
    normalized_or_relative_bound: str
    per_source_cap: str
    ratio_to_per_source_cap: str
    certified_below_cap: bool
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


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def sci(value: Decimal) -> str:
    return f"{value:.18E}"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(phi_tail_path: Path, node_c0_path: Path, interval_path: Path) -> dict[str, str]:
    return {
        "phi_tail_json": phi_tail_path.relative_to(REPO_ROOT).as_posix(),
        "phi_tail_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md",
        "node_c0_json": node_c0_path.relative_to(REPO_ROOT).as_posix(),
        "node_c0_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
        "grid_json": (
            "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def build_tail_source_rows(phi_tail: dict, interval: dict) -> list[TailSourceRow]:
    summary = phi_tail["summary"]
    cap = dec(interval["summary"]["proposed_per_error_source_cap_for_five_sources"])
    channels = [
        (
            "Phi value n>30 tail",
            summary["value_phi_tail_bound"],
            summary["normalized_value_tail_bound_using_c0_proxy"],
        ),
        (
            "Phi prime n>30 derivative-core tail",
            summary["derivative_phi_prime_tail_bound"],
            summary["normalized_derivative_core_tail_bound_using_c0_proxy"],
        ),
        (
            "Phi(0) n>30 denominator tail",
            summary["c0_tail_bound"],
            summary["denominator_relative_tail_bound_using_c0_proxy"],
        ),
    ]
    rows = []
    for channel, tail_bound, normalized in channels:
        ratio = dec(normalized) / cap
        rows.append(
            TailSourceRow(
                channel=channel,
                tail_bound=tail_bound,
                normalized_or_relative_bound=normalized,
                per_source_cap=interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
                ratio_to_per_source_cap=sci(ratio),
                certified_below_cap=bool(ratio < 1),
                proof_boundary=(
                    "Finite-grid Phi-tail source bound only; not a Laguerre node/weight interval, "
                    "not a quadrature-remainder bound, and not a full residual numerator enclosure."
                ),
            )
        )
    return rows


def build_side_condition_diagnostics(phi_tail: dict, node_c0: dict) -> dict:
    node_summary = node_c0["summary"]
    phi_summary = phi_tail["summary"]
    return {
        "padded_tail_range": "0<=x<=1",
        "phi_tail_padded_x_upper": phi_summary["padded_x_range_upper"],
        "certified_worst_x_square_upper_bound": node_summary["worst_x_square_upper_bound"],
        "certified_worst_x_square_upper_bound_decimal": node_summary["worst_x_square_upper_bound_decimal"],
        "node_range_x_le_1_certified": node_summary["node_range_x_le_1_certified"],
        "certified_c0_lower": node_summary["certified_c0_lower"],
        "n1_phi0_term_lower": node_summary["n1_phi0_term_lower"],
        "c0_lower_certified_by_n1_term": node_summary["c0_lower_certified_by_n1_term"],
        "side_conditions_certified": bool(
            node_summary["node_range_x_le_1_certified"] and node_summary["c0_lower_certified_by_n1_term"]
        ),
        "side_condition_source": (
            "Gershgorin/AM-GM certifies every recorded Laguerre node has x<=1, and Arb certifies "
            "Phi(0)>=0.44 from the n=1 term."
        ),
        "proof_boundary": (
            "Side-condition certificate only; individual roots, weights, and quadrature remainders remain open."
        ),
    }


def build_rows(paths: dict[str, str], side_conditions: dict, tail_rows: list[TailSourceRow], summary_bits: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgptgc_01_conditional_tail_scout_import",
            role="analytic_tail_bound",
            readiness="available_conditional",
            claim=(
                "Import the padded-range n>30 majorants for Phi, Phi', and Phi(0) from the tail-bound scout."
            ),
            diagnostics={
                "tail_source_rows": [asdict(row) for row in tail_rows],
                "tail_bounds_below_per_source_cap": summary_bits["tail_sources_below_per_source_cap"],
            },
            source_artifacts=[paths["phi_tail_json"], paths["phi_tail_note"]],
            proof_boundary="Analytic tail-bound import only; not a node, weight, or quadrature certificate.",
        ),
        MatrixRow(
            id="nlrgptgc_02_side_conditions_certified",
            role="finite_side_condition_certificate",
            readiness="available_for_phi_tail_source",
            claim=(
                "The recorded finite grid satisfies the side conditions required by the padded tail scout: "
                "node-induced x<=1 and Phi(0)>=0.44."
            ),
            diagnostics=side_conditions,
            source_artifacts=[paths["node_c0_json"], paths["node_c0_note"]],
            proof_boundary=(
                "Finite-grid side-condition certificate only; not individual node/weight intervals or quadrature error."
            ),
        ),
        MatrixRow(
            id="nlrgptgc_03_finite_grid_tail_source_certificate",
            role="finite_tail_source_certificate",
            readiness="available_for_intervalization",
            claim=(
                "Combining the certified side conditions with the padded majorants certifies the Phi/Phi'/Phi(0) "
                "n-tail source below the per-source intervalization cap for the recorded finite grid."
            ),
            diagnostics=summary_bits,
            source_artifacts=[paths["phi_tail_note"], paths["node_c0_note"], paths["interval_note"]],
            proof_boundary=(
                "Finite-grid tail-source certificate only; not a full Phi/Phi' node evaluation enclosure, "
                "not coefficient propagation, and not quadrature error."
            ),
        ),
        MatrixRow(
            id="nlrgptgc_04_intervalization_obligation_handoff",
            role="intervalization_handoff",
            readiness="available_for_intervalization",
            claim=(
                "The n-tail part of nlrgit_03 may now be treated as certified for the recorded finite grid, "
                "while the non-tail finite-part and quadrature obligations remain open."
            ),
            diagnostics={
                "source_obligation_id": "nlrgit_03_phi_and_c0_interval_tail",
                "certified_component": "n>30 tails in Phi, Phi', and Phi(0) on the recorded finite grid",
                "remaining_obligations": summary_bits["remaining_obligations"],
            },
            source_artifacts=[paths["interval_json"], paths["interval_note"], paths["dependency_graph"]],
            proof_boundary="Handoff only; this does not close the intervalization target or uniform collar theorem.",
        ),
        MatrixRow(
            id="nlrgptgc_05_full_phi_evaluation_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The finite-grid Phi-tail source certificate supplies full interval values for Phi and Phi' at every quadrature node.",
            gap=(
                "It bounds only the omitted n>30 source. Full node evaluations still need individual node "
                "intervals, finite n<=30 interval evaluation, coefficient propagation, and quadrature-radius control."
            ),
            source_artifacts=[paths["interval_note"], paths["uniform_remainder_target"]],
            proof_boundary="Rejected promotion only; not a full residual numerator certificate, RH, or Lambda <= 0.",
        ),
        MatrixRow(
            id="nlrgptgc_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted finite-grid interval proof may use this tail-source certificate only after separately "
                "certifying the finite n<=30 evaluations, coefficient balls, quadrature remainder, and aggregation."
            ),
            source_artifacts=[paths["dependency_graph"]],
            proof_boundary="Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(phi_tail_path: Path, node_c0_path: Path, interval_path: Path) -> dict:
    phi_tail = load_json(phi_tail_path)
    node_c0 = load_json(node_c0_path)
    interval = load_json(interval_path)
    paths = source_paths(phi_tail_path, node_c0_path, interval_path)
    tail_rows = build_tail_source_rows(phi_tail, interval)
    side_conditions = build_side_condition_diagnostics(phi_tail, node_c0)
    remaining_obligations = [
        "individual Laguerre node and weight intervals beyond the coarse x<=1 range",
        "finite n<=30 interval evaluation of Phi and Phi' at certified node intervals",
        "coefficient-ratio ball propagation for the cancellation-reduced polynomial core",
        "generalized Gauss-Laguerre quadrature-remainder or interval adaptive integration",
        "rounding aggregation and finite-grid to full-collar coverage",
    ]
    tail_sources_below_cap = all(row.certified_below_cap for row in tail_rows)
    finite_grid_tail_source_certified = bool(side_conditions["side_conditions_certified"] and tail_sources_below_cap)
    summary_bits = {
        "tail_source_rows": len(tail_rows),
        "side_conditions_certified": side_conditions["side_conditions_certified"],
        "tail_sources_below_per_source_cap": tail_sources_below_cap,
        "finite_grid_tail_source_certified": finite_grid_tail_source_certified,
        "per_source_intervalization_cap": interval["summary"]["proposed_per_error_source_cap_for_five_sources"],
        "value_ratio_to_per_source_cap": tail_rows[0].ratio_to_per_source_cap,
        "derivative_ratio_to_per_source_cap": tail_rows[1].ratio_to_per_source_cap,
        "c0_ratio_to_per_source_cap": tail_rows[2].ratio_to_per_source_cap,
        "remaining_obligations": remaining_obligations,
        "proof_boundary": (
            "Certified finite-grid Phi-tail source only; it does not certify the full cancellation-reduced "
            "numerator or the quadrature integral."
        ),
    }
    rows = build_rows(paths, side_conditions, tail_rows, summary_bits)
    summary = {
        "matrix_rows": len(rows),
        "tail_source_rows": len(tail_rows),
        "certified_side_conditions": 2,
        "ready_to_apply_rows": 0,
        "source_phi_tail_rows": phi_tail["summary"]["matrix_rows"],
        "source_node_c0_rows": node_c0["summary"]["matrix_rows"],
        "source_interval_obligations": interval["summary"]["obligation_rows"],
        "padded_tail_range": side_conditions["padded_tail_range"],
        "node_range_x_le_1_certified": side_conditions["node_range_x_le_1_certified"],
        "certified_c0_lower": side_conditions["certified_c0_lower"],
        "c0_lower_certified_by_n1_term": side_conditions["c0_lower_certified_by_n1_term"],
        "per_source_intervalization_cap": summary_bits["per_source_intervalization_cap"],
        "max_tail_source_ratio_to_cap": sci(max(dec(row.ratio_to_per_source_cap) for row in tail_rows)),
        "value_ratio_to_per_source_cap": summary_bits["value_ratio_to_per_source_cap"],
        "derivative_ratio_to_per_source_cap": summary_bits["derivative_ratio_to_per_source_cap"],
        "c0_ratio_to_per_source_cap": summary_bits["c0_ratio_to_per_source_cap"],
        "tail_sources_below_per_source_cap": tail_sources_below_cap,
        "finite_grid_tail_source_certified": finite_grid_tail_source_certified,
        "remaining_obligation_count": len(remaining_obligations),
        "target_closing": False,
        "main_finding": (
            "The Phi-tail source for the recorded finite grid is now certificate-composed rather than merely "
            "conditional: the padded n>30 majorants are far below the 2.0e-3 per-source cap, while the node-c0 "
            "certificate proves x<=1 and Phi(0)>=0.44 for the same grid. This certifies the omitted-n tail "
            "component only; finite n<=30 node evaluation, weights, quadrature error, coefficient propagation, "
            "rounding, and grid-to-collar coverage remain open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate",
        "date": "2026-07-07",
        "status": "finite-grid Phi-tail source certificate",
        "source_phi_tail_bound_scout": paths["phi_tail_note"],
        "source_phi_tail_bound_json": paths["phi_tail_json"],
        "source_node_c0_range_certificate": paths["node_c0_note"],
        "source_node_c0_range_json": paths["node_c0_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_cancellation_reduced_grid_scout": paths["grid_note"],
        "source_cancellation_reduced_grid_json": paths["grid_json"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py",
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py"
        ),
        "proof_boundary": (
            "Finite-grid Phi-tail source certificate only. It composes the padded n-tail majorants with the "
            "certified x<=1 and Phi(0)>=0.44 side conditions, but it does not provide individual Laguerre "
            "node or weight intervals, does not certify finite n<=30 Phi/Phi' node evaluations, does not "
            "bound quadrature remainder, does not prove a uniform collar theorem, does not prove scaled-curvature "
            "monotonicity, and does not prove RH or Lambda <= 0."
        ),
        "tail_source_rows": [asdict(row) for row in tail_rows],
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "Only the omitted n>30 Phi-tail source is certified for the recorded finite grid.",
            "Finite n<=30 node evaluation, Laguerre weights, quadrature, coefficient propagation, rounding, and grid-to-collar coverage remain open.",
            "The finite grid is not promoted to a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['tail_source_rows']} certified tail sources, "
        f"{summary['certified_side_conditions']} certified side conditions, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Phi-Tail Grid Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite-grid Phi-tail source certificate. This is not a proof",
        "of a finite-grid interval certificate, quadrature-remainder theorem,",
        "uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate`.",
        "",
        "Proof boundary: this artifact certifies only the omitted `n>30`",
        "Phi-tail source for the recorded finite grid by combining the padded",
        "tail majorants with certified `x<=1` and `Phi(0)>=0.44` side",
        "conditions. It does not certify finite `n<=30` node evaluations.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Certified Inputs",
        "",
        "```text",
        f"padded tail range: {summary['padded_tail_range']}",
        f"node range x<=1 certified: {summary['node_range_x_le_1_certified']}",
        f"certified c0 lower: {summary['certified_c0_lower']}",
        f"c0 lower certified by n=1 term: {summary['c0_lower_certified_by_n1_term']}",
        f"per-source intervalization cap: {summary['per_source_intervalization_cap']}",
        f"finite-grid tail source certified: {summary['finite_grid_tail_source_certified']}",
        "```",
        "",
        "Tail source ratios to cap:",
        "",
        "```text",
    ]
    for row in artifact["tail_source_rows"]:
        lines.append(
            f"{row['channel']}: normalized/relative={row['normalized_or_relative_bound']}, "
            f"ratio_to_cap={row['ratio_to_per_source_cap']}, certified={row['certified_below_cap']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Remaining obligations:",
            "",
            "```text",
        ]
    )
    for item in artifact["matrix_rows"][3]["diagnostics"]["remaining_obligations"]:
        lines.append(item)
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_phi_tail_bound_scout"],
            artifact["source_phi_tail_bound_json"],
            artifact["source_node_c0_range_certificate"],
            artifact["source_node_c0_range_json"],
            artifact["source_intervalization_target"],
            artifact["source_intervalization_target_json"],
            artifact["source_cancellation_reduced_grid_scout"],
            artifact["source_cancellation_reduced_grid_json"],
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
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--node-c0-json", type=Path, default=DEFAULT_NODE_C0_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    phi_tail_path = args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    node_c0_path = args.node_c0_json if args.node_c0_json.is_absolute() else REPO_ROOT / args.node_c0_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(phi_tail_path, node_c0_path, interval_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
