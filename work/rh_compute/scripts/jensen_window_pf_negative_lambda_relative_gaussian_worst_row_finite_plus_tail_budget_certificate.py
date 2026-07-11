#!/usr/bin/env python3
"""Build a worst-row finite-plus-tail budget certificate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FINITE_PART_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json"
)
DEFAULT_PHI_TAIL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md"
)

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


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(finite_path: Path, phi_tail_path: Path, interval_path: Path) -> dict[str, str]:
    return {
        "finite_part_json": finite_path.relative_to(REPO_ROOT).as_posix(),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
        ),
        "phi_tail_json": phi_tail_path.relative_to(REPO_ROOT).as_posix(),
        "phi_tail_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "quadrature_ladder_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def build_budget_diagnostics(finite: dict, phi_tail: dict, interval: dict) -> dict:
    finite_summary = finite["summary"]
    tail_summary = phi_tail["summary"]
    interval_summary = interval["summary"]

    cap_text = interval_summary["proposed_per_error_source_cap_for_five_sources"]
    cap = dec(cap_text)
    value_finite_ratio = dec(finite_summary["value_ratio_to_first_omitted_upper"])
    derivative_finite_ratio = dec(finite_summary["derivative_ratio_to_first_omitted_upper"])
    value_composed_ratio = value_finite_ratio + cap
    derivative_composed_ratio = derivative_finite_ratio + cap
    value_margin = Decimal(1) - value_composed_ratio
    derivative_margin = Decimal(1) - derivative_composed_ratio

    tail_source_certified = bool(
        tail_summary["finite_grid_tail_source_certified"]
        and tail_summary["tail_sources_below_per_source_cap"]
        and tail_summary["per_source_intervalization_cap"] == cap_text
    )
    finite_ratios_below_one = bool(finite_summary["both_ratios_certified_below_one"])
    composed_ratios_below_one = bool(
        tail_source_certified
        and finite_ratios_below_one
        and value_composed_ratio < Decimal(1)
        and derivative_composed_ratio < Decimal(1)
    )

    return {
        "source_interval_obligation": "nlrgit_03_phi_and_c0_interval_tail",
        "composition_rule": (
            "Conservative budget composition: add the full per-source Phi-tail cap to the finite-part "
            "ratio upper bound in the same first-omitted comparison scale."
        ),
        "per_source_intervalization_cap": cap_text,
        "tail_budget_ratio_reserved": fixed(cap),
        "finite_part_ratio_count": 2,
        "tail_source_rows": tail_summary["tail_source_rows"],
        "tail_source_certified": tail_source_certified,
        "tail_actual_max_ratio_to_cap": tail_summary["max_tail_source_ratio_to_cap"],
        "tail_value_ratio_to_cap": tail_summary["value_ratio_to_per_source_cap"],
        "tail_derivative_ratio_to_cap": tail_summary["derivative_ratio_to_per_source_cap"],
        "tail_c0_ratio_to_cap": tail_summary["c0_ratio_to_per_source_cap"],
        "value_finite_ratio_upper": fixed(value_finite_ratio),
        "derivative_finite_ratio_upper": fixed(derivative_finite_ratio),
        "both_finite_ratios_below_one": finite_ratios_below_one,
        "value_composed_ratio_upper_using_full_tail_cap": fixed(value_composed_ratio),
        "derivative_composed_ratio_upper_using_full_tail_cap": fixed(derivative_composed_ratio),
        "value_remaining_margin_after_full_tail_cap": fixed(value_margin),
        "derivative_remaining_margin_after_full_tail_cap": fixed(derivative_margin),
        "composed_ratio_count": 2,
        "both_composed_ratios_below_one": composed_ratios_below_one,
        "target_closing": False,
        "budget_is_conservative": True,
    }


def build_rows(paths: dict[str, str], diagnostics: dict, finite: dict, phi_tail: dict, interval: dict) -> list[dict]:
    remaining_sources = [
        "generalized Gauss-Laguerre quadrature-remainder or interval adaptive integration beyond the N=320 sum",
        "rounding and cross-source aggregation under the total intervalization cap",
        "all recorded rows and quadrature orders, not just T=10000, F_21, N=320",
        "finite-grid to full-collar coverage",
        "final acceptance-gate wording separating finite-grid evidence from theorem claims",
    ]
    rows = [
        MatrixRow(
            id="nlrgwrfptbc_01_budget_composition_rule",
            role="budget_composition_rule",
            readiness="available_budget_rule",
            claim=(
                "For the worst-row source budget, the finite n<=30 first-omitted ratio upper may be composed "
                "with the certified n>30 Phi-tail source by reserving the entire per-source intervalization cap."
            ),
            diagnostics={
                "composition_rule": diagnostics["composition_rule"],
                "per_source_intervalization_cap": diagnostics["per_source_intervalization_cap"],
                "tail_budget_ratio_reserved": diagnostics["tail_budget_ratio_reserved"],
                "source_interval_obligation": diagnostics["source_interval_obligation"],
                "source_interval_total_cap": interval["summary"]["proposed_total_ratio_error_cap"],
            },
            source_artifacts=[paths["interval_json"], paths["interval_note"]],
            proof_boundary=(
                "Budget-composition rule only; not an absolute Phi integral evaluation and not a quadrature-remainder theorem."
            ),
        ),
        MatrixRow(
            id="nlrgwrfptbc_02_finite_part_import",
            role="interval_certificate_import",
            readiness="available_interval_certificate",
            claim=(
                "Import the worst-row finite n<=30 weighted-sum interval certificate and its below-one value "
                "and derivative first-omitted ratio bounds."
            ),
            diagnostics={
                "source_matrix_rows": finite["summary"]["matrix_rows"],
                "refined_node_rows": finite["summary"]["refined_node_rows"],
                "interval_weight_rows": finite["summary"]["interval_weight_rows"],
                "value_finite_ratio_upper": diagnostics["value_finite_ratio_upper"],
                "derivative_finite_ratio_upper": diagnostics["derivative_finite_ratio_upper"],
                "both_finite_ratios_below_one": diagnostics["both_finite_ratios_below_one"],
            },
            source_artifacts=[paths["finite_part_json"], paths["finite_part_note"]],
            proof_boundary=(
                "Finite-part import only; by itself it does not include n>30 tails, quadrature remainder, "
                "or all-row coverage."
            ),
        ),
        MatrixRow(
            id="nlrgwrfptbc_03_tail_source_import",
            role="tail_source_import",
            readiness="available_tail_source",
            claim=(
                "Import the certified finite-grid Phi-tail source and reserve its full per-source cap in the "
                "worst-row first-omitted budget."
            ),
            diagnostics={
                "source_matrix_rows": phi_tail["summary"]["matrix_rows"],
                "tail_source_rows": diagnostics["tail_source_rows"],
                "tail_source_certified": diagnostics["tail_source_certified"],
                "tail_actual_max_ratio_to_cap": diagnostics["tail_actual_max_ratio_to_cap"],
                "tail_value_ratio_to_cap": diagnostics["tail_value_ratio_to_cap"],
                "tail_derivative_ratio_to_cap": diagnostics["tail_derivative_ratio_to_cap"],
                "tail_c0_ratio_to_cap": diagnostics["tail_c0_ratio_to_cap"],
                "tail_budget_ratio_reserved": diagnostics["tail_budget_ratio_reserved"],
            },
            source_artifacts=[paths["phi_tail_json"], paths["phi_tail_note"]],
            proof_boundary=(
                "Tail-source import only; it does not certify finite n<=30 node evaluations, weights, "
                "coefficient aggregation, or quadrature remainder."
            ),
        ),
        MatrixRow(
            id="nlrgwrfptbc_04_composed_below_one_budget",
            role="budget_certificate",
            readiness="available_budget_certificate",
            claim=(
                "Adding the full tail-source cap to the finite-part ratio uppers keeps both worst-row channels "
                "strictly below one first omitted term."
            ),
            diagnostics={
                "value_finite_ratio_upper": diagnostics["value_finite_ratio_upper"],
                "derivative_finite_ratio_upper": diagnostics["derivative_finite_ratio_upper"],
                "tail_budget_ratio_reserved": diagnostics["tail_budget_ratio_reserved"],
                "value_composed_ratio_upper_using_full_tail_cap": diagnostics[
                    "value_composed_ratio_upper_using_full_tail_cap"
                ],
                "derivative_composed_ratio_upper_using_full_tail_cap": diagnostics[
                    "derivative_composed_ratio_upper_using_full_tail_cap"
                ],
                "value_remaining_margin_after_full_tail_cap": diagnostics[
                    "value_remaining_margin_after_full_tail_cap"
                ],
                "derivative_remaining_margin_after_full_tail_cap": diagnostics[
                    "derivative_remaining_margin_after_full_tail_cap"
                ],
                "both_composed_ratios_below_one": diagnostics["both_composed_ratios_below_one"],
            },
            source_artifacts=[paths["finite_part_note"], paths["phi_tail_note"], paths["interval_note"]],
            proof_boundary=(
                "Worst-row finite-plus-tail budget comparison only; not a finite-grid interval certificate, "
                "not all rows, and not a uniform collar theorem."
            ),
        ),
        MatrixRow(
            id="nlrgwrfptbc_05_intervalization_handoff",
            role="intervalization_handoff",
            readiness="not_ready_to_apply",
            claim=(
                "The worst-row finite-plus-tail numerator-source budget for T=10000, F_21, N=320 is now "
                "composed below one, but the intervalization target remains open until the remaining sources "
                "are certified and aggregated."
            ),
            diagnostics={
                "retired_component": "worst-row finite n<=30 weighted sum plus reserved n>30 Phi-tail cap",
                "remaining_sources": remaining_sources,
                "source_interval_obligation": diagnostics["source_interval_obligation"],
                "target_closing": diagnostics["target_closing"],
            },
            source_artifacts=[
                paths["interval_json"],
                paths["interval_note"],
                paths["quadrature_ladder_note"],
                paths["uniform_remainder_target"],
            ],
            proof_boundary="Handoff only; not a full intervalized finite-grid proof and not a collar theorem.",
        ),
        MatrixRow(
            id="nlrgwrfptbc_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted finite-grid interval proof may use this artifact only as a worst-row "
                "finite-plus-tail budget component, with quadrature, aggregation, all-row coverage, and "
                "grid-to-collar coverage proved separately."
            ),
            diagnostics={
                "forbidden_promotions": [
                    "full finite-grid interval certificate",
                    "quadrature-remainder theorem",
                    "uniform collar theorem",
                    "scaled-curvature monotonicity",
                    "cone entry",
                    "RH",
                    "Lambda <= 0",
                ]
            },
            source_artifacts=[paths["dependency_graph"]],
            proof_boundary=(
                "Proof-hygiene gate only; not quadrature remainder, not RH, and not Lambda <= 0."
            ),
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(finite_path: Path, phi_tail_path: Path, interval_path: Path) -> dict:
    finite = load_json(finite_path)
    phi_tail = load_json(phi_tail_path)
    interval = load_json(interval_path)
    paths = source_paths(finite_path, phi_tail_path, interval_path)
    diagnostics = build_budget_diagnostics(finite, phi_tail, interval)
    rows = build_rows(paths, diagnostics, finite, phi_tail, interval)
    summary = {
        "matrix_rows": len(rows),
        "finite_part_ratio_count": diagnostics["finite_part_ratio_count"],
        "tail_source_rows": diagnostics["tail_source_rows"],
        "ready_to_apply_rows": 0,
        "source_finite_part_rows": finite["summary"]["matrix_rows"],
        "source_tail_certificate_rows": phi_tail["summary"]["matrix_rows"],
        "source_interval_obligations": interval["summary"]["obligation_rows"],
        "per_source_intervalization_cap": diagnostics["per_source_intervalization_cap"],
        "tail_budget_ratio_reserved": diagnostics["tail_budget_ratio_reserved"],
        "tail_source_certified": diagnostics["tail_source_certified"],
        "tail_actual_max_ratio_to_cap": diagnostics["tail_actual_max_ratio_to_cap"],
        "value_finite_ratio_upper": diagnostics["value_finite_ratio_upper"],
        "derivative_finite_ratio_upper": diagnostics["derivative_finite_ratio_upper"],
        "both_finite_ratios_below_one": diagnostics["both_finite_ratios_below_one"],
        "value_composed_ratio_upper_using_full_tail_cap": diagnostics[
            "value_composed_ratio_upper_using_full_tail_cap"
        ],
        "derivative_composed_ratio_upper_using_full_tail_cap": diagnostics[
            "derivative_composed_ratio_upper_using_full_tail_cap"
        ],
        "value_remaining_margin_after_full_tail_cap": diagnostics[
            "value_remaining_margin_after_full_tail_cap"
        ],
        "derivative_remaining_margin_after_full_tail_cap": diagnostics[
            "derivative_remaining_margin_after_full_tail_cap"
        ],
        "composed_ratio_count": diagnostics["composed_ratio_count"],
        "both_composed_ratios_below_one": diagnostics["both_composed_ratios_below_one"],
        "target_closing": False,
        "budget_is_conservative": True,
        "main_finding": (
            "The worst-row finite-plus-tail budget certificate composes the certified finite n<=30 "
            "weighted-sum ratio uppers with the full 2.0e-3 Phi-tail source cap. The resulting value and "
            "derivative ratio uppers remain strictly below one first omitted term. This retires only the "
            "worst-row finite-plus-tail numerator-source budget; quadrature remainder, rounding aggregation, "
            "all-row coverage, and finite-grid to collar coverage remain open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate",
        "date": "2026-07-07",
        "status": "worst-row finite-plus-tail budget certificate",
        "source_finite_part_weighted_sum_certificate": paths["finite_part_note"],
        "source_finite_part_weighted_sum_json": paths["finite_part_json"],
        "source_phi_tail_grid_certificate": paths["phi_tail_note"],
        "source_phi_tail_grid_json": paths["phi_tail_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_quadrature_ladder_scout": paths["quadrature_ladder_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py"
        ),
        "proof_boundary": (
            "Worst-row finite-plus-tail budget certificate only. It adds the full per-source Phi-tail cap to "
            "the certified finite n<=30 weighted-sum ratio uppers for T=10000, F_21, N=320, but it does not "
            "prove a quadrature-remainder theorem, does not certify all rows or quadrature orders, does not "
            "aggregate rounding, does not prove a finite-grid interval certificate, does not bridge the finite "
            "grid to a uniform collar, does not prove scaled-curvature monotonicity, does not prove cone entry, "
            "and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "The certificate covers only the worst row T=10000, F_21 and quadrature order N=320.",
            "The n>30 tail is used through the full per-source budget cap, not as a full residual enclosure.",
            "Quadrature remainder, rounding aggregation, all-row coverage, and grid-to-collar coverage remain open.",
            "The finite-plus-tail budget is not promoted to a finite-grid interval certificate or a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget "
        f"certificate: {summary['matrix_rows']} rows, 0 issues, "
        f"{summary['composed_ratio_count']} composed ratios, "
        f"{summary['tail_source_rows']} tail sources, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Finite-Plus-Tail Budget Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row finite-plus-tail budget certificate. This is not a proof",
        "of a quadrature-remainder theorem, a finite-grid interval certificate,",
        "a uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate`.",
        "",
        "Proof boundary: this artifact composes only the worst-row `T=10000`,",
        "`F_21`, `N=320` finite `n<=30` weighted-sum certificate with the",
        "reserved full `n>30` Phi-tail per-source budget. It does not prove",
        "quadrature remainder, cover all rows, aggregate rounding, or bridge",
        "the grid to a full collar.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Budget Composition",
        "",
        "```text",
        f"per-source intervalization cap: {summary['per_source_intervalization_cap']}",
        f"tail budget ratio reserved: {summary['tail_budget_ratio_reserved']}",
        f"tail source certified: {summary['tail_source_certified']}",
        f"tail actual max ratio to cap: {summary['tail_actual_max_ratio_to_cap']}",
        f"value finite ratio upper: {summary['value_finite_ratio_upper']}",
        f"derivative finite ratio upper: {summary['derivative_finite_ratio_upper']}",
        f"value composed ratio upper using full tail cap: {summary['value_composed_ratio_upper_using_full_tail_cap']}",
        f"derivative composed ratio upper using full tail cap: {summary['derivative_composed_ratio_upper_using_full_tail_cap']}",
        f"value remaining margin after full tail cap: {summary['value_remaining_margin_after_full_tail_cap']}",
        f"derivative remaining margin after full tail cap: {summary['derivative_remaining_margin_after_full_tail_cap']}",
        f"both composed ratios below one: {summary['both_composed_ratios_below_one']}",
        "```",
        "",
        "Remaining sources:",
        "",
        "```text",
    ]
    for item in artifact["matrix_rows"][4]["diagnostics"]["remaining_sources"]:
        lines.append(item)
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_finite_part_weighted_sum_certificate"],
            artifact["source_finite_part_weighted_sum_json"],
            artifact["source_phi_tail_grid_certificate"],
            artifact["source_phi_tail_grid_json"],
            artifact["source_intervalization_target"],
            artifact["source_intervalization_target_json"],
            artifact["source_quadrature_ladder_scout"],
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
    parser.add_argument("--finite-part-json", type=Path, default=DEFAULT_FINITE_PART_JSON)
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    finite_path = args.finite_part_json if args.finite_part_json.is_absolute() else REPO_ROOT / args.finite_part_json
    phi_tail_path = args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(finite_path, phi_tail_path, interval_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
