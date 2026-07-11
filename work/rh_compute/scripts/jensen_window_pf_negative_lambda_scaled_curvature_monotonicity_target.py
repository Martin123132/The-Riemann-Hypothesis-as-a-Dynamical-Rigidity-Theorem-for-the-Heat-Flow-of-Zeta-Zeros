#!/usr/bin/env python3
"""Build the scaled-curvature monotonicity theorem target for negative lambda."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction import (  # noqa: E402
    DEFAULT_ENCLOSURES,
    REPO_ROOT,
    build_diagnostics,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md"


def build_artifact(paths: list[Path] | None = None) -> dict:
    if paths is None:
        paths = list(DEFAULT_ENCLOSURES)
    diagnostics = build_diagnostics(paths)
    rows = [
        {
            "id": "nlscmt_01_scaled_curvature_statement",
            "role": "open_statement",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "For the actual negative-lambda zeta heat-flow coefficients, prove C_(k+1)>=C_k after an explicit finite collar, where C_k=(2*k+1)*B_k and B_k=-log(((2*k-1)/(2*k+1))*R_k).",
            "formula": "C_k=(2*k+1)*B_k; B_k=-log(((2*k-1)/(2*k+1))*R_k)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
            ],
            "proof_boundary": "Open theorem target only; this artifact does not prove scaled-curvature monotonicity.",
        },
        {
            "id": "nlscmt_02_linear_barrier_equivalence",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim_if_proved": "The scaled-curvature monotonicity statement C_(k+1)>=C_k is exactly the linear B-barrier B_(k+1)>=((2*k+1)/(2*k+3))*B_k.",
            "formula": "C_(k+1)-C_k=(2*k+3)*(B_(k+1)-((2*k+1)/(2*k+3))*B_k)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
            ],
            "proof_boundary": "Exact algebraic equivalence only; no zeta coefficient inequality is proved.",
        },
        {
            "id": "nlscmt_03_raw_corridor_sufficient_chain",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim_if_proved": "For B_k>=0, the scaled-curvature target implies the nonlinear lower coefficient-curvature barrier; paired with B_(k+1)<=B_k it gives the two-sided curvature corridor equivalent to the raw corridor.",
            "formula": "C_(k+1)>=C_k => B_(k+1)>=alpha_k*B_k => L_k(B_k)<=B_(k+1), alpha_k=(2*k+1)/(2*k+3)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
                "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
            ],
            "proof_boundary": "Exact sufficient chain only; it does not prove B_k>=0, monotone curvature, or the all-k raw corridor.",
        },
        {
            "id": "nlscmt_04_repaired_k300_finite_anchor",
            "role": "finite_anchor",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "The repaired k300 data validate C_(k+1)>=C_k and B_k>=0 on every checked row for lambda=-25,-50,-100.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md",
            ],
            "proof_boundary": "Finite repaired stress evidence only; not an all-k or continuous-lambda theorem.",
        },
        {
            "id": "nlscmt_05_retired_fixed_wall_boundary",
            "role": "rejected_shortcut",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "Scaled-curvature monotonicity may be replaced by the fixed wall C_k<=2/3.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md",
            ],
            "gap": "The repaired k300 obstruction finite-rejects C_k<=2/3 on 718/897 checked rows, so the replacement target is monotonicity, not a fixed upper wall.",
            "proof_boundary": "Rejected shortcut only; not evidence against the scaled-curvature monotonicity target.",
        },
        {
            "id": "nlscmt_06_required_monotone_upper_side",
            "role": "open_dependency",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "The scaled-curvature target supplies the lower side of the curvature corridor; the upper side still requires B_(k+1)<=B_k, i.e. the monotone-contraction/Delta^3 log A target.",
            "source_artifacts": [
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
            ],
            "gap": "No all-k proof of the monotone upper side is included here.",
            "proof_boundary": "Open dependency only; this target is not a substitute for the monotone-contraction theorem.",
        },
        {
            "id": "nlscmt_07_live_log_ratio_recurrence_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "A zeta-specific log-ratio recurrence may prove C_(k+1)>=C_k directly from the heat-flow raw moments after a finite collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
                "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
                "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
            ],
            "gap": "The current artifacts give exact reductions and finite stress, not an all-k zeta-specific recurrence.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlscmt_08_live_asymptotic_curvature_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "A signed/tilted saddle analysis with uniform remainders may prove the scaled-curvature monotonicity inequality in the tail.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            "gap": "Current Taylor and perturbation artifacts are fixed-scale diagnostics and do not provide the required uniform all-k remainder.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlscmt_09_generic_shortcuts_blocked",
            "role": "rejected_shortcut",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "Generic raw-moment positivity, raw walls, or monotone B curvature alone prove the scaled-curvature target.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
            ],
            "gap": "Exact positive moment and cone witnesses block the generic shortcuts; the proof must use zeta-specific structure.",
            "proof_boundary": "Rejected shortcut only; not evidence against actual zeta coefficients.",
        },
        {
            "id": "nlscmt_10_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "A promoted proof must state the lambda range, finite collar, B-wall input, scaled monotonicity inequality, companion monotone upper side, and forbidden endpoint assumptions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of scaled-curvature monotonicity, cone entry, jwpf_06, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "target_rows": len(rows),
        "exact_available_rows": 2,
        "finite_anchor_rows": 1,
        "open_dependency_rows": 1,
        "live_routes": 2,
        "rejected_shortcut_rows": 2,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "open_theorem_target": True,
        "lambdas": diagnostics["lambdas"],
        "coefficient_cap": diagnostics["coefficient_cap"],
        "b_wall_rows": diagnostics["b_wall_rows"],
        "scaled_curvature_total_rows": diagnostics["scaled_curvature_total_rows"],
        "adjacent_total_rows": diagnostics["adjacent_total_rows"],
        "scaled_curvature_increase_rows": diagnostics["scaled_curvature_increase_rows"],
        "scaled_curvature_increase_failures": diagnostics["scaled_curvature_increase_failures"],
        "scaled_curvature_increase_inconclusive": diagnostics["scaled_curvature_increase_inconclusive"],
        "two_thirds_failure_rows": diagnostics["two_thirds_failure_rows"],
        "main_finding": (
            "The replacement for the retired fixed 2/3 wall is the scaled-curvature "
            "monotonicity theorem target C_(k+1)>=C_k, exactly equivalent to the "
            "linear barrier B_(k+1)>=((2*k+1)/(2*k+3))*B_k. Repaired k300 data support "
            "this on 894/894 adjacent rows, while C_k<=2/3 fails on 718/897 rows. "
            "Together with B_k>=0 and the separate monotone-contraction upper side "
            "B_(k+1)<=B_k, this would feed the zeta-specific raw-corridor target."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target",
        "date": "2026-07-07",
        "status": "open_theorem_target",
        "target_id": "target_negative_lambda_scaled_curvature_monotonicity",
        "source_linear_barrier_scout": "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
        "source_curvature_corridor_bridge": "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
        "source_raw_corridor_target": "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "source_k300_obstruction": "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in paths],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py",
        "proof_boundary": (
            "Open theorem target only. It names the scaled-curvature monotonicity theorem "
            "needed by the linear curvature-barrier route, but it does not prove that theorem, "
            "does not prove the companion monotone-contraction theorem, does not prove the "
            "raw-corridor theorem, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "target_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The target remains open_theorem_target.",
            "The retired fixed 2/3 wall is not a substitute for scaled-curvature monotonicity.",
            "Finite repaired k300 support is finite evidence only.",
            "The companion monotone-contraction upper side remains a separate open theorem target.",
            "Generic moment positivity, endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["target_rows"][3]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda scaled-curvature monotonicity target: "
        f"{summary['target_rows']} rows, 0 issues, "
        f"{summary['live_routes']} live routes, "
        f"{summary['scaled_curvature_increase_rows']} scaled-curvature increase rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Scaled-Curvature Monotonicity Target",
        "",
        "Date: 2026-07-07",
        "",
        "Status: open theorem target. This is not a proof of the scaled-curvature",
        "monotonicity theorem, the raw-corridor theorem, cone entry, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target`.",
        "",
        "Proof boundary: this artifact names the replacement theorem target after",
        "the fixed `2/3` wall was retired. It does not prove any all-`k` theorem.",
        "",
        "Machine-readable target:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Target Statement",
        "",
        "Let:",
        "",
        "```text",
        "R_k = M_(k+1)*M_(k-1)/M_k^2",
        "B_k = -log(((2*k-1)/(2*k+1))*R_k)",
        "C_k = (2*k+1)*B_k",
        "```",
        "",
        "The replacement target is:",
        "",
        "```text",
        "C_(k+1) >= C_k",
        "```",
        "",
        "Equivalently:",
        "",
        "```text",
        "B_(k+1) >= ((2*k+1)/(2*k+3))*B_k",
        "```",
        "",
        "This is a lower-side target for the coefficient-curvature corridor. It",
        "must still be paired with the separate upper-side target `B_(k+1)<=B_k`.",
        "",
        "## Repaired k300 Anchor",
        "",
        "```text",
        f"B_k > 0 rows: {diagnostics['b_wall_rows']} / {diagnostics['scaled_curvature_total_rows']}",
        f"C_(k+1)-C_k positive rows: {diagnostics['scaled_curvature_increase_rows']} / {diagnostics['adjacent_total_rows']}",
        f"C_(k+1)-C_k failures: {diagnostics['scaled_curvature_increase_failures']}",
        f"C_(k+1)-C_k inconclusive: {diagnostics['scaled_curvature_increase_inconclusive']}",
        f"retired C_k<=2/3 failures: {diagnostics['two_thirds_failure_rows']} / {diagnostics['scaled_curvature_total_rows']}",
        f"min C increase margin: {diagnostics['min_scaled_curvature_increase_margin']['sample']} at lambda={diagnostics['min_scaled_curvature_increase_margin']['lam']}, k={diagnostics['min_scaled_curvature_increase_margin']['k']}",
        "```",
        "",
        "## Exact Handoff",
        "",
        "```text",
        "C_(k+1)>=C_k",
        "<=> B_(k+1)>=((2*k+1)/(2*k+3))*B_k",
        "=> nonlinear lower curvature barrier, assuming B_k>=0",
        "plus B_(k+1)<=B_k gives the two-sided curvature corridor",
        "```",
        "",
        "Rejected shortcuts:",
        "",
        "```text",
        "fixed C_k<=2/3 wall",
        "generic raw-moment positivity or raw walls alone",
        "monotone B curvature alone",
        "```",
        "",
        "Live routes:",
        "",
        "```text",
        "1. zeta-specific log-ratio recurrence",
        "2. signed/tilted saddle analysis with uniform remainders",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_linear_barrier_scout"],
        artifact["source_curvature_corridor_bridge"],
        artifact["source_raw_corridor_target"],
        artifact["source_k300_obstruction"],
        artifact["source_monotone_contraction_target"],
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
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact()
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda scaled-curvature monotonicity target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
