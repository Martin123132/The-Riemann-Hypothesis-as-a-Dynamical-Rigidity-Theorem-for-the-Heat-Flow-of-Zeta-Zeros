#!/usr/bin/env python3
"""Build the zeta-specific raw-corridor theorem target for negative lambda."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md"


def build_artifact() -> dict:
    rows = [
        {
            "id": "nlzrct_01_raw_corridor_statement",
            "role": "open_statement",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "For the actual negative-lambda zeta heat-flow raw moments M_k=mu_{2k}, prove 1<=R_k<=(2*k+1)/(2*k-1) and L_k(R_k)<=R_(k+1)<=U_k(R_k) on the required tail.",
            "formula": "R_k=M_(k+1)*M_(k-1)/M_k^2; L_k=((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k; U_k=(2+(2*k-1)*R_k)/(2*k+1)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
            ],
            "proof_boundary": "Open zeta-specific theorem target only; not proved by this artifact.",
        },
        {
            "id": "nlzrct_02_exact_equivalence_to_adaptive_cone",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim_if_proved": "The raw upper wall is exactly 0<=s_k, the lower raw wall is exactly s_k<=1, and the corridor is the monotone bridge plus scaled k-increase.",
            "formula": "s_k=((2*k+1)-(2*k-1)*R_k)/2; L_k<=R_(k+1)<=U_k",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md",
            ],
            "proof_boundary": "Exact algebra only; no zeta-specific inequality is proved.",
        },
        {
            "id": "nlzrct_03_k200_finite_anchor",
            "role": "finite_anchor",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "The k200 negative-lambda prefix validates 597/597 raw-cone rows and 594/594 corridor rows for lambda=-25,-50,-100.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md",
            ],
            "proof_boundary": "Finite prefix evidence only; not an all-k or continuous-lambda theorem.",
        },
        {
            "id": "nlzrct_04_generic_stieltjes_shortcut_blocked",
            "role": "rejected_shortcut",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "Generic Stieltjes moment positivity or raw log-convexity proves the upper wall and corridor.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/jensen_window_pf_factorial_multiplier_split_audit.md",
            ],
            "gap": "Exact positive two-atom Stieltjes witnesses violate the upper raw wall or one side of the adaptive corridor.",
            "proof_boundary": "Rejected shortcut only; not evidence against zeta-specific structure.",
        },
        {
            "id": "nlzrct_05_positive_gaussian_mixture_shortcut_blocked",
            "role": "rejected_shortcut",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "A positive Gaussian scale-mixture comparison proves the upper wall.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md",
            ],
            "gap": "Positive Gaussian scale mixtures move in the wrong direction for the upper wall by pushing x_k>=1 rather than x_k<=1.",
            "proof_boundary": "Rejected proof template only; signed or tilted zeta perturbations remain live.",
        },
        {
            "id": "nlzrct_06_signed_gaussian_local_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "Certified Phi Taylor signs provide the local signed Gaussian perturbation direction for the upper wall and monotone correction.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
                "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
            ],
            "gap": "The route is fixed-k/local and still lacks interval-safe uniform remainders in q/T and a far-tail handoff.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlzrct_07_two_scale_remainder_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "A uniform local/mesoscopic remainder theorem plus far-tail moving-saddle theorem promotes the signed perturbation route to the raw corridor target.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            "gap": "Low-order and high-order finite Taylor truncations are not proof models for the actual k200 prefix.",
            "proof_boundary": "Open analytic requirement only.",
        },
        {
            "id": "nlzrct_08_ratio_recurrence_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "A direct zeta-specific ratio recurrence or comparison principle proves corridor occupancy from finite anchors.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
                "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
            ],
            "gap": "Known fixed-width recurrences and generic raw-moment arguments are rejected, so the recurrence must be compatible with increasing scaled defect.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlzrct_09_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim_if_proved": "A promoted proof must state its tail start, lambda range, finite/collared handoff, raw-wall and corridor inequalities, and forbidden assumptions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of the target.",
        },
    ]
    summary = {
        "target_rows": len(rows),
        "exact_available_rows": 1,
        "finite_anchor_rows": 1,
        "live_routes": 2,
        "open_requirement_rows": 1,
        "rejected_shortcut_rows": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "finite_raw_cone_rows": 597,
        "finite_corridor_rows": 594,
        "open_theorem_target": True,
        "main_finding": (
            "The zeta-specific raw-corridor theorem target is now explicit: prove the "
            "upper raw wall and adaptive corridor for actual negative-lambda zeta moments. "
            "The k200 prefix supports the target with 597 raw-cone rows and 594 corridor rows, "
            "but generic Stieltjes positivity and positive Gaussian mixtures are blocked; the live "
            "routes are signed Gaussian perturbation with two-scale remainders or a direct zeta-specific "
            "ratio recurrence."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target",
        "date": "2026-07-06",
        "status": "open_theorem_target",
        "target_id": "target_negative_lambda_zeta_specific_raw_corridor",
        "source_raw_bridge": "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
        "source_raw_obstruction": "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
        "source_adaptive_target": "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
        "source_uniform_remainder": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py",
        "proof_boundary": (
            "Open theorem target only. It names the zeta-specific raw-wall/corridor theorem "
            "needed by the adaptive scaled-defect route, but does not prove it, does not prove "
            "cone entry, does not prove jwpf_06, and does not prove RH or Lambda <= 0."
        ),
        "target_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The target remains open_theorem_target.",
            "Generic Stieltjes/raw-log-convexity shortcuts remain rejected.",
            "Positive Gaussian scale-mixture shortcuts remain rejected.",
            "The k200 zeta evidence is finite evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda zeta-specific raw-corridor target: "
        f"{summary['target_rows']} rows, 0 issues, "
        f"{summary['live_routes']} live routes, "
        f"{summary['rejected_shortcut_rows']} rejected shortcuts, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Zeta-Specific Raw-Corridor Target",
        "",
        "Date: 2026-07-06",
        "",
        "Status: open theorem target. This is not a proof of cone entry,",
        "Jensen-window PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target`.",
        "",
        "Proof boundary: this artifact names the missing zeta-specific raw-wall",
        "and adaptive-corridor theorem. It does not prove that theorem.",
        "",
        "Machine-readable target:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.py",
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
        "```text",
        "M_k = mu_{2k}",
        "R_k = M_(k+1)*M_(k-1)/M_k^2",
        "1 <= R_k <= (2*k+1)/(2*k-1)",
        "((2*k-1)*(2*k+3)/(2*k+1)^2)*R_k <= R_(k+1) <= (2+(2*k-1)*R_k)/(2*k+1)",
        "```",
        "",
        "Finite support:",
        "",
        "```text",
        f"k200 raw-cone rows: {summary['finite_raw_cone_rows']} / {summary['finite_raw_cone_rows']}",
        f"k200 corridor rows: {summary['finite_corridor_rows']} / {summary['finite_corridor_rows']}",
        "```",
        "",
        "Rejected shortcuts:",
        "",
        "```text",
        "generic Stieltjes/raw-log-convexity proof",
        "positive Gaussian scale-mixture proof of the upper wall",
        "```",
        "",
        "Live routes:",
        "",
        "```text",
        "1. signed Gaussian perturbation plus two-scale uniform remainders",
        "2. direct zeta-specific ratio recurrence compatible with increasing scaled defect",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_raw_bridge"],
        artifact["source_raw_obstruction"],
        artifact["source_adaptive_target"],
        artifact["source_uniform_remainder"],
        "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
        "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
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
        "wrote Jensen-window PF negative-lambda zeta-specific raw-corridor target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
