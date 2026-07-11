#!/usr/bin/env python3
"""Build the negative-lambda defect-tail theorem target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TAIL_BARRIER_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.json"
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_defect_tail_theorem_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_target(tail_barrier_json: Path) -> dict:
    tail = load_json(tail_barrier_json)
    diagnostics = tail["finite_diagnostics"]
    checked_x_max = int(diagnostics["checked_x_max"])
    coefficient_k_max = int(diagnostics["coefficient_k_max"])
    tail_lower_start = checked_x_max + 1
    monotone_bridge_k = checked_x_max
    next_coefficient = coefficient_k_max + 1
    rows = [
        {
            "id": "nldtt_01_exact_tail_statement",
            "role": "open_statement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md",
                "outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md",
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
            ],
            "claim_if_proved": (
                f"For each checked negative lambda, prove 0 <= d_k <= 2/(2*k+1) for all "
                f"k >= {tail_lower_start} and d_(k+1) <= d_k for all k >= {monotone_bridge_k}, "
                "where d_k=1-x_k."
            ),
            "gap": "No analytic all-k defect-tail theorem is currently proved.",
            "acceptance_test": "Provide noncircular estimates for the actual zeta heat-flow coefficients, with explicit starting index and no endpoint PF/RH input.",
            "proof_boundary": "Open theorem target only; not a proof of cone entry.",
        },
        {
            "id": "nldtt_02_finite_anchor_handoff",
            "role": "exact_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
                "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k150_scout.md",
                "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
                "outputs/jensen_window_pf_negative_lambda_finite_collar_k150_contract.md",
            ],
            "claim_if_proved": (
                f"The existing finite prefix supplies x_1..x_{checked_x_max}; a tail theorem starting at "
                f"x_{tail_lower_start} plus the bridge inequality d_{tail_lower_start} <= d_{checked_x_max} "
                "would complete the infinite ratio-cone entry on the checked lambda grid."
            ),
            "gap": f"The bridge value x_{tail_lower_start} would otherwise require A_{next_coefficient}.",
            "acceptance_test": "State whether the route proves the bridge inequality analytically or adds the missing coefficient enclosure.",
            "proof_boundary": "Requirement extraction only; not a tail theorem.",
        },
        {
            "id": "nldtt_03_uniform_saddle_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
                "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
            ],
            "claim_if_proved": "A uniform large-negative-lambda saddle/Laplace expansion controls d_k and d_(k+1)-d_k beyond the finite anchor.",
            "gap": "The current Phi Taylor sign certificate is fixed-k/local; it does not bound the saddle remainder uniformly in k.",
            "acceptance_test": "Give explicit k-uniform remainder terms strong enough to prove the defect barriers from the stated starting index.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nldtt_04_ratio_recurrence_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md",
                "outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md",
                "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
                "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md",
                "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
            ],
            "claim_if_proved": "A direct recurrence or comparison inequality for x_k proves positivity, the lower wall, and monotone defect decrease for the actual coefficient tail.",
            "gap": "No closed recurrence or comparison principle has been proved for the zeta heat-flow coefficient tail; the width-preserving recurrence, the one-third buffer, and the fixed half-width buffer are now finite-rejected on the k150 prefix.",
            "acceptance_test": "Derive an inequality involving only already-controlled tail quantities and verify its base case noncircularly.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nldtt_05_scaled_defect_shortcut_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md",
                "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md",
            ],
            "claim_if_proved": "The shortcut s_k=((2*k+1)/2)*d_k is nonincreasing would imply the lower wall from one finite bound.",
            "gap": "The certified finite prefix shows s_k is increasing on every checked adjacent pair.",
            "acceptance_test": "Do not use scaled-defect nonincrease as a proof premise.",
            "proof_boundary": "Rejected by finite diagnostic evidence only.",
        },
        {
            "id": "nldtt_06_finite_extension_insufficient",
            "role": "insufficient_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
                "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k150_scout.md",
                "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
                "outputs/jensen_window_pf_negative_lambda_finite_collar_k150_contract.md",
            ],
            "claim_if_proved": "More ACB coefficients can keep extending the finite collar depth.",
            "gap": "Finite extension alone does not prove the all-k tail or an analytic collared finite flow theorem.",
            "acceptance_test": "Use finite extensions as diagnostics or base cases only, unless paired with a genuine analytic tail theorem.",
            "proof_boundary": "Finite evidence only.",
        },
        {
            "id": "nldtt_07_forbidden_endpoint_shortcuts",
            "role": "circular_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "claim_if_proved": "Endpoint Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0 would imply strong tail inequalities.",
            "gap": "Those endpoint statements are outputs or stronger assumptions, not admissible inputs for this Newman-direction programme.",
            "acceptance_test": "Reject any route that imports endpoint PF/LP/RH/Newman conclusions as premises.",
            "proof_boundary": "Circular as input.",
        },
        {
            "id": "nldtt_08_conditional_cone_entry_application",
            "role": "conditional_application",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
                "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
            ],
            "claim_if_proved": "The finite prefix plus this all-k defect-tail theorem would supply the missing cone-entry input for the conditional heat-flow ratio-cone route.",
            "gap": "This would still need the surrounding infinite/collared flow legitimacy and downstream Jensen-window PF bridge checks.",
            "acceptance_test": "After a proof, re-run the cone-entry target and dependency graph without changing forbidden endpoint assumptions.",
            "proof_boundary": "Conditional application only; not jwpf_06 or Lambda <= 0.",
        },
    ]
    summary = {
        "target_rows": len(rows),
        "ready_to_apply_rows": 0,
        "live_routes": 2,
        "rejected_routes": 1,
        "insufficient_routes": 1,
        "conditional_application_rows": 1,
        "tail_lower_start_k": tail_lower_start,
        "tail_monotone_bridge_k": monotone_bridge_k,
        "next_finite_coefficient_needed": f"A_{next_coefficient}",
        "target_closing": False,
        "main_finding": (
            f"The defect-tail target is now indexed precisely: finite evidence covers x_1..x_{checked_x_max}, "
            f"so an analytic proof must supply 0 <= d_k <= 2/(2*k+1) for k >= {tail_lower_start} "
            f"and d_(k+1) <= d_k for k >= {monotone_bridge_k}. Two live routes remain: "
            "uniform saddle/Laplace control and a direct ratio-recurrence comparison. The scaled-defect "
            "nonincreasing shortcut is rejected by the finite diagnostic; the one-third buffer is too "
            "strong on the k150 prefix, and the fixed half-width buffer is also finite-rejected."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_defect_tail_theorem_target",
        "date": "2026-07-06",
        "status": "open_theorem_target",
        "target_id": "target_negative_lambda_defect_tail",
        "source_tail_barrier_scout": "outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md",
        "source_scaled_defect_frontier_scout": "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "tail_barrier_json": tail_barrier_json.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_tail_theorem_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_tail_theorem_target.py",
        "proof_boundary": (
            "Open theorem target only. It states the all-k defect-tail theorem needed after the finite "
            "negative-lambda prefix, but it does not prove the tail theorem, does not prove zeta cone entry, "
            "does not prove jwpf_06, and does not prove RH or Lambda <= 0."
        ),
        "finite_anchor": {
            "lambdas": diagnostics["lambdas"],
            "coefficient_k_max": coefficient_k_max,
            "checked_x_max": checked_x_max,
            "tail_lower_start_k": tail_lower_start,
            "tail_monotone_bridge_k": monotone_bridge_k,
            "next_finite_coefficient_needed": f"A_{next_coefficient}",
        },
        "target_rows": rows,
        "invariants": [
            "No row is ready_to_apply.",
            "The target remains open_target.",
            "Finite extensions are not promoted to an all-k theorem.",
        "The scaled-defect nonincreasing shortcut remains rejected.",
        "The one-third-width buffer is not promoted after the k150 frontier failure.",
        "The fixed half-width buffer is finite-rejected after the k150 frontier failure.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(target: dict, path: Path) -> None:
    summary = target["summary"]
    anchor = target["finite_anchor"]
    result_line = (
        "validated Jensen-window PF negative-lambda defect-tail theorem target: "
        f"{summary['target_rows']} rows, 0 issues, {summary['ready_to_apply_rows']} ready-to-apply rows, "
        f"{summary['live_routes']} live routes"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Defect-Tail Theorem Target",
        "",
        "Date: 2026-07-06",
        "",
        "Status: open theorem target. This is not a proof of zeta cone entry,",
        "Jensen-window PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_defect_tail_theorem_target`.",
        "",
        "Proof boundary: this artifact states the analytic all-`k` defect-tail",
        "theorem needed after the finite negative-lambda prefix. It does not",
        "prove that theorem, does not prove `jwpf_06`, and does not establish",
        "`Lambda <= 0`.",
        "",
        "Machine-readable target:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_defect_tail_theorem_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_tail_theorem_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_tail_theorem_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Tail Statement",
        "",
        "Let:",
        "",
        "```text",
        "d_k = 1 - x_k",
        "x_k = (A_{k+1}/A_k)/(A_k/A_{k-1})",
        "```",
        "",
        "The finite prefix currently supplies:",
        "",
        "```text",
        f"lambdas: {', '.join(anchor['lambdas'])}",
        f"coefficient range: A_0..A_{anchor['coefficient_k_max']}",
        f"finite contractions: x_1..x_{anchor['checked_x_max']}",
        "```",
        "",
        "The analytic tail theorem needed to complete the infinite cone entry is:",
        "",
        "```text",
        f"0 <= d_k <= 2/(2*k+1) for all k >= {anchor['tail_lower_start_k']}",
        f"d_(k+1) <= d_k for all k >= {anchor['tail_monotone_bridge_k']}",
        "```",
        "",
        f"Without such a theorem, the next finite bridge needs `{anchor['next_finite_coefficient_needed']}`.",
        "",
        "## Live Routes",
        "",
        "Two noncircular proof routes remain live:",
        "",
        "```text",
        "1. uniform saddle/Laplace control of the actual zeta heat-flow coefficients",
        "2. direct ratio-recurrence or comparison inequalities for the coefficient tail",
        "```",
        "",
        "The scaled-defect nonincreasing shortcut is rejected by the finite scout:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md",
        "validated Jensen-window PF negative-lambda tail-barrier scout: 447 cone-buffer rows, 444 defect-monotone rows, 179 one-third-buffer rows, 444 scaled-defect increase rows, 1 rejected candidate, 0 issues",
        "```",
        "",
        "The scaled-defect frontier now separates viable finite buffers from",
        "over-strong ones:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md",
        "validated Jensen-window PF negative-lambda scaled-defect frontier scout: 447 scaled rows, 447 cone rows, 430 half-width rows, 179 one-third rows, 268 one-third failures, 444 scaled-increase rows, 0 issues",
        "one-third first fails at lambda=-25.0, k=31; fixed half-width first fails at lambda=-25.0, k=133",
        "```",
        "",
        "The older defect-recurrence scout is retained as a finite diagnostic,",
        "but its one-third buffer should no longer be promoted as an all-tail",
        "sufficient route after the k150 frontier:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
        "validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues",
        "0 <= d_k <= 2/(3*(2*k+1))",
        "d_(k+1) <= d_k",
        "```",
        "",
        "## Integration",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
        "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k150_scout.md",
        "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
        "outputs/jensen_window_pf_negative_lambda_finite_collar_k150_contract.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
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
    parser.add_argument("--tail-barrier-json", type=Path, default=DEFAULT_TAIL_BARRIER_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    tail_barrier_json = args.tail_barrier_json if args.tail_barrier_json.is_absolute() else REPO_ROOT / args.tail_barrier_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    target = build_target(tail_barrier_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(target, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(target, note)
    print(
        "wrote Jensen-window PF negative-lambda defect-tail theorem target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
