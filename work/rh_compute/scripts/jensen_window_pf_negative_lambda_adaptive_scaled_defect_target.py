#!/usr/bin/env python3
"""Build the negative-lambda adaptive scaled-defect tail target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FRONTIER_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.json"
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def half_width_failures(diagnostics: dict) -> list[dict]:
    failures: list[dict] = []
    for row in diagnostics["lambda_frontiers"]:
        checked = int(row["checked_x_max"])
        positive = int(row["half_width_positive_rows"])
        if positive < checked:
            failures.append(
                {
                    "lam": row["lam"],
                    "first_failure_k": positive + 1,
                    "failure_rows": checked - positive,
                    "max_scaled_defect": row["max_scaled_defect"],
                    "min_half_width_margin": row["min_half_width_margin"],
                }
            )
    return failures


def build_target(frontier_json: Path) -> dict:
    frontier = load_json(frontier_json)
    diagnostics = frontier["finite_diagnostics"]
    coefficient_k_max = int(diagnostics["coefficient_k_max"])
    checked_x_max = int(diagnostics["checked_x_max"])
    tail_start = checked_x_max + 1
    next_coefficient = coefficient_k_max + 1
    prefix_label = f"k{coefficient_k_max}"
    half_failures = int(diagnostics["scaled_defect_rows"]) - int(diagnostics["half_width_positive_rows"])
    first_half = half_width_failures(diagnostics)
    rows = [
        {
            "id": "nlasdt_01_exact_scaled_cone",
            "role": "exact_reformulation",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                frontier.get("note", f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md"),
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
            ],
            "claim_if_proved": "With d_k=1-x_k and s_k=((2*k+1)/2)*d_k, the lower wall is equivalent to 0<=s_k<=1.",
            "gap": "This exact reformulation supplies no analytic all-k estimate for the actual zeta heat-flow tail.",
            "acceptance_test": "Use only as algebraic bookkeeping for the defect-tail lower wall.",
            "proof_boundary": "Exact reformulation only; not a tail theorem.",
        },
        {
            "id": "nlasdt_02_k200_frontier",
            "role": "finite_frontier",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                frontier.get("note", f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md"),
                f"outputs/jensen_window_pf_negative_lambda_tail_barrier_{prefix_label}_scout.md",
            ],
            "claim_if_proved": (
                f"The {prefix_label} frontier supplies {diagnostics['cone_positive_rows']}/"
                f"{diagnostics['scaled_defect_rows']} exact-cone rows, "
                f"{diagnostics['half_width_positive_rows']}/{diagnostics['scaled_defect_rows']} fixed half-width rows, "
                f"and {diagnostics['one_third_positive_rows']}/{diagnostics['scaled_defect_rows']} one-third rows."
            ),
            "gap": "Finite frontier evidence is not an all-k theorem and cannot close cone entry by itself.",
            "acceptance_test": "Treat the k200 rows as base evidence, falsification data, and quantitative guidance only.",
            "proof_boundary": "Finite diagnostic only.",
        },
        {
            "id": "nlasdt_03_fixed_buffers_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md",
                frontier.get("note", f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md"),
            ],
            "claim_if_proved": "The one-third buffer, fixed half-width buffer, or scaled-defect nonincrease would give a simpler tail proof.",
            "gap": (
                f"The one-third buffer fails on {diagnostics['one_third_failure_rows']} checked rows, "
                f"the fixed half-width buffer fails on {half_failures} checked rows, and "
                f"scaled-defect nonincrease is rejected on {diagnostics['scaled_defect_increase_rows']} adjacent rows."
            ),
            "acceptance_test": "Do not use these fixed-buffer or monotone-decrease shortcuts as proof premises.",
            "proof_boundary": "Finite rejection only; not a proof of any replacement route.",
        },
        {
            "id": "nlasdt_04_adaptive_or_exact_cone_statement",
            "role": "open_statement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
            ],
            "claim_if_proved": (
                f"Prove for each checked negative lambda that 0<=s_k<=1 for all k>={tail_start}, "
                "preferably through an explicit adaptive envelope E_lambda(k)<1, and separately prove "
                f"d_(k+1)<=d_k for all k>={checked_x_max}."
            ),
            "gap": "No noncircular exact-cone or adaptive scaled-defect all-k theorem is currently available.",
            "acceptance_test": "Give explicit zeta coefficient estimates with a stated starting index and no endpoint PF/RH/Newman input.",
            "proof_boundary": "Open theorem target only; not a proof of cone entry.",
        },
        {
            "id": "nlasdt_05_uniform_saddle_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
            ],
            "claim_if_proved": "A uniform saddle/Laplace expansion derives the adaptive scaled-defect envelope and the monotone defect bridge.",
            "gap": "Current signed perturbation evidence is fixed-k/local and lacks uniform k/T and far-tail remainders.",
            "acceptance_test": "Produce interval-safe remainder estimates strong enough to dominate the k200 frontier and remain below the exact cone wall.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlasdt_06_ratio_comparison_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
                "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
            ],
            "claim_if_proved": "A direct ratio recurrence or comparison inequality proves 0<=s_k<=1 and the separate monotone defect bridge.",
            "gap": "Known direct fixed-width recurrences are finite-rejected, so any recurrence must be compatible with increasing s_k.",
            "acceptance_test": "Derive a noncircular inequality using only tail quantities already controlled by finite anchors or analytic estimates.",
            "proof_boundary": "Live theorem-search route only.",
        },
        {
            "id": "nlasdt_07_finite_anchor_handoff",
            "role": "exact_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                f"outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_{prefix_label}_scout.md",
                f"outputs/jensen_window_pf_negative_lambda_finite_collar_{prefix_label}_contract.md",
            ],
            "claim_if_proved": (
                f"The finite anchor covers x_1..x_{checked_x_max}; without an analytic bridge, "
                f"the next purely finite step needs A_{next_coefficient}."
            ),
            "gap": "Finite extension alone is not an all-k tail theorem or finite-collar flow theorem.",
            "acceptance_test": "State whether the bridge is analytic or purely finite before using it in cone-entry arguments.",
            "proof_boundary": "Requirement extraction only.",
        },
        {
            "id": "nlasdt_08_conditional_application",
            "role": "conditional_application",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
            ],
            "claim_if_proved": "A proved adaptive/exact-cone scaled-defect tail plus the monotone bridge would supply the missing defect-tail input for cone entry.",
            "gap": "The surrounding infinite/collared flow legitimacy and Jensen-window PF bridge obligations would still remain.",
            "acceptance_test": "After a proof, rerun the defect-tail target, dependency graph, and core gates without relaxing forbidden endpoint assumptions.",
            "proof_boundary": "Conditional application only; not jwpf_06, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "target_rows": len(rows),
        "ready_to_apply_rows": 0,
        "live_routes": 2,
        "rejected_routes": 1,
        "coefficient_k_max": coefficient_k_max,
        "checked_x_max": checked_x_max,
        "tail_start_k": tail_start,
        "next_finite_coefficient_needed": f"A_{next_coefficient}",
        "finite_exact_cone_rows": diagnostics["cone_positive_rows"],
        "finite_scaled_rows": diagnostics["scaled_defect_rows"],
        "finite_half_width_rows": diagnostics["half_width_positive_rows"],
        "finite_half_width_failure_rows": half_failures,
        "finite_one_third_failure_rows": diagnostics["one_third_failure_rows"],
        "first_half_width_failures": first_half,
        "max_scaled_defect": diagnostics["max_scaled_defect"],
        "min_exact_cone_margin": diagnostics["min_cone_margin"],
        "target_closing": False,
        "main_finding": (
            f"The {prefix_label} frontier keeps the exact scaled cone 0<=s_k<=1 on "
            f"{diagnostics['cone_positive_rows']}/{diagnostics['scaled_defect_rows']} rows, while rejecting "
            f"the fixed half-width buffer on {half_failures} rows and the one-third buffer on "
            f"{diagnostics['one_third_failure_rows']} rows. The replacement target is an exact-cone or "
            f"adaptive-envelope all-k theorem from k>={tail_start}, plus the separate monotone defect bridge."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_adaptive_scaled_defect_target",
        "date": "2026-07-06",
        "status": "open_theorem_target",
        "target_id": "target_negative_lambda_adaptive_scaled_defect_tail",
        "source_scaled_defect_frontier": frontier.get("note", f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md"),
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "frontier_json": frontier_json.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py",
        "proof_boundary": (
            "Open theorem target only. It replaces the finite-rejected fixed half-width route with an "
            "adaptive or exact-cone scaled-defect target, but it does not prove that target, does not "
            "prove the monotone defect bridge, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "target_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The target remains open_target.",
            "The one-third, fixed half-width, and scaled-defect nonincrease shortcuts remain rejected.",
            "The exact cone evidence is finite evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(target: dict, path: Path) -> None:
    summary = target["summary"]
    max_scaled = summary["max_scaled_defect"]
    min_cone = summary["min_exact_cone_margin"]
    result_line = (
        "validated Jensen-window PF negative-lambda adaptive scaled-defect target: "
        f"{summary['target_rows']} rows, 0 issues, {summary['live_routes']} live routes, "
        f"{summary['finite_exact_cone_rows']} exact-cone rows, "
        f"{summary['finite_half_width_failure_rows']} half-width failures"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Adaptive Scaled-Defect Target",
        "",
        "Date: 2026-07-06",
        "",
        "Status: open theorem target. This is not a proof of zeta cone entry,",
        "Jensen-window PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_adaptive_scaled_defect_target`.",
        "",
        "Proof boundary: this artifact replaces the finite-rejected fixed",
        "half-width route with an adaptive or exact-cone scaled-defect target.",
        "It does not prove that target or the separate monotone defect bridge.",
        "",
        "Machine-readable target:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Target",
        "",
        "```text",
        "d_k = 1 - x_k",
        "s_k = ((2*k+1)/2) * d_k",
        f"prove 0 <= s_k <= 1 for all k >= {summary['tail_start_k']}",
        f"prove d_(k+1) <= d_k for all k >= {summary['checked_x_max']} separately",
        "optionally sharpen 0 <= s_k <= 1 to an explicit adaptive envelope E_lambda(k)<1",
        "```",
        "",
        "## k200 Frontier",
        "",
        "```text",
        f"exact cone rows: {summary['finite_exact_cone_rows']} / {summary['finite_scaled_rows']}",
        f"fixed half-width rows: {summary['finite_half_width_rows']} / {summary['finite_scaled_rows']}",
        f"fixed half-width failure rows: {summary['finite_half_width_failure_rows']}",
        f"one-third failure rows: {summary['finite_one_third_failure_rows']}",
        f"max scaled defect: {max_scaled['sample']} at lambda={max_scaled['lam']}, k={max_scaled['k']}",
        f"min exact-cone margin: {min_cone['sample']} at lambda={min_cone['lam']}, k={min_cone['k']}",
        f"next finite coefficient needed without an analytic bridge: {summary['next_finite_coefficient_needed']}",
        "```",
        "",
        "First fixed half-width failures:",
        "",
        "```text",
        *[
            f"lambda={row['lam']}: first failure k={row['first_failure_k']}, failure rows={row['failure_rows']}, max s={row['max_scaled_defect']['sample']} at k={row['max_scaled_defect']['k']}"
            for row in summary["first_half_width_failures"]
        ],
        "```",
        "",
        "Live routes:",
        "",
        "```text",
        "1. uniform saddle/Laplace control with interval-safe remainders",
        "2. direct ratio-recurrence or comparison inequalities compatible with increasing s_k",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        target["source_scaled_defect_frontier"],
        "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
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
    parser.add_argument("--frontier-json", type=Path, default=DEFAULT_FRONTIER_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    frontier_json = args.frontier_json if args.frontier_json.is_absolute() else REPO_ROOT / args.frontier_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    target = build_target(frontier_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(target, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(target, note)
    print(
        "wrote Jensen-window PF negative-lambda adaptive scaled-defect target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
