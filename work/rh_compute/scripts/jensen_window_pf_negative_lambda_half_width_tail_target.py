#!/usr/bin/env python3
"""Build the negative-lambda half-width scaled-defect tail target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_FRONTIER_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.json"
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_half_width_tail_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_target(frontier_json: Path) -> dict:
    frontier = load_json(frontier_json)
    diagnostics = frontier["finite_diagnostics"]
    coefficient_k_max = int(diagnostics["coefficient_k_max"])
    checked_x_max = int(diagnostics["checked_x_max"])
    prefix_label = f"k{coefficient_k_max}"
    tail_lower_start = checked_x_max + 1
    monotone_bridge_k = checked_x_max
    next_coefficient = coefficient_k_max + 1
    half_margin = diagnostics["min_half_width_margin"]
    max_scaled = diagnostics["max_scaled_defect"]
    one_third = diagnostics["min_one_third_margin"]
    half_width_failures = diagnostics["scaled_defect_rows"] - diagnostics["half_width_positive_rows"]
    half_width_rejected = half_width_failures > 0
    first_half_failures = []
    for row in diagnostics["lambda_frontiers"]:
        local_failures = int(row["checked_x_max"]) - int(row["half_width_positive_rows"])
        if local_failures > 0:
            first_half_failures.append(
                {
                    "lam": row["lam"],
                    "first_failure_k": int(row["half_width_positive_rows"]) + 1,
                    "failure_rows": local_failures,
                    "max_scaled_defect": row["max_scaled_defect"],
                    "min_half_width_margin": row["min_half_width_margin"],
                }
            )

    rows = [
        {
            "id": "nlhwt_01_scaled_defect_definition",
            "role": "exact_reformulation",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md",
                f"outputs/jensen_window_pf_negative_lambda_tail_barrier_{prefix_label}_scout.md",
            ],
            "claim_if_proved": "With d_k=1-x_k and s_k=((2*k+1)/2)*d_k, the lower wall x_k>=(2*k-1)/(2*k+1) is equivalent to s_k<=1.",
            "gap": "This is only a change of variables; it supplies no analytic estimate for the actual zeta heat-flow tail.",
            "acceptance_test": "Use this row only as algebraic bookkeeping for tail inequalities.",
            "proof_boundary": "Exact reformulation only; not a tail theorem.",
        },
        {
            "id": "nlhwt_02_half_width_tail_statement",
            "role": "rejected_route" if half_width_rejected else "open_statement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
                f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md",
            ],
            "claim_if_proved": (
                f"For each checked negative lambda, prove 0 <= s_k <= 1/2 for all k >= {tail_lower_start}; "
                "equivalently 0 <= d_k <= 1/(2*k+1)."
            ),
            "gap": (
                f"Finite stress rejects this half-width statement on {half_width_failures}/"
                f"{diagnostics['scaled_defect_rows']} checked {prefix_label} rows."
                if half_width_rejected
                else "No noncircular analytic proof of the half-width tail bound is currently available."
            ),
            "acceptance_test": (
                "Do not use the half-width bound as a proof premise unless the statement is weakened or the finite failure is explained by a corrected target."
                if half_width_rejected
                else "Provide explicit coefficient estimates for the actual zeta heat-flow tail, with no endpoint PF/RH/Newman input and a stated starting index."
            ),
            "proof_boundary": (
                "Finite-rejected route only; not a live half-width theorem target."
                if half_width_rejected
                else "Open theorem target only; not a proof of the half-width tail bound."
            ),
        },
        {
            "id": "nlhwt_03_monotone_defect_still_required",
            "role": "exact_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
                "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
            ],
            "claim_if_proved": (
                f"The half-width lower-wall bound plus d_(k+1)<=d_k for all k >= {monotone_bridge_k} "
                "would supply the defect-tail inequalities needed by the negative-lambda cone-entry route."
            ),
            "gap": (
                "The half-width bound is already finite-rejected, and in any case would not by itself prove defect monotonicity or the bridge inequality at the finite/tail join."
                if half_width_rejected
                else "The half-width bound does not by itself prove defect monotonicity or the bridge inequality at the finite/tail join."
            ),
            "acceptance_test": "Any application must separately prove d_(k+1)<=d_k from the stated bridge index or explain how the monotone wall is supplied.",
            "proof_boundary": "Requirement extraction only; not a proof of cone entry.",
        },
        {
            "id": "nlhwt_04_finite_anchor",
            "role": "finite_anchor",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md",
                f"outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_{prefix_label}_scout.md",
                f"outputs/jensen_window_pf_negative_lambda_finite_collar_{prefix_label}_contract.md",
            ],
            "claim_if_proved": (
                f"The k{coefficient_k_max} finite prefix supplies checked contractions x_1..x_{checked_x_max} "
                f"with {diagnostics['half_width_positive_rows']}/{diagnostics['scaled_defect_rows']} half-width rows."
            ),
            "gap": f"The next purely finite bridge would require A_{next_coefficient}; finite extension alone is not an all-k theorem.",
            "acceptance_test": "Treat the finite rows as base evidence or falsification data only; do not promote them without an analytic tail proof.",
            "proof_boundary": "Finite evidence only; not an all-k half-width theorem.",
        },
        {
            "id": "nlhwt_05_uniform_saddle_route",
            "role": "repair_route" if half_width_rejected else "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
                "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
            ],
            "claim_if_proved": (
                "A uniform saddle/Laplace expansion bounds the exact scaled defect below a corrected adaptive threshold and controls the monotone defect bridge beyond the finite anchor."
                if half_width_rejected
                else "A uniform saddle/Laplace expansion bounds s_k below 1/2 and controls the monotone defect bridge beyond the finite anchor."
            ),
            "gap": (
                "The fixed half-width threshold is finite-rejected; the saddle route must now derive the actual limiting/adaptive frontier rather than force s_k<=1/2."
                if half_width_rejected
                else "Current Taylor and saddle diagnostics are fixed-k or local/mesoscopic; they do not yet give a k-uniform half-width remainder."
            ),
            "acceptance_test": (
                f"Give explicit error terms whose sign and size keep 0<=s_k<=1 from k>={tail_lower_start} and replace the failed half-width threshold with a valid adaptive bound."
                if half_width_rejected
                else f"Give explicit error terms whose sign and size are strong enough to keep s_k<=1/2 from k>={tail_lower_start} on the checked lambda grid."
            ),
            "proof_boundary": "Repair-route theorem search only; not a proof." if half_width_rejected else "Live theorem-search route only.",
        },
        {
            "id": "nlhwt_06_ratio_recurrence_route",
            "role": "repair_route" if half_width_rejected else "live_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
                f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md",
                "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
            ],
            "claim_if_proved": (
                "A direct recurrence or comparison inequality proves the exact cone 0<=s_k<=1 and the separate monotone defect bridge without the failed half-width threshold."
                if half_width_rejected
                else "A direct recurrence or comparison inequality keeps s_k<=1/2 and proves the separate monotone defect bridge."
            ),
            "gap": (
                "The direct width-preserving recurrence, scaled-defect nonincrease, one-third buffer, and now the fixed half-width buffer are finite-rejected."
                if half_width_rejected
                else "The direct width-preserving recurrence is finite-rejected, and scaled-defect nonincrease is also finite-rejected."
            ),
            "acceptance_test": "Derive a noncircular inequality compatible with increasing scaled defect and the k150 half-width failure.",
            "proof_boundary": "Repair-route theorem search only; not a proof." if half_width_rejected else "Live theorem-search route only.",
        },
        {
            "id": "nlhwt_07_falsification_protocol",
            "role": "stress_protocol",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md",
                "outputs/core_proof_programme_gates.md",
            ],
            "claim_if_proved": (
                "Larger finite prefixes should falsify or refine any proposed fixed/adaptive scaled-defect buffer before proof promotion."
                if half_width_rejected
                else "Larger finite prefixes should keep the half-width margin positive if this route is to remain plausible."
            ),
            "gap": (
                f"The minimum {prefix_label} half-width margin is already negative, so the fixed half-width route is rejected rather than merely unproved."
                if half_width_rejected
                else f"The minimum {prefix_label} half-width margin is finite-only and could still fail at larger k or under a uniform asymptotic estimate."
            ),
            "acceptance_test": f"Stress k beyond {tail_lower_start} and compare the observed scaled-defect frontier with any proposed analytic asymptotic bound.",
            "proof_boundary": "Falsification protocol only; finite stress rejects this threshold but cannot prove a replacement theorem.",
        },
        {
            "id": "nlhwt_08_rejected_shortcuts",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md",
                f"outputs/jensen_window_pf_negative_lambda_tail_barrier_{prefix_label}_scout.md",
            ],
            "claim_if_proved": "The older one-third-width buffer, fixed half-width buffer, or scaled-defect nonincrease shortcuts would simplify the route.",
            "gap": (
                f"The one-third buffer fails on {diagnostics['one_third_failure_rows']}/"
                f"{diagnostics['scaled_defect_rows']} {prefix_label} rows, the half-width buffer fails on "
                f"{half_width_failures}/{diagnostics['scaled_defect_rows']} rows, and s_k is increasing on all "
                f"{diagnostics['scaled_defect_increase_rows']} checked adjacent rows."
            ),
            "acceptance_test": "Do not use one-third-width, fixed half-width, or scaled-defect nonincrease as global proof premises.",
            "proof_boundary": "Rejected by finite diagnostic evidence only.",
        },
        {
            "id": "nlhwt_09_conditional_application",
            "role": "conditional_application",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
            ],
            "claim_if_proved": (
                "A proved corrected scaled-defect buffer plus a proved monotone defect bridge would sharpen the negative-lambda defect-tail theorem target."
                if half_width_rejected
                else "A proved half-width tail bound plus a proved monotone defect bridge would sharpen the negative-lambda defect-tail theorem target."
            ),
            "gap": "Even after a corrected target, the surrounding infinite/collared heat-flow and Jensen-window PF bridge obligations would remain.",
            "acceptance_test": "After a corrected proof, rerun the defect-tail target, cone-entry target, dependency graph, and core gates without changing forbidden endpoint assumptions.",
            "proof_boundary": "Conditional application only; not jwpf_06, RH, or Lambda <= 0.",
        },
    ]

    summary = {
        "target_rows": len(rows),
        "ready_to_apply_rows": 0,
        "live_routes": 0 if half_width_rejected else 2,
        "repair_routes": 2 if half_width_rejected else 0,
        "rejected_routes": 2 if half_width_rejected else 1,
        "stress_protocol_rows": 1,
        "conditional_application_rows": 1,
        "coefficient_k_max": coefficient_k_max,
        "checked_x_max": checked_x_max,
        "tail_lower_start_k": tail_lower_start,
        "tail_monotone_bridge_k": monotone_bridge_k,
        "next_finite_coefficient_needed": f"A_{next_coefficient}",
        "finite_half_width_rows": diagnostics["half_width_positive_rows"],
        "finite_half_width_failure_rows": half_width_failures,
        "finite_scaled_rows": diagnostics["scaled_defect_rows"],
        "finite_one_third_failure_rows": diagnostics["one_third_failure_rows"],
        "first_half_width_failures": first_half_failures,
        "min_half_width_margin": half_margin,
        "max_scaled_defect": max_scaled,
        "min_one_third_margin": one_third,
        "target_rejected": half_width_rejected,
        "target_closing": False,
        "main_finding": (
            (
                f"The k{coefficient_k_max} scaled-defect frontier finite-rejects the fixed half-width "
                f"bound s_k<=1/2: only {diagnostics['half_width_positive_rows']}/"
                f"{diagnostics['scaled_defect_rows']} checked rows pass, with "
                f"{half_width_failures} half-width failures and minimum margin {half_margin['sample']} "
                f"at lambda={half_margin['lam']}, k={half_margin['k']}. The route must be replaced by "
                "an exact-cone or adaptive scaled-defect tail target plus the separate monotone defect bridge."
            )
            if half_width_rejected
            else (
                f"The k{coefficient_k_max} scaled-defect frontier keeps the half-width bound "
                f"s_k<=1/2 positive on {diagnostics['half_width_positive_rows']}/"
                f"{diagnostics['scaled_defect_rows']} checked rows, with minimum finite margin "
                f"{half_margin['sample']} at lambda={half_margin['lam']}, k={half_margin['k']}. "
                "The half-width route is now an explicit open theorem target, but it still needs "
                "a noncircular all-k proof and a separate defect-monotone bridge; the one-third "
                f"buffer is rejected by {diagnostics['one_third_failure_rows']} finite failures."
            )
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_half_width_tail_target",
        "date": "2026-07-06",
        "status": "finite_rejected_target" if half_width_rejected else "open_theorem_target",
        "target_id": "target_negative_lambda_half_width_tail",
        "source_scaled_defect_frontier": frontier.get("note", f"outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_{prefix_label}_scout.md"),
        "source_defect_tail_target": "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
        "source_cone_entry_target": "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
        "frontier_json": frontier_json.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_half_width_tail_target.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_half_width_tail_target.py",
        "proof_boundary": (
            (
                "Finite-rejected target record. It shows the fixed half-width scaled-defect route "
                "fails on the checked prefix, and it does not prove a replacement tail theorem, "
                "does not prove the monotone defect bridge, does not prove zeta cone entry, "
                "does not prove jwpf_06, and does not prove RH or Lambda <= 0."
            )
            if half_width_rejected
            else (
                "Open theorem target only. It isolates the half-width scaled-defect tail route, "
                "but it does not prove the half-width theorem, does not prove the monotone "
                "defect bridge, does not prove zeta cone entry, does not prove jwpf_06, "
                "and does not prove RH or Lambda <= 0."
            )
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
            "The target remains finite_rejected_target." if half_width_rejected else "The target remains open_target.",
            "The fixed half-width buffer is rejected by finite evidence." if half_width_rejected else "The half-width buffer is not promoted from finite evidence to an all-k theorem.",
            "The monotone defect bridge remains a separate required theorem.",
            "The one-third-width, fixed half-width, and scaled-defect nonincrease shortcuts remain rejected." if half_width_rejected else "The one-third-width and scaled-defect nonincrease shortcuts remain rejected.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(target: dict, path: Path) -> None:
    summary = target["summary"]
    anchor = target["finite_anchor"]
    half_margin = summary["min_half_width_margin"]
    max_scaled = summary["max_scaled_defect"]
    one_third = summary["min_one_third_margin"]
    target_rejected = bool(summary.get("target_rejected"))
    failure_suffix = (
        f", {summary['finite_half_width_failure_rows']} half-width failures"
        if target_rejected
        else ""
    )
    result_line = (
        "validated Jensen-window PF negative-lambda half-width tail target: "
        f"{summary['target_rows']} rows, 0 issues, {summary['ready_to_apply_rows']} ready-to-apply rows, "
        f"{summary['live_routes']} live routes, {summary['finite_half_width_rows']} half-width rows"
        f"{failure_suffix}"
    )
    status_line = (
        "Status: finite-rejected target. This is not a proof of a replacement",
        "scaled-defect tail theorem, zeta cone entry, Jensen-window PF-infinity,",
    ) if target_rejected else (
        "Status: open theorem target. This is not a proof of zeta cone entry,",
        "Jensen-window PF-infinity, RH, or `Lambda <= 0`.",
    )
    boundary_lines = (
        "Proof boundary: this artifact records that the fixed half-width",
        "scaled-defect route is rejected by finite stress. It does not prove",
        "a corrected route, does not prove the separate monotone defect bridge,",
        "and does not establish `Lambda <= 0`.",
    ) if target_rejected else (
        "Proof boundary: this artifact isolates the half-width scaled-defect",
        "tail route. It does not prove that route, does not prove the separate",
        "monotone defect bridge, and does not establish `Lambda <= 0`.",
    )
    target_intro = (
        "The fixed half-width tail target was:",
        "",
    ) if target_rejected else (
        "The half-width tail target is:",
        "",
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Half-Width Tail Target",
        "",
        "Date: 2026-07-06",
        "",
        *status_line,
        "RH, or `Lambda <= 0`." if target_rejected else "",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_half_width_tail_target`.",
        "",
        *boundary_lines,
        "",
        "Machine-readable target:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_half_width_tail_target.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_half_width_tail_target.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_half_width_tail_target.py",
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
        "d_k = 1 - x_k",
        "s_k = ((2*k+1)/2) * d_k",
        "```",
        "",
        *target_intro,
        "```text",
        f"0 <= s_k <= 1/2 for all k >= {anchor['tail_lower_start_k']}",
        f"d_(k+1) <= d_k for all k >= {anchor['tail_monotone_bridge_k']} remains separately required",
        "```",
        "",
        (
            "Finite k150 stress rejects this fixed half-width target; it is retained only as a documented failed shortcut."
            if target_rejected
            else "This remains an open target, not a proved tail theorem."
        ),
        "",
        "The finite prefix currently supplies:",
        "",
        "```text",
        f"lambdas: {', '.join(anchor['lambdas'])}",
        f"coefficient range: A_0..A_{anchor['coefficient_k_max']}",
        f"finite contractions: x_1..x_{anchor['checked_x_max']}",
        f"next finite coefficient needed: {anchor['next_finite_coefficient_needed']}",
        "```",
        "",
        "## Finite Evidence",
        "",
        "```text",
        f"half-width rows: {summary['finite_half_width_rows']} / {summary['finite_scaled_rows']}",
        f"half-width failure rows: {summary['finite_half_width_failure_rows']}",
        f"one-third failure rows: {summary['finite_one_third_failure_rows']}",
        f"minimum half-width margin: {half_margin['sample']} at lambda={half_margin['lam']}, k={half_margin['k']}",
        f"maximum scaled defect: {max_scaled['sample']} at lambda={max_scaled['lam']}, k={max_scaled['k']}",
        f"minimum one-third margin: {one_third['sample']} at lambda={one_third['lam']}, k={one_third['k']}",
        "```",
        "",
        *(
            [
                "First half-width failures inferred from the increasing scaled-defect frontier:",
                "",
                "```text",
                *[
                    (
                        f"lambda={row['lam']}: first failure k={row['first_failure_k']}, "
                        f"failure rows={row['failure_rows']}, max s={row['max_scaled_defect']['sample']} at k={row['max_scaled_defect']['k']}"
                    )
                    for row in summary.get("first_half_width_failures", [])
                ],
                "```",
                "",
            ]
            if target_rejected
            else []
        ),
        "## Live Routes",
        "",
        "```text",
        (
            "none for the fixed half-width threshold; repair routes must target the exact cone or an adaptive scaled-defect buffer"
            if target_rejected
            else "1. uniform saddle/Laplace estimates for s_k and the monotone bridge"
        ),
        *([] if target_rejected else ["2. direct ratio-recurrence or comparison inequalities compatible with increasing s_k"]),
        "```",
        "",
        "Rejected shortcuts:",
        "",
        "```text",
        "one-third-width buffer",
        *(("fixed half-width buffer",) if target_rejected else ()),
        "scaled-defect nonincrease",
        "endpoint PF/RH/Newman assumptions",
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
        "wrote Jensen-window PF negative-lambda half-width tail target: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
