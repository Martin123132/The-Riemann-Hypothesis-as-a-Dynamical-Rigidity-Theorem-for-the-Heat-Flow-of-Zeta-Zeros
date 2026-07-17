#!/usr/bin/env python3
"""Build exact obligations for the negative-lambda adaptive envelope route."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.md"


def build_artifact() -> dict:
    rows = [
        {
            "id": "nlaeo_01_scaled_defect_coordinates",
            "role": "exact_reindexing",
            "readiness": "available_exact",
            "claim": "With x_k=(A_{k-1}*A_{k+1})/A_k^2, d_k=1-x_k, and s_k=((2*k+1)/2)*d_k.",
            "formula": "s_k=((2*k+1)/2)*(1-x_k)",
            "proof_boundary": "Exact coordinate reindexing only; not a tail theorem.",
        },
        {
            "id": "nlaeo_02_lower_wall_equivalence",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The scaled upper cone wall s_k<=1 is equivalent to x_k>=(2*k-1)/(2*k+1), which is already supplied by the boundary-threshold lemma.",
            "formula": "1-s_k=((2*k+1)*x_k-(2*k-1))/2",
            "source_artifacts": [
                "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
            ],
            "proof_boundary": "Exact equivalence plus existing exact lower-threshold input; does not prove x_k<=1 or monotonicity.",
        },
        {
            "id": "nlaeo_03_fixed_buffer_thresholds",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The fixed half-width and one-third buffers are stronger thresholds x_k>=2*k/(2*k+1) and x_k>=1-2/(3*(2*k+1)); finite k200 stress rejects them as global routes.",
            "formula": "s_k<=1/2 iff x_k>=2*k/(2*k+1); s_k<=1/3 iff x_k>=1-2/(3*(2*k+1))",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md",
                "outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md",
            ],
            "proof_boundary": "Exact threshold translation plus finite rejection only; not a replacement proof.",
        },
        {
            "id": "nlaeo_04_nonnegative_defect_upper_wall",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "claim": "The lower side 0<=s_k is exactly the adjacent log-concavity upper-wall condition x_k<=1.",
            "formula": "s_k>=0 iff x_k<=1",
            "source_artifacts": [
                "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
            ],
            "proof_boundary": "Exact equivalence only; adjacent log-concavity for the needed zeta tail remains open.",
        },
        {
            "id": "nlaeo_05_monotone_defect_bridge",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "claim": "The separate monotone defect bridge d_(k+1)<=d_k is exactly x_(k+1)>=x_k.",
            "formula": "d_k-d_(k+1)=x_(k+1)-x_k",
            "source_artifacts": [
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
                "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
            ],
            "proof_boundary": "Exact equivalence only; monotone contractions remain an open theorem target.",
        },
        {
            "id": "nlaeo_06_scaled_k_monotone_identity",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The observed k-increase of scaled defect is equivalent to 2+(2*k+1)*x_k-(2*k+3)*x_(k+1)>=0.",
            "formula": "s_(k+1)-s_k=(2+(2*k+1)*x_k-(2*k+3)*x_(k+1))/2",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md",
            ],
            "proof_boundary": "Exact identity only; the all-k inequality remains open.",
        },
        {
            "id": "nlaeo_07_monotone_envelope_implication",
            "role": "conditional_route",
            "readiness": "not_ready_to_apply",
            "claim": "If s_k is eventually nondecreasing and has limsup bounded by an explicit E_lambda<1, then the scaled upper cone wall follows on the tail.",
            "formula": "s_k nondecreasing and limsup s_k<=E_lambda<1 => s_k<=E_lambda<1",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
            ],
            "proof_boundary": "Conditional route only; no limiting envelope or all-k monotonicity is proved.",
        },
        {
            "id": "nlaeo_08_lambda_order_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "claim": "The finite lambda-order pattern would need a continuous negative-lambda monotonicity theorem, such as s_k(-T1)>=s_k(-T2) for 25<=T1<=T2 on the required ray.",
            "formula": "25<=T1<=T2 and k>=K => s_k(-T1)>=s_k(-T2)",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md",
            ],
            "proof_boundary": "Open theorem-search requirement only; the finite lambda grid does not imply continuous lambda monotonicity.",
        },
        {
            "id": "nlaeo_09_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof must combine exact lower-threshold input, x_k<=1, x_(k+1)>=x_k, and a noncircular upper envelope, while preserving fixed-buffer rejection gates.",
            "formula": "lower threshold exact + upper wall + monotone bridge + envelope < 1",
            "proof_boundary": "Proof-hygiene gate only; not a mathematical proof.",
        },
    ]
    summary = {
        "obligation_rows": len(rows),
        "exact_available_rows": 3,
        "open_requirement_rows": 3,
        "conditional_route_rows": 1,
        "rejected_route_rows": 1,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The adaptive scaled-defect target decomposes exactly: the lower threshold "
            "s_k<=1 is already the boundary-threshold lemma, while 0<=s_k is x_k<=1, "
            "the monotone bridge is x_(k+1)>=x_k, and scaled k-monotonicity is "
            "2+(2*k+1)*x_k-(2*k+3)*x_(k+1)>=0. The remaining proof burden is a "
            "noncircular upper-wall/monotone-envelope theorem compatible with the finite "
            "half-width and one-third rejection gates."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_adaptive_envelope_obligations",
        "date": "2026-07-06",
        "status": "exact algebraic obligation diagnostic",
        "source_adaptive_target": "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
        "source_envelope_matrix": "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md",
        "source_boundary_threshold": "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
        "source_monotone_contraction_target": "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py",
        "proof_boundary": (
            "Exact algebraic obligation diagnostic only. It translates the adaptive scaled-defect route "
            "into exact coefficient-ratio inequalities and records which subclaims remain open, but it "
            "does not prove adjacent log-concavity, monotone contractions, the limiting envelope, cone "
            "entry, jwpf_06, RH, or Lambda <= 0."
        ),
        "obligation_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the adaptive target.",
            "The lower threshold is exact input, not a proof of the whole exact cone.",
            "The upper wall x_k<=1 remains open on the needed tail.",
            "The monotone bridge x_(k+1)>=x_k remains open on the needed tail.",
            "The fixed half-width and one-third buffers remain finite-rejected.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    result_line = (
        "validated Jensen-window PF negative-lambda adaptive envelope obligations: "
        f"{summary['obligation_rows']} obligation rows, 0 issues, "
        f"{summary['exact_available_rows']} exact rows, "
        f"{summary['open_requirement_rows']} open requirements, "
        f"{summary['rejected_route_rows']} rejected routes"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Adaptive Envelope Obligations",
        "",
        "Date: 2026-07-06",
        "",
        "Status: exact algebraic obligation diagnostic. This is not a proof of",
        "the adaptive scaled-defect target, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_adaptive_envelope_obligations`.",
        "",
        "Proof boundary: this artifact translates the adaptive route into exact",
        "ratio inequalities and separates exact inputs from open requirements.",
        "It does not prove the open requirements.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## 2026-07-10 Upper-Wall Handoff",
        "",
        "`outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md` proves",
        "`x_k<=1`, so the obligation `0<=s_k` is no longer open. The algebraic ledger",
        "below remains useful as a historical decomposition, while the current hard",
        "requirement is the adjacent-`k` monotone bridge and its scaled-envelope",
        "strengthening.",
        "",
        "## 2026-07-11 Lambda=-100 Closure Handoff",
        "",
        "The full cone-entry and raw-corridor theorems now discharge every obligation",
        "listed below at `lambda=-100`: both pointwise walls, the adjacent monotone",
        "bridge, and scaled-defect increase. The simultaneous all-three-lambda",
        "formulation remains historical rather than proved.",
        "",
        "## Exact Translations",
        "",
        "```text",
        "x_k=(A_{k-1}*A_{k+1})/A_k^2",
        "d_k=1-x_k",
        "s_k=((2*k+1)/2)*d_k",
        "1-s_k=((2*k+1)*x_k-(2*k-1))/2",
        "s_k>=0 iff x_k<=1",
        "d_k-d_(k+1)=x_(k+1)-x_k",
        "s_(k+1)-s_k=(2+(2*k+1)*x_k-(2*k+3)*x_(k+1))/2",
        "```",
        "",
        "## Obligations",
        "",
        "```text",
        "already exact: x_k >= (2*k-1)/(2*k+1) from the boundary-threshold lemma",
        "still open: x_k <= 1 on the needed tail",
        "still open: x_(k+1) >= x_k on the needed tail",
        "still open: limiting/adaptive envelope E_lambda<1",
        "rejected: fixed half-width and one-third buffers as global routes",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_adaptive_target"],
        artifact["source_envelope_matrix"],
        artifact["source_boundary_threshold"],
        artifact["source_monotone_contraction_target"],
        "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md",
        "outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md",
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
        "wrote Jensen-window PF negative-lambda adaptive envelope obligations: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
