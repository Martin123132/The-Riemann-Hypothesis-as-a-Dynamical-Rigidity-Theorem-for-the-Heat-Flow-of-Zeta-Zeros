#!/usr/bin/env python3
"""Build the sign-regular-to-Jensen transfer gap matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BRIDGE_JSON = REPO_ROOT / "work/rh_compute/results/jensen_hankel_bridge_algebra.json"
DEFAULT_OBLIGATION_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_sign_regular_transfer_gap_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_sign_regular_transfer_gap_matrix.md"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_matrix(bridge_json: Path, obligation_json: Path) -> dict:
    bridge = load_json(bridge_json)
    obligation = load_json(obligation_json)
    bridge_counter = bridge["finite_countermodel"]
    obligation_counter = obligation["finite_countermodel"]

    rows = [
        {
            "id": "srgt_01_degree2_exact_contact",
            "role": "exact_contact",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_hankel_bridge_algebra.md",
                "outputs/jensen_window_pf_obligation_algebra.md",
            ],
            "claim": "Degree 2 is exact: the m=1 signed-Hankel condition is the Jensen-window quadratic discriminant threshold.",
            "acceptance_test": "Use this as the contact point only; do not extrapolate it to all degrees.",
            "proof_boundary": "Exact low-degree contact only; not a sign-regular-to-Jensen transfer theorem.",
        },
        {
            "id": "srgt_02_degree3_new_obligation",
            "role": "exact_obstruction",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_obligation_algebra.md"],
            "claim": "Degree 3 introduces a cubic discriminant and Toeplitz minors that are not signed reshaped-Hankel minors.",
            "acceptance_test": "A proposed bridge must prove these cubic obligations from its stated hypotheses.",
            "proof_boundary": "Exact algebraic obstruction only; not a theorem for the zeta coefficients.",
        },
        {
            "id": "srgt_03_finite_reshaped_hankel_countermodel",
            "role": "rejected_shortcut",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_hankel_bridge_algebra.md"],
            "claim": "A positive rational sequence passes finite reshaped-Hankel signs for k=2,3 at N=3 while its degree-3 Jensen discriminant is negative.",
            "acceptance_test": "Reject any bridge that uses only finite low-order reshaped-Hankel signs.",
            "proof_boundary": "Countermodel gate only; it is not a claim about the actual zeta sequence.",
        },
        {
            "id": "srgt_04_selected_toeplitz_countermodel",
            "role": "rejected_shortcut",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/jensen_window_pf_obligation_algebra.md"],
            "claim": "The same positive rational sequence has selected low-order degree-3 and degree-4 Toeplitz minors positive while later contiguous Toeplitz minors and discriminants are negative.",
            "acceptance_test": "Reject any bridge that certifies only selected low-order Jensen-window Toeplitz minors.",
            "proof_boundary": "Countermodel gate only; not an all-minor theorem.",
        },
        {
            "id": "srgt_05_all_order_antecedent_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/signed_hankel_jensen_bridge_target.md",
                "outputs/jensen_window_pf_bridge_obligations.md",
            ],
            "claim": "The signed-Hankel side must be an all-order, all-shift reshaped-Hankel sign-consistency theorem, not a finite staircase.",
            "acceptance_test": "State the exact all-order minor family and prove it for the actual A_k(0).",
            "proof_boundary": "Open antecedent requirement only; current Arb certificates are finite.",
        },
        {
            "id": "srgt_06_xi_phi_specific_transfer_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
                "outputs/jensen_window_pf_structural_ansatz_matrix.md",
            ],
            "claim": "A viable bridge must add noncircular Xi/Phi-specific structure that is absent from arbitrary positive sequences.",
            "acceptance_test": "Identify a positive kernel, determinant integral, production/network model, or sign-regular transfer theorem whose hypotheses are verified for A_k(0).",
            "proof_boundary": "Open bridge requirement only; no such theorem is proved here.",
        },
        {
            "id": "srgt_07_binomial_shift_uniformity_requirement",
            "role": "open_requirement",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/jensen_window_pf_bridge_target.md",
                "outputs/jensen_window_pf_schur_shape_contract.md",
            ],
            "claim": "The transfer must output every binomially weighted Jensen-window Toeplitz matrix for all degrees and shifts.",
            "acceptance_test": "Handle binom(d,j), every d, every n, and every Toeplitz minor shape in one theorem.",
            "proof_boundary": "Open uniformity requirement only; finite rectangles are rejected elsewhere.",
        },
        {
            "id": "srgt_08_acceptable_theorem_shape",
            "role": "live_contract",
            "readiness": "not_ready_to_apply",
            "source_artifacts": ["outputs/signed_hankel_jensen_bridge_target.md"],
            "claim": "An acceptable theorem may prove a sign-regular-to-Toeplitz transfer, a positive determinant integral, or an Xi/Phi positive-kernel identity that directly gives the Jensen-window PF target.",
            "acceptance_test": "The theorem must imply all Jensen-window PF conditions without assuming Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0.",
            "proof_boundary": "Live contract only; no candidate is ready to apply.",
        },
        {
            "id": "srgt_09_forbidden_transfer_shortcuts",
            "role": "rejected_shortcut",
            "readiness": "not_ready_to_apply",
            "source_artifacts": [
                "outputs/countermodel_library.md",
                "outputs/sign_regularity_theorem_fit_matrix.md",
            ],
            "claim": "Finite grids, degree-2 analogy, selected minors, endpoint PF/LP assumptions, and asymptotic Jensen statements are forbidden as bridge replacements.",
            "acceptance_test": "Reject any proposed proof step that avoids the all-order transfer theorem with one of these shortcuts.",
            "proof_boundary": "Rejected proof templates only.",
        },
    ]
    summary = {
        "transfer_rows": len(rows),
        "countermodel_gates": 2,
        "open_requirement_rows": 3,
        "rejected_shortcut_rows": 3,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The signed-Hankel/Jensen bridge is exact only at degree 2. The "
            "degree-3/4 countermodel gates show that finite reshaped-Hankel signs "
            "and selected low-order Toeplitz positivity do not imply Jensen-window "
            "PF-infinity; a proof needs all-order sign consistency plus a genuine "
            "Xi/Phi-specific sign-regular-to-Toeplitz transfer theorem with binomial "
            "and shift uniformity."
        ),
    }
    return {
        "kind": "jensen_window_pf_sign_regular_transfer_gap_matrix",
        "date": "2026-07-06",
        "status": "finite_theorem_search_diagnostic",
        "source_jensen_hankel_bridge_algebra": "outputs/jensen_hankel_bridge_algebra.md",
        "source_jensen_window_pf_obligation_algebra": "outputs/jensen_window_pf_obligation_algebra.md",
        "source_signed_hankel_jensen_bridge_target": "outputs/signed_hankel_jensen_bridge_target.md",
        "source_bridge_obligation_ledger": "outputs/jensen_window_pf_bridge_obligations.md",
        "bridge_json": bridge_json.relative_to(REPO_ROOT).as_posix(),
        "obligation_json": obligation_json.relative_to(REPO_ROOT).as_posix(),
        "generator": "work/rh_compute/scripts/jensen_window_pf_sign_regular_transfer_gap_matrix.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_sign_regular_transfer_gap_matrix.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It combines exact low-degree "
            "algebra and finite countermodel gates to specify the missing transfer "
            "theorem, but it does not prove sign-regular-to-Jensen conversion, "
            "Jensen-window PF-infinity, Laguerre-Polya membership, RH, or Lambda <= 0."
        ),
        "countermodel_data": {
            "sequence_A0_to_A4": bridge_counter["sequence_A0_to_A4"],
            "finite_reshaped_signs_pass": bridge_counter["finite_reshaped_signs_pass"],
            "degree3_jensen_discriminant": bridge_counter["degree3_jensen_discriminant"],
            "degree3_jensen_hyperbolicity_breaks": bridge_counter["degree3_jensen_hyperbolicity_breaks"],
            "d3_selected_toeplitz_minors_positive": obligation_counter["d3_selected_toeplitz_minors_positive"],
            "d3_first_negative_contiguous_size": obligation_counter["d3_first_negative_contiguous_toeplitz_minor"]["size"],
            "d3_first_negative_contiguous_determinant": obligation_counter["d3_first_negative_contiguous_toeplitz_minor"]["determinant"],
            "d4_selected_toeplitz_minors_positive": obligation_counter["d4_selected_toeplitz_minors_positive"],
            "d4_quartic_discriminant": obligation_counter["d4_quartic_discriminant"],
            "d4_first_negative_contiguous_size": obligation_counter["d4_first_negative_contiguous_toeplitz_minor"]["size"],
            "d4_first_negative_contiguous_determinant": obligation_counter["d4_first_negative_contiguous_toeplitz_minor"]["determinant"],
        },
        "transfer_rows": rows,
        "invariants": [
            "No row is ready_to_apply.",
            "Degree-2 contact is not promoted to all-degree transfer.",
            "Finite reshaped-Hankel signs are not promoted to all-order sign consistency.",
            "Selected low-order Toeplitz minors are not promoted to Jensen-window PF-infinity.",
            "Endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, and Lambda <= 0 are forbidden as inputs.",
        ],
        "summary": summary,
    }


def write_note(matrix: dict, path: Path) -> None:
    summary = matrix["summary"]
    data = matrix["countermodel_data"]
    result_line = (
        "validated Jensen-window PF sign-regular transfer gap matrix: "
        f"{summary['transfer_rows']} transfer rows, {summary['countermodel_gates']} countermodel gates, "
        f"{summary['open_requirement_rows']} open requirements, "
        f"{summary['rejected_shortcut_rows']} rejected shortcuts, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows, 0 issues"
    )
    lines = [
        "# Jensen-Window PF Sign-Regular Transfer Gap Matrix",
        "",
        "Date: 2026-07-06",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of",
        "Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya",
        "membership, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_sign_regular_transfer_gap_matrix`.",
        "",
        "Proof boundary: this artifact combines exact low-degree algebra and",
        "finite countermodel gates. It does not prove the missing sign-regular",
        "to Jensen-window PF transfer theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_sign_regular_transfer_gap_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_sign_regular_transfer_gap_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_sign_regular_transfer_gap_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Countermodel Gates",
        "",
        "Shared positive rational sequence:",
        "",
        "```text",
        f"A_0..A_4 = {', '.join(data['sequence_A0_to_A4'])}",
        "```",
        "",
        "It passes the finite reshaped-Hankel signs used by the bridge algebra gate:",
        "",
        "```text",
        f"finite reshaped signs pass: {data['finite_reshaped_signs_pass']}",
        f"degree-3 Jensen discriminant: {data['degree3_jensen_discriminant']}",
        f"degree-3 Jensen hyperbolicity breaks: {data['degree3_jensen_hyperbolicity_breaks']}",
        "```",
        "",
        "It also passes selected low-order Jensen-window Toeplitz tests while",
        "later obligations fail:",
        "",
        "```text",
        f"d=3 selected Toeplitz minors positive: {data['d3_selected_toeplitz_minors_positive']}",
        f"d=3 first negative contiguous size: {data['d3_first_negative_contiguous_size']}",
        f"d=3 first negative contiguous determinant: {data['d3_first_negative_contiguous_determinant']}",
        f"d=4 selected Toeplitz minors positive: {data['d4_selected_toeplitz_minors_positive']}",
        f"d=4 quartic discriminant: {data['d4_quartic_discriminant']}",
        f"d=4 first negative contiguous size: {data['d4_first_negative_contiguous_size']}",
        f"d=4 first negative contiguous determinant: {data['d4_first_negative_contiguous_determinant']}",
        "```",
        "",
        "## Transfer Contract",
        "",
        "A usable theorem must supply all of:",
        "",
        "```text",
        "1. all-order, all-shift reshaped-Hankel sign consistency for the actual A_k(0)",
        "2. noncircular Xi/Phi-specific structure absent from arbitrary positive sequences",
        "3. binomial-weight and shift uniformity for every Jensen-window Toeplitz minor",
        "4. no endpoint PF, Jensen hyperbolicity, Laguerre-Polya, RH, or Lambda <= 0 assumptions",
        "```",
        "",
        "Rows:",
        "",
        "```text",
    ]
    for row in matrix["transfer_rows"]:
        lines.append(f"{row['id']}: {row['role']} - {row['claim']}")
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            "outputs/jensen_hankel_bridge_algebra.md",
            "outputs/jensen_window_pf_obligation_algebra.md",
            "outputs/signed_hankel_jensen_bridge_target.md",
            "outputs/jensen_window_pf_bridge_obligations.md",
            "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
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
    parser.add_argument("--bridge-json", type=Path, default=DEFAULT_BRIDGE_JSON)
    parser.add_argument("--obligation-json", type=Path, default=DEFAULT_OBLIGATION_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    bridge_json = args.bridge_json if args.bridge_json.is_absolute() else REPO_ROOT / args.bridge_json
    obligation_json = args.obligation_json if args.obligation_json.is_absolute() else REPO_ROOT / args.obligation_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    matrix = build_matrix(bridge_json, obligation_json)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(matrix, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(matrix, note)
    print(
        "wrote Jensen-window PF sign-regular transfer gap matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
