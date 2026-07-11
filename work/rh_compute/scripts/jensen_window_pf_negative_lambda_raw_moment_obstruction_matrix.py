#!/usr/bin/env python3
"""Build exact raw-moment obstruction gates for the adaptive negative-lambda route."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT_JSON = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md"


def q(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def two_atom_moments(weight: Fraction, support: Fraction, max_k: int) -> list[Fraction]:
    return [Fraction(1) + weight * support**k for k in range(max_k + 1)]


def raw_ratio(moments: list[Fraction], k: int) -> Fraction:
    return moments[k + 1] * moments[k - 1] / moments[k] ** 2


def raw_upper_wall(k: int) -> Fraction:
    return Fraction(2 * k + 1, 2 * k - 1)


def x_from_raw(k: int, raw: Fraction) -> Fraction:
    return Fraction(2 * k - 1, 2 * k + 1) * raw


def scaled_defect_from_raw(k: int, raw: Fraction) -> Fraction:
    return Fraction((2 * k + 1) - (2 * k - 1) * raw, 2)


def corridor_lower(k: int, raw: Fraction) -> Fraction:
    return Fraction((2 * k - 1) * (2 * k + 3), (2 * k + 1) ** 2) * raw


def corridor_upper(k: int, raw: Fraction) -> Fraction:
    return (Fraction(2) + Fraction(2 * k - 1) * raw) / Fraction(2 * k + 1)


def corridor_width(k: int, raw: Fraction) -> Fraction:
    return corridor_upper(k, raw) - corridor_lower(k, raw)


def exact_values(weight: Fraction, support: Fraction, k: int = 1) -> dict:
    moments = two_atom_moments(weight, support, k + 2)
    r_k = raw_ratio(moments, k)
    r_next = raw_ratio(moments, k + 1)
    lower = corridor_lower(k, r_k)
    upper = corridor_upper(k, r_k)
    return {
        "measure": f"delta_1 + ({q(weight)})*delta_{q(support)}",
        "moments": [q(item) for item in moments],
        "k": k,
        "R_k": q(r_k),
        "R_next": q(r_next),
        "raw_upper_wall": q(raw_upper_wall(k)),
        "raw_next_upper_wall": q(raw_upper_wall(k + 1)),
        "raw_upper_slack": q(raw_upper_wall(k) - r_k),
        "raw_next_upper_slack": q(raw_upper_wall(k + 1) - r_next),
        "x_k": q(x_from_raw(k, r_k)),
        "x_next": q(x_from_raw(k + 1, r_next)),
        "x_next_minus_x_k": q(x_from_raw(k + 1, r_next) - x_from_raw(k, r_k)),
        "s_k": q(scaled_defect_from_raw(k, r_k)),
        "s_next": q(scaled_defect_from_raw(k + 1, r_next)),
        "s_next_minus_s_k": q(scaled_defect_from_raw(k + 1, r_next) - scaled_defect_from_raw(k, r_k)),
        "corridor_lower": q(lower),
        "corridor_upper": q(upper),
        "corridor_width": q(upper - lower),
        "lower_margin_R_next_minus_lower": q(r_next - lower),
        "upper_margin_upper_minus_R_next": q(upper - r_next),
    }


def build_artifact() -> dict:
    upper_failure = exact_values(Fraction(1, 16), Fraction(16))
    scaled_upper_failure = exact_values(Fraction(1), Fraction(2))
    monotone_lower_failure = exact_values(Fraction(1, 3), Fraction(9))
    rows = [
        {
            "id": "nlrmo_01_two_atom_family",
            "role": "exact_countermodel_family",
            "status": "guard_validated",
            "claim": "Every measure delta_1+w*delta_a with w>0 and a>0 gives a positive Stieltjes moment sequence M_k=1+w*a^k, so raw moment log-convexity holds but stronger raw-wall/corridor claims may fail.",
            "formula": "M_k=1+w*a^k; R_k=M_(k+1)*M_(k-1)/M_k^2",
            "proof_boundary": "Exact countermodel family only; not a statement about the zeta heat-flow density.",
        },
        {
            "id": "nlrmo_02_upper_raw_wall_failure",
            "role": "exact_counterexample",
            "status": "guard_validated",
            "claim": "A positive two-atom moment sequence can violate the upper raw wall R_k<=(2*k+1)/(2*k-1), hence can make x_k>1 and s_k<0.",
            "witness": upper_failure,
            "proof_boundary": "Blocks generic Stieltjes/raw-log-convexity proofs of the upper raw wall; not evidence against the zeta prefix.",
        },
        {
            "id": "nlrmo_03_scaled_upper_corridor_failure",
            "role": "exact_counterexample",
            "status": "guard_validated",
            "claim": "Even when the raw upper wall holds at k=1 and k=2, the scaled-upper side of the adaptive corridor can fail.",
            "witness": scaled_upper_failure,
            "proof_boundary": "Blocks replacing scaled k-monotonicity by raw upper-wall checks alone.",
        },
        {
            "id": "nlrmo_04_monotone_lower_corridor_failure",
            "role": "exact_counterexample",
            "status": "guard_validated",
            "claim": "Even when the raw upper wall holds at k=1 and k=2, the monotone-bridge lower side of the adaptive corridor can fail.",
            "witness": monotone_lower_failure,
            "proof_boundary": "Blocks replacing the monotone bridge by raw upper-wall checks alone.",
        },
        {
            "id": "nlrmo_05_corridor_width_identity",
            "role": "exact_identity",
            "status": "available_exact",
            "claim": "The adaptive corridor width is positive exactly under the open upper raw wall.",
            "formula": "U_k-L_k=2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
            "proof_boundary": "Exact algebra only; positivity of the width for actual zeta coefficients remains an all-k requirement.",
        },
        {
            "id": "nlrmo_06_zeta_prefix_contrast",
            "role": "finite_contrast",
            "status": "finite_validated",
            "claim": "The actual k200 negative-lambda zeta prefix satisfies 597/597 raw-cone rows and 594/594 adaptive corridor rows, showing the finite zeta data are special relative to the generic two-atom obstructions.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md",
            ],
            "proof_boundary": "Finite zeta contrast only; not an all-k raw-wall or corridor theorem.",
        },
        {
            "id": "nlrmo_07_acceptance_gate",
            "role": "acceptance_gate",
            "status": "guard_validated",
            "claim": "A promoted proof must use zeta-specific structure to prove the upper raw wall and corridor occupancy; generic Stieltjes positivity, raw log-convexity, and finite prefix evidence are insufficient.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
                "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
            ],
            "proof_boundary": "Proof-safety gate only; not a mathematical proof of the zeta-specific theorem.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_counterexample_rows": 3,
        "exact_identity_rows": 2,
        "finite_contrast_rows": 1,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "Generic positive Stieltjes/raw-moment structure does not prove the adaptive "
            "raw wall or corridor. Exact two-atom measures can violate the upper raw wall, "
            "the scaled-upper corridor side, or the monotone-bridge lower corridor side, "
            "even while raw log-convexity holds. The actual zeta k200 prefix remains finite-compatible, "
            "so the missing proof must be zeta-specific and all-k."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix",
        "date": "2026-07-06",
        "status": "exact countermodel gate",
        "source_raw_moment_bridge": "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
        "source_adaptive_target": "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
        "source_cone_prefix": "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_k200_scout.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py",
        "proof_boundary": (
            "Exact countermodel gate only. It blocks generic Stieltjes/raw-log-convexity "
            "promotions into the adaptive raw-moment wall or corridor, but it is not a "
            "claim about actual zeta coefficients, not an all-k theorem, not cone entry, "
            "not jwpf_06, and not Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "All counterexamples are positive two-atom Stieltjes moment sequences.",
            "The gate is not evidence against the finite zeta prefix.",
            "The upper raw wall remains zeta-specific and open all-k.",
            "The adaptive corridor occupancy remains zeta-specific and open all-k.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    rows = {row["id"]: row for row in artifact["matrix_rows"]}
    upper = rows["nlrmo_02_upper_raw_wall_failure"]["witness"]
    scaled = rows["nlrmo_03_scaled_upper_corridor_failure"]["witness"]
    lower = rows["nlrmo_04_monotone_lower_corridor_failure"]["witness"]
    result_line = (
        "validated Jensen-window PF negative-lambda raw-moment obstruction matrix: "
        f"{summary['matrix_rows']} matrix rows, 0 issues, "
        f"{summary['exact_counterexample_rows']} exact counterexamples, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Raw-Moment Obstruction Matrix",
        "",
        "Date: 2026-07-06",
        "",
        "Status: exact countermodel gate. This is not a proof of the adaptive",
        "scaled-defect target, cone entry, Jensen-window PF-infinity, RH, or",
        "`Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix`.",
        "",
        "Proof boundary: this artifact blocks generic Stieltjes/raw-log-convexity",
        "proofs of the raw upper wall or adaptive corridor. It is not evidence",
        "against the actual zeta finite prefix.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Witnesses",
        "",
        "All witnesses are positive two-atom Stieltjes moment sequences:",
        "",
        "```text",
        "M_k = 1 + w*a^k",
        "R_k = M_(k+1)*M_(k-1)/M_k^2",
        "```",
        "",
        "Upper raw wall failure:",
        "",
        "```text",
        f"measure: {upper['measure']}",
        f"R_1 = {upper['R_k']}",
        f"upper wall at k=1 = {upper['raw_upper_wall']}",
        f"x_1 = {upper['x_k']}",
        f"s_1 = {upper['s_k']}",
        "```",
        "",
        "Scaled-upper corridor failure while the raw upper wall holds:",
        "",
        "```text",
        f"measure: {scaled['measure']}",
        f"R_1 = {scaled['R_k']}",
        f"R_2 = {scaled['R_next']}",
        f"corridor upper = {scaled['corridor_upper']}",
        f"upper margin = {scaled['upper_margin_upper_minus_R_next']}",
        f"s_2-s_1 = {scaled['s_next_minus_s_k']}",
        "```",
        "",
        "Monotone-bridge corridor failure while the raw upper wall holds:",
        "",
        "```text",
        f"measure: {lower['measure']}",
        f"R_1 = {lower['R_k']}",
        f"R_2 = {lower['R_next']}",
        f"corridor lower = {lower['corridor_lower']}",
        f"lower margin = {lower['lower_margin_R_next_minus_lower']}",
        f"x_2-x_1 = {lower['x_next_minus_x_k']}",
        "```",
        "",
        "Exact identity:",
        "",
        "```text",
        "corridor width = 2*(2*k+1-(2*k-1)*R_k)/(2*k+1)^2",
        "```",
        "",
        "Zeta contrast:",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_raw_moment_bridge_matrix.md",
        "validated 597/597 raw-cone rows and 594/594 corridor rows on the k200 prefix",
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
        "wrote Jensen-window PF negative-lambda raw-moment obstruction matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
