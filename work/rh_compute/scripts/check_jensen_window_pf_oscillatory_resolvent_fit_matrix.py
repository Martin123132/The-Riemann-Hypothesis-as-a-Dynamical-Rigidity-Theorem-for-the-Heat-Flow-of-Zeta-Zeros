#!/usr/bin/env python3
"""Validate the Jensen-window PF oscillatory resolvent fit matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_oscillatory_resolvent_fit_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md"

ALLOWED_VERDICTS = {
    "rejected_by_entry_signs",
    "rejected_by_sign_invariance",
    "rejected_by_mu2_gap",
    "route_mismatch",
    "language_only",
    "live_if_transform_proved",
    "live_if_theorem_proved",
    "finite_evidence_only",
}

REQUIRED_ROW_IDS = {
    "or_01_nonnegative_jacobi_or_production_resolvent",
    "or_02_diagonal_similarity_to_nonnegative_matrix",
    "or_03_absolute_value_resolvent_majorant",
    "or_04_classical_oscillatory_matrix_spectrum",
    "or_05_indefinite_or_krein_space_moment_model",
    "or_06_positive_spectral_measure_after_transform",
    "or_07_xi_phi_resolvent_kernel",
    "or_08_finite_signed_jacobi_pattern",
}

LIVE_VERDICTS = {"live_if_transform_proved", "live_if_theorem_proved"}

REQUIRED_CONDITIONS = {
    "Construct an exact finite or operator model whose matrix element, resolvent, or moment generating function is E(t)=1/H(-t) for every degree and shift.",
    "Prove coefficientwise nonnegativity of mu_m for every m, not only spectral reality, oscillation, sign-regularity, or formal continued-fraction existence.",
    "Verify all hypotheses for the actual zeta heat-flow Jensen windows without assuming endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0.",
    "Explain how the theorem survives lambda_n<0, beta_1<0, raw Motzkin cancellation, same-length mixed signs, and the absolute-value sign-lift obstruction.",
    "State whether the theorem proves only the column recurrence or also gives an all-Schur/Toeplitz lift.",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Oscillatory Resolvent Fit Matrix",
    "Status: oscillatory/resolvent theorem-search matrix",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_oscillatory_resolvent_fit_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_oscillatory_resolvent_fit_matrix.py",
    "validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows",
    "or_01_nonnegative_jacobi_or_production_resolvent",
    "or_02_diagonal_similarity_to_nonnegative_matrix",
    "or_03_absolute_value_resolvent_majorant",
    "or_04_classical_oscillatory_matrix_spectrum",
    "or_05_indefinite_or_krein_space_moment_model",
    "or_06_positive_spectral_measure_after_transform",
    "or_07_xi_phi_resolvent_kernel",
    "or_08_finite_signed_jacobi_pattern",
    "mu_m >= 0",
    "lambda_n<0",
    "beta_1<0",
    "positive scalar readout",
    "all-Schur/Toeplitz lift",
    "outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md",
    "validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md",
    "validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",
    "outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md",
    "validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues",
    "outputs/sign_regularity_theorem_map.md",
    "outputs/sign_regularity_theorem_fit_matrix.md",
    "outputs/jensen_window_pf_theorem_machinery_fit_matrix.md",
    "outputs/jensen_window_pf_signed_j_fraction_theorem_target.md",
    "outputs/jensen_window_pf_modified_signed_model_target.md",
    "outputs/jensen_window_pf_positive_readout_theorem_target.md",
    "validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes",
    "outputs/jensen_window_pf_positive_spectral_moment_obstruction.md",
    "validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues",
    "Kill Gates",
)


@dataclass(frozen=True)
class ResolventIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ResolventIssue:
    return ResolventIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[ResolventIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(matrix: dict) -> list[ResolventIssue]:
    issues: list[ResolventIssue] = []
    if matrix.get("kind") != "jensen_window_pf_oscillatory_resolvent_fit_matrix":
        issues.append(issue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    if matrix.get("parent_model_row") != "msm_06_oscillatory_sign_regular_resolvent":
        issues.append(issue("<matrix>", "bad-parent-model-row", repr(matrix.get("parent_model_row"))))
    if matrix.get("parent_route_row") != "rp_09_signed_or_modified_continued_fraction":
        issues.append(issue("<matrix>", "bad-parent-route-row", repr(matrix.get("parent_route_row"))))
    if matrix.get("parent_target") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<matrix>", "bad-parent-target", repr(matrix.get("parent_target"))))
    boundary = str(matrix.get("proof_boundary", "")).lower()
    for required in ("not an oscillatory", "not a positive spectral", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<matrix>", "weak-proof-boundary", required))
    obj = matrix.get("object", {})
    for required in ("E(t)", "mu_m", "formal_jacobi_model", "signed_parameters", "desired_resolvent_output"):
        if required not in obj or not str(obj.get(required, "")).strip():
            issues.append(issue("object", "missing-object", required))
    conditions = set(matrix.get("required_conditions_for_use", []))
    missing = REQUIRED_CONDITIONS - conditions
    if missing:
        issues.append(issue("required_conditions_for_use", "missing-condition", repr(sorted(missing))))
    return issues


def validate_row(row: dict) -> list[ResolventIssue]:
    issues: list[ResolventIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    for key in (
        "id",
        "theorem_family",
        "verdict",
        "readiness",
        "source_artifacts",
        "fit",
        "gap_or_rejection",
        "next_action",
        "proof_boundary",
    ):
        if key not in row:
            issues.append(issue(row_id, "missing-field", key))
    verdict = row.get("verdict")
    if verdict not in ALLOWED_VERDICTS:
        issues.append(issue(row_id, "bad-verdict", repr(verdict)))
    if row.get("readiness") != "not_ready_to_apply":
        issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
    refs = row.get("source_artifacts", [])
    if not isinstance(refs, list) or not refs:
        issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
    else:
        for ref in refs:
            if isinstance(ref, str):
                issues.extend(validate_ref(row_id, ref))
            else:
                issues.append(issue(row_id, "bad-ref", repr(ref)))
    text = " ".join(str(row.get(key, "")) for key in ("fit", "gap_or_rejection", "next_action", "proof_boundary")).lower()
    if verdict in LIVE_VERDICTS:
        if not any(word in text for word in ("no transform", "no xi", "missing", "not currently", "search")):
            issues.append(issue(row_id, "live-row-lacks-gap", text))
    if str(verdict).startswith("rejected") and not any(word in text for word in ("reject", "not", "cannot", "gap", "invariant")):
        issues.append(issue(row_id, "rejected-row-lacks-rejection", text))
    if verdict == "route_mismatch" and "not by themselves" not in text:
        issues.append(issue(row_id, "mismatch-lacks-boundary", text))
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "jwpf_06 is proved", "resolvent theorem is proved"):
        if forbidden in text:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_rows(matrix: dict) -> tuple[list[ResolventIssue], int, int]:
    rows = matrix.get("fit_rows", [])
    issues: list[ResolventIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("fit_rows", "missing-fit-rows", repr(rows))], 0, 0
    seen: set[str] = set()
    ready_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("fit_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(issue(row_id, "duplicate-row", row_id))
        seen.add(row_id)
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        issues.extend(validate_row(row))
    for missing in sorted(REQUIRED_ROW_IDS - seen):
        issues.append(issue(missing, "missing-fit-row", missing))
    return issues, len(rows), ready_count


def validate_kill_gates(matrix: dict) -> list[ResolventIssue]:
    kill_gates = matrix.get("kill_gates", [])
    issues: list[ResolventIssue] = []
    if not isinstance(kill_gates, list) or len(kill_gates) < 6:
        return [issue("kill_gates", "too-few-kill-gates", repr(kill_gates))]
    text = " ".join(str(item) for item in kill_gates).lower()
    for required in ("formal j-fraction", "positive scalar readout", "diagonal", "absolute", "endpoint", "all-schur"):
        if required not in text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))
    return issues


def validate_summary(matrix: dict, rows: int, ready_count: int) -> list[ResolventIssue]:
    issues: list[ResolventIssue] = []
    summary = matrix.get("summary", {})
    expected = {
        "fit_rows": 8,
        "ready_to_apply_rows": 0,
        "rejected_or_mismatch_rows": 4,
        "language_or_finite_rows": 2,
        "live_conditional_rows": 2,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if rows != 8:
        issues.append(issue("summary", "bad-row-count", str(rows)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("positive spectral measure", "xi/phi", "entrywise", "oscillatory", "coefficientwise"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in matrix.get("invariants", [])).lower()
    for required in ("no fit row", "live row", "rejected row", "does not prove"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ResolventIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ResolventIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "oscillatory theorem is proved",
        "jwpf_06 is proved",
        "finite signed jacobi rows prove",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[ResolventIssue], int, int]:
    matrix = load_json(matrix_path)
    issues: list[ResolventIssue] = []
    issues.extend(validate_top_level(matrix))
    row_issues, rows, ready_count = validate_rows(matrix)
    issues.extend(row_issues)
    issues.extend(validate_kill_gates(matrix))
    issues.extend(validate_summary(matrix, rows, ready_count))
    issues.extend(validate_note(note_path))
    return issues, rows, ready_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, rows, ready_count = validate(args.matrix, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "fit_rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"OSCILLATORY-RESOLVENT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF oscillatory resolvent fit matrix: "
            f"{rows} fit rows, {len(issues)} issues, {ready_count} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
