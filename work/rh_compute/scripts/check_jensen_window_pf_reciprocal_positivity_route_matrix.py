#!/usr/bin/env python3
"""Validate the Jensen-window PF reciprocal-positivity route matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md"


ALLOWED_VERDICTS = {
    "endpoint_equivalence_only",
    "live_if_representation_proved",
    "blocked_until_input_pf",
    "route_mismatch",
    "rejected_by_countermodel",
    "finite_evidence_only",
    "circular_if_used_as_input",
}

REQUIRED_IDS = {
    "rp_01_real_negative_zero_or_pf_window_endpoint",
    "rp_02_positive_renewal_or_convolution_kernel",
    "rp_03_positive_stieltjes_or_j_fraction",
    "rp_04_companion_or_production_matrix_total_positivity",
    "rp_05_kaluza_reciprocal_sign_theorem",
    "rp_06_generic_ratio_or_log_concavity_conditions",
    "rp_07_finite_recurrence_stress_grid",
    "rp_08_full_schur_or_toeplitz_lift",
    "rp_09_signed_or_modified_continued_fraction",
}

REQUIRED_FEATURES = {
    "proves coefficient nonnegativity of 1/H(-t) for all m",
    "handles H(t)=sum_j binom(d,j) A_{n+j} t^j / A_n for all d,n",
    "does not assume Jensen-window PF-infinity, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0",
    "does not promote finite recurrence stress rows into an all-order theorem",
    "states whether it proves only column shapes or could lift to all Schur/Toeplitz shapes",
}

LIVE_IDS = {
    "rp_02_positive_renewal_or_convolution_kernel",
    "rp_04_companion_or_production_matrix_total_positivity",
    "rp_09_signed_or_modified_continued_fraction",
}

REJECTED_IDS = {
    "rp_06_generic_ratio_or_log_concavity_conditions",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Reciprocal Positivity Route Matrix",
    "Status: reciprocal-positivity theorem-search matrix",
    "This is not a proof of Schur positivity",
    "E(t) = 1 / H(-t)",
    "work/rh_compute/results/jensen_window_pf_reciprocal_positivity_route_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_reciprocal_positivity_route_matrix.py",
    "validated Jensen-window PF reciprocal positivity route matrix: 9 rows, 0 issues, 0 ready-to-apply rows",
    "rp_01_real_negative_zero_or_pf_window_endpoint",
    "rp_02_positive_renewal_or_convolution_kernel",
    "rp_03_positive_stieltjes_or_j_fraction",
    "rp_04_companion_or_production_matrix_total_positivity",
    "rp_05_kaluza_reciprocal_sign_theorem",
    "rp_06_generic_ratio_or_log_concavity_conditions",
    "rp_07_finite_recurrence_stress_grid",
    "rp_08_full_schur_or_toeplitz_lift",
    "rp_09_signed_or_modified_continued_fraction",
    "12,600",
    "72,600",
    "outputs/jensen_window_pf_reciprocal_fraction_scout.md",
    "validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md",
    "validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md",
    "validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",
    "outputs/jensen_window_pf_signed_j_fraction_theorem_target.md",
    "validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_modified_signed_model_target.md",
    "validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates",
    "outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md",
    "validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues",
    "outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md",
    "validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_positive_readout_theorem_target.md",
    "validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes",
    "outputs/jensen_window_pf_positive_spectral_moment_obstruction.md",
    "validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_coefficient_extended_stress.md",
    "validated Jensen-window PF reciprocal coefficient extended stress: 72600 rows, 0 issues",
    "No row is `ready_to_apply`",
    "no row closes `jwpf_06`",
    "outputs/jensen_window_pf_column_recurrence_contract.md",
    "outputs/jensen_window_pf_column_recurrence_finite_coverage.md",
    "outputs/arb_jensen_window_column_recurrence_stress.md",
    "outputs/jensen_window_pf_ratio_condition_scout.md",
    "outputs/jensen_window_pf_contraction_log_concavity_scout.md",
)


@dataclass(frozen=True)
class RouteIssue:
    row_id: str
    issue: str
    detail: str


def issue(row_id: str, name: str, detail: str) -> RouteIssue:
    return RouteIssue(row_id=row_id, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(row_id: str, ref: str) -> list[RouteIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(row_id, "missing-ref", ref)]


def has_boundary(text: str) -> bool:
    lowered = text.lower()
    return any(
        marker in lowered
        for marker in ("not ", "only", "missing", "rejected", "finite", "candidate", "blocked")
    )


def validate_top_level(matrix: dict) -> list[RouteIssue]:
    issues: list[RouteIssue] = []
    if matrix.get("kind") != "jensen_window_pf_reciprocal_positivity_route_matrix":
        issues.append(issue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    if matrix.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<matrix>", "bad-target-obligation", repr(matrix.get("target_obligation"))))
    target = str(matrix.get("column_recurrence_target", ""))
    for required in ("1/H(-t)", "binom(d,j)", "A_{n+j}", "every m,d,n"):
        if required not in target:
            issues.append(issue("<matrix>", "bad-column-target", required))
    boundary = str(matrix.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<matrix>", "weak-proof-boundary", matrix.get("proof_boundary", "")))
    allowed = set(matrix.get("allowed_verdicts", []))
    if not ALLOWED_VERDICTS.issubset(allowed):
        issues.append(issue("<matrix>", "missing-allowed-verdicts", repr(sorted(ALLOWED_VERDICTS - allowed))))
    features = set(matrix.get("required_route_features", []))
    if not REQUIRED_FEATURES.issubset(features):
        issues.append(issue("<matrix>", "missing-required-features", repr(sorted(REQUIRED_FEATURES - features))))
    return issues


def validate_row(row: dict) -> list[RouteIssue]:
    issues: list[RouteIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    required_fields = (
        "id",
        "route_family",
        "verdict",
        "readiness",
        "hypotheses_verified",
        "would_prove_column_recurrence_if_verified",
        "would_close_jwpf06_if_verified",
        "source_artifacts",
        "mechanism",
        "missing_hypothesis",
        "countermodel_or_gate",
        "next_action",
        "scope",
        "proof_boundary",
    )
    for key in required_fields:
        if key not in row:
            issues.append(issue(row_id, "missing-field", key))

    verdict = row.get("verdict")
    if verdict not in ALLOWED_VERDICTS:
        issues.append(issue(row_id, "bad-verdict", repr(verdict)))
    if row.get("readiness") != "not_ready_to_apply":
        issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
    if row.get("hypotheses_verified") is not False:
        issues.append(issue(row_id, "hypotheses-verified", repr(row.get("hypotheses_verified"))))
    if row.get("would_close_jwpf06_if_verified") is not False:
        issues.append(issue(row_id, "target-closing-row", repr(row.get("would_close_jwpf06_if_verified"))))
    if row_id in LIVE_IDS and verdict != "live_if_representation_proved":
        issues.append(issue(row_id, "live-row-wrong-verdict", repr(verdict)))
    if row_id in LIVE_IDS and "missing" not in str(row.get("missing_hypothesis", "")).lower():
        issues.append(issue(row_id, "live-row-missing-gap", str(row.get("missing_hypothesis", ""))))
    if row_id in REJECTED_IDS and verdict != "rejected_by_countermodel":
        issues.append(issue(row_id, "rejected-row-wrong-verdict", repr(verdict)))

    refs = row.get("source_artifacts")
    if not isinstance(refs, list) or not refs:
        issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
    else:
        for ref in refs:
            if not isinstance(ref, str):
                issues.append(issue(row_id, "bad-ref", repr(ref)))
            else:
                issues.extend(validate_ref(row_id, ref))

    for key in ("mechanism", "missing_hypothesis", "countermodel_or_gate", "next_action", "scope"):
        if not str(row.get(key, "")).strip():
            issues.append(issue(row_id, f"empty-{key}", key))
    if not has_boundary(str(row.get("proof_boundary", ""))):
        issues.append(issue(row_id, "weak-proof-boundary", str(row.get("proof_boundary", ""))))

    combined = " ".join(str(row.get(key, "")) for key in row).lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "jwpf_06 is proved"):
        if forbidden in combined:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_note(path: Path) -> list[RouteIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[RouteIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "the bridge is proved", "jwpf_06 is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[RouteIssue], int, int]:
    matrix = load_json(matrix_path)
    issues: list[RouteIssue] = []
    issues.extend(validate_top_level(matrix))

    rows = matrix.get("rows", [])
    if not isinstance(rows, list) or not rows:
        issues.append(issue("<matrix>", "missing-rows", repr(rows)))
        rows = []

    seen: set[str] = set()
    ready_count = 0
    live_count = 0
    target_closing_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("<matrix>", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(issue(row_id, "duplicate-id", row_id))
        seen.add(row_id)
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        if row.get("verdict") == "live_if_representation_proved":
            live_count += 1
        if row.get("would_close_jwpf06_if_verified") is True:
            target_closing_count += 1
        issues.extend(validate_row(row))

    for missing in sorted(REQUIRED_IDS - seen):
        issues.append(issue(missing, "missing-required-row", missing))

    summary = matrix.get("summary", {})
    expected = {
        "rows": 9,
        "ready_to_apply_rows": 0,
        "live_representation_candidates": 3,
        "target_closing_rows": 0,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("<summary>", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if ready_count != 0:
        issues.append(issue("<matrix>", "ready-to-apply-row-present", str(ready_count)))
    if live_count != 3:
        issues.append(issue("<matrix>", "bad-live-count", str(live_count)))
    if target_closing_count != 0:
        issues.append(issue("<matrix>", "target-closing-present", str(target_closing_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("reciprocal", "circular", "ordinary positive", "signed j-fraction", "finite grids", "generic ratio"):
        if required not in finding:
            issues.append(issue("<summary>", "missing-main-finding-text", required))

    issues.extend(validate_note(note_path))
    return issues, len(rows), ready_count


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
                    "rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-RECIPROCAL-ROUTE {item.row_id} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF reciprocal positivity route matrix: "
            f"{rows} rows, {len(issues)} issues, {ready_count} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
