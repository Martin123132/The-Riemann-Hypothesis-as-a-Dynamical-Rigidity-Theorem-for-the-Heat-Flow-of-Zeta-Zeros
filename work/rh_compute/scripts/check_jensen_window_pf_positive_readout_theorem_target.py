#!/usr/bin/env python3
"""Validate the Jensen-window PF positive-readout theorem target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TARGET = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_positive_readout_theorem_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_positive_readout_theorem_target.md"

ALLOWED_VERDICTS = {
    "live_if_transform_constructed",
    "live_if_kernel_identity_proved",
    "conditional_wrapper_only",
    "language_only",
    "circular_if_endpoint_used",
    "finite_evidence_only",
    "rejected_by_signed_readout",
    "rejected_by_exactness_gap",
}

REQUIRED_IDS = {
    "pr_01_positive_spectral_transform",
    "pr_02_xi_phi_positive_resolvent_kernel",
    "pr_03_abstract_stieltjes_pick_or_hamburger_wrapper",
    "pr_04_indefinite_signed_spectral_readout",
    "pr_05_endpoint_root_factorization_readout",
    "pr_06_finite_quadrature_or_interpolation_fit",
    "pr_07_raw_signed_jacobi_scalar_resolvent",
    "pr_08_absolute_value_or_majorant_readout",
}

LIVE_VERDICTS = {"live_if_transform_constructed", "live_if_kernel_identity_proved"}

REQUIRED_OBLIGATIONS = {
    "Exactness: prove the scalar readout is exactly E(t)=1/H(-t) for every degree d and shift n, not a majorant, approximation, interpolation, or finite fit.",
    "Positivity: prove the readout gives mu_m>=0 for every m by a positive measure, positive kernel, nonnegative path model, or equivalent positive functional.",
    "Noncircularity: verify all hypotheses using Xi/Phi or zeta heat-flow data without assuming endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0.",
    "Signed-obstruction compatibility: explain lambda_n<0, beta_1<0, raw Motzkin cancellation, same-length mixed signs, and the absolute-value sign-lift gap.",
    "Scope: state whether the theorem proves only column recurrence positivity or also supplies an all-Schur/Toeplitz lift.",
    "Uniformity: state the all-degree, all-shift, and all-m coverage and identify any limiting or compactness step separately.",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Positive Readout Theorem Target",
    "Status: positive-readout theorem target",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_positive_readout_theorem_target.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_positive_readout_theorem_target.py",
    "validated Jensen-window PF positive readout theorem target: 8 candidate rows, 0 issues, 0 ready-to-apply rows, 2 live foundational routes",
    "or_06_positive_spectral_measure_after_transform",
    "or_07_xi_phi_resolvent_kernel",
    "E(t)=1/H(-t)",
    "mu_m>=0",
    "positive scalar readout",
    "all-Schur/Toeplitz lift",
    "pr_01_positive_spectral_transform",
    "pr_02_xi_phi_positive_resolvent_kernel",
    "pr_03_abstract_stieltjes_pick_or_hamburger_wrapper",
    "pr_04_indefinite_signed_spectral_readout",
    "pr_05_endpoint_root_factorization_readout",
    "pr_06_finite_quadrature_or_interpolation_fit",
    "pr_07_raw_signed_jacobi_scalar_resolvent",
    "pr_08_absolute_value_or_majorant_readout",
    "outputs/jensen_window_pf_oscillatory_resolvent_fit_matrix.md",
    "validated Jensen-window PF oscillatory resolvent fit matrix: 8 fit rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_signed_j_fraction_theorem_target.md",
    "validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_reciprocal_positivity_route_matrix.md",
    "validated Jensen-window PF reciprocal positivity route matrix: 9 rows, 0 issues, 0 ready-to-apply rows",
    "outputs/jensen_window_pf_modified_signed_model_target.md",
    "validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates",
    "outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md",
    "validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",
    "outputs/jensen_window_pf_positive_spectral_moment_obstruction.md",
    "validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues",
    "outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md",
    "validated Jensen-window PF nonordinary positive transform ansatz matrix: 8 ansatz rows, 0 issues, 0 ready-to-apply rows, 3 live ansatz rows",
    "outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md",
    "validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",
    "outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md",
    "validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues",
    "Kill Gates",
)


@dataclass(frozen=True)
class PositiveReadoutIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> PositiveReadoutIssue:
    return PositiveReadoutIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[PositiveReadoutIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(target: dict) -> list[PositiveReadoutIssue]:
    issues: list[PositiveReadoutIssue] = []
    if target.get("kind") != "jensen_window_pf_positive_readout_theorem_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("parent_fit_matrix") != "jensen_window_pf_oscillatory_resolvent_fit_matrix":
        issues.append(issue("<target>", "bad-parent-fit-matrix", repr(target.get("parent_fit_matrix"))))
    if set(target.get("parent_live_rows", [])) != {
        "or_06_positive_spectral_measure_after_transform",
        "or_07_xi_phi_resolvent_kernel",
    }:
        issues.append(issue("<target>", "bad-parent-live-rows", repr(target.get("parent_live_rows"))))
    if target.get("parent_route_row") != "rp_09_signed_or_modified_continued_fraction":
        issues.append(issue("<target>", "bad-parent-route-row", repr(target.get("parent_route_row"))))
    if target.get("parent_target") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<target>", "bad-parent-target", repr(target.get("parent_target"))))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("not a construction", "not an xi/phi", "not lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    objects = target.get("objects", {})
    for required in ("H(t)", "E(t)", "mu_m", "signed_jacobi_data", "positive_readout"):
        if required not in objects or not str(objects.get(required, "")).strip():
            issues.append(issue("objects", "missing-object", required))
    obligations = set(target.get("required_theorem_obligations", []))
    missing = REQUIRED_OBLIGATIONS - obligations
    if missing:
        issues.append(issue("required_theorem_obligations", "missing-obligation", repr(sorted(missing))))
    return issues


def validate_row(row: dict) -> list[PositiveReadoutIssue]:
    issues: list[PositiveReadoutIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    for key in (
        "id",
        "candidate_family",
        "verdict",
        "readiness",
        "source_artifacts",
        "fit",
        "missing_piece",
        "acceptance_test",
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
    text = " ".join(str(row.get(key, "")) for key in ("fit", "missing_piece", "acceptance_test", "proof_boundary")).lower()
    if verdict in LIVE_VERDICTS and not any(word in text for word in ("no transform", "no kernel", "missing", "not currently")):
        issues.append(issue(row_id, "live-row-lacks-missing-piece", text))
    if verdict == "conditional_wrapper_only" and "after" not in text:
        issues.append(issue(row_id, "wrapper-lacks-dependency", text))
    if verdict == "circular_if_endpoint_used" and "circular" not in text:
        issues.append(issue(row_id, "circular-row-lacks-boundary", text))
    if str(verdict).startswith("rejected") and not any(word in text for word in ("reject", "not exact", "not a positive", "gap", "prevent")):
        issues.append(issue(row_id, "rejected-row-lacks-rejection", text))
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "jwpf_06 is proved", "positive readout is proved"):
        if forbidden in text:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_rows(target: dict) -> tuple[list[PositiveReadoutIssue], int, int, int]:
    rows = target.get("candidate_rows", [])
    issues: list[PositiveReadoutIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("candidate_rows", "missing-candidate-rows", repr(rows))], 0, 0, 0
    seen: set[str] = set()
    ready_count = 0
    live_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("candidate_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(issue(row_id, "duplicate-row", row_id))
        seen.add(row_id)
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        if row.get("verdict") in LIVE_VERDICTS:
            live_count += 1
        issues.extend(validate_row(row))
    for missing in sorted(REQUIRED_IDS - seen):
        issues.append(issue(missing, "missing-candidate-row", missing))
    return issues, len(rows), ready_count, live_count


def validate_kill_gates(target: dict) -> list[PositiveReadoutIssue]:
    gates = target.get("kill_gates", [])
    issues: list[PositiveReadoutIssue] = []
    if not isinstance(gates, list) or len(gates) < 6:
        return [issue("kill_gates", "too-few-kill-gates", repr(gates))]
    text = " ".join(str(item) for item in gates).lower()
    for required in ("majorant", "endpoint", "indefinite", "beta_1", "actual zeta", "all-schur"):
        if required not in text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))
    return issues


def validate_summary(target: dict, rows: int, ready_count: int, live_count: int) -> list[PositiveReadoutIssue]:
    issues: list[PositiveReadoutIssue] = []
    summary = target.get("summary", {})
    expected = {
        "candidate_rows": 8,
        "ready_to_apply_rows": 0,
        "live_foundational_routes": 2,
        "conditional_wrapper_rows": 1,
        "language_or_finite_rows": 2,
        "rejected_or_circular_rows": 3,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if rows != 8:
        issues.append(issue("summary", "bad-row-count", str(rows)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if live_count != 2:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("positive spectral transform", "xi/phi", "wrapper", "endpoint", "finite", "absolute-value"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("no candidate row", "live foundational", "rejected", "does not prove"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[PositiveReadoutIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[PositiveReadoutIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "positive readout theorem is proved",
        "jwpf_06 is proved",
        "finite rows prove",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[PositiveReadoutIssue], int, int, int]:
    target = load_json(target_path)
    issues: list[PositiveReadoutIssue] = []
    issues.extend(validate_top_level(target))
    row_issues, rows, ready_count, live_count = validate_rows(target)
    issues.extend(row_issues)
    issues.extend(validate_summary(target, rows, ready_count, live_count))
    issues.extend(validate_kill_gates(target))
    issues.extend(validate_note(note_path))
    return issues, rows, ready_count, live_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, rows, ready_count, live_count = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "candidate_rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "live_foundational_routes": live_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-POSITIVE-READOUT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF positive readout theorem target: "
            f"{rows} candidate rows, {len(issues)} issues, {ready_count} ready-to-apply rows, "
            f"{live_count} live foundational routes"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
