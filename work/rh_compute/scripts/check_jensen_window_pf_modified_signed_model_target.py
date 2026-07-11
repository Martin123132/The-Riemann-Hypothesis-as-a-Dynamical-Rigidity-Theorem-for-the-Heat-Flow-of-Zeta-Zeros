#!/usr/bin/env python3
"""Validate the Jensen-window PF modified signed-model target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TARGET = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_modified_signed_model_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_modified_signed_model_target.md"

ALLOWED_VERDICTS = {
    "rejected_by_motzkin_obstruction",
    "rejected_by_sign_repair_obstruction",
    "live_if_constructed",
    "live_if_theorem_proved",
    "live_if_identity_proved",
    "finite_evidence_only",
    "circular_if_used_as_input",
}

REQUIRED_MODEL_IDS = {
    "msm_01_raw_ordinary_motzkin_j_fraction",
    "msm_02_diagonal_sign_conjugation",
    "msm_03_global_length_parity_sign",
    "msm_04_state_space_doubled_model",
    "msm_05_modified_production_matrix_or_riordan",
    "msm_06_oscillatory_sign_regular_resolvent",
    "msm_07_cancellation_positive_identity",
    "msm_08_finite_grid_extrapolation",
    "msm_09_endpoint_real_rooted_or_pf_input",
}

LIVE_VERDICTS = {
    "live_if_constructed",
    "live_if_theorem_proved",
    "live_if_identity_proved",
}

REQUIRED_CONDITIONS = {
    "Recover exactly E(t)=1/H(-t) for all degrees d and shifts n of the actual zeta heat-flow Jensen windows.",
    "Prove coefficientwise positivity mu_m>=0 for every m, not only a finite prefix.",
    "Use no endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0 assumption as input.",
    "State whether the result proves only column recurrence shapes or also supplies an all-Schur/Toeplitz lift.",
    "Survive the raw Motzkin path obstruction, beta_1 diagonal obstruction, same-length mixed-sign path obstruction, and finite-prefix promotion gates.",
    "If the model uses cancellation, end with a manifestly nonnegative formula rather than a difference of positive quantities.",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Modified Signed-Model Target",
    "Status: modified signed-model theorem-search audit",
    "This is not a proof of",
    "work/rh_compute/results/jensen_window_pf_modified_signed_model_target.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_modified_signed_model_target.py",
    "validated Jensen-window PF modified signed-model target: 9 model rows, 0 issues, 0 ready-to-apply rows, 4 live modified candidates",
    "msm_01_raw_ordinary_motzkin_j_fraction",
    "msm_02_diagonal_sign_conjugation",
    "msm_03_global_length_parity_sign",
    "msm_04_state_space_doubled_model",
    "msm_05_modified_production_matrix_or_riordan",
    "msm_06_oscillatory_sign_regular_resolvent",
    "msm_07_cancellation_positive_identity",
    "msm_08_finite_grid_extrapolation",
    "msm_09_endpoint_real_rooted_or_pf_input",
    "E(t) = 1 / H(-t)",
    "mu_m >= 0",
    "raw Motzkin path obstruction",
    "same-length mixed-sign path obstruction",
    "all-Schur/Toeplitz lift",
    "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md",
    "validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md",
    "validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues",
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
    "finite grids",
    "endpoint real-rootedness",
    "not merely a difference of positive sums",
)


@dataclass(frozen=True)
class ModelTargetIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ModelTargetIssue:
    return ModelTargetIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[ModelTargetIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(target: dict) -> list[ModelTargetIssue]:
    issues: list[ModelTargetIssue] = []
    if target.get("kind") != "jensen_window_pf_modified_signed_model_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("parent_route_row") != "rp_09_signed_or_modified_continued_fraction":
        issues.append(issue("<target>", "bad-parent-route-row", repr(target.get("parent_route_row"))))
    if target.get("parent_target") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<target>", "bad-parent-target", repr(target.get("parent_target"))))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("not a construction", "not an all-order", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    objects = target.get("objects", {})
    for required in ("H(t)", "E(t)", "mu_m", "raw_motzkin_model", "modified_model"):
        if required not in objects or not str(objects.get(required, "")).strip():
            issues.append(issue("objects", "missing-object", required))
    conditions = set(target.get("required_conditions_for_live_model", []))
    missing = REQUIRED_CONDITIONS - conditions
    if missing:
        issues.append(issue("required_conditions", "missing-condition", repr(sorted(missing))))
    return issues


def validate_row(row: dict) -> list[ModelTargetIssue]:
    issues: list[ModelTargetIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    for key in (
        "id",
        "model_family",
        "verdict",
        "readiness",
        "live_candidate",
        "source_artifacts",
        "fit",
        "gap_or_rejection",
        "required_upgrade",
        "proof_boundary",
    ):
        if key not in row:
            issues.append(issue(row_id, "missing-field", key))
    verdict = row.get("verdict")
    if verdict not in ALLOWED_VERDICTS:
        issues.append(issue(row_id, "bad-verdict", repr(verdict)))
    if row.get("readiness") != "not_ready_to_apply":
        issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
    should_be_live = verdict in LIVE_VERDICTS
    if row.get("live_candidate") is not should_be_live:
        issues.append(issue(row_id, "bad-live-candidate-flag", repr(row.get("live_candidate"))))
    refs = row.get("source_artifacts")
    if not isinstance(refs, list) or not refs:
        issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
    else:
        for ref in refs:
            if not isinstance(ref, str):
                issues.append(issue(row_id, "bad-ref", repr(ref)))
            else:
                issues.extend(validate_ref(row_id, ref))
    text = " ".join(str(row.get(key, "")) for key in ("fit", "gap_or_rejection", "required_upgrade", "proof_boundary")).lower()
    if should_be_live:
        if not any(word in text for word in ("construct", "produce", "theorem", "identity", "missing", "no exact", "not currently")):
            issues.append(issue(row_id, "live-row-lacks-missing-upgrade", text))
    if verdict.startswith("rejected") and not any(word in text for word in ("obstruction", "reject", "cannot", "invariant")):
        issues.append(issue(row_id, "rejected-row-lacks-obstruction", text))
    if row_id == "msm_08_finite_grid_extrapolation" and "finite" not in text:
        issues.append(issue(row_id, "finite-row-lacks-finite-boundary", text))
    if row_id == "msm_09_endpoint_real_rooted_or_pf_input" and "circular" not in text:
        issues.append(issue(row_id, "endpoint-row-lacks-circular-boundary", text))
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "jwpf_06 is proved"):
        if forbidden in text:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_rows(target: dict) -> tuple[list[ModelTargetIssue], int, int, int]:
    rows = target.get("model_rows", [])
    issues: list[ModelTargetIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("model_rows", "missing-model-rows", repr(rows))], 0, 0, 0
    seen: set[str] = set()
    ready_count = 0
    live_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("model_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(issue(row_id, "duplicate-row", row_id))
        seen.add(row_id)
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        if row.get("live_candidate") is True:
            live_count += 1
        issues.extend(validate_row(row))
    for missing in sorted(REQUIRED_MODEL_IDS - seen):
        issues.append(issue(missing, "missing-model-row", missing))
    return issues, len(rows), ready_count, live_count


def validate_summary(target: dict, rows: int, ready_count: int, live_count: int) -> list[ModelTargetIssue]:
    issues: list[ModelTargetIssue] = []
    summary = target.get("summary", {})
    expected = {
        "model_rows": 9,
        "ready_to_apply_rows": 0,
        "live_modified_candidates": 4,
        "rejected_or_dead_repairs": 3,
        "finite_or_circular_rows": 2,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if rows != 9:
        issues.append(issue("summary", "bad-row-count", str(rows)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if live_count != 4:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("raw ordinary motzkin", "length-parity", "four modified routes", "finite grids", "endpoint pf"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("no model row", "live candidate", "finite evidence", "does not prove"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_kill_gates(target: dict) -> list[ModelTargetIssue]:
    kill_gates = target.get("kill_gates", [])
    issues: list[ModelTargetIssue] = []
    if not isinstance(kill_gates, list) or len(kill_gates) < 7:
        return [issue("kill_gates", "too-few-kill-gates", repr(kill_gates))]
    text = " ".join(str(item) for item in kill_gates).lower()
    for required in (
        "raw ordinary motzkin",
        "diagonal sign conjugation",
        "path-length parity",
        "absolute-value sign-state cover",
        "finite",
        "endpoint",
        "formal signed-fraction",
        "all-schur",
    ):
        if required not in text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))
    return issues


def validate_note(path: Path) -> list[ModelTargetIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ModelTargetIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "modified signed model proves",
        "jwpf_06 is proved",
        "finite grids prove",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[ModelTargetIssue], int, int, int]:
    target = load_json(target_path)
    issues: list[ModelTargetIssue] = []
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
                    "model_rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "live_modified_candidates": live_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-MODIFIED-SIGNED-MODEL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF modified signed-model target: "
            f"{rows} model rows, {len(issues)} issues, {ready_count} ready-to-apply rows, "
            f"{live_count} live modified candidates"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
