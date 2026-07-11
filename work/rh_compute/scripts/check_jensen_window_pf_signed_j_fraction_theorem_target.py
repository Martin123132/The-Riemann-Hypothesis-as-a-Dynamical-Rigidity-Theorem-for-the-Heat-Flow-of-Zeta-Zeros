#!/usr/bin/env python3
"""Validate the Jensen-window PF signed J-fraction theorem target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TARGET = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_signed_j_fraction_theorem_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_signed_j_fraction_theorem_target.md"

ALLOWED_VERDICTS = {
    "rejected_by_fraction_sign",
    "open_missing_theorem",
    "live_if_conjugation_theorem_proved",
    "live_if_positive_model_proved",
    "language_only",
    "finite_evidence_only",
    "circular_endpoint_only",
}

REQUIRED_FIT_IDS = {
    "sj_01_standard_positive_stieltjes_j_fraction",
    "sj_02_signed_j_fraction_signature",
    "sj_03_oscillatory_jacobi_matrix_after_sign_conjugation",
    "sj_04_production_matrix_or_lattice_path_model",
    "sj_05_indefinite_moment_problem",
    "sj_06_finite_signature_extrapolation",
    "sj_07_endpoint_real_rooted_model",
}

REQUIRED_OBJECTS = {
    "h_j",
    "g_j",
    "H(t)",
    "E(t)",
    "mu_m",
    "Delta_r",
    "Delta_r^*",
    "Q_r",
    "beta_n",
    "lambda_n",
    "kappa_n",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Signed J-Fraction Theorem Target",
    "Status: signed J-fraction theorem target specification",
    "This is not a proof of",
    "work/rh_compute/results/jensen_window_pf_signed_j_fraction_theorem_target.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_signed_j_fraction_theorem_target.py",
    "validated Jensen-window PF signed J-fraction theorem target: 7 fit rows, 0 issues, 0 ready-to-apply rows",
    "(-1)^(r(r-1)/2) Delta_r > 0",
    "lambda_n = Delta_{n+1} Delta_{n-1} / Delta_n^2",
    "kappa_n = -lambda_n",
    "beta_n = Q_{n+1} - Q_n",
    "outputs/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.md",
    "validated Jensen-window PF reciprocal signed Jacobi beta scout: 3 symbolic rows, 3675 beta rows, 2940 positive rows, 630 negative rows, 105 terminal-zero rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",
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
    "raw ordinary Motzkin/J-fraction",
    "mu_m >= 0 for every m, degree d, and shift n",
    "sj_01_standard_positive_stieltjes_j_fraction",
    "sj_02_signed_j_fraction_signature",
    "sj_03_oscillatory_jacobi_matrix_after_sign_conjugation",
    "sj_04_production_matrix_or_lattice_path_model",
    "sj_05_indefinite_moment_problem",
    "sj_06_finite_signature_extrapolation",
    "sj_07_endpoint_real_rooted_model",
    "outputs/jensen_window_pf_reciprocal_signed_j_fraction_scout.md",
    "validated Jensen-window PF reciprocal signed J-fraction scout: 2 symbolic rows, 3675 signed Hankel rows, 2940 signed-lambda rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_fraction_scout.md",
    "validated Jensen-window PF reciprocal fraction scout: 3 symbolic rows, 735 finite rows, 0 issues",
    "https://epubs.siam.org/doi/10.1137/090781127",
    "https://arxiv.org/pdf/2202.03793",
    "outputs/sign_regularity_theorem_map.md",
    "outputs/sign_regularity_theorem_fit_matrix.md",
    "Acceptance Conditions",
    "Kill Gates",
    "all-Schur/Toeplitz lift",
)


@dataclass(frozen=True)
class SignedJTargetIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> SignedJTargetIssue:
    return SignedJTargetIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[SignedJTargetIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_target(target: dict) -> list[SignedJTargetIssue]:
    issues: list[SignedJTargetIssue] = []
    if target.get("kind") != "jensen_window_pf_signed_j_fraction_theorem_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("target_route_row") != "rp_09_signed_or_modified_continued_fraction":
        issues.append(issue("<target>", "bad-target-route-row", repr(target.get("target_route_row"))))
    if target.get("parent_target") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<target>", "bad-parent-target", repr(target.get("parent_target"))))
    boundary = str(target.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<target>", "weak-proof-boundary", str(target.get("proof_boundary", ""))))

    objects = set(target.get("objects", {}))
    missing_objects = REQUIRED_OBJECTS - objects
    if missing_objects:
        issues.append(issue("objects", "missing-objects", repr(sorted(missing_objects))))
    for key in REQUIRED_OBJECTS:
        if key in target.get("objects", {}) and not str(target["objects"][key]).strip():
            issues.append(issue("objects", f"empty-{key}", key))

    missing_theorem = target.get("missing_theorem", {})
    if missing_theorem.get("id") != "signed_j_fraction_column_recurrence_theorem":
        issues.append(issue("missing_theorem", "bad-id", repr(missing_theorem.get("id"))))
    for key in ("hypotheses", "desired_conclusion", "acceptance_conditions"):
        value = missing_theorem.get(key)
        if not isinstance(value, list) or not value:
            issues.append(issue("missing_theorem", f"missing-{key}", repr(value)))
    theorem_text = " ".join(
        " ".join(missing_theorem.get(key, []))
        for key in ("hypotheses", "desired_conclusion", "acceptance_conditions")
    ).lower()
    for required in ("actual", "signed", "mu_m", "not assume", "all-degree", "all-shift"):
        if required not in theorem_text:
            issues.append(issue("missing_theorem", "missing-theorem-text", required))

    rows = target.get("fit_rows", [])
    if not isinstance(rows, list) or not rows:
        issues.append(issue("fit_rows", "missing-fit-rows", repr(rows)))
        rows = []
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
        if row.get("verdict") not in ALLOWED_VERDICTS:
            issues.append(issue(row_id, "bad-verdict", repr(row.get("verdict"))))
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        for key in ("theorem_family", "fit", "gap_or_rejection", "next_action"):
            if not str(row.get(key, "")).strip():
                issues.append(issue(row_id, f"missing-{key}", key))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                if isinstance(ref, str):
                    issues.extend(validate_ref(row_id, ref))
                else:
                    issues.append(issue(row_id, "bad-ref", repr(ref)))
    for missing in sorted(REQUIRED_FIT_IDS - seen):
        issues.append(issue(missing, "missing-fit-row", missing))
    if ready_count:
        issues.append(issue("fit_rows", "ready-to-apply-row-present", str(ready_count)))

    anchors = target.get("source_anchors", [])
    if not isinstance(anchors, list) or len(anchors) < 3:
        issues.append(issue("source_anchors", "too-few-anchors", repr(anchors)))
    for anchor in anchors:
        if not isinstance(anchor, dict):
            issues.append(issue("source_anchors", "bad-anchor", repr(anchor)))
            continue
        url = anchor.get("url")
        if not isinstance(url, str) or not url:
            issues.append(issue("source_anchors", "missing-url", repr(anchor)))
        else:
            issues.extend(validate_ref(str(anchor.get("id", "anchor")), url))

    kill_gates = target.get("kill_gates", [])
    if not isinstance(kill_gates, list) or len(kill_gates) < 5:
        issues.append(issue("kill_gates", "too-few-kill-gates", repr(kill_gates)))
    kill_text = " ".join(str(item) for item in kill_gates).lower()
    for required in ("finite", "ordinary positive", "endpoint", "formal existence", "all-schur"):
        if required not in kill_text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))

    summary = target.get("summary", {})
    expected_summary = {
        "fit_rows": 7,
        "open_missing_theorem_rows": 1,
        "live_conditional_rows": 2,
        "rejected_or_finite_rows": 3,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    main_finding = str(summary.get("main_finding", "")).lower()
    for required in ("signed j-fraction", "missing theorem", "coefficientwise", "finite extrapolation"):
        if required not in main_finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_note(path: Path) -> list[SignedJTargetIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[SignedJTargetIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "signed j-fraction theorem is proved",
        "jwpf_06 is proved",
        "finite evidence proves",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[SignedJTargetIssue], int, int]:
    target = load_json(target_path)
    issues = validate_target(target)
    issues.extend(validate_note(note_path))
    fit_rows = target.get("summary", {}).get("fit_rows", 0)
    ready_rows = target.get("summary", {}).get("ready_to_apply_rows", 0)
    return issues, int(fit_rows), int(ready_rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, fit_rows, ready_rows = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "fit_rows": fit_rows,
                    "ready_to_apply_rows": ready_rows,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-SIGNED-J-TARGET {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF signed J-fraction theorem target: "
            f"{fit_rows} fit rows, {len(issues)} issues, {ready_rows} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
