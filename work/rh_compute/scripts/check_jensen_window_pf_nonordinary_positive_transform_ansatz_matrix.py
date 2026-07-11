#!/usr/bin/env python3
"""Validate the Jensen-window PF nonordinary positive-transform ansatz matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.md"

ALLOWED_VERDICTS = {
    "rejected_by_delta2_obstruction",
    "rejected_by_exactness_gap",
    "conditional_wrapper_only",
    "live_if_functional_constructed",
    "live_if_kernel_identity_proved",
    "live_if_exact_model_constructed",
    "finite_evidence_only",
    "circular_if_endpoint_used",
}

REQUIRED_IDS = {
    "npt_01_raw_power_moment_measure",
    "npt_02_positive_majorant_or_absolute_value_basis",
    "npt_03_fixed_triangular_basis_change",
    "npt_04_nonpower_positive_functional",
    "npt_05_xi_phi_positive_kernel_identity",
    "npt_06_modified_state_space_transfer_model",
    "npt_07_finite_quadrature_or_prefix_fit",
    "npt_08_endpoint_or_lp_factorization_basis",
}

LIVE_VERDICTS = {
    "live_if_functional_constructed",
    "live_if_kernel_identity_proved",
    "live_if_exact_model_constructed",
}

REQUIRED_ACCEPTANCE_TESTS = {
    "Exactness: prove the proposed readout has Taylor coefficients exactly mu_m=[t^m]1/H(-t) for every m,d,n.",
    "Nonordinary escape: explain explicitly why the raw Hankel obstruction Delta_2(mu)=-g_2 does not apply; the proof must not claim mu_m are ordinary power moments of a positive measure.",
    "Positivity: identify the positive object, such as a positive kernel, positive functional on a non-power basis, nonnegative transfer model, or Xi/Phi integral identity, and prove it gives mu_m>=0.",
    "Noncircularity: verify all hypotheses from Xi/Phi or zeta heat-flow data without endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0.",
    "Signed-data compatibility: survive lambda_n<0, beta_1<0, raw Motzkin cancellation, same-length mixed signs, and the absolute-value sign-lift gap.",
    "Uniformity: state all-degree, all-shift, and all-m coverage and isolate any limiting or compactness argument.",
    "Lift scope: state whether the result is column-only or supplies an all-Schur/Toeplitz lift.",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Nonordinary Positive Transform Ansatz Matrix",
    "Status: nonordinary positive-transform ansatz matrix",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_nonordinary_positive_transform_ansatz_matrix.py",
    "validated Jensen-window PF nonordinary positive transform ansatz matrix: 8 ansatz rows, 0 issues, 0 ready-to-apply rows, 3 live ansatz rows",
    "outputs/jensen_window_pf_positive_readout_theorem_target.md",
    "outputs/jensen_window_pf_positive_spectral_moment_obstruction.md",
    "validated Jensen-window PF positive spectral moment obstruction: 3 symbolic rows, 735 finite Delta2 obstruction rows, 0 issues",
    "outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md",
    "validated Jensen-window PF nonpower functional low-degree scout: 7 scout rows, 0 issues, 0 ready-to-apply rows, 1 live contract rows",
    "npt_01_raw_power_moment_measure",
    "npt_02_positive_majorant_or_absolute_value_basis",
    "npt_03_fixed_triangular_basis_change",
    "npt_04_nonpower_positive_functional",
    "npt_05_xi_phi_positive_kernel_identity",
    "npt_06_modified_state_space_transfer_model",
    "npt_07_finite_quadrature_or_prefix_fit",
    "npt_08_endpoint_or_lp_factorization_basis",
    "Delta_2(mu)=-g_2",
    "all-Schur/Toeplitz lift",
    "Kill Gates",
)


@dataclass(frozen=True)
class TransformAnsatzIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> TransformAnsatzIssue:
    return TransformAnsatzIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[TransformAnsatzIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(matrix: dict) -> list[TransformAnsatzIssue]:
    issues: list[TransformAnsatzIssue] = []
    if matrix.get("kind") != "jensen_window_pf_nonordinary_positive_transform_ansatz_matrix":
        issues.append(issue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    if matrix.get("parent_positive_readout_target") != "jensen_window_pf_positive_readout_theorem_target":
        issues.append(issue("<matrix>", "bad-parent-positive-readout", repr(matrix.get("parent_positive_readout_target"))))
    if set(matrix.get("parent_live_rows", [])) != {
        "pr_01_positive_spectral_transform",
        "pr_02_xi_phi_positive_resolvent_kernel",
    }:
        issues.append(issue("<matrix>", "bad-parent-live-rows", repr(matrix.get("parent_live_rows"))))
    if matrix.get("parent_obstruction") != "jensen_window_pf_positive_spectral_moment_obstruction":
        issues.append(issue("<matrix>", "bad-parent-obstruction", repr(matrix.get("parent_obstruction"))))
    if matrix.get("parent_target") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<matrix>", "bad-parent-target", repr(matrix.get("parent_target"))))
    boundary = str(matrix.get("proof_boundary", "")).lower()
    for required in ("theorem-search hygiene", "not a construction", "not an xi/phi", "not lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<matrix>", "weak-proof-boundary", required))
    objects = matrix.get("objects", {})
    for required in ("H(t)", "E(t)", "mu_m", "ordinary_moment_obstruction", "nonordinary_transform"):
        if required not in objects or not str(objects.get(required, "")).strip():
            issues.append(issue("objects", "missing-object", required))
    tests = set(matrix.get("required_acceptance_tests", []))
    for missing in sorted(REQUIRED_ACCEPTANCE_TESTS - tests):
        issues.append(issue("required_acceptance_tests", "missing-test", missing))
    return issues


def validate_row(row: dict) -> list[TransformAnsatzIssue]:
    issues: list[TransformAnsatzIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    for key in (
        "id",
        "ansatz_family",
        "verdict",
        "readiness",
        "source_artifacts",
        "fit",
        "failure_or_gap",
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
    text = " ".join(
        str(row.get(key, ""))
        for key in ("ansatz_family", "fit", "failure_or_gap", "acceptance_test", "proof_boundary")
    ).lower()
    if verdict in LIVE_VERDICTS and not any(word in text for word in ("no ", "not currently", "unknown", "construct")):
        issues.append(issue(row_id, "live-row-lacks-missing-construction", text))
    if verdict == "conditional_wrapper_only" and "wrapper" not in text:
        issues.append(issue(row_id, "wrapper-row-lacks-wrapper-language", text))
    if verdict == "circular_if_endpoint_used" and "circular" not in text:
        issues.append(issue(row_id, "circular-row-lacks-circular-language", text))
    if str(verdict).startswith("rejected") and not any(word in text for word in ("reject", "fail", "overshoot", "delta_2")):
        issues.append(issue(row_id, "rejected-row-lacks-rejection", text))
    if row_id == "npt_01_raw_power_moment_measure" and "delta_2" not in text:
        issues.append(issue(row_id, "missing-delta2-obstruction", text))
    if row_id == "npt_04_nonpower_positive_functional" and "non-power" not in text:
        issues.append(issue(row_id, "missing-nonpower-language", text))
    if row_id == "npt_05_xi_phi_positive_kernel_identity" and "xi/phi" not in text:
        issues.append(issue(row_id, "missing-xi-phi-language", text))
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "jwpf_06 is proved",
        "positive transform is proved",
        "nonordinary positive transform is proved",
    ):
        if forbidden in text:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_rows(matrix: dict) -> tuple[list[TransformAnsatzIssue], int, int, int]:
    rows = matrix.get("ansatz_rows", [])
    issues: list[TransformAnsatzIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("ansatz_rows", "missing-ansatz-rows", repr(rows))], 0, 0, 0
    seen: set[str] = set()
    ready_count = 0
    live_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("ansatz_rows", "bad-row", repr(row)))
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
        issues.append(issue(missing, "missing-ansatz-row", missing))
    return issues, len(rows), ready_count, live_count


def validate_kill_gates(matrix: dict) -> list[TransformAnsatzIssue]:
    gates = matrix.get("kill_gates", [])
    issues: list[TransformAnsatzIssue] = []
    if not isinstance(gates, list) or len(gates) < 6:
        return [issue("kill_gates", "too-few-kill-gates", repr(gates))]
    text = " ".join(str(item) for item in gates).lower()
    for required in ("ordinary power moments", "majorant", "inverse coefficient", "endpoint", "generic spectral", "all-schur"):
        if required not in text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))
    return issues


def validate_summary(matrix: dict, rows: int, ready_count: int, live_count: int) -> list[TransformAnsatzIssue]:
    issues: list[TransformAnsatzIssue] = []
    summary = matrix.get("summary", {})
    expected = {
        "ansatz_rows": 8,
        "ready_to_apply_rows": 0,
        "live_ansatz_rows": 3,
        "rejected_rows": 3,
        "conditional_wrapper_rows": 1,
        "finite_evidence_rows": 1,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if rows != 8:
        issues.append(issue("summary", "bad-row-count", str(rows)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if live_count != 3:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("delta_2=-g_2", "non-power basis", "xi/phi", "modified exact state-space", "not ready"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in matrix.get("invariants", [])).lower()
    for required in ("no ansatz row", "no finite fit", "live ansatz", "delta_2=-g_2"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[TransformAnsatzIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[TransformAnsatzIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "positive transform is proved",
        "nonordinary positive transform is proved",
        "finite fits prove",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[TransformAnsatzIssue], int, int, int]:
    matrix = load_json(matrix_path)
    issues: list[TransformAnsatzIssue] = []
    issues.extend(validate_top_level(matrix))
    row_issues, rows, ready_count, live_count = validate_rows(matrix)
    issues.extend(row_issues)
    issues.extend(validate_summary(matrix, rows, ready_count, live_count))
    issues.extend(validate_kill_gates(matrix))
    issues.extend(validate_note(note_path))
    return issues, rows, ready_count, live_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, rows, ready_count, live_count = validate(args.matrix, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "ansatz_rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "live_ansatz_rows": live_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-NONORDINARY-TRANSFORM {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF nonordinary positive transform ansatz matrix: "
            f"{rows} ansatz rows, {len(issues)} issues, {ready_count} ready-to-apply rows, "
            f"{live_count} live ansatz rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
