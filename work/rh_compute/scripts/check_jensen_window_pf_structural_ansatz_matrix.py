#!/usr/bin/env python3
"""Validate the Jensen-window PF structural ansatz matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json"
DEFAULT_ALGEBRA = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_obligation_algebra.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_structural_ansatz_matrix.md"


ALLOWED_STATUSES = {
    "live_structural_candidate",
    "blocked_until_input_pf",
    "rejected_by_countermodel",
}

REQUIRED_HARD_TEST_IDS = {
    "hard_01_degree2_contact",
    "hard_02_degree3_2x2_off_hankel",
    "hard_03_degree3_3x3_cubic",
    "hard_04_degree4_banded_coefficients",
    "hard_05_countermodel_kill_gate",
}

REQUIRED_FORMULAS_BY_TEST = {
    "hard_01_degree2_contact": {
        "4*(a1**2 - a0*a2)",
        "a1**2 - a0*a2 = -det([[a0,a1],[a1,a2]])",
    },
    "hard_02_degree3_2x2_off_hankel": {
        "-3*(a0*a2 - 3*a1**2)",
        "-a0*a3 + 9*a1*a2",
        "-3*(a1*a3 - 3*a2**2)",
    },
    "hard_03_degree3_3x3_cubic": {
        "a0**2*a3 - 18*a0*a1*a2 + 27*a1**3",
        "a0*a3**2 - 18*a1*a2*a3 + 27*a2**3",
    },
    "hard_04_degree4_banded_coefficients": {
        "-2*(3*a0*a2 - 8*a1**2)",
        "-4*(a0*a3 - 6*a1*a2)",
        "-a0*a4 + 16*a1*a3",
        "-2*(3*a0*a2*a4 - 8*a0*a3**2 - 8*a1**2*a4 + 96*a1*a2*a3 - 108*a2**3)",
    },
}

REQUIRED_COUNTERMODEL_FIELDS = {
    "d3_first_negative_contiguous_toeplitz_minor.size=8",
    "d4_first_negative_contiguous_toeplitz_minor.size=6",
    "blocked_promotion includes selected low-order Jensen-window Toeplitz minors",
}

REQUIRED_ANSATZ_IDS = {
    "ansatz_01_direct_shifted_hankel_only_transfer",
    "ansatz_02_positive_cauchy_binet_kernel",
    "ansatz_03_planar_network_or_production_matrix",
    "ansatz_04_binomial_multiplier_preserver",
    "ansatz_05_positive_determinant_integral",
    "ansatz_06_finite_grid_interpolation",
}

LIVE_CANDIDATES = {
    "ansatz_02_positive_cauchy_binet_kernel",
    "ansatz_03_planar_network_or_production_matrix",
    "ansatz_05_positive_determinant_integral",
}

REJECTED_ROWS = {
    "ansatz_01_direct_shifted_hankel_only_transfer",
    "ansatz_06_finite_grid_interpolation",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Structural Ansatz Matrix",
    "Status: structural proof-search matrix",
    "This is not a proof of Jensen-window",
    "jwpf_06_sign_regular_to_jensen_pf_conversion",
    "work/rh_compute/results/jensen_window_pf_structural_ansatz_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_structural_ansatz_matrix.py",
    "validated Jensen-window PF structural ansatz matrix: 6 ansatz rows, 0 issues, 0 ready-to-apply rows",
    "hard_01_degree2_contact",
    "hard_02_degree3_2x2_off_hankel",
    "hard_03_degree3_3x3_cubic",
    "hard_04_degree4_banded_coefficients",
    "hard_05_countermodel_kill_gate",
    "ansatz_02_positive_cauchy_binet_kernel",
    "ansatz_03_planar_network_or_production_matrix",
    "ansatz_05_positive_determinant_integral",
    "No row is `ready_to_apply`",
    "Schur Shape Contract",
    "outputs/jensen_window_pf_schur_shape_contract.md",
    "outputs/jensen_window_pf_column_recurrence_contract.md",
    "work/rh_compute/results/jensen_window_pf_schur_shape_contract.json",
    "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json",
    "python work/rh_compute/scripts/jensen_window_pf_schur_shape_contract.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_schur_shape_contract.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_column_recurrence_contract.py",
    "validated Jensen-window PF Schur shape contract: 4 grid rows, 0 issues, 2 frontier rows",
    "validated Jensen-window PF column recurrence contract: 4 degree rows, 0 issues, 2 hard frontier rows",
    "15,709",
    "column-shape Jacobi-Trudi",
    "C_m = h_0^m * e_m",
    "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin path obstruction scout: 3 symbolic rows, 735 mu2 cancellation rows, 630 beta1 diagonal obstruction rows, 0 issues",
    "outputs/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.md",
    "validated Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: 3 symbolic rows, 5145 mixed-sign witness rows, 0 issues",
    "outputs/jensen_window_pf_state_space_sign_lift_obstruction_scout.md",
    "validated Jensen-window PF state-space sign-lift obstruction scout: 3 symbolic rows, 735 mu2 sign-lift obstruction rows, 0 issues",
    "raw ordinary Motzkin/J-fraction",
    "Cauchy-Binet Low-Degree Scout",
    "outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md",
    "work/rh_compute/results/jensen_window_pf_cauchy_binet_low_degree_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_low_degree_scout.py",
    "validated Jensen-window PF Cauchy-Binet low-degree scout: 15 formula rows, 0 issues, 0 kernel identities found",
    "outputs/jensen_window_pf_log_concavity_frontier_scout.md",
    "work/rh_compute/results/jensen_window_pf_log_concavity_frontier_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_log_concavity_frontier_scout.py",
    "validated Jensen-window PF log-concavity frontier scout: 14 contiguous rows, 0 issues",
    "first Bernstein",
    "countermodel negatives",
    "outputs/jensen_window_pf_ratio_condition_scout.md",
    "outputs/jensen_window_pf_contraction_log_concavity_scout.md",
    "work/rh_compute/results/jensen_window_pf_ratio_condition_scout.json",
    "work/rh_compute/results/jensen_window_pf_contraction_log_concavity_scout.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_ratio_condition_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_contraction_log_concavity_scout.py",
    "validated Jensen-window PF ratio-condition scout: 7 candidate rows, 0 issues, 4 rejected by countermodel, 1 rejected by construction",
    "validated Jensen-window PF contraction-log-concavity scout: 1 rejected by construction, 0 issues, 2 negative frontier rows",
    "d3_first_negative_contiguous_toeplitz_minor.size=8",
    "d4_first_negative_contiguous_toeplitz_minor.size=6",
)


@dataclass(frozen=True)
class AnsatzIssue:
    row_id: str
    issue: str
    detail: str


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def issue(row_id: str, name: str, detail: str) -> AnsatzIssue:
    return AnsatzIssue(row_id=row_id, issue=name, detail=detail)


def algebra_formula_pool(algebra: dict) -> set[str]:
    formulas = {
        str(algebra.get("degree2", {}).get("jensen_discriminant", "")),
        str(algebra.get("degree2", {}).get("signed_hankel_relation", "")),
    }
    for section in ("degree3", "degree4"):
        for row in algebra.get(section, {}).get("selected_toeplitz_minors", []):
            if isinstance(row, dict):
                formulas.add(str(row.get("determinant", "")))
    return {formula for formula in formulas if formula}


def validate_ref(row_id: str, ref: str) -> list[AnsatzIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(row_id, "missing-ref", ref)]


def has_boundary(text: str) -> bool:
    lowered = text.lower()
    return any(
        marker in lowered
        for marker in ("not ", "only", "missing", "rejected", "conditional", "open")
    )


def validate_top_level(matrix: dict) -> list[AnsatzIssue]:
    issues: list[AnsatzIssue] = []
    if matrix.get("kind") != "jensen_window_pf_structural_ansatz_matrix":
        issues.append(issue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    if matrix.get("target_obligation") != "jwpf_06_sign_regular_to_jensen_pf_conversion":
        issues.append(issue("<matrix>", "bad-target-obligation", repr(matrix.get("target_obligation"))))
    boundary = str(matrix.get("proof_boundary", "")).lower()
    if "not a proof" not in boundary or "lambda <= 0" not in boundary:
        issues.append(issue("<matrix>", "weak-proof-boundary", matrix.get("proof_boundary", "")))
    statuses = set(matrix.get("status_legend", {}))
    missing_statuses = ALLOWED_STATUSES - statuses
    if missing_statuses:
        issues.append(issue("<matrix>", "missing-status-legend", repr(sorted(missing_statuses))))
    return issues


def validate_hard_tests(matrix: dict, algebra: dict) -> list[AnsatzIssue]:
    issues: list[AnsatzIssue] = []
    tests = matrix.get("hard_low_degree_tests", [])
    if not isinstance(tests, list) or not tests:
        return [issue("hard_tests", "missing-hard-tests", "hard_low_degree_tests must be a nonempty list")]

    by_id: dict[str, dict] = {}
    formula_pool = algebra_formula_pool(algebra)
    for row in tests:
        if not isinstance(row, dict):
            issues.append(issue("hard_tests", "bad-hard-test", repr(row)))
            continue
        test_id = str(row.get("id", ""))
        if not test_id:
            issues.append(issue("hard_tests", "missing-id", repr(row)))
            continue
        if test_id in by_id:
            issues.append(issue(test_id, "duplicate-hard-test", test_id))
        by_id[test_id] = row
        if not str(row.get("source", "")).startswith("work/rh_compute/results/jensen_window_pf_obligation_algebra.json"):
            issues.append(issue(test_id, "bad-source", str(row.get("source", ""))))
        if not str(row.get("implication", "")).strip():
            issues.append(issue(test_id, "missing-implication", test_id))

    for missing in sorted(REQUIRED_HARD_TEST_IDS - set(by_id)):
        issues.append(issue(missing, "missing-required-hard-test", missing))

    for test_id, required_formulas in REQUIRED_FORMULAS_BY_TEST.items():
        row = by_id.get(test_id, {})
        formulas = set(row.get("required_formulas", []))
        missing = required_formulas - formulas
        if missing:
            issues.append(issue(test_id, "missing-required-formulas", repr(sorted(missing))))
        for formula in formulas:
            if formula not in formula_pool:
                issues.append(issue(test_id, "formula-not-in-algebra", formula))

    counter_row = by_id.get("hard_05_countermodel_kill_gate", {})
    counter_fields = set(counter_row.get("required_countermodel_fields", []))
    missing_counter_fields = REQUIRED_COUNTERMODEL_FIELDS - counter_fields
    if missing_counter_fields:
        issues.append(issue("hard_05_countermodel_kill_gate", "missing-countermodel-fields", repr(sorted(missing_counter_fields))))
    counter = algebra.get("finite_countermodel", {})
    if counter.get("d3_first_negative_contiguous_toeplitz_minor", {}).get("size") != 8:
        issues.append(issue("hard_05_countermodel_kill_gate", "bad-algebra-d3-countermodel-size", repr(counter.get("d3_first_negative_contiguous_toeplitz_minor"))))
    if counter.get("d4_first_negative_contiguous_toeplitz_minor", {}).get("size") != 6:
        issues.append(issue("hard_05_countermodel_kill_gate", "bad-algebra-d4-countermodel-size", repr(counter.get("d4_first_negative_contiguous_toeplitz_minor"))))
    if "selected low-order Jensen-window Toeplitz minors" not in str(counter.get("blocked_promotion", "")):
        issues.append(issue("hard_05_countermodel_kill_gate", "missing-blocked-promotion-text", str(counter.get("blocked_promotion", ""))))
    return issues


def validate_ansatz_row(row: dict) -> list[AnsatzIssue]:
    issues: list[AnsatzIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    required_fields = (
        "id",
        "ansatz_family",
        "status",
        "readiness",
        "hoped_mechanism",
        "required_identity",
        "hard_tests_required",
        "failed_hard_tests",
        "low_degree_fit",
        "current_failure",
        "next_action",
        "refs",
        "hypotheses_verified",
        "closes_jwpf06_if_proved",
        "proof_boundary",
    )
    for key in required_fields:
        if key not in row:
            issues.append(issue(row_id, "missing-field", key))

    status = row.get("status")
    if status not in ALLOWED_STATUSES:
        issues.append(issue(row_id, "bad-status", repr(status)))
    if row.get("readiness") == "ready_to_apply":
        issues.append(issue(row_id, "ready-to-apply-present", row_id))
    if row.get("hypotheses_verified") is not False:
        issues.append(issue(row_id, "hypotheses-should-not-be-verified", repr(row.get("hypotheses_verified"))))
    if not isinstance(row.get("closes_jwpf06_if_proved"), bool):
        issues.append(issue(row_id, "bad-closes-flag", repr(row.get("closes_jwpf06_if_proved"))))

    hard_tests = set(row.get("hard_tests_required", []))
    if "hard_05_countermodel_kill_gate" not in hard_tests:
        issues.append(issue(row_id, "missing-countermodel-kill-gate", repr(sorted(hard_tests))))
    if row_id in LIVE_CANDIDATES:
        if row.get("status") != "live_structural_candidate":
            issues.append(issue(row_id, "live-candidate-has-wrong-status", repr(row.get("status"))))
        if row.get("closes_jwpf06_if_proved") is not True:
            issues.append(issue(row_id, "live-candidate-not-marked-conditional-closing", repr(row.get("closes_jwpf06_if_proved"))))
        missing_tests = REQUIRED_HARD_TEST_IDS - hard_tests
        if missing_tests:
            issues.append(issue(row_id, "live-candidate-missing-hard-tests", repr(sorted(missing_tests))))
        combined = f"{row.get('current_failure', '')} {row.get('next_action', '')}".lower()
        if "missing" not in combined and "no " not in combined:
            issues.append(issue(row_id, "missing-gap-language", combined))
    if row_id in REJECTED_ROWS:
        if row.get("status") != "rejected_by_countermodel":
            issues.append(issue(row_id, "rejected-row-has-wrong-status", repr(row.get("status"))))
        if row.get("closes_jwpf06_if_proved") is True:
            issues.append(issue(row_id, "rejected-row-closes-target", row_id))
        if not row.get("failed_hard_tests"):
            issues.append(issue(row_id, "rejected-row-missing-failed-tests", row_id))

    refs = row.get("refs")
    if not isinstance(refs, list) or not refs:
        issues.append(issue(row_id, "missing-refs", "refs must be a nonempty list"))
    else:
        for ref in refs:
            if not isinstance(ref, str):
                issues.append(issue(row_id, "bad-ref", repr(ref)))
            else:
                issues.extend(validate_ref(row_id, ref))

    for key in ("hoped_mechanism", "required_identity", "current_failure", "next_action"):
        if not str(row.get(key, "")).strip():
            issues.append(issue(row_id, f"missing-{key}", row_id))
    if not has_boundary(str(row.get("proof_boundary", ""))):
        issues.append(issue(row_id, "weak-proof-boundary", str(row.get("proof_boundary", ""))))
    combined = " ".join(str(row.get(key, "")) for key in ("hoped_mechanism", "required_identity", "current_failure", "proof_boundary")).lower()
    for forbidden in ("therefore rh", "proves lambda <= 0", "bridge is proved"):
        if forbidden in combined:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_ansatz_rows(matrix: dict) -> tuple[list[AnsatzIssue], int, int]:
    issues: list[AnsatzIssue] = []
    rows = matrix.get("ansatz_rows", [])
    if not isinstance(rows, list) or not rows:
        return [issue("ansatz_rows", "missing-rows", "ansatz_rows must be a nonempty list")], 0, 0

    seen: set[str] = set()
    ready_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("ansatz_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        if row_id in seen:
            issues.append(issue(row_id, "duplicate-id", row_id))
        seen.add(row_id)
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        issues.extend(validate_ansatz_row(row))

    for missing in sorted(REQUIRED_ANSATZ_IDS - seen):
        issues.append(issue(missing, "missing-required-ansatz-row", missing))
    return issues, len(rows), ready_count


def validate_note(path: Path) -> list[AnsatzIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[AnsatzIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in ("therefore rh", "we have proved lambda <= 0", "the bridge is proved"):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, algebra_path: Path, note_path: Path) -> tuple[list[AnsatzIssue], int, int]:
    matrix = load_json(matrix_path)
    algebra = load_json(algebra_path)
    issues: list[AnsatzIssue] = []
    issues.extend(validate_top_level(matrix))
    issues.extend(validate_hard_tests(matrix, algebra))
    row_issues, row_count, ready_count = validate_ansatz_rows(matrix)
    issues.extend(row_issues)
    issues.extend(validate_note(note_path))
    return issues, row_count, ready_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX)
    parser.add_argument("--algebra", type=Path, default=DEFAULT_ALGEBRA)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, row_count, ready_count = validate(args.matrix, args.algebra, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "ansatz_rows": row_count,
                    "ready_to_apply_rows": ready_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-STRUCTURAL-ANSATZ {item.row_id} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF structural ansatz matrix: "
            f"{row_count} ansatz rows, {len(issues)} issues, {ready_count} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
