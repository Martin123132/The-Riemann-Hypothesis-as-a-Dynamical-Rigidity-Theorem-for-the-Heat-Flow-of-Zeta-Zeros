#!/usr/bin/env python3
"""Validate the Jensen-window PF nonpower-functional cone candidate matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_nonpower_functional_cone_candidate_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_nonpower_functional_cone_candidate_matrix.md"

ALLOWED_VERDICTS = {
    "rejected_by_signed_low_degree_polynomials",
    "rejected_as_standalone_bridge",
    "finite_evidence_only",
    "rejected_as_tautological",
    "circular_if_endpoint_used",
    "live_if_kernel_cone_constructed",
    "live_if_integral_identity_constructed",
    "conditional_wrapper_only",
}

REQUIRED_IDS = {
    "npc_01_raw_g_coordinate_cone",
    "npc_02_ratio_log_concavity_cone",
    "npc_03_finite_zeta_prefix_cone",
    "npc_04_tautological_pullback_cone",
    "npc_05_endpoint_pf_or_laguerre_polya_cone",
    "npc_06_xi_phi_kernel_test_function_cone",
    "npc_07_cauchy_binet_determinant_integral_cone",
    "npc_08_order_unit_or_sos_cone",
}

LIVE_VERDICTS = {"live_if_kernel_cone_constructed", "live_if_integral_identity_constructed"}

REQUIRED_CONE_CONTRACT = {
    "Cone definition: specify C without using the unknown target coefficients mu_m.",
    "Membership: prove K_m in C for every m,d,n from Xi/Phi or zeta heat-flow data.",
    "Positive functional: prove L is positive on C independently of the desired conclusion.",
    "Exact readout: prove L(K_m)=mu_m=[t^m]1/H(-t) for the original reciprocal coefficients.",
    "Low-degree forcing: explain how the cone certifies h_1^2-h_0*h_2 and h_1^3-2*h_0*h_1*h_2+h_0^2*h_3.",
    "Noncircularity: do not assume endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0.",
    "Uniformity: state all-m, all-degree, all-shift coverage and any limiting step.",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Nonpower Functional Cone Candidate Matrix",
    "Status: nonpower functional cone-candidate matrix",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_nonpower_functional_cone_candidate_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_nonpower_functional_cone_candidate_matrix.py",
    "validated Jensen-window PF nonpower functional cone candidate matrix: 8 cone rows, 0 issues, 0 ready-to-apply rows, 2 live cone rows",
    "npc_01_raw_g_coordinate_cone",
    "npc_06_xi_phi_kernel_test_function_cone",
    "npc_07_cauchy_binet_determinant_integral_cone",
    "outputs/jensen_window_pf_nonpower_functional_low_degree_scout.md",
    "outputs/jensen_window_pf_ratio_condition_scout.md",
    "outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md",
    "outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md",
    "validated Jensen-window PF Cauchy-Binet cone frontier matrix: 8 frontier rows, 0 issues, 0 ready-to-apply rows, 2 live frontier rows",
    "h_1^2-h_0*h_2",
    "h_1^3-2*h_0*h_1*h_2+h_0^2*h_3",
    "Kill Gates",
)


@dataclass(frozen=True)
class ConeCandidateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ConeCandidateIssue:
    return ConeCandidateIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[ConeCandidateIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(matrix: dict) -> list[ConeCandidateIssue]:
    issues: list[ConeCandidateIssue] = []
    if matrix.get("kind") != "jensen_window_pf_nonpower_functional_cone_candidate_matrix":
        issues.append(issue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    expected = {
        "parent_low_degree_scout": "jensen_window_pf_nonpower_functional_low_degree_scout",
        "parent_ansatz_matrix": "jensen_window_pf_nonordinary_positive_transform_ansatz_matrix",
        "parent_ansatz_row": "npt_04_nonpower_positive_functional",
        "parent_target": "jwpf_06_sign_regular_to_jensen_pf_conversion",
    }
    for key, value in expected.items():
        if matrix.get(key) != value:
            issues.append(issue("<matrix>", f"bad-{key}", f"{matrix.get(key)!r} != {value!r}"))
    boundary = str(matrix.get("proof_boundary", "")).lower()
    for required in ("cone-candidate", "does not construct", "does not prove reciprocal", "does not prove lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<matrix>", "weak-proof-boundary", required))
    objects = matrix.get("objects", {})
    for required in ("candidate_cone_C", "basis_K_m", "mu_m", "low_degree_forcing", "ordinary_moment_blocker"):
        if required not in objects or not str(objects.get(required, "")).strip():
            issues.append(issue("objects", "missing-object", required))
    contract = set(matrix.get("required_cone_contract", []))
    for missing in sorted(REQUIRED_CONE_CONTRACT - contract):
        issues.append(issue("required_cone_contract", "missing-contract-row", missing))
    return issues


def validate_row(row: dict) -> list[ConeCandidateIssue]:
    issues: list[ConeCandidateIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    for key in (
        "id",
        "candidate_cone",
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
        for key in ("candidate_cone", "fit", "failure_or_gap", "acceptance_test", "proof_boundary")
    ).lower()
    if verdict in LIVE_VERDICTS:
        for required in ("no ", "construct", "identity"):
            if required not in text:
                issues.append(issue(row_id, "live-row-missing-gap-text", required))
    if str(verdict).startswith("rejected") and not any(word in text for word in ("reject", "rejected", "insufficient", "nothing")):
        issues.append(issue(row_id, "rejected-row-lacks-rejection", text))
    if verdict == "circular_if_endpoint_used" and "circular" not in text:
        issues.append(issue(row_id, "circular-row-lacks-circular-language", text))
    if verdict == "conditional_wrapper_only" and "wrapper" not in text:
        issues.append(issue(row_id, "wrapper-row-lacks-wrapper-language", text))
    if row_id == "npc_01_raw_g_coordinate_cone" and "mu_2=g1^2-g2" not in text:
        issues.append(issue(row_id, "missing-mu2-obstruction", text))
    if row_id == "npc_06_xi_phi_kernel_test_function_cone" and "xi/phi" not in text:
        issues.append(issue(row_id, "missing-xi-phi-language", text))
    if row_id == "npc_07_cauchy_binet_determinant_integral_cone" and "cauchy-binet" not in text:
        issues.append(issue(row_id, "missing-cauchy-binet-language", text))
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "jwpf_06 is proved",
        "positive cone is constructed",
        "reciprocal coefficient positivity is proved",
    ):
        if forbidden in text:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_rows(matrix: dict) -> tuple[list[ConeCandidateIssue], int, int, int]:
    rows = matrix.get("cone_rows", [])
    issues: list[ConeCandidateIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("cone_rows", "missing-cone-rows", repr(rows))], 0, 0, 0
    seen: set[str] = set()
    ready_count = 0
    live_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("cone_rows", "bad-row", repr(row)))
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
        issues.append(issue(missing, "missing-cone-row", missing))
    return issues, len(rows), ready_count, live_count


def validate_kill_gates(matrix: dict) -> list[ConeCandidateIssue]:
    gates = matrix.get("kill_gates", [])
    issues: list[ConeCandidateIssue] = []
    if not isinstance(gates, list) or len(gates) < 6:
        return [issue("kill_gates", "too-few-kill-gates", repr(gates))]
    text = " ".join(str(item) for item in gates).lower()
    for required in ("g_j", "ratio/log-concavity", "finite", "tautological", "lambda <= 0", "c, k_m, l"):
        if required not in text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))
    return issues


def validate_summary(matrix: dict, rows: int, ready_count: int, live_count: int) -> list[ConeCandidateIssue]:
    issues: list[ConeCandidateIssue] = []
    summary = matrix.get("summary", {})
    expected = {
        "cone_rows": 8,
        "ready_to_apply_rows": 0,
        "live_cone_rows": 2,
        "rejected_rows": 4,
        "finite_evidence_rows": 1,
        "conditional_wrapper_rows": 1,
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
    for required in ("raw g-coordinate", "ratio/log-concavity", "xi/phi", "cauchy-binet", "not constructed"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in matrix.get("invariants", [])).lower()
    for required in ("no cone row", "no finite prefix", "live cone", "does not prove"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ConeCandidateIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ConeCandidateIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "jwpf_06 is proved",
        "positive cone is constructed",
        "reciprocal coefficient positivity is proved",
        "finite rows prove",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[ConeCandidateIssue], int, int, int]:
    matrix = load_json(matrix_path)
    issues: list[ConeCandidateIssue] = []
    issues.extend(validate_top_level(matrix))
    row_issues, rows, ready_count, live_count = validate_rows(matrix)
    issues.extend(row_issues)
    issues.extend(validate_kill_gates(matrix))
    issues.extend(validate_summary(matrix, rows, ready_count, live_count))
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
                    "cone_rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "live_cone_rows": live_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-NONPOWER-CONE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF nonpower functional cone candidate matrix: "
            f"{rows} cone rows, {len(issues)} issues, {ready_count} ready-to-apply rows, "
            f"{live_count} live cone rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
