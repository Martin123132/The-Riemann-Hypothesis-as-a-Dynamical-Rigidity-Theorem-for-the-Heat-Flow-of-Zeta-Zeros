#!/usr/bin/env python3
"""Validate the Jensen-window PF Cauchy-Binet cone frontier matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MATRIX = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_cauchy_binet_cone_frontier_matrix.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cauchy_binet_cone_frontier_matrix.md"

ALLOWED_VERDICTS = {
    "low_degree_certificate_only",
    "rejected_as_standalone_bridge",
    "live_if_frontier_identity_constructed",
    "live_if_all_shape_kernel_constructed",
    "conditional_wrapper_only",
    "finite_evidence_only",
    "circular_if_endpoint_used",
    "language_only",
}

REQUIRED_IDS = {
    "cbcf_01_selected_low_degree_bernstein_certificate",
    "cbcf_02_adjacent_log_concavity_cone",
    "cbcf_03_column_frontier_determinant_integral",
    "cbcf_04_all_shape_cauchy_binet_kernel",
    "cbcf_05_gram_or_hankel_square_kernel",
    "cbcf_06_finite_quadrature_or_minor_fit",
    "cbcf_07_endpoint_factorization_integral",
    "cbcf_08_indefinite_andreief_or_signed_kernel",
}

LIVE_VERDICTS = {"live_if_frontier_identity_constructed", "live_if_all_shape_kernel_constructed"}

REQUIRED_FRONTIER_CONTRACT = {
    "Exact identity: produce a Cauchy-Binet, Andreief, Gram, or determinant-integral identity whose value is the actual Jensen-window Toeplitz/Jacobi-Trudi determinant or reciprocal coefficient.",
    "Positivity: prove every factor, measure, kernel, or minor in the identity is nonnegative from Xi/Phi or zeta heat-flow data.",
    "Frontier coverage: cover the first hard column-shape rows d3_column_recurrence_m8 and d4_column_recurrence_m6, not only the selected low-degree minors.",
    "All-shape scope: state whether the identity proves only column recurrence or every Schur/Toeplitz minor shape.",
    "Noncircularity: do not assume endpoint PF, Jensen hyperbolicity, Laguerre-Polya membership, RH, or Lambda <= 0.",
    "Uniformity: prove all-degree, all-shift, and all-m coverage or isolate any limiting step.",
    "Countermodel gate: explain why the exact rational log-concave countermodel is outside the new zeta-specific hypotheses without smuggling in the target.",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Cauchy-Binet Cone Frontier Matrix",
    "Status: Cauchy-Binet cone frontier matrix",
    "This is not a proof",
    "work/rh_compute/results/jensen_window_pf_cauchy_binet_cone_frontier_matrix.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_cauchy_binet_cone_frontier_matrix.py",
    "validated Jensen-window PF Cauchy-Binet cone frontier matrix: 8 frontier rows, 0 issues, 0 ready-to-apply rows, 2 live frontier rows",
    "npc_07_cauchy_binet_determinant_integral_cone",
    "15 selected low-degree formulas",
    "0 kernel identities found",
    "d3_column_recurrence_m8",
    "d4_column_recurrence_m6",
    "cbcf_03_column_frontier_determinant_integral",
    "cbcf_04_all_shape_cauchy_binet_kernel",
    "Delta_2(mu)=-g_2",
    "outputs/jensen_window_pf_nonpower_functional_cone_candidate_matrix.md",
    "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md",
    "validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
    "validated Jensen-window PF monotone contraction theorem target: 9 rows, 0 issues, 0 ready-to-apply rows, 2 live routes",
    "outputs/jensen_window_pf_monotone_contraction_stress.md",
    "validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues",
    "outputs/jensen_window_pf_cauchy_binet_low_degree_scout.md",
    "outputs/jensen_window_pf_log_concavity_frontier_scout.md",
    "Kill Gates",
)


@dataclass(frozen=True)
class CauchyBinetConeIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> CauchyBinetConeIssue:
    return CauchyBinetConeIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: str) -> list[CauchyBinetConeIssue]:
    if ref.startswith(("http://", "https://")):
        return []
    if (REPO_ROOT / ref).exists():
        return []
    return [issue(section, "missing-ref", ref)]


def validate_top_level(matrix: dict) -> list[CauchyBinetConeIssue]:
    issues: list[CauchyBinetConeIssue] = []
    if matrix.get("kind") != "jensen_window_pf_cauchy_binet_cone_frontier_matrix":
        issues.append(issue("<matrix>", "bad-kind", repr(matrix.get("kind"))))
    expected = {
        "parent_cone_matrix": "jensen_window_pf_nonpower_functional_cone_candidate_matrix",
        "parent_cone_row": "npc_07_cauchy_binet_determinant_integral_cone",
        "parent_structural_ansatz": "ansatz_02_positive_cauchy_binet_kernel",
        "parent_target": "jwpf_06_sign_regular_to_jensen_pf_conversion",
    }
    for key, value in expected.items():
        if matrix.get(key) != value:
            issues.append(issue("<matrix>", f"bad-{key}", f"{matrix.get(key)!r} != {value!r}"))
    boundary = str(matrix.get("proof_boundary", "")).lower()
    for required in ("frontier matrix", "does not construct", "does not prove", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<matrix>", "weak-proof-boundary", required))
    objects = matrix.get("objects", {})
    for required in ("low_degree_scout", "frontier_failures", "column_recurrence_frontier", "cone_contract"):
        if required not in objects or not str(objects.get(required, "")).strip():
            issues.append(issue("objects", "missing-object", required))
    contract = set(matrix.get("required_frontier_contract", []))
    for missing in sorted(REQUIRED_FRONTIER_CONTRACT - contract):
        issues.append(issue("required_frontier_contract", "missing-contract-row", missing))
    return issues


def validate_row(row: dict) -> list[CauchyBinetConeIssue]:
    issues: list[CauchyBinetConeIssue] = []
    row_id = str(row.get("id", "<missing-id>"))
    for key in (
        "id",
        "candidate",
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
        for key in ("candidate", "fit", "failure_or_gap", "acceptance_test", "proof_boundary")
    ).lower()
    if verdict in LIVE_VERDICTS:
        for required in ("no ", "identity", "construct"):
            if required not in text:
                issues.append(issue(row_id, "live-row-missing-gap-text", required))
    if verdict == "rejected_as_standalone_bridge" and "countermodel" not in text:
        issues.append(issue(row_id, "standalone-rejection-lacks-countermodel", text))
    if verdict == "conditional_wrapper_only" and "wrapper" not in text:
        issues.append(issue(row_id, "wrapper-row-lacks-wrapper-language", text))
    if verdict == "finite_evidence_only" and "finite" not in text:
        issues.append(issue(row_id, "finite-row-lacks-finite-language", text))
    if verdict == "circular_if_endpoint_used" and "circular" not in text:
        issues.append(issue(row_id, "circular-row-lacks-circular-language", text))
    if verdict == "language_only" and "does not prove" not in text:
        issues.append(issue(row_id, "language-row-lacks-boundary", text))
    if row_id == "cbcf_01_selected_low_degree_bernstein_certificate" and "kernel_identity_found=false" not in text:
        issues.append(issue(row_id, "missing-kernel-identity-failure", text))
    if row_id == "cbcf_03_column_frontier_determinant_integral" and "d3_column_recurrence_m8" not in text:
        issues.append(issue(row_id, "missing-d3-frontier", text))
    if row_id == "cbcf_04_all_shape_cauchy_binet_kernel" and "all-shape" not in text:
        issues.append(issue(row_id, "missing-all-shape-language", text))
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "kernel identity is proved",
        "cauchy-binet identity is proved",
        "determinant integral is constructed",
    ):
        if forbidden in text:
            issues.append(issue(row_id, "forbidden-overclaim", forbidden))
    return issues


def validate_rows(matrix: dict) -> tuple[list[CauchyBinetConeIssue], int, int, int]:
    rows = matrix.get("frontier_rows", [])
    issues: list[CauchyBinetConeIssue] = []
    if not isinstance(rows, list) or not rows:
        return [issue("frontier_rows", "missing-frontier-rows", repr(rows))], 0, 0, 0
    seen: set[str] = set()
    ready_count = 0
    live_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("frontier_rows", "bad-row", repr(row)))
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
        issues.append(issue(missing, "missing-frontier-row", missing))
    return issues, len(rows), ready_count, live_count


def validate_kill_gates(matrix: dict) -> list[CauchyBinetConeIssue]:
    gates = matrix.get("kill_gates", [])
    issues: list[CauchyBinetConeIssue] = []
    if not isinstance(gates, list) or len(gates) < 6:
        return [issue("kill_gates", "too-few-kill-gates", repr(gates))]
    text = " ".join(str(item) for item in gates).lower()
    for required in ("bernstein", "log-concavity", "finite", "lambda <= 0", "delta_2", "indefinite"):
        if required not in text:
            issues.append(issue("kill_gates", "missing-kill-gate-text", required))
    return issues


def validate_summary(matrix: dict, rows: int, ready_count: int, live_count: int) -> list[CauchyBinetConeIssue]:
    issues: list[CauchyBinetConeIssue] = []
    summary = matrix.get("summary", {})
    expected = {
        "frontier_rows": 8,
        "ready_to_apply_rows": 0,
        "live_frontier_rows": 2,
        "rejected_rows": 3,
        "finite_evidence_rows": 1,
        "conditional_wrapper_rows": 1,
        "language_only_rows": 1,
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
    for required in ("low-degree bernstein", "adjacent log-concavity", "column frontier", "all-shape", "neither identity"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in matrix.get("invariants", [])).lower()
    for required in ("no frontier row", "no kernel identity", "no low-degree", "live frontier", "does not prove"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[CauchyBinetConeIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[CauchyBinetConeIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "kernel identity is proved",
        "cauchy-binet identity is proved",
        "determinant integral is constructed",
        "finite rows prove",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(matrix_path: Path, note_path: Path) -> tuple[list[CauchyBinetConeIssue], int, int, int]:
    matrix = load_json(matrix_path)
    issues: list[CauchyBinetConeIssue] = []
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
                    "frontier_rows": rows,
                    "ready_to_apply_rows": ready_count,
                    "live_frontier_rows": live_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-CB-CONE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF Cauchy-Binet cone frontier matrix: "
            f"{rows} frontier rows, {len(issues)} issues, {ready_count} ready-to-apply rows, "
            f"{live_count} live frontier rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
