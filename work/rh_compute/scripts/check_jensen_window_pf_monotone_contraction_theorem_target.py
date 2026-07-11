#!/usr/bin/env python3
"""Validate the Jensen-window PF monotone-contraction theorem target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_TARGET = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_theorem_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_monotone_contraction_theorem_target.md"

REQUIRED_ROW_IDS = {
    "mct_01_exact_monotone_contraction_statement",
    "mct_02_factorial_normalized_moment_reduction",
    "mct_03_total_positivity_kernel_route",
    "mct_04_heat_flow_differential_inequality_route",
    "mct_05_plain_stieltjes_moment_logconvexity",
    "mct_06_endpoint_pf_or_laguerre_polya_route",
    "mct_07_generic_ratio_conditions",
    "mct_08_finite_stress_support",
    "mct_09_column_frontier_application",
}

ALLOWED_ROLES = {
    "open_statement",
    "exact_reduction",
    "live_route",
    "insufficient_route",
    "circular_route",
    "rejected_route",
    "finite_evidence",
    "conditional_application",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Monotone-Contraction Theorem Target",
    "Status: open theorem target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_monotone_contraction_theorem_target`",
    "work/rh_compute/results/jensen_window_pf_monotone_contraction_theorem_target.json",
    "python work/rh_compute/scripts/check_jensen_window_pf_monotone_contraction_theorem_target.py",
    "validated Jensen-window PF monotone contraction theorem target: 9 rows, 0 issues, 0 ready-to-apply rows, 2 live routes",
    "A_k(lambda) = mu_{2k}(lambda) * k! / (2k)!",
    "A_{k+2}*A_k^3 >= A_{k+1}^3*A_{k-1}",
    "Delta^3 log A_{k-1} >= 0",
    "((2*k-1)*(2*k+3))/(2*k+1)^2",
    "raw moment log-convexity alone is insufficient",
    "outputs/jensen_window_pf_monotone_contraction_stress.md",
    "validated Jensen-window PF monotone contraction stress: 2875 rows, 2875 positive rows, 0 issues",
    "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md",
    "validated Jensen-window PF monotone contraction frontier scout: 2 exact rows, 88 Bernstein coefficients, 210 finite zeta rows, 0 issues",
    "outputs/jensen_window_pf_monotone_contraction_all_m_counterexample.md",
    "validated Jensen-window PF monotone-contraction all-m counterexample: degree 7, m=11, exact full-cone witness, 6 lower walls, negative normalized value, 0 issues",
    "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
    "validated Jensen-window PF heat-flow cone-entry asymptotic target: 8 rows, 0 issues, 0 ready-to-apply rows, 1 live routes",
    "outputs/jensen_window_pf_heat_flow_ratio_cone_invariance_lemma.md",
    "validated Jensen-window PF heat-flow ratio cone invariance lemma: 6 exact rows, 315 lower rows, 315 upper rows, 310 monotone rows, 0 issues",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "validated Jensen-window PF heat-flow boundary threshold lemma: 5 exact rows, 315 strong-threshold rows, 315 heat-threshold rows, 0 issues",
    "outputs/jensen_window_pf_heat_flow_monotone_closure_scout.md",
    "validated Jensen-window PF heat-flow monotone closure scout: 4 exact rows, 315 threshold rows, 305 flow-bracket rows, 0 issues",
)


@dataclass(frozen=True)
class TargetIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> TargetIssue:
    return TargetIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[TargetIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(target: dict) -> list[TargetIssue]:
    issues: list[TargetIssue] = []
    if target.get("kind") != "jensen_window_pf_monotone_contraction_theorem_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") != "open_theorem_target":
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    if target.get("target_id") != "target_monotone_contraction_zeta_coefficients":
        issues.append(issue("<target>", "bad-target-id", repr(target.get("target_id"))))
    if target.get("parent_route") != "jensen_window_pf_cauchy_binet_cone_frontier_matrix":
        issues.append(issue("<target>", "bad-parent-route", repr(target.get("parent_route"))))
    if target.get("parent_frontier_rows") != ["d3_column_recurrence_m8", "d4_column_recurrence_m6"]:
        issues.append(issue("<target>", "bad-frontier-rows", repr(target.get("parent_frontier_rows"))))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("open theorem target", "does not prove", "all-m", "cauchy-binet", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    return issues


def validate_normalization(target: dict) -> list[TargetIssue]:
    normalization = target.get("normalization", {})
    issues: list[TargetIssue] = []
    expected = {
        "A_k_lambda": "mu_2k(lambda) * k! / (2k)!",
        "monotone_contraction": "x_{k+1} >= x_k",
        "log_difference_form": "Delta^2 log A_{k-1} <= 0 and Delta^3 log A_{k-1} >= 0",
        "polynomial_form": "A_{k+2}*A_k^3 >= A_{k+1}^3*A_{k-1}",
        "raw_moment_adjacent_factor": "x_k(A) = ((2*k-1)/(2*k+1)) * (mu_{2k+2}*mu_{2k-2}/mu_{2k}^2)",
        "raw_moment_monotone_equivalent": "mu_{2k+4}*mu_{2k}^3/(mu_{2k+2}^3*mu_{2k-2}) >= ((2*k-1)*(2*k+3))/(2*k+1)^2",
    }
    for key, value in expected.items():
        if normalization.get(key) != value:
            issues.append(issue("normalization", f"bad-{key}", repr(normalization.get(key))))
    return issues


def validate_contract(target: dict) -> list[TargetIssue]:
    rows = target.get("theorem_contract", [])
    issues: list[TargetIssue] = []
    if not isinstance(rows, list) or len(rows) < 7:
        return [issue("theorem_contract", "too-few-contract-rows", repr(rows))]
    text = " ".join(str(row) for row in rows).lower()
    for required in ("a_k", "x_k", "factorial", "noncircular", "countermodel", "column", "lambda <= 0"):
        if required not in text:
            issues.append(issue("theorem_contract", "missing-contract-text", required))
    return issues


def validate_target_rows(target: dict) -> tuple[list[TargetIssue], int, int, int]:
    rows = target.get("target_rows", [])
    issues: list[TargetIssue] = []
    if not isinstance(rows, list):
        return [issue("target_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))
    ready_count = 0
    live_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("target_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "source_artifacts", "claim_if_proved", "gap", "acceptance_test", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("role") == "live_route":
            live_count += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        text = " ".join(str(row.get(key, "")) for key in ("gap", "acceptance_test", "proof_boundary")).lower()
        if row.get("role") in {"finite_evidence", "exact_reduction", "conditional_application"} and "not" not in text:
            issues.append(issue(row_id, "weak-boundary", text))
        if row.get("role") == "circular_route" and "circular" not in text:
            issues.append(issue(row_id, "missing-circular-language", text))
        if row.get("role") == "rejected_route" and "countermodel" not in text:
            issues.append(issue(row_id, "missing-countermodel-language", text))
    return issues, len(rows), ready_count, live_count


def validate_summary(target: dict, row_count: int, ready_count: int, live_count: int) -> list[TargetIssue]:
    summary = target.get("summary", {})
    issues: list[TargetIssue] = []
    expected = {
        "target_rows": 9,
        "ready_to_apply_rows": 0,
        "live_routes": 2,
        "exact_reduction_rows": 1,
        "rejected_or_insufficient_rows": 3,
        "finite_evidence_rows": 1,
        "conditional_application_rows": 1,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if live_count != 2:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("open theorem target", "delta^3", "factorial", "finite stress", "all-m", "no analytic theorem"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("no row", "open_target", "finite stress", "all-m", "endpoint", "factorial"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[TargetIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[TargetIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "analytic monotone-contraction theorem is proved",
        "monotone contractions are proved for all zeta windows",
        "cauchy-binet identity is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[TargetIssue], int, int, int]:
    target = load_json(target_path)
    issues: list[TargetIssue] = []
    issues.extend(validate_top_level(target))
    issues.extend(validate_normalization(target))
    issues.extend(validate_contract(target))
    row_issues, row_count, ready_count, live_count = validate_target_rows(target)
    issues.extend(row_issues)
    issues.extend(validate_summary(target, row_count, ready_count, live_count))
    issues.extend(validate_note(note_path))
    return issues, row_count, ready_count, live_count


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, row_count, ready_count, live_count = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(
            json.dumps(
                {
                    "ok": ok,
                    "rows": row_count,
                    "ready_to_apply_rows": ready_count,
                    "live_routes": live_count,
                    "issues": [asdict(item) for item in issues],
                },
                indent=2,
                sort_keys=True,
            )
        )
    else:
        for item in issues:
            print(f"JWPF-MONOTONE-TARGET {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF monotone contraction theorem target: "
            f"{row_count} rows, {len(issues)} issues, {ready_count} ready-to-apply rows, "
            f"{live_count} live routes"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
