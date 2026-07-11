#!/usr/bin/env python3
"""Validate the relative-Gaussian intervalization target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target import (  # noqa: E402
    DEFAULT_GRID_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgit_01_budget_import",
    "nlrgit_02_common_error_cap",
    "nlrgit_03_interval_sources",
    "nlrgit_04_grid_to_collar_gap",
    "nlrgit_05_floating_grid_promotion_rejected",
    "nlrgit_06_acceptance_gate",
}

REQUIRED_OBLIGATION_IDS = {
    "nlrgit_01_residual_core_identity",
    "nlrgit_02_laguerre_node_weight_intervals",
    "nlrgit_03_phi_and_c0_interval_tail",
    "nlrgit_04_quadrature_remainder_error",
    "nlrgit_05_ratio_and_coefficient_ball_propagation",
    "nlrgit_06_rounding_and_aggregation_budget",
    "nlrgit_07_finite_grid_not_uniform_collar",
    "nlrgit_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "exact_sufficient_condition",
    "open_requirement",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Intervalization Target",
    "Status: open numerical-certification target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian intervalization target: 6 rows, 0 issues, 8 obligations, 5 open requirements, 0 ready-to-apply rows",
    "observed worst value ratio: 0.970710059020329541",
    "observed worst derivative ratio: 0.969356777475803869",
    "proposed total ratio error cap: 1.000000000000000000E-2",
    "closed if total cap met: True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class IntervalizationIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> IntervalizationIssue:
    return IntervalizationIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[IntervalizationIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[IntervalizationIssue]:
    issues: list[IntervalizationIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "open numerical-certification target":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_cancellation_reduced_grid_scout",
        "source_actual_endpoint_scout",
        "source_asymptotic_remainder_target",
        "source_formal_tail_obstruction",
        "source_residual_budget",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "open numerical-certification target",
        "sufficient error budget",
        "does not provide interval",
        "does not prove",
        "uniform collar",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict, grid_path: Path) -> list[IntervalizationIssue]:
    try:
        recomputed = build_artifact(grid_path)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[IntervalizationIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[IntervalizationIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {})
    issues: list[IntervalizationIssue] = []
    source = diagnostics.get("source_grid_summary", {})
    expected_source = {
        "grid_rows": 20,
        "t_grid_count": 5,
        "index_count": 4,
        "quadrature_order_count": 4,
        "selected_quadrature_order": 192,
        "max_value_ratio_to_first_omitted_over_orders": "0.970710059020329541",
        "max_derivative_ratio_to_first_omitted_over_orders": "0.969356777475803869",
        "max_value_ratio_spread_over_orders": "9.63500159925082018e-15",
        "max_derivative_ratio_spread_over_orders": "9.62202755450426673e-15",
    }
    for key, value in expected_source.items():
        if source.get(key) != value:
            issues.append(issue("source_grid_summary", f"bad-{key}", f"{source.get(key)!r} != {value!r}"))
    summary = diagnostics.get("budget_summary", {})
    expected_budget = {
        "value_slack_to_one": "2.928994097967045900E-2",
        "derivative_slack_to_one": "3.064322252419613100E-2",
        "common_slack_to_one": "2.928994097967045900E-2",
        "common_half_slack": "1.464497048983522950E-2",
        "proposed_total_ratio_error_cap": "1.000000000000000000E-2",
        "proposed_per_error_source_cap_for_five_sources": "2.000000000000000000E-3",
        "closed_if_total_cap_met": True,
    }
    for key, value in expected_budget.items():
        if summary.get(key) != value:
            issues.append(issue("budget_summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    budget_rows = diagnostics.get("certification_budget_rows", [])
    if diagnostics.get("certification_budget_row_count") != 2:
        issues.append(issue("diagnostics", "bad-budget-row-count", repr(diagnostics.get("certification_budget_row_count"))))
    if [row.get("channel") for row in budget_rows] != ["value", "derivative"]:
        issues.append(issue("budget_rows", "bad-channel-order", repr([row.get("channel") for row in budget_rows])))
    for row in budget_rows:
        if row.get("closed_if_cap_met") is not True:
            issues.append(issue(str(row.get("channel")), "cap-does-not-close", repr(row)))
        if "not yet proved" not in str(row.get("proof_boundary", "")).lower():
            issues.append(issue(str(row.get("channel")), "weak-budget-boundary", repr(row.get("proof_boundary"))))
    obligations = diagnostics.get("obligation_rows", [])
    if diagnostics.get("obligation_row_count") != 8:
        issues.append(issue("diagnostics", "bad-obligation-count", repr(diagnostics.get("obligation_row_count"))))
    if diagnostics.get("open_requirement_rows") != 5:
        issues.append(issue("diagnostics", "bad-open-requirement-count", repr(diagnostics.get("open_requirement_rows"))))
    if diagnostics.get("ready_to_apply_rows") != 0:
        issues.append(issue("diagnostics", "unexpected-ready-to-apply", repr(diagnostics.get("ready_to_apply_rows"))))
    if {row.get("id") for row in obligations} != REQUIRED_OBLIGATION_IDS:
        issues.append(issue("obligations", "bad-obligation-ids", repr([row.get("id") for row in obligations])))
    for row in obligations:
        row_id = str(row.get("id"))
        if row.get("readiness") == "ready_to_apply":
            issues.append(issue(row_id, "unexpected-ready-to-apply", row_id))
        if row.get("role") == "open_requirement" and row.get("target_error_cap") is None:
            issues.append(issue(row_id, "missing-target-error-cap", row_id))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "rejected")):
            issues.append(issue(row_id, "weak-obligation-boundary", boundary))
    note = str(diagnostics.get("proof_boundary_note", "")).lower()
    for required in ("certification roadmap", "does not", "node/weight", "quadrature", "phi", "continuum-in-t"):
        if required not in note:
            issues.append(issue("diagnostics", "weak-proof-boundary-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[IntervalizationIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[IntervalizationIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_rows = 0
    sufficient_rows = 0
    open_rows = 0
    rejected_rows = 0
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        role = row.get("role")
        if role not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(role)))
        if role == "exact_reduction":
            exact_rows += 1
        if role == "exact_sufficient_condition":
            sufficient_rows += 1
        if role == "open_requirement":
            open_rows += 1
        if role == "rejected_route":
            rejected_rows += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") not in {"available_exact", "not_ready_to_apply"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), exact_rows, sufficient_rows, open_rows, rejected_rows


def validate_summary(
    artifact: dict,
    row_count: int,
    exact_rows: int,
    sufficient_rows: int,
    open_rows: int,
    rejected_rows: int,
) -> list[IntervalizationIssue]:
    summary = artifact.get("summary", {})
    issues: list[IntervalizationIssue] = []
    expected = {
        "matrix_rows": 6,
        "certification_budget_rows": 2,
        "obligation_rows": 8,
        "open_requirement_rows": 5,
        "ready_to_apply_rows": 0,
        "observed_worst_value_ratio": "0.970710059020329541",
        "observed_worst_derivative_ratio": "0.969356777475803869",
        "common_slack_to_one": "2.928994097967045900E-2",
        "proposed_total_ratio_error_cap": "1.000000000000000000E-2",
        "proposed_per_error_source_cap_for_five_sources": "2.000000000000000000E-3",
        "closed_if_total_cap_met": True,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 6:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_rows != 1:
        issues.append(issue("summary", "bad-exact-row-count", str(exact_rows)))
    if sufficient_rows != 1:
        issues.append(issue("summary", "bad-sufficient-row-count", str(sufficient_rows)))
    if open_rows != 2:
        issues.append(issue("summary", "bad-open-row-count", str(open_rows)))
    if rejected_rows != 1:
        issues.append(issue("summary", "bad-rejected-row-count", str(rejected_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "common ratio slack",
        "2.928994097967e-2",
        "total ratio error",
        "1.0e-2",
        "laguerre nodes",
        "phi tails",
        "quadrature error",
        "full collar",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "not promoted", "finite-grid", "sufficient target", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[IntervalizationIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[IntervalizationIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "uniform residual estimate is proved",
        "residual-tail estimates are proved",
        "first-omitted-term theorem is proved",
        "actual remainder theorem is proved",
        "interval certificate is complete",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path, grid_path: Path) -> tuple[list[IntervalizationIssue], dict]:
    artifact = load_json(target_path)
    issues: list[IntervalizationIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, grid_path))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, exact_rows, sufficient_rows, open_rows, rejected_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_rows, sufficient_rows, open_rows, rejected_rows))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    issues, summary = validate(target, note, grid_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-INTERVALIZATION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian intervalization target: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('obligation_rows')} obligations, "
            f"{summary.get('open_requirement_rows')} open requirements, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
