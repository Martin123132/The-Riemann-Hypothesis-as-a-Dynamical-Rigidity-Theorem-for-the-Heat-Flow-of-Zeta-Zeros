#!/usr/bin/env python3
"""Validate the relative-Gaussian next-increment stencil stress diagnostic."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgniss_01_next_increment_coordinate",
    "nlrgniss_02_pointwise_budget_failure",
    "nlrgniss_03_structured_stencil_survival",
    "nlrgniss_04_half_safety_not_respected",
    "nlrgniss_05_degree16_unavailable_frontier",
    "nlrgniss_06_pointwise_triangle_route_rejected_for_current_baselines",
    "nlrgniss_07_live_direct_stencil_tail_route",
    "nlrgniss_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_diagnostic",
    "finite_frontier",
    "rejected_route",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Next-Increment Stencil Stress",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: 8 rows, 0 issues, 2 tested next-increment rows, 2 pointwise budget failures, 2 stencil-sign-preserving rows, 0 ready-to-apply rows",
    "tested next-increment rows: 2",
    "missing next-increment rows: 2",
    "pointwise budget failures: 2",
    "stencil-sign-preserving rows: 2",
    "half-safety stencil failures: 2",
    "worst pointwise over-budget factor: 3.010798908654295615E+3",
    "worst stencil abs over half margin: 3.998025843926772743E+0",
    "nlrgts_M5_T2000: M->6",
    "nlrgts_M6_T2000: M->7",
    "degree-16 Taylor ratio is not available",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md",
    "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
)


@dataclass(frozen=True)
class NextIncrementStressIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> NextIncrementStressIssue:
    return NextIncrementStressIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[NextIncrementStressIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_stencils() -> list[NextIncrementStressIssue]:
    k = sp.symbols("k")
    dm, d0, dp, dpp = sp.symbols("dm d0 dp dpp")
    b_error = 2 * d0 - dm - dp
    next_b_error = 2 * dp - d0 - dpp
    companion = b_error - next_b_error
    weighted_gap = (2 * k + 3) * next_b_error - (2 * k + 1) * b_error
    expected_companion = -dm + 3 * d0 - 3 * dp + dpp
    expected_weighted = (2 * k + 1) * dm - (6 * k + 5) * d0 + (6 * k + 7) * dp - (2 * k + 3) * dpp
    issues: list[NextIncrementStressIssue] = []
    if sp.simplify(companion - expected_companion) != 0:
        issues.append(issue("symbolic", "bad-companion-stencil", str(sp.simplify(companion - expected_companion))))
    if sp.simplify(weighted_gap - expected_weighted) != 0:
        issues.append(issue("symbolic", "bad-weighted-stencil", str(sp.simplify(weighted_gap - expected_weighted))))
    return issues


def validate_top_level(artifact: dict) -> list[NextIncrementStressIssue]:
    issues: list[NextIncrementStressIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_pointwise_tail_budget",
        "source_high_order_taylor_scout",
        "source_uniform_remainder_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite", "does not prove", "infinite", "scaled-curvature", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[NextIncrementStressIssue]:
    params = artifact.get("matrix_rows", [{}])[1].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            int(params.get("max_taylor_degree", 14)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 256)),
            int(params.get("tail_start_k", 22)),
            [int(value) for value in params.get("sample_T_values", [25, 50, 100, 200, 500, 1000, 2000])],
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[NextIncrementStressIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[NextIncrementStressIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[1].get("diagnostics", {})
    expected = {
        "tested_next_increment_rows": 2,
        "missing_next_increment_rows": 2,
        "pointwise_budget_failure_rows": 2,
        "stencil_sign_preserving_rows": 2,
        "half_safety_stencil_failure_rows": 2,
        "worst_pointwise_over_budget": {"sample": "3.010798908654295615E+3", "source_row": "nlrgts_M6_T2000"},
        "worst_stencil_abs_over_half_margin": {
            "sample": "3.998025843926772743E+0",
            "source_row": "nlrgts_M6_T2000",
            "stencil": "companion",
        },
    }
    issues: list[NextIncrementStressIssue] = []
    for key, value in expected.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))
    stress_rows = diagnostics.get("stress_rows", [])
    if [row.get("source_row") for row in stress_rows] != ["nlrgts_M5_T2000", "nlrgts_M6_T2000"]:
        issues.append(issue("diagnostics", "bad-stress-row-order", repr([row.get("source_row") for row in stress_rows])))
    for row in stress_rows:
        row_id = str(row.get("source_row", "<missing-source>"))
        if row.get("pointwise_budget_satisfied") is not False:
            issues.append(issue(row_id, "pointwise-budget-not-failed", repr(row.get("pointwise_budget_satisfied"))))
        if row.get("stencil_signs_preserved") is not True:
            issues.append(issue(row_id, "stencil-sign-not-preserved", repr(row.get("stencil_signs_preserved"))))
        if row.get("half_safety_stencils_satisfied") is not False:
            issues.append(issue(row_id, "half-safety-not-failed", repr(row.get("half_safety_stencils_satisfied"))))
        if len(row.get("stencil_increments", [])) != 3:
            issues.append(issue(row_id, "bad-stencil-increment-count", repr(len(row.get("stencil_increments", [])))))
    missing_rows = diagnostics.get("missing_rows", [])
    if [row.get("source_row") for row in missing_rows] != ["nlrgts_M7_T1000", "nlrgts_M7_T2000"]:
        issues.append(issue("diagnostics", "bad-missing-row-order", repr([row.get("source_row") for row in missing_rows])))
    for row in missing_rows:
        if "degree-16" not in str(row.get("reason", "")):
            issues.append(issue(str(row.get("source_row", "<missing-source>")), "bad-missing-reason", repr(row.get("reason"))))
    return issues


def validate_rows(artifact: dict) -> tuple[list[NextIncrementStressIssue], int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[NextIncrementStressIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    ready_to_apply = 0
    live_routes = 0
    rejected_routes = 0
    finite_diagnostics = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("role") == "finite_diagnostic":
            finite_diagnostics += 1
        if row.get("role") == "live_route":
            live_routes += 1
        if row.get("role") == "rejected_route":
            rejected_routes += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") not in {"available_exact", "not_ready_to_apply"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "live", "hygiene", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), finite_diagnostics, live_routes, rejected_routes if ready_to_apply == 0 else -rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    finite_diagnostics: int,
    live_routes: int,
    rejected_routes: int,
) -> list[NextIncrementStressIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "tested_next_increment_rows": 2,
        "missing_next_increment_rows": 2,
        "pointwise_budget_failure_rows": 2,
        "stencil_sign_preserving_rows": 2,
        "half_safety_stencil_failure_rows": 2,
        "worst_pointwise_over_budget": "3.010798908654295615E+3",
        "worst_stencil_abs_over_half_margin": "3.998025843926772743E+0",
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[NextIncrementStressIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if finite_diagnostics != 3:
        issues.append(issue("summary", "bad-finite-diagnostic-count", str(finite_diagnostics)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    if rejected_routes != 1:
        issues.append(issue("summary", "bad-rejected-route-count", str(rejected_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("pointwise half-safety", "3.010798908654295615e+3", "preserve all tested finite signs", "direct signed stencil-tail"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "not promoted", "pointwise budget", "stencil sign survival", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[NextIncrementStressIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[NextIncrementStressIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "infinite taylor-tail estimate is proved",
        "uniform taylor-tail theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[NextIncrementStressIssue], dict]:
    artifact = load_json(target_path)
    issues: list[NextIncrementStressIssue] = []
    issues.extend(validate_symbolic_stencils())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, finite_diagnostics, live_routes, rejected_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, finite_diagnostics, live_routes, rejected_routes))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    issues, summary = validate(target, note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-NEXT-INCREMENT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('tested_next_increment_rows')} tested next-increment rows, "
            f"{summary.get('pointwise_budget_failure_rows')} pointwise budget failures, "
            f"{summary.get('stencil_sign_preserving_rows')} stencil-sign-preserving rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
