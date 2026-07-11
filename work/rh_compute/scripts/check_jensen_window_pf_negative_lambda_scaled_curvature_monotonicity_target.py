#!/usr/bin/env python3
"""Validate the scaled-curvature monotonicity theorem target."""

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

from jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlscmt_01_scaled_curvature_statement",
    "nlscmt_02_linear_barrier_equivalence",
    "nlscmt_03_raw_corridor_sufficient_chain",
    "nlscmt_04_repaired_k300_finite_anchor",
    "nlscmt_05_retired_fixed_wall_boundary",
    "nlscmt_06_required_monotone_upper_side",
    "nlscmt_07_live_log_ratio_recurrence_route",
    "nlscmt_08_live_asymptotic_curvature_route",
    "nlscmt_09_generic_shortcuts_blocked",
    "nlscmt_10_acceptance_gate",
}

ALLOWED_ROLES = {
    "open_statement",
    "exact_reduction",
    "finite_anchor",
    "rejected_shortcut",
    "open_dependency",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Scaled-Curvature Monotonicity Target",
    "Status: open theorem target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.py",
    "validated Jensen-window PF negative-lambda scaled-curvature monotonicity target: 10 rows, 0 issues, 2 live routes, 894 scaled-curvature increase rows, 0 ready-to-apply rows",
    "C_(k+1) >= C_k",
    "B_(k+1) >= ((2*k+1)/(2*k+3))*B_k",
    "B_k > 0 rows: 897 / 897",
    "C_(k+1)-C_k positive rows: 894 / 894",
    "retired C_k<=2/3 failures: 718 / 897",
    "plus B_(k+1)<=B_k gives the two-sided curvature corridor",
    "zeta-specific log-ratio recurrence",
    "signed/tilted saddle analysis with uniform remainders",
    "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
    "outputs/jensen_window_pf_negative_lambda_coefficient_curvature_corridor_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
    "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
)


@dataclass(frozen=True)
class ScaledCurvatureIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ScaledCurvatureIssue:
    return ScaledCurvatureIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ScaledCurvatureIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identity() -> list[ScaledCurvatureIssue]:
    k, b, bn = sp.symbols("k b bn")
    alpha = (2 * k + 1) / (2 * k + 3)
    c_gap = (2 * k + 3) * bn - (2 * k + 1) * b
    linear_gap = (2 * k + 3) * (bn - alpha * b)
    if sp.simplify(c_gap - linear_gap) != 0:
        return [issue("symbolic", "bad-linear-equivalence", str(sp.simplify(c_gap - linear_gap)))]
    return []


def validate_top_level(artifact: dict) -> list[ScaledCurvatureIssue]:
    issues: list[ScaledCurvatureIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "open_theorem_target":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    if artifact.get("target_id") != "target_negative_lambda_scaled_curvature_monotonicity":
        issues.append(issue("<artifact>", "bad-target-id", repr(artifact.get("target_id"))))
    for key in (
        "source_linear_barrier_scout",
        "source_curvature_corridor_bridge",
        "source_raw_corridor_target",
        "source_k300_obstruction",
        "source_monotone_contraction_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("open theorem target", "scaled-curvature", "does not prove", "monotone-contraction", "raw-corridor", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[ScaledCurvatureIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_artifact([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[ScaledCurvatureIssue] = []
    for key in ("target_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[ScaledCurvatureIssue], int, int, int, int]:
    rows = artifact.get("target_rows", [])
    issues: list[ScaledCurvatureIssue] = []
    if not isinstance(rows, list):
        return [issue("target_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    live_routes = 0
    rejected_shortcuts = 0
    ready_to_apply = 0
    exact_available = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("target_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim_if_proved", "source_artifacts", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        readiness = row.get("readiness")
        if row.get("role") == "exact_reduction":
            exact_available += 1
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-open-readiness", repr(readiness)))
        if readiness == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "live_route":
            live_routes += 1
        if row.get("role") == "rejected_shortcut":
            rejected_shortcuts += 1
            text = f"{row.get('gap', '')} {row.get('proof_boundary', '')}".lower()
            if not any(marker in text for marker in ("reject", "blocked", "finite-reject", "not evidence")):
                issues.append(issue(row_id, "weak-rejection-language", text))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("open", "not", "only", "finite", "exact", "rejected", "live")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), live_routes, rejected_shortcuts, exact_available if ready_to_apply == 0 else -exact_available


def validate_summary(
    artifact: dict,
    row_count: int,
    live_routes: int,
    rejected_shortcuts: int,
    exact_available: int,
) -> list[ScaledCurvatureIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "target_rows": 10,
        "exact_available_rows": 2,
        "finite_anchor_rows": 1,
        "open_dependency_rows": 1,
        "live_routes": 2,
        "rejected_shortcut_rows": 2,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "open_theorem_target": True,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "coefficient_cap": 300,
        "b_wall_rows": 897,
        "scaled_curvature_total_rows": 897,
        "adjacent_total_rows": 894,
        "scaled_curvature_increase_rows": 894,
        "scaled_curvature_increase_failures": 0,
        "scaled_curvature_increase_inconclusive": 0,
        "two_thirds_failure_rows": 718,
    }
    issues: list[ScaledCurvatureIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 10:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if live_routes != 2:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    if rejected_shortcuts != 2:
        issues.append(issue("summary", "bad-rejected-shortcut-count", str(rejected_shortcuts)))
    if exact_available != 2:
        issues.append(issue("summary", "bad-exact-available-count", str(exact_available)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "replacement",
        "fixed 2/3",
        "scaled-curvature monotonicity",
        "c_(k+1)>=c_k",
        "linear barrier",
        "894/894",
        "718/897",
        "monotone-contraction",
        "raw-corridor",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "open_theorem_target", "fixed 2/3", "finite", "monotone-contraction", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ScaledCurvatureIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ScaledCurvatureIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "scaled-curvature monotonicity theorem is proved",
        "raw-corridor theorem is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[ScaledCurvatureIssue], dict]:
    artifact = load_json(target_path)
    issues: list[ScaledCurvatureIssue] = []
    issues.extend(validate_symbolic_identity())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, live_routes, rejected_shortcuts, exact_available = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, live_routes, rejected_shortcuts, exact_available))
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
            print(f"JWPF-NEG-LAMBDA-SCALED-CURVATURE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda scaled-curvature monotonicity target: "
            f"{summary.get('target_rows')} rows, {len(issues)} issues, "
            f"{summary.get('live_routes')} live routes, "
            f"{summary.get('scaled_curvature_increase_rows')} scaled-curvature increase rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
