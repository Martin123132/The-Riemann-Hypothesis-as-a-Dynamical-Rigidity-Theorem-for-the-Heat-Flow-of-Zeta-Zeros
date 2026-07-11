#!/usr/bin/env python3
"""Validate the scaled-curvature log-ceiling bridge."""

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

from jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlsclcb_01_log_curvature_coordinate",
    "nlsclcb_02_affine_ceiling_equivalence",
    "nlsclcb_03_ceiling_dominates_raw_log_upper",
    "nlsclcb_04_repaired_k300_affine_ceiling",
    "nlsclcb_05_sharpness_warning",
    "nlsclcb_06_scaled_log_corridor_target",
    "nlsclcb_07_live_affine_recurrence_route",
    "nlsclcb_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "exact_sufficient_condition",
    "finite_stress",
    "finite_sharpness",
    "open_dependency",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Scaled-Curvature Log-Ceiling Bridge",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.py",
    "validated Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: 8 rows, 0 issues, 894 scaled-ceiling rows, 894 scaled-log-corridor rows, 894 ceiling-dominance rows, 0 ready-to-apply rows",
    "C_(k+1)>=C_k",
    "delta_k <= h_(k+1)-((2*k+1)/(2*k+3))*h_k-2*p_k/(2*k+3)",
    "raw_log_upper(p_k)-scaled_ceiling(p_k)=alpha_k*B_k-L_k(B_k)>=0",
    "p-wall rows: 897 lower and 897 upper / 897",
    "scaled-ceiling rows: 894 / 894",
    "scaled-log-corridor rows: 894 / 894",
    "ceiling-dominance rows: 894 / 894",
    "min raw-upper minus scaled-ceiling slack: 4.569274603858794776E-9 at lambda=-100.0, k=298",
    "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
    "outputs/jensen_window_pf_negative_lambda_raw_log_decrement_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
    "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
)


@dataclass(frozen=True)
class LogCeilingIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> LogCeilingIssue:
    return LogCeilingIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[LogCeilingIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[LogCeilingIssue]:
    k, h, hn, p, d = sp.symbols("k h hn p d")
    alpha = (2 * k + 1) / (2 * k + 3)
    ceiling = hn - alpha * h - 2 * p / (2 * k + 3)
    c_gap = (2 * k + 3) * (hn - p - d) - (2 * k + 1) * (h - p)
    issues: list[LogCeilingIssue] = []
    if sp.simplify(c_gap - (2 * k + 3) * (ceiling - d)) != 0:
        issues.append(issue("symbolic", "bad-affine-ceiling-equivalence", str(sp.simplify(c_gap - (2 * k + 3) * (ceiling - d)))))

    a, t = sp.symbols("a t")
    claimed_derivative = a * (1 - a) * (1 - t) / (1 - a + a * t)
    direct_derivative = a - a * t / (1 - a + a * t)
    if sp.simplify(claimed_derivative - direct_derivative) != 0:
        issues.append(issue("calculus", "bad-dominance-derivative", str(sp.simplify(claimed_derivative - direct_derivative))))
    return issues


def validate_top_level(artifact: dict) -> list[LogCeilingIssue]:
    issues: list[LogCeilingIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_scaled_curvature_target",
        "source_raw_log_bridge",
        "source_linear_barrier_scout",
        "source_raw_corridor_target",
        "source_monotone_contraction_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "finite", "does not prove", "affine recurrence", "monotone-contraction", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[LogCeilingIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_artifact([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[LogCeilingIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[LogCeilingIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[LogCeilingIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_available = 0
    live_routes = 0
    ready_to_apply = 0
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
        readiness = row.get("readiness")
        if row.get("role") in {"exact_reduction", "exact_sufficient_condition"}:
            exact_available += 1
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-open-readiness", repr(readiness)))
        if readiness == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "live_route":
            live_routes += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "open", "live", "hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), exact_available, live_routes if ready_to_apply == 0 else -live_routes


def validate_summary(artifact: dict, row_count: int, exact_available: int, live_routes: int) -> list[LogCeilingIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "exact_available_rows": 3,
        "finite_stress_rows": 1,
        "finite_sharpness_rows": 1,
        "open_dependency_rows": 1,
        "live_routes": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "coefficient_cap": 300,
        "raw_log_total_rows": 897,
        "adjacent_total_rows": 894,
        "p_lower_rows": 897,
        "p_upper_rows": 897,
        "scaled_ceiling_rows": 894,
        "raw_upper_rows": 894,
        "ceiling_dominance_rows": 894,
        "log_lower_rows": 894,
        "scaled_log_corridor_rows": 894,
        "scaled_width_rows": 894,
    }
    issues: list[LogCeilingIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_available != 3:
        issues.append(issue("summary", "bad-exact-count", str(exact_available)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("scaled-curvature", "affine log-ratio ceiling", "delta_k<=", "894/894", "897/897", "nonlinear raw-log upper wall"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "affine log-ratio", "finite stress", "nonlinear raw-log", "monotone-contraction", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[LogCeilingIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[LogCeilingIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "affine recurrence is proved",
        "scaled-curvature monotonicity is proved",
        "raw-corridor theorem is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[LogCeilingIssue], dict]:
    artifact = load_json(target_path)
    issues: list[LogCeilingIssue] = []
    issues.extend(validate_symbolic_identities())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, exact_available, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_available, live_routes))
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
            print(f"JWPF-NEG-LAMBDA-SCALED-LOG-CEILING {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda scaled-curvature log-ceiling bridge: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('scaled_ceiling_rows')} scaled-ceiling rows, "
            f"{summary.get('scaled_log_corridor_rows')} scaled-log-corridor rows, "
            f"{summary.get('ceiling_dominance_rows')} ceiling-dominance rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
