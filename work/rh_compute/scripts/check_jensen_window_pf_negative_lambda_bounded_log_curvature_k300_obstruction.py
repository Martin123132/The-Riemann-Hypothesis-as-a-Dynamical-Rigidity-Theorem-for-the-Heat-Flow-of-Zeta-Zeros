#!/usr/bin/env python3
"""Validate the k300 obstruction to the fixed bounded log-curvature wall."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nllcko_01_exact_scaled_curvature_rewrite",
    "nllcko_02_repaired_k300_two_thirds_obstruction",
    "nllcko_03_first_checked_failure",
    "nllcko_04_worst_checked_failure",
    "nllcko_05_k22_prefix_scope",
    "nllcko_06_replacement_scaled_monotonicity_route",
    "nllcko_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_obstruction",
    "finite_counterexample",
    "finite_extremum",
    "historical_scope",
    "replacement_route_diagnostic",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Bounded Log-Curvature k300 Obstruction",
    "Status: finite obstruction gate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction.py",
    "validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: 7 rows, 0 issues, 718 two-thirds failures, 894 scaled-curvature increase rows, 0 ready-to-apply rows",
    "B_k <= 2/(3*(2*k+1))  iff  C_k <= 2/3",
    "C_k <= 2/3 rows: 179 / 897",
    "C_k > 2/3 rows: 718 / 897",
    "inconclusive rows: 0",
    "lambda=-25.0, k=31",
    "max C_k = 1.144219413064916367E+0 at lambda=-25.0, k=299",
    "C_(k+1)-C_k positive rows: 894 / 894",
    "B_k-B_(k+1) positive rows: 894 / 894",
    "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
    "outputs/jensen_window_pf_negative_lambda_k300_precision_repair_audit.md",
    "outputs/jensen_window_pf_negative_lambda_linear_curvature_barrier_scout.md",
    "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
)


@dataclass(frozen=True)
class ObstructionIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ObstructionIssue:
    return ObstructionIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ObstructionIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[ObstructionIssue]:
    issues: list[ObstructionIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_bounded_log_curvature_k300_obstruction":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite obstruction gate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_bounded_log_curvature_target",
        "source_k300_precision_repair",
        "source_curvature_corridor_bridge",
        "source_linear_barrier_scout",
        "source_raw_corridor_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite obstruction", "fixed 2/3", "does not prove", "raw-corridor", "cone entry", "jwpf_06", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[ObstructionIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_artifact([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[ObstructionIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[ObstructionIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[ObstructionIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
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
        if row.get("role") == "exact_reduction":
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(readiness)))
        if readiness == "ready_to_apply":
            ready_to_apply += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "exact", "not", "only", "hygiene", "scope")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), ready_to_apply


def validate_summary(artifact: dict, row_count: int, ready_to_apply: int) -> list[ObstructionIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 7,
        "exact_reduction_rows": 1,
        "finite_obstruction_rows": 3,
        "historical_scope_rows": 1,
        "replacement_route_rows": 1,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "coefficient_cap": 300,
        "scaled_curvature_total_rows": 897,
        "two_thirds_bound_rows": 179,
        "two_thirds_failure_rows": 718,
        "two_thirds_inconclusive_rows": 0,
        "adjacent_total_rows": 894,
        "scaled_curvature_increase_rows": 894,
        "b_decrease_rows": 894,
        "target_retired": True,
    }
    issues: list[ObstructionIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 7:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ready_to_apply != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_to_apply)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("fixed bounded log-curvature", "c_k=(2*k+1)*b_k<=2/3", "179/897", "718/897", "linear/scaled curvature corridor"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "fixed 2/3", "linear curvature-barrier", "finite", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ObstructionIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ObstructionIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "raw-corridor theorem is proved",
        "linear curvature barrier is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[ObstructionIssue], dict]:
    artifact = load_json(target_path)
    issues: list[ObstructionIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, ready_to_apply = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, ready_to_apply))
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
            print(f"JWPF-NEG-LAMBDA-BOUNDED-CURV-K300-OBSTRUCTION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda bounded log-curvature k300 obstruction: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('two_thirds_failure_rows')} two-thirds failures, "
            f"{summary.get('scaled_curvature_increase_rows')} scaled-curvature increase rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
