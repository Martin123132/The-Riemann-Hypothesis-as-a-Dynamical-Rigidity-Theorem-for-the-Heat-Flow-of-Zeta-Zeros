#!/usr/bin/env python3
"""Validate the negative-lambda defect-tail theorem target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_defect_tail_theorem_target import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_target,
)


REQUIRED_ROW_IDS = {
    "nldtt_01_exact_tail_statement",
    "nldtt_02_finite_anchor_handoff",
    "nldtt_03_uniform_saddle_route",
    "nldtt_04_ratio_recurrence_route",
    "nldtt_05_scaled_defect_shortcut_rejected",
    "nldtt_06_finite_extension_insufficient",
    "nldtt_07_forbidden_endpoint_shortcuts",
    "nldtt_08_conditional_cone_entry_application",
}

ALLOWED_ROLES = {
    "open_statement",
    "exact_requirement",
    "live_route",
    "rejected_route",
    "insufficient_route",
    "circular_route",
    "conditional_application",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Defect-Tail Theorem Target",
    "Status: open theorem target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_defect_tail_theorem_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_defect_tail_theorem_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_tail_theorem_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_tail_theorem_target.py",
    "validated Jensen-window PF negative-lambda defect-tail theorem target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes",
    "0 <= d_k <= 2/(2*k+1) for all k >= 150",
    "d_(k+1) <= d_k for all k >= 149",
    "Without such a theorem, the next finite bridge needs `A_151`",
    "uniform saddle/Laplace control",
    "direct ratio-recurrence",
    "scaled-defect nonincreasing shortcut is rejected",
    "outputs/jensen_window_pf_negative_lambda_tail_barrier_k150_scout.md",
    "validated Jensen-window PF negative-lambda tail-barrier scout: 447 cone-buffer rows, 444 defect-monotone rows, 179 one-third-buffer rows, 444 scaled-defect increase rows, 1 rejected candidate, 0 issues",
    "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md",
    "validated Jensen-window PF negative-lambda scaled-defect frontier scout: 447 scaled rows, 447 cone rows, 430 half-width rows, 179 one-third rows, 268 one-third failures, 444 scaled-increase rows, 0 issues",
    "outputs/jensen_window_pf_negative_lambda_defect_recurrence_scout.md",
    "validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues",
    "0 <= d_k <= 2/(3*(2*k+1))",
    "one-third buffer should no longer be promoted as an all-tail",
    "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
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
    if target.get("kind") != "jensen_window_pf_negative_lambda_defect_tail_theorem_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") != "open_theorem_target":
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    if target.get("target_id") != "target_negative_lambda_defect_tail":
        issues.append(issue("<target>", "bad-target-id", repr(target.get("target_id"))))
    for key in (
        "source_tail_barrier_scout",
        "source_scaled_defect_frontier_scout",
        "source_cone_entry_target",
        "tail_barrier_json",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<target>", target.get(key)))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("open theorem target", "does not prove", "tail theorem", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(target: dict) -> list[TargetIssue]:
    issues: list[TargetIssue] = []
    ref = target.get("tail_barrier_json")
    if not isinstance(ref, str):
        return [issue("recompute", "missing-tail-barrier-ref", repr(ref))]
    try:
        recomputed = build_target(REPO_ROOT / ref)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    for key in ("finite_anchor", "target_rows", "invariants", "summary"):
        if target.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded target differs from recomputed target"))
    return issues


def validate_anchor(target: dict) -> list[TargetIssue]:
    anchor = target.get("finite_anchor", {})
    expected = {
        "coefficient_k_max": 150,
        "checked_x_max": 149,
        "tail_lower_start_k": 150,
        "tail_monotone_bridge_k": 149,
        "next_finite_coefficient_needed": "A_151",
    }
    issues: list[TargetIssue] = []
    for key, value in expected.items():
        if anchor.get(key) != value:
            issues.append(issue("finite_anchor", f"bad-{key}", f"{anchor.get(key)!r} != {value!r}"))
    if anchor.get("lambdas") != ["-100.0", "-50.0", "-25.0"]:
        issues.append(issue("finite_anchor", "bad-lambdas", repr(anchor.get("lambdas"))))
    return issues


def validate_rows(target: dict) -> tuple[list[TargetIssue], int, int, int]:
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
        if row.get("role") in {"exact_requirement", "conditional_application", "insufficient_route"} and "not" not in text:
            issues.append(issue(row_id, "weak-boundary", text))
        if row.get("role") == "rejected_route" and "rejected" not in text and "do not use" not in text:
            issues.append(issue(row_id, "missing-rejection-language", text))
        if row.get("role") == "circular_route" and "circular" not in text:
            issues.append(issue(row_id, "missing-circular-language", text))
    return issues, len(rows), ready_count, live_count


def validate_summary(target: dict, row_count: int, ready_count: int, live_count: int) -> list[TargetIssue]:
    summary = target.get("summary", {})
    expected = {
        "target_rows": 8,
        "ready_to_apply_rows": 0,
        "live_routes": 2,
        "rejected_routes": 1,
        "insufficient_routes": 1,
        "conditional_application_rows": 1,
        "tail_lower_start_k": 150,
        "tail_monotone_bridge_k": 149,
        "next_finite_coefficient_needed": "A_151",
        "target_closing": False,
    }
    issues: list[TargetIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if live_count != 2:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "defect-tail",
        "x_1..x_149",
        "k >= 150",
        "k >= 149",
        "uniform saddle",
        "ratio-recurrence",
        "rejected",
        "one-third",
        "half-width",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("no row", "open_target", "finite extensions", "scaled-defect", "one-third", "half-width", "lambda <= 0"):
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
        "lambda <= 0 is proved",
        "cone entry is proved",
        "all-k defect-tail theorem is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[TargetIssue], dict]:
    target = load_json(target_path)
    issues: list[TargetIssue] = []
    issues.extend(validate_top_level(target))
    issues.extend(validate_recomputed(target))
    issues.extend(validate_anchor(target))
    row_issues, row_count, ready_count, live_count = validate_rows(target)
    issues.extend(row_issues)
    issues.extend(validate_summary(target, row_count, ready_count, live_count))
    issues.extend(validate_note(note_path))
    return issues, target.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, summary = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-DEFECT-TAIL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda defect-tail theorem target: "
            f"{summary.get('target_rows')} rows, "
            f"{len(issues)} issues, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows, "
            f"{summary.get('live_routes')} live routes"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
