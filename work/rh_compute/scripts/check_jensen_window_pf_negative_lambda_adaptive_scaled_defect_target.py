#!/usr/bin/env python3
"""Validate the negative-lambda adaptive scaled-defect tail target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_adaptive_scaled_defect_target import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_target,
)


REQUIRED_ROW_IDS = {
    "nlasdt_01_exact_scaled_cone",
    "nlasdt_02_k200_frontier",
    "nlasdt_03_fixed_buffers_rejected",
    "nlasdt_04_adaptive_or_exact_cone_statement",
    "nlasdt_05_uniform_saddle_route",
    "nlasdt_06_ratio_comparison_route",
    "nlasdt_07_finite_anchor_handoff",
    "nlasdt_08_conditional_application",
}

ALLOWED_ROLES = {
    "exact_reformulation",
    "finite_frontier",
    "rejected_route",
    "open_statement",
    "live_route",
    "exact_requirement",
    "conditional_application",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Adaptive Scaled-Defect Target",
    "Status: open theorem target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_adaptive_scaled_defect_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.py",
    "validated Jensen-window PF negative-lambda adaptive scaled-defect target: 8 rows, 0 issues, 2 live routes, 597 exact-cone rows, 76 half-width failures",
    "prove 0 <= s_k <= 1 for all k >= 200",
    "prove d_(k+1) <= d_k for all k >= 199 separately",
    "exact cone rows: 597 / 597",
    "fixed half-width rows: 521 / 597",
    "fixed half-width failure rows: 76",
    "one-third failure rows: 418",
    "max scaled defect: 5.376643171065356005E-1 at lambda=-25.0, k=199",
    "lambda=-50.0: first failure k=191",
    "lambda=-25.0: first failure k=133",
    "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md",
    "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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
    if target.get("kind") != "jensen_window_pf_negative_lambda_adaptive_scaled_defect_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") != "open_theorem_target":
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    if target.get("target_id") != "target_negative_lambda_adaptive_scaled_defect_tail":
        issues.append(issue("<target>", "bad-target-id", repr(target.get("target_id"))))
    for key in (
        "source_scaled_defect_frontier",
        "source_defect_tail_target",
        "frontier_json",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<target>", target.get(key)))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("open theorem target", "adaptive", "exact-cone", "does not prove", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(target: dict) -> list[TargetIssue]:
    ref = target.get("frontier_json")
    if not isinstance(ref, str):
        return [issue("recompute", "missing-frontier-ref", repr(ref))]
    try:
        recomputed = build_target(REPO_ROOT / ref)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[TargetIssue] = []
    for key in ("target_rows", "summary", "invariants"):
        if target.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded target differs from recomputed target"))
    return issues


def validate_rows(target: dict) -> tuple[list[TargetIssue], int, int]:
    rows = target.get("target_rows", [])
    issues: list[TargetIssue] = []
    if not isinstance(rows, list):
        return [issue("target_rows", "bad-rows", repr(type(rows)))], 0, 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))
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
        if row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("role") == "live_route":
            live_count += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        text = " ".join(str(row.get(key, "")) for key in ("gap", "acceptance_test", "proof_boundary")).lower()
        if row.get("role") == "rejected_route" and "do not use" not in text and "reject" not in text:
            issues.append(issue(row_id, "missing-rejection-language", text))
        if row.get("role") in {"open_statement", "conditional_application", "finite_frontier", "exact_requirement"}:
            if not any(marker in text for marker in ("not", "no ", "only", "finite")):
                issues.append(issue(row_id, "weak-boundary", text))
    return issues, len(rows), live_count


def validate_summary(target: dict, row_count: int, live_count: int) -> list[TargetIssue]:
    summary = target.get("summary", {})
    expected = {
        "target_rows": 8,
        "ready_to_apply_rows": 0,
        "live_routes": 2,
        "rejected_routes": 1,
        "coefficient_k_max": 200,
        "checked_x_max": 199,
        "tail_start_k": 200,
        "next_finite_coefficient_needed": "A_201",
        "finite_exact_cone_rows": 597,
        "finite_scaled_rows": 597,
        "finite_half_width_rows": 521,
        "finite_half_width_failure_rows": 76,
        "finite_one_third_failure_rows": 418,
        "target_closing": False,
    }
    issues: list[TargetIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if live_count != 2:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    first_failures = summary.get("first_half_width_failures", [])
    if not isinstance(first_failures, list) or len(first_failures) != 2:
        issues.append(issue("summary", "bad-first-half-width-failures", repr(first_failures)))
    else:
        expected_first = [("-50.0", 191), ("-25.0", 133)]
        actual = [(row.get("lam"), row.get("first_failure_k")) for row in first_failures]
        if actual != expected_first:
            issues.append(issue("summary", "bad-first-half-width-failures", repr(actual)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("k200", "exact scaled cone", "597/597", "half-width", "76 rows", "one-third", "418 rows", "adaptive-envelope", "monotone defect"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("no row", "open_target", "one-third", "half-width", "nonincrease", "exact cone", "lambda <= 0"):
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
        "adaptive scaled-defect theorem is proved",
        "exact-cone theorem is proved",
        "cone entry is proved",
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
    row_issues, row_count, live_count = validate_rows(target)
    issues.extend(row_issues)
    issues.extend(validate_summary(target, row_count, live_count))
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
            print(f"JWPF-NEG-LAMBDA-ADAPTIVE-SCALED-DEFECT {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda adaptive scaled-defect target: "
            f"{summary.get('target_rows')} rows, {len(issues)} issues, "
            f"{summary.get('live_routes')} live routes, "
            f"{summary.get('finite_exact_cone_rows')} exact-cone rows, "
            f"{summary.get('finite_half_width_failure_rows')} half-width failures"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
