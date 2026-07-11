#!/usr/bin/env python3
"""Validate the negative-lambda half-width scaled-defect tail target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_half_width_tail_target import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_target,
)


REQUIRED_ROW_IDS = {
    "nlhwt_01_scaled_defect_definition",
    "nlhwt_02_half_width_tail_statement",
    "nlhwt_03_monotone_defect_still_required",
    "nlhwt_04_finite_anchor",
    "nlhwt_05_uniform_saddle_route",
    "nlhwt_06_ratio_recurrence_route",
    "nlhwt_07_falsification_protocol",
    "nlhwt_08_rejected_shortcuts",
    "nlhwt_09_conditional_application",
}

ALLOWED_ROLES = {
    "exact_reformulation",
    "open_statement",
    "exact_requirement",
    "finite_anchor",
    "live_route",
    "repair_route",
    "stress_protocol",
    "rejected_route",
    "conditional_application",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Half-Width Tail Target",
    "Status: finite-rejected target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_half_width_tail_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_half_width_tail_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_half_width_tail_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_half_width_tail_target.py",
    "validated Jensen-window PF negative-lambda half-width tail target: 9 rows, 0 issues, 0 ready-to-apply rows, 0 live routes, 430 half-width rows, 17 half-width failures",
    "s_k = ((2*k+1)/2) * d_k",
    "0 <= s_k <= 1/2 for all k >= 150",
    "d_(k+1) <= d_k for all k >= 149 remains separately required",
    "Finite k150 stress rejects this fixed half-width target",
    "half-width rows: 430 / 447",
    "half-width failure rows: 17",
    "one-third failure rows: 268",
    "minimum half-width margin: -1.107359370078869990E-2 at lambda=-25.0, k=149",
    "maximum scaled defect: 5.110735937007886999E-1 at lambda=-25.0, k=149",
    "minimum one-third margin: -1.777402603674553666E-1 at lambda=-25.0, k=149",
    "lambda=-25.0: first failure k=133",
    "none for the fixed half-width threshold",
    "one-third-width buffer",
    "fixed half-width buffer",
    "scaled-defect nonincrease",
    "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k150_scout.md",
    "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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
    if target.get("kind") != "jensen_window_pf_negative_lambda_half_width_tail_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") not in {"open_theorem_target", "finite_rejected_target"}:
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    if target.get("target_id") != "target_negative_lambda_half_width_tail":
        issues.append(issue("<target>", "bad-target-id", repr(target.get("target_id"))))
    for key in (
        "source_scaled_defect_frontier",
        "source_defect_tail_target",
        "source_cone_entry_target",
        "frontier_json",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<target>", target.get(key)))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("does not prove", "half-width", "monotone defect", "cone entry", "lambda <= 0"):
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
        text = " ".join(str(row.get(key, "")) for key in ("claim_if_proved", "gap", "acceptance_test", "proof_boundary")).lower()
        if row.get("role") in {"open_statement", "exact_requirement", "finite_anchor", "conditional_application", "stress_protocol"}:
            if not any(marker in text for marker in ("not", "no ", "only", "separately", "finite")):
                issues.append(issue(row_id, "weak-boundary", text))
        if row.get("role") == "rejected_route" and "rejected" not in text and "do not use" not in text:
            issues.append(issue(row_id, "missing-rejection-language", text))
        if "endpoint" in text and "forbidden" not in " ".join(str(item) for item in target.get("invariants", [])).lower():
            issues.append(issue(row_id, "endpoint-language-without-invariant", text))
    return issues, len(rows), ready_count, live_count


def validate_summary(target: dict, row_count: int, ready_count: int, live_count: int) -> list[TargetIssue]:
    summary = target.get("summary", {})
    expected = {
        "target_rows": 9,
        "ready_to_apply_rows": 0,
        "live_routes": 0,
        "repair_routes": 2,
        "rejected_routes": 2,
        "stress_protocol_rows": 1,
        "conditional_application_rows": 1,
        "coefficient_k_max": 150,
        "checked_x_max": 149,
        "tail_lower_start_k": 150,
        "tail_monotone_bridge_k": 149,
        "next_finite_coefficient_needed": "A_151",
        "finite_half_width_rows": 430,
        "finite_half_width_failure_rows": 17,
        "finite_scaled_rows": 447,
        "finite_one_third_failure_rows": 268,
        "target_rejected": True,
        "target_closing": False,
    }
    issues: list[TargetIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if live_count != 0:
        issues.append(issue("summary", "bad-live-count", str(live_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("k150", "finite-rejects", "half-width", "s_k<=1/2", "430/447", "17 half-width failures", "lambda=-25.0", "k=149", "adaptive", "monotone"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    first_failures = summary.get("first_half_width_failures", [])
    if not isinstance(first_failures, list) or not first_failures:
        issues.append(issue("summary", "missing-first-half-width-failures", repr(first_failures)))
    elif first_failures[0].get("lam") != "-25.0" or first_failures[0].get("first_failure_k") != 133:
        issues.append(issue("summary", "bad-first-half-width-failure", repr(first_failures[0])))
    for key in ("min_half_width_margin", "max_scaled_defect", "min_one_third_margin"):
        extremum = summary.get(key, {})
        if not isinstance(extremum, dict) or not {"sample", "lam", "k"} <= set(extremum):
            issues.append(issue("summary", f"bad-{key}", repr(extremum)))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("no row", "finite_rejected_target", "half-width", "monotone defect", "one-third", "scaled-defect", "lambda <= 0"):
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
        "half-width theorem is proved",
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
            print(f"JWPF-NEG-LAMBDA-HALF-WIDTH-TAIL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda half-width tail target: "
            f"{summary.get('target_rows')} rows, "
            f"{len(issues)} issues, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows, "
            f"{summary.get('live_routes')} live routes, "
            f"{summary.get('finite_half_width_rows')} half-width rows, "
            f"{summary.get('finite_half_width_failure_rows')} half-width failures"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
