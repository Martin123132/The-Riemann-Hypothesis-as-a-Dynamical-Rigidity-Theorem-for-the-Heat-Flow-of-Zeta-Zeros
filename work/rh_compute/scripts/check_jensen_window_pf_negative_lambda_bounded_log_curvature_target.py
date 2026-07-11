#!/usr/bin/env python3
"""Validate the negative-lambda bounded log-curvature theorem target."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_bounded_log_curvature_target import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_diagnostics,
)


REQUIRED_ROW_IDS = {
    "nllct_01_bounded_log_curvature_statement",
    "nllct_02_raw_moment_curvature_equivalence",
    "nllct_03_saddle_laplace_route",
    "nllct_04_moment_ratio_growth_route",
    "nllct_05_plain_moment_logconvexity_insufficient",
    "nllct_06_finite_compatibility",
    "nllct_07_monotone_contraction_interaction",
    "nllct_08_conditional_defect_tail_application",
}

ALLOWED_ROLES = {
    "open_statement",
    "exact_reduction",
    "live_route",
    "insufficient_route",
    "finite_diagnostic",
    "open_dependency",
    "conditional_application",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Bounded Log-Curvature Target",
    "Status: open theorem target",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_bounded_log_curvature_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_bounded_log_curvature_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_bounded_log_curvature_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_bounded_log_curvature_target.py",
    "validated Jensen-window PF negative-lambda bounded log-curvature target: 8 rows, 0 issues, 0 ready-to-apply rows, 2 live routes, 63 raw-threshold rows",
    "B_k = -Delta^2 log A_k = -log(x_k)",
    "A_k = M_k*k!/(2*k)!",
    "0 <= B_k <= 2/(3*(2*k+1))",
    "x_k = ((2*k-1)/(2*k+1)) * (M_(k+1)*M_(k-1)/M_k^2)",
    "Ordinary moment log-convexity gives only the nonnegative left side",
    "uniform saddle/Laplace",
    "raw moment-ratio growth",
    "outputs/jensen_window_pf_negative_lambda_log_curvature_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
)


def issue(section: str, name: str, detail: str) -> dict:
    return {"section": section, "issue": name, "detail": detail}


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[dict]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(target: dict) -> list[dict]:
    issues: list[dict] = []
    if target.get("kind") != "jensen_window_pf_negative_lambda_bounded_log_curvature_target":
        issues.append(issue("<target>", "bad-kind", repr(target.get("kind"))))
    if target.get("status") != "open_theorem_target":
        issues.append(issue("<target>", "bad-status", repr(target.get("status"))))
    if target.get("target_id") != "target_negative_lambda_bounded_log_curvature":
        issues.append(issue("<target>", "bad-target-id", repr(target.get("target_id"))))
    for key in (
        "source_log_curvature_bridge",
        "source_defect_tail_target",
        "source_monotone_contraction_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<target>", target.get(key)))
    for ref in target.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<target>", ref))
    boundary = str(target.get("proof_boundary", "")).lower()
    for required in ("open theorem target", "does not prove", "bounded log-curvature", "defect-tail", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<target>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(target: dict) -> list[dict]:
    issues: list[dict] = []
    paths = [REPO_ROOT / ref for ref in target.get("enclosure_jsonl", [])]
    try:
        diagnostics = asdict(build_diagnostics(paths))
    except Exception as exc:
        return [issue("finite_diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    recorded = target.get("finite_diagnostics", {})
    for key, value in diagnostics.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))
    expected = {
        "coefficient_k_max": 22,
        "checked_x_max": 21,
        "raw_logconvexity_rows": 63,
        "raw_logconvexity_positive_rows": 63,
        "raw_curvature_threshold_rows": 63,
        "raw_curvature_threshold_positive_rows": 63,
        "simple_log_buffer_rows": 63,
        "simple_log_buffer_positive_rows": 63,
        "threshold_positive_rows": 63,
    }
    for key, value in expected.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"unexpected-{key}", f"{recorded.get(key)!r} != {value!r}"))
    return issues


def validate_rows(target: dict) -> tuple[list[dict], int, int, int]:
    rows = target.get("target_rows", [])
    issues: list[dict] = []
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
        if row.get("role") in {"open_statement", "live_route", "open_dependency", "conditional_application"} and "not" not in text:
            issues.append(issue(row_id, "weak-boundary", text))
        if row.get("role") == "insufficient_route" and "insufficient" not in text and "too weak" not in text:
            issues.append(issue(row_id, "missing-insufficient-language", text))
    return issues, len(rows), ready_count, live_count


def validate_summary(target: dict, row_count: int, ready_count: int, live_count: int) -> list[dict]:
    summary = target.get("summary", {})
    expected = {
        "target_rows": 8,
        "ready_to_apply_rows": 0,
        "live_routes": 2,
        "exact_reduction_rows": 2,
        "insufficient_routes": 1,
        "conditional_application_rows": 1,
        "raw_curvature_threshold_rows": 63,
        "simple_log_buffer_rows": 63,
        "target_closing": False,
    }
    issues: list[dict] = []
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
    for required in ("bounded log-curvature", "b_k", "raw moments", "plain moment log-convexity", "too weak"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in target.get("invariants", [])).lower()
    for required in ("no row", "open_theorem_target", "plain raw moment", "finite prefix", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[dict]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[dict] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "defect-tail theorem is proved",
        "bounded log-curvature estimate is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[dict], dict]:
    target = load_json(target_path)
    issues: list[dict] = []
    issues.extend(validate_top_level(target))
    issues.extend(validate_recomputed(target))
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
        print(json.dumps({"ok": ok, "summary": summary, "issues": issues}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-BOUNDED-LOG-CURVATURE {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda bounded log-curvature target: "
            f"{summary.get('target_rows')} rows, {len(issues)} issues, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows, "
            f"{summary.get('live_routes')} live routes, "
            f"{summary.get('raw_curvature_threshold_rows')} raw-threshold rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
