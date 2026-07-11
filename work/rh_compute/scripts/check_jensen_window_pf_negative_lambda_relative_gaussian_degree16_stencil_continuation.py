#!/usr/bin/env python3
"""Validate the relative-Gaussian degree-16 stencil continuation diagnostic."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgd16sc_01_degree16_coefficient",
    "nlrgd16sc_02_all_positive_baselines_tested",
    "nlrgd16sc_03_pointwise_budget_still_fails",
    "nlrgd16sc_04_stencil_survival_filter",
    "nlrgd16sc_05_collar_signal",
    "nlrgd16sc_06_degree16_promotion_rejected",
    "nlrgd16sc_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "finite_certificate",
    "finite_diagnostic",
    "live_route",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Stencil Continuation",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: 7 rows, 0 issues, 4 tested continuation rows, 3 stencil-sign-preserving rows, 1 stencil-sign-failure rows, 0 ready-to-apply rows",
    "degree: 16",
    "sign: negative",
    "tested continuation rows: 4",
    "pointwise budget failures: 4",
    "stencil-sign-preserving rows: 3",
    "stencil-sign-failure rows: 1",
    "degree-14 baseline survivors: 1",
    "degree-14 baseline failures: 1",
    "nlrgts_M7_T1000: M->8",
    "stencil signs preserved=False",
    "failing stencils=['companion']",
    "nlrgts_M7_T2000: M->8",
    "stencil signs preserved=True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
)


@dataclass(frozen=True)
class Degree16ContinuationIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> Degree16ContinuationIssue:
    return Degree16ContinuationIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[Degree16ContinuationIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[Degree16ContinuationIssue]:
    issues: list[Degree16ContinuationIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_pointwise_tail_budget",
        "source_next_increment_stress",
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


def validate_recomputed(artifact: dict) -> list[Degree16ContinuationIssue]:
    params = artifact.get("matrix_rows", [{}])[1].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            int(params.get("baseline_max_taylor_degree", 14)),
            int(params.get("coefficient_max_taylor_degree", 16)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 256)),
            int(params.get("tail_start_k", 22)),
            [int(value) for value in params.get("sample_T_values", [25, 50, 100, 200, 500, 1000, 2000])],
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[Degree16ContinuationIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[Degree16ContinuationIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[1].get("diagnostics", {})
    expected_scalars = {
        "tested_continuation_rows": 4,
        "pointwise_budget_failure_rows": 4,
        "stencil_sign_preserving_rows": 3,
        "stencil_sign_failure_rows": 1,
        "half_safety_stencil_failure_rows": 3,
        "degree14_baseline_rows": 2,
        "degree14_baseline_survivors": 1,
        "degree14_baseline_failures": 1,
        "surviving_rows": ["nlrgts_M5_T2000", "nlrgts_M6_T2000", "nlrgts_M7_T2000"],
        "failed_rows": ["nlrgts_M7_T1000"],
        "worst_pointwise_over_budget": {"sample": "3.010798908654295615E+3", "source_row": "nlrgts_M6_T2000"},
        "worst_stencil_abs_over_half_margin": {
            "sample": "3.998025843926772743E+0",
            "source_row": "nlrgts_M6_T2000",
            "stencil": "companion",
        },
    }
    issues: list[Degree16ContinuationIssue] = []
    for key, value in expected_scalars.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))
    coefficient = diagnostics.get("degree16_coefficient", {})
    if coefficient.get("degree") != 16:
        issues.append(issue("degree16", "bad-degree", repr(coefficient.get("degree"))))
    if coefficient.get("sign") != "negative":
        issues.append(issue("degree16", "bad-sign", repr(coefficient.get("sign"))))
    if "[-237753170.9494156506565064820654116030884" not in str(coefficient.get("ratio_to_c0", "")):
        issues.append(issue("degree16", "bad-ratio-prefix", repr(coefficient.get("ratio_to_c0"))))
    rows = diagnostics.get("continuation_rows", [])
    if [row.get("source_row") for row in rows] != [
        "nlrgts_M5_T2000",
        "nlrgts_M6_T2000",
        "nlrgts_M7_T1000",
        "nlrgts_M7_T2000",
    ]:
        issues.append(issue("diagnostics", "bad-row-order", repr([row.get("source_row") for row in rows])))
    by_id = {row.get("source_row"): row for row in rows}
    failure = by_id.get("nlrgts_M7_T1000", {})
    if failure.get("stencil_signs_preserved") is not False:
        issues.append(issue("nlrgts_M7_T1000", "expected-sign-failure", repr(failure.get("stencil_signs_preserved"))))
    if failure.get("failing_stencils") != ["companion"]:
        issues.append(issue("nlrgts_M7_T1000", "bad-failing-stencils", repr(failure.get("failing_stencils"))))
    survivor = by_id.get("nlrgts_M7_T2000", {})
    if survivor.get("stencil_signs_preserved") is not True:
        issues.append(issue("nlrgts_M7_T2000", "expected-sign-survival", repr(survivor.get("stencil_signs_preserved"))))
    if survivor.get("half_safety_stencils_satisfied") is not True:
        issues.append(issue("nlrgts_M7_T2000", "expected-half-safety-survival", repr(survivor.get("half_safety_stencils_satisfied"))))
    for row in rows:
        row_id = str(row.get("source_row", "<missing-source>"))
        if row.get("pointwise_budget_satisfied") is not False:
            issues.append(issue(row_id, "expected-pointwise-failure", repr(row.get("pointwise_budget_satisfied"))))
        if len(row.get("stencil_increments", [])) != 3:
            issues.append(issue(row_id, "bad-stencil-count", repr(len(row.get("stencil_increments", [])))))
    return issues


def validate_rows(artifact: dict) -> tuple[list[Degree16ContinuationIssue], int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[Degree16ContinuationIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    finite_diagnostics = 0
    live_routes = 0
    rejected_routes = 0
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
        if row.get("role") == "finite_diagnostic":
            finite_diagnostics += 1
        if row.get("role") == "live_route":
            live_routes += 1
        if row.get("role") == "rejected_route":
            rejected_routes += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "not", "only", "live", "hygiene", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), finite_diagnostics, live_routes, rejected_routes if ready_to_apply == 0 else -rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    finite_diagnostics: int,
    live_routes: int,
    rejected_routes: int,
) -> list[Degree16ContinuationIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 7,
        "tested_continuation_rows": 4,
        "pointwise_budget_failure_rows": 4,
        "stencil_sign_preserving_rows": 3,
        "stencil_sign_failure_rows": 1,
        "half_safety_stencil_failure_rows": 3,
        "degree14_baseline_rows": 2,
        "degree14_baseline_survivors": 1,
        "degree14_baseline_failures": 1,
        "worst_pointwise_over_budget": "3.010798908654295615E+3",
        "worst_stencil_abs_over_half_margin": "3.998025843926772743E+0",
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[Degree16ContinuationIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 7:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if finite_diagnostics != 3:
        issues.append(issue("summary", "bad-finite-diagnostic-count", str(finite_diagnostics)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    if rejected_routes != 1:
        issues.append(issue("summary", "bad-rejected-route-count", str(rejected_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("degree-16", "three preserve", "t=1000", "t=2000", "large-t/q-over-t collar"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "not promoted", "t=1000", "t=2000", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[Degree16ContinuationIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[Degree16ContinuationIssue] = []
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


def validate(target_path: Path, note_path: Path) -> tuple[list[Degree16ContinuationIssue], dict]:
    artifact = load_json(target_path)
    issues: list[Degree16ContinuationIssue] = []
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-DEG16-STENCIL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('tested_continuation_rows')} tested continuation rows, "
            f"{summary.get('stencil_sign_preserving_rows')} stencil-sign-preserving rows, "
            f"{summary.get('stencil_sign_failure_rows')} stencil-sign-failure rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
