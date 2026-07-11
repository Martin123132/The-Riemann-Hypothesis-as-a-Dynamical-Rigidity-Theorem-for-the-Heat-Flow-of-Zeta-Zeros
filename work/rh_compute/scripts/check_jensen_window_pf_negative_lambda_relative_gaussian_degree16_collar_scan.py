#!/usr/bin/env python3
"""Validate the relative-Gaussian degree-16 collar scan diagnostic."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgd16cs_01_integer_T_scan",
    "nlrgd16cs_02_baseline_threshold",
    "nlrgd16cs_03_continuation_threshold",
    "nlrgd16cs_04_half_safety_threshold",
    "nlrgd16cs_05_pointwise_route_blocked_on_scan",
    "nlrgd16cs_06_live_collar_theorem_target",
    "nlrgd16cs_07_scan_promotion_rejected",
}

ALLOWED_ROLES = {
    "finite_diagnostic",
    "finite_threshold",
    "finite_obstruction",
    "live_route",
    "rejected_route",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Collar Scan",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: 7 rows, 0 issues, 1301 scan rows, 1045 continuation-positive rows, 718 half-safety rows, 0 ready-to-apply rows",
    "scan rows: 1301",
    "baseline-positive rows: 1283",
    "continuation-positive rows: 1045",
    "half-safety rows: 718",
    "pointwise budget successes: 0",
    "pointwise budget failures: 1283",
    "first baseline-positive T: 918",
    "first continuation-positive T: 1156",
    "first continuation q/T: 1.946366782006920415E-2",
    "first half-safety T: 1483",
    "first half-safety q/T: 1.517194875252865813E-2",
    "worst pointwise over-budget: 2.638097770251728679E+5 at T=918",
    "worst stencil abs over half-margin: 7.209437565352806308E+2 at T=918",
    "T=1155: baseline=True, continuation=False, half_safety=False",
    "T=1156: baseline=True, continuation=True, half_safety=False",
    "T=1482: baseline=True, continuation=True, half_safety=False",
    "T=1483: baseline=True, continuation=True, half_safety=True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class Degree16CollarScanIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> Degree16CollarScanIssue:
    return Degree16CollarScanIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[Degree16CollarScanIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[Degree16CollarScanIssue]:
    issues: list[Degree16CollarScanIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree16_continuation",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite",
        "does not prove",
        "real-t collar",
        "infinite taylor-tail",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[Degree16CollarScanIssue]:
    params = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            int(params.get("coefficient_max_taylor_degree", 16)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 256)),
            int(params.get("tail_start_k", 22)),
            int(params.get("scan_start_T", 900)),
            int(params.get("scan_end_T", 2200)),
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[Degree16CollarScanIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_expected_row(
    section: str,
    row: dict,
    expected: dict,
) -> list[Degree16CollarScanIssue]:
    issues: list[Degree16CollarScanIssue] = []
    for key, value in expected.items():
        if row.get(key) != value:
            issues.append(issue(section, f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
    return issues


def validate_diagnostics(artifact: dict) -> list[Degree16CollarScanIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {})
    expected = {
        "parameters": {
            "baseline_M": 7,
            "coefficient_max_taylor_degree": 16,
            "continuation_M": 8,
            "precision_bits": 256,
            "scan_end_T": 2200,
            "scan_start_T": 900,
            "tail_cutoff_n": 80,
            "tail_start_k": 22,
        },
        "scan_rows": 1301,
        "baseline_positive_rows": 1283,
        "continuation_positive_rows": 1045,
        "half_safety_rows": 718,
        "pointwise_budget_success_rows": 0,
        "pointwise_budget_failure_rows": 1283,
        "worst_pointwise_over_budget": {"sample": "2.638097770251728679E+5", "T": "918"},
        "worst_stencil_abs_over_half_margin": {"sample": "7.209437565352806308E+2", "T": "918"},
    }
    issues: list[Degree16CollarScanIssue] = []
    for key, value in expected.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))

    threshold_expectations = {
        "first_baseline_positive": {
            "T": "918",
            "q_over_T": "2.450980392156862745E-2",
            "baseline_positive": True,
            "continuation_positive": False,
            "half_safety_stencils_satisfied": False,
            "pointwise_budget_satisfied": False,
            "pointwise_over_budget": "2.638097770251728679E+5",
            "worst_stencil_abs_over_half_margin": "7.209437565352806308E+2",
        },
        "first_continuation_positive": {
            "T": "1156",
            "q_over_T": "1.946366782006920415E-2",
            "baseline_positive": True,
            "continuation_positive": True,
            "half_safety_stencils_satisfied": False,
            "pointwise_budget_satisfied": False,
            "companion_margin_M8": "8.851471227548859301E-10",
        },
        "first_half_safety": {
            "T": "1483",
            "q_over_T": "1.517194875252865813E-2",
            "baseline_positive": True,
            "continuation_positive": True,
            "half_safety_stencils_satisfied": True,
            "pointwise_budget_satisfied": False,
            "companion_margin_M8": "1.124347089699779384E-7",
            "worst_stencil_abs_over_half_margin": "9.996039982772186202E-1",
        },
    }
    for key, expected_row in threshold_expectations.items():
        row = diagnostics.get(key, {})
        issues.extend(validate_expected_row(key, row, expected_row))

    selected_rows = diagnostics.get("selected_rows", [])
    if not isinstance(selected_rows, list):
        return issues + [issue("diagnostics", "bad-selected-rows", repr(type(selected_rows)))]
    by_T = {row.get("T"): row for row in selected_rows if isinstance(row, dict)}
    selected_expectations = {
        "1155": {
            "continuation_positive": False,
            "half_safety_stencils_satisfied": False,
            "companion_margin_M8": "-7.079873820087283923E-10",
            "worst_stencil_abs_over_half_margin": "2.001306345621169903E+0",
        },
        "1156": {
            "continuation_positive": True,
            "half_safety_stencils_satisfied": False,
            "companion_margin_M8": "8.851471227548859301E-10",
            "worst_stencil_abs_over_half_margin": "1.998356122693134981E+0",
        },
        "1482": {
            "continuation_positive": True,
            "half_safety_stencils_satisfied": False,
            "companion_margin_M8": "1.125174810934742418E-7",
            "worst_stencil_abs_over_half_margin": "1.002251673746439061E+0",
        },
        "1483": {
            "continuation_positive": True,
            "half_safety_stencils_satisfied": True,
            "companion_margin_M8": "1.124347089699779384E-7",
            "worst_stencil_abs_over_half_margin": "9.996039982772186202E-1",
        },
    }
    for T, expected_row in selected_expectations.items():
        if T not in by_T:
            issues.append(issue("selected_rows", "missing-T", T))
            continue
        if by_T[T].get("baseline_positive") is not True:
            issues.append(issue(f"T={T}", "expected-baseline-positive", repr(by_T[T].get("baseline_positive"))))
        if by_T[T].get("pointwise_budget_satisfied") is not False:
            issues.append(issue(f"T={T}", "expected-pointwise-failure", repr(by_T[T].get("pointwise_budget_satisfied"))))
        issues.extend(validate_expected_row(f"T={T}", by_T[T], expected_row))
    return issues


def validate_rows(artifact: dict) -> tuple[list[Degree16CollarScanIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[Degree16CollarScanIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    finite_diagnostics = 0
    finite_thresholds = 0
    finite_obstructions = 0
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
        role = row.get("role")
        if role not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(role)))
        if role == "finite_diagnostic":
            finite_diagnostics += 1
        if role == "finite_threshold":
            finite_thresholds += 1
        if role == "finite_obstruction":
            finite_obstructions += 1
        if role == "live_route":
            live_routes += 1
        if role == "rejected_route":
            rejected_routes += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "not", "only", "live", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), finite_diagnostics, finite_thresholds, finite_obstructions, live_routes + rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    finite_diagnostics: int,
    finite_thresholds: int,
    finite_obstructions: int,
    route_rows: int,
) -> list[Degree16CollarScanIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 7,
        "scan_rows": 1301,
        "baseline_positive_rows": 1283,
        "continuation_positive_rows": 1045,
        "half_safety_rows": 718,
        "pointwise_budget_success_rows": 0,
        "pointwise_budget_failure_rows": 1283,
        "first_baseline_positive_T": "918",
        "first_continuation_positive_T": "1156",
        "first_half_safety_T": "1483",
        "first_continuation_q_over_T": "1.946366782006920415E-2",
        "first_half_safety_q_over_T": "1.517194875252865813E-2",
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[Degree16CollarScanIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 7:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if finite_diagnostics != 1:
        issues.append(issue("summary", "bad-finite-diagnostic-count", str(finite_diagnostics)))
    if finite_thresholds != 3:
        issues.append(issue("summary", "bad-finite-threshold-count", str(finite_thresholds)))
    if finite_obstructions != 1:
        issues.append(issue("summary", "bad-finite-obstruction-count", str(finite_obstructions)))
    if route_rows != 2:
        issues.append(issue("summary", "bad-route-row-count", str(route_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "t=900..2200",
        "t=918",
        "t=1156",
        "1.946366782006920415e-2",
        "t=1483",
        "1.517194875252865813e-2",
        "pointwise budget fails",
        "direct signed stencil-tail collar",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "integer-t", "pointwise budget", "k=22", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[Degree16CollarScanIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[Degree16CollarScanIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "analytic collar theorem is proved",
        "real-t collar theorem is proved",
        "infinite taylor-tail estimate is proved",
        "uniform taylor-tail theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[Degree16CollarScanIssue], dict]:
    artifact = load_json(target_path)
    issues: list[Degree16CollarScanIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, finite_diagnostics, finite_thresholds, finite_obstructions, route_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(
        validate_summary(
            artifact,
            row_count,
            finite_diagnostics,
            finite_thresholds,
            finite_obstructions,
            route_rows,
        )
    )
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-DEG16-COLLAR {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('scan_rows')} scan rows, "
            f"{summary.get('continuation_positive_rows')} continuation-positive rows, "
            f"{summary.get('half_safety_rows')} half-safety rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
