#!/usr/bin/env python3
"""Validate the relative-Gaussian asymptotic remainder target artifact."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target import (  # noqa: E402
    DEFAULT_FORMAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_RESIDUAL_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgart_01_budget_import",
    "nlrgart_02_first_omitted_multiplier_target",
    "nlrgart_03_safe_1000x_first_omitted_budget",
    "nlrgart_04_optimized_window_target",
    "nlrgart_05_bad_infinite_sum_rejected",
    "nlrgart_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "exact_sufficient_condition",
    "finite_diagnostic",
    "conditional_handoff",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Asymptotic Remainder Target",
    "Status: exact theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: 6 rows, 0 issues, 4 first-omitted rows, 4 optimized-window rows, 0 ready-to-apply rows",
    "common multiplier limit: [1419.93907869608 +/- 1.36e-12]",
    "safe common multiplier tested: 1000",
    "safe value budget fraction: [0.67090090649120 +/- 5.01e-15]",
    "safe derivative budget fraction: [0.704255566315065 +/- 6.49e-16]",
    "F_24: value term=[0.000361133847314933 +/- 4.01e-19]",
    "common least-term multiplier limit after window: [4.86406709155425e+29 +/- 6.04e+14]",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class AsymptoticTargetIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> AsymptoticTargetIssue:
    return AsymptoticTargetIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[AsymptoticTargetIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[AsymptoticTargetIssue]:
    issues: list[AsymptoticTargetIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree40_residual_tail_budget",
        "source_formal_tail_obstruction",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "exact theorem-search diagnostic",
        "calibrates sufficient",
        "does not prove",
        "remainder estimates",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict, residual_path: Path, formal_path: Path) -> list[AsymptoticTargetIssue]:
    try:
        recomputed = build_artifact(residual_path, formal_path)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[AsymptoticTargetIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[AsymptoticTargetIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[0].get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    issues: list[AsymptoticTargetIssue] = []
    expected_params = {
        "tail_start_k": 22,
        "indices": [21, 22, 23, 24],
        "collar_interval_u": "[0, 1/1156]",
        "degree40_residual_first_j": 21,
        "formal_window_j_range": "21..120",
        "least_term_j": 103,
        "value_budget_A": "5.382819486765314521E-01",
        "derivative_budget_B": "9.315354075509573936E-03",
        "safe_common_first_omitted_multiplier": "1000",
    }
    if params != expected_params:
        issues.append(issue("parameters", "bad-parameters", repr(params)))
    expected_scalars = {
        "first_omitted_row_count": 4,
        "optimized_window_row_count": 4,
        "common_first_omitted_multiplier_limit": "[1419.93907869608 +/- 1.36e-12]",
        "safe_common_multiplier_value_fraction": "[0.67090090649120 +/- 5.01e-15]",
        "safe_common_multiplier_derivative_fraction": "[0.704255566315065 +/- 6.49e-16]",
        "safe_common_multiplier_closes_first_omitted_target": True,
        "common_least_multiplier_limit_after_window": "[4.86406709155425e+29 +/- 6.04e+14]",
        "max_window_value_budget_fraction": "[0.00093155697477913 +/- 2.65e-18]",
        "max_window_derivative_budget_fraction": "[0.00099588669254223 +/- 1.83e-18]",
    }
    for key, value in expected_scalars.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))
    first_rows = diagnostics.get("first_omitted_multiplier_rows", [])
    opt_rows = diagnostics.get("optimized_window_multiplier_rows", [])
    if [row.get("index") for row in first_rows] != [21, 22, 23, 24]:
        issues.append(issue("first_omitted_rows", "bad-index-order", repr([row.get("index") for row in first_rows])))
    if [row.get("index") for row in opt_rows] != [21, 22, 23, 24]:
        issues.append(issue("optimized_window_rows", "bad-index-order", repr([row.get("index") for row in opt_rows])))
    for row in first_rows:
        label = f"F_{row.get('index')}"
        if row.get("limiting_channel") != "derivative":
            issues.append(issue(label, "bad-limiting-channel", repr(row.get("limiting_channel"))))
        if "not proved" not in str(row.get("proof_boundary", "")).lower():
            issues.append(issue(label, "weak-first-row-boundary", repr(row.get("proof_boundary"))))
    f24 = {row.get("index"): row for row in first_rows}.get(24, {})
    if f24.get("common_budget_multiplier_limit") != "[1419.939078696080 +/- 9.63e-13]":
        issues.append(issue("F_24", "bad-common-budget-limit", repr(f24.get("common_budget_multiplier_limit"))))
    for row in opt_rows:
        label = f"F_{row.get('index')}"
        if row.get("least_value_term_j") != 103 or row.get("least_derivative_term_j") != 103:
            issues.append(issue(label, "bad-least-term-j", repr(row)))
        if "conditional" not in str(row.get("proof_boundary", "")).lower():
            issues.append(issue(label, "weak-optimized-row-boundary", repr(row.get("proof_boundary"))))
    note = str(diagnostics.get("proof_boundary_note", "")).lower()
    for required in ("sufficient targets", "do not prove", "actual zeta heat-flow", "remainder"):
        if required not in note:
            issues.append(issue("diagnostics", "weak-proof-boundary-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[AsymptoticTargetIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[AsymptoticTargetIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_rows = 0
    finite_rows = 0
    route_rows = 0
    acceptance_rows = 0
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
        if role in {"exact_reduction", "exact_sufficient_condition"}:
            exact_rows += 1
        if role in {"finite_diagnostic", "conditional_handoff"}:
            finite_rows += 1
        if role == "rejected_route":
            route_rows += 1
        if role == "acceptance_gate":
            acceptance_rows += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") not in {"not_ready_to_apply", "available_exact"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "conditional", "rejected", "proof-hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), exact_rows, finite_rows, route_rows, acceptance_rows


def validate_summary(
    artifact: dict,
    row_count: int,
    exact_rows: int,
    finite_rows: int,
    route_rows: int,
    acceptance_rows: int,
) -> list[AsymptoticTargetIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 6,
        "first_omitted_rows": 4,
        "optimized_window_rows": 4,
        "common_first_omitted_multiplier_limit": "[1419.93907869608 +/- 1.36e-12]",
        "safe_common_multiplier": "1000",
        "safe_common_multiplier_value_fraction": "[0.67090090649120 +/- 5.01e-15]",
        "safe_common_multiplier_derivative_fraction": "[0.704255566315065 +/- 6.49e-16]",
        "safe_common_multiplier_closes_first_omitted_target": True,
        "common_least_multiplier_limit_after_window": "[4.86406709155425e+29 +/- 6.04e+14]",
        "max_window_value_budget_fraction": "[0.00093155697477913 +/- 2.65e-18]",
        "max_window_derivative_budget_fraction": "[0.00099588669254223 +/- 1.83e-18]",
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[AsymptoticTargetIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 6:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_rows != 2:
        issues.append(issue("summary", "bad-exact-row-count", str(exact_rows)))
    if finite_rows != 2:
        issues.append(issue("summary", "bad-finite-row-count", str(finite_rows)))
    if route_rows != 1:
        issues.append(issue("summary", "bad-route-row-count", str(route_rows)))
    if acceptance_rows != 1:
        issues.append(issue("summary", "bad-acceptance-row-count", str(acceptance_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "degree-40 residual budget",
        "1000x first-omitted-term",
        "half-safety",
        "f_21..f_24",
        "optimized-window",
        "real analytic remainder theorem",
        "rather than an infinite formal-tail sum",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "1000x", "conditional", "infinite formal", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[AsymptoticTargetIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[AsymptoticTargetIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "residual estimate is proved",
        "first-omitted-term theorem is proved",
        "actual remainder theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path, residual_path: Path, formal_path: Path) -> tuple[list[AsymptoticTargetIssue], dict]:
    artifact = load_json(target_path)
    issues: list[AsymptoticTargetIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, residual_path, formal_path))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, exact_rows, finite_rows, route_rows, acceptance_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_rows, finite_rows, route_rows, acceptance_rows))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--residual-json", type=Path, default=DEFAULT_RESIDUAL_JSON)
    parser.add_argument("--formal-json", type=Path, default=DEFAULT_FORMAL_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    residual_path = args.residual_json if args.residual_json.is_absolute() else REPO_ROOT / args.residual_json
    formal_path = args.formal_json if args.formal_json.is_absolute() else REPO_ROOT / args.formal_json
    issues, summary = validate(target, note, residual_path, formal_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-ASYMPTOTIC-REMAINDER {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian asymptotic remainder target: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('first_omitted_rows')} first-omitted rows, "
            f"{summary.get('optimized_window_rows')} optimized-window rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
