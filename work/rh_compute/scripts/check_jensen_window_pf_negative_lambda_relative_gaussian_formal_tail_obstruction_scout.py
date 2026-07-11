#!/usr/bin/env python3
"""Validate the relative-Gaussian formal-tail obstruction scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgftos_01_formal_residual_terms",
    "nlrgftos_02_degree240_coefficient_profile",
    "nlrgftos_03_finite_window_decay_until_least_term",
    "nlrgftos_04_turnaround_after_j103",
    "nlrgftos_05_monotone_geometric_tail_rejected",
    "nlrgftos_06_fixed_radius_cauchy_termwise_rejected",
    "nlrgftos_07_required_remainder_theorem",
    "nlrgftos_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_diagnostic",
    "finite_obstruction",
    "rejected_route",
    "exact_obstruction",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Formal-Tail Obstruction Scout",
    "Status: finite theorem-search obstruction",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout`",
    "templates. It does not reject or prove an actual contour",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: 8 rows, 0 issues, 4 profile rows, 4 formal-tail turnaround rows, 0 ready-to-apply rows",
    "max Taylor degree: 240",
    "tested j range: 21..120",
    "value budget A: 5.382819486765314521E-01",
    "derivative budget B: 9.315354075509573936E-03",
    "F_24: least value j=103",
    "first value growth after j>=80: j=103 ratio=[4.49075370734479227 +/- 8.55e-19]",
    "first derivative growth after j>=80: j=103 ratio=[4.53435324974143468 +/- 2.91e-18]",
    "A monotone decreasing or fixed-ratio geometric formal tail from",
    "fixed-radius Cauchy coefficient estimate",
    "This does not close the degree-40 residual budget.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class FormalTailObstructionIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> FormalTailObstructionIssue:
    return FormalTailObstructionIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[FormalTailObstructionIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[FormalTailObstructionIssue]:
    issues: list[FormalTailObstructionIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search obstruction":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree40_residual_tail_budget",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite",
        "rejects naive formal-tail",
        "does not bound",
        "actual residual",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[FormalTailObstructionIssue]:
    params = artifact.get("matrix_rows", [{}, {}])[1].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            int(params.get("max_taylor_degree", 240)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 256)),
            int(params.get("tail_start_k", 22)),
            int(params.get("collar_start_T", 1156)),
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[FormalTailObstructionIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[FormalTailObstructionIssue]:
    diagnostics = artifact.get("matrix_rows", [{}, {}])[1].get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    issues: list[FormalTailObstructionIssue] = []
    expected_params = {
        "max_taylor_degree": 240,
        "last_tested_j": 120,
        "tail_cutoff_n": 80,
        "precision_bits": 256,
        "tail_start_k": 22,
        "collar_start_T": 1156,
        "real_interval_u": "[0, 1/1156]",
        "profile_indices": [21, 22, 23, 24],
        "formal_residual_first_j": 21,
        "value_budget_A": "5.382819486765314521E-01",
        "derivative_budget_B": "9.315354075509573936E-03",
    }
    if params != expected_params:
        issues.append(issue("parameters", "bad-parameters", repr(params)))
    expected_scalars = {
        "ratio_ball_rows": 121,
        "highest_ratio_degree": 240,
        "profile_row_count": 4,
        "formal_tail_turnaround_rows": 4,
        "all_first_value_growth_from_j80_at_j103": True,
        "all_first_derivative_growth_from_j80_at_j103": True,
        "max_first_value_growth_ratio_from_j80": "[4.49075370734479227 +/- 8.56e-19]",
        "max_first_derivative_growth_ratio_from_j80": "[4.53435324974143468 +/- 2.92e-18]",
    }
    for key, value in expected_scalars.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))
    rows = diagnostics.get("profile_rows", [])
    if [row.get("index") for row in rows] != [21, 22, 23, 24]:
        issues.append(issue("profile_rows", "bad-index-order", repr([row.get("index") for row in rows])))
    for row in rows:
        label = f"F_{row.get('index')}"
        if row.get("term_j_range") != "21..120":
            issues.append(issue(label, "bad-term-j-range", repr(row.get("term_j_range"))))
        for key in (
            "least_value_term_j",
            "least_derivative_term_j",
            "first_value_growth_j_from_80",
            "first_derivative_growth_j_from_80",
        ):
            if row.get(key) != 103:
                issues.append(issue(label, f"bad-{key}", repr(row.get(key))))
        if row.get("formal_monotone_tail_from_j80_rejected") is not True:
            issues.append(issue(label, "missing-formal-tail-rejection", repr(row.get("formal_monotone_tail_from_j80_rejected"))))
        for key in (
            "value_sum_j21_to_j40",
            "value_sum_j41_to_j80",
            "value_sum_j81_to_j120",
            "derivative_sum_j21_to_j40",
            "derivative_sum_j41_to_j80",
            "derivative_sum_j81_to_j120",
        ):
            if not isinstance(row.get(key), str) or "+/-" not in row.get(key):
                issues.append(issue(label, f"bad-{key}", repr(row.get(key))))
        boundary = str(row.get("proof_boundary", "")).lower()
        if "finite formal-term diagnostic" not in boundary or "does not bound" not in boundary:
            issues.append(issue(label, "weak-row-boundary", boundary))
    obstruction = str(diagnostics.get("fixed_radius_cauchy_obstruction", "")).lower()
    for required in ("fixed-radius cauchy", "(a+j)", "tends to infinity", "cannot be summed"):
        if required not in obstruction:
            issues.append(issue("fixed_radius_cauchy_obstruction", "weak-obstruction-text", required))
    note = str(diagnostics.get("proof_boundary_note", "")).lower()
    for required in ("rejects naive", "does not reject", "actual", "remainder theorem"):
        if required not in note:
            issues.append(issue("diagnostics", "weak-proof-boundary-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[FormalTailObstructionIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[FormalTailObstructionIssue] = []
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
        if role in {"exact_reduction", "exact_obstruction"}:
            exact_rows += 1
        if role in {"finite_diagnostic", "finite_obstruction"}:
            finite_rows += 1
        if role in {"live_route", "rejected_route"}:
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
        if not any(marker in boundary for marker in ("only", "not", "rejected", "proof-hygiene", "method obstruction")):
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
) -> list[FormalTailObstructionIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "profile_rows": 4,
        "formal_tail_turnaround_rows": 4,
        "last_tested_j": 120,
        "max_taylor_degree": 240,
        "max_first_value_growth_ratio_from_j80": "[4.49075370734479227 +/- 8.56e-19]",
        "max_first_derivative_growth_ratio_from_j80": "[4.53435324974143468 +/- 2.92e-18]",
        "least_value_term_js": {"F_21": 103, "F_22": 103, "F_23": 103, "F_24": 103},
        "least_derivative_term_js": {"F_21": 103, "F_22": 103, "F_23": 103, "F_24": 103},
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[FormalTailObstructionIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_rows != 2:
        issues.append(issue("summary", "bad-exact-row-count", str(exact_rows)))
    if finite_rows != 3:
        issues.append(issue("summary", "bad-finite-row-count", str(finite_rows)))
    if route_rows != 2:
        issues.append(issue("summary", "bad-route-row-count", str(route_rows)))
    if acceptance_rows != 1:
        issues.append(issue("summary", "bad-acceptance-row-count", str(acceptance_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "degree-240",
        "naive infinite formal-tail",
        "degree-40 residual budget",
        "j=103",
        "j=104",
        "fixed-radius cauchy",
        "actual asymptotic remainder theorem",
        "not an infinite absolute sum",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "formal asymptotic", "geometric", "degree-40 residual budget", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[FormalTailObstructionIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[FormalTailObstructionIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "residual estimate is proved",
        "residual estimates are proved",
        "actual remainder theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[FormalTailObstructionIssue], dict]:
    artifact = load_json(target_path)
    issues: list[FormalTailObstructionIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-FORMAL-TAIL-OBSTRUCTION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('profile_rows')} profile rows, "
            f"{summary.get('formal_tail_turnaround_rows')} formal-tail turnaround rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
