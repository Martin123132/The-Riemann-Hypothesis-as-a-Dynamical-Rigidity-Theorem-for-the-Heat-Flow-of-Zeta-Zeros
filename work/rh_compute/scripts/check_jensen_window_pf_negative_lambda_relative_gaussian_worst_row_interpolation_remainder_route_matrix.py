#!/usr/bin/env python3
"""Validate the worst-row interpolation-remainder route matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix import (  # noqa: E402
    DEFAULT_ARB_SCOUT_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_QUADRATURE_ROUTE_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwrirm_01_arb_interpolant_import",
    "nlrgwrirm_02_panel_mass_allocation",
    "nlrgwrirm_03_bernstein_sufficient_condition",
    "nlrgwrirm_04_heaviest_panel_budget_table",
    "nlrgwrirm_05_minimal_degree_table",
    "nlrgwrirm_06_endpoint_and_branch_condition",
    "nlrgwrirm_07_arb_cauchy_delta_promotion_rejected",
    "nlrgwrirm_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "scope_reduction",
    "budget_diagnostic",
    "sufficient_condition",
    "route_budget",
    "route_constraint",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "diagnostic_only",
    "open_requirement",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Interpolation-Remainder Route Matrix",
    "Status: worst-row interpolation-remainder route matrix",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian worst-row interpolation-remainder route matrix: 8 rows, 0 issues, 6 panel masses, 20 Bernstein budgets, 16 minimal-degree rows, 0 ready-to-apply rows",
    "heaviest panel: 20<=y<=50",
    "heaviest panel mass upper: 0.60216159332489531950482867922076390781009134801886",
    "rho=2 degree=128 value sup budget: 1.596890001110721912195980431156E-2",
    "rho=2 degree=160 value sup budget: 6.858590330079954287832319494469E+7",
    "panel_error <= gamma_panel_mass * 4*rho^(-N)/(rho-1) * M",
    "Arb Cauchy deltas do not prove the interpolation remainder.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class InterpolationRouteIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> InterpolationRouteIssue:
    return InterpolationRouteIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[InterpolationRouteIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[InterpolationRouteIssue]:
    issues: list[InterpolationRouteIssue] = []
    if (
        artifact.get("kind")
        != "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix"
    ):
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row interpolation-remainder route matrix":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_arb_chebyshev_interpolant_moment_scout",
        "source_arb_chebyshev_interpolant_moment_scout_json",
        "source_quadrature_remainder_route_matrix",
        "source_quadrature_remainder_route_json",
        "source_floating_chebyshev_panel_moment_scout",
        "source_compact_interval_integration_scout",
        "source_far_tail_split_certificate",
        "source_finite_part_weighted_sum_interval_certificate",
        "source_intervalization_target",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "route matrix only",
        "berstein",
        "does not prove an analytic",
        "does not prove a sup-norm",
        "does not prove a taylor-model",
        "does not prove the true compact integral",
        "finite-grid interval certificate",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            # Keep the common misspelling path off the artifact while still requiring the word below.
            if required == "berstein" and "bernstein" in boundary:
                continue
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[InterpolationRouteIssue]:
    try:
        recomputed = build_artifact(DEFAULT_ARB_SCOUT_JSON, DEFAULT_QUADRATURE_ROUTE_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[InterpolationRouteIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[InterpolationRouteIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[InterpolationRouteIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    ready_to_apply = 0
    open_requirements = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary", "source_artifacts"):
            if key not in row or row[key] in (None, ""):
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") not in ALLOWED_READINESS:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        if row.get("readiness") == "open_requirement":
            open_requirements += 1
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "does not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    rejected = by_id.get("nlrgwrirm_07_arb_cauchy_delta_promotion_rejected", {})
    if rejected.get("diagnostics", {}).get("arb_cauchy_delta_is_proof") is not False:
        issues.append(issue("nlrgwrirm_07_arb_cauchy_delta_promotion_rejected", "bad-rejection-flag", repr(rejected)))
    if "do not prove an interpolation remainder" not in str(rejected.get("gap", "")).lower():
        issues.append(issue("nlrgwrirm_07_arb_cauchy_delta_promotion_rejected", "weak-gap", repr(rejected.get("gap"))))
    if open_requirements < 2:
        issues.append(issue("matrix_rows", "missing-open-requirements", str(open_requirements)))
    return issues, len(rows), open_requirements


def validate_summary(artifact: dict, row_count: int) -> list[InterpolationRouteIssue]:
    summary = artifact.get("summary", {})
    issues: list[InterpolationRouteIssue] = []
    expected = {
        "matrix_rows": 8,
        "T": 10000,
        "index": 21,
        "alpha": "20.5",
        "compact_interval": "0<=y<=200",
        "panel_count": 6,
        "heaviest_panel": "20<=y<=50",
        "heaviest_panel_mass_upper": "0.60216159332489531950482867922076390781009134801886",
        "degree_count": 5,
        "rho_count": 4,
        "bernstein_budget_row_count": 20,
        "minimal_degree_row_count": 16,
        "value_unscaled_expectation_error_cap": "6.782032247872604818E-40",
        "derivative_unscaled_expectation_error_cap": "1.424226772053247012E-38",
        "rho_2_degree_128_value_sup_budget": "1.596890001110721912195980431156E-2",
        "rho_2_degree_160_value_sup_budget": "6.858590330079954287832319494469E+7",
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "arb_cauchy_delta_is_proof": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if dec(summary.get("heaviest_panel_mass_upper", "0")) <= Decimal("0.6"):
        issues.append(issue("summary", "bad-heaviest-mass", repr(summary.get("heaviest_panel_mass_upper"))))
    if dec(summary.get("rho_2_degree_128_value_sup_budget", "0")) >= Decimal("1"):
        issues.append(issue("summary", "degree-128-budget-too-large", repr(summary.get("rho_2_degree_128_value_sup_budget"))))
    if dec(summary.get("rho_2_degree_160_value_sup_budget", "0")) <= Decimal("1e6"):
        issues.append(issue("summary", "degree-160-budget-too-small", repr(summary.get("rho_2_degree_160_value_sup_budget"))))
    budget_rows = None
    minimal_rows = None
    for row in artifact.get("matrix_rows", []):
        if isinstance(row, dict) and row.get("id") == "nlrgwrirm_04_heaviest_panel_budget_table":
            budget_rows = row.get("diagnostics", {}).get("bernstein_budget_rows")
        if isinstance(row, dict) and row.get("id") == "nlrgwrirm_05_minimal_degree_table":
            minimal_rows = row.get("diagnostics", {}).get("minimal_degree_rows")
    if not isinstance(budget_rows, list) or len(budget_rows) != 20:
        issues.append(issue("summary", "bad-budget-rows", repr(budget_rows)))
        budget_rows = []
    if not isinstance(minimal_rows, list) or len(minimal_rows) != 16:
        issues.append(issue("summary", "bad-minimal-degree-rows", repr(minimal_rows)))
        minimal_rows = []
    by_min = {
        (row.get("assumed_joint_sup_norm"), row.get("rho")): row
        for row in minimal_rows
        if isinstance(row, dict)
    }
    expected_minimal = {
        ("1e-12", "2.0"): 95,
        ("1e-6", "2.0"): 115,
        ("1", "2.0"): 134,
        ("1e6", "2.0"): 154,
        ("1e6", "3.0"): 97,
    }
    for key, expected_degree in expected_minimal.items():
        if by_min.get(key, {}).get("minimal_degree_in_16_to_256") != expected_degree:
            issues.append(issue("summary", "bad-minimal-degree", f"{key}: {by_min.get(key)!r}"))
    none_125 = [
        row
        for row in minimal_rows
        if isinstance(row, dict)
        and row.get("rho") == "1.25"
        and row.get("minimal_degree_in_16_to_256") is not None
    ]
    if none_125:
        issues.append(issue("summary", "rho-1.25-unexpectedly-closed", repr(none_125)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "interpolation-remainder route matrix",
        "heaviest gamma panel",
        "bernstein route formula",
        "degree 128",
        "degree 160",
        "route matrix only",
        "no analytic-domain",
        "true interpolation-remainder",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    required_input = str(summary.get("required_missing_input", "")).lower()
    for required in ("analytic bernstein-ellipse domain", "sup norm", "value and derivative cores"):
        if required not in required_input:
            issues.append(issue("summary", "weak-required-missing-input", required))
    endpoint = str(summary.get("endpoint_route_condition", "")).lower()
    for required in ("y=0", "sqrt", "x-variable", "even analytic y-core"):
        if required not in endpoint:
            issues.append(issue("summary", "weak-endpoint-condition", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in (
        "ready_to_apply",
        "sufficient-condition",
        "cauchy deltas",
        "endpoint y=0",
        "compact interval",
        "lambda <= 0",
    ):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[InterpolationRouteIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[InterpolationRouteIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "interpolation-remainder theorem is proved",
        "compact interval-integration certificate is complete",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
        "rh is proved",
        "lambda <= 0 is proved",
        "arb cauchy deltas prove the interpolation remainder",
        "analytic-domain bound is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-promotion-language", forbidden))
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact_path = args.artifact if args.artifact.is_absolute() else REPO_ROOT / args.artifact
    note_path = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = load_json(artifact_path)
    issues: list[InterpolationRouteIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, _open_requirements = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"INTERPOLATION-REMAINDER-ROUTE {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
