#!/usr/bin/env python3
"""Validate the relative-Gaussian quadrature-remainder route matrix."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix import (  # noqa: E402
    DEFAULT_FINITE_PLUS_TAIL_JSON,
    DEFAULT_FIRST_OMITTED_JSON,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_QUADRATURE_LADDER_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgqrrm_01_gauss_laguerre_remainder_formula",
    "nlrgqrrm_02_ratio_cap_to_derivative_sup_caps",
    "nlrgqrrm_03_finite_plus_tail_plus_quadrature_budget",
    "nlrgqrrm_04_cauchy_derivative_route",
    "nlrgqrrm_05_interval_adaptive_integration_route",
    "nlrgqrrm_06_order_spread_promotion_rejected",
    "nlrgqrrm_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_formula_target",
    "sufficient_condition",
    "budget_composition",
    "live_route",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_exact_formula",
    "open_requirement",
    "available_budget_arithmetic",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Quadrature-Remainder Route Matrix",
    "Status: quadrature-remainder route matrix",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian quadrature-remainder route matrix: 7 rows, 0 issues, derivative order 640, 2 derivative-sup caps, 0 ready-to-apply rows",
    "Laguerre error prefactor: 2.791244149880588131514918942024363250922057058755595558776957E-159",
    "value required 640th derivative sup cap: 2.4297524271256407631984588349963763207485E+119",
    "derivative required 640th derivative sup cap: 5.1024800969638456035049425453208614326121E+120",
    "value ratio upper after quadrature cap: 0.9853967992836557769419015895036210773888",
    "derivative ratio upper after quadrature cap: 0.9714065674762067320093698741711260875260",
    "both ratios still below one after quadrature cap: True",
    "Rejected route: do not promote floating N=96..320 order-spread into a remainder theorem.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class QuadratureRouteIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> QuadratureRouteIssue:
    return QuadratureRouteIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[QuadratureRouteIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[QuadratureRouteIssue]:
    issues: list[QuadratureRouteIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "quadrature-remainder route matrix":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_first_omitted_denominator_certificate",
        "source_first_omitted_denominator_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_quadrature_ladder_scout",
        "source_quadrature_ladder_json",
        "source_finite_plus_tail_budget_certificate",
        "source_finite_plus_tail_budget_json",
        "source_cancellation_reduced_grid_scout",
        "source_phi_tail_grid_certificate",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "route matrix only",
        "does not prove any derivative",
        "does not produce interval adaptive integration",
        "finite-grid interval certificate",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[QuadratureRouteIssue]:
    try:
        recomputed = build_artifact(
            DEFAULT_FIRST_OMITTED_JSON,
            DEFAULT_INTERVAL_JSON,
            DEFAULT_QUADRATURE_LADDER_JSON,
            DEFAULT_FINITE_PLUS_TAIL_JSON,
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[QuadratureRouteIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[QuadratureRouteIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[QuadratureRouteIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    live_routes = 0
    rejected_routes = 0
    ready_to_apply = 0
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
        if row.get("role") == "live_route":
            live_routes += 1
        if row.get("role") == "rejected_route":
            rejected_routes += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
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
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    rejected = rows_by_id.get("nlrgqrrm_06_order_spread_promotion_rejected", {}).get("diagnostics", {})
    if rejected.get("ladder_order_spread_is_proof") is not False:
        issues.append(issue("nlrgqrrm_06_order_spread_promotion_rejected", "bad-rejection-flag", repr(rejected)))
    return issues, len(rows), live_routes, rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    live_routes: int,
    rejected_routes: int,
) -> list[QuadratureRouteIssue]:
    summary = artifact.get("summary", {})
    issues: list[QuadratureRouteIssue] = []
    expected = {
        "matrix_rows": 7,
        "T": 10000,
        "index": 21,
        "quadrature_order": 320,
        "derivative_order": 640,
        "laguerre_error_prefactor": "2.791244149880588131514918942024363250922057058755595558776957E-159",
        "quadrature_ratio_radius_cap": "0.0000010",
        "intervalization_per_source_cap": "2.000000000000000000E-3",
        "quadrature_cap_below_per_source_cap": True,
        "value_unscaled_expectation_error_cap": "6.782032247872604818E-40",
        "derivative_unscaled_expectation_error_cap": "1.424226772053247012E-38",
        "value_required_640th_derivative_supremum_cap": "2.4297524271256407631984588349963763207485E+119",
        "derivative_required_640th_derivative_supremum_cap": "5.1024800969638456035049425453208614326121E+120",
        "finite_plus_tail_value_ratio_upper": "0.9853957992836557769419015895036210773888",
        "finite_plus_tail_derivative_ratio_upper": "0.9714055674762067320093698741711260875260",
        "value_ratio_upper_after_quadrature_cap": "0.9853967992836557769419015895036210773888",
        "derivative_ratio_upper_after_quadrature_cap": "0.9714065674762067320093698741711260875260",
        "both_ratios_still_below_one_after_quadrature_cap": True,
        "live_route_rows": 2,
        "rejected_route_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if live_routes != summary.get("live_route_rows"):
        issues.append(issue("summary", "live-route-count-mismatch", f"{live_routes} != {summary.get('live_route_rows')!r}"))
    if rejected_routes != summary.get("rejected_route_rows"):
        issues.append(
            issue("summary", "rejected-route-count-mismatch", f"{rejected_routes} != {summary.get('rejected_route_rows')!r}")
        )
    if not (dec(summary.get("laguerre_error_prefactor", "1")) < Decimal("1e-150")):
        issues.append(issue("summary", "prefactor-not-small", repr(summary.get("laguerre_error_prefactor"))))
    if not (dec(summary.get("value_ratio_upper_after_quadrature_cap", "2")) < Decimal(1)):
        issues.append(issue("summary", "value-ratio-not-below-one", repr(summary.get("value_ratio_upper_after_quadrature_cap"))))
    if not (dec(summary.get("derivative_ratio_upper_after_quadrature_cap", "2")) < Decimal(1)):
        issues.append(
            issue("summary", "derivative-ratio-not-below-one", repr(summary.get("derivative_ratio_upper_after_quadrature_cap")))
        )
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "classical n=320",
        "640th-derivative",
        "2.79e-159",
        "1e-6",
        "finite-plus-tail",
        "route matrix only",
        "no derivative",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "order-spread", "sufficient conditions", "node-only", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[QuadratureRouteIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[QuadratureRouteIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "quadrature-remainder theorem is proved",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
        "rh is proved",
        "lambda <= 0 is proved",
        "order-spread proves the quadrature remainder",
        "derivative bound is proved",
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
    issues: list[QuadratureRouteIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, live_routes, rejected_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, live_routes, rejected_routes))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"QUADRATURE-ROUTE {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
