#!/usr/bin/env python3
"""Validate the cancellation-reduced actual-remainder grid scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout import (  # noqa: E402
    DEFAULT_INDICES,
    DEFAULT_MPMATH_DPS,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_PHI_TERM_COUNT,
    DEFAULT_POLYNOMIAL_M,
    DEFAULT_PRECISION_BITS,
    DEFAULT_QUADRATURE_ORDERS,
    DEFAULT_RATIO_CUTOFF_N,
    DEFAULT_RESIDUAL_JSON,
    DEFAULT_SELECTED_ORDER,
    DEFAULT_T_GRID,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgcrrgs_01_cancellation_reduced_identity",
    "nlrgcrrgs_02_floating_grid_quadrature",
    "nlrgcrrgs_03_first_omitted_grid_comparison",
    "nlrgcrrgs_04_far_T_shape_evidence",
    "nlrgcrrgs_05_uniform_promotion_rejected",
    "nlrgcrrgs_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "floating_diagnostic",
    "finite_diagnostic",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Cancellation-Reduced Remainder Grid Scout",
    "Status: floating cancellation-reduced theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: 6 rows, 0 issues, 20 grid rows, 5 T values, 0 ready-to-apply rows",
    "all grid ratios below one: True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class CancellationReducedGridIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> CancellationReducedGridIssue:
    return CancellationReducedGridIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[CancellationReducedGridIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[CancellationReducedGridIssue]:
    issues: list[CancellationReducedGridIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "floating cancellation-reduced theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree40_residual_tail_budget",
        "source_formal_tail_obstruction",
        "source_asymptotic_remainder_target",
        "source_actual_endpoint_scout",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "floating cancellation-reduced",
        "floating generalized laguerre",
        "does not provide interval",
        "does not prove",
        "uniform collar",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict, residual_path: Path) -> list[CancellationReducedGridIssue]:
    try:
        recomputed = build_artifact(
            residual_path,
            DEFAULT_T_GRID,
            DEFAULT_QUADRATURE_ORDERS,
            DEFAULT_SELECTED_ORDER,
            DEFAULT_INDICES,
            DEFAULT_POLYNOMIAL_M,
            DEFAULT_RATIO_CUTOFF_N,
            DEFAULT_PRECISION_BITS,
            DEFAULT_MPMATH_DPS,
            DEFAULT_PHI_TERM_COUNT,
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[CancellationReducedGridIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[CancellationReducedGridIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[1].get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    expected_params = {
        "t_grid": [1156, 1500, 2000, 5000, 10000],
        "indices": [21, 22, 23, 24],
        "quadrature_rule": "generalized Gauss-Laguerre for Gamma(i+1/2, 1)",
        "quadrature_orders": [64, 96, 128, 192],
        "selected_quadrature_order": 192,
        "polynomial_degree": 40,
        "polynomial_M": 20,
        "residual_first_j": 21,
        "ratio_cutoff_n": 80,
        "precision_bits_for_ratio_build": 384,
        "mpmath_dps": 80,
        "phi_n_terms": 30,
        "value_budget_A": "5.382819486765314521E-01",
        "derivative_budget_B": "9.315354075509573936E-03",
        "value_core": "Phi(sqrt(v))/Phi(0)-sum_{j=0}^{20} r_j*v^j",
        "derivative_core": "sqrt(v)*Phi'(sqrt(v))/(2*Phi(0))-sum_{j=1}^{20} j*r_j*v^j",
        "scaled_value_residual": "|E[value_core(uY)]|/u^3",
        "scaled_derivative_residual": "|E[derivative_core(uY)]|/u^2",
    }
    issues: list[CancellationReducedGridIssue] = []
    for key, value in expected_params.items():
        if params.get(key) != value:
            issues.append(issue("parameters", f"bad-{key}", f"{params.get(key)!r} != {value!r}"))
    if "normalizing_phi0_mpmath" not in params:
        issues.append(issue("parameters", "missing-normalizing-phi0", repr(params)))
    grid_rows = diagnostics.get("grid_rows", [])
    if diagnostics.get("grid_row_count") != 20:
        issues.append(issue("diagnostics", "bad-grid-row-count", repr(diagnostics.get("grid_row_count"))))
    if diagnostics.get("t_grid_count") != 5:
        issues.append(issue("diagnostics", "bad-t-grid-count", repr(diagnostics.get("t_grid_count"))))
    if diagnostics.get("index_count") != 4:
        issues.append(issue("diagnostics", "bad-index-count", repr(diagnostics.get("index_count"))))
    if diagnostics.get("quadrature_order_count") != 4:
        issues.append(issue("diagnostics", "bad-quadrature-order-count", repr(diagnostics.get("quadrature_order_count"))))
    if diagnostics.get("all_grid_ratios_below_one") is not True:
        issues.append(issue("diagnostics", "grid-ratios-not-below-one", "expected true"))
    thresholds = {
        "max_value_ratio_to_first_omitted_over_orders": 0.98,
        "max_derivative_ratio_to_first_omitted_over_orders": 0.98,
        "max_value_ratio_spread_over_orders": 1.0e-12,
        "max_derivative_ratio_spread_over_orders": 1.0e-12,
        "max_selected_value_budget_fraction": 1.0e-3,
        "max_selected_derivative_budget_fraction": 1.0e-3,
    }
    for key, upper in thresholds.items():
        try:
            value = float(diagnostics.get(key))
        except (TypeError, ValueError):
            issues.append(issue("diagnostics", f"bad-{key}", repr(diagnostics.get(key))))
            continue
        if not value < upper:
            issues.append(issue("diagnostics", f"threshold-{key}", f"{value} >= {upper}"))
    if diagnostics.get("max_value_ratio_location") != {"T": 10000, "index": 21}:
        issues.append(issue("diagnostics", "bad-max-value-location", repr(diagnostics.get("max_value_ratio_location"))))
    if diagnostics.get("max_derivative_ratio_location") != {"T": 10000, "index": 21}:
        issues.append(
            issue("diagnostics", "bad-max-derivative-location", repr(diagnostics.get("max_derivative_ratio_location")))
        )
    expected_pairs = [(T, i) for T in [1156, 1500, 2000, 5000, 10000] for i in [21, 22, 23, 24]]
    actual_pairs = [(row.get("T"), row.get("index")) for row in grid_rows]
    if actual_pairs != expected_pairs:
        issues.append(issue("grid_rows", "bad-grid-order", repr(actual_pairs)))
    for row in grid_rows:
        label = f"T={row.get('T')},F_{row.get('index')}"
        for key in (
            "selected_value_ratio_to_first_omitted",
            "selected_derivative_ratio_to_first_omitted",
            "max_value_ratio_to_first_omitted_over_orders",
            "max_derivative_ratio_to_first_omitted_over_orders",
        ):
            try:
                value = float(row.get(key))
            except (TypeError, ValueError):
                issues.append(issue(label, f"bad-{key}", repr(row.get(key))))
                continue
            if not 0.0 < value < 1.0:
                issues.append(issue(label, f"ratio-outside-unit-interval-{key}", str(value)))
        if "floating quadrature" not in str(row.get("proof_boundary", "")).lower():
            issues.append(issue(label, "weak-row-boundary", repr(row.get("proof_boundary"))))
    note = str(diagnostics.get("cancellation_reduction_note", "")).lower()
    for required in (
        "subtracted inside",
        "gamma expectation",
        "far-t double-precision cancellation",
        "floating-node diagnostic",
        "not an interval enclosure",
    ):
        if required not in note:
            issues.append(issue("diagnostics", "weak-cancellation-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[CancellationReducedGridIssue], int, int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[CancellationReducedGridIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_rows = 0
    floating_rows = 0
    finite_rows = 0
    rejected_rows = 0
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
        if role == "exact_reduction":
            exact_rows += 1
        if role == "floating_diagnostic":
            floating_rows += 1
        if role == "finite_diagnostic":
            finite_rows += 1
        if role == "rejected_route":
            rejected_rows += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") not in {"available_exact", "not_ready_to_apply"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), exact_rows, floating_rows, finite_rows, rejected_rows


def validate_summary(
    artifact: dict,
    row_count: int,
    exact_rows: int,
    floating_rows: int,
    finite_rows: int,
    rejected_rows: int,
) -> list[CancellationReducedGridIssue]:
    summary = artifact.get("summary", {})
    issues: list[CancellationReducedGridIssue] = []
    expected_exact = {
        "matrix_rows": 6,
        "grid_rows": 20,
        "t_grid_count": 5,
        "index_count": 4,
        "quadrature_order_count": 4,
        "selected_quadrature_order": 192,
        "max_value_ratio_location": {"T": 10000, "index": 21},
        "max_derivative_ratio_location": {"T": 10000, "index": 21},
        "all_grid_ratios_below_one": True,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected_exact.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 6:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_rows != 1:
        issues.append(issue("summary", "bad-exact-row-count", str(exact_rows)))
    if floating_rows != 1:
        issues.append(issue("summary", "bad-floating-row-count", str(floating_rows)))
    if finite_rows != 2:
        issues.append(issue("summary", "bad-finite-row-count", str(finite_rows)))
    if rejected_rows != 1:
        issues.append(issue("summary", "bad-rejected-row-count", str(rejected_rows)))
    thresholds = {
        "max_value_ratio_to_first_omitted_over_orders": 0.98,
        "max_derivative_ratio_to_first_omitted_over_orders": 0.98,
        "max_selected_value_budget_fraction": 1.0e-3,
        "max_selected_derivative_budget_fraction": 1.0e-3,
    }
    for key, upper in thresholds.items():
        try:
            value = float(summary.get(key))
        except (TypeError, ValueError):
            issues.append(issue("summary", f"bad-{key}", repr(summary.get(key))))
            continue
        if not value < upper:
            issues.append(issue("summary", f"threshold-{key}", f"{value} >= {upper}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "subtracting the degree-40 polynomial inside",
        "far-t catastrophic cancellation",
        "t=1156,1500,2000,5000,10000",
        "below one first omitted",
        "0.9707100590",
        "0.9693567775",
        "finite floating evidence",
        "interval-certified",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "not promoted", "not interval", "finite grid", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[CancellationReducedGridIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[CancellationReducedGridIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "uniform residual estimate is proved",
        "residual-tail estimates are proved",
        "first-omitted-term theorem is proved",
        "actual remainder theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path, residual_path: Path) -> tuple[list[CancellationReducedGridIssue], dict]:
    artifact = load_json(target_path)
    issues: list[CancellationReducedGridIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, residual_path))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, exact_rows, floating_rows, finite_rows, rejected_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_rows, floating_rows, finite_rows, rejected_rows))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--residual-json", type=Path, default=DEFAULT_RESIDUAL_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    residual_path = args.residual_json if args.residual_json.is_absolute() else REPO_ROOT / args.residual_json
    issues, summary = validate(target, note, residual_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-CANCELLATION-GRID {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian cancellation-reduced remainder grid scout: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('grid_rows')} grid rows, "
            f"{summary.get('t_grid_count')} T values, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
