#!/usr/bin/env python3
"""Validate the relative-Gaussian quadrature ladder scout."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout import (  # noqa: E402
    DEFAULT_GRID_JSON,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_RESIDUAL_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgqls_01_worst_row_import",
    "nlrgqls_02_high_order_ladder_stability",
    "nlrgqls_03_quadrature_radius_target",
    "nlrgqls_04_floating_overflow_boundary",
    "nlrgqls_05_acceptance_gate",
}

ALLOWED_ROLES = {
    "finite_diagnostic",
    "floating_diagnostic",
    "open_requirement",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Quadrature Ladder Scout",
    "Status: high-order floating quadrature ladder scout",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: 5 rows, 0 issues, 7 ladder rows, 320 reference order, 0 ready-to-apply rows",
    "ladder orders: [96, 128, 160, 192, 224, 256, 320]",
    "reference order: 320",
    "max value ratio: 0.9707100590203351233111029",
    "max derivative ratio: 0.9693567774758094418100653",
    "max value ratio spread: 7.767821744285352977404633e-15",
    "max derivative ratio spread: 7.756680887393827019918227e-15",
    "proposed quadrature ratio radius cap: 1.0e-6",
    "cap keeps worst ladder below one: True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class QuadratureLadderIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> QuadratureLadderIssue:
    return QuadratureLadderIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("E+00", "E+0"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[QuadratureLadderIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[QuadratureLadderIssue]:
    issues: list[QuadratureLadderIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "high-order floating quadrature ladder scout":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_cancellation_reduced_grid_scout",
        "source_cancellation_reduced_grid_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_residual_budget",
        "source_node_c0_certificate",
        "source_phi_tail_bound_scout",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "floating quadrature ladder scout only",
        "does not provide interval nodes or weights",
        "does not prove a quadrature-remainder theorem",
        "uniform collar",
        "scaled-curvature",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(
    artifact: dict,
    grid_path: Path,
    interval_path: Path,
    residual_path: Path,
) -> list[QuadratureLadderIssue]:
    try:
        recomputed = build_artifact(grid_path, interval_path, residual_path)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[QuadratureLadderIssue] = []
    for key in ("location_diagnostics", "matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_location_diagnostics(artifact: dict) -> list[QuadratureLadderIssue]:
    locations = artifact.get("location_diagnostics", [])
    issues: list[QuadratureLadderIssue] = []
    if not isinstance(locations, list):
        return [issue("location_diagnostics", "bad-locations", repr(type(locations)))]
    if len(locations) != 1:
        issues.append(issue("location_diagnostics", "bad-location-count", str(len(locations))))
    if not locations:
        return issues
    location = locations[0]
    expected = {
        "T": 10000,
        "index": 21,
        "ladder_orders": [96, 128, 160, 192, 224, 256, 320],
        "reference_order": 320,
        "ladder_row_count": 7,
        "max_value_ratio": "0.9707100590203351233111029",
        "max_derivative_ratio": "0.9693567774758094418100653",
        "value_ratio_spread": "7.767821744285352977404633e-15",
        "derivative_ratio_spread": "7.756680887393827019918227e-15",
        "all_ladder_ratios_below_one": True,
        "cap_keeps_worst_ladder_below_one": True,
        "proposed_quadrature_ratio_radius_cap": "1.0e-6",
    }
    for key, value in expected.items():
        if location.get(key) != value:
            issues.append(issue("location_diagnostics", f"bad-{key}", f"{location.get(key)!r} != {value!r}"))
    if dec(location.get("value_ratio_spread")) >= Decimal("1e-12"):
        issues.append(issue("location_diagnostics", "value-spread-too-large", repr(location.get("value_ratio_spread"))))
    if dec(location.get("derivative_ratio_spread")) >= Decimal("1e-12"):
        issues.append(
            issue("location_diagnostics", "derivative-spread-too-large", repr(location.get("derivative_ratio_spread")))
        )
    rows = location.get("ladder_rows", [])
    if len(rows) != 7:
        issues.append(issue("ladder_rows", "bad-row-count", str(len(rows))))
    if [row.get("quadrature_order") for row in rows] != [96, 128, 160, 192, 224, 256, 320]:
        issues.append(issue("ladder_rows", "bad-order-sequence", repr([row.get("quadrature_order") for row in rows])))
    for row in rows:
        row_id = f"N={row.get('quadrature_order')}"
        for key in (
            "quadrature_order",
            "value_ratio_to_first_omitted",
            "derivative_ratio_to_first_omitted",
            "value_delta_from_reference_ratio",
            "derivative_delta_from_reference_ratio",
            "proof_boundary",
        ):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if not Decimal("0.96") < dec(row.get("value_ratio_to_first_omitted")) < Decimal("0.98"):
            issues.append(issue(row_id, "bad-value-ratio", repr(row.get("value_ratio_to_first_omitted"))))
        if not Decimal("0.96") < dec(row.get("derivative_ratio_to_first_omitted")) < Decimal("0.98"):
            issues.append(issue(row_id, "bad-derivative-ratio", repr(row.get("derivative_ratio_to_first_omitted"))))
        if "floating high-order quadrature" not in str(row.get("proof_boundary", "")).lower():
            issues.append(issue(row_id, "weak-row-boundary", repr(row.get("proof_boundary"))))
    return issues


def validate_summary(artifact: dict) -> list[QuadratureLadderIssue]:
    summary = artifact.get("summary", {})
    issues: list[QuadratureLadderIssue] = []
    expected = {
        "matrix_rows": 5,
        "location_count": 1,
        "ladder_orders": [96, 128, 160, 192, 224, 256, 320],
        "reference_order": 320,
        "total_ladder_rows": 7,
        "max_value_ratio": "0.9707100590203351233111029",
        "max_derivative_ratio": "0.9693567774758094418100653",
        "max_value_ratio_spread": "7.767821744285352977404633e-15",
        "max_derivative_ratio_spread": "7.756680887393827019918227e-15",
        "all_ladder_ratios_below_one": True,
        "proposed_quadrature_ratio_radius_cap": "1.0e-6",
        "intervalization_per_source_cap": "2.000000000000000000E-3",
        "cap_below_per_source_cap": True,
        "cap_keeps_worst_ladder_below_one": True,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if dec(summary.get("max_value_ratio")) + dec(summary.get("proposed_quadrature_ratio_radius_cap")) >= Decimal(1):
        issues.append(issue("summary", "value-cap-does-not-close", repr(summary.get("max_value_ratio"))))
    if dec(summary.get("max_derivative_ratio")) + dec(summary.get("proposed_quadrature_ratio_radius_cap")) >= Decimal(1):
        issues.append(issue("summary", "derivative-cap-does-not-close", repr(summary.get("max_derivative_ratio"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("worst", "n=96..320", "below one", "order-spread below 1e-14", "1e-6", "floating evidence"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "floating evidence", "target, not a proved remainder", "nodes and weights", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_rows(artifact: dict) -> list[QuadratureLadderIssue]:
    rows = artifact.get("matrix_rows", [])
    issues: list[QuadratureLadderIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))]
    if len(rows) != 5:
        issues.append(issue("matrix_rows", "bad-row-count", str(len(rows))))
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary", "source_artifacts"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    target = rows_by_id.get("nlrgqls_03_quadrature_radius_target", {}).get("diagnostics", {})
    if target.get("cap_below_per_source_cap") is not True or target.get("cap_keeps_worst_ladder_below_one") is not True:
        issues.append(issue("nlrgqls_03_quadrature_radius_target", "cap-not-sufficient", repr(target)))
    return issues


def validate_note(path: Path) -> list[QuadratureLadderIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[QuadratureLadderIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "quadrature-remainder theorem is proved",
        "finite-grid interval certificate is complete",
        "uniform residual estimate is proved",
        "laguerre weight intervals are certified",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path, grid_path: Path, interval_path: Path, residual_path: Path) -> tuple[list[QuadratureLadderIssue], dict]:
    artifact = load_json(target_path)
    issues: list[QuadratureLadderIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, grid_path, interval_path, residual_path))
    issues.extend(validate_location_diagnostics(artifact))
    issues.extend(validate_summary(artifact))
    issues.extend(validate_rows(artifact))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--residual-json", type=Path, default=DEFAULT_RESIDUAL_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    residual_path = args.residual_json if args.residual_json.is_absolute() else REPO_ROOT / args.residual_json
    issues, summary = validate(target, note, grid_path, interval_path, residual_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-QUAD-LADDER {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian quadrature ladder scout: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('total_ladder_rows')} ladder rows, "
            f"{summary.get('reference_order')} reference order, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
