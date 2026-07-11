#!/usr/bin/env python3
"""Validate the relative-Gaussian coefficient-core propagation certificate."""

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

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate import (  # noqa: E402
    DEFAULT_INDICES,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_POLYNOMIAL_M,
    DEFAULT_RATIO_CUTOFF_N,
    DEFAULT_T_GRID,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgccc_01_ratio_ball_source",
    "nlrgccc_02_value_core_gamma_propagation",
    "nlrgccc_03_derivative_core_gamma_propagation",
    "nlrgccc_04_intervalization_handoff",
    "nlrgccc_05_coefficient_only_promotion_rejected",
    "nlrgccc_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "arb_coefficient_certificate",
    "coefficient_radius_propagation",
    "exact_budget_handoff",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_arb",
    "available_for_intervalization",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Coefficient-Core Propagation Certificate",
    "Status: coefficient-core propagation certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian coefficient-core propagation certificate: 6 rows, 0 issues, 22 coefficient rows, 20 propagation rows, 2 intervalization rows, 0 ready-to-apply rows",
    "sum_{j=0}^{20} rad(r_j)*(i+1/2)_j*u^(j-3)",
    "sum_{j=1}^{20} j*rad(r_j)*(i+1/2)_j*u^(j-2)",
    "maximum value coefficient ratio: 8.610446518945492184E-81 at T=10000, F_21",
    "maximum derivative coefficient ratio: 5.523600431279558041E-83 at T=10000, F_21",
    "ratio cap: 1.000000000000000000E-6",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class CoefficientCoreIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> CoefficientCoreIssue:
    return CoefficientCoreIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[CoefficientCoreIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[CoefficientCoreIssue]:
    issues: list[CoefficientCoreIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "coefficient-core propagation certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree16_certificate",
        "source_cancellation_reduced_grid_scout",
        "source_intervalization_target",
        "source_first_omitted_denominator_certificate",
        "source_quadrature_ladder_scout",
        "source_node_c0_range_certificate",
        "source_phi_tail_grid_certificate",
        "source_uniform_remainder_target",
        "source_formal_tail_obstruction",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "coefficient-core propagation certificate only",
        "does not enclose",
        "laguerre node/weight intervals",
        "quadrature-remainder",
        "uniform collar",
        "scaled-curvature",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[CoefficientCoreIssue]:
    try:
        recomputed = build_artifact()
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[CoefficientCoreIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[CoefficientCoreIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[CoefficientCoreIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    available_for_intervalization = 0
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
        if row.get("readiness") == "available_for_intervalization":
            available_for_intervalization += 1
        if row.get("readiness") == "ready_to_apply":
            issues.append(issue(row_id, "forbidden-ready-to-apply", row_id))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    return issues, len(rows), available_for_intervalization


def validate_coefficients(artifact: dict) -> list[CoefficientCoreIssue]:
    rows_by_id = {row.get("id"): row for row in artifact.get("matrix_rows", []) if isinstance(row, dict)}
    diagnostics = rows_by_id.get("nlrgccc_01_ratio_ball_source", {}).get("diagnostics", {})
    issues: list[CoefficientCoreIssue] = []
    coefficient_rows = diagnostics.get("coefficient_rows", [])
    if diagnostics.get("coefficient_rows_count") != 22 or len(coefficient_rows) != 22:
        issues.append(issue("coefficient_rows", "bad-row-count", repr(diagnostics.get("coefficient_rows_count"))))
    for expected_index, row in enumerate(coefficient_rows):
        if not isinstance(row, dict):
            issues.append(issue("coefficient_rows", "bad-row", repr(row)))
            continue
        row_id = f"r_{expected_index}"
        if row.get("coefficient_index") != expected_index:
            issues.append(issue(row_id, "bad-index", repr(row.get("coefficient_index"))))
        if row.get("degree") != 2 * expected_index:
            issues.append(issue(row_id, "bad-degree", repr(row.get("degree"))))
        if row.get("sign") not in {"positive", "negative"}:
            issues.append(issue(row_id, "bad-sign", repr(row.get("sign"))))
        if dec(row.get("radius_upper", "0")) <= 0:
            issues.append(issue(row_id, "nonpositive-radius", repr(row.get("radius_upper"))))
    if len(coefficient_rows) == 22 and coefficient_rows[21].get("sign") != "negative":
        issues.append(issue("r_21", "bad-sign", repr(coefficient_rows[21].get("sign"))))
    if diagnostics.get("maximum_ratio_ball_radius_index") != 21:
        issues.append(issue("coefficient_rows", "bad-max-radius-index", repr(diagnostics.get("maximum_ratio_ball_radius_index"))))
    if dec(diagnostics.get("maximum_ratio_ball_radius_upper", "1")) >= Decimal("1E-70"):
        issues.append(issue("coefficient_rows", "max-radius-too-large", repr(diagnostics.get("maximum_ratio_ball_radius_upper"))))
    pointwise = diagnostics.get("pointwise_radius_envelopes", {})
    expected_pointwise_keys = (
        "value_core_radius_on_v_le_809_1156_upper",
        "derivative_core_radius_on_v_le_809_1156_upper",
        "value_core_radius_on_v_le_1_upper",
        "derivative_core_radius_on_v_le_1_upper",
    )
    for key in expected_pointwise_keys:
        if dec(pointwise.get(key, "1")) >= Decimal("1E-78"):
            issues.append(issue("pointwise_radius_envelopes", "radius-too-large", f"{key}={pointwise.get(key)!r}"))
    return issues


def validate_propagation(artifact: dict) -> list[CoefficientCoreIssue]:
    rows_by_id = {row.get("id"): row for row in artifact.get("matrix_rows", []) if isinstance(row, dict)}
    diagnostics = rows_by_id.get("nlrgccc_02_value_core_gamma_propagation", {}).get("diagnostics", {})
    issues: list[CoefficientCoreIssue] = []
    propagation_rows = diagnostics.get("propagation_rows", [])
    summary = diagnostics.get("propagation_summary", {})
    if len(propagation_rows) != len(DEFAULT_T_GRID) * len(DEFAULT_INDICES):
        issues.append(issue("propagation_rows", "bad-row-count", str(len(propagation_rows))))
    seen = {(row.get("T"), row.get("index")) for row in propagation_rows if isinstance(row, dict)}
    expected_seen = {(T, index) for T in DEFAULT_T_GRID for index in DEFAULT_INDICES}
    if seen != expected_seen:
        issues.append(issue("propagation_rows", "bad-grid", repr(sorted(seen))))
    for row in propagation_rows:
        if not isinstance(row, dict):
            issues.append(issue("propagation_rows", "bad-row", repr(row)))
            continue
        row_id = f"T={row.get('T')},F={row.get('index')}"
        if row.get("polynomial_M") != DEFAULT_POLYNOMIAL_M or row.get("first_j") != DEFAULT_POLYNOMIAL_M + 1:
            issues.append(issue(row_id, "bad-polynomial-parameters", repr(row)))
        for key in (
            "value_coefficient_scaled_radius_upper",
            "value_first_omitted_denominator_lower",
            "value_coefficient_ratio_to_first_omitted_upper",
            "derivative_coefficient_scaled_radius_upper",
            "derivative_first_omitted_denominator_lower",
            "derivative_coefficient_ratio_to_first_omitted_upper",
        ):
            if dec(row.get(key, "0")) <= 0:
                issues.append(issue(row_id, "nonpositive-field", f"{key}={row.get(key)!r}"))
        if dec(row.get("value_coefficient_ratio_to_first_omitted_upper", "1")) >= Decimal("1E-6"):
            issues.append(issue(row_id, "value-ratio-above-cap", repr(row.get("value_coefficient_ratio_to_first_omitted_upper"))))
        if dec(row.get("derivative_coefficient_ratio_to_first_omitted_upper", "1")) >= Decimal("1E-6"):
            issues.append(
                issue(row_id, "derivative-ratio-above-cap", repr(row.get("derivative_coefficient_ratio_to_first_omitted_upper")))
            )
        boundary = str(row.get("proof_boundary", "")).lower()
        for required in ("coefficient-radius", "does not", "quadrature"):
            if required not in boundary:
                issues.append(issue(row_id, "weak-row-boundary", required))
    expected_summary = {
        "propagation_rows": 20,
        "t_grid": list(DEFAULT_T_GRID),
        "indices": list(DEFAULT_INDICES),
        "polynomial_M": 20,
        "polynomial_degree": 40,
        "first_j": 21,
        "maximum_value_coefficient_ratio_to_first_omitted_upper": "8.610446518945492184E-81",
        "maximum_value_coefficient_ratio_location": "T=10000, F_21",
        "maximum_derivative_coefficient_ratio_to_first_omitted_upper": "5.523600431279558041E-83",
        "maximum_derivative_coefficient_ratio_location": "T=10000, F_21",
        "all_coefficient_ratios_below_ratio_cap": True,
        "all_coefficient_ratios_below_intervalization_source_cap": True,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("propagation_summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))
    return issues


def validate_summary(artifact: dict, row_count: int, available_for_intervalization: int) -> list[CoefficientCoreIssue]:
    summary = artifact.get("summary", {})
    issues: list[CoefficientCoreIssue] = []
    expected = {
        "matrix_rows": 6,
        "coefficient_rows": 22,
        "propagation_rows": 20,
        "available_for_intervalization_rows": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "all_coefficient_ratios_below_ratio_cap": True,
        "all_coefficient_ratios_below_intervalization_source_cap": True,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if available_for_intervalization != summary.get("available_for_intervalization_rows"):
        issues.append(
            issue(
                "summary",
                "available-count-mismatch",
                f"{available_for_intervalization} != {summary.get('available_for_intervalization_rows')!r}",
            )
        )
    if dec(summary.get("maximum_value_coefficient_ratio_to_first_omitted_upper", "1")) >= Decimal("1E-70"):
        issues.append(issue("summary", "value-max-too-large", repr(summary.get("maximum_value_coefficient_ratio_to_first_omitted_upper"))))
    if dec(summary.get("maximum_derivative_coefficient_ratio_to_first_omitted_upper", "1")) >= Decimal("1E-72"):
        issues.append(
            issue(
                "summary",
                "derivative-max-too-large",
                repr(summary.get("maximum_derivative_coefficient_ratio_to_first_omitted_upper")),
            )
        )
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("coefficient", "finite-grid", "not phi", "quadrature", "grid-to-collar"):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    return issues


def validate_note(path: Path) -> list[CoefficientCoreIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[CoefficientCoreIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    forbidden_phrases = (
        "first-omitted residual theorem is proved",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
        "rh is proved",
        "lambda <= 0 is proved",
    )
    lowered = text.lower()
    for phrase in forbidden_phrases:
        if phrase in lowered:
            issues.append(issue("note", "forbidden-promotion-language", phrase))
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
    issues: list[CoefficientCoreIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, available_for_intervalization = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_coefficients(artifact))
    issues.extend(validate_propagation(artifact))
    issues.extend(validate_summary(artifact, row_count, available_for_intervalization))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"COEFFICIENT-CORE {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
