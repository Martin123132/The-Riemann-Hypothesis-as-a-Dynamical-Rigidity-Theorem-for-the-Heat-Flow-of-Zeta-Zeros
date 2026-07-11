#!/usr/bin/env python3
"""Validate the worst-row Christoffel-weight interval certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate import (  # noqa: E402
    DEFAULT_INTERVAL_JSON,
    DEFAULT_MIDPOINT_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_ROOT_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwcwic_01_taylor_denominator_formula",
    "nlrgwcwic_02_taylor_denominator_intervals",
    "nlrgwcwic_03_christoffel_weight_intervals",
    "nlrgwcwic_04_underflow_repair_promoted_to_interval",
    "nlrgwcwic_05_mass_sum_cross_check",
    "nlrgwcwic_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_formula",
    "arb_interval_certificate",
    "finite_consistency_check",
    "acceptance_gate",
}

ALLOWED_READINESS = {"available_exact", "available_interval_certificate", "not_ready_to_apply"}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Christoffel-Weight Interval Certificate",
    "Status: worst-row Christoffel-weight interval certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Christoffel-weight interval certificate: 6 rows, 0 issues, 320 interval weights, 0 Taylor denominator obstructions, 30 repaired floating underflows, 0 ready-to-apply rows",
    "Taylor denominator obstructions: 0",
    "direct recurrence denominator obstructions: 320",
    "repaired floating underflows: 30",
    "mass Gamma(43/2) contained: True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_laguerre_root_bracket_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_midpoint_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class WeightIntervalIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> WeightIntervalIssue:
    return WeightIntervalIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[WeightIntervalIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[WeightIntervalIssue]:
    issues: list[WeightIntervalIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate"
    ):
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row Christoffel-weight interval certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_worst_row_laguerre_root_bracket_certificate",
        "source_worst_row_laguerre_root_bracket_json",
        "source_worst_row_christoffel_weight_midpoint_scout",
        "source_worst_row_christoffel_weight_midpoint_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_quadrature_ladder_scout",
        "source_node_c0_certificate",
        "source_phi_tail_grid_certificate",
        "source_coefficient_core_certificate",
        "source_first_omitted_denominator_certificate",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "weight interval certificate only",
        "does not evaluate phi",
        "does not prove a quadrature-remainder",
        "all recorded rows",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[WeightIntervalIssue]:
    try:
        recomputed = build_artifact(DEFAULT_ROOT_JSON, DEFAULT_MIDPOINT_JSON, DEFAULT_INTERVAL_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[WeightIntervalIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[WeightIntervalIssue], int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[WeightIntervalIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
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
            issues.append(issue(row_id, "forbidden-ready-to-apply", row_id))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "does not")):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    return issues, len(rows)


def validate_weight_rows(artifact: dict) -> list[WeightIntervalIssue]:
    rows_by_id = {row.get("id"): row for row in artifact.get("matrix_rows", []) if isinstance(row, dict)}
    diagnostics = rows_by_id.get("nlrgwcwic_03_christoffel_weight_intervals", {}).get("diagnostics", {})
    weight_rows = diagnostics.get("weight_interval_rows", [])
    summary = diagnostics.get("weight_summary", {})
    issues: list[WeightIntervalIssue] = []
    if len(weight_rows) != 320:
        issues.append(issue("weight_rows", "bad-row-count", str(len(weight_rows))))
        return issues
    underflows = 0
    min_weight: tuple[Decimal, int] | None = None
    max_weight: tuple[Decimal, int] | None = None
    max_relative_width: tuple[Decimal, int] | None = None
    for expected_index, row in enumerate(weight_rows, start=1):
        row_id = f"weight_{expected_index}"
        if row.get("root_index") != expected_index:
            issues.append(issue(row_id, "bad-root-index", repr(row.get("root_index"))))
        denominator_lower = dec(row.get("denominator_lower", "0"))
        denominator_upper = dec(row.get("denominator_upper", "0"))
        if not (denominator_upper < 0 or denominator_lower > 0):
            issues.append(issue(row_id, "denominator-contains-zero", f"{denominator_lower}, {denominator_upper}"))
        weight_lower = dec(row.get("weight_lower", "0"))
        weight_upper = dec(row.get("weight_upper", "0"))
        if weight_lower <= 0 or weight_upper <= 0 or weight_lower > weight_upper:
            issues.append(issue(row_id, "bad-positive-weight-interval", f"{weight_lower}, {weight_upper}"))
        if bool(row.get("scipy_underflowed_to_zero")):
            underflows += 1
            if row.get("scipy_float_weight") != "0.0":
                issues.append(issue(row_id, "bad-underflow-flag", repr(row.get("scipy_float_weight"))))
        relative_width = dec(row.get("relative_weight_width", "0"))
        if relative_width <= 0:
            issues.append(issue(row_id, "bad-relative-width", repr(row.get("relative_weight_width"))))
        if min_weight is None or weight_lower < min_weight[0]:
            min_weight = (weight_lower, expected_index)
        if max_weight is None or weight_upper > max_weight[0]:
            max_weight = (weight_upper, expected_index)
        if max_relative_width is None or relative_width > max_relative_width[0]:
            max_relative_width = (relative_width, expected_index)
        boundary = str(row.get("proof_boundary", "")).lower()
        for required in ("worst-row", "not phi", "lambda <= 0"):
            if required not in boundary:
                issues.append(issue(row_id, "weak-proof-boundary", required))
    if underflows != 30:
        issues.append(issue("weight_rows", "bad-underflow-count", str(underflows)))
    if min_weight is None or min_weight[1] != 320:
        issues.append(issue("weight_rows", "bad-min-weight-root", repr(min_weight)))
    if max_weight is None or max_weight[1] != 43:
        issues.append(issue("weight_rows", "bad-max-weight-root", repr(max_weight)))
    if max_relative_width is None or max_relative_width[1] != 1:
        issues.append(issue("weight_rows", "bad-max-relative-width-root", repr(max_relative_width)))
    expected_summary = {
        "weight_interval_rows": 320,
        "quadrature_order": 320,
        "index": 21,
        "T": 10000,
        "alpha": "41/2",
        "all_taylor_denominators_separated": True,
        "taylor_denominator_contains_zero_rows": 0,
        "direct_interval_denominator_contains_zero_rows": 320,
        "all_weight_intervals_positive": True,
        "zero_scipy_float_weights_repaired": 30,
        "minimum_weight_root_index": 320,
        "maximum_weight_root_index": 43,
        "maximum_relative_weight_width_root_index": 1,
        "mass_contained_by_weight_sum_interval": True,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("weight_summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))
    for key, snippet in (
        ("minimum_weight_interval", "2.5663925426991066e-492"),
        ("maximum_weight_interval", "7.762348957842933e+17"),
        ("maximum_relative_weight_width", "6.645705603047399351898944020607455890895e-16"),
        ("relative_weight_sum_mass_error_interval", "+/- 1.26e-17"),
    ):
        if snippet not in str(summary.get(key, "")):
            issues.append(issue("weight_summary", f"bad-{key}", f"missing {snippet!r} in {summary.get(key)!r}"))
    return issues


def validate_summary(artifact: dict, row_count: int) -> list[WeightIntervalIssue]:
    summary = artifact.get("summary", {})
    issues: list[WeightIntervalIssue] = []
    expected = {
        "matrix_rows": 6,
        "weight_interval_rows": 320,
        "quadrature_order": 320,
        "index": 21,
        "T": 10000,
        "precision_bits": 1024,
        "taylor_denominator_contains_zero_rows": 0,
        "direct_interval_denominator_contains_zero_rows": 320,
        "zero_scipy_float_weights_repaired": 30,
        "minimum_weight_root_index": 320,
        "maximum_weight_root_index": 43,
        "maximum_relative_weight_width_root_index": 1,
        "mass_contained_by_weight_sum_interval": True,
        "available_interval_certificate_rows": 3,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected in expected.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "taylor-centered arb enclosure",
        "all 320 denominator intervals avoid zero",
        "all 320 christoffel-weight intervals are positive",
        "remain open",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    return issues


def validate_note(path: Path) -> list[WeightIntervalIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[WeightIntervalIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    forbidden_phrases = (
        "phi node intervals are certified",
        "quadrature-remainder theorem is proved",
        "finite-grid interval certificate is complete",
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
    issues: list[WeightIntervalIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_weight_rows(artifact))
    issues.extend(validate_summary(artifact, row_count))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"WEIGHT-INTERVAL {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
