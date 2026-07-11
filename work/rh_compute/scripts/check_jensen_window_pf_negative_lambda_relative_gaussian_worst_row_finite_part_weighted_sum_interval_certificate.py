#!/usr/bin/env python3
"""Validate the worst-row finite-part weighted-sum interval certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate import (  # noqa: E402
    DEFAULT_COEFFICIENT_JSON,
    DEFAULT_FIRST_OMITTED_JSON,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_PHI_TAIL_JSON,
    DEFAULT_WEIGHT_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwrfpwsic_01_exact_moment_rewrite",
    "nlrgwrfpwsic_02_refined_nodes_and_weights",
    "nlrgwrfpwsic_03_finite_phi_weighted_sum",
    "nlrgwrfpwsic_04_first_omitted_comparison",
    "nlrgwrfpwsic_05_intervalization_handoff",
    "nlrgwrfpwsic_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "arb_interval_certificate",
    "intervalization_handoff",
    "acceptance_gate",
}

ALLOWED_READINESS = {"available_exact", "available_interval_certificate", "not_ready_to_apply"}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Finite-Part Weighted-Sum Interval Certificate",
    "Status: worst-row finite-part weighted-sum interval certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-part weighted-sum interval certificate: 6 rows, 0 issues, 320 refined nodes, 320 interval weights, 2 below-one ratios, 0 ready-to-apply rows",
    "value ratio upper: 0.9833957992836557769419015895036210773888",
    "derivative ratio upper: 0.9694055674762067320093698741711260875260",
    "both ratios certified below one: True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_christoffel_weight_interval_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_coefficient_core_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class FinitePartIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> FinitePartIssue:
    return FinitePartIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[FinitePartIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[FinitePartIssue]:
    issues: list[FinitePartIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate"
    ):
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row finite-part weighted-sum interval certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_weight_interval_certificate",
        "source_weight_interval_json",
        "source_root_bracket_certificate",
        "source_phi_tail_grid_certificate",
        "source_phi_tail_grid_json",
        "source_coefficient_core_certificate",
        "source_coefficient_core_json",
        "source_first_omitted_denominator_certificate",
        "source_first_omitted_denominator_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_quadrature_ladder_scout",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite-part weighted-sum interval certificate only",
        "does not compose the n>30 tail",
        "does not prove a quadrature-remainder",
        "all recorded rows",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[FinitePartIssue]:
    try:
        recomputed = build_artifact(
            DEFAULT_WEIGHT_JSON,
            DEFAULT_PHI_TAIL_JSON,
            DEFAULT_COEFFICIENT_JSON,
            DEFAULT_FIRST_OMITTED_JSON,
            DEFAULT_INTERVAL_JSON,
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[FinitePartIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[FinitePartIssue], int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[FinitePartIssue] = []
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


def validate_node_rows(artifact: dict) -> list[FinitePartIssue]:
    rows_by_id = {row.get("id"): row for row in artifact.get("matrix_rows", []) if isinstance(row, dict)}
    diagnostics = rows_by_id.get("nlrgwrfpwsic_03_finite_phi_weighted_sum", {}).get("diagnostics", {})
    node_rows = diagnostics.get("node_evaluation_rows", [])
    issues: list[FinitePartIssue] = []
    if len(node_rows) != 320:
        issues.append(issue("node_rows", "bad-row-count", str(len(node_rows))))
        return issues
    max_relative_weight_width: tuple[Decimal, int] | None = None
    max_node_width: tuple[Decimal, int] | None = None
    for expected_index, row in enumerate(node_rows, start=1):
        row_id = f"node_{expected_index}"
        if row.get("root_index") != expected_index:
            issues.append(issue(row_id, "bad-root-index", repr(row.get("root_index"))))
        node_width = dec(row.get("node_width", "0"))
        if node_width <= 0:
            issues.append(issue(row_id, "bad-node-width", repr(row.get("node_width"))))
        relative_weight_width = dec(row.get("relative_weight_width", "0"))
        if relative_weight_width <= 0:
            issues.append(issue(row_id, "bad-relative-weight-width", repr(row.get("relative_weight_width"))))
        if max_relative_weight_width is None or relative_weight_width > max_relative_weight_width[0]:
            max_relative_weight_width = (relative_weight_width, expected_index)
        if max_node_width is None or node_width > max_node_width[0]:
            max_node_width = (node_width, expected_index)
        boundary = str(row.get("proof_boundary", "")).lower()
        for required in ("refined", "not quadrature remainder", "lambda <= 0"):
            if required not in boundary:
                issues.append(issue(row_id, "weak-proof-boundary", required))
    if max_relative_weight_width is None or max_relative_weight_width[1] != 1:
        issues.append(issue("node_rows", "bad-max-relative-weight-width-root", repr(max_relative_weight_width)))
    if max_node_width is None or max_node_width[1] != 320:
        issues.append(issue("node_rows", "bad-widest-node-root", repr(max_node_width)))
    return issues


def validate_summary(artifact: dict, row_count: int) -> list[FinitePartIssue]:
    summary = artifact.get("summary", {})
    issues: list[FinitePartIssue] = []
    expected = {
        "matrix_rows": 6,
        "refined_node_rows": 320,
        "interval_weight_rows": 320,
        "quadrature_order": 320,
        "index": 21,
        "T": 10000,
        "refinement_steps": 120,
        "precision_bits": 4096,
        "phi_term_count": 30,
        "value_residual_separated_from_zero": True,
        "derivative_residual_separated_from_zero": True,
        "below_one_ratio_count": 2,
        "both_ratios_certified_below_one": True,
        "available_interval_certificate_rows": 3,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    value_ratio = dec(summary.get("value_ratio_to_first_omitted_upper", "2"))
    derivative_ratio = dec(summary.get("derivative_ratio_to_first_omitted_upper", "2"))
    if not (Decimal("0") < value_ratio < Decimal("1")):
        issues.append(issue("summary", "value-ratio-not-below-one", str(value_ratio)))
    if not (Decimal("0") < derivative_ratio < Decimal("1")):
        issues.append(issue("summary", "derivative-ratio-not-below-one", str(derivative_ratio)))
    for key, snippet in (
        ("widest_refined_node_width", "3.595897829114163907E-35"),
        ("maximum_relative_weight_width", "5.764230762018462413243017966627116467124e-34"),
        ("value_ratio_to_first_omitted_upper", "0.9833957992836557769419015895036210773888"),
        ("derivative_ratio_to_first_omitted_upper", "0.9694055674762067320093698741711260875260"),
    ):
        if snippet not in str(summary.get(key, "")):
            issues.append(issue("summary", f"bad-{key}", f"missing {snippet!r} in {summary.get(key)!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "finite n<=30 phi",
        "exact gamma moments",
        "certified below one first omitted",
        "remain open",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    return issues


def validate_note(path: Path) -> list[FinitePartIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[FinitePartIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    forbidden_phrases = (
        "quadrature-remainder theorem is proved",
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
    issues: list[FinitePartIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_node_rows(artifact))
    issues.extend(validate_summary(artifact, row_count))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"FINITE-PART {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
