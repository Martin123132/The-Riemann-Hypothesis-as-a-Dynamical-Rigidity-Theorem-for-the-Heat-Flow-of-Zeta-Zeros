#!/usr/bin/env python3
"""Validate the worst-row finite-plus-tail budget certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate import (  # noqa: E402
    DEFAULT_FINITE_PART_JSON,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_PHI_TAIL_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwrfptbc_01_budget_composition_rule",
    "nlrgwrfptbc_02_finite_part_import",
    "nlrgwrfptbc_03_tail_source_import",
    "nlrgwrfptbc_04_composed_below_one_budget",
    "nlrgwrfptbc_05_intervalization_handoff",
    "nlrgwrfptbc_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "budget_composition_rule",
    "interval_certificate_import",
    "tail_source_import",
    "budget_certificate",
    "intervalization_handoff",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_budget_rule",
    "available_interval_certificate",
    "available_tail_source",
    "available_budget_certificate",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Finite-Plus-Tail Budget Certificate",
    "Status: worst-row finite-plus-tail budget certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian worst-row finite-plus-tail budget certificate: 6 rows, 0 issues, 2 composed ratios, 3 tail sources, 0 ready-to-apply rows",
    "value finite ratio upper: 0.9833957992836557769419015895036210773888",
    "derivative finite ratio upper: 0.9694055674762067320093698741711260875260",
    "value composed ratio upper using full tail cap: 0.9853957992836557769419015895036210773888",
    "derivative composed ratio upper using full tail cap: 0.9714055674762067320093698741711260875260",
    "both composed ratios below one: True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class FinitePlusTailIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> FinitePlusTailIssue:
    return FinitePlusTailIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[FinitePlusTailIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[FinitePlusTailIssue]:
    issues: list[FinitePlusTailIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_plus_tail_budget_certificate"
    ):
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row finite-plus-tail budget certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_finite_part_weighted_sum_certificate",
        "source_finite_part_weighted_sum_json",
        "source_phi_tail_grid_certificate",
        "source_phi_tail_grid_json",
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
        "finite-plus-tail budget certificate only",
        "does not prove a quadrature-remainder",
        "all rows",
        "finite-grid interval certificate",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[FinitePlusTailIssue]:
    try:
        recomputed = build_artifact(DEFAULT_FINITE_PART_JSON, DEFAULT_PHI_TAIL_JSON, DEFAULT_INTERVAL_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[FinitePlusTailIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[FinitePlusTailIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[FinitePlusTailIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    budget_certificate_rows = 0
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
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "budget_certificate":
            budget_certificate_rows += 1
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "does not")):
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
    handoff = rows_by_id.get("nlrgwrfptbc_05_intervalization_handoff", {}).get("diagnostics", {})
    if handoff.get("target_closing") is not False:
        issues.append(issue("nlrgwrfptbc_05_intervalization_handoff", "target-closing", repr(handoff)))
    if len(handoff.get("remaining_sources", [])) != 5:
        issues.append(issue("nlrgwrfptbc_05_intervalization_handoff", "bad-remaining-sources", repr(handoff)))
    return issues, len(rows), budget_certificate_rows


def validate_summary(artifact: dict, row_count: int, budget_certificate_rows: int) -> list[FinitePlusTailIssue]:
    summary = artifact.get("summary", {})
    issues: list[FinitePlusTailIssue] = []
    expected = {
        "matrix_rows": 6,
        "finite_part_ratio_count": 2,
        "tail_source_rows": 3,
        "ready_to_apply_rows": 0,
        "source_finite_part_rows": 6,
        "source_tail_certificate_rows": 6,
        "source_interval_obligations": 8,
        "per_source_intervalization_cap": "2.000000000000000000E-3",
        "tail_budget_ratio_reserved": "0.002000000000000000000",
        "tail_source_certified": True,
        "tail_actual_max_ratio_to_cap": "3.778846084314790602E-1292",
        "value_finite_ratio_upper": "0.9833957992836557769419015895036210773888",
        "derivative_finite_ratio_upper": "0.9694055674762067320093698741711260875260",
        "both_finite_ratios_below_one": True,
        "value_composed_ratio_upper_using_full_tail_cap": "0.9853957992836557769419015895036210773888",
        "derivative_composed_ratio_upper_using_full_tail_cap": "0.9714055674762067320093698741711260875260",
        "value_remaining_margin_after_full_tail_cap": "0.0146042007163442230580984104963789226112",
        "derivative_remaining_margin_after_full_tail_cap": "0.0285944325237932679906301258288739124740",
        "composed_ratio_count": 2,
        "both_composed_ratios_below_one": True,
        "target_closing": False,
        "budget_is_conservative": True,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if budget_certificate_rows != 1:
        issues.append(issue("summary", "bad-budget-certificate-row-count", str(budget_certificate_rows)))
    if not (dec(summary.get("value_composed_ratio_upper_using_full_tail_cap", "2")) < Decimal(1)):
        issues.append(issue("summary", "value-composed-ratio-not-below-one", repr(summary.get("value_composed_ratio_upper_using_full_tail_cap"))))
    if not (dec(summary.get("derivative_composed_ratio_upper_using_full_tail_cap", "2")) < Decimal(1)):
        issues.append(issue("summary", "derivative-composed-ratio-not-below-one", repr(summary.get("derivative_composed_ratio_upper_using_full_tail_cap"))))
    if dec(summary.get("tail_actual_max_ratio_to_cap", "1")) >= Decimal("1e-1000"):
        issues.append(issue("summary", "tail-ratio-not-tiny", repr(summary.get("tail_actual_max_ratio_to_cap"))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "finite n<=30",
        "full 2.0e-3",
        "strictly below one first omitted",
        "worst-row finite-plus-tail numerator-source budget",
        "quadrature remainder",
        "all-row coverage",
        "remain open",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "worst row", "n>30 tail", "quadrature", "uniform collar", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[FinitePlusTailIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[FinitePlusTailIssue] = []
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
        "all rows are certified",
        "full residual numerator certificate",
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
    issues: list[FinitePlusTailIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, budget_certificate_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, budget_certificate_rows))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"FINITE-PLUS-TAIL {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
