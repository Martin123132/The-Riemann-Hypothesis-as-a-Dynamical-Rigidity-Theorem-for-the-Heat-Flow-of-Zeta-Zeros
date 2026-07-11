#!/usr/bin/env python3
"""Validate the worst-row far-tail split certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate import (  # noqa: E402
    DEFAULT_FIRST_OMITTED_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_QUADRATURE_ROUTE_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgwrftsc_01_split_and_monotonicity",
    "nlrgwrftsc_02_finite_phi_tail_bound",
    "nlrgwrftsc_03_polynomial_tail_moments",
    "nlrgwrftsc_04_tail_vs_quadrature_cap",
    "nlrgwrftsc_05_interval_integration_handoff",
    "nlrgwrftsc_06_far_tail_promotion_rejected",
    "nlrgwrftsc_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "arb_tail_certificate",
    "intervalization_handoff",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_tail_certificate",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Far-Tail Split Certificate",
    "Status: worst-row far-tail split certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian worst-row far-tail split certificate: 7 rows, 0 issues, split y=200, 2 tail ratios below quadrature cap, 0 ready-to-apply rows",
    "split y: 200",
    "remaining compact interval: 0<=y<=200",
    "value tail ratio below quadrature cap: True",
    "derivative tail ratio below quadrature cap: True",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class FarTailIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> FarTailIssue:
    return FarTailIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[FarTailIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[FarTailIssue]:
    issues: list[FarTailIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "worst-row far-tail split certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_first_omitted_denominator_certificate",
        "source_first_omitted_denominator_json",
        "source_quadrature_remainder_route_matrix",
        "source_quadrature_remainder_route_json",
        "source_node_c0_range_certificate",
        "source_phi_tail_grid_certificate",
        "source_cancellation_reduced_grid_scout",
        "source_intervalization_target",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "far-tail split certificate only",
        "does not integrate the compact interval",
        "all rows",
        "finite-grid interval certificate",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[FarTailIssue]:
    try:
        recomputed = build_artifact(DEFAULT_FIRST_OMITTED_JSON, DEFAULT_QUADRATURE_ROUTE_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[FarTailIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[FarTailIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[FarTailIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    tail_certificate_rows = 0
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
        if row.get("role") == "arb_tail_certificate":
            tail_certificate_rows += 1
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not ", "does not", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", repr(row.get("proof_boundary"))))
        refs = row.get("source_artifacts", [])
        if not isinstance(refs, list) or not refs:
            issues.append(issue(row_id, "missing-source-artifacts", repr(refs)))
        else:
            for ref in refs:
                issues.extend(validate_ref(row_id, ref))
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    handoff = rows_by_id.get("nlrgwrftsc_05_interval_integration_handoff", {}).get("diagnostics", {})
    if handoff.get("remaining_compact_interval") != "0<=y<=200":
        issues.append(issue("nlrgwrftsc_05_interval_integration_handoff", "bad-compact-interval", repr(handoff)))
    rejected = rows_by_id.get("nlrgwrftsc_06_far_tail_promotion_rejected", {})
    if "compact interval" not in str(rejected.get("gap", "")).lower():
        issues.append(issue("nlrgwrftsc_06_far_tail_promotion_rejected", "weak-rejection-gap", repr(rejected.get("gap"))))
    return issues, len(rows), tail_certificate_rows


def validate_summary(artifact: dict, row_count: int, tail_certificate_rows: int) -> list[FarTailIssue]:
    summary = artifact.get("summary", {})
    issues: list[FarTailIssue] = []
    expected = {
        "matrix_rows": 7,
        "T": 10000,
        "index": 21,
        "split_y": 200,
        "phi_term_count": 30,
        "polynomial_M": 20,
        "precision_bits": 512,
        "quadrature_ratio_radius_cap": "0.0000010",
        "value_tail_ratio_below_quadrature_cap": True,
        "derivative_tail_ratio_below_quadrature_cap": True,
        "value_unscaled_tail_below_quadrature_cap": True,
        "derivative_unscaled_tail_below_quadrature_cap": True,
        "remaining_compact_interval": "0<=y<=200",
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if tail_certificate_rows != 4:
        issues.append(issue("summary", "bad-tail-certificate-row-count", str(tail_certificate_rows)))
    if dec(summary.get("value_ratio_to_first_omitted_upper", "1")) >= Decimal("1e-20"):
        issues.append(issue("summary", "value-tail-ratio-too-large", repr(summary.get("value_ratio_to_first_omitted_upper"))))
    if dec(summary.get("derivative_ratio_to_first_omitted_upper", "1")) >= Decimal("1e-20"):
        issues.append(
            issue("summary", "derivative-tail-ratio-too-large", repr(summary.get("derivative_ratio_to_first_omitted_upper")))
        )
    if dec(summary.get("value_tail_to_unscaled_cap_ratio_upper", "1")) >= Decimal("1e-15"):
        issues.append(issue("summary", "value-tail-cap-ratio-too-large", repr(summary.get("value_tail_to_unscaled_cap_ratio_upper"))))
    if dec(summary.get("derivative_tail_to_unscaled_cap_ratio_upper", "1")) >= Decimal("1e-15"):
        issues.append(
            issue("summary", "derivative-tail-cap-ratio-too-large", repr(summary.get("derivative_tail_to_unscaled_cap_ratio_upper")))
        )
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "arb-certified far-tail split",
        "y>=200",
        "incomplete-gamma",
        "quadrature route matrix caps",
        "compact integration",
        "remain open",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "far tail", "compact interval", "all-row", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[FarTailIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[FarTailIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "compact interval-integration certificate is complete",
        "quadrature-remainder theorem is proved",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
        "rh is proved",
        "lambda <= 0 is proved",
        "far-tail certificate proves the quadrature-remainder source",
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
    issues: list[FarTailIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, tail_certificate_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, tail_certificate_rows))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"FAR-TAIL {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
