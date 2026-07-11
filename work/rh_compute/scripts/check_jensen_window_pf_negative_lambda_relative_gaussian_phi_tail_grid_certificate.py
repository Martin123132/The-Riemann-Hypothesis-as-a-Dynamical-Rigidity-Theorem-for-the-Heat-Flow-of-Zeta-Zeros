#!/usr/bin/env python3
"""Validate the relative-Gaussian finite-grid Phi-tail source certificate."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate import (  # noqa: E402
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NODE_C0_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_PHI_TAIL_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgptgc_01_conditional_tail_scout_import",
    "nlrgptgc_02_side_conditions_certified",
    "nlrgptgc_03_finite_grid_tail_source_certificate",
    "nlrgptgc_04_intervalization_obligation_handoff",
    "nlrgptgc_05_full_phi_evaluation_promotion_rejected",
    "nlrgptgc_06_acceptance_gate",
}

ALLOWED_ROLES = {
    "analytic_tail_bound",
    "finite_side_condition_certificate",
    "finite_tail_source_certificate",
    "intervalization_handoff",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "available_conditional",
    "available_for_phi_tail_source",
    "available_for_intervalization",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Phi-Tail Grid Certificate",
    "Status: finite-grid Phi-tail source certificate",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: 6 rows, 0 issues, 3 certified tail sources, 2 certified side conditions, 0 ready-to-apply rows",
    "node range x<=1 certified: True",
    "certified c0 lower: 0.44",
    "finite-grid tail source certified: True",
    "Phi value n>30 tail",
    "Phi prime n>30 derivative-core tail",
    "Phi(0) n>30 denominator tail",
    "individual Laguerre node and weight intervals beyond the coarse x<=1 range",
    "finite n<=30 interval evaluation of Phi and Phi' at certified node intervals",
    "generalized Gauss-Laguerre quadrature-remainder or interval adaptive integration",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_node_c0_range_certificate.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class PhiTailGridIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> PhiTailGridIssue:
    return PhiTailGridIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[PhiTailGridIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[PhiTailGridIssue]:
    issues: list[PhiTailGridIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite-grid Phi-tail source certificate":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_phi_tail_bound_scout",
        "source_phi_tail_bound_json",
        "source_node_c0_range_certificate",
        "source_node_c0_range_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_cancellation_reduced_grid_scout",
        "source_cancellation_reduced_grid_json",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "phi-tail source certificate only",
        "does not provide individual laguerre",
        "does not certify finite n<=30",
        "does not bound quadrature",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(
    artifact: dict,
    phi_tail_path: Path,
    node_c0_path: Path,
    interval_path: Path,
) -> list[PhiTailGridIssue]:
    try:
        recomputed = build_artifact(phi_tail_path, node_c0_path, interval_path)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[PhiTailGridIssue] = []
    for key in ("tail_source_rows", "matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_tail_rows(artifact: dict) -> list[PhiTailGridIssue]:
    rows = artifact.get("tail_source_rows", [])
    issues: list[PhiTailGridIssue] = []
    if not isinstance(rows, list):
        return [issue("tail_source_rows", "bad-rows", repr(type(rows)))]
    expected_channels = [
        "Phi value n>30 tail",
        "Phi prime n>30 derivative-core tail",
        "Phi(0) n>30 denominator tail",
    ]
    if [row.get("channel") for row in rows if isinstance(row, dict)] != expected_channels:
        issues.append(issue("tail_source_rows", "bad-channel-order", repr([row.get("channel") for row in rows])))
    if len(rows) != 3:
        issues.append(issue("tail_source_rows", "bad-row-count", str(len(rows))))
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("tail_source_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("channel", "<missing-channel>"))
        for key in (
            "tail_bound",
            "normalized_or_relative_bound",
            "per_source_cap",
            "ratio_to_per_source_cap",
            "certified_below_cap",
            "proof_boundary",
        ):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("per_source_cap") != "2.000000000000000000E-3":
            issues.append(issue(row_id, "bad-per-source-cap", repr(row.get("per_source_cap"))))
        if row.get("certified_below_cap") is not True:
            issues.append(issue(row_id, "not-certified-below-cap", repr(row)))
        if not dec(row.get("ratio_to_per_source_cap", "1")) < Decimal("1e-1000"):
            issues.append(issue(row_id, "ratio-not-tiny", repr(row.get("ratio_to_per_source_cap"))))
        boundary = str(row.get("proof_boundary", "")).lower()
        for required in ("tail source", "not a laguerre", "not a quadrature", "not a full residual"):
            if required not in boundary:
                issues.append(issue(row_id, "weak-row-boundary", required))
    return issues


def validate_summary(artifact: dict) -> list[PhiTailGridIssue]:
    summary = artifact.get("summary", {})
    issues: list[PhiTailGridIssue] = []
    expected = {
        "matrix_rows": 6,
        "tail_source_rows": 3,
        "certified_side_conditions": 2,
        "ready_to_apply_rows": 0,
        "source_phi_tail_rows": 6,
        "source_node_c0_rows": 5,
        "source_interval_obligations": 8,
        "padded_tail_range": "0<=x<=1",
        "node_range_x_le_1_certified": True,
        "certified_c0_lower": "0.44",
        "c0_lower_certified_by_n1_term": True,
        "per_source_intervalization_cap": "2.000000000000000000E-3",
        "max_tail_source_ratio_to_cap": "3.778846084314790602E-1292",
        "value_ratio_to_per_source_cap": "1.146230845771967159E-1297",
        "derivative_ratio_to_per_source_cap": "3.778846084314790602E-1292",
        "c0_ratio_to_per_source_cap": "1.415251175348174758E-1301",
        "tail_sources_below_per_source_cap": True,
        "finite_grid_tail_source_certified": True,
        "remaining_obligation_count": 5,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    for key in ("max_tail_source_ratio_to_cap", "value_ratio_to_per_source_cap", "derivative_ratio_to_per_source_cap", "c0_ratio_to_per_source_cap"):
        if not dec(summary.get(key, "1")) < Decimal("1e-1000"):
            issues.append(issue("summary", f"{key}-not-tiny", repr(summary.get(key))))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "certificate-composed",
        "n>30",
        "2.0e-3",
        "x<=1",
        "phi(0)>=0.44",
        "finite n<=30",
        "quadrature error",
        "grid-to-collar",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "n>30", "finite n<=30", "quadrature", "uniform collar", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_matrix_rows(artifact: dict) -> tuple[list[PhiTailGridIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[PhiTailGridIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    ready_to_apply = 0
    available_for_intervalization = 0
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
        if row.get("readiness") not in ALLOWED_READINESS:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        if row.get("readiness") == "available_for_intervalization":
            available_for_intervalization += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "rejected", "handoff")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    side = rows_by_id.get("nlrgptgc_02_side_conditions_certified", {}).get("diagnostics", {})
    if side.get("side_conditions_certified") is not True:
        issues.append(issue("nlrgptgc_02_side_conditions_certified", "side-conditions-not-certified", repr(side)))
    if side.get("certified_worst_x_square_upper_bound") != "809/1156":
        issues.append(issue("nlrgptgc_02_side_conditions_certified", "bad-x-square-bound", repr(side)))
    handoff = rows_by_id.get("nlrgptgc_04_intervalization_obligation_handoff", {}).get("diagnostics", {})
    if handoff.get("source_obligation_id") != "nlrgit_03_phi_and_c0_interval_tail":
        issues.append(issue("nlrgptgc_04_intervalization_obligation_handoff", "bad-source-obligation", repr(handoff)))
    if len(handoff.get("remaining_obligations", [])) != 5:
        issues.append(issue("nlrgptgc_04_intervalization_obligation_handoff", "bad-remaining-obligations", repr(handoff)))
    return issues, len(rows), available_for_intervalization


def validate_note(path: Path) -> list[PhiTailGridIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[PhiTailGridIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "finite-grid interval certificate is complete",
        "quadrature-remainder theorem is proved",
        "uniform collar theorem is proved",
        "full residual numerator certificate",
        "finite n<=30 node evaluations are certified",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(
    target_path: Path,
    note_path: Path,
    phi_tail_path: Path,
    node_c0_path: Path,
    interval_path: Path,
) -> tuple[list[PhiTailGridIssue], dict]:
    artifact = load_json(target_path)
    issues: list[PhiTailGridIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, phi_tail_path, node_c0_path, interval_path))
    issues.extend(validate_tail_rows(artifact))
    issues.extend(validate_summary(artifact))
    row_issues, row_count, available_for_intervalization = validate_matrix_rows(artifact)
    issues.extend(row_issues)
    if row_count != 6:
        issues.append(issue("matrix_rows", "bad-row-count", str(row_count)))
    if available_for_intervalization != 2:
        issues.append(issue("matrix_rows", "bad-available-for-intervalization-count", str(available_for_intervalization)))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--node-c0-json", type=Path, default=DEFAULT_NODE_C0_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    phi_tail_path = args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    node_c0_path = args.node_c0_json if args.node_c0_json.is_absolute() else REPO_ROOT / args.node_c0_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    issues, summary = validate(target, note, phi_tail_path, node_c0_path, interval_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-PHI-TAIL-GRID {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian Phi-tail grid certificate: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('tail_source_rows')} certified tail sources, "
            f"{summary.get('certified_side_conditions')} certified side conditions, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
