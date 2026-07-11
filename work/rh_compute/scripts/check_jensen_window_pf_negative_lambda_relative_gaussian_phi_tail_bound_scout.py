#!/usr/bin/env python3
"""Validate the relative-Gaussian Phi-tail bound scout."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout import (  # noqa: E402
    DEFAULT_GRID_JSON,
    DEFAULT_INTERVAL_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgtb_01_padded_node_range_scout",
    "nlrgtb_02_value_phi_tail_majorant",
    "nlrgtb_03_derivative_phi_prime_tail_majorant",
    "nlrgtb_04_c0_tail_and_lower_proxy",
    "nlrgtb_05_intervalization_handoff",
    "nlrgtb_06_full_obligation_promotion_rejected",
}

ALLOWED_ROLES = {
    "finite_node_range_scout",
    "analytic_tail_bound",
    "denominator_tail_scout",
    "conditional_budget_handoff",
    "rejected_route",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Phi Tail Bound Scout",
    "Status: analytic padded-range tail scout",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: 6 rows, 0 issues, 3 tail bounds below 1e-1000, 2 conditional requirements, 0 ready-to-apply rows",
    "max observed node x: 8.224376518739603403E-01",
    "padded x range: 1.000000000000000000E+00",
    "tail start n: 31",
    "tail bounds below per-source cap: True",
    "Replace the floating SciPy node range by interval enclosures proving x<=1 for all grid nodes.",
    "Replace the floating c0 lower proxy by an interval-certified lower bound Phi(0)>=0.44.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_actual_endpoint_remainder_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_asymptotic_remainder_target.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class PhiTailIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> PhiTailIssue:
    return PhiTailIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("E+00", "E+0"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[PhiTailIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[PhiTailIssue]:
    issues: list[PhiTailIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_bound_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "analytic padded-range tail scout":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_cancellation_reduced_grid_scout",
        "source_cancellation_reduced_grid_json",
        "source_intervalization_target",
        "source_intervalization_target_json",
        "source_actual_endpoint_scout",
        "source_asymptotic_remainder_target",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "padded-range tail scout only",
        "conditionally",
        "does not interval-certify",
        "does not bound quadrature",
        "uniform collar",
        "scaled-curvature",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict, grid_path: Path, interval_path: Path) -> list[PhiTailIssue]:
    try:
        recomputed = build_artifact(grid_path, interval_path)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[PhiTailIssue] = []
    for key in ("matrix_rows", "summary", "tail_diagnostics", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_summary(artifact: dict) -> list[PhiTailIssue]:
    summary = artifact.get("summary", {})
    issues: list[PhiTailIssue] = []
    expected = {
        "matrix_rows": 6,
        "tail_bound_rows": 3,
        "conditional_requirement_rows": 2,
        "ready_to_apply_rows": 0,
        "source_grid_rows": 20,
        "source_t_grid_count": 5,
        "source_index_count": 4,
        "source_selected_quadrature_order": 192,
        "node_range_row_count": 80,
        "max_observed_x": "8.224376518739603403E-01",
        "padded_x_range_upper": "1.000000000000000000E+00",
        "all_observed_nodes_inside_padded_range": True,
        "phi_truncation_n": 30,
        "tail_start_n": 31,
        "per_source_intervalization_cap": "2.000000000000000000E-3",
        "tail_bounds_below_per_source_cap": True,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))

    numeric_bounds = {
        "value_phi_tail_bound": Decimal("1e-1000"),
        "derivative_phi_prime_tail_bound": Decimal("1e-1000"),
        "c0_tail_bound": Decimal("1e-1000"),
        "normalized_value_tail_bound_using_c0_proxy": Decimal("1e-1000"),
        "normalized_derivative_core_tail_bound_using_c0_proxy": Decimal("1e-1000"),
        "denominator_relative_tail_bound_using_c0_proxy": Decimal("1e-1000"),
    }
    for key, limit in numeric_bounds.items():
        try:
            value = dec(summary.get(key))
        except Exception as exc:
            issues.append(issue("summary", f"bad-decimal-{key}", f"{summary.get(key)!r}: {exc}"))
            continue
        if not value < limit:
            issues.append(issue("summary", f"{key}-not-small", f"{value} >= {limit}"))

    if not (Decimal("0.82") < dec(summary.get("max_observed_x")) < Decimal("0.83")):
        issues.append(issue("summary", "bad-max-observed-x-range", repr(summary.get("max_observed_x"))))
    if not (Decimal("0.17") < dec(summary.get("observed_slack_to_padded_range")) < Decimal("0.18")):
        issues.append(issue("summary", "bad-observed-slack", repr(summary.get("observed_slack_to_padded_range"))))

    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "padded range",
        "n>30",
        "2.0e-3",
        "c0 lower proxy 0.44",
        "interval-prove",
        "phi(0)>=0.44",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    return issues


def validate_tail_diagnostics(artifact: dict) -> list[PhiTailIssue]:
    tail = artifact.get("tail_diagnostics", {})
    issues: list[PhiTailIssue] = []
    expected = {
        "phi_truncation_n": 30,
        "tail_start_n": 31,
        "x_range_upper": "1.000000000000000000E+00",
        "value_tail_majorant_degree": 4,
        "derivative_tail_majorant_degree": 6,
        "c0_lower_proxy": "4.400000000000000000E-01",
    }
    for key, value in expected.items():
        if tail.get(key) != value:
            issues.append(issue("tail_diagnostics", f"bad-{key}", f"{tail.get(key)!r} != {value!r}"))
    for key in (
        "value_tail_geometric_ratio_bound",
        "derivative_tail_geometric_ratio_bound",
        "c0_tail_geometric_ratio_bound",
    ):
        if not dec(tail.get(key)) < Decimal("1e-80"):
            issues.append(issue("tail_diagnostics", f"{key}-too-large", repr(tail.get(key))))
    if not dec(tail.get("c0_truncated_through_n30")) > Decimal("0.446"):
        issues.append(issue("tail_diagnostics", "c0-truncation-too-small", repr(tail.get("c0_truncated_through_n30"))))
    if not dec(tail.get("c0_truncated_margin_over_proxy")) > Decimal("0.006"):
        issues.append(issue("tail_diagnostics", "c0-proxy-margin-too-small", repr(tail.get("c0_truncated_margin_over_proxy"))))
    boundary = str(tail.get("proof_boundary", "")).lower()
    for required in ("analytic", "0<=x<=1", "floating proxy", "interval promotion"):
        if required not in boundary:
            issues.append(issue("tail_diagnostics", "weak-proof-boundary", required))
    return issues


def validate_rows(artifact: dict) -> list[PhiTailIssue]:
    rows = artifact.get("matrix_rows", [])
    issues: list[PhiTailIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))]
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    if len(rows) != 6:
        issues.append(issue("matrix_rows", "bad-row-count", str(len(rows))))
    ready_to_apply = 0
    role_counts: dict[str, int] = {}
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary", "source_artifacts"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        role = row.get("role")
        role_counts[str(role)] = role_counts.get(str(role), 0) + 1
        if role not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(role)))
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "conditional", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    if role_counts.get("analytic_tail_bound") != 2:
        issues.append(issue("matrix_rows", "bad-analytic-tail-row-count", repr(role_counts)))
    handoff = rows_by_id.get("nlrgtb_05_intervalization_handoff", {}).get("diagnostics", {})
    if handoff.get("source_obligation_id") != "nlrgit_03_phi_and_c0_interval_tail":
        issues.append(issue("nlrgtb_05_intervalization_handoff", "bad-source-obligation", repr(handoff)))
    if len(handoff.get("conditional_requirements", [])) != 2:
        issues.append(issue("nlrgtb_05_intervalization_handoff", "bad-conditional-requirements", repr(handoff)))
    if len(handoff.get("remaining_non_tail_obligations", [])) != 5:
        issues.append(issue("nlrgtb_05_intervalization_handoff", "bad-remaining-obligations", repr(handoff)))
    return issues


def validate_note(path: Path) -> list[PhiTailIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[PhiTailIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "finite-grid interval certificate is complete",
        "uniform residual estimate is proved",
        "quadrature-remainder theorem is proved",
        "node interval certificate is complete",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path, grid_path: Path, interval_path: Path) -> tuple[list[PhiTailIssue], dict]:
    artifact = load_json(target_path)
    issues: list[PhiTailIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact, grid_path, interval_path))
    issues.extend(validate_summary(artifact))
    issues.extend(validate_tail_diagnostics(artifact))
    issues.extend(validate_rows(artifact))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    issues, summary = validate(target, note, grid_path, interval_path)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-PHI-TAIL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian Phi tail bound scout: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('tail_bound_rows')} tail bounds below 1e-1000, "
            f"{summary.get('conditional_requirement_rows')} conditional requirements, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
