#!/usr/bin/env python3
"""Validate the relative-Gaussian endpoint parity-repair matrix."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix import (  # noqa: E402
    DEFAULT_INTERPOLATION_ROUTE_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_PHI_TAIL_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgeprm_01_endpoint_branch_obligation",
    "nlrgeprm_02_finite_truncation_odd_taylor_audit",
    "nlrgeprm_03_low_order_branch_size_scout",
    "nlrgeprm_04_infinite_even_core_tail_charge_route",
    "nlrgeprm_05_x_variable_endpoint_route",
    "nlrgeprm_06_near_evenness_promotion_rejected",
    "nlrgeprm_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "route_constraint",
    "arb_series_diagnostic",
    "live_route",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "open_requirement",
    "diagnostic_only",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint Parity-Repair Matrix",
    "Status: endpoint parity-repair route matrix",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix: 7 rows, 0 issues, 8 odd Taylor rows, order 15, 0 ready-to-apply rows",
    "first odd coefficient abs upper: 1.5013918057273630212637957520005850517148059822056E-1300",
    "max low-order odd degree: 15",
    "max low-order odd abs upper: 1.5394049948039318674157069296860766261840406690307E-1255",
    "Near-evenness and low-order cancellation are not an endpoint analyticity proof.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class EndpointParityIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> EndpointParityIssue:
    return EndpointParityIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[EndpointParityIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[EndpointParityIssue]:
    issues: list[EndpointParityIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "endpoint parity-repair route matrix":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_interpolation_remainder_route_matrix",
        "source_interpolation_remainder_route_json",
        "source_arb_chebyshev_interpolant_moment_scout",
        "source_finite_part_weighted_sum_interval_certificate",
        "source_phi_tail_grid_certificate",
        "source_phi_tail_grid_json",
        "source_compact_interval_integration_scout",
        "source_intervalization_target",
        "source_uniform_remainder_target",
        "source_dependency_graph",
        "source_formal_tail_obstruction",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "route matrix only",
        "does not prove exact evenness",
        "does not certify the finite-truncation tail charge",
        "does not produce an x-variable",
        "interpolation remainder theorem",
        "finite-grid interval certificate",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[EndpointParityIssue]:
    try:
        recomputed = build_artifact(DEFAULT_INTERPOLATION_ROUTE_JSON, DEFAULT_PHI_TAIL_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[EndpointParityIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[EndpointParityIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[EndpointParityIssue] = []
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
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    rejected = by_id.get("nlrgeprm_06_near_evenness_promotion_rejected", {})
    if rejected.get("diagnostics", {}).get("near_evenness_is_proof") is not False:
        issues.append(issue("nlrgeprm_06_near_evenness_promotion_rejected", "bad-rejection-flag", repr(rejected)))
    if "near-evenness" not in str(rejected.get("gap", "")).lower():
        issues.append(issue("nlrgeprm_06_near_evenness_promotion_rejected", "weak-gap", repr(rejected.get("gap"))))
    return issues, len(rows), live_routes, rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    live_routes: int,
    rejected_routes: int,
) -> list[EndpointParityIssue]:
    summary = artifact.get("summary", {})
    issues: list[EndpointParityIssue] = []
    expected = {
        "matrix_rows": 7,
        "T": 10000,
        "compact_interval": "0<=y<=200",
        "heaviest_panel": "20<=y<=50",
        "phi_term_count": 30,
        "series_order": 15,
        "precision_bits": 8192,
        "odd_taylor_row_count": 8,
        "first_odd_coefficient_abs_upper": "1.5013918057273630212637957520005850517148059822056E-1300",
        "max_low_order_odd_degree": 15,
        "max_low_order_odd_abs_upper": "1.5394049948039318674157069296860766261840406690307E-1255",
        "live_route_rows": 2,
        "rejected_route_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "near_evenness_is_proof": False,
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
    if dec(summary.get("first_odd_coefficient_abs_upper", "1")) >= Decimal("1e-1200"):
        issues.append(issue("summary", "first-odd-too-large", repr(summary.get("first_odd_coefficient_abs_upper"))))
    if dec(summary.get("max_low_order_odd_abs_upper", "1")) >= Decimal("1e-1200"):
        issues.append(issue("summary", "max-odd-too-large", repr(summary.get("max_low_order_odd_abs_upper"))))
    odd_rows = None
    for row in artifact.get("matrix_rows", []):
        if isinstance(row, dict) and row.get("id") == "nlrgeprm_02_finite_truncation_odd_taylor_audit":
            odd_rows = row.get("diagnostics", {}).get("odd_taylor_rows")
    if not isinstance(odd_rows, list) or len(odd_rows) != 8:
        issues.append(issue("summary", "bad-odd-rows", repr(odd_rows)))
        odd_rows = []
    degrees = [row.get("degree") for row in odd_rows if isinstance(row, dict)]
    if degrees != [1, 3, 5, 7, 9, 11, 13, 15]:
        issues.append(issue("summary", "bad-odd-degrees", repr(degrees)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "endpoint parity-repair matrix",
        "first-panel branch obligation",
        "finite n<=30",
        "odd coefficient",
        "tail-sized",
        "does not prove",
        "exact evenness",
        "x variable",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in (
        "ready_to_apply",
        "near-evenness",
        "finite phi truncation",
        "y=0 endpoint branch",
        "x-variable",
        "lambda <= 0",
    ):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[EndpointParityIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[EndpointParityIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "endpoint analyticity is proved",
        "interpolation-remainder theorem is proved",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
        "rh is proved",
        "lambda <= 0 is proved",
        "near-evenness proves endpoint analyticity",
        "finite phi truncation is exactly even",
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
    issues: list[EndpointParityIssue] = []
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
            print(f"ENDPOINT-PARITY {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
