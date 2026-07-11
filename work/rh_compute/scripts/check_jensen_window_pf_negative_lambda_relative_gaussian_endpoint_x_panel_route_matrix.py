#!/usr/bin/env python3
"""Validate the relative-Gaussian endpoint x-panel route matrix."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix import (  # noqa: E402
    DEFAULT_ENDPOINT_PARITY_JSON,
    DEFAULT_INTERPOLATION_ROUTE_JSON,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
    result_line,
)


REQUIRED_ROW_IDS = {
    "nlrgexrm_01_endpoint_route_import",
    "nlrgexrm_02_exact_x_change_of_variables",
    "nlrgexrm_03_first_panel_mass_budget",
    "nlrgexrm_04_x_bernstein_budget_table",
    "nlrgexrm_05_x_moment_taylor_obligation",
    "nlrgexrm_06_tiny_mass_promotion_rejected",
    "nlrgexrm_07_acceptance_gate",
}

ALLOWED_ROLES = {
    "scope_reduction",
    "change_of_variables",
    "budget_diagnostic",
    "route_budget",
    "open_requirement",
    "rejected_route",
    "acceptance_gate",
}

ALLOWED_READINESS = {
    "diagnostic_only",
    "open_requirement",
    "not_ready_to_apply",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint X-Panel Route Matrix",
    "Status: endpoint x-panel route matrix",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian endpoint x-panel route matrix: 7 rows, 0 issues, x interval 0<=x<=0.01, 18 Bernstein budgets, 0 ready-to-apply rows",
    "first panel mass upper: 1.6155560239464103501984383230311864446832165667162E-21",
    "first/heaviest mass ratio: 2.682927708866247294881825729652E-21",
    "first-panel value sup budget without Bernstein: 6.996592450057070652146864508995E-20",
    "rho=3 degree=32 value sup budget: 6.482413531562058989909765298732E-5",
    "2*T^(alpha+1)*x^(2*alpha+1)*exp(-T*x^2)/Gamma(alpha+1)",
    "T^(-k/2)*lower_gamma(alpha+1+k/2,1)/Gamma(alpha+1)",
    "Tiny first-panel Gamma mass is not a proof of the endpoint remainder.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class EndpointXPanelIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> EndpointXPanelIssue:
    return EndpointXPanelIssue(section=section, issue=name, detail=detail)


def dec(value: object) -> Decimal:
    return Decimal(str(value).replace("e", "E"))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[EndpointXPanelIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[EndpointXPanelIssue]:
    issues: list[EndpointXPanelIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_endpoint_x_panel_route_matrix":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "endpoint x-panel route matrix":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_interpolation_remainder_route_matrix",
        "source_interpolation_remainder_route_json",
        "source_endpoint_parity_repair_matrix",
        "source_endpoint_parity_repair_json",
        "source_arb_chebyshev_interpolant_moment_scout",
        "source_compact_interval_integration_scout",
        "source_quadrature_remainder_route_matrix",
        "source_finite_part_weighted_sum_interval_certificate",
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
        "does not prove an x-domain",
        "does not prove a sup-norm",
        "does not implement exact x moments",
        "interpolation remainder",
        "compact interval certificate",
        "finite-grid interval certificate",
        "uniform collar",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[EndpointXPanelIssue]:
    try:
        recomputed = build_artifact(DEFAULT_INTERPOLATION_ROUTE_JSON, DEFAULT_ENDPOINT_PARITY_JSON)
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[EndpointXPanelIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[EndpointXPanelIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[EndpointXPanelIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "not-list", repr(type(rows)))], 0, 0
    seen = {row.get("id") for row in rows if isinstance(row, dict)}
    missing = REQUIRED_ROW_IDS - seen
    extra = seen - REQUIRED_ROW_IDS
    if missing:
        issues.append(issue("matrix_rows", "missing-row-ids", repr(sorted(missing))))
    if extra:
        issues.append(issue("matrix_rows", "extra-row-ids", repr(sorted(extra))))
    ready_to_apply = 0
    open_requirements = 0
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
        if row.get("readiness") == "open_requirement":
            open_requirements += 1
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
    if open_requirements < 2:
        issues.append(issue("matrix_rows", "missing-open-requirements", str(open_requirements)))
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    rejected = by_id.get("nlrgexrm_06_tiny_mass_promotion_rejected", {})
    if rejected.get("diagnostics", {}).get("tiny_mass_is_proof") is not False:
        issues.append(issue("nlrgexrm_06_tiny_mass_promotion_rejected", "bad-rejection-flag", repr(rejected)))
    if "does not supply" not in str(rejected.get("gap", "")).lower():
        issues.append(issue("nlrgexrm_06_tiny_mass_promotion_rejected", "weak-gap", repr(rejected.get("gap"))))
    return issues, len(rows), ready_to_apply


def validate_summary(artifact: dict, row_count: int) -> list[EndpointXPanelIssue]:
    summary = artifact.get("summary", {})
    issues: list[EndpointXPanelIssue] = []
    expected = {
        "matrix_rows": 7,
        "T": 10000,
        "index": 21,
        "alpha": "20.5",
        "x_interval": "0<=x<=0.01",
        "y_interval": "0<=y<=1",
        "x_right_endpoint": "0.01",
        "compact_interval": "0<=y<=200",
        "panel_count": 6,
        "first_panel_mass_upper": "1.6155560239464103501984383230311864446832165667162E-21",
        "heaviest_panel": "20<=y<=50",
        "heaviest_panel_mass_upper": "0.60216159332489531950482867922076390781009134801886",
        "first_to_heaviest_mass_ratio": "2.682927708866247294881825729652E-21",
        "value_unscaled_expectation_error_cap": "6.782032247872604818E-40",
        "derivative_unscaled_expectation_error_cap": "1.424226772053247012E-38",
        "first_panel_value_sup_budget_without_bernstein": "6.996592450057070652146864508995E-20",
        "first_panel_derivative_sup_budget_without_bernstein": "1.469284414511984837177801591129E-18",
        "x_weight_power": 42,
        "degree_count": 6,
        "rho_count": 3,
        "bernstein_x_budget_row_count": 18,
        "rho_2_degree_32_value_sup_budget": "7.512533939108907946133043813769E-11",
        "rho_2_degree_32_derivative_sup_budget": "1.577632127212870668931635692769E-9",
        "rho_3_degree_32_value_sup_budget": "6.482413531562058989909765298732E-5",
        "rho_3_degree_32_derivative_sup_budget": "1.361306841628032388091331484753E-3",
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "tiny_mass_is_proof": False,
        "y_plane_branch_resolved": False,
        "x_panel_certificate_produced": False,
    }
    for key, expected_value in expected.items():
        if summary.get(key) != expected_value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {expected_value!r}"))
    if row_count != summary.get("matrix_rows"):
        issues.append(issue("summary", "row-count-mismatch", f"{row_count} != {summary.get('matrix_rows')!r}"))
    if dec(summary.get("first_panel_mass_upper", "1")) >= Decimal("1e-20"):
        issues.append(issue("summary", "first-panel-mass-too-large", repr(summary.get("first_panel_mass_upper"))))
    if dec(summary.get("first_to_heaviest_mass_ratio", "1")) >= Decimal("1e-20"):
        issues.append(issue("summary", "mass-ratio-too-large", repr(summary.get("first_to_heaviest_mass_ratio"))))
    if dec(summary.get("rho_3_degree_32_value_sup_budget", "0")) <= Decimal("1e-5"):
        issues.append(issue("summary", "rho3-degree32-budget-too-small", repr(summary.get("rho_3_degree_32_value_sup_budget"))))
    budget_rows = None
    for row in artifact.get("matrix_rows", []):
        if isinstance(row, dict) and row.get("id") == "nlrgexrm_04_x_bernstein_budget_table":
            budget_rows = row.get("diagnostics", {}).get("bernstein_x_budget_rows")
    if not isinstance(budget_rows, list) or len(budget_rows) != 18:
        issues.append(issue("summary", "bad-budget-rows", repr(budget_rows)))
        budget_rows = []
    pairs = {(row.get("degree"), row.get("rho")) for row in budget_rows if isinstance(row, dict)}
    for expected_pair in ((8, "1.5"), (16, "2.0"), (32, "3.0")):
        if expected_pair not in pairs:
            issues.append(issue("summary", "missing-budget-pair", repr(expected_pair)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "endpoint x-panel route matrix",
        "y=t*x^2",
        "0<=x<=0.01",
        "x^42",
        "numerically lightweight",
        "route matrix only",
        "no x-domain",
        "true x-panel interpolation-remainder",
    ):
        if required not in finding:
            issues.append(issue("summary", "weak-main-finding", required))
    required_input = str(summary.get("required_missing_input", "")).lower()
    for required in ("first-panel", "x-domain", "sup-norm", "taylor-model", "exact transformed moments", "splice"):
        if required not in required_input:
            issues.append(issue("summary", "weak-required-missing-input", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in (
        "ready_to_apply",
        "mass is not promoted",
        "y=0 branch",
        "bernstein x-panel budgets",
        "x moment formula",
        "lambda <= 0",
    ):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[EndpointXPanelIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[EndpointXPanelIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-required-string", required))
    lowered = text.lower()
    for forbidden in (
        "x-panel certificate is complete",
        "interpolation-remainder theorem is proved",
        "compact interval certificate is complete",
        "finite-grid interval certificate is complete",
        "uniform collar theorem is proved",
        "rh is proved",
        "lambda <= 0 is proved",
        "tiny first-panel gamma mass proves",
        "y=0 branch is resolved",
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
    issues: list[EndpointXPanelIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, _ready_to_apply = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count))
    issues.extend(validate_note(note_path))

    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"ENDPOINT-X-PANEL {item.section} [{item.issue}] {item.detail}")
        print(result_line(artifact).replace("0 issues", f"{len(issues)} issues"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
