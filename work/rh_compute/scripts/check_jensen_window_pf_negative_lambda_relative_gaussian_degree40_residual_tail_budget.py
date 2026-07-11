#!/usr/bin/env python3
"""Validate the degree-40 residual-tail budget artifact."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgd40rtb_01_residual_coordinate",
    "nlrgd40rtb_02_degree40_margin_extraction",
    "nlrgd40rtb_03_value_residual_sufficient_budget",
    "nlrgd40rtb_04_derivative_residual_sufficient_budget",
    "nlrgd40rtb_05_finite_tail_profile_through_degree80",
    "nlrgd40rtb_06_live_residual_tail_theorem",
    "nlrgd40rtb_07_finite_profile_promotion_rejected",
    "nlrgd40rtb_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_interval_diagnostic",
    "exact_sufficient_condition",
    "finite_diagnostic",
    "live_route",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-40 Residual Tail Budget",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget`",
    "prove the residual estimates.",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: 8 rows, 0 issues, 5 budget inequalities, 4 finite tail profile rows, 0 ready-to-apply rows",
    "F_i(u) = P_i^(40)(u) + R_i(u),  i in {21,22,23,24}",
    "|R_i(u)|  <= 5.382819486765314521E-01 * u^3",
    "|R_i'(u)| <= 9.315354075509573936E-03 * u",
    "The value-tail budget is limited by the companion product.",
    "companion_product_value_residual: allocated=8.612511187827506021E+00",
    "weighted_gap_derivative_residual: allocated=6.856100606742001169E+00",
    "F_24, degrees 42..80: value fraction=[0.000931556974617909435 +/- 2.12e-22], derivative fraction=[0.000995886692207194084 +/- 4.49e-22]",
    "The finite profile is a plausibility diagnostic only.",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class ResidualTailBudgetIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ResidualTailBudgetIssue:
    return ResidualTailBudgetIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ResidualTailBudgetIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[ResidualTailBudgetIssue]:
    issues: list[ResidualTailBudgetIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_degree40_ladder_stress",
        "source_uniform_remainder_target",
        "source_stencil_remainder_obligations",
        "source_dependency_graph",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite",
        "sufficient",
        "residual-tail targets",
        "does not prove",
        "infinite tail",
        "scaled-curvature",
        "cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[ResidualTailBudgetIssue]:
    params = artifact.get("matrix_rows", [{}, {}])[1].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            int(params.get("finite_degree", 40)),
            int(params.get("profile_max_degree", 80)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 384)),
            int(params.get("tail_start_k", 22)),
            int(params.get("collar_start_T", 1156)),
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[ResidualTailBudgetIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_diagnostics(artifact: dict) -> list[ResidualTailBudgetIssue]:
    diagnostics = artifact.get("matrix_rows", [{}, {}])[1].get("diagnostics", {})
    params = diagnostics.get("parameters", {})
    issues: list[ResidualTailBudgetIssue] = []
    expected_params = {
        "finite_degree": 40,
        "profile_max_degree": 80,
        "tail_cutoff_n": 80,
        "precision_bits": 384,
        "tail_start_k": 22,
        "collar_start_T": 1156,
        "real_interval_u": "[0, 1/1156]",
        "real_interval_T": "[1156, infinity)",
        "residual_first_j": 21,
        "finite_tail_profile_last_j": 40,
        "weighted_gap_abs_weight_sum": 368,
    }
    if params != expected_params:
        issues.append(issue("parameters", "bad-parameters", repr(params)))
    expected_margins = {
        "normalizer_min_lower": "[0.439627934927551258358085649181 +/- 2.62e-31]",
        "B_product_min_lower": "[37.2033572903651204156102209610 +/- 2.68e-29]",
        "companion_product_min_lower": "[17.2250223756550133023672457518 +/- 2.11e-29]",
        "weighted_gap_derivative_min_lower": "[27.4244024269680055723947152492 +/- 3.52e-29]",
        "finite_normalizer_abs_upper": "[1.00000000000000000000000000000 +/- 1e-34]",
        "finite_derivative_abs_upper": "[917.618341545662737252398795656 +/- 4.88e-28]",
    }
    if diagnostics.get("finite_collar_margins") != expected_margins:
        issues.append(issue("finite_collar_margins", "bad-margins", repr(diagnostics.get("finite_collar_margins"))))
    expected_scalars = {
        "value_residual_half_safety_budget_A": "5.382819486765314521E-01",
        "derivative_residual_half_safety_budget_B": "9.315354075509573936E-03",
        "limiting_value_budget": "companion_product",
        "budget_inequality_row_count": 5,
        "finite_tail_profile_row_count": 4,
    }
    for key, value in expected_scalars.items():
        if diagnostics.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{diagnostics.get(key)!r} != {value!r}"))
    expected_raw_thresholds = {
        "normalizer": "3.395695876365208626E+08",
        "B_product": "5.375875774509763687E+03",
        "companion_product": "1.076563897353062904E+00",
        "weighted_gap_value_part": "9.044023628251490976E+00",
    }
    if diagnostics.get("raw_value_thresholds") != expected_raw_thresholds:
        issues.append(issue("raw_value_thresholds", "bad-thresholds", repr(diagnostics.get("raw_value_thresholds"))))
    budget_names = [row.get("name") for row in diagnostics.get("budget_inequality_rows", [])]
    if budget_names != [
        "normalizer_value_residual",
        "B_product_value_residual",
        "companion_product_value_residual",
        "weighted_gap_value_residual",
        "weighted_gap_derivative_residual",
    ]:
        issues.append(issue("budget_inequality_rows", "bad-row-order", repr(budget_names)))
    limiting_names = [row.get("name") for row in diagnostics.get("budget_inequality_rows", []) if row.get("limiting") is True]
    if limiting_names != ["companion_product_value_residual", "weighted_gap_derivative_residual"]:
        issues.append(issue("budget_inequality_rows", "bad-limiting-rows", repr(limiting_names)))
    profile_rows = diagnostics.get("finite_tail_profile_rows", [])
    if [row.get("index") for row in profile_rows] != [21, 22, 23, 24]:
        issues.append(issue("finite_tail_profile_rows", "bad-index-order", repr([row.get("index") for row in profile_rows])))
    for row in profile_rows:
        if row.get("term_degree_range") != "42..80":
            issues.append(issue(f"F_{row.get('index')}", "bad-term-degree-range", repr(row.get("term_degree_range"))))
        if row.get("largest_value_term_degree") != 42 or row.get("largest_derivative_term_degree") != 42:
            issues.append(issue(f"F_{row.get('index')}", "bad-largest-term-degree", repr(row)))
    if diagnostics.get("max_finite_profile_value_budget_fraction") >= 0.001:
        issues.append(issue("finite_tail_profile", "value-profile-too-large", repr(diagnostics.get("max_finite_profile_value_budget_fraction"))))
    if diagnostics.get("max_finite_profile_derivative_budget_fraction") >= 0.001:
        issues.append(issue("finite_tail_profile", "derivative-profile-too-large", repr(diagnostics.get("max_finite_profile_derivative_budget_fraction"))))
    boundary = str(diagnostics.get("proof_boundary_note", "")).lower()
    for required in ("sufficient", "do not prove", "infinite residual tail"):
        if required not in boundary:
            issues.append(issue("diagnostics", "weak-proof-boundary-note", required))
    return issues


def validate_rows(artifact: dict) -> tuple[list[ResidualTailBudgetIssue], int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[ResidualTailBudgetIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_or_sufficient = 0
    finite_diagnostics = 0
    route_rows = 0
    ready_to_apply = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        role = row.get("role")
        if role not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(role)))
        if role in {"exact_reduction", "exact_sufficient_condition"}:
            exact_or_sufficient += 1
        if role in {"finite_interval_diagnostic", "finite_diagnostic"}:
            finite_diagnostics += 1
        if role in {"live_route", "rejected_route"}:
            route_rows += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") not in {"not_ready_to_apply", "available_exact"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("not", "only", "sufficient", "rejected", "proof-hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    if ready_to_apply:
        issues.append(issue("matrix_rows", "unexpected-ready-to-apply", str(ready_to_apply)))
    return issues, len(rows), exact_or_sufficient, finite_diagnostics, route_rows


def validate_summary(
    artifact: dict,
    row_count: int,
    exact_or_sufficient: int,
    finite_diagnostics: int,
    route_rows: int,
) -> list[ResidualTailBudgetIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "budget_inequality_rows": 5,
        "finite_tail_profile_rows": 4,
        "value_residual_half_safety_budget_A": "5.382819486765314521E-01",
        "derivative_residual_half_safety_budget_B": "9.315354075509573936E-03",
        "limiting_value_budget": "companion_product",
        "max_finite_profile_value_budget_fraction": "9.315569746179094314E-04",
        "max_finite_profile_derivative_budget_fraction": "9.958866922071939881E-04",
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[ResidualTailBudgetIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_or_sufficient != 3:
        issues.append(issue("summary", "bad-exact-sufficient-count", str(exact_or_sufficient)))
    if finite_diagnostics != 2:
        issues.append(issue("summary", "bad-finite-diagnostic-count", str(finite_diagnostics)))
    if route_rows != 2:
        issues.append(issue("summary", "bad-route-row-count", str(route_rows)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in (
        "degree-40 arb collar margins",
        "0<=u<=1/1156",
        "half-safety a",
        "half-safety b",
        "i=21..24",
        "less than 0.1%",
        "no analytic majorant",
        "infinite residual tail",
    ):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "not proved", "degree-42..80", "k=22", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ResidualTailBudgetIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ResidualTailBudgetIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "residual estimates are proved",
        "infinite taylor-tail theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[ResidualTailBudgetIssue], dict]:
    artifact = load_json(target_path)
    issues: list[ResidualTailBudgetIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_diagnostics(artifact))
    row_issues, row_count, exact_or_sufficient, finite_diagnostics, route_rows = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_or_sufficient, finite_diagnostics, route_rows))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    target = args.target if args.target.is_absolute() else REPO_ROOT / args.target
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    issues, summary = validate(target, note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": [asdict(item) for item in issues]}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-DEG40-RESIDUAL-TAIL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('budget_inequality_rows')} budget inequalities, "
            f"{summary.get('finite_tail_profile_rows')} finite tail profile rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
