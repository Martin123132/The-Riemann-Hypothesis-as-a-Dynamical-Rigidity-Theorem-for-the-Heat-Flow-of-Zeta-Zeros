#!/usr/bin/env python3
"""Validate the relative-Gaussian pointwise tail budget diagnostic."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
    decimal_format,
    relative_tail_ratio_from_log_bound,
)


REQUIRED_ROW_IDS = {
    "nlrgptb_01_pointwise_log_tail_envelope",
    "nlrgptb_02_stencil_weight_sums",
    "nlrgptb_03_uniform_eta_sufficient_condition",
    "nlrgptb_04_half_safety_eta_budget",
    "nlrgptb_05_relative_tail_ratio_conversion",
    "nlrgptb_06_bottleneck_diagnostic",
    "nlrgptb_07_live_pointwise_tail_theorem",
    "nlrgptb_08_finite_budget_promotion_rejected",
    "nlrgptb_09_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "exact_sufficient_condition",
    "finite_diagnostic",
    "live_route",
    "rejected_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Pointwise Tail Budget",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: 9 rows, 0 issues, 4 positive baseline rows, 4 budget rows, 0 ready-to-apply rows",
    "|E_B| <= eta_(k-1)+2*eta_k+eta_(k+1)",
    "|E_U| <= eta_(k-1)+3*eta_k+3*eta_(k+1)+eta_(k+2)",
    "|E_C| <= (2*k+1)*eta_(k-1)+(6*k+5)*eta_k+(6*k+7)*eta_(k+1)+(2*k+3)*eta_(k+2)",
    "weighted-gap sum: 16*k+16",
    "at k=22, weighted-gap sum: 368",
    "limiting stencil counts: {'companion': 3, 'weighted_gap': 1}",
    "weakest half-safety eta: 1.458113205526978052E-9 at nlrgts_M6_T2000",
    "relative tail ratio bound: 1.458113204463930993E-9",
    "rho <= 1-exp(-eta)",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
)


@dataclass(frozen=True)
class PointwiseTailBudgetIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> PointwiseTailBudgetIssue:
    return PointwiseTailBudgetIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[PointwiseTailBudgetIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_weight_sums() -> list[PointwiseTailBudgetIssue]:
    k = sp.symbols("k", integer=True, positive=True)
    b_weights = [-1, 2, -1]
    companion_weights = [-1, 3, -3, 1]
    weighted_gap_weights = [2 * k + 1, -(6 * k + 5), 6 * k + 7, -(2 * k + 3)]
    issues: list[PointwiseTailBudgetIssue] = []
    if sum(abs(value) for value in b_weights) != 4:
        issues.append(issue("symbolic", "bad-B-weight-sum", repr(b_weights)))
    if sum(abs(value) for value in companion_weights) != 8:
        issues.append(issue("symbolic", "bad-companion-weight-sum", repr(companion_weights)))
    if sp.simplify(sum(abs(value) for value in weighted_gap_weights) - (16 * k + 16)) != 0:
        issues.append(issue("symbolic", "bad-weighted-gap-weight-sum", repr(weighted_gap_weights)))
    if (16 * 22 + 16) != 368:
        issues.append(issue("symbolic", "bad-k22-weighted-sum", str(16 * 22 + 16)))
    return issues


def validate_top_level(artifact: dict) -> list[PointwiseTailBudgetIssue]:
    issues: list[PointwiseTailBudgetIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_stencil_remainder_obligations",
        "source_uniform_remainder_target",
        "source_taylor_moment_budget",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "finite", "does not prove", "analytic tail", "scaled-curvature", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[PointwiseTailBudgetIssue]:
    params = artifact.get("matrix_rows", [{}])[3].get("diagnostics", {}).get("parameters", {})
    try:
        recomputed = build_artifact(
            int(params.get("max_taylor_degree", 14)),
            int(params.get("tail_cutoff_n", 80)),
            int(params.get("precision_bits", 256)),
            int(params.get("tail_start_k", 22)),
            [int(value) for value in params.get("sample_T_values", [25, 50, 100, 200, 500, 1000, 2000])],
        )
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[PointwiseTailBudgetIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_budget_rows(artifact: dict) -> list[PointwiseTailBudgetIssue]:
    diagnostics = artifact.get("matrix_rows", [{}])[3].get("diagnostics", {})
    budget_rows = diagnostics.get("budget_rows", [])
    issues: list[PointwiseTailBudgetIssue] = []
    if len(budget_rows) != 4:
        issues.append(issue("budget_rows", "bad-budget-row-count", str(len(budget_rows))))
    limiting_counts: dict[str, int] = {}
    weakest_half: tuple[Decimal, str, str] | None = None
    for row in budget_rows:
        row_id = str(row.get("source_row", "<missing-source>"))
        try:
            k = int(row["k"])
            b_margin = Decimal(row["B_margin"])
            companion_margin = Decimal(row["companion_margin"])
            weighted_gap_margin = Decimal(row["weighted_gap_margin"])
            b_eta = b_margin / Decimal(4)
            companion_eta = companion_margin / Decimal(8)
            weighted_eta = weighted_gap_margin / Decimal(16 * k + 16)
            limits = {"B": b_eta, "companion": companion_eta, "weighted_gap": weighted_eta}
            limiting, uniform_eta = min(limits.items(), key=lambda item: item[1])
            half_eta = uniform_eta / Decimal(2)
            rho = relative_tail_ratio_from_log_bound(half_eta)
        except Exception as exc:
            issues.append(issue(row_id, "bad-budget-row", f"{type(exc).__name__}: {exc}"))
            continue
        expected = {
            "B_weight_sum": decimal_format(Decimal(4)),
            "companion_weight_sum": decimal_format(Decimal(8)),
            "weighted_gap_weight_sum": decimal_format(Decimal(16 * k + 16)),
            "B_eta_bound": decimal_format(b_eta),
            "companion_eta_bound": decimal_format(companion_eta),
            "weighted_gap_eta_bound": decimal_format(weighted_eta),
            "uniform_eta_bound": decimal_format(uniform_eta),
            "half_safety_eta_bound": decimal_format(half_eta),
            "half_safety_relative_tail_ratio_bound": decimal_format(rho),
        }
        for key, value in expected.items():
            if row.get(key) != value:
                issues.append(issue(row_id, f"bad-{key}", f"{row.get(key)!r} != {value!r}"))
        if row.get("limiting_stencil") != limiting:
            issues.append(issue(row_id, "bad-limiting-stencil", f"{row.get('limiting_stencil')!r} != {limiting!r}"))
        limiting_counts[limiting] = limiting_counts.get(limiting, 0) + 1
        if weakest_half is None or half_eta < weakest_half[0]:
            weakest_half = (half_eta, decimal_format(half_eta), row_id)
    if limiting_counts != {"companion": 3, "weighted_gap": 1}:
        issues.append(issue("budget_rows", "bad-limiting-counts", repr(limiting_counts)))
    if weakest_half != (Decimal("1.4581132055269780525E-9"), "1.458113205526978052E-9", "nlrgts_M6_T2000"):
        issues.append(issue("budget_rows", "bad-weakest-half", repr(weakest_half)))
    return issues


def validate_rows(artifact: dict) -> tuple[list[PointwiseTailBudgetIssue], int, int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[PointwiseTailBudgetIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_sufficient = 0
    ready_to_apply = 0
    live_routes = 0
    rejected_routes = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("role") == "exact_sufficient_condition":
            exact_sufficient += 1
        if row.get("role") == "live_route":
            live_routes += 1
        if row.get("role") == "rejected_route":
            rejected_routes += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") not in {"available_exact", "not_ready_to_apply"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "live", "hygiene", "rejected")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
        if row.get("role") in {"live_route", "rejected_route"} and "prove" in str(row.get("proof_boundary", "")).lower() and "not" not in str(row.get("proof_boundary", "")).lower():
            issues.append(issue(row_id, "promotional-boundary", str(row.get("proof_boundary", ""))))
    return issues, len(rows), exact_sufficient, live_routes, rejected_routes if ready_to_apply == 0 else -rejected_routes


def validate_summary(
    artifact: dict,
    row_count: int,
    exact_sufficient: int,
    live_routes: int,
    rejected_routes: int,
) -> list[PointwiseTailBudgetIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 9,
        "positive_baseline_rows": 4,
        "blocked_baseline_rows": 31,
        "budget_rows": 4,
        "exact_sufficient_rows": 2,
        "weakest_half_safety_eta_bound": "1.458113205526978052E-9",
        "weakest_half_safety_relative_tail_ratio_bound": "1.458113204463930993E-9",
        "limiting_stencil_counts": {"companion": 3, "weighted_gap": 1},
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[PointwiseTailBudgetIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_sufficient != 2:
        issues.append(issue("summary", "bad-exact-sufficient-count", str(exact_sufficient)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    if rejected_routes != 1:
        issues.append(issue("summary", "bad-rejected-route-count", str(rejected_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("pointwise log-tail target", "4, 8, and 368", "1.458113205526978052e-9", "1.458113204463930993e-9", "does not prove"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "required tolerances", "not promoted", "companion", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[PointwiseTailBudgetIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[PointwiseTailBudgetIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "analytic tail estimates are proved",
        "uniform taylor-tail remainder theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[PointwiseTailBudgetIssue], dict]:
    artifact = load_json(target_path)
    issues: list[PointwiseTailBudgetIssue] = []
    issues.extend(validate_symbolic_weight_sums())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_budget_rows(artifact))
    row_issues, row_count, exact_sufficient, live_routes, rejected_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_sufficient, live_routes, rejected_routes))
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-POINTWISE-TAIL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('positive_baseline_rows')} positive baseline rows, "
            f"{summary.get('budget_rows')} budget rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
