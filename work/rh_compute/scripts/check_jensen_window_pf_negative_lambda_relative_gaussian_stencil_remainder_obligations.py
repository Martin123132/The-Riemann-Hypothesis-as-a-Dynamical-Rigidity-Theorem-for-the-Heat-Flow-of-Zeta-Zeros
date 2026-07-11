#!/usr/bin/env python3
"""Validate the relative-Gaussian stencil remainder obligations."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgsro_01_log_tail_error_coordinate",
    "nlrgsro_02_B_error_stencil",
    "nlrgsro_03_companion_error_stencil",
    "nlrgsro_04_weighted_gap_error_stencil",
    "nlrgsro_05_positive_baseline_margin_budget",
    "nlrgsro_06_blocked_baseline_rows",
    "nlrgsro_07_sufficient_remainder_budget",
    "nlrgsro_08_live_uniform_stencil_bound_route",
    "nlrgsro_09_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_diagnostic",
    "finite_obstruction",
    "exact_sufficient_condition",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Stencil Remainder Obligations",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: 9 rows, 0 issues, 4 positive baseline rows, 31 blocked baseline rows, 4 exact stencil rows, 0 ready-to-apply rows",
    "E_B = 2*epsilon_k-epsilon_(k-1)-epsilon_(k+1)",
    "E_U = -epsilon_(k-1)+3*epsilon_k-3*epsilon_(k+1)+epsilon_(k+2)",
    "E_C = (2*k+1)*epsilon_(k-1)-(6*k+5)*epsilon_k+(6*k+7)*epsilon_(k+1)-(2*k+3)*epsilon_(k+2)",
    "positive baseline rows: 4",
    "blocked baseline rows: 31",
    "weakest half-margin budget: 1.166490564421582442E-8 at nlrgts_M6_T2000",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
)


@dataclass(frozen=True)
class StencilRemainderIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> StencilRemainderIssue:
    return StencilRemainderIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[StencilRemainderIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_stencils() -> list[StencilRemainderIssue]:
    k = sp.symbols("k")
    em, e0, ep, epp = sp.symbols("em e0 ep epp")
    B_error = 2 * e0 - em - ep
    next_B_error = 2 * ep - e0 - epp
    companion_error = B_error - next_B_error
    weighted_error = (2 * k + 3) * next_B_error - (2 * k + 1) * B_error
    expected_companion = -em + 3 * e0 - 3 * ep + epp
    expected_weighted = (2 * k + 1) * em - (6 * k + 5) * e0 + (6 * k + 7) * ep - (2 * k + 3) * epp
    issues: list[StencilRemainderIssue] = []
    if sp.simplify(companion_error - expected_companion) != 0:
        issues.append(issue("symbolic", "bad-companion-error", str(sp.simplify(companion_error - expected_companion))))
    if sp.simplify(weighted_error - expected_weighted) != 0:
        issues.append(issue("symbolic", "bad-weighted-error", str(sp.simplify(weighted_error - expected_weighted))))
    return issues


def validate_top_level(artifact: dict) -> list[StencilRemainderIssue]:
    issues: list[StencilRemainderIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_taylor_stencil_scout",
        "source_relative_gaussian_bridge",
        "source_uniform_remainder_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "finite", "does not prove", "uniform", "scaled-curvature", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[StencilRemainderIssue]:
    params = artifact.get("matrix_rows", [{}])[4].get("diagnostics", {}).get("parameters", {})
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
    issues: list[StencilRemainderIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[StencilRemainderIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[StencilRemainderIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_rows = 0
    ready_to_apply = 0
    live_routes = 0
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
        if row.get("role") in {"exact_reduction", "exact_sufficient_condition"}:
            exact_rows += 1
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") not in {"available_exact", "not_ready_to_apply"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("role") == "live_route":
            live_routes += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "live", "hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), exact_rows, live_routes if ready_to_apply == 0 else -live_routes


def validate_summary(artifact: dict, row_count: int, exact_rows: int, live_routes: int) -> list[StencilRemainderIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 9,
        "exact_stencil_rows": 4,
        "truncation_rows": 35,
        "positive_baseline_rows": 4,
        "blocked_baseline_rows": 31,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[StencilRemainderIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_rows != 5:
        issues.append(issue("summary", "bad-exact-count", str(exact_rows)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("three exact epsilon stencils", "4/35", "1.166490564421582442e-8", "uniform epsilon-stencil"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "finite positive baseline", "uniform epsilon-stencil", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[StencilRemainderIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[StencilRemainderIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "uniform epsilon-stencil remainder theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[StencilRemainderIssue], dict]:
    artifact = load_json(target_path)
    issues: list[StencilRemainderIssue] = []
    issues.extend(validate_symbolic_stencils())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, exact_rows, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_rows, live_routes))
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-STENCIL-REMAINDER {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('positive_baseline_rows')} positive baseline rows, "
            f"{summary.get('blocked_baseline_rows')} blocked baseline rows, "
            f"{summary.get('exact_stencil_rows')} exact stencil rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
