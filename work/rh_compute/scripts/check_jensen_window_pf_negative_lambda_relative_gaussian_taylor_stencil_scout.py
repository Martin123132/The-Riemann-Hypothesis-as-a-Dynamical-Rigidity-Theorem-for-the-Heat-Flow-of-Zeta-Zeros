#!/usr/bin/env python3
"""Validate the relative-Gaussian Taylor stencil scout."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgts_01_relative_gaussian_taylor_multiplier",
    "nlrgts_02_B_leading_sign",
    "nlrgts_03_companion_leading_sign",
    "nlrgts_04_weighted_gap_leading_sign",
    "nlrgts_05_finite_truncation_stencil_matrix",
    "nlrgts_06_finite_truncation_promotion_rejected",
    "nlrgts_07_live_uniform_stencil_remainder_route",
    "nlrgts_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "formal_model",
    "certified_local_sign",
    "finite_diagnostic",
    "rejected_route",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Taylor Stencil Scout",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: 8 rows, 0 issues, 3 certified leading-sign rows, 35 truncation rows, 4 all-positive stencil rows, 0 ready-to-apply rows",
    "B_k = (a^2-2*b)/T^2 + O_k(T^-3)",
    "B_k-B_(k+1) = 2*(a^3-3*a*b+3*c)/T^3 + O_k(T^-4)",
    "C_(k+1)-C_k = 2*(a^2-2*b)/T^2 + O_k(T^-3)",
    "all-positive stencil rows: 4",
    "invalid normalizers: 9",
    "weighted-gap failure rows: 12",
    "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md",
    "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
    "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
)


@dataclass(frozen=True)
class TaylorStencilIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> TaylorStencilIssue:
    return TaylorStencilIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[TaylorStencilIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_expansions() -> list[TaylorStencilIssue]:
    k, e, a, b, c = sp.symbols("k e a b c")

    def f_at(offset: int) -> sp.Expr:
        q = k + sp.Rational(1, 2) + offset
        F = 1 + a * q * e + b * q * (q + 1) * e**2 + c * q * (q + 1) * (q + 2) * e**3
        return sp.series(sp.log(F), e, 0, 5).removeO().expand()

    f_m = f_at(-1)
    f_0 = f_at(0)
    f_p = f_at(1)
    f_pp = f_at(2)
    B = sp.expand(2 * f_0 - f_m - f_p)
    B_next = sp.expand(2 * f_p - f_0 - f_pp)
    companion = sp.expand(B - B_next)
    weighted = sp.expand((2 * k + 3) * B_next - (2 * k + 1) * B)
    checks = {
        "B-e2": (B.coeff(e, 2), a**2 - 2 * b),
        "companion-e3": (companion.coeff(e, 3), 2 * (a**3 - 3 * a * b + 3 * c)),
        "weighted-e2": (weighted.coeff(e, 2), 2 * (a**2 - 2 * b)),
    }
    issues: list[TaylorStencilIssue] = []
    for name, (actual, expected) in checks.items():
        if sp.simplify(actual - expected) != 0:
            issues.append(issue("symbolic", f"bad-{name}", str(sp.simplify(actual - expected))))
    return issues


def validate_top_level(artifact: dict) -> list[TaylorStencilIssue]:
    issues: list[TaylorStencilIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_relative_gaussian_bridge",
        "source_taylor_moment_budget",
        "source_high_order_taylor_scout",
        "source_uniform_remainder_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite", "does not prove", "uniform", "scaled-curvature", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[TaylorStencilIssue]:
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
    issues: list[TaylorStencilIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[TaylorStencilIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[TaylorStencilIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    leading_rows = 0
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
        if row.get("readiness") == "ready_to_apply":
            ready_to_apply += 1
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("role") == "certified_local_sign":
            leading_rows += 1
        if row.get("role") == "live_route":
            live_routes += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "fixed-k", "not", "only", "live", "hygiene", "rejected", "formal")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), leading_rows, live_routes if ready_to_apply == 0 else -live_routes


def validate_summary(artifact: dict, row_count: int, leading_rows: int, live_routes: int) -> list[TaylorStencilIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "certified_leading_sign_rows": 3,
        "truncation_rows": 35,
        "invalid_normalizer_rows": 9,
        "positive_normalizer_rows": 26,
        "b_positive_rows": 24,
        "b_decrease_rows": 18,
        "c_increase_rows": 12,
        "all_positive_rows": 4,
        "upper_wall_violation_rows": 2,
        "companion_failure_rows": 8,
        "weighted_gap_failure_rows": 12,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[TaylorStencilIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if leading_rows != 3:
        issues.append(issue("summary", "bad-leading-count", str(leading_rows)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("fixed-k", "positive leading signs", "4/35", "uniform taylor-tail remainder"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "fixed-k", "finite taylor", "weighted four-point", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[TaylorStencilIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[TaylorStencilIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "uniform taylor-tail remainder theorem is proved",
        "scaled-curvature monotonicity is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[TaylorStencilIssue], dict]:
    artifact = load_json(target_path)
    issues: list[TaylorStencilIssue] = []
    issues.extend(validate_symbolic_expansions())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, leading_rows, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, leading_rows, live_routes))
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-TAYLOR-STENCIL {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('certified_leading_sign_rows')} certified leading-sign rows, "
            f"{summary.get('truncation_rows')} truncation rows, "
            f"{summary.get('all_positive_rows')} all-positive stencil rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
