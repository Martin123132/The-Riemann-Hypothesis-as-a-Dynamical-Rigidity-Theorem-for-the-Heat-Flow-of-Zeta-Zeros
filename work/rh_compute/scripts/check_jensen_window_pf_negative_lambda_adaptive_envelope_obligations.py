#!/usr/bin/env python3
"""Validate exact obligations for the negative-lambda adaptive envelope route."""

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

from jensen_window_pf_negative_lambda_adaptive_envelope_obligations import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlaeo_01_scaled_defect_coordinates",
    "nlaeo_02_lower_wall_equivalence",
    "nlaeo_03_fixed_buffer_thresholds",
    "nlaeo_04_nonnegative_defect_upper_wall",
    "nlaeo_05_monotone_defect_bridge",
    "nlaeo_06_scaled_k_monotone_identity",
    "nlaeo_07_monotone_envelope_implication",
    "nlaeo_08_lambda_order_requirement",
    "nlaeo_09_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reindexing",
    "exact_reduction",
    "rejected_route",
    "open_requirement",
    "conditional_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Adaptive Envelope Obligations",
    "Status: exact algebraic obligation diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_adaptive_envelope_obligations`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_obligations.py",
    "validated Jensen-window PF negative-lambda adaptive envelope obligations: 9 obligation rows, 0 issues, 3 exact rows, 3 open requirements, 1 rejected routes",
    "x_k=(A_{k-1}*A_{k+1})/A_k^2",
    "d_k=1-x_k",
    "s_k=((2*k+1)/2)*d_k",
    "1-s_k=((2*k+1)*x_k-(2*k-1))/2",
    "s_k>=0 iff x_k<=1",
    "d_k-d_(k+1)=x_(k+1)-x_k",
    "s_(k+1)-s_k=(2+(2*k+1)*x_k-(2*k+3)*x_(k+1))/2",
    "already exact: x_k >= (2*k-1)/(2*k+1) from the boundary-threshold lemma",
    "still open: x_k <= 1 on the needed tail",
    "still open: x_(k+1) >= x_k on the needed tail",
    "still open: limiting/adaptive envelope E_lambda<1",
    "rejected: fixed half-width and one-third buffers as global routes",
    "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
    "outputs/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.md",
    "outputs/jensen_window_pf_heat_flow_boundary_threshold_lemma.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
    "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md",
    "outputs/jensen_window_pf_negative_lambda_half_width_tail_target.md",
)


@dataclass(frozen=True)
class ObligationIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> ObligationIssue:
    return ObligationIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[ObligationIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[ObligationIssue]:
    k, x, y = sp.symbols("k x y")
    s = (2 * k + 1) * (1 - x) / 2
    s_next = (2 * k + 3) * (1 - y) / 2
    d = 1 - x
    d_next = 1 - y
    checks = {
        "lower-wall": sp.simplify((1 - s) - (((2 * k + 1) * x - (2 * k - 1)) / 2)),
        "half-width": sp.simplify((sp.Rational(1, 2) - s) - (((2 * k + 1) * x - 2 * k) / 2)),
        "one-third": sp.simplify((sp.Rational(1, 3) - s) - ((3 * (2 * k + 1) * x - (6 * k + 1)) / 6)),
        "defect-monotone": sp.simplify((d - d_next) - (y - x)),
        "scaled-k-increase": sp.simplify((s_next - s) - ((2 + (2 * k + 1) * x - (2 * k + 3) * y) / 2)),
    }
    issues: list[ObligationIssue] = []
    for name, value in checks.items():
        if sp.simplify(value) != 0:
            issues.append(issue("symbolic", f"bad-{name}", str(value)))
    half_margin = sp.Rational(1, 2) - s
    half_threshold_margin = sp.simplify(half_margin.subs(x, 2 * k / (2 * k + 1)))
    if half_threshold_margin != 0:
        issues.append(issue("symbolic", "bad-half-threshold", str(half_threshold_margin)))
    one_third_margin = sp.Rational(1, 3) - s
    one_third_threshold_margin = sp.simplify(one_third_margin.subs(x, 1 - sp.Rational(2, 1) / (3 * (2 * k + 1))))
    if one_third_threshold_margin != 0:
        issues.append(issue("symbolic", "bad-one-third-threshold-margin", str(one_third_threshold_margin)))
    one_third_threshold = sp.simplify(
        (6 * k + 1) / (6 * k + 3) - (1 - sp.Rational(2, 1) / (3 * (2 * k + 1)))
    )
    if one_third_threshold != 0:
        issues.append(issue("symbolic", "bad-one-third-threshold", str(one_third_threshold)))
    return issues


def validate_top_level(artifact: dict) -> list[ObligationIssue]:
    issues: list[ObligationIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_adaptive_envelope_obligations":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact algebraic obligation diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_adaptive_target",
        "source_envelope_matrix",
        "source_boundary_threshold",
        "source_monotone_contraction_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "open", "does not prove", "adjacent log-concavity", "monotone contractions", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[ObligationIssue]:
    recomputed = build_artifact()
    issues: list[ObligationIssue] = []
    for key in ("obligation_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[ObligationIssue], int]:
    rows = artifact.get("obligation_rows", [])
    issues: list[ObligationIssue] = []
    if not isinstance(rows, list):
        return [issue("obligation_rows", "bad-rows", repr(type(rows)))], 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("obligation_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "claim", "formula", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") not in {"available_exact", "not_ready_to_apply"}:
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("role") in {"open_requirement", "conditional_route", "rejected_route", "acceptance_gate"} and row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-open-readiness", repr(row.get("readiness"))))
        if row.get("role") in {"exact_reindexing", "exact_reduction"} and row.get("readiness") != "available_exact":
            issues.append(issue(row_id, "bad-exact-readiness", repr(row.get("readiness"))))
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "not", "open", "only", "finite")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows)


def validate_summary(artifact: dict, row_count: int) -> list[ObligationIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "obligation_rows": 9,
        "exact_available_rows": 3,
        "open_requirement_rows": 3,
        "conditional_route_rows": 1,
        "rejected_route_rows": 1,
        "acceptance_gate_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    issues: list[ObligationIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("decomposes exactly", "boundary-threshold", "x_k<=1", "x_(k+1)>=x_k", "scaled k-monotonicity", "halF-width".lower(), "one-third"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("no row", "lower threshold", "upper wall", "monotone bridge", "half-width", "one-third", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[ObligationIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[ObligationIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "adaptive scaled-defect target is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[ObligationIssue], dict]:
    artifact = load_json(target_path)
    issues: list[ObligationIssue] = []
    issues.extend(validate_symbolic_identities())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count))
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
            print(f"JWPF-NEG-LAMBDA-ADAPTIVE-ENVELOPE-OBLIGATION {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda adaptive envelope obligations: "
            f"{summary.get('obligation_rows')} obligation rows, {len(issues)} issues, "
            f"{summary.get('exact_available_rows')} exact rows, "
            f"{summary.get('open_requirement_rows')} open requirements, "
            f"{summary.get('rejected_route_rows')} rejected routes"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
