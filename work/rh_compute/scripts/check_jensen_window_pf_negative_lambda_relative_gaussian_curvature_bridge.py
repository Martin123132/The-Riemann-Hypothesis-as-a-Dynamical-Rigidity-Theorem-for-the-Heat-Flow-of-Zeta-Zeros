#!/usr/bin/env python3
"""Validate the relative-Gaussian curvature bridge."""

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

from jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "nlrgcb_01_relative_gaussian_log_coordinate",
    "nlrgcb_02_deficit_as_negative_second_difference",
    "nlrgcb_03_scaled_curvature_weighted_third_gap",
    "nlrgcb_04_companion_upper_side_third_difference",
    "nlrgcb_05_repaired_k300_curvature_ladder",
    "nlrgcb_06_lambda_order_finite_pattern",
    "nlrgcb_07_live_tail_profile_route",
    "nlrgcb_08_acceptance_gate",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "finite_stress",
    "finite_order_pattern",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Relative-Gaussian Curvature Bridge",
    "Status: exact finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.py",
    "validated Jensen-window PF negative-lambda relative-Gaussian curvature bridge: 8 rows, 0 issues, 897 B-positive rows, 894 B-decrease rows, 894 C-increase rows, 598 C-lambda-order rows, 0 ready-to-apply rows",
    "B_k = h_k-log(R_k) = 2*f_k-f_(k-1)-f_(k+1)",
    "C_(k+1)-C_k",
    "(2*k+1)*f_(k-1)-(6*k+5)*f_k+(6*k+7)*f_(k+1)-(2*k+3)*f_(k+2)",
    "B_k-B_(k+1) = -f_(k-1)+3*f_k-3*f_(k+1)+f_(k+2) >= 0",
    "B-positive rows: 897 / 897",
    "B-decrease rows: 894 / 894",
    "C-increase rows: 894 / 894",
    "C lambda-order rows: 598 / 598",
    "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md",
    "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
    "outputs/jensen_window_pf_negative_lambda_scaled_curvature_log_ceiling_bridge.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
    "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
)


@dataclass(frozen=True)
class RelativeGaussianIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> RelativeGaussianIssue:
    return RelativeGaussianIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[RelativeGaussianIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_symbolic_identities() -> list[RelativeGaussianIssue]:
    k = sp.symbols("k")
    fm1, f0, fp1, fp2 = sp.symbols("fm1 f0 fp1 fp2")
    b0 = 2 * f0 - fm1 - fp1
    b1 = 2 * fp1 - f0 - fp2
    c_gap = (2 * k + 3) * b1 - (2 * k + 1) * b0
    weighted_gap = (2 * k + 1) * fm1 - (6 * k + 5) * f0 + (6 * k + 7) * fp1 - (2 * k + 3) * fp2
    third_side = -fm1 + 3 * f0 - 3 * fp1 + fp2
    issues: list[RelativeGaussianIssue] = []
    if sp.simplify(c_gap - weighted_gap) != 0:
        issues.append(issue("symbolic", "bad-weighted-gap", str(sp.simplify(c_gap - weighted_gap))))
    if sp.simplify((b0 - b1) - third_side) != 0:
        issues.append(issue("symbolic", "bad-third-difference", str(sp.simplify((b0 - b1) - third_side))))
    if sp.simplify(b0 + (fp1 + fm1 - 2 * f0)) != 0:
        issues.append(issue("symbolic", "bad-negative-second-difference", str(sp.simplify(b0 + (fp1 + fm1 - 2 * f0)))))
    return issues


def validate_top_level(artifact: dict) -> list[RelativeGaussianIssue]:
    issues: list[RelativeGaussianIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "exact finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_gaussian_curvature_matrix",
        "source_scaled_curvature_target",
        "source_scaled_log_ceiling_bridge",
        "source_monotone_contraction_target",
        "source_raw_corridor_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("exact", "finite", "does not prove", "weighted four-point", "monotone-contraction", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[RelativeGaussianIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_artifact([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[RelativeGaussianIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[RelativeGaussianIssue], int, int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[RelativeGaussianIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    rows_by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(rows_by_id)):
        issues.append(issue(missing, "missing-row", missing))
    exact_available = 0
    live_routes = 0
    ready_to_apply = 0
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
        readiness = row.get("readiness")
        if row.get("role") == "exact_reduction":
            exact_available += 1
            if readiness != "available_exact":
                issues.append(issue(row_id, "bad-exact-readiness", repr(readiness)))
        elif readiness != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-open-readiness", repr(readiness)))
        if readiness == "ready_to_apply":
            ready_to_apply += 1
        if row.get("role") == "live_route":
            live_routes += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("exact", "finite", "not", "only", "open", "live", "hygiene")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), exact_available, live_routes if ready_to_apply == 0 else -live_routes


def validate_summary(artifact: dict, row_count: int, exact_available: int, live_routes: int) -> list[RelativeGaussianIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "exact_available_rows": 4,
        "finite_stress_rows": 1,
        "finite_order_rows": 1,
        "live_routes": 1,
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "coefficient_cap": 300,
        "raw_log_total_rows": 897,
        "adjacent_total_rows": 894,
        "lambda_order_total_rows": 598,
        "b_positive_rows": 897,
        "b_decrease_rows": 894,
        "c_increase_rows": 894,
        "c_lambda_order_rows": 598,
    }
    issues: list[RelativeGaussianIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 8:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if exact_available != 4:
        issues.append(issue("summary", "bad-exact-count", str(exact_available)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("weighted four-point", "relative to the gaussian", "897/897", "894/894", "598/598"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "weighted four-point", "finite stress", "monotone-contraction", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[RelativeGaussianIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[RelativeGaussianIssue] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "weighted four-point inequality is proved",
        "scaled-curvature monotonicity is proved",
        "raw-corridor theorem is proved",
        "cone entry is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[RelativeGaussianIssue], dict]:
    artifact = load_json(target_path)
    issues: list[RelativeGaussianIssue] = []
    issues.extend(validate_symbolic_identities())
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, exact_available, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, exact_available, live_routes))
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
            print(f"JWPF-NEG-LAMBDA-REL-GAUSS-CURV {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda relative-Gaussian curvature bridge: "
            f"{summary.get('matrix_rows')} rows, {len(issues)} issues, "
            f"{summary.get('b_positive_rows')} B-positive rows, "
            f"{summary.get('b_decrease_rows')} B-decrease rows, "
            f"{summary.get('c_increase_rows')} C-increase rows, "
            f"{summary.get('c_lambda_order_rows')} C-lambda-order rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
