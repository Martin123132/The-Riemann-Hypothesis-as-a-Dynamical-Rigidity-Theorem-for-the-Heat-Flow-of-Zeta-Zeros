#!/usr/bin/env python3
"""Validate the negative-lambda Taylor-moment budget scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_taylor_moment_budget import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_diagnostics,
)


REQUIRED_ROW_IDS = {
    "nltmb_01_gaussian_moment_factorization",
    "nltmb_02_truncated_multiplier_formula",
    "nltmb_03_tail_start_stability_samples",
    "nltmb_04_first_correction_budget",
    "nltmb_05_normalizer_positivity_requirement",
    "nltmb_06_log_curvature_remainder_requirement",
    "nltmb_07_monotone_remainder_requirement",
    "nltmb_08_low_order_finite_T_promotion_rejected",
    "nltmb_09_local_mesoscopic_theorem_handoff",
}

ALLOWED_ROLES = {
    "exact_reduction",
    "formal_model",
    "finite_diagnostic",
    "scale_budget",
    "open_requirement",
    "rejected_route",
    "conditional_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Taylor Moment Budget",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_taylor_moment_budget`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_taylor_moment_budget.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_taylor_moment_budget.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_taylor_moment_budget.py",
    "validated Jensen-window PF negative-lambda Taylor moment budget: 9 budget rows, 7 tail-start samples, 4 invalid truncation rows, 2 bounded truncation rows, 0 ready-to-apply rows, 0 issues",
    "q = k+1/2",
    "F_k = 1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3+R_k^(>=8)",
    "x_k = F_(k+1)*F_(k-1)/F_k^2",
    "k=22 enters only after T >=",
    "low-order Taylor truncation is not a finite proof",
    "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
    "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
    "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
    "outputs/jensen_window_pf_negative_lambda_bounded_log_curvature_target.md",
)


def issue(section: str, name: str, detail: str) -> dict:
    return {"section": section, "issue": name, "detail": detail}


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[dict]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_taylor_moment_budget":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite_theorem_search_diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_uniform_remainder_target",
        "source_signed_gaussian_matrix",
        "source_phi_taylor_sign_scout",
        "source_bounded_log_curvature_target",
        "sign_scout_json",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite theorem-search diagnostic",
        "does not prove a taylor remainder",
        "does not prove bounded log-curvature",
        "does not prove cone entry",
        "lambda <= 0",
    ):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    formulae = " ".join(str(value) for value in artifact.get("formulae", {}).values()).lower()
    for required in ("q=k+1/2", "r_k^(>=8)", "x_k", "b_k=-log"):
        if required not in formulae:
            issues.append(issue("formulae", "missing-formula-text", required))
    return issues


def validate_recomputed(artifact: dict) -> list[dict]:
    ref = artifact.get("sign_scout_json")
    if not isinstance(ref, str):
        return [issue("diagnostics", "missing-sign-scout-ref", repr(ref))]
    try:
        diagnostics = asdict(build_diagnostics(REPO_ROOT / ref, 22, [25, 50, 100, 200, 500, 1000, 2000]))
    except Exception as exc:
        return [issue("diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    recorded = artifact.get("diagnostics", {})
    issues: list[dict] = []
    for key, value in diagnostics.items():
        if recorded.get(key) != value:
            issues.append(issue("diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))
    statuses = [row.get("truncation_status") for row in recorded.get("tail_start_samples", [])]
    expected_statuses = [
        "invalid_truncation_normalizer",
        "invalid_truncation_normalizer",
        "invalid_truncation_normalizer",
        "invalid_truncation_normalizer",
        "overbound_positive_truncation",
        "bounded_positive_truncation",
        "bounded_positive_truncation",
    ]
    if statuses != expected_statuses:
        issues.append(issue("diagnostics", "bad-sample-statuses", repr(statuses)))
    budget = recorded.get("first_correction_budget", {})
    if budget.get("threshold") != "abs(a)*(k+1/2)/T <= 1/2":
        issues.append(issue("first_correction_budget", "bad-threshold", repr(budget.get("threshold"))))
    return issues


def validate_rows(artifact: dict) -> tuple[list[dict], int, int, int]:
    rows = artifact.get("budget_rows", [])
    issues: list[dict] = []
    if not isinstance(rows, list):
        return [issue("budget_rows", "bad-rows", repr(type(rows)))], 0, 0, 0
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        issues.append(issue("budget_rows", "bad-row-ids", repr(sorted(ids))))
    ready_count = 0
    open_requirement_count = 0
    for row in rows:
        if not isinstance(row, dict):
            issues.append(issue("budget_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "readiness", "source_artifacts", "claim", "acceptance_test", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        if row.get("role") not in ALLOWED_ROLES:
            issues.append(issue(row_id, "bad-role", repr(row.get("role"))))
        if row.get("readiness") == "ready_to_apply":
            ready_count += 1
            issues.append(issue(row_id, "ready-to-apply-overclaim", row_id))
        elif row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("role") == "open_requirement":
            open_requirement_count += 1
        for ref in row.get("source_artifacts", []):
            issues.extend(validate_ref(row_id, ref))
        text = " ".join(str(row.get(key, "")) for key in ("claim", "acceptance_test", "proof_boundary")).lower()
        if row.get("role") in {"open_requirement", "conditional_handoff"} and "not" not in text:
            issues.append(issue(row_id, "weak-boundary", text))
        if row.get("role") == "rejected_route" and "reject" not in text:
            issues.append(issue(row_id, "missing-rejection-language", text))
    return issues, len(rows), ready_count, open_requirement_count


def validate_summary(artifact: dict, row_count: int, ready_count: int, open_requirement_count: int) -> list[dict]:
    summary = artifact.get("summary", {})
    issues: list[dict] = []
    expected = {
        "budget_rows": 9,
        "tail_start_samples": 7,
        "invalid_truncation_rows": 4,
        "positive_truncation_rows": 3,
        "overbound_truncation_rows": 1,
        "bounded_truncation_rows": 2,
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 9:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if ready_count != 0:
        issues.append(issue("summary", "ready-row-present", str(ready_count)))
    if open_requirement_count != 3:
        issues.append(issue("summary", "bad-open-requirement-count", str(open_requirement_count)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("gaussian-moment", "q/t", "low-order taylor", "actual finite prefix", "remainder"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("no row", "not the actual zeta", "r_k^(>=8)", "q/t constants", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[dict]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[dict] = []
    for required in REQUIRED_NOTE_STRINGS:
        if required not in text:
            issues.append(issue("note", "missing-text", required))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "lambda <= 0 is proved",
        "bounded log-curvature theorem is proved",
        "cone entry is proved",
        "taylor-tail remainder theorem is proved",
        "low-order taylor truncation proves",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(artifact_path: Path, note_path: Path) -> tuple[list[dict], dict]:
    artifact = load_json(artifact_path)
    issues: list[dict] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, ready_count, open_requirement_count = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, ready_count, open_requirement_count))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, summary = validate(args.artifact, args.note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": issues}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-TAYLOR-MOMENT {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda Taylor moment budget: "
            f"{summary.get('budget_rows')} budget rows, "
            f"{summary.get('tail_start_samples')} tail-start samples, "
            f"{summary.get('invalid_truncation_rows')} invalid truncation rows, "
            f"{summary.get('bounded_truncation_rows')} bounded truncation rows, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows, {len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
