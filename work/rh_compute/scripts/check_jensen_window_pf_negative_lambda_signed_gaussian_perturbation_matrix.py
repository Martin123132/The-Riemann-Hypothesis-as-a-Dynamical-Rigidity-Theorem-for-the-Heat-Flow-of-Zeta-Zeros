#!/usr/bin/env python3
"""Validate the negative-lambda signed Gaussian perturbation matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import json
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_SIGN_SCOUT_JSON,
    REPO_ROOT,
    build_diagnostics,
)


REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Signed Gaussian Perturbation Matrix",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.py",
    "validated Jensen-window PF negative-lambda signed Gaussian perturbation matrix: 8 matrix rows, 2 certified Taylor signs, 1 fixed-k activation estimates, 0 ready-to-apply rows, 0 issues",
    "log x_k = (2*b-a^2)/T^2 + O_k(T^-3)",
    "B_k = -log(x_k) = (a^2-2*b)/T^2 + O_k(T^-3)",
    "T >=",
    "ignores remainders",
    "outputs/jensen_window_pf_phi_taylor_cone_entry_sign_scout.md",
    "outputs/jensen_window_pf_negative_lambda_gaussian_curvature_matrix.md",
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
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_phi_taylor_sign_scout",
        "source_gaussian_curvature_matrix",
        "source_bounded_log_curvature_target",
        "sign_scout_json",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite theorem-search", "does not prove", "uniform", "bounded log-curvature", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    ref = artifact.get("sign_scout_json")
    if not isinstance(ref, str):
        return [issue("finite_diagnostics", "missing-sign-scout-ref", repr(ref))]
    try:
        diagnostics = asdict(build_diagnostics(REPO_ROOT / ref, artifact.get("finite_diagnostics", {}).get("fixed_k_activation_depth", 21)))
    except Exception as exc:
        return [issue("finite_diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    recorded = artifact.get("finite_diagnostics", {})
    for key, value in diagnostics.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))
    if recorded.get("certified_taylor_signs") != 2:
        issues.append(issue("finite_diagnostics", "bad-certified-signs", repr(recorded.get("certified_taylor_signs"))))
    if recorded.get("fixed_k_activation_depth") != 21:
        issues.append(issue("finite_diagnostics", "bad-activation-depth", repr(recorded.get("fixed_k_activation_depth"))))
    return issues


def validate_rows_and_summary(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    rows = artifact.get("matrix_rows", [])
    if not isinstance(rows, list) or len(rows) != 8:
        issues.append(issue("matrix_rows", "bad-row-count", repr(len(rows) if isinstance(rows, list) else type(rows))))
    roles = {row.get("role") for row in rows if isinstance(row, dict)}
    for required in ("exact_baseline", "formal_model", "certified_local_sign", "finite_depth_estimate", "open_gap", "rejected_template", "live_route"):
        if required not in roles:
            issues.append(issue("matrix_rows", "missing-role", required))
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            issues.append(issue("matrix_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("not", "only", "open", "reject", "live")):
            issues.append(issue(row_id, "weak-boundary", boundary))

    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 8,
        "ready_to_apply_rows": 0,
        "certified_taylor_signs": 2,
        "fixed_k_activation_estimates": 1,
        "live_routes": 2,
        "rejected_templates": 1,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("signed", "positive leading deficit", "monotone correction", "uniform remainder", "not a positive gaussian mixture"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("no row", "fixed-k", "positive gaussian", "bounded log-curvature", "lambda <= 0"):
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
        "uniform remainder is proved",
        "cone entry is proved",
        "jwpf_06 is proved",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-text", forbidden))
    return issues


def validate(target_path: Path, note_path: Path) -> tuple[list[dict], dict]:
    artifact = load_json(target_path)
    issues: list[dict] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    issues.extend(validate_rows_and_summary(artifact))
    issues.extend(validate_note(note_path))
    return issues, artifact.get("summary", {})


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--sign-scout-json", type=Path, default=DEFAULT_SIGN_SCOUT_JSON)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    issues, summary = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": issues}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-SIGNED-GAUSSIAN {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda signed Gaussian perturbation matrix: "
            f"{summary.get('matrix_rows')} matrix rows, "
            f"{summary.get('certified_taylor_signs')} certified Taylor signs, "
            f"{summary.get('fixed_k_activation_estimates')} fixed-k activation estimates, "
            f"{summary.get('ready_to_apply_rows')} ready-to-apply rows, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
