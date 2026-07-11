#!/usr/bin/env python3
"""Validate the negative-lambda defect-recurrence scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_defect_recurrence_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_diagnostics,
)


REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Defect-Recurrence Scout",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_defect_recurrence_scout`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_defect_recurrence_scout.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_defect_recurrence_scout.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_defect_recurrence_scout.py",
    "validated Jensen-window PF negative-lambda defect-recurrence scout: 63 buffered rows, 60 defect-monotone rows, 60 width-recurrence rejections, 1 live sufficient routes, 0 issues",
    "0 <= d_k <= 2/(3*(2*k+1))",
    "d_(k+1) <= d_k",
    "width-preserving rejected rows: 60 / 60",
    "outputs/jensen_window_pf_negative_lambda_tail_barrier_scout.md",
    "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
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
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_defect_recurrence_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in ("source_tail_barrier_scout", "source_defect_tail_target", "generator", "checker"):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite theorem-search", "does not prove", "all-k tail", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    paths = [REPO_ROOT / ref for ref in artifact.get("enclosure_jsonl", [])]
    try:
        diagnostics = asdict(build_diagnostics(paths))
    except Exception as exc:
        return [issue("finite_diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    recorded = artifact.get("finite_diagnostics", {})
    for key, value in diagnostics.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))
    expected = {
        "coefficient_k_max": 22,
        "checked_x_max": 21,
        "buffered_sufficient_rows": 63,
        "buffered_sufficient_positive_rows": 63,
        "defect_monotone_rows": 60,
        "defect_monotone_positive_rows": 60,
        "width_preserving_rows": 60,
        "width_preserving_positive_rows": 0,
        "width_preserving_rejected_rows": 60,
    }
    for key, value in expected.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"unexpected-{key}", f"{recorded.get(key)!r} != {value!r}"))
    return issues


def validate_rows_and_summary(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    rows = artifact.get("recurrence_rows", [])
    if not isinstance(rows, list) or len(rows) != 4:
        issues.append(issue("recurrence_rows", "bad-row-count", repr(len(rows) if isinstance(rows, list) else type(rows))))
    roles = {row.get("role") for row in rows if isinstance(row, dict)}
    for required in ("live_sufficient_condition", "insufficient_condition", "rejected_candidate"):
        if required not in roles:
            issues.append(issue("recurrence_rows", "missing-role", required))
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            issues.append(issue("recurrence_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("not", "only", "insufficient", "rejected")):
            issues.append(issue(row_id, "weak-boundary", boundary))

    summary = artifact.get("summary", {})
    expected = {
        "buffered_sufficient_rows": 63,
        "defect_monotone_rows": 60,
        "width_preserving_rows": 60,
        "width_preserving_rejected_rows": 60,
        "ready_to_apply_rows": 0,
        "live_sufficient_routes": 1,
        "rejected_routes": 1,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("buffered sufficient", "0<=d_k", "k>=22", "k>=21", "width-preserving", "rejected"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "buffered", "width-preserving", "scaled-defect", "lambda <= 0"):
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
        "defect-tail theorem is proved",
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
            print(f"JWPF-NEG-LAMBDA-DEFECT-RECURRENCE {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda defect-recurrence scout: "
            f"{summary.get('buffered_sufficient_rows')} buffered rows, "
            f"{summary.get('defect_monotone_rows')} defect-monotone rows, "
            f"{summary.get('width_preserving_rejected_rows')} width-recurrence rejections, "
            f"{summary.get('live_sufficient_routes')} live sufficient routes, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
