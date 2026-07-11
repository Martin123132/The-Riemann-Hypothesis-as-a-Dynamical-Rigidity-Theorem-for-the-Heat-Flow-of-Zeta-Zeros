#!/usr/bin/env python3
"""Validate the negative-lambda scaled-defect frontier scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_scaled_defect_frontier_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_diagnostics,
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


def repo_ref(path: Path) -> str:
    abs_path = path if path.is_absolute() else REPO_ROOT / path
    return abs_path.relative_to(REPO_ROOT).as_posix()


def validate_top_level(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_scaled_defect_frontier_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in ("source_tail_barrier_scout", "source_defect_tail_target", "generator", "checker"):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for key in ("result_json", "note"):
        if key in artifact:
            issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite", "exact cone", "half-width", "one-third", "does not prove", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    paths = [REPO_ROOT / ref for ref in artifact.get("enclosure_jsonl", [])]
    try:
        recomputed = asdict(build_diagnostics(paths))
    except Exception as exc:
        return [issue("finite_diagnostics", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    recorded = artifact.get("finite_diagnostics", {})
    for key, value in recomputed.items():
        if recorded.get(key) != value:
            issues.append(issue("finite_diagnostics", f"bad-{key}", f"{recorded.get(key)!r} != {value!r}"))

    rows = int(recorded.get("scaled_defect_rows", 0))
    if recorded.get("cone_positive_rows") != rows:
        issues.append(issue("finite_diagnostics", "cone-not-all-positive", repr(recorded.get("cone_positive_rows"))))
    if not (0 < int(recorded.get("one_third_failure_rows", 0)) < rows):
        issues.append(issue("finite_diagnostics", "bad-one-third-frontier", repr(recorded.get("one_third_failure_rows"))))
    if recorded.get("scaled_defect_increase_positive_rows") != recorded.get("scaled_defect_increase_rows"):
        issues.append(
            issue(
                "finite_diagnostics",
                "scaled-defect-not-all-increasing",
                f"{recorded.get('scaled_defect_increase_positive_rows')!r} != {recorded.get('scaled_defect_increase_rows')!r}",
            )
        )
    return issues


def validate_rows_and_summary(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    rows = artifact.get("frontier_rows", [])
    if not isinstance(rows, list) or len(rows) != 4:
        issues.append(issue("frontier_rows", "bad-row-count", repr(len(rows) if isinstance(rows, list) else type(rows))))
    roles = {row.get("role") for row in rows if isinstance(row, dict)}
    for required in ("finite_barrier_certificate", "frontier_obstruction", "rejected_candidate"):
        if required not in roles:
            issues.append(issue("frontier_rows", "missing-role", required))
    half_width_rows = [row for row in rows if isinstance(row, dict) and row.get("id") == "nlsdf_02_half_width_buffer"]
    if half_width_rows and half_width_rows[0].get("role") not in {"live_sufficient_buffer", "frontier_obstruction"}:
        issues.append(issue("frontier_rows", "bad-half-width-role", repr(half_width_rows[0].get("role"))))
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            issues.append(issue("frontier_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "not", "only")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))

    diagnostics = artifact.get("finite_diagnostics", {})
    summary = artifact.get("summary", {})
    expected = {
        "scaled_defect_rows": diagnostics.get("scaled_defect_rows"),
        "cone_positive_rows": diagnostics.get("cone_positive_rows"),
        "half_width_positive_rows": diagnostics.get("half_width_positive_rows"),
        "one_third_positive_rows": diagnostics.get("one_third_positive_rows"),
        "one_third_failure_rows": diagnostics.get("one_third_failure_rows"),
        "scaled_defect_increase_rows": diagnostics.get("scaled_defect_increase_rows"),
        "scaled_defect_increase_positive_rows": diagnostics.get("scaled_defect_increase_positive_rows"),
        "ready_to_apply_rows": 0,
        "target_closing": False,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("exact", "half-width", "one-third", "too strong", "analytic tail theorem"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("finite", "half-width", "one-third", "nonincrease", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path, target_path: Path, artifact: dict) -> list[dict]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[dict] = []
    summary = artifact.get("summary", {})
    diagnostics = artifact.get("finite_diagnostics", {})
    result_json = artifact.get("result_json") or repo_ref(target_path)
    note_ref = artifact.get("note") or repo_ref(path)
    required_strings = [
        "# Jensen-Window PF Negative-Lambda Scaled-Defect Frontier Scout",
        "Status: finite theorem-search diagnostic",
        "This is not a proof",
        "Artifact kind: `jensen_window_pf_negative_lambda_scaled_defect_frontier_scout`",
        result_json,
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py",
        (
            "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_defect_frontier_scout.py "
            f"--target {result_json} --note {note_ref}"
        ),
        (
            "validated Jensen-window PF negative-lambda scaled-defect frontier scout: "
            f"{summary.get('scaled_defect_rows')} scaled rows, "
            f"{summary.get('cone_positive_rows')} cone rows, "
            f"{summary.get('half_width_positive_rows')} half-width rows, "
            f"{summary.get('one_third_positive_rows')} one-third rows, "
            f"{summary.get('one_third_failure_rows')} one-third failures, "
            f"{summary.get('scaled_defect_increase_positive_rows')} scaled-increase rows, 0 issues"
        ),
        "s_k = ((2*k+1)/2) * d_k",
        "0 <= s_k <= 1",
        "s_k <= 1/3",
        "s_k <= 1/2",
        f"coefficient range: A_0..A_{diagnostics.get('coefficient_k_max')}",
        f"checked contractions: x_1..x_{diagnostics.get('checked_x_max')}",
        f"exact cone rows: {diagnostics.get('cone_positive_rows')} / {diagnostics.get('scaled_defect_rows')}",
        f"half-width rows: {diagnostics.get('half_width_positive_rows')} / {diagnostics.get('scaled_defect_rows')}",
        f"one-third rows: {diagnostics.get('one_third_positive_rows')} / {diagnostics.get('scaled_defect_rows')}",
        f"one-third failure rows: {diagnostics.get('one_third_failure_rows')}",
        str(artifact.get("source_tail_barrier_scout")),
        "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
    ]
    required_strings.extend(str(ref) for ref in artifact.get("enclosure_jsonl", []))
    for required in required_strings:
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
    issues.extend(validate_note(note_path, target_path, artifact))
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
        print(json.dumps({"ok": ok, "summary": summary, "issues": issues}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-SCALED-DEFECT-FRONTIER {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda scaled-defect frontier scout: "
            f"{summary.get('scaled_defect_rows')} scaled rows, "
            f"{summary.get('cone_positive_rows')} cone rows, "
            f"{summary.get('half_width_positive_rows')} half-width rows, "
            f"{summary.get('one_third_positive_rows')} one-third rows, "
            f"{summary.get('one_third_failure_rows')} one-third failures, "
            f"{summary.get('scaled_defect_increase_positive_rows')} scaled-increase rows, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
