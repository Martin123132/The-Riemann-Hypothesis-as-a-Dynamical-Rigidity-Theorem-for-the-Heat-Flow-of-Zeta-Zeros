#!/usr/bin/env python3
"""Validate the negative-lambda tail-barrier scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from jensen_window_pf_negative_lambda_tail_barrier_scout import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_diagnostics,
)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def issue(section: str, name: str, detail: str) -> dict:
    return {"section": section, "issue": name, "detail": detail}


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
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_tail_barrier_scout":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in (
        "source_prefix_scout",
        "source_finite_collar_contract",
        "source_cone_entry_target",
        "generator",
        "checker",
    ):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    for key in ("result_json", "note"):
        if key in artifact:
            issues.extend(validate_ref("<artifact>", artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in (
        "finite theorem-search",
        "defect inequalities",
        "does not prove an all-k tail theorem",
        "collared finite flow theorem",
        "cone entry",
        "jwpf_06",
        "lambda <= 0 unsettled",
    ):
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

    for row_key, positive_key in (
        ("cone_buffer_rows", "cone_buffer_positive_rows"),
        ("defect_monotone_rows", "defect_monotone_positive_rows"),
        ("scaled_defect_increase_rows", "scaled_defect_increase_positive_rows"),
    ):
        if recorded.get(positive_key) != recorded.get(row_key):
            issues.append(
                issue(
                    "finite_diagnostics",
                    f"not-all-positive-{positive_key}",
                    f"{recorded.get(positive_key)!r} != {recorded.get(row_key)!r}",
                )
            )
    if recorded.get("rejected_candidate_count") != 1:
        issues.append(issue("finite_diagnostics", "bad-rejected-candidate-count", repr(recorded.get("rejected_candidate_count"))))
    return issues


def validate_rows_and_summary(artifact: dict) -> list[dict]:
    issues: list[dict] = []
    rows = artifact.get("tail_barrier_rows", [])
    if not isinstance(rows, list) or len(rows) != 5:
        issues.append(issue("tail_barrier_rows", "bad-row-count", repr(len(rows) if isinstance(rows, list) else type(rows))))
    rejected = 0
    for row in rows if isinstance(rows, list) else []:
        if not isinstance(row, dict):
            issues.append(issue("tail_barrier_rows", "bad-row", repr(row)))
            continue
        row_id = str(row.get("id", "<missing-id>"))
        for key in ("id", "role", "claim", "proof_boundary"):
            if key not in row:
                issues.append(issue(row_id, "missing-field", key))
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("only", "not", "finite", "target", "rejects")):
            issues.append(issue(row_id, "weak-row-boundary", boundary))
        if row.get("role") == "rejected_candidate":
            rejected += 1
    if rejected != 1:
        issues.append(issue("tail_barrier_rows", "bad-rejected-count", repr(rejected)))

    summary = artifact.get("summary", {})
    diagnostics = artifact.get("finite_diagnostics", {})
    expected_summary = {
        "checked_x_max": diagnostics.get("checked_x_max"),
        "cone_buffer_rows": diagnostics.get("cone_buffer_rows"),
        "defect_monotone_rows": diagnostics.get("defect_monotone_rows"),
        "one_third_buffer_rows": diagnostics.get("scaled_defect_one_third_rows"),
        "scaled_defect_increase_rows": diagnostics.get("scaled_defect_increase_rows"),
        "rejected_candidate_count": diagnostics.get("rejected_candidate_count"),
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if summary.get("ready_to_apply_rows") != 0:
        issues.append(issue("summary", "bad-ready_to_apply_rows", repr(summary.get("ready_to_apply_rows"))))
    if summary.get("target_closing") is not False:
        issues.append(issue("summary", "bad-target_closing", repr(summary.get("target_closing"))))
    finding = str(summary.get("main_finding", "")).lower()
    next_collar_x = int(diagnostics.get("checked_x_max", 0)) + 1
    next_coefficient = int(diagnostics.get("coefficient_k_max", 0)) + 1
    for required in ("one-third", "defect", "scaled defect", "rejected", f"x_{next_collar_x}", f"a_{next_coefficient}", "tail theorem"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("ready_to_apply", "one-third", "rejected", "tail", "lambda <= 0"):
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
    active_depth = int(diagnostics.get("checked_x_max", 0)) - 2
    next_active_depth = active_depth + 1
    next_collar_x = int(diagnostics.get("checked_x_max", 0)) + 1
    next_coefficient = int(diagnostics.get("coefficient_k_max", 0)) + 1
    required_strings = [
        "# Jensen-Window PF Negative-Lambda Tail-Barrier Scout",
        "Status: finite theorem-search diagnostic",
        "This is not a proof",
        "Artifact kind: `jensen_window_pf_negative_lambda_tail_barrier_scout`",
        artifact.get("result_json") or repo_ref(target_path),
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_tail_barrier_scout.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_tail_barrier_scout.py",
        (
            "validated Jensen-window PF negative-lambda tail-barrier scout: "
            f"{summary.get('cone_buffer_rows')} cone-buffer rows, "
            f"{summary.get('defect_monotone_rows')} defect-monotone rows, "
            f"{summary.get('one_third_buffer_rows')} one-third-buffer rows, "
            f"{summary.get('scaled_defect_increase_rows')} scaled-defect increase rows, "
            f"{summary.get('rejected_candidate_count')} rejected candidate, 0 issues"
        ),
        "0 <= d_k <= 2/(2*k+1)",
        "d_{k+1} <= d_k",
        f"coefficient range: A_0..A_{diagnostics.get('coefficient_k_max')}",
        f"checked contractions: x_1..x_{diagnostics.get('checked_x_max')}",
        f"cone-buffer rows: {diagnostics.get('cone_buffer_positive_rows')} / {diagnostics.get('cone_buffer_rows')}",
        (
            "one-third width buffer rows: "
            f"{diagnostics.get('scaled_defect_one_third_rows')} / {diagnostics.get('scaled_defect_rows')}"
        ),
        f"defect monotone rows: {diagnostics.get('defect_monotone_positive_rows')} / {diagnostics.get('defect_monotone_rows')}",
        (
            "scaled-defect increase rows: "
            f"{diagnostics.get('scaled_defect_increase_positive_rows')} / {diagnostics.get('scaled_defect_increase_rows')}"
        ),
        "Rejected Shortcut",
        f"purely finite upgrade to `K={next_active_depth}` needs",
        f"`x_{next_collar_x}`, hence `A_{next_coefficient}`",
        "outputs/jensen_window_pf_negative_lambda_cone_entry_prefix_scout.md",
        "outputs/jensen_window_pf_negative_lambda_finite_collar_contract.md",
        "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
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
        "cone entry is proved",
        "all-k tail theorem is proved",
        "collared finite flow theorem is proved",
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
    issues, summary = validate(args.target, args.note)
    ok = not issues
    if args.json:
        print(json.dumps({"ok": ok, "summary": summary, "issues": issues}, indent=2, sort_keys=True))
    else:
        for item in issues:
            print(f"JWPF-NEG-LAMBDA-TAIL-BARRIER {item['section']} [{item['issue']}] {item['detail']}")
        print(
            "validated Jensen-window PF negative-lambda tail-barrier scout: "
            f"{summary.get('cone_buffer_rows')} cone-buffer rows, "
            f"{summary.get('defect_monotone_rows')} defect-monotone rows, "
            f"{summary.get('one_third_buffer_rows')} one-third-buffer rows, "
            f"{summary.get('scaled_defect_increase_rows')} scaled-defect increase rows, "
            f"{summary.get('rejected_candidate_count')} rejected candidate, "
            f"{len(issues)} issues"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
