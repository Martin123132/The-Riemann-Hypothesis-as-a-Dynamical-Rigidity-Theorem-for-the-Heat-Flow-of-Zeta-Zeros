#!/usr/bin/env python3
"""Validate the negative-lambda adaptive scaled-defect envelope matrix."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_adaptive_envelope_matrix import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_matrix,
)


REQUIRED_ROW_IDS = {
    "nlaem_01_exact_cone_prefix",
    "nlaem_02_k_monotone_frontier",
    "nlaem_03_lambda_magnitude_order",
    "nlaem_04_corner_dominance",
    "nlaem_05_fixed_buffer_rejections",
    "nlaem_06_monotone_envelope_route",
    "nlaem_07_acceptance_tests",
}

ALLOWED_ROLES = {
    "finite_certificate",
    "finite_pattern",
    "finite_frontier",
    "rejected_route",
    "live_route",
    "acceptance_gate",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda Adaptive Envelope Matrix",
    "Status: finite theorem-search diagnostic",
    "This is not a proof",
    "Artifact kind: `jensen_window_pf_negative_lambda_adaptive_envelope_matrix`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.json",
    "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m25_m50_m100_k0_k200.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_adaptive_envelope_matrix.py",
    "validated Jensen-window PF negative-lambda adaptive envelope matrix: 7 matrix rows, 0 issues, 594 k-increase rows, 398 lambda-order rows, 76 half-width failures",
    "lambdas: -25.0, -50.0, -100.0",
    "coefficient range: A_0..A_200",
    "checked contractions: x_1..x_199",
    "exact cone rows: 597 / 597",
    "k-increase rows: 594 / 594",
    "lambda-order rows: 398 / 398",
    "half-width failure rows: 76",
    "one-third failure rows: 418",
    "-25.0>=-50.0: 199 / 199 positive",
    "-50.0>=-100.0: 199 / 199 positive",
    "max scaled defect: 5.376643171065356005E-1 at lambda=-25.0, k=199",
    "min upper cone slack: 4.623356828934643995E-1 at lambda=-25.0, k=199",
    "min exact-cone margin: 1.346067912947116963E-2 at lambda=-100.0, k=1",
    "min k-increase: 4.449655594664470276E-4 at lambda=-25.0, k=198",
    "min lambda-order gap: 1.663609726220394924E-2 for -50.0>=-100.0, k=1",
    "prove s_k(lambda) increases in k without crossing 1",
    "prove s_k(lambda) decreases as |lambda| increases",
    "outputs/jensen_window_pf_negative_lambda_scaled_defect_frontier_k200_scout.md",
    "outputs/jensen_window_pf_negative_lambda_adaptive_scaled_defect_target.md",
    "outputs/jensen_window_pf_negative_lambda_defect_tail_theorem_target.md",
)


@dataclass(frozen=True)
class MatrixIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: str) -> MatrixIssue:
    return MatrixIssue(section=section, issue=name, detail=detail)


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(section: str, ref: object) -> list[MatrixIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue(section, "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue(section, "missing-ref", ref)]
    return []


def validate_top_level(artifact: dict) -> list[MatrixIssue]:
    issues: list[MatrixIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_adaptive_envelope_matrix":
        issues.append(issue("<artifact>", "bad-kind", repr(artifact.get("kind"))))
    if artifact.get("status") != "finite theorem-search diagnostic":
        issues.append(issue("<artifact>", "bad-status", repr(artifact.get("status"))))
    for key in ("source_scaled_defect_frontier", "source_adaptive_target", "generator", "checker"):
        issues.extend(validate_ref("<artifact>", artifact.get(key)))
    for ref in artifact.get("enclosure_jsonl", []):
        issues.extend(validate_ref("<artifact>", ref))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for required in ("finite", "does not prove", "continuous lambda", "adaptive", "cone entry", "lambda <= 0"):
        if required not in boundary:
            issues.append(issue("<artifact>", "weak-proof-boundary", required))
    return issues


def validate_recomputed(artifact: dict) -> list[MatrixIssue]:
    refs = artifact.get("enclosure_jsonl", [])
    if not isinstance(refs, list) or not refs:
        return [issue("recompute", "missing-enclosures", repr(refs))]
    try:
        recomputed = build_matrix([REPO_ROOT / ref for ref in refs])
    except Exception as exc:
        return [issue("recompute", "recompute-failed", f"{type(exc).__name__}: {exc}")]
    issues: list[MatrixIssue] = []
    for key in ("matrix_rows", "summary", "invariants"):
        if artifact.get(key) != recomputed.get(key):
            issues.append(issue("recompute", f"bad-{key}", "recorded artifact differs from recomputed artifact"))
    return issues


def validate_rows(artifact: dict) -> tuple[list[MatrixIssue], int, int]:
    rows = artifact.get("matrix_rows", [])
    issues: list[MatrixIssue] = []
    if not isinstance(rows, list):
        return [issue("matrix_rows", "bad-rows", repr(type(rows)))], 0, 0
    by_id = {row.get("id"): row for row in rows if isinstance(row, dict)}
    for missing in sorted(REQUIRED_ROW_IDS - set(by_id)):
        issues.append(issue(missing, "missing-row", missing))
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
        if row.get("readiness") != "not_ready_to_apply":
            issues.append(issue(row_id, "bad-readiness", repr(row.get("readiness"))))
        if row.get("role") == "live_route":
            live_routes += 1
        boundary = str(row.get("proof_boundary", "")).lower()
        if not any(marker in boundary for marker in ("finite", "not", "only", "theorem-search")):
            issues.append(issue(row_id, "weak-proof-boundary", boundary))
    return issues, len(rows), live_routes


def validate_summary(artifact: dict, row_count: int, live_routes: int) -> list[MatrixIssue]:
    summary = artifact.get("summary", {})
    expected = {
        "matrix_rows": 7,
        "ready_to_apply_rows": 0,
        "coefficient_k_max": 200,
        "checked_x_max": 199,
        "lambdas": ["-25.0", "-50.0", "-100.0"],
        "scaled_defect_rows": 597,
        "cone_positive_rows": 597,
        "half_width_positive_rows": 521,
        "half_width_failure_rows": 76,
        "one_third_positive_rows": 179,
        "one_third_failure_rows": 418,
        "k_increase_rows": 594,
        "k_increase_positive_rows": 594,
        "lambda_order_rows": 398,
        "lambda_order_positive_rows": 398,
        "target_closing": False,
    }
    issues: list[MatrixIssue] = []
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(issue("summary", f"bad-{key}", f"{summary.get(key)!r} != {value!r}"))
    if row_count != 7:
        issues.append(issue("summary", "bad-row-count", str(row_count)))
    if live_routes != 1:
        issues.append(issue("summary", "bad-live-route-count", str(live_routes)))
    expected_pairs = [
        ("-25.0>=-50.0", 199, "2.185649379581140840E-2"),
        ("-50.0>=-100.0", 199, "1.663609726220394924E-2"),
    ]
    actual_pairs = summary.get("lambda_pairs", [])
    if not isinstance(actual_pairs, list) or len(actual_pairs) != 2:
        issues.append(issue("summary", "bad-lambda-pairs", repr(actual_pairs)))
    else:
        actual = [(row.get("pair"), row.get("positive_rows"), row.get("min_gap", {}).get("sample")) for row in actual_pairs]
        if actual != expected_pairs:
            issues.append(issue("summary", "bad-lambda-pairs", repr(actual)))
    extrema = {
        "max_scaled_defect": ("-25.0", 199, "5.376643171065356005E-1"),
        "min_upper_slack": ("-25.0", 199, "4.623356828934643995E-1"),
        "min_cone_margin": ("-100.0", 1, "1.346067912947116963E-2"),
        "min_k_increase": ("-25.0", 198, "4.449655594664470276E-4"),
    }
    for key, (lam, k, sample) in extrema.items():
        row = summary.get(key, {})
        if (row.get("lam"), row.get("k"), row.get("sample")) != (lam, k, sample):
            issues.append(issue("summary", f"bad-{key}", repr(row)))
    lambda_gap = summary.get("min_lambda_gap", {})
    if (lambda_gap.get("pair"), lambda_gap.get("k"), lambda_gap.get("sample")) != (
        "-50.0>=-100.0",
        1,
        "1.663609726220394924E-2",
    ):
        issues.append(issue("summary", "bad-min-lambda-gap", repr(lambda_gap)))
    finding = str(summary.get("main_finding", "")).lower()
    for required in ("k200", "interval-increasing", "594/594", "398/398", "upper cone slack", "no all-k tail theorem"):
        if required not in finding:
            issues.append(issue("summary", "missing-main-finding-text", required))
    invariants = " ".join(str(item) for item in artifact.get("invariants", [])).lower()
    for required in ("no row", "finite", "k-monotone", "lambda-magnitude", "one-third", "half-width", "nonincrease", "lambda <= 0"):
        if required not in invariants:
            issues.append(issue("invariants", "missing-invariant-text", required))
    return issues


def validate_note(path: Path) -> list[MatrixIssue]:
    if not path.exists():
        return [issue("note", "missing-note", str(path))]
    text = path.read_text(encoding="utf-8")
    issues: list[MatrixIssue] = []
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


def validate(target_path: Path, note_path: Path) -> tuple[list[MatrixIssue], dict]:
    artifact = load_json(target_path)
    issues: list[MatrixIssue] = []
    issues.extend(validate_top_level(artifact))
    issues.extend(validate_recomputed(artifact))
    row_issues, row_count, live_routes = validate_rows(artifact)
    issues.extend(row_issues)
    issues.extend(validate_summary(artifact, row_count, live_routes))
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
            print(f"JWPF-NEG-LAMBDA-ADAPTIVE-ENVELOPE {item.section} [{item.issue}] {item.detail}")
        print(
            "validated Jensen-window PF negative-lambda adaptive envelope matrix: "
            f"{summary.get('matrix_rows')} matrix rows, {len(issues)} issues, "
            f"{summary.get('k_increase_positive_rows')} k-increase rows, "
            f"{summary.get('lambda_order_positive_rows')} lambda-order rows, "
            f"{summary.get('half_width_failure_rows')} half-width failures"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
