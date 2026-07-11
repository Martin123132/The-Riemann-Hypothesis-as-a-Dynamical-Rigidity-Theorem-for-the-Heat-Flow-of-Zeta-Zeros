#!/usr/bin/env python3
"""Validate the first-summand saddle-wall theorem target and scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_first_summand_saddle_wall_target import (  # noqa: E402
    DEFAULT_FULL_SOURCE_JSONL,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "fsswt_01_exact_first_moment_wall",
    "fsswt_02_exact_saddle_geometry",
    "fsswt_03_all_k_saddle_bracket",
    "fsswt_04_high_precision_scout",
    "fsswt_05_positive_scaled_profile",
    "fsswt_06_full_kernel_overlap",
    "fsswt_07_formal_saddle_scale",
    "fsswt_08_quantitative_wall_target",
    "fsswt_09_full_wall_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda First-Summand Saddle-Wall Target",
    "Status: exact saddle geometry with analytic wall closure",
    "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_saddle_wall_target`",
    "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.json",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.py",
    "validated Jensen-window PF negative-lambda first-summand saddle-wall closure: 9 rows, 0 issues, 9 positive samples, 9 quarter-k2 samples, 9 bracketed saddles, 0 open requirements, 2 ready-to-apply rows",
    "S_k''(u)=-2*k/u^2-200-16*q-96*q/(2*q-3)^2<0",
    "(log(k)+22/25)/8 < s_k < log(k)/4",
    "L_k^(1)>=1/(4*k^2), k>=319",
    "0.9264<=u<=5",
    "u>=5",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
    "|L_k-L_k^(1)|<=16/(k-1)^6",
    "finite upper truncation is not certified",
    "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
    "outputs/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class TargetIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> TargetIssue:
    return TargetIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(ref: object) -> list[TargetIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict, source_path: Path) -> list[TargetIssue]:
    findings: list[TargetIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_first_summand_saddle_wall_target":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "exact saddle geometry with analytic wall closure":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_dominance_certificate",
        "source_collar_extension",
        "source_paired_remainder_certificate",
        "source_paired_ray_certificate",
        "source_full_enclosure_jsonl",
        "source_cone_entry_target",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("exact", "finite", "does not prove", "all-k", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected = {
        "target_rows": 9,
        "sample_rows": 9,
        "positive_sample_rows": 9,
        "target_sample_rows": 9,
        "bracketed_saddle_rows": 9,
        "minimum_scaled_k2_log_gap_at_k": 300,
        "scaled_profile_strictly_increasing": True,
        "open_requirement_rows": 0,
        "ready_to_apply_rows": 2,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        minimum = Decimal(summary["minimum_scaled_k2_log_gap"])
        if not Decimal("0.365") < minimum < Decimal("0.366"):
            findings.append(issue("summary", "bad-minimum-scaled-gap", minimum))
    except Exception as exc:
        findings.append(issue("summary", "bad-decimal", exc))

    diagnostics = artifact.get("diagnostics", {})
    rows = diagnostics.get("sample_rows", [])
    if len(rows) != 9:
        findings.append(issue("diagnostics", "bad-sample-count", len(rows)))
    if not all(row.get("positive_log_gap") and row.get("above_quarter_k2_target") for row in rows):
        findings.append(issue("diagnostics", "bad-sample-sign", rows))
    if not all(row.get("saddle_inside_bracket") for row in rows):
        findings.append(issue("diagnostics", "bad-saddle-bracket", rows))
    cross = diagnostics.get("full_kernel_k300_crosscheck", {})
    try:
        if not Decimal(cross["absolute_midpoint_difference"]) < Decimal("1e-30"):
            findings.append(issue("diagnostics", "weak-k300-crosscheck", cross))
        if not Decimal(diagnostics["exact_geometry"]["wall_transfer_margin_at_k319_lower"]) > Decimal("1e-6"):
            findings.append(issue("diagnostics", "weak-transfer-margin", diagnostics["exact_geometry"]))
    except Exception as exc:
        findings.append(issue("diagnostics", "bad-crosscheck-decimal", exc))

    target_rows = artifact.get("rows", [])
    ids = {row.get("id") for row in target_rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    if len(target_rows) != 9:
        findings.append(issue("rows", "bad-row-count", len(target_rows)))
    if sum(row.get("role") == "open_requirement" for row in target_rows) != 0:
        findings.append(issue("rows", "bad-open-requirement-count", target_rows))
    try:
        recomputed = build_artifact(source_path)
        for key in ("diagnostics", "summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[TargetIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue("note", "missing-text", value) for value in REQUIRED_NOTE_STRINGS if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--full-source-jsonl", type=Path, default=DEFAULT_FULL_SOURCE_JSONL)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings: list[TargetIssue] = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(validate_artifact(artifact, args.full_source_jsonl))
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF negative-lambda first-summand saddle-wall closure: "
            f"9 rows, {len(findings)} issues, 9 positive samples, 9 quarter-k2 samples, "
            "9 bracketed saddles, 0 open requirements, 2 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
