#!/usr/bin/env python3
"""Validate the lambda=-100 k320 cone-collar extension certificate."""

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

from jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_SOURCE_JSONL,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "m100k320_01_repaired_source",
    "m100k320_02_coefficient_positivity",
    "m100k320_03_cone_walls",
    "m100k320_04_adjacent_wall",
    "m100k320_05_new_collar_extension",
    "m100k320_06_tail_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Negative-Lambda -100 k320 Collar Extension Certificate",
    "Status: finite Arb cone-collar extension certificate",
    "Artifact kind: `jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate`",
    "work/rh_compute/results/acb_enclosures_negative_lambda_cone_entry_lam_m100_k245_k320_dps220.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate.py",
    "validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: 6 rows, 0 issues, 76 positive coefficients, 74 cone rows, 73 adjacent-wall rows, 19 new extension rows, 0 ready-to-apply rows",
    "(2*k-1)/(2*k+1) < x_k < 1",
    "x_(k+1)>x_k",
    "k=300..318",
    "analytic tail",
    "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class CertificateIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> CertificateIssue:
    return CertificateIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(ref: object) -> list[CertificateIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict, source_path: Path) -> list[CertificateIssue]:
    findings: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_m100_k320_collar_extension_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "finite Arb cone-collar extension certificate":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_enclosure_jsonl",
        "source_enclosure_summary",
        "source_k300_audit",
        "source_dominance_certificate",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("finite", "does not prove", "eventual-k", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected = {
        "certificate_rows": 6,
        "coefficient_rows": 76,
        "cone_rows": 74,
        "adjacent_wall_rows": 73,
        "new_extension_rows": 19,
        "new_collar_end_k": 318,
        "minimum_extension_log_gap_at_k": 318,
        "ready_to_apply_rows": 0,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        min_gap = Decimal(summary["minimum_extension_log_gap_lower"])
        if not Decimal("3.70e-6") < min_gap < Decimal("3.72e-6"):
            findings.append(issue("summary", "bad-min-gap", min_gap))
    except Exception as exc:
        findings.append(issue("summary", "bad-decimal", exc))

    diagnostics = artifact.get("diagnostics", {})
    expected_flags = {
        "all_source_coefficients_positive": True,
        "all_cone_rows_certified": True,
        "all_adjacent_wall_rows_certified": True,
    }
    for key, value in expected_flags.items():
        if diagnostics.get(key) is not value:
            findings.append(issue("diagnostics", f"bad-{key}", diagnostics.get(key)))
    if len(diagnostics.get("cone_detail_rows", [])) != 74:
        findings.append(issue("diagnostics", "bad-cone-detail-count", len(diagnostics.get("cone_detail_rows", []))))
    if len(diagnostics.get("wall_detail_rows", [])) != 73:
        findings.append(issue("diagnostics", "bad-wall-detail-count", len(diagnostics.get("wall_detail_rows", []))))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    if len(rows) != 6:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    if sum(row.get("readiness") == "not_ready_to_apply" for row in rows) != 1:
        findings.append(issue("rows", "bad-open-row-count", rows))
    try:
        recomputed = build_artifact(source_path)
        for key in ("diagnostics", "summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[CertificateIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue("note", "missing-text", value) for value in REQUIRED_NOTE_STRINGS if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--source-jsonl", type=Path, default=DEFAULT_SOURCE_JSONL)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings: list[CertificateIssue] = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(validate_artifact(artifact, args.source_jsonl))
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF negative-lambda -100 k320 collar extension certificate: "
            f"6 rows, {len(findings)} issues, 76 positive coefficients, 74 cone rows, "
            "73 adjacent-wall rows, 19 new extension rows, 0 ready-to-apply rows"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
