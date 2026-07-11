#!/usr/bin/env python3
"""Validate the T=1156 zeta monotone-wall counterexample certificate."""

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

from jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    DEFAULT_SOURCE_JSONL,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "t1156mwc_01_rigorous_coefficient_source",
    "t1156mwc_02_exact_contraction_reduction",
    "t1156mwc_03_upper_wall_occupancy",
    "t1156mwc_04_strict_monotone_wall_violation",
    "t1156mwc_05_polynomial_witness",
    "t1156mwc_06_fixed_k_promotion_gate",
    "t1156mwc_07_revised_tail_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF T=1156 Monotone-Wall Counterexample Certificate",
    "Status: interval zeta-kernel counterexample certificate",
    "Artifact kind: `jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate`",
    "work/rh_compute/results/acb_enclosures_lambda_m1156_k119_k122_dps250.jsonl",
    "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate.py",
    "validated Jensen-window PF T=1156 monotone-wall counterexample certificate: 7 rows, 0 issues, 4 coefficient enclosures, 1 zeta monotone-wall violation",
    "x_120=A_121*A_119/A_120^2",
    "x_121=A_122*A_120/A_121^2",
    "failure of",
    "the adjacent-k wall `x_(k+1)>=x_k`",
    "fixed-k=22 certificate",
    "choose a moderate fixed negative lambda",
    "outputs/jensen_window_pf_heat_flow_cone_entry_asymptotic_target.md",
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
    if artifact.get("kind") != "jensen_window_pf_negative_lambda_t1156_monotone_wall_counterexample_certificate":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "interval zeta-kernel counterexample certificate":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in (
        "source_enclosure_jsonl",
        "source_enclosure_summary",
        "source_enclosure_generator",
        "source_tail_bounds",
        "source_fixed_k_certificate",
        "source_cone_entry_target",
        "generator",
        "checker",
    ):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("one actual zeta-kernel", "another lambda", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))

    summary = artifact.get("summary", {})
    expected = {
        "certificate_rows": 7,
        "source_coefficient_rows": 4,
        "upper_wall_contractions": 2,
        "zeta_monotone_wall_violations": 1,
        "violated_wall_index": 120,
        "fixed_k_t1156_all_k_promotion_blocked": True,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        if not Decimal(summary["gap_upper"]) < Decimal("-1.68e-8"):
            findings.append(issue("summary", "weak-gap", summary.get("gap_upper")))
        if not Decimal(summary["log_ratio_upper"]) < Decimal("-1.68e-8"):
            findings.append(issue("summary", "weak-log-ratio", summary.get("log_ratio_upper")))
    except Exception as exc:
        findings.append(issue("summary", "bad-decimal", exc))

    diagnostics = artifact.get("diagnostics", {})
    expected_diagnostics = {
        "lambda": "-1156",
        "T": 1156,
        "violated_wall_index": 120,
        "coefficient_indices": [119, 120, 121, 122],
        "all_source_coefficients_positive": True,
        "both_contractions_inside_upper_wall": True,
        "monotone_wall_strictly_violated": True,
    }
    for key, value in expected_diagnostics.items():
        if diagnostics.get(key) != value:
            findings.append(issue("diagnostics", f"bad-{key}", diagnostics.get(key)))

    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    if len(rows) != 7:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    for row in rows:
        if not isinstance(row, dict) or not row.get("claim") or not row.get("proof_boundary"):
            findings.append(issue("rows", "incomplete-row", row))

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
            "validated Jensen-window PF T=1156 monotone-wall counterexample certificate: "
            f"7 rows, {len(findings)} issues, 4 coefficient enclosures, 1 zeta monotone-wall violation"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
