#!/usr/bin/env python3
"""Validate the log-concave Mellin monotone-wall countermodel."""

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

from jensen_window_pf_log_concave_mellin_monotone_wall_countermodel import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    REPO_ROOT,
    build_artifact,
)


REQUIRED_ROW_IDS = {
    "lcmmc_01_log_concave_density",
    "lcmmc_02_normalized_mellin_closed_form",
    "lcmmc_03_first_contraction",
    "lcmmc_04_second_contraction",
    "lcmmc_05_monotone_wall_violation",
    "lcmmc_06_zeta_specific_handoff",
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Log-Concave Mellin Monotone-Wall Countermodel",
    "Status: exact interval countermodel gate. This is not a proof or disproof",
    "Artifact kind: `jensen_window_pf_log_concave_mellin_monotone_wall_countermodel`",
    "work/rh_compute/results/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.json",
    "python work/rh_compute/scripts/jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py",
    "python work/rh_compute/scripts/check_jensen_window_pf_log_concave_mellin_monotone_wall_countermodel.py",
    "validated Jensen-window PF log-concave Mellin monotone-wall countermodel: 6 rows, 0 issues, 2 upper-wall contractions, 1 monotone-wall violation",
    "f(y)=exp(-5*y)*1_[0,1](y)",
    "H(p)=integral_0^1 y^(p-1)*exp(-5*y)dy/Gamma(p)",
    "=gamma_lower(p,5)/(5^p*Gamma(p))",
    "Both contractions satisfy the Berwald-Borell upper wall",
    "but `x_2<x_1`",
    "outputs/jensen_window_pf_kernel_mellin_upper_wall_certificate.md",
    "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
    "outputs/signed_hankel_jensen_dependency_graph.md",
)


@dataclass(frozen=True)
class CountermodelIssue:
    section: str
    issue: str
    detail: str


def issue(section: str, name: str, detail: object) -> CountermodelIssue:
    return CountermodelIssue(section=section, issue=name, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate_ref(ref: object) -> list[CountermodelIssue]:
    if not isinstance(ref, str) or not ref:
        return [issue("artifact", "bad-ref", repr(ref))]
    if not (REPO_ROOT / ref).exists():
        return [issue("artifact", "missing-ref", ref)]
    return []


def validate_artifact(artifact: dict) -> list[CountermodelIssue]:
    findings = []
    if artifact.get("kind") != "jensen_window_pf_log_concave_mellin_monotone_wall_countermodel":
        findings.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("status") != "exact interval countermodel gate":
        findings.append(issue("artifact", "bad-status", artifact.get("status")))
    if artifact.get("date") != "2026-07-10":
        findings.append(issue("artifact", "bad-date", artifact.get("date")))
    for key in ("source_upper_wall_certificate", "source_monotone_target", "generator", "checker"):
        findings.extend(validate_ref(artifact.get(key)))
    boundary = str(artifact.get("proof_boundary", "")).lower()
    for text in ("countermodel", "not the zeta kernel", "does not disprove", "rh", "lambda <= 0"):
        if text not in boundary:
            findings.append(issue("artifact", "weak-proof-boundary", text))
    summary = artifact.get("summary", {})
    expected = {
        "countermodel_rows": 6,
        "log_concave_density_rows": 1,
        "upper_wall_rows": 2,
        "monotone_wall_violations": 1,
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            findings.append(issue("summary", f"bad-{key}", summary.get(key)))
    try:
        x1_low = Decimal(summary["x_1_lower"])
        x1_up = Decimal(summary["x_1_upper"])
        x2_low = Decimal(summary["x_2_lower"])
        x2_up = Decimal(summary["x_2_upper"])
        gap_up = Decimal(summary["x_2_minus_x_1_upper"])
        if not (Decimal(0) < x1_low <= x1_up < Decimal(1)):
            findings.append(issue("summary", "bad-x1", f"{x1_low}, {x1_up}"))
        if not (Decimal(0) < x2_low <= x2_up < Decimal(1)):
            findings.append(issue("summary", "bad-x2", f"{x2_low}, {x2_up}"))
        if not (gap_up < Decimal("-0.02")):
            findings.append(issue("summary", "weak-gap", gap_up))
    except Exception as exc:
        findings.append(issue("summary", "bad-decimal", exc))
    rows = artifact.get("rows", [])
    ids = {row.get("id") for row in rows if isinstance(row, dict)}
    if ids != REQUIRED_ROW_IDS:
        findings.append(issue("rows", "bad-row-ids", sorted(str(item) for item in ids)))
    if len(rows) != 6:
        findings.append(issue("rows", "bad-row-count", len(rows)))
    for row in rows:
        if not isinstance(row, dict):
            findings.append(issue("rows", "bad-row", row))
            continue
        if not row.get("claim") or not row.get("proof_boundary"):
            findings.append(issue(str(row.get("id")), "missing-boundary", row))
    try:
        recomputed = build_artifact()
        for key in ("summary", "rows"):
            if artifact.get(key) != recomputed.get(key):
                findings.append(issue("recompute", f"mismatch-{key}", "stored differs"))
    except Exception as exc:
        findings.append(issue("recompute", "exception", exc))
    return findings


def validate_note(path: Path) -> list[CountermodelIssue]:
    if not path.exists():
        return [issue("note", "missing", path)]
    text = path.read_text(encoding="utf-8")
    return [issue("note", "missing-text", value) for value in REQUIRED_NOTE_STRINGS if value not in text]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = []
    try:
        artifact = load_json(args.artifact)
    except Exception as exc:
        artifact = {}
        findings.append(issue("artifact", "load-failed", exc))
    if artifact:
        findings.extend(validate_artifact(artifact))
    findings.extend(validate_note(args.note))
    ok = not findings
    if args.json:
        print(json.dumps({"ok": ok, "issues": [asdict(item) for item in findings]}, indent=2, sort_keys=True))
    else:
        for finding in findings:
            print(f"ISSUE {finding.section} [{finding.issue}] {finding.detail}")
        print(
            "validated Jensen-window PF log-concave Mellin monotone-wall countermodel: "
            f"6 rows, {len(findings)} issues, 2 upper-wall contractions, 1 monotone-wall violation"
        )
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
