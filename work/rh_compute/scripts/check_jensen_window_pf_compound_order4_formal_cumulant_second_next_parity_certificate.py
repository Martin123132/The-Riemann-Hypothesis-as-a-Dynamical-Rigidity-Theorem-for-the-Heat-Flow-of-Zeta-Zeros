#!/usr/bin/env python3
"""Validate the formal-cumulant second-next parity certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


EXPECTED_SUMMARY = (
    "validated order-four formal cumulant second-next parity: 5 rows, 0 issues, "
    "3 exact rows, 56 epsilon-eight audits, 7 second-next coefficients, "
    "2 open analytic rows"
)


@dataclass(frozen=True)
class Finding:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> Finding:
    return Finding(section=section, issue=issue, detail=str(detail))


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def validate(artifact_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = load_json(artifact_path)
    try:
        rebuilt = build_artifact()
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))
    expected_summary = {
        "rows": 5,
        "exact_rows": 3,
        "open_rows": 2,
        "cumulant_orders": 7,
        "epsilon_eight_coefficient_comparisons": 56,
        "second_next_coefficients": 7,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    if artifact.get("scaled_hierarchy") != {
        "q^-4": [2, 3, 4],
        "q^-3": [5, 6],
        "q^-2": [7, 8],
    }:
        findings.append(finding("hierarchy", "mismatch", artifact.get("scaled_hierarchy")))
    coefficient_rows = artifact.get("coefficient_rows", {})
    for order_text, row in coefficient_rows.items():
        expected_terms = 42 if int(order_text) % 2 == 0 else 30
        if row.get("raw_coefficient_terms") != expected_terms:
            findings.append(finding("coefficients", "bad-term-count", (order_text, row)))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(finding("artifact", "missing-source", ref))
    generator = artifact.get("generator")
    if not isinstance(generator, str) or not (REPO_ROOT / generator).exists():
        findings.append(finding("artifact", "missing-generator", generator))
    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: exact formal epsilon-ten and second-next parity theorem.",
        "All 56 coefficients",
        "q^-4*D_r",
        "q^-3*D_r",
        "q^-2*D_r",
        "not a claim that",
        "beyond-epsilon-ten central and two-tail estimates remain open",
        "not a proof of the exact cumulant corridors",
    )
    for marker in required:
        if marker not in text:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = validate(args.artifact, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-four formal cumulant second-next parity: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
