#!/usr/bin/env python3
"""Validate the formal-cumulant next-parity certificate."""

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

from jensen_window_pf_compound_order4_formal_cumulant_next_parity_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
)


EXPECTED_SUMMARY = (
    "validated order-four formal cumulant next parity: 6 rows, 0 issues, "
    "4 exact rows, 42 epsilon-six audits, 7 next-parity coefficients, "
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
        "rows": 6,
        "exact_rows": 4,
        "open_rows": 2,
        "cumulant_orders": 7,
        "epsilon_six_coefficient_comparisons": 42,
        "next_parity_coefficients": 7,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    if artifact.get("scaled_hierarchy") != {
        "q^-3": [2, 3, 4],
        "q^-2": [5, 6],
        "q^-1": [7, 8],
    }:
        findings.append(finding("hierarchy", "mismatch", artifact.get("scaled_hierarchy")))
    coefficient_rows = artifact.get("coefficient_rows", {})
    if set(coefficient_rows) != {str(order) for order in range(2, 9)}:
        findings.append(finding("coefficients", "bad-orders", sorted(coefficient_rows)))
    for order_text, row in coefficient_rows.items():
        if not row.get("raw_coefficient") or not row.get("scaled_coefficient"):
            findings.append(finding("coefficients", "missing-polynomial", order_text))

    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(finding("artifact", "missing-source", ref))
    generator = artifact.get("generator")
    if not isinstance(generator, str) or not (REPO_ROOT / generator).exists():
        findings.append(finding("artifact", "missing-generator", generator))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: exact formal epsilon-eight and next-parity coefficient theorem.",
        "All 42",
        "[epsilon^j] kappa_r=0",
        "q^-3*C_r",
        "q^-2*C_r",
        "q^-1*C_r",
        "No exact-density remainder has been promoted.",
        "not a proof of the exact cumulant ray",
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
        print(f"order-four formal cumulant next parity: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
