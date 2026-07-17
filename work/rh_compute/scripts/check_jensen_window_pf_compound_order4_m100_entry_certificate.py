#!/usr/bin/env python3
"""Validate the contiguous order-four lambda=-100 entry certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_m100_entry_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated Jensen-window PF compound order-four lambda=-100 entry: "
    "10 rows, 0 issues, 9 exact rows, 317 prefix margins, 1 analytic tail, "
    "1 all-shift entry theorem, 1 open forward handoff"
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
        "rows": 10,
        "exact_rows": 9,
        "open_analytic_rows": 1,
        "prefix_order_four_margins": 317,
        "analytic_order_four_tails": 1,
        "global_curvature_theorems": 1,
        "full_kernel_transfers": 1,
        "all_shift_order_four_entry_theorems": 1,
        "open_forward_handoffs": 1,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    if artifact.get("exact", {}).get("all_shift_entry") != (
        "H_(4,n)(-100)>0, every n>=0"
    ):
        findings.append(finding("entry", "missing-all-shift-theorem", artifact.get("exact")))

    arithmetic = artifact.get("tail_arithmetic", {})
    if arithmetic.get("buffer_endpoint_numerator") != 76_466_200:
        findings.append(finding("tail", "wrong-buffer-endpoint", arithmetic))
    if arithmetic.get("buffer_shifted_coefficients_descending") != [
        753,
        479_920,
        76_466_200,
    ]:
        findings.append(finding("tail", "wrong-buffer-polynomial", arithmetic))
    if arithmetic.get("perturbation_shifted_coefficients_descending") != [
        1,
        1902,
        1_482_055,
        604_691_300,
        135_889_991_935,
        15_877_422_036_942,
        748_002_501_678_169,
    ]:
        findings.append(finding("tail", "wrong-perturbation-polynomial", arithmetic))
    if artifact.get("prefix", {}).get("n_range") != [0, 316]:
        findings.append(finding("prefix", "wrong-range", artifact.get("prefix")))

    for ref, expected_hash in artifact.get("source_hashes", {}).items():
        path = REPO_ROOT / ref
        if not path.exists():
            findings.append(finding("source", "missing", ref))
        elif sha256(path) != expected_hash:
            findings.append(finding("source", "sha256-mismatch", ref))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(finding("artifact", "missing-source", ref))
    generator = artifact.get("generator")
    if not isinstance(generator, str) or not (REPO_ROOT / generator).exists():
        findings.append(finding("artifact", "missing-generator", generator))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: all-shift contiguous order-four entry theorem",
        "H_(4,n)(-100)>0, 0<=n<=316",
        "K_1(t)<=7/(2*t^2) for every real t>=319",
        "P_n<=4/k^2",
        "76466200",
        "H_(4,n)(-100)>0 for every integer n>=0",
        "closes the contiguous order-four entry problem",
        "next live theorem is forward propagation",
        "not\na proof of forward order-four invariance",
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
        print(f"compound order-four lambda=-100 entry: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
