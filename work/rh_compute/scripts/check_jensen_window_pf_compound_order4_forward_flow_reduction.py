#!/usr/bin/env python3
"""Validate the contiguous order-four forward-flow reduction."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_forward_flow_reduction import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


EXPECTED_SUMMARY = (
    "validated Jensen-window PF compound order-four forward flow: 8 rows, "
    "0 issues, 7 exact rows, 3 exact identities, 2 cooperative flow lemmas, "
    "1 maximum-principle reduction, 1 open spatial-tail bound"
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
        "rows": 8,
        "exact_rows": 7,
        "exact_identities": 3,
        "cooperative_flow_lemmas": 2,
        "stable_gap_factorizations": 1,
        "maximum_principle_reductions": 1,
        "open_spatial_tail_bounds": 1,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    exact = artifact.get("exact_flow", {})
    if exact.get("a_n") != "(4*n+26)*H_(3,n)/H_(3,n+1)":
        findings.append(finding("flow", "wrong-off-diagonal", exact))
    stable = artifact.get("stable_gap", {})
    if stable.get("factorization") != "H_(4,n)=S_n*F_n":
        findings.append(finding("stable", "wrong-factorization", stable))
    if stable.get("alpha_positive") is not True:
        findings.append(finding("stable", "nonpositive-alpha", stable))
    maximum = artifact.get("maximum_principle", {})
    if "sup_(-100<=lambda<=L,n>=0)" not in maximum.get(
        "sufficient_coefficient_bound", ""
    ):
        findings.append(finding("maximum", "missing-uniform-target", maximum))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: exact cooperative order-four flow reduction",
        "This is not a proof of forward",
        "A_j'=(4*j+2)*A_(j+1)",
        "Q_n'=a_n*Q_(n+1)+b_n*Q_n",
        "order five",
        "H_(4,n)=S_n*F_n",
        "|F_n|<=1",
        "z_n=F_n/(n+1)",
        "sole forward-propagation blocker",
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
        print(f"compound order-four forward flow: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
