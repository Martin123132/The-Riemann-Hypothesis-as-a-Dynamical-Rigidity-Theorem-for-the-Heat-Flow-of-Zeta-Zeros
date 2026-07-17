#!/usr/bin/env python3
"""Validate the lambda-zero contiguous order-four eventual theorem."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order4_lambda0_eventual_positivity_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact,
)


EXPECTED_SUMMARY = (
    "validated Jensen-window PF compound order-four lambda-zero eventual "
    "positivity: 8 rows, 0 issues, 7 exact rows, 7 symbolic coefficients, "
    "501 finite prefix rows, 1 eventual theorem, 1 open effective splice"
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
        "symbolic_coefficients_checked": 7,
        "finite_prefix_rows": 501,
        "eventual_positivity_theorems": 1,
        "open_effective_splices": 1,
        "all_shift_lambda0_order4_proved": False,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))

    normalization = artifact.get("normalization", {})
    if normalization.get("coefficient_identity") != "A_k(0)=gamma(k)/4^(k+1)":
        findings.append(finding("normalization", "bad-coefficient-identity", normalization))
    if normalization.get("determinant_identity") != (
        "H_(4,n)[A(0)]=4^(-4*n-16)*H_(4,n)[gamma]"
    ):
        findings.append(finding("normalization", "bad-determinant-scale", normalization))

    published = artifact.get("published_ratio_input", {})
    if published.get("needed_limit") != "G_2(M)->1":
        findings.append(finding("published", "missing-G2-limit", published))
    if "1910.01227" not in published.get("arxiv", ""):
        findings.append(finding("published", "bad-source", published))

    cancellation = artifact.get("symbolic_cancellation", {})
    if cancellation.get("coefficients_h0_through_h6") != [
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        "768*G2**6",
    ]:
        findings.append(finding("cancellation", "bad-coefficients", cancellation))
    if cancellation.get("permutations_checked") != 24:
        findings.append(finding("cancellation", "bad-permutation-count", cancellation))

    eventual = artifact.get("eventual_theorem", {})
    if eventual.get("strict_sign") != (
        "there exists N_H4 such that H_(4,n)(0)>0 for every n>=N_H4"
    ):
        findings.append(finding("eventual", "bad-sign-conclusion", eventual))
    if eventual.get("threshold_effective_here") is not False:
        findings.append(finding("eventual", "threshold-overclaim", eventual))

    prefix = artifact.get("finite_prefix", {})
    if prefix.get("n_range") != [0, 500]:
        findings.append(finding("prefix", "bad-range", prefix.get("n_range")))
    if len(prefix.get("rows", [])) != 501:
        findings.append(finding("prefix", "bad-row-count", len(prefix.get("rows", []))))
    if prefix.get("all_H4_positive") is not True:
        findings.append(finding("prefix", "H4-not-closed", prefix))
    if prefix.get("all_stable_margins_positive") is not True:
        findings.append(finding("prefix", "margin-not-closed", prefix))
    if prefix.get("minimum_margin_n") != 500:
        findings.append(finding("prefix", "unexpected-minimum", prefix.get("minimum_margin_n")))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: exact eventual lambda-zero order-four positivity theorem",
        "not a proof of all-shift order four",
        "A_k(0)=gamma(k)/4^(k+1)",
        "Theorem 2.1",
        "No claim that",
        "768*G2^6*h^6",
        "H_(4,n)(0)>0 for every 0<=n<=500",
        "open splice",
        "make N_H4 explicit",
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
        print(f"lambda-zero order-four eventual positivity: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
