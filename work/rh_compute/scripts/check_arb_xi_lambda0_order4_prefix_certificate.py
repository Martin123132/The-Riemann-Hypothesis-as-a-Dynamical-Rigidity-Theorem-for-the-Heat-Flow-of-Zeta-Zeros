#!/usr/bin/env python3
"""Validate the cached Arb Xi lambda-zero order-four prefix."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from arb_xi_lambda0_order4_prefix_certificate import (  # noqa: E402
    DEFAULT_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    build_artifact_from_cache,
)


EXPECTED_SUMMARY = (
    "validated Arb Xi lambda-zero order-four prefix: 507 coefficients, "
    "501 positive H4 rows, 501 positive stable margins, 0 inconclusive, 0 issues"
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


def validate(artifact_path: Path, cache_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = load_json(artifact_path)
    try:
        rebuilt = build_artifact_from_cache(cache_path)
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))

    expected = {
        "coefficient_rows": 507,
        "prefix_rows": 501,
        "positive_coefficients": 507,
        "positive_H4_rows": 501,
        "positive_stable_margin_rows": 501,
        "inconclusive_rows": 0,
    }
    if artifact.get("summary") != expected:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    method = artifact.get("method", {})
    if method.get("precision_bits") != 24576:
        findings.append(finding("method", "bad-precision", method))
    if method.get("series_order") != 1016:
        findings.append(finding("method", "bad-series-order", method))
    if "xi(1/2+z)" not in method.get("method", ""):
        findings.append(finding("method", "bad-series-method", method))
    finite = artifact.get("finite", {})
    if finite.get("n_range") != [0, 500]:
        findings.append(finding("finite", "bad-n-range", finite.get("n_range")))
    if finite.get("coefficient_range") != [0, 506]:
        findings.append(finding("finite", "bad-k-range", finite.get("coefficient_range")))
    if finite.get("minimum_margin_n") != 500:
        findings.append(finding("finite", "unexpected-minimum", finite.get("minimum_margin_n")))
    if len(finite.get("rows", [])) != 501:
        findings.append(finding("finite", "bad-row-count", len(finite.get("rows", []))))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous direct-Xi-series finite certificate",
        "This is not a proof of an all-shift order-four theorem",
        "24576",
        "A_k(0)>0 for every 0<=k<=506",
        "H_(4,n)(0)>0 for every 0<=n<=500",
        "extends the previous rigorous endpoint prefix",
        "eventual Xi-asymptotic theorem",
    )
    for marker in required:
        if marker not in text:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = validate(args.artifact, args.cache, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"Arb Xi lambda-zero order-four prefix: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
