#!/usr/bin/env python3
"""Validate the exact-cumulant exact-tail certificate."""

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

from jensen_window_pf_compound_order4_exact_cumulant_exact_tail_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four exact cumulant exact tails: 9 rows, 0 issues, "
    "8 exact rows, 2 positive-coefficient polynomials, 2 closed exact tails, "
    "1 open central residual"
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
        "rows": 9,
        "exact_rows": 8,
        "open_analytic_rows": 1,
        "positive_coefficient_polynomials": 2,
        "local_curvature_ratio": "59319/100000",
        "formal_tails_inherited": 2,
        "exact_tails_closed": 2,
        "open_exact_components": 1,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    curvature = artifact.get("curvature_certificate", {})
    if curvature.get("central_upper_minimum_coefficient") != "16":
        findings.append(finding("curvature", "bad-upper-polynomial", curvature))
    if curvature.get("shifted_lower_minimum_coefficient") != "8":
        findings.append(finding("curvature", "bad-lower-polynomial", curvature))
    if artifact.get("adaptive_collar", {}).get("curvature_ratio") != "59319/100000":
        findings.append(finding("collar", "bad-curvature-ratio", artifact.get("adaptive_collar")))
    if artifact.get("tail_budgets", {}).get("finite", {}).get("target") != (
        "each exact tail <1/(500000*q^3)"
    ):
        findings.append(finding("tails", "bad-finite-target", artifact.get("tail_budgets")))
    if artifact.get("tail_budgets", {}).get("ray", {}).get("target") != (
        "each exact tail <1/(300000*u*q^3)"
    ):
        findings.append(finding("tails", "bad-ray-target", artifact.get("tail_budgets")))

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
        "Status: exact-density two-tail theorem with one open central residual.",
        "positive-coefficient substitutions",
        "u_-/u>39/40",
        "q_-/q>4/5",
        "W_u''(y)>59319/100000>1/2",
        "W_u(+-Y)>=Y^2/4",
        "each exact tail <1/(500000*q^3)",
        "each exact tail <1/(300000*u*q^3)",
        "all four tails",
        "Only the central exact-minus-formal residual remains",
        "open target",
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
        print(f"order-four exact cumulant exact tails: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
