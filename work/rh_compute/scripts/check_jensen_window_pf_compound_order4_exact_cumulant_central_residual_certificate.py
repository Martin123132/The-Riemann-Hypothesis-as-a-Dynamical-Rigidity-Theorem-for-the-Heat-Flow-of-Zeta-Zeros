#!/usr/bin/env python3
"""Validate the exact-cumulant central-residual certificate."""

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

from jensen_window_pf_compound_order4_exact_cumulant_central_residual_certificate import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four exact cumulant central residual: 11 rows, 0 issues, "
    "11 exact rows, 222 Bell terms, 2 central regimes, 4 inherited tails, "
    "0 open partition components"
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
        "rows": 11,
        "exact_rows": 11,
        "open_analytic_rows": 0,
        "partition_extension_orders": 4,
        "bell_polynomial_terms": 222,
        "formal_tails_closed": 2,
        "exact_tails_closed": 2,
        "central_residuals_closed": 2,
        "global_partition_residual_closed": True,
        "exact_minus_epsilon_ten_cumulant_budgets_closed": True,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    finite = artifact.get("finite", {})
    ray = artifact.get("ray", {})
    if finite.get("bell", {}).get("integer_cap") != 452_921_174_468:
        findings.append(finding("finite", "bad-bell-cap", finite.get("bell")))
    if finite.get("budget", {}).get("proved_bound") != (
        "central residual <1/(500000*q^3), 2<=u<=20"
    ):
        findings.append(finding("finite", "bad-central-bound", finite.get("budget")))
    if ray.get("budget", {}).get("proved_bound") != (
        "central residual <1/(300000*u*q^3), u>=20"
    ):
        findings.append(finding("ray", "bad-central-bound", ray.get("budget")))
    if float(finite.get("budget", {}).get("budget_ratio_upper", "inf")) >= 0.364:
        findings.append(finding("finite", "budget-ratio-too-large", finite.get("budget")))
    if float(ray.get("budget", {}).get("budget_ratio_upper", "inf")) >= 1e-62:
        findings.append(finding("ray", "budget-ratio-too-large", ray.get("budget")))

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
        "Status: global exact central residual and unit-disk partition theorem.",
        "Z_11,...,Z_14",
        "epsilon-fifteen Bell-polynomial",
        "seventeenth-order potential Taylor remainder",
        "central residual <1/(500000*q^3)",
        "central residual <1/(300000*u*q^3)",
        "|A_u-P_u^[10]|<1/(100000*q^3)",
        "|A_u-P_u^[10]|<1/(20000*u*q^3)",
        "9/1000",
        "1/(100u)",
        "exact-density remainder is no longer open",
        "not a proof of order-four entry",
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
        print(f"order-four exact cumulant central residual: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
