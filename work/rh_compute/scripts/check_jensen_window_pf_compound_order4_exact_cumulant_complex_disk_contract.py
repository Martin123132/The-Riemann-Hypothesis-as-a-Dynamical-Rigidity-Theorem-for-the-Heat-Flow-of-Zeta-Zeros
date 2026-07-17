#!/usr/bin/env python3
"""Validate the exact-cumulant complex-disk contract."""

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

from jensen_window_pf_compound_order4_exact_cumulant_complex_disk_contract import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four exact cumulant complex-disk contract: "
    "10 rows, 0 issues, 9 exact rows, 70 cumulant audits, "
    "2 sufficient partition targets, 1 open central-tail theorem"
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
        "partition_degrees": 10,
        "cumulant_coefficient_audits": 70,
        "finite_partition_target": "1/(100000*q^3)",
        "ray_partition_target": "1/(20000*u*q^3)",
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    audit = artifact.get("formal_partition", {}).get("cumulant_audit", {})
    if audit.get("coefficient_comparisons") != 70:
        findings.append(finding("audit", "bad-comparison-count", audit))
    if len(artifact.get("formal_partition", {}).get("partition_rows", {})) != 10:
        findings.append(finding("partition", "bad-row-count", artifact.get("formal_partition")))
    finite = artifact.get("cauchy_budgets", {}).get("finite", {})
    ray = artifact.get("cauchy_budgets", {}).get("ray", {})
    if finite.get("scaled_cumulant_cap") != "532/61875":
        findings.append(finding("finite", "bad-scaled-cap", finite))
    if ray.get("scaled_cumulant_cap") != "7693/1237500/u":
        findings.append(finding("ray", "bad-scaled-cap", ray))

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
        "Status: exact complex-disk reduction with open partition residual.",
        "audits all 70 cumulant coefficients",
        "|P_u^[10](z)-1|<1/100",
        "|log P^[10]-S^[10]|<1/(7500*q^3)",
        "|log P^[10]-S^[10]|<1/(100000*u*q^3)",
        "|A_u(z)-P_u^[10](z)|<1/(100000*q^3)",
        "|A_u(z)-P_u^[10](z)|<1/(20000*u*q^3)",
        "not yet proved",
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
        print(f"order-four exact cumulant complex-disk contract: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
