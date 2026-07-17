#!/usr/bin/env python3
"""Validate the global exact order-four cumulant corridor theorem."""

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

from jensen_window_pf_compound_order4_exact_cumulant_corridor_theorem import (  # noqa: E402
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four exact cumulant corridors: 8 rows, 0 issues, "
    "7 exact rows, 7 global exact corridors, 2 strict reserve regimes, "
    "1 open localized-curvature composition"
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
        "open_analytic_rows": 1,
        "exact_corridors": 7,
        "finite_corridors_closed": 7,
        "ray_corridors_closed": 7,
        "global_exact_corridors_closed": True,
        "open_localized_curvature_compositions": 1,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    if artifact.get("finite_reserve") != "79999/10000000":
        findings.append(finding("reserve", "bad-finite", artifact.get("finite_reserve")))
    if artifact.get("ray_reserve") != "29/(1000u)":
        findings.append(finding("reserve", "bad-ray", artifact.get("ray_reserve")))
    for ref, expected_hash in artifact.get("source_hashes", {}).items():
        path = REPO_ROOT / ref
        if not path.exists():
            findings.append(finding("source", "missing", ref))
        elif sha256(path) != expected_hash:
            findings.append(finding("source", "sha256-mismatch", ref))
    for ref in artifact.get("sources", []):
        if not isinstance(ref, str) or not (REPO_ROOT / ref).exists():
            findings.append(finding("artifact", "missing-source", ref))
    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: global exact cumulant corridor theorem",
        "all seven candidate corridors",
        "79999/10000000",
        "29/(1000u)",
        "exact-density theorem is closed",
        "curvature ray is not yet",
        "samples are not",
        "not a proof of the curvature ray",
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
        print(f"order-four exact cumulant corridors: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
