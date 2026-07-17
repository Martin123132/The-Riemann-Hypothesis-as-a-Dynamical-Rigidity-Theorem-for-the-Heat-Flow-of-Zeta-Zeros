#!/usr/bin/env python3
"""Validate the finite exact-corridor localized-curvature certificate."""

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

from jensen_window_pf_compound_order4_exact_corridor_localized_curvature_finite_certificate import (  # noqa: E402
    DEFAULT_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    REPO_ROOT,
    build_artifact,
    certify_block,
    deterministic_tasks,
    initialize_worker,
    load_cache,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four exact-corridor finite curvature: 7 rows, 0 issues, "
    "6 exact rows, 20700 mode blocks, 41400 t-collar gates, "
    "20700 positive localized blocks, 1 open asymptotic ray"
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


def validate(artifact_path: Path, note_path: Path, cache_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = load_json(artifact_path)
    tasks = deterministic_tasks()
    try:
        records = load_cache(cache_path, tasks)
        rebuilt = build_artifact(records, cache_path)
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))
    expected_summary = {
        "rows": 7,
        "exact_rows": 6,
        "open_analytic_rows": 1,
        "mode_blocks": 20_700,
        "chunks": 1_035,
        "corridors_used": 7,
        "t_collar_gates": 41_400,
        "positive_localized_blocks": 20_700,
        "finite_corridor_to_curvature_closed": True,
        "open_asymptotic_rays": 1,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))
    if cache_path.exists() and sha256(cache_path) != artifact.get("cache", {}).get("sha256"):
        findings.append(finding("cache", "sha256-mismatch", cache_path))
    extrema = artifact.get("extrema", {})
    if extrema.get("maximum_scaled_block") != 3000:
        findings.append(finding("extrema", "unexpected-maximum-block", extrema))
    if float(extrema.get("maximum_scaled_localized_upper", "inf")) >= 3.5:
        findings.append(finding("extrema", "scaled-ceiling-failed", extrema))
    if float(extrema.get("minimum_t_pad_lower", "-inf")) <= 2:
        findings.append(finding("extrema", "t-pad-failed", extrema))

    initialize_worker()
    sample_blocks = (0, 1, 2999, 3000, 3199, 3200, 3699, 3700, 4699, 4700, 5699, 5700, 10699, 10700, 20698, 20699)
    for block in sample_blocks:
        row = certify_block(block)
        if row.get("passed") is not True:
            findings.append(finding("sample", "block-failed", row))

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
        "Status: rigorous exact-corridor localized-curvature theorem",
        "width `10^-4`",
        "width `10^-3`",
        "t+-2",
        "20700",
        "41400",
        "t^2 U(t)<7/2",
        "2<=u<=20",
        "t=V'(20)",
        "only curvature segment left is `u>=20`",
        "not a proof of the remaining asymptotic ray",
    )
    for marker in required:
        if marker not in text:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = validate(args.artifact, args.note, args.cache)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-four exact-corridor finite curvature: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
