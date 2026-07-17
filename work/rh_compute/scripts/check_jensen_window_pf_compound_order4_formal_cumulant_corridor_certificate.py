#!/usr/bin/env python3
"""Validate the finite formal-cumulant corridor certificate."""

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

from jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate import (  # noqa: E402
    DEFAULT_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MODE_END,
    MODE_START,
    SOURCE_TARGET,
    build_artifact,
    cache_sha256,
    deterministic_tasks,
    load_cache,
)


EXPECTED_SUMMARY = (
    "validated order-four formal cumulant corridor certificate: "
    "1800000 blocks, 1800 chunks, 7 positive corridors, 0 issues"
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
    tasks = deterministic_tasks()
    try:
        records = load_cache(cache_path, tasks)
    except Exception as exc:
        return [finding("cache", "load-failed", repr(exc))]
    if len(records) != 1800:
        findings.append(finding("cache", "bad-row-count", len(records)))
    if records:
        if records[0].get("mode_left") != str(MODE_START):
            findings.append(finding("cache", "bad-start", records[0].get("mode_left")))
        if records[-1].get("mode_right") != str(MODE_END):
            findings.append(finding("cache", "bad-end", records[-1].get("mode_right")))
        for previous, current in zip(records, records[1:]):
            if previous.get("mode_right") != current.get("mode_left"):
                findings.append(finding("cache", "mode-gap", current.get("chunk_index")))
                break

    cache_section = artifact.get("cache", {})
    actual_cache_hash = cache_sha256(cache_path)
    if cache_section.get("sha256") != actual_cache_hash:
        findings.append(
            finding("cache", "sha256-mismatch", (cache_section.get("sha256"), actual_cache_hash))
        )
    source_section = artifact.get("source_formal_target", {})
    actual_source_hash = cache_sha256(SOURCE_TARGET)
    if source_section.get("sha256") != actual_source_hash:
        findings.append(
            finding("source", "sha256-mismatch", (source_section.get("sha256"), actual_source_hash))
        )

    if not findings:
        rebuilt = build_artifact(records, cache_path)
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))

    parameters = artifact.get("parameters", {})
    if parameters.get("block_count") != 1_800_000 or parameters.get("chunk_count") != 1_800:
        findings.append(finding("summary", "bad-counts", parameters))
    quantities = artifact.get("quantities", {})
    if len(quantities) != 7:
        findings.append(finding("summary", "bad-quantity-count", len(quantities)))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous interval certificate for the epsilon-six formal cumulant",
        "`1800000` adjacent Arb blocks",
        "cover `2<=u<=20`",
        "lies strictly inside all",
        "This certificate concerns the exact truncated formal polynomial, not the",
        "prove the formal corridor",
        "not a proof of the exact cumulant corridors",
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
        print(f"order-four formal cumulant corridor certificate: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
