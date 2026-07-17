#!/usr/bin/env python3
"""Validate the second-next finite coefficient certificate."""

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

from jensen_window_pf_compound_order4_formal_cumulant_second_next_parity_finite_certificate import (  # noqa: E402
    COEFFICIENT_BOUNDS,
    DEFAULT_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    MODE_END,
    MODE_START,
    NORMALIZED_JET_CAPS,
    SOURCE_SECOND_NEXT,
    TAYLOR_DEGREE,
    build_artifact,
    deterministic_tasks,
    load_cache,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four second-next finite bounds: "
    "3600 Taylor blocks, 360 chunks, 7 signed coefficients, 0 issues"
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
    if len(records) != 360:
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
    actual_cache_hash = sha256(cache_path)
    if cache_section.get("sha256") != actual_cache_hash:
        findings.append(finding("cache", "sha256-mismatch", actual_cache_hash))
    source_section = artifact.get("source_second_next", {})
    actual_source_hash = sha256(SOURCE_SECOND_NEXT)
    if source_section.get("sha256") != actual_source_hash:
        findings.append(finding("source", "sha256-mismatch", actual_source_hash))
    if not findings:
        rebuilt = build_artifact(records, cache_path)
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))

    parameters = artifact.get("parameters", {})
    if (
        parameters.get("block_count") != 3_600
        or parameters.get("chunk_count") != 360
        or parameters.get("taylor_degree") != TAYLOR_DEGREE
        or parameters.get("maximum_potential_order") != 12
    ):
        findings.append(finding("summary", "bad-parameters", parameters))
    coefficients = artifact.get("coefficients", {})
    if set(coefficients) != {str(order) for order in COEFFICIENT_BOUNDS}:
        findings.append(finding("summary", "bad-coefficient-orders", sorted(coefficients)))
    normalized_jets = artifact.get("normalized_jets", {})
    if set(normalized_jets) != {str(order) for order in NORMALIZED_JET_CAPS}:
        findings.append(finding("summary", "bad-normalized-jet-orders", sorted(normalized_jets)))
    if artifact.get("normalized_jet_caps", {}).get("12") != "71":
        findings.append(finding("summary", "bad-L12-cap", artifact.get("normalized_jet_caps")))

    text = note_path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous finite interval theorem for the second-next coefficient",
        "`3600` adjacent",
        "blocks of width `1/200`",
        "centered sixth-order Arb Taylor model",
        "full-block seventh-derivative",
        "through `V^(12)`",
        "same Taylor cells also prove the normalized potential-jet caps",
        "0 < L_12(u) < 71",
        "not a point sample or floating-point fit",
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
        print(f"order-four second-next finite certificate: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
