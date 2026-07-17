#!/usr/bin/env python3
"""Validate the compact localized order-four curvature certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

from jensen_window_pf_compound_order4_localized_curvature_compact_certificate import (  # noqa: E402
    DEFAULT_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT_JSON,
    OUTER_END,
    OUTER_START,
    TILE_WIDTH,
    cache_sha256,
    compact_certificate,
    fraction_grid,
    load_cache,
)


EXPECTED_SUMMARY = (
    "validated localized order-four compact curvature certificate: "
    "107452 H tiles, 1073 positive blocks, 1 open analytic ray, 0 issues"
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


def validate_note(path: Path) -> list[Finding]:
    findings: list[Finding] = []
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous compact interval certificate with an open analytic ray.",
        "K_1(t)<=7/(2*t^2), 319<=t<=V'(2).",
        "The deterministic cache contains `107452` adjacent mode tiles",
        "The sole remaining part of the first-summand curvature theorem",
        "K_1(t)<=7/(2*t^2) on the mode ray u>=2.",
        "not promoted to that theorem",
        "not a proof of complete order-four entry, PF-infinity, RH, or `Lambda <= 0`",
    )
    lowered = text.lower()
    for marker in required:
        if marker.lower() not in lowered:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def validate_cache_shape(records: list[dict]) -> list[Finding]:
    findings: list[Finding] = []
    expected_parameters = {
        "precision_bits": 256,
        "panels": 200,
        "window_y": 6,
        "eighth_envelope": "1/50000",
    }
    for index, record in enumerate(records):
        if record.get("index") != index:
            findings.append(finding("cache", "bad-index", (index, record.get("index"))))
            break
        if record.get("passed") is not True:
            findings.append(finding("cache", "failed-tile", index))
            break
        if record.get("parameters") != expected_parameters:
            findings.append(finding("cache", "bad-parameters", (index, record.get("parameters"))))
            break
        derivatives = record.get("h_derivatives", {})
        if sorted(derivatives) != [str(order) for order in range(2, 9)]:
            findings.append(finding("cache", "bad-derivative-orders", (index, sorted(derivatives))))
            break
        if index and records[index - 1].get("mode_right") != record.get("mode_left"):
            findings.append(finding("cache", "mode-gap", index))
            break
    if records:
        if records[0].get("mode_left") != str(OUTER_START):
            findings.append(finding("cache", "bad-start", records[0].get("mode_left")))
        if records[-1].get("mode_right") != str(OUTER_END):
            findings.append(finding("cache", "bad-end", records[-1].get("mode_right")))
    return findings


def validate(
    artifact_path: Path,
    cache_path: Path,
    note_path: Path,
) -> tuple[list[Finding], dict]:
    findings: list[Finding] = []
    artifact = load_json(artifact_path)
    tasks = [
        (index, left, right)
        for index, (left, right) in enumerate(
            fraction_grid(OUTER_START, OUTER_END, TILE_WIDTH)
        )
    ]
    try:
        records = load_cache(cache_path, tasks)
    except Exception as exc:  # pragma: no cover - command-line corruption path
        return [finding("cache", "load-failed", repr(exc))], artifact

    if len(records) != 107452:
        findings.append(finding("cache", "bad-row-count", len(records)))
    findings.extend(validate_cache_shape(records))
    cache_section = artifact.get("cache", {})
    actual_hash = cache_sha256(cache_path)
    if cache_section.get("sha256") != actual_hash:
        findings.append(
            finding("cache", "sha256-mismatch", (cache_section.get("sha256"), actual_hash))
        )

    if not findings:
        try:
            rebuilt = compact_certificate(records, cache_path, progress=False)
        except Exception as exc:  # pragma: no cover - command-line corruption path
            findings.append(finding("rebuild", "failed", repr(exc)))
        else:
            if rebuilt != artifact:
                findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))

    summary = artifact.get("summary", {})
    if summary.get("accepted_central_blocks") != 1073:
        findings.append(finding("summary", "bad-block-count", summary.get("accepted_central_blocks")))
    if summary.get("all_blocks_passed") is not True:
        findings.append(finding("summary", "blocks-not-positive", summary.get("all_blocks_passed")))
    if artifact.get("remaining_target") != "Prove K_1(t)<=7/(2t^2) on the analytic mode ray u>=2.":
        findings.append(finding("status", "bad-open-ray", artifact.get("remaining_target")))
    findings.extend(validate_note(note_path))
    return findings, artifact


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings, _artifact = validate(args.artifact, args.cache, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"localized order-four compact curvature certificate: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
