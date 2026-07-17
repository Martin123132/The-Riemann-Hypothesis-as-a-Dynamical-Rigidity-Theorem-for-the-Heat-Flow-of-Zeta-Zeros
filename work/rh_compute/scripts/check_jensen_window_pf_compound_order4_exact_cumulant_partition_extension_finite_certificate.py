#!/usr/bin/env python3
"""Validate the finite epsilon-11 through epsilon-14 partition extension."""

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

from jensen_window_pf_compound_order4_exact_cumulant_partition_extension_finite_certificate import (  # noqa: E402
    DEFAULT_JET_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    DEFAULT_PARTITION_CACHE,
    MODE_END,
    MODE_START,
    PARTITION_BOUNDS,
    REPO_ROOT,
    SHIFTED_END,
    SHIFTED_START,
    build_artifact,
    deterministic_tasks,
    formal_partition_extension,
    initialize_worker,
    jet_chunk_task,
    load_cache,
    partition_chunk_task,
    sha256,
)


EXPECTED_SUMMARY = (
    "validated order-four finite partition extension: 8 rows, 0 issues, "
    "7 exact rows, 4 partition orders, 78 scalar functions, "
    "5400 partition blocks, 5430 shifted-jet blocks, 5 new jet caps"
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


def validate(
    artifact_path: Path,
    note_path: Path,
    partition_cache: Path,
    jet_cache: Path,
) -> list[Finding]:
    findings: list[Finding] = []
    artifact = load_json(artifact_path)
    expected_summary = {
        "rows": 8,
        "exact_rows": 7,
        "open_analytic_rows": 1,
        "partition_orders": 4,
        "partition_scalar_functions": 78,
        "partition_blocks": 5400,
        "shifted_jet_blocks": 5430,
        "new_normalized_jet_caps": 5,
        "finite_partition_extension_closed": True,
    }
    if artifact.get("summary") != expected_summary:
        findings.append(finding("summary", "mismatch", artifact.get("summary")))

    partition_tasks = deterministic_tasks(MODE_START, MODE_END)
    jet_tasks = deterministic_tasks(SHIFTED_START, SHIFTED_END)
    try:
        partition_records = load_cache(partition_cache, partition_tasks, "partition")
        jet_records = load_cache(jet_cache, jet_tasks, "shifted_jet17")
    except Exception as exc:
        findings.append(finding("cache", "load-failed", repr(exc)))
        partition_records = []
        jet_records = []
    cache = artifact.get("cache", {})
    if partition_cache.exists() and sha256(partition_cache) != cache.get(
        "partition_sha256"
    ):
        findings.append(finding("cache", "partition-sha256-mismatch", partition_cache))
    if jet_cache.exists() and sha256(jet_cache) != cache.get("jet17_sha256"):
        findings.append(finding("cache", "jet17-sha256-mismatch", jet_cache))

    try:
        algebra = formal_partition_extension()
        coefficients = algebra.pop("compiled")
        metadata = algebra.pop("metadata")
        rebuilt = build_artifact(
            algebra,
            partition_records,
            jet_records,
            partition_cache,
            jet_cache,
        )
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))
        initialize_worker(coefficients, metadata)
        partition_samples = (0, 1, 539, 1800, 2700, 5399)
        for block in partition_samples:
            row = partition_chunk_task((block, block, 1))
            if row.get("passed") is not True:
                findings.append(finding("sample", "partition-failed", row))
        jet_samples = (0, 1, 542, 2715, 5428, 5429)
        for block in jet_samples:
            row = jet_chunk_task((block, block, 1))
            if row.get("passed") is not True:
                findings.append(finding("sample", "jet17-failed", row))

    if artifact.get("partition_bounds") != {
        str(degree): str(bound) for degree, bound in PARTITION_BOUNDS.items()
    }:
        findings.append(finding("bounds", "partition-mismatch", artifact.get("partition_bounds")))
    if artifact.get("shifted_jet17_cap") != "4000":
        findings.append(finding("bounds", "jet17-cap-mismatch", artifact.get("shifted_jet17_cap")))
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
        "Status: rigorous finite partition extension and high-jet theorem.",
        "78 scalar coefficient functions",
        "||Z_11||_1<2",
        "||Z_12||_1<2",
        "||Z_13||_1<21/10",
        "||Z_14||_1<12/5",
        "L_17(v)<4000",
        "5400",
        "5430",
        "epsilon-fifteen Bell remainder",
        "No exact cumulant",
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
    parser.add_argument("--partition-cache", type=Path, default=DEFAULT_PARTITION_CACHE)
    parser.add_argument("--jet-cache", type=Path, default=DEFAULT_JET_CACHE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    findings = validate(
        args.artifact, args.note, args.partition_cache, args.jet_cache
    )
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-four finite partition extension: {len(findings)} issues")
        return 1
    print(EXPECTED_SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
