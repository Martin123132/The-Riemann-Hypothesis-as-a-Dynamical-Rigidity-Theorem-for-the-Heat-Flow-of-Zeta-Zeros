#!/usr/bin/env python3
"""Validate the order-eight shifted-jet t=699..999 certificate."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate import (  # noqa: E402
    DEFAULT_BLOCK_CACHE,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    END_T,
    START_T,
    build_artifact,
    deterministic_tasks,
    load_block_cache,
    load_source,
)


@dataclass(frozen=True)
class CertificateIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> CertificateIssue:
    return CertificateIssue(section, code, str(detail))


def validate_artifact(path: Path, cache: Path) -> tuple[dict, list[CertificateIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[CertificateIssue] = []
    if artifact.get("kind") != "jensen_window_pf_compound_order8_shifted_jet_t699_t999_certificate":
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "699<=t<=999" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    source = load_source()
    tasks = deterministic_tasks(source)
    try:
        blocks = load_block_cache(cache, tasks)
    except Exception as exc:
        issues.append(issue("cache", "load-failure", exc))
        blocks = []
    if len(blocks) != len(tasks):
        issues.append(issue("cache", "incomplete", f"{len(blocks)}/{len(tasks)}"))
    else:
        if blocks[0].get("anchor") != str(START_T):
            issues.append(issue("cache", "bad-start", blocks[0].get("anchor")))
        if blocks[-1].get("right") != str(END_T):
            issues.append(issue("cache", "bad-end", blocks[-1].get("right")))
        for previous, current in zip(blocks, blocks[1:]):
            if previous.get("right") != current.get("anchor"):
                issues.append(issue("cache", "coverage-gap", current.get("index")))
                break
        for block in blocks:
            if Decimal(block.get("scaled_curvature_upper", "Infinity")) >= Decimal("2000"):
                issues.append(issue("cache", "curvature-failure", block.get("index")))
                break
            if Decimal(
                block.get("common_coordinates", {})
                .get("U", {})
                .get("dimensionless_common_lower", "-Infinity")
            ) <= 0:
                issues.append(issue("cache", "nonpositive-U", block.get("index")))
                break
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 6,
        "ready_rows": 5,
        "open_rows": 1,
        "blocks": len(tasks),
        "passed_blocks": len(tasks),
        "point_shift_count": 13,
        "stable_layers": 5,
        "continuous_curvature_theorems": 1,
        "open_ray_targets": 1,
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    rows = artifact.get("rows", [])
    if len(rows) != 6 or len({row.get("id") for row in rows}) != 6:
        issues.append(issue("rows", "bad-rows", rows))
    if len(blocks) == len(tasks):
        try:
            canonical = build_artifact(blocks, source)
        except Exception as exc:
            issues.append(issue("rebuild", "exception", exc))
        else:
            if artifact != canonical:
                issues.append(issue("rebuild", "artifact-drift", path))
    return artifact, issues


def validate_note(path: Path) -> list[CertificateIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous continuous first-summand curvature theorem",
        "This is not a proof",
        "common t+-6 collar",
        "thirteen rigorous shifted quadratures",
        "Five staged localized tent",
        "s_1''(t)<=2000/t^2 for every real 699<=t<=999",
        "prove s_1''(t)<=4000/t^2 for every real t>=999",
        "outputs/jensen_window_pf_compound_order8_first_summand_curvature_bridge.md",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-eight entry is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_BLOCK_CACHE)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact, issues = validate_artifact(args.artifact, args.cache)
    issues.extend(validate_note(args.note))
    if issues:
        for finding in issues:
            print(
                f"ORDER8-SHIFTED-JET {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-eight shifted-jet certificate: {len(issues)} issues")
        return 1
    summary = artifact["summary"]
    print(
        "validated order-eight shifted-jet certificate: "
        f"{summary['blocks']} blocks, 0 issues, "
        f"{summary['continuous_curvature_theorems']} continuous theorem, "
        f"largest scaled upper {summary['largest_scaled_upper']}, "
        f"{summary['open_ray_targets']} open ray"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
