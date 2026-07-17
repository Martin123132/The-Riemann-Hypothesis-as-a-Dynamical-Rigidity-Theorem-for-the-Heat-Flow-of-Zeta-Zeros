#!/usr/bin/env python3
"""Validate the localized order-nine lower-bridge certificate."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order9_localized_lower_bridge_certificate as bridge  # noqa: E402


@dataclass(frozen=True)
class BridgeIssue:
    section: str
    code: str
    detail: str


def issue(section: str, code: str, detail: object) -> BridgeIssue:
    return BridgeIssue(section, code, str(detail))


def validate_artifact(path: Path) -> tuple[dict, list[BridgeIssue]]:
    if not path.exists():
        return {}, [issue("artifact", "missing-file", path)]
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[BridgeIssue] = []
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order9_localized_lower_bridge_certificate"
    ):
        issues.append(issue("artifact", "bad-kind", artifact.get("kind")))
    if "1250<=t<=5700" not in str(artifact.get("status", "")):
        issues.append(issue("artifact", "bad-status", artifact.get("status")))
    theorem = "w_1''(t)<=4200/t^2 for every real 1250<=t<=5700"
    if artifact.get("theorem") != theorem:
        issues.append(issue("artifact", "bad-theorem", artifact.get("theorem")))

    coverage = artifact.get("coverage", {})
    if coverage.get("t_range") != ["1250", "5700"]:
        issues.append(issue("coverage", "bad-range", coverage.get("t_range")))
    if coverage.get("root_segment_count") != len(
        bridge.deterministic_segments()
    ):
        issues.append(
            issue(
                "coverage",
                "bad-root-count",
                coverage.get("root_segment_count"),
            )
        )
    if coverage.get("all_blocks_passed") is not True:
        issues.append(issue("coverage", "not-all-passed", None))
    try:
        largest = Decimal(
            coverage.get("largest_scaled_curvature_upper", "Infinity")
        )
        weakest = Decimal(
            coverage.get("weakest_common_coordinate_lower", "-Infinity")
        )
    except Exception as exc:
        issues.append(issue("coverage", "bad-decimal", exc))
    else:
        if largest >= Decimal(bridge.CURVATURE_CONSTANT):
            issues.append(issue("coverage", "curvature-failure", largest))
        if weakest <= 0:
            issues.append(issue("coverage", "nonpositive-coordinate", weakest))
    width_counts = coverage.get("width_counts", {})
    try:
        counted_blocks = sum(int(count) for count in width_counts.values())
    except Exception as exc:
        issues.append(issue("coverage", "bad-width-counts", exc))
    else:
        if counted_blocks != coverage.get("accepted_block_count"):
            issues.append(
                issue("coverage", "width-count-mismatch", counted_blocks)
            )
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 4,
        "ready_rows": 4,
        "open_rows": 0,
        "curvature_theorems": 1,
        "root_segments": coverage.get("root_segment_count"),
        "accepted_blocks": coverage.get("accepted_block_count"),
    }
    for key, expected in expected_summary.items():
        if summary.get(key) != expected:
            issues.append(issue("summary", f"bad-{key}", summary.get(key)))
    return artifact, issues


def validate_bound_sources(artifact: dict) -> list[BridgeIssue]:
    issues: list[BridgeIssue] = []
    source_contract = artifact.get("source_contract", {})
    sources = source_contract.get("sources", [])
    expected_roles = {"exact_t_h2_h22", "exact_point_h0_h8"}
    if {source.get("role") for source in sources} != expected_roles:
        issues.append(issue("sources", "bad-source-roles", sources))
    for source in sources:
        role = source.get("role", "unknown")
        for path_key, hash_key in (
            ("cache", "cache_sha256"),
            ("manifest", "manifest_sha256"),
        ):
            try:
                path = bridge.REPO_ROOT / source[path_key]
            except Exception as exc:
                issues.append(issue("sources", f"bad-{role}-{path_key}", exc))
                continue
            if not path.exists():
                issues.append(issue("sources", f"missing-{role}-{path_key}", path))
                continue
            if bridge.sha256(path) != source.get(hash_key):
                issues.append(issue("sources", f"hash-{role}-{path_key}", path))
        cache_path = bridge.REPO_ROOT / str(source.get("cache", ""))
        if cache_path.exists():
            rows = bridge.jsonl_row_count(cache_path)
            if rows != source.get("rows"):
                issues.append(issue("sources", f"rows-{role}", rows))

    for path_key, hash_key in (
        ("run_contract", "run_contract_sha256"),
        ("segment_cache", "segment_cache_sha256"),
    ):
        try:
            path = bridge.REPO_ROOT / source_contract[path_key]
        except Exception as exc:
            issues.append(issue("sources", f"bad-{path_key}", exc))
            continue
        if not path.exists():
            issues.append(issue("sources", f"missing-{path_key}", path))
            continue
        if bridge.sha256(path) != source_contract.get(hash_key):
            issues.append(issue("sources", f"hash-{path_key}", path))
    return issues


def validate_cache_and_canonical(artifact: dict) -> list[BridgeIssue]:
    issues: list[BridgeIssue] = []
    source_contract = artifact.get("source_contract", {})
    try:
        segment_path = bridge.REPO_ROOT / source_contract["segment_cache"]
        contract_path = bridge.REPO_ROOT / source_contract["run_contract"]
        contract = json.loads(contract_path.read_text(encoding="utf-8"))
        tasks = bridge.deterministic_segments()
        segments = bridge.load_segment_cache(segment_path, tasks)
        if len(segments) != len(tasks):
            raise RuntimeError(
                f"segment cache has {len(segments)}/{len(tasks)} rows"
            )
        canonical = bridge.build_artifact(
            segments,
            run_contract=contract,
            run_contract_path=contract_path,
            segment_cache_path=segment_path,
        )
    except Exception as exc:
        issues.append(issue("cache", "exception", exc))
    else:
        if artifact != canonical:
            issues.append(issue("cache", "artifact-drift", bridge.DEFAULT_OUT))
    return issues


def validate_interval_rebuild(
    artifact: dict,
    *,
    workers: int,
) -> list[BridgeIssue]:
    issues: list[BridgeIssue] = []
    source_contract = artifact.get("source_contract", {})
    try:
        source_by_role = {
            source["role"]: source for source in source_contract["sources"]
        }
        h_cache = (
            bridge.REPO_ROOT / source_by_role["exact_t_h2_h22"]["cache"]
        )
        point_cache = (
            bridge.REPO_ROOT / source_by_role["exact_point_h0_h8"]["cache"]
        )
        segment_path = bridge.REPO_ROOT / source_contract["segment_cache"]
        tasks = bridge.deterministic_segments()
        stored = bridge.load_segment_cache(segment_path, tasks)
        if workers <= 1:
            bridge.initialize_worker(str(h_cache), str(point_cache))
            rebuilt = map(bridge.segment_task, tasks)
            executor = None
        else:
            executor = ProcessPoolExecutor(
                max_workers=workers,
                initializer=bridge.initialize_worker,
                initargs=(str(h_cache), str(point_cache)),
            )
            rebuilt = executor.map(bridge.segment_task, tasks, chunksize=1)
        try:
            for index, (actual, expected) in enumerate(zip(rebuilt, stored)):
                if actual != expected:
                    issues.append(
                        issue("rebuild", "segment-drift", index)
                    )
                    break
        finally:
            if executor is not None:
                executor.shutdown(wait=True, cancel_futures=True)
    except Exception as exc:
        issues.append(issue("rebuild", "exception", exc))
    return issues


def validate_note(path: Path) -> list[BridgeIssue]:
    if not path.exists():
        return [issue("note", "missing-file", path)]
    text = path.read_text(encoding="utf-8")
    required = (
        "Status: rigorous first-summand curvature theorem on `1250<=t<=5700`",
        "This is not a proof",
        "H2-H22",
        "H0-H8",
        "`t+-7`",
        "B,J,R,S,T,U,V",
        "w_1''(t)<=4200/t^2 for every real 1250<=t<=5700",
        "above `t=5700`",
        "umbrella certificate",
    )
    issues = []
    for marker in required:
        if marker not in text:
            issues.append(issue("note", "missing-marker", marker))
    lowered = text.lower()
    for forbidden in (
        "therefore rh",
        "we have proved lambda <= 0",
        "order-nine entry is proved",
        "pf-infinity follows",
    ):
        if forbidden in lowered:
            issues.append(issue("note", "forbidden-overclaim", forbidden))
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=bridge.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=bridge.DEFAULT_NOTE)
    parser.add_argument("--rebuild", action="store_true")
    parser.add_argument("--workers", type=int, default=bridge.DEFAULT_WORKERS)
    args = parser.parse_args()

    artifact, issues = validate_artifact(args.artifact)
    issues.extend(validate_note(args.note))
    if artifact:
        issues.extend(validate_bound_sources(artifact))
        issues.extend(validate_cache_and_canonical(artifact))
        if args.rebuild:
            issues.extend(
                validate_interval_rebuild(
                    artifact,
                    workers=max(1, args.workers),
                )
            )
    if issues:
        for finding in issues:
            print(
                f"ORDER9-LOWER-BRIDGE {finding.section} "
                f"[{finding.code}] {finding.detail}"
            )
        print(f"order-nine localized lower bridge: {len(issues)} issues")
        return 1
    coverage = artifact["coverage"]
    rebuild = " with full interval rebuild" if args.rebuild else ""
    print(
        "validated order-nine localized lower bridge"
        f"{rebuild}: {coverage['root_segment_count']} root segments, "
        f"{coverage['accepted_block_count']} accepted blocks, 0 issues, "
        "largest scaled upper "
        f"{coverage['largest_scaled_curvature_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
