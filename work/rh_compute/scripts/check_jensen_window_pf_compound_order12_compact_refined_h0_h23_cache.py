#!/usr/bin/env python3
"""Validate the refined order-twelve compact H0-H23 source cache."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_compact_refined_h0_h23_cache as cache  # noqa: E402


def manifest_issues(
    manifest: dict,
    records: list[dict],
    expected: list[tuple],
    cache_path: Path,
) -> list[str]:
    issues = []
    payload = manifest.get("cache", {})
    latest = manifest.get("latest_run", {})
    complete = len(records) == len(expected)
    next_task = expected[len(records)] if not complete else None
    lower_rows = sum(
        record.get("profile") == cache.base.LOWER_PROFILE_NAME
        for record in records
    )
    upper_rows = sum(
        record.get("profile") == cache.base.UPPER_PROFILE_NAME
        for record in records
    )
    expected_parameters = {
        "start_t": str(cache.START_T),
        "end_t": str(cache.END_T),
        "step_t": str(cache.STEP_T),
        "required_rows": len(expected),
        "profiles": cache.PROFILES,
        "max_moment": cache.MAX_MOMENT,
        "window_y": cache.WINDOW_Y,
        "taylor_order": cache.TAYLOR_ORDER,
    }
    if (
        manifest.get("kind")
        != "jensen_window_pf_compound_order12_compact_refined_h0_h23_cache"
    ):
        issues.append("manifest kind changed")
    if manifest.get("source_contract") != cache.SOURCE_CONTRACT:
        issues.append("source contract changed")
    if manifest.get("parameters") != expected_parameters:
        issues.append("manifest parameters changed")
    if payload.get("path") != cache.relative(cache_path):
        issues.append("cache path changed")
    if payload.get("sha256") != cache.sha256(cache_path):
        issues.append("cache hash changed")
    if payload.get("row_count") != len(records):
        issues.append("manifest row count changed")
    if payload.get("required_row_count") != len(expected):
        issues.append("required row count changed")
    if payload.get("all_rows_passed") is not True:
        issues.append("manifest pass flag changed")
    if payload.get("complete") is not complete:
        issues.append("manifest completion flag changed")
    if payload.get("lower_profile_rows") != lower_rows:
        issues.append("lower-profile row count changed")
    if payload.get("upper_profile_rows") != upper_rows:
        issues.append("upper-profile row count changed")
    if lower_rows + upper_rows != len(records):
        issues.append("cache contains an uncounted profile row")
    expected_next = (
        (next_task[0], str(next_task[1]), next_task[2])
        if next_task is not None
        else (None, None, None)
    )
    actual_next = (
        latest.get("next_index"),
        latest.get("next_target_t"),
        latest.get("next_profile"),
    )
    if actual_next != expected_next:
        issues.append("manifest resume cursor changed")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-complete", action="store_true")
    parser.add_argument("--skip-live-rebuilds", action="store_true")
    parser.add_argument("--rebuild-index", action="append", type=int, default=[])
    args = parser.parse_args()
    try:
        expected = cache.deterministic_tasks()
        records = cache.load_cache(cache.DEFAULT_CACHE, expected)
        manifest = json.loads(cache.DEFAULT_MANIFEST.read_text(encoding="utf-8"))
        issues = manifest_issues(manifest, records, expected, cache.DEFAULT_CACHE)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"refined compact exact source: source failure: {exc}")
        return 1
    complete = len(records) == len(expected)
    if args.require_complete and not complete:
        issues.append(f"cache is incomplete: {len(records)}/{len(expected)} rows")
    rebuild_indices = set(args.rebuild_index)
    if records and not args.skip_live_rebuilds:
        rebuild_indices.update({0, len(records) // 2, len(records) - 1})
    for index in sorted(rebuild_indices):
        if not 0 <= index < len(records):
            issues.append(f"live rebuild index is outside cache: {index}")
        elif cache.exact_task(expected[index]) != records[index]:
            issues.append(f"live rebuild mismatch at row {index}")
    if issues:
        print(f"refined compact exact source: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated refined compact exact source: "
        f"{len(records)}/{len(expected)} rows, complete={complete}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
