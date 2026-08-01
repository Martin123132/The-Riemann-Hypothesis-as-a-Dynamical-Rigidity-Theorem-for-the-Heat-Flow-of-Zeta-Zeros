#!/usr/bin/env python3
"""Validate the upper-compact step-four source-profile pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_upper_profile_lattice_pilot as pilot  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-live-source-rebuilds", action="store_true")
    args = parser.parse_args()
    try:
        stored = json.loads(pilot.DEFAULT_OUT.read_text(encoding="utf-8"))
        profile_names = tuple(stored["parameters"]["evaluated_profiles"])
        rebuilt = pilot.build_artifact(profile_names, generate=False)
    except (KeyError, OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"upper compact profile lattice: source failure: {exc}")
        return 1
    issues = []
    if stored != rebuilt:
        issues.append("stored artifact differs from the live interval rebuild")
    expected_rows = pilot.ROWS_PER_PROFILE * len(profile_names)
    if rebuilt.get("summary", {}).get("exact_source_rows") != expected_rows:
        issues.append("exact source-row count changed")
    if not args.skip_live_source_rebuilds:
        records = pilot.load_source_cache()[:expected_rows]
        expected = pilot.all_tasks()[:expected_rows]
        for index in sorted({0, expected_rows // 2, expected_rows - 1}):
            if pilot.exact_task(expected[index]) != records[index]:
                issues.append(f"live exact source rebuild mismatch at row {index}")
    if issues:
        print(f"upper compact profile lattice: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    summary = rebuilt["summary"]
    print(
        "validated upper compact profile lattice: "
        f"profiles={summary['evaluated_profiles']}, "
        f"passing={summary['passing_profiles']}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
