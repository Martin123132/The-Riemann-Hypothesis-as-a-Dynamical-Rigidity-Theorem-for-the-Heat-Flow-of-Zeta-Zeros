#!/usr/bin/env python3
"""Validate the rigorous compact endpoint lattice pilot."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_endpoint_lattice_pilot as pilot  # noqa: E402


def main() -> int:
    try:
        stored = json.loads(pilot.DEFAULT_OUT.read_text(encoding="utf-8"))
        rebuilt = pilot.build_artifact(generate=False)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"compact endpoint lattice: source failure: {exc}")
        return 1
    issues = []
    if stored != rebuilt:
        issues.append("stored artifact differs from the live interval rebuild")
    if rebuilt.get("summary", {}).get("quarter_blocks") != 18:
        issues.append("endpoint quarter-block count changed")
    records = pilot.load_source_cache()
    expected = pilot.tasks()
    for index in (0, len(expected) // 2, len(expected) - 1):
        if pilot.exact_task(expected[index]) != records[index]:
            issues.append(f"live exact source rebuild mismatch at row {index}")
    if issues:
        print(f"compact endpoint lattice: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    summary = rebuilt["summary"]
    print(
        "validated compact endpoint lattice: "
        f"passed={summary['all_endpoints_passed']}, "
        f"blocks={summary['quarter_blocks']}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
