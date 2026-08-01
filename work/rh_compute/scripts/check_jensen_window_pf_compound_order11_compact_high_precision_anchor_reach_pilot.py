#!/usr/bin/env python3
"""Rebuild the high-precision compact single-anchor reach pilot."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_high_precision_anchor_reach_pilot as pilot  # noqa: E402


def main() -> int:
    try:
        stored = json.loads(pilot.DEFAULT_OUT.read_text(encoding="utf-8"))
        rebuilt = pilot.build_artifact()
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"compact high-precision anchor reach: source failure: {exc}")
        return 1
    issues = []
    if stored != rebuilt:
        issues.append("stored artifact differs from the live rebuild")
    summary = rebuilt.get("summary", {})
    if summary.get("expansion_rows") != len(pilot.offsets()):
        issues.append("expansion-row count changed")
    if summary.get("quarter_blocks") != 2 * len(pilot.offsets()):
        issues.append("quarter-block count changed")
    if issues:
        print(f"compact high-precision anchor reach: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated compact high-precision anchor reach: "
        f"{summary['passing_expansion_rows']}/{summary['expansion_rows']} expansions, "
        "0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
