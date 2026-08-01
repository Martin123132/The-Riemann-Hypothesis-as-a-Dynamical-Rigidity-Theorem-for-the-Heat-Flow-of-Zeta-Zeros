#!/usr/bin/env python3
"""Validate the compact unit-grid pilot anchor artifact."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_unit_pilot_anchors as pilot  # noqa: E402
import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as source  # noqa: E402


def main() -> int:
    try:
        artifact = json.loads(pilot.DEFAULT_OUT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"compact unit pilot anchors: source failure: {exc}")
        return 1
    rows = artifact.get("rows", [])
    issues = []
    if len(rows) != len(pilot.TARGETS) or any(row.get("passed") is not True for row in rows):
        issues.append("unit-pilot row count or pass flags changed")
    if [row.get("target_t") for row in rows] != [str(t) for t in pilot.TARGETS]:
        issues.append("unit-pilot target order changed")
    for index in (0, len(pilot.TARGETS) // 2, len(pilot.TARGETS) - 1):
        if source.exact_task(pilot.tasks()[index]) != {
            **rows[index],
            "kind": "order10_compact_sparse_point_h0_h23_jet",
        }:
            issues.append(f"live quadrature rebuild mismatch at pilot row {index}")
    if issues:
        print(f"compact unit pilot anchors: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print("validated compact unit-grid pilot anchors: 14 rows, 3 live rebuilds, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
