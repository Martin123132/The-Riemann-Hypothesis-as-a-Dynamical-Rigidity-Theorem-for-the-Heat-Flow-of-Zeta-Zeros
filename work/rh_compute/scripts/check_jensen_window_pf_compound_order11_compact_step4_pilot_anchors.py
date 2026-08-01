#!/usr/bin/env python3
"""Validate the compact step-four pilot H0-H23 anchors."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_step4_pilot_anchors as pilot  # noqa: E402


def main() -> int:
    try:
        actual = json.loads(pilot.DEFAULT_OUT.read_text(encoding="utf-8"))
        expected = pilot.build_artifact()
    except (OSError, json.JSONDecodeError, RuntimeError) as exc:
        print(f"compact step-four pilot anchors: source failure: {exc}")
        return 1
    issues = []
    if actual != expected:
        issues.append("artifact differs from four live quadrature rebuilds")
    if actual.get("summary", {}).get("targets") != ["20000", "20008", "38016", "38024"]:
        issues.append("pilot targets changed")
    if any(row.get("passed") is not True for row in actual.get("rows", [])):
        issues.append("a pilot row is not passing")
    if issues:
        print(f"compact step-four pilot anchors: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print("validated compact step-four pilot anchors: 4 live rebuilds, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
