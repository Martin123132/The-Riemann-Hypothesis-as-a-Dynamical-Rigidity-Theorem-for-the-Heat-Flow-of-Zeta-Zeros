#!/usr/bin/env python3
"""Validate the order-twelve compact H2-H24 boundary extension."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_compact_h2_h24_boundary_extension as source  # noqa: E402


def main() -> int:
    try:
        artifact = json.loads(source.DEFAULT_OUT.read_text(encoding="utf-8"))
        rebuilt = source.build_artifact()
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"order-twelve compact H boundary: source failure: {exc}")
        return 1
    issues = []
    if artifact != rebuilt:
        issues.append("live two-tile rebuild changed")
    if artifact.get("summary") != {
        "tiles": 2,
        "failed_tiles": 0,
        "covered_tiles": [["5690", "5691"], ["38029", "38030"]],
    }:
        issues.append("boundary summary changed")
    if issues:
        print(f"order-twelve compact H boundary: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print("validated order-twelve compact H boundary: 2 tiles, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
