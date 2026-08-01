#!/usr/bin/env python3
"""Validate the order-eleven compact H2-H24 boundary extension."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension as extension  # noqa: E402


def main() -> int:
    try:
        actual = json.loads(extension.DEFAULT_OUT.read_text(encoding="utf-8"))
        expected = extension.build_artifact()
    except (OSError, json.JSONDecodeError, RuntimeError) as exc:
        print(f"order-eleven compact H boundary extension: source failure: {exc}")
        return 1
    issues = []
    if actual != expected:
        issues.append("artifact differs from two independent live tile rebuilds")
    rows = actual.get("rows", [])
    if (
        len(rows) != 2
        or [row.get("target_t_left") for row in rows] != ["5691", "38028"]
        or [row.get("target_t_right") for row in rows] != ["5692", "38029"]
        or any(row.get("passed") is not True for row in rows)
    ):
        issues.append("boundary tile identity or coverage changed")
    if issues:
        print(f"order-eleven compact H boundary extension: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print("validated order-eleven compact H boundary extension: 2 live rebuilds, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
