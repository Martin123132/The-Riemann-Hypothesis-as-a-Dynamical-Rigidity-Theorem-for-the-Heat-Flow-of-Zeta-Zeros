#!/usr/bin/env python3
"""Validate the high-precision order-eleven compact stencil pilot."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_high_precision_stencil_pilot as pilot  # noqa: E402


def main() -> int:
    try:
        artifact = json.loads(pilot.DEFAULT_OUT.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"high-precision compact stencil: source failure: {exc}")
        return 1
    rows = artifact.get("rows", [])
    issues = []
    if len(rows) != 19 or any(row.get("passed") is not True for row in rows):
        issues.append("stencil row count or pass flags changed")
    if [row.get("target_t") for row in rows] != [str(t) for t in pilot.TARGETS]:
        issues.append("stencil target order changed")
    for index in (0, 9, 18):
        if pilot.stress_task(pilot.tasks()[index]) != rows[index]:
            issues.append(f"live high-precision rebuild mismatch at row {index}")
    if issues:
        print(f"high-precision compact stencil: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print("validated high-precision compact stencil: 19 rows, 3 live rebuilds, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
