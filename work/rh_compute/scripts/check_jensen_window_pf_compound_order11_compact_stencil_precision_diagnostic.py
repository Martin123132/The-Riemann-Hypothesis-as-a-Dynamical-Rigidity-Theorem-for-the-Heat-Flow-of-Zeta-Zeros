#!/usr/bin/env python3
"""Rebuild the order-eleven compact source-precision diagnostic."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_stencil_precision_diagnostic as diagnostic  # noqa: E402


def main() -> int:
    try:
        stored = json.loads(diagnostic.DEFAULT_OUT.read_text(encoding="utf-8"))
        rebuilt = diagnostic.build_artifact()
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"compact stencil precision diagnostic: source failure: {exc}")
        return 1
    issues = []
    if stored != rebuilt:
        issues.append("stored artifact differs from the live rebuild")
    comparison = rebuilt.get("comparison", {})
    if not comparison.get("ordinary_scaled_curvature_upper"):
        issues.append("ordinary comparison bound is missing")
    if not comparison.get("high_precision_scaled_curvature_upper"):
        issues.append("high-precision comparison bound is missing")
    if issues:
        print(f"compact stencil precision diagnostic: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated compact stencil precision diagnostic: "
        f"high_passed={comparison['high_precision_passed_cap']}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
