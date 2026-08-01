#!/usr/bin/env python3
"""Validate the exact order-eleven compact H0-H23 right extension."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_sparse_h0_h23_right_extension as extension  # noqa: E402


def main() -> int:
    try:
        actual = json.loads(extension.DEFAULT_OUT.read_text(encoding="utf-8"))
        expected = extension.build_artifact()
    except (OSError, json.JSONDecodeError, RuntimeError) as exc:
        print(f"order-eleven compact right extension: source failure: {exc}")
        return 1
    row = actual.get("row", {})
    issues = []
    if actual != expected:
        issues.append("artifact differs from independent live quadrature rebuild")
    if (
        row.get("target_t") != "38028"
        or row.get("profile") != "far"
        or row.get("passed") is not True
        or set(row.get("h_derivatives", {})) != {str(order) for order in range(24)}
    ):
        issues.append("right-anchor identity or H0-H23 payload changed")
    if issues:
        print(f"order-eleven compact right extension: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print("validated order-eleven compact H0-H23 right extension: 1 live rebuild, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
