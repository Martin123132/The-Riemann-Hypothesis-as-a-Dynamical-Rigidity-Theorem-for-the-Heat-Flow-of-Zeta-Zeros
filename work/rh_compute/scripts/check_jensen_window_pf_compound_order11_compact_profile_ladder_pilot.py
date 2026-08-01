#!/usr/bin/env python3
"""Rebuild the compact exact-source profile ladder pilot."""

from __future__ import annotations

import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_compact_profile_ladder_pilot as pilot  # noqa: E402


def main() -> int:
    try:
        stored = json.loads(pilot.DEFAULT_OUT.read_text(encoding="utf-8"))
        rebuilt = pilot.build_artifact()
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as exc:
        print(f"compact profile ladder: source failure: {exc}")
        return 1
    issues = []
    if stored != rebuilt:
        issues.append("stored artifact differs from the live rebuild")
    summary = rebuilt.get("summary", {})
    if summary.get("profiles") != len(pilot.PROFILES):
        issues.append("profile count changed")
    if set(summary.get("scaled_curvature_uppers", {})) != set(pilot.PROFILES):
        issues.append("profile summary keys changed")
    if issues:
        print(f"compact profile ladder: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated compact profile ladder: "
        f"passing={summary['passing_profiles']}, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
