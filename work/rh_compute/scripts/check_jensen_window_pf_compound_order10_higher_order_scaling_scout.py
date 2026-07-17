#!/usr/bin/env python3
"""Validate the order-ten higher-order scaling scout."""

from __future__ import annotations

import argparse
from decimal import Decimal
from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_higher_order_scaling_scout as target  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    issues: list[str] = []
    if not args.artifact.exists():
        issues.append(f"missing artifact: {args.artifact}")
        artifact = {}
    else:
        artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    if artifact:
        if artifact.get("kind") != (
            "jensen_window_pf_compound_order10_higher_order_scaling_scout"
        ):
            issues.append("bad artifact kind")
        expected_transitions = {"1/8": "1350", "1/4": "1500", "1/2": "1800"}
        if artifact.get("transitions") != expected_transitions:
            issues.append(f"unexpected transitions: {artifact.get('transitions')!r}")
        summary = artifact.get("summary", {})
        expected_summary = {
            "rows": 4,
            "ready_rows": 4,
            "open_rows": 0,
            "selected_anchors": 15,
            "selected_blocks": 45,
            "passing_blocks": 33,
            "failing_blocks": 12,
            "continuum_bridge_theorems": 0,
            "full_kernel_theorems": 0,
        }
        if summary != expected_summary:
            issues.append(f"bad summary: {summary!r}")
        blocks = artifact.get("blocks", [])
        if len(blocks) != 45:
            issues.append(f"bad block count: {len(blocks)}")
        for index, row in enumerate(blocks):
            left = Fraction(row.get("anchor", "0"))
            right = Fraction(row.get("right", "0"))
            width = Fraction(row.get("scout_width", "0"))
            if right - left != width or width not in target.WIDTHS:
                issues.append(f"bad selected block {index}: {left}..{right}")
            try:
                scaled = Decimal(row.get("scaled_curvature_upper", "nan"))
            except Exception as exc:
                issues.append(f"bad scaled upper in block {index}: {exc}")
            else:
                if row.get("passed") is not (scaled < Decimal(5500)):
                    issues.append(f"classification mismatch in block {index}")
        try:
            target.validate_sources()
        except Exception as exc:
            issues.append(f"source validation failed: {exc}")
        if args.rebuild:
            try:
                canonical = target.build_artifact()
            except Exception as exc:
                issues.append(f"canonical rebuild failed: {exc}")
            else:
                if artifact != canonical:
                    issues.append("artifact differs from canonical rebuild")
        rows = artifact.get("rows", [])
        if len(rows) != 4 or any(
            row.get("readiness") != "ready_to_apply" for row in rows
        ):
            issues.append("bad ledger rows")

    if not args.note.exists():
        issues.append(f"missing note: {args.note}")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in (
            "rigorous selected real-interval blocks",
            "width 1/8: first selected passing anchor 1350",
            "width 1/4: first selected passing anchor 1500",
            "width 1/2: first selected passing anchor 1800",
            "not proved threshold locations",
        ):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")
    if issues:
        print(f"order-ten higher-order scaling scout: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-ten higher-order scaling scout: "
        "45 selected blocks, 33 passed, 12 failed, 0 continuum claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
