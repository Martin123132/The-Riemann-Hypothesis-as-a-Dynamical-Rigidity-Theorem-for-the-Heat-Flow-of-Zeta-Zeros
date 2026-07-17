#!/usr/bin/env python3
"""Validate the order-ten localized lower-bridge certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_localized_lower_bridge_certificate as target  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    if not args.artifact.exists():
        issues.append(f"missing artifact: {args.artifact}")
        artifact = {}
    else:
        artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    if artifact:
        if artifact.get("kind") != (
            "jensen_window_pf_compound_order10_localized_lower_bridge_certificate"
        ):
            issues.append("bad artifact kind")
        if artifact.get("theorem") != (
            "z_1''(t)<=4200/t^2 for every real 1251<=t<=5700"
        ):
            issues.append("lower-bridge theorem changed")
        expected_summary = {
            "rows": 5,
            "ready_rows": 5,
            "open_rows": 0,
            "segment_rows": 279,
            "accepted_blocks": 9996,
            "near_blocks": 996,
            "middle_blocks": 1200,
            "far_blocks": 7800,
            "lower_bridge_curvature_theorems": 1,
            "full_half_line_theorems": 0,
            "full_kernel_theorems": 0,
        }
        if artifact.get("summary") != expected_summary:
            issues.append(f"bad summary: {artifact.get('summary')!r}")
        composition = artifact.get("composition", {})
        maximum = composition.get("global_maximum", {})
        minimum = composition.get("global_minimum_margin", {})
        try:
            scaled = Decimal(maximum.get("scaled_curvature_upper", "nan"))
            margin = Decimal(minimum.get("curvature_margin_lower", "nan"))
        except Exception as exc:
            issues.append(f"bad global diagnostics: {exc}")
        else:
            if not (
                Decimal(0)
                < scaled
                < Decimal(target.THEOREM_CURVATURE_CONSTANT)
            ):
                issues.append(f"bad global scaled upper: {scaled}")
            if margin <= 0:
                issues.append(f"bad global margin: {margin}")
        source = artifact.get("source_contract", {})
        segment_record = source.get("segment_cache", {})
        segment_path = target.REPO_ROOT / segment_record.get("path", "")
        if (
            not segment_path.exists()
            or target.segments.sha256(segment_path) != segment_record.get("sha256")
            or segment_record.get("row_count") != 279
        ):
            issues.append("segment cache hash or row count changed")
        try:
            canonical = target.build_artifact()
        except Exception as exc:
            issues.append(f"canonical rebuild failed: {exc}")
        else:
            if artifact != canonical:
                issues.append("artifact differs from canonical composition")
        rows = artifact.get("rows", [])
        if len(rows) != 5 or any(
            row.get("readiness") != "ready_to_apply" for row in rows
        ):
            issues.append("bad theorem-ledger rows")

    if not args.note.exists():
        issues.append(f"missing note: {args.note}")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in (
            "rigorous first-summand continuous curvature theorem",
            "z_1''(t)<=4200/t^2 for every real 1251<=t<=5700",
            "9,996",
            "not yet a full-kernel or RH theorem",
            "first-summand-to-full-kernel transfer",
        ):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")
    if issues:
        print(f"order-ten localized lower bridge: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-ten localized lower bridge: "
        "279 segments, 9,996 contiguous blocks, 0 issues, "
        "1 first-summand theorem, 0 full-kernel claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
