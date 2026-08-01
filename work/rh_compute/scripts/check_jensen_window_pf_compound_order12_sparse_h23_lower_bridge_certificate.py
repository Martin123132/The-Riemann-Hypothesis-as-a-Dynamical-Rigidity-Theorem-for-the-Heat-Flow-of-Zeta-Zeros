#!/usr/bin/env python3
"""Validate the complete order-twelve sparse-H23 physical lower certificate."""

from __future__ import annotations

import argparse
from decimal import Decimal
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate as certificate  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=certificate.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=certificate.DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []
    try:
        actual = json.loads(args.artifact.read_text(encoding="utf-8"))
        expected = certificate.build_artifact()
    except (OSError, json.JSONDecodeError, RuntimeError, KeyError) as exc:
        print(f"order-twelve lower-bridge certificate: source failure: {exc}")
        return 1
    if actual != expected:
        issues.append("artifact differs from deterministic live-source rebuild")
    summary = actual.get("summary", {})
    if (
        actual.get("kind")
        != "jensen_window_pf_compound_order12_sparse_h23_lower_bridge_certificate"
        or actual.get("theorem") != certificate.THEOREM
        or summary.get("segments") != 263
        or summary.get("required_segments") != 263
        or summary.get("quarter_blocks") != 16788
        or summary.get("required_quarter_blocks") != 16788
        or summary.get("covered_interval") != ["1503", "5700"]
        or summary.get("global_first_summand_theorems") != 0
        or summary.get("rh_claims") != 0
    ):
        issues.append("identity, proof boundary, or complete-coverage summary changed")
    try:
        if Decimal(summary["largest_scaled_curvature_upper"]) >= Decimal(8000):
            issues.append("scaled curvature cap failed")
        if Decimal(summary["smallest_margin_lower"]) <= 0:
            issues.append("nonpositive global margin")
    except (KeyError, ValueError):
        issues.append("numeric summary parse failure")
    if not args.note.exists():
        issues.append("missing note")
    else:
        note = args.note.read_text(encoding="utf-8")
        for marker in ("263/263", "16788/16788", "This is not a proof of RH"):
            if marker not in note:
                issues.append(f"note missing marker: {marker}")
    if issues:
        print(f"order-twelve lower-bridge certificate: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated complete order-twelve lower bridge: 263 segments, "
        "16788 quarter blocks, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
