#!/usr/bin/env python3
"""Validate the order-ten seventh-nested first-summand point scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order10_first_summand_point_scout as target  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    args = parser.parse_args()
    issues: list[str] = []

    if not args.artifact.exists():
        issues.append(f"missing artifact: {args.artifact}")
    else:
        artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
        expected = target.build_artifact()
        if artifact != expected:
            issues.append("artifact does not match a fresh deterministic build")
        summary = artifact.get("summary", {})
        expected_summary = {
            "rows": 5,
            "ready_rows": 4,
            "open_rows": 1,
            "selected_points": 6,
            "positive_coordinate_balls": 48,
            "point_curvature_cap_certificates": 6,
            "continuous_curvature_theorems": 0,
            "full_kernel_transfer_theorems": 0,
        }
        if summary != expected_summary:
            issues.append(f"bad summary: {summary!r}")
        selected = artifact.get("diagnostics", {}).get("selected_points", [])
        if [row.get("t") for row in selected] != [
            "1251",
            "1300",
            "1500",
            "2000",
            "3000",
            "4000",
        ]:
            issues.append("selected-point grid changed")
        for row in selected:
            if row.get("passed") is not True:
                issues.append(f"selected point did not pass: {row.get('t')}")
            if set(row.get("coordinate_balls", {})) != {
                "B",
                "J",
                "R",
                "S",
                "T",
                "U",
                "V",
                "W",
            }:
                issues.append(f"coordinate set changed at t={row.get('t')}")
        rows = artifact.get("rows", [])
        open_rows = [row for row in rows if row.get("readiness") == "not_ready_to_apply"]
        if len(open_rows) != 1 or open_rows[0].get("id") != "co10fsp_05_continuum_handoff":
            issues.append("unexpected continuum-handoff contract")

    if not args.note.exists():
        issues.append(f"missing note: {args.note}")
    else:
        note = args.note.read_text(encoding="utf-8")
        for required in (
            "t^2*z_1''(t)<250<5500",
            "isolated points do not control a real",
            "higher-summand",
            "not a",
        ):
            if required not in note:
                issues.append(f"note missing required text: {required!r}")

    if issues:
        print(f"order-ten first-summand point scout: {len(issues)} issues")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated order-ten first-summand point scout: "
        "6 selected points, 48 positive coordinates, 6 point curvature caps, "
        "0 continuum theorems, 1 open handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
