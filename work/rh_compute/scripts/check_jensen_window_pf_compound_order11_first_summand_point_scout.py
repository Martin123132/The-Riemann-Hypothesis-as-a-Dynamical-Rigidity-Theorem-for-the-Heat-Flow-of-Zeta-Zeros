#!/usr/bin/env python3
"""Validate the order-eleven first-summand point scout."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_compound_order11_first_summand_point_scout as target  # noqa: E402


@dataclass(frozen=True)
class Finding:
    section: str
    issue: str
    detail: str


def finding(section: str, issue: str, detail: object) -> Finding:
    return Finding(section, issue, str(detail))


def validate(artifact_path: Path, note_path: Path) -> list[Finding]:
    findings: list[Finding] = []
    artifact = json.loads(artifact_path.read_text(encoding="utf-8"))
    try:
        rebuilt = target.build_artifact()
    except Exception as exc:
        findings.append(finding("rebuild", "failed", repr(exc)))
    else:
        if rebuilt != artifact:
            findings.append(finding("rebuild", "artifact-mismatch", "rebuilt JSON differs"))
    if artifact.get("kind") != "jensen_window_pf_compound_order11_first_summand_point_scout":
        findings.append(finding("artifact", "bad-kind", artifact.get("kind")))
    if artifact.get("exact") != {
        "generic_stage": "G_p=p*B-Delta^2*l_(p-1); l_p=2*l_(p-1)-l_(p-2)+log(1-exp(-G_p))",
        "eighth_gap": "X=9*B-Delta^2*z",
        "order10_log_coordinate": "y=2*z-w+log(1-exp(-X))",
        "scaled_curvature": "t^2*y_1''(t)",
    }:
        findings.append(finding("artifact", "exact-mismatch", artifact.get("exact")))
    summary = artifact.get("summary", {})
    expected_summary = {
        "rows": 6,
        "ready_rows": 5,
        "open_rows": 1,
        "selected_points": 8,
        "positive_coordinate_balls": 72,
        "point_curvature_enclosures": 8,
        "continuous_curvature_theorems": 0,
        "full_kernel_transfer_theorems": 0,
    }
    if summary != expected_summary:
        findings.append(finding("summary", "mismatch", summary))
    points = artifact.get("diagnostics", {}).get("selected_points", [])
    if [row.get("t") for row in points] != [str(value) for value in target.SELECTED_ANCHORS]:
        findings.append(finding("points", "anchor-mismatch", [row.get("t") for row in points]))
    for row in points:
        if row.get("passed") is not True:
            findings.append(finding("points", "failed-row", row.get("t")))
        if set(row.get("coordinate_lower_bounds", {})) != {"B", "J", "R", "S", "T", "U", "V", "W", "X"}:
            findings.append(finding("points", "coordinate-set-mismatch", row.get("t")))
    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "X=9*B-Delta^2*z",
        "y=2*z-w+log(1-exp(-X))",
        "exact-point Arb enclosures, not a continuous ray theorem",
        "not PF-infinity",
        "RH, or `Lambda<=0`",
    ):
        if marker not in note:
            findings.append(finding("note", "missing-marker", marker))
    return findings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=target.DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=target.DEFAULT_NOTE)
    args = parser.parse_args()
    findings = validate(args.artifact, args.note)
    if findings:
        for row in findings:
            print(f"{row.section}: {row.issue}: {row.detail}")
        print(f"order-eleven first-summand point scout: {len(findings)} issues")
        return 1
    maximum = json.loads(args.artifact.read_text(encoding="utf-8"))["diagnostics"]["maximum_scaled_curvature_upper"]
    print(f"validated order-eleven first-summand point scout: 8 points, 72 positive coordinates, maximum scaled upper {maximum}, 0 issues")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
