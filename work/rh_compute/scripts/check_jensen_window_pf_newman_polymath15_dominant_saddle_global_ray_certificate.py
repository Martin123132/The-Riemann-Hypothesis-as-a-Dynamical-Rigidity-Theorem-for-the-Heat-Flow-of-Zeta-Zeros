#!/usr/bin/env python3
"""Validate the complete dominant-saddle exact-H ray certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate import (
    DEFAULT_OUT,
    arb_transition_data,
    build_exact,
)


EXPECTED_IDS = [
    f"np15dsgrc_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "boundary_spacing"),
        (2, "fixed_center_cutoff"),
        (3, "adjacent_block"),
        (4, "arb_jump"),
        (5, "global_collar"),
        (6, "global_exact_h_ray"),
        (7, "cutoff_gap_closed"),
        (8, "live_handoff"),
    )
]


def validate(path: Path, rebuild: bool) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    build_exact()
    if artifact.get("parameters", {}).get("tL_min") != "25/1":
        issues.append("global ray threshold mismatch")
    if artifact.get("arb", {}).get("raw_adjacent_jump_lt") != "10^-30":
        issues.append("adjacent jump cap mismatch")
    if artifact.get("arb", {}).get("global_collar_lt") != "1/8000":
        issues.append("global collar cap mismatch")
    curvature = artifact.get("rational_budget", {}).get("curvature", {})
    if curvature.get("full_margin_over_L2") != "3/40":
        issues.append("exact-H margin mismatch")
    rows = {row.get("id"): row for row in artifact.get("rows", [])}
    if rows.get("np15dsgrc_07_cutoff_gap_closed", {}).get("readiness") != "closed":
        issues.append("cutoff closure row not closed")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "no cutoff-distance restriction",
        "does not cover",
        "0<t*log(x/(4*pi))<25",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    if rebuild and not issues:
        fresh = arb_transition_data()
        for key in (
            "adjacent_index_lower_ball",
            "cutoff_spacing_ball",
            "raw_adjacent_jump_ball",
            "global_collar_ball",
        ):
            if fresh[key] != artifact["arb"].get(key):
                issues.append(f"Arb rebuild mismatch for {key}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--rebuild", action="store_true")
    args = parser.parse_args()
    issues = validate(args.artifact, args.rebuild)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman Polymath-15 dominant-saddle global ray certificate: "
        f"{len(artifact['rows'])} rows, 0 issues, 1 adjacent-cutoff Arb bound, "
        "1 global collar, 1 strict exact-H ray theorem, 1 remaining boundary layer"
    )


if __name__ == "__main__":
    main()
