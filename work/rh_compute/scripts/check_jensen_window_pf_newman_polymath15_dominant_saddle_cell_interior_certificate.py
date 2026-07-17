#!/usr/bin/env python3
"""Validate the dominant-saddle exact-H cell-interior certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate import (
    DEFAULT_OUT,
    arb_collar_data,
    build_exact,
    rational_budget,
)


EXPECTED_IDS = [
    f"np15dscic_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "fixed_n_remainder"),
        (2, "schwarz_collar"),
        (3, "scalar_collar"),
        (4, "cauchy_raw_jets"),
        (5, "normalized_jets"),
        (6, "main_jets"),
        (7, "remainder_curvature"),
        (8, "exact_h_cell_theorem"),
        (9, "cutoff_collar_gap"),
        (10, "boundary_layer_gap"),
    )
]


def validate(path: Path, rebuild: bool) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_dominant_saddle_cell_interior_certificate"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    build_exact()
    params = artifact.get("parameters", {})
    expected_params = {
        "rho": "1/4",
        "cell_distance": "1/2",
        "L_min": "50/1",
        "tL_min": "25/1",
        "delta": "1/8000",
    }
    for key, expected in expected_params.items():
        if params.get(key) != expected:
            issues.append(f"parameter mismatch for {key}")
    curvature = artifact.get("rational_budget", {}).get("curvature", {})
    expected_curvature = {
        "finite_main_margin_over_L2": "19/250",
        "remainder_error_over_L2": "1/1000",
        "full_margin_over_L2": "3/40",
    }
    for key, expected in expected_curvature.items():
        if curvature.get(key) != expected:
            issues.append(f"curvature budget mismatch for {key}")
    arb_data = artifact.get("arb", {})
    if arb_data.get("combined_raw_collar_lt") != "1/8000":
        issues.append("raw collar cap mismatch")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not cover",
        "transition collars",
        "t*log(x/(4*pi))<25",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    if rebuild and not issues:
        fresh_arb = arb_collar_data()
        fresh_budget = rational_budget()
        if fresh_budget != artifact.get("rational_budget"):
            issues.append("rational budget rebuild mismatch")
        for key in (
            "positive_eC_correction_ball",
            "eAB_argument_ball",
            "raw_eAB_over_center_amplitude_ball",
            "raw_eC_over_center_amplitude_ball",
            "combined_raw_collar_ball",
        ):
            if fresh_arb[key] != artifact["arb"].get(key):
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
        "validated Newman Polymath-15 dominant-saddle cell-interior certificate: "
        f"{len(artifact['rows'])} rows, 0 issues, 3 Arb collar bounds, "
        "3 Cauchy jet bounds, 1 strict exact-H theorem, 2 open regimes"
    )


if __name__ == "__main__":
    main()
