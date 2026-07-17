#!/usr/bin/env python3
"""Validate the dominant-saddle arithmetic ray certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate import (
    DEFAULT_OUT,
    arb_moment_data,
    build_exact,
    rational_budget,
)


EXPECTED_IDS = [
    f"np15dsarc_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "alpha_geometry"),
        (2, "ray_geometry"),
        (3, "coefficient_envelope"),
        (4, "arb_moments"),
        (5, "tail_jets"),
        (6, "tail_budget"),
        (7, "phase_floor"),
        (8, "arithmetic_ray"),
        (9, "remainder_gap"),
        (10, "live_handoff"),
    )
]


def validate(path: Path, rebuild: bool) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_dominant_saddle_arithmetic_ray_certificate"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    build_exact()
    if artifact.get("parameters", {}).get("p_exponent") != "349/100":
        issues.append("power exponent mismatch")
    if artifact.get("parameters", {}).get("tL_min") != "24/1":
        issues.append("dominant-ray threshold mismatch")
    budget = artifact.get("rational_budget", {})
    curvature = budget.get("curvature_budget", {})
    if curvature.get("phase_floor_over_L2") != "23/100":
        issues.append("phase floor mismatch")
    if curvature.get("error_over_L2") != "77/500":
        issues.append("error cap mismatch")
    if curvature.get("strict_margin_over_L2") != "19/250":
        issues.append("strict margin mismatch")
    moments = artifact.get("arb", {}).get("moments", {})
    expected_moment_caps = {
        "S0": "0.1280000000000000000000000000000000000000000000000000000000000000000000",
        "S1": "0.1150000000000000000000000000000000000000000000000000000000000000000000",
        "S2": "0.1210000000000000000000000000000000000000000000000000000000000000000000",
    }
    for key, prefix in expected_moment_caps.items():
        if not moments.get(key, {}).get("strict_upper", "").startswith(
            f"[{prefix}"
        ):
            issues.append(f"Arb moment cap mismatch for {key}")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not",
        "analytic remainder",
        "strict Laguerre positivity",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    if rebuild and not issues:
        fresh_arb = arb_moment_data()
        fresh_budget = rational_budget()
        if fresh_budget != artifact.get("rational_budget"):
            issues.append("rational budget rebuild mismatch")
        for key in ("S0", "S1", "S2"):
            if fresh_arb["moments"][key] != artifact["arb"]["moments"][key]:
                issues.append(f"Arb moment rebuild mismatch for {key}")
        if fresh_arb["x_at_L48_ball"] != artifact["arb"].get("x_at_L48_ball"):
            issues.append("L=48 scale rebuild mismatch")
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
        "validated Newman Polymath-15 dominant-saddle arithmetic ray certificate: "
        f"{len(artifact['rows'])} rows, 0 issues, 3 Arb moment bounds, "
        "3 rational tail budgets, 1 strict finite-main ray theorem, 2 open transfer rows"
    )


if __name__ == "__main__":
    main()
