#!/usr/bin/env python3
"""Validate the corrected critical-scaled coercivity scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout import (
    C_VALUES,
    DEFAULT_OUT,
    X_VALUES,
    compute_rows,
)


def validate(path: Path, rebuild: bool) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_critical_scaled_coercivity_scout"
    ):
        issues.append("artifact kind mismatch")
    rows = artifact.get("rows", [])
    expected_keys = [(x, c) for x in X_VALUES for c in C_VALUES]
    actual_keys = [(row.get("x"), row.get("c")) for row in rows]
    if actual_keys != expected_keys:
        issues.append(f"row key/order mismatch: {actual_keys}")
    if any(mp.mpf(row.get("corrected_curvature", "nan")) <= 0 for row in rows):
        issues.append("nonpositive corrected diagnostic")
    exact_rows = [row for row in rows if row.get("c") == "0"]
    if len(exact_rows) != len(X_VALUES):
        issues.append("wrong exact xi row count")
    for row in exact_rows:
        if mp.mpf(row.get("relative_corrected_to_exact_delta", "1")) >= mp.mpf("1e-3"):
            issues.append(f"exact xi discrepancy too large at x={row.get('x')}")
    convergence = artifact.get("convergence", {})
    if convergence.get("stable_corrected_signs") is not True:
        issues.append("coarse/fine corrected sign instability")
    if mp.mpf(convergence.get("max_relative_corrected_curvature_delta", "1")) >= mp.mpf("1e-30"):
        issues.append("coarse/fine corrected curvature drift")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "point diagnostics",
        "do not certify intervals",
        "corrected coercivity target",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    if rebuild and not issues:
        fresh = compute_rows(artifact["parameters"]["fine_dps"])
        for stored, rebuilt in zip(rows, fresh, strict=True):
            left = mp.mpf(stored["corrected_curvature"])
            right = mp.mpf(rebuilt["corrected_curvature"])
            if abs(left - right) / abs(right) >= mp.mpf("1e-30"):
                issues.append(
                    f"rebuild corrected curvature drift at {(stored['x'], stored['c'])}"
                )
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
        "validated Newman Polymath-15 critical scaled coercivity scout: "
        f"{len(artifact['rows'])} rows, 0 issues, 24 positive corrected diagnostics, "
        "4 exact xi cross-checks, 1 nonpromotion boundary"
    )


if __name__ == "__main__":
    main()
