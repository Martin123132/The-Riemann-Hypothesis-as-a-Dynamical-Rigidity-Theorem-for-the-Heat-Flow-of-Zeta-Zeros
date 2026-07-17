#!/usr/bin/env python3
"""Validate the modular-blend high-frequency cancellation scout."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from jensen_window_pf_newman_theta_modular_blend_high_frequency_scout import (
    DEFAULT_OUT,
    compute_rows,
)


def validate(path: Path, rebuild: bool) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != "jensen_window_pf_newman_theta_modular_blend_high_frequency_scout":
        issues.append("artifact kind mismatch")
    rows = artifact.get("rows", [])
    expected_keys = [(t, x) for t in ("0.0", "0.50000000") for x in (80, 100, 120, 150, 200)]
    actual_keys = [(row.get("t"), row.get("x")) for row in rows]
    if actual_keys != expected_keys:
        issues.append(f"row key/order mismatch: {actual_keys}")
    for row in rows:
        partial = mp.mpf(row["partial_laguerre"])
        full = mp.mpf(row["full_laguerre"])
        ratio = mp.mpf(row["full_to_partial_ratio"])
        if partial <= 0 or full <= 0 or ratio <= 0:
            issues.append(f"nonpositive diagnostic row: {(row['t'], row['x'])}")
        if row["t"] == "0.0":
            if mp.mpf(row["max_relative_xi_jet_delta"]) >= mp.mpf("1e-25"):
                issues.append(f"xi jet cross-check drift: {row['x']}")
            if mp.mpf(row["relative_xi_laguerre_delta"]) >= mp.mpf("1e-25"):
                issues.append(f"xi Laguerre cross-check drift: {row['x']}")
    by_key = {(row["t"], row["x"]): row for row in rows}
    for t in ("0.0", "0.50000000"):
        for x, minimum_digits in ((120, 5), (150, 11), (200, 20)):
            if mp.mpf(by_key[(t, x)]["cancellation_digits"]) <= minimum_digits:
                issues.append(f"high-frequency cancellation guard failed at {(t, x)}")
    convergence = artifact.get("convergence", {})
    if mp.mpf(convergence.get("max_relative_partial_laguerre_delta", "1")) >= mp.mpf("1e-30"):
        issues.append("partial coarse/fine convergence failed")
    if mp.mpf(convergence.get("max_relative_full_laguerre_delta", "1")) >= mp.mpf("1e-25"):
        issues.append("full coarse/fine convergence failed")
    boundary = artifact.get("proof_boundary", "")
    for marker in ("do not certify", "Lambda<=0", "RH"):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    if rebuild and not issues:
        params = artifact["parameters"]
        rebuilt = compute_rows(
            params["fine_nodes_per_panel"],
            params["dps"],
            params["theta_terms"],
        )
        for stored, fresh in zip(rows, rebuilt, strict=True):
            for field in ("partial_laguerre", "full_laguerre"):
                left = mp.mpf(stored[field])
                right = mp.mpf(fresh[field])
                if abs(left - right) / abs(right) >= mp.mpf("1e-40"):
                    issues.append(
                        f"rebuild drift in {field} at {(stored['t'], stored['x'])}"
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
        "validated Newman modular-blend high-frequency scout: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        "2 times, 5 frequencies, 3 fixed blocks, 1 cancellation non-promotion gate"
    )


if __name__ == "__main__":
    main()
