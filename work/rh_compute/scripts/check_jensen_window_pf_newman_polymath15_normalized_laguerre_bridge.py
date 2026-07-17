#!/usr/bin/env python3
"""Validate the Polymath-15 normalized Laguerre bridge."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from jensen_window_pf_newman_polymath15_normalized_laguerre_bridge import (
    DEFAULT_OUT,
    SOURCE_URL,
    TIME_VALUES,
    X_VALUES,
    build_exact,
    phase_diagnostics,
)


EXPECTED_IDS = [
    f"np15nlb_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "published_input"),
        (2, "real_axis_collapse"),
        (3, "normalized_curvature"),
        (4, "phase_potential"),
        (5, "cauchy_jet_transfer"),
        (6, "curvature_error"),
        (7, "single_saddle"),
        (8, "phase_scout"),
        (9, "scale_match"),
        (10, "transition_collar"),
        (11, "boundary_layer"),
        (12, "live_handoff"),
    )
]


def validate(path: Path, rebuild: bool) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if (
        artifact.get("kind")
        != "jensen_window_pf_newman_polymath15_normalized_laguerre_bridge"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    if artifact.get("source", {}).get("url") != SOURCE_URL:
        issues.append("published source URL mismatch")
    build_exact()
    exact = artifact.get("exact", {})
    if exact.get("published_input", {}).get("cutoff") != (
        "N=floor(sqrt(x/(4*pi)+t/16))"
    ):
        issues.append("published cutoff mismatch")
    if "t*log(x)=O(1)" not in exact.get("boundary_layer", ""):
        issues.append("boundary-layer reduction missing")

    diagnostics = artifact.get("diagnostics", [])
    expected_keys = [(t, x) for t in TIME_VALUES for x in X_VALUES]
    actual_keys = [(row.get("t"), row.get("x")) for row in diagnostics]
    if actual_keys != expected_keys:
        issues.append(f"diagnostic key/order mismatch: {actual_keys}")
    for row in diagnostics:
        key = (row.get("t"), row.get("x"))
        if mp.mpf(row.get("potential", "nan")) <= 0:
            issues.append(f"nonpositive sampled potential at {key}")
        if mp.mpf(row.get("single_saddle_lower_without_potential", "nan")) <= 0:
            issues.append(f"nonpositive sampled phase floor at {key}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove",
        "uniform finite-main curvature margin",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    if rebuild and not issues:
        fresh = phase_diagnostics(artifact["parameters"]["diagnostic_dps"])
        for stored, rebuilt in zip(diagnostics, fresh, strict=True):
            for field in (
                "beta_first",
                "beta_second",
                "potential",
                "single_saddle_lower_without_potential",
            ):
                left = mp.mpf(stored[field])
                right = mp.mpf(rebuilt[field])
                if abs(left - right) > mp.mpf("1e-32") * max(1, abs(right)):
                    issues.append(
                        f"rebuild drift for {field} at {(stored['t'], stored['x'])}"
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
        "validated Newman Polymath-15 normalized Laguerre bridge: "
        f"{len(artifact['rows'])} rows, 0 issues, 4 exact identities, "
        "2 exact transfer lemmas, 16 phase diagnostics, 2 open tail regimes"
    )


if __name__ == "__main__":
    main()
