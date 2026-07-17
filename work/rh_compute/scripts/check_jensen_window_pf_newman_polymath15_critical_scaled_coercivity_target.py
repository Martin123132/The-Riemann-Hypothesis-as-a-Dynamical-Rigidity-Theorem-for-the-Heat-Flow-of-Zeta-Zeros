#!/usr/bin/env python3
"""Validate the critical scaled corrected-coercivity target."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target import (
    DEFAULT_OUT,
    build_exact,
)


EXPECTED_IDS = [
    f"np15csct_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "critical_region"),
        (2, "uncorrected_coercivity"),
        (3, "phase_free_cross"),
        (4, "rs_correction"),
        (5, "corrected_main"),
        (6, "correction_curvature"),
        (7, "refined_remainder"),
        (8, "shrinking_collar"),
        (9, "live_coercivity"),
        (10, "nonpromotion"),
    )
]


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    exact = build_exact()
    for marker in (
        "B^2(X^2+Y^2)",
        "Im(conj(D)*D')",
        "J_(N,t)=P_(N,t)-Q_(N,t)",
        "tilde_epsilon",
    ):
        if not any(marker in str(value) for value in exact.values()):
            issues.append(f"exact target marker missing: {marker}")
    rows = {row.get("id"): row for row in artifact.get("rows", [])}
    if rows.get("np15csct_09_live_coercivity", {}).get("readiness") != "not_ready_to_apply":
        issues.append("live coercivity target was promoted")
    if rows.get("np15csct_10_nonpromotion", {}).get("readiness") != "guard_validated":
        issues.append("nonpromotion guard missing")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the corrected coercivity inequality",
        "critical layer",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman Polymath-15 critical scaled coercivity target: "
        f"{len(artifact['rows'])} rows, 0 issues, 2 exact curvature identities, "
        "1 published endpoint correction, 1 refined remainder, 1 open coercivity target"
    )


if __name__ == "__main__":
    main()
