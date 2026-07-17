#!/usr/bin/env python3
"""Validate the scaled Newman boundary-layer asymptotic theorem."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem import (
    DEFAULT_OUT,
    build_exact,
)


EXPECTED_IDS = [
    f"np15sblat_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "scaled_regime"),
        (2, "power_envelope"),
        (3, "zeta_moments"),
        (4, "euler_floor"),
        (5, "c2_convergence"),
        (6, "phase_amplitude_identity"),
        (7, "main_asymptotic_curvature"),
        (8, "remainder_transfer"),
        (9, "exact_h_asymptotic"),
        (10, "live_handoff"),
    )
]


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem"
    ):
        issues.append("artifact kind mismatch")
    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")
    exact = build_exact()
    if exact.get("critical_layer") != "0<t*log(x/(4*pi))<=4+o(1)":
        issues.append("critical scaled layer mismatch")
    if "p_epsilon=1+epsilon/16" not in exact.get("coefficient_envelope", ""):
        issues.append("summable envelope mismatch")
    if "1/zeta(3/2+epsilon/8)" not in exact.get("zeta_floor", ""):
        issues.append("Euler-product floor mismatch")
    if "For every epsilon>0" not in exact.get("theorem", ""):
        issues.append("epsilon quantifier missing")
    if "c_*=4911678521/1933561194" not in exact.get(
        "current_global_boundary", ""
    ):
        issues.append("superseding global boundary missing")
    rows = {row.get("id"): row for row in artifact.get("rows", [])}
    for row_id in (
        "np15sblat_07_main_asymptotic_curvature",
        "np15sblat_08_remainder_transfer",
        "np15sblat_09_exact_h_asymptotic",
    ):
        if rows.get(row_id, {}).get("readiness") != "ready_to_apply":
            issues.append(f"proved asymptotic row not ready: {row_id}")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "For every fixed epsilon>0",
        "does not supply a practical numerical L_epsilon",
        "t*log(x/(4*pi))<=4+o(1)",
        "historical boundary of the absolute-value method",
        "c_*=4911678521/1933561194",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    sources = set(artifact.get("sources", []))
    if (
        "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md"
        not in sources
    ):
        issues.append("superseding oscillatory theorem source missing")
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
        "validated Newman Polymath-15 scaled boundary-layer asymptotic theorem: "
        f"{len(artifact['rows'])} rows, 0 issues, 1 summable envelope, "
        "1 Euler-product floor, 1 exact phase identity, 1 exact-H asymptotic theorem"
    )


if __name__ == "__main__":
    main()
