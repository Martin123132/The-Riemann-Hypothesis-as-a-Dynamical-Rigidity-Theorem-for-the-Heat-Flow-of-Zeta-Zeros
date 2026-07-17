#!/usr/bin/env python3
"""Validate the Newman oscillatory zeta handoff theorem artifact."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem import (
    DEFAULT_OUT,
    build_exact,
    exponent_pair,
    phase_exponent,
    required_scaled_time,
    transition_radius,
)


EXPECTED_IDS = [
    f"np15ozht_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "scaled_geometry"),
        (2, "exponent_pair_input"),
        (3, "exact_envelope"),
        (4, "low_absolute_range"),
        (5, "high_oscillatory_range"),
        (6, "zeta_jet_handoff"),
        (7, "euler_floor"),
        (8, "phase_curvature_identity"),
        (9, "main_curvature"),
        (10, "remainder_transfer"),
        (11, "exact_h_theorem"),
        (12, "live_handoff"),
    )
]
EXPECTED_THRESHOLD = Fraction(4911678521, 1933561194)
EXPECTED_RADIUS = Fraction(125662, 155153)


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != (
        "jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem"
    ):
        issues.append("artifact kind mismatch")

    ids = [row.get("id") for row in artifact.get("rows", [])]
    if ids != EXPECTED_IDS:
        issues.append(f"row id/order mismatch: {ids}")

    exact = build_exact()
    stored_exact = artifact.get("exact", {})
    if stored_exact != exact:
        issues.append("stored exact block differs from independent rebuild")

    envelope = exact["envelope"]
    threshold = Fraction(envelope["critical_scaled_time"]["exact"])
    radius = Fraction(envelope["critical_radius"]["exact"])
    if threshold != EXPECTED_THRESHOLD:
        issues.append(f"critical threshold mismatch: {threshold}")
    if radius != EXPECTED_RADIUS:
        issues.append(f"critical radius mismatch: {radius}")
    if transition_radius(2, 1) != EXPECTED_RADIUS:
        issues.append("pair 2/1 transition identity failed")
    if required_scaled_time(2, EXPECTED_RADIUS) != EXPECTED_THRESHOLD:
        issues.append("pair 2/1 threshold identity failed")
    if phase_exponent(2, EXPECTED_RADIUS) != phase_exponent(1, EXPECTED_RADIUS):
        issues.append("pair 2/1 phase lines do not meet")

    transition_rows = envelope.get("transition_rows", [])
    if len(transition_rows) != 10:
        issues.append(f"transition row count mismatch: {len(transition_rows)}")
    maxima = [row for row in transition_rows if row.get("is_critical_maximum")]
    if len(maxima) != 1:
        issues.append(f"critical transition count mismatch: {len(maxima)}")
    elif (maxima[0].get("upper_index"), maxima[0].get("lower_index")) != (2, 1):
        issues.append("critical transition is not pair 2 to pair 1")

    interval_rows = envelope.get("interval_rows", [])
    if len(interval_rows) != 11:
        issues.append(f"interval row count mismatch: {len(interval_rows)}")
    if [row.get("index") for row in interval_rows] != list(range(10, -1, -1)):
        issues.append("active exponent-pair order mismatch")

    endpoint_values: list[Fraction] = []
    for row in interval_rows:
        index = row.get("index")
        if not isinstance(index, int):
            issues.append("noninteger exponent-pair index")
            continue
        kappa, lam = exponent_pair(index)
        if row.get("kappa", {}).get("exact") != f"{kappa.numerator}/{kappa.denominator}":
            issues.append(f"kappa serialization mismatch for pair {index}")
        if row.get("lambda", {}).get("exact") != f"{lam.numerator}/{lam.denominator}":
            issues.append(f"lambda serialization mismatch for pair {index}")
        left = Fraction(row["left_radius"]["exact"])
        right = Fraction(row["right_radius"]["exact"])
        if not (0 < left < right <= 1):
            issues.append(f"invalid active interval for pair {index}")
        midpoint = (left + right) / 2
        for test_radius in (left, midpoint, right):
            chosen = phase_exponent(index, test_radius)
            finite_minimum = min(
                phase_exponent(other, test_radius) for other in range(11)
            )
            if chosen != finite_minimum:
                issues.append(
                    f"pair {index} is not the finite-envelope minimum at "
                    f"r={test_radius}"
                )
        derivative_offset = lam - kappa - Fraction(1, 2)
        if derivative_offset < 0:
            issues.append(f"negative endpoint derivative offset for pair {index}")
        left_required = required_scaled_time(index, left)
        right_required = required_scaled_time(index, right)
        endpoint_values.extend((left_required, right_required))
        if row["left_required_scaled_time"]["exact"] != (
            f"{left_required.numerator}/{left_required.denominator}"
        ):
            issues.append(f"left threshold mismatch for pair {index}")
        if row["right_required_scaled_time"]["exact"] != (
            f"{right_required.numerator}/{right_required.denominator}"
        ):
            issues.append(f"right threshold mismatch for pair {index}")

    if endpoint_values and max(endpoint_values) != EXPECTED_THRESHOLD:
        issues.append("finite-envelope endpoint maximum mismatch")

    rows = {row.get("id"): row for row in artifact.get("rows", [])}
    for row_id in (
        "np15ozht_03_exact_envelope",
        "np15ozht_04_low_absolute_range",
        "np15ozht_05_high_oscillatory_range",
        "np15ozht_06_zeta_jet_handoff",
        "np15ozht_09_main_curvature",
        "np15ozht_10_remainder_transfer",
        "np15ozht_11_exact_h_theorem",
    ):
        if rows.get(row_id, {}).get("readiness") != "ready_to_apply":
            issues.append(f"proved theorem row not ready: {row_id}")
    if rows.get("np15ozht_12_live_handoff", {}).get("readiness") != "not_ready_to_apply":
        issues.append("open handoff was promoted")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "For every fixed epsilon>0",
        "4911678521/1933561194+epsilon",
        "does not provide a practical L_epsilon",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")
    if "not a proof of Lambda<=0 or RH" not in artifact.get("status", ""):
        issues.append("status nonpromotion marker missing")
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
        "validated Newman Polymath-15 oscillatory zeta handoff theorem: "
        f"{len(artifact['rows'])} rows, 0 issues, 11 exponent pairs, "
        "10 exact transitions, threshold 4911678521/1933561194, "
        "1 zeta-jet handoff, 1 exact-H asymptotic theorem"
    )


if __name__ == "__main__":
    main()
