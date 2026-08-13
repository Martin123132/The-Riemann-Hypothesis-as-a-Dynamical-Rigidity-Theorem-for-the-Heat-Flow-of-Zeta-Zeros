#!/usr/bin/env python3
"""Validate exact low-mode projection of the physical Hardy kernels."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from jensen_window_pf_newman_c1_hardy_exact_low_mode_projected_kernel_gate import (
    ACTIVE_LOW_MODES,
    DEFAULT_NOTE,
    DEFAULT_OUT,
    HEIGHT_KEYS,
    KIND,
    build_artifact,
)


EXPECTED_CONCEPTUAL_IDS = [
    f"helpk_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "saved_kernels"),
        (2, "integer_modes"),
        (3, "rational_projection"),
        (4, "exact_annihilation"),
        (5, "fit_replay"),
        (6, "derivative_budget"),
        (7, "scalar_guard"),
        (8, "physical_guard"),
    )
]


def validate(path: Path, note_path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != KIND:
        issues.append("artifact kind mismatch")
    if artifact.get("status") != (
        "exact_low_mode_projection_with_diagnostic_kernel_replay"
    ):
        issues.append("artifact status mismatch")
    conceptual_ids = [
        row.get("id") for row in artifact.get("conceptual_rows", [])
    ]
    if conceptual_ids != EXPECTED_CONCEPTUAL_IDS:
        issues.append(f"conceptual row id/order mismatch: {conceptual_ids}")

    rebuilt = build_artifact()
    if artifact != rebuilt:
        issues.append("artifact differs from deterministic rebuild")

    rows = artifact.get("rows", [])
    if len(rows) != 24:
        issues.append(f"physical row count drifted: {len(rows)}")
    annihilations = 0
    for row in rows:
        key = (row.get("family"), row.get("name"))
        if key not in ACTIVE_LOW_MODES:
            issues.append(f"unregistered physical row: {key}")
            continue
        active = tuple(row.get("active_low_modes", []))
        if active != ACTIVE_LOW_MODES[key]:
            issues.append(f"active-mode registry drifted for {key}")
        inactive = tuple(row.get("exact_annihilated_low_modes", []))
        expected_inactive = tuple(index for index in range(6) if index not in active)
        if inactive != expected_inactive:
            issues.append(f"inactive-mode registry drifted for {key}")
        annihilations += len(inactive)

        coefficients = [
            Fraction(value) for value in row.get("exact_projected_coefficients", [])
        ]
        moments = [Fraction(value) for value in row.get("exact_low_mode_moments", [])]
        if len(coefficients) != 15 or len(moments) != 6:
            issues.append(f"exact vector dimensions drifted for {key}")
            continue
        for mode_index in inactive:
            if moments[mode_index] != 0:
                issues.append(f"inactive mode {mode_index} is nonzero for {key}")
        for value in row.get("projected_coefficients_decimal", []):
            if not math.isfinite(float(value)):
                issues.append(f"nonfinite projected coefficient for {key}")
        if row.get("projected_held_out_relative_error", 1.0) >= 1.0e-9:
            issues.append(f"held-out projection error too large for {key}")
        if len(row.get("derivative_budget_multipliers", [])) != 7:
            issues.append(f"derivative budget dimension drifted for {key}")
        if len(row.get("normalized_derivative_budget_multipliers", [])) != 7:
            issues.append(f"normalized derivative budget dimension drifted for {key}")
        projections = row.get("residual_projections", {})
        if tuple(projections) != HEIGHT_KEYS:
            issues.append(f"residual height registry drifted for {key}")

    summary = artifact.get("summary", {})
    if annihilations != 84 or summary.get("exact_low_mode_annihilations") != 84:
        issues.append("exact annihilation count drifted")
    if summary.get("physical_rows") != 24:
        issues.append("summary physical row count drifted")
    if summary.get("maximum_projected_held_out_relative_error", 1.0) >= 1.0e-9:
        issues.append("summary held-out error threshold failed")
    if summary.get("maximum_coefficient_change_l2", 1.0) >= 1.0e-6:
        issues.append("coefficient projection changed a vector too far")
    if summary.get("increased_projection_rows") != 21:
        issues.append("two-height scalar obstruction drifted")
    derivative_maxima = summary.get("normalized_derivative_budget_maxima", [])
    if [row.get("derivative_order") for row in derivative_maxima] != list(range(7)):
        issues.append("derivative maximum registry drifted")
    expected_ceilings = (0.485, 0.020, 0.0010, 3.2e-5, 6.6e-7, 2.1e-8, 5.3e-8)
    for row, ceiling in zip(derivative_maxima, expected_ceilings, strict=True):
        value = row.get("maximum_normalized_multiplier", float("inf"))
        if not (0.0 <= value < ceiling):
            issues.append(
                f"derivative order {row.get('derivative_order')} exceeds ceiling"
            )

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the external evaluator derivative theorem",
        "interval kernel approximation",
        "physical carrier value",
        "Lambda<=0",
        "RH",
        "prize-level conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = note_path.read_text(encoding="utf-8")
    for marker in (
        "Exact Low-Mode Projected-Kernel Gate",
        "84` exact annihilations",
        "180-digit phase precision",
        "Derivative sensitivity",
        "21/24",
        "does not revive the rejected scalar-extrapolation route",
        "No new `pi` is introduced",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    issues = validate(args.artifact, args.note)
    if issues:
        for issue in issues:
            print(f"ISSUE: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman C1 Hardy exact low-mode projected-kernel gate: "
        "24 rows, 84 exact annihilations, 6 low modes, "
        f"max held-out error {artifact['summary']['maximum_projected_held_out_relative_error']:.12g}, "
        "21/24 projected errors still increase, 0 physical values and 1 open derivative theorem"
    )


if __name__ == "__main__":
    main()
