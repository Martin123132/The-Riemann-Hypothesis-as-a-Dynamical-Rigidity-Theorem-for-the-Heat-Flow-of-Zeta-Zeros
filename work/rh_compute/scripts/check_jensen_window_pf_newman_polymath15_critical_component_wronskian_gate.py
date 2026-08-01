#!/usr/bin/env python3
"""Validate the corrected-component Wronskian and resonance gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

import jensen_window_pf_newman_polymath15_critical_component_wronskian_gate as gate


EXPECTED_IDS = [
    "np15cwg_01_component_split",
    "np15cwg_02_pairwise_wronskian",
    "np15cwg_03_dirichlet_rates",
    "np15cwg_04_reference_rate",
    "np15cwg_05_crossing_system",
    "np15cwg_06_ordered_speed_guard",
    "np15cwg_07_corrected_diagnostics",
    "np15cwg_08_cancellation_guard",
    "np15cwg_09_small_ball_target",
    "np15cwg_10_proof_boundary",
]


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    expected_kind = (
        "jensen_window_pf_newman_polymath15_"
        "critical_component_wronskian_gate"
    )
    if artifact.get("kind") != expected_kind:
        issues.append("artifact kind mismatch")
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    exact = gate.build_exact()
    if artifact.get("exact") != exact:
        issues.append("exact payload drifted")
    role_counts = {
        role: sum(row.get("role") == role for row in rows)
        for role in {
            "exact_definition",
            "exact_identity",
            "exact_reduction",
            "exact_countermodel",
            "finite_diagnostics",
            "nonpromotion_gate",
            "open_theorem_target",
            "proof_guard",
        }
    }
    expected_counts = {
        "exact_definition": 1,
        "exact_identity": 3,
        "exact_reduction": 1,
        "exact_countermodel": 1,
        "finite_diagnostics": 1,
        "nonpromotion_gate": 1,
        "open_theorem_target": 1,
        "proof_guard": 1,
    }
    if role_counts != expected_counts:
        issues.append(f"role counts drifted: {role_counts}")
    open_row = next(
        (row for row in rows if row.get("role") == "open_theorem_target"),
        {},
    )
    if open_row.get("readiness") != "not_ready_to_apply":
        issues.append("small-ball target was promoted")
    guard = next(
        (row for row in rows if row.get("role") == "exact_countermodel"),
        {},
    )
    if guard.get("readiness") != "guard_validated":
        issues.append("ordered-speed guard is not active")

    fresh_coarse = gate.diagnostics(gate.COARSE_DPS)
    fresh_fine = gate.diagnostics(gate.FINE_DPS)
    fresh_convergence = gate.compare_diagnostics(fresh_coarse, fresh_fine)
    if artifact.get("diagnostics") != fresh_fine:
        issues.append("stored corrected diagnostics drifted")
    if artifact.get("convergence") != fresh_convergence:
        issues.append("stored convergence audit drifted")
    diagnostic_rows = fresh_fine.get("rows", [])
    if len(diagnostic_rows) != 4:
        issues.append("expected four corrected crossing diagnostics")
    if not all(
        row.get("all_component_phase_speeds_negative")
        and row.get("dirichlet_phase_speeds_strictly_increasing")
        and row.get("endpoint_is_maximum_phase_speed")
        for row in diagnostic_rows
    ):
        issues.append("ordered negative component-speed diagnostic failed")
    wronskians = [
        mp.mpf(row["direct_wronskian"]) for row in diagnostic_rows
    ]
    if wronskians and not (min(wronskians) < 0 < max(wronskians)):
        issues.append("aggregate Wronskian signs do not straddle zero")
    if any(
        mp.mpf(row["pairwise_identity_relative_error"]) >= mp.mpf("1e-40")
        for row in diagnostic_rows
    ):
        issues.append("pairwise numerical identity lost precision")
    if mp.mpf(
        fresh_fine.get("maximum_absolute_term_cancellation_factor", "0")
    ) <= 500:
        issues.append("Lehmer cancellation stress missing")

    status = artifact.get("status", "")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "small-ball theorem remains open",
        "countermodel",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "does not prove",
        "positive-time simplicity",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = gate.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pairwise Identity",
        "rank(K)<=2",
        "inertia (1,1,m-2)",
        "Ordered-Speed Countermodel",
        "Corrected Diagnostics",
        "Live Target",
        "not a proof of `Lambda <= 0` or RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=gate.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman Polymath-15 critical component Wronskian gate: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        "5 exact identities/reductions, 1 exact ordered-speed countermodel, "
        "4 corrected diagnostics, 1 open arithmetic small-ball target"
    )


if __name__ == "__main__":
    main()
