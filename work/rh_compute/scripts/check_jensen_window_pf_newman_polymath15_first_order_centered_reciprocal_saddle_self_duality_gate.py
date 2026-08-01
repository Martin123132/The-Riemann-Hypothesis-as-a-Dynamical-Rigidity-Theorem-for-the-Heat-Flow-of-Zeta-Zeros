#!/usr/bin/env python3
"""Validate the reciprocal-saddle heat-amplitude self-duality gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_saddle_self_duality_gate as builder


EXPECTED_IDS = [
    "rssd_00_pi_provenance",
    "rssd_01_model_phase",
    "rssd_02_dual_phase",
    "rssd_03_amplitude_identity",
    "rssd_04_actual_xi_defect",
    "rssd_05_phase_saddle_comparison",
    "rssd_06_uniform_defect_bound",
    "rssd_07_cutoff_guard",
    "rssd_08_critical_block_map",
    "rssd_09_tail_interpretation",
    "rssd_10_hypothetical_guard",
    "rssd_11_correction_guard",
    "rssd_12_route_decision",
    "rssd_13_open_theorem",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    t, a_log, z, sigma = sp.symbols("t A z sigma", real=True)
    y = 2 * a_log - z
    lhs = (
        t * y**2 / 4
        - sigma * y
        + y
        - a_log
        - (t * z**2 / 4 - sigma * z)
    )
    rhs = (2 * sigma - 1 - t * a_log) * (z - a_log)
    if sp.simplify(sp.expand(lhs - rhs)) != 0:
        issues.append("independent reciprocal amplitude identity failed")
    omega = sp.symbols("omega", positive=True)
    physical_phase = omega * (2 * a_log - z) - omega
    reciprocal_phase = -omega * z + omega * (2 * a_log - 1)
    if sp.simplify(physical_phase - reciprocal_phase) != 0:
        issues.append("independent physical reciprocal phase failed")

    x = sp.symbols("x", positive=True)
    s = sp.Rational(1, 2) - sp.I * x / 2
    rational = sp.simplify(
        sp.im(sp.expand_complex(1 / (2 * s) + 1 / (s - 1)))
    )
    if rational != 3 * x / (1 + x**2):
        issues.append("independent alpha imaginary rational part failed")

    r_star = Fraction(125_662, 155_153)
    if 2 - r_star != Fraction(184_644, 155_153):
        issues.append("independent critical reciprocal radius failed")
    if not 2 - r_star > 1:
        issues.append("critical reciprocal block did not leave cutoff")


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    if artifact.get("source_sha256") != builder.source_hashes():
        issues.append("source hash chain drifted")
    if artifact.get("source_audit") != builder.source_audit():
        issues.append("source audit drifted")
    if artifact.get("exact") != builder.build_exact():
        issues.append("exact payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if len(rows) != 14:
        issues.append("row count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 1:
        issues.append("open-theorem count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") in {"countermodel", "nonpromotion_guard"}
        ]
    ) != 3:
        issues.append("guard count drifted")

    independent_symbolic_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "a_omega^2=omega/(2*pi)",
        "physical phase exp[i*omega log(u)]",
        "u_nu=a_omega^2/nu",
        "Delta=2*sigma-1-t*log(a_omega)",
        "Delta_Xi=t*(Re(alpha)-log(a_omega))",
        "0<T_0-omega<7t/(4x)",
        "|Delta_Xi|<t/(2x)",
        "exp[-13/(2x)]",
        "2-r_*=184644/155153",
        "3133668399/48144906818",
        "self-duality alone proves neither cancellation",
        "last saddle and endpoint",
        "fixed-c<2 zero-free wall",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove a discrete Poisson/B-process remainder theorem",
        "reciprocal block cancellation",
        "improved c_* threshold",
        "Abel-scalar gap",
        "contact exclusion",
        "Lambda<=0",
        "PF-infinity",
        "RH",
        "Clay-prize conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pi Provenance",
        "Reciprocal Saddle",
        "Exact Heat Self-Duality",
        "Xi Defect Bound",
        "Critical Block",
        "Hypothetical Stress Test",
        "Correction Guard",
        "Route Decision",
        "Open Theorem",
        "Boundary",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=builder.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman reciprocal-saddle self-duality gate: "
        "14 rows, 1 exact stationary map, 1 exact heat-amplitude "
        "self-duality, 1 Xi defect bound, 1 critical reciprocal-tail "
        "map, 2 nonpromotion guards, 1 open theorem"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
