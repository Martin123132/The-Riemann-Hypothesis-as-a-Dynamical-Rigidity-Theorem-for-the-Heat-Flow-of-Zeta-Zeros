#!/usr/bin/env python3
"""Validate the reciprocal normalizer-phase reinforcement guard."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_reciprocal_normalizer_phase_reinforcement_guard as builder


EXPECTED_IDS = [
    "rnpr_00_source_coordinate",
    "rnpr_01_pi_provenance",
    "rnpr_02_zero_time_ratio",
    "rnpr_03_zero_time_bound",
    "rnpr_04_heat_phase_transport",
    "rnpr_05_heat_phase_bound",
    "rnpr_06_complex_lock",
    "rnpr_07_reinforcement_guard",
    "rnpr_08_projection_guard",
    "rnpr_09_exponent_guard",
    "rnpr_10_endpoint_guard",
    "rnpr_11_route_rejection",
    "rnpr_12_replacement_target",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    t, alpha_r, alpha_i = sp.symbols("t alpha_R alpha_I", real=True)
    h = -t * alpha_i / 2
    if sp.simplify(t * alpha_r * alpha_i + 2 * h * alpha_r) != 0:
        issues.append("independent heat phase cancellation failed")

    pi_upper = Fraction(22, 7)
    coefficient = pi_upper**2 / 64 + pi_upper / 16
    if coefficient != Fraction(275, 784):
        issues.append("independent rational pi coefficient failed")
    if not coefficient < Fraction(3, 8):
        issues.append("independent heat phase bound failed")

    if not Fraction(15, 16) < 1:
        issues.append("final phase coefficient failed")
    if not builder.C_TWO_DEFICIT == Fraction(
        3_133_668_399, 48_144_906_818
    ):
        issues.append("c=2 exponent deficit drifted")


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
    if artifact.get("source_urls", {}).get("polymath15_primary") != builder.POLYMATH_SOURCE_URL:
        issues.append("primary source URL drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if len(rows) != 13:
        issues.append("row count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 1:
        issues.append("open-target count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") in {"countermodel", "nonpromotion_guard"}
        ]
    ) != 3:
        issues.append("nonpromotion-guard count drifted")

    independent_symbolic_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "gamma_0/Q(T)=exp(i*delta_0(T))",
        "0<delta_0(T)<3/(8T)=3/(4x)",
        "psi_t=-delta_0(T)+integral_T^omega log(v/T_0)dv",
        "|psi_t|<15/(16x)<1/x",
        "|R-1|<8/x",
        "modulus greater than 2-8/x",
        "z+R*conj(z)=2Re(z)+O(|z|/x)",
        "N^(-d_2-eta)",
        "raw primal-plus-reciprocal destructive interference",
        "weighted-minus-unweighted tail",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "rejects only a power gain from raw reciprocal pairing",
        "does not evaluate the signed weighted-minus-unweighted tail kernel",
        "prove cancellation",
        "improve c_*",
        "Abel-scalar gap",
        "contact",
        "Lambda<=0",
        "PF-infinity",
        "RH",
        "Clay-prize conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Exact Normalizer",
        "Pi Provenance",
        "Zero-Time Phase",
        "Heat Transport",
        "Conjugate Lock",
        "Reinforcement Guard",
        "Endpoint Scope",
        "Route Decision",
        "Replacement Target",
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
        "validated Newman reciprocal normalizer-phase reinforcement guard: "
        "13 rows, 2 exact phase identities, 1 uniform conjugate lock, "
        "3 nonpromotion guards, 1 rejected raw-pair route, "
        "1 open signed-kernel target"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
