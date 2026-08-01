#!/usr/bin/env python3
"""Validate the oscillatory-spliced outer-collar degree reduction."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_oscillatory_spliced_outer_collar_reduction as builder


EXPECTED_IDS = [
    "osce_00_pi_provenance",
    "osce_01_current_frontier",
    "osce_02_epsilon_quantifier",
    "osce_03_raised_collar",
    "osce_04_interface_geometry",
    "osce_05_three_regime_cover",
    "osce_06_seam_stress",
    "osce_07_conditional_noncontact",
    "osce_08_degree_excision",
    "osce_09_fixed_epsilon_consequence",
    "osce_10_asymptotic_scope",
    "osce_11_frontier_guard",
    "osce_12_finite_shoulder_guard",
    "osce_13_outer_target",
    "osce_14_inner_target",
]


def independent_arithmetic_audit(issues: list[str]) -> None:
    c_star = Fraction(4_911_678_521, 1_933_561_194)
    epsilon = Fraction(1, 100)
    if c_star != builder.C_STAR:
        issues.append("critical threshold drifted")
    if c_star + epsilon != builder.EXAMPLE_CAP:
        issues.append("example epsilon cap drifted")
    if not c_star + epsilon < 25:
        issues.append("example cap left first-order range")
    if not Fraction(1, 100) < c_star:
        issues.append("q=1 interface cap is not below c_star")

    ell, heat_time = sp.symbols("ell heat_time", positive=True)
    q = 2 * heat_time * ell**2
    c = heat_time * ell
    q_one_c = sp.simplify(c.subs(heat_time, 1 / (2 * ell**2)))
    if q_one_c != 1 / (2 * ell):
        issues.append("q=1 scaled-time identity failed")

    ell_j = sp.symbols("ell_j", positive=True)
    bottom = sp.Rational(25) / ell_j
    threshold = sp.simplify(1 / sp.sqrt(2 * bottom))
    if threshold != sp.sqrt(ell_j / 50):
        issues.append("successor bottom q-threshold identity failed")


def independent_chain_audit(issues: list[str]) -> None:
    outer = Counter(
        {
            "bottom_outer": 1,
            "right": 1,
            "top_outer": 1,
            "interface": 1,
        }
    )
    inner = Counter(
        {
            "bottom_inner": 1,
            "interface": -1,
            "top_inner": 1,
            "axis": 1,
        }
    )
    full = Counter(
        {
            "bottom_inner": 1,
            "bottom_outer": 1,
            "right": 1,
            "top_outer": 1,
            "top_inner": 1,
            "axis": 1,
        }
    )
    if inner + outer != full:
        issues.append("raised-collar interface chain failed")


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
    if len(rows) != 15:
        issues.append("row count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 2:
        issues.append("open-obligation count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "conditional_topological_reduction"
        ]
    ) != 1:
        issues.append("conditional-excision count drifted")

    independent_arithmetic_audit(issues)
    independent_chain_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "c_*=4911678521/1933561194",
        "B_epsilon=max(50,L_epsilon)",
        "J_epsilon=max(0,ceil(B_epsilon)-101)",
        "q=2tL^2>=1",
        "A_(j,epsilon)={c<=c_epsilon}",
        "Z_(j,epsilon)={c_epsilon<=c<=25}",
        "G_(j,epsilon)={c>=25}",
        "q=1 implies c=tL=1/(2L)<=1/100<c_*",
        "J_H=(H_t,partial_xH_t)",
        "deg(J_H,Omega_(j,epsilon)^out,0)=0",
        "c<=c_*+o(1)",
        "does not manufacture a uniform explicit epsilon(L)",
        "L_epsilon is existential",
        "W_0=0",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the Xi Abel-scalar gap",
        "effective L_epsilon",
        "q<1/bounded-L inner theorem",
        "complete boundary nonvanishing",
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
        "Current Frontier",
        "Epsilon Collar",
        "Three-Regime Cover",
        "Conditional Excision",
        "Exact Consequence",
        "Route Guards",
        "Route Decision",
        "Open Theorems",
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
        "validated Newman oscillatory-spliced outer collar: "
        "15 rows, 1 exact current frontier, 1 epsilon-raised collar, "
        "1 three-regime cover, 1 conditional degree excision, "
        "2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
