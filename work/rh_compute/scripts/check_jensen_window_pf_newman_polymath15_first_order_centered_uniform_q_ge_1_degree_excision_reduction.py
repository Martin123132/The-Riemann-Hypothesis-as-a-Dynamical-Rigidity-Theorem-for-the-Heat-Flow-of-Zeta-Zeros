#!/usr/bin/env python3
"""Validate the uniform-q>=1 outer-collar degree-excision reduction."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_uniform_q_ge_1_degree_excision_reduction as builder


EXPECTED_IDS = [
    "uqde_00_pi_provenance",
    "uqde_01_successor_geometry",
    "uqde_02_outer_interface",
    "uqde_03_endpoint_thresholds",
    "uqde_04_interface_regularity",
    "uqde_05_two_regime_cover",
    "uqde_06_interface_abel_side",
    "uqde_07_conditional_noncontact",
    "uqde_08_chart_guard",
    "uqde_09_oriented_chain",
    "uqde_10_degree_excision",
    "uqde_11_many_turn_guard",
    "uqde_12_boundary_only_guard",
    "uqde_13_route_fork",
    "uqde_14_outer_gap_target",
    "uqde_15_inner_degree_target",
]


def independent_geometry_audit(issues: list[str]) -> None:
    ell = sp.symbols("L", real=True, positive=True)
    t_bottom = sp.Rational(25) / ell
    t_top = sp.Rational(25) / (ell - 1)
    if sp.simplify(
        1 / sp.sqrt(2 * t_bottom) - sp.sqrt(ell / 50)
    ) != 0:
        issues.append("bottom q=1 threshold identity failed")
    if sp.simplify(
        1 / (2 * t_top) - (ell - 1) / 50
    ) != 0:
        issues.append("top q=1 threshold identity failed")
    if sp.Rational(1, 2 * 50**2) != sp.Rational(1, 5000):
        issues.append("interface kink time failed")
    if not 101**2 > sp.Rational(101, 50):
        issues.append("bottom threshold/right-edge inequality failed")
    if not 101**2 > sp.Rational(100, 50):
        issues.append("top threshold/right-edge inequality failed")


def independent_chain_audit(issues: list[str]) -> None:
    outer = Counter(
        {
            "bottom_outer": 1,
            "right": 1,
            "top_outer": 1,
            "interface_outer": 1,
        }
    )
    inner = Counter(
        {
            "bottom_inner": 1,
            "interface_outer": -1,
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
        issues.append("independent interface chain cancellation failed")


def independent_heat_guards(issues: list[str]) -> None:
    x, t, m = sp.symbols("x t m", real=True, positive=True)
    field = sp.exp(m**2 * t) * sp.sin(m * x + sp.pi / 4)
    slope = sp.diff(field, x)
    if sp.simplify(sp.diff(field, t) + sp.diff(field, x, 2)) != 0:
        issues.append("many-turn backward-heat identity failed")
    jet_norm = sp.trigsimp(field**2 + slope**2)
    expected_norm = sp.exp(2 * m**2 * t) * (
        sp.sin(m * x + sp.pi / 4) ** 2
        + m**2 * sp.cos(m * x + sp.pi / 4) ** 2
    )
    if sp.simplify(jet_norm - expected_norm) != 0:
        issues.append("many-turn jet norm identity failed")
    current = sp.trigsimp(
        field * sp.diff(slope, x) - slope * sp.diff(field, x)
    )
    if sp.simplify(current + m**2 * sp.exp(2 * m**2 * t)) != 0:
        issues.append("many-turn argument current failed")

    f = x**2 - 2 * t
    if sp.simplify(sp.diff(f, t) + sp.diff(f, x, 2)) != 0:
        issues.append("boundary-only heat identity failed")
    jacobian = sp.Matrix(
        [
            [sp.diff(f, x), sp.diff(f, t)],
            [sp.diff(sp.diff(f, x), x), sp.diff(sp.diff(f, x), t)],
        ]
    )
    if jacobian.det().subs({x: 0, t: 0}) != 4:
        issues.append("boundary-only local degree determinant failed")


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
    if len(rows) != 16:
        issues.append("row count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 2:
        issues.append("open-obligation count drifted")
    if len([row for row in rows if row.get("role") == "countermodel"]) != 2:
        issues.append("route-guard count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "conditional_topological_reduction"
        ]
    ) != 1:
        issues.append("conditional-excision count drifted")

    independent_geometry_audit(issues)
    independent_chain_audit(issues)
    independent_heat_guards(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "L_*(t)=max(50,(2t)^(-1/2))",
        "q=2tL^2>=1",
        "t=1/5000",
        "A_j={tL<=25}",
        "G_j={tL>=25}",
        "meet at tL=25",
        "J_H=(H_t,partial_x H_t)",
        "partial D_j=partial Omega_j^in+partial Omega_j^out",
        "deg(J_H,Omega_j^out,0)",
        "wind(J_H(partial Omega_j^in),0)",
        "forward horizontal path winds -m",
        "Jacobian determinant in (x,t) coordinates is +4",
        "H_j<3*pi/2",
        "whole closed collar",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the uniform Xi Abel-scalar gap",
        "exact-H noncontact consequence",
        "inner q<1/bounded-L successor theorem",
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
        "Successor Collar",
        "Two-Regime Cover",
        "Conditional Noncontact",
        "Oriented Excision",
        "Route Guards",
        "Route Fork",
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
        "validated Newman uniform-q>=1 degree excision: "
        "16 rows, 1 exact interface chain, "
        "1 conditional outer-degree excision, "
        "1 paired many-turn cancellation guard, "
        "1 boundary-only contact guard, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
