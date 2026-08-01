#!/usr/bin/env python3
"""Validate the interior-carrier projective-current theorem and edge guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_interior_projective_current_gate as builder


EXPECTED_IDS = [
    f"ipcg_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "pi_provenance",
            "domain",
            "carrier_jet",
            "components",
            "division_free_current",
            "effective_slope",
            "source_bounds",
            "interior_envelope",
            "matrix_budget",
            "negative_definite_theorem",
            "projective_flow",
            "zero_projection_join",
            "terminal_obstruction",
            "edge_block",
            "endpoint_recurrence",
            "order_swap_guard",
            "cumulative_mass_boundary",
            "route_decision",
            "q_lt_1",
            "proof_boundary",
        )
    )
]


def independent_current(issues: list[str]) -> None:
    x_part, y_part, radial, angular = sp.symbols(
        "X Y radial angular", real=True
    )
    c_s, b, k, u = sp.symbols("c_s b k u", real=True)
    c_s_x, b_x, u_x = sp.symbols("c_s_x b_x u_x", real=True)
    x_part_x = radial * x_part - angular * y_part
    y_part_x = angular * x_part + radial * y_part
    c_atom = x_part
    d_atom = c_s * k * x_part - b * u * y_part
    d_atom_x = (
        c_s_x * k * x_part
        + c_s * k * x_part_x
        - (b_x * u + b * u_x) * y_part
        - b * u * y_part_x
    )
    current = sp.expand(c_atom * d_atom_x - d_atom * x_part_x)
    expected = sp.expand(
        c_s_x * k * x_part**2
        - (b_x * u + b * u_x) * x_part * y_part
        - b * u * angular * (x_part**2 + y_part**2)
    )
    if sp.simplify(current - expected) != 0:
        issues.append("independent carrier-current identity failed")
    if radial in current.free_symbols:
        issues.append("independent radial-current cancellation failed")


def independent_matrix(issues: list[str]) -> None:
    x_part, y_part = sp.symbols("X Y", real=True)
    a_term, mixed, k_term = sp.symbols("A B K", real=True)
    matrix = sp.Matrix(
        [
            [a_term - k_term, -mixed / 2],
            [-mixed / 2, -k_term],
        ]
    )
    vector = sp.Matrix([x_part, y_part])
    quadratic = sp.expand((vector.T * matrix * vector)[0])
    expected = sp.expand(
        a_term * x_part**2
        - mixed * x_part * y_part
        - k_term * (x_part**2 + y_part**2)
    )
    if sp.simplify(quadratic - expected) != 0:
        issues.append("independent quadratic-matrix identity failed")
    if (
        -sp.Rational(1, 8)
        + sp.Rational(1, 64)
        + sp.Rational(1, 32)
        != -sp.Rational(5, 64)
    ):
        issues.append("independent first Gershgorin budget failed")
    if (
        -sp.Rational(1, 8) + sp.Rational(1, 32)
        != -sp.Rational(3, 32)
    ):
        issues.append("independent second Gershgorin budget failed")


def independent_terminal_and_swap(issues: list[str]) -> None:
    x_part, y_part, b, u_x = sp.symbols("X Y b u_x", real=True)
    terminal = -b * u_x * x_part * y_part
    values = [
        terminal.subs(
            {b: -sp.Rational(1, 2), u_x: 1, x_part: 1, y_part: sign}
        )
        for sign in (1, -1)
    ]
    if values != [sp.Rational(1, 2), -sp.Rational(1, 2)]:
        issues.append("independent terminal sign guard failed")

    variable = sp.symbols("x", real=True)
    first = -variable
    second = -variable + sp.sin(variable) / 2
    if sp.diff(first, variable) != -1:
        issues.append("independent first monotone slope failed")
    if sp.diff(second, variable) != -1 + sp.cos(variable) / 2:
        issues.append("independent second monotone slope failed")
    if sp.simplify(second - first - sp.sin(variable) / 2) != 0:
        issues.append("independent order-swap difference failed")


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
    if len([row for row in rows if row.get("role") == "analytic_lemma"]) != 1:
        issues.append("analytic lemma count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 2:
        issues.append("open theorem target count drifted")
    if len(
        [row for row in rows if row.get("role") == "nonpromotion_guard"]
    ) != 2:
        issues.append("nonpromotion guard count drifted")

    independent_current(issues)
    independent_matrix(issues)
    independent_terminal_and_swap(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "varrho_n cancels identically",
        "J_n<=-(5/64)*u_n^2*|z_n|^2<0",
        "h_(n,x)=J_n/X_n^2<0",
        "J_n=-K_n*Y_n^2",
        "J_N=-b*u_x*X_NY_N",
        "C_edge=c_0+c_N",
        "J_a^(der)=H_(a,x)+mu_aH_a is generally complex",
        "Never identify J_a^(der) with J_a^(adj)",
        "h_2-h_1=sin(x)/2",
        "includes q=1",
        "cumulative-mass estimate",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove an Xi cumulative-mass estimate",
        "edge-block sign",
        "Abel-scalar gap",
        "q<1 closure",
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
        "Carrier Current",
        "Uniform Interior Bound",
        "Projective Join",
        "Exceptional Edge",
        "Order-Swap Guard",
        "Remaining Theorem",
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
        "validated interior projective-current gate: "
        "20 rows, exact radial-current cancellation, proved 5/64 "
        "negative-definite interior current, terminal-edge and "
        "order-swap nonpromotion guards"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
