#!/usr/bin/env python3
"""Validate the endpoint first-pivot odd-small-ball reduction and guard."""

from __future__ import annotations

import argparse
import cmath
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_endpoint_first_pivot_odd_small_ball_guard as builder


EXPECTED_IDS = [
    "efposbg_00_pi_provenance",
    "efposbg_01_terminal_singleton",
    "efposbg_02_odd_numerator",
    "efposbg_03_endpoint_numerator",
    "efposbg_04_denominator_cancellation",
    "efposbg_05_small_ball",
    "efposbg_06_correlation",
    "efposbg_07_endpoint_phase",
    "efposbg_08_sufficient_routes",
    "efposbg_09_linked_model",
    "efposbg_10_density",
    "efposbg_11_triangle",
    "efposbg_12_linked_pivot_failure",
    "efposbg_13_endpoint_neighborhood",
    "efposbg_14_scope_guard",
    "efposbg_15_adjacent_boundary",
    "efposbg_16_route_decision",
    "efposbg_17_open_targets",
]


def norm_square(value: sp.Expr) -> sp.Expr:
    return sp.expand_complex(value * sp.conjugate(value))


def independent_denominator(issues: list[str]) -> None:
    x_0, x_1, y_0, y_1 = sp.symbols(
        "x_0 x_1 y_0 y_1", real=True
    )
    c_0, c_1, d_0, d_1 = sp.symbols(
        "c_0 c_1 d_0 d_1", real=True
    )
    value = x_0 + sp.I * x_1
    endpoint = y_0 + sp.I * y_1
    terminal = c_0 + sp.I * c_1
    denominator = d_0 + sp.I * d_1
    lhs = (
        norm_square((value + endpoint) / denominator)
        - norm_square(terminal / denominator)
    ) * norm_square(denominator)
    rhs = norm_square(value + endpoint) - norm_square(terminal)
    if sp.simplify(lhs - rhs) != 0:
        issues.append("independent common-denominator cancellation failed")

    expanded = (
        norm_square(value + endpoint)
        - norm_square(value)
        - norm_square(endpoint)
        - 2 * sp.re(value * sp.conjugate(endpoint))
    )
    if sp.simplify(expanded) != 0:
        issues.append("independent correlation expansion failed")


def independent_terminal_layer(issues: list[str]) -> None:
    for n_value in range(1, 8193):
        degree = n_value.bit_length() - 1
        if n_value // (1 << degree) != 1:
            issues.append(
                f"terminal singleton failed at N={n_value}"
            )
            break


def independent_triangle(issues: list[str]) -> None:
    if not 30**2 > (7**2) * 15:
        issues.append("independent a+b>1 proof failed")
    if not sp.Rational(1, 3) > sp.Rational(1, 5):
        issues.append("independent a>b proof failed")
    delta = 1 / math.sqrt(3) + 1 / math.sqrt(5) - 1
    if not 0.024 < delta < 0.025:
        issues.append("independent endpoint radius failed")

    omega = 62643 / 100
    value = (
        1
        + cmath.exp(1j * omega * math.log(3)) / math.sqrt(3)
        + cmath.exp(1j * omega * math.log(5)) / math.sqrt(5)
    )
    if abs(value) >= 0.001:
        issues.append("independent linked-phase diagnostic failed")


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
    countermodels = [
        row
        for row in rows
        if str(row.get("role", "")).startswith("countermodel")
    ]
    if len(countermodels) != 2:
        issues.append(f"countermodel count drifted: {len(countermodels)}")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 1:
        issues.append(f"open-target row count drifted: {len(open_rows)}")

    independent_denominator(issues)
    independent_terminal_layer(issues)
    independent_triangle(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "M_K=floor(N/2^K)=1",
        "B_(0,0)=O_N/(1+d_1)",
        "r_0=R_N/(1+d_1)",
        "normalization denominator affects neither",
        "d_1 remains inside O_N",
        "Dbar(-R_N,rho_K)",
        "2Re(O_N*conj(R_N))",
        "3^u*5^v=1",
        "inf_omega|O_5(omega)|=0",
        "|e|<delta_0",
        "not a counterexample",
        "actual linked-point nonvanishing",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "No uniform actual Xi small-ball exclusion",
        "physical pivot failure",
        "later Schur pivot",
        "full endpoint disk stability",
        "linked-point lower bound",
        "signed crossing theorem",
        "q<1 closure",
        "contact exclusion",
        "Lambda<=0",
        "PF-infinity",
        "RH proof",
        "Clay-prize conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pi Provenance",
        "Terminal Singleton",
        "Common Numerators",
        "Exact First Pivot",
        "Endpoint Phase",
        "Linked-Phase Countermodel",
        "Small-Endpoint Extension",
        "Adjacent-Chart Boundary",
        "Route Decision",
        "Live Theorems",
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
        "validated Newman endpoint first-pivot odd small-ball guard: "
        "18 rows, 1 exact denominator-cancelled pivot, "
        "1 linked-logarithmic-phase countermodel, "
        "1 small-endpoint countermodel family, 2 open Xi routes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
