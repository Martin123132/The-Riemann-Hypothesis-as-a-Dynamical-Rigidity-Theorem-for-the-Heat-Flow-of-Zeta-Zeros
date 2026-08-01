#!/usr/bin/env python3
"""Validate the endpoint odd-fibre correlation feasibility guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_endpoint_odd_fibre_correlation_feasibility_guard as builder


EXPECTED_IDS = [
    "eofcfg_00_pi_provenance",
    "eofcfg_01_physical_coefficients",
    "eofcfg_02_odd_fibre",
    "eofcfg_03_endpoint",
    "eofcfg_04_terminal",
    "eofcfg_05_physical_pivot",
    "eofcfg_06_correlation",
    "eofcfg_07_cartesian_pivot",
    "eofcfg_08_endpoint_rotation",
    "eofcfg_09_total_decomposition",
    "eofcfg_10_moment_null_family",
    "eofcfg_11_five_current_guard",
    "eofcfg_12_linked_two_jet_guard",
    "eofcfg_13_nonpromotion",
    "eofcfg_14_feasibility_verdict",
    "eofcfg_15_retained_schur_role",
    "eofcfg_16_signed_primary_target",
    "eofcfg_17_q_lt_1_target",
]


def independent_correlation(issues: list[str]) -> None:
    x_value, y_value, gamma, time_0 = sp.symbols(
        "x y gamma T", real=True
    )
    subtotal = x_value + sp.I * y_value
    endpoint = gamma * (time_0 + sp.I)
    norm = sp.expand_complex(
        (subtotal + endpoint) * sp.conjugate(subtotal + endpoint)
    )
    expected = (
        (x_value + gamma * time_0) ** 2
        + (y_value + gamma) ** 2
    )
    if sp.simplify(norm - expected) != 0:
        issues.append("independent Cartesian pivot failed")
    correlation = sp.expand_complex(
        sp.re(subtotal * sp.conjugate(endpoint))
    )
    if sp.simplify(
        correlation - gamma * (time_0 * x_value + y_value)
    ) != 0:
        issues.append("independent endpoint correlation failed")


def independent_moment_null(issues: list[str]) -> None:
    ell, h = sp.symbols("ell h", real=True)
    weights = [1, -3, 3, -1]
    nodes = [ell + index * h for index in range(4)]
    moments = [
        sp.expand(
            sum(
                weight * node**order
                for weight, node in zip(weights, nodes)
            )
        )
        for order in range(4)
    ]
    if moments[:3] != [0, 0, 0]:
        issues.append("independent moment-null family failed")
    if sp.simplify(moments[3] + 6 * h**3) != 0:
        issues.append("independent third-moment audit failed")
    support = [3, 6, 12, 24]
    if [value for value in support if value % 2] != [3]:
        issues.append("independent odd-support audit failed")
    if 32 in support:
        issues.append("independent terminal-support audit failed")


def independent_linked_jet(issues: list[str]) -> None:
    w = sp.symbols("w")
    for lam, expected_pivot in (
        (sp.Rational(1, 4), sp.Rational(1, 2)),
        (sp.Integer(1), sp.Integer(-1)),
    ):
        polynomial = sp.expand(1 + lam * (w - 1) ** 3)
        if polynomial.subs(w, 1) != 1:
            issues.append("independent linked value failed")
        if sp.diff(polynomial, w).subs(w, 1) != 0:
            issues.append("independent linked first jet failed")
        if sp.diff(polynomial, w, 2).subs(w, 1) != 0:
            issues.append("independent linked second jet failed")
        constant = polynomial.coeff(w, 0)
        leading = polynomial.coeff(w, 3)
        if constant**2 - leading**2 != expected_pivot:
            issues.append("independent linked pivot failed")


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
        if row.get("role") == "countermodel"
    ]
    if len(countermodels) != 2:
        issues.append(f"countermodel count drifted: {len(countermodels)}")
    countermodel_consequences = [
        row for row in rows if row.get("role") == "countermodel_consequence"
    ]
    if len(countermodel_consequences) != 1:
        issues.append(
            "countermodel-consequence count drifted: "
            f"{len(countermodel_consequences)}"
        )
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 2:
        issues.append(f"open-route count drifted: {len(open_rows)}")

    independent_correlation(issues)
    independent_moment_null(issues)
    independent_linked_jet(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "S_odd=phi*O_N",
        "g_0=gamma_N*(T_0+i)",
        "|f_(p_K)|=|C_K|=rho_K",
        "physical odd-index subtotal plus endpoint",
        "gamma_N*(T_0*X_odd+Y_odd)",
        "E_[1]-S_even",
        "q_3,q_6,q_12,q_24",
        "H_0,H_1,H_2,D_0,D_1",
        "F_lambda(w)=1+lambda*(w-1)^3",
        "falsification coordinate",
        "|mathcal_C_N|>A_L+epsilon_term",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "no-go for the current Schur inputs",
        "not a proof that no future Xi-specific",
        "No uniform actual Xi pivot",
        "physical pivot failure",
        "signed prefix lower bound",
        "crossing count",
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
        "Physical Rephasing",
        "First Pivot",
        "Total Versus Odd Fibre",
        "Five-Current Null Family",
        "Linked Two-Jet Guard",
        "Feasibility Verdict",
        "Retained Schur Role",
        "Primary Target",
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
        "validated Newman endpoint odd-fibre correlation feasibility guard: "
        "18 rows, 1 exact physical odd-fibre pivot, "
        "1 endpoint-projection normal form, "
        "1 five-current null-family guard, 1 linked-two-jet guard, "
        "1 Schur route downgrade, 2 open routes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
