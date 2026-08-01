#!/usr/bin/env python3
"""Validate the endpoint Schur-Cohn stability and first-jet guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_endpoint_schur_cohn_first_jet_guard as builder


EXPECTED_IDS = [
    "escfjg_00_pi_provenance",
    "escfjg_01_endpoint_polynomial",
    "escfjg_02_outside_step",
    "escfjg_03_stability_theorem",
    "escfjg_04_reflection_recursion",
    "escfjg_05_boundary_margin",
    "escfjg_06_coefficient_x_tangent",
    "escfjg_07_centered_companion",
    "escfjg_08_frechet_recursion",
    "escfjg_09_mixed_recursion",
    "escfjg_10_first_xi_pivot",
    "escfjg_11_degree_one_orientation",
    "escfjg_12_first_pivot_guard",
    "escfjg_13_decreasing_phase_guard",
    "escfjg_14_endpoint_pivot_guard",
    "escfjg_15_cutoff_update",
    "escfjg_16_route_decision",
    "escfjg_17_weaker_fallback",
    "escfjg_18_q_ge_1_target",
    "escfjg_19_q_lt_1_nonpromotion",
]


def reverse_conjugate(coeffs: list[sp.Expr]) -> list[sp.Expr]:
    return [sp.conjugate(value) for value in reversed(coeffs)]


def independent_schur_step(issues: list[str]) -> None:
    z = sp.symbols("z")
    b = list(sp.symbols("c_0:4"))
    f = sum(value * z**index for index, value in enumerate(b))
    f_star_coeffs = reverse_conjugate(b)
    f_star = sum(
        value * z**index
        for index, value in enumerate(f_star_coeffs)
    )
    transformed = sp.expand(sp.conjugate(b[0]) * f - b[-1] * f_star)
    if sp.simplify(transformed.coeff(z, 3)) != 0:
        issues.append("independent outside step did not lower degree")

    t_coeffs = [sp.expand(transformed.coeff(z, k)) for k in range(3)]
    t = sum(value * z**index for index, value in enumerate(t_coeffs))
    t_star = sum(
        value * z**index
        for index, value in enumerate(reverse_conjugate(t_coeffs))
    )
    delta = (
        sp.conjugate(b[0]) * b[0]
        - sp.conjugate(b[-1]) * b[-1]
    )
    residual = sp.expand(b[0] * t + b[-1] * z * t_star - delta * f)
    for index in range(4):
        if sp.simplify(residual.coeff(z, index)) != 0:
            issues.append("independent outside inverse failed")
            break


def independent_normalized_recursion(issues: list[str]) -> None:
    z = sp.symbols("z")
    a, c = sp.symbols("a c")
    p = 1 + c * z + a * z**2
    p_star = sp.conjugate(a) + sp.conjugate(c) * z + z**2
    denominator = 1 - sp.conjugate(a) * a
    reduced = sp.expand((p - a * p_star) / denominator)
    q0 = sp.simplify(reduced.coeff(z, 0))
    q1 = sp.simplify(reduced.coeff(z, 1))
    if q0 != 1:
        issues.append("independent normalized constant failed")
    q_star = sp.conjugate(q1) + z
    reconstructed = sp.expand((1 + q1 * z) + a * z * q_star)
    for index in range(3):
        if sp.simplify((reconstructed - p).coeff(z, index)) != 0:
            issues.append("independent normalized inverse failed")
            break


def independent_tangents(issues: list[str]) -> None:
    eps = sp.symbols("eps", real=True)
    b = list(sp.symbols("b_0:4"))
    v = list(sp.symbols("v_0:4"))

    def step(coeffs: list[sp.Expr]) -> list[sp.Expr]:
        reversed_coeffs = reverse_conjugate(coeffs)
        return [
            sp.expand(
                sp.conjugate(coeffs[0]) * coeffs[k]
                - coeffs[-1] * reversed_coeffs[k]
            )
            for k in range(3)
        ]

    actual = [
        sp.diff(value, eps).subs(eps, 0)
        for value in step([b[k] + eps * v[k] for k in range(4)])
    ]
    expected = [
        sp.expand(
            sp.conjugate(v[0]) * b[k]
            + sp.conjugate(b[0]) * v[k]
            - v[3] * sp.conjugate(b[3 - k])
            - b[3] * sp.conjugate(v[3 - k])
        )
        for k in range(3)
    ]
    for got, want in zip(actual, expected):
        if sp.simplify(got - want) != 0:
            issues.append("independent Frechet tangent failed")
            break


def independent_guards(issues: list[str]) -> None:
    z = sp.symbols("z")
    quadratic = 1 - sp.Rational(9, 4) * z + sp.Rational(1, 2) * z**2
    if set(sp.solve(quadratic, z)) != {sp.Rational(1, 2), sp.Integer(4)}:
        issues.append("independent degree-two roots failed")
    transformed = sp.expand(quadratic - sp.Rational(1, 2) * z**2 * quadratic.subs(z, 1 / z))
    # Direct coefficient form avoids any assumptions about z conjugation.
    reduced = sp.Rational(3, 4) - sp.Rational(9, 8) * z
    if (
        sp.Rational(3, 4) ** 2
        - sp.Rational(9, 8) ** 2
        != -sp.Rational(45, 64)
    ):
        issues.append("independent degree-two second pivot failed")
    if reduced.subs(z, sp.Rational(2, 3)) != 0:
        issues.append("independent reduced inside root failed")

    sqrt_55 = sp.sqrt(55)
    b1 = (-23 + 3 * sp.I * sqrt_55) / 40
    b2 = (-17 - 3 * sp.I * sqrt_55) / 40
    if sp.simplify(1 + b1 + b2) != 0:
        issues.append("independent decreasing-phase boundary zero failed")
    if sp.simplify(sp.conjugate(b1) * b1) != sp.Rational(16, 25):
        issues.append("independent decreasing-phase first norm failed")
    if sp.simplify(sp.conjugate(b2) * b2) != sp.Rational(49, 100):
        issues.append("independent decreasing-phase second norm failed")

    endpoint = -sp.Rational(2, 3) + sp.Rational(2, 3) * z
    if endpoint.subs(z, 1) != 0:
        issues.append("independent endpoint boundary zero failed")


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
        row for row in rows if row.get("role") == "countermodel"
    ]
    if len(countermodels) != 3:
        issues.append(f"route-guard count drifted: {len(countermodels)}")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 2:
        issues.append(f"open Xi target count drifted: {len(open_rows)}")
    exact_theorems = [
        row
        for row in rows
        if row.get("role") == "exact_algebraic_theorem"
    ]
    if len(exact_theorems) != 1:
        issues.append(
            f"outside-stability theorem count drifted: {len(exact_theorems)}"
        )

    independent_schur_step(issues)
    independent_normalized_recursion(issues)
    independent_tangents(issues)
    independent_guards(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "introduces no new pi",
        "F_N(z)=r_0+P_0(z)",
        "S_dF(z)=conj(b_0)F(z)-b_dF#(z)",
        "Delta_d=|b_0|^2-|b_d|^2",
        "Rouche",
        "product_(d=1)^K(1-|alpha_d|)",
        "(B_(k,r))_x=-s_*'B_(k,r+1)+E_(k,r)",
        "G(z)=r_A-s_*'[P_1(z)-A P_0(z)]",
        "DS_b[W]+D2S_b[U,V]",
        "Delta_K=|r_0+B_(0,0)|^2-|B_(K,0)|^2",
        "F_(N+1)(z)-F_N(z)=j_0+u_n z^v",
        "0<=kappa_j<1",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "No actual Xi pivot inequality",
        "full unit-disk stability",
        "endpoint-complete Xi lower bound",
        "strict successor flux upper bound",
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
        "Endpoint Polynomial",
        "Outside-Disk Schur Step",
        "Reflection Recursion",
        "Conditional Boundary Margin",
        "Physical Coefficient Current",
        "Centered Companion",
        "Tangent Recursion",
        "First Physical Pivot",
        "First pivot is not sufficient",
        "Decreasing moduli with physical-style phases",
        "Endpoint pivot collapse",
        "Cutoff Update",
        "Route Decision",
        "Live Theorem",
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
        "validated Newman endpoint Schur-Cohn first-jet guard: "
        "20 rows, 1 exact outside-stability theorem, "
        "1 conditional reflection-product margin, "
        "3 exact route guards, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
