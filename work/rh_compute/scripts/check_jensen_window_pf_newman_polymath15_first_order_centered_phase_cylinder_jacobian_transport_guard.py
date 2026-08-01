#!/usr/bin/env python3
"""Validate the phase-cylinder Jacobian and transport guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_phase_cylinder_jacobian_transport_guard as builder


EXPECTED_IDS = [
    "pcjtg_00_pi_provenance",
    "pcjtg_01_phase_lift",
    "pcjtg_02_phase_partials",
    "pcjtg_03_pair_jacobian",
    "pcjtg_04_linked_current",
    "pcjtg_05_winding_transport",
    "pcjtg_06_xi_pair_sign_guard",
    "pcjtg_07_linked_monotonicity_guard",
    "pcjtg_08_correction_cylinder",
    "pcjtg_09_endpoint_cylinder",
    "pcjtg_10_endpoint_guard",
    "pcjtg_11_three_cylinders",
    "pcjtg_12_endpoint_first_jet",
    "pcjtg_13_route_decision",
    "pcjtg_14_q_ge_1_target",
    "pcjtg_15_q_lt_1_nonpromotion",
]


def independent_pair_kernel(issues: list[str]) -> None:
    k_1, k_2, l_1, l_2, a_1, a_2 = sp.symbols(
        "k_1 k_2 l_1 l_2 a_1 a_2", real=True
    )
    phi_1, phi_2 = sp.symbols("phi_1 phi_2", real=True)
    u_1 = a_1 * sp.exp(sp.I * phi_1)
    u_2 = a_2 * sp.exp(sp.I * phi_2)
    theta_vector = sp.I * (k_1 * u_1 + k_2 * u_2)
    xi_vector = sp.I * (l_1 * u_1 + l_2 * u_2)
    jacobian = sp.expand_complex(
        sp.conjugate(theta_vector) * xi_vector
    ).as_real_imag()[1]
    expected = (
        a_1
        * a_2
        * (k_1 * l_2 - k_2 * l_1)
        * sp.sin(phi_2 - phi_1)
    )
    if sp.trigsimp(jacobian - expected) != 0:
        issues.append("independent phase pair kernel failed")


def independent_linked_current(issues: list[str]) -> None:
    z, xi, h, log_m = sp.symbols("z xi h log_m")
    k = sp.symbols("k", integer=True, nonnegative=True)
    monomial = z**k * sp.exp(sp.I * xi * log_m)
    p_theta = sp.I * z * sp.diff(monomial, z)
    p_xi = sp.diff(monomial, xi)
    expected = sp.I * (h * k + log_m) * monomial
    if sp.simplify(h * p_theta + p_xi - expected) != 0:
        issues.append("independent linked derivative failed")

    p_r, p_i, h_r, h_i = sp.symbols(
        "p_r p_i h_r h_i", real=True
    )
    p_value = p_r + sp.I * p_i
    h_value = h_r + sp.I * h_i
    lhs = sp.expand_complex(
        sp.conjugate(p_value) * sp.I * h_value
    ).as_real_imag()[1]
    rhs = sp.expand_complex(
        sp.conjugate(p_value) * h_value
    ).as_real_imag()[0]
    if sp.expand(lhs - rhs) != 0:
        issues.append("independent linked argument current failed")


def independent_guards(issues: list[str]) -> None:
    h = sp.log(2)
    linked_frequency = sp.simplify(sp.log(3) - h)
    if linked_frequency != sp.log(sp.Rational(3, 2)):
        issues.append("independent n=2,3 linked frequency failed")

    r = sp.sqrt(2) / 2
    if not sp.simplify(r**2 - r) < 0:
        issues.append("independent linked monotonicity guard failed")
    if not 1 - r > 0:
        issues.append("independent linked nonvanishing guard failed")

    z = sp.symbols("z")
    prefix = 1 + sp.Rational(2, 3) * z
    if sp.solve(prefix, z)[0] != -sp.Rational(3, 2):
        issues.append("independent endpoint prefix root failed")
    if sp.simplify(prefix.subs(z, 1) - sp.Rational(5, 3)) != 0:
        issues.append("independent endpoint augmentation failed")


def independent_cylinders(issues: list[str]) -> None:
    mu, d_n, d_1 = sp.symbols("mu d_n d_1")
    ratio = (1 + mu * d_n) / (1 + mu * d_1)
    expected = (d_n - d_1) / (1 + mu * d_1) ** 2
    if sp.simplify(sp.diff(ratio, mu) - expected) != 0:
        issues.append("independent correction-cylinder derivative failed")

    theta = sp.symbols("theta", real=True)
    r_prime = sp.symbols("r_prime", positive=True, real=True)
    p_theta = sp.I * sp.exp(sp.I * theta)
    p_xi = -r_prime
    jacobian = sp.expand_complex(
        sp.conjugate(p_theta.subs(theta, 0)) * p_xi
    ).as_real_imag()[1]
    if jacobian != r_prime:
        issues.append("independent winding orientation failed")


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
    if artifact.get("constants") != {"L_min": 50}:
        issues.append("constant payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 2:
        issues.append(f"open Xi target count drifted: {len(open_rows)}")
    countermodels = [
        row for row in rows if row.get("role") == "countermodel"
    ]
    if len(countermodels) != 3:
        issues.append(f"route-guard count drifted: {len(countermodels)}")
    topology = [
        row
        for row in rows
        if row.get("role") == "exact_topological_theorem"
    ]
    if len(topology) != 1:
        issues.append(
            f"winding-transport theorem count drifted: {len(topology)}"
        )

    independent_pair_kernel(issues)
    independent_linked_current(issues)
    independent_guards(issues)
    independent_cylinders(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "do not define pi",
        "P_theta=i*sum_alpha",
        "J_(theta,xi)=Im(conj(P_theta)*P_xi)",
        "k_alpha*l_beta-k_beta*l_alpha",
        "dP/dxi=h*P_theta+P_xi=i*mathcal D_hP=i*H_1",
        "wind_theta P(.,xi_1)-wind_theta P(.,xi_0)=",
        "J_23=w_2*w_3*log(3)*sin(xi*log(3/2))",
        "1/2-1/sqrt(2)",
        "J_(theta,mu)=Im(conj(P_theta)P_mu)",
        "J_(theta,tau)=Im(conj(P_theta)r_0)",
        '"augmented_zero": "P(1)+e=0"',
        "0<=kappa_j<1",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "No sign or nonzero lower bound for the full Xi Jacobian",
        "no three-cylinder crossing budget",
        "no joined Xi lower bound",
        "no strict successor flux upper bound",
        "no q<1 closure",
        "no contact exclusion",
        "no Lambda<=0",
        "no RH proof",
        "no Clay-prize conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pi Provenance",
        "Phase Cylinder",
        "Pair-Kernel Jacobian",
        "Physical Linked Direction",
        "Winding Transport",
        "Genuine dyadic/odd pair",
        "Zero-free but nonmonotone linked block",
        "Correction Cylinder",
        "Endpoint Cylinder",
        "Three-Cylinder Programme",
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
        "validated Newman phase-cylinder Jacobian transport guard: "
        "16 rows, 1 exact pair-kernel Jacobian, "
        "1 cylinder winding-transport theorem, "
        "3 exact route guards, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
