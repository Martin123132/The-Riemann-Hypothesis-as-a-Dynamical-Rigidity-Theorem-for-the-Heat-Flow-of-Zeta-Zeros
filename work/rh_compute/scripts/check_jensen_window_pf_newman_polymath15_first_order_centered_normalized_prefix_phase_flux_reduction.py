#!/usr/bin/env python3
"""Validate the normalized-prefix and first-jet phase-flux reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_normalized_prefix_phase_flux_reduction as builder


EXPECTED_IDS = [
    "nfnpfr_00_pi_provenance",
    "nfnpfr_01_normalized_prefix",
    "nfnpfr_02_normalized_slope",
    "nfnpfr_03_coefficient_current",
    "nfnpfr_04_inward_current",
    "nfnpfr_05_ordered_rotation",
    "nfnpfr_06_prefix_currents",
    "nfnpfr_07_recrossing_countermodel",
    "nfnpfr_07a_dyadic_completion",
    "nfnpfr_08_signed_flux_guard",
    "nfnpfr_09_first_jet_flux",
    "nfnpfr_10_flux_integrand",
    "nfnpfr_11_q_ge_1_target",
    "nfnpfr_12_q_lt_1_target",
    "nfnpfr_13_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    theta = sp.symbols("theta", real=True)
    core = sp.cos(theta) + sp.cos(3 * theta) / 2
    expected = sp.cos(theta) * (
        2 * sp.cos(theta) ** 2 - sp.Rational(1, 2)
    )
    if sp.trigsimp(sp.expand_trig(core) - expected) != 0:
        issues.append("independent recrossing factorization failed")
    roots = [
        sp.pi / 3,
        sp.pi / 2,
        2 * sp.pi / 3,
        4 * sp.pi / 3,
        3 * sp.pi / 2,
        5 * sp.pi / 3,
    ]
    derivative = sp.diff(core, theta)
    values = [sp.simplify(derivative.subs(theta, root)) for root in roots]
    if any(sp.simplify(core.subs(theta, root)) != 0 for root in roots):
        issues.append("independent recrossing root failed")
    if len([value for value in values if value.is_positive]) != 3:
        issues.append("independent upward count failed")
    if len([value for value in values if value.is_negative]) != 3:
        issues.append("independent downward count failed")

    c = sp.symbols("c", real=True)
    r = sp.sqrt(2) / 2
    completed = (
        sp.chebyshevt(1, c)
        + r * sp.chebyshevt(2, c)
        + r**2 * sp.chebyshevt(3, c)
    )
    polynomial = sp.expand(2 * sp.sqrt(2) * completed)
    if sp.discriminant(polynomial, c) != -1696:
        issues.append("independent dyadic discriminant failed")
    if not sp.simplify(polynomial.subs(c, -1)) < 0:
        issues.append("independent dyadic left sign failed")
    if not sp.simplify(polynomial.subs(c, 1)) > 0:
        issues.append("independent dyadic right sign failed")

    rho1, rho2, nu1, nu2, delta = sp.symbols(
        "rho1 rho2 nu1 nu2 delta", real=True
    )
    r1, r2 = sp.symbols("r1 r2", positive=True, real=True)
    q1 = r1 * (sp.cos(delta) + sp.I * sp.sin(delta))
    q2 = r2
    gamma1 = rho1 + sp.I * nu1
    gamma2 = rho2 + sp.I * nu2
    product = sp.expand_complex(
        (gamma1 * q1 + gamma2 * q2) * sp.conjugate(q1 + q2)
    )
    radial = (
        rho1 * r1**2
        + rho2 * r2**2
        + r1
        * r2
        * (
            (rho1 + rho2) * sp.cos(delta)
            - (nu1 - nu2) * sp.sin(delta)
        )
    )
    angular = (
        nu1 * r1**2
        + nu2 * r2**2
        + r1
        * r2
        * (
            (rho1 - rho2) * sp.sin(delta)
            + (nu1 + nu2) * sp.cos(delta)
        )
    )
    if sp.trigsimp(sp.re(product) - radial) != 0:
        issues.append("independent prefix radial kernel failed")
    if sp.trigsimp(sp.im(product) - angular) != 0:
        issues.append("independent prefix angular kernel failed")


def independent_numeric_audit(issues: list[str]) -> None:
    l_value = float(builder.L_MIN)
    x_min = 4 * math.pi * math.exp(l_value)
    radial_floor = math.log(2) / (16 * l_value**2 * x_min)
    error = 4 * builder.D_X_CONSTANT / x_min**2
    if not error < radial_floor:
        issues.append("independent inward-current margin failed")

    t_max = 25 / l_value
    a_max = math.sqrt(math.exp(l_value) + t_max / 16)
    n_value = math.floor(a_max)
    width = 4 * math.pi * (2 * n_value + 1)
    lower = width * (1 / (2 * n_value) - error)
    if not lower > 4 * math.pi:
        issues.append("independent adjacent rotation margin failed")


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
    if artifact.get("constants") != {
        "L_min": 50,
        "d_x": 4_223,
        "d": 2_189,
    }:
        issues.append("constant payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row
        for row in rows
        if row.get("readiness") == "not_ready_to_apply"
    ]
    if len(open_rows) != 2:
        issues.append(f"open target count drifted: {len(open_rows)}")
    countermodels = [
        row for row in rows if row.get("role") == "countermodel"
    ]
    if len(countermodels) != 1:
        issues.append(f"countermodel count drifted: {len(countermodels)}")

    independent_symbolic_audit(issues)
    independent_numeric_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2",
        "No circle or prefix polygon is selected",
        "G_k=conj(eta)*F_k=sum_(n=1)^k q_n",
        "G_1=1",
        "Z_A=s_*'*u_N*Z_0+B_N",
        "mathcal_C_N=mathsf_A",
        "gamma_n=q_(n,x)/q_n",
        "rho_n<-log(2)/(16L^2*x)<0",
        "every adjacent relative phase advances by more than 4*pi",
        "Im(G_(k,x)*conj(G_k))",
        "zero_count",
        "upward_count",
        "first_jet_winding",
        "dyadic_completion_scout",
        "4sqrt(2)c^3+4c^2-sqrt(2)c-2",
        "\"discriminant\": \"-1696\"",
        "N_up-N_down",
        "partial_x arg(Gamma_ell)",
        "u_(N,x)=1/(8*pi*a^2)=1/(4*T_0)",
        "0<=kappa_j=(2*pi)^(-1)*Delta_arg(Gamma_j)<1",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the Xi prefix lower bound",
        "strict composed successor flux bound",
        "q<1",
        "finite shoulder",
        "Lambda<=0",
        "PF-infinity",
        "RH",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pi Provenance",
        "Normalized Polygon",
        "Carrier Currents",
        "Recrossing Guard",
        "multiplicative chain repairs",
        "First-Jet Flux",
        "Live Theorem",
        "not an Xi",
        "Clay-prize conclusion",
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
        "validated Newman normalized-prefix phase-flux reduction: "
        "15 rows, 2 certified spiral bounds, 1 exact recrossing guard, "
        "1 dyadic completion diagnostic, "
        "1 O(N) first-jet flux, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
