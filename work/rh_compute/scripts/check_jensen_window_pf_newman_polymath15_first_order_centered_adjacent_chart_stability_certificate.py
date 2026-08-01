#!/usr/bin/env python3
"""Validate the centered adjacent-chart stability certificate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_adjacent_chart_stability_certificate as builder


EXPECTED_IDS = [
    "nfocacs_01_domain",
    "nfocacs_02_ratio",
    "nfocacs_03_scaled_blocks",
    "nfocacs_04_value_budget",
    "nfocacs_05_w_derivative",
    "nfocacs_06_epsilon_derivative",
    "nfocacs_07_ratio_bounds",
    "nfocacs_08_endpoint_rate",
    "nfocacs_09_adjacent_real_jet",
    "nfocacs_10_nuisances",
    "nfocacs_11_scalar_stability",
    "nfocacs_12_absorption",
    "nfocacs_13_bulk_handoff",
    "nfocacs_14_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    t, delta, r1 = sp.symbols("t delta r1")
    h = sp.Rational(1, 2) - sp.I * sp.pi * t / 8
    if sp.simplify(
        -h * delta
        + sp.I * sp.pi * t * (r1 - delta) / 8
        + delta / 2
        - sp.I * sp.pi * t * r1 / 8
    ) != 0:
        issues.append("independent heat/log cancellation failed")

    epsilon, w = sp.symbols(
        "epsilon w", positive=True, real=True
    )
    delta_series = sp.series(
        sp.log(1 + epsilon * w), epsilon, 0, 3
    ).removeO()
    psi_series = sp.series(
        (
            2 * sp.log(1 + epsilon * w)
            - 2 * epsilon * w
            + epsilon**2 * w**2
        )
        / epsilon**2,
        epsilon,
        0,
        3,
    ).removeO()
    b1 = -w / 2 - 2 * sp.I * sp.pi * w**3 / 3
    geometric = (
        -delta_series / 2
        - sp.I * sp.pi * psi_series
        - sp.series(
            sp.log(1 + epsilon * b1), epsilon, 0, 3
        ).removeO()
    )
    if sp.simplify(sp.expand(geometric).coeff(epsilon, 1)) != 0:
        issues.append("independent geometric first coefficient failed")

    a = sp.symbols("a", positive=True, real=True)
    a_x = 1 / (8 * sp.pi * a)
    eps_expr = 1 / a
    if sp.simplify(
        sp.diff(eps_expr, a) * a_x
        + eps_expr**3 / (8 * sp.pi)
    ) != 0:
        issues.append("independent epsilon derivative failed")

    alpha_value, alpha_one, log_a, k_rate = sp.symbols(
        "alpha alpha_one log_a K_rate"
    )
    s_x = -sp.I / 2 - sp.I * t * alpha_one / 4
    m_x = -sp.I * (
        alpha_value + t * alpha_value * alpha_one / 2
    ) / 2
    mu = k_rate - m_x + s_x * log_a
    chi = alpha_value - log_a
    expected = (
        k_rate
        + sp.I * chi / 2
        + sp.I * t * alpha_one * chi / 4
    )
    if sp.simplify(mu - expected) != 0:
        issues.append("independent mu cancellation failed")


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    if artifact.get("source_sha256") != builder.source_hashes():
        issues.append("source hash chain drifted")
    if artifact.get("source_audit") != builder.source_audit():
        issues.append("source audit drifted")
    if artifact.get("exact") != builder.symbolic_audit():
        issues.append("exact payload drifted")
    if artifact.get("interval") != builder.interval_budget():
        issues.append("interval payload drifted")
    expected_constants = {
        "log_ratio": 8,
        "log_ratio_w": 30,
        "log_ratio_epsilon": 100,
        "log_ratio_x": 2,
        "ratio_value": 11,
        "ratio_x": 3,
        "endpoint": 100,
        "endpoint_rate": 2,
        "adjacent_jet": 2500,
        "explicit_d": 1500,
        "frame": 5,
        "scalar_chart": 5000,
    }
    if artifact.get("constants") != expected_constants:
        issues.append("constant payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row
        for row in rows
        if row.get("readiness") == "not_ready_to_apply"
    ]
    if len(open_rows) != 1:
        issues.append(f"open obligation count drifted: {len(open_rows)}")

    independent_symbolic_audit(issues)
    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "epsilon<=exp(-25)",
        "rho=exp(Omega)*(1+d_+)/B",
        "Omega=h*r0+R_M-delta/2",
        "r0=epsilon^2/(4*pi*i)",
        "w_x=-epsilon/(8*pi)",
        "mu_a=K_x/K+i*(alpha-log(a))/2",
        "|Z|<8*epsilon^2",
        "|rho-1|<11*epsilon^2",
        "|Delta U|<2500*exp(-7L/4)",
        "|A_(a,N+1)-A_(a,N)|<5000*exp(-7L/4)",
        "<10^-7*exp(-5L/4)",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove",
        "q>=1",
        "q<1",
        "finite phase cells",
        "one-sided successor winding",
        "contact exclusion",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Scaled Domain",
        "Exact Corrected Ratio",
        "Value Bound",
        "Derivative Bound",
        "Centered Endpoint Rate",
        "Chart-Stability Theorem",
        "bounding an apparent",
        "`1/epsilon` derivative",
        "not the bulk scalar lower bound",
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
        "validated Newman first-order centered adjacent-chart stability: "
        "14 rows, |rho-1|<11/a^2, |rho_x|<3/a^3, "
        "|Delta A|<5000*e^-7L/4, 1 open bulk Xi obligation"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
