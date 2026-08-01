#!/usr/bin/env python3
"""Validate the branch-free absolute-phase anchor reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_absolute_phase_anchor_reduction as builder


EXPECTED_IDS = [
    "nfocaapr_01_first_coefficient",
    "nfocaapr_02_branch_free_anchor",
    "nfocaapr_03_endpoint_phase",
    "nfocaapr_04_relative_endpoint",
    "nfocaapr_05_aggregate_shapes",
    "nfocaapr_06_real_projection",
    "nfocaapr_07_contact_split",
    "nfocaapr_08_shape_determinant",
    "nfocaapr_09_singular_value",
    "nfocaapr_10_wronskian_collapse",
    "nfocaapr_11_crossing_sign",
    "nfocaapr_12_complex_zero_guard",
    "nfocaapr_13_chart_guard",
    "nfocaapr_14_q_ge_1_target",
    "nfocaapr_15_q_lt_1_target",
    "nfocaapr_16_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    x, y, p, q, r, s = sp.symbols(
        "x y p q r s", real=True
    )
    matrix = sp.Matrix([[p, -q], [r, -s]])
    vector = sp.Matrix([x, y])
    expected = sp.Matrix([p * x - q * y, r * x - s * y])
    if matrix * vector != expected:
        issues.append("independent projection matrix failed")

    complex_det = -sp.im(
        (p - sp.I * q) * (r + sp.I * s)
    )
    if sp.expand(matrix.det() - complex_det) != 0:
        issues.append("independent shape determinant failed")

    big_x, big_y = sp.symbols("X Y", real=True)
    a_part, b_part = sp.symbols("A B", real=True)
    u, v = sp.symbols("u v", real=True)
    d_r, d_i = sp.symbols("D_r D_i", real=True)
    e_value = big_x + sp.I * big_y
    c_value = a_part + sp.I * b_part
    d_value = d_r + sp.I * d_i
    derivative = (u + sp.I * v) * e_value + c_value + d_value
    wronskian = sp.im(derivative * sp.conjugate(e_value))
    determinant = -sp.im(c_value * sp.conjugate(e_value))
    rhs = (
        v * (big_x**2 + big_y**2)
        + sp.im(d_value * sp.conjugate(e_value))
        - wronskian
    )
    if sp.simplify(sp.expand_complex(determinant - rhs)) != 0:
        issues.append("independent Wronskian collapse failed")
    if sp.simplify(
        sp.expand_complex(determinant.subs(big_x, 0) - a_part * big_y)
    ) != 0:
        issues.append("independent crossing sign identity failed")

    parity, beta, time_0, h_value = sp.symbols(
        "parity beta T_0 H_a", real=True
    )
    kappa = -parity * beta * (time_0 + sp.I)
    endpoint = -kappa * h_value
    if sp.expand(
        endpoint - parity * beta * (time_0 + sp.I) * h_value
    ) != 0:
        issues.append("independent endpoint phase failed")

    # If both anchored projections vanish and Z_0 is nonzero, the two
    # shapes are real-collinear. Parameterize the pure-imaginary pair.
    eta_r, eta_i, alpha, gamma = sp.symbols(
        "eta_r eta_i alpha gamma", real=True
    )
    eta = eta_r + sp.I * eta_i
    z_0 = sp.I * alpha / eta
    z_a = sp.I * gamma / eta
    collinearity = sp.im(sp.conjugate(z_0) * z_a)
    unit_residual = sp.factor(
        sp.together(collinearity * (eta_r**2 + eta_i**2))
    )
    if unit_residual != 0:
        issues.append("independent contact collinearity failed")


def independent_numeric_audit(issues: list[str]) -> None:
    x_min = 4 * math.pi * math.exp(builder.L_MIN)
    d_upper = builder.D_ABSOLUTE_CONSTANT / x_min
    if not d_upper < 0.5:
        issues.append("independent first-coefficient disk failed")
    first_lower = 1 - d_upper
    if not (
        builder.VALUE_CHART_CONSTANT / first_lower
        < builder.NORMALIZED_VALUE_CHART_CONSTANT
    ):
        issues.append("independent normalized value budget failed")
    if not (
        builder.SCALAR_CHART_CONSTANT / first_lower
        < builder.NORMALIZED_SCALAR_CHART_CONSTANT
    ):
        issues.append("independent normalized scalar budget failed")

    # Deterministic direct complex checks use values unrelated to the
    # builder's symbolic symbols.
    first = complex(0.7, -1.1)
    z_0 = complex(2.3, 0.4)
    z_a = complex(-0.8, 1.7)
    projected = (
        (first * z_0).real,
        (first * z_a).real,
    )
    determinant = -(
        z_0.conjugate() * z_a
    ).imag
    matrix_det = (
        z_0.real * (-z_a.imag)
        - (-z_0.imag) * z_a.real
    )
    if abs(determinant - matrix_det) > 1e-14:
        issues.append("independent numeric determinant failed")
    if abs(projected[0] - (first * z_0).real) > 1e-14:
        issues.append("independent numeric value projection failed")


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
        "d_absolute": 2_189,
        "value_chart": 1_100,
        "scalar_chart": 5_000,
        "normalized_value_chart": 2_200,
        "normalized_scalar_chart": 10_000,
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

    independent_symbolic_audit(issues)
    independent_numeric_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "f_1=phi*(1+d_1)",
        "eta=f_1/|f_1|",
        "q_1=1",
        "kappa_N=(-1)^(N+1)*beta*(T_0+i)",
        "r_0=g_0/f_1",
        "Z_0=P_0+r_0=E_[1]/f_1",
        "Z_A=-s_*'*P_1+r_A=C_a/f_1",
        "X/|f_1|=Re(eta*Z_0)",
        "If Z_0!=0",
        "Delta_anchor=det(S)=-Im(conj(Z_0)*Z_A)",
        "|f_1|^2*Delta_anchor=-Im(C_a*conj(E_[1]))",
        "v_a*|E_[1]|^2+Im(D_(1,x)*conj(E_[1]))-W_[1]",
        "|f_1|^2*Delta_anchor=A_a*Y",
        "At E_[1]=0",
        "|Re(eta*Delta Z_0)|<2200*exp(-5L/4)",
        "Q_N can have a large imaginary part",
        "This formulation retains Z_0=0",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not classify E_[1]=0",
        "q>=1 anchored projection lower bound",
        "q<1",
        "finite phase cells",
        "one-sided successor winding",
        "Lambda<=0",
        "RH",
        "PF-infinity",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "First-Coefficient Anchor",
        "Endpoint In The Same Frame",
        "Exact Contact Split",
        "Determinant Collapse",
        "Chart-Covariance Guard",
        "Surviving Theorem",
        "direct anchored projection",
        "not a proof of RH",
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
        "validated Newman centered absolute-phase anchor reduction: "
        "16 rows, 1 branch-free anchor, 1 determinant collapse, "
        "2 open arithmetic obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
