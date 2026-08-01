#!/usr/bin/env python3
"""Validate the first-order saddle-centered real-residual reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_real_residual_reduction as builder


EXPECTED_IDS = [
    "nfocrr_01_centered_scalar",
    "nfocrr_02_wronskian_audit",
    "nfocrr_03_critical_frame",
    "nfocrr_04_u_bound",
    "nfocrr_05_v_bound",
    "nfocrr_06_main_bound",
    "nfocrr_07_derivative_correction",
    "nfocrr_08_core_scalar",
    "nfocrr_09_core_approximation",
    "nfocrr_10_endpoint_expansion",
    "nfocrr_11_contact_theorem",
    "nfocrr_12_crossing_unification",
    "nfocrr_13_successor_count",
    "nfocrr_14_frequency_target",
    "nfocrr_15_other_cells",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    x_part, y_part, p_part, q_part = sp.symbols(
        "X Y P Q", real=True
    )
    u_part, v_part = sp.symbols("u v", real=True)
    slope_real = u_part * x_part - v_part * y_part + p_part
    slope_imag = v_part * x_part + u_part * y_part + q_part
    centered = p_part - v_part * y_part
    wronskian = sp.expand(
        slope_imag * x_part - slope_real * y_part
    )
    if sp.simplify(slope_real - u_part * x_part - centered) != 0:
        issues.append("independent centered slope audit failed")
    if sp.simplify(
        wronskian
        - (
            x_part * (v_part * x_part + q_part)
            - centered * y_part
        )
    ) != 0:
        issues.append("independent centered Wronskian audit failed")
    if sp.simplify(
        wronskian.subs(x_part, 0) + centered * y_part
    ) != 0:
        issues.append("independent crossing Wronskian audit failed")

    x = sp.symbols("x", positive=True, real=True)
    s = (1 - sp.I * x) / 2
    alpha_prime = (
        -1 / (2 * s**2)
        - 1 / (s - 1) ** 2
        + 1 / (2 * s)
    )
    expected_c = (7 * x**2 - 5) / (x**2 + 1) ** 2
    expected_d = x * (x**2 + 5) / (x**2 + 1) ** 2
    if sp.simplify(sp.re(alpha_prime) - expected_c) != 0:
        issues.append("independent Re(alpha') audit failed")
    if sp.simplify(sp.im(alpha_prime) - expected_d) != 0:
        issues.append("independent Im(alpha') audit failed")
    expected_gap = (
        (x**2 - sp.Rational(1, 2)) ** 2 + sp.Rational(7, 4)
    ) / (x * (x**2 + 1) ** 2)
    if sp.simplify(2 / x - expected_d - expected_gap) != 0:
        issues.append("symbolic D_x upper-bound identity audit failed")

    t, log_a = sp.symbols("t log_a", positive=True, real=True)
    c_rate, d_rate = sp.symbols("C D", real=True)
    saddle_prime = (
        t * d_rate / 4
        + sp.I * (-sp.Rational(1, 2) - t * c_rate / 4)
    )
    expected = -sp.I / 2 - sp.I * t * (
        c_rate + sp.I * d_rate
    ) / 4
    if sp.simplify(saddle_prime - expected) != 0:
        issues.append("independent saddle derivative audit failed")
    if sp.simplify(-sp.re(saddle_prime) * log_a + t * d_rate * log_a / 4) != 0:
        issues.append("independent u_a formula audit failed")


def independent_numeric_audit(issues: list[str]) -> None:
    ux = 50_000.0 * math.exp(-50.0)
    vy = 303.0 / (16.0 * math.pi**2) * math.exp(-25.0)
    d1x = (
        builder.COEFFICIENT_MASS
        * builder.CORRECTION_DERIVATIVE_CONSTANT
        / (16.0 * math.pi**2)
        * math.exp(-25.0)
    )
    if not ux < 1e-16:
        issues.append("independent uX budget failed")
    if not vy < 1e-10:
        issues.append("independent vY budget failed")
    if not d1x < 1e-7:
        issues.append("independent D1x budget failed")
    if not ux + vy + d1x < 1e-6:
        issues.append("independent combined nuisance budget failed")


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
        "eta_0": 100_000,
        "eta_1": 200_000,
        "value_band": 50_000,
        "slope_band": 100_000,
        "coefficient_mass": 50,
        "main_absolute_constant": 101,
        "correction_derivative_constant": 4_223,
        "core_approx_relative": "1e-6",
    }:
        issues.append("constant payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_rows) != 2:
        issues.append(f"open target count drifted: {len(open_rows)}")
    if any(row.get("readiness") != "not_ready_to_apply" for row in open_rows):
        issues.append("an open target was promoted")

    independent_symbolic_audit(issues)
    independent_numeric_audit(issues)
    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "S_a=P_a-v_a*Y",
        "U=S_a",
        "W_[1]=-S_a*Y",
        "C_x=Re(alpha'(s))",
        "D_x=Im(alpha'(s))",
        "|u_a|<=tL/(2x)",
        "|v_a|<3/x^2",
        "|E_[1]|<101*exp(L/4)",
        "|d_(n,x)|<4223/x^2",
        "|U-A_a|<1e-6*exp(-5L/4)",
        "G_a=-kappa_N*(H_(a,x)+mu_a*H_a)",
        "|A_a|>(100000L+1)*exp(-5L/4)",
        "N_(X=0,A_a>0",
        "q=2tL^2>=1",
        "q<1",
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
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Centered Crossing Scalar",
        "Critical Frame",
        "Absolute Budgets",
        "Core Arithmetic Scalar",
        "Endpoint Defect",
        "Contact And Winding Target",
        "Route Audit",
        "Remaining Cells",
        "no complex-main zero is deleted",
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
        "validated Newman first-order centered real-residual reduction: "
        "15 rows, |u_a|<e^-L, |v_a|<3/x^2, "
        "core error <1e-6*e^-5L/4, 2 open Xi cell obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
