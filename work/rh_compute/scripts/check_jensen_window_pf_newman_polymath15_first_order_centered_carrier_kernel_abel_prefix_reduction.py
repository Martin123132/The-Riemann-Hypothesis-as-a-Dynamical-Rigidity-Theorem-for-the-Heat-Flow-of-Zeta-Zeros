#!/usr/bin/env python3
"""Validate the centered carrier-kernel and Abel-prefix reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_carrier_kernel_abel_prefix_reduction as builder


EXPECTED_IDS = [
    "nfockapr_00_pi_provenance",
    "nfockapr_01_coordinates",
    "nfockapr_02_bulk_pairwise",
    "nfockapr_03_endpoint_pairwise",
    "nfockapr_04_abel_prefix",
    "nfockapr_05_endpoint_complete_prefix",
    "nfockapr_06_prefix_currents",
    "nfockapr_07_contact_scalar",
    "nfockapr_08_terminal_band",
    "nfockapr_09_all_fiber_sufficient",
    "nfockapr_10_ordinary_connection",
    "nfockapr_11_exceptional_derivative",
    "nfockapr_12_ordinary_countermodel",
    "nfockapr_13_zero_fiber_countermodel",
    "nfockapr_14_prefix_sector_guard",
    "nfockapr_15_relative_carrier_rotation",
    "nfockapr_16_q_ge_1_target",
    "nfockapr_17_q_lt_1_target",
    "nfockapr_18_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    x_height, heat_time, cutoff_index = sp.symbols(
        "x heat_time cutoff_index", real=True
    )
    shifted_height = x_height / 2 + sp.pi * heat_time / 8
    saddle_squared = shifted_height / (2 * sp.pi)
    if sp.simplify(
        saddle_squared - (x_height / (4 * sp.pi) + heat_time / 16)
    ) != 0:
        issues.append("independent saddle pi provenance failed")
    cell_left = 4 * sp.pi * (cutoff_index**2 - heat_time / 16)
    cell_right = 4 * sp.pi * (
        (cutoff_index + 1) ** 2 - heat_time / 16
    )
    if sp.simplify(
        cell_right - cell_left
        - 4 * sp.pi * (2 * cutoff_index + 1)
    ) != 0:
        issues.append("independent cutoff-cell pi provenance failed")

    z1, z2, z3, u3, h1, h2 = sp.symbols(
        "z1 z2 z3 u3 h1 h2"
    )
    u1 = u3 + h1 + h2
    u2 = u3 + h2
    total = z1 + z2 + z3
    moment = u1 * z1 + u2 * z2 + u3 * z3
    prefix = u3 * total + h1 * z1 + h2 * (z1 + z2)
    if sp.expand(moment - prefix) != 0:
        issues.append("independent Abel-prefix identity failed")

    c, b, u, x, y, vr, vi, fr, fi, h = sp.symbols(
        "c b u X Y V_R V_I F_R F_I h", real=True
    )
    s_prime = c + sp.I * b
    w = x + sp.I * y
    v = vr + sp.I * vi
    f = fr + sp.I * fi
    wa = s_prime * u * w + v + s_prime * h * f
    scalar = vr - b * u * y + h * (c * fr - b * fi)
    if sp.simplify(
        sp.expand_complex(sp.re(wa) - scalar - c * u * x)
    ) != 0:
        issues.append("independent contact-scalar identity failed")

    product = sp.expand_complex(wa * sp.conjugate(w))
    fw = sp.expand_complex(f * sp.conjugate(w))
    vw = sp.expand_complex(v * sp.conjugate(w))
    radius = x**2 + y**2
    p_expected = (
        c * u * radius
        + sp.re(vw)
        + h * (c * sp.re(fw) - b * sp.im(fw))
    )
    q_expected = (
        b * u * radius
        + sp.im(vw)
        + h * (b * sp.re(fw) + c * sp.im(fw))
    )
    if sp.simplify(sp.re(product) - p_expected) != 0:
        issues.append("independent prefix P identity failed")
    if sp.simplify(sp.im(product) - q_expected) != 0:
        issues.append("independent prefix Q identity failed")

    rho_minus_u, v_a, d_real = sp.symbols(
        "rho_minus_u v_a d_real", real=True
    )
    derivative = (
        sp.re(wa)
        - rho_minus_u * x
        - v_a * y
        + d_real
    )
    expected = (
        scalar
        + (c * u - rho_minus_u) * x
        - v_a * y
        + d_real
    )
    if sp.simplify(sp.expand_complex(derivative - expected)) != 0:
        issues.append("independent exceptional derivative failed")

    cosine = sp.Rational(-35, 37)
    sine = sp.Rational(-12, 37)
    unit = cosine + sp.I * sine
    carriers = (
        sp.Integer(1),
        sp.Rational(4, 5) * unit,
        sp.Rational(7, 10) * unit,
    )
    distances = (sp.Integer(3), sp.Integer(2), sp.Integer(1))
    s_value = sum(carriers)
    u_value = sum(
        distance * carrier
        for distance, carrier in zip(distances, carriers)
    )
    ordinary_s_prime = sp.Rational(13, 1056) - sp.I / 2
    q_value = sp.simplify(
        sp.im(
            ordinary_s_prime
            * u_value
            * sp.conjugate(s_value)
        )
    )
    r_squared = sp.simplify(s_value * sp.conjugate(s_value))
    if q_value != 0 or r_squared != sp.Rational(61, 148):
        issues.append("independent ordinary countermodel failed")

    s_prime_value = sp.Rational(1, 16) - sp.I / 2
    omega = (-8 + sp.I) / sp.sqrt(65)
    zero_carriers = (
        omega,
        -sp.Rational(3, 5) * omega,
        -sp.Rational(2, 5) * omega,
    )
    zero_sum = sp.simplify(sum(zero_carriers))
    zero_moment = sp.simplify(
        sum(
            distance * carrier
            for distance, carrier in zip(distances, zero_carriers)
        )
    )
    zero_slope = sp.simplify(s_prime_value * zero_moment)
    if zero_sum != 0 or sp.re(zero_slope) != 0 or zero_slope == 0:
        issues.append("independent zero-fiber countermodel failed")


def independent_numeric_audit(issues: list[str]) -> None:
    x_min = 4 * math.pi * math.exp(builder.L_MIN)
    rotation_lower = (
        0.5 * math.log(2)
        - 4 * builder.D_X_CONSTANT / x_min**2
    )
    if not rotation_lower > 1 / 3:
        issues.append("independent carrier rotation bound failed")
    terminal_ratio = math.exp(-1.5 * builder.L_MIN) / (
        4 * builder.L_MIN
    )
    if not terminal_ratio < builder.TERMINAL_RELATIVE_AT_L_MIN:
        issues.append("independent terminal ratio failed")
    if not terminal_ratio < 2e-35:
        issues.append("independent L=50 terminal margin failed")

    ordinary = builder.ordinary_countermodel()
    if ordinary.get("aligned_Q") != "-53/8":
        issues.append("ordinary aligned sign drifted")
    if ordinary.get("opposed_Q") != "7/40":
        issues.append("ordinary opposed sign drifted")
    if ordinary.get("kernel_zero") != "0":
        issues.append("ordinary exact zero drifted")
    if ordinary.get("radius_squared_at_zero") != "61/148":
        issues.append("ordinary nonzero radius drifted")

    zero_fiber = builder.zero_fiber_countermodel()
    if zero_fiber.get("W_0") != "0":
        issues.append("zero-fiber value drifted")
    if zero_fiber.get("Re_W_A") != "0":
        issues.append("zero-fiber real slope drifted")
    if zero_fiber.get("abs_W_A_squared") != "637/1280":
        issues.append("zero-fiber nonzero slope drifted")


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
        "terminal_band_constant": 25_000,
        "terminal_relative_at_L_min": 2e-35,
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
    kernel_rows = [
        row for row in rows if row.get("role") == "exact_kernel"
    ]
    if len(kernel_rows) != 3:
        issues.append(f"exact kernel count drifted: {len(kernel_rows)}")
    countermodels = [
        row for row in rows if row.get("role") == "countermodel"
    ]
    if len(countermodels) != 2:
        issues.append(f"countermodel count drifted: {len(countermodels)}")

    independent_symbolic_audit(issues)
    independent_numeric_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2",
        "a^2=T_0/(2*pi)=x/(4*pi)+t/16",
        "x_N=4*pi*(N^2-t/16)",
        "exp(i*(theta+2*pi))=exp(i*theta)",
        "prefix polygon",
        "W_0=e+S, W_A=g+s_*'*U",
        "P_B=Re(s_*'*U*conj(S))",
        "Q_B=Im(s_*'*U*conj(S))",
        "g*conj(e)=B_0^2*(T_0^2+1)*H_a*J_a",
        "U=u_N*S+sum_(k=1)^(N-1)h_k*F_k",
        "V_N=g-s_*'*u_N*e",
        "W_A=s_*'*u_N*W_0+V_N",
        "P=c*u_N*R^2",
        "Q=b*u_N*R^2",
        "mathcal_C_N=Re(V_N)-b*u_N*mathsf_Y",
        "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X",
        "Q=-mathsf_Y*mathcal_C_N",
        "25000*exp(-11L/4)/|f_1|",
        "<2e-35*A_L",
        "|mathcal_C_N|>A_L+epsilon_term",
        "partial_x mathsf_X=mathcal_C_N",
        "radius_squared_at_zero",
        "abs_W_A_squared",
        "weighted prefix-sector",
        "vartheta_(n,x)=Im(chi_(n,x)*conj(chi_n))",
        "vartheta_(n,x)>log(2)/2-16892/x^2>1/3",
        "Delta x=4*pi*(2N+1)",
        "takes every unit-circle value",
        "This single division-free prefix theorem includes W_0=0",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the q>=1 Xi prefix lower bound",
        "signed successor count",
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
        "Pi Provenance",
        "No particular circle is",
        "Pairwise Kernels",
        "Abel Prefix Form",
        "Continuous Crossing Scalar",
        "Exact Route Guards",
        "actual Xi relative carriers satisfy",
        "Live Theorem",
        "|W_0|^2=61/148",
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
        "validated Newman centered carrier-kernel Abel reduction: "
        "19 rows, 3 exact kernel forms, 1 all-fiber prefix scalar, "
        "2 exact countermodels, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
