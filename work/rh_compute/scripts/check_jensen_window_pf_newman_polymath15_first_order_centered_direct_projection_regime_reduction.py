#!/usr/bin/env python3
"""Validate the direct anchored-projection regime reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_direct_projection_regime_reduction as builder


EXPECTED_IDS = [
    "nfocdprr_01_carriers",
    "nfocdprr_02_saddle_amplitude",
    "nfocdprr_03_endpoint_projection",
    "nfocdprr_04_direct_pair",
    "nfocdprr_05_relative_derivative",
    "nfocdprr_06_anchored_derivative",
    "nfocdprr_07_frame_bound",
    "nfocdprr_08_ordinary_identity",
    "nfocdprr_09_ordinary_condition",
    "nfocdprr_10_exceptional_condition",
    "nfocdprr_11_zero_fiber",
    "nfocdprr_12_heat_guard",
    "nfocdprr_13_chart_guard",
    "nfocdprr_14_q_ge_1_target",
    "nfocdprr_15_q_lt_1_target",
    "nfocdprr_16_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    log_a, u, delta, time = sp.symbols(
        "log_a u delta time", real=True
    )
    y = log_a - u
    alpha_r = log_a - delta
    sigma = sp.Rational(1, 2) + time * alpha_r / 2
    direct = time * y**2 / 4 - sigma * y
    edge = time * log_a**2 / 4 - sigma * log_a
    centered = (
        edge
        + (sp.Rational(1, 2) - time * delta / 2) * u
        + time * u**2 / 4
    )
    if sp.expand(direct - centered) != 0:
        issues.append("independent saddle centering failed")

    z_0, z_0x, nu, lam, d_term = sp.symbols(
        "z_0 z_0x nu lam d_term"
    )
    z_a = z_0x + (nu - lam) * z_0 - d_term
    lhs = nu * z_0 + z_0x
    rhs = lam * z_0 + z_a + d_term
    if sp.expand(lhs - rhs) != 0:
        issues.append("independent relative derivative identity failed")

    x_part, y_part, p_part, q_part = sp.symbols(
        "X Y p q", real=True
    )
    ratio = p_part + sp.I * q_part
    value = x_part + sp.I * y_part
    if sp.expand_complex(
        sp.re(ratio * value)
        - (p_part * x_part - q_part * y_part)
    ) != 0:
        issues.append("independent ratio projection failed")

    c, b, ell, cosine, sine = sp.symbols(
        "c b ell cosine sine", real=True
    )
    projected = sp.re(
        -(c + sp.I * b) * ell * (cosine + sp.I * sine)
    )
    if sp.expand_complex(
        projected - ell * (-c * cosine + b * sine)
    ) != 0:
        issues.append("independent carrier slope projection failed")

    x, t = sp.symbols("x t", real=True)
    model = t - x**2 / 2 + sp.I * x
    if sp.simplify(sp.diff(model, t) + sp.diff(model, x, 2)) != 0:
        issues.append("independent heat-model PDE failed")
    jacobian = sp.Matrix(
        [
            [sp.diff(sp.re(model), x), sp.diff(sp.re(model), t)],
            [sp.diff(sp.im(model), x), sp.diff(sp.im(model), t)],
        ]
    )
    at_origin = jacobian.subs({x: 0, t: 0})
    if at_origin.det() != -1:
        issues.append("independent heat-model Jacobian failed")
    if sp.re(sp.diff(model, x)).subs({x: 0, t: 0}) != 0:
        issues.append("independent directional-zero guard failed")


def independent_numeric_audit(issues: list[str]) -> None:
    x_min = 4 * math.pi * math.exp(builder.L_MIN)
    d_upper = builder.D_ABSOLUTE_CONSTANT / x_min
    first_lower = 1 - d_upper
    rho_upper = builder.D_X_CONSTANT / x_min**2 / first_lower
    frame = rho_upper + math.exp(-builder.L_MIN) + 3 / x_min**2
    if not frame < 2 * math.exp(-builder.L_MIN):
        issues.append("independent frame bound failed")
    if not 1e-7 / first_lower < 2e-7:
        issues.append("independent normalized D bound failed")
    split_ratio = 100_000 * math.exp(-builder.L_MIN / 4) / builder.L_MIN
    if not split_ratio < 0.01:
        issues.append("independent explicit split width failed")
    if not math.sqrt(1 - split_ratio**2) > 0.9999:
        issues.append("independent explicit split square-root failed")

    tau = 3.0
    delta = 0.4
    kappa = 0.7
    cap = 1.2
    lower = kappa * math.sqrt(tau**2 - delta**2) - cap * delta
    y_abs = math.sqrt(tau**2 - delta**2)
    worst_direct = kappa * y_abs - cap * delta
    if abs(lower - worst_direct) > 1e-14:
        issues.append("independent ordinary lower-bound arithmetic failed")

    ratio = complex(cap, kappa)
    value = complex(delta, -y_abs)
    direct_projection = (ratio * value).real
    if abs(direct_projection - (cap * delta + kappa * y_abs)) > 1e-14:
        issues.append("independent complex ratio evaluation failed")


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
        "d_x": 4_223,
        "normalized_D": 2e-7,
        "frame": 2,
        "split_ratio_denominator": 100,
        "split_sqrt": 0.9999,
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
        "eta*q_n=r_n*zeta_n",
        "u_n=log(a/n)=-ell_n",
        "B_a=1/2-t*delta_a/2>49/100",
        "Re(eta*r_0)=(-1)^N*B_0*T_0*H_a",
        "mathsf_X:=Re(eta*Z_0)",
        "mathsf_A:=Re(eta*Z_A)",
        "Z_A=Z_(0,x)+(nu_1-lambda_a)*Z_0-D_(1,x)/f_1",
        "W_A=W_(0,x)+(rho_1-lambda_a)*W_0",
        "|rho_1-lambda_a|<2*exp(-L)",
        "mathsf_A=p*mathsf_X-q*mathsf_Y",
        "q=-Delta_anchor/R^2",
        "h=W_(0,x)/W_0+rho_1-lambda_a-D_(1,x)/E_[1]",
        "theta_x=Im(W_(0,x)*conj(W_0))/|W_0|^2",
        "|D_(1,x)/E_[1]|<2e-7*exp(-5L/4)/tau",
        "kappa*sqrt(tau^2-delta^2)-K*delta",
        "Choose tau_L=L*exp(-L)",
        "delta_L/tau_L<100000*exp(-L/4)/L<1/100",
        "sqrt(tau_L^2-delta_L^2)>0.9999*tau_L",
        "|mathsf_A-partial_x mathsf_X|",
        "At Z_0=0",
        "W_*(x,t)=t-x^2/2+i*x",
        "det D_(x,t)(Re W_*,Im W_*)=-1",
        "directional real-slope bound",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the q>=1 ratio",
        "directional-transversality",
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
        "Direct Unit Carriers",
        "Derivative Identity",
        "Ordinary Regime",
        "Exceptional Regime",
        "Live Theorem",
        "directional real transversality",
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
        "validated Newman centered direct-projection regime reduction: "
        "16 rows, 2 exact projection normal forms, "
        "1 heat-direction countermodel, 2 open Xi regimes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
