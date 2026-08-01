#!/usr/bin/env python3
"""Validate the prime-power heat-block and composition guard."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_prime_power_heat_block_composition_guard as builder


EXPECTED_IDS = [
    "pphbcg_00_pi_provenance",
    "pphbcg_01_absolute_rate",
    "pphbcg_02_terminal_collar",
    "pphbcg_03_prime_power_amplitude",
    "pphbcg_04_corrected_contraction",
    "pphbcg_05_one_block_winding",
    "pphbcg_06_gaussian_mixture",
    "pphbcg_07_tilted_mass_guard",
    "pphbcg_08_joined_p_free",
    "pphbcg_09_two_block_countermodel",
    "pphbcg_10_composition_guard",
    "pphbcg_11_q_ge_1_target",
    "pphbcg_12_q_lt_1_target",
    "pphbcg_13_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    omega, b, log_a, log_n, correction = sp.symbols(
        "omega b log_a log_n correction", real=True
    )
    u_n = log_a - log_n
    v_a = omega - b * log_a
    if sp.expand(
        omega - b * log_n + correction
        - (v_a + b * u_n + correction)
    ) != 0:
        issues.append("independent absolute-rate identity failed")

    b_amp, time, u_value, h = sp.symbols(
        "B time u h", real=True
    )
    current = b_amp * u_value + time * u_value**2 / 4
    successor = (
        b_amp * (u_value - h)
        + time * (u_value - h) ** 2 / 4
    )
    expected = (
        -b_amp * h
        - time * u_value * h / 2
        + time * h**2 / 4
    )
    if sp.expand(successor - current - expected) != 0:
        issues.append("independent prime-power ratio failed")

    z = sp.symbols("z")
    r = sp.sqrt(2) / 2
    long_block = z * (1 + r * z + r**2 * z**2)
    short_block = -sp.Rational(4, 5) * z * (1 + r * z)
    joined_inner = sp.expand(5 * (long_block + short_block) / z)
    if sp.discriminant(joined_inner, z) != -sp.Rational(19, 2):
        issues.append("independent joined discriminant failed")
    joined_poly = sp.Poly(joined_inner, z)
    joined_product = sp.simplify(
        joined_poly.all_coeffs()[-1] / joined_poly.all_coeffs()[0]
    )
    if joined_product != sp.Rational(2, 5):
        issues.append("independent joined root modulus failed")
    if not sp.discriminant(joined_inner, z) < 0:
        issues.append("independent joined conjugate-pair audit failed")

    long_inner = sp.Poly(1 + r * z + r**2 * z**2, z)
    long_product = sp.simplify(
        long_inner.all_coeffs()[-1] / long_inner.all_coeffs()[0]
    )
    if long_product != 2:
        issues.append("independent long-block exterior roots failed")
    if sp.solve(1 + r * z, z)[0] != -sp.sqrt(2):
        issues.append("independent short-block exterior root failed")

    values: list[int] = []
    n_max = 40
    for base in range(1, n_max + 1, 2):
        value = base
        while value <= n_max:
            values.append(value)
            value *= 2
    if sorted(values) != list(range(1, n_max + 1)):
        issues.append("independent joined 2-free decomposition failed")


def independent_numeric_audit(issues: list[str]) -> None:
    x_min = 4 * math.pi * math.exp(builder.L_MIN)
    d_upper = builder.D_ABSOLUTE_CONSTANT / x_min
    log_derivative = (
        builder.D_X_CONSTANT / x_min**2 / (1 - d_upper)
    )
    rate_error_scaled = 3 + log_derivative * x_min**2
    if not rate_error_scaled < builder.TERMINAL_RATE_CONSTANT / 2:
        issues.append("independent terminal-rate threshold failed")
    if not x_min**1.5 > builder.TERMINAL_RATE_CONSTANT:
        issues.append("independent nonterminal separation failed")

    ratio = (
        2 ** (-49 / 100) * (1 + d_upper) / (1 - d_upper)
    )
    if not ratio < 2 ** (-12 / 25):
        issues.append("independent corrected contraction failed")

    collar_scaled = 50_700 * (1 + math.pi / (8 * x_min))
    if not collar_scaled < builder.TERMINAL_WIDTH_CONSTANT:
        issues.append("independent terminal collar failed")

    tail = math.erfc(1 / math.sqrt(2)) / 2
    if not abs(tail - 0.15865525393145707) < 1e-15:
        issues.append("independent Gaussian-tail calibration failed")


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
        "terminal_rate": 16_900,
        "terminal_width": 50_701,
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
        "No circle, prime-power block, or prefix polygon defines",
        "nu_n=partial_x arg(z_n)=",
        "u_n>16900/x^2",
        "0<=x-x_N=4*pi*N^2*(exp(2u_N)-1)<50701/x",
        "p^(-12/25)",
        "Enestrom-Kakeya",
        "wind(Q_R,0)=0",
        "wind(P_R,0)=1",
        "external m-phase is factored out",
        "exp(t*log(n)^2/4)=E exp(sigma*Z*log(n))",
        "barPhi(1)=0.158655",
        "O_M^(p)(s_Z)",
        "\"joined_winding\": 3",
        "\"joined_quadratic_discriminant\": \"-19/2\"",
        "winding is not additive under vector addition",
        "0<=kappa_j<1",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove a uniform corrected-block homotopy",
        "joined Xi prefix lower bound",
        "strict composed successor flux bound",
        "q<1",
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
        "Absolute Carrier Rate",
        "Prime-Power Blocks",
        "One-Block Theorem",
        "Gaussian Heat Identity",
        "Two-Block Guard",
        "Joined Structure",
        "Live Theorem",
        "not an Xi counterexample",
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
        "validated Newman prime-power heat-block guard: "
        "14 rows, 1 absolute-rate collar, "
        "1 exact normalized-block zero-free/shifted-winding theorem, "
        "1 Gaussian-mixture audit, "
        "1 joined p-free decomposition, "
        "1 exact two-block winding-3 countermodel, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
