#!/usr/bin/env python3
"""Validate the joined dyadic odd-prefix first-jet reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_joined_dyadic_odd_prefix_first_jet_reduction as builder


EXPECTED_IDS = [
    "jdopfj_00_pi_provenance",
    "jdopfj_01_coefficient_chain",
    "jdopfj_02_heat_shift",
    "jdopfj_03_joined_odd_prefix",
    "jdopfj_04_correction_current",
    "jdopfj_05_phase_lift",
    "jdopfj_06_five_current_closure",
    "jdopfj_07_endpoint_jet",
    "jdopfj_08_cutoff_jump",
    "jdopfj_09_phase_frozen_theorem",
    "jdopfj_10_phase_frozen_margin",
    "jdopfj_11_homotopy_target",
    "jdopfj_12_averaging_guard",
    "jdopfj_13_second_moment_guard",
    "jdopfj_14_route_decision",
    "jdopfj_15_q_ge_1_target",
    "jdopfj_16_q_lt_1_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    time, h, k, log_m, spectral = sp.symbols(
        "time h k log_m spectral"
    )
    lhs = time * (k * h + log_m) ** 2 / 4 - spectral * (
        k * h + log_m
    )
    rhs = (
        time * (k * h) ** 2 / 4
        - spectral * k * h
        + time * log_m**2 / 4
        - (spectral - time * k * h / 2) * log_m
    )
    if sp.expand(lhs - rhs) != 0:
        issues.append("independent deterministic heat shift failed")

    z, xi = sp.symbols("z xi")
    integer_k = sp.symbols("integer_k", integer=True, nonnegative=True)
    monomial = z**integer_k * sp.exp(sp.I * xi * log_m)

    def operator(value: sp.Expr) -> sp.Expr:
        return h * z * sp.diff(value, z) - sp.I * sp.diff(value, xi)

    expected_log = integer_k * h + log_m
    if sp.simplify(operator(monomial) - expected_log * monomial) != 0:
        issues.append("independent first moment operator failed")
    if sp.simplify(
        operator(operator(monomial)) - expected_log**2 * monomial
    ) != 0:
        issues.append("independent second moment operator failed")

    s_1, log_a, u_x = sp.symbols("s_1 log_a u_x")
    h_0, h_1, h_2, d_0, d_1 = sp.symbols(
        "H_0 H_1 H_2 D_0 D_1"
    )
    h_0_x = -s_1 * h_1 + d_0
    h_1_x = -s_1 * h_2 + d_1
    direct = h_1_x - u_x * h_0 - log_a * h_0_x
    reduced = (
        -s_1 * (h_2 - log_a * h_1)
        + d_1
        - log_a * d_0
        - u_x * h_0
    )
    if sp.expand(direct - reduced) != 0:
        issues.append("independent five-current closure failed")

    mu, d_n, d_anchor = sp.symbols("mu d_n d_anchor")
    ratio = (1 + mu * d_n) / (1 + mu * d_anchor)
    expected_mu = (d_n - d_anchor) / (1 + mu * d_anchor) ** 2
    if sp.simplify(sp.diff(ratio, mu) - expected_mu) != 0:
        issues.append("independent correction homotopy failed")

    moment_weights = (sp.Integer(1), sp.Integer(-2), sp.Integer(1))
    moment_logs = (sp.Integer(0), h, 2 * h)
    moment_values = [
        sp.expand(sum(
            weight * log_value**order
            for weight, log_value in zip(
                moment_weights, moment_logs, strict=True
            )
        ))
        for order in range(3)
    ]
    if moment_values != [0, 0, 2 * h**2]:
        issues.append("independent second-moment guard failed")

    r = sp.sqrt(2) / 2
    first_fibre = z * (1 + r * z + r**2 * z**2)
    second_fibre = -sp.Rational(4, 5) * z * (1 + r * z)
    mean_inner = sp.expand(10 * (first_fibre + second_fibre) / (2 * z))
    if sp.discriminant(mean_inner, z) != -sp.Rational(19, 2):
        issues.append("independent two-fibre discriminant failed")
    coefficients = sp.Poly(mean_inner, z).all_coeffs()
    if sp.simplify(coefficients[-1] / coefficients[0]) != sp.Rational(
        2, 5
    ):
        issues.append("independent two-fibre root modulus failed")


def independent_valuation_audit(issues: list[str]) -> None:
    n_max = 128
    represented: list[int] = []
    for layer in range(8):
        bound = n_max // (2**layer)
        represented.extend(
            (2**layer) * odd
            for odd in range(1, bound + 1, 2)
        )
    if sorted(represented) != list(range(1, n_max + 1)):
        issues.append("independent dyadic valuation reconstruction failed")

    for n in range(2, n_max + 1):
        layer = 0
        odd = n
        while odd % 2 == 0:
            odd //= 2
            layer += 1
        previous = n - 1
        old_bound = previous // (2**layer)
        new_bound = n // (2**layer)
        if not (
            new_bound == odd
            and new_bound == old_bound + 1
            and odd % 2 == 1
        ):
            issues.append(f"cutoff layer entry failed at n={n}")
            break


def independent_phase_frozen_audit(issues: list[str]) -> None:
    z = sp.symbols("z")
    coefficients = [
        sp.Integer(1),
        sp.Rational(1, 2),
        sp.Rational(1, 4),
        sp.Rational(1, 8),
    ]
    polynomial = sum(
        coefficient * z**index
        for index, coefficient in enumerate(coefficients)
    )
    roots = sp.nroots(polynomial, n=40, maxsteps=100)
    minimum_modulus = min(abs(complex(root)) for root in roots)
    if not minimum_modulus > 1.999999999999:
        issues.append("independent phase-frozen root-radius audit failed")
    if not all(
        coefficients[index] > coefficients[index + 1]
        for index in range(len(coefficients) - 1)
    ):
        issues.append("independent decreasing-coefficient audit failed")


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
        "uncorrected_dyadic_ratio": "2^(-49/100)",
        "corrected_dyadic_ratio": "2^(-12/25)",
    }:
        issues.append("constant payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_xi_rows = [
        row for row in rows if row.get("role") == "open_theorem_target"
    ]
    if len(open_xi_rows) != 2:
        issues.append(f"open Xi target count drifted: {len(open_xi_rows)}")
    countermodels = [
        row for row in rows if row.get("role") == "countermodel"
    ]
    if len(countermodels) != 2:
        issues.append(f"route-guard count drifted: {len(countermodels)}")
    hypothetical = [
        row
        for row in rows
        if row.get("role") == "hypothetical_exact_theorem"
    ]
    if len(hypothetical) != 1:
        issues.append(
            f"phase-frozen theorem count drifted: {len(hypothetical)}"
        )

    independent_symbolic_audit(issues)
    independent_valuation_audit(issues)
    independent_phase_frozen_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "do not define pi",
        "A_t(2^k*m;s)=A_t(2^k;s)",
        "M_k=floor(N/2^k)",
        "H_(r,x)=-s_*'*H_(r+1)+D_r",
        "H_0,H_1,H_2,D_0,D_1",
        "D_h=h*z*partial_z-i*partial_xi",
        "P_N(z_*,omega_*,1)=H_0",
        "Z_(A,x)=r_(A,x)-s_*''P_1-s_*'P_(1,x)",
        "Delta Z_0=q_n+j_0=Q_N/f_1",
        "0<B_(k+1)<2^(-49/100)B_k",
        "B_K*(R_0-1)^K",
        '"mean_winding": 3',
        '"moments": "H_0=0, H_1=0, H_2=2h^2"',
        "0<=kappa_j<1",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "actual odd-phase/correction homotopy",
        "joined Xi lower bound",
        "strict successor flux upper bound",
        "q<1",
        "contact exclusion",
        "Lambda<=0",
        "PF-infinity",
        "RH",
        "Clay-prize",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pi Provenance",
        "Deterministic Heat Shift",
        "Joined Odd Prefixes",
        "Phase Lift",
        "Five-Current First Jet",
        "Cutoff And Endpoint",
        "Phase-Synchronized Hypothetical",
        "Homotopy Target",
        "Averaging does not preserve winding",
        "second logarithmic moment is indispensable",
        "not Xi counterexamples",
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
        "validated Newman joined dyadic odd-prefix first-jet reduction: "
        "17 rows, 1 exact heat-shift factorization, "
        "1 five-current closure, "
        "1 phase-frozen joined zero-free theorem, "
        "2 exact route guards, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
