#!/usr/bin/env python3
"""Independently validate the Newman theta-summand spectral-square gate."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_theta_summand_spectral_square_gate as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_theta_summand_spectral_square_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md"
)

EXPECTED_IDS = [
    "ntssg_01_summand_shift",
    "ntssg_02_differential_profile",
    "ntssg_03_theta_boundary",
    "ntssg_04_half_transform",
    "ntssg_04b_deformed_half_transform",
    "ntssg_05_mellin_transform",
    "ntssg_06_xi_reconstruction",
    "ntssg_07_laguerre_pair_expansion",
    "ntssg_08_finite_boundary_jet",
    "ntssg_09_finite_laguerre_obstruction",
    "ntssg_10_summand_sign_witnesses",
    "ntssg_11_infinite_cancellation_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]
    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored theta-summand payload differs from reconstruction")
    if stored.get("status") != (
        "exact theta differential/Mellin reduction with finite-summand spectral obstruction"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 12:
        issues.append("expected 12 rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    for role, expected in [
        ("exact_identity", 7),
        ("nonpromotion_gate", 1),
        ("exact_theorem", 2),
        ("numerical_counter_witness", 1),
        ("open_handoff", 1),
    ]:
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    exact = stored.get("exact", {})

    # Re-derive phi_1=(D^2-1)h/8 directly.
    u = sp.symbols("u", real=True)
    h = sp.exp(u - sp.pi * sp.exp(4 * u))
    phi = (
        2 * sp.pi**2 * sp.exp(9 * u) - 3 * sp.pi * sp.exp(5 * u)
    ) * sp.exp(-sp.pi * sp.exp(4 * u))
    if sp.simplify((sp.diff(h, u, 2) - h) / 8 - phi) != 0:
        issues.append("differential profile failed")
    if exact.get("differential_profile", {}).get("identity") != "phi_1(u)=(h''(u)-h(u))/8":
        issues.append("stored differential profile drifted")

    # Re-derive R'(0) from theta'(1)=-theta(1)/4.
    theta, theta_prime, s_value = sp.symbols("theta theta_prime S", real=True)
    theta_prime_value = -theta / 4
    r_prime = s_value + 2 * theta_prime_value
    if sp.simplify(r_prime.subs(theta, 1 + 2 * s_value) + sp.Rational(1, 2)) != 0:
        issues.append("theta primitive boundary failed")
    boundary = exact.get("theta_boundary", {})
    if boundary.get("primitive_boundary") != "R'(0)=-1/2":
        issues.append("stored R boundary drifted")
    if "every j>=0" not in boundary.get("evenness_cancellation", ""):
        issues.append("infinite odd-jet cancellation missing")

    # Verify the half-transform boundary normalization.
    if exact.get("half_transform", {}).get("identity") != (
        "H_0(x)=1/16-(1+x^2)*C(x)/8"
    ):
        issues.append("half-transform formula drifted")
    if "-R'(0)-x^2*C(x)" not in exact.get("half_transform", {}).get("boundary_term", ""):
        issues.append("half-transform integration-by-parts term missing")

    # Re-derive the all-t endpoint-subtracted transform.
    t, x = sp.symbols("t x", real=True)
    w = sp.exp(t * u**2) * sp.cos(x * u)
    expected_wpp_minus_w = sp.exp(t * u**2) * (
        (4 * t**2 * u**2 + 2 * t - 1 - x**2) * sp.cos(x * u)
        - 4 * t * u * x * sp.sin(x * u)
    )
    if sp.simplify(sp.diff(w, u, 2) - w - expected_wpp_minus_w) != 0:
        issues.append("deformed integration-by-parts multiplier failed")
    deformed = exact.get("deformed_half_transform", {})
    if deformed.get("identity") != "H_t(x)=1/16+D_t[C_t](x)/8":
        issues.append("deformed half-transform identity drifted")
    if "-4*t^2*C_t''(x)" not in deformed.get("expanded", ""):
        issues.append("deformed second-derivative coefficient missing")
    if "4*t*x*C_t'(x)" not in deformed.get("expanded", ""):
        issues.append("deformed first-derivative coefficient missing")
    a0, ap, app = sp.symbols("A Ap App", real=True)
    h0 = sp.Rational(1, 16) + a0 / 8
    hp = ap / 8
    hpp = app / 8
    reduced_laguerre = sp.simplify(hp**2 - h0 * hpp)
    expected_reduction = (ap**2 - (a0 + sp.Rational(1, 2)) * app) / 64
    if sp.simplify(reduced_laguerre - expected_reduction) != 0:
        issues.append("deformed curvature reduction failed")
    if "A_t(c)=-1/2" not in deformed.get("multiple_contact", ""):
        issues.append("deformed multiple-contact criterion missing")

    # Re-derive the Gamma factor and completed-xi normalization algebraically.
    z = sp.symbols("z")
    alpha = (5 + sp.I * z) / 4
    raw_transform = sp.pi ** (-(1 + sp.I * z) / 4) * (
        sp.Rational(1, 2) * sp.gamma(alpha + 1)
        - sp.Rational(3, 4) * sp.gamma(alpha)
    )
    reduced = -(
        1 + z**2
    ) / 32 * sp.pi ** (-(1 + sp.I * z) / 4) * sp.gamma((1 + sp.I * z) / 4)
    if sp.simplify(sp.expand_func(raw_transform) - sp.expand_func(reduced)) != 0:
        issues.append("bilateral Gamma transform failed")
    spectral_s = (1 + sp.I * z) / 2
    xi_prefactor = sp.Rational(1, 2) * spectral_s * (spectral_s - 1)
    if sp.simplify(xi_prefactor / 4 + (1 + z**2) / 32) != 0:
        issues.append("xi/4 prefactor failed")
    xi_reconstruction = exact.get("xi_reconstruction", {})
    if xi_reconstruction.get("identity") != "Qhat(z)*zeta(s)=xi(s)/4":
        issues.append("stored xi reconstruction drifted")
    if "not a new zero-free factor" not in xi_reconstruction.get("scope", ""):
        issues.append("xi reconstruction non-promotion boundary missing")

    # Re-derive the quadratic self/cross expansion for two abstract functions.
    f, fp, fpp, g, gp, gpp = sp.symbols("f fp fpp g gp gpp")
    direct = (fp + gp) ** 2 - (f + g) * (fpp + gpp)
    split = (fp**2 - f * fpp) + (gp**2 - g * gpp) + (
        2 * fp * gp - f * gpp - g * fpp
    )
    if sp.expand(direct - split) != 0:
        issues.append("Laguerre self/cross expansion failed")

    # Prove the finite endpoint jet sign pattern.
    q = sp.symbols("q", real=True)
    derivative_polynomial = -8 * q**2 + 30 * q - 15
    roots = sp.solve(derivative_polynomial, q)
    lower = (15 - sp.sqrt(105)) / 8
    upper = (15 + sp.sqrt(105)) / 8
    if roots != [lower, upper]:
        issues.append("summand derivative roots failed")
    if not (sp.N(lower) < sp.N(sp.pi) < sp.N(upper) < sp.N(4 * sp.pi)):
        issues.append("summand derivative sign ordering failed")
    defect = exact.get("finite_theta_boundary_defect", {})
    for phrase in ("phi_1'(0)>0", "phi_n'(0)<0", "Phi'(0)=0"):
        if phrase not in defect.get("signs", ""):
            issues.append(f"finite endpoint sign proof missing: {phrase}")
    if ">0" not in defect.get("partial_boundary_jet", ""):
        issues.append("positive finite endpoint jet missing")

    # Re-derive the leading negative Laguerre tail.
    x, a, b = sp.symbols("x a b", positive=True)
    asymptotic_h = -a / x**2 + b / x**4
    asymptotic_l = sp.expand(
        sp.diff(asymptotic_h, x) ** 2
        - asymptotic_h * sp.diff(asymptotic_h, x, 2)
    )
    if asymptotic_l.coeff(x, -6) != -2 * a**2:
        issues.append("finite theta Laguerre asymptotic failed")
    obstruction = exact.get("finite_theta_laguerre_obstruction", {})
    if "-2*A_N^2/x^6" not in obstruction.get("laguerre_tail", ""):
        issues.append("stored finite theta obstruction drifted")
    if "Every finite theta truncation" not in obstruction.get("conclusion", ""):
        issues.append("finite truncation conclusion missing")

    # Compare fast quadrature to independent 60-digit references.
    numerics = stored.get("numerics", {})
    reference = numerics.get("independent_mpmath_60_digit_references", {})
    by_key: dict[str, mp.mpf] = {}
    for row in numerics.get("rows", []):
        t_key = "0" if row.get("t") == "0" else "1/2"
        x_key = row.get("x")
        if row.get("kind") == "cross_1_2":
            key = f"t={t_key},x={x_key},cross12"
            value = mp.mpf(row.get("cross_laguerre", "nan"))
        else:
            key = f"t={t_key},x={x_key},self1"
            value = mp.mpf(row.get("self_laguerre", "nan"))
        by_key[key] = value
        full_key = f"t={t_key},x={x_key},full"
        if full_key in reference:
            by_key[full_key] = mp.mpf(row.get("full_laguerre", "nan"))
        if value >= 0:
            issues.append(f"negative summand witness failed: {key}")
        if mp.mpf(row.get("full_laguerre", "nan")) <= 0:
            issues.append(f"full sampled Laguerre sign failed: {full_key}")
    for key, expected_text in reference.items():
        if key not in by_key:
            issues.append(f"missing numerical witness key: {key}")
            continue
        expected = mp.mpf(expected_text)
        if abs(by_key[key] - expected) > mp.mpf("5e-8") * abs(expected):
            issues.append(f"independent numerical mismatch: {key}")
    if mp.mpf(numerics.get("max_relative_coarse_fine_delta_negative_witnesses", "inf")) >= mp.mpf("1e-8"):
        issues.append("summand quadrature convergence margin failed")
    if "Numerical witnesses only" not in numerics.get("scope", ""):
        issues.append("numerical scope boundary missing")

    decision = exact.get("nonpromotion_decision", "")
    for phrase in ("finite positive theta grouping", "infinite theta/modular"):
        if phrase not in decision:
            issues.append(f"nonpromotion decision missing: {phrase}")
    handoff = exact.get("open_handoff", "")
    for phrase in ("A_t'(x)^2", "odd-endpoint cancellation", "noncircular curvature", "closed routes"):
        if phrase not in handoff:
            issues.append(f"open handoff missing: {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "Theta Primitive",
        "Mellin Audit",
        "Qhat(z)*zeta(s)=xi(s)/4",
        "Finite-Sum Theorem",
        "every finite truncation fails",
        "entrywise sign is false",
        "infinite odd-endpoint cancellation",
    ]
    for phrase in required:
        if phrase not in note:
            issues.append(f"note missing required phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Jensen-window PF Newman theta-summand spectral-square gate: "
        "12 rows, 0 issues, 7 exact transform identities, "
        "1 xi-reconstruction non-promotion gate, 2 exact finite-truncation theorems, "
        "1 numerical sign diagnostic, 1 infinite-cancellation handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
