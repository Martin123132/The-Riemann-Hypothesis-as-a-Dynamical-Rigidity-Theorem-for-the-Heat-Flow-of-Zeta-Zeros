#!/usr/bin/env python3
"""Independently validate the Laguerre scale-mixture gate."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_laguerre_scale_mixture_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_laguerre_scale_mixture_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_laguerre_scale_mixture_gate.md"


EXPECTED_IDS = [
    "lsmg_01_jensen_moment_integral",
    "lsmg_02_kummer_kernel",
    "lsmg_03_laguerre_kernel",
    "lsmg_04_single_scale_hyperbolicity",
    "lsmg_05_positive_mixture_guard",
    "lsmg_06_gamma_reduction",
    "lsmg_07_half_integer_factorization",
    "lsmg_08_half_integer_all_degree_theorem",
    "lsmg_09_log_concavity_guard",
    "lsmg_10_xi_gamma_laguerre_handoff",
]


def independent_kernel(d: int, w: sp.Symbol, y: sp.Expr) -> sp.Expr:
    return sp.expand(
        sum(
            sp.factorial(d)
            / sp.factorial(d - j)
            * (w * y) ** j
            / sp.factorial(2 * j)
            for j in range(d + 1)
        )
    )


def independent_gamma(
    d: int, alpha: sp.Expr, theta: sp.Expr, w: sp.Symbol
) -> sp.Expr:
    return sp.expand(
        sum(
            sp.factorial(d)
            / sp.factorial(d - j)
            * sp.rf(alpha, j)
            * (theta * w) ** j
            / sp.factorial(2 * j)
            for j in range(d + 1)
        )
    )


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = gate.build_payload()
    if stored != rebuilt:
        issues.append("stored Laguerre scale-mixture payload differs from reconstruction")
    if stored.get("status") != "exact Laguerre scale-mixture reduction with sharp preservation guards":
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 gate rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    if sum(row.get("role") == "exact_countermodel" for row in rows) != 2:
        issues.append("expected two exact preservation countermodels")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("expected exactly one open Xi-specific handoff")

    w, y = sp.symbols("w y", real=True)
    for d in range(1, 11):
        direct = independent_kernel(d, w, y)
        kummer = sp.hyperexpand(
            sp.hyper([-d], [sp.Rational(1, 2)], -w * y / 4)
        )
        laguerre = (
            sp.factorial(d)
            / sp.rf(sp.Rational(1, 2), d)
            * sp.assoc_laguerre(d, sp.Rational(-1, 2), -w * y / 4)
        )
        if sp.simplify(direct - kummer) != 0:
            issues.append(f"independent Kummer identity failed at D={d}")
        if sp.simplify(direct - laguerre) != 0:
            issues.append(f"independent Laguerre identity failed at D={d}")

        fixed_scale = sp.Poly(direct.subs(y, sp.Rational(3, 2)), w, domain=sp.QQ)
        if fixed_scale.count_roots(-sp.oo, 0) != d:
            issues.append(f"fixed-scale negative-root count failed at D={d}")
        if sp.discriminant(fixed_scale.as_expr(), w) == 0:
            issues.append(f"fixed-scale roots were not simple at D={d}")

    # Recompute the minimal positive-mixture failure from its moments.
    mass = sp.Rational(9, 10)
    a, b = sp.Rational(1, 4), sp.Rational(5, 2)
    mixed = sp.expand(
        mass * independent_kernel(2, w, a)
        + (1 - mass) * independent_kernel(2, w, b)
    )
    if mixed != 1 + sp.Rational(19, 40) * w + sp.Rational(109, 1920) * w**2:
        issues.append("independent two-atom polynomial failed")
    if sp.factor(sp.discriminant(mixed, w)) != -sp.Rational(7, 4800):
        issues.append("independent positive-mixture discriminant failed")

    # Gamma moments and the terminating 2F1 expansion agree independently.
    theta = sp.symbols("theta", positive=True)
    alpha = sp.Rational(11, 7)
    for d in range(1, 8):
        direct = independent_gamma(d, alpha, theta, w)
        hyper = sp.hyperexpand(
            sp.hyper([-d, alpha], [sp.Rational(1, 2)], -theta * w / 4)
        )
        if sp.simplify(direct - hyper) != 0:
            issues.append(f"independent Gamma reduction failed at D={d}")

    # Check the two Jacobi branches and exact real-root locations over a wider grid.
    xvar = -w / 4
    for m in range(7):
        alpha_half = sp.Rational(2 * m + 1, 2)
        for d in range(1, 13):
            polynomial = independent_gamma(d, alpha_half, 1, w)
            if d >= m:
                cofactor_x = sp.expand(
                    sum(
                        sp.rf(-m, j)
                        * sp.rf(sp.Rational(2 * d + 1, 2), j)
                        * xvar**j
                        / (sp.rf(sp.Rational(1, 2), j) * sp.factorial(j))
                        for j in range(m + 1)
                    )
                )
                expected = sp.expand((1 - xvar) ** (d - m) * cofactor_x)
                endpoint_power = d - m
                if sp.simplify(polynomial - expected) != 0:
                    issues.append(f"Euler-Jacobi factorization failed at m={m}, D={d}")
                quotient = sp.cancel(polynomial / (w + 4) ** endpoint_power)
                quotient_poly = sp.Poly(quotient, w, domain=sp.QQ)
                if quotient_poly.degree() != m:
                    issues.append(f"half-integer cofactor degree failed at m={m}, D={d}")
                if m and quotient_poly.count_roots(-4, 0) != m:
                    issues.append(f"half-integer interior root count failed at m={m}, D={d}")
                if m and sp.discriminant(quotient_poly.as_expr(), w) == 0:
                    issues.append(f"half-integer cofactor simplicity failed at m={m}, D={d}")
            else:
                expected = (
                    sp.factorial(d)
                    / sp.rf(sp.Rational(1, 2), d)
                    * sp.jacobi(
                        d,
                        sp.Rational(-1, 2),
                        m - d,
                        1 - 2 * xvar,
                    )
                )
                if sp.simplify(polynomial - expected) != 0:
                    issues.append(f"direct Jacobi factorization failed at m={m}, D={d}")
                direct_poly = sp.Poly(polynomial, w, domain=sp.QQ)
                if direct_poly.count_roots(-4, 0) != d:
                    issues.append(f"direct Jacobi root count failed at m={m}, D={d}")
                if sp.discriminant(direct_poly.as_expr(), w) == 0:
                    issues.append(f"direct Jacobi root simplicity failed at m={m}, D={d}")

    # The exponential law is both log-concave and an exact D=3 failure.
    exponential = independent_gamma(3, 1, 1, w)
    expected_exponential = 1 + sp.Rational(3, 2) * w + w**2 / 2 + w**3 / 20
    if exponential != expected_exponential:
        issues.append("independent exponential-mixing polynomial failed")
    if sp.factor(sp.discriminant(exponential, w)) != -sp.Rational(1, 200):
        issues.append("independent log-concavity countermodel discriminant failed")

    checks = stored.get("exact", {}).get("checks", {})
    if len(checks.get("kernel_checks", [])) != 8:
        issues.append("stored kernel-check count drifted")
    if len(checks.get("gamma_checks", [])) != 6:
        issues.append("stored Gamma-check count drifted")
    if len(checks.get("half_integer_checks", [])) != 60:
        issues.append("stored half-integer-check count drifted")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of PF-infinity, RH, or `Lambda <= 0`",
        "every scale entering the integral is individually hyperbolic",
        "counterexample to positive-mixture preservation",
        "log-concavity is insufficient",
        "genuine all-degree hyperbolic benchmark family",
        "positive mixing of these blocks is not enough",
        "Xi-specific theorem",
        "https://dlmf.nist.gov/18.16",
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
        "validated Jensen-window PF Laguerre scale-mixture gate: "
        "10 rows, 0 issues, 3 exact kernel identities, "
        "1 individual-kernel hyperbolicity theorem, 1 positive-mixture countermodel, "
        "1 Gamma reduction, 1 half-integer all-degree theorem, "
        "1 log-concavity countermodel, 1 open Xi-specific handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
