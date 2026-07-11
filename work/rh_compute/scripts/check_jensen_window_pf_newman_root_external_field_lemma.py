#!/usr/bin/env python3
"""Independently validate the Newman root external-field lemma."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_root_external_field_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_root_external_field_lemma.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_root_external_field_lemma.md"


EXPECTED_IDS = [
    "nref_01_canonical_product",
    "nref_02_squared_external_field",
    "nref_03_squared_stiffness",
    "nref_04_newman_field_conversion",
    "nref_05_newman_stiffness_conversion",
    "nref_06_pair_flow_reduction",
    "nref_07_gap_stiffness_expansion",
    "nref_08_external_field_sign_guard",
    "nref_09_cosine_equilibrium",
    "nref_10_xi_balance_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = lemma.build_payload()
    if stored != rebuilt:
        issues.append("stored external-field payload differs from reconstruction")
    if stored.get("status") != "exact Newman root external-field reduction with sign countermodels":
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 theorem rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("expected exactly one open handoff")

    # Rebuild the local F(s) and H(z)=F(-z^2) jets independently.
    s, z, x = sp.symbols("s z x", positive=True)
    a0, a1, a2 = sp.symbols("a0 a1 a2", nonzero=True)
    rho = -x**2
    delta = s - rho
    f = sp.expand(delta**2 * (a0 + a1 * delta + a2 * delta**2))
    h = sp.expand(f.subs(s, -z**2))
    squared_field = a1 / a0
    signed_field = 1 / x - 2 * x * squared_field
    if sp.simplify(sp.diff(f, s, 3).subs(s, rho) / sp.diff(f, s, 2).subs(s, rho) - 3 * squared_field) != 0:
        issues.append("independent squared-variable jet identity failed")
    if sp.simplify(sp.diff(h, z, 3).subs(z, x) / sp.diff(h, z, 2).subs(z, x) - 3 * signed_field) != 0:
        issues.append("independent signed-variable jet identity failed")

    # Verify the field and stiffness conversions on several exact root systems.
    for xv, roots, multiplicities in [
        (sp.Rational(1), [sp.Rational(2), sp.Rational(4)], [1, 3]),
        (sp.Rational(3, 2), [sp.Rational(1, 2), sp.Rational(5, 2)], [2, 1]),
        (sp.Rational(2), [sp.Rational(1), sp.Rational(3), sp.Rational(5)], [1, 2, 1]),
    ]:
        e_s = sp.factor(sum(m / (r**2 - xv**2) for r, m in zip(roots, multiplicities)))
        k_s = sp.factor(sum(m / (r**2 - xv**2) ** 2 for r, m in zip(roots, multiplicities)))
        signed: list[tuple[sp.Rational, int]] = [(-xv, 2)]
        for root, multiplicity in zip(roots, multiplicities):
            signed.extend([(root, multiplicity), (-root, multiplicity)])
        b_x = sp.factor(sum(m / (xv - r) for r, m in signed))
        k_x = sp.factor(sum(m / (xv - r) ** 2 for r, m in signed))
        if sp.simplify(b_x - (1 / xv - 2 * xv * e_s)) != 0:
            issues.append(f"field conversion failed independently at x={xv}")
        if sp.simplify(k_x - (1 / (2 * xv**2) + 2 * e_s + 4 * xv**2 * k_s)) != 0:
            issues.append(f"stiffness conversion failed independently at x={xv}")
        if k_x <= 0:
            issues.append(f"signed stiffness was not positive at x={xv}")

    # Derive center and squared-gap equations directly from the zero ODE.
    c, g = sp.symbols("c g", real=True, nonzero=True)
    outsiders = [sp.Rational(-7, 3), sp.Rational(1, 4), sp.Rational(11, 3)]
    xp, xm = c + g / 2, c - g / 2
    xp_dot = 2 / g + 2 * sum(1 / (xp - root) for root in outsiders)
    xm_dot = -2 / g + 2 * sum(1 / (xm - root) for root in outsiders)
    center_dot = sp.factor((xp_dot + xm_dot) / 2)
    gap_dot = sp.factor(xp_dot - xm_dot)
    pair_sum = sum(1 / ((xp - root) * (xm - root)) for root in outsiders)
    if sp.simplify(gap_dot - (4 / g - 2 * g * pair_sum)) != 0:
        issues.append("independent gap equation failed")
    if sp.simplify(center_dot.subs(g, 0) - 2 * sum(1 / (c - root) for root in outsiders)) != 0:
        issues.append("independent center-field limit failed")

    tau, stiffness = sp.symbols("tau stiffness", positive=True)
    q_trial = 8 * tau - 16 * stiffness * tau**2
    # q'=8-4*K*q+O(q^2,tau*q) fixes the first two coefficients.
    residual = sp.expand(sp.diff(q_trial, tau) - (8 - 4 * stiffness * q_trial))
    if sp.expand(residual).coeff(tau, 0) != 0 or sp.expand(residual).coeff(tau, 1) != 0:
        issues.append("independent squared-gap stiffness expansion failed")

    # Both signs occur inside finite positive-coefficient LP products.
    for polynomial, expected_field, expected_shift in [
        ((1 + s) ** 2 * (1 + s / 2), sp.Rational(1), sp.Rational(-3, 64)),
        ((1 + s) ** 2 * (1 + 2 * s), sp.Rational(-2), sp.Rational(9, 64)),
    ]:
        u = sp.cancel(polynomial / (s + 1) ** 2)
        field = sp.simplify(sp.diff(u, s).subs(s, -1) / u.subs(s, -1))
        shift = sp.simplify(sp.Rational(1, 64) - field / 16)
        if field != expected_field or shift != expected_shift:
            issues.append(f"external-field sign model failed for {sp.expand(polynomial)}")

    # The cosine lattice is the exact zero-center-field equilibrium.
    y = sp.symbols("y", real=True)
    kernel = sp.cosh(sp.sqrt(y))
    ode = sp.simplify(4 * y * sp.diff(kernel, y, 2) + 2 * sp.diff(kernel, y) - kernel)
    if ode != 0:
        issues.append("independent cosine-kernel ODE failed")
    rho_symbol = sp.symbols("rho_symbol", negative=True)
    equilibrium_field = -1 / (2 * rho_symbol)
    x_symbol = sp.sqrt(-rho_symbol)
    if sp.simplify(1 / x_symbol - 2 * x_symbol * equilibrium_field) != 0:
        issues.append("independent cosine equilibrium field failed")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "universal `8*(t-t_*)` square-root birth",
        "neither LP membership nor coefficient positivity fixes the field",
        "exact arithmetic-lattice equilibrium",
        "escape to infinite height",
        "https://arxiv.org/abs/1801.05914",
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
        "validated Jensen-window PF Newman root external-field lemma: "
        "10 rows, 0 issues, 5 exact canonical-product identities, "
        "1 pair-flow reduction, 1 gap-stiffness expansion, 2 sign countermodels, "
        "1 cosine equilibrium benchmark, 1 open Xi-balance handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
