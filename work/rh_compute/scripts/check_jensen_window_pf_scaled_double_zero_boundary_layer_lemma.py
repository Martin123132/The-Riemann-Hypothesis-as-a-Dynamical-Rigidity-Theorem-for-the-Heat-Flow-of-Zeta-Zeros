#!/usr/bin/env python3
"""Independently validate the scaled double-zero boundary-layer lemma."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_scaled_double_zero_boundary_layer_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_scaled_double_zero_boundary_layer_lemma.md"


EXPECTED_IDS = [
    "sdzb_01_scaled_jensen_identity",
    "sdzb_02_first_degree_correction",
    "sdzb_03_coefficient_heat_pde",
    "sdzb_04_double_zero_transversality",
    "sdzb_05_universal_boundary_layer",
    "sdzb_06_root_gap_law",
    "sdzb_07_finite_collision_shift",
    "sdzb_08_exact_toy_validation",
    "sdzb_09_eventual_threshold_exhaustion",
    "sdzb_10_uniform_layer_handoff",
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
        issues.append("stored boundary-layer payload differs from reconstruction")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 theorem rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")

    exact = stored.get("exact", {})
    checks = exact.get("checks", {})
    if checks.get("falling_factor_orders") != [0, 12]:
        issues.append("falling-factor order range drifted")
    if len(checks.get("falling_factor_checks", [])) != 13:
        issues.append("falling-factor check count drifted")

    h = sp.symbols("h")
    for order in range(17):
        product = sp.expand(sp.prod(1 - offset * h for offset in range(order)))
        if product.coeff(h, 1) != -sp.Rational(order * (order - 1), 2):
            issues.append(f"independent falling correction failed at j={order}")
        expected_quadratic = sp.Rational(
            order * (order - 1) * (order - 2) * (3 * order - 1), 24
        )
        if product.coeff(h, 2) != expected_quadratic:
            issues.append(f"independent second falling correction failed at j={order}")

    lam, z = sp.symbols("lambda z", real=True)
    toy = z**2 + (2 + 12 * lam) * z + 1 + 4 * lam + 12 * lam**2
    if sp.factor(sp.diff(toy, lam) - 2 * sp.diff(toy, z) - 4 * z * sp.diff(toy, z, 2)) != 0:
        issues.append("independent toy PDE failed")

    D = sp.symbols("D", positive=True)
    scaled = 1 + 4 * lam + 12 * lam**2 + (2 + 12 * lam) * z + (1 - 1 / D) * z**2
    disc = sp.factor(sp.discriminant(scaled, z))
    near = (
        sp.sqrt(2) * sp.sqrt(D - 1) - sp.sqrt(2 * D + 1)
    ) / (6 * sp.sqrt(2 * D + 1))
    if sp.simplify(disc.subs(lam, near)) != 0:
        issues.append("independent near-threshold check failed")
    series = sp.series(near.subs(D, 1 / h), h, 0, 4).removeO()
    if sp.simplify(series - (-h / 8 + h**2 / 64 - 5 * h**3 / 256)) != 0:
        issues.append("independent threshold expansion failed")

    rho, tau, eta = sp.symbols("rho tau eta", nonzero=True, real=True)
    quadratic = eta**2 + 8 * rho * tau - rho**2
    if sp.simplify(
        sp.discriminant(quadratic, eta) - 4 * (rho**2 - 8 * rho * tau)
    ) != 0:
        issues.append("independent local discriminant failed")
    if sp.simplify(quadratic.subs({tau: rho / 8, eta: 0})) != 0:
        issues.append("independent rho/(8D) collision check failed")

    n, f2, f3 = sp.symbols("n f2 f3", nonzero=True, real=True)
    expected_second = -rho * (6 * n + 1) / 64 - rho**2 * f3 / (48 * f2)
    recorded_second = exact.get("checks", {}).get("second_parameter_shift", "")
    parsed_second = sp.sympify(
        recorded_second, locals={"rho": rho, "n": n, "f2": f2, "f3": f3}
    )
    if sp.simplify(parsed_second - expected_second) != 0:
        issues.append("independent external-field collision correction failed")

    thresholds = exact.get("eventual_thresholds", {})
    if thresholds.get("exhaustion") != "sup_D Theta_D=Lambda":
        issues.append("threshold-exhaustion statement drifted")
    if "full forward ray" not in thresholds.get("warning", ""):
        issues.append("forward-ray warning missing")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of PF-infinity, RH, or",
        "finite-degree collision lies on the bad side",
        "global root external field first enters",
        "A single endpoint test",
        "degree- and",
        "height-uniform bound",
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
        "validated Jensen-window PF scaled double-zero boundary-layer lemma: "
        "10 rows, 0 issues, 3 exact scaling identities, 1 heat PDE, "
        "1 double-zero transversality, 1 universal boundary layer, "
        "1 D^(-3/2) gap law, 1 external-field D^(-2) collision law, 1 exact toy family, "
        "1 threshold-exhaustion theorem, 1 open uniform handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
