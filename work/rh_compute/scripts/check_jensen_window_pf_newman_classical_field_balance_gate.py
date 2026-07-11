#!/usr/bin/env python3
"""Independently validate the Newman classical-field balance gate."""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_classical_field_balance_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_classical_field_balance_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_classical_field_balance_gate.md"


EXPECTED_IDS = [
    "ncfb_01_principal_value_field",
    "ncfb_02_arithmetic_pair_equilibrium",
    "ncfb_03_weighted_perturbation_identity",
    "ncfb_04_classical_continuum_field",
    "ncfb_05_positive_time_quantile_drift",
    "ncfb_06_published_high_zero_localization",
    "ncfb_07_positive_boundary_compactness",
    "ncfb_08_macroscopic_location_guard",
    "ncfb_09_required_local_upgrade",
    "ncfb_10_lambda_uniform_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = gate.build_payload()
    if stored != rebuilt:
        issues.append("stored classical-field payload differs from reconstruction")
    if stored.get("status") != "exact classical-field balance reduction with fixed-time compactness gate":
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 gate rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    if sum(row.get("role") == "open_handoff" for row in rows) != 1:
        issues.append("expected exactly one open lambda-uniform handoff")
    if sum(row.get("role") == "exact_countermodel" for row in rows) != 1:
        issues.append("expected one countermodel row containing both field signs")

    # Recheck the one-root perturbation identity independently.
    c0, dc, root, dr = sp.symbols("c0 dc root dr", real=True)
    lhs = 1 / (c0 + dc - root - dr) - 1 / (c0 - root)
    rhs = (dr - dc) / ((c0 + dc - root - dr) * (c0 - root))
    if sp.simplify(lhs - rhs) != 0:
        issues.append("independent field perturbation identity failed")
    for values in [
        (sp.Rational(2), sp.Rational(1, 7), sp.Rational(-3), sp.Rational(1, 5)),
        (sp.Rational(-1), sp.Rational(2, 9), sp.Rational(4), sp.Rational(-1, 6)),
    ]:
        substitution = dict(zip((c0, dc, root, dr), values))
        if sp.simplify((lhs - rhs).subs(substitution)) != 0:
            issues.append(f"rational perturbation check failed at {values}")

    # Every finite reflection-symmetric arithmetic truncation cancels exactly.
    h = sp.Rational(7, 5)
    for cutoff in range(1, 31):
        distances = [sp.Rational(2 * k + 1, 2) * h for k in range(1, cutoff + 1)]
        field = sum(1 / (-distance) + 1 / distance for distance in distances)
        if sp.simplify(field) != 0:
            issues.append(f"arithmetic pair equilibrium failed at cutoff={cutoff}")

    # Verify the principal-value constant through its two half-line series.
    k = sp.symbols("k", integer=True, nonnegative=True)
    odd_sum = sp.summation(1 / (2 * k + 1) ** 2, (k, 0, sp.oo))
    if sp.simplify(odd_sum - sp.pi**2 / 8) != 0:
        issues.append("independent odd-square sum failed")
    mp.mp.dps = 60
    half_integral = mp.quad(lambda u: mp.log(u) / (1 - u * u), [0, 1])
    if abs(half_integral + mp.pi**2 / 8) > mp.mpf("1e-45"):
        issues.append("independent logarithmic half-integral failed")

    # Check the exact continuum formula and its first asymptotic correction.
    def continuum_formula(a_value: mp.mpf, x_value: mp.mpf) -> mp.mpf:
        r_value = a_value / x_value
        j_value = mp.mpf("0.5") * (
            mp.log(r_value) * mp.log((1 + r_value) / (1 - r_value))
            + mp.polylog(2, -r_value)
            - mp.polylog(2, r_value)
        )
        return (
            -mp.log(x_value / (4 * mp.pi)) * mp.atanh(r_value)
            - mp.pi**2 / 4
            - j_value
        ) / (2 * mp.pi)

    a_value = mp.mpf(20)
    expected_scaled = (1 - mp.log(a_value / (4 * mp.pi))) / (2 * mp.pi)
    previous_error = None
    for x_value in (mp.mpf("1e4"), mp.mpf("1e6"), mp.mpf("1e8")):
        field = continuum_formula(a_value, x_value)
        scaled = (field + mp.pi / 8) * x_value / a_value
        error = abs(scaled - expected_scaled)
        if previous_error is not None and error >= previous_error:
            issues.append("continuum first-correction convergence was not improving")
        previous_error = error
    if previous_error is None or previous_error > mp.mpf("1e-12"):
        issues.append("continuum -pi/8 asymptotic check failed")

    # Recompute the positive-time quantile velocity.
    level, time = sp.symbols("level time", positive=True)
    log_density = sp.log(level / (4 * sp.pi))
    velocity = sp.factor(
        -(log_density / 16)
        / (log_density / (4 * sp.pi) + time / (16 * level))
    )
    expected_velocity = -sp.pi * level * log_density / (
        4 * level * log_density + sp.pi * time
    )
    if sp.simplify(velocity - expected_velocity) != 0:
        issues.append("independent quantile velocity identity failed")
    if sp.limit(velocity, level, sp.oo) != -sp.pi / 4:
        issues.append("independent quantile velocity limit failed")

    # Bounded even displacements can force either sign of an unbounded field.
    eps = sp.symbols("eps", positive=True)
    n = sp.Rational(5)
    baseline = 1 / (2 * n)
    upper = baseline - 1 / eps + 1 + 1 / (2 * n + eps) - 1 / (2 * n + 1)
    lower = baseline + 1 / eps - 1 + 1 / (2 * n - eps) - 1 / (2 * n - 1)
    if sp.limit(eps * upper, eps, 0, dir="+") != -1:
        issues.append("independent negative-field sensitivity model failed")
    if sp.limit(eps * lower, eps, 0, dir="+") != 1:
        issues.append("independent positive-field sensitivity model failed")

    published = stored.get("exact", {}).get("published_fixed_time_localization", {})
    if published.get("source") != "https://arxiv.org/abs/1904.12438":
        issues.append("published localization source drifted")
    if "exp(C/t)" not in published.get("theorem", ""):
        issues.append("fixed-time localization scale missing")
    if "diverges exponentially" not in published.get("uniformity_warning", ""):
        issues.append("fixed-time nonuniformity warning missing")

    checks = stored.get("exact", {}).get("checks", {})
    if len(checks.get("lattice_checks", [])) != 20:
        issues.append("stored arithmetic lattice check count drifted")
    if checks.get("odd_square_sum") != "pi**2/8":
        issues.append("stored odd-square sum drifted")

    note = NOTE.read_text(encoding="utf-8")
    required = [
        "not a proof of RH or `Lambda <= 0`",
        "B_rho(x)=-pi/8+O_a(1/x)",
        "reference field is therefore -pi/8",
        "scale `exp(C/Lambda)`",
        "threshold diverges as `Lambda->0+`",
        "move every zero by less than one",
        "macroscopic counting control, cannot determine the collision field",
        "fixed-lambda asymptotics do not prove `Lambda<=0`",
        "https://arxiv.org/abs/1904.12438",
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
        "validated Jensen-window PF Newman classical-field balance gate: "
        "10 rows, 0 issues, 3 exact field identities, 1 arithmetic equilibrium, "
        "1 continuum -pi/8 benchmark, 1 quantile-drift match, "
        "1 published fixed-time localization theorem, 2 exact sensitivity countermodels, "
        "1 compactness reduction, 1 open lambda-uniform handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
