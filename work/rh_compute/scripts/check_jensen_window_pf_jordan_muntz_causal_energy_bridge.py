#!/usr/bin/env python3
"""Validate the Jordan-Muntz causal-energy and mollifier bridge."""

from __future__ import annotations

import argparse
from itertools import combinations
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_jordan_muntz_causal_energy_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_jordan_muntz_causal_energy_bridge.md"
)

REQUIRED_IDS = {
    f"jmce_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "dirichlet_series"),
        (2, "jordan_error"),
        (3, "power_remainder"),
        (4, "full_muntz_identity"),
        (5, "tail_absolute"),
        (6, "remainder_hilbert_space"),
        (7, "remainder_mellin"),
        (8, "unitary_dilations"),
        (9, "partial_mobius_sum"),
        (10, "partial_mellin"),
        (11, "partial_plancherel"),
        (12, "mollified_right_line"),
        (13, "pointwise_limit"),
        (14, "fatou_handoff"),
        (15, "positive_time_energy"),
        (16, "positive_laplace_transform"),
        (17, "boundary_energy_automatic"),
        (18, "l2_implies_fixed_shift"),
        (19, "fixed_shift_implies_l2"),
        (20, "fixed_shift_equivalence"),
        (21, "cofinal_error_criterion"),
        (22, "rh_natural_convergence"),
        (23, "cofinal_partial_norm"),
        (24, "finite_energy_formula"),
        (25, "signed_measure_energy"),
        (26, "max_kernel_tn"),
        (27, "tn_nonclosure_guard"),
        (28, "causality_countermodel"),
        (29, "l2_coefficients_guard"),
        (30, "burnol_fit"),
        (31, "nyman_nontransfer"),
        (32, "open_mollifier_gate"),
    )
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Jordan-Muntz Causal-Energy Bridge",
    "E_omega(x)",
    "sum_(d>=1)mu(d)d^(-omega)R_omega(x/d)",
    "R_omega in H",
    "M_N(u)=sum_(d<=N)mu(d)d^(-u)",
    "Mellin-Plancherel gives the exact finite identity",
    "W_omega(t)",
    "sup_N ||E_(omega,N)||_H < infinity",
    "G_omega(s)",
    "J_omega<infinity",
    "iff D(omega)",
    "Finite Brownian Energy",
    "1/max(1,u,v)-1/X",
    "cumulative-matrix factorization",
    "same boundary squared norm `1/(2a)`",
    "compact support",
    "all-height finite inequality",
    "does not",
    "RH, or `Lambda <= 0`",
)


def mobius(n: int) -> int:
    return int(sp.mobius(n))


def coefficient(omega: mp.mpf, n: int) -> mp.mpf:
    return mp.fsum(
        mobius(int(d))
        * mp.power(int(d), -omega)
        * mp.power(n // int(d), omega)
        for d in sp.divisors(n)
    )


def independent_muntz_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    a_omega = 1 / ((1 + omega) * mp.zeta(1 + 2 * omega))
    for x in (mp.mpf(7), mp.mpf("12.5"), mp.mpf(31)):
        upper = int(mp.floor(x))
        cumulative = mp.fsum(
            coefficient(omega, n) for n in range(1, upper + 1)
        )
        direct = cumulative - a_omega * mp.power(x, 1 + omega)

        local = mp.mpf("0")
        reciprocal_partial = mp.mpf("0")
        for d in range(1, upper + 1):
            mu_d = mobius(d)
            reciprocal_partial += mu_d * mp.power(
                d, -(1 + 2 * omega)
            )
            y = x / d
            power_sum = mp.fsum(
                mp.power(m, omega)
                for m in range(1, int(mp.floor(y)) + 1)
            )
            remainder = (
                power_sum
                - mp.power(y, 1 + omega) / (1 + omega)
            )
            local += mu_d * mp.power(d, -omega) * remainder

        reciprocal_tail = (
            1 / mp.zeta(1 + 2 * omega) - reciprocal_partial
        )
        full_muntz = (
            local
            - mp.power(x, 1 + omega)
            / (1 + omega)
            * reciprocal_tail
        )
        if abs(direct - full_muntz) > mp.mpf("1e-55"):
            issues.append(
                f"full Muntz identity failed at x={x}: "
                f"{abs(direct-full_muntz)}"
            )
    return issues


def independent_finite_energy_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    x_max = mp.mpf("10.5")
    upper = int(mp.floor(x_max))
    coefficients = [mp.mpf("0")] + [
        coefficient(omega, n) for n in range(1, upper + 1)
    ]
    cumulative = []
    running = mp.mpf("0")
    for n in range(1, upper + 1):
        running += coefficients[n]
        cumulative.append(running)
    a_omega = 1 / ((1 + omega) * mp.zeta(1 + 2 * omega))

    direct = mp.mpf("0")
    for n in range(1, upper + 1):
        right = min(mp.mpf(n + 1), x_max)
        if right <= n:
            continue
        c_value = cumulative[n - 1]
        direct += mp.quad(
            lambda x: (
                c_value - a_omega * mp.power(x, 1 + omega)
            )
            ** 2
            / (x * x),
            [mp.mpf(n), right],
        )

    discrete = mp.fsum(
        coefficients[n]
        * coefficients[m]
        * (
            1 / mp.mpf(max(n, m))
            - 1 / x_max
        )
        for n in range(1, upper + 1)
        for m in range(1, upper + 1)
    )
    cross = (
        -2
        * a_omega
        / omega
        * mp.fsum(
            coefficients[n]
            * (
                mp.power(x_max, omega)
                - mp.power(n, omega)
            )
            for n in range(1, upper + 1)
        )
    )
    continuum = (
        a_omega
        * a_omega
        * (
            mp.power(x_max, 1 + 2 * omega) - 1
        )
        / (1 + 2 * omega)
    )
    closed = discrete + cross + continuum
    if abs(direct - closed) > mp.mpf("1e-50"):
        issues.append(
            "finite energy formula mismatch: "
            f"{abs(direct-closed)}"
        )
    if direct < 0:
        issues.append("finite energy became negative")
    return issues


def independent_tn_check() -> list[str]:
    issues: list[str] = []
    points = (1, 2, 4, 7)
    x_max = sp.Rational(10)
    matrix = sp.Matrix(
        [
            [
                sp.Rational(1, max(u, v)) - 1 / x_max
                for v in points
            ]
            for u in points
        ]
    )
    for order in range(1, len(points) + 1):
        for row_indices in combinations(range(len(points)), order):
            for column_indices in combinations(
                range(len(points)), order
            ):
                minor = sp.factor(
                    matrix.extract(
                        row_indices, column_indices
                    ).det()
                )
                if minor.is_nonnegative is not True:
                    issues.append(
                        "max-kernel TN check failed: "
                        f"rows={row_indices}, cols={column_indices}, "
                        f"minor={minor}"
                    )
    return issues


def chi_factor(value: mp.mpc) -> mp.mpc:
    return (
        mp.power(mp.pi, value - mp.mpf("0.5"))
        * mp.gamma((1 - value) / 2)
        / mp.gamma(value / 2)
    )


def independent_functional_equation_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4
    for height in (mp.mpf("3.25"), mp.mpf("11.75")):
        s = mp.mpf("0.5") + 1j * height
        mollifier = mp.fsum(
            mobius(n) * mp.power(n, -(s + omega))
            for n in range(1, 18)
        )
        left = abs(
            mp.zeta(s - omega) / s * mollifier
        ) ** 2
        right = (
            abs(chi_factor(s - omega)) ** 2
            / abs(s) ** 2
            * abs(mp.zeta(s + omega) * mollifier) ** 2
        )
        relative = abs(left - right) / max(
            mp.mpf("1e-60"), abs(left), abs(right)
        )
        if relative > mp.mpf("1e-55"):
            issues.append(
                "functional-equation mollifier rewrite failed at "
                f"t={height}: {relative}"
            )
    return issues


def completed_xi(value: mp.mpc) -> mp.mpc:
    return (
        mp.mpf("0.5")
        * value
        * (value - 1)
        * mp.power(mp.pi, -value / 2)
        * mp.gamma(value / 2)
        * mp.zeta(value)
    )


def independent_archimedean_factor_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 70
    omega = mp.mpf(1) / 4

    def archimedean_factor(value: mp.mpc) -> mp.mpc:
        return (
            (value - omega)
            * (value - omega - 1)
            / ((value + omega) * (value + omega - 1))
            * mp.power(mp.pi, omega)
            * mp.gamma((value - omega) / 2)
            / mp.gamma((value + omega) / 2)
        )

    for value in (
        mp.mpc("0.9", "2.75"),
        mp.mpc("1.35", "6.25"),
    ):
        quotient = completed_xi(
            value - omega
        ) / completed_xi(value + omega)
        factored = (
            archimedean_factor(value)
            * mp.zeta(value - omega)
            / mp.zeta(value + omega)
        )
        relative = abs(quotient - factored) / max(
            mp.mpf("1e-60"), abs(quotient), abs(factored)
        )
        if relative > mp.mpf("1e-55"):
            issues.append(
                "archimedean quotient factorization failed at "
                f"s={value}: {relative}"
            )

    s_zero = 1 + omega
    residue_l = (
        (1 + 2 * omega)
        * (2 * omega)
        / s_zero
        * mp.power(mp.pi, -omega)
        * mp.gamma(mp.mpf("0.5") + omega)
        / mp.gamma(mp.mpf("0.5"))
    )
    quotient_at_pole = (
        mp.mpf("0.5") / completed_xi(1 + 2 * omega)
    )
    a_omega = 1 / (
        (1 + omega) * mp.zeta(1 + 2 * omega)
    )
    if abs(residue_l * quotient_at_pole - a_omega) > mp.mpf(
        "1e-55"
    ):
        issues.append(
            "archimedean pole residue does not match A_omega"
        )
    return issues


def independent_countermodel_checks() -> list[str]:
    issues: list[str] = []
    a, z, y, t = sp.symbols(
        "a z y t", positive=True, real=True
    )
    boundary_energy = sp.integrate(
        1 / (a * a + y * y), (y, -sp.oo, sp.oo)
    ) / (2 * sp.pi)
    if sp.simplify(boundary_energy - 1 / (2 * a)) != 0:
        issues.append("causality boundary-energy identity failed")

    good_laplace = sp.integrate(
        sp.exp(-a * t) * sp.exp(-z * t),
        (t, 0, sp.oo),
    )
    bad_laplace = sp.integrate(
        -sp.exp(a * t) * sp.exp(-z * t),
        (t, 0, sp.oo),
        conds="none",
    )
    if sp.simplify(good_laplace - 1 / (a + z)) != 0:
        issues.append("good causal transform failed")
    if sp.simplify(bad_laplace - 1 / (a - z)) != 0:
        issues.append("bad positive-time transform failed")

    n = sp.symbols("n", integer=True, positive=True)
    square_sum = sp.summation(1 / n**2, (n, 1, sp.oo))
    if sp.simplify(square_sum - sp.pi**2 / 6) != 0:
        issues.append("l2 coefficient countermodel square sum failed")
    harmonic_limit = sp.limit(
        sp.harmonic(n), n, sp.oo
    )
    if harmonic_limit != sp.oo:
        issues.append("l2 coefficient vector-series divergence failed")
    return issues


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []

    if payload.get("kind") != (
        "jensen_window_pf_jordan_muntz_causal_energy_bridge"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    ids = {item.get("id") for item in rows}
    if len(rows) != 32:
        issues.append(f"expected 32 rows, found {len(rows)}")
    if ids != REQUIRED_IDS:
        issues.append(
            f"row ids differ: missing={sorted(REQUIRED_IDS-ids)}, "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )
    if len(ids) != len(rows):
        issues.append("duplicate row id")
    for item in rows:
        if not item.get("proof_boundary"):
            issues.append(f"{item.get('id')}: missing proof boundary")

    expected_statuses = {
        "jmce_17_boundary_energy_automatic": "guard_validated",
        "jmce_18_l2_implies_fixed_shift": "theorem_candidate",
        "jmce_19_fixed_shift_implies_l2": "theorem_candidate",
        "jmce_20_fixed_shift_equivalence": "theorem_candidate",
        "jmce_21_cofinal_error_criterion": "theorem_candidate",
        "jmce_22_rh_natural_convergence": "conditional_on_rh",
        "jmce_23_cofinal_partial_norm": "theorem_candidate",
        "jmce_27_tn_nonclosure_guard": "guard_validated",
        "jmce_28_causality_countermodel": "guard_validated",
        "jmce_29_l2_coefficients_guard": "guard_validated",
        "jmce_31_nyman_nontransfer": "guard_validated",
        "jmce_32_open_mollifier_gate": "open_target",
    }
    by_id = {item.get("id"): item for item in rows}
    for row_id, status in expected_statuses.items():
        if by_id.get(row_id, {}).get("status") != status:
            issues.append(f"{row_id}: expected status {status}")

    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 32,
        "theorem_candidate_count": 5,
        "conditional_estimate_count": 1,
        "guard_count": 5,
        "open_arithmetic_gate_count": 1,
        "finite_energy_formula_checked": True,
        "muntz_identity_checked": True,
        "max_kernel_tn_checked": True,
        "uniform_mollifier_bound_proved": False,
        "rh_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key}: expected {expected!r}, "
                f"found {audit.get(key)!r}"
            )

    for marker in REQUIRED_NOTE:
        if marker not in note:
            issues.append(f"note missing: {marker}")

    issues.extend(independent_muntz_check())
    issues.extend(independent_finite_energy_check())
    issues.extend(independent_tn_check())
    issues.extend(independent_functional_equation_check())
    issues.extend(independent_archimedean_factor_check())
    issues.extend(independent_countermodel_checks())

    lowered = note.lower()
    for forbidden in (
        "this proves rh",
        "we have proved rh",
        "the uniform estimate is proved",
        "total nonnegativity proves the bound",
        "square-summable coefficients imply convergence",
        "boundary energy proves causality",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    if not math.isclose(
        float(audit.get("row_count", 0)), 32.0
    ):
        issues.append("non-numeric row-count audit")

    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        print(
            "failed Jordan-Muntz causal-energy bridge: "
            f"{len(issues)} issues"
        )
        return 1

    print(
        "validated Jordan-Muntz causal-energy bridge: "
        "32 rows, 0 issues, 5 theorem candidates, "
        "5 proof guards, 1 open all-height mollifier gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
