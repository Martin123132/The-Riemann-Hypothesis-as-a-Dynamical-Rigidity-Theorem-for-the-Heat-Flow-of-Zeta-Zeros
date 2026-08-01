#!/usr/bin/env python3
"""Validate the Suzuki cofinal L2-hierarchy theorem audit."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_suzuki_cofinal_l2_hierarchy.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_suzuki_cofinal_l2_hierarchy.md"
)

REQUIRED_IDS = {
    f"sl2_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "shifted_quotient"),
        (2, "mellin_laplace_identity"),
        (3, "central_taylor_polynomial"),
        (4, "log_polynomial"),
        (5, "residual_transform"),
        (6, "l2_to_hardy"),
        (7, "hardy_to_pole_free"),
        (8, "pole_free_to_inner"),
        (9, "inner_remainder_hardy"),
        (10, "hardy_to_l2_residual"),
        (11, "fixed_shift_equivalence"),
        (12, "cofinal_reverse"),
        (13, "rh_forward"),
        (14, "cofinal_equivalence"),
        (15, "jordan_error_form"),
        (16, "level_one_target"),
        (17, "higher_level_target"),
        (18, "boundary_trace_guard"),
        (19, "fixed_shift_guard"),
        (20, "finite_energy_guard"),
        (21, "riesz_comparison_guard"),
        (22, "open_energy_target"),
    )
}

REQUIRED_NOTE_STRINGS = (
    "# Jensen-Window PF Suzuki Cofinal L2 Hierarchy",
    "This is not a proof of an `L2` residual, RH",
    "Q_omega(z)=xi(1/2+z-omega)/xi(1/2+z+omega)",
    "P_(omega,k)(t)",
    "[Q_omega(z)-T_(omega,k-1)(z)]/z^k",
    "r_(omega,k) in L2(0,infinity)",
    "iff D(omega)",
    "there exist omega_j->0 and integers k_j>=1",
    "integral_1^infinity |H_(omega,1)(x)-1|^2 dx/x",
    "q_(omega,1)",
    "-2*xi'(1/2+omega)/xi(1/2+omega)",
    "normalized cumulative",
    "(SL2.7)",
    "Q_good(z)=(a-z)/(a+z)",
    "Q_bad(z)=(a+z)/(a-z)",
    "boundary unimodularity gives",
    "regularized boundary trace",
    "positive-time support, not finite boundary",
    "anti-causal decaying function",
    "Riesz-Smoothing Comparison",
    "https://arxiv.org/abs/1204.1823",
    "https://arxiv.org/abs/1204.1827",
    "https://arxiv.org/abs/2207.07722",
    "does not prove RH or `Lambda <= 0`",
)


def evaluate_polynomial(coefficients: list[Fraction], z: Fraction) -> Fraction:
    value = Fraction(0)
    for coefficient in reversed(coefficients):
        value = value * z + coefficient
    return value


def independent_rational_inner_check() -> list[str]:
    issues: list[str] = []
    a = Fraction(3, 2)
    z = Fraction(1, 7)

    q_good = lambda j: (
        Fraction(1) if j == 0 else Fraction(2) * ((-1) ** j) / (a**j)
    )
    q_bad = lambda j: (
        Fraction(1) if j == 0 else Fraction(2) / (a**j)
    )
    good_value = (a - z) / (a + z)
    bad_value = (a + z) / (a - z)

    for k in range(1, 7):
        good_coefficients = [q_good(j) for j in range(k)]
        bad_coefficients = [q_bad(j) for j in range(k)]
        good_remainder = (
            good_value - evaluate_polynomial(good_coefficients, z)
        ) / (z**k)
        bad_remainder = (
            bad_value - evaluate_polynomial(bad_coefficients, z)
        ) / (z**k)
        expected_good = (
            Fraction(2) * ((-1) ** k) / (a ** (k - 1) * (a + z))
        )
        expected_bad = Fraction(2) / (a ** (k - 1) * (a - z))
        if good_remainder != expected_good:
            issues.append(f"k={k}: exact good remainder identity failed")
        if bad_remainder != expected_bad:
            issues.append(f"k={k}: exact bad remainder identity failed")

        # The inverse Laplace image of q_j/z^(k-j) is the matching
        # coefficient of t^(k-j-1)/(k-j-1)!.
        laplace_polynomial = sum(
            good_coefficients[j] / (z ** (k - j))
            for j in range(k)
        )
        expected_polynomial = (
            evaluate_polynomial(good_coefficients, z) / (z**k)
        )
        if laplace_polynomial != expected_polynomial:
            issues.append(f"k={k}: polynomial Laplace identity failed")
    return issues


def independent_energy_and_boundary_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 60
    a = mp.mpf("1.5")

    for t in (mp.mpf("0.25"), mp.mpf("1"), mp.mpf("7")):
        z = 1j * t
        good = (a - z) / (a + z)
        bad = (a + z) / (a - z)
        if abs(abs(good) - 1) > mp.mpf("1e-50"):
            issues.append(f"good boundary modulus failed at t={t}")
        if abs(abs(bad) - 1) > mp.mpf("1e-50"):
            issues.append(f"bad boundary modulus failed at t={t}")

    for k in range(1, 6):
        amplitude = 2 * ((-1) ** k) / (a ** (k - 1))
        numerical_energy = mp.quad(
            lambda t: abs(amplitude * mp.exp(-a * t)) ** 2,
            [0, mp.inf],
        )
        exact_energy = 2 / (a ** (2 * k - 1))
        if abs(numerical_energy - exact_energy) > mp.mpf("1e-45"):
            issues.append(f"k={k}: decaying residual energy failed")

        bad_energy_1 = (
            abs(amplitude) ** 2 * (mp.exp(2 * a) - 1) / (2 * a)
        )
        bad_energy_2 = (
            abs(amplitude) ** 2 * (mp.exp(4 * a) - 1) / (2 * a)
        )
        if not bad_energy_2 > 10 * bad_energy_1:
            issues.append(f"k={k}: growing residual guard failed")

        bad_amplitude = 2 / (a ** (k - 1))
        boundary_energy = mp.quad(
            lambda y: abs(bad_amplitude / (a - 1j * y)) ** 2,
            [-mp.inf, mp.inf],
        ) / (2 * mp.pi)
        if abs(boundary_energy - exact_energy) > mp.mpf("1e-45"):
            issues.append(f"k={k}: finite bad boundary energy failed")

        for y in (mp.mpf("0.5"), mp.mpf("3")):
            anticausal_transform = mp.quad(
                lambda t: bad_amplitude
                * mp.exp(a * t)
                * mp.exp(-1j * y * t),
                [-mp.inf, 0],
            )
            expected_transform = bad_amplitude / (a - 1j * y)
            if abs(anticausal_transform - expected_transform) > mp.mpf(
                "1e-45"
            ):
                issues.append(
                    f"k={k}, y={y}: anti-causal transform failed"
                )
    return issues


def independent_central_derivative_check() -> list[str]:
    issues: list[str] = []
    mp.mp.dps = 60
    omega = mp.mpf(1) / 4

    def xi(value):
        return (
            value
            * (value - 1)
            * mp.power(mp.pi, -value / 2)
            * mp.gamma(value / 2)
            * mp.zeta(value)
        )

    def quotient(z):
        return xi(mp.mpf("0.5") + z - omega) / xi(
            mp.mpf("0.5") + z + omega
        )

    direct_q1 = mp.diff(quotient, 0)
    point = mp.mpf("0.5") + omega
    expected_q1 = -2 * mp.diff(lambda s: mp.log(xi(s)), point)
    if abs(direct_q1 - expected_q1) > mp.mpf("1e-50"):
        issues.append(
            "central derivative identity failed: "
            f"error={abs(direct_q1 - expected_q1)}"
        )
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

    if payload.get("kind") != "jensen_window_pf_suzuki_cofinal_l2_hierarchy":
        issues.append("bad kind")

    rows = payload.get("rows", [])
    by_id = {item.get("id"): item for item in rows}
    if len(rows) != 22:
        issues.append(f"expected 22 rows, found {len(rows)}")
    missing = REQUIRED_IDS - set(by_id)
    if missing:
        issues.append(f"missing row ids: {sorted(missing)}")
    if len(by_id) != len(rows):
        issues.append("duplicate row id")

    for item in rows:
        if not item.get("role") or not item.get("status"):
            issues.append(f"{item.get('id')}: missing role or status")
        if not str(item.get("proof_boundary", "")).strip():
            issues.append(f"{item.get('id')}: missing proof boundary")

    for row_id in (
        "sl2_18_boundary_trace_guard",
        "sl2_19_fixed_shift_guard",
        "sl2_20_finite_energy_guard",
    ):
        item = by_id.get(row_id, {})
        if item.get("role") != "countermodel_gate":
            issues.append(f"{row_id}: bad guard role")
        if item.get("status") != "guard_validated":
            issues.append(f"{row_id}: bad guard status")

    literature_guard = by_id.get("sl2_21_riesz_comparison_guard", {})
    if literature_guard.get("role") != "literature_fit_guard":
        issues.append("Riesz comparison has wrong role")
    if literature_guard.get("status") != "guard_validated":
        issues.append("Riesz comparison has wrong status")

    target = by_id.get("sl2_22_open_energy_target", {})
    if target.get("role") != "open_arithmetic_gate":
        issues.append("open target has wrong role")
    if target.get("status") != "open_target":
        issues.append("open target has wrong status")

    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 22,
        "countermodel_gate_count": 3,
        "literature_fit_guard_count": 1,
        "open_arithmetic_gate_count": 1,
        "fixed_shift_equivalence_candidate_count": 1,
        "l2_residual_proved": False,
        "rh_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key}: expected {expected!r}, "
                f"found {audit.get(key)!r}"
            )

    source_urls = {item.get("url") for item in payload.get("sources", [])}
    for required in (
        "https://arxiv.org/abs/1204.1823",
        "https://arxiv.org/abs/1204.1827",
        "https://arxiv.org/abs/2207.07722",
        "outputs/jensen_window_pf_suzuki_fixed_omega_phase_diagram.md",
        "outputs/jensen_window_pf_suzuki_jordan_error_kernel_reduction.md",
    ):
        if required not in source_urls:
            issues.append(f"missing source: {required}")

    for needle in REQUIRED_NOTE_STRINGS:
        if needle not in note:
            issues.append(f"note missing: {needle}")

    issues.extend(independent_rational_inner_check())
    issues.extend(independent_energy_and_boundary_check())
    issues.extend(independent_central_derivative_check())

    lowered = note.lower()
    for forbidden in (
        "we have proved rh",
        "this proves rh",
        "therefore rh is true",
        "finite energy proves",
        "boundary modulus proves innerness",
    ):
        if forbidden in lowered:
            issues.append(f"forbidden promotion language: {forbidden}")

    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        print(f"failed Suzuki cofinal L2 hierarchy: {len(issues)} issues")
        return 1

    print(
        "validated Suzuki cofinal L2 hierarchy: "
        "22 rows, 0 issues, 3 countermodel guards, "
        "1 literature-fit guard, 1 cofinal equivalence candidate, "
        "1 open arithmetic energy gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
