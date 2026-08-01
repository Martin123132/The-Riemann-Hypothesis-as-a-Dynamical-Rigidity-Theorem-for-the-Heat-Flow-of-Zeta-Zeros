#!/usr/bin/env python3
"""Validate the weighted fractional-part autocorrelation Gram bridge."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_weighted_fractional_autocorrelation_gram_bridge.md"
)

REQUIRED_IDS = {
    f"wfagb_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "parameters"),
        (2, "weighted_autocorrelation"),
        (3, "scaling_law"),
        (4, "log_stationary_kernel"),
        (5, "fractional_part_mellin"),
        (6, "autocorrelation_mellin"),
        (7, "fourier_density"),
        (8, "positive_definiteness"),
        (9, "burnol_reciprocal_norm"),
        (10, "cross_kernel"),
        (11, "stationary_gram_form"),
        (12, "spectral_gram_form"),
        (13, "mellin_plancherel_match"),
        (14, "diagonal_split"),
        (15, "uniform_diagonal_bound"),
        (16, "signed_off_diagonal"),
        (17, "absolute_value_barrier"),
        (18, "source_scope"),
        (19, "nonpromotion_guard"),
        (20, "open_off_diagonal_gate"),
    )
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Weighted Fractional Autocorrelation Gram Bridge",
    "A_alpha(lambda)",
    "C_alpha(-u)=C_alpha(u)",
    "alpha-1<Re(s)<0",
    "|zeta(a+it)|^2 / (a^2+t^2)",
    "Phi_alpha>=0",
    "I_alpha(d,e)",
    "Gamma_(alpha,N)",
    "beta=1/2+omega",
    "Diag_(alpha,N)",
    "Off_(alpha,N)",
    "O(N^(1-alpha))",
    "Weighted signed off-diagonal gate",
    "https://arxiv.org/abs/math/0306251",
    "weighted signed off-diagonal gate is open",
)


def fractional_breakpoints(
    numerator: int,
    denominator: int,
    x_max: Fraction,
) -> set[Fraction]:
    rate = Fraction(numerator, denominator)
    points = {Fraction(0), x_max}
    for integer in range(1, int(x_max) + 1):
        points.add(Fraction(integer))
    max_second = int(rate * x_max)
    for integer in range(1, max_second + 1):
        point = Fraction(integer, 1) / rate
        if point <= x_max:
            points.add(point)
    return points


def mp_fraction(value: Fraction) -> mp.mpf:
    return mp.mpf(value.numerator) / value.denominator


def monomial_integral(
    coefficient: mp.mpf,
    exponent: mp.mpf,
    left: Fraction,
    right: Fraction,
) -> mp.mpf:
    if coefficient == 0:
        return mp.mpf("0")
    left_mp = mp_fraction(left)
    right_mp = mp_fraction(right)
    return coefficient * (
        mp.power(right_mp, exponent + 1)
        - mp.power(left_mp, exponent + 1)
    ) / (exponent + 1)


def truncated_autocorrelation(
    alpha: mp.mpf,
    numerator: int,
    denominator: int,
    x_max: Fraction,
) -> mp.mpf:
    """Integrate A_alpha(p/q) to x_max by exact cell antiderivatives."""
    rate = Fraction(numerator, denominator)
    ordered = sorted(
        fractional_breakpoints(numerator, denominator, x_max)
    )
    rate_mp = mp_fraction(rate)
    total = mp.mpf("0")
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        first_floor = midpoint.numerator // midpoint.denominator
        second_value = rate * midpoint
        second_floor = (
            second_value.numerator // second_value.denominator
        )
        linear_coefficient = -(
            mp.mpf(second_floor) + rate_mp * first_floor
        )
        constant_coefficient = mp.mpf(first_floor * second_floor)
        total += monomial_integral(
            rate_mp,
            alpha,
            left,
            right,
        )
        total += monomial_integral(
            linear_coefficient,
            alpha - 1,
            left,
            right,
        )
        total += monomial_integral(
            constant_coefficient,
            alpha - 2,
            left,
            right,
        )
    return total


def direct_two_term_burnol_integral(
    alpha: mp.mpf,
    x_max: Fraction,
) -> mp.mpf:
    coefficient = mp.power(2, -alpha)
    points = fractional_breakpoints(1, 1, x_max)
    points.update(fractional_breakpoints(1, 2, x_max))
    ordered = sorted(points)
    total = mp.mpf("0")
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        floor_one = midpoint.numerator // midpoint.denominator
        half_midpoint = midpoint / 2
        floor_half = (
            half_midpoint.numerator // half_midpoint.denominator
        )
        slope = 1 - coefficient / 2
        constant = -floor_one + coefficient * floor_half
        total += monomial_integral(
            slope**2,
            alpha,
            left,
            right,
        )
        total += monomial_integral(
            2 * slope * constant,
            alpha - 1,
            left,
            right,
        )
        total += monomial_integral(
            constant**2,
            alpha - 2,
            left,
            right,
        )
    return total


def independent_parameter_check() -> list[str]:
    issues: list[str] = []
    for alpha in (mp.mpf("0.2"), mp.mpf("0.5"), mp.mpf("0.8")):
        a_value = (1 - alpha) / 2
        beta = (1 + alpha) / 2
        if abs(a_value + beta - 1) > mp.mpf("1e-60"):
            issues.append(f"a+beta normalization failed at alpha={alpha}")
        if not (0 < a_value < mp.mpf("0.5") < beta < 1):
            issues.append(f"parameter ranges failed at alpha={alpha}")
    return issues


def independent_mellin_density_check() -> list[str]:
    issues: list[str] = []
    alpha = mp.mpf("0.5")
    a_value = (1 - alpha) / 2
    for t_value in (mp.mpf("0"), mp.mpf("0.7"), mp.mpf("3.2")):
        s_value = -a_value + 1j * t_value
        mellin = (
            -mp.zeta(-s_value)
            * mp.zeta(1 - alpha + s_value)
            / (s_value * (s_value + 1 - alpha))
        )
        density = (
            abs(mp.zeta(a_value + 1j * t_value)) ** 2
            / (a_value**2 + t_value**2)
        )
        if abs(mellin - density) > mp.mpf("1e-55"):
            issues.append(f"Mellin-density identity failed at t={t_value}")
        if density < 0:
            issues.append(f"negative spectral density at t={t_value}")
    return issues


def independent_kernel_scaling_check() -> list[str]:
    issues: list[str] = []
    alpha = mp.mpf("0.5")
    a_value = (1 - alpha) / 2
    a_two_thirds = truncated_autocorrelation(
        alpha, 2, 3, Fraction(12)
    )
    a_three_halves = truncated_autocorrelation(
        alpha, 3, 2, Fraction(8)
    )

    scaling_error = abs(
        a_two_thirds
        - mp.power(mp.mpf(2) / 3, 1 - alpha) * a_three_halves
    )
    if scaling_error > mp.mpf("1e-50"):
        issues.append(
            "matched finite autocorrelation scaling failed: "
            f"error={mp.nstr(scaling_error, 8)}"
        )

    cross_23 = mp.power(2, alpha - 1) * a_two_thirds
    cross_32 = mp.power(3, alpha - 1) * a_three_halves
    if abs(cross_23 - cross_32) > mp.mpf("1e-50"):
        issues.append("cross-kernel symmetry failed for d=2,e=3")

    a_two = truncated_autocorrelation(
        alpha, 2, 1, Fraction(10)
    )
    a_half = truncated_autocorrelation(
        alpha, 1, 2, Fraction(20)
    )
    c_log_two = mp.power(2, -a_value) * a_two
    c_minus_log_two = mp.power(mp.mpf("0.5"), -a_value) * a_half
    if abs(c_log_two - c_minus_log_two) > mp.mpf("1e-50"):
        issues.append("matched finite log-kernel evenness failed")
    return issues


def independent_gram_and_diagonal_check() -> list[str]:
    issues: list[str] = []
    alpha = mp.mpf("0.5")
    x_max = Fraction(24)
    j_11 = truncated_autocorrelation(alpha, 1, 1, x_max)
    j_22 = (
        mp.power(2, alpha - 1)
        * truncated_autocorrelation(alpha, 1, 1, x_max / 2)
    )
    j_12 = truncated_autocorrelation(alpha, 1, 2, x_max)
    coefficient_two = -mp.power(2, -alpha)
    pair_expansion = (
        j_11
        + coefficient_two**2 * j_22
        + 2 * coefficient_two * j_12
    )
    direct_integral = direct_two_term_burnol_integral(alpha, x_max)
    if abs(pair_expansion - direct_integral) > mp.mpf("1e-50"):
        issues.append("direct/pair-expanded finite Burnol norm failed")
    if direct_integral < -mp.mpf("1e-50"):
        issues.append("finite reciprocal norm is negative")

    n_max = 50
    diagonal_sum = mp.fsum(
        int(sp.mobius(d)) ** 2 * mp.power(d, -(1 + alpha))
        for d in range(1, n_max + 1)
    )
    if not (0 <= diagonal_sum <= mp.zeta(1 + alpha)):
        issues.append("uniform diagonal bound failed")
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
        "jensen_window_pf_weighted_fractional_"
        "autocorrelation_gram_bridge"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    ids = {item.get("id") for item in rows}
    if len(rows) != 20:
        issues.append(f"expected 20 rows, found {len(rows)}")
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
        "wfagb_13_mellin_plancherel_match": "available_exact",
        "wfagb_18_source_scope": "source_backed_guard",
        "wfagb_19_nonpromotion_guard": "guarded",
        "wfagb_20_open_off_diagonal_gate": "open_target",
    }
    by_id = {item.get("id"): item for item in rows}
    for row_id, expected in expected_statuses.items():
        if by_id.get(row_id, {}).get("status") != expected:
            issues.append(f"{row_id}: expected status {expected}")

    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 20,
        "exact_identity_count": 12,
        "normalization_crosscheck_count": 1,
        "proof_guard_count": 3,
        "literature_guard_count": 1,
        "open_signed_off_diagonal_gate_count": 1,
        "uniform_gram_bound_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key}: expected {expected}, found {audit.get(key)}"
            )

    for required in REQUIRED_NOTE:
        if required not in note:
            issues.append(f"note missing required text: {required!r}")

    forbidden_promotions = (
        "therefore RH is proved",
        "proves PF-infinity",
        "establishes Lambda <= 0",
        "the off-diagonal gate is closed",
    )
    lowered = note.lower()
    for forbidden in forbidden_promotions:
        if forbidden.lower() in lowered:
            issues.append(f"forbidden promotion: {forbidden!r}")

    mp.mp.dps = 60
    issues.extend(independent_parameter_check())
    issues.extend(independent_mellin_density_check())
    issues.extend(independent_kernel_scaling_check())
    issues.extend(independent_gram_and_diagonal_check())

    print(
        "validated weighted fractional autocorrelation Gram bridge: "
        f"{len(rows)} rows, {len(issues)} issues, "
        f"{audit.get('exact_identity_count')} exact identities, "
        f"{audit.get('proof_guard_count')} proof guards, "
        f"{audit.get('open_signed_off_diagonal_gate_count')} "
        "open signed off-diagonal gate"
    )
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
