#!/usr/bin/env python3
"""Validate the weighted-autocorrelation OU tail-energy reduction."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_weighted_autocorrelation_ou_tail_energy_reduction.md"
)

REQUIRED_IDS = {
    f"waoter_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "parameters"),
        (2, "integrable_amplitude"),
        (3, "period_mean"),
        (4, "bernoulli_primitive"),
        (5, "quantitative_period_averaging"),
        (6, "large_ratio_asymptotic"),
        (7, "ou_kernel_split"),
        (8, "fourier_density_split"),
        (9, "remainder_sign_witness"),
        (10, "main_gram"),
        (11, "ou_factorization"),
        (12, "continuous_tail_energy"),
        (13, "discrete_tail_energy"),
        (14, "pnt_tail_limit"),
        (15, "fatou_weighted_l2"),
        (16, "reciprocal_zeta_series"),
        (17, "fixed_alpha_zero_free"),
        (18, "rh_implies_main_energy"),
        (19, "cofinal_rh_equivalence"),
        (20, "full_gram_noncomparison"),
        (21, "hilbert_large_sieve_barrier"),
        (22, "open_comparison_gate"),
    )
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Weighted Autocorrelation OU Tail-Energy Reduction",
    "L_alpha",
    "K_alpha=1/(2alpha)+zeta(2-alpha)/6",
    "Psi(y)=B_2({y})/2",
    "R_alpha(u)",
    "-0.0687554899394677",
    "not positive definite",
    "G_(alpha,N)",
    "sum_(d=k)^N mu(d)/d",
    "F_alpha(s)",
    "zeta(s)!=0 for Re(s)>(1+alpha)/2",
    "G_(1/(j+1),N)<infinity for every integer j>=1",
    "Equation (WAOTER.6)",
    "O(N^(1-alpha))",
    "lambda_n=log n",
    "skew kernel",
    "product kernel",
    "https://arxiv.org/abs/math/0306251",
    "https://arxiv.org/abs/2203.14950",
    "open target",
    "does not prove",
)


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


def breakpoints_for_rate(rate: int, x_max: int) -> list[Fraction]:
    points = {Fraction(0), Fraction(x_max)}
    points.update(Fraction(integer) for integer in range(1, x_max + 1))
    points.update(
        Fraction(integer, rate)
        for integer in range(1, rate * x_max + 1)
    )
    return sorted(points)


def truncated_autocorrelation(
    alpha: mp.mpf,
    rate: int,
    x_max: int,
) -> mp.mpf:
    total = mp.mpf("0")
    for left, right in zip(
        breakpoints_for_rate(rate, x_max),
        breakpoints_for_rate(rate, x_max)[1:],
    ):
        midpoint = (left + right) / 2
        first_floor = midpoint.numerator // midpoint.denominator
        second = rate * midpoint
        second_floor = second.numerator // second.denominator
        total += monomial_integral(
            mp.mpf(rate), alpha, left, right
        )
        total += monomial_integral(
            -(second_floor + rate * first_floor),
            alpha - 1,
            left,
            right,
        )
        total += monomial_integral(
            mp.mpf(first_floor * second_floor),
            alpha - 2,
            left,
            right,
        )
    return total


def truncated_single_amplitude(
    alpha: mp.mpf,
    x_max: int,
) -> mp.mpf:
    total = mp.mpf("0")
    for integer in range(x_max):
        left = Fraction(integer)
        right = Fraction(integer + 1)
        total += monomial_integral(
            mp.mpf("1"), alpha - 1, left, right
        )
        total += monomial_integral(
            -mp.mpf(integer), alpha - 2, left, right
        )
    return total


def independent_period_averaging_check() -> list[str]:
    issues: list[str] = []
    alpha = mp.mpf("0.5")
    k_alpha = 1 / (2 * alpha) + mp.zeta(2 - alpha) / 6
    x_max = 16
    mean_integral = truncated_single_amplitude(alpha, x_max) / 2
    for rate in (4, 8, 16):
        actual = truncated_autocorrelation(alpha, rate, x_max)
        error = abs(actual - mean_integral)
        bound = k_alpha * mp.power(rate, -alpha)
        if error > bound + mp.mpf("1e-50"):
            issues.append(
                f"period-averaging bound failed at lambda={rate}"
            )
    return issues


def independent_ou_energy_check() -> list[str]:
    issues: list[str] = []
    alpha = mp.mpf("0.5")
    a_value = (1 - alpha) / 2
    beta = (1 + alpha) / 2
    n_max = 12

    direct = mp.mpf("0")
    for d in range(1, n_max + 1):
        mu_d = int(sp.mobius(d))
        for e in range(1, n_max + 1):
            mu_e = int(sp.mobius(e))
            ratio = mp.mpf(min(d, e)) / max(d, e)
            direct += (
                mu_d
                * mu_e
                * mp.power(d * e, -beta)
                * mp.power(ratio, a_value)
            )

    discrete = mp.mpf("0")
    for k in range(1, n_max + 1):
        weight = (
            mp.power(k, 1 - alpha)
            - mp.power(k - 1, 1 - alpha)
        )
        tail = mp.fsum(
            int(sp.mobius(d)) / mp.mpf(d)
            for d in range(k, n_max + 1)
        )
        discrete += weight * tail**2

    if abs(direct - discrete) > mp.mpf("1e-55"):
        issues.append("OU Gram/discrete tail-energy identity failed")
    if direct < -mp.mpf("1e-55"):
        issues.append("OU Gram form is negative")
    return issues


def independent_remainder_sign_check() -> list[str]:
    issues: list[str] = []
    alpha = mp.mpf("0.5")
    a_value = (1 - alpha) / 2
    l_alpha = -mp.zeta(1 - alpha) / (2 * (1 - alpha))
    numerator = abs(mp.zeta(a_value)) ** 2 - 2 * a_value * l_alpha
    expected = mp.mpf(
        "-0.06875548993946772397902047370185327984625"
    )
    if abs(numerator - expected) > mp.mpf("1e-40"):
        issues.append("remainder sign-witness value changed")
    if numerator >= 0:
        issues.append("remainder sign witness is not negative")
    return issues


def independent_threshold_and_spacing_check() -> list[str]:
    issues: list[str] = []
    alpha = mp.mpf("0.5")
    beta = (1 + alpha) / 2
    if abs(2 * beta - alpha - 1) > mp.mpf("1e-60"):
        issues.append("Cauchy-Schwarz half-plane threshold failed")

    for n_value in (2, 5, 20, 100):
        delta = mp.log(1 + mp.mpf(1) / n_value)
        if not (
            mp.mpf(1) / (n_value + 1)
            <= delta
            <= mp.mpf(1) / n_value
        ):
            issues.append(f"log-spacing bounds failed at n={n_value}")

    n_max = 200
    spacing_cost = mp.fsum(
        int(sp.mobius(n)) ** 2
        * mp.power(n, -2 * beta)
        / mp.log(1 + mp.mpf(1) / n)
        for n in range(2, n_max + 1)
    )
    critical_cost = mp.fsum(
        int(sp.mobius(n)) ** 2 * mp.power(n, -alpha)
        for n in range(2, n_max + 1)
    )
    ratio = spacing_cost / critical_cost
    if not (mp.mpf("0.9") < ratio < mp.mpf("1.1")):
        issues.append("weighted Hilbert spacing cost scale failed")
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
        "jensen_window_pf_weighted_autocorrelation_"
        "ou_tail_energy_reduction"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    ids = {item.get("id") for item in rows}
    if len(rows) != 22:
        issues.append(f"expected 22 rows, found {len(rows)}")
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

    by_id = {item.get("id"): item for item in rows}
    expected_statuses = {
        "waoter_09_remainder_sign_witness": "guard_validated",
        "waoter_14_pnt_tail_limit": "source_backed",
        "waoter_19_cofinal_rh_equivalence": "available_exact",
        "waoter_20_full_gram_noncomparison": "guarded",
        "waoter_21_hilbert_large_sieve_barrier": "source_backed_guard",
        "waoter_22_open_comparison_gate": "open_target",
    }
    for row_id, expected in expected_statuses.items():
        if by_id.get(row_id, {}).get("status") != expected:
            issues.append(f"{row_id}: expected status {expected}")

    audit = payload.get("audit", {})
    expected_audit = {
        "row_count": 22,
        "exact_identity_count": 13,
        "classical_theorem_step_count": 2,
        "proof_guard_count": 4,
        "literature_guard_count": 2,
        "open_comparison_gate_count": 1,
        "cofinal_ou_energy_rh_equivalence_proved": True,
        "full_gram_comparison_proved": False,
        "full_burnol_bound_proved": False,
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
        "the full Burnol bound follows",
        "proves PF-infinity",
        "establishes Lambda <= 0",
        "the comparison gate is closed",
    )
    lowered = note.lower()
    for forbidden in forbidden_promotions:
        if forbidden.lower() in lowered:
            issues.append(f"forbidden promotion: {forbidden!r}")

    mp.mp.dps = 60
    issues.extend(independent_period_averaging_check())
    issues.extend(independent_ou_energy_check())
    issues.extend(independent_remainder_sign_check())
    issues.extend(independent_threshold_and_spacing_check())

    print(
        "validated weighted-autocorrelation OU tail-energy reduction: "
        f"{len(rows)} rows, {len(issues)} issues, "
        f"{audit.get('exact_identity_count')} exact identities, "
        f"{audit.get('proof_guard_count')} proof guards, "
        f"{audit.get('open_comparison_gate_count')} "
        "open comparison gate"
    )
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
