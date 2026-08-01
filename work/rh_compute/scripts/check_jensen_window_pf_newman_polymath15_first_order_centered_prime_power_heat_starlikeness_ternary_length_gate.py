#!/usr/bin/env python3
"""Validate the ternary all-length prime-power heat-starlikeness gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_ternary_length_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
LENGTH_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
BASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_base_certificate"
)
PHASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_logarithmic_phase_flow_gate"
)
SOURCES = {
    "length_propagation": (
        REPO_ROOT / "work/rh_compute/results" / f"{LENGTH_STEM}.json"
    ),
    "base_certificate": (
        REPO_ROOT / "work/rh_compute/results" / f"{BASE_STEM}.json"
    ),
    "phase_flow": (
        REPO_ROOT / "work/rh_compute/results" / f"{PHASE_STEM}.json"
    ),
}
Q0 = Fraction(17, 20)
ALPHA = Fraction(577351, 1_000_000)
WIDE_A = Fraction(29, 50)
EXPECTED_ROWS = [
    "pphstl_01_domain",
    "pphstl_02_fourier",
    "pphstl_03_common_groups",
    "pphstl_04_terminal_cluster",
    "pphstl_05_cluster_slope",
    "pphstl_06_small_d",
    "pphstl_07_large_d",
    "pphstl_08_unique_defect",
    "pphstl_09_zero_reserve",
    "pphstl_10_finite_energy",
    "pphstl_11_tail_implication",
    "pphstl_12_tail_decay",
    "pphstl_13_fejer_absorption",
    "pphstl_14_ideal_theorem",
    "pphstl_15_actual_sums",
    "pphstl_16_actual_current",
    "pphstl_17_pi_provenance",
    "pphstl_18_boundary",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        issues.append(f"cannot load {label}: {error}")
        return {}


def put(
    terms: dict[tuple[int, int], Fraction],
    q_power: int,
    a_power: int,
    coefficient: int | Fraction,
) -> None:
    key = (q_power, a_power)
    terms[key] = terms.get(key, Fraction()) + Fraction(coefficient)
    if terms[key] == 0:
        terms.pop(key)


def q_bernstein_vector(power: int, degree: int) -> list[Fraction]:
    denominator = comb(degree, power)
    result = []
    for index in range(degree + 1):
        first = max(0, power - degree + index)
        last = min(power, index)
        total = Fraction()
        for chosen_upper in range(first, last + 1):
            total += (
                Fraction(
                    comb(index, chosen_upper)
                    * comb(degree - index, power - chosen_upper),
                    denominator,
                )
                * Q0 ** (power - chosen_upper)
            )
        result.append(total)
    return result


def a_bernstein_vector(
    power: int, degree: int, upper: Fraction
) -> list[Fraction]:
    denominator = comb(degree, power)
    return [
        (
            upper**power * Fraction(comb(index, power), denominator)
            if index >= power
            else Fraction()
        )
        for index in range(degree + 1)
    ]


def rebuild_certificate(
    source: dict[tuple[int, int], Fraction],
    upper: Fraction,
    target: Fraction,
) -> dict:
    terms = dict(source)
    put(terms, 0, 0, -target)
    q_degree = max(key[0] for key in terms)
    a_degree = max(key[1] for key in terms)
    q_vectors = {
        power: q_bernstein_vector(power, q_degree)
        for power in {key[0] for key in terms}
    }
    a_vectors = {
        power: a_bernstein_vector(power, a_degree, upper)
        for power in {key[1] for key in terms}
    }
    minimum: Fraction | None = None
    argminimum = None
    for qi in range(q_degree + 1):
        for ai in range(a_degree + 1):
            coefficient = Fraction()
            for (q_power, a_power), scalar in terms.items():
                coefficient += (
                    scalar
                    * q_vectors[q_power][qi]
                    * a_vectors[a_power][ai]
                )
            if minimum is None or coefficient < minimum:
                minimum = coefficient
                argminimum = [qi, ai]
    if minimum is None or argminimum is None:
        raise RuntimeError("empty reconstructed certificate")
    return {
        "power_degrees": [q_degree, a_degree],
        "bernstein_shape": [q_degree + 1, a_degree + 1],
        "claimed_lower": f"{target.numerator}/{target.denominator}",
        "minimum_shifted_coefficient": (
            f"{minimum.numerator}/{minimum.denominator}"
        ),
        "argminimum": argminimum,
        "positive": minimum > 0,
    }


def cluster_polynomial(d: int) -> dict[tuple[int, int], Fraction]:
    B = d + 2
    result: dict[tuple[int, int], Fraction] = {}
    for q_power, a_power, scalar in (
        (0, 0, B),
        (3, 1, -2 * (B + 1)),
        (4, 2, B + 2),
        (2 * d + 6, 2, B + 2),
        (2 * d + 7, 3, -2 * (B + 3)),
        (4 * d + 8, 4, B + 4),
    ):
        put(result, q_power, a_power, scalar)
    return result


def four_group_polynomial() -> dict[tuple[int, int], Fraction]:
    result: dict[tuple[int, int], Fraction] = {}
    for k in range(4):
        A = 2 * k + 2
        base_q = 14 * k - 2 * k * k
        base_b = 2 * k
        pieces = [
            (A, 0, 0),
            (-2 * (A + 1), 1, 6 - 2 * k),
            (A + 2, 2, 10 - 4 * k),
        ]
        for scalar, b_increment, q_increment in pieces:
            b_power = base_b + b_increment
            put(
                result,
                base_q + q_increment + 3 * b_power,
                b_power,
                Fraction(scalar) * ALPHA**b_power,
            )
    return result


def first_difference_polynomial(
    M: int,
) -> dict[tuple[int, int], Fraction]:
    exponent = [k * (2 * M - 2 - k) for k in range(M)]
    result: dict[tuple[int, int], Fraction] = {}
    for k in range(M):
        put(result, 2 * exponent[k], 2 * k, 2 * (k + 1))
    for k in range(M - 1):
        put(
            result,
            exponent[k] + exponent[k + 1],
            2 * k + 1,
            -2 * (2 * k + 3),
        )
    for k in range(M - 2):
        put(
            result,
            exponent[k] + exponent[k + 2],
            2 * k + 2,
            2 * k + 4,
        )
    return result


def terminal_polynomial(M: int) -> dict[tuple[int, int], Fraction]:
    result: dict[tuple[int, int], Fraction] = {}
    q_power = M * (M - 2)
    a_power = M - 2
    put(result, q_power, a_power, M)
    put(result, q_power + 1, a_power + 1, -2 * (M + 1))
    put(result, q_power + 2 * M - 2, a_power + 2, M + 2)
    return result


def finite_energy_polynomial(
    M: int,
) -> dict[tuple[int, int], Fraction]:
    result = first_difference_polynomial(M)
    for key, coefficient in terminal_polynomial(M).items():
        put(result, key[0], key[1], (M - 1) ** 2 * coefficient)
    return result


def compare_certificate(
    stored: dict,
    rebuilt: dict,
    label: str,
    issues: list[str],
) -> None:
    for key in (
        "power_degrees",
        "bernstein_shape",
        "claimed_lower",
        "minimum_shifted_coefficient",
        "argminimum",
        "positive",
    ):
        if stored.get(key) != rebuilt.get(key):
            issues.append(f"{label} certificate drifted at {key}")


def validate_symbolic_algebra(issues: list[str]) -> None:
    a, q = sp.symbols("a q", positive=True)
    for M in (5, 6, 9):
        coefficients = [
            a**k * q ** (k * (2 * M - 2 - k)) for k in range(M)
        ]
        A = [sp.Integer(0)] * (M + 2)
        A[0] = 2 * sum(
            (k + 1) * coefficients[k] ** 2 for k in range(M)
        )
        for d in range(1, M):
            A[d] = sum(
                (2 * k + d + 2)
                * coefficients[k]
                * coefficients[k + d]
                for k in range(M - d)
            )
        terminal = sp.expand(A[M - 2] - 2 * A[M - 1])
        expected = sp.expand(
            a ** (M - 2)
            * q ** (M * (M - 2))
            * (
                M
                - 2 * (M + 1) * a * q
                + (M + 2) * a**2 * q ** (2 * M - 2)
            )
        )
        if sp.expand(terminal - expected) != 0:
            issues.append(f"terminal formula failed at M={M}")

    B, d = sp.symbols("B d", integer=True, positive=True)
    t = a**2 * q ** (2 * d + 6)
    interior = B - 2 * (B + 1) * a * q**3 + (B + 2) * a**2 * q**4
    boundary = (
        B
        + 2
        - 2 * (B + 3) * a * q
        + (B + 4) * a**2 * q ** (2 * d + 2)
    )
    slope = sp.diff(sp.expand(interior + t * boundary), B)
    expected_slope = (
        1
        - 2 * a * q**3
        + a**2 * q**4
        + t * (1 - 2 * a * q + a**2 * q ** (2 * d + 2))
    )
    if sp.expand(slope - expected_slope) != 0:
        issues.append("terminal-cluster slope formula failed")


def validate_exact_inequalities(issues: list[str]) -> None:
    if not 3 * ALPHA**2 > 1:
        issues.append("rational reciprocal-sqrt(3) upper failed")
    if not ALPHA < WIDE_A:
        issues.append("wide a-box does not contain alpha")
    if not 25 < 27:
        issues.append("1/sqrt(3)<3/5 comparison failed")
    if not 75 > 64:
        issues.append("cluster-slope radical comparison failed")
    if not 44**2 * 3 > 76**2:
        issues.append("d=6 analytic cluster endpoint failed")

    M = 17
    # R_M<3/5 is equivalent to
    # 25(M+1)^2<3(3M+1)^2 after squaring positive sides.
    if not 25 * (M + 1) ** 2 < 3 * (3 * M + 1) ** 2:
        issues.append("terminal R_17<3/5 failed")
    # h_M<M/3 follows from 3(M+1)^2<4M^2.
    if not 3 * (M + 1) ** 2 < 4 * M**2:
        issues.append("terminal h_17<M/3 failed")
    # sqrt(5)>2 gives the rational M=17 endpoint.
    if not 4352 * 100 < 468750:
        issues.append("terminal M=17 decay endpoint failed")
    if not M**2 - 5 * M + 2 > 0:
        issues.append("terminal majorant monotonicity failed")


def validate_certificates(stored: dict, issues: list[str]) -> None:
    exact = stored.get("exact", {})
    stored_clusters = exact.get("nonterminal", {}).get("certificates", [])
    if len(stored_clusters) != 5:
        issues.append("wrong stored small-d certificate count")
    else:
        for d, item in enumerate(stored_clusters, start=1):
            rebuilt = rebuild_certificate(
                cluster_polynomial(d), WIDE_A, Fraction(17, 100)
            )
            compare_certificate(
                item, rebuilt, f"nonterminal d={d}", issues
            )
            if item.get("d") != d:
                issues.append(f"nonterminal d={d} label drifted")

    zero = exact.get("zero_reserve", {})
    compare_certificate(
        zero.get("four_group_certificate", {}),
        rebuild_certificate(
            four_group_polynomial(), Fraction(1), Fraction(11, 250)
        ),
        "four-group zero reserve",
        issues,
    )
    compare_certificate(
        zero.get("m5_certificate", {}),
        rebuild_certificate(
            first_difference_polynomial(5), ALPHA, Fraction(1, 25)
        ),
        "M=5 zero reserve",
        issues,
    )

    stored_energy = exact.get("terminal_defect", {}).get(
        "finite_certificates", []
    )
    if len(stored_energy) != 12:
        issues.append("wrong stored finite-energy certificate count")
    else:
        for M, item in zip(range(5, 17), stored_energy):
            rebuilt = rebuild_certificate(
                finite_energy_polynomial(M), ALPHA, Fraction(1, 25)
            )
            compare_certificate(
                item, rebuilt, f"finite energy M={M}", issues
            )
            if item.get("M") != M:
                issues.append(f"finite energy M={M} label drifted")


def validate_actual_transfer(issues: list[str]) -> None:
    x_lower = 12 * 2**50
    relative = Fraction(18000, x_lower)
    coefficient_sum = Fraction(5, 2)
    moment_sum = Fraction(15, 4)
    current_mass = coefficient_sum * (
        coefficient_sum + moment_sum
    )
    if not current_mass == Fraction(125, 8) < 16:
        issues.append("ternary current-mass bound failed")
    if not 16 * (2 * relative + relative**2) < 48 * relative:
        issues.append("ternary frozen-current bound failed")
    if not Fraction(8446) * 5 == 42230:
        issues.append("ternary epsilon mass failed")
    ray_error = (
        Fraction(864020, x_lower)
        + Fraction(633450, x_lower**2)
    )
    if not ray_error < Fraction(1, 1000):
        issues.append("ternary actual-ray transfer failed")
    if not Fraction(17, 1000) - Fraction(1, 1000) == Fraction(
        2, 125
    ):
        issues.append("ternary actual-ray margin failed")


def validate_sources(stored: dict, issues: list[str]) -> None:
    expected = {
        key: file_hash(path)
        for key, path in SOURCES.items()
        if path.is_file()
    }
    if stored.get("source_audit", {}).get("source_sha256") != expected:
        issues.append("source hashes drifted")
    length = load_json(
        SOURCES["length_propagation"], "length source", issues
    )
    if length.get("summary", {}).get(
        "analytic_all_length_large_prime_families"
    ) != 1:
        issues.append("length source summary drifted")
    base = load_json(SOURCES["base_certificate"], "base source", issues)
    if base.get("exact", {}).get("actual_transfer", {}).get(
        "coefficient_error"
    ) != "|R_k-1|<r_x=18000/x":
        issues.append("base actual-transfer source drifted")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "result", issues)
    if not stored:
        return issues
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")

    summary = stored.get("summary", {})
    expected_summary = {
        "rows": 18,
        "exact_small_d_cluster_certificates": 5,
        "exact_zero_reserve_certificates": 2,
        "exact_finite_energy_certificates": 12,
        "analytic_terminal_tail_families": 1,
        "ideal_all_length_ternary_families": 1,
        "actual_all_length_ternary_families": 1,
        "open_small_prime_length_families": 1,
        "joined_abel_gaps": 0,
        "successor_winding_bounds": 0,
    }
    if summary != expected_summary:
        issues.append("summary drifted")
    row_ids = [row.get("id") for row in stored.get("rows", [])]
    if row_ids != EXPECTED_ROWS:
        issues.append("row IDs drifted")

    validate_sources(stored, issues)
    validate_symbolic_algebra(issues)
    validate_exact_inequalities(issues)
    validate_certificates(stored, issues)
    validate_actual_transfer(issues)

    exact = stored.get("exact", {})
    if exact.get("ideal_theorem", {}).get(
        "uniform"
    ) != "p=3,M>=5: J_0>17/1000":
        issues.append("ideal ternary theorem drifted")
    if exact.get("actual_transfer", {}).get(
        "actual_theorem"
    ) != "p=3,M>=5: J_ray>2/125":
        issues.append("actual ternary theorem drifted")
    boundary = stored.get("proof_boundary", "")
    for phrase in (
        "p=2,M>=10",
        "p=3 lengths 2..3",
        "Xi Abel gap",
        "Lambda<=0",
        "RH",
        "prize-level proof",
    ):
        if phrase not in boundary:
            issues.append(f"proof boundary missing: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    print(
        "validated ternary all-length prime-power heat-starlikeness gate: "
        "18 rows, 5 small-d clusters, 12 finite energies, "
        "1 analytic terminal tail, 1 ideal all-length ternary family, "
        "1 actual all-length ternary family, "
        "1 open small-prime length family, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
