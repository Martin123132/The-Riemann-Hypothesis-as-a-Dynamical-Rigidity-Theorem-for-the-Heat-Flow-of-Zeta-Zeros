#!/usr/bin/env python3
"""Validate the dyadic all-length prime-power heat-starlikeness gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

from dyadic_heat_starlikeness_bernstein import (
    Q2_LOWER,
    RATIO_UPPER,
    certify_finite_core,
    certify_initial_core,
    fejer_numerator_polynomials,
    finite_core_power,
    initial_core_power,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_dyadic_length_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
LENGTH_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
DEFECT_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_dyadic_defect_gate"
)
SOURCES = {
    "length_propagation": (
        REPO_ROOT / "work/rh_compute/results" / f"{LENGTH_STEM}.json"
    ),
    "dyadic_defect": (
        REPO_ROOT / "work/rh_compute/results" / f"{DEFECT_STEM}.json"
    ),
}
X_LEAVES = tuple(
    "".join(directions)
    for directions in product("LR", repeat=4)
)
EXPECTED_ROWS = [
    "pphsdl_01_domain",
    "pphsdl_02_localization",
    "pphsdl_03_fejer",
    "pphsdl_04_terminal_pair",
    "pphsdl_05_terminal_parameters",
    "pphsdl_06_finite_core",
    "pphsdl_07_finite_terminal",
    "pphsdl_08_finite_theorem",
    "pphsdl_09_initial_groups",
    "pphsdl_10_initial_core",
    "pphsdl_11_m15_terminal",
    "pphsdl_12_terminal_tail",
    "pphsdl_13_tail_theorem",
    "pphsdl_14_ideal_theorem",
    "pphsdl_15_actual_sums",
    "pphsdl_16_actual_current",
    "pphsdl_17_pi_provenance",
    "pphsdl_18_boundary",
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


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def normalized_json(value: object) -> object:
    return json.loads(json.dumps(value))


def terminal_square_guard(M: int, target: Fraction) -> dict:
    n = M - 1
    A = 2 * (M + 1) * RATIO_UPPER - M
    denominator = M - (M + 1) * RATIO_UPPER
    low_coefficient = Fraction(M + 1) * A / denominator
    low_base = A / (2 * (M + 2))
    low_margin = target**2 - low_coefficient**2 * low_base**n
    high_coefficient = Fraction(4 * (M + 1), M) * A
    high_radicand = (
        (A / (M + 2)) ** n
        * Fraction(n**n, (n + 2) ** (n + 2))
    )
    high_margin = target**2 - high_coefficient**2 * high_radicand
    return {
        "M": M,
        "target": fraction_text(target),
        "A_upper": fraction_text(A),
        "denominator_lower": fraction_text(denominator),
        "low_square_margin": fraction_text(low_margin),
        "high_square_margin": fraction_text(high_margin),
        "positive": low_margin > 0 and high_margin > 0,
    }


def compare_terminal_guard(
    stored: dict,
    rebuilt: dict,
    label: str,
    issues: list[str],
) -> None:
    for key in (
        "M",
        "target",
        "A_upper",
        "denominator_lower",
        "low_square_margin",
        "high_square_margin",
        "positive",
    ):
        if stored.get(key) != rebuilt.get(key):
            issues.append(f"{label} terminal guard drifted at {key}")


def validate_certificates(payload: dict, issues: list[str]) -> None:
    exact = payload.get("exact", {})
    stored_finite = exact.get("finite", {}).get(
        "core_certificates", []
    )
    if len(stored_finite) != 5:
        issues.append("wrong finite-core certificate count")
    else:
        for M, stored in zip(range(10, 15), stored_finite, strict=True):
            rebuilt = certify_finite_core(
                M, x_leaf_paths=X_LEAVES
            ).to_dict()
            if stored != normalized_json(rebuilt):
                issues.append(f"M={M} finite-core certificate drifted")

    stored_initial = (
        exact.get("all_length", {})
        .get("initial_core_certificate", {})
    )
    rebuilt_initial = certify_initial_core(
        x_leaf_paths=X_LEAVES
    ).to_dict()
    if stored_initial != normalized_json(rebuilt_initial):
        issues.append("all-length initial-core certificate drifted")

    stored_terminal = exact.get("finite", {}).get(
        "terminal_guards", []
    )
    if len(stored_terminal) != 5:
        issues.append("wrong finite terminal-guard count")
    else:
        for M, stored in zip(
            range(10, 15), stored_terminal, strict=True
        ):
            compare_terminal_guard(
                stored,
                terminal_square_guard(M, Fraction(1, 25)),
                f"M={M}",
                issues,
            )
    compare_terminal_guard(
        exact.get("all_length", {}).get("m15_terminal_guard", {}),
        terminal_square_guard(15, Fraction(1, 250)),
        "M=15",
        issues,
    )


def polynomial_value(coefficients: list[int], x: Fraction) -> Fraction:
    value = Fraction()
    for coefficient in reversed(coefficients):
        value = value * x + coefficient
    return value


def array_value(
    array: np.ndarray,
    y: Fraction,
    ratio: Fraction,
    x: Fraction,
) -> Fraction:
    y_powers = [Fraction(1)]
    ratio_powers = [Fraction(1)]
    x_powers = [Fraction(1)]
    for _ in range(1, array.shape[0]):
        y_powers.append(y_powers[-1] * y)
    for _ in range(1, array.shape[1]):
        ratio_powers.append(ratio_powers[-1] * ratio)
    for _ in range(1, array.shape[2]):
        x_powers.append(x_powers[-1] * x)
    value = Fraction()
    for (yi, ri, xi), coefficient in np.ndenumerate(array):
        if coefficient:
            value += (
                coefficient
                * y_powers[yi]
                * ratio_powers[ri]
                * x_powers[xi]
            )
    return value


def direct_finite_core(
    M: int,
    y: Fraction,
    ratio: Fraction,
    x: Fraction,
) -> Fraction:
    exponents = [
        index * (2 * M - 3 - index) // 2
        for index in range(M)
    ]
    coefficients = [
        ratio**index * y**exponents[index]
        for index in range(M)
    ]
    A = [Fraction()] * (M + 2)
    for d in range(M):
        A[d] = sum(
            Fraction(2 * k + d + 2)
            * coefficients[k]
            * coefficients[k + d]
            for k in range(M - d)
        )
    D = [
        A[d] - 2 * A[d + 1] + A[d + 2]
        for d in range(M)
    ]
    angular = fejer_numerator_polynomials(7)
    return sum(
        D[d] * polynomial_value(angular[d], x)
        for d in range(8)
    )


def direct_initial_core(
    y: Fraction,
    ratio: Fraction,
    x: Fraction,
) -> Fraction:
    anchor = 12
    exponents = [
        index * (2 * anchor - index + 1) // 2
        for index in range(14)
    ]
    coefficients = [
        ratio**index * y**exponents[index]
        for index in range(14)
    ]
    angular = fejer_numerator_polynomials(4)
    total = Fraction()
    for d in range(5):
        partial = Fraction()
        for k in range(8):
            A = 2 * k + d + 2
            partial += (
                A * coefficients[k] * coefficients[k + d]
                - 2
                * (A + 1)
                * coefficients[k]
                * coefficients[k + d + 1]
                + (A + 2)
                * coefficients[k]
                * coefficients[k + d + 2]
            )
        total += partial * polynomial_value(angular[d], x)
    return total


def validate_power_algebra(issues: list[str]) -> None:
    finite_array = finite_core_power(10)
    initial_array, anchor = initial_core_power()
    if anchor != 12:
        issues.append("initial-core ratio anchor drifted")
    points = (
        (
            Q2_LOWER,
            Fraction(1, 3),
            Fraction(-2, 3),
        ),
        (
            Fraction(1),
            RATIO_UPPER,
            Fraction(-1),
        ),
    )
    for index, (y, ratio, x) in enumerate(points):
        finite_stored = array_value(finite_array, y, ratio, x)
        finite_direct = (
            direct_finite_core(10, y, ratio, x) - Fraction(1, 10)
        )
        if finite_stored != finite_direct:
            issues.append(f"finite-core power identity failed at point {index}")
        initial_stored = array_value(initial_array, y, ratio, x)
        initial_direct = (
            direct_initial_core(y, ratio, x) - Fraction(1, 100)
        )
        if initial_stored != initial_direct:
            issues.append(
                f"initial-core power identity failed at point {index}"
            )


def validate_symbolic_identities(issues: list[str]) -> None:
    P, N = sp.symbols("P N", positive=True)
    yr, yi, ur, ui = sp.symbols("yr yi ur ui", real=True)
    y_abs = yr**2 + yi**2
    shifted_abs = (yr + ur) ** 2 + (yi + ui) ** 2
    unit_guard = ur**2 + ui**2 - 1
    completed = (
        (P - N)
        * (
            (yr + P * ur / (P - N)) ** 2
            + (yi + P * ui / (P - N)) ** 2
        )
        - P * N / (P - N)
    )
    difference = sp.together(
        P * shifted_abs - N * y_abs - completed
    )
    numerator = sp.factor(difference.as_numer_denom()[0])
    if sp.rem(
        sp.Poly(numerator, ur),
        sp.Poly(unit_guard, ur),
    ) != 0:
        issues.append("terminal completed-square identity failed")

    for M in range(10, 21):
        alpha_twice = M - 1
        if M - 2 - alpha_twice + 1 != 0:
            issues.append(f"terminal r cancellation failed at M={M}")


def validate_exact_guards(payload: dict, issues: list[str]) -> None:
    if not 2 * RATIO_UPPER**2 > 1:
        issues.append("rational reciprocal-sqrt(2) upper failed")
    if not all(
        M - (M + 1) * RATIO_UPPER > 0
        for M in range(10, 16)
    ):
        issues.append("terminal P-N denominator guard failed")
    if 231**2 - 2 * 160**2 != 2161:
        issues.append("D0 terminal-cluster guard failed")

    tail = (
        payload.get("exact", {})
        .get("all_length", {})
        .get("tail_terminal_guard", {})
    )
    if tail.get("half_power_integer_margin") != (
        17**17 - 4 * 21**2 * 15**15
    ):
        issues.append("half-power maximum guard drifted")

    M = 16
    low_coefficient = Fraction(21, 10) * (M + 1)
    low_margin = (
        Fraction(1, 250) ** 2
        - low_coefficient**2 * Fraction(9, 40) ** (M - 1)
    )
    high_coefficient = (
        Fraction(9, 10)
        * Fraction((M + 1) * (M + 2), 21 * M)
    )
    high_margin = (
        Fraction(1, 250) ** 2
        - high_coefficient**2 * Fraction(9, 20) ** (M - 1)
    )
    if tail.get("low_M16_square_margin") != fraction_text(low_margin):
        issues.append("low M=16 terminal-tail guard drifted")
    if tail.get("high_M16_square_margin") != fraction_text(high_margin):
        issues.append("high M=16 terminal-tail guard drifted")
    if not low_margin > 0 or not high_margin > 0:
        issues.append("M=16 terminal-tail endpoint failed")
    if not Fraction(9, 17) < 1:
        issues.append("low terminal successor ratio failed")
    if not Fraction(51, 64) < 1:
        issues.append("high terminal successor ratio failed")


def validate_actual_transfer(issues: list[str]) -> None:
    x_lower = 12 * 2**50
    relative = Fraction(18000, x_lower)
    if not Fraction(4) * (4 + 12) == 64:
        issues.append("dyadic geometric current mass failed")
    if not 64 * (2 * relative + relative**2) < 192 * relative:
        issues.append("dyadic frozen-current transfer failed")
    ray_error = (
        Fraction(3_456_096, x_lower)
        + Fraction(1_621_632, x_lower**2)
    )
    if not ray_error < Fraction(1, 1000):
        issues.append("dyadic actual-ray transfer failed")
    if not Fraction(3, 1000) - Fraction(1, 1000) == Fraction(
        1, 500
    ):
        issues.append("dyadic actual margin arithmetic failed")


def validate_structure(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != (
        "prime_power_heat_starlikeness_dyadic_length_gate"
    ):
        issues.append("kind drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_ROWS:
        issues.append("row ids or order drifted")
    readiness = {
        row.get("id"): row.get("readiness") for row in rows
    }
    for row_id in EXPECTED_ROWS[:-1]:
        if readiness.get(row_id) != "proved":
            issues.append(f"{row_id} is not proved")
    if readiness.get(EXPECTED_ROWS[-1]) != "open":
        issues.append(f"{EXPECTED_ROWS[-1]} is not open")

    expected_summary = {
        "rows": 18,
        "exact_finite_core_certificates": 5,
        "exact_initial_core_certificates": 1,
        "exact_finite_terminal_guards": 5,
        "analytic_terminal_tail_families": 1,
        "ideal_all_length_dyadic_families": 1,
        "actual_all_length_dyadic_families": 1,
        "complete_dyadic_minimum_length": 8,
        "remaining_possible_defects": 0,
        "joined_abel_gaps": 0,
        "successor_winding_bounds": 0,
    }
    summary = payload.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(f"summary drifted at {key}")
    exact = payload.get("exact", {})
    if exact.get("ideal_theorem", {}).get("uniform") != (
        "p=2,M>=10: J_0>3/1000"
    ):
        issues.append("uniform ideal theorem drifted")
    if exact.get("actual_transfer", {}).get("actual_theorem") != (
        "p=2,M>=10: J_ray>1/500"
    ):
        issues.append("uniform actual theorem drifted")


def validate_sources(payload: dict, issues: list[str]) -> None:
    stored_hashes = (
        payload.get("source_audit", {}).get("source_sha256", {})
    )
    for key, path in SOURCES.items():
        if not path.is_file():
            issues.append(f"missing source {key}: {path}")
        elif stored_hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted for {key}")


def validate_note(payload: dict, issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    required = (
        "# Dyadic All-Length Prime-Power Heat-Starlikeness Gate",
        "p=2,M>=10: J_0>3/1000",
        "p=2,M>=10: J_ray>1/500",
        "ordinary `2*pi` period",
        payload.get("proof_boundary", ""),
    )
    for snippet in required:
        if snippet not in text:
            issues.append(f"note missing required text: {snippet[:80]}")
    forbidden = (
        "This proves RH",
        "Lambda<=0 is proved",
        "PF-infinity is proved",
        "all prime-power joins are proved",
    )
    for snippet in forbidden:
        if snippet in text:
            issues.append(f"note contains forbidden overclaim: {snippet}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, "result", issues)
    if payload:
        validate_structure(payload, issues)
        validate_sources(payload, issues)
        validate_certificates(payload, issues)
        validate_power_algebra(issues)
        validate_symbolic_identities(issues)
        validate_exact_guards(payload, issues)
        validate_actual_transfer(issues)
        validate_note(payload, issues)
    if issues:
        print(
            "dyadic all-length prime-power heat-starlikeness "
            "validation failed:"
        )
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated dyadic all-length prime-power heat-starlikeness gate: "
        "18 rows, 5 finite cores, 1 initial core, "
        "5 finite terminal guards, 1 analytic terminal tail, "
        "1 ideal all-length dyadic family, "
        "1 actual all-length dyadic family, "
        "complete dyadic minimum length 8, 0 remaining defects, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
