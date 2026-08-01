#!/usr/bin/env python3
"""Validate the dyadic all-length Fourier-defect localization gate."""

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
    "prime_power_heat_starlikeness_dyadic_defect_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
LENGTH_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
TERNARY_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_ternary_length_gate"
)
SOURCES = {
    "length_propagation": (
        REPO_ROOT / "work/rh_compute/results" / f"{LENGTH_STEM}.json"
    ),
    "ternary_length": (
        REPO_ROOT / "work/rh_compute/results" / f"{TERNARY_STEM}.json"
    ),
}
Q0 = Fraction(47, 50)
U = Fraction(707107, 1_000_000)
EXPECTED_ROWS = [
    "pphsdd_01_domain",
    "pphsdd_02_fourier_groups",
    "pphsdd_03_common_threshold",
    "pphsdd_04_cluster_slope",
    "pphsdd_05_d1_initial",
    "pphsdd_06_d1_m10",
    "pphsdd_07_d2_initial",
    "pphsdd_08_low_terminal",
    "pphsdd_09_middle_terminal",
    "pphsdd_10_terminal_tail",
    "pphsdd_11_localization",
    "pphsdd_12_two_defects",
    "pphsdd_13_boundary",
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
    second_power: int,
    coefficient: int | Fraction,
) -> None:
    if q_power < 0 or second_power < 0:
        raise ValueError("negative polynomial degree")
    key = (q_power, second_power)
    terms[key] = terms.get(key, Fraction()) + Fraction(coefficient)
    if not terms[key]:
        terms.pop(key)


def monomial_vector(
    power: int,
    degree: int,
    lower: Fraction,
    upper: Fraction,
) -> list[Fraction]:
    width = upper - lower
    result = []
    for index in range(degree + 1):
        value = Fraction()
        for j in range(min(power, index) + 1):
            value += (
                Fraction(comb(power, j))
                * lower ** (power - j)
                * width**j
                * Fraction(comb(index, j), comb(degree, j))
            )
        result.append(value)
    return result


def rebuild_certificate(
    source: dict[tuple[int, int], Fraction],
    second_interval: tuple[Fraction, Fraction],
    target: Fraction,
) -> dict:
    terms = dict(source)
    put(terms, 0, 0, -target)
    q_degree = max(q_power for q_power, _ in terms)
    second_degree = max(second_power for _, second_power in terms)
    q_vectors = {
        power: monomial_vector(power, q_degree, Q0, Fraction(1))
        for power in {q_power for q_power, _ in terms}
    }
    second_vectors = {
        power: monomial_vector(
            power,
            second_degree,
            second_interval[0],
            second_interval[1],
        )
        for power in {second_power for _, second_power in terms}
    }
    minimum: Fraction | None = None
    argminimum = None
    for qi in range(q_degree + 1):
        for si in range(second_degree + 1):
            value = Fraction()
            for (q_power, second_power), scalar in terms.items():
                value += (
                    scalar
                    * q_vectors[q_power][qi]
                    * second_vectors[second_power][si]
                )
            if minimum is None or value < minimum:
                minimum = value
                argminimum = [qi, si]
    if minimum is None or argminimum is None:
        raise RuntimeError("empty reconstructed certificate")
    return {
        "q_interval": ["47/50", "1"],
        "second_interval": [
            f"{second_interval[0].numerator}/{second_interval[0].denominator}",
            f"{second_interval[1].numerator}/{second_interval[1].denominator}",
        ],
        "power_degrees": [q_degree, second_degree],
        "bernstein_shape": [q_degree + 1, second_degree + 1],
        "claimed_lower": f"{target.numerator}/{target.denominator}",
        "minimum_shifted_coefficient": (
            f"{minimum.numerator}/{minimum.denominator}"
        ),
        "argminimum": argminimum,
        "positive": minimum > 0,
    }


def terminal_polynomial(
    B: int, d: int
) -> dict[tuple[int, int], Fraction]:
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


def initial_polynomial(
    d: int, count: int
) -> tuple[dict[tuple[int, int], Fraction], int]:
    entries: list[tuple[int, int, int]] = []
    minimum = 0
    for k in range(count):
        A = 2 * k + d + 2
        group = (
            (2 * k * (d + 1 - k), 2 * k, A),
            (2 * k * (d - k), 2 * k + 1, -2 * (A + 1)),
            (
                2 * k * (d - 1 - k) - 2,
                2 * k + 2,
                A + 2,
            ),
        )
        entries.extend(group)
        minimum = min(minimum, *(item[0] for item in group))
    shift = -minimum
    result: dict[tuple[int, int], Fraction] = {}
    for q_power, ratio_power, scalar in entries:
        put(result, q_power + shift, ratio_power, scalar)
    return result, shift


def normalized_difference_polynomial(
    M: int, d: int
) -> dict[tuple[int, int], Fraction]:
    exponent = [k * (2 * M - 2 - k) for k in range(M)]
    result: dict[tuple[int, int], Fraction] = {}
    for increment, scale in ((0, 1), (1, -2), (2, 1)):
        offset = d + increment
        for k in range(M - offset):
            put(
                result,
                exponent[k] + exponent[k + offset] - exponent[d],
                2 * k + increment,
                scale * (2 * k + offset + 2),
            )
    return result


def compare_certificate(
    stored: dict,
    rebuilt: dict,
    label: str,
    issues: list[str],
) -> None:
    for key in (
        "q_interval",
        "second_interval",
        "power_degrees",
        "bernstein_shape",
        "claimed_lower",
        "minimum_shifted_coefficient",
        "argminimum",
        "positive",
    ):
        if stored.get(key) != rebuilt.get(key):
            issues.append(f"{label} certificate drifted at {key}")


def validate_certificates(payload: dict, issues: list[str]) -> None:
    low = payload.get("exact", {}).get("low_offsets", {})
    d1_terms, d1_shift = initial_polynomial(1, 7)
    d2_terms, d2_shift = initial_polynomial(2, 3)
    d1 = low.get("d1_initial_certificate", {})
    d2 = low.get("d2_initial_certificate", {})
    if d1.get("q_shift") != d1_shift:
        issues.append("d=1 initial-cluster shift drifted")
    if d2.get("q_shift") != d2_shift:
        issues.append("d=2 initial-cluster shift drifted")
    compare_certificate(
        d1,
        rebuild_certificate(
            d1_terms, (Fraction(3, 5), U), Fraction(1, 200)
        ),
        "d=1 initial",
        issues,
    )
    compare_certificate(
        d2,
        rebuild_certificate(
            d2_terms, (Fraction(2, 3), U), Fraction(1, 25)
        ),
        "d=2 initial",
        issues,
    )
    compare_certificate(
        low.get("d1_m10_certificate", {}),
        rebuild_certificate(
            normalized_difference_polynomial(10, 1),
            (Fraction(0), U),
            Fraction(1, 100),
        ),
        "d=1 M=10",
        issues,
    )
    for label, B, d in (
        ("d=1 terminal", 15, 1),
        ("d=2 terminal", 14, 2),
    ):
        compare_certificate(
            low.get(
                "d1_terminal_certificate"
                if d == 1
                else "d2_terminal_certificate",
                {},
            ),
            rebuild_certificate(
                terminal_polynomial(B, d),
                (Fraction(0), U),
                Fraction(1),
            ),
            label,
            issues,
        )

    middle = (
        payload.get("exact", {})
        .get("middle_offsets", {})
        .get("certificates", [])
    )
    if len(middle) != 12:
        issues.append("wrong middle-terminal certificate count")
    else:
        for d, stored in zip(range(3, 15), middle, strict=True):
            if stored.get("d") != d:
                issues.append(f"middle-terminal d={d} label drifted")
            compare_certificate(
                stored,
                rebuild_certificate(
                    terminal_polynomial(d + 2, d),
                    (Fraction(0), U),
                    Fraction(1, 25),
                ),
                f"middle-terminal d={d}",
                issues,
            )


def validate_symbolic_algebra(issues: list[str]) -> None:
    a, q, R = sp.symbols("a q R", positive=True)
    M = 10
    coefficients = [
        a**k * q ** (k * (2 * M - 2 - k)) for k in range(M)
    ]
    A = [sp.Integer(0)] * (M + 2)
    for d in range(M):
        A[d] = sum(
            (2 * k + d + 2)
            * coefficients[k]
            * coefficients[k + d]
            for k in range(M - d)
        )
    D1 = sp.expand(A[1] - 2 * A[2] + A[3])
    rebuilt = sum(
        sp.Rational(value.numerator, value.denominator)
        * q**q_power
        * a**a_power
        for (q_power, a_power), value in
        normalized_difference_polynomial(M, 1).items()
    )
    if sp.expand(D1 / coefficients[1] - rebuilt) != 0:
        issues.append("normalized M=10 D1 identity failed")

    for d, count in ((1, 7), (2, 3)):
        direct = sp.Integer(0)
        for k in range(count):
            A_value = 2 * k + d + 2
            weight = R ** (2 * k) * q ** (2 * k * (d + 1 - k))
            ratio = R * q ** (-2 * k)
            direct += weight * (
                A_value
                - 2 * (A_value + 1) * ratio
                + (A_value + 2) * ratio**2 * q**-2
            )
        terms, shift = initial_polynomial(d, count)
        polynomial = sum(
            sp.Rational(value.numerator, value.denominator)
            * q**q_power
            * R**ratio_power
            for (q_power, ratio_power), value in terms.items()
        )
        if sp.expand(q**shift * direct - polynomial) != 0:
            issues.append(f"d={d} initial-cluster identity failed")

    B, d_symbol = sp.symbols("B d", integer=True, positive=True)
    cluster = (
        B
        - 2 * (B + 1) * a * q**3
        + (B + 2) * a**2 * q**4
        + a**2
        * q ** (2 * d_symbol + 6)
        * (
            B
            + 2
            - 2 * (B + 3) * a * q
            + (B + 4) * a**2 * q ** (2 * d_symbol + 2)
        )
    )
    slope = sp.diff(cluster, B)
    expected = (
        1
        - 2 * a * q**3
        + a**2 * q**4
        + a**2
        * q ** (2 * d_symbol + 6)
        * (1 - 2 * a * q + a**2 * q ** (2 * d_symbol + 2))
    )
    if sp.expand(slope - expected) != 0:
        issues.append("terminal-cluster slope identity failed")

    x = sp.symbols("x", nonnegative=True)
    minimal = (
        d_symbol
        + 2
        - 2 * (d_symbol + 3) * a * q**3
        + (d_symbol + 4) * a**2 * q**4
        + a**2
        * q**4
        * x
        * (
            d_symbol
            + 4
            - 2 * (d_symbol + 5) * a * q
            + (d_symbol + 6) * a**2 * x
        )
    )
    E = (
        1
        - 2 * a * q**3
        + a**2 * q**4
        + a**2 * q**4 * x * (1 - 2 * a * q + a**2 * x)
    )
    K = (
        2
        - 6 * a * q**3
        + 4 * a**2 * q**4
        + a**2 * q**4 * x * (4 - 10 * a * q + 6 * a**2 * x)
    )
    if sp.expand(minimal - d_symbol * E - K) != 0:
        issues.append("terminal dE+K decomposition failed")


def validate_exact_guards(issues: list[str]) -> None:
    if not 2 * U**2 > 1:
        issues.append("rational reciprocal-sqrt(2) upper failed")
    if not 7 * U < 5:
        issues.append("common-group 1/sqrt(2)<5/7 guard failed")
    if not Fraction(71, 25) ** 2 > 8:
        issues.append("slope delta>1/25 guard failed")
    if 75**2 - 2 * 53**2 != 7:
        issues.append("analytic tail endpoint identity failed")
    if not 75**2 > 2 * 53**2:
        issues.append("analytic tail endpoint sign failed")


def validate_structure(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != (
        "prime_power_heat_starlikeness_dyadic_defect_gate"
    ):
        issues.append("kind drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_ROWS:
        issues.append("row ids or order drifted")
    readiness = {
        row.get("id"): row.get("readiness") for row in rows
    }
    for row_id in EXPECTED_ROWS[:11]:
        if readiness.get(row_id) != "proved":
            issues.append(f"{row_id} is not proved")
    for row_id in EXPECTED_ROWS[11:]:
        if readiness.get(row_id) != "open":
            issues.append(f"{row_id} is not open")

    summary = payload.get("summary", {})
    expected_summary = {
        "rows": 13,
        "exact_initial_cluster_certificates": 2,
        "exact_shortest_length_certificates": 1,
        "exact_low_terminal_certificates": 2,
        "exact_middle_terminal_certificates": 12,
        "analytic_terminal_tail_families": 1,
        "proved_positive_offset_families": 1,
        "remaining_possible_defects": 2,
        "ideal_all_length_dyadic_families": 0,
        "actual_all_length_dyadic_families": 0,
        "joined_abel_gaps": 0,
        "successor_winding_bounds": 0,
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(f"summary drifted at {key}")

    localization = payload.get("exact", {}).get("localization", {})
    if localization.get("possible_negative_indices") != ["0", "M-2"]:
        issues.append("possible-defect list drifted")
    if localization.get("theorem") != (
        "p=2,M>=10: D_d>0 for 1<=d<=M-3; D_(M-1)>0"
    ):
        issues.append("localization theorem drifted")


def validate_sources(payload: dict, issues: list[str]) -> None:
    audit = payload.get("source_audit", {})
    stored_hashes = audit.get("source_sha256", {})
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
        "# Dyadic All-Length Fourier-Defect Localization Gate",
        "possible negative indices: D_0 and D_(M-2)",
        "75^2-2*53^2=7>0",
        payload.get("proof_boundary", ""),
    )
    for snippet in required:
        if snippet not in text:
            issues.append(f"note missing required text: {snippet[:80]}")
    forbidden = (
        "This proves RH",
        "Lambda<=0 is proved",
        "PF-infinity is proved",
        "uniform dyadic heat-starlikeness is proved",
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
        validate_symbolic_algebra(issues)
        validate_exact_guards(issues)
        validate_note(payload, issues)
    if issues:
        print(
            "dyadic all-length Fourier-defect localization validation failed:"
        )
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated dyadic all-length Fourier-defect localization gate: "
        "13 rows, 2 initial clusters, 1 shortest-length certificate, "
        "2 low terminal anchors, 12 middle terminal clusters, "
        "1 analytic terminal tail, 1 proved positive-offset family, "
        "2 possible defects, 0 ideal all-length dyadic families, "
        "0 actual all-length dyadic families, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
