#!/usr/bin/env python3
"""Check exact short small-prime heat-starlikeness counter-witnesses."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_short_family_counter_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_short_family_counter_gate.md"
)
SOURCES = {
    "base": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_base_certificate.json"
    ),
    "geometric": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "geometric_prime_power_phase_monotonicity_gate.json"
    ),
    "occupation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "complete_prime_power_chain_occupation_gate.json"
    ),
}
Q_STAR = Fraction(9999, 10000)
TARGET = Fraction(1, 100)
Pair = tuple[Fraction, Fraction]
EXPECTED = (
    (2, 2, Fraction(-99, 100)),
    (2, 3, Fraction(-1, 2)),
    (2, 4, Fraction(0)),
    (2, 5, Fraction(1, 3)),
    (2, 6, Fraction(-1, 2)),
    (2, 7, Fraction(-1, 4)),
    (3, 2, Fraction(-99, 100)),
    (3, 3, Fraction(-1, 2)),
)
EXPECTED_ROWS = [
    f"pphssf_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "current",
            "domain",
            "dyadic",
            "ternary",
            "exact_sign",
            "sharpness",
            "continuity",
            "route",
            "pi",
            "boundary",
        ),
        start=1,
    )
]


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/")
    return Fraction(int(numerator), int(denominator))


def parse_pair(payload: dict) -> Pair:
    return (
        parse_fraction(payload["rational"]),
        parse_fraction(payload["sqrt_prime_coefficient"]),
    )


def add(left: Pair, right: Pair) -> Pair:
    return left[0] + right[0], left[1] + right[1]


def multiply(left: Pair, right: Pair, prime: int) -> Pair:
    return (
        left[0] * right[0] + prime * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def scale(value: Pair, scalar: Fraction) -> Pair:
    return value[0] * scalar, value[1] * scalar


def coefficient(
    prime: int,
    length: int,
    index: int,
    q: Fraction,
    s: Fraction,
) -> Pair:
    scalar = q ** (index * (2 * length - 2 - index)) * s**index
    if index % 2 == 0:
        return scalar / prime ** (index // 2), Fraction(0)
    return Fraction(0), scalar / prime ** ((index + 1) // 2)


def chebyshev(maximum: int, x_value: Fraction) -> list[Fraction]:
    values = [Fraction(1)]
    if maximum:
        values.append(x_value)
    for _ in range(2, maximum + 1):
        values.append(2 * x_value * values[-1] - values[-2])
    return values


def direct_current(
    prime: int,
    length: int,
    q: Fraction,
    s: Fraction,
    x_value: Fraction,
) -> Pair:
    c = [coefficient(prime, length, k, q, s) for k in range(length)]
    t = chebyshev(length - 1, x_value)
    result = (Fraction(0), Fraction(0))
    for k in range(length):
        result = add(
            result,
            scale(multiply(c[k], c[k], prime), Fraction(k + 1)),
        )
    for k in range(length):
        for ell in range(k + 1, length):
            result = add(
                result,
                scale(
                    multiply(c[k], c[ell], prime),
                    Fraction(k + ell + 2) * t[ell - k],
                ),
            )
    return result


def is_negative(value: Pair, prime: int) -> bool:
    rational, radical = value
    if radical == 0:
        return rational < 0
    if radical < 0 and rational <= 0:
        return True
    if radical < 0 and rational > 0:
        return prime * radical**2 > rational**2
    if radical > 0 and rational < 0:
        return rational**2 > prime * radical**2
    return False


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_sources(payload: dict, issues: list[str]) -> None:
    stored = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in SOURCES.items():
        if not path.is_file():
            issues.append(f"missing source {key}")
        elif stored.get(key) != file_hash(path):
            issues.append(f"source hash drifted for {key}")


def validate_structure(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != (
        "prime_power_heat_starlikeness_short_family_counter_gate"
    ):
        issues.append("kind drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_ROWS:
        issues.append("row ids or order drifted")
    if any(
        row.get("readiness") != "proved" for row in rows[:-1]
    ):
        issues.append("a proved row drifted")
    if not rows or rows[-1].get("readiness") != "open":
        issues.append("boundary row is not open")
    expected_summary = {
        "rows": 10,
        "exact_interior_negative_witnesses": 8,
        "short_dyadic_lengths_rejected": 6,
        "short_ternary_lengths_rejected": 2,
        "uniform_short_block_positive_families": 0,
        "joined_phase_handoffs": 1,
        "joined_abel_gaps": 0,
        "successor_winding_bounds": 0,
    }
    summary = payload.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(f"summary drifted at {key}")


def validate_witnesses(payload: dict, issues: list[str]) -> None:
    rows = payload.get("witnesses", [])
    if len(rows) != len(EXPECTED):
        issues.append("wrong witness count")
        return
    boxes = {
        2: (Fraction(47, 50), Fraction(22, 25)),
        3: (Fraction(17, 20), Fraction(73, 100)),
    }
    for stored, expected in zip(rows, EXPECTED, strict=True):
        prime, length, x_value = expected
        if (
            stored.get("prime"),
            stored.get("block_length"),
            parse_fraction(stored.get("x_ang", "0/1")),
        ) != expected:
            issues.append(f"witness identity drifted at p={prime},M={length}")
            continue
        q_value = parse_fraction(stored["q"])
        s_value = parse_fraction(stored["s"])
        q_lower, s_lower = boxes[prime]
        if not q_lower < q_value < 1:
            issues.append(f"q is not interior at p={prime},M={length}")
        if not s_lower < s_value < 1:
            issues.append(f"s is not interior at p={prime},M={length}")
        if not -1 < x_value < 1:
            issues.append(f"x is not interior at p={prime},M={length}")
        rebuilt = direct_current(
            prime, length, q_value, s_value, x_value
        )
        if parse_pair(stored["current_pair"]) != rebuilt:
            issues.append(f"current pair drifted at p={prime},M={length}")
        shifted = add(rebuilt, (TARGET, Fraction(0)))
        if parse_pair(stored["shifted_pair"]) != shifted:
            issues.append(f"shifted pair drifted at p={prime},M={length}")
        if not is_negative(shifted, prime):
            issues.append(f"J_0<-1/100 failed at p={prime},M={length}")
        boundary = direct_current(
            prime,
            length,
            Fraction(1),
            Fraction(1),
            x_value,
        )
        if parse_pair(stored["boundary_q_s_one_pair"]) != boundary:
            issues.append(f"boundary pair drifted at p={prime},M={length}")
        if not is_negative(boundary, prime):
            issues.append(f"boundary sign failed at p={prime},M={length}")
    if not Q_STAR == parse_fraction(
        payload.get("assumptions", {}).get("q_star", "0/1")
    ):
        issues.append("q star drifted")
    if not Q_STAR == parse_fraction(
        payload.get("assumptions", {}).get("s_star", "0/1")
    ):
        issues.append("s star drifted")


def validate_note(payload: dict, issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append("missing note")
        return
    text = NOTE.read_text(encoding="utf-8")
    required = (
        "# Short Small-Prime Heat-Starlikeness Counter-Gate",
        "q=s=9999/10000",
        "p=2,M=2..7: universal blockwise J_0>0 is false",
        "p=3,M=2..3: universal blockwise J_0>0 is false",
        "ordinary `2*pi` period",
        payload.get("proof_boundary", ""),
    )
    for snippet in required:
        if snippet not in text:
            issues.append(f"note missing required text: {snippet[:70]}")
    for forbidden in (
        "This proves RH",
        "Lambda<=0 is proved",
        "an Xi ray attains every witness",
    ):
        if forbidden in text:
            issues.append(f"note contains forbidden overclaim: {forbidden}")


def main() -> int:
    issues: list[str] = []
    if not RESULT.is_file():
        issues.append("missing result")
        payload = {}
    else:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
    if payload:
        validate_structure(payload, issues)
        validate_sources(payload, issues)
        validate_witnesses(payload, issues)
        validate_note(payload, issues)
    if issues:
        print("short small-prime heat-starlikeness counter-gate failed:")
        for issue in issues:
            print(f"- {issue}")
        return 1
    print(
        "validated short small-prime heat-starlikeness counter-gate: "
        "10 rows, 8 exact interior negative witnesses, "
        "6 short dyadic lengths rejected, "
        "2 short ternary lengths rejected, "
        "0 uniform short-block positive families, "
        "1 joined-phase handoff, 0 joined Abel gaps, "
        "0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
