#!/usr/bin/env python3
"""Validate the prime-power heat-starlikeness length-propagation gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from itertools import product
from pathlib import Path

import sympy as sp

from prime_power_heat_bernstein import certify_family


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_length_propagation_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
BASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_heat_starlikeness_base_certificate"
)
PHASE_STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "prime_power_logarithmic_phase_flow_gate"
)
SOURCES = {
    "base_certificate": (
        REPO_ROOT / "work/rh_compute/results" / f"{BASE_STEM}.json"
    ),
    "phase_flow": (
        REPO_ROOT / "work/rh_compute/results" / f"{PHASE_STEM}.json"
    ),
}
EXPECTED_IDS = [
    "pphslp_01_domain",
    "pphslp_02_append",
    "pphslp_03_current_recurrence",
    "pphslp_04_contracted_old",
    "pphslp_05_dyadic_m9_box",
    "pphslp_06_dyadic_m9_bernstein",
    "pphslp_07_fourier",
    "pphslp_08_ratios",
    "pphslp_09_monotone",
    "pphslp_10_interior_convexity",
    "pphslp_11_boundary_reduction",
    "pphslp_12_boundary_minimization",
    "pphslp_13_finite_boundary",
    "pphslp_14_tail_boundary",
    "pphslp_15_fejer",
    "pphslp_16_first_difference",
    "pphslp_17_large_prime_ideal",
    "pphslp_18_uniform_transfer",
    "pphslp_19_actual_margins",
    "pphslp_20_boundary",
]


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def rebuild_m9_certificate() -> dict:
    paths = tuple(
        "".join(directions)
        for directions in product("LR", repeat=5)
    )
    return certify_family(
        prime=2,
        block_length=9,
        q_lower=Fraction(47, 50),
        s_lower=Fraction(22, 25),
        x_leaf_paths=paths,
        claimed_lower=Fraction(1, 20),
    ).to_dict()


def validate_append_identity(issues: list[str]) -> None:
    length = 5
    z = sp.symbols("z", nonzero=True)
    q, d = sp.symbols("q d", positive=True, real=True)
    coefficients = sp.symbols(
        f"c0:{length}", positive=True, real=True
    )
    old_scaled = [
        coefficients[index] * q ** (2 * index)
        for index in range(length)
    ]
    q_old = sum(
        coefficient * z**index
        for index, coefficient in enumerate(old_scaled)
    )
    h_old = sum(
        index * coefficient * z**index
        for index, coefficient in enumerate(old_scaled)
    )
    q_new = q_old + d * z**length
    h_new = h_old + length * d * z**length

    def conjugate_laurent(expression: sp.Expr) -> sp.Expr:
        return sp.expand(expression.xreplace({z: 1 / z}))

    direct = sp.expand(
        q_new * conjugate_laurent(q_new)
        + (
            h_new * conjugate_laurent(q_new)
            + q_new * conjugate_laurent(h_new)
        )
        / 2
    )
    old_current = sp.expand(
        q_old * conjugate_laurent(q_old)
        + (
            h_old * conjugate_laurent(q_old)
            + q_old * conjugate_laurent(h_old)
        )
        / 2
    )
    cross = sum(
        d
        * (length + index + 2)
        * old_scaled[index]
        * (
            z ** (length - index)
            + z ** (index - length)
        )
        / 2
        for index in range(length)
    )
    expected = sp.expand(
        old_current + (length + 1) * d**2 + cross
    )
    if sp.expand(direct - expected) != 0:
        issues.append("append-one current identity failed")

    p, s = sp.symbols("p s", positive=True)
    for index in range(length):
        old = (
            p ** (-sp.Rational(index, 2))
            * q ** (index * (2 * length - 2 - index))
            * s**index
        )
        new = (
            p ** (-sp.Rational(index, 2))
            * q ** (index * (2 * (length + 1) - 2 - index))
            * s**index
        )
        if sp.simplify(new / old - q ** (2 * index)) != 0:
            issues.append(f"append coefficient failed at k={index}")
    new_term = (
        p ** (-sp.Rational(length, 2))
        * q
        ** (
            length
            * (2 * (length + 1) - 2 - length)
        )
        * s**length
    )
    expected_term = (
        p ** (-sp.Rational(length, 2))
        * q ** (length**2)
        * s**length
    )
    if sp.simplify(new_term - expected_term) != 0:
        issues.append("new append coefficient failed")


def validate_fejer_identity(issues: list[str]) -> None:
    degree = 6
    z = sp.symbols("z", nonzero=True)
    values = sp.symbols(f"A0:{degree + 2}", real=True)
    substitutions = {
        values[degree]: sp.Integer(0),
        values[degree + 1]: sp.Integer(0),
    }
    target = values[0] / 2 + sum(
        values[index] * (z**index + z ** (-index)) / 2
        for index in range(1, degree)
    )
    reconstructed = 0
    for index in range(degree):
        second = (
            values[index]
            - 2 * values[index + 1]
            + values[index + 2]
        ).subs(substitutions)
        geometric = sum(z**power for power in range(index + 1))
        kernel = sp.expand(
            geometric
            * sum(z ** (-power) for power in range(index + 1))
            / (index + 1)
        )
        reconstructed += (
            sp.Rational(index + 1, 2) * second * kernel
        )
    difference = sp.expand(
        (reconstructed - target).subs(substitutions)
    )
    if difference != 0:
        issues.append("Fejer-kernel reconstruction failed")


def validate_large_prime_algebra(issues: list[str]) -> None:
    # 1/sqrt(5)<1/2 and 2/sqrt(5)<9/10.
    if not Fraction(1, 5) < Fraction(1, 4):
        issues.append("large-prime half-ratio bound failed")
    if not Fraction(4, 5) < Fraction(81, 100):
        issues.append("large-prime Bernstein radical bound failed")

    a_value, ratio, next_ratio = sp.symbols(
        "A r r_plus", positive=True
    )
    common = (
        a_value
        - 2 * (a_value + 1) * ratio
        + (a_value + 2) * ratio * next_ratio
    )
    lower = (1 - ratio) * (
        a_value - (a_value + 2) * ratio
    )
    if sp.expand(common.subs(next_ratio, ratio) - lower) != 0:
        issues.append("interior convexity factorization failed")

    # The a-derivative is negative because
    # (d+4)^2<5(d+3)^2 for every d>=0.
    d = sp.symbols("d", nonnegative=True, integer=True)
    derivative_gap = sp.expand(
        5 * (d + 3) ** 2 - (d + 4) ** 2
    )
    if derivative_gap != 4 * d**2 + 22 * d + 29:
        issues.append("terminal a-monotonicity polynomial failed")

    # For d<=6, certify every ordinary degree-n Bernstein coefficient
    # with the rational upper bound 2/sqrt(5)<9/10.
    for d_value in range(7):
        degree = 2 * d_value + 2
        for index in range(degree):
            lower_coefficient = Fraction(d_value + 2) - (
                Fraction(9, 10)
                * (d_value + 3)
                * Fraction(index, degree)
            )
            if not lower_coefficient > 0:
                issues.append(
                    "small-d Bernstein coefficient failed: "
                    f"d={d_value}, index={index}"
                )
        final_lower = Fraction(3 * d_value + 1, 10)
        if not final_lower > 0:
            issues.append(
                f"small-d final Bernstein coefficient failed: d={d_value}"
            )

    # For d>=7, the dropped-positive-term lower bound is
    # (d-7)/10, with strictness at d=7 inherited from 2/sqrt(5)<9/10.
    if not Fraction(7 + 2) - Fraction(9, 10) * (7 + 3) == 0:
        issues.append("large-d boundary threshold failed")

    # 14/5-6/sqrt(5)>1/10 is equivalent to 729>720.
    if not Fraction(729, 100) > Fraction(36, 5):
        issues.append("uniform first-second-difference margin failed")


def validate_actual_transfer(issues: list[str]) -> None:
    x_lower = 12 * 2**50
    relative = Fraction(18000, x_lower)
    if not relative < 1:
        issues.append("relative perturbation is not contractive")

    # k log(p)<=log(N), t log(N)<51/4, and delta_a<1/(4x)
    # give eta<51/(32x)<2/x.
    if not Fraction(51, 32) < 2:
        issues.append("uniform saddle-offset exponent failed")

    # For p>=5 use 1/sqrt(5)<9/20.
    coefficient_sum = Fraction(20, 11)
    moment_sum = Fraction(180, 121)
    absolute_current_sum = coefficient_sum * (
        coefficient_sum + moment_sum
    )
    if not absolute_current_sum < 7:
        issues.append("large-prime geometric current sum failed")
    if not 7 * (2 * relative + relative**2) < 21 * relative:
        issues.append("large-prime frozen-current transfer failed")

    large_prime_error = (
        Fraction(378_006, x_lower)
        + Fraction(405_408, x_lower**2)
    )
    if not large_prime_error < Fraction(1, 1000):
        issues.append("large-prime fixed-ray transfer failed")

    # For p=2,M=9 use 1/sqrt(2)<3/4, hence sum c_k<4
    # and sum k c_k<12.
    dyadic_current_sum = Fraction(4) * (4 + 12)
    if dyadic_current_sum != 64:
        issues.append("dyadic M=9 geometric current sum failed")
    if not 64 * (2 * relative + relative**2) < 192 * relative:
        issues.append("dyadic M=9 frozen-current transfer failed")
    dyadic_error = (
        Fraction(3_456_096, x_lower)
        + Fraction(1_621_632, x_lower**2)
    )
    if not dyadic_error < Fraction(1, 1000):
        issues.append("dyadic M=9 fixed-ray transfer failed")

    if not Fraction(1, 20) - Fraction(1, 1000) == Fraction(
        49, 1000
    ):
        issues.append("propagated actual margin failed")


def validate_sources(
    stored: dict, sources: dict[str, dict], issues: list[str]
) -> None:
    expected_hashes = {
        key: file_hash(path)
        for key, path in SOURCES.items()
        if path.is_file()
    }
    if stored.get("source_audit", {}).get("source_sha256") != expected_hashes:
        issues.append("source hashes drifted")
    if sources["base_certificate"].get("summary", {}).get(
        "length_propagation_theorems"
    ) != 0:
        issues.append("base source boundary drifted")
    if "J_heat=|Q_heat|^2" not in sources["phase_flow"].get(
        "exact", {}
    ).get("replacement_target", {}).get("current", ""):
        issues.append("phase-flow source target drifted")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "gate result", issues)
    sources = {
        key: load_json(path, f"{key} source", issues)
        for key, path in SOURCES.items()
    }
    if issues:
        return issues

    if stored.get("kind") != STEM:
        issues.append("result kind drifted")
    if stored.get("date") != "2026-07-30":
        issues.append("result date drifted")
    if [row.get("id") for row in stored.get("rows", [])] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    expected_summary = {
        "rows": 20,
        "append_one_identities": 1,
        "exact_dyadic_next_length_certificates": 1,
        "analytic_all_length_large_prime_families": 1,
        "actual_ray_positive_propagated_families": 2,
        "open_small_prime_length_families": 2,
        "proved_joined_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    rebuilt = rebuild_m9_certificate()
    rebuilt_json = json.loads(json.dumps(rebuilt))
    certificate = stored.get("exact", {}).get(
        "dyadic_m9", {}
    ).get("certificate", {})
    if certificate != rebuilt_json:
        issues.append("stored dyadic M=9 certificate drifted")
    if (
        not rebuilt.get("positive")
        or rebuilt.get("claimed_rational_lower") != "1/20"
        or len(rebuilt.get("x_leaf_paths", ())) != 32
    ):
        issues.append("dyadic M=9 exact certificate failed")

    exact = stored.get("exact", {})
    if "Q_(M+1)(z)=Q_M(q^2 z)" not in exact.get(
        "append_one", {}
    ).get("polynomial", ""):
        issues.append("stored append-one identity drifted")
    if "p>=5,M>=2: J_0>1/20" not in exact.get(
        "large_prime_fejer", {}
    ).get("theorem", ""):
        issues.append("stored large-prime theorem drifted")
    if "J_ray>49/1000" not in exact.get(
        "uniform_actual_transfer", {}
    ).get("margins", ""):
        issues.append("stored actual margin drifted")

    boundary = stored.get("proof_boundary", "")
    for required in (
        "p=2,M>=10",
        "p=3,M>=5",
        "p-free",
        "Xi Abel gap",
        "successor winding",
        "Lambda<=0",
        "RH",
    ):
        if required not in boundary:
            issues.append(f"proof boundary missing: {required}")

    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
    else:
        note = NOTE.read_text(encoding="utf-8")
        for required in (
            "# Prime-Power Heat-Starlikeness Length Propagation",
            "## Append-One Identity",
            "## Dyadic M=9",
            "## All Large-Prime Lengths",
            "## Actual Fixed-Ray Transfer",
            "1 analytic all-length p>=5 family",
            "2 open small-prime length families",
        ):
            if required not in note:
                issues.append(f"note marker missing: {required}")

    validate_append_identity(issues)
    validate_fejer_identity(issues)
    validate_large_prime_algebra(issues)
    validate_actual_transfer(issues)
    validate_sources(stored, sources, issues)
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"ERROR: {issue}")
    if issues:
        return 1
    print(
        "validated prime-power heat-starlikeness length propagation gate: "
        "20 rows, 1 append-one identity, "
        "1 exact dyadic next-length certificate, "
        "1 analytic all-length p>=5 family, "
        "2 actual ray-positive propagated families, "
        "2 open small-prime length families, "
        "0 joined Abel gaps, 0 successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
