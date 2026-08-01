#!/usr/bin/env python3
"""Check the corrected four-moment Mangoldt-Abel contact gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "mangoldt_abel_contact_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "first_order_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_signed_contact_reduction.json"
    ),
    "absolute_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "legacy_symmetry": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "legacy_prime_curvature_symmetry_audit.json"
    ),
}
EXPECTED_ROWS = [
    "macg_01_correction",
    "macg_02_anchor",
    "macg_03_four_moments",
    "macg_04_mangoldt",
    "macg_05_heat_factor",
    "macg_06_heat_gram",
    "macg_07_bilinear",
    "macg_08_contact",
    "macg_09_band",
    "macg_10_minor",
    "macg_11_legacy",
    "macg_12_route",
]
ComplexQ = tuple[Fraction, Fraction]


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def parse_complex(payload: dict) -> ComplexQ:
    return parse_fraction(payload["real"]), parse_fraction(payload["imag"])


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


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
    if payload.get("kind") != "corrected_mangoldt_abel_contact_gate":
        issues.append("kind drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_ROWS:
        issues.append("row ids or order drifted")
    expected_readiness = [
        "proved",
        "proved",
        "proved",
        "proved",
        "proved",
        "proved",
        "proved",
        "proved",
        "conditional",
        "proved",
        "proved",
        "open",
    ]
    if [row.get("readiness") for row in rows] != expected_readiness:
        issues.append("row readiness drifted")
    expected_summary = {
        "rows": 12,
        "exact_correction_polynomials": 1,
        "logarithmic_moments": 4,
        "symmetric_mangoldt_moments": 3,
        "prime_edge_witnesses": 8,
        "indefinite_prime_edge_minors": 24,
        "proved_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    summary = payload.get("summary", {})
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            issues.append(f"summary drifted at {key}")


def validate_symbolic(issues: list[str]) -> None:
    ell, alpha, alpha_prime, heat_time, s = sp.symbols(
        "ell alpha alpha_prime heat_time s"
    )
    quadratic = alpha_prime * heat_time**2 / 8
    a_0 = (
        1
        + 1 / (6 * s)
        + alpha_prime * heat_time / 4
        + quadratic * alpha**2
    )
    a_1 = -2 * quadratic * alpha
    a_2 = quadratic
    direct = (
        1
        + 1 / (6 * s)
        + alpha_prime
        * (
            heat_time / 4
            + heat_time**2 * (alpha - ell) ** 2 / 8
        )
    )
    if sp.expand(direct - a_0 - a_1 * ell - a_2 * ell**2) != 0:
        issues.append("correction polynomial failed")

    varrho_1, varrho_2, log_a = sp.symbols(
        "varrho_1 varrho_2 log_a"
    )
    direct_abel = sp.expand(
        (log_a - ell) * (1 + varrho_1 * ell + varrho_2 * ell**2)
    )
    collapsed_abel = (
        log_a
        + (varrho_1 * log_a - 1) * ell
        + (varrho_2 * log_a - varrho_1) * ell**2
        - varrho_2 * ell**3
    )
    if sp.expand(direct_abel - collapsed_abel) != 0:
        issues.append("four-moment polynomial failed")

    log_d, log_m, s_star = sp.symbols("log_d log_m s_star")
    left = sp.exp(
        heat_time * (log_d + log_m) ** 2 / 4
        - s_star * (log_d + log_m)
    )
    right = (
        sp.exp(heat_time * log_d**2 / 4 - s_star * log_d)
        * sp.exp(heat_time * log_m**2 / 4 - s_star * log_m)
        * sp.exp(heat_time * log_d * log_m / 2)
    )
    if sp.simplify(left - right) != 0:
        issues.append("heat factorization failed")


def factor_vector(number: int) -> dict[int, Fraction]:
    return {
        int(prime): Fraction(int(exponent))
        for prime, exponent in sp.factorint(number).items()
    }


def lambda_vector(number: int) -> dict[int, Fraction]:
    factors = sp.factorint(number)
    if number <= 1 or len(factors) != 1:
        return {}
    return {int(next(iter(factors))): Fraction(1)}


def vector_add(
    target: dict[int, Fraction],
    source: dict[int, Fraction],
    scalar: Fraction,
) -> None:
    for prime, value in source.items():
        target[prime] = target.get(prime, Fraction(0)) + scalar * value
        if target[prime] == 0:
            del target[prime]


def validate_mangoldt(payload: dict, issues: list[str]) -> None:
    audit = payload.get("witnesses", {}).get("mangoldt", {})
    limit = int(audit.get("limit", 0))
    if limit != 128:
        issues.append("Mangoldt audit limit drifted")
        return
    terms_seen = 0
    for number in range(1, limit + 1):
        one_sided: dict[int, Fraction] = {}
        symmetric: dict[int, Fraction] = {}
        for divisor_value in sp.divisors(number):
            divisor = int(divisor_value)
            other = number // divisor
            left = lambda_vector(divisor)
            right = lambda_vector(other)
            if left:
                terms_seen += 1
            vector_add(one_sided, left, Fraction(1))
            vector_add(symmetric, left, Fraction(1, 2))
            vector_add(symmetric, right, Fraction(1, 2))
        expected = factor_vector(number)
        if one_sided != expected:
            issues.append(f"one-sided Mangoldt identity failed at n={number}")
        if symmetric != expected:
            issues.append(f"symmetric Mangoldt identity failed at n={number}")
    if audit.get("coefficient_audits") != limit:
        issues.append("Mangoldt coefficient count drifted")
    if audit.get("moment_families") != 3:
        issues.append("Mangoldt moment-family count drifted")
    if audit.get("prime_power_terms_seen") != terms_seen:
        issues.append("Mangoldt prime-power count drifted")
    if audit.get("mismatches") != 0:
        issues.append("stored Mangoldt mismatches are nonzero")


def validate_rational_moment(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("rational_moment", {})
    logarithms = [parse_fraction(value) for value in witness["logarithms"]]
    base = [parse_complex(value) for value in witness["base"]]
    varrho_1 = parse_complex(witness["varrho_1"])
    varrho_2 = parse_complex(witness["varrho_2"])
    eta = parse_complex(witness["eta"])
    log_a = parse_fraction(witness["log_a"])

    moments: list[ComplexQ] = []
    for order in range(4):
        total = (Fraction(0), Fraction(0))
        for ell, value in zip(logarithms, base, strict=True):
            total = add(total, scale(value, ell**order))
        moments.append(total)
    if [parse_complex(value) for value in witness["moments"]] != moments:
        issues.append("stored logarithmic moments drifted")

    direct_s = (Fraction(0), Fraction(0))
    direct_u = (Fraction(0), Fraction(0))
    for ell, value in zip(logarithms, base, strict=True):
        polynomial = add(
            (Fraction(1), Fraction(0)),
            add(scale(varrho_1, ell), scale(varrho_2, ell**2)),
        )
        carrier = multiply(eta, multiply(value, polynomial))
        direct_s = add(direct_s, carrier)
        direct_u = add(direct_u, scale(carrier, log_a - ell))

    collapsed_s = multiply(
        eta,
        add(
            moments[0],
            add(
                multiply(varrho_1, moments[1]),
                multiply(varrho_2, moments[2]),
            ),
        ),
    )
    coefficient_1 = add(
        scale(varrho_1, log_a), (Fraction(-1), Fraction(0))
    )
    coefficient_2 = add(
        scale(varrho_2, log_a), scale(varrho_1, Fraction(-1))
    )
    collapsed_u = multiply(
        eta,
        add(
            scale(moments[0], log_a),
            add(
                multiply(coefficient_1, moments[1]),
                add(
                    multiply(coefficient_2, moments[2]),
                    multiply(scale(varrho_2, Fraction(-1)), moments[3]),
                ),
            ),
        ),
    )
    if direct_s != collapsed_s:
        issues.append("rational S collapse failed")
    if direct_u != collapsed_u:
        issues.append("rational U collapse failed")
    if parse_complex(witness["direct_S"]) != direct_s:
        issues.append("stored direct S drifted")
    if parse_complex(witness["collapsed_S"]) != collapsed_s:
        issues.append("stored collapsed S drifted")
    if parse_complex(witness["direct_U"]) != direct_u:
        issues.append("stored direct U drifted")
    if parse_complex(witness["collapsed_U"]) != collapsed_u:
        issues.append("stored collapsed U drifted")

    endpoint = parse_complex(witness["endpoint"])
    endpoint_slope = parse_complex(witness["endpoint_slope"])
    s_prime = parse_complex(witness["s_prime"])
    u_n = parse_fraction(witness["u_N"])
    x_value = add(endpoint, collapsed_s)[0]
    contact = add(endpoint_slope, multiply(s_prime, collapsed_u))[0]
    contact -= s_prime[0] * u_n * x_value
    if parse_fraction(witness["mathsf_X"]) != x_value:
        issues.append("stored mathsf_X drifted")
    if parse_fraction(witness["mathcal_C_N"]) != contact:
        issues.append("stored mathcal_C_N drifted")


def validate_prime_edges(payload: dict, issues: list[str]) -> None:
    rows = payload.get("witnesses", {}).get("prime_edges", [])
    expected_cutoffs = [2, 3, 4, 5, 10, 25, 64, 127]
    if [row.get("N") for row in rows] != expected_cutoffs:
        issues.append("prime-edge cutoffs drifted")
        return
    for row in rows:
        cutoff = int(row["N"])
        prime = int(row["p"])
        if not sp.isprime(prime):
            issues.append(f"prime edge is composite at N={cutoff}")
        if not (prime <= cutoff < prime * prime):
            issues.append(f"prime edge misses hyperbolic boundary at N={cutoff}")
        expected = {
            str(moment): f"-log({prime})^{2 * moment}/4"
            for moment in (1, 2, 3)
        }
        if row.get("determinants") != expected:
            issues.append(f"prime-edge determinant drifted at N={cutoff}")


def validate_note(payload: dict, issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append("missing note")
        return
    text = NOTE.read_text(encoding="utf-8")
    required = (
        "# Corrected Mangoldt-Abel Contact Gate",
        "four finite logarithmic moments",
        "Lambda(d)+Lambda(m)",
        "dm<=N",
        "complex bilinear form",
        "negative prime-edge minor",
        "Type-I/II",
        "not a proof",
        payload.get("proof_boundary", ""),
    )
    for needle in required:
        if needle not in text:
            issues.append(f"note missing required text: {needle[:60]}")


def main() -> int:
    if not RESULT.is_file():
        print(f"missing result: {RESULT}", file=sys.stderr)
        return 1
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    issues: list[str] = []
    validate_sources(payload, issues)
    validate_structure(payload, issues)
    validate_symbolic(issues)
    validate_mangoldt(payload, issues)
    validate_rational_moment(payload, issues)
    validate_prime_edges(payload, issues)
    validate_note(payload, issues)
    boundary = payload.get("proof_boundary", "")
    for marker in (
        "no signed endpoint-plus-moment lower bound",
        "Abel-scalar gap",
        "contact exclusion",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof boundary missing marker: {marker}")

    for issue in issues:
        print(f"ISSUE: {issue}")
    if issues:
        return 1
    print(
        "validated corrected Mangoldt-Abel contact gate: "
        "12 rows, 1 correction polynomial, 4 logarithmic moments, "
        "3 symmetric Mangoldt moments, 8 prime-edge witnesses, "
        "24 indefinite minors, 0 Abel gaps, 0 winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
