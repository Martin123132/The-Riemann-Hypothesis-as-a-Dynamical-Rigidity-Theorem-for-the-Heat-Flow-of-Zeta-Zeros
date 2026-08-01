#!/usr/bin/env python3
"""Check the centered Mangoldt adjacent-cutoff transport gate."""

from __future__ import annotations

import argparse
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
    "mangoldt_adjacent_cutoff_transport_gate"
)
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
SOURCE_PATHS = {
    "contact_centering": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_contact_centering_gate.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
    "chart_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_chart_stability_certificate.json"
    ),
}

ComplexQ = tuple[Fraction, Fraction]
LogVector = dict[int, ComplexQ]


def parse_fraction(value: str) -> Fraction:
    numerator, denominator = value.split("/", maxsplit=1)
    return Fraction(int(numerator), int(denominator))


def parse_complex(value: dict[str, str]) -> ComplexQ:
    return parse_fraction(value["real"]), parse_fraction(value["imag"])


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def negate(value: ComplexQ) -> ComplexQ:
    return -value[0], -value[1]


def scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def vector_add(left: LogVector, right: LogVector) -> LogVector:
    result = dict(left)
    for prime, value in right.items():
        result[prime] = add(
            result.get(prime, (Fraction(0), Fraction(0))),
            value,
        )
        if result[prime] == (Fraction(0), Fraction(0)):
            del result[prime]
    return result


def vector_subtract(left: LogVector, right: LogVector) -> LogVector:
    return vector_add(
        left,
        {prime: negate(value) for prime, value in right.items()},
    )


def validate_sources(payload: dict, issues: list[str]) -> None:
    stored = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if stored.get(key) != actual:
            issues.append(f"source hash drifted: {key}")


def validate_symbolic(issues: list[str]) -> None:
    delta_e, delta_g, z_n = sp.symbols("delta_e delta_g z_n")
    s_prime = sp.symbols("s_prime")
    log_a, log_n = sp.symbols("log_a log_n")

    delta_endpoint = delta_g - s_prime * log_a * delta_e
    delta_mangoldt = log_n * z_n
    delta_value = z_n + delta_e
    centered = sp.expand(
        delta_endpoint
        - s_prime * delta_mangoldt
        + s_prime * log_a * delta_value
    )
    reduced = delta_g + s_prime * (log_a - log_n) * z_n
    if sp.expand(centered - reduced) != 0:
        issues.append("symbolic centered endpoint cancellation failed")

    ell, f_n, modulus = sp.symbols(
        "ell f_n modulus",
        nonzero=True,
    )
    recurrence_slope = sp.symbols("recurrence_slope")
    physical = (
        recurrence_slope - s_prime * ell * f_n
    ) / modulus
    substituted = reduced.subs(
        {
            delta_g: recurrence_slope / modulus,
            z_n: f_n / modulus,
            log_a - log_n: -ell,
        }
    )
    if sp.simplify(substituted - physical) != 0:
        issues.append("symbolic physical recurrence match failed")


def validate_rational_jump(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("rational_jump", {})
    try:
        kappa = parse_complex(witness["kappa_N"])
        j_value = parse_complex(witness["J_a"])
        j_x = parse_complex(witness["J_a_x"])
        mu = parse_complex(witness["mu_a"])
        f_n = parse_complex(witness["f_n"])
        modulus = parse_fraction(witness["abs_f_1"])
        s_prime = parse_complex(witness["s_prime"])
        log_a = parse_fraction(witness["log_a"])
        log_n = parse_fraction(witness["log_n"])
        ell = parse_fraction(witness["ell"])
    except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
        issues.append(f"rational jump witness malformed: {exc}")
        return

    if modulus <= 0:
        issues.append("rational jump anchor modulus is not positive")
        return
    if ell != log_n - log_a:
        issues.append("rational jump logarithmic displacement drifted")

    inverse_modulus = Fraction(1, 1) / modulus
    recurrent_bracket = add(j_x, multiply(mu, j_value))
    delta_e = scale(multiply(kappa, j_value), inverse_modulus)
    delta_g = scale(
        multiply(kappa, recurrent_bracket),
        inverse_modulus,
    )
    z_n = scale(f_n, inverse_modulus)
    delta_z_lambda = scale(z_n, log_n)
    delta_w_0 = add(z_n, delta_e)
    delta_endpoint = add(
        delta_g,
        negate(multiply(s_prime, scale(delta_e, log_a))),
    )
    centered = add(
        delta_endpoint,
        add(
            negate(multiply(s_prime, delta_z_lambda)),
            multiply(s_prime, scale(delta_w_0, log_a)),
        ),
    )
    reduced = add(
        delta_g,
        multiply(s_prime, scale(z_n, log_a - log_n)),
    )
    physical = scale(
        add(
            negate(multiply(s_prime, scale(f_n, ell))),
            multiply(kappa, recurrent_bracket),
        ),
        inverse_modulus,
    )

    if centered != reduced:
        issues.append("rational centered cancellation failed")
    if reduced != physical:
        issues.append("rational physical recurrence match failed")

    expected = {
        "Delta_e": delta_e,
        "Delta_g": delta_g,
        "z_n": z_n,
        "Delta_Z_Lambda": delta_z_lambda,
        "Delta_W_0": delta_w_0,
        "Delta_E_N": delta_endpoint,
        "Delta_mathfrak_B_direct": centered,
        "Delta_mathfrak_B_reduced": reduced,
        "Delta_mathfrak_B_physical": physical,
    }
    for key, value in expected.items():
        try:
            stored = parse_complex(witness[key])
        except (KeyError, TypeError, ValueError, ZeroDivisionError) as exc:
            issues.append(f"stored rational jump value malformed: {key}: {exc}")
            continue
        if stored != value:
            issues.append(f"stored rational jump value drifted: {key}")


def prime_factorization(number: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    remainder = number
    candidate = 2
    while candidate * candidate <= remainder:
        while remainder % candidate == 0:
            factors[candidate] = factors.get(candidate, 0) + 1
            remainder //= candidate
        candidate = 3 if candidate == 2 else candidate + 2
    if remainder > 1:
        factors[remainder] = factors.get(remainder, 0) + 1
    return factors


def lambda_vector(number: int) -> dict[int, int]:
    factors = prime_factorization(number)
    if number <= 1 or len(factors) != 1:
        return {}
    return {next(iter(factors)): 1}


def deterministic_carrier(number: int) -> ComplexQ:
    return (
        Fraction((7 * number) % 19 - 9, number + 5),
        Fraction((11 * number) % 23 - 11, number + 7),
    )


def add_log_coefficients(
    target: LogVector,
    logarithm: dict[int, int],
    value: ComplexQ,
    scalar: Fraction = Fraction(1),
) -> None:
    for prime, exponent in logarithm.items():
        contribution = scale(value, scalar * exponent)
        target[prime] = add(
            target.get(prime, (Fraction(0), Fraction(0))),
            contribution,
        )
        if target[prime] == (Fraction(0), Fraction(0)):
            del target[prime]


def direct_moment(cutoff: int) -> LogVector:
    result: LogVector = {}
    for number in range(1, cutoff + 1):
        add_log_coefficients(
            result,
            prime_factorization(number),
            deterministic_carrier(number),
        )
    return result


def balanced_parts(cutoff: int) -> tuple[LogVector, LogVector]:
    split = math.isqrt(cutoff)
    square: LogVector = {}
    wing: LogVector = {}
    for left in range(1, split + 1):
        for right in range(1, split + 1):
            value = deterministic_carrier(left * right)
            add_log_coefficients(
                square,
                lambda_vector(left),
                value,
                Fraction(1, 2),
            )
            add_log_coefficients(
                square,
                lambda_vector(right),
                value,
                Fraction(1, 2),
            )
        for right in range(split + 1, cutoff // left + 1):
            value = deterministic_carrier(left * right)
            add_log_coefficients(
                wing,
                lambda_vector(left),
                value,
            )
            add_log_coefficients(
                wing,
                lambda_vector(right),
                value,
            )
    return square, wing


def direct_boundary(number: int) -> LogVector:
    result: LogVector = {}
    add_log_coefficients(
        result,
        prime_factorization(number),
        deterministic_carrier(number),
    )
    return result


def square_components(
    number: int,
) -> tuple[LogVector, LogVector, LogVector]:
    root = math.isqrt(number)
    if root * root != number:
        raise ValueError("square components require a perfect square")

    moved: LogVector = {}
    diagonal: LogVector = {}
    boundary: LogVector = {}
    for divisor in range(1, root):
        value = deterministic_carrier(divisor * root)
        add_log_coefficients(moved, lambda_vector(divisor), value)
        add_log_coefficients(moved, lambda_vector(root), value)

    add_log_coefficients(
        diagonal,
        lambda_vector(root),
        deterministic_carrier(number),
    )
    for divisor in range(1, root):
        if number % divisor:
            continue
        other = number // divisor
        value = deterministic_carrier(number)
        add_log_coefficients(boundary, lambda_vector(divisor), value)
        add_log_coefficients(boundary, lambda_vector(other), value)
    return moved, diagonal, boundary


def validate_transitions(payload: dict, issues: list[str]) -> None:
    rows = payload.get("witnesses", {}).get("cutoff_transitions", [])
    expected_cutoffs = list(range(2, 81))
    if [row.get("N") for row in rows] != expected_cutoffs:
        issues.append("cutoff transition registry drifted")
        return

    ordinary = 0
    squares = 0
    mismatches = 0
    for row in rows:
        cutoff = int(row["N"])
        next_cutoff = cutoff + 1
        old_square, old_wing = balanced_parts(cutoff)
        new_square, new_wing = balanced_parts(next_cutoff)
        old_total = vector_add(old_square, old_wing)
        new_total = vector_add(new_square, new_wing)
        if old_total != direct_moment(cutoff):
            issues.append(f"balanced level failed at N={cutoff}")
            mismatches += 1
        if new_total != direct_moment(next_cutoff):
            issues.append(f"balanced level failed at N={next_cutoff}")
            mismatches += 1

        delta_square = vector_subtract(new_square, old_square)
        delta_wing = vector_subtract(new_wing, old_wing)
        total_jump = vector_add(delta_square, delta_wing)
        expected = direct_boundary(next_cutoff)
        is_square = math.isqrt(next_cutoff) ** 2 == next_cutoff

        if total_jump != expected:
            issues.append(f"total transition failed at N={cutoff}")
            mismatches += 1

        if is_square:
            squares += 1
            moved, diagonal, boundary = square_components(next_cutoff)
            if delta_square != vector_add(moved, diagonal):
                issues.append(
                    f"perfect-square central transition failed at N={cutoff}"
                )
                mismatches += 1
            if delta_wing != vector_subtract(boundary, moved):
                issues.append(
                    f"perfect-square wing transition failed at N={cutoff}"
                )
                mismatches += 1
            if vector_add(diagonal, boundary) != expected:
                issues.append(
                    f"perfect-square recombination failed at N={cutoff}"
                )
                mismatches += 1
            transition = "perfect_square"
            transfer_coordinates = len(moved)
        else:
            ordinary += 1
            if delta_square:
                issues.append(
                    f"ordinary square component changed at N={cutoff}"
                )
                mismatches += 1
            if delta_wing != expected:
                issues.append(
                    f"ordinary wing insertion failed at N={cutoff}"
                )
                mismatches += 1
            transition = "ordinary"
            transfer_coordinates = 0

        expected_metadata = {
            "n": next_cutoff,
            "transition": transition,
            "D_before": math.isqrt(cutoff),
            "D_after": math.isqrt(next_cutoff),
            "delta_square_coordinates": len(delta_square),
            "delta_wing_coordinates": len(delta_wing),
            "row_transfer_coordinates": transfer_coordinates,
            "boundary_coordinates": len(expected),
            "mismatches": 0,
        }
        for key, value in expected_metadata.items():
            if row.get(key) != value:
                issues.append(
                    f"stored transition metadata drifted at N={cutoff}: {key}"
                )

    if ordinary != 71 or squares != 8:
        issues.append("independent transition class counts drifted")
    if mismatches != 0:
        issues.append(
            f"independent transition audit has {mismatches} mismatches"
        )


def validate_exact_statements(payload: dict, issues: list[str]) -> None:
    exact = payload.get("exact", {})
    required = {
        "anchored_component_jumps": (
            "Delta W_0=z_n+Delta e",
            "Delta Z_Lambda=log(n)z_n",
        ),
        "centered_jet_jump": (
            "endpoint-value jump cancels exactly",
            "Delta g+s_*'log(a/n)z_n",
        ),
        "physical_recurrence_match": (
            "Delta W_A",
            "Delta A_a/|f_1|",
        ),
        "real_projection_bound": (
            "|Re(Delta W_0)|<2200exp(-5L/4)",
            "|Re(Delta mathfrak B_N)|<10000exp(-7L/4)",
        ),
        "complex_boundary": (
            "No complex absolute-value bound",
            "real projection",
        ),
        "ordinary_split_jump": (
            "Delta Z_square=0",
            "Delta Z_wing=log(n)z_n",
        ),
        "square_split_jump": (
            "Delta Z_square=R_r+D_r",
            "Delta Z_wing=A_n-R_r",
        ),
        "square_recombination": (
            "D_r+A_n=log(n)z_n",
            "no residual term",
        ),
    }
    for key, markers in required.items():
        text = str(exact.get(key, ""))
        for marker in markers:
            if marker not in text:
                issues.append(f"exact statement marker missing: {key}: {marker}")


def validate_rows(payload: dict, issues: list[str]) -> None:
    rows = payload.get("rows", [])
    suffixes = (
        "components",
        "endpoint",
        "cancellation",
        "recurrence",
        "projection",
        "complex_guard",
        "ordinary",
        "square",
        "square_rejoin",
        "transition_audit",
        "band",
        "route",
    )
    expected_ids = [
        f"mactg_{index:02d}_{suffix}"
        for index, suffix in enumerate(suffixes, start=1)
    ]
    if [row.get("id") for row in rows] != expected_ids:
        issues.append("row registry drifted")

    for row in rows:
        for field in (
            "claim",
            "readiness",
            "exact_statement",
            "implication",
            "boundary",
        ):
            if not row.get(field):
                issues.append(f"empty row field {row.get('id')}: {field}")

    readiness = {row.get("id"): row.get("readiness") for row in rows}
    if readiness.get("mactg_06_complex_guard") != "guard_validated":
        issues.append("complex guard readiness drifted")
    if readiness.get("mactg_10_transition_audit") != "finite_certificate":
        issues.append("finite transition audit readiness drifted")
    if readiness.get("mactg_12_route") != "open":
        issues.append("route row must remain open")
    for row_id in expected_ids:
        if row_id in {
            "mactg_06_complex_guard",
            "mactg_10_transition_audit",
            "mactg_12_route",
        }:
            continue
        if readiness.get(row_id) != "proved":
            issues.append(f"proved row readiness drifted: {row_id}")


def validate_summary(payload: dict, issues: list[str]) -> None:
    expected = {
        "rows": 12,
        "centered_endpoint_cancellations": 1,
        "cutoff_transitions": 79,
        "ordinary_transitions": 71,
        "perfect_square_transitions": 8,
        "transition_mismatches": 0,
        "proved_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    summary = payload.get("summary", {})
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(f"summary drifted: {key}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", nargs="?", type=Path, default=DEFAULT_RESULT)
    args = parser.parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    issues: list[str] = []

    if payload.get("kind") != "mangoldt_adjacent_cutoff_transport_gate":
        issues.append("kind drifted")
    if payload.get("date") != "2026-07-30":
        issues.append("date drifted")
    status = str(payload.get("status", "")).lower()
    boundary = str(payload.get("proof_boundary", "")).lower()
    if "signed lower bound remains open" not in status:
        issues.append("status boundary missing: signed lower bound remains open")
    for token in (
        "no signed lower bound",
        "abel gap",
        "lambda<=0",
        "rh",
        "prize-level",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing: {token}")

    validate_sources(payload, issues)
    validate_symbolic(issues)
    validate_rational_jump(payload, issues)
    validate_transitions(payload, issues)
    validate_exact_statements(payload, issues)
    validate_rows(payload, issues)
    validate_summary(payload, issues)

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1

    summary = payload["summary"]
    print(
        "validated centered Mangoldt cutoff transport gate: "
        f"{summary['rows']} rows, "
        f"{summary['centered_endpoint_cancellations']} endpoint cancellation, "
        f"{summary['cutoff_transitions']} transitions, "
        f"{summary['ordinary_transitions']} ordinary, "
        f"{summary['perfect_square_transitions']} square, "
        f"{summary['transition_mismatches']} mismatches, "
        f"{summary['proved_abel_gaps']} Abel gaps, "
        f"{summary['proved_successor_winding_bounds']} winding bounds"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
