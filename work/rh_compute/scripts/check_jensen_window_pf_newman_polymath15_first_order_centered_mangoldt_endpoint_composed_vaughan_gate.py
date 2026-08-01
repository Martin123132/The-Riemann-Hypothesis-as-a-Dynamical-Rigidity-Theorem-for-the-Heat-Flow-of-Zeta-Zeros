#!/usr/bin/env python3
"""Check the centered Mangoldt endpoint-composed Vaughan gate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "mangoldt_endpoint_composed_vaughan_gate"
)
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
SOURCE_PATHS = {
    "contact_centering": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_contact_centering_gate.json"
    ),
    "cutoff_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_adjacent_cutoff_transport_gate.json"
    ),
    "corrected_mangoldt": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_abel_contact_gate.json"
    ),
}

ComplexQ = tuple[Fraction, Fraction]
LogCoeff = dict[int, int]
LogMoment = dict[int, ComplexQ]


def parse_fraction(value: str) -> Fraction:
    numerator, denominator = value.split("/", maxsplit=1)
    return Fraction(int(numerator), int(denominator))


def parse_complex(value: dict[str, str]) -> ComplexQ:
    return parse_fraction(value["real"]), parse_fraction(value["imag"])


def parse_moment(value: dict[str, dict[str, str]]) -> LogMoment:
    return {
        int(prime): parse_complex(coefficient)
        for prime, coefficient in value.items()
    }


def parse_real_log(value: dict[str, str]) -> dict[int, Fraction]:
    return {
        int(prime): parse_fraction(coefficient)
        for prime, coefficient in value.items()
    }


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def negate(value: ComplexQ) -> ComplexQ:
    return -value[0], -value[1]


def scale(value: ComplexQ, scalar: Fraction | int) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def moment_add(left: LogMoment, right: LogMoment) -> LogMoment:
    result = dict(left)
    for prime, value in right.items():
        result[prime] = add(
            result.get(prime, (Fraction(0), Fraction(0))),
            value,
        )
        if result[prime] == (Fraction(0), Fraction(0)):
            del result[prime]
    return result


def moment_subtract(left: LogMoment, right: LogMoment) -> LogMoment:
    return moment_add(
        left,
        {prime: negate(value) for prime, value in right.items()},
    )


def moment_multiply(value: LogMoment, scalar: ComplexQ) -> LogMoment:
    result: LogMoment = {}
    for prime, coefficient in value.items():
        product = multiply(scalar, coefficient)
        if product != (Fraction(0), Fraction(0)):
            result[prime] = product
    return result


def coeff_add(
    left: LogCoeff,
    right: LogCoeff,
    scalar: int = 1,
) -> LogCoeff:
    result = dict(left)
    for prime, exponent in right.items():
        result[prime] = result.get(prime, 0) + scalar * exponent
        if result[prime] == 0:
            del result[prime]
    return result


def prime_factorization(number: int) -> LogCoeff:
    factors: LogCoeff = {}
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


def lambda_vector(number: int) -> LogCoeff:
    factors = prime_factorization(number)
    if number <= 1 or len(factors) != 1:
        return {}
    return {next(iter(factors)): 1}


def mobius_value(number: int) -> int:
    factors = prime_factorization(number)
    if any(exponent > 1 for exponent in factors.values()):
        return 0
    return -1 if len(factors) % 2 else 1


def divisors(number: int) -> list[int]:
    result = [1]
    for prime, exponent in prime_factorization(number).items():
        current = list(result)
        power = 1
        for _ in range(exponent):
            power *= prime
            result.extend(value * power for value in current)
    return sorted(result)


def divisor_count(number: int) -> int:
    count = 1
    for exponent in prime_factorization(number).values():
        count *= exponent + 1
    return count


def scaled_coeff_add(
    target: LogCoeff,
    source: LogCoeff,
    scalar: int,
) -> None:
    for prime, exponent in source.items():
        target[prime] = target.get(prime, 0) + scalar * exponent
        if target[prime] == 0:
            del target[prime]


def direct_component_coefficients(
    number: int,
    cutoff_lambda: int,
    cutoff_mu: int,
) -> dict[str, LogCoeff]:
    low: LogCoeff = {}
    type_i_main: LogCoeff = {}
    type_i_correction: LogCoeff = {}
    type_ii: LogCoeff = {}

    for lambda_index in divisors(number):
        if lambda_index <= cutoff_lambda:
            scaled_coeff_add(
                low,
                lambda_vector(lambda_index),
                1,
            )

    for mu_index in divisors(number):
        if mu_index > cutoff_mu:
            continue
        mu_coefficient = mobius_value(mu_index)
        if not mu_coefficient:
            continue
        remaining = number // mu_index
        for log_index in divisors(remaining):
            scaled_coeff_add(
                type_i_main,
                prime_factorization(log_index),
                mu_coefficient,
            )

    for lambda_index in divisors(number):
        lambda_coefficient = lambda_vector(lambda_index)
        if (
            lambda_index > cutoff_lambda
            or not lambda_coefficient
        ):
            continue
        after_lambda = number // lambda_index
        for mu_index in divisors(after_lambda):
            if mu_index > cutoff_mu:
                continue
            mu_coefficient = mobius_value(mu_index)
            if not mu_coefficient:
                continue
            remaining = after_lambda // mu_index
            scaled_coeff_add(
                type_i_correction,
                lambda_coefficient,
                mu_coefficient * divisor_count(remaining),
            )

    for lambda_index in divisors(number):
        lambda_coefficient = lambda_vector(lambda_index)
        if (
            lambda_index <= cutoff_lambda
            or not lambda_coefficient
        ):
            continue
        after_lambda = number // lambda_index
        for mu_index in divisors(after_lambda):
            if mu_index <= cutoff_mu:
                continue
            mu_coefficient = mobius_value(mu_index)
            if not mu_coefficient:
                continue
            remaining = after_lambda // mu_index
            scaled_coeff_add(
                type_ii,
                lambda_coefficient,
                mu_coefficient * divisor_count(remaining),
            )

    type_i = coeff_add(type_i_main, type_i_correction, -1)
    rebuilt = coeff_add(coeff_add(low, type_i), type_ii)
    logarithm = prime_factorization(number)
    if rebuilt != logarithm:
        raise ArithmeticError(
            f"direct Vaughan identity failed at n={number}, "
            f"U={cutoff_lambda}, V={cutoff_mu}"
        )
    return {
        "log": logarithm,
        "low": low,
        "type_i": type_i,
        "type_ii": type_ii,
    }


def deterministic_carrier(number: int) -> ComplexQ:
    return (
        Fraction((13 * number) % 29 - 14, number + 9),
        Fraction((17 * number) % 31 - 15, number + 11),
    )


def component_cache(
    maximum: int,
    split: int,
) -> dict[str, list[LogCoeff]]:
    result = {
        component: [{} for _ in range(maximum + 1)]
        for component in ("log", "low", "type_i", "type_ii")
    }
    for number in range(1, maximum + 1):
        row = direct_component_coefficients(number, split, split)
        for component, coefficient in row.items():
            result[component][number] = coefficient
    return result


def weighted_moment(
    coefficients: list[LogCoeff],
    cutoff: int,
    carriers: dict[int, ComplexQ] | None = None,
) -> LogMoment:
    result: LogMoment = {}
    for number in range(1, cutoff + 1):
        value = (
            carriers.get(number, (Fraction(0), Fraction(0)))
            if carriers is not None
            else deterministic_carrier(number)
        )
        if value == (Fraction(0), Fraction(0)):
            continue
        for prime, exponent in coefficients[number].items():
            contribution = scale(value, exponent)
            result[prime] = add(
                result.get(prime, (Fraction(0), Fraction(0))),
                contribution,
            )
            if result[prime] == (Fraction(0), Fraction(0)):
                del result[prime]
    return result


def entry_moment(
    coefficient: LogCoeff,
    value: ComplexQ,
) -> LogMoment:
    result: LogMoment = {}
    for prime, exponent in coefficient.items():
        contribution = scale(value, exponent)
        if contribution != (Fraction(0), Fraction(0)):
            result[prime] = contribution
    return result


def integer_cube_root(number: int) -> int:
    lower = 0
    upper = number + 1
    while upper - lower > 1:
        middle = (lower + upper) // 2
        if middle**3 <= number:
            lower = middle
        else:
            upper = middle
    return lower


def validate_sources(payload: dict, issues: list[str]) -> None:
    stored = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if stored.get(key) != actual:
            issues.append(f"source hash drifted: {key}")


def validate_convolution_identity(issues: list[str]) -> None:
    for cutoff_lambda, cutoff_mu in (
        (1, 1),
        (2, 3),
        (4, 4),
        (5, 7),
        (8, 6),
    ):
        for number in range(1, 257):
            try:
                direct_component_coefficients(
                    number,
                    cutoff_lambda,
                    cutoff_mu,
                )
            except ArithmeticError as exc:
                issues.append(str(exc))
                return


def validate_transitions(payload: dict, issues: list[str]) -> None:
    rows = payload.get("witnesses", {}).get("cutoff_transitions", [])
    expected_cutoffs = list(range(2, 216))
    if [row.get("N") for row in rows] != expected_cutoffs:
        issues.append("cutoff transition registry drifted")
        return

    maximum = 216
    roots = {
        max(1, integer_cube_root(number))
        for number in range(2, maximum + 1)
    }
    caches = {
        root: component_cache(maximum, root)
        for root in sorted(roots)
    }
    components = ("low", "type_i", "type_ii")
    ordinary = 0
    cubes = 0
    squares = 0
    overlaps = 0
    mismatches = 0

    for row in rows:
        cutoff = int(row["N"])
        next_cutoff = cutoff + 1
        root_before = max(1, integer_cube_root(cutoff))
        root_after = max(1, integer_cube_root(next_cutoff))
        old_coefficients = caches[root_before]
        new_coefficients = caches[root_after]
        old_moments = {
            component: weighted_moment(
                old_coefficients[component],
                cutoff,
            )
            for component in components
        }
        new_moments = {
            component: weighted_moment(
                new_coefficients[component],
                next_cutoff,
            )
            for component in components
        }
        deltas = {
            component: moment_subtract(
                new_moments[component],
                old_moments[component],
            )
            for component in components
        }
        total_delta: LogMoment = {}
        for component in components:
            total_delta = moment_add(total_delta, deltas[component])
        expected = entry_moment(
            prime_factorization(next_cutoff),
            deterministic_carrier(next_cutoff),
        )
        if total_delta != expected:
            issues.append(f"total transition failed at N={cutoff}")
            mismatches += 1

        cube_transition = root_after != root_before
        transfers: dict[str, LogMoment] = {}
        entries: dict[str, LogMoment] = {}
        for component in components:
            transfer_coefficients = [{} for _ in range(maximum + 1)]
            if cube_transition:
                for number in range(1, cutoff + 1):
                    transfer_coefficients[number] = coeff_add(
                        new_coefficients[component][number],
                        old_coefficients[component][number],
                        -1,
                    )
            transfers[component] = weighted_moment(
                transfer_coefficients,
                cutoff,
            )
            entries[component] = entry_moment(
                new_coefficients[component][next_cutoff],
                deterministic_carrier(next_cutoff),
            )
            if deltas[component] != moment_add(
                transfers[component],
                entries[component],
            ):
                issues.append(
                    f"component transition failed at N={cutoff}: "
                    f"{component}"
                )
                mismatches += 1

        total_transfer: LogMoment = {}
        for component in components:
            total_transfer = moment_add(
                total_transfer,
                transfers[component],
            )
        if total_transfer:
            issues.append(f"cube transfer did not cancel at N={cutoff}")
            mismatches += 1

        square_transition = math.isqrt(next_cutoff) ** 2 == next_cutoff
        overlap = cube_transition and square_transition
        if cube_transition:
            cubes += 1
            transition = "perfect_cube"
        else:
            ordinary += 1
            transition = "ordinary"
        squares += int(square_transition)
        overlaps += int(overlap)

        expected_metadata = {
            "n": next_cutoff,
            "K_before": root_before,
            "K_after": root_after,
            "transition": transition,
            "square_transition": square_transition,
            "sixth_power_overlap": overlap,
            "delta_low_coordinates": len(deltas["low"]),
            "delta_type_i_coordinates": len(deltas["type_i"]),
            "delta_type_ii_coordinates": len(deltas["type_ii"]),
            "transfer_low_coordinates": len(transfers["low"]),
            "transfer_type_i_coordinates": len(
                transfers["type_i"]
            ),
            "transfer_type_ii_coordinates": len(
                transfers["type_ii"]
            ),
            "boundary_coordinates": len(expected),
            "mismatches": 0,
        }
        for key, value in expected_metadata.items():
            if row.get(key) != value:
                issues.append(
                    f"stored transition metadata drifted at "
                    f"N={cutoff}: {key}"
                )

    if (ordinary, cubes, squares, overlaps) != (209, 5, 13, 1):
        issues.append("independent transition class counts drifted")
    if mismatches:
        issues.append(
            f"independent transition audit has {mismatches} mismatches"
        )


def sum_carriers(carriers: dict[int, ComplexQ]) -> ComplexQ:
    total = (Fraction(0), Fraction(0))
    for value in carriers.values():
        total = add(total, value)
    return total


def validate_prime_power_countermodel(
    payload: dict,
    issues: list[str],
) -> None:
    witness = payload.get("witnesses", {}).get("prime_power_chain", {})
    cutoff = int(witness.get("N", 0))
    split_u = int(witness.get("U", 0))
    split_v = int(witness.get("V", 0))
    if (cutoff, split_u, split_v) != (64, 4, 4):
        issues.append("prime-power witness parameters drifted")
        return
    carriers = {
        int(number): parse_complex(value)
        for number, value in witness.get("carriers", {}).items()
    }
    expected_carriers = {
        1: (Fraction(1), Fraction(0)),
        5: (Fraction(-2), Fraction(0)),
        25: (Fraction(1), Fraction(0)),
    }
    if carriers != expected_carriers:
        issues.append("prime-power carriers drifted")

    cache = component_cache(cutoff, split_u)
    low = weighted_moment(cache["low"], cutoff, carriers)
    type_i = weighted_moment(cache["type_i"], cutoff, carriers)
    type_ii = weighted_moment(cache["type_ii"], cutoff, carriers)
    total = moment_add(moment_add(low, type_i), type_ii)
    if sum_carriers(carriers) != (Fraction(0), Fraction(0)):
        issues.append("prime-power W_0 cancellation failed")
    if low or type_i != {5: (Fraction(1), Fraction(0))}:
        issues.append("prime-power low/Type-I reconstruction failed")
    if type_ii != {5: (Fraction(-1), Fraction(0))} or total:
        issues.append("prime-power Type-II cancellation failed")

    s_prime = parse_complex(witness["s_prime"])
    projected_i = {
        prime: value[0]
        for prime, value in moment_multiply(type_i, s_prime).items()
        if value[0]
    }
    projected_ii = {
        prime: value[0]
        for prime, value in moment_multiply(type_ii, s_prime).items()
        if value[0]
    }
    penalty = {
        prime: abs(projected_i.get(prime, Fraction(0)))
        + abs(projected_ii.get(prime, Fraction(0)))
        for prime in set(projected_i) | set(projected_ii)
    }
    if penalty != {5: Fraction(1, 8)}:
        issues.append("prime-power separate absolute penalty failed")

    expected_moments = {
        "Z_low": low,
        "Z_type_i": type_i,
        "Z_type_ii": type_ii,
        "Z_Lambda": total,
    }
    for key, value in expected_moments.items():
        if parse_moment(witness.get(key, {})) != value:
            issues.append(f"stored prime-power moment drifted: {key}")
    if parse_real_log(
        witness.get("separate_absolute_penalty", {})
    ) != penalty:
        issues.append("stored prime-power penalty drifted")
    dominance = {
        int(prime): parse_fraction(value)
        for prime, value in witness.get(
            "Delta_endpoint_dominance", {}
        ).items()
    }
    if dominance != {5: Fraction(-1, 8)}:
        issues.append("stored endpoint-dominance margin drifted")


def validate_prime_edge_countermodel(
    payload: dict,
    issues: list[str],
) -> None:
    witness = payload.get("witnesses", {}).get(
        "prime_edge_endpoint", {}
    )
    cutoff = int(witness.get("N", 0))
    split_u = int(witness.get("U", 0))
    split_v = int(witness.get("V", 0))
    if (cutoff, split_u, split_v) != (11, 2, 2):
        issues.append("prime-edge witness parameters drifted")
        return
    carriers = {
        int(number): parse_complex(value)
        for number, value in witness.get("carriers", {}).items()
    }
    if carriers != {
        1: (Fraction(-1), Fraction(0)),
        11: (Fraction(1), Fraction(0)),
    }:
        issues.append("prime-edge carriers drifted")

    cache = component_cache(cutoff, split_u)
    low = weighted_moment(cache["low"], cutoff, carriers)
    type_i = weighted_moment(cache["type_i"], cutoff, carriers)
    type_ii = weighted_moment(cache["type_ii"], cutoff, carriers)
    total = moment_add(moment_add(low, type_i), type_ii)
    s_prime = parse_complex(witness["s_prime"])
    endpoint_g = moment_multiply(total, s_prime)
    mathfrak_b = moment_subtract(
        endpoint_g,
        moment_multiply(total, s_prime),
    )

    if sum_carriers(carriers) != (Fraction(0), Fraction(0)):
        issues.append("prime-edge W_0 cancellation failed")
    if low or type_i != {11: (Fraction(1), Fraction(0))}:
        issues.append("prime-edge low/Type-I reconstruction failed")
    if type_ii or total != type_i:
        issues.append("prime-edge Type-II reconstruction failed")
    if mathfrak_b:
        issues.append("prime-edge endpoint cancellation failed")
    if not (
        witness.get("prime_edge") == 11
        and witness.get("prime_edge_strict") is True
        and 11 > math.isqrt(cutoff)
    ):
        issues.append("prime-edge condition drifted")

    expected_moments = {
        "endpoint_g_log_coordinates": endpoint_g,
        "Z_low": low,
        "Z_type_i": type_i,
        "Z_type_ii": type_ii,
        "Z_Lambda": total,
        "mathfrak_B_log_coordinates": mathfrak_b,
    }
    for key, value in expected_moments.items():
        if parse_moment(witness.get(key, {})) != value:
            issues.append(f"stored prime-edge moment drifted: {key}")


def validate_exact_statements(payload: dict, issues: list[str]) -> None:
    exact = payload.get("exact", {})
    required = {
        "convolution_identity": (
            "log=Lambda*1",
            "mu*1=delta",
            "mu_1*lambda_1*1",
        ),
        "weighted_vaughan": (
            "C_0+C_I+C_II=log",
            "Z_Lambda=Z_0+Z_I+Z_II",
        ),
        "endpoint_composition": (
            "mathfrak B_N=R_(0;U,V)-s_*'Z_I-s_*'Z_II",
            "before any real projection",
        ),
        "contact_fibres": (
            "At W_0=0",
            "At mathsf_X=0",
        ),
        "dominance_criterion": (
            "|Re mathfrak B_N|>=|R_0|-|R_I|-|R_II|",
            "sufficient for the Abel gap",
        ),
        "cube_transport": (
            "sum_j tau_j(m)=0",
            "log(r^3)",
        ),
        "countermodel_consequence": (
            "Z_I=log(5)",
            "Z_II=-log(5)",
        ),
        "prime_edge_consequence": (
            "Z_I=log(11)",
            "mathfrak B_N=0",
        ),
    }
    for key, markers in required.items():
        text = str(exact.get(key, ""))
        for marker in markers:
            if marker not in text:
                issues.append(
                    f"exact statement marker missing: {key}: {marker}"
                )


def validate_rows(payload: dict, issues: list[str]) -> None:
    rows = payload.get("rows", [])
    suffixes = (
        "convolution",
        "weighted",
        "endpoint",
        "fibres",
        "criterion",
        "joint",
        "cubic",
        "ordinary",
        "cube",
        "overlap",
        "chain_guard",
        "edge_guard",
        "fit",
        "route",
    )
    expected_ids = [
        f"mecvg_{index:02d}_{suffix}"
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
    expected_readiness = {
        **{
            row_id: "proved"
            for row_id in expected_ids
        },
        "mecvg_05_criterion": "proved_sufficient",
        "mecvg_10_overlap": "finite_certificate",
        "mecvg_11_chain_guard": "guard_validated",
        "mecvg_12_edge_guard": "guard_validated",
        "mecvg_13_fit": "open",
        "mecvg_14_route": "open",
    }
    for row_id, expected in expected_readiness.items():
        if readiness.get(row_id) != expected:
            issues.append(f"row readiness drifted: {row_id}")


def validate_summary(payload: dict, issues: list[str]) -> None:
    expected = {
        "rows": 14,
        "convolution_identities": 1,
        "vaughan_components": 3,
        "cutoff_transitions": 214,
        "ordinary_transitions": 209,
        "perfect_cube_transitions": 5,
        "square_transitions": 13,
        "sixth_power_overlaps": 1,
        "transition_mismatches": 0,
        "exact_countermodels": 2,
        "proved_componentwise_dominance_theorems": 0,
        "proved_signed_lower_bounds": 0,
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

    if payload.get("kind") != "mangoldt_endpoint_composed_vaughan_gate":
        issues.append("kind drifted")
    if payload.get("date") != "2026-07-31":
        issues.append("date drifted")
    status = str(payload.get("status", "")).lower()
    boundary = str(payload.get("proof_boundary", "")).lower()
    if "signed lower bound remains open" not in status:
        issues.append("status boundary missing")
    for token in (
        "no xi-specific",
        "signed lower bound",
        "abel gap",
        "lambda<=0",
        "rh",
        "prize-level",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing: {token}")

    validate_sources(payload, issues)
    validate_convolution_identity(issues)
    validate_transitions(payload, issues)
    validate_prime_power_countermodel(payload, issues)
    validate_prime_edge_countermodel(payload, issues)
    validate_exact_statements(payload, issues)
    validate_rows(payload, issues)
    validate_summary(payload, issues)

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1

    summary = payload["summary"]
    print(
        "validated centered Mangoldt endpoint-composed Vaughan gate: "
        f"{summary['rows']} rows, "
        f"{summary['convolution_identities']} convolution identity, "
        f"{summary['cutoff_transitions']} transitions, "
        f"{summary['perfect_cube_transitions']} cube transfers, "
        f"{summary['sixth_power_overlaps']} square/cube overlap, "
        f"{summary['exact_countermodels']} exact countermodels, "
        f"{summary['proved_signed_lower_bounds']} signed lower bounds, "
        f"{summary['proved_abel_gaps']} Abel gaps, "
        f"{summary['proved_successor_winding_bounds']} winding bounds"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
