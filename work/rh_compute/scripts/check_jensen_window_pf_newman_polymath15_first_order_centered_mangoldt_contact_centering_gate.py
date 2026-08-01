#!/usr/bin/env python3
"""Check the contact-centered Mangoldt/Abel equivalence gate."""

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
    "mangoldt_contact_centering_gate"
)
DEFAULT_RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
SOURCE_PATHS = {
    "mangoldt_abel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_abel_contact_gate.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "absolute_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
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


def sum_complex(values: list[ComplexQ]) -> ComplexQ:
    total = (Fraction(0), Fraction(0))
    for value in values:
        total = add(total, value)
    return total


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
    h_0, h_1, h_2, h_3 = sp.symbols("H_0 H_1 H_2 H_3")
    varrho_1, varrho_2, eta = sp.symbols(
        "varrho_1 varrho_2 eta"
    )
    log_a = sp.symbols("log_a", real=True)
    bulk = eta * (h_0 + varrho_1 * h_1 + varrho_2 * h_2)
    abel = eta * (
        log_a * h_0
        + (varrho_1 * log_a - 1) * h_1
        + (varrho_2 * log_a - varrho_1) * h_2
        - varrho_2 * h_3
    )
    mangoldt = eta * (h_1 + varrho_1 * h_2 + varrho_2 * h_3)
    if sp.expand(abel - log_a * bulk + mangoldt) != 0:
        issues.append("symbolic contact centering failed")

    c, b, x_value, y_value, log_n = sp.symbols(
        "c b mathsf_X mathsf_Y log_N", real=True
    )
    s_prime = c + sp.I * b
    w_0 = x_value + sp.I * y_value
    u_n = log_a - log_n
    left = sp.re(sp.expand(s_prime * log_a * w_0)) - c * u_n * x_value
    right = c * log_n * x_value - b * log_a * y_value
    if sp.simplify(left - right) != 0:
        issues.append("symbolic Cartesian centering failed")


def validate_centering_witness(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("rational_centering", {})
    h_values = [parse_complex(value) for value in witness["H"]]
    h_0, h_1, h_2, h_3 = h_values
    varrho_1 = parse_complex(witness["varrho_1"])
    varrho_2 = parse_complex(witness["varrho_2"])
    eta = parse_complex(witness["eta"])
    e = parse_complex(witness["e"])
    g = parse_complex(witness["g"])
    s_prime = parse_complex(witness["s_prime"])
    log_a = parse_fraction(witness["log_a"])
    log_n = parse_fraction(witness["log_N"])
    u_n = log_a - log_n

    s_bulk = multiply(
        eta,
        add(
            h_0,
            add(
                multiply(varrho_1, h_1),
                multiply(varrho_2, h_2),
            ),
        ),
    )
    u_bulk = multiply(
        eta,
        add(
            scale(h_0, log_a),
            add(
                multiply(
                    add(
                        scale(varrho_1, log_a),
                        (Fraction(-1), Fraction(0)),
                    ),
                    h_1,
                ),
                add(
                    multiply(
                        add(
                            scale(varrho_2, log_a),
                            negate(varrho_1),
                        ),
                        h_2,
                    ),
                    multiply(negate(varrho_2), h_3),
                ),
            ),
        ),
    )
    z_lambda = multiply(
        eta,
        add(
            h_1,
            add(
                multiply(varrho_1, h_2),
                multiply(varrho_2, h_3),
            ),
        ),
    )
    if u_bulk != add(scale(s_bulk, log_a), negate(z_lambda)):
        issues.append("rational U centering failed")

    w_0 = add(e, s_bulk)
    endpoint_center = add(
        g, negate(multiply(s_prime, scale(e, log_a)))
    )
    centered_jet = add(
        endpoint_center,
        add(
            negate(multiply(s_prime, z_lambda)),
            multiply(s_prime, scale(w_0, log_a)),
        ),
    )
    direct_jet = add(g, multiply(s_prime, u_bulk))
    if centered_jet != direct_jet:
        issues.append("rational centered jet failed")
    contact = direct_jet[0] - s_prime[0] * u_n * w_0[0]
    expanded = (
        add(endpoint_center, negate(multiply(s_prime, z_lambda)))[0]
        + s_prime[0] * log_n * w_0[0]
        - s_prime[1] * log_a * w_0[1]
    )
    if contact != expanded:
        issues.append("rational contact scalar failed")

    zero_e = negate(s_bulk)
    zero_endpoint = add(
        g, negate(multiply(s_prime, scale(zero_e, log_a)))
    )
    zero_jet = add(zero_endpoint, negate(multiply(s_prime, z_lambda)))
    if zero_jet != direct_jet:
        issues.append("rational zero-fibre jet failed")

    expected = {
        "u_N": u_n,
        "S": s_bulk,
        "U": u_bulk,
        "Z_Lambda": z_lambda,
        "W_0": w_0,
        "direct_jet": direct_jet,
        "centered_jet": centered_jet,
        "mathcal_C_N": contact,
        "zero_fibre_e": zero_e,
        "zero_fibre_jet": zero_jet,
    }
    for key, value in expected.items():
        if isinstance(value, tuple):
            stored = parse_complex(witness[key])
        else:
            stored = parse_fraction(witness[key])
        if stored != value:
            issues.append(f"stored centering witness drifted: {key}")


def validate_abel_duality(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("abel_duality", {})
    z_values = [
        (Fraction(1, 3), Fraction(-1, 7)),
        (Fraction(-2, 5), Fraction(1, 11)),
        (Fraction(3, 8), Fraction(2, 13)),
        (Fraction(-1, 4), Fraction(-3, 17)),
        (Fraction(2, 9), Fraction(1, 6)),
    ]
    increments = [parse_fraction(value) for value in witness["increments"]]
    log_n = parse_fraction(witness["log_N"])
    log_a = parse_fraction(witness["log_a"])
    prefixes: list[ComplexQ] = []
    running = (Fraction(0), Fraction(0))
    for value in z_values:
        running = add(running, value)
        prefixes.append(running)
    logarithms = [
        log_n - sum(increments[index:], Fraction(0))
        for index in range(len(z_values))
    ]
    direct_z = sum_complex(
        [
            scale(value, logarithm)
            for value, logarithm in zip(
                z_values, logarithms, strict=True
            )
        ]
    )
    prefix_sum = sum_complex(
        [
            scale(prefixes[index], increments[index])
            for index in range(len(increments))
        ]
    )
    s_bulk = prefixes[-1]
    prefix_z = add(scale(s_bulk, log_n), negate(prefix_sum))
    direct_u = sum_complex(
        [
            scale(value, log_a - logarithm)
            for value, logarithm in zip(
                z_values, logarithms, strict=True
            )
        ]
    )
    abel_u = add(scale(s_bulk, log_a - log_n), prefix_sum)
    if direct_z != prefix_z or direct_u != abel_u:
        issues.append("independent Abel duality failed")
    stored = {
        "S": s_bulk,
        "Z_Lambda_direct": direct_z,
        "Z_Lambda_prefix": prefix_z,
        "U_direct": direct_u,
        "U_abel": abel_u,
    }
    for key, value in stored.items():
        if parse_complex(witness[key]) != value:
            issues.append(f"stored Abel witness drifted: {key}")
    if witness.get("terms") != 5 or witness.get("mismatches") != 0:
        issues.append("Abel witness summary drifted")


def deterministic_carrier(number: int) -> ComplexQ:
    return (
        Fraction((7 * number) % 19 - 9, number + 5),
        Fraction((11 * number) % 23 - 11, number + 7),
    )


def factor_vector(number: int) -> dict[int, int]:
    return {
        int(prime): int(exponent)
        for prime, exponent in sp.factorint(number).items()
    }


def lambda_vector(number: int) -> dict[int, int]:
    factors = factor_vector(number)
    if number <= 1 or len(factors) != 1:
        return {}
    return {next(iter(factors)): 1}


def add_log_vector(
    target: LogVector,
    source: dict[int, int],
    value: ComplexQ,
    scalar: Fraction = Fraction(1),
) -> None:
    for prime, exponent in source.items():
        contribution = scale(value, scalar * exponent)
        target[prime] = add(
            target.get(prime, (Fraction(0), Fraction(0))),
            contribution,
        )
        if target[prime] == (Fraction(0), Fraction(0)):
            del target[prime]


def direct_vector(cutoff: int) -> LogVector:
    result: LogVector = {}
    for number in range(1, cutoff + 1):
        add_log_vector(
            result,
            factor_vector(number),
            deterministic_carrier(number),
        )
    return result


def one_sided_vector(cutoff: int) -> LogVector:
    result: LogVector = {}
    for divisor in range(1, cutoff + 1):
        lam = lambda_vector(divisor)
        for other in range(1, cutoff // divisor + 1):
            add_log_vector(
                result,
                lam,
                deterministic_carrier(divisor * other),
            )
    return result


def balanced_vector(cutoff: int) -> tuple[LogVector, int, int]:
    split = math.isqrt(cutoff)
    result: LogVector = {}
    central_pairs = 0
    wing_pairs = 0
    for divisor in range(1, split + 1):
        for other in range(1, split + 1):
            value = deterministic_carrier(divisor * other)
            add_log_vector(
                result,
                lambda_vector(divisor),
                value,
                Fraction(1, 2),
            )
            add_log_vector(
                result,
                lambda_vector(other),
                value,
                Fraction(1, 2),
            )
            central_pairs += 1
        for other in range(split + 1, cutoff // divisor + 1):
            value = deterministic_carrier(divisor * other)
            add_log_vector(result, lambda_vector(divisor), value)
            add_log_vector(result, lambda_vector(other), value)
            wing_pairs += 1
    return result, central_pairs, wing_pairs


def validate_hyperbola(payload: dict, issues: list[str]) -> None:
    rows = payload.get("witnesses", {}).get("hyperbola", [])
    expected_cutoffs = (2, 3, 4, 5, 10, 17, 31, 64)
    if [row.get("N") for row in rows] != list(expected_cutoffs):
        issues.append("hyperbola cutoff registry drifted")
        return
    for row in rows:
        cutoff = int(row["N"])
        direct = direct_vector(cutoff)
        one_sided = one_sided_vector(cutoff)
        balanced, central_pairs, wing_pairs = balanced_vector(cutoff)
        if direct != one_sided or direct != balanced:
            issues.append(f"balanced hyperbola failed at N={cutoff}")
        if row.get("D") != math.isqrt(cutoff):
            issues.append(f"stored split drifted at N={cutoff}")
        if row.get("log_prime_coordinates") != len(direct):
            issues.append(f"stored coordinate count drifted at N={cutoff}")
        if row.get("central_pairs") != central_pairs:
            issues.append(f"stored central count drifted at N={cutoff}")
        if row.get("wing_pairs") != wing_pairs:
            issues.append(f"stored wing count drifted at N={cutoff}")
        if row.get("mismatches") != 0:
            issues.append(f"stored mismatch at N={cutoff}")


def validate_rows(payload: dict, issues: list[str]) -> None:
    rows = payload.get("rows", [])
    expected_ids = [f"mccg_{index:02d}_{suffix}" for index, suffix in (
        (1, "log_moment"),
        (2, "abel_centering"),
        (3, "centered_jet"),
        (4, "contact_scalar"),
        (5, "zero_fibre"),
        (6, "real_crossing"),
        (7, "abel_duality"),
        (8, "endpoint_factor"),
        (9, "hyperbola"),
        (10, "band"),
        (11, "notation"),
        (12, "route"),
    )]
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
    if readiness.get("mccg_10_band") != "conditional":
        issues.append("band row must remain conditional")
    if readiness.get("mccg_12_route") != "open":
        issues.append("route row must remain open")


def validate_summary(payload: dict, issues: list[str]) -> None:
    expected = {
        "rows": 12,
        "physical_mangoldt_moments": 1,
        "contact_fibres": 2,
        "abel_duality_witnesses": 1,
        "balanced_hyperbola_audits": 8,
        "hyperbola_mismatches": 0,
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

    if payload.get("kind") != "mangoldt_contact_centering_gate":
        issues.append("kind drifted")
    if payload.get("date") != "2026-07-30":
        issues.append("date drifted")
    status = str(payload.get("status", "")).lower()
    boundary = str(payload.get("proof_boundary", "")).lower()
    for token in ("signed lower bound remains open",):
        if token not in status:
            issues.append(f"status boundary missing: {token}")
    for token in ("no signed lower bound", "rh", "prize-level"):
        if token not in boundary:
            issues.append(f"proof boundary missing: {token}")

    validate_sources(payload, issues)
    validate_symbolic(issues)
    validate_centering_witness(payload, issues)
    validate_abel_duality(payload, issues)
    validate_hyperbola(payload, issues)
    validate_rows(payload, issues)
    validate_summary(payload, issues)

    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1

    summary = payload["summary"]
    print(
        "validated contact-centered Mangoldt gate: "
        f"{summary['rows']} rows, "
        f"{summary['physical_mangoldt_moments']} physical moment, "
        f"{summary['contact_fibres']} contact fibres, "
        f"{summary['abel_duality_witnesses']} Abel duality, "
        f"{summary['balanced_hyperbola_audits']} balanced hyperbola audits, "
        f"{summary['hyperbola_mismatches']} mismatches, "
        f"{summary['proved_abel_gaps']} Abel gaps, "
        f"{summary['proved_successor_winding_bounds']} winding bounds"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
