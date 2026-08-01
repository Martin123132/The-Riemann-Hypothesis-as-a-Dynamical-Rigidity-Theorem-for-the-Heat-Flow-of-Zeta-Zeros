#!/usr/bin/env python3
"""Check the physical pair kernel and p-free quadratic-rejoin audit."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "joined_pair_kernel_pfree_rejoin_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "polarization": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_phase_current_polarization_gate.json"
    ),
    "ray_bottom": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "ray_bottom_logarithmic_flow_reduction.json"
    ),
    "joined_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "complete_chain": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "complete_prime_power_chain_occupation_gate.json"
    ),
}
EXPECTED_ROWS = [
    "jpkr_01_rate",
    "jpkr_02_pair",
    "jpkr_03_moment",
    "jpkr_04_heat_shift",
    "jpkr_05_quadratic",
    "jpkr_06_guard",
    "jpkr_07_correction",
    "jpkr_08_route",
    "jpkr_09_pi",
    "jpkr_10_boundary",
]
ComplexQ = tuple[Fraction, Fraction]


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def parse_complex(payload: dict) -> ComplexQ:
    return parse_fraction(payload["real"]), parse_fraction(payload["imag"])


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def subtract(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] - right[0], left[1] - right[1]


def scale(value: ComplexQ, scalar: Fraction) -> ComplexQ:
    return value[0] * scalar, value[1] * scalar


def conjugate(value: ComplexQ) -> ComplexQ:
    return value[0], -value[1]


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def norm_squared(value: ComplexQ) -> Fraction:
    return value[0] ** 2 + value[1] ** 2


def current(derivative: ComplexQ, value: ComplexQ) -> Fraction:
    return multiply(derivative, conjugate(value))[1]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def validate_sources(payload: dict, issues: list[str]) -> None:
    stored = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in SOURCES.items():
        if not path.is_file():
            issues.append(f"missing source {key}")
        elif stored.get(key) != file_hash(path):
            issues.append(f"source hash drifted for {key}")


def validate_structure(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != "joined_pair_kernel_pfree_rejoin_gate":
        issues.append("kind drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_ROWS:
        issues.append("row ids or order drifted")
    if any(row.get("readiness") != "proved" for row in rows[:-1]):
        issues.append("a proved row drifted")
    if not rows or rows[-1].get("readiness") != "open":
        issues.append("boundary row is not open")
    expected = {
        "rows": 10,
        "exact_pair_kernel_witnesses": 1,
        "exact_moment_collapses": 1,
        "heat_shifted_pfree_telescopes": 1,
        "quadratic_rejoin_witnesses": 1,
        "strict_cross_cancellation_guards": 1,
        "proved_cross_current_signs": 0,
        "proved_joined_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    summary = payload.get("summary", {})
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(f"summary drifted at {key}")


def validate_physical_kernel(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("physical_kernel", {})
    parameters = witness["parameters"]
    x_value = parse_fraction(parameters["x"])
    c_value = parse_fraction(parameters["c"])
    b_value = parse_fraction(parameters["b"])
    omega = parse_fraction(parameters["Omega_eta"])
    logarithms = [
        parse_fraction(value) for value in witness["logarithmic_nodes"]
    ]
    carriers = [parse_complex(value) for value in witness["carriers"]]
    corrections = [
        parse_complex(value) for value in witness["corrections"]
    ]
    rates: list[ComplexQ] = []
    derivatives: list[ComplexQ] = []
    h_0 = (Fraction(0), Fraction(0))
    h_1 = (Fraction(0), Fraction(0))
    d_0 = (Fraction(0), Fraction(0))
    for ell, carrier, correction in zip(
        logarithms, carriers, corrections, strict=True
    ):
        rate = (
            x_value * (-c_value * ell + correction[0]),
            omega + x_value * (-b_value * ell + correction[1]),
        )
        rates.append(rate)
        derivatives.append(multiply(rate, carrier))
        h_0 = add(h_0, carrier)
        h_1 = add(h_1, scale(carrier, ell))
        d_0 = add(d_0, multiply(correction, carrier))
    if [parse_complex(value) for value in witness["rates"]] != rates:
        issues.append("stored physical rates drifted")
    if parse_complex(witness["H_0"]) != h_0:
        issues.append("H_0 drifted")
    if parse_complex(witness["H_1"]) != h_1:
        issues.append("H_1 drifted")
    if parse_complex(witness["D_0"]) != d_0:
        issues.append("D_0 drifted")

    total_derivative = (Fraction(0), Fraction(0))
    diagonal = Fraction(0)
    for rate, derivative, carrier in zip(
        rates, derivatives, carriers, strict=True
    ):
        total_derivative = add(total_derivative, derivative)
        diagonal += rate[1] * norm_squared(carrier)
    pair_sum = Fraction(0)
    for left in range(len(carriers)):
        for right in range(left + 1, len(carriers)):
            correlation = multiply(
                carriers[left], conjugate(carriers[right])
            )
            formula = (
                (rates[left][0] - rates[right][0]) * correlation[1]
                + (rates[left][1] + rates[right][1]) * correlation[0]
            )
            direct = current(derivatives[left], carriers[right])
            direct += current(derivatives[right], carriers[left])
            if formula != direct:
                issues.append(f"pair kernel failed at {left},{right}")
            pair_sum += formula
    direct_total = current(total_derivative, h_0)
    s_prime = (c_value, b_value)
    bracket = add(multiply(scale(s_prime, -1), h_1), d_0)
    moment = omega * norm_squared(h_0)
    moment += x_value * multiply(bracket, conjugate(h_0))[1]
    if not direct_total == diagonal + pair_sum == moment:
        issues.append("bulk moment collapse failed")
    if parse_fraction(witness["direct_current"]) != direct_total:
        issues.append("stored direct current drifted")
    if parse_fraction(witness["diagonal_sum"]) != diagonal:
        issues.append("stored diagonal sum drifted")
    if parse_fraction(witness["pair_sum"]) != pair_sum:
        issues.append("stored pair sum drifted")
    if parse_fraction(witness["moment_current"]) != moment:
        issues.append("stored moment current drifted")


def validate_heat_shift(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("heat_shift", {})
    rows = witness.get("coefficient_rows", [])
    if len(rows) != 9:
        issues.append("heat-shift row count drifted")
        return
    for row in rows:
        index = row["k"]
        expected_t = Fraction((index + 1) ** 2, 4)
        expected_s = Fraction(-(index + 1))
        if parse_fraction(row["left_t_h2_coefficient"]) != expected_t:
            issues.append(f"left heat coefficient failed at k={index}")
        if parse_fraction(row["right_t_h2_coefficient"]) != expected_t:
            issues.append(f"right heat coefficient failed at k={index}")
        if parse_fraction(row["left_s_h_coefficient"]) != expected_s:
            issues.append(f"left Mellin coefficient failed at k={index}")
        if parse_fraction(row["right_s_h_coefficient"]) != expected_s:
            issues.append(f"right Mellin coefficient failed at k={index}")


def validate_quadratic_rejoin(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("quadratic_rejoin", {})
    t_values = [parse_complex(value) for value in witness["T_values"]]
    t_derivatives = [
        parse_complex(value) for value in witness["T_derivatives"]
    ]
    layers = (
        subtract(t_values[0], t_values[1]),
        subtract(t_values[1], t_values[2]),
    )
    derivatives = (
        subtract(t_derivatives[0], t_derivatives[1]),
        subtract(t_derivatives[1], t_derivatives[2]),
    )
    if [parse_complex(value) for value in witness["layers"]] != list(
        layers
    ):
        issues.append("stored telescoping layers drifted")
    diagonal = current(derivatives[0], layers[0])
    diagonal += current(derivatives[1], layers[1])
    cross = current(derivatives[0], layers[1])
    cross += current(derivatives[1], layers[0])
    joined_value = add(layers[0], layers[1])
    joined_derivative = add(derivatives[0], derivatives[1])
    joined = current(joined_derivative, joined_value)
    if joined_value != t_values[0]:
        issues.append("linear value telescope failed")
    if joined_derivative != t_derivatives[0]:
        issues.append("linear derivative telescope failed")
    if joined != diagonal + cross:
        issues.append("quadratic rejoin failed")
    if diagonal != Fraction(4, 3):
        issues.append("diagonal guard is not 4/3")
    if cross != Fraction(-4, 3):
        issues.append("cross guard is not -4/3")
    if joined != 0:
        issues.append("joined guard current is not zero")


def validate_note(payload: dict, issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append("missing note")
        return
    text = NOTE.read_text(encoding="utf-8")
    required = (
        "# Joined Pair Kernel And P-Free Rejoin Gate",
        "K_(n,m)",
        "H_0,H_1,D_0",
        "L_k=T_k-T_(k+1)",
        "linear telescope",
        "does not create a",
        "division-free Abel-scalar gap",
        payload.get("proof_boundary", ""),
    )
    for needle in required:
        if needle not in text:
            issues.append(f"note missing required text: {needle[:60]}")


def validate_content_hash(payload: dict, issues: list[str]) -> None:
    stored = payload.get("content_sha256")
    copy = dict(payload)
    copy.pop("content_sha256", None)
    rebuilt = hashlib.sha256(
        canonical_json(copy).encode("utf-8")
    ).hexdigest()
    if stored != rebuilt:
        issues.append("content hash drifted")


def main() -> int:
    issues: list[str] = []
    if not RESULT.is_file():
        print(f"missing result: {RESULT}", file=sys.stderr)
        return 1
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        validate_sources(payload, issues)
        validate_structure(payload, issues)
        validate_physical_kernel(payload, issues)
        validate_heat_shift(payload, issues)
        validate_quadratic_rejoin(payload, issues)
        validate_note(payload, issues)
        validate_content_hash(payload, issues)
    except Exception as exc:
        issues.append(f"checker exception: {exc}")
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    summary = payload["summary"]
    print(
        "validated joined pair-kernel p-free rejoin gate: "
        f"{summary['rows']} rows, "
        f"{summary['exact_pair_kernel_witnesses']} exact pair kernel, "
        f"{summary['exact_moment_collapses']} moment collapse, "
        f"{summary['heat_shifted_pfree_telescopes']} heat telescope, "
        f"{summary['quadratic_rejoin_witnesses']} quadratic rejoin, "
        f"{summary['strict_cross_cancellation_guards']} cross-cancellation guard, "
        f"{summary['proved_cross_current_signs']} cross-current signs, "
        f"{summary['proved_joined_abel_gaps']} joined Abel gaps, "
        f"{summary['proved_successor_winding_bounds']} successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
