#!/usr/bin/env python3
"""Check the exact joined phase-current polarization gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "joined_phase_current_polarization_gate"
)
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "ray_bottom": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "ray_bottom_logarithmic_flow_reduction.json"
    ),
    "phase_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_logarithmic_phase_flow_gate.json"
    ),
    "large_prime": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_length_propagation_gate.json"
    ),
    "ternary": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_ternary_length_gate.json"
    ),
    "dyadic": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_dyadic_length_gate.json"
    ),
    "short": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_starlikeness_short_family_counter_gate.json"
    ),
}
EXPECTED_ROWS = [
    "jpcp_01_polarization",
    "jpcp_02_contact",
    "jpcp_03_chain",
    "jpcp_04_long",
    "jpcp_05_cubic",
    "jpcp_06_negative",
    "jpcp_07_zero_fibre",
    "jpcp_08_route",
    "jpcp_09_pi",
    "jpcp_10_boundary",
]
ComplexQ = tuple[Fraction, Fraction]


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def parse_complex(payload: dict) -> ComplexQ:
    return parse_fraction(payload["real"]), parse_fraction(payload["imag"])


def add(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return left[0] + right[0], left[1] + right[1]


def conjugate(value: ComplexQ) -> ComplexQ:
    return value[0], -value[1]


def multiply(left: ComplexQ, right: ComplexQ) -> ComplexQ:
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


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
    if payload.get("kind") != "joined_phase_current_polarization_gate":
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
        "exact_polarization_witnesses": 1,
        "exact_chain_insertion_witnesses": 1,
        "positive_self_current_cubic_contacts": 1,
        "negative_joined_current_witnesses": 1,
        "long_family_diagonal_inputs": 3,
        "proved_joined_cross_current_bounds": 0,
        "proved_joined_abel_gaps": 0,
        "proved_successor_winding_bounds": 0,
    }
    summary = payload.get("summary", {})
    for key, value in expected.items():
        if summary.get(key) != value:
            issues.append(f"summary drifted at {key}")


def validate_polarization(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("polarization", {})
    values = [parse_complex(row) for row in witness.get("values", [])]
    derivatives = [
        parse_complex(row) for row in witness.get("derivatives", [])
    ]
    if len(values) != 3 or len(derivatives) != 3:
        issues.append("polarization witness size drifted")
        return
    total_value = (Fraction(0), Fraction(0))
    total_derivative = (Fraction(0), Fraction(0))
    diagonal = Fraction(0)
    for value, derivative in zip(values, derivatives, strict=True):
        total_value = add(total_value, value)
        total_derivative = add(total_derivative, derivative)
        diagonal += current(derivative, value)
    cross = Fraction(0)
    for left in range(3):
        for right in range(left + 1, 3):
            cross += current(derivatives[left], values[right])
            cross += current(derivatives[right], values[left])
    joined = current(total_derivative, total_value)
    if joined != diagonal + cross:
        issues.append("exact polarization reconstruction failed")
    if parse_complex(witness["total_value"]) != total_value:
        issues.append("stored total value drifted")
    if parse_complex(witness["total_derivative"]) != total_derivative:
        issues.append("stored total derivative drifted")
    if parse_fraction(witness["joined_current"]) != joined:
        issues.append("stored joined current drifted")


def validate_chain_insertion(payload: dict, issues: list[str]) -> None:
    witness = payload.get("witnesses", {}).get("chain_insertion", {})
    q_value = parse_complex(witness["Q"])
    q_norm = q_value[0] ** 2 + q_value[1] ** 2
    j_ray = parse_fraction(witness["J_ray"])
    tau = parse_fraction(witness["tau"])
    projector = parse_fraction(witness["projector_rate"])
    external = parse_fraction(witness["external_rate"])
    inserted = tau * j_ray + (projector + external - tau) * q_norm
    direct = parse_fraction(witness["direct_self_current"])
    if q_norm != parse_fraction(witness["Q_norm_squared"]):
        issues.append("chain Q norm drifted")
    if inserted != direct:
        issues.append("chain insertion identity failed")
    if inserted != parse_fraction(witness["inserted_self_current"]):
        issues.append("stored inserted current drifted")


def validate_guards(payload: dict, issues: list[str]) -> None:
    guards = payload.get("guards", {})
    cubic = guards.get("cubic_contact", {})
    negative = guards.get("negative_current", {})
    for name, guard, amplitude in (
        ("cubic", cubic, Fraction(1, 3)),
        ("negative", negative, Fraction(1, 2)),
    ):
        if parse_fraction(guard.get("amplitude", "0")) != amplitude:
            issues.append(f"{name} amplitude drifted")
            continue
        first = parse_complex(guard["first_value"])
        second = parse_complex(guard["second_value"])
        joined = add(first, second)
        joined_derivative = parse_complex(guard["joined_derivative"])
        exact = current(joined_derivative, joined)
        formula = (1 - amplitude) * (1 - 3 * amplitude)
        if parse_complex(guard["joined_value"]) != joined:
            issues.append(f"{name} joined value drifted")
        if exact != formula:
            issues.append(f"{name} factorized current failed")
        if parse_fraction(guard["joined_current"]) != exact:
            issues.append(f"{name} stored current drifted")
        if parse_fraction(
            guard["real_projection_linear_coefficient"]
        ) != 1 - 3 * amplitude:
            issues.append(f"{name} projection linear coefficient drifted")
        if parse_fraction(
            guard["real_projection_cubic_coefficient"]
        ) != 4 * amplitude:
            issues.append(f"{name} projection cubic coefficient drifted")
        if parse_fraction(guard["first_self_current"]) <= 0:
            issues.append(f"{name} first self current is not positive")
        if parse_fraction(guard["second_self_current"]) <= 0:
            issues.append(f"{name} second self current is not positive")
    if parse_complex(cubic["joined_value"]) != (Fraction(0), Fraction(2, 3)):
        issues.append("cubic joined value drifted")
    if parse_complex(cubic["joined_derivative"]) != (
        Fraction(0),
        Fraction(0),
    ):
        issues.append("cubic first jet does not vanish")
    if parse_fraction(cubic["joined_current"]) != 0:
        issues.append("cubic joined current is not zero")
    if parse_fraction(
        cubic["real_projection_linear_coefficient"]
    ) != 0:
        issues.append("cubic projection retained a linear term")
    if parse_fraction(
        cubic["real_projection_cubic_coefficient"]
    ) != Fraction(4, 3):
        issues.append("cubic projection coefficient is not 4/3")
    if parse_fraction(negative["joined_current"]) != Fraction(-1, 4):
        issues.append("negative joined current is not -1/4")


def validate_note(payload: dict, issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append("missing note")
        return
    text = NOTE.read_text(encoding="utf-8")
    required = (
        "# Joined Phase-Current Polarization Gate",
        "K(W)=Im((D_j W) conjugate(W))",
        "cos(theta)+(1/3)cos(3theta)=(4/3)cos(theta)^3",
        "p=2,M>=8",
        "p=3,M>=4",
        "p>=5,M>=2",
        "division-free Abel scalar",
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
        validate_polarization(payload, issues)
        validate_chain_insertion(payload, issues)
        validate_guards(payload, issues)
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
        "validated joined phase-current polarization gate: "
        f"{summary['rows']} rows, "
        f"{summary['exact_polarization_witnesses']} exact polarization, "
        f"{summary['exact_chain_insertion_witnesses']} exact chain insertion, "
        f"{summary['positive_self_current_cubic_contacts']} cubic contact, "
        f"{summary['negative_joined_current_witnesses']} negative joined current, "
        f"{summary['long_family_diagonal_inputs']} long-family inputs, "
        f"{summary['proved_joined_cross_current_bounds']} joined cross-current bounds, "
        f"{summary['proved_joined_abel_gaps']} joined Abel gaps, "
        f"{summary['proved_successor_winding_bounds']} successor winding bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
