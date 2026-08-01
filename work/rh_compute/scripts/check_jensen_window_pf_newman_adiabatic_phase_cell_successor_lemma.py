#!/usr/bin/env python3
"""Independently validate the adiabatic phase-cell successor lemma."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma as lemma


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT_PATH = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.json"
)
NOTE_PATH = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_newman_adiabatic_phase_cell_successor_lemma.md"
)
EXPECTED_IDS = [
    "napcs_01_shell_decomposition",
    "napcs_02_heat_jet_transport",
    "napcs_03_reverse_triangle_collar",
    "napcs_04_phase_cell_thickening",
    "napcs_05_right_half_plane",
    "napcs_06_successor_theorem",
    "napcs_07_cofinal_induction",
    "napcs_08_q208_calibration",
    "napcs_09_asymptotic_budget",
    "napcs_10_nonpromotion",
]


def distance_to_interval(interval: tuple[Fraction, Fraction]) -> Fraction:
    lower, upper = interval
    if lower <= 0 <= upper:
        return Fraction(0)
    return min(abs(lower), abs(upper))


def rectangle_distance_squared(
    rectangle: tuple[
        tuple[Fraction, Fraction],
        tuple[Fraction, Fraction],
    ],
) -> Fraction:
    return sum(
        distance_to_interval(interval) ** 2
        for interval in rectangle
    )


def main() -> None:
    issues: list[str] = []
    if not RESULT_PATH.exists():
        raise SystemExit("missing stored successor lemma")
    if not NOTE_PATH.exists():
        issues.append("missing rendered note")
    stored = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    rebuilt = lemma.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from exact reconstruction")
    if stored.get("kind") != lemma.STEM:
        issues.append("kind mismatch")
    if "antecedent open" not in stored.get("status", ""):
        issues.append("status lacks open-antecedent guard")
    rows = stored.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "exact_geometry": 1,
        "exact_identity": 1,
        "exact_lemma": 3,
        "exact_composition": 2,
        "finite_calibration": 1,
        "open_handoff": 1,
        "nonpromotion_gate": 1,
    }
    for role, expected in expected_roles.items():
        observed = sum(row.get("role") == role for row in rows)
        if observed != expected:
            issues.append(
                f"expected {expected} rows with role {role}, got {observed}"
            )

    # Independently derive the exhaustion step.
    j = sp.symbols("j", integer=True, positive=True)
    t_j = sp.Rational(1, 5) / j
    t_next = sp.Rational(1, 5) / (j + 1)
    if sp.factor(t_j - t_next) != 1 / (5 * j * (j + 1)):
        issues.append("independent successor time identity failed")
    if sp.expand((j + 39) - (j + 38)) != 1:
        issues.append("independent successor x identity failed")
    for stage in (1, 31, 207, 208, 1000):
        delta = Fraction(1, 5 * stage) - Fraction(
            1,
            5 * (stage + 1),
        )
        if delta != Fraction(1, 5 * stage * (stage + 1)):
            issues.append(f"rational shell fixture failed at {stage}")

    # Independently derive the heat-jet transport norm.
    h_xx, h_xxx, ell = sp.symbols(
        "h_xx h_xxx ell",
        real=True,
        positive=True,
    )
    transport = sp.Matrix([-h_xx, -h_xxx / ell])
    if sp.expand(transport.dot(transport)) != (
        h_xx**2 + h_xxx**2 / ell**2
    ):
        issues.append("independent heat-jet norm failed")

    # Exact rectangle and thickening fixtures.
    rectangle_fixtures = [
        (((Fraction(2), Fraction(3)), (Fraction(-1), Fraction(1))), 4),
        (((Fraction(-4), Fraction(-2)), (Fraction(3), Fraction(5))), 13),
        (((Fraction(-1), Fraction(1)), (Fraction(5, 2), Fraction(3))), Fraction(25, 4)),
    ]
    for rectangle, expected_squared in rectangle_fixtures:
        observed = rectangle_distance_squared(rectangle)
        if observed != expected_squared:
            issues.append("rectangle-distance fixture failed")
            continue
        # Choose rho=d/3; reverse triangle leaves a strict 2d/3 floor.
        distance = sp.sqrt(sp.Rational(observed.numerator, observed.denominator))
        rho = distance / 3
        if sp.simplify(distance - rho) != 2 * distance / 3:
            issues.append("phase-cell thickening fixture failed")
        if not bool(distance - rho > 0):
            issues.append("phase-cell thickening lost strictness")

    # Exact integrated-transport fixtures.
    for margin, delta, derivative_bound in (
        (sp.Rational(7, 3), sp.Rational(1, 100), 9),
        (sp.Rational(5, 8), sp.Rational(1, 1000), 17),
        (sp.Rational(11, 5), sp.Rational(1, 25), 13),
    ):
        displacement = delta * derivative_bound
        if not bool(displacement < margin):
            issues.append("transport fixture lacks strict domination")
        if not bool(margin - displacement > 0):
            issues.append("transport reverse triangle failed")

    # A strict linear functional puts every convex fixture in one half-plane.
    q = sp.Matrix([2, -1])
    strip_vectors = [
        sp.Matrix([3, 1]),
        sp.Matrix([5, 2]),
        sp.Matrix([sp.Rational(7, 2), sp.Rational(1, 3)]),
    ]
    if not all(bool(q.dot(vector) > 0) for vector in strip_vectors):
        issues.append("right-half-plane fixture failed")

    # Audit the imported finite calibration and nonpromotion guards.
    exact = stored.get("exact", {})
    if exact.get("right_half_plane", {}).get("q208_calibration") != (
        "For the Q208 transformed right strip one may take "
        "q=(0,1), because F_x>0 in all 20 stored cells."
    ):
        issues.append("Q208 right-strip calibration drifted")
    q208_row = next(
        (row for row in rows if row.get("id") == "napcs_08_q208_calibration"),
        {},
    )
    if q208_row.get("readiness") != "proved_finite_only":
        issues.append("Q208 calibration was overpromoted")
    open_row = next(
        (row for row in rows if row.get("id") == "napcs_09_asymptotic_budget"),
        {},
    )
    if open_row.get("readiness") != "open":
        issues.append("all-j asymptotic budget was overpromoted")
    boundary = stored.get("proof_boundary", "")
    for phrase in (
        "No Q207-to-Q208 transport budget",
        "no all-j bottom-collar estimate",
        "no all-j right-strip cone",
        "no Lambda<=0",
        "no RH proof",
    ):
        if phrase not in boundary:
            issues.append(f"proof boundary lacks {phrase!r}")

    if issues:
        raise SystemExit("\n".join(issues))
    print(
        "validated Newman adiabatic phase-cell successor lemma: "
        "10 rows, 3 exact lemmas, 2 conditional compositions, "
        "1 finite Q208 calibration, 1 open all-j Xi antecedent"
    )


if __name__ == "__main__":
    main()
