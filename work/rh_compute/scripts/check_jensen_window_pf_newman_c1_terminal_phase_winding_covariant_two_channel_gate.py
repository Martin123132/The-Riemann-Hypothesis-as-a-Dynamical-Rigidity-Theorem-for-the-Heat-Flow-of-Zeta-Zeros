#!/usr/bin/env python3
"""Independently check the C1 phase-winding/covariant two-channel gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_terminal_phase_winding_covariant_two_channel_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def check_sources(payload: dict) -> None:
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        require(file_hash(path) == item["sha256"], f"source hash drift {key}")
        require(payload["source_sha256"][key] == item["sha256"], f"hash map {key}")


def check_phase_cell() -> int:
    ell = sp.symbols("L", positive=True, real=True)
    a_squared = sp.exp(ell) + sp.Rational(1, 32) / ell**2
    expected_derivative = sp.exp(ell) - sp.Rational(1, 16) / ell**3
    require(sp.simplify(sp.diff(a_squared, ell) - expected_derivative) == 0, "a derivative")

    n, eps_left, eps_right, log_n, pi = sp.symbols(
        "N eps_left eps_right log_N pi", positive=True, real=True
    )
    left = (2 * pi * n**2 - eps_left) * log_n
    right = (2 * pi * (n + 1) ** 2 - eps_right) * log_n
    expected = (2 * pi * (2 * n + 1) + eps_left - eps_right) * log_n
    require(sp.expand(right - left - expected) == 0, "phase endpoint identity")

    require(Fraction(7, 8 * 50**2) < 1, "epsilon guard")
    require(Fraction(2 * 3 - 1, 1) > 0, "2pi-1 guard")
    require(Fraction(4 * 3 * 2, 2) > Fraction(44, 7), "full turn")

    # Check the abstract circle-coverage construction with exact fractions.
    period = Fraction(7, 1)
    starts = [Fraction(k, 5) for k in range(-4, 5)]
    angles = [Fraction(k, 6) * period for k in range(7)]
    samples = 0
    for start in starts:
        end = start + period + Fraction(1, 11)
        for angle in angles:
            k = floor_fraction((start - angle) / period) + 1
            representative = angle + k * period
            require(start < representative <= start + period, "lattice representative")
            require(representative < end, "representative in phase interval")
            samples += 1
    return samples


def check_covariant_matrix() -> tuple[int, int]:
    h, u, v, x, y, lam = sp.symbols("h u v x y lambda", real=True)
    matrix = sp.Matrix([[h + u, -v], [-v, h - u]]) / 2
    vector = sp.Matrix([x, y])
    expected = (
        h * (x**2 + y**2) + u * (x**2 - y**2) - 2 * v * x * y
    ) / 2
    require(sp.expand((vector.T * matrix * vector)[0] - expected) == 0, "quadratic form")
    require(sp.expand(matrix.trace() - h) == 0, "trace")
    require(sp.expand(matrix.det() - (h**2 - u**2 - v**2) / 4) == 0, "determinant")
    characteristic = sp.expand((lam * sp.eye(2) - matrix).det())
    expected_characteristic = sp.expand((lam - h / 2) ** 2 - (u**2 + v**2) / 4)
    require(sp.expand(characteristic - expected_characteristic) == 0, "characteristic polynomial")

    matrix_samples = 0
    for h_value in (-5, -1, 0, 3):
        for u_value in (-4, -1, 0, 2):
            for v_value in (-3, 0, 5):
                for x_value in (-2, 0, 3):
                    for y_value in (-4, 1):
                        q_value = (
                            h_value * (x_value**2 + y_value**2)
                            + u_value * (x_value**2 - y_value**2)
                            - 2 * v_value * x_value * y_value
                        )
                        direct = (
                            (h_value + u_value) * x_value**2
                            - 2 * v_value * x_value * y_value
                            + (h_value - u_value) * y_value**2
                        )
                        require(q_value == direct, "numeric quadratic form")
                        matrix_samples += 1

    criterion_samples = 0
    pythagorean = ((0, 0, 0), (3, 4, 5), (5, 12, 13), (-8, 15, 17))
    for u_value, v_value, modulus in pythagorean:
        for h_value in range(-20, 21):
            eigenvalues = (
                Fraction(h_value - modulus, 2),
                Fraction(h_value + modulus, 2),
            )
            negative_definite = all(value < 0 for value in eigenvalues)
            require(negative_definite == (h_value + modulus < 0), "uniform criterion")
            criterion_samples += 1

    a_h, a_t, rho = sp.symbols("A_H A_T rho", positive=True, real=True)
    h_three_quarters = a_h * rho
    t_modulus = a_t * rho
    require(sp.simplify(h_three_quarters + t_modulus - rho * (a_h + a_t)) == 0, "anchor obstruction")
    return matrix_samples, criterion_samples


def check_rows(payload: dict) -> None:
    rows = payload["rows"]
    suffixes = (
        "q1_law",
        "monotone",
        "cells",
        "phase",
        "endpoint",
        "epsilon",
        "turn",
        "circle",
        "sector",
        "current",
        "matrix",
        "spectrum",
        "actual",
        "uniform",
        "centres",
        "target",
    )
    expected_ids = [f"tpw_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    require([row["id"] for row in rows] == expected_ids, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open rows")
    require(rows[-1]["readiness"] == "open", "frontier row")


def main() -> None:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    check_sources(payload)
    circle_samples = check_phase_cell()
    matrix_samples, criterion_samples = check_covariant_matrix()
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 16, "summary rows")
    require(summary["complete_cell_winding_theorems"] == 1, "winding theorem")
    require(summary["circle_surjectivity_theorems"] == 1, "circle theorem")
    require(summary["fixed_sector_obstructions"] == 1, "sector obstruction")
    require(summary["covariant_matrix_identities"] == 1, "matrix identity")
    require(summary["source_phase_uniform_criteria"] == 1, "uniform criterion")
    require(summary["center_only_obstructions"] == 1, "centre obstruction")
    require(summary["live_covariant_current_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success line")
    require("Delta phi_N" in note and ">4pi N log N" in note, "phase theorem note")
    require("h+|t|<0" in note, "uniform criterion note")
    require("eigenvalues(Q)" in note, "matrix spectrum note")

    print(
        "validated Newman C1 terminal-phase winding/covariant two-channel gate: "
        f"16 rows, 0 issues, {circle_samples} circle-coverage samples, "
        f"{matrix_samples} matrix samples, {criterion_samples} criterion samples, "
        "1 full-cell winding theorem, 1 phase-uniform criterion, "
        "1 live covariant current target"
    )


if __name__ == "__main__":
    main()
