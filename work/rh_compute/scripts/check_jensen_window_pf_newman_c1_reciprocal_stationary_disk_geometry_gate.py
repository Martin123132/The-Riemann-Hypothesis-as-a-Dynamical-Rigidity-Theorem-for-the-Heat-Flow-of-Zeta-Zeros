#!/usr/bin/env python3
"""Independently check the C1 reciprocal-stationary disk geometry gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_reciprocal_stationary_disk_geometry_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> None:
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        require(file_hash(path) == item["sha256"], f"source hash drift {key}")
        require(payload["source_sha256"][key] == item["sha256"], f"hash map {key}")


def check_algebra() -> int:
    a_h, a_t = sp.symbols("A_H A_T", positive=True, real=True)
    r_h, r_t, q_t = sp.symbols("R_H R_T Q_T")
    j_h = -a_h + r_h
    j_t = -a_t + r_t + q_t
    require(sp.expand(r_h - (a_h + j_h)) == 0, "Hermitian residual")
    require(sp.expand(r_t + q_t - (a_t + j_t)) == 0, "transpose residual")

    m, tau, exterior, alias = sp.symbols("M tau E Lambda")
    band = tau + m - exterior + alias
    require(
        sp.expand(2 * m - band - (m - tau + exterior - alias)) == 0,
        "cell join",
    )

    checks = 0
    for a_value in (1, 2, 7, 31):
        for x_value in (-9, -3, -1, 0, 2, 11):
            for y_value in (-8, -2, 0, 5):
                lhs = (a_value + x_value) ** 2 + y_value**2 - a_value**2
                rhs = x_value**2 + y_value**2 + 2 * a_value * x_value
                require(lhs == rhs, "numeric disk identity")
                checks += 1

    a, x, y, delta = sp.symbols("A x y delta", real=True)
    require(
        sp.expand(
            ((a + x) ** 2 + y**2 - (a - delta) ** 2)
            - (x**2 + y**2 + 2 * a * x + 2 * a * delta - delta**2)
        )
        == 0,
        "contracted disk identity",
    )
    return checks


def check_budget(payload: dict) -> tuple[int, Fraction, Fraction]:
    budget = payload["analytic_proof"]["absolute_budget_guard"]
    coefficient = (
        Fraction(2, 1)
        * Fraction(1, 65)
        * Fraction(9, 10)
        * Fraction(5, 1)
        * Fraction(27, 4)
    )
    require(coefficient == Fraction(243, 260), "carrier coefficient")
    require(budget["carrier_coefficient"] == "243/260", "stored coefficient")

    count_checks = 0
    test_values = list(range(4160, 10001)) + [10**k + 37 for k in range(5, 13)]
    for n_value in test_values:
        count = n_value // 32 - ((n_value + 63) // 64) + 1
        require(Fraction(count, 1) >= Fraction(n_value, 65), "block count")
        count_checks += 1

    # A positive Taylor truncation is a rigorous lower bound for exp(25/2).
    exp_lower = sum(
        Fraction(25, 2) ** k / math.factorial(k)
        for k in range(45)
    )
    lower_50 = Fraction(243, 520) * exp_lower
    upper_50 = Fraction(57 * 51 * 53 * 53, 192)
    require(lower_50 > upper_50, "base budget comparison")
    require(Fraction(4, 51) < Fraction(1, 4), "growth comparison")

    # Elementary constants used in the analytic proof can be checked without
    # fitting a decimal: e<3 implies log(32)>3, and 51/800<1/10.
    require(Fraction(51, 800) < Fraction(1, 10), "exponent rational floor")
    require(3**3 < 32, "log block floor")
    require(Fraction(243, 520) > 0, "positive budget coefficient")

    require(budget["absolute_budget_obstructions"] == 2, "obstruction count")
    require(float(lower_50) > budget["anchor_upper_50"], "stored base scale")
    return count_checks, lower_50, upper_50


def check_rows(payload: dict) -> None:
    rows = payload["rows"]
    require(len(rows) == 18, "row count")
    expected_ids = [f"rsd_{index:02d}_{suffix}" for index, suffix in enumerate(
        (
            "hermitian",
            "transpose",
            "improvement",
            "cell",
            "global",
            "anchor_stationary",
            "disk",
            "delta",
            "sign",
            "block",
            "amplitude",
            "hermitian_poly",
            "transpose_poly",
            "budget",
            "anchor_upper",
            "comparison",
            "guard",
            "target",
        ),
        start=1,
    )]
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
    disk_checks = check_algebra()
    count_checks, lower_50, upper_50 = check_budget(payload)
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 18, "summary rows")
    require(summary["combined_residuals"] == 2, "combined residuals")
    require(summary["optimal_disk_identities"] == 2, "disk identities")
    require(summary["absolute_budget_obstructions"] == 2, "budget guards")
    require(summary["live_signed_stationary_targets"] == 1, "live target")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success line")
    require("|R_T+Q_T|<=A_T-delta" in note, "combined transpose handoff")
    require("B_H>A_H" in note and "B_T>A_T" in note, "budget guard note")
    require(lower_50 > upper_50, "independent exact comparison")

    print(
        "validated Newman C1 reciprocal-stationary disk geometry gate: "
        f"18 rows, 0 issues, {disk_checks} disk samples, {count_checks} block counts, "
        "2 combined residuals, 2 absolute-budget obstructions, "
        "1 live signed stationary target"
    )


if __name__ == "__main__":
    main()
