#!/usr/bin/env python3
"""Independently check current/transfer ownership and the quadratic primitive."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_current_transfer_quadratic_primitive_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 6, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")
    return len(payload["source_audit"])


def check_symbolics() -> int:
    audits = 0
    x = sp.symbols("x0:4", real=True)
    y = sp.symbols("y0:4", real=True)
    xv = sp.Matrix(x)
    yv = sp.Matrix(y)
    j = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    determinant = x[0] * x[1] - x[2] * x[3]
    require(sp.expand((xv.T * j * xv)[0]) == determinant, "J determinant")
    audits += 1
    expected = x[0] * y[1] + x[1] * y[0] - x[2] * y[3] - x[3] * y[2]
    derivative = sum(sp.diff(determinant, x[index]) * y[index] for index in range(4))
    require(sp.expand(derivative) == expected, "determinant derivative")
    audits += 1
    require(sp.expand(2 * (xv.T * j * yv)[0]) == expected, "J derivative")
    audits += 1

    h, r, k = sp.symbols("h r k", positive=True, real=True)
    p_star = r**2 / 2
    s0 = sp.Rational(16, 7) * h ** sp.Rational(-1, 2)
    t1 = 4 * h ** sp.Rational(3, 2) * k**2
    require(sp.simplify(8 * p_star**2 * s0 * t1) == sp.Rational(128, 7) * h * r**4 * k**2, "carrier")
    audits += 1
    require(
        sp.simplify(sp.Rational(24, 5) * p_star * t1)
        == sp.Rational(48, 5) * h ** sp.Rational(3, 2) * r**2 * k**2,
        "edge",
    )
    audits += 1

    vr, vi, nr, ni, ar, ai, qr, qi = sp.symbols("vr vi nr ni ar ai qr qi", real=True)
    zv = vr + sp.I * vi
    zn = nr + sp.I * ni
    za = ar + sp.I * ai
    zq = qr + sp.I * qi
    polarized = sp.re(zv * sp.conjugate(zn) - za * sp.conjugate(zq) + zv * zn - za * zq).expand(complex=True)
    require(sp.simplify(polarized - 2 * (vr * nr - ar * qr)) == 0, "polarization")
    audits += 1
    return audits


def check_rationals(payload: dict) -> int:
    audits = 0
    h = Fraction(1, 72_000_000_000)
    l = Fraction(50)
    r = Fraction(53, 2)
    exp25_lower = sum(Fraction(25) ** index / math.factorial(index) for index in range(61))
    require(exp25_lower > 72_000_000_000, "exp floor")
    audits += 1
    require(72_000_000_000**7 > 13_000**18, "fractional-power floor")
    audits += 1

    current_main = Fraction(422_400, 7 * 329 * 14_300) * r**4 / (l * 13_000)
    current_edge = Fraction(2_928, 5 * 100 * 14_300) * h * r**2
    require(current_main + current_edge < Fraction(1, 100), "pointwise budget")
    audits += 1

    exact = payload["exact"]
    require(Fraction(exact["current_main_numerator"], exact["current_main_denominator"]) == current_main, "stored main")
    require(Fraction(exact["current_edge_numerator"], exact["current_edge_denominator"]) == current_edge, "stored edge")
    audits += 1

    closed = Fraction(280_759, 858_000) + Fraction(1, 100)
    reserve = Fraction(98, 100) - closed
    sufficient = Fraction(98, 100) + closed
    require(closed == Fraction(289_339, 858_000), "closed")
    require(reserve == Fraction(551_501, 858_000), "reserve")
    require(sufficient == Fraction(1_130_179, 858_000), "sufficient")
    require(Fraction(exact["closed_numerator"], exact["closed_denominator"]) == closed, "stored closed")
    require(Fraction(exact["reserve_numerator"], exact["reserve_denominator"]) == reserve, "stored reserve")
    require(Fraction(exact["sufficient_numerator"], exact["sufficient_denominator"]) == sufficient, "stored sufficient")
    audits += 1

    require(-Fraction(7, 36) + Fraction(4, 53) - Fraction(1, 50) < 0, "main monotonicity")
    require(-Fraction(25, 36) + Fraction(2, 53) < 0, "edge monotonicity")
    audits += 1
    return audits


def check_rows(payload: dict) -> None:
    suffixes = (
        "domain",
        "prior",
        "current",
        "bound",
        "ratio",
        "monotone",
        "endpoint",
        "move",
        "closed",
        "four",
        "primitive",
        "phase",
        "transfer",
        "exclusion",
        "reclass",
        "handoff",
        "boundary",
    )
    expected = [f"ctq_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    require([row["id"] for row in payload["rows"]] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in payload["rows"]) == 1, "open rows")
    require(payload["rows"][-2]["readiness"] == "open", "frontier row")


def main() -> None:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    source_audits = check_sources(payload)
    symbolic_audits = check_symbolics()
    rational_audits = check_rationals(payload)
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 17, "summary rows")
    require(summary["ownership_corrections"] == 1, "ownership")
    require(summary["valid_transfer_budgets_retained"] == 1, "transfer")
    require(summary["pointwise_moving_budget_denominator"] == 100, "moving budget")
    require(summary["quadratic_constituents"] == 4, "quadratic constituents")
    require(summary["quadratic_primitives"] == 1, "primitive")
    require(summary["corrected_closed_fraction"] == "289339/858000", "closed fraction")
    require(summary["necessary_signed_fraction"] == "551501/858000", "necessary fraction")
    require(summary["sufficient_signed_fraction"] == "1130179/858000", "sufficient fraction")
    require(summary["open_rows"] == 1, "open summary")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("E_move^cur" in note, "current notation")
    require("E_move^tr" in note, "transfer notation")
    require("determinant" in note, "primitive note")
    require("not a proof" in note, "proof boundary")

    print(
        "validated Newman C1 current/transfer quadratic-primitive gate: "
        f"17 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{rational_audits} rational audits, {source_audits} source audits, "
        "pointwise |E_move^cur|<rho*A_T/100, 1 determinant primitive, "
        "1 signed-main reclassification, 1 open theorem row"
    )


if __name__ == "__main__":
    main()
