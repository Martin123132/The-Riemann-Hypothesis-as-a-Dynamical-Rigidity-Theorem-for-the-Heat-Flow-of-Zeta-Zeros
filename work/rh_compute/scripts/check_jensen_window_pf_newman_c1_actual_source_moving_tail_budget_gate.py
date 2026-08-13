#!/usr/bin/env python3
"""Independently check the actual-source moving-tail common-unit budget."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_moving_tail_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 8, "source count")
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
    require(sp.expand((xv.T * j * xv)[0]) == x[0] * x[1] - x[2] * x[3], "J quadratic")
    audits += 1
    require(
        sp.expand(2 * (xv.T * j * yv)[0])
        == x[0] * y[1] + x[1] * y[0] - x[2] * y[3] - x[3] * y[2],
        "J bilinear",
    )
    audits += 1

    h, l, r, k = sp.symbols("h l r k", positive=True, real=True)
    eps = sp.Rational(7, 8) * h**2 / l**2
    p_star = r**2 / 2
    s0 = sp.Rational(16, 7) * h ** sp.Rational(-1, 2)
    s1 = sp.Rational(16, 7) * r * h ** sp.Rational(-1, 2)
    t1 = 4 * h ** sp.Rational(3, 2) * k**2
    require(sp.simplify(8 * eps * p_star**2 * s0 * t1) == 16 * h**3 * r**4 * k**2 / l**2, "main term")
    audits += 1
    require(sp.simplify(8 * eps**2 * p_star**2 * t1 * s1) == 14 * h**5 * r**5 * k**2 / l**4, "displacement term")
    audits += 1
    require(
        sp.simplify(sp.Rational(24, 5) * eps * p_star * t1)
        == sp.Rational(42, 5) * h ** sp.Rational(7, 2) * r**2 * k**2 / l**2,
        "edge term",
    )
    audits += 1
    return audits


def check_rationals(payload: dict) -> int:
    audits = 0
    h = Fraction(1, 72_000_000_000)
    l = Fraction(50)
    r = Fraction(53, 2)

    require(Fraction(4) < Fraction(120**2, 1000), "delta absorption")
    audits += 1
    require(Fraction(1001, 2000) < Fraction(3, 5), "A envelope")
    audits += 1
    q_upper = (r / 2 + Fraction(1, 1_000_000)) * Fraction(1001, 1000) + Fraction(1, 1_000_000)
    require(q_upper < Fraction(3, 5) * r, "Q envelope")
    audits += 1
    n_upper = r * Fraction(3, 5) * r / 2 + Fraction(1, 500_000) * Fraction(1001, 1000)
    require(n_upper < r**2 / 2, "N envelope")
    audits += 1
    require(Fraction(1001, 1000) < r**2 / 2, "V envelope")
    audits += 1

    exp25_lower = sum(Fraction(25) ** k / math.factorial(k) for k in range(61))
    require(exp25_lower > 72_000_000_000, "exp(25) floor")
    audits += 1

    main = Fraction(52800, 329 * 14300) * h**2 * r**4 / l**3
    displacement = Fraction(46200, 329 * 14300) * h**4 * r**5 / l**5
    edge = Fraction(2562, 500 * 14300) * h**3 * r**2 / l**2
    require(main + displacement + edge < Fraction(1, 1_000_000), "endpoint budget")
    audits += 1

    exact = payload["exact"]
    require(Fraction(exact["endpoint_main_numerator"], exact["endpoint_main_denominator"]) == main, "stored main")
    require(
        Fraction(exact["endpoint_displacement_numerator"], exact["endpoint_displacement_denominator"])
        == displacement,
        "stored displacement",
    )
    require(Fraction(exact["endpoint_edge_numerator"], exact["endpoint_edge_denominator"]) == edge, "stored edge")
    audits += 1

    require(Fraction(198, 329) < Fraction(61, 100), "source ratio")
    audits += 1
    edge_source = Fraction(42, 5) * Fraction(61, 100)
    require(edge_source == Fraction(2562, 500), "edge source product")
    audits += 1

    closed = Fraction(280759, 858000) + Fraction(1, 1_000_000)
    reserve = Fraction(98, 100) - closed
    require(closed == Fraction(140379929, 429000000), "closed budget")
    require(reserve == Fraction(280040071, 429000000), "reserve")
    require(Fraction(exact["closed_numerator"], exact["closed_denominator"]) == closed, "stored closed")
    require(Fraction(exact["reserve_numerator"], exact["reserve_denominator"]) == reserve, "stored reserve")
    audits += 1

    require(-Fraction(43, 36) + Fraction(4, 53) < 0, "main monotonicity")
    require(-Fraction(79, 36) < 0, "displacement monotonicity")
    require(-Fraction(61, 36) < 0, "edge monotonicity")
    audits += 1
    return audits


def check_rows(payload: dict) -> None:
    suffixes = (
        "ownership",
        "definition",
        "weight",
        "observation",
        "mass",
        "moment",
        "terminal",
        "displacement",
        "bilinear",
        "pointwise",
        "integrated",
        "main",
        "delta",
        "edge",
        "close",
        "target",
    )
    expected = [f"mtb_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
    rows = payload["rows"]
    require([row["id"] for row in rows] == expected, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open rows")
    require(rows[-1]["readiness"] == "open", "frontier row")


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
    require(summary["rows"] == 16, "summary rows")
    require(summary["exact_moving_tail_definitions"] == 1, "definitions")
    require(summary["source_free_observation_envelopes"] == 4, "observation envelopes")
    require(summary["source_free_moment_envelopes"] == 3, "moment envelopes")
    require(summary["monotone_normalized_terms"] == 3, "normalized terms")
    require(summary["moving_tail_budget_denominator"] == 1_000_000, "moving budget")
    require(summary["closed_secondary_fraction"] == "140379929/429000000", "closed fraction")
    require(summary["remaining_reserve_fraction"] == "280040071/429000000", "reserve fraction")
    require(summary["remaining_absolute_secondary_packages"] == 1, "remaining packages")
    require(summary["live_signed_near_affine_targets"] == 1, "live targets")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("rho A_T/1000000" in note, "moving budget note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source moving-tail budget gate: "
        f"16 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{rational_audits} rational audits, {source_audits} source audits, "
        "|E_move|<rho*A_T/1000000, 1 remaining absolute secondary package, "
        "1 live signed near-affine target"
    )


if __name__ == "__main__":
    main()
