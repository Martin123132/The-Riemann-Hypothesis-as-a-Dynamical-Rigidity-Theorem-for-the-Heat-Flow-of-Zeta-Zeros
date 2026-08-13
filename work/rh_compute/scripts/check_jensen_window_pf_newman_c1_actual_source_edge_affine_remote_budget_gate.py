#!/usr/bin/env python3
"""Independently check the actual-source edge-affine remote budget."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_edge_affine_remote_budget_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 7, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")
    return len(payload["source_audit"])


def check_symbolics() -> int:
    L, X = sp.symbols("L X", positive=True, real=True)
    t = 1 / (2 * L**2)
    T = X / 2 + sp.pi / (16 * L**2)
    exponent = t * sp.pi**2 / 64 - sp.pi * T / 4 + sp.pi * X / 8
    require(sp.simplify(exponent + sp.pi**2 / (128 * L**2)) == 0, "source exponent")

    x, v, a, q, n, ux = sp.symbols("x v a q n ux", real=True)
    p = -sp.I * v * x**3 / 4 + (a + q) * x**2 / 2 + (v * ux / 2 + sp.I * n) * x
    require(sp.diff(p, x, 3) == -3 * sp.I * v / 2, "ideal cubic derivative")

    u, sigma, heat = sp.symbols("u sigma heat", positive=True, real=True)
    lam = sp.symbols("lambda", real=True)
    polynomial = sp.Function("P")(lam)
    phase = heat * lam**2 / 4 - sigma * lam
    g = sp.diff(phase, lam)
    amplitude = (sp.exp(phase) * polynomial).subs(lam, sp.log(u))
    transformed = sp.simplify(
        sp.diff(amplitude, u, 2)
        * u**2
        / sp.exp(phase.subs(lam, sp.log(u)))
    )
    expected = (
        sp.diff(polynomial, lam, 2)
        + (2 * g - 1) * sp.diff(polynomial, lam)
        + (g**2 - g + heat / 2) * polynomial
    )
    require(
        sp.simplify(transformed - expected.subs(lam, sp.log(u))) == 0,
        "amplitude second derivative",
    )

    r0, dr, rx0, drx, c, dc, delta, d = sp.symbols("r0 dr rx0 drx c dc delta d")
    r = r0 + dr
    cc = c + dc
    aa = r * cc
    qq = (r + delta) * cc + d
    nn = r * qq + (rx0 + drx) * cc
    require(sp.expand(aa - r0 * c) == sp.expand(dr * cc + r0 * dc), "A drift")
    require(
        sp.expand(qq - r0 * c) == sp.expand(dr * cc + r0 * dc + delta * cc + d),
        "Q drift",
    )
    require(
        sp.expand(nn - (r0**2 + rx0) * c)
        == sp.expand(dr * qq + r0 * (qq - r0 * c) + drx * cc + rx0 * dc),
        "N drift",
    )
    return 4


def check_rationals() -> int:
    audits = 0

    require(Fraction(22, 7) < Fraction(81, 25), "sqrt pi surrogate")
    audits += 1
    require(Fraction(51, 100) ** 3 < Fraction(9, 64), "T0 ratio")
    audits += 1
    require(Fraction(9, 5) ** 4 < 12, "frequency fourth root")
    audits += 1

    rho_coefficient = Fraction(7, 44) * Fraction(47, 50) * Fraction(2, 3)
    require(rho_coefficient == Fraction(329, 3300), "rho coefficient")
    audits += 1
    require(Fraction(198, 329) < Fraction(61, 100), "source ratio")
    audits += 1
    require(Fraction(101, 100) * Fraction(51, 100) < Fraction(3, 5), "edge value")
    audits += 1

    h = Fraction(1, 72_000_000_000)
    R = Fraction(53, 2)
    low = Fraction(321, 80)
    require(Fraction(3, 20) + 4 * h / R + low * h**2 / R**2 < Fraction(1, 6), "P bound")
    audits += 1
    require(Fraction(9, 20) + 8 * h / R + low * h**2 / R**2 < Fraction(1, 2), "P prime")
    audits += 1
    require(Fraction(9, 10) + 8 * h / R < 1, "P second")
    audits += 1
    require(Fraction(1, 6) + Fraction(1, 2) / R < Fraction(1, 5), "A prime")
    audits += 1
    require(Fraction(1, 1) / R**2 + Fraction(3, 2) / R + Fraction(3, 8) < Fraction(1, 2), "A second")
    audits += 1

    L, y = sp.symbols("L y", real=True, nonnegative=True)
    anchor_lower = sp.Rational(7, 2816) * L**2 * (L - 2) ** 2
    ideal_upper = (L + 3) ** 3 / 81 + 1
    difference = sp.together(sp.Rational(13, 100) * anchor_lower - ideal_upper)
    shifted = sp.Poly(sp.expand(sp.numer(difference).subs(L, y + 50)), y)
    require(all(coefficient > 0 for coefficient in shifted.all_coeffs()), "ideal anchor ratio")
    audits += 1

    correction_l50 = h * 120**2 * R**3 * 301 / 6
    require(correction_l50 < Fraction(1, 5), "correction remote")
    audits += 1
    correction_log_derivative = -Fraction(1, 2) + Fraction(2, 120) + Fraction(3, 53) + Fraction(6, 301)
    require(correction_log_derivative < 0, "correction monotonicity")
    audits += 1

    affine_remote = 2 * Fraction(61, 100) * Fraction(131, 1000)
    require(affine_remote == Fraction(7991, 50000), "remote product")
    require(affine_remote < Fraction(4, 25), "remote close")
    audits += 1
    closed = Fraction(1, 6) + Fraction(1, 17160) + Fraction(1, 2000) + Fraction(4, 25)
    require(closed == Fraction(280759, 858000), "closed fraction")
    require(Fraction(98, 100) - closed == Fraction(560081, 858000), "reserve fraction")
    audits += 1
    return audits


def check_rows(payload: dict) -> None:
    suffixes = (
        "source",
        "terminal",
        "ratio",
        "edge",
        "split",
        "polynomial",
        "amplitude",
        "pairing",
        "ideal",
        "drift",
        "drift_remote",
        "actual",
        "budget",
        "near",
        "close",
        "target",
    )
    expected = [f"ear_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
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
    rational_audits = check_rationals()
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 16, "summary rows")
    require(summary["source_scale_ratios"] == 1, "source ratios")
    require(summary["genuine_edge_coefficient_bounds"] == 4, "edge bounds")
    require(summary["paired_denominator_cancellations"] == 1, "pair cancellation")
    require(summary["absolute_affine_remote_budgets"] == 1, "remote budgets")
    require(summary["affine_remote_fraction"] == "4/25", "remote fraction")
    require(summary["closed_secondary_fraction"] == "280759/858000", "closed fraction")
    require(summary["remaining_reserve_fraction"] == "560081/858000", "reserve fraction")
    require(summary["remaining_absolute_secondary_packages"] == 2, "remaining packages")
    require(summary["live_signed_near_affine_targets"] == 1, "live targets")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("4rho A_T/25" in note, "remote budget note")
    require("not a proof of RH" in note, "proof boundary note")

    print(
        "validated Newman C1 actual-source edge-affine remote budget gate: "
        f"16 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{rational_audits} rational audits, {source_audits} source audits, "
        "B_0T_0/rho<61/100, |E_aff^C|<4rho*A_T/25, "
        "2 remaining absolute secondary packages, 1 live signed near-affine target"
    )


if __name__ == "__main__":
    main()
