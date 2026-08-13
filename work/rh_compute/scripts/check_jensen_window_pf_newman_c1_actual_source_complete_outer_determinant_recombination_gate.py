#!/usr/bin/env python3
"""Independently check the complete outer-determinant recombination gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_complete_outer_determinant_recombination_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def real_vector(prefix: str) -> sp.Matrix:
    return sp.Matrix(sp.symbols(f"{prefix}_V {prefix}_N {prefix}_A {prefix}_Q", real=True))


def bilinear(left: sp.Matrix, J: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand((left.T * J * right)[0])


def increment_current(
    base: sp.Matrix,
    increment: sp.Matrix,
    base_dot: sp.Matrix,
    increment_dot: sp.Matrix,
    J: sp.Matrix,
) -> sp.Expr:
    return sp.expand(
        2 * bilinear(base, J, increment_dot)
        + 2 * bilinear(base_dot, J, increment)
        + 2 * bilinear(increment, J, increment_dot)
    )


def rational_vector(values: tuple[Fraction, Fraction, Fraction, Fraction]) -> sp.Matrix:
    return sp.Matrix([sp.Rational(value.numerator, value.denominator) for value in values])


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 11, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"source hash drift {key}")
        require(digest == payload["source_sha256"][key], f"source map drift {key}")
    return len(payload["source_audit"])


def check_symbolics() -> int:
    J = sp.Rational(1, 2) * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    p, a, o, near, remote, edge = (
        real_vector(name) for name in ("cp", "ca", "co", "cn", "cc", "cedge")
    )
    dp, da, do, dnear, dremote = (
        real_vector(name) for name in ("cdp", "cda", "cdo", "cdn", "cdc")
    )
    e, de = o - a, do - da
    audits = 0

    require(sp.simplify(a + e - o) == sp.zeros(4, 1), "owner sum")
    audits += 1
    r_nr = 2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e)
    r_quad = 2 * bilinear(e, J, de)
    r_mixed = 2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a)
    r_aa = 2 * bilinear(a, J, da)
    outer = increment_current(p, o, dp, do, J)
    require(sp.simplify(r_nr + r_quad + r_mixed + r_aa - outer) == 0, "outer recombination")
    audits += 1

    t_edge = 2 * J * edge
    affine = sp.expand(t_edge.dot(do))
    require(sp.simplify(affine - 2 * bilinear(edge, J, do)) == 0, "edge row")
    audits += 1
    base = edge + p
    geometric = sp.expand(outer + affine)
    shifted = increment_current(base, o, dp, do, J)
    require(sp.simplify(geometric - shifted) == 0, "edge shift")
    audits += 1
    primitive = 2 * (bilinear(base + o, J, dp + do) - bilinear(base, J, dp))
    require(sp.simplify(geometric - primitive) == 0, "primitive")
    audits += 1

    split = {o[index]: near[index] + remote[index] for index in range(4)}
    split.update({do[index]: dnear[index] + dremote[index] for index in range(4)})
    geometric_split = sp.expand(geometric.subs(split, simultaneous=True))
    near_current = increment_current(base, near, dp, dnear, J)
    remote_current = increment_current(base + near, remote, dp + dnear, dremote, J)
    require(sp.simplify(geometric_split - near_current - remote_current) == 0, "near remote")
    audits += 1

    defect, ddefect = real_vector("cd"), real_vector("cdd")
    finite = {o[index]: -defect[index] for index in range(4)}
    finite.update({do[index]: -ddefect[index] for index in range(4)})
    geometric_finite = sp.expand(geometric.subs(finite, simultaneous=True))
    finite_current = increment_current(base, -defect, dp, -ddefect, J)
    require(sp.simplify(geometric_finite - finite_current) == 0, "finite defect")
    audits += 1

    band, dband = real_vector("cf"), real_vector("cdf")
    defect_to_band = {defect[index]: band[index] - p[index] for index in range(4)}
    defect_to_band.update(
        {ddefect[index]: dband[index] - dp[index] for index in range(4)}
    )
    finite_band = sp.expand(finite_current.subs(defect_to_band, simultaneous=True))
    expected_band = 2 * (
        bilinear(edge + 2 * p - band, J, 2 * dp - dband)
        - bilinear(edge + p, J, dp)
    )
    require(sp.simplify(finite_band - expected_band) == 0, "finite band")
    audits += 1

    reverse = increment_current(base, remote, dp, dremote, J)
    reverse += increment_current(base + remote, near, dp + dremote, dnear, J)
    require(sp.simplify(near_current + remote_current - reverse) == 0, "allocation order")
    audits += 1
    return audits


def check_finite(payload: dict) -> int:
    J = sp.Rational(1, 2) * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    p = rational_vector((Fraction(1, 2), Fraction(-2, 3), Fraction(3, 5), Fraction(-4, 7)))
    dp = rational_vector((Fraction(5, 11), Fraction(-6, 13), Fraction(7, 17), Fraction(-8, 19)))
    a = rational_vector((Fraction(1, 3), Fraction(2, 5), Fraction(-3, 7), Fraction(4, 9)))
    da = rational_vector((Fraction(-5, 11), Fraction(7, 13), Fraction(-11, 17), Fraction(13, 19)))
    near = rational_vector((Fraction(-2, 5), Fraction(3, 7), Fraction(-5, 9), Fraction(7, 11)))
    dnear = rational_vector((Fraction(11, 13), Fraction(-13, 17), Fraction(17, 19), Fraction(-19, 23)))
    remote = rational_vector((Fraction(3, 11), Fraction(-5, 13), Fraction(7, 17), Fraction(-11, 19)))
    dremote = rational_vector((Fraction(-13, 23), Fraction(17, 29), Fraction(-19, 31), Fraction(23, 37)))
    edge = rational_vector((Fraction(2, 3), Fraction(-3, 5), Fraction(5, 7), Fraction(-7, 11)))
    o, do = near + remote, dnear + dremote
    e, de = o - a, do - da
    base = edge + p

    r_nr = 2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e)
    r_quad = 2 * bilinear(e, J, de)
    r_mixed = 2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a)
    r_aa = 2 * bilinear(a, J, da)
    outer = increment_current(p, o, dp, do, J)
    affine = 2 * bilinear(edge, J, do)
    geometric = sp.expand(outer + affine)
    near_current = increment_current(base, near, dp, dnear, J)
    remote_current = increment_current(base + near, remote, dp + dnear, dremote, J)
    defect, ddefect = -o, -do
    finite_current = increment_current(base, -defect, dp, -ddefect, J)
    band, dband = p + defect, dp + ddefect
    band_current = 2 * (
        bilinear(edge + 2 * p - band, J, 2 * dp - dband)
        - bilinear(edge + p, J, dp)
    )

    require(sp.simplify(r_nr + r_quad + r_mixed + r_aa - outer) == 0, "finite packages")
    require(sp.simplify(outer + affine - geometric) == 0, "finite edge")
    require(sp.simplify(geometric - near_current - remote_current) == 0, "finite near remote")
    require(sp.simplify(geometric - finite_current) == 0, "finite defect")
    require(sp.simplify(geometric - band_current) == 0, "finite band")
    require(geometric != 0, "finite nonzero")

    stored = payload["finite_rational"]
    require(stored["outer"] == str(outer), "stored outer")
    require(stored["affine"] == str(affine), "stored affine")
    require(stored["geometric"] == str(geometric), "stored geometric")
    require(stored["near_current"] == str(near_current), "stored near")
    require(stored["remote_current"] == str(remote_current), "stored remote")
    require(stored["finite_defect_current"] == str(finite_current), "stored defect")
    require(stored["finite_band_current"] == str(band_current), "stored band")
    return 6


def check_rows(payload: dict) -> None:
    suffixes = (
        "domain",
        "owners",
        "terminal",
        "outer",
        "packages",
        "edge",
        "geom",
        "collapse",
        "outer_split",
        "near",
        "remote",
        "order",
        "defect",
        "band",
        "ties",
        "threshold",
        "target",
        "boundary",
    )
    expected = [f"odr_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
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
    finite_audits = check_finite(payload)
    check_rows(payload)

    exact = payload["exact"]
    require(exact["collapsed_packages"] == 8, "collapsed packages")
    require(exact["remaining_pointwise_packages"] == 2, "remaining packages")
    require((exact["necessary_numerator"], exact["necessary_denominator"]) == (97, 100), "necessary threshold")
    require((exact["sufficient_numerator"], exact["sufficient_denominator"]) == (99, 100), "sufficient threshold")
    require("B=-49/100" in exact["current_criterion"], "calibrated level condition")

    summary = payload["summary"]
    require(summary["rows"] == 18, "summary rows")
    require(summary["symbolic_audits"] == 9, "summary symbolics")
    require(summary["finite_rational_audits"] == 6, "summary finite")
    require(summary["source_audits"] == 11, "summary sources")
    require(summary["collapsed_packages"] == 8, "summary collapse")
    require(summary["remaining_pointwise_packages"] == 2, "summary remaining")
    require(summary["near_remote_determinant_splits"] == 1, "near remote split")
    require(summary["finite_cell_determinant_forms"] == 1, "finite-cell form")
    require(summary["necessary_geometric_fraction"] == "97/100", "necessary fraction")
    require(summary["sufficient_geometric_fraction"] == "99/100", "sufficient fraction")
    require(summary["open_rows"] == 1, "open summary")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("E_geom^cur+E_move^cur" in note, "two-package ledger note")
    require("97/100" in note and "99/100" in note, "threshold note")
    require("B=-49/100" in note, "calibrated level note")
    require("common-cutoff paired-remote" in note, "remote ownership note")
    require("not a proof of RH" in note, "RH boundary")
    require("prize-level conclusion" in note, "prize boundary")

    print(
        "validated Newman C1 complete outer-determinant recombination gate: "
        f"18 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{finite_audits} exact-rational audits, {source_audits} source audits, "
        "8 packages collapsed, 2 pointwise packages remain, "
        "1 finite-cell determinant form, 1 open arithmetic row"
    )


if __name__ == "__main__":
    main()
