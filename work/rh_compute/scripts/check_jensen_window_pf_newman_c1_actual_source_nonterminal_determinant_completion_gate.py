#!/usr/bin/env python3
"""Independently check the nonterminal determinant-completion gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_nonterminal_determinant_completion_gate"
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


def rational_vector(values: tuple[Fraction, Fraction, Fraction, Fraction]) -> sp.Matrix:
    return sp.Matrix([sp.Rational(value.numerator, value.denominator) for value in values])


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
    J = sp.Rational(1, 2) * sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    p, e, a = (real_vector(name) for name in ("check_p", "check_e", "check_a"))
    dp, de, da = (real_vector(name) for name in ("check_dp", "check_de", "check_da"))
    t_e = real_vector("check_tE")

    raw = sp.expand(
        2 * bilinear(p, J, da + de)
        + 2 * bilinear(dp, J, a + e)
        + 2 * bilinear(a + e, J, da + de)
        + t_e.dot(da + de)
    )
    r_nr = sp.expand(2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e))
    r_quad = sp.expand(2 * bilinear(e, J, de))
    r_mixed = sp.expand(2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a))
    r_aa = sp.expand(2 * bilinear(a, J, da))
    r_aff = sp.expand(t_e.dot(da + de))

    audits = 0
    require(sp.simplify(raw - r_nr - r_quad - r_mixed - r_aa - r_aff) == 0, "five-block split")
    audits += 1
    require(sp.simplify(raw - (r_quad + r_mixed + r_aa + r_aff) - r_nr) == 0, "ledger defect")
    audits += 1
    require(sp.simplify(r_nr - 2 * bilinear(dp, J, e) - 2 * bilinear(p, J, de)) == 0, "cross primitive")
    audits += 1
    completed = 2 * (bilinear(p + e, J, dp + de) - bilinear(p, J, dp))
    require(sp.simplify(completed - r_nr - r_quad) == 0, "determinant completion")
    audits += 1

    alpha_p = sp.Matrix([dp[1], dp[0], -dp[3], -dp[2]])
    beta_p = sp.Matrix([p[1], p[0], -p[3], -p[2]])
    require(sp.simplify(alpha_p - 2 * J * dp) == sp.zeros(4, 1), "alpha row")
    audits += 1
    require(sp.simplify(beta_p - 2 * J * p) == sp.zeros(4, 1), "beta row")
    audits += 1
    require(sp.simplify(alpha_p.dot(e) + beta_p.dot(de) - r_nr) == 0, "P_nr projection")
    audits += 1

    beta_edge = real_vector("check_betaE")
    beta_full = beta_p + beta_edge
    require(sp.simplify(beta_full - beta_edge - beta_p) == sp.zeros(4, 1), "beta split")
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
    e = rational_vector((Fraction(-2, 5), Fraction(3, 7), Fraction(-5, 9), Fraction(7, 11)))
    de = rational_vector((Fraction(11, 13), Fraction(-13, 17), Fraction(17, 19), Fraction(-19, 23)))
    a = rational_vector((Fraction(1, 3), Fraction(2, 5), Fraction(-3, 7), Fraction(4, 9)))
    da = rational_vector((Fraction(-5, 11), Fraction(7, 13), Fraction(-11, 17), Fraction(13, 19)))
    t_e = rational_vector((Fraction(2, 3), Fraction(-3, 5), Fraction(5, 7), Fraction(-7, 11)))

    raw = sp.expand(
        2 * bilinear(p, J, da + de)
        + 2 * bilinear(dp, J, a + e)
        + 2 * bilinear(a + e, J, da + de)
        + t_e.dot(da + de)
    )
    r_nr = sp.expand(2 * bilinear(p, J, de) + 2 * bilinear(dp, J, e))
    r_quad = sp.expand(2 * bilinear(e, J, de))
    r_mixed = sp.expand(2 * bilinear(p + e, J, da) + 2 * bilinear(dp + de, J, a))
    r_aa = sp.expand(2 * bilinear(a, J, da))
    r_aff = sp.expand(t_e.dot(da + de))
    represented_old = sp.expand(r_quad + r_mixed + r_aa + r_aff)
    completed = sp.expand(2 * (bilinear(p + e, J, dp + de) - bilinear(p, J, dp)))
    projection = sp.expand((2 * J * dp).dot(e) + (2 * J * p).dot(de))

    require(sp.simplify(raw - represented_old - r_nr) == 0, "finite split")
    require(sp.simplify(raw - represented_old) == r_nr, "finite defect")
    require(sp.simplify(projection - r_nr) == 0, "finite projection")
    require(sp.simplify(completed - r_nr - r_quad) == 0, "finite completion")
    require(r_nr != 0 and sp.simplify((r_nr + r_quad) - r_quad) == r_nr, "finite failed join")

    stored = payload["finite_rational"]
    require(stored["raw"] == str(raw), "stored raw")
    require(stored["represented_old"] == str(represented_old), "stored old")
    require(stored["missing_retained_nonterminal"] == str(r_nr), "stored missing")
    require(stored["quadratic"] == str(r_quad), "stored quadratic")
    require(stored["completed_nonterminal"] == str(completed), "stored completed")
    require(stored["polynomial_projection"] == str(projection), "stored projection")
    require(stored["zero_terminal_edge_proposed"] == str(r_quad), "stored proposed")
    require(stored["zero_terminal_edge_complete"] == str(r_nr + r_quad), "stored complete")
    return 5


def check_rows(payload: dict) -> None:
    suffixes = (
        "domain",
        "raw",
        "split",
        "terminal",
        "affine",
        "quadratic",
        "missing",
        "rows",
        "polynomial",
        "functional",
        "primitive",
        "completion",
        "falsification",
        "ledger",
        "main",
        "target",
        "boundary",
    )
    expected = [f"ndc_{index:02d}_{suffix}" for index, suffix in enumerate(suffixes, start=1)]
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
    require(exact["corrected_package_count"] == 9, "package count")
    require(exact["missing_linear_packages"] == 1, "missing package count")
    require(exact["determinant_completions"] == 1, "completion count")
    require((exact["necessary_numerator"], exact["necessary_denominator"]) == (551501, 858000), "necessary threshold")
    require((exact["sufficient_numerator"], exact["sufficient_denominator"]) == (1130179, 858000), "sufficient threshold")

    summary = payload["summary"]
    require(summary["rows"] == 17, "summary rows")
    require(summary["symbolic_audits"] == 8, "summary symbolics")
    require(summary["finite_rational_audits"] == 5, "summary finite")
    require(summary["source_audits"] == 7, "summary sources")
    require(summary["prior_package_count"] == 8, "prior package count")
    require(summary["corrected_package_count"] == 9, "corrected package count")
    require(summary["missing_linear_packages"] == 1, "summary missing")
    require(summary["determinant_completions"] == 1, "summary completion")
    require(summary["necessary_signed_fraction"] == "551501/858000", "necessary fraction")
    require(summary["sufficient_signed_fraction"] == "1130179/858000", "sufficient fraction")
    require(summary["open_rows"] == 1, "open summary")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("E_FAN+E_quad" in note and "not a determinant-current identity" in note, "failed-test note")
    require("E_nr+E_quad" in note, "completed determinant note")
    require("not a proof of RH" in note, "RH boundary")
    require("prize-level conclusion" in note, "prize boundary")

    print(
        "validated Newman C1 nonterminal determinant-completion gate: "
        f"17 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{finite_audits} exact-rational audits, {source_audits} source audits, "
        "1 missing linear package identified, 1 corrected determinant completion, "
        "1 open arithmetic row"
    )


if __name__ == "__main__":
    main()
