#!/usr/bin/env python3
"""Independently check the finite determinant reciprocal-block symmetry gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_c1_actual_source_finite_determinant_reciprocal_block_symmetry_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, message: str) -> None:
    require(sp.simplify(sp.expand(expression)) == 0, message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def real_vector(prefix: str) -> sp.Matrix:
    return sp.Matrix(sp.symbols(f"{prefix}_V {prefix}_N {prefix}_A {prefix}_Q", real=True))


def observation_matrix(vector: sp.Matrix) -> sp.Matrix:
    return sp.Matrix([[vector[0], vector[2]], [vector[3], vector[1]]])


def bilinear(left: sp.Matrix, J: sp.Matrix, right: sp.Matrix) -> sp.Expr:
    return sp.expand((left.T * J * right)[0])


def quadratic(vector: sp.Matrix, J: sp.Matrix) -> sp.Expr:
    return bilinear(vector, J, vector)


def primitive_increment(base: sp.Matrix, defect: sp.Matrix, J: sp.Matrix) -> sp.Expr:
    return sp.expand(quadratic(base - defect, J) - quadratic(base, J))


def current_increment(
    base: sp.Matrix,
    defect: sp.Matrix,
    base_dot: sp.Matrix,
    defect_dot: sp.Matrix,
    J: sp.Matrix,
) -> sp.Expr:
    return sp.expand(
        -2 * bilinear(base, J, defect_dot)
        - 2 * bilinear(base_dot, J, defect)
        + 2 * bilinear(defect, J, defect_dot)
    )


def rational_vector(values: tuple[Fraction, Fraction, Fraction, Fraction]) -> sp.Matrix:
    return sp.Matrix([sp.Rational(value.numerator, value.denominator) for value in values])


def check_sources(payload: dict) -> int:
    require(set(payload["source_sha256"]) == set(payload["source_audit"]), "source keys")
    require(len(payload["source_audit"]) == 3, "source count")
    for key, item in payload["source_audit"].items():
        path = REPO_ROOT / item["path"]
        require(path.is_file(), f"missing source {key}")
        digest = file_hash(path)
        require(digest == item["sha256"], f"stored source hash {key}")
        require(digest == payload["source_sha256"][key], f"source hash map {key}")
    return 3


def check_symbolic(payload: dict) -> int:
    J = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    base = real_vector("cb")
    pieces = [real_vector(f"ce{index}") for index in range(3)]
    dbase = real_vector("cdb")
    dpieces = [real_vector(f"cde{index}") for index in range(3)]
    defect = sum(pieces, sp.zeros(4, 1))
    ddefect = sum(dpieces, sp.zeros(4, 1))
    audits = 0

    B, D = observation_matrix(base), observation_matrix(defect)
    require_zero(B.det() - quadratic(base, J), "matrix determinant")
    audits += 1
    require_zero(sp.trace(B.adjugate() * D) - 2 * bilinear(base, J, defect), "adjugate")
    audits += 1
    primitive = primitive_increment(base, defect, J)
    require_zero(primitive - D.det() + sp.trace(B.adjugate() * D), "division-free increment")
    audits += 1

    sequential_primitive = 0
    running_base = base
    for piece in pieces:
        sequential_primitive += primitive_increment(running_base, piece, J)
        running_base -= piece
    require_zero(primitive - sequential_primitive, "primitive telescope")
    audits += 1

    generic_base, generic_defect = real_vector("cg"), real_vector("cr")
    generic_dbase, generic_ddefect = real_vector("cdg"), real_vector("cdr")
    current = current_increment(generic_base, generic_defect, generic_dbase, generic_ddefect, J)
    direct = 2 * bilinear(
        generic_base - generic_defect,
        J,
        generic_dbase - generic_ddefect,
    ) - 2 * bilinear(generic_base, J, generic_dbase)
    require_zero(current - direct, "current formula")
    audits += 1

    total_current = current_increment(base, defect, dbase, ddefect, J)
    sequential_current = 0
    running_base, running_dot = base, dbase
    for piece, dpiece in zip(pieces, dpieces, strict=True):
        sequential_current += current_increment(running_base, piece, running_dot, dpiece, J)
        running_base -= piece
        running_dot -= dpiece
    require_zero(total_current - sequential_current, "current telescope")
    audits += 1

    naive = sum(primitive_increment(base, piece, J) for piece in pieces)
    cross = 2 * sum(
        bilinear(pieces[left], J, pieces[right])
        for left in range(3)
        for right in range(left + 1, 3)
    )
    require_zero(primitive - naive - cross, "cross ownership")
    audits += 1

    l11, l12, l21, l22 = sp.symbols("cl11 cl12 cl21 cl22", real=True)
    r11, r12, r21, r22 = sp.symbols("cr11 cr12 cr21 cr22", real=True)
    L = sp.Matrix([[l11, l12], [l21, l22]])
    R = sp.Matrix([[r11, r12], [r21, r22]])
    require_zero((L * B * R).det() - L.det() * R.det() * B.det(), "group covariance")
    audits += 1

    u, a, q, n, s = sp.symbols("cu ca cq cn cs", real=True)
    C = sp.Matrix([[u, a], [q, n]])
    plus = sp.diag(s, s)
    minus = sp.diag(s, -s)
    require_zero((plus - C).det() - plus.det() - C.det() + s * (u + n), "plus gauge")
    audits += 1
    require_zero((minus - C).det() - minus.det() - C.det() + s * (-u + n), "minus gauge")
    audits += 1

    require(audits == 10, "symbolic count")
    require(payload["symbolic"]["symbolic_audits"] == audits, "stored symbolic count")
    return audits


def check_rational(payload: dict) -> int:
    J = sp.Rational(1, 2) * sp.Matrix(
        [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]]
    )
    base = rational_vector((Fraction(-7, 9), Fraction(11, 12), Fraction(5, 14), Fraction(-13, 15)))
    pieces = [
        rational_vector((Fraction(2, 5), Fraction(-3, 7), Fraction(4, 9), Fraction(5, 11))),
        rational_vector((Fraction(-5, 8), Fraction(7, 10), Fraction(-9, 13), Fraction(11, 16))),
        rational_vector((Fraction(3, 17), Fraction(4, 19), Fraction(5, 21), Fraction(-7, 23))),
    ]
    dbase = rational_vector((Fraction(3, 8), Fraction(5, 9), Fraction(-7, 12), Fraction(11, 13)))
    dpieces = [
        rational_vector((Fraction(-2, 7), Fraction(3, 10), Fraction(5, 12), Fraction(-7, 15))),
        rational_vector((Fraction(4, 11), Fraction(-6, 13), Fraction(8, 17), Fraction(9, 20))),
        rational_vector((Fraction(-5, 14), Fraction(7, 18), Fraction(-11, 22), Fraction(13, 25))),
    ]
    defect = sum(pieces, sp.zeros(4, 1))
    ddefect = sum(dpieces, sp.zeros(4, 1))
    audits = 0

    primitive = primitive_increment(base, defect, J)
    running_base = base
    sequential = 0
    for piece in pieces:
        sequential += primitive_increment(running_base, piece, J)
        running_base -= piece
    require_zero(primitive - sequential, "rational primitive")
    audits += 1

    current = current_increment(base, defect, dbase, ddefect, J)
    running_base, running_dot = base, dbase
    sequential_current = 0
    for piece, dpiece in zip(pieces, dpieces, strict=True):
        sequential_current += current_increment(running_base, piece, running_dot, dpiece, J)
        running_base -= piece
        running_dot -= dpiece
    require_zero(current - sequential_current, "rational current")
    audits += 1

    naive = sum(primitive_increment(base, piece, J) for piece in pieces)
    require(naive != primitive, "independent cross witness")
    audits += 1

    B, D = observation_matrix(base), observation_matrix(defect)
    L = sp.Matrix([[1, -3], [0, 1]])
    R = sp.Matrix([[1, 0], [2, 1]])
    transformed = (L * (B - D) * R).det() - (L * B * R).det()
    require_zero(transformed - primitive, "rational symmetry")
    audits += 1

    Bp = sp.Matrix([[3, 2], [0, 3]])
    Lp = 3 * Bp.inv()
    require_zero(Lp.det() - 1, "positive canonical determinant")
    require(Lp * Bp == 3 * sp.eye(2), "positive canonical form")
    audits += 1

    Bm = sp.Matrix([[3, 2], [0, -3]])
    Lm = 3 * Bm.inv()
    Rm = sp.diag(1, -1)
    require_zero(Lm.det() * Rm.det() - 1, "negative canonical determinant")
    require(Lm * Bm * Rm == sp.diag(3, -3), "negative canonical form")
    audits += 1

    require(audits == 6, "rational count")
    require(payload["finite_rational"]["exact_rational_audits"] == audits, "stored rational count")
    return audits


def check_rows(payload: dict) -> None:
    rows = payload["rows"]
    require(len(rows) == 20, "row count")
    require(len({row["id"] for row in rows}) == 20, "row ids")
    require(sum(row["readiness"] == "open" for row in rows) == 1, "open rows")
    require(rows[-1]["readiness"] == "guard_validated", "boundary row")
    require(all(row["claim"] and row["certificate"] and row["proof_boundary"] for row in rows), "row fields")


def main() -> int:
    require(RESULT_PATH.is_file(), "missing result")
    require(NOTE_PATH.is_file(), "missing note")
    payload = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(payload["kind"] == KIND, "kind")
    require(payload["schema_version"] == 1, "schema")
    source_audits = check_sources(payload)
    symbolic_audits = check_symbolic(payload)
    rational_audits = check_rational(payload)
    check_rows(payload)

    summary = payload["summary"]
    require(summary["rows"] == 20 and summary["issues"] == 0, "summary rows")
    require(summary["symbolic_audits"] == 10, "summary symbolic")
    require(summary["exact_rational_audits"] == 6, "summary rational")
    require(summary["source_audits"] == 3, "summary sources")
    require(summary["reciprocal_block_classes"] == 3, "summary blocks")
    require(summary["determinant_symmetries"] == 1, "summary symmetry")
    require(summary["singular_base_guards"] == 1, "summary singular")
    require(summary["open_rows"] == 1, "summary open")

    note = NOTE_PATH.read_text(encoding="utf-8")
    require(payload["success"] in note, "note success")
    require("SL(2,R) x SL(2,R)" in note, "symmetry note")
    require("diagonal" in note and "adjacent" in note and "far" in note, "block note")
    require("singular" in note and "division-free" in note, "singular note")
    require("not a proof of RH" in note, "RH boundary")
    require("prize-level conclusion" in note, "prize boundary")

    success = (
        "validated Newman C1 finite determinant reciprocal-block symmetry gate: "
        f"20 rows, 0 issues, {symbolic_audits} symbolic audits, "
        f"{rational_audits} exact-rational audits, {source_audits} source audits, "
        "3 reciprocal block classes, 1 determinant symmetry, 1 singular-base guard, "
        "1 open arithmetic row"
    )
    print(success)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
