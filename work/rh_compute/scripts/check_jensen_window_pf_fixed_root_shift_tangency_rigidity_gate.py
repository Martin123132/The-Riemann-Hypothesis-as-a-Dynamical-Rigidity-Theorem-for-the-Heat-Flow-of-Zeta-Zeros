#!/usr/bin/env python3
"""Independently check the fixed-root shift-tangency rigidity gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_fixed_root_shift_tangency_rigidity_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 32,
    "derivative_checks": 105,
    "adjacent_shift_checks": 24,
    "same_root_transfer_checks": 40,
    "classification_checks": 249,
    "recurrence_checks": 108,
    "persistent_witness_checks": 67,
    "hankel_checks": 25,
    "m100_applications": 1,
    "heat_interval_applications": 1,
    "newton_identity_checks": 55,
    "newton_bound_checks": 55,
    "root_drift_factor_checks": 21,
    "zeta_near_chain_bounds": 0,
    "rh_conclusions": 0,
}

EXPECTED_ROWS = [
    "frt_01_sources", "frt_02_derivative", "frt_03_adjacent", "frt_04_transfer",
    "frt_05_polar_dependency", "frt_06_zero_root", "frt_07_recurrence",
    "frt_08_indexing", "frt_09_tail", "frt_10_differential",
    "frt_11_classification", "frt_12_converse", "frt_13_ogf",
    "frt_14_normalization", "frt_15_persistent", "frt_16_persistent_root",
    "frt_17_persistent_ogf", "frt_18_persistent_egf", "frt_19_hierarchy",
    "frt_20_source_scope", "frt_21_hankel_vector", "frt_22_rank",
    "frt_23_quartic_hankel", "frt_24_heat_interval", "frt_25_conditional_transfer",
    "frt_26_normalized_defect", "frt_27_newton", "frt_28_near_bound",
    "frt_29_hankel_margin", "frt_30_root_drift", "frt_31_zeta_near",
    "frt_32_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_sources(artifact: dict) -> None:
    expected = {
        "heat_hierarchy", "polar_contact", "coefficient_pf", "degree_cascade",
        "order4_m100", "order4_forward",
    }
    source = artifact.get("source_audit", {})
    require(set(source) == expected, "source keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def jensen(values: list[sp.Expr], degree: int, shift: int, z: sp.Symbol) -> sp.Expr:
    return sp.expand(sum(sp.binomial(degree, j) * values[shift + j] * z**j for j in range(degree + 1)))


def check_derivative_and_adjacent() -> tuple[int, int]:
    z = sp.symbols("z")
    values = [sp.Rational(5 * k**2 + 2 * k + 9, 3 * k + 4) for k in range(24)]
    derivative_checks = 0
    for degree in range(2, 9):
        for shift in (0, 2):
            polynomial = jensen(values, degree, shift, z)
            for order in range(degree + 1):
                expected = sp.factorial(degree) / sp.factorial(degree - order) * sum(
                    sp.binomial(degree - order, j) * z**j * values[shift + order + j]
                    for j in range(degree - order + 1)
                )
                require_zero(sp.diff(polynomial, z, order) - expected, "independent derivative")
                derivative_checks += 1

    adjacent_checks = 0
    for degree in range(1, 8):
        for shift in (0, 1, 3):
            require_zero(
                jensen(values, degree + 1, shift, z)
                - jensen(values, degree, shift, z)
                - z * jensen(values, degree, shift + 1, z),
                "independent adjacent identity",
            )
            adjacent_checks += 1
    return derivative_checks, adjacent_checks


def check_transfer() -> int:
    z = sp.symbols("z")
    checks = 0
    for root in (sp.Rational(-5, 2), sp.Rational(-7, 3)):
        for multiplicity in range(1, 7):
            lower = (z - root) ** multiplicity * (-root + 3 * z + z**3)
            upper = (z - root) ** (multiplicity + 1) * (1 - 2 * z + z**2)
            shifted, remainder = sp.div(sp.expand(upper - lower), z)
            require_zero(remainder, "independent transfer divisibility")
            for order in range(multiplicity):
                require_zero(sp.diff(shifted, z, order).subs(z, root), "independent transfer jet")
                checks += 1
    return checks


def check_persistent() -> tuple[int, int]:
    n, z = sp.symbols("n z", integer=True, nonnegative=True)
    root = -10
    value = lambda k: (33 * k**2 + 3 * k + 4) / (
        sp.Integer(4) * sp.Integer(10) ** k
    )
    recurrence = value(n) + 3 * root * value(n + 1) + 3 * root**2 * value(n + 2) + root**3 * value(n + 3)
    require_zero(recurrence, "independent persistent recurrence")
    root_checks = 0
    for shift in range(8):
        polynomial = sum(sp.binomial(4, j) * value(shift + j) * z**j for j in range(5))
        require_zero(polynomial.subs(z, root), "independent persistent root")
        require_zero(sp.diff(polynomial, z).subs(z, root), "independent persistent derivative")
        require(sp.diff(polynomial, z, 2).subs(z, root) != 0, "independent persistent multiplicity")
        root_checks += 3

    egf = (1 + sp.Rational(9, 10) * z + sp.Rational(33, 400) * z**2) * sp.exp(z / 10)
    for index in range(9):
        require_zero(sp.diff(egf, z, index).subs(z, 0) - value(index), "independent EGF")
        root_checks += 1

    hankel_checks = 0
    c = sp.Matrix([1, 3 * root, 3 * root**2, root**3])
    for shift in range(5):
        matrix = sp.Matrix(4, 4, lambda i, j: value(shift + 1 + i + j))
        require(matrix * c == sp.zeros(4, 1), "independent quartic Hankel kernel")
        require_zero(matrix.det(), "independent quartic determinant")
        hankel_checks += 2
    return root_checks, hankel_checks


def finite_difference(values: list[sp.Expr], order: int, index: int) -> sp.Expr:
    return sp.simplify(sum((-1) ** (order - j) * sp.binomial(order, j) * values[index + j] for j in range(order + 1)))


def check_newton() -> int:
    values = [sp.Rational(2 * k**5 + 7 * k**2 - 3 * k + 11, k + 3) for k in range(16)]
    checks = 0
    for order in range(1, 7):
        for k in range(order, 16):
            base = sum(sp.binomial(k, s) * finite_difference(values, s, 0) for s in range(order))
            remainder = sum(
                sp.binomial(k - 1 - j, order - 1) * finite_difference(values, order, j)
                for j in range(k - order + 1)
            )
            require_zero(values[k] - base - remainder, "independent Newton identity")
            checks += 1
    return checks


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["frt_30_root_drift", "frt_31_zeta_near"],
        "open rows drifted",
    )


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH", "## Adjacent Transfer", "J_A^(d+1,n)",
        "## Fixed-Root Classification", "exp(qz)Q", "A_(n+1)+3rA_(n+2)",
        "## Persistent Rational Witness", "33/2", "exp(z/10)",
        "## Finite Hankel Bridge", "H_(4,n+1)", "lambda=-100", "[-100,0]",
        "## Quantitative Near-Rigidity", "binom(k,3)", "sigma_min",
        "## Pi Provenance", "## Proof Boundary",
    ]
    for marker in markers:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "counts drifted")
    require("finite Hankel heat-interval obstruction" in artifact.get("status", ""), "status drifted")
    require("no heat-uniform drifting-root exclusion" in artifact.get("proof_boundary", ""), "boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {"transfer", "fixed_root", "persistent", "hankel", "near_rigidity", "handoff"},
        "certificate keys drifted",
    )
    check_sources(artifact)
    derivative_checks, adjacent_checks = check_derivative_and_adjacent()
    transfer_checks = check_transfer()
    persistent_checks, hankel_checks = check_persistent()
    newton_checks = check_newton()
    check_rows(artifact)
    check_note(note)
    print(
        "validated fixed-root shift-tangency rigidity gate: "
        f"{EXPECTED_COUNTS['rows']} rows, {derivative_checks} independent derivatives, "
        f"{adjacent_checks} adjacent shifts, {transfer_checks} transfer jets, "
        f"{persistent_checks} persistent checks, {hankel_checks} Hankel checks, "
        f"{newton_checks} Newton identities, "
        f"{EXPECTED_COUNTS['m100_applications']} lambda=-100 anchor, "
        f"{EXPECTED_COUNTS['heat_interval_applications']} heat-interval application, "
        f"{EXPECTED_COUNTS['zeta_near_chain_bounds']} zeta near-chain bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
