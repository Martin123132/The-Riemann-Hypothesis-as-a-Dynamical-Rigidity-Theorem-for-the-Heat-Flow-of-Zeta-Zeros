#!/usr/bin/env python3
"""Validate the explicit five-moment real quadratic reduction for Phi_B."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPTS = REPO_ROOT / "work/rh_compute/scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import jensen_window_pf_newman_polymath15_first_order_centered_contiguous_terminal_tail_anchored_five_moment_explicit_phi_quadratic_reduction as gate  # noqa: E402


STEM = gate.STEM
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    f"ttqp_{index:02d}_{suffix}"
    for index, suffix in enumerate(
        (
            "domain",
            "base",
            "derivatives",
            "observations",
            "phi",
            "quadratic",
            "rank",
            "inertia",
            "kernel",
            "value_fibre",
            "contact",
            "target",
            "handoff",
        ),
        start=1,
    )
]
EXPECTED_COUNTS = {
    "rows": 13,
    "complex_correction_free_moments": 5,
    "real_moment_coordinates": 10,
    "complex_base_coefficient_vectors": 5,
    "real_bulk_observations": 4,
    "explicit_phi_quadratic_forms": 1,
    "quadratic_rank_upper_bound": 4,
    "generic_positive_directions": 2,
    "generic_negative_directions": 2,
    "generic_null_directions": 6,
    "division_free_value_fibres": 1,
    "double_contact_obstructions": 1,
    "signed_phi_targets": 1,
    "signed_phi_bounds": 0,
    "xi_level_current_theorems": 0,
}
SUMMARY = (
    "validated explicit five-moment Phi_B quadratic reduction: "
    "13 rows, 4 real observations, rank <=4, generic inertia (2,2,6), "
    "1 double-contact obstruction, 0 signed Phi_B bounds"
)


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_json(path: Path, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing result: {path}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid result: {exc}")
        return {}


def real_row(vector: sp.Matrix) -> sp.Matrix:
    real = [sp.expand(entry).as_real_imag()[0] for entry in vector]
    imag = [sp.expand(entry).as_real_imag()[1] for entry in vector]
    return sp.Matrix(real + [-entry for entry in imag])


def independent_coefficient_audit(issues: list[str]) -> None:
    l_n = sp.symbols("ell", real=True)
    r_1, r_2, r_1_x, r_2_x = sp.symbols("r_1 r_2 r_1_x r_2_x")
    g = sp.Matrix(sp.symbols("g_0:5"))

    h_0 = g[0] + r_1 * g[1] + r_2 * g[2]
    h_1 = g[1] + r_1 * g[2] + r_2 * g[3]
    h_2 = g[2] + r_1 * g[3] + r_2 * g[4]
    d_0 = r_1_x * g[1] + r_2_x * g[2]
    d_1 = r_1_x * g[2] + r_2_x * g[3]

    expected = (
        l_n * h_0 - h_1,
        l_n**2 * h_0 - 2 * l_n * h_1 + h_2,
        l_n * d_0 - d_1,
    )
    vectors = (
        sp.Matrix(
            [l_n, l_n * r_1 - 1, l_n * r_2 - r_1, -r_2, 0]
        ),
        sp.Matrix(
            [
                l_n**2,
                l_n**2 * r_1 - 2 * l_n,
                l_n**2 * r_2 - 2 * l_n * r_1 + 1,
                r_1 - 2 * l_n * r_2,
                r_2,
            ]
        ),
        sp.Matrix(
            [0, l_n * r_1_x, l_n * r_2_x - r_1_x, -r_2_x, 0]
        ),
    )
    for index, (target, vector) in enumerate(zip(expected, vectors)):
        if sp.expand(target - (vector.T * g)[0]) != 0:
            issues.append(f"independent base vector failed: {index}")


def independent_quadratic_audit(issues: list[str]) -> None:
    # Use four independent generic real rows. This verifies the exact
    # realification and determinant polarization without importing the gate's
    # symbolic construction.
    entries = sp.symbols("u_0:40", real=True)
    v = sp.Matrix(entries[0:10])
    q = sp.Matrix(entries[10:20])
    a = sp.Matrix(entries[20:30])
    n = sp.Matrix(entries[30:40])
    g = sp.Matrix(sp.symbols("x_0:10", real=True))
    x_t, a_t, x_t_x, a_t_x = sp.symbols(
        "x_t a_t x_t_x a_t_x", real=True
    )

    value = (v.T * g)[0]
    value_x = (q.T * g)[0]
    scalar = (a.T * g)[0]
    scalar_x = (n.T * g)[0]
    tail = x_t * a_t_x - a_t * x_t_x
    total = (
        (x_t + value) * (a_t_x + scalar_x)
        - (a_t + scalar) * (x_t_x + value_x)
    )
    phi = sp.expand(total - tail)

    matrix = (
        v * n.T + n * v.T - a * q.T - q * a.T
    ) / 2
    linear = x_t * n + a_t_x * v - a_t * q - x_t_x * a
    reduced = sp.expand((g.T * matrix * g)[0] + (linear.T * g)[0])
    if sp.expand(phi - reduced) != 0:
        issues.append("independent real quadratic identity failed")
    if matrix != matrix.T:
        issues.append("independent quadratic symmetry failed")

    factor = sp.Matrix.hstack(v, n, a, q)
    core = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    if matrix != factor * core * factor.T:
        issues.append("independent rank factorization failed")

    contact = sp.expand(phi.subs({x_t: -value, a_t: -scalar}))
    negative_tail = sp.expand(
        -tail.subs({x_t: -value, a_t: -scalar})
    )
    if sp.expand(contact - negative_tail) != 0:
        issues.append("independent double-contact identity failed")


def independent_rank_witness(issues: list[str]) -> None:
    observation = sp.Matrix(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 1, -sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, 0, 0, sp.Rational(3, 2), -sp.Rational(1, 2), 0, 0, 0],
            [
                -sp.Rational(3, 2),
                sp.Rational(5, 4),
                -sp.Rational(1, 4),
                0,
                0,
                2,
                -sp.Rational(1, 2),
                0,
                0,
                0,
            ],
        ]
    )
    if observation.rank() != 4:
        issues.append("independent witness observation rank failed")
    if observation[:, [0, 1, 5, 6]].det() != sp.Rational(5, 16):
        issues.append("independent witness minor failed")

    # A full-column factor preserves the inertia of the nonsingular core on
    # its four-dimensional image. The two off-diagonal blocks each have one
    # positive and one negative eigenvalue.
    core = sp.Matrix(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, -1],
            [0, 0, -1, 0],
        ]
    )
    if core.eigenvals() != {sp.Integer(-1): 2, sp.Integer(1): 2}:
        issues.append("independent core inertia failed")


def validate(payload: dict, issues: list[str]) -> None:
    if payload.get("kind") != STEM:
        issues.append("unexpected artifact kind")
    if payload.get("date") != "2026-07-31":
        issues.append("artifact date drifted")
    if payload.get("counts") != EXPECTED_COUNTS:
        issues.append("counts drifted")
    rows = payload.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or order drifted")
    if sum(row.get("readiness") == "open" for row in rows) != 2:
        issues.append("open-row count drifted")

    status = payload.get("status", "")
    for token in (
        "exact terminal-tail five-moment",
        "signed Phi_B bound open",
        "no retained aggregate sign",
        "or RH",
    ):
        if token not in status:
            issues.append(f"status missing token: {token}")
    boundary = payload.get("proof_boundary", "")
    for token in (
        "no upper bound for Phi_B",
        "contact exclusion",
        "Lambda<=0",
        "RH",
    ):
        if token not in boundary:
            issues.append(f"proof boundary missing token: {token}")

    hashes = payload.get("source_audit", {}).get("source_sha256", {})
    for key, path in gate.SOURCE_PATHS.items():
        if not path.is_file():
            issues.append(f"missing source: {path}")
        elif hashes.get(key) != file_hash(path):
            issues.append(f"source hash drifted: {key}")

    try:
        sources = gate.load_sources()
        if payload.get("source_audit") != gate.source_audit(sources):
            issues.append("source audit recomputation drifted")
        if payload.get("symbolic_certificate") != gate.symbolic_certificate():
            issues.append("symbolic certificate drifted")
        if payload.get("exact") != gate.exact_payload():
            issues.append("exact payload drifted")
    except Exception as exc:  # noqa: BLE001
        issues.append(f"certificate recomputation failed: {exc}")

    if payload.get("counts", {}).get("signed_phi_bounds") != 0:
        issues.append("unexpected signed Phi_B bound")

    independent_coefficient_audit(issues)
    independent_quadratic_audit(issues)
    independent_rank_witness(issues)


def validate_note(issues: list[str]) -> None:
    if not NOTE.is_file():
        issues.append(f"missing note: {NOTE}")
        return
    text = NOTE.read_text(encoding="utf-8")
    for token in (
        "Terminal-Tail Explicit Five-Moment Phi Quadratic Reduction",
        "Date: 2026-07-31",
        "0 Phi_B bounds",
        "h^2*Phi_B=g^T M g+ell_T^T g",
        "rank(M)<=4",
        "inertia (2,2,6)",
        "Phi_B=-P_T>1/400",
        "<=h^2/400",
    ):
        if token not in text:
            issues.append(f"note marker missing: {token}")


def main() -> int:
    issues: list[str] = []
    payload = load_json(RESULT, issues)
    if payload:
        validate(payload, issues)
    validate_note(issues)
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1
    print(SUMMARY)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
