#!/usr/bin/env python3
"""Independently check the ideal-cubic Fresnel reciprocal gate."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
    "ideal_cubic_reciprocal_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "rows": 20,
    "off_saddle_vector_identities": 2,
    "exact_weighted_mode_transport_identities": 4,
    "active_roster_reindexings": 1,
    "reciprocal_phase_identities": 6,
    "second_order_abel_identities": 1,
    "reciprocal_absolute_family_barriers": 1,
    "grouped_ideal_cubic_bounds": 0,
    "endpoint_composed_h2_bounds": 0,
    "quadratic_remainder_bounds": 0,
    "signed_flow_bounds": 0,
    "phi_b_bounds": 0,
}

EXPECTED_ROWS = [f"icrg_{index:02d}_{suffix}" for index, suffix in enumerate(
    [
        "weighted",
        "kernel",
        "weighted_transport",
        "cocycle",
        "operator",
        "first",
        "second",
        "saddle",
        "active",
        "fubini",
        "gauge",
        "critical",
        "hessian",
        "curvature",
        "self_dual",
        "coefficient",
        "abel",
        "barrier",
        "handoff",
        "boundary",
    ],
    start=1,
)]


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr | sp.MatrixBase, label: str) -> None:
    if isinstance(expression, sp.MatrixBase):
        for item in expression:
            simplified = sp.simplify(item)
            if simplified != 0:
                raise RuntimeError(f"{label}: {simplified}")
        return
    simplified = sp.simplify(expression)
    if simplified != 0:
        raise RuntimeError(f"{label}: {simplified}")


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(set(source) == {
        "aggregate_scaling",
        "basis_conjugation",
        "endpoint_abel",
        "endpoint_composition",
        "nilpotent_transport",
    }, "source-audit keys drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source for {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted for {key}")


def check_symbolics() -> None:
    i = sp.I
    x, x_0, u_x = sp.symbols("x x_0 u_x", real=True)
    t, delta, d = sp.symbols("t delta d", nonnegative=True, real=True)
    d_1, d_2 = sp.symbols("d_1 d_2", nonnegative=True, real=True)
    a_z, b_z = sp.symbols("A_z B_z", real=True)
    identity_4 = sp.eye(4)
    identity_8 = sp.eye(8)
    d_0 = sp.Matrix(
        [
            [0, 0, 0, 0],
            [0, 0, i / 2, i / 2],
            [i / 2, 0, 0, 0],
            [i / 2, 0, 0, 0],
        ]
    )
    d_hat = d_0.row_join(sp.zeros(4)).col_join(identity_4.row_join(d_0))

    def y_at(variable: sp.Expr) -> sp.Matrix:
        f = sp.Matrix(
            [
                1,
                -variable**2 / 4 - i * u_x / 2,
                i * variable / 2,
                i * variable / 2,
            ]
        )
        return f.col_join(variable * f)

    def eps(variable: sp.Expr) -> sp.Expr:
        return t * (variable + delta) / 2

    def exponent(variable: sp.Expr) -> sp.Expr:
        return t * variable**2 / 4 - (sp.Rational(1, 2) - t * delta / 2) * variable

    def e_matrix(step: sp.Expr) -> sp.Matrix:
        return sum(
            (
                (-step) ** power * d_hat**power / sp.factorial(power)
                for power in range(4)
            ),
            sp.zeros(8),
        )

    def p_matrix(step: sp.Expr, epsilon: sp.Expr) -> sp.Matrix:
        return sp.exp(
            -step * (sp.Rational(1, 2) + epsilon) + t * step**2 / 4
        ) * e_matrix(step)

    y = y_at(x)
    y_0 = y_at(x_0)
    z = sp.exp(exponent(x)) * y
    z_0 = sp.exp(exponent(x_0)) * y_0
    require_zero(sp.diff(y, x) - d_hat * y, "independent ideal jet")
    require_zero(
        sp.diff(z, x)
        - (d_hat + (-sp.Rational(1, 2) + eps(x)) * identity_8) * z,
        "independent weighted jet",
    )
    require_zero(
        p_matrix(d, eps(x)) * z - sp.exp(-d) * z.subs(x, x - d),
        "independent weighted shift",
    )

    epsilon_symbol = sp.symbols("epsilon", real=True)
    require_zero(
        p_matrix(d_2, epsilon_symbol - t * d_1 / 2)
        * p_matrix(d_1, epsilon_symbol)
        - p_matrix(d_1 + d_2, epsilon_symbol),
        "independent cocycle",
    )

    def kernel(epsilon: sp.Expr) -> sp.Matrix:
        return a_z * (
            d_hat + (-sp.Rational(1, 2) + epsilon) * identity_8
        ) + b_z * identity_8

    direct = a_z * sp.diff(z, x) + b_z * z + z_0
    require_zero(kernel(eps(x)) * z + z_0 - direct, "independent kernel")

    state = z.col_join(z_0)
    operator_0 = kernel(eps(x)).row_join(identity_8)

    def shifted_operator(step: sp.Expr) -> sp.Matrix:
        return (
            kernel(eps(x) - t * step / 2) * p_matrix(step, eps(x))
        ).row_join(p_matrix(step, eps(x_0)))

    require_zero(
        shifted_operator(d) * state
        - sp.exp(-d)
        * (
            kernel(eps(x) - t * d / 2) * z.subs(x, x - d)
            + z_0.subs(x_0, x_0 - d)
        ),
        "independent off-saddle shift",
    )
    first = shifted_operator(d) - operator_0
    second = shifted_operator(d_1 + d_2) - 2 * shifted_operator(d_1) + operator_0
    require(any(sp.simplify(item) != 0 for item in first), "first transport collapsed")
    require(any(sp.simplify(item) != 0 for item in second), "second transport collapsed")

    def saddle_matrix(epsilon: sp.Expr) -> sp.Matrix:
        return (
            d_hat**2
            - identity_8 / 12
            + 2 * epsilon * d_hat
            + (epsilon**2 + t / 2) * identity_8
        )

    saddle_0 = saddle_matrix(eps(x_0))
    saddle_1 = saddle_matrix(eps(x_0) - t * d / 2) * p_matrix(d, eps(x_0))
    require(
        any(sp.simplify(item) != 0 for item in saddle_1 - saddle_0),
        "saddle transport collapsed",
    )

    alpha, r, u, q = sp.symbols("alpha r u q", positive=True, real=True)
    psi = alpha * sp.log(u) - r * (u - q)
    critical = {r: alpha / q, u: q}
    require_zero(sp.diff(psi, r).subs(critical), "independent r critical")
    require_zero(sp.diff(psi, u).subs(critical), "independent u critical")
    require_zero(psi.subs(critical) - alpha * sp.log(q), "independent phase value")
    hessian = sp.hessian(psi, (r, u)).subs(critical)
    require_zero(hessian.det() + 1, "independent Hessian")
    require(hessian.det() < 0, "joint Hessian is not indefinite")
    require_zero(
        (-alpha / q**2) * (q**2 / alpha) + 1,
        "independent sequential curvatures",
    )
    require_zero(
        q / (2 * alpha ** sp.Rational(3, 2)) * sp.sqrt(alpha) / q
        - 1 / (2 * alpha),
        "independent reciprocal scale",
    )

    z_values = sp.symbols("a:6")
    w_values = sp.symbols("b:6")
    partial = [sum(z_values[: index + 1]) for index in range(6)]
    double_partial = [sum(partial[: index + 1]) for index in range(6)]
    lhs = sum(w_values[index] * z_values[index] for index in range(6))
    rhs = (
        w_values[5] * partial[5]
        - (w_values[5] - w_values[4]) * double_partial[4]
        + sum(
            (w_values[index + 2] - 2 * w_values[index + 1] + w_values[index])
            * double_partial[index]
            for index in range(4)
        )
    )
    require_zero(lhs - rhs, "independent second Abel identity")


def check_active_roster() -> None:
    alpha = 120
    b_value = 6
    roster = list(range(3, 241))
    for v in [
        Fraction(1, 5),
        Fraction(1, 2),
        Fraction(1, 1),
        Fraction(7, 3),
        Fraction(11, 2),
        Fraction(10, 1),
    ]:
        direct = [
            r
            for r in roster
            if Fraction(alpha) * v / b_value <= r <= Fraction(alpha) * v
        ]
        lower = ceil_fraction(Fraction(alpha) * v / b_value)
        upper = (Fraction(alpha) * v).numerator // (Fraction(alpha) * v).denominator
        interval = [r for r in roster if lower <= r <= upper]
        require(direct == interval, f"active-roster audit failed at v={v}")


def check_bounds(artifact: dict) -> None:
    values = artifact["bound_certificate"]["exact_rational_values"]
    margin = Fraction(values["core_margin"])
    require(margin == Fraction(1997, 24000), "core margin drifted")
    require(margin > Fraction(1, 13), "core margin target failed")
    exponent = Fraction(values["exponent"])
    growth = Fraction(values["growth_exponent"])
    coefficient = Fraction(values["family_coefficient"])
    require(exponent == Fraction(201, 400), "exponent drifted")
    require(growth == Fraction(199, 400), "growth exponent drifted")
    require(coefficient == Fraction(200, 2587), "family coefficient drifted")
    require(coefficient == Fraction(1, 26) / growth, "integral coefficient failed")


def check_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_ROWS, "row order drifted")
    for row in rows:
        for field in ["role", "readiness", "claim", "certificate", "proof_boundary"]:
            require(bool(row.get(field)), f"empty {field} in {row.get('id')}")
    require(
        [row["id"] for row in rows if row["readiness"] == "open"]
        == ["icrg_19_handoff"],
        "open-row boundary drifted",
    )
    barrier = next(row for row in rows if row["id"] == "icrg_18_barrier")
    require("does not lower-bound" in barrier["proof_boundary"], "barrier overclaim")


def check_note(note: str) -> None:
    required = [
        "This is not a proof",
        "## Weighted Off-Saddle Jet",
        "## Exact Roster Transport",
        "## Active Roster",
        "## Reciprocal Self-Duality",
        "determinant -1 is exact",
        "## Second-Order Abel And Barrier",
        "lower bound only on that chosen nonnegative scale budget",
        "## Proof Boundary",
    ]
    for marker in required:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "count block drifted")
    require("no grouped h2" in artifact.get("status", ""), "status boundary drifted")
    require("It proves no grouped" in artifact.get("proof_boundary", ""), "proof boundary drifted")
    check_sources(artifact)
    check_symbolics()
    check_active_roster()
    check_bounds(artifact)
    check_rows(artifact)
    check_note(note)
    counts = artifact["counts"]
    print(
        "validated ideal-cubic Fresnel reciprocal gate: "
        f"{counts['rows']} rows, "
        f"{counts['exact_weighted_mode_transport_identities']} weighted transports, "
        f"{counts['active_roster_reindexings']} Fubini reindexing, "
        f"{counts['reciprocal_phase_identities']} reciprocal phase identities, "
        f"{counts['reciprocal_absolute_family_barriers']} absolute-family barrier, "
        f"{counts['grouped_ideal_cubic_bounds']} grouped bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
