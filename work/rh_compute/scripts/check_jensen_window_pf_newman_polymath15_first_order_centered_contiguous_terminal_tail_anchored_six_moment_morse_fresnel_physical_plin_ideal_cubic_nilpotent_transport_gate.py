#!/usr/bin/env python3
"""Independently check the ideal-cubic nilpotent transport gate."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "physical_plin_ideal_cubic_nilpotent_transport_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"


def fail(message: str) -> None:
    raise RuntimeError(message)


def require_zero(expression: sp.Expr | sp.MatrixBase, label: str) -> None:
    if isinstance(expression, sp.MatrixBase):
        if any(sp.simplify(sp.expand(item)) != 0 for item in expression):
            fail(label)
        return
    if sp.simplify(sp.expand(expression)) != 0:
        fail(label)


def matrix_one_norm(matrix: sp.MatrixBase) -> sp.Expr:
    return max(
        sp.simplify(sum(abs(matrix[i, j]) for i in range(matrix.rows)))
        for j in range(matrix.cols)
    )


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/", 1)
    return Fraction(int(numerator), int(denominator))


def check_sources(artifact: dict) -> None:
    audit = artifact.get("source_audit", {})
    if set(audit) != {
        "basis_conjugation",
        "endpoint_abel",
        "q1_saddle",
        "physical_plin",
    }:
        fail("source audit roster drifted")
    for key, source in audit.items():
        path = REPO_ROOT / source.get("path", "")
        if not path.is_file():
            fail(f"missing source: {key}")
        actual = sha256(path.read_bytes()).hexdigest()
        if actual != source.get("sha256"):
            fail(f"source hash drifted: {key}")


def build_matrices() -> tuple[sp.Matrix, sp.Matrix]:
    i = sp.I
    d_0 = sp.Matrix(
        [
            [0, 0, 0, 0],
            [0, 0, i / 2, i / 2],
            [i / 2, 0, 0, 0],
            [i / 2, 0, 0, 0],
        ]
    )
    d_hat = d_0.row_join(sp.zeros(4)).col_join(sp.eye(4).row_join(d_0))
    return d_0, d_hat


def check_jets_and_cubic() -> None:
    i = sp.I
    x, u_x = sp.symbols("x u_x", real=True)
    d_0, d_hat = build_matrices()
    field = sp.Matrix([1, -x**2 / 4 - i * u_x / 2, i * x / 2, i * x / 2])
    augmented = field.col_join(x * field)
    require_zero(d_0**3, "independent D_0 nilpotence failed")
    if d_0**2 == sp.zeros(4):
        fail("D_0 nilpotence order collapsed")
    require_zero(sp.diff(field, x) - d_0 * field, "independent field jet")
    require_zero(d_hat**4, "independent Dhat nilpotence failed")
    if d_hat**3 == sp.zeros(8):
        fail("Dhat nilpotence order collapsed")
    require_zero(sp.diff(augmented, x) - d_hat * augmented, "independent joined jet")
    if [(d_hat**k).rank() for k in range(1, 4)] != [5, 2, 1]:
        fail("augmented power ranks drifted")

    alpha = sp.symbols("alpha_V alpha_N alpha_A alpha_Q", real=True)
    beta = sp.symbols("beta_V beta_N beta_A beta_Q", real=True)
    row = sp.Matrix([[*alpha, *[i * value for value in beta]]])
    polynomial = sp.expand((row * augmented)[0])
    coefficients = [
        sp.Poly(polynomial, x).coeff_monomial(x**j) for j in range(4)
    ]
    expected = [
        alpha[0] - i * u_x * alpha[1] / 2,
        i * (alpha[2] + alpha[3]) / 2 + i * beta[0] + u_x * beta[1] / 2,
        -alpha[1] / 4 - (beta[2] + beta[3]) / 2,
        -i * beta[1] / 4,
    ]
    for index in range(4):
        require_zero(coefficients[index] - expected[index], f"cubic row {index}")


def check_saddle_matrix() -> None:
    x, t, delta_a = sp.symbols("x t delta_a", real=True)
    _, d_hat = build_matrices()
    identity = sp.eye(8)
    epsilon = t * (x + delta_a) / 2
    g = -sp.Rational(1, 2) + epsilon
    field_i = sp.I
    u_x = sp.symbols("u_x", real=True)
    field = sp.Matrix(
        [1, -x**2 / 4 - field_i * u_x / 2, field_i * x / 2, field_i * x / 2]
    )
    augmented = field.col_join(x * field)
    matrix = (
        d_hat**2
        - identity / 12
        + 2 * epsilon * d_hat
        + (epsilon**2 + t / 2) * identity
    )
    h_scalar = g**2 + g + t / 2 + sp.Rational(1, 6)
    direct = (
        sp.diff(augmented, x, 2)
        + (2 * g + 1) * sp.diff(augmented, x)
        + h_scalar * augmented
    )
    require_zero(direct - matrix * augmented, "independent saddle matrix")
    core = d_hat**2 - identity / 12
    if sp.factor(core.det()) != sp.Rational(1, 12**8):
        fail("independent core determinant failed")
    norms = [matrix_one_norm(d_hat**k) for k in range(1, 4)]
    if norms != [2, sp.Rational(5, 2), sp.Rational(3, 2)]:
        fail("independent augmented norms drifted")


def check_transport() -> None:
    x, d_1, d_2 = sp.symbols("x d_1 d_2", real=True)
    u_x = sp.symbols("u_x", real=True)
    i = sp.I
    _, d_hat = build_matrices()
    identity = sp.eye(8)
    field = sp.Matrix([1, -x**2 / 4 - i * u_x / 2, i * x / 2, i * x / 2])
    augmented = field.col_join(x * field)

    def transport(step: sp.Expr) -> sp.Matrix:
        return sum(
            (
                ((-step) ** k / sp.factorial(k)) * d_hat**k
                for k in range(4)
            ),
            sp.zeros(8),
        )

    require_zero(
        augmented.subs(x, x - d_1) - transport(d_1) * augmented,
        "independent one-step transport",
    )
    require_zero(
        transport(d_2) * transport(d_1) - transport(d_1 + d_2),
        "independent transport semigroup",
    )
    second = transport(d_1 + d_2) - 2 * transport(d_1) + identity
    expected = (
        (d_1 - d_2) * d_hat
        + ((d_1 + d_2) ** 2 / 2 - d_1**2) * d_hat**2
        + (2 * d_1**3 - (d_1 + d_2) ** 3) * d_hat**3 / 6
    )
    require_zero(second - expected, "independent second transport")

    t, epsilon = sp.symbols("t epsilon", real=True)

    def mode_matrix(e_value: sp.Expr) -> sp.Matrix:
        return (
            d_hat**2
            - identity / 12
            + 2 * e_value * d_hat
            + (e_value**2 + t / 2) * identity
        )

    first_saddle = mode_matrix(epsilon - t * d_1 / 2) * transport(
        d_1
    ) - mode_matrix(epsilon)
    second_saddle = (
        mode_matrix(epsilon - t * (d_1 + d_2) / 2)
        * transport(d_1 + d_2)
        - 2
        * mode_matrix(epsilon - t * d_1 / 2)
        * transport(d_1)
        + mode_matrix(epsilon)
    )
    if not any(item != 0 for item in first_saddle):
        fail("physical first saddle transport collapsed")
    if not any(item != 0 for item in second_saddle):
        fail("physical second saddle transport collapsed")


def check_bounds(artifact: dict) -> None:
    alpha_min = 4096
    first = (
        Fraction(20, 9)
        + Fraction(5, 4) * Fraction(100, 81 * alpha_min)
        + Fraction(1, 4) * Fraction(1000, 729 * alpha_min**2)
    )
    second = (
        Fraction(5, 2)
        + Fraction(500, 81)
        + Fraction(5, 2) * Fraction(1000, 729 * alpha_min)
    )
    if not first < Fraction(9, 4):
        fail("independent first collar bound failed")
    if not second < 9:
        fail("independent second collar bound failed")
    perturbation = Fraction(161, 8000)
    diagnostics = artifact.get("bound_certificate", {}).get(
        "rational_diagnostics", {}
    )
    expected = {
        "first_difference_scaled_at_alpha_4096": first,
        "second_difference_scaled_at_alpha_4096": second,
        "operator_perturbation": perturbation,
    }
    for key, value in expected.items():
        if parse_fraction(diagnostics.get(key, "0/1")) != value:
            fail(f"artifact rational diagnostic drifted: {key}")


def check_artifact_and_note(artifact: dict, note: str) -> None:
    if artifact.get("kind") != STEM:
        fail("artifact kind drifted")
    expected_counts = {
        "rows": 18,
        "four_vector_nilpotence_order": 3,
        "augmented_nilpotence_order": 4,
        "ideal_cubic_coefficients": 4,
        "invisible_ideal_row_directions": 2,
        "exact_mode_transport_identities": 5,
        "collar_variation_bounds": 2,
        "grouped_ideal_cubic_bounds": 0,
        "numerical_observation_row_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_flow_bounds": 0,
        "phi_b_bounds": 0,
    }
    if artifact.get("counts") != expected_counts:
        fail("artifact counts drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 18 or len({row.get("id") for row in rows}) != 18:
        fail("artifact rows drifted")
    tokens = (
        "F_0(x)=(1,-x^2/4-i*u_(N,x)/2,i*x/2,i*x/2)^T",
        "Dhat^4=0 but Dhat^3!=0",
        "det(M_0)=12^(-8)",
        "d_r-d_(r+1)=Delta^2G/alpha",
        "(9/alpha^2)||Y_r||_1",
        "This is not a proof",
        "grouped ideal-cubic c-prime or roster estimate",
    )
    for token in tokens:
        if token not in note:
            fail(f"note token missing: {token}")


def main() -> int:
    if not RESULT_PATH.is_file() or not NOTE_PATH.is_file():
        fail("nilpotent transport result or note missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    check_sources(artifact)
    check_jets_and_cubic()
    check_saddle_matrix()
    check_transport()
    check_bounds(artifact)
    check_artifact_and_note(artifact, note)
    counts = artifact["counts"]
    print(
        "validated ideal-cubic nilpotent transport gate: "
        f"{counts['rows']} rows, nilpotence orders "
        f"{counts['four_vector_nilpotence_order']}/"
        f"{counts['augmented_nilpotence_order']}, "
        f"{counts['ideal_cubic_coefficients']} cubic coefficients, "
        f"{counts['exact_mode_transport_identities']} transport identities, "
        f"{counts['collar_variation_bounds']} collar bounds, "
        f"{counts['grouped_ideal_cubic_bounds']} grouped bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
