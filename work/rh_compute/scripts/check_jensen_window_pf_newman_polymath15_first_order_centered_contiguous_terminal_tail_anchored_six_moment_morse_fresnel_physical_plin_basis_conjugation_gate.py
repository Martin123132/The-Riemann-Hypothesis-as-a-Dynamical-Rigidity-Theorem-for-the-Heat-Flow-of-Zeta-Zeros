#!/usr/bin/env python3
"""Independently check the centered-basis Morse conjugation gate."""

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
    "physical_plin_basis_conjugation_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{STEM}.md"
K_VALUES = (6, 14, 55, 336, 2738, 27936)
H_DENOMINATOR = 72_000_000_000


def fail(message: str) -> None:
    raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.simplify(sp.expand(expression)) != 0:
        fail(label)


def parse_fraction(text: str) -> Fraction:
    numerator, denominator = text.split("/", 1)
    return Fraction(int(numerator), int(denominator))


def check_source_hashes(artifact: dict) -> None:
    audit = artifact.get("source_audit", {})
    if set(audit) != {
        "critical_ray",
        "q1_saddle",
        "observation_image",
        "physical_plin",
        "joined_envelope",
    }:
        fail("source audit keys drifted")
    for key, source in audit.items():
        path = REPO_ROOT / source.get("path", "")
        if not path.is_file():
            fail(f"missing source for {key}: {path}")
        actual = sha256(path.read_bytes()).hexdigest()
        if actual != source.get("sha256"):
            fail(f"source hash drifted for {key}")


def check_basis_and_weights(artifact: dict) -> None:
    lam, ell, x = sp.symbols("lambda ell x", real=True)
    p = sp.symbols("p_0:6")
    q = [
        sp.expand(
            sum(
                sp.binomial(k, j) * (-ell) ** (k - j) * p[k]
                for k in range(j, 6)
            )
        )
        for j in range(6)
    ]
    centered = sum(p[k] * x**k for k in range(6))
    uncentered = sum(q[j] * lam**j for j in range(6))
    require_zero(
        uncentered.subs(lam, ell + x) - centered,
        "independent forward basis check failed",
    )

    inverse = [
        sp.expand(
            sum(
                sp.binomial(k, j) * ell ** (k - j) * q[k]
                for k in range(j, 6)
            )
        )
        for j in range(6)
    ]
    for j in range(6):
        require_zero(inverse[j] - p[j], f"independent inverse row {j}")

    weights = [
        sp.expand(
            sum(
                sp.binomial(k, j) * K_VALUES[j] * ell ** (k - j)
                for j in range(k + 1)
            )
        )
        for k in range(6)
    ]
    expected = [
        6,
        6 * ell + 14,
        6 * ell**2 + 28 * ell + 55,
        6 * ell**3 + 42 * ell**2 + 165 * ell + 336,
        6 * ell**4 + 56 * ell**3 + 330 * ell**2 + 1344 * ell + 2738,
        6 * ell**5
        + 70 * ell**4
        + 550 * ell**3
        + 3360 * ell**2
        + 13690 * ell
        + 27936,
    ]
    for j in range(6):
        require_zero(weights[j] - expected[j], f"translated weight row {j}")

    rows = artifact["symbolic_certificate"]["translated_tail_weights"]["rows"]
    if len(rows) != 6:
        fail("artifact translated-weight count drifted")
    for j, expression in enumerate(weights):
        if sp.sstr(expression) not in rows[j]:
            fail(f"artifact translated weight {j} drifted")


def check_physical_conjugation() -> None:
    x, ell, lam = sp.symbols("x ell lambda", real=True)
    t, sigma, b_a = sp.symbols("t sigma B_a", real=True)
    p = sp.symbols("p_0:6")
    polynomial = sum(p[k] * x**k for k in range(6))
    S = t * lam**2 / 4 - sigma * lam
    centered_S = t * x**2 / 4 - b_a * x
    physical_b = sigma - t * ell / 2
    require_zero(
        (S.subs(lam, ell + x) - S.subs(lam, ell) - centered_S).subs(
            b_a, physical_b
        ),
        "independent exponent conjugation failed",
    )

    alpha_real, delta_a = sp.symbols("alpha_real delta_a", real=True)
    sigma_physical = sp.Rational(1, 2) + t * alpha_real / 2
    b_physical = sp.Rational(1, 2) - t * delta_a / 2
    require_zero(
        (b_physical - (sigma_physical - t * ell / 2)).subs(
            delta_a, ell - alpha_real
        ),
        "independent B_a source identity failed",
    )

    g = t * x / 2 - b_a
    operator = (
        sp.diff(polynomial, x, 2)
        + (2 * g + 1) * sp.diff(polynomial, x)
        + (g**2 + g + t / 2 + sp.Rational(1, 6)) * polynomial
    )
    conjugated = sp.exp(-centered_S) * (
        sp.diff(sp.exp(centered_S) * polynomial, x, 2)
        + sp.diff(sp.exp(centered_S) * polynomial, x)
        + sp.exp(centered_S) * polynomial / 6
    )
    require_zero(conjugated - operator, "independent saddle conjugation failed")


def check_joined_rows() -> None:
    x, t, b_a = sp.symbols("x t B_a", real=True)
    g = t * x / 2 - b_a
    alpha = sp.symbols("alpha_0:4", real=True)
    beta = sp.symbols("beta_0:4", real=True)
    fields = [sp.Function(f"G_{j}")(x) for j in range(4)]
    f_alpha = sum(alpha[j] * fields[j] for j in range(4))
    f_beta = sum(beta[j] * fields[j] for j in range(4))
    joined = f_alpha + sp.I * x * f_beta

    def first(expression: sp.Expr) -> sp.Expr:
        return sp.diff(expression, x) + g * expression

    def second(expression: sp.Expr) -> sp.Expr:
        return (
            sp.diff(expression, x, 2)
            + (2 * g + 1) * sp.diff(expression, x)
            + (g**2 + g + t / 2 + sp.Rational(1, 6)) * expression
        )

    first_expected = sum(
        alpha[j] * first(fields[j])
        + sp.I * beta[j] * (x * first(fields[j]) + fields[j])
        for j in range(4)
    )
    require_zero(first(joined) - first_expected, "joined first-order identity")

    second_expected = sum(
        alpha[j] * second(fields[j])
        + sp.I
        * beta[j]
        * (
            x * second(fields[j])
            + 2 * sp.diff(fields[j], x)
            + (2 * g + 1) * fields[j]
        )
        for j in range(4)
    )
    require_zero(second(joined) - second_expected, "joined saddle identity")

    a_z, b_z, e_x, e_0, x_0 = sp.symbols("A_z B_z E_x E_0 x_0")
    off = e_x * (a_z * first(joined) + b_z * joined) + e_0 * joined.subs(
        x, x_0
    )
    expected = 0
    for j, field in enumerate(fields):
        field_0 = field.subs(x, x_0)
        k_alpha = e_x * (a_z * first(field) + b_z * field) + e_0 * field_0
        k_beta = e_x * (
            a_z * (x * first(field) + field) + b_z * x * field
        ) + e_0 * x_0 * field_0
        expected += alpha[j] * k_alpha + sp.I * beta[j] * k_beta
    require_zero(off - expected, "joined off-saddle identity")


def check_tail_rounding(artifact: dict) -> None:
    ell = sp.symbols("ell", nonnegative=True)
    weights = [
        sp.Poly(
            sum(
                sp.binomial(k, j) * K_VALUES[j] * ell ** (k - j)
                for j in range(k + 1)
            ),
            ell,
        )
        for k in range(6)
    ]
    at_50 = [int(poly.eval(50)) for poly in weights]
    if at_50 != [6, 314, 16455, 863586, 45394938, 2390362436]:
        fail("independent L=50 weights drifted")

    c_50 = (
        8764 * at_50[0]
        + 8765 * at_50[1]
        + 8762 * at_50[2]
        + 8760 * at_50[3]
        + Fraction(2, 3) * at_50[4]
        + Fraction(1, 63) * at_50[5]
    )
    ratio = 2 * c_50 / H_DENOMINATOR**2
    if not ratio < Fraction(1, 300_000_000_000):
        fail("independent perturbation ratio failed")

    certificate = artifact.get("bound_certificate", {})
    if parse_fraction(certificate.get("perturbation_ratio", "0/1")) != ratio:
        fail("artifact perturbation ratio drifted")
    if certificate.get("weights_at_50", {}) != {
        f"Khat_{j}": value for j, value in enumerate(at_50)
    }:
        fail("artifact L=50 table drifted")

    coarse = certificate.get("coarse_joined_bounds", {})
    required = {
        "p_0": "|p_0-bar(gamma)_0|<8764h^2R",
        "p_1": "|p_1-bar(gamma)_1|<8765h^2R",
        "p_2": "|p_2-bar(gamma)_2|<8762h^2R",
        "p_3": "|p_3-bar(gamma)_3|<8760h^2R",
        "p_4": "|p_4|<(2/3)h^2R",
        "p_5": "|p_5|<h^2R/63",
    }
    for key, value in required.items():
        if coarse.get(key) != value:
            fail(f"coarse joined bound drifted: {key}")


def check_artifact_and_note(artifact: dict, note: str) -> None:
    if artifact.get("kind") != STEM:
        fail("artifact kind drifted")
    counts = artifact.get("counts", {})
    expected = {
        "rows": 18,
        "basis_coefficients": 6,
        "inverse_basis_coefficients": 6,
        "translated_tail_weights": 6,
        "exact_morse_conjugations": 3,
        "joined_operator_identities": 3,
        "tail_perturbation_bounds": 1,
        "numerical_observation_row_bounds": 0,
        "grouped_ideal_cubic_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_flow_bounds": 0,
        "phi_b_bounds": 0,
    }
    if counts != expected:
        fail("artifact counts drifted")
    rows = artifact.get("rows", [])
    if len(rows) != 18 or len({row.get("id") for row in rows}) != 18:
        fail("artifact row roster drifted")
    if "no signed flow" not in artifact.get("status", ""):
        fail("artifact status boundary drifted")

    tokens = (
        "q_j=sum_(k=j)^5 binom(k,j)(-ell)^(k-j)p_k",
        "B_a=sigma-(t/2)ell",
        "L_B[P]=exp(-S_tilde)*(D_x^2+D_x+1/6)",
        "alpha dot K_alpha+i beta dot K_beta",
        "|T[P_lin-Pbar]|<h^2 R/300000000000",
        "This is not a proof",
        "grouped ideal-cubic c-prime sum",
    )
    for token in tokens:
        if token not in note:
            fail(f"note token missing: {token}")


def main() -> int:
    if not RESULT_PATH.is_file() or not NOTE_PATH.is_file():
        fail("basis-conjugation result or note is missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    check_source_hashes(artifact)
    check_basis_and_weights(artifact)
    check_physical_conjugation()
    check_joined_rows()
    check_tail_rounding(artifact)
    check_artifact_and_note(artifact, note)
    counts = artifact["counts"]
    print(
        "validated physical P_lin basis-conjugation gate: "
        f"{counts['rows']} rows, {counts['basis_coefficients']} basis "
        f"coefficients, {counts['translated_tail_weights']} translated "
        f"weights, {counts['joined_operator_identities']} joined identities, "
        f"{counts['grouped_ideal_cubic_bounds']} ideal-cubic bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
