#!/usr/bin/env python3
"""Independently check the complete first-Xi-jet Jensen layer gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_first_xi_jet_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression_hash(expression: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(expression)).encode("utf-8")).hexdigest()


def independent_polynomial(m: int, y: sp.Symbol, q: sp.Expr) -> sp.Expr:
    source = y**m
    result = sp.Integer(0)
    for total_order in range(m + 1):
        for cubic_order in range(total_order + 1):
            heat_order = total_order - cubic_order
            derivative_order = 2 * heat_order + 3 * cubic_order
            if derivative_order > m:
                continue
            result += (
                q**heat_order
                * sp.diff(source, y, derivative_order)
                / (
                    sp.factorial(heat_order)
                    * sp.factorial(cubic_order)
                    * 12**cubic_order
                )
            )
    return sp.expand(result)


def independent_tertiary(
    polynomial: sp.Expr,
    n: int | sp.Expr,
    y: sp.Symbol,
    r: sp.Expr,
) -> sp.Expr:
    shift = sp.sympify(2 * n + 1)
    return sp.expand(
        (r - y / 2) * sp.diff(polynomial, y, 2)
        + shift * sp.diff(polynomial, y) / 4
    )


def independent_background(
    polynomial: sp.Expr,
    n: int | sp.Expr,
    y: sp.Symbol,
    q: sp.Expr,
    r: sp.Expr,
    v: sp.Expr,
) -> sp.Expr:
    tertiary = independent_tertiary(polynomial, n, y, r)
    shift = sp.sympify(2 * n + 1)
    return sp.expand(
        v * sp.diff(polynomial, y, 2)
        + independent_tertiary(tertiary, n, y, r) / 2
        + q * shift * sp.diff(polynomial, y) / 2
        + (q * y + shift / 8) * sp.diff(polynomial, y, 2)
        + (q**2 + r / 2) * sp.diff(polynomial, y, 3)
        + q * sp.diff(polynomial, y, 4) / 4
        + sp.diff(polynomial, y, 5) / 80
    )


def finite_model_coefficients(m: int, n: int, unit_slope: sp.Rational) -> dict[str, sp.Expr]:
    eps, x, z, y = sp.symbols("eps x z y")
    rho = sp.Rational(-5, 3)
    q = sp.Rational(-3, 7)
    r = sp.Rational(2, 11)
    v = sp.Rational(-5, 13)
    delta = (
        rho * eps**6 / 8
        + rho * q * eps**8 / 4
        + rho * r * eps**10 / 4
        + rho * v * eps**12 / 4
    )
    source = (x - rho) ** m * (1 + unit_slope * (x - rho))

    def generator(polynomial: sp.Expr) -> sp.Expr:
        return sp.expand(
            4 * x * sp.diff(polynomial, x, 2)
            + (4 * n + 2) * sp.diff(polynomial, x)
        )

    evolved = sp.Integer(0)
    iterate = source
    for order in range(m + 2):
        evolved += delta**order * iterate / sp.factorial(order)
        iterate = generator(iterate)

    transformed = sp.Integer(0)
    for (degree,), coefficient in sp.Poly(sp.expand(evolved), x).terms():
        ratio = sp.prod(1 - index * eps**6 for index in range(degree))
        transformed += coefficient * ratio * z**degree
    local = sp.expand(transformed.subs(z, rho + rho * y * eps**4))
    return {
        "lead": sp.expand(local.coeff(eps, 4 * m) / rho**m),
        "odd_one": sp.expand(local.coeff(eps, 4 * m + 1) / rho**m),
        "tertiary": sp.expand(local.coeff(eps, 4 * m + 2) / rho**m),
        "odd_three": sp.expand(local.coeff(eps, 4 * m + 3) / rho**m),
        "quaternary": sp.expand(local.coeff(eps, 4 * m + 4) / rho**m),
        "odd_five": sp.expand(local.coeff(eps, 4 * m + 5) / rho**m),
    }


def main() -> None:
    require(RESULT_PATH.exists(), "missing first-Xi-jet result")
    require(NOTE_PATH.exists(), "missing first-Xi-jet note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(artifact.get("kind") == KIND, "kind drift")

    expected_counts = {
        "rows": 16,
        "sources": 3,
        "multiplicities": 14,
        "operator_checks": 28,
        "exact_finite_models": 12,
        "shift_values": 2,
        "appell_identities": 28,
        "full_xi_jet_checks": 12,
        "exact_m3_contact_corrections": 2,
        "canonical_product_identities": 2,
        "sign_countermodels": 2,
        "open_handoffs": 2,
        "degree_uniform_bounds": 0,
        "rh_conclusions": 0,
    }
    require(artifact.get("counts") == expected_counts, "count drift")

    source_audit = artifact.get("source_audit", {})
    require(len(source_audit) == 3, "source count drift")
    for source in source_audit.values():
        path = REPO_ROOT / source["path"]
        require(path.exists(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")

    expected_ids = [
        "fxj_01_sources",
        "fxj_02_scaling",
        "fxj_03_expansion",
        "fxj_04_tertiary",
        "fxj_05_operator",
        "fxj_06_normal",
        "fxj_07_unit",
        "fxj_08_appell",
        "fxj_09_finite",
        "fxj_10_contact",
        "fxj_11_newman",
        "fxj_12_product",
        "fxj_13_countermodels",
        "fxj_14_actual_xi",
        "fxj_15_uniform",
        "fxj_16_boundary",
    ]
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == expected_ids, "row id drift")
    require(sum(row.get("readiness") == "open" for row in rows) == 2, "open-row drift")

    y, q, r, v = sp.symbols("y q r v")
    operator_rows = artifact.get("operator_audit", [])
    require(len(operator_rows) == 14, "operator-row drift")
    operator_checks = 0
    appell_checks = 0
    for expected_m, row in zip(range(3, 17), operator_rows, strict=True):
        require(row.get("multiplicity") == expected_m, f"multiplicity drift m={expected_m}")
        primary = independent_polynomial(expected_m, y, q)
        require(expression_hash(primary) == row["primary_sha256"], f"primary hash drift m={expected_m}")
        require(
            sp.expand(sp.diff(primary, q) - sp.diff(primary, y, 2)) == 0,
            f"heat Appell drift m={expected_m}",
        )
        require(
            sp.expand(
                sp.diff(independent_polynomial(expected_m + 1, y, q), y)
                - (expected_m + 1) * primary
            )
            == 0,
            f"degree Appell drift m={expected_m}",
        )
        shift_index = {entry["shift"]: entry for entry in row["shift_audit"]}
        for n in (0, 3):
            tertiary = independent_tertiary(primary, n, y, r)
            background = independent_background(primary, n, y, q, r, v)
            stored = shift_index[n]
            require(expression_hash(tertiary) == stored["tertiary_sha256"], f"tertiary hash drift m={expected_m}, n={n}")
            require(expression_hash(background) == stored["quaternary_background_sha256"], f"background hash drift m={expected_m}, n={n}")
            operator_checks += 1
        appell_checks += 2
    require(operator_checks == 28, "operator count drift")
    require(appell_checks == 28, "Appell count drift")

    finite_rows = artifact.get("exact_finite_model_audit", [])
    require(len(finite_rows) == 12, "stored finite-model count drift")
    require(all(row.get("matched") is True for row in finite_rows), "stored finite-model match drift")
    rho_value = sp.Rational(-5, 3)
    q_value = sp.Rational(-3, 7)
    r_value = sp.Rational(2, 11)
    v_value = sp.Rational(-5, 13)
    unit_value = sp.Rational(7, 5)
    finite_checks = 0
    unit_checks = 0
    for n in (0, 3):
        for m in range(3, 9):
            with_unit = finite_model_coefficients(m, n, unit_value)
            without_unit = finite_model_coefficients(m, n, sp.Rational(0))
            primary = independent_polynomial(m, y, q_value)
            tertiary = independent_tertiary(primary, n, y, r_value)
            background = independent_background(primary, n, y, q_value, r_value, v_value)
            unit_target = sp.expand(rho_value * independent_polynomial(m + 1, y, q_value))
            complete = sp.expand(background + unit_value * unit_target)
            require(sp.expand(with_unit["lead"] - primary) == 0, f"finite lead drift m={m}, n={n}")
            require(with_unit["odd_one"] == 0, f"epsilon-one drift m={m}, n={n}")
            require(sp.expand(with_unit["tertiary"] - tertiary) == 0, f"finite tertiary drift m={m}, n={n}")
            require(with_unit["odd_three"] == 0, f"epsilon-three drift m={m}, n={n}")
            require(sp.expand(with_unit["quaternary"] - complete) == 0, f"finite quaternary drift m={m}, n={n}")
            require(with_unit["odd_five"] == 0, f"epsilon-five drift m={m}, n={n}")
            isolated = sp.expand((with_unit["quaternary"] - without_unit["quaternary"]) / unit_value)
            require(sp.expand(isolated - unit_target) == 0, f"unit isolation drift m={m}, n={n}")
            finite_checks += 1
            unit_checks += 1
    require(finite_checks == 12, "finite-model count drift")
    require(unit_checks == 12, "unit-isolation count drift")

    n, ell, t = sp.symbols("n ell t")
    p3 = independent_polynomial(3, y, q)
    p4 = independent_polynomial(4, y, q)
    q3 = -2 ** sp.Rational(-7, 3)
    y3 = 2 ** sp.Rational(-2, 3)
    r3 = y3 / 2
    s3 = (1 - 2 * n) / 4
    tertiary3 = independent_tertiary(p3, n, y, r)
    fourth3 = independent_background(p3, n, y, q, r, v) + ell * p4
    value_equation = sp.simplify(
        s3**2 * sp.diff(p3, y, 2) / 2
        + s3 * sp.diff(tertiary3, y)
        + fourth3
    ).subs({q: q3, y: y3, r: r3})
    v3 = sp.solve(sp.simplify(value_equation), v)[0]
    derivative_equation = sp.simplify(
        t * sp.diff(p3, y, 2)
        + s3**2 * sp.diff(p3, y, 3) / 2
        + s3 * sp.diff(tertiary3, y, 2)
        + sp.diff(fourth3, y)
    ).subs({q: q3, y: y3, r: r3, v: v3})
    t3 = sp.solve(sp.simplify(derivative_equation), t)[0]
    require(sp.simplify(v3 + (2 * ell + 3 * n + 2) / 8) == 0, "m3 heat correction drift")
    require(sp.simplify(t3 - 2 ** sp.Rational(-7, 3) * (2 * ell + n + 4)) == 0, "m3 root correction drift")

    s = sp.symbols("s")
    model_rows = artifact.get("canonical_product_audit", {}).get("models", [])
    require(len(model_rows) == 2, "sign-model count drift")
    expected_ells = []
    for other_root in (sp.Integer(-4), sp.Rational(-1, 4)):
        source = sp.expand((s + 1) ** 3 * (s - other_root))
        jet = sp.simplify(
            sp.diff(source, s, 4).subs(s, -1)
            / (4 * sp.diff(source, s, 3).subs(s, -1))
        )
        expected_ells.append(sp.simplify(-jet))
        require(all(coefficient > 0 for coefficient in sp.Poly(source, s).all_coeffs()), "positive coefficient witness drift")
    require(expected_ells == [sp.Rational(-1, 3), sp.Rational(4, 3)], "ell sign witnesses drift")

    note = NOTE_PATH.read_text(encoding="utf-8")
    for marker in (
        "P'''''/80",
        "v_3=-(2*ell+3*n+2)/8",
        "t_3=2^(-7/3)*(2*ell+n+4)",
        "ell=-1/3",
        "ell=4/3",
        "No pi enters this local operator calculation.",
        "This gate proves no Xi regular-field sign",
    ):
        require(marker in note, f"note marker missing: {marker}")

    print(
        "validated first Xi-jet Jensen layer gate: "
        "16 rows, 3 sources, 14 multiplicities, 28 independent operators, "
        "12 independent finite models, 28 Appell identities, "
        "12 full Xi-jet isolations, 2 exact m3 corrections, "
        "2 sign countermodels, 2 open handoffs, 0 degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
