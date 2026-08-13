#!/usr/bin/env python3
"""Independently check the tertiary universal Jensen correction gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_tertiary_universal_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_polynomial(m: int, y: sp.Symbol, q: sp.Expr) -> sp.Expr:
    """Expand the two commuting exponentials directly on y**m."""
    source = y**m
    result = sp.Integer(0)
    for heat_order in range(m // 2 + 1):
        after_heat = sp.diff(source, y, 2 * heat_order)
        if after_heat == 0:
            continue
        for cubic_order in range((m - 2 * heat_order) // 3 + 1):
            result += (
                q**heat_order
                * sp.diff(after_heat, y, 3 * cubic_order)
                / (
                    sp.factorial(heat_order)
                    * sp.factorial(cubic_order)
                    * 12**cubic_order
                )
            )
    return sp.expand(result)


def independent_tertiary(
    m: int,
    n: int,
    y: sp.Symbol,
    q: sp.Expr,
    r: sp.Expr,
) -> sp.Expr:
    primary = independent_polynomial(m, y, q)
    return sp.expand(
        (r - y / 2) * sp.diff(primary, y, 2)
        + (2 * n + 1) * sp.diff(primary, y) / 4
    )


def finite_model_coefficients(m: int, n: int, unit_slope: sp.Rational) -> dict[str, sp.Expr]:
    """Recompute one exact radial-heat/Jensen model at fixed rational data."""
    eps, x, z, y = sp.symbols("eps x z y")
    rho = sp.Rational(-3, 2)
    q = sp.Rational(-2, 5)
    r = sp.Rational(1, 7)
    delta = rho * eps**6 / 8 + rho * q * eps**8 / 4 + rho * r * eps**10 / 4
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
        falling_ratio = sp.prod(1 - index * eps**6 for index in range(degree))
        transformed += coefficient * falling_ratio * z**degree
    local = sp.expand(transformed.subs(z, rho + rho * y * eps**4))
    return {
        "lead": sp.expand(local.coeff(eps, 4 * m) / rho**m),
        "odd_one": sp.expand(local.coeff(eps, 4 * m + 1) / rho**m),
        "tertiary": sp.expand(local.coeff(eps, 4 * m + 2) / rho**m),
        "odd_three": sp.expand(local.coeff(eps, 4 * m + 3) / rho**m),
        "fourth": sp.expand(local.coeff(eps, 4 * m + 4) / rho**m),
    }


def main() -> None:
    require(RESULT_PATH.exists(), "missing tertiary universal result")
    require(NOTE_PATH.exists(), "missing tertiary universal note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(artifact.get("kind") == KIND, "kind drift")

    expected_counts = {
        "rows": 15,
        "sources": 3,
        "multiplicities": 14,
        "operator_checks": 14,
        "exact_finite_models": 12,
        "shift_values": 2,
        "appell_identities": 28,
        "unit_isolation_checks": 12,
        "exact_m3_contact_corrections": 2,
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
        "tul_01_sources",
        "tul_02_scaling",
        "tul_03_expansion",
        "tul_04_operator",
        "tul_05_cancellation",
        "tul_06_shift",
        "tul_07_appell",
        "tul_08_unit",
        "tul_09_contact",
        "tul_10_m3",
        "tul_11_newman",
        "tul_12_double",
        "tul_13_xi",
        "tul_14_uniform",
        "tul_15_boundary",
    ]
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == expected_ids, "row id drift")
    require(sum(row.get("readiness") == "open" for row in rows) == 2, "open-row drift")

    y, q, r = sp.symbols("y q r")
    operator_rows = artifact.get("operator_audit", [])
    require(len(operator_rows) == 14, "operator-row drift")
    appell_checks = 0
    for expected_m, row in zip(range(3, 17), operator_rows, strict=True):
        require(row.get("multiplicity") == expected_m, f"multiplicity drift m={expected_m}")
        primary = independent_polynomial(expected_m, y, q)
        observed_primary = sp.sympify(row["secondary_polynomial"])
        require(sp.expand(observed_primary - primary) == 0, f"P_m drift m={expected_m}")
        for n, key in ((0, "tertiary_shift_zero"), (3, "tertiary_shift_three")):
            observed = sp.sympify(row[key])
            target = independent_tertiary(expected_m, n, y, q, r)
            require(sp.expand(observed - target) == 0, f"Q_mn drift m={expected_m}, n={n}")
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
        appell_checks += 2
    require(appell_checks == 28, "Appell count drift")

    finite_rows = artifact.get("exact_finite_model_audit", [])
    finite_index = {(row["multiplicity"], row["shift"]): row for row in finite_rows}
    require(len(finite_index) == 12, "finite-model row drift")
    rho_value = sp.Rational(-3, 2)
    q_value = sp.Rational(-2, 5)
    r_value = sp.Rational(1, 7)
    unit_value = sp.Rational(7, 3)
    finite_checks = 0
    unit_checks = 0
    for n in (0, 3):
        for m in range(3, 9):
            with_unit = finite_model_coefficients(m, n, unit_value)
            without_unit = finite_model_coefficients(m, n, sp.Rational(0))
            target_primary = independent_polynomial(m, y, q_value)
            target_tertiary = independent_tertiary(m, n, y, q_value, r_value)
            require(sp.expand(with_unit["lead"] - target_primary) == 0, f"finite lead drift m={m}, n={n}")
            require(with_unit["odd_one"] == 0, f"epsilon-one drift m={m}, n={n}")
            require(sp.expand(with_unit["tertiary"] - target_tertiary) == 0, f"finite Q drift m={m}, n={n}")
            require(with_unit["odd_three"] == 0, f"epsilon-three drift m={m}, n={n}")
            isolated_unit = sp.expand((with_unit["fourth"] - without_unit["fourth"]) / unit_value)
            target_unit = sp.expand(rho_value * independent_polynomial(m + 1, y, q_value))
            require(sp.expand(isolated_unit - target_unit) == 0, f"unit isolation drift m={m}, n={n}")

            stored = finite_index[(m, n)]
            stored_primary = sp.sympify(stored["secondary_limit"]).subs(q, q_value)
            stored_tertiary = sp.sympify(stored["tertiary_correction"]).subs({q: q_value, r: r_value})
            stored_unit = sp.sympify(stored["first_local_unit_coefficient"]).subs({q: q_value, sp.Symbol("rho"): rho_value})
            require(sp.expand(stored_primary - target_primary) == 0, f"stored lead drift m={m}, n={n}")
            require(sp.expand(stored_tertiary - target_tertiary) == 0, f"stored Q drift m={m}, n={n}")
            require(sp.expand(stored_unit - target_unit) == 0, f"stored unit drift m={m}, n={n}")
            finite_checks += 1
            unit_checks += 1
    require(finite_checks == 12, "finite-model count drift")
    require(unit_checks == 12, "unit-isolation count drift")

    p3 = independent_polynomial(3, y, q)
    q3 = -2 ** sp.Rational(-7, 3)
    y3 = 2 ** sp.Rational(-2, 3)
    n, s = sp.symbols("n s")
    tertiary3 = independent_tertiary(3, n, y, q, r)
    value_equation = sp.simplify(tertiary3.subs({q: q3, y: y3}))
    r3 = sp.solve(value_equation, r)[0]
    derivative_equation = sp.simplify(
        s * sp.diff(p3, y, 2).subs({q: q3, y: y3})
        + sp.diff(tertiary3, y).subs({q: q3, y: y3, r: r3})
    )
    s3 = sp.solve(derivative_equation, s)[0]
    require(sp.simplify(r3 - 2 ** sp.Rational(-5, 3)) == 0, "m=3 r correction drift")
    require(sp.simplify(s3 - (1 - 2 * n) / 4) == 0, "m=3 root correction drift")

    certificate = artifact.get("symbolic_certificate", {})
    require("no independent partial_y^4" in certificate.get("cancellation", ""), "cancellation marker missing")
    require("rho*[U'(rho)/U(rho)]*P_(m+1)" in certificate.get("unit_delay", ""), "Xi-unit marker missing")
    require("still-uncomputed universal background" in certificate.get("xi_handoff", ""), "open-background marker missing")

    note = NOTE_PATH.read_text(encoding="utf-8")
    for marker in (
        "Q_(m,n)=((r-y/2)*partial_y^2+(2n+1)*partial_y/4)P_m.",
        "r_3=2^(-5/3)",
        "rho*[U'(rho)/U(rho)]*P_(m+1)(y,q)",
        "No pi enters this local operator calculation.",
        "This gate proves no Xi-specific collision exclusion",
    ):
        require(marker in note, f"note marker missing: {marker}")

    print(
        "validated tertiary universal Jensen layer gate: "
        "15 rows, 3 sources, 14 multiplicities, 14 independent operators, "
        "12 independent finite models, 28 Appell identities, "
        "12 Xi-unit isolations, 2 exact m3 corrections, "
        "2 open handoffs, 0 degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
