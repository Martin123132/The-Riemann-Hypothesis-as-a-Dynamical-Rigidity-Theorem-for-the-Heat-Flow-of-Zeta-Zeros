#!/usr/bin/env python3
"""Independently check the secondary cubic Jensen boundary layer gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_arbitrary_multiplicity_secondary_cubic_layer_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_polynomial(m: int, y: sp.Symbol, q: sp.Symbol) -> sp.Expr:
    result = sp.Integer(0)
    source = y**m
    for total_order in range(m + 1):
        for cubic_order in range(total_order + 1):
            heat_order = total_order - cubic_order
            derivative_order = 2 * heat_order + 3 * cubic_order
            if derivative_order > m:
                continue
            result += (
                q**heat_order
                * sp.diff(source, y, derivative_order)
                / (sp.factorial(heat_order) * sp.factorial(cubic_order) * 12**cubic_order)
            )
    return sp.expand(result)


def independent_numeric_model(m: int, n: int) -> sp.Expr:
    eps, x, z, y = sp.symbols("eps x z y")
    rho = sp.Rational(-3, 2)
    q = sp.Rational(-2, 5)
    u1 = sp.Rational(7, 3)
    delta = rho * eps**6 / 8 + rho * q * eps**8 / 4
    source = (x - rho) ** m * (1 + u1 * (x - rho))

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
        ratio = sp.prod(1 - j * eps**6 for j in range(degree))
        transformed += coefficient * ratio * z**degree
    local = sp.expand(transformed.subs(z, rho + rho * y * eps**4))
    for power in range(4 * m):
        require(local.coeff(eps, power) == 0, f"numeric valuation drift m={m}, n={n}")
    return sp.expand(local.coeff(eps, 4 * m))


def main() -> None:
    require(RESULT_PATH.exists(), "missing secondary cubic result")
    require(NOTE_PATH.exists(), "missing secondary cubic note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(artifact.get("kind") == KIND, "kind drift")

    expected_counts = {
        "rows": 16,
        "sources": 4,
        "multiplicities": 14,
        "universal_coefficient_checks": 200,
        "exact_finite_models": 12,
        "fixed_shifts": 2,
        "open_handoffs": 2,
        "xi_specific_signs": 0,
        "degree_uniform_bounds": 0,
        "rh_conclusions": 0,
    }
    require(artifact.get("counts") == expected_counts, "count drift")

    source_audit = artifact.get("source_audit", {})
    require(len(source_audit) == 4, "source count drift")
    for source in source_audit.values():
        path = REPO_ROOT / source["path"]
        require(path.exists(), f"missing source: {path}")
        require(file_hash(path) == source["sha256"], f"source hash drift: {path}")

    rows = artifact.get("rows", [])
    expected_ids = [f"scb_{index:02d}_{suffix}" for index, suffix in [
        (1, "sources"),
        (2, "scaling"),
        (3, "limit"),
        (4, "coefficients"),
        (5, "cubic_origin"),
        (6, "universality"),
        (7, "heat_ray"),
        (8, "threshold"),
        (9, "collision"),
        (10, "m3"),
        (11, "newman"),
        (12, "double"),
        (13, "first_jet"),
        (14, "tertiary"),
        (15, "uniform"),
        (16, "boundary"),
    ]]
    require([row.get("id") for row in rows] == expected_ids, "row id drift")
    require(sum(row.get("readiness") == "open" for row in rows) == 2, "open-row drift")

    y, q = sp.symbols("y q")
    coefficient_rows = artifact.get("coefficient_audit", [])
    require(len(coefficient_rows) == 14, "coefficient-row drift")
    independent_checks = 0
    for row in coefficient_rows:
        m = row["multiplicity"]
        observed = sp.sympify(row["polynomial"])
        target = independent_polynomial(m, y, q)
        require(sp.expand(observed - target) == 0, f"independent polynomial drift m={m}")
        require(observed.coeff(y, m - 1) == 0, f"sum-root coefficient drift m={m}")
        at_zero = sp.expand(observed.subs(q, 0))
        require(at_zero.coeff(y, m - 2) == 0, f"pair-sum coefficient drift m={m}")
        require(at_zero.coeff(y, m - 3) != 0, f"cubic witness missing m={m}")
        independent_checks += len(row["coefficients"])
    require(independent_checks == 200, "independent coefficient count drift")

    rho = sp.symbols("rho")
    tau0 = rho / 8
    jensen_cubic = rho**3 / 3
    heat_cubic = -16 * rho * tau0**2
    require(sp.simplify(jensen_cubic + heat_cubic - rho**3 / 12) == 0, "cubic residual drift")

    numeric_checks = 0
    for n in (0, 3):
        for m in range(3, 9):
            observed = independent_numeric_model(m, n)
            target = sp.expand(
                sp.Rational(-3, 2) ** m
                * independent_polynomial(m, y, q).subs(q, sp.Rational(-2, 5))
            )
            require(sp.expand(observed - target) == 0, f"numeric finite model drift m={m}, n={n}")
            numeric_checks += 1
    require(numeric_checks == 12, "numeric model count drift")

    p3 = independent_polynomial(3, y, q)
    require(p3 == y**3 + 6 * q * y + sp.Rational(1, 2), "P3 drift")
    disc3 = sp.discriminant(p3, y)
    require(
        sp.expand(disc3 + sp.Rational(27, 4) * (128 * q**3 + 1)) == 0,
        "P3 discriminant drift",
    )
    q3 = -2 ** sp.Rational(-7, 3)
    y3 = 2 ** sp.Rational(-2, 3)
    require(sp.simplify(p3.subs({q: q3, y: y3})) == 0, "P3 contact root drift")
    require(sp.simplify(sp.diff(p3, y).subs({q: q3, y: y3})) == 0, "P3 contact derivative drift")

    certificate = artifact.get("symbolic_certificate", {})
    require("D^(-4/3)" in certificate.get("collision_refinement", ""), "parameter scale missing")
    require("D^(-5/3)" in certificate.get("collision_refinement", ""), "root scale missing")
    require("D^(-2/3)" in certificate.get("open_target", ""), "Xi-unit order missing")

    note = NOTE_PATH.read_text(encoding="utf-8")
    for marker in (
        "exp(q*partial_y^2+partial_y^3/12)y^m",
        "q_m<0",
        "No pi enters this universal local calculation.",
        "It proves no Xi-specific exclusion",
    ):
        require(marker in note, f"note marker missing: {marker}")

    print(
        "validated secondary cubic Jensen layer gate: "
        "16 rows, 4 sources, 14 multiplicities, 200 universal coefficients, "
        "12 independent finite models, 1 cubic commutator, 1 exact m=3 threshold, "
        "2 open handoffs, 0 Xi-specific signs, 0 degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
