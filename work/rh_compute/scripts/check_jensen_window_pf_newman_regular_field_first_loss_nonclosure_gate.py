#!/usr/bin/env python3
"""Independently check the arbitrary-field first-loss nonclosure gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_newman_regular_field_first_loss_nonclosure_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def expression_hash(expression: sp.Expr) -> str:
    return hashlib.sha256(sp.srepr(sp.expand(expression)).encode("utf-8")).hexdigest()


def independent_heat_monomial(m: int, z: sp.Symbol, sign: int) -> sp.Expr:
    polynomial = z**m
    total = sp.Integer(0)
    iterate = polynomial
    order = 0
    while iterate != 0:
        total += sign**order * iterate / sp.factorial(order)
        iterate = sp.diff(iterate, z, 2)
        order += 1
    return sp.expand(total)


def independent_evolution(polynomial: sp.Expr, x: sp.Symbol, tau: sp.Rational) -> sp.Expr:
    total = sp.Integer(0)
    iterate = polynomial
    order = 0
    while iterate != 0:
        total += (-tau) ** order * iterate / sp.factorial(order)
        iterate = sp.diff(iterate, x, 2)
        order += 1
    return sp.expand(total)


def main() -> None:
    require(RESULT_PATH.exists(), "missing regular-field result")
    require(NOTE_PATH.exists(), "missing regular-field note")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    require(artifact.get("kind") == KIND, "kind drift")

    expected_counts = {
        "rows": 16,
        "sources": 3,
        "field_models": 10,
        "multiplicities": 15,
        "positive_coefficient_checks": 10,
        "negative_root_checks": 10,
        "signed_field_checks": 10,
        "positive_time_samples": 10,
        "forward_hermite_checks": 15,
        "backward_nonreal_checks": 15,
        "open_handoffs": 2,
        "actual_xi_field_bounds": 0,
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
        "rfn_01_sources",
        "rfn_02_parameters",
        "rfn_03_product",
        "rfn_04_coefficients",
        "rfn_05_field",
        "rfn_06_even",
        "rfn_07_conversion",
        "rfn_08_preserver",
        "rfn_09_forward",
        "rfn_10_drift",
        "rfn_11_backward",
        "rfn_12_boundary",
        "rfn_13_multiplicity",
        "rfn_14_nonclosure",
        "rfn_15_xi",
        "rfn_16_uniform",
    ]
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == expected_ids, "row id drift")
    require(sum(row.get("readiness") == "open" for row in rows) == 2, "open-row drift")

    ell, s = sp.symbols("ell s", real=True)
    c = sp.symbols("c", positive=True)
    A = ell**2 + 2
    B = ell**2 - ell + 2
    alpha = sp.cancel((ell**2 + 1) / A)
    beta = sp.cancel((ell**2 - ell + 3) / B)
    rho = -c**2
    unit = (s + alpha * c**2) * (s + beta * c**2)
    field = sp.factor(rho * sp.diff(unit, s).subs(s, rho) / unit.subs(s, rho))
    require(sp.simplify(field - ell) == 0, "independent arbitrary-field identity failed")
    require(sp.expand(B - ((ell - sp.Rational(1, 2)) ** 2 + sp.Rational(7, 4))) == 0, "B positivity decomposition failed")

    stored_models = artifact.get("model_audit", [])
    require(len(stored_models) == 10, "stored model count drift")
    s_var, x = sp.symbols("s x")
    for model in stored_models:
        field_value = sp.Rational(model["field"])
        m = int(model["multiplicity"])
        c_value = sp.Rational(model["c"])
        A_value = field_value**2 + 2
        B_value = field_value**2 - field_value + 2
        alpha_value = 1 - 1 / A_value
        beta_value = 1 + 1 / B_value
        source = sp.expand(
            (s_var + c_value**2) ** m
            * (s_var + alpha_value * c_value**2)
            * (s_var + beta_value * c_value**2)
        )
        lifted = sp.expand(source.subs(s_var, -x**2))
        unit_value = (s_var + alpha_value * c_value**2) * (s_var + beta_value * c_value**2)
        rho_value = -c_value**2
        recovered = sp.factor(
            rho_value
            * sp.diff(unit_value, s_var).subs(s_var, rho_value)
            / unit_value.subs(s_var, rho_value)
        )
        regular = sp.cancel(lifted / (x - c_value) ** m)
        signed = sp.factor(sp.diff(regular, x).subs(x, c_value) / regular.subs(x, c_value))
        coefficients = sp.Poly(source, s_var).all_coeffs()
        require(0 < alpha_value < 1 < beta_value, f"independent ordering failed ell={field_value}")
        require(all(coefficient > 0 for coefficient in coefficients), f"independent coefficient sign failed ell={field_value}")
        require(recovered == field_value, f"independent field failed ell={field_value}")
        require(signed == (m + 4 * field_value) / (2 * c_value), f"independent signed field failed ell={field_value}")
        require(expression_hash(source) == model["source_sha256"], f"source hash drift ell={field_value}")
        tau = sp.Rational(model["positive_time"])
        evolved = independent_evolution(lifted, x, tau)
        degree = sp.Poly(evolved, x).degree()
        real_count = sp.Poly(evolved, x).count_roots(-sp.oo, sp.oo)
        require(degree == model["positive_time_degree"], f"heat degree drift ell={field_value}")
        require(real_count == degree, f"independent positive-time root count failed ell={field_value}")
        require(real_count == model["positive_time_real_roots"], f"stored root count drift ell={field_value}")
        require(expression_hash(evolved) == model["positive_time_sha256"], f"heat hash drift ell={field_value}")

    z = sp.symbols("z")
    stored_hermite = artifact.get("hermite_audit", [])
    require(len(stored_hermite) == 15, "stored Hermite count drift")
    for m, stored in zip(range(2, 17), stored_hermite, strict=True):
        forward = independent_heat_monomial(m, z, -1)
        backward = independent_heat_monomial(m, z, 1)
        next_forward = independent_heat_monomial(m + 1, z, -1)
        require(stored["multiplicity"] == m, f"Hermite multiplicity drift m={m}")
        require(expression_hash(forward) == stored["forward_sha256"], f"forward hash drift m={m}")
        require(expression_hash(backward) == stored["backward_sha256"], f"backward hash drift m={m}")
        require(sp.expand(next_forward - z * forward + 2 * sp.diff(forward, z)) == 0, f"independent recurrence failed m={m}")
        forward_real = sp.Poly(forward, z).count_roots(-sp.oo, sp.oo)
        backward_real = sp.Poly(backward, z).count_roots(-sp.oo, sp.oo)
        require(forward_real == m, f"independent forward roots failed m={m}")
        require(backward_real == m % 2, f"independent backward roots failed m={m}")
        require(stored["forward_real_roots"] == forward_real, f"stored forward roots drift m={m}")
        require(stored["backward_real_roots"] == backward_real, f"stored backward roots drift m={m}")

    note = NOTE_PATH.read_text(encoding="utf-8")
    required_phrases = [
        "every finite real regular field",
        "exact Newman-style threshold",
        "cannot change the real-versus-nonreal side",
        "genuine nonclosure theorem",
        "actual Xi field is arbitrary",
        "No `pi` enters this construction",
        "not a proof of RH",
    ]
    for phrase in required_phrases:
        require(phrase in note, f"note missing required phrase: {phrase}")

    print(
        "validated regular-field first-loss nonclosure gate: "
        "16 rows, 3 sources, 10 arbitrary-field models, 15 multiplicities, "
        "10 positive-time samples, 15 forward Hermite checks, "
        "15 backward nonreal checks, 2 open handoffs, 0 Xi field bounds, "
        "0 degree-uniform bounds"
    )


if __name__ == "__main__":
    main()
