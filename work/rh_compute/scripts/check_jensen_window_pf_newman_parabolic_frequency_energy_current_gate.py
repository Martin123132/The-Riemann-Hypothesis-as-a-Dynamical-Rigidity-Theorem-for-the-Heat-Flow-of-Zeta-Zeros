#!/usr/bin/env python3
"""Independently validate the parabolic-frequency energy/current audit."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_parabolic_frequency_energy_current_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_parabolic_frequency_energy_current_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "pfec_01_heat_setup",
    "pfec_02_direct_derivative",
    "pfec_03_flux_derivative",
    "pfec_04_completed_balance",
    "pfec_05_scale_derivatives",
    "pfec_06_coefficient_factorization",
    "pfec_07_outer_ratio_bound",
    "pfec_08_parabolic_chart",
    "pfec_09_frequency_chart",
    "pfec_10_fixed_slice_balance",
    "pfec_11_spacetime_balance",
    "pfec_12_spatial_gradient_bridge",
    "pfec_13_conditional_no_contact",
    "pfec_14_contact_cancellation",
    "pfec_15_quadratic_countermodel",
    "pfec_16_integral_nonpromotion",
    "pfec_17_conditional_utility",
    "pfec_18_route_handoff",
]


def load_json(path: Path, label: str, issues: list[str]) -> dict:
    if not path.is_file():
        issues.append(f"missing {label}")
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        issues.append(f"invalid {label}: {exc}")
        return {}


def validate_general_balance(issues: list[str]) -> None:
    H, Hx, Hxx, Hxxx = sp.symbols("H Hx Hxx Hxxx", real=True)
    s, st, sx = sp.symbols("s st sx", real=True)
    energy_t = -2 * H * Hxx + 2 * s * st * Hx**2 - 2 * s**2 * Hx * Hxxx
    flux_x = (
        Hx**2
        + H * Hxx
        + 2 * s * sx * Hx * Hxx
        + s**2 * Hxx**2
        + s**2 * Hx * Hxxx
    )
    completed = (
        2 * (s * Hxx + sx * Hx) ** 2
        + 2 * (1 + s * st - sx**2) * Hx**2
    )
    if sp.expand(energy_t + 2 * flux_x - completed) != 0:
        issues.append("independent general balance failed")


def validate_scale(issues: list[str]) -> None:
    t, x = sp.symbols("t x", positive=True)
    L = sp.log(x / (4 * sp.pi))
    q = 2 * t * L**2
    s = (L**2 + 1 / (2 * t)) ** sp.Rational(-1, 2)
    targets = {
        "time": 1 / (2 * t * (1 + q)),
        "space": -2 * t * L / (x * (1 + q)),
        "sst": 1 / (1 + q) ** 2,
        "sx2": 4 * q * t**2 / (x**2 * (1 + q) ** 3),
    }
    values = {
        "time": sp.diff(sp.log(s), t),
        "space": sp.diff(sp.log(s), x),
        "sst": s * sp.diff(s, t),
        "sx2": sp.diff(s, x) ** 2,
    }
    for key, value in values.items():
        if sp.simplify(value - targets[key]) != 0:
            issues.append(f"independent scale derivative failed: {key}")


def validate_coefficient(issues: list[str]) -> None:
    q, t, x = sp.symbols("q t x", positive=True)
    sst = 1 / (1 + q) ** 2
    sx2 = 4 * q * t**2 / (x**2 * (1 + q) ** 3)
    ratio = 4 * q * t**2 / (x**2 * (1 + q))
    coefficient = 1 + sst - sx2
    factorized = 1 + (1 - ratio) / (1 + q) ** 2
    if sp.simplify(coefficient - factorized) != 0:
        issues.append("independent coefficient factorization failed")
    if 4 * Fraction(1, 4) ** 2 / 38**2 != Fraction(1, 5776):
        issues.append("independent outer ratio cap failed")
    q_le_one_bound = 1 + (
        1 - Fraction(1, 5776)
    ) * Fraction(1, 4)
    if q_le_one_bound != Fraction(28879, 23104):
        issues.append("independent q<=1 coefficient bound failed")


def validate_contact_model(issues: list[str]) -> None:
    x, t, x0, t0, s0 = sp.symbols(
        "x t x0 t0 s0", real=True, positive=True
    )
    H = (x - x0) ** 2 - 2 * (t - t0)
    if sp.simplify(sp.diff(H, t) + sp.diff(H, x, 2)) != 0:
        issues.append("independent quadratic heat equation failed")
    Hx = sp.diff(H, x)
    Hxx = sp.diff(H, x, 2)
    E = H**2 + s0**2 * Hx**2
    J = H * Hx + s0**2 * Hx * Hxx
    rhs = 2 * s0**2 * Hxx**2 + 2 * Hx**2
    if sp.expand(sp.diff(E, t) + 2 * sp.diff(J, x) - rhs) != 0:
        issues.append("independent quadratic balance failed")
    contact = {x: x0, t: t0}
    expected = {
        "E": 0,
        "Et": 0,
        "J": 0,
        "Jx": 4 * s0**2,
        "rhs": 8 * s0**2,
    }
    values = {
        "E": E.subs(contact),
        "Et": sp.diff(E, t).subs(contact),
        "J": J.subs(contact),
        "Jx": sp.diff(J, x).subs(contact),
        "rhs": rhs.subs(contact),
    }
    for key, value in values.items():
        if sp.simplify(value - expected[key]) != 0:
            issues.append(f"independent contact value failed: {key}")


def validate_spatial_bridge(issues: list[str]) -> None:
    Hx, D, A = sp.symbols("Hx D A", real=True)
    bulk = D**2 + A * Hx**2
    jet_gradient_squared = D**2 + Hx**2
    if sp.expand(
        bulk - jet_gradient_squared - (A - 1) * Hx**2
    ) != 0:
        issues.append("independent spatial-gradient domination failed")

    u, v, left, right = sp.symbols(
        "u v left right", positive=True
    )
    two_side_cost = u**2 / left + v**2 / right
    endpoint_floor = (u + v) ** 2 / (left + right)
    square_remainder = (
        (u * right - v * left) ** 2
        / (left * right * (left + right))
    )
    if sp.simplify(
        two_side_cost - endpoint_floor - square_remainder
    ) != 0:
        issues.append("independent two-endpoint contact floor failed")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "stored energy/current result", issues)
    if not NOTE.is_file():
        issues.append("missing rendered note")
    if not stored:
        return issues

    rebuilt = gate.build_payload()
    if stored != rebuilt:
        issues.append("stored payload differs from deterministic reconstruction")
    if stored.get("kind") != STEM:
        issues.append("kind drifted")
    if stored.get("status") != (
        "exact parabolic-frequency energy-current balance with "
        "positive bulk and pointwise nonpromotion guard"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_summary = {
        "rows": 18,
        "exact_identity_rows": 8,
        "coefficient_bound_rows": 3,
        "pointwise_bridge_rows": 2,
        "nonpromotion_guard_rows": 3,
        "open_xi_targets": 1,
        "pointwise_contact_exclusion": False,
        "q209_certified": False,
        "cofinal_descendant_theorem": False,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    validate_general_balance(issues)
    validate_scale(issues)
    validate_coefficient(issues)
    validate_contact_model(issues)
    validate_spatial_bridge(issues)

    exact = stored.get("exact", {})
    coefficient = exact.get("bulk_coefficient", {})
    for phrase in [
        "partial_t E_s+2 partial_x J_s=",
        "2(s H_xx+s_x H_x)^2",
        "2(1+s s_t-s_x^2)H_x^2",
    ]:
        if phrase not in exact.get("local_balance", ""):
            issues.append(f"local balance missing phrase: {phrase}")
    for phrase in [
        "r(t,x)=4q t^2/[x^2(1+q)]",
        "1/5776",
        "A_pf>=1+5775/[5776(1+q)^2]>1",
        "28879/23104",
    ]:
        if phrase not in json.dumps(coefficient, sort_keys=True):
            issues.append(f"coefficient audit missing phrase: {phrase}")

    cancellation = exact.get("contact_cancellation", {})
    for phrase in [
        "partial_t E_pf=0",
        "partial_x J_pf=s_pf^2 H_xx^2",
        "exactly canceled",
    ]:
        if phrase not in json.dumps(cancellation, sort_keys=True):
            issues.append(f"contact cancellation missing phrase: {phrase}")
    bridge = exact.get("spatial_gradient_bridge", {})
    for phrase in [
        "partial_x V_pf",
        "B_pf=|partial_x V_pf|^2",
        ">=|partial_x V_pf|^2",
    ]:
        if phrase not in json.dumps(bridge, sort_keys=True):
            issues.append(f"spatial bridge missing phrase: {phrase}")
    for phrase in [
        "[sqrt(E_pf(t,a))+sqrt(E_pf(t,b))]^2/(b-a)",
        "then V_pf(t,x) has no zero",
        "Cauchy-Schwarz",
    ]:
        if phrase not in json.dumps(exact, sort_keys=True):
            issues.append(f"contact criterion missing phrase: {phrase}")
    countermodel = exact.get("quadratic_countermodel", {})
    for phrase in [
        "(x-x_0)^2-2(t-t_0)",
        "H_t=-H_xx",
        "J_x=4s_0^2",
        "cannot exclude",
    ]:
        if phrase not in json.dumps(countermodel, sort_keys=True):
            issues.append(f"countermodel missing phrase: {phrase}")

    nonpromotion = exact.get("nonpromotion", {})
    for phrase in [
        "isolated point",
        "not a closed scalar parabolic inequality",
        "No kappa is defined by dividing",
    ]:
        if phrase not in json.dumps(nonpromotion, sort_keys=True):
            issues.append(f"nonpromotion guard missing phrase: {phrase}")

    ray = load_json(gate.SOURCES["ray_aligned"], "ray source", issues)
    contact = load_json(
        gate.SOURCES["contact_normal"], "contact source", issues
    )
    shell = load_json(gate.SOURCES["q209_shell"], "Q209 source", issues)
    if ray.get("summary", {}).get("old_collar_open_antecedents") != 1:
        issues.append("ray source outer-collar count failed")
    if ray.get("summary", {}).get("new_strip_open_antecedents") != 0:
        issues.append("ray source new-strip count failed")
    if "quadratic_countermodel" not in contact.get("exact", {}):
        issues.append("contact source countermodel failed")
    if shell.get("summary", {}).get("q209_certified") is not False:
        issues.append("Q209 source noncertification failed")

    proof_boundary = stored.get("proof_boundary", "")
    for phrase in [
        "does not exclude a first-jet contact",
        "does not prove Xi endpoint",
        "strict bulk-budget inequality",
        "certify Q209",
        "prove Lambda<=0",
        "prove RH",
        "Clay-prize conclusion",
    ]:
        if phrase not in proof_boundary:
            issues.append(f"proof boundary missing phrase: {phrase}")

    if NOTE.is_file():
        note = NOTE.read_text(encoding="utf-8")
        for phrase in [
            "No term involving `s_t` or `s_x` has been discarded.",
            "28879/23104",
            "Pi provenance:",
            "## Exact Pointwise Bridge",
            "Cauchy-Schwarz",
            "There is no contradiction",
            "positive or monotone spatial integral",
            "0 pointwise contact exclusions",
        ]:
            if phrase not in note:
                issues.append(f"note missing phrase: {phrase}")
    return issues


def main() -> int:
    issues = validate()
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman parabolic-frequency energy/current gate: "
        "18 rows, 8 exact identities, 3 coefficient bounds, "
        "2 pointwise bridge rows, 3 nonpromotion guards, "
        "1 open Xi bulk-margin target, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
