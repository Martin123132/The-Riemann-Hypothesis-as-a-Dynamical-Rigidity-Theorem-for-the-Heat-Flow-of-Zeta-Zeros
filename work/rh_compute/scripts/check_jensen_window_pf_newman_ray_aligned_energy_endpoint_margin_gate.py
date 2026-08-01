#!/usr/bin/env python3
"""Independently validate the ray-aligned energy endpoint-margin gate."""

from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate as gate


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_ray_aligned_energy_endpoint_margin_gate"
RESULT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
EXPECTED_IDS = [
    "raem_01_active_schedule",
    "raem_02_compact_edge_inventory",
    "raem_03_compact_ratio",
    "raem_04_left_full_margin",
    "raem_05_left_energy_conversion",
    "raem_06_left_stage_margin",
    "raem_07_right_curvature",
    "raem_08_right_hxx_cap",
    "raem_09_right_beta_floor",
    "raem_10_right_energy_margin",
    "raem_11_two_endpoint_floor",
    "raem_12_strict_budget_target",
    "raem_13_floor_conditioning",
    "raem_14_dominant_coverage",
    "raem_15_interior_numerator_gap",
    "raem_16_tautology_guard",
    "raem_17_normalizer_and_finite_guards",
    "raem_18_localization_handoff",
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


def validate_compact_edge(issues: list[str]) -> None:
    compact = load_json(gate.SOURCES["compact"], "compact source", issues)
    records = compact.get("certificate", {}).get("records", [])
    edge = [row for row in records if row.get("x_high") == "38"]
    if len(edge) != 10:
        issues.append("independent compact edge count failed")
        return
    if any(row.get("branch") != "derivative" for row in edge):
        issues.append("independent compact branch audit failed")
    if any(
        gate.ball_lower(row.get("certified_ratio_lower", ""))
        <= Decimal("1.75")
        for row in edge
    ):
        issues.append("independent compact ratio floor failed")
    intervals = sorted(
        (Fraction(row["t_low"]), Fraction(row["t_high"]))
        for row in edge
    )
    expected = [
        (Fraction(k, 50), Fraction(k + 1, 50))
        for k in range(10)
    ]
    if intervals != expected:
        issues.append("independent compact edge tiling failed")


def validate_left_margin(issues: list[str]) -> None:
    x = 38
    ratio_excess = Fraction(7, 4) - 1
    delta_0 = Fraction(1, 2800)
    margin = ratio_excess * 4 * x**3 * delta_0
    if margin != Fraction(20577, 350):
        issues.append("independent left raw margin failed")

    H, Hx, s = sp.symbols("H Hx s", real=True, positive=True)
    functional = 64 * x**3 * H + 16 * x**4 * Hx
    energy = H**2 + s**2 * Hx**2
    norm_squared = (64 * x**3) ** 2 + (16 * x**4 / s) ** 2
    u = sp.Matrix([H, s * Hx])
    c = sp.Matrix([64 * x**3, 16 * x**4 / s])
    if sp.expand(functional - (c.dot(u))) != 0:
        issues.append("independent left functional representation failed")
    gram_residual = sp.expand(
        norm_squared * energy - functional**2
    )
    if sp.factor(gram_residual) != (
        16 * x**3 * H / s - 64 * x**3 * s * Hx
    ) ** 2:
        # Compare algebraically because SymPy may choose another square form.
        if sp.simplify(
            gram_residual
            - (
                (64 * x**3) * s * Hx
                - (16 * x**4 / s) * H
            ) ** 2
        ) != 0:
            issues.append("independent weighted Cauchy inequality failed")


def validate_right_margin(issues: list[str]) -> None:
    z0 = Fraction(282, 125) + Fraction(1, 8000)
    z1 = Fraction(59, 100) + Fraction(7, 160000)
    z2 = Fraction(77, 500) + Fraction(13, 800000)
    ell = Fraction(27, 100)
    cap = z2 + 2 * ell * z1 + (
        ell**2 + Fraction(1, 10**23)
    ) * z0
    if cap >= Fraction(13, 20):
        issues.append("independent normalized Hxx cap failed")
    beta2 = Fraction(6300, 6301)
    y = Fraction(99, 1000)
    upper = y**2 / beta2 + Fraction(13, 20) * y
    if Fraction(3, 40) - upper != Fraction(593211, 700000000):
        issues.append("independent right margin gap failed")
    if upper >= Fraction(3, 40):
        issues.append("independent right margin is not strict")


def validate_asymptotic(issues: list[str]) -> None:
    x = sp.Integer(38)
    ell = sp.symbols("ell", positive=True)
    margin = sp.Rational(20577, 350)
    l0 = sp.log(x / (4 * sp.pi))
    left_energy = margin**2 / (
        4096 * x**6
        + 256 * x**8 * (l0**2 + ell / 50)
    )
    width = 4 * sp.pi * sp.exp(ell) - 38
    limit = sp.simplify(
        sp.limit(ell * sp.exp(ell) * left_energy / width, ell, sp.oo)
    )
    expected = sp.Rational(9, 231853260800) / sp.pi
    if sp.simplify(limit - expected) != 0:
        issues.append("independent floor asymptotic failed")


def validate_diagnostics(stored: dict, issues: list[str]) -> None:
    rows = stored.get("diagnostics", [])
    if [row.get("stage") for row in rows] != [25, 100, 1000, 10000]:
        issues.append("diagnostic stages drifted")
        return
    x = 38.0
    l0 = math.log(x / (4 * math.pi))
    margin = 3 * x**3 / 2800
    for row in rows:
        ell = 101 + row["stage"]
        t = 25 / ell
        inv_s2 = l0**2 + 1 / (2 * t)
        energy = margin**2 / (
            4096 * x**6 + 256 * x**8 * inv_s2
        )
        actual = float(row["certified_left_energy_lower"])
        if not math.isclose(actual, energy, rel_tol=1e-14, abs_tol=0):
            issues.append(f"diagnostic left energy failed: {row['stage']}")


def validate() -> list[str]:
    issues: list[str] = []
    stored = load_json(RESULT, "stored endpoint-margin result", issues)
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
        "exact ray-aligned endpoint margins and whole-collar "
        "bulk-budget conditioning audit"
    ):
        issues.append("status drifted")
    rows = stored.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_summary = {
        "rows": 18,
        "compact_edge_records": 10,
        "left_raw_margins": 1,
        "right_normalized_margins": 1,
        "exact_contact_floors": 1,
        "open_full_collar_bulk_budgets": 1,
        "nonpromotion_guards": 2,
        "pointwise_contact_exclusions": 0,
        "cofinal_descendant_theorem": False,
    }
    if stored.get("summary") != expected_summary:
        issues.append("summary drifted")

    validate_compact_edge(issues)
    validate_left_margin(issues)
    validate_right_margin(issues)
    validate_asymptotic(issues)
    validate_diagnostics(stored, issues)

    exact = stored.get("exact", {})
    for phrase in [
        "|J_t'(38)|>3*38^3/2800=20577/350",
        "E_pf(t,38)>M_38^2/",
        "sqrt(E_pf(t,R_j))>(99/1000)A_t(R_j)",
        "593211/700000000",
        "exp(-L_j)/L_j",
    ]:
        if phrase not in json.dumps(exact, sort_keys=True):
            issues.append(f"exact endpoint audit missing phrase: {phrase}")
    for phrase in [
        "x_dom(t_(j+1))=R_j",
        "x_dom(t_j)=R_(j-1)",
        "[38,x_dom(t))",
    ]:
        if phrase not in json.dumps(
            exact.get("dominant_coverage", {}), sort_keys=True
        ):
            issues.append(f"dominant coverage missing phrase: {phrase}")
    for phrase in [
        "identity, not the strict",
        "A_t cannot be dropped",
        "finite evaluation",
    ]:
        if phrase not in json.dumps(
            exact.get("nonpromotion", {}), sort_keys=True
        ):
            issues.append(f"nonpromotion guard missing phrase: {phrase}")

    energy = load_json(
        gate.SOURCES["energy_current"], "energy source", issues
    )
    global_ray = load_json(
        gate.SOURCES["dominant_global"], "dominant source", issues
    )
    if energy.get("summary", {}).get("pointwise_bridge_rows") != 2:
        issues.append("energy source bridge count failed")
    if "3/40*L^2" not in global_ray.get("exact", {}).get("theorem", ""):
        issues.append("dominant source curvature failed")

    proof_boundary = stored.get("proof_boundary", "")
    for phrase in [
        "does not prove a full-collar Xi bulk upper bound",
        "strict no-contact budget",
        "Q209",
        "Lambda<=0",
        "RH",
        "Clay-prize conclusion",
    ]:
        if phrase not in proof_boundary:
            issues.append(f"proof boundary missing phrase: {phrase}")

    if NOTE.is_file():
        note = NOTE.read_text(encoding="utf-8")
        for phrase in [
            "20577/350",
            "99/1000",
            "A_t(R_j)",
            "That upper budget is not currently proved.",
            "conditioning diagnostic, not a",
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
        "validated ray-aligned energy endpoint-margin gate: 18 rows, "
        "10 compact edge records, 1 left raw margin, "
        "1 right normalized margin, 1 exact contact floor, "
        "1 open full-collar bulk budget, 2 nonpromotion guards, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
