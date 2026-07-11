#!/usr/bin/env python3
"""Independently validate the Newman Gasper residual two-block gate."""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_gasper_residual_two_block_gate as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_gasper_residual_two_block_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_gasper_residual_two_block_gate.md"

EXPECTED_IDS = [
    "grtb_01_one_block_residual_positive",
    "grtb_02_two_block_tail_parameter",
    "grtb_03_laguerre_beta_quadratic",
    "grtb_04_acb_cauchy_derivatives",
    "grtb_05_beta_uniform_negative_laguerre",
    "grtb_06_two_block_obstruction",
    "grtb_07_multiplier_discriminant_guard",
    "grtb_08_signed_multiblock_handoff",
]


def validate() -> list[str]:
    issues: list[str] = []
    if not RESULT.exists():
        return ["missing stored result"]
    if not NOTE.exists():
        return ["missing rendered note"]

    stored = json.loads(RESULT.read_text(encoding="utf-8"))
    rebuilt = target.build_payload()
    if stored != rebuilt:
        issues.append("stored Gasper residual payload differs from reconstruction")
    if stored.get("status") != (
        "exact and interval-certified positive two-block Gasper obstruction"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 8:
        issues.append("expected 8 rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "exact_theorem": 2,
        "exact_identity": 1,
        "interval_method": 1,
        "interval_theorem": 1,
        "exact_interval_composition": 1,
        "exact_route_obstruction": 1,
        "open_handoff": 1,
    }
    for role, expected in expected_roles.items():
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    exact = stored.get("exact", {})

    # Verify the exact beta-quadratic Laguerre identity.
    beta = sp.symbols("beta", real=True)
    f, fp, fpp, g, gp, gpp = sp.symbols("f fp fpp g gp gpp", real=True)
    direct = (fp - beta * gp) ** 2 - (f - beta * g) * (fpp - beta * gpp)
    split = (
        fp**2
        - f * fpp
        - beta * (2 * fp * gp - f * gpp - g * fpp)
        + beta**2 * (gp**2 - g * gpp)
    )
    if sp.expand(direct - split) != 0:
        issues.append("two-block Laguerre quadratic failed")
    laguerre = exact.get("laguerre_quadratic", {})
    if "L[R_beta]=L[E_0]-beta*B" not in laguerre.get("identity", ""):
        issues.append("stored Laguerre quadratic drifted")

    # Check the pointwise first-residual proof and tail coefficient.
    one_block = exact.get("one_block_positive_residual", {})
    for phrase in ("q=exp(-4u)", "pi-a-3/2>1", "Phi(u)>phi_1(u)>Psi_9(u)"):
        text = " ".join(str(value) for value in one_block.values())
        if phrase not in text:
            issues.append(f"one-block positivity proof missing: {phrase}")
    if not sp.N(sp.pi - sp.Rational(3, 2) / sp.pi - sp.Rational(3, 2)) > 1:
        issues.append("one-block elementary coefficient check failed")
    family = exact.get("two_block_family", {})
    if family.get("beta_match") != "beta_tail=pi-3/(2*pi)":
        issues.append("tail-matched beta drifted")
    if "pi-3/(2*pi)-beta" not in family.get("residual_tail", ""):
        issues.append("two-block residual tail coefficient missing")

    # Re-derive the multiplier quartic and its discriminant.
    y = sp.symbols("y")
    quartic = (
        256 * y**4
        - 576 * y**3
        + (432 + 16 * beta) * y**2
        - (120 + 20 * beta) * y
        + 9
        + 5 * beta
    )
    discriminant = sp.factor(sp.discriminant(quartic, y))
    expected_discriminant = -sp.Integer(2) ** 24 * (
        20 * beta**5
        - 64 * beta**4
        + 184 * beta**3
        - 432 * beta**2
        + 972 * beta
        - 729
    )
    if sp.expand(discriminant - expected_discriminant) != 0:
        issues.append("multiplier quartic discriminant failed")
    f_poly = -discriminant / sp.Integer(2) ** 24
    if sp.simplify(f_poly.subs(beta, sp.Rational(11, 10)) - sp.Rational(4459, 5000)) != 0:
        issues.append("multiplier discriminant threshold value failed")
    guard = exact.get("multiplier_guard", {})
    bernstein = [sp.Rational(value) for value in guard.get("bernstein_coefficients", [])]
    if len(bernstein) != 5 or not all(value > 0 for value in bernstein):
        issues.append("multiplier derivative Bernstein certificate failed")
    if "off-imaginary zeros" not in guard.get("conclusion", ""):
        issues.append("universal-multiplier rejection missing")

    # Validate the Acb/Arb Cauchy derivative and beta-cover certificates.
    interval = stored.get("interval_certificate", {})
    if interval.get("rigorous") is not True:
        issues.append("interval certificate is not marked rigorous")
    if interval.get("arithmetic") != "python-flint Acb/Arb balls at 110 decimal digits":
        issues.append("interval arithmetic metadata drifted")
    target.ctx.dps = 110
    step = target.arb(interval.get("step", "nan"))
    cauchy_radius = target.arb(interval.get("cauchy_radius", "nan"))
    box_radius = target.arb(interval.get("acb_box_radius", "nan"))
    if not (step > 0 and cauchy_radius > 0 and step + cauchy_radius < box_radius):
        issues.append("Cauchy circle is not contained in the Acb box")
    method = interval.get("derivative_method", "")
    for phrase in ("h^2*M/r^3", "2*h^2*M/r^4", "Cauchy"):
        if phrase not in method:
            issues.append(f"Cauchy derivative method missing: {phrase}")

    certificate_rows = interval.get("rows", [])
    if [row.get("x") for row in certificate_rows] != ["66", "50"]:
        issues.append("interval spectral points drifted")
    if len(certificate_rows) != 2:
        issues.append("expected two beta-cover interval rows")
    parsed_intervals: list[tuple[target.arb, target.arb]] = []
    for row in certificate_rows:
        try:
            left_beta = target.arb(row["beta_interval"][0])
            right_beta = target.arb(row["beta_interval"][1])
            c0 = target.arb(row["laguerre_quadratic"]["c0"])
            c1 = target.arb(row["laguerre_quadratic"]["c1"])
            c2 = target.arb(row["laguerre_quadratic"]["c2"])
            left_value = target.arb(row["left_endpoint_value"])
            right_value = target.arb(row["right_endpoint_value"])
        except (KeyError, TypeError, ValueError):
            issues.append(f"unparseable interval certificate row at x={row.get('x')}")
            continue
        parsed_intervals.append((left_beta, right_beta))
        recomputed_left = c0 + left_beta * c1 + left_beta**2 * c2
        recomputed_right = c0 + right_beta * c1 + right_beta**2 * c2
        if not (c2 > 0 and left_value < 0 and right_value < 0):
            issues.append(f"beta-uniform convexity certificate failed at x={row.get('x')}")
        if not (recomputed_left < 0 and recomputed_right < 0):
            issues.append(f"recomputed beta endpoints failed at x={row.get('x')}")
        if not all(row.get("certified", {}).values()):
            issues.append(f"stored interval flags failed at x={row.get('x')}")

    if len(parsed_intervals) == 2:
        first, second = parsed_intervals
        beta_match = target.arb(interval.get("beta_match", "nan"))
        if not (first[0] == 0 and first[1].overlaps(second[0]) and second[1].overlaps(beta_match)):
            issues.append("beta intervals do not cover [0,beta_match]")

    obstruction = exact.get("two_block_obstruction", {})
    for phrase in ("every beta>=0", "eventually negative", "L[R_beta](66)<0"):
        joined = " ".join(str(value) for value in obstruction.values())
        if phrase not in joined:
            issues.append(f"two-block composition missing: {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    required_note = (
        "Phi(u)>phi_1(u)>Psi_9(u)",
        "x=66",
        "x=50",
        "negative throughout its closed beta interval",
        "signed, coupled-square problem",
        "This is not a proof or disproof of RH",
    )
    for phrase in required_note:
        if phrase not in note:
            issues.append(f"rendered note missing: {phrase}")
    for forbidden in ("therefore RH", "we have disproved RH", "Lambda <= 0 is proved"):
        if forbidden.lower() in note.lower():
            issues.append(f"rendered note contains forbidden promotion: {forbidden}")
    return issues


def main() -> int:
    issues = validate()
    for issue in issues:
        print(f"GASPER-TWO-BLOCK [{issue}]")
    print(
        "validated Jensen-window PF Newman Gasper residual two-block gate: "
        f"8 rows, {len(issues)} issues, 2 exact kernel theorems, "
        "1 exact Laguerre quadratic, 2 Acb derivative certificates, "
        "2 beta intervals covered, 1 exhaustive positive-residual obstruction, "
        "1 multiplier guard, 1 signed handoff"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
