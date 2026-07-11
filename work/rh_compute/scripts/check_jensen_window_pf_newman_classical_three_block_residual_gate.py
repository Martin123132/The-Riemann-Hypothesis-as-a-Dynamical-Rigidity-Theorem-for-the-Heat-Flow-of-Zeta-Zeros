#!/usr/bin/env python3
"""Independently validate the classical three-block Gasper residual gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_classical_three_block_residual_gate as target


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_classical_three_block_residual_gate.json"
)
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_classical_three_block_residual_gate.md"

EXPECTED_IDS = [
    "gctb_01_classical_normalizations",
    "gctb_02_polya_p2_positive_residual",
    "gctb_03_de_bruijn_positive_residual",
    "gctb_04_tail_origin_compactness",
    "gctb_05_bivariate_laguerre_identity",
    "gctb_06_acb_cauchy_jets",
    "gctb_07_positive_three_block_mesh",
    "gctb_08_classical_residual_obstructions",
    "gctb_09_gasper_square_scope",
    "gctb_10_coupled_signed_handoff",
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
        issues.append("stored classical three-block payload differs from reconstruction")
    if stored.get("status") != (
        "exact and interval-certified classical three-block residual obstruction"
    ):
        issues.append("status drifted")

    rows = stored.get("rows", [])
    if len(rows) != 10:
        issues.append("expected 10 rows")
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row ids or ordering drifted")
    expected_roles = {
        "established_benchmark": 1,
        "exact_kernel_theorem": 1,
        "interval_kernel_theorem": 1,
        "exact_interval_reduction": 1,
        "exact_identity": 1,
        "interval_method": 1,
        "interval_theorem": 1,
        "interval_composition": 1,
        "source_scope_guard": 1,
        "open_handoff": 1,
    }
    for role, expected in expected_roles.items():
        if sum(row.get("role") == role for row in rows) != expected:
            issues.append(f"expected {expected} row(s) with role {role}")

    required_sources = {
        "https://arxiv.org/abs/0801.2996",
        "https://arxiv.org/abs/1502.06844",
        "https://doi.org/10.1215/S0012-7094-50-01720-0",
        "https://doi.org/10.1007/BF02565336",
    }
    if not required_sources.issubset(set(stored.get("sources", []))):
        issues.append("classical primary-source set drifted")

    exact = stored.get("exact", {})
    benchmarks = exact.get("classical_benchmarks", {})
    if benchmarks.get("polya_p2", {}).get("transform") != (
        "B_P2=P_(9/4)-a*P_(5/4)"
    ):
        issues.append("Polya P2 normalization drifted")
    if benchmarks.get("de_bruijn", {}).get("transform") != (
        "B_dB=P_(9/4)+b*P_(5/4)+P_(1/4)"
    ):
        issues.append("de Bruijn normalization drifted")

    # Re-derive the P2 lower bound in q=r^2.
    r = sp.symbols("r", positive=True)
    q = r**2
    pi = sp.pi
    a = sp.Rational(3, 2) / pi
    b = pi - a
    p2_block = 1 + r**9 - a * q * (1 + r**5)
    p2_lower = sp.factor((1 + pi * q) * (1 - a * q) - p2_block)
    p2_expected = q * (pi - sp.Rational(3, 2) * q - r**7 + a * r**5)
    if sp.expand(p2_lower - p2_expected) != 0:
        issues.append("Polya P2 pointwise lower bound failed")
    p2 = exact.get("polya_p2_positive_residual", {})
    if "q*(pi-5/2)>0" not in p2.get("proof", ""):
        issues.append("Polya P2 positivity proof drifted")

    # Independently rebuild the de Bruijn Bernstein coefficients.
    db_block = 1 + r**9 + b * q * (1 + r**5) + q**2 * (1 + r)
    taylor = sum((pi * q) ** k / sp.factorial(k) for k in range(5))
    db_poly = sp.Poly(sp.cancel(((1 - a * q) * taylor - db_block) / q**2), r)
    db_bernstein = target.bernstein_coefficients(db_poly)
    db = exact.get("de_bruijn_positive_residual", {})
    if [sp.sstr(value) for value in db_bernstein] != db.get(
        "bernstein_coefficients", []
    ):
        issues.append("de Bruijn Bernstein coefficients drifted")
    try:
        db_balls = [target.arb(value) for value in db.get("bernstein_balls", [])]
    except (TypeError, ValueError):
        db_balls = []
        issues.append("de Bruijn Bernstein balls are unparseable")
    if len(db_balls) != 7 or not all(value > 0 for value in db_balls):
        issues.append("de Bruijn Bernstein positivity failed")

    # Re-derive the exact bivariate Laguerre identity.
    beta, gamma = sp.symbols("beta gamma", real=True)
    e, ep, epp = sp.symbols("e ep epp", real=True)
    p, pp, ppp = sp.symbols("p pp ppp", real=True)
    s, sp1, spp = sp.symbols("s sp spp", real=True)
    direct = (ep - beta * pp - gamma * sp1) ** 2 - (
        e - beta * p - gamma * s
    ) * (epp - beta * ppp - gamma * spp)
    split = (
        ep**2
        - e * epp
        - beta * (2 * ep * pp - e * ppp - p * epp)
        - gamma * (2 * ep * sp1 - e * spp - s * epp)
        + beta**2 * (pp**2 - p * ppp)
        + beta * gamma * (2 * pp * sp1 - p * spp - s * ppp)
        + gamma**2 * (sp1**2 - s * spp)
    )
    if sp.expand(direct - split) != 0:
        issues.append("bivariate Laguerre identity failed")
    laguerre = exact.get("laguerre_polynomial", {})
    for phrase in ("-beta*B[E,P5]", "beta*gamma*B[P5,P1]", "gamma^2*L[P1]"):
        if phrase not in laguerre.get("identity", ""):
            issues.append(f"stored bivariate identity missing: {phrase}")

    origin = stored.get("origin_certificate", {})
    if origin.get("rigorous") is not True:
        issues.append("origin certificate is not marked rigorous")
    if origin.get("partial_cutoff") != 8:
        issues.append("origin partial cutoff drifted")
    try:
        phi_upper = target.arb(origin.get("phi_upper", "nan"))
        comparison_upper = target.arb(origin.get("comparison_upper", "nan"))
        tail_ratio = target.arb(origin.get("tail_ratio", "nan"))
    except (TypeError, ValueError):
        phi_upper = comparison_upper = tail_ratio = target.arb("nan")
        issues.append("origin certificate balls are unparseable")
    if not (tail_ratio > 0 and tail_ratio < 1 and phi_upper < comparison_upper):
        issues.append("origin compactness inequality failed")
    parameter_region = origin.get("parameter_region", "")
    for phrase in ("0<=beta<27/10", "beta+gamma<51/10"):
        if phrase not in parameter_region:
            issues.append(f"parameter region missing: {phrase}")

    interval = stored.get("interval_certificate", {})
    if interval.get("rigorous") is not True:
        issues.append("spectral interval certificate is not marked rigorous")
    if interval.get("arithmetic") != (
        "python-flint Acb/Arb balls at 130 decimal digits"
    ):
        issues.append("spectral arithmetic metadata drifted")
    if interval.get("spectral_points") != ["48", "52", "86"]:
        issues.append("spectral points drifted")
    target.ctx.dps = 130
    step = target.arb(interval.get("step", "nan"))
    cauchy_radius = target.arb(interval.get("cauchy_radius", "nan"))
    box_radius = target.arb(interval.get("acb_box_radius", "nan"))
    if not (
        step > 0
        and cauchy_radius > 0
        and step + cauchy_radius < box_radius
    ):
        issues.append("Cauchy circles are not contained in the Acb boxes")
    for phrase in ("h^2*M/r^3", "2*h^2*M/r^4", "Cauchy"):
        if phrase not in interval.get("derivative_method", ""):
            issues.append(f"Cauchy derivative method missing: {phrase}")

    coefficient_names = (
        "c00",
        "c_beta",
        "c_gamma",
        "c_beta2",
        "c_beta_gamma",
        "c_gamma2",
    )
    coefficients_by_x: dict[int, tuple[target.arb, ...]] = {}
    certificate_rows = interval.get("rows", [])
    if [row.get("x") for row in certificate_rows] != ["48", "52", "86"]:
        issues.append("spectral coefficient rows drifted")
    for row in certificate_rows:
        try:
            x = int(row["x"])
            jets = {
                key: [target.arb(value) for value in row["jets"][key]]
                for key in ("E", "P5", "P1")
            }
            coefficients = tuple(
                target.arb(row["laguerre_coefficients"][name])
                for name in coefficient_names
            )
        except (KeyError, TypeError, ValueError):
            issues.append(f"unparseable spectral row at x={row.get('x')}")
            continue
        coefficients_by_x[x] = coefficients
        recomputed = target.laguerre_coefficients(
            jets["E"], jets["P5"], jets["P1"]
        )
        if not all(left.overlaps(right) for left, right in zip(coefficients, recomputed)):
            issues.append(f"Laguerre coefficients do not match stored jets at x={x}")

    mesh = interval.get("mesh", {})
    if mesh.get("denominator") != 80:
        issues.append("mesh denominator drifted")
    if mesh.get("beta_cell_count") != 216 or mesh.get("sum_cell_count") != 408:
        issues.append("mesh dimensions drifted")
    if mesh.get("total_boxes") != 64908:
        issues.append("mesh total drifted")
    if mesh.get("assignment_order") != ["48", "52", "86"]:
        issues.append("mesh assignment order drifted")
    if mesh.get("assignment_counts") != {"48": 12929, "52": 5833, "86": 46146}:
        issues.append("stored mesh assignment-count regression")

    if set(coefficients_by_x) == {48, 52, 86}:
        denominator = 80
        radius = target.arb(1) / (2 * denominator)
        assignment_counts = {48: 0, 52: 0, 86: 0}
        digest = hashlib.sha256()
        total = 0
        for beta_index in range(216):
            beta_midpoint = target.arb(2 * beta_index + 1) / (2 * denominator)
            beta_box = target.arb(beta_midpoint, radius)
            for gamma_index in range(408 - beta_index):
                gamma_midpoint = target.arb(2 * gamma_index + 1) / (2 * denominator)
                gamma_box = target.arb(gamma_midpoint, radius)
                total += 1
                for x in (48, 52, 86):
                    value = target.evaluate_laguerre_polynomial(
                        coefficients_by_x[x], beta_box, gamma_box
                    )
                    if value < 0:
                        assignment_counts[x] += 1
                        digest.update(
                            f"{beta_index},{gamma_index}:{x}\n".encode("ascii")
                        )
                        break
                else:
                    issues.append(
                        f"uncovered parameter box beta={beta_index}, gamma={gamma_index}"
                    )
                    break
            if issues and issues[-1].startswith("uncovered parameter box"):
                break
        if total != 64908:
            issues.append("rechecked mesh did not traverse 64908 boxes")
        stored_counts = {
            int(key): value for key, value in mesh.get("assignment_counts", {}).items()
        }
        if assignment_counts != stored_counts:
            issues.append("mesh assignment counts drifted")
        if digest.hexdigest() != mesh.get("assignment_sha256"):
            issues.append("mesh assignment digest drifted")

        a_ball = 3 / (2 * target.arb.pi())
        b_ball = target.arb.pi() - a_ball
        p2_value = target.evaluate_laguerre_polynomial(
            coefficients_by_x[86], -a_ball, target.arb(0)
        )
        db_value = target.evaluate_laguerre_polynomial(
            coefficients_by_x[52], b_ball, target.arb(1)
        )
        if not (p2_value < 0 and db_value < 0):
            issues.append("classical residual Laguerre witnesses failed")

    witnesses = interval.get("classical_witnesses", {})
    if witnesses.get("polya_p2", {}).get("x") != "86":
        issues.append("Polya P2 witness point drifted")
    if witnesses.get("de_bruijn", {}).get("x") != "52":
        issues.append("de Bruijn witness point drifted")
    for name in ("polya_p2", "de_bruijn"):
        try:
            value = target.arb(witnesses[name]["laguerre_value"])
        except (KeyError, TypeError, ValueError):
            issues.append(f"unparseable classical witness: {name}")
            continue
        if not value < 0 or witnesses[name].get("strictly_negative") is not True:
            issues.append(f"classical witness is not strictly negative: {name}")

    gasper = exact.get("gasper_square_scope", {})
    for phrase in ("two Bessel moduli", "mixed products", "single-shift"):
        joined = " ".join(str(value) for value in gasper.values())
        if phrase not in joined:
            issues.append(f"Gasper square-scope audit missing: {phrase}")

    note = NOTE.read_text(encoding="utf-8")
    required_note = (
        "Phi(u)-Psi_P2(u)>0",
        "Phi(u)-Psi_dB(u)>0",
        "64908",
        "x=48 assignments=12929",
        "x=52 assignments=5833",
        "x=86 assignments=46146",
        "mixed products",
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
        print(f"CLASSICAL-THREE-BLOCK [{issue}]")
    print(
        "validated Jensen-window PF Newman classical three-block residual gate: "
        f"10 rows, {len(issues)} issues, 2 established classical real-zero "
        "benchmarks, 2 positive-kernel residual theorems, 1 compact parameter "
        "reduction, 1 exact bivariate Laguerre identity, 3 Acb spectral "
        "certificates, 64908 parameter boxes covered, 2 classical residual "
        "obstructions, 1 Gasper square-scope guard, 1 coupled handoff"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
