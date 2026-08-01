#!/usr/bin/env python3
"""Validate the exact one-sided Newman phase/moment bridge gate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_one_sided_phase_moment_bridge_gate as gate


EXPECTED_IDS = [
    "nospmb_01_imported_shape",
    "nospmb_02_score_probability",
    "nospmb_03_zero_free_lift",
    "nospmb_04_phase_contact",
    "nospmb_05_scaled_first_jet",
    "nospmb_06_signed_hankel_bridge",
    "nospmb_06a_score_beta_abel_measure",
    "nospmb_06b_full_scale_interlacing_guard",
    "nospmb_06c_abel_kernel_sign_guard",
    "nospmb_07_score_flow",
    "nospmb_08_normalizer_transfer",
    "nospmb_09_phase_conditioning",
    "nospmb_10_shape_countermodel",
    "nospmb_11_xi_diagnostics",
    "nospmb_12_live_target",
    "nospmb_13_nonpromotion",
]


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if (
        artifact.get("kind")
        != "jensen_window_pf_newman_one_sided_phase_moment_bridge_gate"
    ):
        issues.append("artifact kind mismatch")
    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    exact = gate.build_exact()
    if artifact.get("exact") != exact:
        issues.append("exact payload drifted")

    roles = {
        role: sum(row.get("role") == role for row in rows)
        for role in {
            "imported_exact_theorem",
            "exact_probability_identity",
            "exact_complex_identity",
            "exact_contact_equivalence",
            "exact_target_reduction",
            "exact_coefficient_bridge",
            "exact_abel_measure_bridge",
            "exact_interlacing_countergate",
            "exact_abel_kernel_countergate",
            "exact_dynamical_identity",
            "exact_quantitative_transfer",
            "exact_asymptotic_guard",
            "exact_countermodel",
            "finite_diagnostics",
            "open_theorem_target",
            "nonpromotion_gate",
        }
    }
    if any(count != 1 for count in roles.values()):
        issues.append(f"role counts drifted: {roles}")

    open_row = next(
        (row for row in rows if row.get("role") == "open_theorem_target"),
        {},
    )
    if open_row.get("readiness") != "not_ready_to_apply":
        issues.append("Xi joint-avoidance target was promoted")
    countermodel = next(
        (row for row in rows if row.get("role") == "exact_countermodel"),
        {},
    )
    if countermodel.get("readiness") != "guard_validated":
        issues.append("shape countermodel is not active")

    xi = sp.symbols("xi", real=True)
    transform = (
        8
        * sp.sqrt(2 * sp.pi)
        * sp.exp(-2 * xi**2)
        * sp.sin(xi / 2) ** 2
        / xi**2
    )
    point = 2 * sp.pi
    if sp.simplify(transform.subs(xi, point)) != 0:
        issues.append("countermodel transform does not vanish at 2*pi")
    if sp.simplify(sp.diff(transform, xi).subs(xi, point)) != 0:
        issues.append("countermodel transform is not double at 2*pi")
    quadratic = sp.simplify(
        sp.limit(transform / (xi - point) ** 2, xi, point)
    )
    expected_quadratic = (
        sp.sqrt(2 * sp.pi) * sp.exp(-8 * sp.pi**2) / (2 * sp.pi**2)
    )
    if sp.simplify(quadratic - expected_quadratic) != 0:
        issues.append("countermodel quadratic coefficient drifted")

    w = sp.symbols("w")
    for degree in range(1, 5):
        for shift in range(4):
            moments = sp.symbols(
                f"M0:{shift + degree + 1}", positive=True
            )
            direct = sp.Add(
                *[
                    sp.binomial(degree, index)
                    * moments[shift + index]
                    * w**index
                    / sp.factorial(shift + index)
                    for index in range(degree + 1)
                ]
            )
            laguerre = sp.expand(
                sp.assoc_laguerre(degree, shift, -w)
            )
            mixed = sp.Add(
                *[
                    laguerre.coeff(w, power)
                    * moments[shift + power]
                    * w**power
                    for power in range(degree + 1)
                ]
            )
            mixed = sp.expand(
                sp.factorial(degree)
                * mixed
                / sp.factorial(degree + shift)
            )
            if sp.simplify(direct - mixed) != 0:
                issues.append(
                    "score-Beta Abel shifted-Laguerre identity failed "
                    f"at D={degree}, n={shift}"
                )
    for index in range(5):
        next_moment = sp.symbols(f"N_{index + 1}", positive=True)
        flow_moment = (
            4 * next_moment
            - 2 * next_moment / (index + 1)
        )
        expected_flow = (
            2
            * (2 * index + 1)
            * next_moment
            / (index + 1)
        )
        if sp.simplify(flow_moment - expected_flow) != 0:
            issues.append(
                f"score-Beta Abel flow moment failed at k={index}"
            )
    for shift in range(6):
        lower_root = shift + 2 - sp.sqrt(shift + 2)
        upper_root = shift + 2 + sp.sqrt(shift + 2)
        if sp.ask(sp.Q.positive(10 * lower_root - upper_root)) is not True:
            issues.append(
                "two-scale degree-two interlacing guard failed "
                f"at n={shift}"
            )
    abel_negative_minor = sp.Rational(1, 2) - 1 / sp.sqrt(3)
    if sp.ask(sp.Q.negative(abel_negative_minor)) is not True:
        issues.append("negative Abel-kernel minor guard failed")
    abel_triangular_minor = sp.det(
        sp.Matrix([[1, 1 / sp.sqrt(3)], [0, 1]])
    )
    if abel_triangular_minor != 1:
        issues.append("positive Abel-kernel minor guard failed")

    fresh_coarse = gate.diagnostics(gate.COARSE_DPS)
    fresh_fine = gate.diagnostics(gate.FINE_DPS)
    fresh_convergence = gate.compare_diagnostics(
        fresh_coarse, fresh_fine
    )
    if artifact.get("diagnostics") != fresh_fine:
        issues.append("stored one-sided diagnostics drifted")
    if artifact.get("convergence") != fresh_convergence:
        issues.append("stored convergence audit drifted")
    diagnostic_rows = fresh_fine.get("rows", [])
    if len(diagnostic_rows) != 3:
        issues.append("expected three one-sided diagnostic rows")
    if any(mp.mpf(row["F_imag"]) <= 0 for row in diagnostic_rows):
        issues.append("one-sided imaginary lift is not strictly positive")
    residual_fields = (
        "complex_lift_abs_residual",
        "H_characteristic_abs_residual",
        "H_first_characteristic_abs_residual",
        "phase_formula_abs_residual",
    )
    if any(
        mp.mpf(row[field]) >= mp.mpf("1e-30")
        for row in diagnostic_rows
        for field in residual_fields
    ):
        issues.append("one-sided identity residual lost precision")
    first_crossing = next(
        (
            row
            for row in diagnostic_rows
            if row.get("label") == "first_xi_crossing"
        ),
        {},
    )
    if not first_crossing:
        issues.append("first Xi crossing diagnostic missing")
    elif abs(mp.mpf(first_crossing["B"])) >= mp.mpf("1e-25"):
        issues.append("first Xi crossing characteristic residual too large")

    status = artifact.get("status", "")
    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "zero-free phase lift",
        "signed-Hankel",
        "joint-avoidance theorem remains open",
    ):
        if marker not in status:
            issues.append(f"status marker missing: {marker}")
    for marker in (
        "does not prove",
        "PF-infinity",
        "Lambda<=0",
        "RH",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = gate.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Score Probability",
        "Zero-Free Complex Lift",
        "Exact Contact System",
        "Signed-Hankel Coordinate",
        "Score-Beta Abel Measure",
        "Interlacing Guard",
        "bare Abel operator",
        "Closed Score Flow",
        "Normalizer-Compatible Jet",
        "Conditioning Guard",
        "Shape Countermodel",
        "Live Target",
        "not a proof of `Lambda <= 0` or RH",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact", type=Path, default=gate.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        raise SystemExit(1)
    artifact = json.loads(args.artifact.read_text(encoding="utf-8"))
    print(
        "validated Newman one-sided phase/moment bridge gate: "
        f"{len(artifact['rows'])} rows, 0 issues, "
        "1 zero-free complex lift, 1 score probability, "
        "1 exact phase-contact system, 1 signed-Hankel moment bridge, "
        "1 score/Beta Abel measure, 1 full-scale interlacing guard, "
        "1 Abel-kernel sign guard, 1 closed score flow, "
        "1 normalizer transfer, 1 conditioning guard, "
        "1 shape countermodel, "
        "3 precision-stable diagnostics, 1 open Xi joint-avoidance target"
    )


if __name__ == "__main__":
    main()
