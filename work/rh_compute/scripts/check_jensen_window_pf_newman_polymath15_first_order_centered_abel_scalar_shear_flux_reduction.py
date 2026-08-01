#!/usr/bin/env python3
"""Validate the endpoint-complete Abel-scalar shear-flux reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_abel_scalar_shear_flux_reduction as builder


EXPECTED_IDS = [
    "assfr_00_pi_provenance",
    "assfr_01_slope_shear",
    "assfr_02_gap_nonvanishing",
    "assfr_03_shear_homotopy",
    "assfr_04_positive_scale",
    "assfr_05_reduced_proxy",
    "assfr_06_general_flux",
    "assfr_07_crossing_speed",
    "assfr_08_orientation_handoff",
    "assfr_09_zero_fibre",
    "assfr_10_chart_boundary",
    "assfr_11_heat_family",
    "assfr_12_heat_gap",
    "assfr_13_heat_winding",
    "assfr_14_nonpromotion",
    "assfr_15_pointwise_target",
    "assfr_16_flux_target",
    "assfr_17_q_lt_1_target",
]


def independent_shear_flux(issues: list[str]) -> None:
    x_value, scalar, alpha, ell, homotopy = sp.symbols(
        "x c alpha ell s", real=True
    )
    x_rate, scalar_rate, ell_rate = sp.symbols(
        "x_dot c_dot ell_dot", real=True
    )
    matrix = sp.Matrix(
        [
            [1, 0],
            [homotopy * alpha / ell, 1],
        ]
    )
    if sp.simplify(matrix.det() - 1) != 0:
        issues.append("independent shear determinant failed")

    proxy = x_value + sp.I * scalar / ell
    proxy_rate = x_rate + sp.I * (
        scalar_rate / ell - scalar * ell_rate / ell**2
    )
    current = sp.simplify(
        sp.im(sp.conjugate(proxy) * proxy_rate)
        / sp.re(sp.conjugate(proxy) * proxy)
    )
    expected = (
        ell * x_value * scalar_rate
        - x_value * scalar * ell_rate
        - ell * scalar * x_rate
    ) / (ell**2 * x_value**2 + scalar**2)
    if sp.simplify(current - expected) != 0:
        issues.append("independent Abel phase current failed")
    if sp.simplify(
        current.subs(x_value, 0) + ell * x_rate / scalar
    ) != 0:
        issues.append("independent crossing speed failed")


def independent_orientation_budget(issues: list[str]) -> None:
    level = 50.0
    terms = [
        math.exp(-1.5 * level) / (4.0 * level),
        math.exp(-level) / level,
        101.0
        * math.exp(-0.5 * level)
        / (100000.0 * level + 1.0),
        1.0e-7 / (100000.0 * level + 1.0),
    ]
    total = sum(terms)
    if not total < 2.03e-14:
        issues.append("independent orientation total failed")
    for term in terms:
        if not term > 0:
            issues.append("orientation term is not positive")
    for larger_level in (51.0, 75.0, 100.0):
        larger_terms = [
            math.exp(-1.5 * larger_level) / (4.0 * larger_level),
            math.exp(-larger_level) / larger_level,
            101.0
            * math.exp(-0.5 * larger_level)
            / (100000.0 * larger_level + 1.0),
            1.0e-7 / (100000.0 * larger_level + 1.0),
        ]
        if not sum(larger_terms) < total:
            issues.append("orientation bound is not decreasing")


def independent_heat_guard(issues: list[str]) -> None:
    x_value, time = sp.symbols("x t", real=True)
    for frequency in (1, 2, 5, 9):
        field = sp.exp(frequency**2 * time) * sp.sin(
            frequency * x_value + sp.pi / 4
        )
        slope = sp.diff(field, x_value)
        if sp.simplify(
            sp.diff(field, time) + sp.diff(field, x_value, 2)
        ) != 0:
            issues.append(
                f"independent heat equation failed at m={frequency}"
            )
        numerator = sp.trigsimp(
            field * sp.diff(slope, x_value)
            - slope * sp.diff(field, x_value)
        )
        expected = -frequency**2 * sp.exp(
            2 * frequency**2 * time
        )
        if sp.simplify(numerator - expected) != 0:
            issues.append(
                f"independent heat flux failed at m={frequency}"
            )

        zeros = []
        upward = []
        for index in range(-2, 2 * frequency + 3):
            location = (
                index * math.pi - math.pi / 4
            ) / frequency
            if 0.0 <= location < 2.0 * math.pi:
                zeros.append(location)
                if index % 2 == 0:
                    upward.append(location)
        if len(zeros) != 2 * frequency:
            issues.append(
                f"independent zero count failed at m={frequency}"
            )
        if len(upward) != frequency:
            issues.append(
                f"independent upward count failed at m={frequency}"
            )


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    if artifact.get("source_sha256") != builder.source_hashes():
        issues.append("source hash chain drifted")
    if artifact.get("source_audit") != builder.source_audit():
        issues.append("source audit drifted")
    if artifact.get("exact") != builder.build_exact():
        issues.append("exact payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if len([row for row in rows if row.get("role") == "countermodel"]) != 1:
        issues.append("countermodel-family count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "countermodel_property"
        ]
    ) != 2:
        issues.append("countermodel-property count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "countermodel_consequence"
        ]
    ) != 1:
        issues.append("countermodel-consequence count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "open_theorem_target"
        ]
    ) != 3:
        issues.append("open-obligation count drifted")

    independent_shear_flux(issues)
    independent_orientation_budget(issues)
    independent_heat_guard(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "mathsf_A=mathcal_C_N+alpha*mathsf_X",
        "Gamma_(s,ell)",
        "determinant 1",
        "Psi_ell=mathsf_X+i*mathcal_C_N/ell",
        "(mathcal_C_N)_x=(mathsf_A)_x-alpha_x*mathsf_X",
        "H_0,H_1,H_2,D_0,D_1",
        "partial_s arg(Psi_ell)",
        "less than 2.03e-14",
        "sign(partial_x mathsf_X)=sign(mathcal_C_N)",
        "At W_0=0",
        "exp(m^2*t)*sin(m*x+pi/4)",
        "wind[X_m+i*mathcal_C_m]=-m",
        "do not bound the one-sided crossing count",
        "Delta_(partial D_j)arg(Psi_j)<2*pi",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "pointwise Xi Abel-scalar gap",
        "independent open theorems",
        "No q<1 closure",
        "complete boundary composition",
        "contact exclusion",
        "Lambda<=0",
        "PF-infinity",
        "RH proof",
        "Clay-prize conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pi Provenance",
        "Exact Shear",
        "Reduced Proxy",
        "Crossing Orientation",
        "Zero Fibre",
        "Backward-Heat Guard",
        "Minimal Linked-Point Route",
        "Boundary",
    ):
        if marker not in note:
            issues.append(f"note marker missing: {marker}")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact", type=Path, default=builder.DEFAULT_OUT)
    args = parser.parse_args()
    issues = validate(args.artifact)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}")
        return 1
    print(
        "validated Newman Abel-scalar shear-flux reduction: "
        "18 rows, 1 exact shear/scale homotopy, "
        "1 reduced O(N) physical flux, "
        "1 certified crossing-orientation handoff, "
        "1 backward-heat many-crossing guard, 3 open obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
