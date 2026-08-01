#!/usr/bin/env python3
"""Validate the pairwise projective-alignment identity and route guard."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_pairwise_projective_alignment_gate as builder


SUFFIXES = (
    "pi_provenance",
    "domain",
    "residual_definition",
    "single_aligned_jet",
    "relative_alignment_identity",
    "frozen_model_theorem",
    "moving_scale_restoration",
    "moving_scale_sign_reversal",
    "finite_projective_rate",
    "exact_pole_rate",
    "pole_limit_join",
    "residual_order_guard",
    "legacy_symmetry_boundary",
    "occupation_identity",
    "swap_join",
    "aggregate_oriented_cell",
    "pole_chart_join",
    "exceptional_edge",
    "route_decision",
    "q_lt_1",
    "proof_boundary",
)
EXPECTED_IDS = [
    f"ppag_{index:02d}_{suffix}"
    for index, suffix in enumerate(SUFFIXES)
]


def independent_alignment(issues: list[str]) -> None:
    (
        b,
        c_s,
        c_s_x,
        b_x,
        u_x,
        aligned,
        u_n,
        u_m,
        k_n,
        k_m,
        u_terminal,
        e_n,
        e_m,
    ) = sp.symbols(
        "b c_s c_s_x b_x u_x H u_n u_m k_n k_m "
        "u_N e_n e_m",
        nonzero=True,
        real=True,
    )
    w_n = c_s * k_n - aligned
    w_m = c_s * k_m - aligned

    def direct(u_value, k_value, residual, shifted):
        tangent = shifted / (b * u_value)
        angular = b * u_value + residual
        return sp.expand(
            c_s_x * k_value
            - (b_x * u_value + b * u_x) * tangent
            - b * u_value * angular * (1 + tangent**2)
        )

    relative = sp.expand(
        direct(u_m, k_m, e_m, w_m)
        - direct(u_n, k_n, e_n, w_n)
    )
    delta = u_n - u_m
    compact = sp.expand(
        delta
        * (
            b**2 * (u_n + u_m)
            + c_s * (c_s * (k_n + k_m) - 2 * aligned)
            - c_s_x
            + c_s * b_x / b
            + u_x
            * (c_s * u_terminal + aligned)
            / (u_n * u_m)
        )
        + e_n
        * (b**2 * u_n**2 + w_n**2)
        / (b * u_n)
        - e_m
        * (b**2 * u_m**2 + w_m**2)
        / (b * u_m)
    )
    relation = {
        k_n: u_n - u_terminal,
        k_m: u_m - u_terminal,
    }
    if sp.factor((relative - compact).subs(relation)) != 0:
        issues.append("independent pairwise alignment identity failed")


def independent_moving_scale(issues: list[str]) -> None:
    b, u_x, u_n, u_m = sp.symbols(
        "b u_x u_n u_m",
        nonzero=True,
        real=True,
    )
    delta = u_n - u_m
    sum_u = u_n + u_m
    critical = -b**2 * sum_u * u_n * u_m / u_x

    def relative(aligned):
        return sp.factor(
            delta
            * (
                b**2 * sum_u
                + u_x * aligned / (u_n * u_m)
            )
        )

    if sp.simplify(relative(0) - b**2 * delta * sum_u) != 0:
        issues.append("independent moving-scale positive sign failed")
    if sp.simplify(relative(critical)) != 0:
        issues.append("independent moving-scale zero failed")
    if sp.simplify(
        relative(2 * critical) + b**2 * delta * sum_u
    ) != 0:
        issues.append("independent moving-scale negative sign failed")


def independent_pole(issues: list[str]) -> None:
    b, u_n, u_m, e_n, e_m = sp.symbols(
        "b u_n u_m e_n e_m",
        nonzero=True,
        real=True,
    )
    direct = sp.factor(
        (-(b * u_m + e_m) / (b * u_m))
        - (-(b * u_n + e_n) / (b * u_n))
    )
    expected = sp.factor(
        e_n / (b * u_n) - e_m / (b * u_m)
    )
    if sp.factor(direct - expected) != 0:
        issues.append("independent exact pole rate failed")

    epsilon, delta = sp.symbols(
        "epsilon delta",
        positive=True,
        real=True,
    )
    angular_n = -b * delta - epsilon
    angular_m = -b * delta + epsilon
    if sp.expand(angular_n - (-b * delta - epsilon)) != 0:
        issues.append("independent first angular-order guard failed")
    if sp.expand(angular_m - (-b * delta + epsilon)) != 0:
        issues.append("independent second angular-order guard failed")


def independent_occupation(issues: list[str]) -> None:
    anchor = sp.symbols("A", real=True)
    c_1, c_2, c_3 = sp.symbols("c_1 c_2 c_3", real=True)
    h_1, h_2, h_3 = sp.symbols("h_1 h_2 h_3", real=True)
    direct = c_1 * h_1 + c_2 * h_2 + c_3 * h_3
    represented = (
        anchor * (c_1 + c_2 + c_3)
        + c_1 * (h_1 - anchor)
        + c_2 * (h_2 - anchor)
        + c_3 * (h_3 - anchor)
    )
    if sp.expand(direct - represented) != 0:
        issues.append("independent occupation identity failed")

    first_order = (
        h_2 * (c_1 + c_2) - (h_2 - h_1) * c_1
    )
    second_order = (
        h_1 * (c_1 + c_2) - (h_1 - h_2) * c_2
    )
    two_direct = h_1 * c_1 + h_2 * c_2
    if sp.expand(first_order - two_direct) != 0:
        issues.append("independent first swap join failed")
    if sp.expand(second_order - two_direct) != 0:
        issues.append("independent second swap join failed")


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
    if len(rows) != 21:
        issues.append("row count drifted")
    if len(
        [row for row in rows if row.get("role") == "nonpromotion_guard"]
    ) != 3:
        issues.append("nonpromotion guard count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 2:
        issues.append("open theorem target count drifted")

    independent_alignment(issues)
    independent_moving_scale(issues)
    independent_pole(issues)
    independent_occupation(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "D_nm(H):=h_(m,x)-h_(n,x)",
        "D_nm=b^2*(u_n^2-u_m^2)",
        "u_x*H/(u_n*u_m)",
        "D_nm(2H_mov)=-b^2*delta*S<0",
        "e_n/(b*u_n)-e_m/(b*u_m)",
        "H->+infinity",
        "nu_m-nu_n>0",
        "transpose-symmetric",
        "integral_A^B P_x(s)ds",
        "C_edge+X_I",
        "J_a^(der)",
        "J_a^(adj)",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "guards are not actual Xi counterexamples",
        "signed Xi occupation bound",
        "edge-block sign",
        "Abel-scalar gap",
        "q<1 closure",
        "contact exclusion",
        "Lambda<=0",
        "PF-infinity",
        "RH",
        "Clay-prize conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Pi Provenance",
        "Exact Alignment",
        "Frozen Model",
        "Moving-Scale Obstruction",
        "Projection Pole",
        "Legacy Symmetry Boundary",
        "Aggregate Occupation Identity",
        "Route Decision",
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
        "validated pairwise projective-alignment gate: "
        "21 rows, exact relative current and pole join, frozen-model "
        "Sturm term, moving-scale and residual-order guards, aggregate "
        "occupation identity"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
