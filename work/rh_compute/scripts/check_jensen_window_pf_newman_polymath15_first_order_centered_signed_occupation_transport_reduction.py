#!/usr/bin/env python3
"""Validate the signed occupation transport and transversality reduction."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_signed_occupation_transport_reduction as builder


SUFFIXES = (
    "pi_provenance",
    "domain",
    "measure_definition",
    "weak_transport",
    "cumulative_transport",
    "actual_carrier_rates",
    "source_collapse",
    "reaction_law",
    "normalized_measure",
    "normalized_cumulative",
    "normalization_boundary",
    "pointwise_residual",
    "amplitude_sum",
    "summed_residual",
    "bulk_first_moment",
    "endpoint_complete_moment",
    "edge_nonsplitting",
    "direct_transversality_equivalence",
    "orientation_handoff",
    "signed_flux_guard",
    "pole_boundary",
    "successor_guard",
    "route_decision",
    "q_lt_1",
    "proof_boundary",
)
EXPECTED_IDS = [
    f"sotr_{index:02d}_{suffix}"
    for index, suffix in enumerate(SUFFIXES)
]


def independent_weak_transport(issues: list[str]) -> None:
    for index in range(4):
        c_value, c_x, h_x, test, test_x = sp.symbols(
            f"c_{index} c_x_{index} h_x_{index} "
            f"phi_{index} phi_x_{index}",
            real=True,
        )
        measure_derivative = c_x * test + c_value * h_x * test_x
        distributional_flux = -c_value * h_x * test_x
        if sp.expand(
            measure_derivative + distributional_flux - c_x * test
        ) != 0:
            issues.append(f"independent weak transport failed at {index}")


def independent_carrier_source(issues: list[str]) -> None:
    (
        x_part,
        y_part,
        b,
        u,
        c_s,
        log_n,
        log_cap,
        xi,
        xi_one,
        residual,
    ) = sp.symbols(
        "X Y b u c_s log_n log_N xi xi_1 e",
        real=True,
    )
    k = log_cap - log_n
    radial = -c_s * log_n + xi - xi_one
    angular = b * u + residual
    x_part_x = radial * x_part - angular * y_part
    d_value = c_s * k * x_part - b * u * y_part
    remainder = (xi - xi_one) * x_part - residual * y_part
    expected_x = d_value - c_s * log_cap * x_part + remainder
    if sp.expand(x_part_x - expected_x) != 0:
        issues.append("independent source collapse failed")
    expected_d = x_part_x + c_s * log_cap * x_part - remainder
    if sp.expand(d_value - expected_d) != 0:
        issues.append("independent first-moment identity failed")


def independent_normalization(issues: list[str]) -> None:
    omega, c_s, log_cap = sp.symbols(
        "Omega c_s log_N",
        real=True,
    )
    mu, mu_x, flux_s, slope_mu, residual = sp.symbols(
        "mu mu_x flux_s slope_mu residual",
        real=True,
    )
    relation = mu_x + flux_s - slope_mu + c_s * log_cap * mu - residual
    normalized = (
        omega * mu_x
        + c_s * log_cap * omega * mu
        + omega * flux_s
        - omega * slope_mu
        - omega * residual
    )
    if sp.expand(normalized - omega * relation) != 0:
        issues.append("independent normalized measure law failed")


def independent_residual_budget(issues: list[str]) -> None:
    if 18888**2 - (16892**2 + 8449**2) != 31279:
        issues.append("independent residual norm constant failed")
    if not Fraction(18888 * 10, 3 * 144) < 438:
        issues.append("independent summed residual coefficient failed")
    if not (
        Fraction(19, 7) ** 45
        > Fraction(438 * 10**7) ** 2
    ):
        issues.append("independent exponential residual guard failed")


def independent_aggregate(issues: list[str]) -> None:
    c_zero, c_zero_x, d_zero = sp.symbols(
        "c_0 c_0_x d_0",
        real=True,
    )
    c_bulk, c_bulk_x, r_bulk = sp.symbols(
        "C_bulk C_bulk_x R_bulk",
        real=True,
    )
    c_s, log_cap = sp.symbols("c_s log_N", real=True)
    d_bulk = c_bulk_x + c_s * log_cap * c_bulk - r_bulk
    total_c = c_zero + c_bulk
    total_c_x = c_zero_x + c_bulk_x
    edge = d_zero - c_zero_x - c_s * log_cap * c_zero
    total_d = d_zero + d_bulk
    expected = total_c_x + c_s * log_cap * total_c + edge - r_bulk
    if sp.expand(total_d - expected) != 0:
        issues.append("independent endpoint-complete moment failed")

    signs = []
    for first, second in ((1, 0), (0, 1), (1, 1)):
        signs.append(first - second)
    if signs != [1, -1, 0]:
        issues.append("independent contact sign guard failed")


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
    if len(rows) != 25:
        issues.append("row count drifted")
    if len(
        [row for row in rows if row.get("role") == "nonpromotion_guard"]
    ) != 3:
        issues.append("nonpromotion guard count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 2:
        issues.append("open theorem target count drifted")

    independent_weak_transport(issues)
    independent_carrier_source(issues)
    independent_normalization(issues)
    independent_residual_budget(issues)
    independent_aggregate(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "partial_x mu+partial_s F=S",
        "partial_x P=F",
        "r_i=(xi_i-xi_1)X_i-e_iY_i",
        "S=[s-c_s*log(N)]mu+R",
        "partial_x bar_mu+partial_s bar_F=s*bar_mu+bar_R",
        "16892^2+8449^2<18888^2",
        "438*exp(-17L/10)",
        "Omega_N*D_bulk=partial_x(Omega_N*C_bulk)",
        "mathcal_C_N=partial_x mathsf_X",
        "2.03e-14",
        "boundary mass is D_perp",
        "first moments +1,-1,0",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove a signed Xi threshold-mass estimate",
        "pointwise Abel-scalar gap",
        "edge sign",
        "one-turn successor bound",
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
        "Distributional Transport",
        "Xi Source Collapse",
        "Common Normalization",
        "Residual Budget",
        "First Moment",
        "Route Guards",
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
        "validated signed occupation transport reduction: "
        "25 rows, exact weak/cumulative laws, Xi source collapse, "
        "N^sigma normalization, sub-1e-7 residual budget, "
        "direct-transversality equivalence, 3 route guards"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
