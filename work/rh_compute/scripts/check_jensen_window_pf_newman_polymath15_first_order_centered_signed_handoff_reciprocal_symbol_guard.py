#!/usr/bin/env python3
"""Validate the signed zeta-handoff reciprocal-symbol no-gain guard."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_signed_handoff_reciprocal_symbol_guard as builder


EXPECTED_IDS = [
    "shrs_00_pi_provenance",
    "shrs_01_signed_handoff",
    "shrs_02_stationary_term",
    "shrs_03_interior_symbol",
    "shrs_04_normalized_symbol",
    "shrs_05_active_scale",
    "shrs_06_c2_audit",
    "shrs_07_cstar_audit",
    "shrs_08_kernel_derivative",
    "shrs_09_full_term_derivative",
    "shrs_10_cutoff_guard",
    "shrs_11_endpoint_guard",
    "shrs_12_route_decision",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    t, a_log, z, sigma = sp.symbols("t A z sigma", real=True)
    u_log = a_log - z
    nu_log = a_log + z
    q_u = t * u_log**2 / 4
    q_nu = t * nu_log**2 / 4
    delta = 2 * sigma - 1 - t * a_log
    if sp.simplify(
        q_u + (2 * sigma - 1) * z - q_nu - delta * z
    ) != 0:
        issues.append("independent weighted symbol identity failed")

    signed = sp.exp((2 * sigma - 1) * z) * (sp.exp(q_u) - 1)
    weighted = sp.exp(q_nu + delta * z)
    if sp.simplify(signed / weighted - (1 - sp.exp(-q_u))) != 0:
        issues.append("independent normalized signed symbol failed")

    r = Fraction(125_662, 155_153)
    dual = Fraction(184_644, 155_153)
    weighted_exponent = dual**2 / 4
    ordinary_exponent = 1 - r
    gap = r**2 / 4
    if weighted_exponent - ordinary_exponent != gap:
        issues.append("independent c=2 exponent split failed")
    if gap != Fraction(3_947_734_561, 24_072_453_409):
        issues.append("independent active relative gap drifted")

    c_star = Fraction(4_911_678_521, 1_933_561_194)
    if c_star * r**2 / 8 != Fraction(
        1_989_040_967, 9_549_356_844
    ):
        issues.append("independent c_* relative gap drifted")

    d_two = Fraction(3_133_668_399, 48_144_906_818)
    if d_two + Fraction(3, 4) != Fraction(
        78_484_697_025, 96_289_813_636
    ):
        issues.append("independent endpoint separation drifted")


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

    urls = artifact.get("source_urls", {})
    if urls.get("polymath15_primary") != builder.POLYMATH_SOURCE_URL:
        issues.append("Polymath-15 source URL drifted")
    if urls.get("exponent_pairs_primary") != builder.EXPONENT_PAIR_SOURCE_URL:
        issues.append("exponent-pair source URL drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    if len(rows) != 13:
        issues.append("row count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "exact_first_jet_identity"
        ]
    ) != 2:
        issues.append("first-jet identity count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "exact_rational_certificate"
        ]
    ) != 2:
        issues.append("rational-audit count drifted")
    if len(
        [
            row
            for row in rows
            if row.get("role") == "nonpromotion_guard"
        ]
    ) != 2:
        issues.append("endpoint/cutoff guard count drifted")

    independent_symbolic_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "F_N(u)=w_t(u)*1_(u<=N)-1",
        "K_N(nu)=y^(2sigma-1)F_N(u_nu)",
        "Sigma_N(nu)=1-exp[-t*log(u_nu)^2/4]",
        "Sigma_N=1-N^(-c*r_*^2/8+o(1))",
        "3947734561/24072453409",
        "1989040967/9549356844",
        "partial_x log(Sigma_N)=tU*A_x/(exp(q)-1)",
        "29491/155153",
        "78484697025/96289813636",
        "Retire signed primal/tail pairwise cancellation",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "continuous signed reciprocal symbol",
        "fixed-chart first-x derivative",
        "does not prove a discrete endpoint-complete",
        "does not",
        "dual sum",
        "Abel-scalar gap",
        "contact exclusion",
        "RH proof",
        "prize-level conclusion",
    ):
        if marker not in boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Exact Signed Handoff",
        "Pi Provenance",
        "Reciprocal Symbol",
        "Active Radius",
        "First X Derivative",
        "Endpoint And Cutoff",
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
        "validated Newman signed-handoff reciprocal-symbol guard: "
        "13 rows, 3 exact symbol identities, 2 first-jet identities, "
        "2 rational exponent audits, 2 endpoint/cutoff guards, "
        "1 retired pairwise route"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
