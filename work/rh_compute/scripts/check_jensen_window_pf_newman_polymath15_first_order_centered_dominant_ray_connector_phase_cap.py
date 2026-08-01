#!/usr/bin/env python3
"""Validate the dominant-ray connector phase cap and phase ledger."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_dominant_ray_connector_phase_cap as builder


EXPECTED_IDS = [
    "drcpc_00_pi_provenance",
    "drcpc_01_connector_geometry",
    "drcpc_02_normalized_vector",
    "drcpc_03_normalizer_phase",
    "drcpc_04_leading_ellipse",
    "drcpc_05_uniform_error",
    "drcpc_06_uniform_cone",
    "drcpc_07_leading_phase",
    "drcpc_08_connector_cap",
    "drcpc_09_axis_phase",
    "drcpc_10_composite_ledger",
    "drcpc_11_horizontal_target",
    "drcpc_12_seam_boundary",
    "drcpc_13_endpoint_guard",
    "drcpc_14_pointwise_target",
    "drcpc_15_horizontal_budget",
    "drcpc_16_small_q_and_finite",
]


def independent_rational_audit(issues: list[str]) -> None:
    value_error = Fraction(32, 125) + Fraction(1, 8000)
    slope_error = Fraction(6901, 100000) + Fraction(7, 160000)
    if value_error != Fraction(2049, 8000):
        issues.append("independent value-error fraction failed")
    if slope_error != Fraction(55243, 800000):
        issues.append("independent slope-error fraction failed")
    if not value_error**2 + slope_error**2 < Fraction(4, 15) ** 2:
        issues.append("independent Euclidean error cap failed")
    if Fraction(4, 15) / Fraction(12, 25) != Fraction(5, 9):
        issues.append("independent relative cone cap failed")
    sine_lower = Fraction(3, 5) - Fraction(3, 5) ** 3 / 6
    if sine_lower != Fraction(141, 250) or not sine_lower > Fraction(5, 9):
        issues.append("independent sine cone bound failed")
    lead = Fraction(25, 384) * Fraction(22, 7) + Fraction(1, 24)
    if lead != Fraction(331, 1344):
        issues.append("independent leading phase cap failed")
    total = Fraction(6, 5) + lead
    if total != Fraction(9719, 6720) or not total < Fraction(3, 2):
        issues.append("independent connector phase cap failed")


def independent_ellipse_current(issues: list[str]) -> None:
    beta, lam = sp.symbols("beta lam", real=True, positive=True)
    x_value = sp.cos(beta)
    y_value = lam * sp.sin(beta)
    beta_numerator = sp.simplify(
        x_value * sp.diff(y_value, beta)
        - y_value * sp.diff(x_value, beta)
    )
    lambda_numerator = sp.simplify(
        x_value * sp.diff(y_value, lam)
        - y_value * sp.diff(x_value, lam)
    )
    if sp.simplify(beta_numerator - lam) != 0:
        issues.append("ellipse beta-current identity failed")
    if sp.simplify(
        lambda_numerator - sp.sin(beta) * sp.cos(beta)
    ) != 0:
        issues.append("ellipse lambda-current identity failed")


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
    if len(rows) != 17:
        issues.append("row count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 3:
        issues.append("open-obligation count drifted")
    if len([row for row in rows if row.get("role") == "countermodel"]) != 1:
        issues.append("countermodel count drifted")
    if len(
        [row for row in rows if row.get("role") == "proved_boundary_lemma"]
    ) != 2:
        issues.append("proved-boundary-lemma count drifted")

    independent_rational_audit(issues)
    independent_ellipse_current(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "I_L=[25/L,25/(L-1)]",
        "partial_t beta=Im(alpha^2)/4",
        "V_0(t)=(2*cos(beta),2*lambda*sin(beta))",
        "|V-V_0|/|V_0|<5/9",
        "|delta(t)|<3/5",
        "331/1344",
        "9719/6720<3/2<pi/2",
        "|Delta_(right D_j)arg V|<pi/2",
        "2*pi*kappa_j=H_j+C_j",
        "H_j<3*pi/2",
        "G_m(u)=exp(2*pi*i*m*u)",
        "may not be discarded or charged twice",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the Xi Abel-scalar gap",
        "3*pi/2 horizontal-composite inequality",
        "q<1 closure",
        "finite shoulders or endpoint tracks",
        "complete boundary composition",
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
        "Ray Geometry",
        "One-Saddle Ellipse",
        "Uniform Cone",
        "Connector Cap",
        "Phase Ledger",
        "Seam Guard",
        "Endpoint-Only Guard",
        "Remaining Target",
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
        "validated Newman dominant-ray connector phase cap: "
        "17 rows, 1 uniform cone, 1 strict quarter-turn connector cap, "
        "1 reduced 3*pi/2 horizontal ledger, "
        "1 endpoint-only winding guard, 3 open obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
