#!/usr/bin/env python3
"""Validate the endpoint-complete contact signed-transport reduction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_contact_signed_transport_reduction as builder


EXPECTED_IDS = [f"cstr_{index:02d}_{suffix}" for index, suffix in enumerate(
    (
        "pi_provenance",
        "domain",
        "coordinates",
        "centered_components",
        "endpoint_atom",
        "signed_partition",
        "band_transport",
        "contact_covariance",
        "ordered_transport",
        "conditional_margin",
        "gram_rank",
        "endpoint_null_guard",
        "first_jet",
        "cutoff",
        "zero_projection",
        "legacy_symmetry",
        "source_inequality",
        "route_decision",
        "q_lt_1",
        "boundary",
    )
)]


def independent_centering(issues: list[str]) -> None:
    c, b, u, u_n = sp.symbols("c b u u_N", real=True)
    e_r, g_r, x, y = sp.symbols("e g x y", real=True)
    alpha = c * u_n
    value = e_r + x
    slope = g_r + u * (c * x - b * y)
    centered = (
        g_r - alpha * e_r
        + c * (u - u_n) * x
        - b * u * y
    )
    if sp.simplify(centered - (slope - alpha * value)) != 0:
        issues.append("independent centered-component audit failed")


def independent_transport(issues: list[str]) -> None:
    mp, mm, hp, hm, d0 = sp.symbols(
        "Mp Mm hp hm d0", real=True
    )
    x_value = mp - mm
    mass = (mp + mm) / 2
    lhs = mp * hp - mm * hm + d0
    rhs = mass * (hp - hm) + x_value * (hp + hm) / 2 + d0
    if sp.simplify(lhs - rhs) != 0:
        issues.append("independent mass-transport audit failed")

    h = sp.symbols("h1:6", real=True)
    weights = sp.symbols("w1:6", real=True)
    total = sum(weights)
    prefixes = [sum(weights[:index]) for index in range(1, 5)]
    abel = h[-1] * total - sum(
        (h[index + 1] - h[index]) * prefixes[index]
        for index in range(4)
    )
    direct = sum(h_j * w_j for h_j, w_j in zip(h, weights))
    if sp.simplify(direct - abel) != 0:
        issues.append("independent slope-order Abel audit failed")


def independent_endpoint_guard(issues: list[str]) -> None:
    kappa, h_a, j_re, j_im, time_0, alpha = sp.symbols(
        "kappa H_a J_re J_im T_0 alpha", real=True
    )
    j_a = j_re + sp.I * j_im
    endpoint_value = kappa * h_a * (time_0 + sp.I)
    endpoint_slope = kappa * j_a * (time_0 + sp.I)
    c_0 = sp.re(endpoint_value)
    d_0 = sp.re(endpoint_slope) - alpha * c_0
    expected_d = kappa * (
        time_0 * j_re - j_im - alpha * time_0 * h_a
    )
    effective_slope = (
        j_re / h_a - j_im / (time_0 * h_a) - alpha
    )
    if sp.simplify(c_0 - kappa * time_0 * h_a) != 0:
        issues.append("independent complex-endpoint value failed")
    if sp.simplify(d_0 - expected_d) != 0:
        issues.append("independent complex-endpoint slope failed")
    if sp.simplify(d_0 - c_0 * effective_slope) != 0:
        issues.append("independent complex-endpoint quotient failed")

    delta = sp.Rational(1, 100)
    root = sp.sqrt(65)
    omega = (-8 + sp.I) / root
    s_prime = sp.Rational(1, 16) - sp.I / 2
    endpoint = delta * (8 + sp.I)
    carriers = [
        omega,
        -sp.Rational(3, 5) * omega,
        -sp.Rational(2, 5) * omega,
        -endpoint,
    ]
    distances = [4, 3, 2, 1]
    value = sp.simplify(endpoint + sum(carriers))
    slope = sp.simplify(
        endpoint / 8
        + s_prime * sum(
            distance * carrier
            for distance, carrier in zip(distances, carriers)
        )
    )
    expected = sp.I * (
        sp.Rational(7, 80) * root + sp.Rational(13, 320)
    )
    if value != 0:
        issues.append("independent endpoint-null value failed")
    if sp.simplify(slope - expected) != 0:
        issues.append("independent endpoint-null slope failed")
    amplitudes = [float(sp.N(sp.Abs(value), 30)) for value in carriers]
    if not all(
        amplitudes[index] > amplitudes[index + 1]
        for index in range(3)
    ):
        issues.append("independent endpoint-null amplitude order failed")


def independent_gram_and_cutoff(issues: list[str]) -> None:
    p = sp.Matrix([1, 1, 1])
    q = sp.Matrix([0, 1, 2])
    null = sp.Matrix([1, -2, 1])
    if (p.dot(null), q.dot(null)) != (0, 0):
        issues.append("independent rank-three null vector failed")
    gram = p * p.T + q * q.T
    if gram.rank() != 2:
        issues.append("independent Gram rank failed")

    old, new, x_value, d_x, d_a, a_value = sp.symbols(
        "old new X dX dA A", real=True
    )
    c_old = a_value - old * x_value
    c_new = a_value + d_a - new * (x_value + d_x)
    expected = d_a - new * d_x + (old - new) * x_value
    if sp.simplify(c_new - c_old - expected) != 0:
        issues.append("independent cutoff jump failed")


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
        issues.append("countermodel count drifted")
    if len(
        [row for row in rows if row.get("role") == "open_theorem_target"]
    ) != 2:
        issues.append("open theorem target count drifted")
    if len(
        [row for row in rows if row.get("role") == "route_decision"]
    ) != 2:
        issues.append("route decision count drifted")

    independent_centering(issues)
    independent_transport(issues)
    independent_endpoint_guard(issues)
    independent_gram_and_cutoff(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "sum_(j=0)^N d_j=mathcal_C_N",
        "Re(J_a)/H_a-Im(J_a)/(T_0*H_a)-c*u_N",
        "d_0=kappa*[T_0*Re(J_a)-Im(J_a)]",
        "M*(h_+-h_-)",
        "h_(k+1)-h_k",
        "rank at most one",
        "W_0=0",
        "H_0,H_1,H_2,D_0,D_1",
        "Delta Z_0=q_(N+1)+j_0",
        "q=2*t*L^2>=1",
        "cumulative-real-mass",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    boundary = artifact.get("proof_boundary", "")
    for marker in (
        "not an actual Xi counterexample",
        "No Xi cumulative-mass estimate",
        "Abel-scalar gap",
        "q<1 closure",
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
        "Centered Components",
        "Contact Transport",
        "Slope-Ordered Form",
        "Conditional Margin",
        "Gram Rank Guard",
        "Endpoint-Shaped Null Guard",
        "First Jet",
        "Cutoff Law",
        "Legacy Symmetry Verdict",
        "Open Source Inequality",
        "Route Decision",
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
        "validated Newman contact signed-transport reduction: "
        "20 rows, exact centered endpoint/carrier transport, "
        "slope-order Abel identity, conditional margin, "
        "Gram-rank guard, endpoint-shaped null guard, "
        "five-current first jet, and adjacent-cutoff law"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
