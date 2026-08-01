#!/usr/bin/env python3
"""Validate the centered first-order adjacent-saddle recurrence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp
import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_adjacent_saddle_recurrence as builder


EXPECTED_IDS = [
    "nfocasr_01_recurrence_derivatives",
    "nfocasr_02_endpoint_sum",
    "nfocasr_03_endpoint_derivative",
    "nfocasr_04_scalar_jump",
    "nfocasr_05_mismatch_identity",
    "nfocasr_06_real_projection",
    "nfocasr_07_first_coefficient",
    "nfocasr_08_formal_scale",
    "nfocasr_09_selected_audit",
    "nfocasr_10_route_audit",
    "nfocasr_11_uniform_target",
    "nfocasr_12_bulk_handoff",
    "nfocasr_13_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    p, a, r = sp.symbols("p a r", positive=True, real=True)
    y = p + 1
    phase = sp.exp(
        sp.pi
        * sp.I
        * (p**2 / 2 + p + sp.Rational(3, 8))
    )
    third = sp.diff(phase, p, 3)
    fourth = sp.diff(phase, p, 4)
    if sp.simplify(
        third / phase
        + 3 * sp.pi**2 * y
        + sp.I * sp.pi**3 * y**3
    ) != 0:
        issues.append("independent R''' identity failed")
    if sp.simplify(
        fourth / phase
        - (
            sp.pi**4 * y**4
            - 6 * sp.I * sp.pi**3 * y**2
            - 3 * sp.pi**2
        )
    ) != 0:
        issues.append("independent R'''' identity failed")

    endpoint = phase + third / (12 * sp.pi**2 * a)
    expected = phase * (
        1 - r / 2 - 2 * sp.I * sp.pi * a**2 * r**3 / 3
    )
    if sp.simplify(
        endpoint.subs(p, 2 * a * r - 1)
        - expected.subs(p, 2 * a * r - 1)
    ) != 0:
        issues.append("independent endpoint recurrence failed")

    time, w = sp.symbols("time w", real=True)
    shift = sp.Rational(1, 2) - sp.I * sp.pi * time / 8
    log_ratio_first = (
        -shift * w
        - sp.I * sp.pi * time * w / 8
        - 2 * sp.I * sp.pi * w**3 / 3
    )
    endpoint_first = -w / 2 - 2 * sp.I * sp.pi * w**3 / 3
    if sp.simplify(log_ratio_first - endpoint_first) != 0:
        issues.append("independent first-coefficient cancellation failed")

    u, v = sp.symbols("u v", real=True)
    dx, dy, du = sp.symbols("dX dY dU", real=True)
    correction_real = sp.symbols("D", real=True)
    projected = du - u * dx + v * dy - correction_real
    complex_form = sp.re(
        (du + sp.I * sp.symbols("dV", real=True))
        - (u + sp.I * v) * (dx + sp.I * dy)
    ) - correction_real
    if sp.simplify(projected - complex_form) != 0:
        issues.append("independent real projection failed")


def validate(path: Path) -> list[str]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    issues: list[str] = []
    if artifact.get("kind") != builder.STEM:
        issues.append("artifact kind mismatch")
    if artifact.get("source_sha256") != builder.source_hashes():
        issues.append("source hash chain drifted")
    if artifact.get("source_audit") != builder.source_audit():
        issues.append("source audit drifted")
    if artifact.get("exact") != builder.symbolic_audit():
        issues.append("exact payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row
        for row in rows
        if row.get("readiness") == "not_ready_to_apply"
    ]
    if len(open_rows) != 2:
        issues.append(f"open obligation count drifted: {len(open_rows)}")

    diagnostics = artifact.get("diagnostics", {})
    if diagnostics.get("row_count") != 18:
        issues.append("diagnostic row count drifted")
    if mp.mpf(diagnostics.get("max_relative_delta", "1")) >= mp.mpf(
        "1e-45"
    ):
        issues.append("cross-precision drift exceeded 1e-45")
    if mp.mpf(
        diagnostics.get("max_leading_scaled_e3", "1")
    ) >= mp.mpf("0.24"):
        issues.append("selected edge scale exceeded 0.24")
    if mp.mpf(
        diagnostics.get("max_scalar_jump_scaled_e7", "2")
    ) >= 2:
        issues.append("selected scalar scale exceeded 2")
    diag_rows = diagnostics.get("rows", [])
    if len(diag_rows) != 18:
        issues.append("stored diagnostic rows drifted")
    for row in diag_rows:
        if mp.mpf(row.get("L", "0")) < 50:
            issues.append("diagnostic below L=50")
            break
        finite = mp.mpf(row.get("finite_real_scaled_e3", "0"))
        endpoint = mp.mpf(row.get("endpoint_real_scaled_e3", "0"))
        if finite * endpoint >= 0:
            issues.append("selected leading pieces are not opposite-signed")
            break

    independent_symbolic_audit(issues)
    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "R'''=",
        "R''''=",
        "J_a=H_a(p+2)+H_a(p)",
        "J_(a,x)/R=",
        "Delta A_a=Re[-s_*'*ell*f_(N+1)",
        "Q_N=f_(N+1)+kappa_N*J_a",
        "Delta U-u_a*Delta X+v_a*Delta Y",
        "[epsilon]log(main_sharp/endpoint_C1)=0",
        "O(exp(-7L/4))",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "do not prove the uniform differentiated ratio bound",
        "q>=1 scalar lower bound",
        "q<1",
        "finite phase cells",
        "one-sided successor winding",
        "Lambda<=0",
        "RH",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Exact Adjacent Recurrence",
        "Centered Scalar Jump",
        "Leading Cancellation",
        "Selected High-Height Audit",
        "Route Decision",
        "complex half itself is chart dependent",
        "not an",
        "RH-level object",
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
        "validated Newman first-order centered adjacent-saddle recurrence: "
        "13 rows, exact Delta A identity, vanished a^-1 coefficient, "
        "18 selected L>=50 rows with e^-7L/4 scalar scaling, "
        "2 open scalar obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
