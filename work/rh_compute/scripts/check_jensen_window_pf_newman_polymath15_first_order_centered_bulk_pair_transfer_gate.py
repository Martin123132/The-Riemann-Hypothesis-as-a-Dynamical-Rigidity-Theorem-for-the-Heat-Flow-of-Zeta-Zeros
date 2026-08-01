#!/usr/bin/env python3
"""Validate the centered bulk adjacent-pair transfer gate."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import sympy as sp

import jensen_window_pf_newman_polymath15_first_order_centered_bulk_pair_transfer_gate as builder


EXPECTED_IDS = [
    "nfocbptg_01_reflected_basis",
    "nfocbptg_02_abel_basis",
    "nfocbptg_03_component_transfer",
    "nfocbptg_04_ratio_lock",
    "nfocbptg_05_amplitude_decay",
    "nfocbptg_06_pair_determinant",
    "nfocbptg_07_positive_pair",
    "nfocbptg_08_pair_separation",
    "nfocbptg_09_edge_composition",
    "nfocbptg_10_component_countermodel",
    "nfocbptg_11_pair_countermodel",
    "nfocbptg_12_xi_handoff",
    "nfocbptg_13_scalar_target",
    "nfocbptg_14_nonpromotion",
]


def independent_symbolic_audit(issues: list[str]) -> None:
    c, b, ell, next_ell = sp.symbols(
        "c b ell next_ell", real=True
    )
    rho, theta = sp.symbols("rho theta", positive=True, real=True)
    transfer = sp.Matrix([[1, 0], [-c * ell, b * ell]])
    next_transfer = sp.Matrix(
        [[1, 0], [-c * next_ell, b * next_ell]]
    )
    rotation = sp.Matrix(
        [
            [sp.cos(theta), -sp.sin(theta)],
            [sp.sin(theta), sp.cos(theta)],
        ]
    )
    if sp.factor(transfer.det() - b * ell) != 0:
        issues.append("independent component determinant failed")

    pair = transfer + rho * next_transfer * rotation
    expected = (
        b * (ell + rho**2 * next_ell)
        + rho
        * (
            b * (ell + next_ell) * sp.cos(theta)
            + c * (next_ell - ell) * sp.sin(theta)
        )
    )
    if sp.trigsimp(pair.det() - expected) != 0:
        issues.append("independent pair determinant failed")

    big_l, big_m = sp.symbols("L_0 M_0", positive=True)
    factor = (
        (big_l + rho**2 * big_m) ** 2
        - rho**2 * (big_l + big_m) ** 2
    )
    expected_factor = (
        (1 - rho**2)
        * (big_l - rho * big_m)
        * (big_l + rho * big_m)
    )
    if sp.factor(factor - expected_factor) != 0:
        issues.append("independent determinant factorization failed")

    alpha_n, alpha_prime, time, h = sp.symbols(
        "alpha_n alpha_prime time h"
    )
    difference = (
        alpha_prime
        * time**2
        * ((alpha_n - h) ** 2 - alpha_n**2)
        / 8
    )
    expected_difference = (
        alpha_prime * time**2 * (h**2 - 2 * h * alpha_n) / 8
    )
    if sp.expand(difference - expected_difference) != 0:
        issues.append("independent d-difference failed")

    terms = sp.symbols("f_1:6")
    weights = sp.symbols("e_1:6")
    partials = []
    running = 0
    for term in terms:
        running += term
        partials.append(running)
    abel = weights[-1] * partials[-1]
    for index in range(len(terms) - 1):
        abel -= (weights[index + 1] - weights[index]) * partials[index]
    direct = sum(weight * term for weight, term in zip(weights, terms))
    if sp.expand(abel - direct) != 0:
        issues.append("independent five-term Abel identity failed")

    target_vectors = [
        sp.Matrix([1, 0]),
        sp.Matrix([0, 1]),
        sp.Matrix([-1, -1]),
    ]
    ell_values = sp.symbols("q_1:4", nonzero=True, real=True)
    mapped = []
    for target, ell_value in zip(target_vectors, ell_values):
        coefficient = sp.Matrix(
            [
                target[0],
                (target[1] / ell_value + c * target[0]) / b,
            ]
        )
        local_transfer = sp.Matrix(
            [[1, 0], [-c * ell_value, b * ell_value]]
        )
        mapped.append(local_transfer * coefficient)
    if sp.simplify(sum(mapped, sp.zeros(2, 1))) != sp.zeros(2, 1):
        issues.append("independent aggregate countermodel failed")

    p_1, q_1, p_2, q_2 = sp.symbols(
        "p_1 q_1 p_2 q_2", positive=True, real=True
    )
    phase_x, phase_y = sp.symbols("phase_x phase_y", real=True)
    first_block = sp.Matrix([p_1 * phase_x, q_1 * phase_y])
    second_block = sp.Matrix([p_2 * phase_x, q_2 * phase_y])
    turning = sp.det(sp.Matrix.hstack(first_block, second_block))
    if sp.factor(
        turning
        - (p_1 * q_2 - q_1 * p_2) * phase_x * phase_y
    ) != 0:
        issues.append("independent common-phase turning guard failed")


def independent_numeric_audit(issues: list[str]) -> None:
    x_min = 4 * math.pi * math.exp(50)
    if not builder.D_ABSOLUTE_CONSTANT / x_min < 0.5:
        issues.append("independent d-disk bound failed")

    corrected = (
        -33 / 80
        + 1 / (16 * x_min)
        + 2 * builder.D_DIFFERENCE_CONSTANT / x_min
    )
    if not corrected < -2 / 5:
        issues.append("independent amplitude-decay margin failed")

    h_max = math.log(2)
    derivative_at_h_max = (
        (2 / 5) * math.exp(-2 * h_max / 5) - 1 / 4
    )
    if not derivative_at_h_max > 0:
        issues.append("independent one-minus-ratio derivative failed")
    if not 1 - math.exp(-2 * h_max / 5) > h_max / 4:
        issues.append("independent one-minus-ratio endpoint failed")

    time_max = 0.5
    saddle = math.sqrt(x_min / (4 * math.pi) + time_max / 16)
    if not saddle < x_min / 4:
        issues.append("independent saddle/frame comparison failed")

    # The determinant proof uses only these four scalar inequalities.
    h_test = math.log(2)
    rho_test = math.exp(-2 * h_test / 5)
    frame_b = 0.5
    frame_c = h_test / 16
    lower = h_test * (
        frame_b * (1 - rho_test) - frame_c * rho_test
    )
    if not lower > h_test**2 / 16:
        issues.append("independent determinant endpoint margin failed")


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
    if artifact.get("constants") != {
        "L_min": 50,
        "tL_max": 25,
        "d_absolute": 2_189,
        "d_difference": 44,
        "amplitude_decay": "2/5",
        "pair_determinant_denominator": 16,
        "pair_separation_denominator": 48,
    }:
        issues.append("constant payload drifted")

    rows = artifact.get("rows", [])
    if [row.get("id") for row in rows] != EXPECTED_IDS:
        issues.append("row id/order mismatch")
    open_rows = [
        row
        for row in rows
        if row.get("readiness") == "not_ready_to_apply"
    ]
    if len(open_rows) != 2:
        issues.append(f"open target count drifted: {len(open_rows)}")

    independent_symbolic_audit(issues)
    independent_numeric_audit(issues)

    exact_text = json.dumps(artifact.get("exact", {}), sort_keys=True)
    for marker in (
        "2*(X,A_a)^T=",
        "det(T_ell)=b*ell",
        "M_(1,a)=ell_N*F_N",
        "f_(n+1)/f_n=rho_n*exp(i*delta_n)",
        "d_(n+1)-d_n=",
        "log(rho_n)<-(2/5)h_n",
        "det(B_n)=b*(ell_n+rho_n^2*ell_(n+1))",
        ">=h_n^2/16",
        "h_n^2*|f_n|/(48L)",
        "(1,0),(0,1),(-1,-1)",
        "Cross-block Xi ratios are essential",
        "u_k^T*B_(2k-1)^T*J*B_(2k+1)",
        "det(W_1,W_2)=(p_1*q_2-q_1*p_2)*x*y",
        "normalizer/endpoint phase anchor",
    ):
        if marker not in exact_text:
            issues.append(f"exact marker missing: {marker}")

    proof_boundary = artifact.get("proof_boundary", "")
    for marker in (
        "does not prove the q>=1 bulk scalar lower bound",
        "q<1",
        "finite phase cells",
        "one-sided successor winding",
        "Lambda<=0",
        "RH",
        "PF-infinity",
        "Clay-prize",
    ):
        if marker not in proof_boundary:
            issues.append(f"proof-boundary marker missing: {marker}")

    note = builder.DEFAULT_NOTE.read_text(encoding="utf-8")
    for marker in (
        "Reflected Bulk Coordinates",
        "Exact Ratio Lock",
        "Locked-Pair Determinant",
        "Falsification Gate",
        "Xi-Specific Handoff",
        "cancellation problem remains open",
        "not a proof of RH",
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
        "validated Newman centered bulk pair-transfer gate: "
        "14 rows, rho_n<exp(-2h_n/5), det(B_n)>=h_n^2/16, "
        "1 exact aggregate countermodel family, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
