#!/usr/bin/env python3
"""Validate the local Mertens-path cosine transference reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_local_path_cosine_transference.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_local_path_cosine_transference.md"
)

REQUIRED_IDS = {
    f"mlpath_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "tail_definitions"),
        (2, "dyadic_tail_vector"),
        (3, "weighted_local_path"),
        (4, "affine_tail_path"),
        (5, "path_dct"),
        (6, "affine_dct_transference"),
        (7, "summation_by_parts"),
        (8, "block_parseval"),
        (9, "centered_path_variance"),
        (10, "deleted_dct_gram"),
        (11, "deleted_dct_spectrum"),
        (12, "low_mode_path_criterion"),
        (13, "cofinal_rh_criterion"),
        (14, "two_component_split"),
        (15, "local_mertens_path"),
        (16, "abel_forward"),
        (17, "abel_inverse"),
        (18, "path_norm_upper"),
        (19, "path_norm_lower"),
        (20, "projection_scope_guard"),
        (21, "generic_countermodel"),
        (22, "countermodel_envelope"),
        (23, "countermodel_zero_means"),
        (24, "countermodel_first_mode"),
        (25, "countermodel_divergence"),
        (26, "generic_method_barrier"),
        (27, "linear_phase_sources"),
        (28, "three_quarter_guard"),
        (29, "rational_zero_density_guard"),
        (30, "vaughan_handoff_guard"),
        (31, "scope_guard"),
        (32, "open_two_component_gate"),
    )
}

ROLE_STATUS = {
    "exact_definition": "available_exact",
    "exact_identity": "available_exact",
    "exact_inequality": "available_exact",
    "exact_equivalence": "available_exact",
    "exact_consequence": "available_exact",
    "literature_guard": "source_backed",
    "proof_guard": "guard_validated",
    "open_target": "open_target",
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Mertens Local-Path Cosine Transference",
    "T_(K,m):=sum_(ell=1)^m q_(K+ell)",
    "c_(K,0)=sqrt(K)r_K-tau_(K,0)",
    "The factor `1/sin(pi*r/(2K))`",
    "centered path variance",
    "U_K^*U_K=I_(K-1)-K^(-1)11^T",
    "RH-Equivalent Two-Component Gate",
    "a_m:=(K+m)^(-1-alpha)",
    "[2^alpha(3+alpha)]^(-1)",
    "does not commute with the low-mode projector",
    "Zero-Mean First-Mode Countermodel",
    "T*_(K,m):=A_K sin(2pi*m/K)",
    "r*_K=0",
    "2^(-4-2alpha)pi^(-4)K^(1-alpha)",
    "deliberately not a Mobius model",
    "3/4+epsilon",
    "zero-density terms",
    "one fixed cofinal sequence",
    "RH-equivalent",
    "not a proof",
)


def close(
    left: complex | float,
    right: complex | float,
    tolerance: float = 5e-11,
) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def norm(vector: list[float]) -> float:
    return math.sqrt(dot(vector, vector))


def cosine_basis(size: int) -> list[list[float]]:
    basis = [[1.0 / math.sqrt(size)] * size]
    for frequency in range(1, size):
        basis.append(
            [
                math.sqrt(2.0 / size)
                * math.cos(
                    math.pi * frequency * (index + 0.5) / size
                )
                for index in range(size)
            ]
        )
    return basis


def check_schema(payload: dict, note: str, issues: list[str]) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_local_path_cosine_transference"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    if len(rows) != 32:
        issues.append(f"expected 32 rows, found {len(rows)}")
    ids = {item.get("id") for item in rows}
    if ids != REQUIRED_IDS:
        issues.append(
            "row id mismatch: "
            f"missing={sorted(REQUIRED_IDS-ids)}, "
            f"extra={sorted(ids-REQUIRED_IDS)}"
        )
    for item in rows:
        role = item.get("role")
        expected_status = ROLE_STATUS.get(role)
        if expected_status is None:
            issues.append(f"{item.get('id')}: bad role {role!r}")
        elif item.get("status") != expected_status:
            issues.append(
                f"{item.get('id')}: role/status mismatch "
                f"{role!r}/{item.get('status')!r}"
            )
        for field in ("statement", "proof_boundary"):
            if not item.get(field):
                issues.append(f"{item.get('id')}: missing {field}")

    expected_audit = {
        "row_count": 32,
        "exact_reduction_count": 24,
        "literature_guard_count": 2,
        "proof_guard_count": 5,
        "open_two_component_gate_count": 1,
        "affine_dct_transference_proved": True,
        "centered_path_variance_proved": True,
        "local_mertens_abel_pair_proved": True,
        "full_path_norm_equivalence_proved": True,
        "low_projection_norm_equivalence_proved": False,
        "generic_zero_mean_first_mode_barrier_proved": True,
        "mobius_low_mode_gain_proved": False,
        "full_burnol_bound_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    audit = payload.get("audit", {})
    for key, value in expected_audit.items():
        if audit.get(key) != value:
            issues.append(
                f"audit {key}: expected {value!r}, "
                f"found {audit.get(key)!r}"
            )

    for needle in REQUIRED_NOTE:
        if needle not in note:
            issues.append(f"note missing {needle!r}")

    anchors = set(payload.get("source_anchors", []))
    for required in (
        "https://doi.org/10.1093/qmath/os-8.1.313",
        "https://doi.org/10.1112/jlms/s2-43.2.193",
        "https://doi.org/10.46298/hrj.2005.151",
        "https://arxiv.org/abs/2204.04613",
    ):
        if required not in anchors:
            issues.append(f"missing source anchor {required}")


def check_affine_transference(issues: list[str]) -> None:
    for size in (2, 3, 4, 8, 17, 32):
        basis = cosine_basis(size)
        path = [0.0] + [
            (
                (((29 * index + 7) % 31) - 15) / 19.0
                + 0.17 * math.sin(0.4 * index)
            )
            for index in range(1, size)
        ]
        incoming_tail = 0.37 - 0.011 * size
        tails = [incoming_tail - value for value in path]
        tau = [dot(path, mode) for mode in basis]
        coefficients = [dot(tails, mode) for mode in basis]

        if not close(
            coefficients[0],
            math.sqrt(size) * incoming_tail - tau[0],
            tolerance=3e-12,
        ):
            issues.append(f"affine mean transference failed at K={size}")
        for frequency in range(1, size):
            if not close(
                coefficients[frequency],
                -tau[frequency],
                tolerance=3e-12,
            ):
                issues.append(
                    f"nonconstant transference failed at "
                    f"K={size}, r={frequency}"
                )
                return
            formula = -math.sqrt(2.0 / size) * sum(
                path[index]
                * math.cos(
                    math.pi
                    * frequency
                    * (index + 0.5)
                    / size
                )
                for index in range(1, size)
            )
            if not close(
                coefficients[frequency],
                formula,
                tolerance=4e-12,
            ):
                issues.append(
                    f"path mode formula failed at "
                    f"K={size}, r={frequency}"
                )
                return

        left = dot(tails, tails)
        right = (
            abs(math.sqrt(size) * incoming_tail - tau[0]) ** 2
            + sum(value * value for value in tau[1:])
        )
        if not close(left, right, tolerance=4e-12):
            issues.append(f"affine Parseval failed at K={size}")

        centered = (
            dot(path, path)
            - sum(path) ** 2 / size
        )
        spectral = sum(value * value for value in tau[1:])
        if not close(centered, spectral, tolerance=4e-12):
            issues.append(f"centered variance failed at K={size}")


def check_deleted_dct_gram(issues: list[str]) -> None:
    for size in (2, 3, 4, 9, 16, 31):
        basis = cosine_basis(size)
        for left in range(1, size):
            for right in range(1, size):
                gram = sum(
                    basis[frequency][left]
                    * basis[frequency][right]
                    for frequency in range(1, size)
                )
                expected = (
                    (1.0 if left == right else 0.0)
                    - 1.0 / size
                )
                if not close(gram, expected, tolerance=4e-12):
                    issues.append(
                        f"deleted DCT Gram failed at "
                        f"K={size}, i={left}, j={right}"
                    )
                    return

        all_ones = [1.0] * (size - 1)
        transformed = [
            sum(
                basis[frequency][column] * all_ones[column - 1]
                for column in range(1, size)
            )
            for frequency in range(1, size)
        ]
        ratio_sq = dot(transformed, transformed) / dot(
            all_ones, all_ones
        )
        if not close(ratio_sq, 1.0 / size, tolerance=4e-12):
            issues.append(
                f"weak singular direction failed at K={size}"
            )

        if size > 2:
            zero_sum = [0.0] * (size - 1)
            zero_sum[0] = 1.0
            zero_sum[-1] = -1.0
            transformed = [
                sum(
                    basis[frequency][column]
                    * zero_sum[column - 1]
                    for column in range(1, size)
                )
                for frequency in range(1, size)
            ]
            ratio_sq = dot(transformed, transformed) / dot(
                zero_sum, zero_sum
            )
            if not close(ratio_sq, 1.0, tolerance=4e-12):
                issues.append(
                    f"unit singular direction failed at K={size}"
                )


def check_abel_pair_and_norms(issues: list[str]) -> None:
    for alpha in (0.07, 0.31, 0.74, 0.93):
        for size in (4, 8, 17, 32, 65):
            increments = [
                (
                    (((43 * index + 5) % 37) - 18) / 11.0
                    + 0.2 * math.cos(0.31 * index)
                )
                for index in range(1, size + 1)
            ]
            partial = []
            running = 0.0
            for value in increments:
                running += value
                partial.append(running)

            weights = [
                (size + index) ** (-1.0 - alpha)
                for index in range(1, size + 1)
            ]
            weighted_path = []
            running = 0.0
            for value, weight in zip(increments, weights):
                running += value * weight
                weighted_path.append(running)

            for index in range(size):
                forward = weights[index] * partial[index]
                forward += sum(
                    partial[ell]
                    * (weights[ell] - weights[ell + 1])
                    for ell in range(index)
                )
                if not close(
                    weighted_path[index],
                    forward,
                    tolerance=5e-12,
                ):
                    issues.append(
                        f"Abel forward failed at alpha={alpha}, "
                        f"K={size}, m={index + 1}"
                    )
                    return

                inverse = weighted_path[index] / weights[index]
                inverse += sum(
                    weighted_path[ell]
                    * (
                        1.0 / weights[ell]
                        - 1.0 / weights[ell + 1]
                    )
                    for ell in range(index)
                )
                if not close(
                    partial[index],
                    inverse,
                    tolerance=7e-11,
                ):
                    issues.append(
                        f"Abel inverse failed at alpha={alpha}, "
                        f"K={size}, m={index + 1}"
                    )
                    return

            s_vector = partial[:-1]
            t_vector = weighted_path[:-1]
            s_norm = norm(s_vector)
            t_norm = norm(t_vector)
            upper = (
                (2.0 + alpha)
                * size ** (-1.0 - alpha)
                * s_norm
            )
            lower = (
                size ** (-1.0 - alpha)
                * s_norm
                / (2.0**alpha * (3.0 + alpha))
            )
            if t_norm > upper * (1.0 + 2e-12):
                issues.append(
                    f"path norm upper failed at alpha={alpha}, K={size}"
                )
            if t_norm < lower * (1.0 - 2e-12):
                issues.append(
                    f"path norm lower failed at alpha={alpha}, K={size}"
                )


def check_countermodel(issues: list[str]) -> None:
    for alpha in (0.08, 0.37, 0.69, 0.91):
        for size in (4, 8, 16, 32, 128, 512):
            amplitude = (
                2.0 ** (-2.0 - alpha)
                / math.pi
                * size ** (-alpha)
            )
            path = [
                amplitude
                * math.sin(2.0 * math.pi * index / size)
                for index in range(size + 1)
            ]
            coefficients = [
                path[index] - path[index - 1]
                for index in range(1, size + 1)
            ]
            for index, value in enumerate(coefficients, start=1):
                envelope = (size + index) ** (-1.0 - alpha)
                if abs(value) > envelope * (1.0 + 2e-12):
                    issues.append(
                        f"countermodel envelope failed at "
                        f"alpha={alpha}, K={size}, m={index}"
                    )
                    return
            if not close(sum(coefficients), 0.0, tolerance=2e-12):
                issues.append(
                    f"countermodel block total failed at "
                    f"alpha={alpha}, K={size}"
                )
            if sum(abs(value) for value in coefficients) > (
                2.0 * math.pi * amplitude * (1.0 + 2e-12)
            ):
                issues.append(
                    f"countermodel l1 bound failed at "
                    f"alpha={alpha}, K={size}"
                )

            basis = cosine_basis(size)
            dct_path = path[:-1]
            tau_zero = dot(dct_path, basis[0])
            if not close(tau_zero, 0.0, tolerance=3e-12):
                issues.append(
                    f"countermodel mean failed at "
                    f"alpha={alpha}, K={size}"
                )
            first_mode = abs(-dot(dct_path, basis[1]))
            angle = math.pi / size
            exact = (
                amplitude
                * math.sqrt(2.0 / size)
                * math.cos(angle)
                / 2.0
                * (
                    1.0 / math.sin(angle / 2.0)
                    + 1.0 / math.sin(3.0 * angle / 2.0)
                )
            )
            if not close(first_mode, exact, tolerance=5e-12):
                issues.append(
                    f"countermodel first-mode formula failed at "
                    f"alpha={alpha}, K={size}"
                )
            lower_mode = amplitude * math.sqrt(size) / math.pi
            if first_mode < lower_mode * (1.0 - 3e-12):
                issues.append(
                    f"countermodel first-mode lower bound failed at "
                    f"alpha={alpha}, K={size}"
                )
            weighted = size**alpha * first_mode**2
            weighted_lower = (
                2.0 ** (-4.0 - 2.0 * alpha)
                * math.pi ** (-4.0)
                * size ** (1.0 - alpha)
            )
            if weighted < weighted_lower * (1.0 - 4e-12):
                issues.append(
                    f"countermodel weighted divergence scale failed at "
                    f"alpha={alpha}, K={size}"
                )


def check_exponent_guard(issues: list[str]) -> None:
    for epsilon in (0.0, 0.01, 0.07):
        sigma = 0.75 + epsilon
        critical_alpha = 2.0 * sigma - 1.0
        if not close(
            critical_alpha,
            0.5 + 2.0 * epsilon,
            tolerance=1e-14,
        ):
            issues.append(
                f"three-quarter threshold failed at epsilon={epsilon}"
            )
        alpha = critical_alpha + 0.05
        exponent = 2.0 * sigma - 1.0 - alpha
        if exponent >= 0.0:
            issues.append(
                f"three-quarter summable side failed at epsilon={epsilon}"
            )
        alpha = max(0.01, critical_alpha - 0.05)
        exponent = 2.0 * sigma - 1.0 - alpha
        if exponent <= 0.0:
            issues.append(
                f"three-quarter nonsummable side failed at "
                f"epsilon={epsilon}"
            )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = json.loads(args.result.read_text(encoding="utf-8"))
    note = args.note.read_text(encoding="utf-8")
    issues: list[str] = []
    check_schema(payload, note, issues)
    check_affine_transference(issues)
    check_deleted_dct_gram(issues)
    check_abel_pair_and_norms(issues)
    check_countermodel(issues)
    check_exponent_guard(issues)
    for issue in issues:
        print(f"MERTENS-LOCAL-PATH {issue}")
    audit = payload.get("audit", {})
    print(
        "validated Mertens local-path cosine transference: "
        f"{len(payload.get('rows', []))} rows, {len(issues)} issues, "
        f"{audit.get('exact_reduction_count', 0)} exact reductions, "
        f"{audit.get('proof_guard_count', 0)} proof guards, "
        f"{audit.get('open_two_component_gate_count', 0)} open gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
