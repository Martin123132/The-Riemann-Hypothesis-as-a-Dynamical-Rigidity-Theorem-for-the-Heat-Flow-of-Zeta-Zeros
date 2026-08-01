#!/usr/bin/env python3
"""Validate the joined planar boundary-flux and anchor handoff."""

from __future__ import annotations

import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.py"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"


def close(left: float, right: float, tolerance: float = 3e-9) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def sinc(value: float) -> float:
    if value == 0.0:
        return 1.0
    return math.sin(value) / value


def mobius_values(limit: int) -> list[int]:
    values = [0] * (limit + 1)
    values[1] = 1
    primes: list[int] = []
    composite = [False] * (limit + 1)
    for value in range(2, limit + 1):
        if not composite[value]:
            primes.append(value)
            values[value] = -1
        for prime in primes:
            product = value * prime
            if product > limit:
                break
            composite[product] = True
            if value % prime == 0:
                values[product] = 0
                break
            values[product] = -values[value]
    return values


def q_value(mode: int) -> int:
    return 2 * mode - 1


def theta(size: int, mode: int) -> float:
    return q_value(mode) * math.pi / (2 * size + 1)


def phi(size: int, mode: int, index: int) -> float:
    return (
        2.0
        / math.sqrt(2 * size + 1)
        * math.sin(index * theta(size, mode))
    )


def future_weight(alpha: float, size: int, value: int) -> float:
    if value >= 4 * size:
        return 0.0
    if value < 2 * size:
        raise ValueError("future weight requested before 2K")
    return (
        (2.0 * size / value) ** (1.0 + alpha)
        * (4.0 * size - value)
        / (2.0 * size)
    )


def decrement(alpha: float, size: int, value: int) -> float:
    return (
        future_weight(alpha, size, value)
        - future_weight(alpha, size, value + 1)
    )


def feature(
    alpha: float,
    size: int,
    mode: int,
    value: int,
) -> float:
    if size < value < 2 * size:
        return phi(size, mode, value - size)
    if 2 * size <= value < 4 * size:
        return future_weight(alpha, size, value) * phi(size, mode, size)
    return 0.0


def kernel(size: int, rank: int, left: int, right: int) -> float:
    return sum(
        phi(size, mode, left)
        * phi(size, mode, right)
        / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )


def ordinary_kernel(
    alpha: float,
    size: int,
    rank: int,
    left: int,
    right: int,
) -> float:
    return sum(
        feature(alpha, size, mode, left)
        * feature(alpha, size, mode, right)
        / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )


def mixed_difference(
    alpha: float,
    size: int,
    rank: int,
    left: int,
    shift: int,
) -> float:
    right = left + shift
    return (
        ordinary_kernel(alpha, size, rank, left, right)
        - ordinary_kernel(alpha, size, rank, left + 1, right + 1)
        - ordinary_kernel(alpha, size, rank, left, right + 1)
        + ordinary_kernel(alpha, size, rank, left + 1, right + 2)
    )


def prefix_current(
    current: list[float],
    endpoint: int,
    width: int,
) -> float:
    return sum(
        current[base - 1] * current[base + gap - 1]
        for base in range(1, endpoint + 1)
        for gap in range(1, width + 1)
    )


def check_schema(
    payload: dict,
    note_text: str,
    formal_text: str,
    issues: list[str],
) -> None:
    expected_kind = (
        "jensen_window_pf_mertens_planar_"
        "boundary_flux_anchor_handoff"
    )
    if payload.get("kind") != expected_kind:
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    if payload.get("status") != (
        "exact joined boundary-flux and rank-one anchor reduction "
        "with one open Mobius energy gate"
    ):
        issues.append("result status mismatch")

    rows = payload.get("rows", [])
    summary = payload.get("summary", {})
    if len(rows) != 35 or summary.get("row_count") != 35:
        issues.append("expected exactly 35 ledger rows")
    row_ids = [item.get("id") for item in rows]
    if len(row_ids) != len(set(row_ids)):
        issues.append("duplicate row id")
    for index in range(1, 36):
        prefix = f"pbfah_{index:02d}_"
        if not any(str(row_id).startswith(prefix) for row_id in row_ids):
            issues.append(f"missing row prefix {prefix}")
    for key, expected in (
        ("exact_reduction_count", 29),
        ("conditional_calibration_count", 1),
        ("proof_guard_count", 2),
        ("countermodel_count", 1),
        ("finite_validation_count", 1),
        ("open_joint_energy_gate_count", 1),
    ):
        if summary.get(key) != expected:
            issues.append(f"{key} mismatch")
    for key in (
        "global_flux_proved",
        "all_remaining_cells_projected",
        "interface_joins_proved",
        "shoulder_join_proved",
        "block_compression_proved",
        "boundary_cancellation_witness_proved",
        "orthogonal_anchor_completion_proved",
    ):
        if summary.get(key) is not True:
            issues.append(f"{key} must be true")
    for key in (
        "interior_projection_estimate_proved",
        "joint_mobius_energy_proved",
        "full_reciprocal_tail_proved",
        "rh_proved",
        "pf_infinity_proved",
        "lambda_le_zero_proved",
    ):
        if summary.get(key) is not False:
            issues.append(f"{key} must be false")

    role_status = {
        item.get("id"): (item.get("role"), item.get("status"))
        for item in rows
    }
    if role_status.get("pbfah_31_joined_energy_gate") != (
        "theorem_target",
        "open",
    ):
        issues.append("joined energy gate must remain open")
    if role_status.get("pbfah_32_boundary_cancellation_witness") != (
        "countermodel_gate",
        "guard_validated",
    ):
        issues.append("boundary witness status mismatch")

    for source in payload.get("source", {}).values():
        if not (REPO_ROOT / source).exists():
            issues.append(f"missing source {source}")

    for marker in (
        "Status: exact joined boundary-flux",
        "This is not a proof of RH",
        "One Modewise Flux Law",
        "(PBFA.1)",
        "(PBFA.6)",
        "(PBFA.10)",
        "(PBFA.13)",
        "Exact Shoulder-Omission Witness",
        "B_CC=198pi^2/19^2",
        "A_CF=198pi^2/19^2",
        "Finite Mobius Diagnostic",
        "2.086e-5",
        "shoulder",
        "No estimate for (PBFA.14)",
    ):
        if marker not in note_text:
            issues.append(f"note missing marker {marker!r}")

    for marker in (
        "Corollary 11.22Z.8",
        "(11.22Z.120)",
        "(11.22Z.121)",
        "(11.22Z.122)",
        "(11.22Z.123)",
        "(11.22Z.124)",
        "(11.22Z.125)",
        "(11.22Z.126)",
        "(11.22Z.127)",
        "(11.22Z.128)",
        "(11.22Z.129)",
        "jensen_window_pf_mertens_planar_boundary_flux_anchor_handoff.md",
    ):
        if marker not in formal_text:
            issues.append(f"formal core missing marker {marker!r}")


def check_flux_cells(issues: list[str]) -> int:
    audits = 0
    for alpha, size in ((0.19, 3), (0.47, 5), (0.81, 8)):
        for rank in sorted({1, max(1, size // 2), size}):
            for left in range(size + 1, 4 * size - 1):
                for right in range(left + 1, 4 * size):
                    shift = right - left
                    direct = mixed_difference(
                        alpha,
                        size,
                        rank,
                        left,
                        shift,
                    )
                    flux = 0.0
                    for mode in range(1, rank + 1):
                        delta_right = (
                            feature(alpha, size, mode, right)
                            - feature(alpha, size, mode, right + 1)
                        )
                        delta_next = (
                            feature(alpha, size, mode, right + 1)
                            - feature(alpha, size, mode, right + 2)
                        )
                        flux += (
                            feature(alpha, size, mode, left)
                            * delta_right
                            - feature(alpha, size, mode, left + 1)
                            * delta_next
                        ) / q_value(mode) ** 2
                    if not close(direct, flux, 2e-11):
                        issues.append(
                            f"global flux mismatch alpha={alpha}, "
                            f"K={size}, R={rank}, n={left}, m={right}"
                        )
                    audits += 1

                    if not (size < left < 2 * size):
                        continue
                    local = left - size
                    expected = None
                    if right == 2 * size - 1:
                        expected = sum(
                            (
                                phi(size, mode, local)
                                * (
                                    phi(size, mode, size - 1)
                                    - phi(size, mode, size)
                                )
                                - phi(size, mode, local + 1)
                                * phi(size, mode, size)
                                * decrement(alpha, size, 2 * size)
                            )
                            / q_value(mode) ** 2
                            for mode in range(1, rank + 1)
                        )
                    elif left <= 2 * size - 2 and right >= 2 * size:
                        expected = sum(
                            phi(size, mode, size)
                            * (
                                (
                                    phi(size, mode, local)
                                    - phi(size, mode, local + 1)
                                )
                                * decrement(alpha, size, right)
                                + phi(size, mode, local + 1)
                                * (
                                    decrement(alpha, size, right)
                                    - decrement(alpha, size, right + 1)
                                )
                            )
                            / q_value(mode) ** 2
                            for mode in range(1, rank + 1)
                        )
                    elif left == 2 * size - 1 and right >= 2 * size:
                        expected = sum(
                            (
                                phi(size, mode, size)
                                * (
                                    phi(size, mode, size - 1)
                                    - phi(size, mode, size)
                                )
                                * decrement(alpha, size, right)
                                + phi(size, mode, size) ** 2
                                * (
                                    decrement(alpha, size, right)
                                    - decrement(alpha, size, right + 1)
                                )
                            )
                            / q_value(mode) ** 2
                            for mode in range(1, rank + 1)
                        )
                    if expected is not None:
                        if not close(direct, expected, 2e-10):
                            issues.append(
                                f"cell specialization mismatch K={size}, "
                                f"R={rank}, n={left}, m={right}"
                            )
                        audits += 1

            for left in range(2 * size, 4 * size - 1):
                for right in range(left + 1, 4 * size):
                    direct = mixed_difference(
                        alpha,
                        size,
                        rank,
                        left,
                        right - left,
                    )
                    expected = sum(
                        phi(size, mode, size) ** 2
                        * (
                            future_weight(alpha, size, left)
                            * (
                                decrement(alpha, size, right)
                                - decrement(alpha, size, right + 1)
                            )
                            + decrement(alpha, size, left)
                            * decrement(alpha, size, right + 1)
                        )
                        / q_value(mode) ** 2
                        for mode in range(1, rank + 1)
                    )
                    if not close(direct, expected, 2e-10):
                        issues.append(
                            f"future cell mismatch K={size}, "
                            f"R={rank}, n={left}, m={right}"
                        )
                    audits += 1
    return audits


def check_endpoint_identities(issues: list[str]) -> int:
    audits = 0
    for size in range(2, 25):
        normalizer = 2 * size + 1
        for mode in range(1, size + 1):
            frequency = q_value(mode)
            angle = theta(size, mode)
            endpoint = phi(size, mode, size)
            if not close(phi(size, mode, size + 1), endpoint, 2e-12):
                issues.append(f"endpoint symmetry failed K={size}, r={mode}")
            recurrence = (2 * math.cos(angle) - 1) * endpoint
            if not close(phi(size, mode, size - 1), recurrence, 2e-12):
                issues.append(f"endpoint recurrence failed K={size}, r={mode}")
            jump = phi(size, mode, size - 1) - endpoint
            smoothed = (
                -math.pi**2
                * frequency**2
                / normalizer**2
                * sinc(angle / 2.0) ** 2
                * endpoint
            )
            if not close(jump, smoothed, 2e-12):
                issues.append(f"endpoint smoothing failed K={size}, r={mode}")
            audits += 3
    return audits


def deterministic_vectors(size: int, mu: list[int]) -> list[list[float]]:
    return [
        [float(mu[size + index]) for index in range(1, 4 * size - size)],
        [
            float(((7 * index + 3 * size) % 5) - 2)
            for index in range(1, 3 * size)
        ],
    ]


def check_block_case(
    alpha: float,
    size: int,
    rank: int,
    values: list[float],
    issues: list[str],
) -> int:
    current = values[: size - 1]
    future = values[size - 1 :]
    p_values = [
        kernel(size, rank, index, size)
        for index in range(1, size)
    ]
    c_value = kernel(size, rank, size, size)
    b_values = [
        future_weight(alpha, size, value)
        for value in range(2 * size, 4 * size)
    ]

    current_direct = sum(
        current[left - 1]
        * current[right - 1]
        * kernel(size, rank, left, right)
        for left in range(1, size)
        for right in range(left + 1, size)
    )
    cross_direct = sum(
        current[index - 1]
        * future[offset]
        * p_values[index - 1]
        * b_values[offset]
        for index in range(1, size)
        for offset in range(2 * size)
    )
    future_direct = sum(
        future[left]
        * future[right]
        * b_values[left]
        * b_values[right]
        * c_value
        for left in range(2 * size)
        for right in range(left + 1, 2 * size)
    )

    mode_x = [
        sum(
            current[index - 1] * phi(size, mode, index)
            for index in range(1, size)
        )
        for mode in range(1, rank + 1)
    ]
    mode_diagonal = [
        sum(
            current[index - 1] ** 2 * phi(size, mode, index) ** 2
            for index in range(1, size)
        )
        for mode in range(1, rank + 1)
    ]
    beta = sum(
        future[offset] * b_values[offset]
        for offset in range(2 * size)
    )
    q_future = sum(
        future[offset] ** 2 * b_values[offset] ** 2
        for offset in range(2 * size)
    )
    gamma = sum(
        current[index - 1] * p_values[index - 1]
        for index in range(1, size)
    )
    current_spectral = 0.5 * sum(
        (
            mode_x[mode - 1] ** 2
            - mode_diagonal[mode - 1]
        )
        / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )
    cross_compressed = beta * gamma
    future_compressed = 0.5 * c_value * (beta**2 - q_future)

    for name, direct, compressed in (
        ("CC", current_direct, current_spectral),
        ("CF", cross_direct, cross_compressed),
        ("FF", future_direct, future_compressed),
    ):
        if not close(direct, compressed, 3e-10):
            issues.append(
                f"{name} compression failed alpha={alpha}, "
                f"K={size}, R={rank}"
            )

    boundary_prefixes = [
        prefix_current(current, index, size - 1 - index)
        for index in range(1, size - 1)
    ]
    boundary = sum(
        boundary_prefixes[index - 1]
        * kernel(size, rank, index, size - 1)
        for index in range(1, size - 1)
    )
    shoulder_prefixes = [
        prefix_current(current, index, size - 2 - index)
        for index in range(1, size - 2)
    ]
    shoulder = sum(
        shoulder_prefixes[index - 1] * p_values[index]
        for index in range(1, size - 2)
    )
    interior = 0.0
    for index in range(1, size - 1):
        for width in range(1, size - 1 - index):
            right = index + width
            curvature = (
                kernel(size, rank, index, right)
                - kernel(size, rank, index + 1, right + 1)
                - kernel(size, rank, index, right + 1)
                + kernel(size, rank, index + 1, right + 2)
            )
            interior += (
                prefix_current(current, index, width) * curvature
            )
    if not close(
        current_direct,
        interior + boundary - shoulder,
        4e-10,
    ):
        issues.append(
            f"current Abel reconstruction failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )

    mode_u = [
        sum(
            boundary_prefixes[index - 1]
            * phi(size, mode, index)
            for index in range(1, size - 1)
        )
        for mode in range(1, rank + 1)
    ]
    mode_v = [
        sum(
            shoulder_prefixes[index - 1]
            * phi(size, mode, index + 1)
            for index in range(1, size - 2)
        )
        for mode in range(1, rank + 1)
    ]
    boundary_spectral = sum(
        phi(size, mode, size - 1)
        * mode_u[mode - 1]
        / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )
    if not close(boundary, boundary_spectral, 3e-10):
        issues.append(
            f"boundary spectrum failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )
    shoulder_spectral = sum(
        phi(size, mode, size)
        * mode_v[mode - 1]
        / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )
    if not close(shoulder, shoulder_spectral, 3e-10):
        issues.append(
            f"shoulder spectrum failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )

    z_values = [
        (
            phi(size, mode, size - 1) * mode_u[mode - 1]
            - phi(size, mode, size) * mode_v[mode - 1]
            + beta * phi(size, mode, size) * mode_x[mode - 1]
            + 0.5
            * phi(size, mode, size) ** 2
            * (beta**2 - q_future)
        )
        for mode in range(1, rank + 1)
    ]
    remainder = boundary - shoulder + cross_direct + future_direct
    remainder_spectral = sum(
        z_values[mode - 1] / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )
    if not close(remainder, remainder_spectral, 4e-10):
        issues.append(
            f"remainder spectrum failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )
    total_direct = current_direct + cross_direct + future_direct
    if not close(total_direct, interior + remainder_spectral, 4e-10):
        issues.append(
            f"joined direct decomposition failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )

    cauchy_right = sum(
        1.0 / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    ) * sum(
        z_values[mode - 1] ** 2 / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )
    if remainder_spectral**2 > cauchy_right + 3e-10:
        issues.append(
            f"remainder Cauchy failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )

    current_energy = sum(
        mode_x[mode - 1] ** 2 / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )
    low_energy = sum(
        (
            mode_x[mode - 1]
            + beta * phi(size, mode, size)
        )
        ** 2
        / q_value(mode) ** 2
        for mode in range(1, rank + 1)
    )
    perpendicular = current_energy - gamma**2 / c_value
    completed = perpendicular + c_value * (beta + gamma / c_value) ** 2
    if perpendicular < -3e-10:
        issues.append(
            f"negative perpendicular energy alpha={alpha}, "
            f"K={size}, R={rank}"
        )
    if not close(low_energy, completed, 4e-10):
        issues.append(
            f"orthogonal completion failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )

    diagonal = sum(
        current[index - 1] ** 2
        * kernel(size, rank, index, index)
        for index in range(1, size)
    ) + c_value * q_future
    if not close(low_energy, diagonal + 2.0 * total_direct, 4e-10):
        issues.append(
            f"energy/offdiagonal reconciliation failed alpha={alpha}, "
            f"K={size}, R={rank}"
        )
    return 12


def check_block_compression(issues: list[str]) -> int:
    audits = 0
    mu = mobius_values(100)
    for alpha, size in ((0.23, 4), (0.51, 6), (0.86, 9)):
        for rank in sorted({1, max(1, size // 2), size}):
            for values in deterministic_vectors(size, mu):
                audits += check_block_case(
                    alpha,
                    size,
                    rank,
                    values,
                    issues,
                )
    return audits


def check_boundary_witness(issues: list[str]) -> int:
    size = 9
    current = [-1.0] * 6 + [1.0, 1.0]
    scaled_current = sum(
        current[left - 1] * current[right - 1] * left
        for left in range(1, size)
        for right in range(left + 1, size)
    )
    prefixes = [
        int(prefix_current(current, index, size - 1 - index))
        for index in range(1, size - 1)
    ]
    scaled_boundary = sum(
        index * prefixes[index - 1]
        for index in range(1, size - 1)
    )
    shoulder_prefixes = [
        int(prefix_current(current, index, size - 2 - index))
        for index in range(1, size - 2)
    ]
    scaled_shoulder = sum(
        (index + 1) * shoulder_prefixes[index - 1]
        for index in range(1, size - 2)
    )
    if scaled_current != 0:
        issues.append("K=9 witness current offdiagonal is not zero")
    if prefixes != [3, 6, 9, 10, 9, 6, 5]:
        issues.append(f"K=9 witness prefixes mismatch: {prefixes}")
    if scaled_boundary != 198:
        issues.append(
            f"K=9 witness boundary mismatch: {scaled_boundary}"
        )
    if shoulder_prefixes != [4, 8, 10, 10, 8, 4]:
        issues.append(
            f"K=9 witness shoulder prefixes mismatch: "
            f"{shoulder_prefixes}"
        )
    if scaled_shoulder != 198:
        issues.append(
            f"K=9 witness shoulder mismatch: {scaled_shoulder}"
        )
    if scaled_current - scaled_boundary + scaled_shoulder != 0:
        issues.append("K=9 witness joined interior mismatch")
    return 6


def recompute_mobius_row(alpha: float, size: int) -> dict[str, float | int]:
    rank = math.ceil(size ** (1.0 - alpha / 2.0))
    modes = list(range(1, rank + 1))
    mu = mobius_values(4 * size + 2)
    current = [float(mu[size + index]) for index in range(1, size)]
    future = [
        float(mu[value]) for value in range(2 * size, 4 * size)
    ]
    mode_phi = [
        [phi(size, mode, index) for mode in modes]
        for index in range(1, size)
    ]
    endpoint = [phi(size, mode, size) for mode in modes]
    penultimate = [phi(size, mode, size - 1) for mode in modes]
    current_prefix = [0.0]
    for value in current:
        current_prefix.append(current_prefix[-1] + value)

    def fast_prefix(endpoint_index: int, width: int) -> float:
        return sum(
            current[base - 1]
            * (
                current_prefix[base + width]
                - current_prefix[base]
            )
            for base in range(1, endpoint_index + 1)
        )

    mode_x = [
        sum(
            current[index - 1] * mode_phi[index - 1][mode - 1]
            for index in range(1, size)
        )
        for mode in modes
    ]
    mode_diagonal = [
        sum(
            current[index - 1] ** 2
            * mode_phi[index - 1][mode - 1] ** 2
            for index in range(1, size)
        )
        for mode in modes
    ]
    current_offdiagonal = 0.5 * sum(
        (
            mode_x[mode - 1] ** 2
            - mode_diagonal[mode - 1]
        )
        / q_value(mode) ** 2
        for mode in modes
    )
    boundary_column = [
        sum(
            mode_phi[index - 1][mode - 1]
            * penultimate[mode - 1]
            / q_value(mode) ** 2
            for mode in modes
        )
        for index in range(1, size)
    ]
    boundary = sum(
        fast_prefix(index, size - 1 - index)
        * boundary_column[index - 1]
        for index in range(1, size - 1)
    )
    p_values = [
        sum(
            mode_phi[index - 1][mode - 1]
            * endpoint[mode - 1]
            / q_value(mode) ** 2
            for mode in modes
        )
        for index in range(1, size)
    ]
    shoulder = sum(
        fast_prefix(index, size - 2 - index)
        * p_values[index]
        for index in range(1, size - 2)
    )
    joined_boundary = boundary - shoulder
    interior = current_offdiagonal - joined_boundary
    b_values = [
        future_weight(alpha, size, value)
        for value in range(2 * size, 4 * size)
    ]
    beta = sum(
        future[offset] * b_values[offset]
        for offset in range(2 * size)
    )
    q_future = sum(
        future[offset] ** 2 * b_values[offset] ** 2
        for offset in range(2 * size)
    )
    gamma = sum(
        current[index - 1] * p_values[index - 1]
        for index in range(1, size)
    )
    c_value = sum(
        endpoint[mode - 1] ** 2 / q_value(mode) ** 2
        for mode in modes
    )
    anchor = beta * gamma + 0.5 * c_value * (beta**2 - q_future)
    return {
        "K": size,
        "R": rank,
        "O_CC": current_offdiagonal,
        "B_CC": boundary,
        "A_CF": shoulder,
        "joined_boundary": joined_boundary,
        "O_CC_int": interior,
        "anchor_correction": anchor,
        "O_total": current_offdiagonal + anchor,
    }


def check_mobius_table(payload: dict, issues: list[str]) -> int:
    audit = payload.get("mobius_audit", {})
    if audit.get("alpha") != 0.5:
        issues.append("Mobius audit alpha mismatch")
        return 0
    rows = audit.get("rows", [])
    if [item.get("K") for item in rows] != [
        16,
        32,
        64,
        128,
        256,
        512,
        1024,
    ]:
        issues.append("Mobius audit K grid mismatch")
        return 0
    audits = 0
    for stored in rows:
        computed = recompute_mobius_row(0.5, int(stored["K"]))
        for key in (
            "R",
            "O_CC",
            "B_CC",
            "A_CF",
            "joined_boundary",
            "O_CC_int",
            "anchor_correction",
            "O_total",
        ):
            if not close(
                float(stored[key]),
                float(computed[key]),
                2e-10,
            ):
                issues.append(
                    f"Mobius audit mismatch K={stored['K']}, key={key}"
                )
            audits += 1
        if not close(
            float(stored["O_CC"]),
            float(stored["O_CC_int"])
            + float(stored["B_CC"])
            - float(stored["A_CF"]),
            1e-9,
        ):
            issues.append(
                f"stored current reconstruction mismatch K={stored['K']}"
            )
        audits += 1
        if not close(
            float(stored["joined_boundary"]),
            float(stored["B_CC"]) - float(stored["A_CF"]),
            1e-9,
        ):
            issues.append(
                f"stored boundary join mismatch K={stored['K']}"
            )
        audits += 1
    return audits


def check_reproduction(issues: list[str]) -> None:
    with tempfile.TemporaryDirectory() as temporary:
        root = Path(temporary)
        result = root / "result.json"
        note = root / "note.md"
        process = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--output",
                str(result),
                "--note",
                str(note),
            ],
            cwd=REPO_ROOT,
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        if process.returncode:
            issues.append(
                "builder reproduction failed: "
                + (process.stderr.strip() or process.stdout.strip())
            )
            return
        if result.read_bytes() != RESULT.read_bytes():
            issues.append("builder result reproduction mismatch")
        if note.read_bytes() != NOTE.read_bytes():
            issues.append("builder note reproduction mismatch")


def main() -> int:
    issues: list[str] = []
    for path in (RESULT, NOTE, BUILDER, FORMAL_CORE):
        if not path.exists():
            issues.append(f"missing file {path}")
    if issues:
        for issue in issues:
            print(f"PBFAH [error] {issue}", file=sys.stderr)
        return 1
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        note_text = NOTE.read_text(encoding="utf-8")
        formal_text = FORMAL_CORE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as error:
        print(f"PBFAH [load-failure] {error}", file=sys.stderr)
        return 1

    check_schema(payload, note_text, formal_text, issues)
    flux_audits = check_flux_cells(issues)
    endpoint_audits = check_endpoint_identities(issues)
    block_audits = check_block_compression(issues)
    witness_audits = check_boundary_witness(issues)
    mobius_audits = check_mobius_table(payload, issues)
    check_reproduction(issues)

    if issues:
        for issue in issues:
            print(f"PBFAH [error] {issue}", file=sys.stderr)
        return 1
    print(
        "validated planar boundary-flux/anchor handoff: "
        f"{len(payload['rows'])} rows, 0 issues, "
        f"{flux_audits} flux audits, "
        f"{endpoint_audits} endpoint audits, "
        f"{block_audits} block audits, "
        f"{witness_audits} witness audits, "
        f"{mobius_audits} Mobius audits, "
        "1 open joined energy gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
