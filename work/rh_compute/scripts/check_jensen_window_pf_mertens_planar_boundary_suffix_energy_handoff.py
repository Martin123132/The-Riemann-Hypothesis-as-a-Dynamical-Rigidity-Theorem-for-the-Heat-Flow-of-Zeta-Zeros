#!/usr/bin/env python3
"""Validate the planar boundary suffix-energy handoff."""

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
    "jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.py"
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


def kernel(size: int, rank: int, left: int, right: int) -> float:
    return sum(
        phi(size, mode, left)
        * phi(size, mode, right)
        / q_value(mode) ** 2
        for mode in range(1, rank + 1)
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


def suffixes(current: list[float]) -> list[float]:
    size = len(current) + 1
    result = [0.0] * (size + 1)
    for index in range(size - 1, 0, -1):
        result[index] = result[index + 1] + current[index - 1]
    return result


def check_schema(
    payload: dict,
    note_text: str,
    formal_text: str,
    issues: list[str],
) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_planar_"
        "boundary_suffix_energy_handoff"
    ):
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    if payload.get("status") != (
        "exact positive suffix-energy reduction for the complete "
        "joined boundary/future remainder with one open anchored "
        "Mobius-energy gate"
    ):
        issues.append("result status mismatch")

    rows = payload.get("rows", [])
    summary = payload.get("summary", {})
    if len(rows) != 26 or summary.get("row_count") != 26:
        issues.append("expected exactly 26 ledger rows")
    row_ids = [item.get("id") for item in rows]
    if len(row_ids) != len(set(row_ids)):
        issues.append("duplicate row id")
    for index in range(1, 27):
        prefix = f"pbseh_{index:02d}_"
        if not any(str(row_id).startswith(prefix) for row_id in row_ids):
            issues.append(f"missing row prefix {prefix}")
    for key, expected in (
        ("exact_reduction_count", 23),
        ("conditional_calibration_count", 0),
        ("proof_guard_count", 1),
        ("finite_validation_count", 1),
        ("open_energy_gate_count", 1),
    ):
        if summary.get(key) != expected:
            issues.append(f"{key} mismatch")
    for key in (
        "boundary_suffix_split_proved",
        "endpoint_correction_logarithmic",
        "positive_suffix_energy_proved",
        "subpower_equivalence_proved",
    ):
        if summary.get(key) is not True:
            issues.append(f"{key} must be true")
    for key in (
        "anchored_mobius_energy_proved",
        "interior_projection_estimate_proved",
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
    if role_status.get("pbseh_25_open_energy_gate") != (
        "theorem_target",
        "open",
    ):
        issues.append("anchored energy gate must remain open")
    if role_status.get("pbseh_22_block_mean_square_sufficient") != (
        "exact_reduction",
        "available_exact",
    ):
        issues.append("affine-tent domination status mismatch")

    for source in payload.get("source", {}).values():
        if not (REPO_ROOT / source).exists():
            issues.append(f"missing source {source}")

    for marker in (
        "Status: exact positive suffix-energy reduction",
        "This is not a proof of RH",
        "Boundary And Shoulder As One Edge Kernel",
        "(PBSE.2)",
        "(PBSE.8)",
        "(PBSE.13)",
        "(PBSE.15)",
        "Positive Anchored Suffix Energy",
        "R_B=O_epsilon(K^epsilon)",
        "E_B=O_epsilon(K^epsilon)",
        "(PBSE.18)",
        "(MATBH.39)",
        "0.097890613",
        "No estimate for `H_K`",
    ):
        if marker not in note_text:
            issues.append(f"note missing marker {marker!r}")

    for marker in (
        "Corollary 11.22Z.9",
        "(11.22Z.130)",
        "(11.22Z.131)",
        "(11.22Z.132)",
        "(11.22Z.133)",
        "(11.22Z.134)",
        "(11.22Z.135)",
        "(11.22Z.136)",
        "(11.22Z.137)",
        "(11.22Z.138)",
        "(11.22Z.139)",
        "jensen_window_pf_mertens_planar_boundary_suffix_energy_handoff.md",
    ):
        if marker not in formal_text:
            issues.append(f"formal core missing marker {marker!r}")


def deterministic_vectors(size: int, mu: list[int]) -> list[list[float]]:
    return [
        [float(mu[size + index]) for index in range(1, 4 * size - size)],
        [
            float(((5 * index + 2 * size) % 3) - 1)
            for index in range(1, 3 * size)
        ],
    ]


def check_case(
    alpha: float,
    size: int,
    rank: int,
    values: list[float],
    issues: list[str],
) -> tuple[int, int, int]:
    current = values[: size - 1]
    future = values[size - 1 :]
    normalizer = 2 * size + 1
    odd_harmonic = sum(
        1.0 / q_value(mode)
        for mode in range(1, rank + 1)
    )
    w = [0.0] + [
        kernel(size, rank, index, size - 1)
        for index in range(1, size)
    ]
    p = [0.0] + [
        kernel(size, rank, index, size)
        for index in range(1, size + 1)
    ]
    c_value = p[size]
    delta = [0.0] + [
        w[index] - p[index]
        for index in range(1, size)
    ]
    current_suffix = suffixes(current)
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
        current[index - 1] * p[index]
        for index in range(1, size)
    )

    boundary_prefix = [0.0] * size
    for index in range(1, size - 1):
        boundary_prefix[index] = prefix_current(
            current,
            index,
            size - 1 - index,
        )
    shoulder_prefix = [0.0] * size
    for index in range(1, size - 2):
        shoulder_prefix[index] = prefix_current(
            current,
            index,
            size - 2 - index,
        )

    algebra_audits = 0
    for index in range(1, size - 1):
        shifted = (
            shoulder_prefix[index - 1]
            if index >= 2
            else 0.0
        )
        expected = (
            shifted
            + current[index - 1] * current_suffix[index + 1]
        )
        if not close(boundary_prefix[index], expected, 2e-11):
            issues.append(
                f"shifted prefix mismatch K={size}, R={rank}, i={index}"
            )
        algebra_audits += 1

    boundary = sum(
        boundary_prefix[index] * w[index]
        for index in range(1, size - 1)
    )
    shoulder = sum(
        shoulder_prefix[index] * p[index + 1]
        for index in range(1, size - 2)
    )
    main_suffix = sum(
        w[index]
        * current[index - 1]
        * current_suffix[index + 1]
        for index in range(1, size - 1)
    )
    shoulder_delta = sum(
        delta[index] * shoulder_prefix[index - 1]
        for index in range(2, size - 1)
    )
    if not close(
        boundary - shoulder,
        main_suffix + shoulder_delta,
        3e-10,
    ):
        issues.append(
            f"boundary split mismatch K={size}, R={rank}"
        )
    algebra_audits += 1

    edge_audits = 0
    for base in range(1, size - 1):
        for gap in range(1, size - base):
            occurrence = sum(
                w[index]
                for index in range(base, size - gap)
            ) - sum(
                p[index + 1]
                for index in range(base, size - 1 - gap)
            )
            reduced = w[base] + sum(
                delta[index]
                for index in range(base + 1, size - gap)
            )
            spectral = 0.0
            for mode in range(1, rank + 1):
                endpoint = phi(size, mode, size)
                penultimate = phi(size, mode, size - 1)
                tail = sum(
                    phi(size, mode, index)
                    for index in range(base + 1, size - gap)
                )
                spectral += (
                    penultimate * phi(size, mode, base)
                    + (penultimate - endpoint) * tail
                ) / q_value(mode) ** 2
            if not close(occurrence, reduced, 3e-10):
                issues.append(
                    f"edge reindex mismatch K={size}, R={rank}, "
                    f"a={base}, g={gap}"
                )
            if not close(reduced, spectral, 3e-10):
                issues.append(
                    f"edge spectrum mismatch K={size}, R={rank}, "
                    f"a={base}, g={gap}"
                )
            edge_audits += 2

    bound_audits = 0
    for index in range(1, size):
        exact_delta = sum(
            (
                -math.pi**2
                / normalizer**2
                * phi(size, mode, index)
                * phi(size, mode, size)
                * sinc(theta(size, mode) / 2.0) ** 2
            )
            for mode in range(1, rank + 1)
        )
        if not close(delta[index], exact_delta, 2e-11):
            issues.append(
                f"endpoint delta mismatch K={size}, R={rank}, i={index}"
            )
        point_bound = (
            4.0 * math.pi**2 * rank / normalizer**3
        )
        if abs(delta[index]) > point_bound + 2e-12:
            issues.append(
                f"endpoint point bound failed K={size}, R={rank}, "
                f"i={index}"
            )
        bound_audits += 2

    interval_bound = (
        8.0 * math.pi * odd_harmonic / normalizer**2
    )
    for left in range(1, size):
        running = 0.0
        for right in range(left, size):
            running += delta[right]
            if abs(running) > interval_bound + 3e-11:
                issues.append(
                    f"delta interval bound failed K={size}, R={rank}, "
                    f"I=[{left},{right}]"
                )
            bound_audits += 1

    extended_w = w + [c_value]
    for index in range(1, size + 1):
        increment = extended_w[index] - extended_w[index - 1]
        if increment <= -2e-12:
            issues.append(
                f"nonpositive weight increment K={size}, R={rank}, "
                f"i={index}"
            )
        if index < size:
            upper = (
                4.0 * math.pi * odd_harmonic / normalizer**2
            )
        else:
            upper = 4.0 * math.pi**2 / normalizer**2
        if increment > upper + 3e-11:
            issues.append(
                f"weight increment bound failed K={size}, R={rank}, "
                f"i={index}"
            )
        uniform_upper = (
            4.0
            * math.pi
            * (0.5 + 2.0 * math.pi)
            / normalizer**2
        )
        if increment > uniform_upper + 3e-11:
            issues.append(
                f"uniform increment bound failed K={size}, R={rank}, "
                f"i={index}"
            )
        bound_audits += 2

    energy = sum(
        (w[index] - w[index - 1])
        * (current_suffix[index] + beta) ** 2
        for index in range(1, size)
    ) + (c_value - w[size - 1]) * beta**2
    current_diagonal = sum(
        w[index] * current[index - 1] ** 2
        for index in range(1, size)
    )
    diagonal = current_diagonal + c_value * q_future
    correction = shoulder_delta - beta * sum(
        current[index - 1] * delta[index]
        for index in range(1, size)
    )
    min_gram_left = main_suffix + beta * sum(
        w[index] * current[index - 1]
        for index in range(1, size)
    )
    min_gram_right = 0.5 * (
        energy - current_diagonal - c_value * beta**2
    )
    if not close(min_gram_left, min_gram_right, 4e-10):
        issues.append(
            f"min-kernel Gram mismatch K={size}, R={rank}"
        )
    algebra_audits += 1

    remainder = (
        boundary
        - shoulder
        + beta * gamma
        + 0.5 * c_value * (beta**2 - q_future)
    )
    reduced_remainder = 0.5 * (energy - diagonal) + correction
    if not close(remainder, reduced_remainder, 4e-10):
        issues.append(
            f"joined remainder mismatch K={size}, R={rank}"
        )
    algebra_audits += 1

    if abs(shoulder_delta) > (
        math.pi * odd_harmonic + 3e-10
    ):
        issues.append(
            f"shoulder correction bound failed K={size}, R={rank}"
        )
    anchor_delta = beta * sum(
        current[index - 1] * delta[index]
        for index in range(1, size)
    )
    if abs(anchor_delta) > math.pi**2 / 2.0 + 3e-10:
        issues.append(
            f"anchor delta bound failed K={size}, R={rank}"
        )
    if current_diagonal > math.pi**2 / 4.0 + 3e-10:
        issues.append(
            f"current diagonal bound failed K={size}, R={rank}"
        )
    if c_value * q_future > math.pi**2 / 6.0 + 3e-10:
        issues.append(
            f"future diagonal bound failed K={size}, R={rank}"
        )
    if diagonal > 5.0 * math.pi**2 / 12.0 + 3e-10:
        issues.append(
            f"complete diagonal bound failed K={size}, R={rank}"
        )
    error_bound = (
        math.pi * odd_harmonic + 17.0 * math.pi**2 / 24.0
    )
    if abs(remainder - 0.5 * energy) > error_bound + 3e-10:
        issues.append(
            f"logarithmic approximation failed K={size}, R={rank}"
        )
    affine_tent_energy = (
        sum(
            (current_suffix[index] + beta) ** 2
            for index in range(1, size)
        )
        + beta**2
    )
    energy_upper = (
        4.0
        * math.pi
        * (0.5 + 2.0 * math.pi)
        / normalizer**2
        * affine_tent_energy
    )
    if energy > energy_upper + 4e-10:
        issues.append(
            f"block mean-square calibration failed K={size}, R={rank}"
        )
    if abs(beta) > normalizer / 2.0 + 3e-10:
        issues.append(
            f"future scalar bound failed K={size}, R={rank}"
        )
    bound_audits += 8
    return algebra_audits, edge_audits, bound_audits


def check_finite_algebra(issues: list[str]) -> tuple[int, int, int]:
    algebra_audits = 0
    edge_audits = 0
    bound_audits = 0
    mu = mobius_values(100)
    for alpha, size in ((0.21, 4), (0.53, 6), (0.83, 9)):
        for rank in sorted({1, max(1, size // 2), size}):
            for values in deterministic_vectors(size, mu):
                algebra, edges, bounds = check_case(
                    alpha,
                    size,
                    rank,
                    values,
                    issues,
                )
                algebra_audits += algebra
                edge_audits += edges
                bound_audits += bounds
    return algebra_audits, edge_audits, bound_audits


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
    w = [0.0] + [
        sum(
            mode_phi[index - 1][mode - 1]
            * penultimate[mode - 1]
            / q_value(mode) ** 2
            for mode in modes
        )
        for index in range(1, size)
    ]
    p = [0.0] + [
        sum(
            mode_phi[index - 1][mode - 1]
            * endpoint[mode - 1]
            / q_value(mode) ** 2
            for mode in modes
        )
        for index in range(1, size)
    ]
    c_value = sum(
        endpoint[mode - 1] ** 2 / q_value(mode) ** 2
        for mode in modes
    )
    delta = [0.0] + [
        w[index] - p[index]
        for index in range(1, size)
    ]
    current_suffix = suffixes(current)
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

    boundary = sum(
        fast_prefix(index, size - 1 - index) * w[index]
        for index in range(1, size - 1)
    )
    shoulder = sum(
        fast_prefix(index, size - 2 - index) * p[index + 1]
        for index in range(1, size - 2)
    )
    shoulder_delta = sum(
        delta[index]
        * fast_prefix(index - 1, size - 1 - index)
        for index in range(2, size - 1)
    )
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
        current[index - 1] * p[index]
        for index in range(1, size)
    )
    energy = sum(
        (w[index] - w[index - 1])
        * (current_suffix[index] + beta) ** 2
        for index in range(1, size)
    ) + (c_value - w[size - 1]) * beta**2
    diagonal = sum(
        w[index] * current[index - 1] ** 2
        for index in range(1, size)
    ) + c_value * q_future
    correction = shoulder_delta - beta * sum(
        current[index - 1] * delta[index]
        for index in range(1, size)
    )
    remainder = (
        boundary
        - shoulder
        + beta * gamma
        + 0.5 * c_value * (beta**2 - q_future)
    )
    return {
        "K": size,
        "R": rank,
        "E_B": energy,
        "D_B": diagonal,
        "C_delta": correction,
        "R_B": remainder,
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
        for key in ("R", "E_B", "D_B", "C_delta", "R_B"):
            if not close(
                float(stored[key]),
                float(computed[key]),
                3e-8,
            ):
                issues.append(
                    f"Mobius audit mismatch K={stored['K']}, key={key}"
                )
            audits += 1
        reconstructed = (
            0.5
            * (float(stored["E_B"]) - float(stored["D_B"]))
            + float(stored["C_delta"])
        )
        if not close(float(stored["R_B"]), reconstructed, 3e-8):
            issues.append(
                f"stored suffix identity mismatch K={stored['K']}"
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
            print(f"PBSEH [error] {issue}", file=sys.stderr)
        return 1
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        note_text = NOTE.read_text(encoding="utf-8")
        formal_text = FORMAL_CORE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as error:
        print(f"PBSEH [load-failure] {error}", file=sys.stderr)
        return 1

    check_schema(payload, note_text, formal_text, issues)
    algebra_audits, edge_audits, bound_audits = check_finite_algebra(
        issues
    )
    mobius_audits = check_mobius_table(payload, issues)
    check_reproduction(issues)

    if issues:
        for issue in issues:
            print(f"PBSEH [error] {issue}", file=sys.stderr)
        return 1
    print(
        "validated planar boundary suffix-energy handoff: "
        f"{payload['summary']['row_count']} rows, 0 issues, "
        f"{algebra_audits} algebra audits, "
        f"{edge_audits} edge audits, "
        f"{bound_audits} bound audits, "
        f"{mobius_audits} Mobius audits, "
        f"{payload['summary']['open_energy_gate_count']} "
        "open anchored energy gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
