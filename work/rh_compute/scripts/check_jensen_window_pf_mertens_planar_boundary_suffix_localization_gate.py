#!/usr/bin/env python3
"""Validate the planar boundary suffix-localization route guard."""

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
    "jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.py"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"


def close(left: float, right: float, tolerance: float = 3e-10) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def odd(mode: int) -> int:
    return 2 * mode - 1


def sine_sum(rank: int, value: float) -> float:
    return sum(
        math.sin(odd(mode) * value) / odd(mode)
        for mode in range(1, rank + 1)
    )


def sine_integral(rank: int, left: float, right: float) -> float:
    return sum(
        (
            math.cos(odd(mode) * left)
            - math.cos(odd(mode) * right)
        )
        / odd(mode) ** 2
        for mode in range(1, rank + 1)
    )


def spectral_increments(size: int, rank: int) -> list[float]:
    normalizer = 2 * size + 1
    result = [0.0] * (size + 1)
    for mode in range(1, rank + 1):
        frequency = odd(mode)
        angle = frequency * math.pi / normalizer
        scale = 4.0 / (normalizer * frequency**2)
        for index in range(1, size):
            result[index] += (
                scale
                * (
                    math.sin(index * angle)
                    - math.sin((index - 1) * angle)
                )
                * math.sin((size - 1) * angle)
            )
        result[size] += scale * (
            math.sin(size * angle) ** 2
            - math.sin((size - 1) * angle) ** 2
        )
    return result


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


def future_weight(alpha: float, size: int, value: int) -> float:
    return (
        (2.0 * size / value) ** (1.0 + alpha)
        * (4.0 * size - value)
        / (2.0 * size)
    )


def check_schema(
    payload: dict,
    note_text: str,
    formal_text: str,
    issues: list[str],
) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_planar_"
        "boundary_suffix_localization_gate"
    ):
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    rows = payload.get("rows", [])
    summary = payload.get("summary", {})
    if len(rows) != 26 or summary.get("row_count") != 26:
        issues.append("expected exactly 26 ledger rows")
    row_ids = [item.get("id") for item in rows]
    if len(row_ids) != len(set(row_ids)):
        issues.append("duplicate row id")
    for index in range(1, 27):
        prefix = f"pbslg_{index:02d}_"
        if not any(str(row_id).startswith(prefix) for row_id in row_ids):
            issues.append(f"missing row prefix {prefix}")
    for key, expected in (
        ("exact_result_count", 22),
        ("finite_validation_count", 1),
        ("guard_count", 3),
        ("open_energy_gate_count", 1),
    ):
        if summary.get(key) != expected:
            issues.append(f"{key} mismatch")
    for key in (
        "bulk_localization_proved",
        "terminal_collar_recovery_proved",
        "cofinal_rh_equivalence_proved",
    ):
        if summary.get(key) is not True:
            issues.append(f"{key} must be true")
    for key in (
        "uniform_lower_norm_comparison_proved",
        "anchored_mobius_energy_estimate_proved",
        "interior_projection_estimate_proved",
        "rh_proved",
        "pf_infinity_proved",
        "lambda_le_zero_proved",
    ):
        if summary.get(key) is not False:
            issues.append(f"{key} must be false")
    roles = {
        item.get("id"): (item.get("role"), item.get("status"))
        for item in rows
    }
    if roles.get("pbslg_25_open_energy_gate") != (
        "theorem_target",
        "open",
    ):
        issues.append("energy gate must remain open")
    if roles.get("pbslg_22_no_uniform_lower_comparison") != (
        "countermodel_guard",
        "guard_active",
    ):
        issues.append("countermodel guard status mismatch")
    for source in payload.get("source", {}).values():
        if not (REPO_ROOT / source).exists():
            issues.append(f"missing source {source}")
    for marker in (
        "Status: exact bulk-localization",
        "(PBSL.3)",
        "J:=min(K,K+2-L)",
        "(PBSL.8)",
        "(PBSL.11)",
        "(PBSL.15)",
        "Strict Norm-Separation Countermodel",
        "(PBSL.19)",
        "not an easier independent",
        "No `E_B` estimate",
    ):
        if marker not in note_text:
            issues.append(f"note missing marker {marker!r}")
    for marker in (
        "Corollary 11.22Z.10",
        "(11.22Z.140)",
        "(11.22Z.141)",
        "(11.22Z.142)",
        "(11.22Z.143)",
        "(11.22Z.144)",
        "(11.22Z.145)",
        "(11.22Z.146)",
        "(11.22Z.147)",
        "(11.22Z.148)",
        "jensen_window_pf_mertens_planar_boundary_suffix_localization_gate.md",
    ):
        if marker not in formal_text:
            issues.append(f"formal core missing marker {marker!r}")


def check_odd_sine_algebra(issues: list[str]) -> int:
    audits = 0
    for denominator in (37, 83, 191):
        for numerator in range(1, denominator):
            value = math.pi * numerator / denominator
            for first, last in ((1, 1), (1, 7), (3, 11), (9, 23)):
                direct = sum(
                    math.sin(odd(mode) * value)
                    for mode in range(first, last + 1)
                )
                formula = (
                    math.sin((last - first + 1) * value)
                    * math.sin((first + last - 1) * value)
                    / math.sin(value)
                )
                if not close(direct, formula, 2e-11):
                    issues.append(
                        f"odd-sine partial sum mismatch u={value}, "
                        f"I=[{first},{last}]"
                    )
                audits += 1
    for rank in range(1, 41):
        for index in range(1, 1000):
            value = math.pi * index / 1000.0
            error = abs(sine_sum(rank, value) - math.pi / 4.0)
            bound = 1.0 / ((2 * rank + 1) * math.sin(value))
            if error > bound + 2e-12:
                issues.append(
                    f"square-wave tail bound failed R={rank}, u={value}"
                )
            audits += 1
    return audits


def check_slivers_and_bulk(issues: list[str]) -> tuple[int, int]:
    identity_audits = 0
    bound_audits = 0
    for size in range(4, 65):
        normalizer = 2 * size + 1
        angle_unit = math.pi / normalizer
        ranks = sorted(
            {
                1,
                max(1, size // 4),
                max(1, size // 2),
                size,
                math.ceil(size ** 0.875),
            }
        )
        for rank in ranks:
            delta = spectral_increments(size, rank)
            for index in range(1, size):
                integral = 2.0 / normalizer * (
                    sine_integral(
                        rank,
                        (size - 1 - index) * angle_unit,
                        (size - index) * angle_unit,
                    )
                    + sine_integral(
                        rank,
                        (size - 2 + index) * angle_unit,
                        (size - 1 + index) * angle_unit,
                    )
                )
                if not close(delta[index], integral, 3e-11):
                    issues.append(
                        f"current sliver mismatch K={size}, R={rank}, "
                        f"i={index}"
                    )
                identity_audits += 1
            endpoint_integral = 2.0 / normalizer * sine_integral(
                rank,
                angle_unit,
                3.0 * angle_unit,
            )
            if not close(delta[size], endpoint_integral, 3e-11):
                issues.append(
                    f"endpoint sliver mismatch K={size}, R={rank}"
                )
            identity_audits += 1
            eta = math.asin(8.0 / (math.pi * (2 * rank + 1)))
            collar = math.ceil(normalizer * eta / math.pi - 1e-14)
            bulk_end = min(size, size + 2 - collar)
            if bulk_end < 1:
                issues.append(
                    f"empty bulk K={size}, R={rank}, L={collar}"
                )
                continue
            for index in range(1, bulk_end + 1):
                if delta[index] + 3e-12 < (
                    math.pi**2 / (4.0 * normalizer**2)
                ):
                    issues.append(
                        f"bulk lower bound failed K={size}, R={rank}, "
                        f"L={collar}, J={bulk_end}, i={index}"
                    )
                bound_audits += 1
            if collar > (
                1.0 + 2.0 * normalizer / (math.pi * rank) + 1e-12
            ):
                issues.append(
                    f"collar length bound failed K={size}, R={rank}"
                )
            bound_audits += 1
    return identity_audits, bound_audits


def deterministic_currents(size: int, mu: list[int]) -> list[list[float]]:
    return [
        [float(mu[size + index]) for index in range(1, size)],
        [
            float(((7 * index + 2 * size) % 3) - 1)
            for index in range(1, size)
        ],
        [(-1.0 if index == size - 1 else 0.0) for index in range(1, size)],
    ]


def check_collar_recovery(issues: list[str]) -> int:
    audits = 0
    mu = mobius_values(400)
    for size in (4, 7, 12, 23, 48, 79):
        normalizer = 2 * size + 1
        for rank in sorted({1, max(1, size // 3), size}):
            delta = spectral_increments(size, rank)
            eta = math.asin(8.0 / (math.pi * (2 * rank + 1)))
            collar = math.ceil(normalizer * eta / math.pi - 1e-14)
            bulk_end = min(size, size + 2 - collar)
            missing = size - bulk_end
            for current in deterministic_currents(size, mu):
                for beta in (-1.75, 0.0, 2.25):
                    suffix = 0.0
                    b_values = [0.0] * (size + 1)
                    b_values[size] = beta
                    for index in range(size - 1, 0, -1):
                        suffix += current[index - 1]
                        b_values[index] = suffix + beta
                    energy = sum(
                        delta[index] * b_values[index] ** 2
                        for index in range(1, size + 1)
                    )
                    total = sum(
                        b_values[index] ** 2
                        for index in range(1, size + 1)
                    )
                    bulk = sum(
                        b_values[index] ** 2
                        for index in range(1, bulk_end + 1)
                    )
                    exact_upper = (
                        (1 + 2 * missing) * bulk
                        + missing
                        * (missing + 1)
                        * (2 * missing + 1)
                        / 3.0
                    )
                    coarse_upper = (
                        4.0
                        * normalizer**2
                        * (2 * collar + 1)
                        / math.pi**2
                        * energy
                        + 2.0 * collar**3
                    )
                    if total > exact_upper + 2e-9:
                        issues.append(
                            f"exact collar recovery failed K={size}, "
                            f"R={rank}, beta={beta}"
                        )
                    if total > coarse_upper + 2e-8:
                        issues.append(
                            f"coarse collar recovery failed K={size}, "
                            f"R={rank}, beta={beta}"
                        )
                    audits += 2
    return audits


def recompute_diagnostic(alpha: float, size: int, mu: list[int]) -> dict:
    rank = math.ceil(size ** (1.0 - alpha / 2.0))
    normalizer = 2 * size + 1
    delta = spectral_increments(size, rank)
    eta = math.asin(8.0 / (math.pi * (2 * rank + 1)))
    collar = math.ceil(normalizer * eta / math.pi - 1e-14)
    bulk_end = min(size, size + 2 - collar)
    beta = sum(
        mu[value] * future_weight(alpha, size, value)
        for value in range(2 * size, 4 * size)
    )
    suffix = 0.0
    values = [0.0] * (size + 1)
    values[size] = beta
    for index in range(size - 1, 0, -1):
        suffix += mu[size + index]
        values[index] = suffix + beta
    energy = sum(
        delta[index] * values[index] ** 2
        for index in range(1, size + 1)
    )
    affine = sum(
        values[index] ** 2 for index in range(1, size + 1)
    )
    return {
        "K": size,
        "R": rank,
        "L": collar,
        "J": bulk_end,
        "min_bulk_N2_Delta": min(
            normalizer**2 * delta[index]
            for index in range(1, bulk_end + 1)
        ),
        "terminal_K2_Delta": size**2 * delta[size],
        "mobius_K2_E_over_H": size**2 * energy / affine,
    }


def check_countermodel_and_table(payload: dict, issues: list[str]) -> int:
    audits = 0
    mu = mobius_values(4098)
    stored_rows = payload.get("diagnostic", {}).get("rows", [])
    expected_sizes = [16, 32, 64, 128, 256, 512, 1024]
    if [item.get("K") for item in stored_rows] != expected_sizes:
        issues.append("diagnostic K grid mismatch")
        return audits
    for stored in stored_rows:
        size = int(stored["K"])
        computed = recompute_diagnostic(0.5, size, mu)
        for key in (
            "R",
            "L",
            "J",
            "min_bulk_N2_Delta",
            "terminal_K2_Delta",
            "mobius_K2_E_over_H",
        ):
            if not close(float(stored[key]), float(computed[key]), 3e-8):
                issues.append(
                    f"diagnostic mismatch K={size}, key={key}"
                )
            audits += 1
        rank = int(computed["R"])
        normalizer = 2 * size + 1
        endpoint = spectral_increments(size, rank)[size]
        endpoint_bound = 8.0 * math.pi**2 * rank / normalizer**3
        if endpoint <= 0.0 or endpoint > endpoint_bound + 2e-12:
            issues.append(f"terminal countermodel bound failed K={size}")
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
            print(f"PBSLG [error] {issue}", file=sys.stderr)
        return 1
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        note_text = NOTE.read_text(encoding="utf-8")
        formal_text = FORMAL_CORE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as error:
        print(f"PBSLG [load-failure] {error}", file=sys.stderr)
        return 1
    check_schema(payload, note_text, formal_text, issues)
    sine_audits = check_odd_sine_algebra(issues)
    identity_audits, bulk_audits = check_slivers_and_bulk(issues)
    collar_audits = check_collar_recovery(issues)
    diagnostic_audits = check_countermodel_and_table(payload, issues)
    check_reproduction(issues)
    if issues:
        for issue in issues:
            print(f"PBSLG [error] {issue}", file=sys.stderr)
        return 1
    print(
        "validated planar boundary suffix-localization gate: "
        f"{payload['summary']['row_count']} rows, 0 issues, "
        f"{sine_audits} odd-sine audits, "
        f"{identity_audits} sliver audits, "
        f"{bulk_audits} bulk audits, "
        f"{collar_audits} collar audits, "
        f"{diagnostic_audits} diagnostic/countermodel audits, "
        f"{payload['summary']['open_energy_gate_count']} "
        "open RH-equivalent energy gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
