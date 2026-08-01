#!/usr/bin/env python3
"""Validate the planar joined-energy equivalence route guard."""

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
    "jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.py"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"


def close(left: float, right: float, tolerance: float = 4e-9) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


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


def prefix_current(
    current: list[float],
    prefix: list[float],
    endpoint: int,
    width: int,
) -> float:
    return sum(
        current[base - 1]
        * (prefix[base + width] - prefix[base])
        for base in range(1, endpoint + 1)
    )


def compute_case(
    alpha: float,
    size: int,
    rank: int,
    current: list[float],
    future: list[float],
) -> dict[str, float]:
    normalizer = 2 * size + 1
    frequencies = [2 * mode - 1 for mode in range(1, rank + 1)]
    mode_phi = [
        [
            2.0
            / math.sqrt(normalizer)
            * math.sin(index * frequency * math.pi / normalizer)
            for frequency in frequencies
        ]
        for index in range(size + 1)
    ]
    endpoint = mode_phi[size]
    penultimate = mode_phi[size - 1]
    w = [0.0] + [
        sum(
            mode_phi[index][mode]
            * penultimate[mode]
            / frequencies[mode] ** 2
            for mode in range(rank)
        )
        for index in range(1, size)
    ]
    p = [0.0] + [
        sum(
            mode_phi[index][mode]
            * endpoint[mode]
            / frequencies[mode] ** 2
            for mode in range(rank)
        )
        for index in range(1, size + 1)
    ]
    c_value = p[size]
    delta = [0.0] + [
        w[index] - p[index]
        for index in range(1, size)
    ]
    beta = sum(
        future[offset] * future_weight(alpha, size, 2 * size + offset)
        for offset in range(2 * size)
    )
    q_future = sum(
        future[offset] ** 2
        * future_weight(alpha, size, 2 * size + offset) ** 2
        for offset in range(2 * size)
    )
    x_modes = [
        sum(
            current[index - 1] * mode_phi[index][mode]
            for index in range(1, size)
        )
        for mode in range(rank)
    ]
    current_mode_diagonal = [
        sum(
            current[index - 1] ** 2 * mode_phi[index][mode] ** 2
            for index in range(1, size)
        )
        for mode in range(rank)
    ]
    lossless = sum(
        (x_modes[mode] + beta * endpoint[mode]) ** 2
        / frequencies[mode] ** 2
        for mode in range(rank)
    )
    lossless_diagonal = (
        sum(
            current_mode_diagonal[mode] / frequencies[mode] ** 2
            for mode in range(rank)
        )
        + c_value * q_future
    )
    total_offdiagonal = 0.5 * (lossless - lossless_diagonal)
    current_offdiagonal = 0.5 * sum(
        (
            x_modes[mode] ** 2
            - current_mode_diagonal[mode]
        )
        / frequencies[mode] ** 2
        for mode in range(rank)
    )
    prefix = [0.0]
    for value in current:
        prefix.append(prefix[-1] + value)
    boundary = sum(
        prefix_current(
            current,
            prefix,
            index,
            size - 1 - index,
        )
        * w[index]
        for index in range(1, size - 1)
    )
    shoulder = sum(
        prefix_current(
            current,
            prefix,
            index,
            size - 2 - index,
        )
        * p[index + 1]
        for index in range(1, size - 2)
    )
    current_interior = current_offdiagonal - boundary + shoulder
    suffix = [0.0] * (size + 1)
    for index in range(size - 1, 0, -1):
        suffix[index] = suffix[index + 1] + current[index - 1]
    boundary_energy = sum(
        (w[index] - w[index - 1])
        * (suffix[index] + beta) ** 2
        for index in range(1, size)
    ) + (c_value - w[size - 1]) * beta**2
    boundary_diagonal = sum(
        w[index] * current[index - 1] ** 2
        for index in range(1, size)
    ) + c_value * q_future
    shoulder_delta = sum(
        delta[index]
        * prefix_current(
            current,
            prefix,
            index - 1,
            size - 1 - index,
        )
        for index in range(2, size - 1)
    )
    correction = shoulder_delta - beta * sum(
        current[index - 1] * delta[index]
        for index in range(1, size)
    )
    remainder = (
        0.5 * (boundary_energy - boundary_diagonal)
        + correction
    )
    joined = current_interior + 0.5 * boundary_energy
    return {
        "L": lossless,
        "D_L": lossless_diagonal,
        "O": total_offdiagonal,
        "O_direct": (
            current_offdiagonal
            + beta
            * sum(
                current[index - 1] * p[index]
                for index in range(1, size)
            )
            + 0.5 * c_value * (beta**2 - q_future)
        ),
        "E_B": boundary_energy,
        "D_B": boundary_diagonal,
        "C_delta": correction,
        "R_B": remainder,
        "O_CC_int": current_interior,
        "J_B": joined,
        "J_rhs": (
            0.5 * lossless
            + 0.5 * (boundary_diagonal - lossless_diagonal)
            - correction
        ),
    }


def check_schema(
    payload: dict,
    note_text: str,
    formal_text: str,
    issues: list[str],
) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_planar_"
        "joined_energy_equivalence_gate"
    ):
        issues.append("result kind mismatch")
    rows = payload.get("rows", [])
    summary = payload.get("summary", {})
    if len(rows) != 20 or summary.get("row_count") != 20:
        issues.append("expected exactly 20 rows")
    row_ids = [item.get("id") for item in rows]
    if len(row_ids) != len(set(row_ids)):
        issues.append("duplicate row id")
    for index in range(1, 21):
        prefix = f"pjeg_{index:02d}_"
        if not any(str(row_id).startswith(prefix) for row_id in row_ids):
            issues.append(f"missing row prefix {prefix}")
    for key, expected in (
        ("exact_result_count", 16),
        ("finite_validation_count", 1),
        ("guard_count", 3),
        ("open_energy_gate_count", 1),
    ):
        if summary.get(key) != expected:
            issues.append(f"{key} mismatch")
    for key in (
        "joined_identity_proved",
        "weighted_dyadic_equivalence_proved",
    ):
        if summary.get(key) is not True:
            issues.append(f"{key} must be true")
    for key in (
        "joined_route_distinct",
        "lossless_energy_estimate_proved",
        "boundary_energy_estimate_proved",
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
    if roles.get("pjeg_18_open_lossless_gate") != (
        "theorem_target",
        "open",
    ):
        issues.append("lossless energy gate must remain open")
    for source in payload.get("source", {}).values():
        if not (REPO_ROOT / source).exists():
            issues.append(f"missing source {source}")
    for marker in (
        "Status: exact equivalence",
        "(PJE.2)",
        "J_B:=O_(CC,int)+E_B/2",
        "(PJE.6)",
        "41pi^2/48",
        "(PJE.9)",
        "does not create a new weaker",
        "No estimate for `L`",
    ):
        if marker not in note_text:
            issues.append(f"note missing marker {marker!r}")
    for marker in (
        "Corollary 11.22Z.11",
        "(11.22Z.149)",
        "(11.22Z.150)",
        "(11.22Z.151)",
        "(11.22Z.152)",
        "(11.22Z.153)",
        "(11.22Z.154)",
        "jensen_window_pf_mertens_planar_joined_energy_equivalence_gate.md",
    ):
        if marker not in formal_text:
            issues.append(f"formal core missing marker {marker!r}")


def deterministic_vectors(
    size: int,
    mu: list[int],
) -> list[tuple[list[float], list[float]]]:
    mobius_current = [
        float(mu[size + index]) for index in range(1, size)
    ]
    mobius_future = [
        float(mu[value]) for value in range(2 * size, 4 * size)
    ]
    patterned_current = [
        float(((5 * index + size) % 3) - 1)
        for index in range(1, size)
    ]
    patterned_future = [
        float(((7 * index + size) % 3) - 1)
        for index in range(2 * size)
    ]
    terminal_current = [
        -1.0 if index == size - 1 else 0.0
        for index in range(1, size)
    ]
    terminal_future = [1.0] + [0.0] * (2 * size - 1)
    return [
        (mobius_current, mobius_future),
        (patterned_current, patterned_future),
        (terminal_current, terminal_future),
    ]


def check_arbitrary_algebra(issues: list[str]) -> int:
    audits = 0
    mu = mobius_values(200)
    for alpha, size in ((0.21, 4), (0.47, 7), (0.83, 11)):
        for rank in sorted({1, max(1, size // 2), size}):
            odd_harmonic = sum(
                1.0 / (2 * mode - 1)
                for mode in range(1, rank + 1)
            )
            for current, future in deterministic_vectors(size, mu):
                values = compute_case(
                    alpha,
                    size,
                    rank,
                    current,
                    future,
                )
                for left, right, label in (
                    (values["O"], values["O_direct"], "block"),
                    (
                        values["O"],
                        values["O_CC_int"] + values["R_B"],
                        "remainder",
                    ),
                    (values["J_B"], values["J_rhs"], "joined"),
                ):
                    if not close(left, right, 5e-10):
                        issues.append(
                            f"{label} identity failed K={size}, R={rank}"
                        )
                    audits += 1
                if values["D_L"] < -2e-12 or values["D_L"] > (
                    7.0 * math.pi**2 / 24.0 + 2e-10
                ):
                    issues.append(
                        f"lossless diagonal bound failed K={size}, R={rank}"
                    )
                if values["D_B"] < -2e-12 or values["D_B"] > (
                    5.0 * math.pi**2 / 12.0 + 2e-10
                ):
                    issues.append(
                        f"boundary diagonal bound failed K={size}, R={rank}"
                    )
                correction_bound = (
                    math.pi * odd_harmonic + math.pi**2 / 2.0
                )
                if abs(values["C_delta"]) > correction_bound + 2e-10:
                    issues.append(
                        f"defect bound failed K={size}, R={rank}"
                    )
                joined_bound = (
                    math.pi * odd_harmonic
                    + 41.0 * math.pi**2 / 48.0
                )
                if abs(values["J_B"] - values["L"] / 2.0) > (
                    joined_bound + 2e-10
                ):
                    issues.append(
                        f"joined error bound failed K={size}, R={rank}"
                    )
                audits += 4
    return audits


def check_diagnostics(payload: dict, issues: list[str]) -> int:
    audits = 0
    mu = mobius_values(4098)
    rows = payload.get("diagnostic", {}).get("rows", [])
    expected_sizes = [16, 32, 64, 128, 256, 512, 1024]
    if [item.get("K") for item in rows] != expected_sizes:
        issues.append("diagnostic K grid mismatch")
        return audits
    for stored in rows:
        size = int(stored["K"])
        rank = math.ceil(size ** 0.75)
        current = [
            float(mu[size + index]) for index in range(1, size)
        ]
        future = [
            float(mu[value]) for value in range(2 * size, 4 * size)
        ]
        computed = compute_case(
            0.5,
            size,
            rank,
            current,
            future,
        )
        comparisons = {
            "R": rank,
            "L": computed["L"],
            "D_L": computed["D_L"],
            "O": computed["O"],
            "E_B": computed["E_B"],
            "O_CC_int": computed["O_CC_int"],
            "J_B": computed["J_B"],
            "J_minus_L_over_2": (
                computed["J_B"] - computed["L"] / 2.0
            ),
        }
        for key, value in comparisons.items():
            if not close(float(stored[key]), float(value), 4e-8):
                issues.append(
                    f"diagnostic mismatch K={size}, key={key}"
                )
            audits += 1
        if not close(computed["J_B"], computed["J_rhs"], 5e-10):
            issues.append(f"joined diagnostic identity failed K={size}")
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
            print(f"PJEG [error] {issue}", file=sys.stderr)
        return 1
    try:
        payload = json.loads(RESULT.read_text(encoding="utf-8"))
        note_text = NOTE.read_text(encoding="utf-8")
        formal_text = FORMAL_CORE.read_text(encoding="utf-8")
    except (OSError, json.JSONDecodeError) as error:
        print(f"PJEG [load-failure] {error}", file=sys.stderr)
        return 1
    check_schema(payload, note_text, formal_text, issues)
    algebra_audits = check_arbitrary_algebra(issues)
    diagnostic_audits = check_diagnostics(payload, issues)
    check_reproduction(issues)
    if issues:
        for issue in issues:
            print(f"PJEG [error] {issue}", file=sys.stderr)
        return 1
    print(
        "validated planar joined-energy equivalence gate: "
        f"{payload['summary']['row_count']} rows, 0 issues, "
        f"{algebra_audits} arbitrary-vector audits, "
        f"{diagnostic_audits} Mobius audits, "
        f"{payload['summary']['open_energy_gate_count']} "
        "open lossless energy gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
