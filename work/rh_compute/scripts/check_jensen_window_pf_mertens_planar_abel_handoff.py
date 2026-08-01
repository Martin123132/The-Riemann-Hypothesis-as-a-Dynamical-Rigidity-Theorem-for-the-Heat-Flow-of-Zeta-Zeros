#!/usr/bin/env python3
"""Validate the Mertens planar Abel and mixed-variation handoff."""

from __future__ import annotations

import importlib.util
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
    "jensen_window_pf_mertens_planar_abel_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_mertens_planar_abel_handoff.md"
)
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_planar_abel_handoff.py"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"


def close(
    left: float,
    right: float,
    tolerance: float = 2e-10,
) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def mobius_values(limit: int) -> list[int]:
    values = [1] * (limit + 1)
    prime = [True] * (limit + 1)
    values[0] = 0
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        for multiple in range(p, limit + 1, p):
            prime[multiple] = False
            values[multiple] *= -1
        square = p * p
        for multiple in range(square, limit + 1, square):
            values[multiple] = 0
    return values


def cutoff(alpha: float, size: int) -> int:
    return min(size, math.ceil(size ** (1.0 - alpha / 2.0)))


def angle(size: int, mode: int) -> float:
    return (2 * mode - 1) * math.pi / (2 * size + 1)


def eigenvector_value(size: int, mode: int, index: int) -> float:
    return (
        2.0
        / math.sqrt(2 * size + 1)
        * math.sin(index * angle(size, mode))
    )


def kernel(size: int, rank: int, left: int, right: int) -> float:
    return sum(
        eigenvector_value(size, mode, left)
        * eigenvector_value(size, mode, right)
        / (2 * mode - 1) ** 2
        for mode in range(1, rank + 1)
    )


def cosine_sum(rank: int, value: float) -> float:
    return sum(
        math.cos((2 * mode - 1) * value) / (2 * mode - 1) ** 2
        for mode in range(1, rank + 1)
    )


def dirichlet_sum(rank: int, value: float) -> float:
    return sum(
        math.cos((2 * mode - 1) * value)
        for mode in range(1, rank + 1)
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


def ordinary_position(
    alpha: float,
    size: int,
    value: int,
) -> tuple[int, float] | None:
    if not size < value < 4 * size:
        return None
    if value < 2 * size:
        return value - size, 1.0
    return size, future_weight(alpha, size, value)


def ordinary_kernel(
    alpha: float,
    size: int,
    rank: int,
    left: int,
    right: int,
) -> float:
    if not left < right:
        return 0.0
    left_position = ordinary_position(alpha, size, left)
    right_position = ordinary_position(alpha, size, right)
    if left_position is None or right_position is None:
        return 0.0
    left_index, left_weight = left_position
    right_index, right_weight = right_position
    return (
        left_weight
        * right_weight
        * kernel(size, rank, left_index, right_index)
    )


def extended_weight(
    alpha: float,
    size: int,
    rank: int,
    n: int,
    shift: int,
) -> float:
    if shift < 1 or n + shift > 4 * size - 1:
        return 0.0
    return ordinary_kernel(alpha, size, rank, n, n + shift)


def mixed_difference(
    alpha: float,
    size: int,
    rank: int,
    n: int,
    shift: int,
) -> float:
    return (
        extended_weight(alpha, size, rank, n, shift)
        - extended_weight(alpha, size, rank, n + 1, shift)
        - extended_weight(alpha, size, rank, n, shift + 1)
        + extended_weight(alpha, size, rank, n + 1, shift + 1)
    )


def planar_arrays(
    alpha: float,
    size: int,
) -> tuple[
    int,
    list[list[int]],
    list[list[float]],
    list[list[int]],
]:
    rank = cutoff(alpha, size)
    lower = size + 1
    side = 3 * size - 2
    mu = mobius_values(4 * size + 2)
    arithmetic = [[0 for _ in range(side)] for _ in range(side)]
    weights = [[0.0 for _ in range(side)] for _ in range(side)]
    prefixes = [[0 for _ in range(side)] for _ in range(side)]
    for i in range(side):
        n = lower + i
        for j in range(side):
            shift = j + 1
            if n + shift <= 4 * size - 1:
                arithmetic[i][j] = mu[n] * mu[n + shift]
            weights[i][j] = extended_weight(
                alpha,
                size,
                rank,
                n,
                shift,
            )
            prefixes[i][j] = (
                arithmetic[i][j]
                + (prefixes[i - 1][j] if i else 0)
                + (prefixes[i][j - 1] if j else 0)
                - (prefixes[i - 1][j - 1] if i and j else 0)
            )
    return rank, arithmetic, weights, prefixes


def check_schema(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    if payload.get("kind") != "jensen_window_pf_mertens_planar_abel_handoff":
        issues.append("result kind mismatch")
    if payload.get("date") != "2026-07-24":
        issues.append("result date mismatch")
    if payload.get("status") != (
        "exact planar Abel identity and logarithmic mixed-kernel "
        "variation with one open arithmetic energy gate"
    ):
        issues.append("result status mismatch")
    rows = payload.get("rows", [])
    if len(rows) != 36:
        issues.append("expected exactly 36 rows")
    ids = [entry.get("id", "") for entry in rows]
    if len(set(ids)) != len(ids):
        issues.append("row ids are not unique")
    for index in range(1, 37):
        prefix = f"pah_{index:02d}_"
        if not any(row_id.startswith(prefix) for row_id in ids):
            issues.append(f"missing row prefix {prefix}")
    exact_count = sum(
        entry.get("status") == "available_exact" for entry in rows
    )
    if exact_count != 28:
        issues.append(f"expected 28 exact rows, found {exact_count}")
    role_counts = {
        role: sum(entry.get("role") == role for entry in rows)
        for role in (
            "conditional_calibration",
            "literature_guard",
            "proof_guard",
            "countermodel_guard",
            "finite_validation",
            "open_gate",
        )
    }
    expected_role_counts = {
        "conditional_calibration": 2,
        "literature_guard": 1,
        "proof_guard": 2,
        "countermodel_guard": 1,
        "finite_validation": 1,
        "open_gate": 1,
    }
    if role_counts != expected_role_counts:
        issues.append(
            f"role counts mismatch: {role_counts} != {expected_role_counts}"
        )
    expected_audit = {
        "row_count": 36,
        "exact_reduction_count": 28,
        "conditional_calibration_count": 2,
        "literature_guard_count": 1,
        "proof_guard_count": 3,
        "finite_validation_count": 1,
        "open_planar_gate_count": 1,
        "double_abel_proved": True,
        "mixed_variation_bound_proved": True,
        "curvature_cauchy_proved": True,
        "axiswise_countermodel_proved": True,
        "curvature_energy_proved": False,
        "joint_mobius_gain_proved": False,
        "full_burnol_bound_proved": False,
        "rh_proved": False,
        "pf_infinity_proved": False,
        "lambda_le_zero_proved": False,
    }
    audit = payload.get("audit", {})
    for key, expected in expected_audit.items():
        if audit.get(key) != expected:
            issues.append(f"audit mismatch for {key}")
    for marker in (
        "(PAH.4)",
        "(PAH.17)",
        "(PAH.20)",
        "(PAH.21)",
        "(PAH.23)",
        "(PAH.25)",
        "false `O(1)`",
        "remain open",
    ):
        if marker not in note_text:
            issues.append(f"note missing marker: {marker}")
    for source in (
        "https://doi.org/10.2140/ant.2015.9.2167",
        "https://doi.org/10.1017/fmp.2023.28",
        "https://doi.org/10.5802/aif.2401",
    ):
        if source not in payload.get("source_anchors", []):
            issues.append(f"missing primary source {source}")


def check_builder_reproduction(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    spec = importlib.util.spec_from_file_location("planar_abel_builder", BUILDER)
    if spec is None or spec.loader is None:
        issues.append("could not load builder")
        return
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.build_payload() != payload:
        issues.append("builder payload does not reproduce result")
    if module.render_note() != note_text:
        issues.append("builder note does not reproduce output")
    with tempfile.TemporaryDirectory() as temporary:
        result = Path(temporary) / "result.json"
        note = Path(temporary) / "note.md"
        completed = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--result",
                str(result),
                "--note",
                str(note),
            ],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        if completed.returncode != 0:
            issues.append(
                f"builder subprocess failed: {completed.stderr.strip()}"
            )
        elif (
            json.loads(result.read_text(encoding="utf-8")) != payload
            or note.read_text(encoding="utf-8") != note_text
        ):
            issues.append("builder subprocess output mismatch")


def check_dirichlet_identity(issues: list[str]) -> None:
    for rank in (1, 2, 3, 5, 8, 13):
        for numerator in (1, 2, 5, 9, 17):
            value = numerator * math.pi / 37.0
            direct = dirichlet_sum(rank, value)
            closed = math.sin(2 * rank * value) / (2 * math.sin(value))
            if not close(direct, closed, 2e-12):
                issues.append(
                    f"Dirichlet identity failed at R={rank}, u={value}"
                )
                return
        steps = 12000
        width = math.pi / steps
        integral = sum(
            abs(dirichlet_sum(rank, (index + 0.5) * width))
            for index in range(steps)
        ) * width
        bound = math.pi / 2.0 * (1.0 + math.log(2.0 * rank))
        if integral > bound + 2e-4:
            issues.append(f"Dirichlet L1 bound failed at R={rank}")
            return


def check_future_geometry(issues: list[str]) -> None:
    for alpha in (0.1, 0.37, 0.75, 0.95):
        for size in (3, 5, 8, 13, 21):
            values = [
                future_weight(alpha, size, value)
                for value in range(2 * size, 4 * size + 1)
            ]
            decrements = [
                values[index] - values[index + 1]
                for index in range(len(values) - 1)
            ] + [0.0]
            if any(value < -1e-13 for value in decrements):
                issues.append(
                    f"future decrement negative at alpha={alpha}, K={size}"
                )
                return
            if any(
                decrements[index] + 1e-12 < decrements[index + 1]
                for index in range(len(decrements) - 1)
            ):
                issues.append(
                    f"future decrements not decreasing at alpha={alpha}, "
                    f"K={size}"
                )
                return
            if not close(sum(decrements), 1.0, 2e-12):
                issues.append(
                    f"future decrement telescope failed at K={size}"
                )
                return
            initial_bound = (2.0 + alpha) / (2.0 * size)
            if decrements[0] > initial_bound + 1e-12:
                issues.append(
                    f"initial future decrement bound failed at K={size}"
                )
                return


def check_planar_case(
    alpha: float,
    size: int,
    issues: list[str],
) -> None:
    rank, arithmetic, weights, prefixes = planar_arrays(alpha, size)
    lower = size + 1
    side = 3 * size - 2
    rho = math.pi**2 * size / (2 * size + 1) ** 2
    direct = 0.0
    abel = 0.0
    variation = 0.0
    energy = 0.0
    maximum = 0
    pieces = {
        "cc": 0.0,
        "second_transition": 0.0,
        "cf": 0.0,
        "first_transition": 0.0,
        "ff": 0.0,
    }
    for i in range(side):
        n = lower + i
        for j in range(side):
            shift = j + 1
            m = n + shift
            delta = mixed_difference(alpha, size, rank, n, shift)
            direct += arithmetic[i][j] * weights[i][j]
            abel += prefixes[i][j] * delta
            variation += abs(delta)
            energy += abs(delta) * prefixes[i][j] ** 2
            maximum = max(maximum, abs(prefixes[i][j]))
            if m > 4 * size - 1:
                if abs(delta) > 2e-12:
                    issues.append(
                        f"nonzero curvature beyond support at "
                        f"alpha={alpha}, K={size}, n={n}, h={shift}"
                    )
                    return
                continue
            if n <= 2 * size - 2 and m <= 2 * size - 2:
                piece = "cc"
                local_i = n - size
                s = 2 * local_i + shift
                step = math.pi / (2 * size + 1)
                expected = 2.0 / (2 * size + 1) * (
                    -cosine_sum(rank, s * step)
                    + cosine_sum(rank, (s + 2) * step)
                    + cosine_sum(rank, (s + 1) * step)
                    - cosine_sum(rank, (s + 3) * step)
                )
            elif n <= 2 * size - 2 and m == 2 * size - 1:
                piece = "second_transition"
                local_i = n - size
                p_i = kernel(size, rank, local_i, size)
                p_next = kernel(size, rank, local_i + 1, size)
                decrement = (
                    future_weight(alpha, size, 2 * size)
                    - future_weight(alpha, size, 2 * size + 1)
                )
                expected = (
                    kernel(size, rank, local_i, size - 1)
                    - p_i
                    - p_next * decrement
                )
            elif n <= 2 * size - 2 and m >= 2 * size:
                piece = "cf"
                local_i = n - size
                p_i = kernel(size, rank, local_i, size)
                p_next = kernel(size, rank, local_i + 1, size)
                d_m = (
                    future_weight(alpha, size, m)
                    - future_weight(alpha, size, m + 1)
                )
                d_next = (
                    future_weight(alpha, size, m + 1)
                    - future_weight(alpha, size, m + 2)
                )
                expected = p_i * d_m - p_next * d_next
            elif n == 2 * size - 1 and m >= 2 * size:
                piece = "first_transition"
                p_last = kernel(size, rank, size - 1, size)
                diagonal = kernel(size, rank, size, size)
                d_m = (
                    future_weight(alpha, size, m)
                    - future_weight(alpha, size, m + 1)
                )
                d_next = (
                    future_weight(alpha, size, m + 1)
                    - future_weight(alpha, size, m + 2)
                )
                expected = p_last * d_m - diagonal * d_next
            else:
                piece = "ff"
                diagonal = kernel(size, rank, size, size)
                b_n = future_weight(alpha, size, n)
                b_next = future_weight(alpha, size, n + 1)
                d_n = b_n - b_next
                d_m = (
                    future_weight(alpha, size, m)
                    - future_weight(alpha, size, m + 1)
                )
                d_next = (
                    future_weight(alpha, size, m + 1)
                    - future_weight(alpha, size, m + 2)
                )
                expected = diagonal * (
                    b_n * (d_m - d_next) + d_n * d_next
                )
                if delta < -2e-12:
                    issues.append(
                        f"future curvature sign failed at K={size}"
                    )
                    return
            if not close(delta, expected, 2e-9):
                issues.append(
                    f"{piece} curvature formula failed at alpha={alpha}, "
                    f"K={size}, n={n}, h={shift}: {delta} != {expected}"
                )
                return
            pieces[piece] += abs(delta)
    if not close(direct, abel, 3e-9):
        issues.append(
            f"double Abel identity failed at alpha={alpha}, K={size}: "
            f"{direct} != {abel}"
        )
        return
    pair_direct = 0.0
    mu = mobius_values(4 * size + 2)
    for n in range(size + 1, 4 * size):
        for m in range(n + 1, 4 * size):
            pair_direct += (
                mu[n]
                * mu[m]
                * ordinary_kernel(alpha, size, rank, n, m)
            )
    if not close(direct, pair_direct, 3e-9):
        issues.append(f"pair/rectangle agreement failed at K={size}")
        return
    l1_bound = math.pi / 2.0 * (1.0 + math.log(2.0 * rank))
    piece_bounds = {
        "cc": 4.0
        * math.pi
        * size
        / (2 * size + 1) ** 2
        * l1_bound,
        "second_transition": math.pi**2 / (2.0 * (2 * size + 1))
        + 1.5 * rho,
        "cf": 2.5 * rho,
        "first_transition": 2.0 * rho,
        "ff": 2.0 * rho,
    }
    for piece, value in pieces.items():
        if value > piece_bounds[piece] + 2e-10:
            issues.append(
                f"{piece} variation bound failed at alpha={alpha}, "
                f"K={size}: {value} > {piece_bounds[piece]}"
            )
            return
    if not close(variation, sum(pieces.values()), 2e-9):
        issues.append(f"variation partition failed at K={size}")
        return
    intermediate_bound = (
        4.0
        * math.pi
        * size
        / (2 * size + 1) ** 2
        * l1_bound
        + math.pi**2 / (2.0 * (2 * size + 1))
        + 8.0 * rho
    )
    final_bound = (
        3.0 * math.pi**2 * (1.0 + math.log(2.0 * rank)) / size
    )
    if variation > intermediate_bound + 2e-10:
        issues.append(f"intermediate variation bound failed at K={size}")
        return
    if variation > final_bound + 2e-10:
        issues.append(f"global variation bound failed at K={size}")
        return
    if abs(direct) > variation * maximum + 2e-10:
        issues.append(f"planar maximum inequality failed at K={size}")
        return
    if direct * direct > variation * energy + 2e-10:
        issues.append(f"curvature Cauchy inequality failed at K={size}")
        return
    p_values = [
        kernel(size, rank, index, size)
        for index in range(0, size + 1)
    ]
    if any(
        p_values[index] > p_values[index + 1] + 2e-12
        for index in range(size)
    ):
        issues.append(f"endpoint column monotonicity failed at K={size}")
        return
    strip = sum(
        abs(
            kernel(size, rank, index, size)
            - kernel(size, rank, index, size - 1)
        )
        for index in range(1, max(1, size - 1))
    )
    if strip > math.pi**2 / (2.0 * (2 * size + 1)) + 2e-10:
        issues.append(f"endpoint strip bound failed at K={size}")
        return
    terminal = prefixes[-1][-1]
    block_sum = sum(mu[n] for n in range(size + 1, 4 * size))
    squarefree_count = sum(
        mu[n] * mu[n] for n in range(size + 1, 4 * size)
    )
    expected_terminal = (block_sum * block_sum - squarefree_count) // 2
    if terminal != expected_terminal:
        issues.append(f"full triangle identity failed at K={size}")


def check_countermodel(issues: list[str]) -> None:
    for length in range(2, 9):
        size = length * length
        matrix = [
            [
                int(row_index // length == column_index // length)
                for column_index in range(size)
            ]
            for row_index in range(size)
        ]
        row_prefix = max(
            sum(row[: endpoint + 1])
            for row in matrix
            for endpoint in range(size)
        )
        column_prefix = max(
            sum(matrix[row][column] for row in range(endpoint + 1))
            for column in range(size)
            for endpoint in range(size)
        )
        total = sum(sum(row) for row in matrix)
        if (
            row_prefix != length
            or column_prefix != length
            or total != length**3
        ):
            issues.append(f"axiswise countermodel failed at L={length}")
            return


def check_vaughan_algebra(issues: list[str]) -> None:
    size = 17
    first = [((5 * index + 2) % 11) - 5 for index in range(size)]
    second = [((7 * index + 1) % 13) - 6 for index in range(size)]
    combined = [
        -first[index] + second[index] for index in range(size)
    ]
    for endpoint in (4, 9, 16):
        for width in (1, 3, 7):
            def form(left: list[int], right: list[int]) -> int:
                return sum(
                    left[n] * right[n + shift]
                    for n in range(endpoint + 1)
                    for shift in range(1, width + 1)
                    if n + shift < size
                )

            direct = form(combined, combined)
            expanded = (
                form(first, first)
                + form(second, second)
                - form(first, second)
                - form(second, first)
            )
            if direct != expanded:
                issues.append("planar Vaughan expansion failed")
                return


def check_formal_core(issues: list[str]) -> None:
    text = FORMAL_CORE.read_text(encoding="utf-8")
    for marker in (
        "Corollary 11.22Z.5",
        "(11.22Z.86)",
        "(11.22Z.99)",
        "jensen_window_pf_mertens_planar_abel_handoff.md",
    ):
        if marker not in text:
            issues.append(f"formal core missing marker: {marker}")


def main() -> int:
    issues: list[str] = []
    for path in (RESULT, NOTE, BUILDER, FORMAL_CORE):
        if not path.exists():
            issues.append(f"missing file {path}")
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    note_text = NOTE.read_text(encoding="utf-8")
    check_schema(payload, note_text, issues)
    check_builder_reproduction(payload, note_text, issues)
    check_dirichlet_identity(issues)
    check_future_geometry(issues)
    for alpha, size in (
        (0.12, 3),
        (0.31, 4),
        (0.49, 6),
        (0.67, 9),
        (0.83, 13),
        (0.94, 18),
    ):
        check_planar_case(alpha, size, issues)
    check_countermodel(issues)
    check_vaughan_algebra(issues)
    check_formal_core(issues)
    if issues:
        for issue in issues:
            print(f"ERROR: {issue}", file=sys.stderr)
        return 1
    print(
        "validated Mertens planar Abel handoff: "
        "36 rows, 0 issues, 28 exact reductions, "
        "3 proof guards, 1 open gate"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
