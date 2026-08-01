#!/usr/bin/env python3
"""Check the truncated half-odd kernel and signed Vaughan handoff."""

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
    "jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.md"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_truncated_half_odd_kernel_handoff.py"
)

EXPECTED_KIND = (
    "jensen_window_pf_mertens_truncated_half_odd_kernel_handoff"
)
EXPECTED_STATUS = (
    "exact truncated half-odd kernel, uniformly summable diagonal, "
    "and signed off-diagonal/Vaughan reduction with one open "
    "arithmetic gate"
)
EXPECTED_ROW_COUNT = 45
EXPECTED_EXACT_COUNT = 39
EXPECTED_GUARD_COUNT = 2
EXPECTED_OPEN_COUNT = 1

REQUIRED_NOTE_MARKERS = (
    "P_(K,R)(i,j)",
    "C_R((i-j)pi/N)-C_R((i+j)pi/N)",
    "sin(2R*u)/(2sin u)",
    "S_R(u)>0",
    "P_(K,infinity)(i,j)",
    "pi^2/N^2*min(i,j)",
    "0<=P_(K,R)<=pi^2/N^2*C_K",
    "tr P_(K,R)",
    "G_(alpha,K,R)(n,m)",
    "0<=D_(alpha,K)<7pi^2/24",
    "sup_J sum_(K dyadic<=2^J)",
    "0=p_0<p_1<...<p_K",
    "A_(K,t)(p_t-p_(t-1))",
    "mu(n)=-u_(I,K)(n)+u_(II,K)(n)",
    "-2<u_I,G u_II>",
    "Mobius-specific off-diagonal gain",
)

REQUIRED_FORMAL_MARKERS = (
    "Corollary 11.22Z.3: Truncated Half-Odd Kernel Handoff",
    "P_(K,R)(i,j)",
    "P_(K,infinity)(i,j)",
    "0<=D_(alpha,K)<7pi^2/24",
    "R_alpha<infinity",
    "p_i:=P_(K,R)(i,K)",
    "mu(n)=-u_(I,K)(n)+u_(II,K)(n)",
    "Only the signed low-kernel off-diagonal estimate remains open",
)


def close(left: float, right: float, tolerance: float = 1e-11) -> bool:
    return abs(left - right) <= tolerance * max(
        1.0,
        abs(left),
        abs(right),
    )


def mobius_values(limit: int) -> list[int]:
    mu = [0] * (limit + 1)
    mu[1] = 1
    primes: list[int] = []
    composite = [False] * (limit + 1)
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            mu[n] = -1
        for prime in primes:
            value = n * prime
            if value > limit:
                break
            composite[value] = True
            if n % prime == 0:
                mu[value] = 0
                break
            mu[value] = -mu[n]
    return mu


def angle(size: int, mode: int) -> float:
    return (2 * mode - 1) * math.pi / (2 * size + 1)


def eigenvector_value(size: int, mode: int, index: int) -> float:
    return (
        2.0
        / math.sqrt(2 * size + 1)
        * math.sin(index * angle(size, mode))
    )


def kernel(size: int, cutoff: int, left: int, right: int) -> float:
    return sum(
        eigenvector_value(size, mode, left)
        * eigenvector_value(size, mode, right)
        / (2 * mode - 1) ** 2
        for mode in range(1, cutoff + 1)
    )


def cosine_kernel(
    size: int,
    cutoff: int,
    left: int,
    right: int,
) -> float:
    ambient = 2 * size + 1

    def cosine_sum(value: float) -> float:
        return sum(
            math.cos((2 * mode - 1) * value)
            / (2 * mode - 1) ** 2
            for mode in range(1, cutoff + 1)
        )

    return 2.0 / ambient * (
        cosine_sum((left - right) * math.pi / ambient)
        - cosine_sum((left + right) * math.pi / ambient)
    )


def sine_primitive(cutoff: int, value: float) -> float:
    return sum(
        math.sin((2 * mode - 1) * value) / (2 * mode - 1)
        for mode in range(1, cutoff + 1)
    )


def anchor_weight(alpha: float, size: int, n: int) -> float:
    if 2 * size <= n < 4 * size:
        return (
            (2.0 * size / n) ** (1.0 + alpha)
            * (4 * size - n)
            / (2.0 * size)
        )
    return 0.0


def position_weight(
    alpha: float,
    size: int,
    n: int,
) -> tuple[int, float]:
    if size < n < 2 * size:
        return n - size, 1.0
    if 2 * size <= n < 4 * size:
        return size, anchor_weight(alpha, size, n)
    return size, 0.0


def spectral_feature(
    alpha: float,
    size: int,
    mode: int,
    n: int,
) -> float:
    position, weight = position_weight(alpha, size, n)
    return weight * eigenvector_value(size, mode, position)


def ordinary_kernel(
    alpha: float,
    size: int,
    cutoff: int,
    left: int,
    right: int,
) -> float:
    left_position, left_weight = position_weight(alpha, size, left)
    right_position, right_weight = position_weight(alpha, size, right)
    return (
        left_weight
        * right_weight
        * kernel(size, cutoff, left_position, right_position)
    )


def check_schema(payload: dict, note_text: str, issues: list[str]) -> None:
    if payload.get("kind") != EXPECTED_KIND:
        issues.append("unexpected result kind")
    if payload.get("date") != "2026-07-24":
        issues.append("unexpected result date")
    if payload.get("status") != EXPECTED_STATUS:
        issues.append("unexpected result status")

    rows = payload.get("rows")
    if not isinstance(rows, list):
        issues.append("rows must be a list")
        return
    if len(rows) != EXPECTED_ROW_COUNT:
        issues.append(
            f"expected {EXPECTED_ROW_COUNT} rows, found {len(rows)}"
        )
    row_ids = [entry.get("id") for entry in rows if isinstance(entry, dict)]
    if len(row_ids) != len(set(row_ids)):
        issues.append("row ids are not unique")
    for index, entry in enumerate(rows):
        if not isinstance(entry, dict):
            issues.append(f"row {index} is not an object")
            continue
        for field in ("id", "role", "status", "statement", "proof_boundary"):
            value = entry.get(field)
            if not isinstance(value, str) or not value.strip():
                issues.append(f"row {index} has invalid {field}")

    audit = payload.get("audit", {})
    expected_scalars = {
        "row_count": EXPECTED_ROW_COUNT,
        "exact_reduction_count": EXPECTED_EXACT_COUNT,
        "literature_guard_count": 1,
        "proof_guard_count": EXPECTED_GUARD_COUNT,
        "conditional_calibration_count": 1,
        "finite_validation_count": 1,
        "open_signed_gate_count": EXPECTED_OPEN_COUNT,
    }
    for key, expected in expected_scalars.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key} expected {expected}, found {audit.get(key)}"
            )

    expected_true = (
        "kernel_formula_proved",
        "entrywise_positivity_proved",
        "min_kernel_completion_proved",
        "uniform_diagonal_bound_proved",
        "endpoint_abel_identity_proved",
        "signed_vaughan_gram_identity_proved",
    )
    for key in expected_true:
        if audit.get(key) is not True:
            issues.append(f"audit flag {key} must be true")

    expected_false = (
        "offdiagonal_mobius_gain_proved",
        "full_burnol_bound_proved",
        "rh_proved",
        "pf_infinity_proved",
        "lambda_le_zero_proved",
    )
    for key in expected_false:
        if audit.get(key) is not False:
            issues.append(f"audit flag {key} must be false")

    parameters = payload.get("parameters", {})
    expected_parameters = {
        "alpha_range": "0<alpha<1",
        "dyadic_blocks": "K=2^j",
        "cutoff": "R_(alpha,K)=ceil(K^(1-alpha/2))",
        "ambient_size": "N=2K+1",
        "mode_numbers": "q_r=2r-1",
        "uniform_diagonal_bound": "7pi^2/24",
    }
    for key, expected in expected_parameters.items():
        if parameters.get(key) != expected:
            issues.append(
                f"parameter {key} expected {expected}, "
                f"found {parameters.get(key)}"
            )

    for marker in REQUIRED_NOTE_MARKERS:
        if marker not in note_text:
            issues.append(f"note missing marker: {marker}")


def check_builder_reproduction(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    with tempfile.TemporaryDirectory(prefix="thok_check_") as temp_dir:
        result_path = Path(temp_dir) / "result.json"
        note_path = Path(temp_dir) / "note.md"
        process = subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--result",
                str(result_path),
                "--note",
                str(note_path),
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        if process.returncode != 0:
            issues.append(
                "builder reproduction failed: "
                + (process.stderr.strip() or process.stdout.strip())
            )
            return
        rebuilt = json.loads(result_path.read_text(encoding="utf-8"))
        rebuilt_note = note_path.read_text(encoding="utf-8")
        if rebuilt != payload:
            issues.append("builder reproduction changed the JSON payload")
        if rebuilt_note != note_text:
            issues.append("builder reproduction changed the note")


def cutoff_values(size: int) -> tuple[int, ...]:
    return tuple(
        sorted(
            {
                1,
                max(1, math.isqrt(size)),
                max(1, math.ceil(size**0.73)),
                size,
            }
        )
    )


def check_kernel_algebra(issues: list[str]) -> None:
    for size in (1, 2, 3, 5, 8, 16, 31):
        ambient = 2 * size + 1
        for cutoff in cutoff_values(size):
            trace = 0.0
            endpoint_values = [0.0]
            for left in range(1, size + 1):
                diagonal = kernel(size, cutoff, left, left)
                trace += diagonal
                if diagonal <= 0.0:
                    issues.append(
                        f"nonpositive diagonal at K={size}, R={cutoff}, "
                        f"i={left}"
                    )
                diagonal_bound = math.pi**2 * left / ambient**2
                if diagonal > diagonal_bound + 3e-12:
                    issues.append(
                        f"diagonal envelope failed at K={size}, "
                        f"R={cutoff}, i={left}"
                    )
                endpoint_values.append(kernel(size, cutoff, left, size))
                for right in range(1, size + 1):
                    direct = kernel(size, cutoff, left, right)
                    trigonometric = cosine_kernel(
                        size,
                        cutoff,
                        left,
                        right,
                    )
                    if not close(direct, trigonometric, tolerance=4e-12):
                        issues.append(
                            f"cosine kernel failed at K={size}, "
                            f"R={cutoff}, i={left}, j={right}"
                        )
                    if direct <= -3e-13:
                        issues.append(
                            f"entrywise positivity failed at K={size}, "
                            f"R={cutoff}, i={left}, j={right}"
                        )

            expected_trace = sum(
                1.0 / (2 * mode - 1) ** 2
                for mode in range(1, cutoff + 1)
            )
            if not close(trace, expected_trace, tolerance=5e-12):
                issues.append(
                    f"trace identity failed at K={size}, R={cutoff}"
                )
            for index in range(size):
                if (
                    endpoint_values[index + 1]
                    <= endpoint_values[index] - 3e-13
                ):
                    issues.append(
                        f"endpoint monotonicity failed at K={size}, "
                        f"R={cutoff}, i={index + 1}"
                    )

            if cutoff < size:
                for left in range(1, size + 1):
                    for right in range(1, size + 1):
                        difference = (
                            kernel(size, cutoff + 1, left, right)
                            - kernel(size, cutoff, left, right)
                        )
                        expected = (
                            eigenvector_value(size, cutoff + 1, left)
                            * eigenvector_value(
                                size,
                                cutoff + 1,
                                right,
                            )
                            / (2 * cutoff + 1) ** 2
                        )
                        if not close(
                            difference,
                            expected,
                            tolerance=5e-12,
                        ):
                            issues.append(
                                f"rank-one nesting failed at K={size}, "
                                f"R={cutoff}, i={left}, j={right}"
                            )


def check_sine_positivity(issues: list[str]) -> None:
    for cutoff in (1, 2, 3, 5, 8, 16, 31, 64):
        for index in range(1, 2048):
            value = math.pi * index / 2048.0
            sine_value = sine_primitive(cutoff, value)
            if sine_value <= -2e-13:
                issues.append(
                    f"sine positivity failed at R={cutoff}, "
                    f"sample={index}"
                )
                break
            symmetric = sine_primitive(cutoff, math.pi - value)
            if not close(sine_value, symmetric, tolerance=4e-12):
                issues.append(
                    f"sine symmetry failed at R={cutoff}, sample={index}"
                )
                break


def check_min_kernel_completion(issues: list[str]) -> None:
    for size in (1, 2, 3, 5, 8, 16, 31, 64):
        ambient = 2 * size + 1
        for left in range(1, size + 1):
            for right in range(1, size + 1):
                first = abs(left - right) * math.pi / ambient
                second = (left + right) * math.pi / ambient
                cosine_infinite = lambda value: (
                    math.pi**2 / 8.0 - math.pi * value / 4.0
                )
                completed = 2.0 / ambient * (
                    cosine_infinite(first) - cosine_infinite(second)
                )
                expected = (
                    math.pi**2
                    / ambient**2
                    * min(left, right)
                )
                if not close(completed, expected, tolerance=4e-12):
                    issues.append(
                        f"min-kernel completion failed at K={size}, "
                        f"i={left}, j={right}"
                    )

        for mode in range(1, size + 1):
            theta = angle(size, mode)
            gram_eigenvalue = 1.0 / (
                4.0 * math.sin(theta / 2.0) ** 2
            )
            kernel_eigenvalue = 1.0 / (2 * mode - 1) ** 2
            ratio = kernel_eigenvalue / gram_eigenvalue
            if ratio + 3e-14 < 4.0 / ambient**2:
                issues.append(
                    f"full lower comparison failed at K={size}, "
                    f"r={mode}"
                )
            if ratio > math.pi**2 / ambient**2 + 3e-14:
                issues.append(
                    f"full upper comparison failed at K={size}, "
                    f"r={mode}"
                )


def check_pullback_and_diagonal(issues: list[str]) -> None:
    mu = mobius_values(4096)
    diagonal_bound = 7.0 * math.pi**2 / 24.0
    for alpha in (0.09, 0.47, 0.91):
        for size in (1, 2, 3, 5, 8, 15, 24):
            support = list(range(size + 1, 4 * size))
            for cutoff in cutoff_values(size):
                spectral_energy = 0.0
                for mode in range(1, cutoff + 1):
                    coefficient = sum(
                        mu[n] * spectral_feature(
                            alpha,
                            size,
                            mode,
                            n,
                        )
                        for n in support
                    )
                    spectral_energy += (
                        coefficient**2 / (2 * mode - 1) ** 2
                    )

                diagonal = 0.0
                off_diagonal = 0.0
                ordinary_energy = 0.0
                for left_index, left in enumerate(support):
                    for right_index, right in enumerate(support):
                        entry = ordinary_kernel(
                            alpha,
                            size,
                            cutoff,
                            left,
                            right,
                        )
                        contribution = mu[left] * mu[right] * entry
                        ordinary_energy += contribution
                        if left_index == right_index:
                            diagonal += contribution
                        elif left_index < right_index:
                            off_diagonal += contribution

                if not close(
                    spectral_energy,
                    ordinary_energy,
                    tolerance=2e-10,
                ):
                    issues.append(
                        f"ordinary pullback failed at alpha={alpha}, "
                        f"K={size}, R={cutoff}"
                    )
                if not close(
                    ordinary_energy,
                    diagonal + 2.0 * off_diagonal,
                    tolerance=2e-10,
                ):
                    issues.append(
                        f"diagonal split failed at alpha={alpha}, "
                        f"K={size}, R={cutoff}"
                    )
                if diagonal < -2e-13 or diagonal > diagonal_bound + 2e-12:
                    issues.append(
                        f"uniform diagonal failed at alpha={alpha}, "
                        f"K={size}, R={cutoff}"
                    )

            future_square = sum(
                anchor_weight(alpha, size, n) ** 2
                for n in range(2 * size, 4 * size)
            )
            tent_bound = (
                (2 * size + 1)
                * (4 * size + 1)
                / (12.0 * size)
            )
            if future_square > tent_bound + 2e-13:
                issues.append(
                    f"future tent-square bound failed at alpha={alpha}, "
                    f"K={size}"
                )


def check_endpoint_abel(issues: list[str]) -> None:
    mu = mobius_values(4096)
    for size in (2, 3, 5, 8, 16, 31, 64):
        ambient = 2 * size + 1
        for cutoff in cutoff_values(size):
            endpoint = [
                0.0,
                *[
                    kernel(size, cutoff, index, size)
                    for index in range(1, size + 1)
                ],
            ]
            direct = sum(
                mu[size + index] * endpoint[index]
                for index in range(1, size)
            )
            suffix = [
                sum(
                    mu[size + index]
                    for index in range(start, size)
                )
                for start in range(1, size)
            ]
            abel = sum(
                suffix[start - 1]
                * (endpoint[start] - endpoint[start - 1])
                for start in range(1, size)
            )
            if not close(direct, abel, tolerance=5e-12):
                issues.append(
                    f"endpoint Abel identity failed at K={size}, "
                    f"R={cutoff}"
                )
            suffix_maximum = max((abs(value) for value in suffix), default=0)
            bound = (
                math.pi**2
                * size
                / ambient**2
                * suffix_maximum
            )
            if abs(direct) > bound + 3e-12:
                issues.append(
                    f"endpoint Abel bound failed at K={size}, "
                    f"R={cutoff}"
                )


def vaughan_coefficients(
    mu: list[int],
    base: int,
) -> tuple[dict[int, int], dict[int, int]]:
    u_cut = max(1, math.isqrt(base))
    v_cut = u_cut
    a_coefficients = {
        d: sum(
            mu[b] * mu[c]
            for b in range(1, u_cut + 1)
            for c in range(1, v_cut + 1)
            if b * c == d
        )
        for d in range(1, u_cut * v_cut + 1)
    }
    d_coefficients = {
        d: sum(
            mu[c]
            for c in range(v_cut + 1, d + 1)
            if d % c == 0
        )
        for d in range(v_cut + 1, (2 * base) // u_cut + 1)
    }

    type_i: dict[int, int] = {}
    type_ii: dict[int, int] = {}
    for n in range(base + 1, 2 * base + 1):
        type_i[n] = sum(
            coefficient
            for d, coefficient in a_coefficients.items()
            if n % d == 0
        )
        type_ii[n] = sum(
            coefficient * mu[n // d]
            for d, coefficient in d_coefficients.items()
            if n % d == 0 and n // d > u_cut
        )
    return type_i, type_ii


def inner_product(
    alpha: float,
    size: int,
    cutoff: int,
    left: dict[int, int],
    right: dict[int, int],
) -> float:
    return sum(
        left_value
        * right_value
        * ordinary_kernel(
            alpha,
            size,
            cutoff,
            left_n,
            right_n,
        )
        for left_n, left_value in left.items()
        if size < left_n < 4 * size
        for right_n, right_value in right.items()
        if size < right_n < 4 * size
    )


def check_signed_vaughan_gram(issues: list[str]) -> None:
    mu = mobius_values(4096)
    for alpha in (0.15, 0.57, 0.87):
        for size in (2, 3, 5, 8, 15, 24):
            type_i: dict[int, int] = {}
            type_ii: dict[int, int] = {}
            for base in (size, 2 * size):
                base_i, base_ii = vaughan_coefficients(mu, base)
                for n in range(base + 1, 2 * base + 1):
                    if mu[n] != -base_i[n] + base_ii[n]:
                        issues.append(
                            f"Vaughan coefficient identity failed at "
                            f"K={size}, X={base}, n={n}"
                        )
                        break
                type_i.update(base_i)
                type_ii.update(base_ii)

            support = range(size + 1, 4 * size)
            for cutoff in cutoff_values(size):
                direct = 0.0
                for mode in range(1, cutoff + 1):
                    coefficient = sum(
                        mu[n]
                        * spectral_feature(alpha, size, mode, n)
                        for n in support
                    )
                    direct += coefficient**2 / (2 * mode - 1) ** 2

                type_i_energy = inner_product(
                    alpha,
                    size,
                    cutoff,
                    type_i,
                    type_i,
                )
                type_ii_energy = inner_product(
                    alpha,
                    size,
                    cutoff,
                    type_ii,
                    type_ii,
                )
                cross = inner_product(
                    alpha,
                    size,
                    cutoff,
                    type_i,
                    type_ii,
                )
                reconstructed = (
                    type_i_energy + type_ii_energy - 2.0 * cross
                )
                if not close(direct, reconstructed, tolerance=3e-10):
                    issues.append(
                        f"signed Vaughan Gram failed at alpha={alpha}, "
                        f"K={size}, R={cutoff}"
                    )


def check_formal_core(issues: list[str]) -> None:
    if not FORMAL_CORE.exists():
        issues.append("formal core is missing")
        return
    text = FORMAL_CORE.read_text(encoding="utf-8")
    for marker in REQUIRED_FORMAL_MARKERS:
        if marker not in text:
            issues.append(f"formal core missing marker: {marker}")


def main() -> int:
    issues: list[str] = []
    if not RESULT.exists():
        issues.append(f"missing result: {RESULT}")
    if not NOTE.exists():
        issues.append(f"missing note: {NOTE}")
    if not BUILDER.exists():
        issues.append(f"missing builder: {BUILDER}")
    if issues:
        for issue in issues:
            print(f"ISSUE {issue}")
        return 1

    payload = json.loads(RESULT.read_text(encoding="utf-8"))
    note_text = NOTE.read_text(encoding="utf-8")
    check_schema(payload, note_text, issues)
    check_builder_reproduction(payload, note_text, issues)
    check_kernel_algebra(issues)
    check_sine_positivity(issues)
    check_min_kernel_completion(issues)
    check_pullback_and_diagonal(issues)
    check_endpoint_abel(issues)
    check_signed_vaughan_gram(issues)
    check_formal_core(issues)

    audit = payload.get("audit", {})
    print(
        "validated Mertens truncated half-odd kernel handoff: "
        f"{len(payload.get('rows', []))} rows, "
        f"{len(issues)} issues, "
        f"{audit.get('exact_reduction_count')} exact reductions, "
        f"{audit.get('proof_guard_count')} proof guards, "
        f"{audit.get('open_signed_gate_count')} open gate"
    )
    for issue in issues:
        print(f"ISSUE {issue}")
    return 1 if issues else 0


if __name__ == "__main__":
    raise SystemExit(main())
