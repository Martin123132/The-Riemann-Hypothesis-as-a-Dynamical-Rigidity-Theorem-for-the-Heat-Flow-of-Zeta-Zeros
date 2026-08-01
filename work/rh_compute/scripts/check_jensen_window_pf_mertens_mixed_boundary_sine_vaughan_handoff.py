#!/usr/bin/env python3
"""Check the mixed-boundary sine and modewise Vaughan handoff."""

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
    "jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.md"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff.py"
)

EXPECTED_KIND = (
    "jensen_window_pf_mertens_mixed_boundary_sine_vaughan_handoff"
)
EXPECTED_STATUS = (
    "exact mixed-boundary sine and modewise Vaughan reduction "
    "with one open additive-twist gate"
)
EXPECTED_ROW_COUNT = 35
EXPECTED_EXACT_COUNT = 27
EXPECTED_GUARD_COUNT = 3
EXPECTED_OPEN_COUNT = 1

REQUIRED_NOTE_MARKERS = (
    "C_K(i,j)=min(i,j)",
    "v_0=0",
    "v_(K+1)=v_K",
    "theta_(K,r):=(2r-1)pi/(2K+1)",
    "lambda_(K,r)",
    "1/[4sin^2(theta_(K,r)/2)]",
    "X_(K,r):=<x_K,phi_(K,r)>",
    "g_(K,r)(K+i):=phi_(K,r)(i)",
    "4/[pi^2(2r-1)^2]",
    "(4/pi^2)S_alpha<=E_alpha<=(9/4)S_alpha",
    "K^(1-alpha)log^(-2A)K",
    "K^(1/2+epsilon)",
    "X_(K,r)=-T_(I,K,r)+T_(II,K,r)",
    "(MBSV.21)",
    "all remain open",
)

REQUIRED_FORMAL_MARKERS = (
    "Corollary 11.22Z.1: Mixed-Boundary Sine Vaughan Handoff",
    "theta_(K,r):=(2r-1)pi/(2K+1)",
    "C_K phi_(K,r)=lambda_(K,r)phi_(K,r)",
    "(4/pi^2)S_alpha<=E_alpha<=(9/4)S_alpha",
    "X_(K,r)=-T_(I,K,r)+T_(II,K,r)",
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


def anchor_weight(alpha: float, size: int, n: int) -> float:
    if 2 * size <= n < 4 * size:
        return (
            (2.0 * size / n) ** (1.0 + alpha)
            * (4 * size - n)
            / (2.0 * size)
        )
    return 0.0


def angle(size: int, mode: int) -> float:
    return (2 * mode - 1) * math.pi / (2 * size + 1)


def eigenvector(size: int, mode: int) -> list[float]:
    theta = angle(size, mode)
    factor = 2.0 / math.sqrt(2 * size + 1)
    return [
        factor * math.sin(index * theta)
        for index in range(1, size + 1)
    ]


def eigenvalue(size: int, mode: int) -> float:
    return 1.0 / (
        4.0 * math.sin(angle(size, mode) / 2.0) ** 2
    )


def covariance(size: int) -> list[list[float]]:
    return [
        [float(min(i, j)) for j in range(1, size + 1)]
        for i in range(1, size + 1)
    ]


def inverse_laplacian(size: int) -> list[list[float]]:
    matrix = [[0.0] * size for _ in range(size)]
    for index in range(size):
        matrix[index][index] = 1.0 if index == size - 1 else 2.0
        if index + 1 < size:
            matrix[index][index + 1] = -1.0
            matrix[index + 1][index] = -1.0
    return matrix


def matrix_vector(
    matrix: list[list[float]],
    vector: list[float],
) -> list[float]:
    return [
        sum(value * entry for value, entry in zip(row, vector))
        for row in matrix
    ]


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def spectral_feature(
    alpha: float,
    size: int,
    mode: int,
    n: int,
) -> float:
    vector = eigenvector(size, mode)
    if size < n < 2 * size:
        return vector[n - size - 1]
    if 2 * size <= n < 4 * size:
        return anchor_weight(alpha, size, n) * vector[-1]
    return 0.0


def check_schema(payload: dict, note_text: str, issues: list[str]) -> None:
    if payload.get("kind") != EXPECTED_KIND:
        issues.append("unexpected result kind")
    if payload.get("date") != "2026-07-23":
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
        "literature_guard_count": 2,
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
        "mixed_boundary_spectrum_proved",
        "spectral_mobius_features_proved",
        "odd_frequency_equivalence_proved",
        "modewise_vaughan_proved",
    )
    for key in expected_true:
        if audit.get(key) is not True:
            issues.append(f"audit flag {key} must be true")
    expected_false = (
        "square_root_additive_twist_proved",
        "signed_type_i_ii_gain_proved",
        "full_burnol_bound_proved",
        "rh_proved",
        "pf_infinity_proved",
        "lambda_le_zero_proved",
    )
    for key in expected_false:
        if audit.get(key) is not False:
            issues.append(f"audit flag {key} must be false")

    for marker in REQUIRED_NOTE_MARKERS:
        if marker not in note_text:
            issues.append(f"note missing marker: {marker}")


def check_builder_reproduction(
    payload: dict,
    note_text: str,
    issues: list[str],
) -> None:
    with tempfile.TemporaryDirectory(prefix="mbsv_check_") as temp_dir:
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


def check_mixed_spectrum(issues: list[str]) -> None:
    for size in (1, 2, 3, 7, 16, 31, 64):
        cov = covariance(size)
        lap = inverse_laplacian(size)
        for row in range(size):
            for column in range(size):
                product = sum(
                    cov[row][middle] * lap[middle][column]
                    for middle in range(size)
                )
                expected = 1.0 if row == column else 0.0
                if not close(product, expected, tolerance=2e-12):
                    issues.append(
                        f"CQ inverse failed at K={size}, "
                        f"row={row + 1}, column={column + 1}"
                    )
                    break

        vectors = [
            eigenvector(size, mode)
            for mode in range(1, size + 1)
        ]
        for left in range(size):
            for right in range(size):
                expected = 1.0 if left == right else 0.0
                if not close(
                    dot(vectors[left], vectors[right]),
                    expected,
                    tolerance=3e-12,
                ):
                    issues.append(
                        f"sine orthogonality failed at K={size}, "
                        f"r={left + 1}, q={right + 1}"
                    )
                    break

        for mode, vector in enumerate(vectors, start=1):
            cov_vector = matrix_vector(cov, vector)
            target = [
                eigenvalue(size, mode) * value
                for value in vector
            ]
            for index, (left, right) in enumerate(
                zip(cov_vector, target)
            ):
                if not close(left, right, tolerance=3e-11):
                    issues.append(
                        f"covariance eigenpair failed at K={size}, "
                        f"r={mode}, i={index + 1}"
                    )
                    break

            odd = 2 * mode - 1
            normalized = eigenvalue(size, mode) / (size * size)
            lower = 4.0 / (math.pi * math.pi * odd * odd)
            upper = 9.0 / (4.0 * odd * odd)
            if normalized + 2e-13 < lower:
                issues.append(
                    f"eigenvalue lower bound failed at K={size}, "
                    f"r={mode}"
                )
            if normalized > upper + 2e-13:
                issues.append(
                    f"eigenvalue upper bound failed at K={size}, "
                    f"r={mode}"
                )

        arbitrary = [
            ((17 * index + 5 * index * index) % 23 - 11) / 7.0
            for index in range(1, size + 1)
        ]
        direct = dot(arbitrary, matrix_vector(cov, arbitrary))
        spectral = sum(
            eigenvalue(size, mode)
            * dot(arbitrary, vectors[mode - 1]) ** 2
            for mode in range(1, size + 1)
        )
        if not close(direct, spectral, tolerance=4e-11):
            issues.append(f"arbitrary spectral energy failed at K={size}")


def check_mobius_features(issues: list[str]) -> None:
    mu = mobius_values(2048)
    for alpha in (0.09, 0.51, 0.91):
        for size in (1, 2, 4, 9, 24, 47):
            beta = sum(
                mu[n] * anchor_weight(alpha, size, n)
                for n in range(2 * size, 4 * size)
            )
            coefficients = [
                *[
                    float(mu[size + index])
                    for index in range(1, size)
                ],
                beta,
            ]
            suffix = [
                sum(coefficients[t - 1 :])
                for t in range(1, size + 1)
            ]
            energy = sum(value * value for value in suffix)
            spectral_energy = 0.0
            odd_square = 0.0
            for mode in range(1, size + 1):
                vector = eigenvector(size, mode)
                coefficient_mode = dot(coefficients, vector)
                feature_mode = sum(
                    mu[n] * spectral_feature(alpha, size, mode, n)
                    for n in range(size + 1, 4 * size)
                )
                if not close(
                    coefficient_mode,
                    feature_mode,
                    tolerance=3e-11,
                ):
                    issues.append(
                        f"spectral feature unfolding failed at "
                        f"alpha={alpha}, K={size}, r={mode}"
                    )
                    break
                envelope = 2.0 / math.sqrt(2 * size + 1)
                for n in range(size + 1, 4 * size):
                    if (
                        abs(spectral_feature(alpha, size, mode, n))
                        > envelope * (1.0 + 2e-12)
                    ):
                        issues.append(
                            f"spectral feature envelope failed at "
                            f"alpha={alpha}, K={size}, r={mode}"
                        )
                        break
                spectral_energy += (
                    eigenvalue(size, mode) * coefficient_mode**2
                )
                odd_square += (
                    coefficient_mode**2 / (2 * mode - 1) ** 2
                )

            if not close(energy, spectral_energy, tolerance=4e-11):
                issues.append(
                    f"Mobius spectral reconstruction failed at "
                    f"alpha={alpha}, K={size}"
                )
            normalized_energy = energy / (size * size)
            if normalized_energy + 3e-11 < 4.0 / math.pi**2 * odd_square:
                issues.append(
                    f"odd-square lower comparison failed at "
                    f"alpha={alpha}, K={size}"
                )
            if normalized_energy > 2.25 * odd_square + 3e-11:
                issues.append(
                    f"odd-square upper comparison failed at "
                    f"alpha={alpha}, K={size}"
                )


def vaughan_interval_mode(
    mu: list[int],
    alpha: float,
    size: int,
    mode: int,
    base: int,
) -> tuple[float, float, float]:
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

    direct = sum(
        mu[n] * spectral_feature(alpha, size, mode, n)
        for n in range(base + 1, 2 * base + 1)
    )
    type_i = sum(
        coefficient
        * sum(
            spectral_feature(alpha, size, mode, d * w)
            for w in range(
                base // d + 1,
                (2 * base) // d + 1,
            )
        )
        for d, coefficient in a_coefficients.items()
    )
    type_ii = sum(
        coefficient
        * sum(
            mu[w] * spectral_feature(alpha, size, mode, d * w)
            for w in range(
                max(u_cut, base // d) + 1,
                (2 * base) // d + 1,
            )
        )
        for d, coefficient in d_coefficients.items()
    )
    return direct, type_i, type_ii


def check_modewise_vaughan(issues: list[str]) -> None:
    mu = mobius_values(2048)
    for alpha in (0.15, 0.57, 0.87):
        for size in (2, 4, 8, 15, 24, 39):
            for mode in range(1, size + 1):
                full_direct = 0.0
                full_type_i = 0.0
                full_type_ii = 0.0
                for base in (size, 2 * size):
                    direct, type_i, type_ii = vaughan_interval_mode(
                        mu,
                        alpha,
                        size,
                        mode,
                        base,
                    )
                    if not close(
                        direct,
                        -type_i + type_ii,
                        tolerance=4e-11,
                    ):
                        issues.append(
                            f"mode Vaughan failed at alpha={alpha}, "
                            f"K={size}, r={mode}, X={base}"
                        )
                        break
                    full_direct += direct
                    full_type_i += type_i
                    full_type_ii += type_ii

                full_feature = sum(
                    mu[n] * spectral_feature(
                        alpha,
                        size,
                        mode,
                        n,
                    )
                    for n in range(size + 1, 4 * size)
                )
                if not close(
                    full_direct,
                    full_feature,
                    tolerance=4e-11,
                ):
                    issues.append(
                        f"mode interval cover failed at alpha={alpha}, "
                        f"K={size}, r={mode}"
                    )
                if not close(
                    full_feature,
                    -full_type_i + full_type_ii,
                    tolerance=4e-11,
                ):
                    issues.append(
                        f"full mode Vaughan failed at alpha={alpha}, "
                        f"K={size}, r={mode}"
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
    check_mixed_spectrum(issues)
    check_mobius_features(issues)
    check_modewise_vaughan(issues)
    check_formal_core(issues)

    audit = payload.get("audit", {})
    print(
        "validated Mertens mixed-boundary sine/Vaughan handoff: "
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
