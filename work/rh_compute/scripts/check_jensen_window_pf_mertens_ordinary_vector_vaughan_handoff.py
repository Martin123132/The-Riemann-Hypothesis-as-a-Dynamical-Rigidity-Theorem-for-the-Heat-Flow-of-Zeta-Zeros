#!/usr/bin/env python3
"""Check the ordinary-Mobius feature-Gram and vector Vaughan handoff."""

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
    "jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.md"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_ordinary_vector_vaughan_handoff.py"
)

EXPECTED_KIND = "jensen_window_pf_mertens_ordinary_vector_vaughan_handoff"
EXPECTED_STATUS = (
    "exact ordinary-Mobius feature Gram and pre-square vector Vaughan "
    "reduction with one open signed gate"
)
EXPECTED_ROW_COUNT = 34
EXPECTED_EXACT_COUNT = 28
EXPECTED_GUARD_COUNT = 3
EXPECTED_OPEN_COUNT = 1

REQUIRED_NOTE_MARKERS = (
    "a_(K,t)(K+i):=1_(t<=i)",
    "B_K:=sum_(K<n<4K)mu(n)a_K(n)",
    "H_K=||B_K||_2^2",
    "L_K(K+i,K+j)=min(i,j)",
    "L_K(n,m)=K*b_K(n)b_K(m)",
    "(14K^2+1)/12",
    "K^(-2-alpha)C_(vec,K)",
    "V_(I,K,X)",
    "V_(II,K,X)",
    "B_K=-V_(I,K)+V_(II,K)",
    "Vaughan before squaring",
    "(MOVH.15)",
    "K^3*log^(-2A)K",
    "one full power",
    "is not a proof of RH",
)

REQUIRED_FORMAL_MARKERS = (
    "Lemma 11.22Z: Ordinary-Mobius Vector Vaughan Handoff",
    "a_(K,t)(K+i):=1_(t<=i)",
    "L_K(K+i,K+j)=min(i,j)",
    "(14K^2+1)/12",
    "B_K=-V_(I,K)+V_(II,K)",
    "Vaughan before squaring",
    "K^3*log^(-2A)K",
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


def feature(alpha: float, size: int, n: int) -> list[float]:
    if size < n < 2 * size:
        endpoint = n - size
        return [
            1.0 if t <= endpoint else 0.0
            for t in range(1, size + 1)
        ]
    if 2 * size <= n < 4 * size:
        value = anchor_weight(alpha, size, n)
        return [value] * size
    return [0.0] * size


def add_scaled(
    target: list[float],
    source: list[float],
    scale: float,
) -> None:
    for index, value in enumerate(source):
        target[index] += scale * value


def squared_norm(vector: list[float]) -> float:
    return sum(value * value for value in vector)


def expected_kernel(
    alpha: float,
    size: int,
    left: int,
    right: int,
) -> float:
    left_current = size < left < 2 * size
    right_current = size < right < 2 * size
    left_future = 2 * size <= left < 4 * size
    right_future = 2 * size <= right < 4 * size
    if left_current and right_current:
        return float(min(left - size, right - size))
    if left_current and right_future:
        return (left - size) * anchor_weight(alpha, size, right)
    if left_future and right_current:
        return (right - size) * anchor_weight(alpha, size, left)
    if left_future and right_future:
        return (
            size
            * anchor_weight(alpha, size, left)
            * anchor_weight(alpha, size, right)
        )
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
        "finite_validation_count": 1,
        "open_signed_gate_count": EXPECTED_OPEN_COUNT,
    }
    for key, expected in expected_scalars.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key} expected {expected}, found {audit.get(key)}"
            )
    expected_true = (
        "ordinary_feature_vector_proved",
        "explicit_gram_kernel_proved",
        "diagonal_summable",
        "componentwise_vaughan_proved",
        "pre_square_signed_difference_preserved",
    )
    for key in expected_true:
        if audit.get(key) is not True:
            issues.append(f"audit flag {key} must be true")
    expected_false = (
        "vector_type_i_ii_gain_proved",
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
    with tempfile.TemporaryDirectory(prefix="movh_check_") as temp_dir:
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


def check_feature_gram(issues: list[str]) -> None:
    mu = mobius_values(1024)
    for alpha in (0.07, 0.43, 0.93):
        for size in (1, 2, 4, 7, 16, 31):
            beta = sum(
                mu[n] * anchor_weight(alpha, size, n)
                for n in range(2 * size, 4 * size)
            )
            direct = [0.0] * size
            for n in range(size + 1, 4 * size):
                add_scaled(direct, feature(alpha, size, n), mu[n])
            suffix = [
                beta
                + sum(
                    mu[size + index]
                    for index in range(t, size)
                )
                for t in range(1, size + 1)
            ]
            for index, (left, right) in enumerate(zip(direct, suffix)):
                if not close(left, right, tolerance=2e-12):
                    issues.append(
                        f"feature suffix failed at alpha={alpha}, "
                        f"K={size}, coordinate={index + 1}"
                    )
                    break

            gram_sum = 0.0
            diagonal = 0.0
            offdiagonal = 0.0
            for left in range(size + 1, 4 * size):
                left_feature = feature(alpha, size, left)
                for right in range(left, 4 * size):
                    right_feature = feature(alpha, size, right)
                    actual = sum(
                        a * b
                        for a, b in zip(left_feature, right_feature)
                    )
                    expected = expected_kernel(
                        alpha,
                        size,
                        left,
                        right,
                    )
                    if not close(actual, expected, tolerance=2e-12):
                        issues.append(
                            f"explicit Gram kernel failed at "
                            f"alpha={alpha}, K={size}, "
                            f"n={left}, m={right}"
                        )
                        break
                    if actual < -1e-14 or actual > size * (1.0 + 1e-12):
                        issues.append(
                            f"Gram envelope failed at alpha={alpha}, "
                            f"K={size}, n={left}, m={right}"
                        )
                    contribution = mu[left] * mu[right] * actual
                    if left == right:
                        diagonal += contribution
                        gram_sum += contribution
                    else:
                        offdiagonal += contribution
                        gram_sum += 2.0 * contribution

            energy = squared_norm(direct)
            if not close(energy, gram_sum, tolerance=3e-11):
                issues.append(
                    f"Gram expansion failed at alpha={alpha}, K={size}"
                )
            if not close(
                energy,
                diagonal + 2.0 * offdiagonal,
                tolerance=3e-11,
            ):
                issues.append(
                    f"signed split failed at alpha={alpha}, K={size}"
                )
            bound = 1.25 * size * size
            if diagonal < -1e-14 or diagonal > bound * (1.0 + 1e-12):
                issues.append(
                    f"diagonal bound failed at alpha={alpha}, K={size}"
                )

            envelope_sum = (
                (2 * size + 1)
                * (4 * size + 1)
                / 12.0
            )
            actual_future_envelope = size * sum(
                ((4 * size - n) / (2.0 * size)) ** 2
                for n in range(2 * size, 4 * size)
            )
            if not close(
                actual_future_envelope,
                envelope_sum,
                tolerance=2e-12,
            ):
                issues.append(
                    f"future diagonal sum failed at alpha={alpha}, "
                    f"K={size}"
                )

            arbitrary = {
                n: ((11 * n + n * n) % 13 - 6) / 5.0
                for n in range(size + 1, 4 * size)
            }
            arbitrary_vector = [0.0] * size
            for n, coefficient in arbitrary.items():
                add_scaled(
                    arbitrary_vector,
                    feature(alpha, size, n),
                    coefficient,
                )
            arbitrary_gram = sum(
                arbitrary[left]
                * arbitrary[right]
                * expected_kernel(alpha, size, left, right)
                for left in arbitrary
                for right in arbitrary
            )
            if not close(
                squared_norm(arbitrary_vector),
                arbitrary_gram,
                tolerance=3e-11,
            ):
                issues.append(
                    f"arbitrary-vector Gram failed at alpha={alpha}, "
                    f"K={size}"
                )


def vaughan_interval_vectors(
    mu: list[int],
    alpha: float,
    size: int,
    base: int,
) -> tuple[list[float], list[float], list[float]]:
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

    direct = [0.0] * size
    for n in range(base + 1, 2 * base + 1):
        add_scaled(direct, feature(alpha, size, n), mu[n])

    type_i = [0.0] * size
    for d, coefficient in a_coefficients.items():
        for w in range(base // d + 1, (2 * base) // d + 1):
            add_scaled(
                type_i,
                feature(alpha, size, d * w),
                coefficient,
            )

    type_ii = [0.0] * size
    for d, coefficient in d_coefficients.items():
        for w in range(
            max(u_cut, base // d) + 1,
            (2 * base) // d + 1,
        ):
            add_scaled(
                type_ii,
                feature(alpha, size, d * w),
                coefficient * mu[w],
            )
    return direct, type_i, type_ii


def check_vector_vaughan(issues: list[str]) -> None:
    mu = mobius_values(2048)
    for alpha in (0.17, 0.59, 0.89):
        for size in (2, 4, 8, 15, 24, 39):
            full_direct = [0.0] * size
            full_type_i = [0.0] * size
            full_type_ii = [0.0] * size
            for base in (size, 2 * size):
                direct, type_i, type_ii = vaughan_interval_vectors(
                    mu,
                    alpha,
                    size,
                    base,
                )
                for coordinate in range(size):
                    reconstructed = (
                        -type_i[coordinate] + type_ii[coordinate]
                    )
                    if not close(
                        direct[coordinate],
                        reconstructed,
                        tolerance=3e-11,
                    ):
                        issues.append(
                            f"componentwise Vaughan failed at "
                            f"alpha={alpha}, K={size}, X={base}, "
                            f"t={coordinate + 1}"
                        )
                        break
                    full_direct[coordinate] += direct[coordinate]
                    full_type_i[coordinate] += type_i[coordinate]
                    full_type_ii[coordinate] += type_ii[coordinate]

            beta = sum(
                mu[n] * anchor_weight(alpha, size, n)
                for n in range(2 * size, 4 * size)
            )
            suffix = [
                beta
                + sum(
                    mu[size + index]
                    for index in range(t, size)
                )
                for t in range(1, size + 1)
            ]
            reconstructed_full = [
                -left + right
                for left, right in zip(full_type_i, full_type_ii)
            ]
            for coordinate in range(size):
                if not close(
                    full_direct[coordinate],
                    suffix[coordinate],
                    tolerance=3e-11,
                ):
                    issues.append(
                        f"full feature cover failed at alpha={alpha}, "
                        f"K={size}, t={coordinate + 1}"
                    )
                    break
                if not close(
                    reconstructed_full[coordinate],
                    suffix[coordinate],
                    tolerance=3e-11,
                ):
                    issues.append(
                        f"full vector Vaughan failed at alpha={alpha}, "
                        f"K={size}, t={coordinate + 1}"
                    )
                    break
            if not close(
                squared_norm(suffix),
                squared_norm(reconstructed_full),
                tolerance=3e-11,
            ):
                issues.append(
                    f"pre-square Vaughan norm failed at alpha={alpha}, "
                    f"K={size}"
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
    check_feature_gram(issues)
    check_vector_vaughan(issues)
    check_formal_core(issues)

    audit = payload.get("audit", {})
    print(
        "validated Mertens ordinary-vector Vaughan handoff: "
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
