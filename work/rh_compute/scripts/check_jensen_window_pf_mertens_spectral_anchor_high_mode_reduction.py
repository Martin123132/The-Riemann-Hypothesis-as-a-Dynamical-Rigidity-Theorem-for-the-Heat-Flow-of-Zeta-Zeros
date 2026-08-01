#!/usr/bin/env python3
"""Check the spectral anchor/high-mode reduction."""

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
    "jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.json"
)
NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.md"
)
FORMAL_CORE = REPO_ROOT / "outputs/formal_core.md"
BUILDER = (
    REPO_ROOT
    / "work/rh_compute/scripts/"
    "jensen_window_pf_mertens_spectral_anchor_high_mode_reduction.py"
)

EXPECTED_KIND = (
    "jensen_window_pf_mertens_spectral_anchor_high_mode_reduction"
)
EXPECTED_STATUS = (
    "exact spectral anchor separation and summable high-mode "
    "reduction with one open low-mode gate"
)
EXPECTED_ROW_COUNT = 24
EXPECTED_EXACT_COUNT = 16
EXPECTED_GUARD_COUNT = 2
EXPECTED_OPEN_COUNT = 1

REQUIRED_NOTE_MARKERS = (
    "X_(K,r)=Y_(K,r)+beta_K*h_(K,r)",
    "K*theta_(K,r)",
    "cos(theta_(K,r)/2)",
    "sum_(r=1)^K|Y_(K,r)|^2",
    "sum_(r>R)1/(2r-1)^2<=1/(2R+1)",
    "<=2K/(2R+1)^2",
    "+4|beta_K|^2/[K(2R+1)]",
    "R_(alpha,K):=ceil(K^(1-alpha/2))",
    "<=2K^(-1)+8K^(-alpha/2)",
    "O(K^(1-alpha/2))` modes remain",
    "|beta_K|=O_epsilon(K^(1/2+epsilon))",
    "O(K^((1-alpha)/2+eta))",
    "beta_K>=2K/9",
    "Low-mode square-root cancellation",
)

REQUIRED_FORMAL_MARKERS = (
    "Corollary 11.22Z.2: Spectral Anchor/High-Mode Reduction",
    "X_(K,r)=Y_(K,r)+beta_K h_(K,r)",
    "h_(K,r)=phi_(K,r)(K)",
    "sum_(r>R)|X_(K,r)|^2/(2r-1)^2",
    "R_(alpha,K):=ceil(K^(1-alpha/2))",
    "2K^(-1)+8K^(-alpha/2)",
    "Only the reduced low-mode square remains open",
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


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


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
        "literature_guard_count": 1,
        "proof_guard_count": EXPECTED_GUARD_COUNT,
        "conditional_calibration_count": 3,
        "finite_validation_count": 1,
        "open_signed_gate_count": EXPECTED_OPEN_COUNT,
    }
    for key, expected in expected_scalars.items():
        if audit.get(key) != expected:
            issues.append(
                f"audit {key} expected {expected}, found {audit.get(key)}"
            )

    expected_true = (
        "current_anchor_split_proved",
        "unconditional_high_modes_summable",
        "unconditional_low_mode_equivalence_proved",
    )
    for key in expected_true:
        if audit.get(key) is not True:
            issues.append(f"audit flag {key} must be true")

    expected_false = (
        "square_root_anchor_proved",
        "low_mode_gain_proved",
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
        "unconditional_cutoff": "ceil(K^(1-alpha/2))",
        "current_tail_decay": "R^(-2)",
        "anchor_tail_decay": "R^(-1)",
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
    with tempfile.TemporaryDirectory(prefix="sahr_check_") as temp_dir:
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


def check_exact_split(issues: list[str]) -> None:
    mu = mobius_values(4096)
    for alpha in (0.09, 0.47, 0.91):
        for size in (1, 2, 4, 9, 24, 47):
            beta = sum(
                mu[n] * anchor_weight(alpha, size, n)
                for n in range(2 * size, 4 * size)
            )
            current = [
                float(mu[size + index])
                for index in range(1, size)
            ]
            current_energy = sum(value * value for value in current)
            spectral_current_energy = 0.0

            arbitrary = [
                ((13 * index + 7 * index * index) % 29 - 14) / 9.0
                for index in range(1, size)
            ]
            arbitrary_beta = ((11 * size + 3) % 17 - 8) / 5.0

            for mode in range(1, size + 1):
                vector = eigenvector(size, mode)
                theta = angle(size, mode)
                endpoint = vector[-1]
                expected_endpoint = (
                    2.0
                    * (-1.0) ** (mode - 1)
                    / math.sqrt(2 * size + 1)
                    * math.cos(theta / 2.0)
                )
                if not close(endpoint, expected_endpoint, tolerance=3e-12):
                    issues.append(
                        f"endpoint phase failed at alpha={alpha}, "
                        f"K={size}, r={mode}"
                    )

                endpoint_bound = 2.0 / size
                if endpoint * endpoint > endpoint_bound + 3e-13:
                    issues.append(
                        f"endpoint envelope failed at K={size}, r={mode}"
                    )

                current_mode = dot(current, vector[:-1])
                split_mode = current_mode + beta * endpoint
                direct_mode = dot([*current, beta], vector)
                if not close(split_mode, direct_mode, tolerance=3e-12):
                    issues.append(
                        f"Mobius anchor split failed at alpha={alpha}, "
                        f"K={size}, r={mode}"
                    )
                spectral_current_energy += current_mode * current_mode

                arbitrary_current = dot(arbitrary, vector[:-1])
                arbitrary_split = (
                    arbitrary_current + arbitrary_beta * endpoint
                )
                arbitrary_direct = dot(
                    [*arbitrary, arbitrary_beta],
                    vector,
                )
                if not close(
                    arbitrary_split,
                    arbitrary_direct,
                    tolerance=3e-12,
                ):
                    issues.append(
                        f"arbitrary anchor split failed at K={size}, "
                        f"r={mode}"
                    )

            if not close(
                spectral_current_energy,
                current_energy,
                tolerance=5e-12,
            ):
                issues.append(
                    f"current Parseval failed at alpha={alpha}, K={size}"
                )
            if current_energy > size + 1e-12:
                issues.append(
                    f"Mobius current envelope failed at K={size}"
                )


def check_odd_tail(issues: list[str]) -> None:
    for cutoff in (1, 2, 3, 7, 16, 47, 128):
        telescoping = 1.0 / (4.0 * cutoff)
        claimed = 1.0 / (2.0 * cutoff + 1.0)
        if telescoping > claimed + 1e-15:
            issues.append(f"odd-tail comparison failed at R={cutoff}")
        for mode in range(cutoff + 1, cutoff + 257):
            odd_term = 1.0 / (2 * mode - 1) ** 2
            telescoping_term = 1.0 / (
                (2 * mode - 2) * (2 * mode)
            )
            if odd_term > telescoping_term + 1e-16:
                issues.append(
                    f"odd-tail termwise bound failed at "
                    f"R={cutoff}, r={mode}"
                )
                break


def cutoff_values(size: int) -> tuple[int, ...]:
    values = {
        1,
        max(1, math.ceil(math.sqrt(size))),
        max(1, size - 1),
        size,
    }
    return tuple(sorted(value for value in values if value <= size))


def check_high_mode_bounds(issues: list[str]) -> None:
    mu = mobius_values(4096)
    for alpha in (0.09, 0.47, 0.91):
        for size in (2, 4, 9, 24, 47):
            beta = sum(
                mu[n] * anchor_weight(alpha, size, n)
                for n in range(2 * size, 4 * size)
            )
            current_modes: list[float] = []
            anchor_modes: list[float] = []
            total_modes: list[float] = []
            for mode in range(1, size + 1):
                vector = eigenvector(size, mode)
                current_mode = sum(
                    mu[size + index] * vector[index - 1]
                    for index in range(1, size)
                )
                anchor_mode = beta * vector[-1]
                current_modes.append(current_mode)
                anchor_modes.append(anchor_mode)
                total_modes.append(current_mode + anchor_mode)

            for cutoff in cutoff_values(size):
                modes = range(cutoff + 1, size + 1)
                current_tail = sum(
                    current_modes[mode - 1] ** 2 / (2 * mode - 1) ** 2
                    for mode in modes
                )
                anchor_tail = sum(
                    anchor_modes[mode - 1] ** 2 / (2 * mode - 1) ** 2
                    for mode in modes
                )
                total_tail = sum(
                    total_modes[mode - 1] ** 2 / (2 * mode - 1) ** 2
                    for mode in modes
                )
                current_bound = size / (2 * cutoff + 1) ** 2
                anchor_bound = (
                    2.0 * beta * beta
                    / (size * (2 * cutoff + 1))
                )
                total_bound = (
                    2.0 * size / (2 * cutoff + 1) ** 2
                    + 4.0 * beta * beta
                    / (size * (2 * cutoff + 1))
                )
                label = (
                    f"alpha={alpha}, K={size}, R={cutoff}"
                )
                if current_tail > current_bound + 2e-12:
                    issues.append(f"current tail failed at {label}")
                if anchor_tail > anchor_bound + 2e-12:
                    issues.append(f"anchor tail failed at {label}")
                if total_tail > total_bound + 3e-12:
                    issues.append(f"total tail failed at {label}")


def check_unconditional_cutoff(issues: list[str]) -> None:
    mu = mobius_values(8192)
    for alpha in (0.09, 0.47, 0.91):
        for size in (2, 4, 9, 24, 47, 128, 257):
            cutoff = math.ceil(size ** (1.0 - alpha / 2.0))
            if not 1 <= cutoff <= size:
                issues.append(
                    f"unconditional cutoff out of range at "
                    f"alpha={alpha}, K={size}"
                )
                continue
            beta = sum(
                mu[n] * anchor_weight(alpha, size, n)
                for n in range(2 * size, 4 * size)
            )
            if abs(beta) > 2.0 * size + 2e-12:
                issues.append(
                    f"trivial anchor bound failed at alpha={alpha}, "
                    f"K={size}"
                )
            total_tail = 0.0
            for mode in range(cutoff + 1, size + 1):
                vector = eigenvector(size, mode)
                current_mode = sum(
                    mu[size + index] * vector[index - 1]
                    for index in range(1, size)
                )
                total_mode = current_mode + beta * vector[-1]
                total_tail += total_mode**2 / (2 * mode - 1) ** 2
            normalized = size ** (-alpha) * total_tail
            claimed = 2.0 * size ** (-1.0) + 8.0 * size ** (
                -alpha / 2.0
            )
            if normalized > claimed + 3e-12:
                issues.append(
                    f"unconditional high bound failed at "
                    f"alpha={alpha}, K={size}"
                )


def check_conditional_cutoff(issues: list[str]) -> None:
    for alpha in (0.09, 0.47, 0.91):
        epsilon = alpha / 4.0
        eta = (1.0 + alpha) / 8.0
        exponent = (1.0 - alpha) / 2.0 + eta
        if not 2.0 * epsilon < alpha:
            issues.append(f"conditional epsilon guard failed at {alpha}")
        if not 0.0 < exponent < 1.0:
            issues.append(f"conditional cutoff exponent failed at {alpha}")
        for size in (2, 4, 9, 24, 47, 128, 257, 1024):
            cutoff = math.ceil(size**exponent)
            beta = size ** (0.5 + epsilon)
            working_bound = (
                2.0
                * size ** (1.0 - alpha)
                / (2 * cutoff + 1) ** 2
                + 4.0
                * size ** (-1.0 - alpha)
                * beta**2
                / (2 * cutoff + 1)
            )
            calibration = (
                2.0 * size ** (-2.0 * eta)
                + 4.0 * size ** (2.0 * epsilon - alpha)
            )
            if working_bound > calibration + 2e-14:
                issues.append(
                    f"conditional calibration failed at "
                    f"alpha={alpha}, K={size}"
                )


def check_generic_anchor_guard(issues: list[str]) -> None:
    for alpha in (0.09, 0.47, 0.91):
        for size in (1, 2, 4, 9, 24, 47, 128):
            synthetic_beta = 0.0
            for n in range(2 * size, 3 * size + 1):
                weight = anchor_weight(alpha, size, n)
                if weight + 2e-14 < 2.0 / 9.0:
                    issues.append(
                        f"synthetic weight floor failed at "
                        f"alpha={alpha}, K={size}, n={n}"
                    )
                synthetic_beta += weight
            if synthetic_beta + 2e-13 < 2.0 * size / 9.0:
                issues.append(
                    f"generic anchor lower bound failed at "
                    f"alpha={alpha}, K={size}"
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
    check_exact_split(issues)
    check_odd_tail(issues)
    check_high_mode_bounds(issues)
    check_unconditional_cutoff(issues)
    check_conditional_cutoff(issues)
    check_generic_anchor_guard(issues)
    check_formal_core(issues)

    audit = payload.get("audit", {})
    print(
        "validated Mertens spectral anchor/high-mode reduction: "
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
