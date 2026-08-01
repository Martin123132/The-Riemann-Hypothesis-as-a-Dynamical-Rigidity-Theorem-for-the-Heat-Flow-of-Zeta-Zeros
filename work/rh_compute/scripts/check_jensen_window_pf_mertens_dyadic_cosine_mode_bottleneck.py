#!/usr/bin/env python3
"""Validate the dyadic reciprocal-tail cosine-mode bottleneck reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck.md"
)

REQUIRED_IDS = {
    f"mdcos_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "tail_definitions"),
        (2, "dyadic_block"),
        (3, "weight_comparison"),
        (4, "tail_difference"),
        (5, "neumann_cosine_basis"),
        (6, "parseval"),
        (7, "gradient_parseval"),
        (8, "constant_mode"),
        (9, "nonconstant_modes"),
        (10, "high_mode_bound"),
        (11, "summable_threshold"),
        (12, "low_mode_equivalence"),
        (13, "cofinal_rh_criterion"),
        (14, "dimension_reduction"),
        (15, "lowest_oscillatory_mode"),
        (16, "davenport_linear_phase"),
        (17, "davenport_partial_summation"),
        (18, "davenport_power_deficit_guard"),
        (19, "general_exponent_threshold"),
        (20, "square_root_sufficient"),
        (21, "square_root_scope_guard"),
        (22, "translation_average_guard"),
        (23, "full_burnol_guard"),
        (24, "literature_context"),
        (25, "focused_source_search_guard"),
        (26, "open_low_mode_gate"),
    )
}

ROLE_STATUS = {
    "exact_definition": "available_exact",
    "exact_identity": "available_exact",
    "exact_inequality": "available_exact",
    "exact_equivalence": "available_exact",
    "exact_consequence": "available_exact",
    "classical_theorem_step": "source_backed",
    "conditional_consequence": "conditional_exact",
    "literature_guard": "source_backed",
    "proof_guard": "guard_validated",
    "open_target": "open_target",
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Mertens Dyadic Cosine-Mode Bottleneck",
    "y_m^(K):=r_(K+m-1)",
    "D_K^*D_K phi_r=nu_r phi_r",
    "nu_r:=4sin^2(pi*r/(2K))",
    "gradient Parseval",
    "sqrt(K)r_(2K)",
    "K^(1-alpha)/(4R^2)",
    "R_K=ceil(sqrt(K))",
    "automatically summable",
    "one-hump anchored Mobius average",
    "Davenport's theorem",
    "K^(1-alpha)(log K)^(-2A)",
    "sigma<(1+alpha)/2",
    "only a calibration of the missing exponent",
    "not claimed here",
    "RH-equivalent",
    "not a proof",
)


def close(
    left: complex | float,
    right: complex | float,
    tolerance: float = 4e-11,
) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def mobius_values(limit: int) -> list[int]:
    mu = [0] * (limit + 1)
    composite = [False] * (limit + 1)
    primes: list[int] = []
    mu[1] = 1
    for n in range(2, limit + 1):
        if not composite[n]:
            primes.append(n)
            mu[n] = -1
        for prime in primes:
            product = n * prime
            if product > limit:
                break
            composite[product] = True
            if n % prime == 0:
                mu[product] = 0
                break
            mu[product] = -mu[n]
    return mu


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


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def check_schema(payload: dict, note: str, issues: list[str]) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_dyadic_cosine_mode_bottleneck"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    if len(rows) != 26:
        issues.append(f"expected 26 rows, found {len(rows)}")
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
        "row_count": 26,
        "exact_reduction_count": 15,
        "classical_theorem_step_count": 2,
        "conditional_consequence_count": 2,
        "literature_guard_count": 2,
        "proof_guard_count": 4,
        "open_low_mode_gate_count": 1,
        "cosine_diagonalization_proved": True,
        "high_modes_automatically_summable": True,
        "low_mode_cofinal_rh_equivalence_proved": True,
        "davenport_closes_low_modes": False,
        "square_root_twist_bound_proved": False,
        "low_mode_gain_proved": False,
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
        "https://doi.org/10.5802/aif.2401",
    ):
        if required not in anchors:
            issues.append(f"missing source anchor {required}")


def check_cosine_spectrum(issues: list[str]) -> None:
    for size in (1, 2, 4, 8, 17, 32):
        basis = cosine_basis(size)
        for left in range(size):
            for right in range(size):
                expected = 1.0 if left == right else 0.0
                if not close(
                    dot(basis[left], basis[right]),
                    expected,
                    tolerance=2e-12,
                ):
                    issues.append(
                        f"cosine orthogonality failed at "
                        f"K={size}, r={left}, s={right}"
                    )
                    return

        for frequency, vector in enumerate(basis):
            laplacian = [0.0] * size
            if size > 1:
                laplacian[0] = vector[0] - vector[1]
                for index in range(1, size - 1):
                    laplacian[index] = (
                        2.0 * vector[index]
                        - vector[index - 1]
                        - vector[index + 1]
                    )
                laplacian[-1] = vector[-1] - vector[-2]
            eigenvalue = 4.0 * math.sin(
                math.pi * frequency / (2.0 * size)
            ) ** 2
            for index in range(size):
                if not close(
                    laplacian[index],
                    eigenvalue * vector[index],
                    tolerance=3e-12,
                ):
                    issues.append(
                        f"Neumann eigenvector failed at "
                        f"K={size}, r={frequency}, m={index + 1}"
                    )
                    return


def check_parseval_and_modes(issues: list[str]) -> None:
    for size in (2, 4, 8, 16, 31):
        vector = [
            (
                (((37 * (index + 1) + 11) % 41) - 20) / 17.0
                + 0.13 * math.sin((index + 1) * 0.7)
            )
            for index in range(size)
        ]
        basis = cosine_basis(size)
        coefficients = [dot(vector, mode) for mode in basis]
        norm_sq = dot(vector, vector)
        if not close(
            norm_sq,
            sum(value * value for value in coefficients),
            tolerance=3e-12,
        ):
            issues.append(f"cosine Parseval failed at K={size}")

        differences = [
            vector[index] - vector[index + 1]
            for index in range(size - 1)
        ]
        gradient_sq = dot(differences, differences)
        spectral_gradient = sum(
            (
                4.0
                * math.sin(math.pi * frequency / (2.0 * size)) ** 2
                * coefficients[frequency] ** 2
            )
            for frequency in range(1, size)
        )
        if not close(
            gradient_sq,
            spectral_gradient,
            tolerance=4e-12,
        ):
            issues.append(f"gradient Parseval failed at K={size}")

        for frequency in range(1, size):
            sine_sum = sum(
                differences[index]
                * math.sin(
                    math.pi * frequency * (index + 1) / size
                )
                for index in range(size - 1)
            )
            formula = (
                math.sqrt(2.0 / size)
                / (
                    2.0
                    * math.sin(
                        math.pi * frequency / (2.0 * size)
                    )
                )
                * sine_sum
            )
            if not close(
                coefficients[frequency],
                formula,
                tolerance=5e-12,
            ):
                issues.append(
                    f"nonconstant mode formula failed at "
                    f"K={size}, r={frequency}"
                )
                return


def check_mobius_tail_blocks(issues: list[str]) -> None:
    limit = 2048
    mu = mobius_values(limit)
    for alpha in (0.13, 0.47, 0.82):
        q = [0.0] + [
            mu[n] * n ** (-1.0 - alpha)
            for n in range(1, limit + 1)
        ]
        tails = [0.0] * (limit + 1)
        for n in range(limit - 1, -1, -1):
            tails[n] = tails[n + 1] + q[n + 1]

        for size in (2, 4, 8, 16, 32, 64, 128):
            vector = [
                tails[size + index] for index in range(size)
            ]
            basis = cosine_basis(size)
            coefficients = [dot(vector, mode) for mode in basis]

            for index in range(size - 1):
                if not close(
                    vector[index] - vector[index + 1],
                    q[size + index + 1],
                    tolerance=3e-12,
                ):
                    issues.append(
                        f"Mobius tail difference failed at "
                        f"alpha={alpha}, K={size}, m={index + 1}"
                    )
                    return

            constant_formula = (
                math.sqrt(size) * tails[2 * size]
                + sum(
                    ell * q[size + ell]
                    for ell in range(1, size + 1)
                )
                / math.sqrt(size)
            )
            if not close(
                coefficients[0],
                constant_formula,
                tolerance=5e-12,
            ):
                issues.append(
                    f"constant mode failed at alpha={alpha}, K={size}"
                )

            block = sum(
                (size + index) ** alpha * vector[index] ** 2
                for index in range(size)
            )
            unweighted = sum(value * value for value in vector)
            lower = size**alpha * unweighted
            upper = (2.0 * size) ** alpha * unweighted
            if block < lower * (1.0 - 2e-13):
                issues.append(
                    f"block lower comparison failed at alpha={alpha}, "
                    f"K={size}"
                )
            if block > upper * (1.0 + 2e-13):
                issues.append(
                    f"block upper comparison failed at alpha={alpha}, "
                    f"K={size}"
                )

            cutoff = math.ceil(math.sqrt(size))
            high = (
                size**alpha
                * sum(
                    coefficients[frequency] ** 2
                    for frequency in range(cutoff, size)
                )
            )
            high_bound = 0.25 * size ** (-alpha)
            if high > high_bound * (1.0 + 3e-12):
                issues.append(
                    f"universal high-mode bound failed at "
                    f"alpha={alpha}, K={size}: {high}>{high_bound}"
                )

            for cutoff in (1, max(1, size // 4), max(1, size // 2)):
                if cutoff >= size:
                    continue
                high = (
                    size**alpha
                    * sum(
                        coefficients[frequency] ** 2
                        for frequency in range(cutoff, size)
                    )
                )
                high_bound = (
                    size ** (1.0 - alpha)
                    / (4.0 * cutoff * cutoff)
                )
                if high > high_bound * (1.0 + 3e-12):
                    issues.append(
                        f"general high-mode bound failed at "
                        f"alpha={alpha}, K={size}, R={cutoff}"
                    )


def check_threshold_bookkeeping(issues: list[str]) -> None:
    for alpha in (0.08, 0.31, 0.73):
        beta = 0.5
        exponent = 1.0 - alpha - 2.0 * beta
        if not close(exponent, -alpha, tolerance=1e-14):
            issues.append(f"sqrt cutoff exponent failed at alpha={alpha}")

        beta = (1.0 - alpha) / 2.0 + alpha / 4.0
        exponent = 1.0 - alpha - 2.0 * beta
        if not close(exponent, -alpha / 2.0, tolerance=1e-14):
            issues.append(f"general cutoff exponent failed at alpha={alpha}")

        critical_sigma = (1.0 + alpha) / 2.0
        if not close(
            2.0 * critical_sigma - 1.0 - alpha,
            0.0,
            tolerance=1e-14,
        ):
            issues.append(
                f"additive exponent threshold failed at alpha={alpha}"
            )
        epsilon = alpha / 4.0
        square_root_exponent = 2.0 * epsilon - alpha
        if square_root_exponent >= 0.0:
            issues.append(
                f"square-root summability exponent failed at alpha={alpha}"
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
    check_cosine_spectrum(issues)
    check_parseval_and_modes(issues)
    check_mobius_tail_blocks(issues)
    check_threshold_bookkeeping(issues)
    for issue in issues:
        print(f"MERTENS-DYADIC-COSINE {issue}")
    audit = payload.get("audit", {})
    print(
        "validated Mertens dyadic cosine-mode bottleneck: "
        f"{len(payload.get('rows', []))} rows, {len(issues)} issues, "
        f"{audit.get('exact_reduction_count', 0)} exact reductions, "
        f"{audit.get('proof_guard_count', 0)} proof guards, "
        f"{audit.get('open_low_mode_gate_count', 0)} open low-mode gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
