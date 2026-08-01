#!/usr/bin/env python3
"""Validate the centered Mertens-bridge and Vaughan handoff reduction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RESULT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_mertens_centered_bridge_vaughan_handoff.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_mertens_centered_bridge_vaughan_handoff.md"
)

REQUIRED_IDS = {
    f"mcbvh_{index:02d}_{suffix}"
    for index, suffix in (
        (1, "tail_setup"),
        (2, "dyadic_dct"),
        (3, "common_kernels"),
        (4, "kernel_identity"),
        (5, "analysis_operator"),
        (6, "cosine_energy_equivalence"),
        (7, "mean_row_bound"),
        (8, "nonconstant_row_bound"),
        (9, "finite_trace"),
        (10, "signed_gram_expansion"),
        (11, "offdiagonal_criterion"),
        (12, "trace_class_interpretation"),
        (13, "linf_operator_guard"),
        (14, "local_mertens_path"),
        (15, "quotient_norm"),
        (16, "extended_abel_operator"),
        (17, "constant_preservation"),
        (18, "abel_operator_upper"),
        (19, "abel_inverse"),
        (20, "abel_inverse_upper"),
        (21, "centered_abel_equivalence"),
        (22, "centered_dct_identity"),
        (23, "high_modes_summable"),
        (24, "centered_bridge_equivalence"),
        (25, "two_component_criterion"),
        (26, "cofinal_rh_criterion"),
        (27, "pairwise_variance"),
        (28, "brownian_bridge_kernel"),
        (29, "bridge_diagonal_offdiagonal"),
        (30, "bridge_diagonal_bound"),
        (31, "bridge_offdiagonal_criterion"),
        (32, "vaughan_source"),
        (33, "shifted_bridge_test"),
        (34, "vaughan_coefficients"),
        (35, "exact_vaughan_handoff"),
        (36, "signed_type_i_ii_criterion"),
        (37, "separate_bound_guard"),
        (38, "averaged_chowla_guard"),
        (39, "mean_and_scope_guard"),
        (40, "open_two_component_gate"),
    )
}

ROLE_STATUS = {
    "exact_definition": "available_exact",
    "exact_identity": "available_exact",
    "exact_inequality": "available_exact",
    "exact_equivalence": "available_exact",
    "exact_consequence": "available_exact",
    "literature_guard": "source_backed",
    "proof_guard": "guard_validated",
    "open_target": "open_target",
}

REQUIRED_NOTE = (
    "# Jensen-Window PF Mertens Centered-Bridge Vaughan Handoff",
    "h_(K,0)(n)",
    "c_(K,r)",
    "Hilbert-Schmidt",
    "purely signed and off-diagonal",
    "A_K(c1)=a_1c1",
    "Centered Abel Quotient",
    "[2^alpha(3+alpha)]^(-1)",
    "No commutator estimate is needed",
    "Discrete Brownian-Bridge Kernel",
    "B_K(i,j)",
    "min(i,j)-ij/K",
    "Delta_K+2C_K",
    "(K^2-1)/6",
    "Pre-Collapse Vaughan Handoff",
    "C_K=-TI_K+TII_K",
    "|TI_K|+|TII_K|",
    "one full power",
    "two noninterchangeable pieces",
    "not a proof",
)


def close(
    left: complex | float,
    right: complex | float,
    tolerance: float = 6e-11,
) -> bool:
    scale = max(1.0, abs(left), abs(right))
    return abs(left - right) <= tolerance * scale


def dot(left: list[float], right: list[float]) -> float:
    return sum(a * b for a, b in zip(left, right))


def norm(vector: list[float]) -> float:
    return math.sqrt(dot(vector, vector))


def centered_energy(vector: list[float]) -> float:
    mean = sum(vector) / len(vector)
    return sum((value - mean) ** 2 for value in vector)


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


def kernel_h(size: int, frequency: int, n: int) -> float:
    if frequency == 0:
        return min(max(n - size, 0), size) / math.sqrt(size)
    offset = n - size
    if not (1 <= offset < size):
        return 0.0
    return (
        math.sqrt(2.0 / size)
        * math.sin(math.pi * frequency * offset / size)
        / (
            2.0
            * math.sin(math.pi * frequency / (2.0 * size))
        )
    )


def check_schema(payload: dict, note: str, issues: list[str]) -> None:
    if payload.get("kind") != (
        "jensen_window_pf_mertens_centered_bridge_vaughan_handoff"
    ):
        issues.append("bad kind")
    rows = payload.get("rows", [])
    if len(rows) != 40:
        issues.append(f"expected 40 rows, found {len(rows)}")
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
        "row_count": 40,
        "exact_reduction_count": 34,
        "literature_guard_count": 2,
        "proof_guard_count": 3,
        "open_two_component_gate_count": 1,
        "global_gram_diagonal_summable": True,
        "analysis_operator_hilbert_schmidt_on_l2": True,
        "analysis_operator_bounded_linf_to_l2": False,
        "centered_abel_quotient_equivalence_proved": True,
        "brownian_bridge_kernel_proved": True,
        "bridge_vaughan_handoff_proved": True,
        "averaged_chowla_closes_bridge": False,
        "signed_type_i_ii_gain_proved": False,
        "affine_mean_series_proved": False,
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
        "https://doi.org/10.5802/aif.2401",
        "https://arxiv.org/abs/1503.05121",
        "https://doi.org/10.2140/ant.2015.9.2167",
    ):
        if required not in anchors:
            issues.append(f"missing source anchor {required}")


def check_common_kernels_and_trace(issues: list[str]) -> None:
    limit = 1024
    mu = mobius_values(limit)
    for alpha in (0.09, 0.37, 0.81):
        q = [0.0] + [
            mu[n] * n ** (-1.0 - alpha)
            for n in range(1, limit + 1)
        ]
        tails = [0.0] * (limit + 1)
        for n in range(limit - 1, -1, -1):
            tails[n] = tails[n + 1] + q[n + 1]

        trace_partial = 0.0
        for size in (2, 4, 8, 16, 32, 64):
            basis = cosine_basis(size)
            tail_vector = [
                tails[size + index] for index in range(size)
            ]
            dct = [dot(tail_vector, mode) for mode in basis]
            for frequency in range(size):
                formula = sum(
                    q[n] * kernel_h(size, frequency, n)
                    for n in range(1, limit + 1)
                )
                if not close(
                    dct[frequency],
                    formula,
                    tolerance=7e-12,
                ):
                    issues.append(
                        f"common kernel failed at alpha={alpha}, "
                        f"K={size}, r={frequency}"
                    )
                    return

                row_norm = sum(
                    (
                        n ** (-1.0 - alpha)
                        * kernel_h(size, frequency, n)
                    )
                    ** 2
                    for n in range(1, limit + 1)
                )
                trace_partial += size**alpha * row_norm
                if frequency == 0:
                    bound = (
                        1.0
                        + 2.0 ** (-1.0 - 2.0 * alpha)
                        / (1.0 + 2.0 * alpha)
                    ) * size ** (-2.0 * alpha)
                else:
                    bound = (
                        size ** (-2.0 * alpha)
                        / (4.0 * frequency * frequency)
                    )
                if row_norm > bound * (1.0 + 4e-12):
                    issues.append(
                        f"row norm bound failed at alpha={alpha}, "
                        f"K={size}, r={frequency}"
                    )

        constant = (
            1.0
            + 2.0 ** (-1.0 - 2.0 * alpha)
            / (1.0 + 2.0 * alpha)
            + math.pi * math.pi / 24.0
        )
        geometric_partial = sum(
            size ** (-alpha)
            for size in (2, 4, 8, 16, 32, 64)
        )
        if trace_partial > (
            constant * geometric_partial * (1.0 + 5e-12)
        ):
            issues.append(f"trace bound failed at alpha={alpha}")

        sizes = (2, 4, 8, 16)
        coefficient_rows: list[tuple[float, list[float]]] = []
        energy = 0.0
        for size in sizes:
            basis = cosine_basis(size)
            tail_vector = [
                tails[size + index] for index in range(size)
            ]
            for frequency, mode in enumerate(basis):
                coefficient = dot(tail_vector, mode)
                energy += size**alpha * coefficient * coefficient
                row_vector = [
                    (
                        n ** (-1.0 - alpha)
                        * kernel_h(size, frequency, n)
                    )
                    for n in range(1, limit + 1)
                ]
                coefficient_rows.append((size**alpha, row_vector))

        diagonal = 0.0
        offdiagonal = 0.0
        for weight, row_vector in coefficient_rows:
            diagonal += weight * sum(
                mu[n] * mu[n] * row_vector[n - 1] ** 2
                for n in range(1, limit + 1)
            )
            offdiagonal += weight * sum(
                mu[m]
                * mu[n]
                * row_vector[m - 1]
                * row_vector[n - 1]
                for m in range(1, limit)
                for n in range(m + 1, limit + 1)
            )
        if not close(
            energy,
            diagonal + 2.0 * offdiagonal,
            tolerance=2e-10,
        ):
            issues.append(f"signed Gram expansion failed at alpha={alpha}")


def apply_abel(
    vector: list[float],
    size: int,
    alpha: float,
) -> list[float]:
    weights = [
        0.0
    ] + [
        (size + index) ** (-1.0 - alpha)
        for index in range(1, size)
    ]
    result = [weights[1] * vector[0]]
    running = result[0]
    for index in range(1, size):
        running += weights[index] * (
            vector[index] - vector[index - 1]
        )
        result.append(running)
    return result


def apply_abel_inverse(
    vector: list[float],
    size: int,
    alpha: float,
) -> list[float]:
    inverse_weights = [
        0.0
    ] + [
        (size + index) ** (1.0 + alpha)
        for index in range(1, size)
    ]
    result = [inverse_weights[1] * vector[0]]
    running = result[0]
    for index in range(1, size):
        running += inverse_weights[index] * (
            vector[index] - vector[index - 1]
        )
        result.append(running)
    return result


def check_centered_abel(issues: list[str]) -> None:
    for alpha in (0.07, 0.42, 0.88):
        for size in (2, 3, 4, 8, 17, 32, 65):
            vector = [
                (
                    (((37 * (index + 1) + 9) % 43) - 21) / 13.0
                    + 0.11 * math.sin(0.7 * index)
                )
                for index in range(size)
            ]
            image = apply_abel(vector, size, alpha)
            recovered = apply_abel_inverse(image, size, alpha)
            for index in range(size):
                if not close(
                    vector[index],
                    recovered[index],
                    tolerance=1e-10,
                ):
                    issues.append(
                        f"extended Abel inverse failed at "
                        f"alpha={alpha}, K={size}, m={index}"
                    )
                    return

            constant = [0.31] * size
            constant_image = apply_abel(constant, size, alpha)
            expected_constant = (
                (size + 1) ** (-1.0 - alpha) * 0.31
            )
            if any(
                not close(value, expected_constant, tolerance=3e-12)
                for value in constant_image
            ):
                issues.append(
                    f"constant preservation failed at alpha={alpha}, "
                    f"K={size}"
                )

            image_norm = norm(image)
            upper = (
                (2.0 + alpha)
                * size ** (-1.0 - alpha)
                * norm(vector)
            )
            if image_norm > upper * (1.0 + 3e-12):
                issues.append(
                    f"Abel operator upper failed at alpha={alpha}, "
                    f"K={size}"
                )
            inverse_upper = (
                2.0**alpha
                * (3.0 + alpha)
                * size ** (1.0 + alpha)
                * image_norm
            )
            if norm(vector) > inverse_upper * (1.0 + 3e-12):
                issues.append(
                    f"Abel inverse upper failed at alpha={alpha}, "
                    f"K={size}"
                )

            source_centered = centered_energy(vector)
            image_centered = centered_energy(image)
            lower_centered = (
                size ** (-2.0 - 2.0 * alpha)
                * source_centered
                / (
                    2.0 ** (2.0 * alpha)
                    * (3.0 + alpha) ** 2
                )
            )
            upper_centered = (
                (2.0 + alpha) ** 2
                * size ** (-2.0 - 2.0 * alpha)
                * source_centered
            )
            if image_centered < lower_centered * (1.0 - 5e-12):
                issues.append(
                    f"centered Abel lower failed at alpha={alpha}, "
                    f"K={size}"
                )
            if image_centered > upper_centered * (1.0 + 5e-12):
                issues.append(
                    f"centered Abel upper failed at alpha={alpha}, "
                    f"K={size}"
                )

            basis = cosine_basis(size)
            coefficients = [dot(image, mode) for mode in basis]
            spectral_centered = sum(
                value * value for value in coefficients[1:]
            )
            if not close(
                image_centered,
                spectral_centered,
                tolerance=5e-12,
            ):
                issues.append(
                    f"centered DCT identity failed at alpha={alpha}, "
                    f"K={size}"
                )


def check_brownian_bridge_and_vaughan(issues: list[str]) -> None:
    limit = 1024
    mu = mobius_values(limit)
    for size in (4, 8, 16, 31, 64, 96):
        increments = [mu[size + index] for index in range(1, size)]
        path = [0.0]
        for value in increments:
            path.append(path[-1] + value)
        variance = centered_energy(path)
        pairwise = sum(
            (path[right] - path[left]) ** 2
            for left in range(size)
            for right in range(left + 1, size)
        ) / size
        if not close(variance, pairwise, tolerance=4e-12):
            issues.append(f"pairwise variance failed at K={size}")

        bridge = sum(
            mu[size + left]
            * mu[size + right]
            * (min(left, right) - left * right / size)
            for left in range(1, size)
            for right in range(1, size)
        )
        if not close(variance, bridge, tolerance=5e-12):
            issues.append(f"Brownian bridge kernel failed at K={size}")

        diagonal = sum(
            mu[size + index] ** 2
            * index
            * (size - index)
            / size
            for index in range(1, size)
        )
        offdiagonal = sum(
            mu[size + index]
            * mu[size + index + shift]
            * index
            * (size - index - shift)
            / size
            for shift in range(1, size - 1)
            for index in range(1, size - shift)
        )
        if not close(
            variance,
            diagonal + 2.0 * offdiagonal,
            tolerance=5e-12,
        ):
            issues.append(
                f"bridge diagonal/offdiagonal failed at K={size}"
            )
        diagonal_bound = (size * size - 1.0) / 6.0
        if diagonal > diagonal_bound * (1.0 + 3e-12):
            issues.append(f"bridge diagonal bound failed at K={size}")

        weight_total = sum(
            index * (size - index - shift) / size
            for shift in range(1, size - 1)
            for index in range(1, size - shift)
        )
        if weight_total > size**3 / 8.0 * (1.0 + 3e-12):
            issues.append(f"absolute bridge scale failed at K={size}")

        u_cut = max(1, math.isqrt(size))
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
        b_coefficients = {
            d: sum(
                mu[c]
                for c in range(v_cut + 1, d + 1)
                if d % c == 0
            )
            for d in range(
                v_cut + 1,
                (2 * size) // u_cut + 1,
            )
        }
        type_i = 0.0
        type_ii = 0.0
        for shift in range(1, size - 1):
            def test(n: int) -> float:
                index = n - size
                if not (1 <= index <= size - 1 - shift):
                    return 0.0
                return (
                    mu[n + shift]
                    * index
                    * (size - index - shift)
                    / size
                )

            for d, coefficient in a_coefficients.items():
                type_i += coefficient * sum(
                    test(d * w)
                    for w in range(
                        size // d + 1,
                        (2 * size) // d + 1,
                    )
                )
            for d, coefficient in b_coefficients.items():
                type_ii += coefficient * sum(
                    mu[w] * test(d * w)
                    for w in range(
                        max(u_cut, size // d) + 1,
                        (2 * size) // d + 1,
                    )
                )
        if not close(
            offdiagonal,
            -type_i + type_ii,
            tolerance=2e-10,
        ):
            issues.append(
                f"bridge Vaughan handoff failed at K={size}: "
                f"C={offdiagonal}, rhs={-type_i + type_ii}"
            )


def check_countermodel_and_exponents(issues: list[str]) -> None:
    for alpha in (0.11, 0.47, 0.83):
        for size in (4, 8, 32, 128):
            amplitude = (
                2.0 ** (-2.0 - alpha)
                / math.pi
                * size ** (-alpha)
            )
            path = [
                amplitude
                * math.sin(2.0 * math.pi * index / size)
                for index in range(size + 1)
            ]
            q = [
                path[index] - path[index - 1]
                for index in range(1, size + 1)
            ]
            bounded_input = max(
                abs(q[index - 1])
                * (size + index) ** (1.0 + alpha)
                for index in range(1, size + 1)
            )
            if bounded_input > 1.0 + 3e-12:
                issues.append(
                    f"linf countermodel input failed at "
                    f"alpha={alpha}, K={size}"
                )
            first_mode = abs(
                -dot(path[:-1], cosine_basis(size)[1])
            )
            weighted_lower = (
                2.0 ** (-4.0 - 2.0 * alpha)
                * math.pi ** (-4.0)
                * size ** (1.0 - alpha)
            )
            if size**alpha * first_mode**2 < (
                weighted_lower * (1.0 - 4e-12)
            ):
                issues.append(
                    f"linf countermodel output failed at "
                    f"alpha={alpha}, K={size}"
                )

        epsilon = alpha / 2.0
        exponent = epsilon - alpha
        if exponent >= 0.0:
            issues.append(
                f"bridge calibration exponent failed at alpha={alpha}"
            )
        absolute_exponent = 1.0 - alpha
        if absolute_exponent <= 0.0:
            issues.append(
                f"absolute bridge barrier failed at alpha={alpha}"
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
    check_common_kernels_and_trace(issues)
    check_centered_abel(issues)
    check_brownian_bridge_and_vaughan(issues)
    check_countermodel_and_exponents(issues)
    for issue in issues:
        print(f"MERTENS-CENTERED-BRIDGE {issue}")
    audit = payload.get("audit", {})
    print(
        "validated Mertens centered-bridge Vaughan handoff: "
        f"{len(payload.get('rows', []))} rows, {len(issues)} issues, "
        f"{audit.get('exact_reduction_count', 0)} exact reductions, "
        f"{audit.get('proof_guard_count', 0)} proof guards, "
        f"{audit.get('open_two_component_gate_count', 0)} open gate"
    )
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
