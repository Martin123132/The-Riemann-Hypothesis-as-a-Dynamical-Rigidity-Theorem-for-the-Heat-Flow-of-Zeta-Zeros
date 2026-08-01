"""Exact Bernstein certificates for the dyadic heat-starlikeness core."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction

import numpy as np

from prime_power_heat_bernstein import (
    affine_power_axis,
    decimal_text,
    leaf_minimum_lowers,
    power_to_bernstein_axis,
    zero_array,
)


Q2_LOWER = Fraction(47, 50) ** 2
RATIO_UPPER = Fraction(707107, 1_000_000)


@dataclass(frozen=True)
class CoreCertificate:
    certificate_kind: str
    block_length: int | None
    angular_cutoff: int
    group_cutoff: int | None
    ratio_anchor: int | None
    q2_interval: tuple[str, str]
    ratio_interval: tuple[str, str]
    x_interval: tuple[str, str]
    power_degrees: tuple[int, int, int]
    bernstein_shape: tuple[int, int, int]
    x_leaf_paths: tuple[str, ...]
    leaf_minimum_lowers: tuple[str, ...]
    global_minimum_lower: str
    global_minimum_fraction: str
    claimed_rational_lower: str
    positive: bool

    def to_dict(self) -> dict:
        return asdict(self)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def polynomial_add(left: list[int], right: list[int]) -> list[int]:
    size = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else 0)
        + (right[index] if index < len(right) else 0)
        for index in range(size)
    ]


def polynomial_scale(values: list[int], scale: int) -> list[int]:
    return [scale * value for value in values]


def chebyshev_polynomials(degree: int) -> list[list[int]]:
    values = [[1]]
    if degree:
        values.append([0, 1])
    for _ in range(2, degree + 1):
        shifted = [0] + polynomial_scale(values[-1], 2)
        values.append(
            polynomial_add(shifted, polynomial_scale(values[-2], -1))
        )
    return values


def fejer_numerator_polynomials(degree: int) -> list[list[int]]:
    chebyshev = chebyshev_polynomials(degree)
    values = []
    for d in range(degree + 1):
        polynomial = [d + 1]
        for frequency in range(1, d + 1):
            polynomial = polynomial_add(
                polynomial,
                polynomial_scale(
                    chebyshev[frequency],
                    2 * (d + 1 - frequency),
                ),
            )
        values.append(polynomial)
    return values


def add_term(
    terms: dict[tuple[int, int, int], Fraction],
    key: tuple[int, int, int],
    coefficient: int | Fraction,
) -> None:
    terms[key] = terms.get(key, Fraction()) + Fraction(coefficient)
    if not terms[key]:
        terms.pop(key)


def terms_to_array(
    terms: dict[tuple[int, int, int], Fraction]
) -> np.ndarray:
    shape = tuple(
        max(key[axis] for key in terms) + 1 for axis in range(3)
    )
    array = zero_array(shape)
    for key, coefficient in terms.items():
        array[key] += coefficient
    return array


def tail_ratio_q2_exponents(M: int) -> list[int]:
    values = []
    for index in range(M):
        exponent = index * (2 * M - 3 - index)
        if exponent % 2:
            raise RuntimeError("tail-ratio q exponent is not even")
        values.append(exponent // 2)
    return values


def finite_core_power(
    M: int,
    angular_cutoff: int = 7,
    claimed_lower: Fraction = Fraction(1, 10),
) -> np.ndarray:
    if M < angular_cutoff + 3:
        raise ValueError("block is too short for the terminal pair")
    angular = fejer_numerator_polynomials(angular_cutoff)
    exponents = tail_ratio_q2_exponents(M)
    terms: dict[tuple[int, int, int], Fraction] = {}
    for d in range(angular_cutoff + 1):
        for increment, scale in ((0, 1), (1, -2), (2, 1)):
            offset = d + increment
            for k in range(M - offset):
                scalar = scale * (2 * k + offset + 2)
                q2_power = exponents[k] + exponents[k + offset]
                ratio_power = 2 * k + offset
                for x_power, angular_coefficient in enumerate(angular[d]):
                    if angular_coefficient:
                        add_term(
                            terms,
                            (q2_power, ratio_power, x_power),
                            scalar * angular_coefficient,
                        )
    add_term(terms, (0, 0, 0), -claimed_lower)
    return terms_to_array(terms)


def anchor_q2_exponent(anchor: int, index: int) -> int:
    exponent = index * (2 * anchor - index + 1)
    if exponent % 2:
        raise RuntimeError("anchor q exponent is not even")
    return exponent // 2


def initial_core_power(
    group_cutoff: int = 7,
    angular_cutoff: int = 4,
    claimed_lower: Fraction = Fraction(1, 100),
) -> tuple[np.ndarray, int]:
    anchor = group_cutoff + angular_cutoff + 1
    angular = fejer_numerator_polynomials(angular_cutoff)
    terms: dict[tuple[int, int, int], Fraction] = {}
    for d in range(angular_cutoff + 1):
        for k in range(group_cutoff + 1):
            A = 2 * k + d + 2
            for offset, scalar in (
                (d, A),
                (d + 1, -2 * (A + 1)),
                (d + 2, A + 2),
            ):
                q2_power = (
                    anchor_q2_exponent(anchor, k)
                    + anchor_q2_exponent(anchor, k + offset)
                )
                ratio_power = 2 * k + offset
                for x_power, angular_coefficient in enumerate(angular[d]):
                    if angular_coefficient:
                        add_term(
                            terms,
                            (q2_power, ratio_power, x_power),
                            scalar * angular_coefficient,
                        )
    add_term(terms, (0, 0, 0), -claimed_lower)
    return terms_to_array(terms), anchor


def certify_power_array(
    power: np.ndarray,
    *,
    certificate_kind: str,
    block_length: int | None,
    angular_cutoff: int,
    group_cutoff: int | None,
    ratio_anchor: int | None,
    ratio_lower: Fraction,
    claimed_lower: Fraction,
    x_leaf_paths: tuple[str, ...],
) -> CoreCertificate:
    power_degrees = tuple(size - 1 for size in power.shape)
    intervals = (
        (Q2_LOWER, Fraction(1)),
        (ratio_lower, RATIO_UPPER),
        (Fraction(-1), Fraction(1)),
    )
    bernstein = power
    for axis, interval in enumerate(intervals):
        bernstein = affine_power_axis(bernstein, axis, *interval)
    for axis in range(3):
        bernstein = power_to_bernstein_axis(bernstein, axis)
    zero = zero_array(bernstein.shape)
    minimum_map = leaf_minimum_lowers(
        bernstein,
        zero,
        2,
        x_leaf_paths,
        axis=2,
    )
    minima = [minimum_map[path] for path in x_leaf_paths]
    minimum = min(minima)
    return CoreCertificate(
        certificate_kind=certificate_kind,
        block_length=block_length,
        angular_cutoff=angular_cutoff,
        group_cutoff=group_cutoff,
        ratio_anchor=ratio_anchor,
        q2_interval=(fraction_text(Q2_LOWER), "1/1"),
        ratio_interval=(
            fraction_text(ratio_lower),
            fraction_text(RATIO_UPPER),
        ),
        x_interval=("-1/1", "1/1"),
        power_degrees=power_degrees,
        bernstein_shape=tuple(size for size in bernstein.shape),
        x_leaf_paths=x_leaf_paths,
        leaf_minimum_lowers=tuple(
            decimal_text(value) for value in minima
        ),
        global_minimum_lower=decimal_text(minimum),
        global_minimum_fraction=fraction_text(minimum),
        claimed_rational_lower=fraction_text(claimed_lower),
        positive=minimum > 0,
    )


def certify_finite_core(
    M: int,
    *,
    x_leaf_paths: tuple[str, ...],
) -> CoreCertificate:
    target = Fraction(1, 10)
    return certify_power_array(
        finite_core_power(M, 7, target),
        certificate_kind="finite_seven_kernel_core",
        block_length=M,
        angular_cutoff=7,
        group_cutoff=None,
        ratio_anchor=M - 2,
        ratio_lower=Fraction(0),
        claimed_lower=target,
        x_leaf_paths=x_leaf_paths,
    )


def certify_initial_core(
    *,
    x_leaf_paths: tuple[str, ...],
) -> CoreCertificate:
    target = Fraction(1, 100)
    power, anchor = initial_core_power(7, 4, target)
    return certify_power_array(
        power,
        certificate_kind="all_length_initial_four_kernel_core",
        block_length=None,
        angular_cutoff=4,
        group_cutoff=7,
        ratio_anchor=anchor,
        ratio_lower=Fraction(0),
        claimed_lower=target,
        x_leaf_paths=x_leaf_paths,
    )
