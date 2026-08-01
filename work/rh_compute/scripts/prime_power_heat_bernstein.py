"""Exact Bernstein tools for finite prime-power heat-block currents."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from decimal import Decimal, ROUND_FLOOR, localcontext
from fractions import Fraction
from math import comb, isqrt

import numpy as np


@dataclass(frozen=True)
class FamilyCertificate:
    prime: int
    block_length: int
    q_lower: str
    s_lower: str
    power_degrees: tuple[int, int, int]
    bernstein_shape: tuple[int, int, int]
    x_leaf_paths: tuple[str, ...]
    leaf_minimum_lowers: tuple[str, ...]
    global_minimum_lower: str
    claimed_rational_lower: str
    positive: bool

    def to_dict(self) -> dict:
        return asdict(self)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def decimal_text(value: Fraction, digits: int = 40) -> str:
    with localcontext() as context:
        context.prec = digits
        context.rounding = ROUND_FLOOR
        return format(
            Decimal(value.numerator) / Decimal(value.denominator),
            "e",
        )


def zero_array(shape: tuple[int, ...]) -> np.ndarray:
    array = np.empty(shape, dtype=object)
    array.fill(Fraction(0))
    return array


def chebyshev_power_coefficients(degree: int) -> list[list[int]]:
    values: list[list[int]] = [[1]]
    if degree >= 1:
        values.append([0, 1])
    for _ in range(2, degree + 1):
        doubled_shift = [0] + [2 * value for value in values[-1]]
        previous = values[-2] + [0] * (
            len(doubled_shift) - len(values[-2])
        )
        values.append(
            [
                doubled_shift[index] - previous[index]
                for index in range(len(doubled_shift))
            ]
        )
    return values


def build_current_power_pairs(
    prime: int, block_length: int
) -> tuple[np.ndarray, np.ndarray]:
    """Build J=|Q|^2+Re(H*conj(Q)) in q,s,x power coordinates.

    The coefficient pair (r, u) represents r+u*sqrt(prime).
    """
    exponents = [
        index * (2 * block_length - 2 - index)
        for index in range(block_length)
    ]
    shape = (
        2 * max(exponents) + 1,
        2 * block_length - 1,
        block_length,
    )
    rational = zero_array(shape)
    radical = zero_array(shape)
    chebyshev = chebyshev_power_coefficients(block_length - 1)

    for index in range(block_length):
        rational[2 * exponents[index], 2 * index, 0] += Fraction(
            index + 1, prime**index
        )

    for left in range(block_length):
        for right in range(left + 1, block_length):
            total_index = left + right
            q_degree = exponents[left] + exponents[right]
            multiplier = left + right + 2
            for x_degree, coefficient in enumerate(
                chebyshev[right - left]
            ):
                numerator = multiplier * coefficient
                if total_index % 2 == 0:
                    rational[q_degree, total_index, x_degree] += (
                        Fraction(
                            numerator,
                            prime ** (total_index // 2),
                        )
                    )
                else:
                    radical[q_degree, total_index, x_degree] += (
                        Fraction(
                            numerator,
                            prime ** ((total_index + 1) // 2),
                        )
                    )
    return rational, radical


def affine_power_axis(
    coefficients: np.ndarray,
    axis: int,
    lower: Fraction,
    upper: Fraction,
) -> np.ndarray:
    source = np.moveaxis(coefficients, axis, 0)
    degree = source.shape[0] - 1
    target = zero_array(source.shape)
    width = upper - lower
    for old_degree in range(degree + 1):
        if not any(source[old_degree].flat):
            continue
        for new_degree in range(old_degree + 1):
            factor = (
                Fraction(comb(old_degree, new_degree))
                * lower ** (old_degree - new_degree)
                * width**new_degree
            )
            target[new_degree] += source[old_degree] * factor
    return np.moveaxis(target, 0, axis)


def power_to_bernstein_axis(
    coefficients: np.ndarray, axis: int
) -> np.ndarray:
    source = np.moveaxis(coefficients, axis, 0)
    degree = source.shape[0] - 1
    target = zero_array(source.shape)
    for bernstein_index in range(degree + 1):
        for power_index in range(bernstein_index + 1):
            target[bernstein_index] += source[power_index] * Fraction(
                comb(bernstein_index, power_index),
                comb(degree, power_index),
            )
    return np.moveaxis(target, 0, axis)


def split_bernstein_axis(
    coefficients: np.ndarray, axis: int
) -> tuple[np.ndarray, np.ndarray]:
    source = np.moveaxis(coefficients, axis, 0)
    degree = source.shape[0] - 1
    work = source.copy()
    left = zero_array(source.shape)
    right = zero_array(source.shape)
    left[0] = work[0]
    right[degree] = work[degree]
    for level in range(1, degree + 1):
        work = (work[:-1] + work[1:]) / 2
        left[level] = work[0]
        right[degree - level] = work[-1]
    return np.moveaxis(left, 0, axis), np.moveaxis(right, 0, axis)


def sqrt_bounds(
    radicand: int, digits: int = 90
) -> tuple[Fraction, Fraction]:
    scale = 10**digits
    lower_integer = isqrt(radicand * scale * scale)
    return (
        Fraction(lower_integer, scale),
        Fraction(lower_integer + 1, scale),
    )


def coefficient_lower(
    rational: Fraction,
    radical: Fraction,
    sqrt_lower: Fraction,
    sqrt_upper: Fraction,
) -> Fraction:
    root = sqrt_lower if radical >= 0 else sqrt_upper
    return rational + radical * root


def minimum_coefficient_lower(
    rational: np.ndarray, radical: np.ndarray, prime: int
) -> Fraction:
    sqrt_lower, sqrt_upper = sqrt_bounds(prime)
    minimum: Fraction | None = None
    for rational_value, radical_value in zip(
        rational.flat, radical.flat
    ):
        lower = coefficient_lower(
            rational_value,
            radical_value,
            sqrt_lower,
            sqrt_upper,
        )
        if minimum is None or lower < minimum:
            minimum = lower
    if minimum is None:
        raise RuntimeError("empty Bernstein coefficient array")
    return minimum


def follow_leaf_path(
    rational: np.ndarray,
    radical: np.ndarray,
    path: str,
    axis: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    current_rational = rational
    current_radical = radical
    for direction in path:
        rational_left, rational_right = split_bernstein_axis(
            current_rational, axis
        )
        radical_left, radical_right = split_bernstein_axis(
            current_radical, axis
        )
        if direction == "L":
            current_rational, current_radical = (
                rational_left,
                radical_left,
            )
        elif direction == "R":
            current_rational, current_radical = (
                rational_right,
                radical_right,
            )
        else:
            raise ValueError(f"bad leaf direction: {direction!r}")
    return current_rational, current_radical


def leaf_minimum_lowers(
    rational: np.ndarray,
    radical: np.ndarray,
    prime: int,
    paths: tuple[str, ...],
    axis: int = 2,
) -> dict[str, Fraction]:
    requested = set(paths)
    minima: dict[str, Fraction] = {}

    def visit(
        current_rational: np.ndarray,
        current_radical: np.ndarray,
        prefix: str,
    ) -> None:
        if prefix in requested:
            minima[prefix] = minimum_coefficient_lower(
                current_rational, current_radical, prime
            )
            return
        descendants = [
            path for path in paths if path.startswith(prefix)
        ]
        if not descendants:
            return
        rational_left, rational_right = split_bernstein_axis(
            current_rational, axis
        )
        radical_left, radical_right = split_bernstein_axis(
            current_radical, axis
        )
        visit(rational_left, radical_left, prefix + "L")
        visit(rational_right, radical_right, prefix + "R")

    visit(rational, radical, "")
    if set(minima) != requested:
        missing = sorted(requested - set(minima))
        raise RuntimeError(f"unreached Bernstein leaves: {missing}")
    return minima


def certify_family(
    *,
    prime: int,
    block_length: int,
    q_lower: Fraction,
    s_lower: Fraction,
    x_leaf_paths: tuple[str, ...],
    claimed_lower: Fraction,
) -> FamilyCertificate:
    rational, radical = build_current_power_pairs(prime, block_length)
    power_degrees = tuple(size - 1 for size in rational.shape)
    boxes = (
        (q_lower, Fraction(1)),
        (s_lower, Fraction(1)),
        (Fraction(-1), Fraction(1)),
    )
    for axis, (lower, upper) in enumerate(boxes):
        rational = affine_power_axis(rational, axis, lower, upper)
        radical = affine_power_axis(radical, axis, lower, upper)
    for axis in range(3):
        rational = power_to_bernstein_axis(rational, axis)
        radical = power_to_bernstein_axis(radical, axis)

    minimum_map = leaf_minimum_lowers(
        rational, radical, prime, x_leaf_paths
    )
    minima = [minimum_map[path] for path in x_leaf_paths]
    global_minimum = min(minima)
    return FamilyCertificate(
        prime=prime,
        block_length=block_length,
        q_lower=fraction_text(q_lower),
        s_lower=fraction_text(s_lower),
        power_degrees=power_degrees,
        bernstein_shape=tuple(int(size) for size in rational.shape),
        x_leaf_paths=x_leaf_paths,
        leaf_minimum_lowers=tuple(decimal_text(value) for value in minima),
        global_minimum_lower=decimal_text(global_minimum),
        claimed_rational_lower=fraction_text(claimed_lower),
        positive=global_minimum > claimed_lower > 0,
    )


FAMILY_SPECS = (
    {
        "prime": 2,
        "block_length": 8,
        "q_lower": Fraction(47, 50),
        "s_lower": Fraction(22, 25),
        "x_leaf_paths": ("LLL", "LLR", "LR", "RL", "RRL", "RRR"),
        "claimed_lower": Fraction(1, 300),
    },
    {
        "prime": 3,
        "block_length": 4,
        "q_lower": Fraction(17, 20),
        "s_lower": Fraction(73, 100),
        "x_leaf_paths": ("L", "R"),
        "claimed_lower": Fraction(1, 25),
    },
    {
        "prime": 5,
        "block_length": 2,
        "q_lower": Fraction(18, 25),
        "s_lower": Fraction(13, 25),
        "x_leaf_paths": ("",),
        "claimed_lower": Fraction(1, 20),
    },
)


def build_base_certificates() -> list[FamilyCertificate]:
    return [certify_family(**specification) for specification in FAMILY_SPECS]
