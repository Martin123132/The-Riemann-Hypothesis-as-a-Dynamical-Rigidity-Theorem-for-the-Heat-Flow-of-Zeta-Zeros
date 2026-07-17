#!/usr/bin/env python3
"""Rigorous H9-H16 point jets from exact half-grid H8 samples and H24 walls."""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_localized_final_gap_interval_core as order10  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_rational,
    arb_upper_text,
    upper_absolute,
)


BASE_H_ORDER = 8
MAXIMUM_H_ORDER = 16
NODE_COUNT = 16
TAYLOR_DEGREE = NODE_COUNT - 1
REMAINDER_H_ORDER = BASE_H_ORDER + TAYLOR_DEGREE + 1


def _invert_fraction_matrix(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    size = len(matrix)
    if size == 0 or any(len(row) != size for row in matrix):
        raise ValueError("finite-difference matrix must be nonempty and square")
    augmented = [
        list(row)
        + [Fraction(int(column == index)) for column in range(size)]
        for index, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if augmented[row][column]),
            None,
        )
        if pivot is None:
            raise RuntimeError("singular finite-difference Vandermonde matrix")
        if pivot != column:
            augmented[column], augmented[pivot] = (
                augmented[pivot],
                augmented[column],
            )
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    value - factor * pivot_value
                    for value, pivot_value in zip(
                        augmented[row], augmented[column]
                    )
                ]
    return [row[size:] for row in augmented]


@lru_cache(maxsize=64)
def stencil_weights(
    offsets: tuple[Fraction, ...],
) -> tuple[tuple[Fraction, ...], ...]:
    """Return exact weights for derivatives one through eight."""
    if len(offsets) != NODE_COUNT or len(set(offsets)) != NODE_COUNT:
        raise ValueError("order-eleven stencil requires 16 distinct nodes")
    vandermonde = [
        [offset**power for offset in offsets]
        for power in range(NODE_COUNT)
    ]
    inverse = _invert_fraction_matrix(vandermonde)
    rows = []
    for derivative in range(1, MAXIMUM_H_ORDER - BASE_H_ORDER + 1):
        weights = tuple(
            Fraction(math.factorial(derivative))
            * inverse[node][derivative]
            for node in range(NODE_COUNT)
        )
        for power in range(NODE_COUNT):
            moment = sum(
                weight * offset**power
                for weight, offset in zip(weights, offsets)
            )
            expected = (
                Fraction(math.factorial(derivative))
                if power == derivative
                else Fraction(0)
            )
            if moment != expected:
                raise RuntimeError(
                    f"bad exact stencil moment d={derivative}, p={power}"
                )
        rows.append(weights)
    return tuple(rows)


def stencil_moment_contract(offsets: tuple[Fraction, ...]) -> dict:
    weights = stencil_weights(offsets)
    return {
        "node_count": len(offsets),
        "taylor_degree": TAYLOR_DEGREE,
        "derivative_orders": [BASE_H_ORDER + order for order in range(1, 9)],
        "all_exact_moments_passed": all(
            sum(
                weight * offset**power
                for weight, offset in zip(row, offsets)
            )
            == (
                Fraction(math.factorial(derivative))
                if power == derivative
                else Fraction(0)
            )
            for derivative, row in enumerate(weights, start=1)
            for power in range(NODE_COUNT)
        ),
    }


def higher_derivative_stencil_jet(
    target: Fraction,
    stencil_nodes: list[Fraction],
    h_rows: list[dict],
    *,
    point_h_source: dict[Fraction, tuple[list[flint.arb], dict]],
) -> tuple[list[flint.arb], dict]:
    """Append rigorous H9-H16 coefficients to one exact H0-H8 point jet."""
    if len(stencil_nodes) != NODE_COUNT or target not in stencil_nodes:
        raise ValueError("invalid order-eleven high-derivative stencil")
    stencil_nodes = sorted(stencil_nodes)
    offsets = tuple(node - target for node in stencil_nodes)
    weights = stencil_weights(offsets)
    try:
        target_series, target_diagnostics = point_h_source[target]
    except KeyError as exc:
        raise RuntimeError(f"point H source misses stencil target {target}") from exc
    if len(target_series) != BASE_H_ORDER + 1:
        raise RuntimeError("stencil target is not an exact H0-H8 jet")

    base_values = []
    for node in stencil_nodes:
        try:
            series, _ = point_h_source[node]
        except KeyError as exc:
            raise RuntimeError(f"point H source misses stencil node {node}") from exc
        if len(series) != BASE_H_ORDER + 1:
            raise RuntimeError(f"stencil node {node} is not an exact H0-H8 jet")
        base_values.append(series[BASE_H_ORDER] * math.factorial(BASE_H_ORDER))

    support_left = min(target, stencil_nodes[0])
    support_right = max(target, stencil_nodes[-1])
    support_rows = [
        row
        for row in h_rows
        if row["target_t_right"] > support_left
        and row["target_t_left"] < support_right
    ]
    if not support_rows:
        raise RuntimeError("H24 stencil remainder has no support rows")
    support_rows.sort(key=lambda row: row["target_t_left"])
    if (
        support_rows[0]["target_t_left"] > support_left
        or support_rows[-1]["target_t_right"] < support_right
    ):
        raise RuntimeError("H24 stencil remainder has incomplete support")
    if any(REMAINDER_H_ORDER not in row["H"] for row in support_rows):
        raise RuntimeError("stencil source misses H24")
    h24_upper = max(
        (upper_absolute(row["H"][REMAINDER_H_ORDER]) for row in support_rows),
        key=float,
    )

    augmented = list(target_series)
    error_uppers = {}
    for derivative, row_weights in enumerate(weights, start=1):
        estimate = flint.arb(0)
        for weight, value in zip(row_weights, base_values):
            estimate += arb_rational(weight) * value
        remainder_factor = sum(
            abs(weight) * abs(offset) ** (TAYLOR_DEGREE + 1)
            for weight, offset in zip(row_weights, offsets)
        ) / math.factorial(TAYLOR_DEGREE + 1)
        error = arb_rational(remainder_factor) * h24_upper
        h_order = BASE_H_ORDER + derivative
        augmented.append(
            (estimate + order10._symmetric(error)) / math.factorial(h_order)
        )
        error_uppers[str(h_order)] = arb_upper_text(error)

    diagnostics = {
        **target_diagnostics,
        "higher_derivative_stencil": {
            "base_h_order": BASE_H_ORDER,
            "derived_h_orders": [BASE_H_ORDER + 1, MAXIMUM_H_ORDER],
            "node_count": NODE_COUNT,
            "node_left": str(stencil_nodes[0]),
            "node_right": str(stencil_nodes[-1]),
            "offsets": [str(offset) for offset in offsets],
            "remainder_h_order": REMAINDER_H_ORDER,
            "h24_upper": arb_upper_text(h24_upper),
            "derivative_error_upper": error_uppers,
        },
    }
    return augmented, diagnostics
