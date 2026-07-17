#!/usr/bin/env python3
"""Propagate sparse exact H0-H23 anchors to rigorous half-grid H0-H16 jets."""

from __future__ import annotations

from fractions import Fraction
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


ANCHOR_MAXIMUM_H_ORDER = 23
REMAINDER_H_ORDER = 24
OUTPUT_MAXIMUM_H_ORDER = 16
MAXIMUM_PROPAGATION_DISTANCE = Fraction(1)


def nearest_anchor(
    target: Fraction,
    anchors: list[Fraction],
) -> Fraction:
    if not anchors:
        raise ValueError("empty sparse H anchor source")
    anchor = min(anchors, key=lambda value: (abs(value - target), value))
    if abs(anchor - target) > MAXIMUM_PROPAGATION_DISTANCE:
        raise RuntimeError(
            f"sparse H anchor is too far from {target}: {anchor}"
        )
    return anchor


def propagated_h0_h16_jet(
    target: Fraction,
    anchor: Fraction,
    anchor_series: list[flint.arb],
    h_rows: list[dict],
    *,
    point_diagnostics: dict,
) -> tuple[list[flint.arb], dict]:
    """Taylor-propagate one exact H0-H23 anchor using an H24 wall."""
    if len(anchor_series) != ANCHOR_MAXIMUM_H_ORDER + 1:
        raise RuntimeError("sparse anchor is not an exact H0-H23 jet")
    displacement = target - anchor
    if abs(displacement) > MAXIMUM_PROPAGATION_DISTANCE:
        raise ValueError("H23 propagation distance exceeds one")

    if displacement:
        support_left = min(anchor, target)
        support_right = max(anchor, target)
        support_rows = [
            row
            for row in h_rows
            if row["target_t_right"] > support_left
            and row["target_t_left"] < support_right
        ]
        if not support_rows:
            raise RuntimeError("H24 propagation has no support rows")
        support_rows.sort(key=lambda row: row["target_t_left"])
        if (
            support_rows[0]["target_t_left"] > support_left
            or support_rows[-1]["target_t_right"] < support_right
        ):
            raise RuntimeError("H24 propagation support is incomplete")
        if any(REMAINDER_H_ORDER not in row["H"] for row in support_rows):
            raise RuntimeError("H24 propagation source misses H24")
        h24_upper = max(
            (upper_absolute(row["H"][REMAINDER_H_ORDER]) for row in support_rows),
            key=float,
        )
    else:
        support_rows = []
        h24_upper = flint.arb(0)

    displacement_arb = arb_rational(displacement)
    radius_arb = arb_rational(abs(displacement))
    result = []
    error_uppers = {}
    for derivative in range(OUTPUT_MAXIMUM_H_ORDER + 1):
        coefficient = flint.arb(0)
        for anchor_order in range(derivative, ANCHOR_MAXIMUM_H_ORDER + 1):
            coefficient += (
                anchor_series[anchor_order]
                * math.comb(anchor_order, derivative)
                * displacement_arb ** (anchor_order - derivative)
            )
        if displacement:
            error = (
                h24_upper
                * radius_arb ** (REMAINDER_H_ORDER - derivative)
                / math.factorial(REMAINDER_H_ORDER - derivative)
                / math.factorial(derivative)
            )
            coefficient += order10._symmetric(error)
        else:
            error = flint.arb(0)
        result.append(coefficient)
        error_uppers[str(derivative)] = arb_upper_text(error)

    diagnostics = {
        **point_diagnostics,
        "sparse_h23_propagation": {
            "anchor_t": str(anchor),
            "target_t": str(target),
            "displacement": str(displacement),
            "anchor_h_orders": [0, ANCHOR_MAXIMUM_H_ORDER],
            "output_h_orders": [0, OUTPUT_MAXIMUM_H_ORDER],
            "remainder_h_order": REMAINDER_H_ORDER,
            "support_rows": len(support_rows),
            "h24_upper": arb_upper_text(h24_upper),
            "coefficient_error_upper": error_uppers,
        },
    }
    return result, diagnostics


def propagated_point_source(
    targets: list[Fraction],
    anchor_source: dict[Fraction, list[flint.arb]],
    h_rows: list[dict],
    *,
    anchor_diagnostic_source: dict[Fraction, dict],
) -> dict[Fraction, tuple[list[flint.arb], dict]]:
    anchors = sorted(anchor_source)
    result = {}
    for target in targets:
        anchor = nearest_anchor(target, anchors)
        try:
            diagnostics = anchor_diagnostic_source[anchor]
        except KeyError as exc:
            raise RuntimeError(
                f"sparse H diagnostic source misses anchor {anchor}"
            ) from exc
        result[target] = propagated_h0_h16_jet(
            target,
            anchor,
            anchor_source[anchor],
            h_rows,
            point_diagnostics=diagnostics,
        )
    return result
