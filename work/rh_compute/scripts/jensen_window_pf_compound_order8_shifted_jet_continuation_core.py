#!/usr/bin/env python3
"""Shifted-jet Taylor-continuation core for order-eight first-summand curvature."""

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

from jensen_window_pf_compound_order4_localized_curvature_interval_core import (  # noqa: E402
    hull,
)
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
    stable_log_series,
)
from jensen_window_pf_compound_order7_shifted_jet_continuation_core import (  # noqa: E402
    ExactTileHullConvolver,
    centered_second_difference,
    point_h_series,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    upper_absolute,
)


PRECISION_BITS = 384
POINT_ORDER = 9
REMAINDER_ORDER = 10
CURVATURE_CONSTANT = 2000


def point_hierarchy(anchor: Fraction, order: int = POINT_ORDER) -> dict:
    """Build the five-layer stable jet at one exact t anchor."""
    h_source: dict[int, list] = {}
    quadrature = []
    for shift in range(-6, 7):
        h_source[shift], diagnostics = point_h_series(anchor + shift, order)
        quadrature.append(diagnostics)

    b = {
        shift: centered_second_difference(h_source, shift)
        for shift in range(-5, 6)
    }
    ell = {
        shift: stable_log_series(b[shift], order)
        for shift in range(-5, 6)
    }
    j_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 2),
            centered_second_difference(ell, shift),
        )
        for shift in range(-4, 5)
    }
    h_stable = {
        shift: series_add(
            series_scale(ell[shift], 2),
            stable_log_series(j_coordinate[shift], order),
        )
        for shift in range(-4, 5)
    }
    r_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 3),
            centered_second_difference(h_stable, shift),
        )
        for shift in range(-3, 4)
    }
    q = {
        shift: series_add(
            series_sub(series_scale(h_stable[shift], 2), ell[shift]),
            stable_log_series(r_coordinate[shift], order),
        )
        for shift in range(-3, 4)
    }
    s_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 4),
            centered_second_difference(q, shift),
        )
        for shift in range(-2, 3)
    }
    p = {
        shift: series_add(
            series_sub(series_scale(q[shift], 2), h_stable[shift]),
            stable_log_series(s_coordinate[shift], order),
        )
        for shift in range(-2, 3)
    }
    t_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 5),
            centered_second_difference(p, shift),
        )
        for shift in range(-1, 2)
    }
    r_final = {
        shift: series_add(
            series_sub(series_scale(p[shift], 2), q[shift]),
            stable_log_series(t_coordinate[shift], order),
        )
        for shift in range(-1, 2)
    }
    u_coordinate = series_sub(
        series_scale(b[0], 6),
        centered_second_difference(r_final, 0),
    )
    s_final = series_add(
        series_sub(series_scale(r_final[0], 2), p[0]),
        stable_log_series(u_coordinate, order),
    )
    return {
        "B": b[0],
        "J": j_coordinate[0],
        "R": r_coordinate[0],
        "S": s_coordinate[0],
        "T": t_coordinate[0],
        "U": u_coordinate,
        "s": s_final,
        "quadrature": quadrature,
    }


def dimensionless_common_hierarchy(
    scale_t: Fraction,
    h_bounds: dict[int, flint.arb],
    order: int = REMAINDER_ORDER,
) -> dict:
    """Enclose normalized five-layer stable jets on one six-shift H collar."""
    if scale_t <= 0:
        raise ValueError("dimensionless scale must be positive")
    highest_h_order = order + 12
    missing = [
        degree
        for degree in range(2, highest_h_order + 1)
        if degree not in h_bounds
    ]
    if missing:
        raise ValueError(f"missing common H derivative bounds: {missing}")

    scale = arb_rational(scale_t)
    inverse_scale = 1 / scale

    def stable_gap_log(gap: list, output_order: int) -> list:
        physical_gap = [
            coefficient * inverse_scale
            for coefficient in gap[: output_order + 1]
        ]
        return stable_log_series(physical_gap, output_order)

    def centered_log_difference(log_jet: list, degree: int):
        return (
            (degree + 1)
            * (degree + 2)
            * log_jet[degree + 2]
            * inverse_scale
        )

    b = [
        scale ** (degree + 1)
        * h_bounds[degree + 2]
        / math.factorial(degree)
        for degree in range(order + 11)
    ]
    ell = stable_gap_log(b, order + 10)

    j_coordinate = [
        2 * b[degree] - centered_log_difference(ell, degree)
        for degree in range(order + 9)
    ]
    log_j = stable_gap_log(j_coordinate, order + 8)
    h_stable = [
        2 * ell[degree] + log_j[degree]
        for degree in range(order + 9)
    ]

    r_coordinate = [
        3 * b[degree] - centered_log_difference(h_stable, degree)
        for degree in range(order + 7)
    ]
    log_r = stable_gap_log(r_coordinate, order + 6)
    q = [
        2 * h_stable[degree] - ell[degree] + log_r[degree]
        for degree in range(order + 7)
    ]

    s_coordinate = [
        4 * b[degree] - centered_log_difference(q, degree)
        for degree in range(order + 5)
    ]
    log_s = stable_gap_log(s_coordinate, order + 4)
    p = [
        2 * q[degree] - h_stable[degree] + log_s[degree]
        for degree in range(order + 5)
    ]

    t_coordinate = [
        5 * b[degree] - centered_log_difference(p, degree)
        for degree in range(order + 3)
    ]
    log_t = stable_gap_log(t_coordinate, order + 2)
    r_final = [
        2 * p[degree] - q[degree] + log_t[degree]
        for degree in range(order + 3)
    ]

    u_coordinate = [
        6 * b[degree] - centered_log_difference(r_final, degree)
        for degree in range(order + 1)
    ]
    log_u = stable_gap_log(u_coordinate, order)
    s_final = [
        2 * r_final[degree] - p[degree] + log_u[degree]
        for degree in range(order + 1)
    ]
    return {
        "B": b,
        "ell": ell,
        "J": j_coordinate,
        "h": h_stable,
        "R": r_coordinate,
        "q": q,
        "S": s_coordinate,
        "p": p,
        "T": t_coordinate,
        "r": r_final,
        "U": u_coordinate,
        "s": s_final,
    }


def dimensionless_continuation_row(
    anchor: Fraction,
    right: Fraction,
    h_bounds: dict[int, flint.arb],
    *,
    remainder_order: int = REMAINDER_ORDER,
) -> dict:
    """Certify one rightward block with a normalized common-collar remainder."""
    flint.ctx.prec = PRECISION_BITS
    if not anchor < right:
        raise ValueError("continuation block must have positive width")
    if remainder_order != POINT_ORDER + 1:
        raise ValueError("remainder order must be one above the point-jet order")

    point = point_hierarchy(anchor, POINT_ORDER)
    common = dimensionless_common_hierarchy(anchor, h_bounds, remainder_order)
    coordinate_rows = {}
    for name in ("B", "J", "R", "S", "T", "U"):
        if not bool(point[name][0] > 0):
            raise RuntimeError(f"nonpositive point {name} at t={anchor}")
        if not bool(common[name][0] > 0):
            raise RuntimeError(
                f"nonpositive common-collar {name} on {anchor}..{right}"
            )
        coordinate_rows[name] = {
            "point_lower": arb_lower_text(point[name][0]),
            "dimensionless_common_lower": arb_lower_text(common[name][0]),
            "dimensionless_common_upper": arb_upper_text(common[name][0]),
        }

    width = arb_rational(right - anchor)
    ratio = width / arb_rational(anchor)
    curvature = 2 * point["s"][2]
    for degree in range(3, POINT_ORDER + 1):
        curvature += (
            upper_absolute(point["s"][degree])
            * degree
            * (degree - 1)
            * width ** (degree - 2)
        )

    normalized_remainder_coefficient = upper_absolute(
        common["s"][remainder_order]
    )
    curvature_remainder = (
        remainder_order
        * (remainder_order - 1)
        * normalized_remainder_coefficient
        / arb_rational(anchor) ** 2
        * ratio ** (remainder_order - 2)
    )
    curvature += curvature_remainder
    scaled = arb_rational(right) ** 2 * curvature
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    if not bool(margin > 0):
        raise RuntimeError(
            f"dimensionless order-eight curvature continuation failed on "
            f"{anchor}..{right}: {scaled}"
        )

    quadrature = point["quadrature"]
    return {
        "anchor": str(anchor),
        "right": str(right),
        "width": str(right - anchor),
        "relative_width": str((right - anchor) / anchor),
        "point_coordinates": {
            name: arb_lower_text(point[name][0])
            for name in ("B", "J", "R", "S", "T", "U")
        },
        "point_scaled_curvature": arb_upper_text(
            arb_rational(anchor) ** 2 * 2 * point["s"][2]
        ),
        "common_coordinates": coordinate_rows,
        "point_polynomial_curvature_upper": arb_upper_text(
            curvature - curvature_remainder
        ),
        "normalized_s_remainder_coefficient_upper": arb_upper_text(
            normalized_remainder_coefficient
        ),
        "curvature_remainder_upper": arb_upper_text(curvature_remainder),
        "scaled_curvature_remainder_upper": arb_upper_text(
            arb_rational(right) ** 2 * curvature_remainder
        ),
        "curvature_upper": arb_upper_text(curvature),
        "scaled_curvature_upper": arb_upper_text(scaled),
        "curvature_margin_lower": arb_lower_text(margin),
        "quadrature": {
            "shift_count": len(quadrature),
            "maximum_panel_error_upper": max(
                (row["maximum_panel_error_upper"] for row in quadrature),
                key=float,
            ),
            "maximum_tail_moment_upper": max(
                (row["maximum_tail_moment_upper"] for row in quadrature),
                key=float,
            ),
            "mode_brackets": [row["mode_bracket"] for row in quadrature],
        },
        "passed": True,
    }


def dimensionless_localized_hierarchy(
    scale_t: Fraction,
    h_rows: list[dict],
    order: int = REMAINDER_ORDER,
) -> dict[str, list[dict]]:
    """Propagate normalized jets through five localized tent convolutions."""
    if scale_t <= 0:
        raise ValueError("dimensionless scale must be positive")
    if not h_rows:
        raise ValueError("empty localized H collar")
    highest_h_order = order + 12
    for row in h_rows:
        missing = [
            degree
            for degree in range(2, highest_h_order + 1)
            if degree not in row["H"]
        ]
        if missing:
            raise ValueError(f"missing localized H derivatives: {missing}")

    scale = arb_rational(scale_t)
    inverse_scale = 1 / scale
    h_convolvers = {
        degree: ExactTileHullConvolver(
            h_rows,
            [row["H"][degree] for row in h_rows],
        )
        for degree in range(2, highest_h_order + 1)
    }
    b_rows = []
    for row in h_rows:
        averages = [
            h_convolvers[degree + 2].average(
                row["target_t_left"],
                row["target_t_right"],
            )
            for degree in range(order + 11)
        ]
        if any(value is None for value in averages):
            continue
        b = [
            scale ** (degree + 1) * value / math.factorial(degree)
            for degree, value in enumerate(averages)
            if value is not None
        ]
        if not bool(b[0] > 0):
            raise RuntimeError(
                f"nonpositive localized dimensionless B on "
                f"{row.get('target_t_left')}"
            )
        ell = stable_log_series(
            [coefficient * inverse_scale for coefficient in b],
            order + 10,
        )
        b_rows.append({**row, "B": b, "ell": ell})

    def stable_stage(
        rows: list[dict],
        *,
        coordinate_name: str,
        output_name: str,
        current_name: str,
        previous_name: str | None,
        leading_factor: int,
        output_order: int,
    ) -> list[dict]:
        convolvers = {
            degree: ExactTileHullConvolver(
                rows,
                [row[current_name][degree] for row in rows],
            )
            for degree in range(2, output_order + 3)
        }
        output = []
        for row in rows:
            averages = [
                convolvers[degree + 2].average(
                    row["target_t_left"],
                    row["target_t_right"],
                )
                for degree in range(output_order + 1)
            ]
            if any(value is None for value in averages):
                continue
            coordinate = [
                leading_factor * row["B"][degree]
                - (degree + 1)
                * (degree + 2)
                * value
                * inverse_scale
                for degree, value in enumerate(averages)
                if value is not None
            ]
            if not bool(coordinate[0] > 0):
                raise RuntimeError(
                    f"nonpositive localized dimensionless {coordinate_name} "
                    f"on {row.get('target_t_left')}: {coordinate[0]}"
                )
            stable = stable_log_series(
                [coefficient * inverse_scale for coefficient in coordinate],
                output_order,
            )
            base = series_scale(row[current_name][: output_order + 1], 2)
            if previous_name is not None:
                base = series_sub(
                    base,
                    row[previous_name][: output_order + 1],
                )
            output.append(
                {
                    **row,
                    coordinate_name: coordinate,
                    output_name: series_add(base, stable),
                }
            )
        return output

    j_rows = stable_stage(
        b_rows,
        coordinate_name="J",
        output_name="h",
        current_name="ell",
        previous_name=None,
        leading_factor=2,
        output_order=order + 8,
    )
    r_rows = stable_stage(
        j_rows,
        coordinate_name="R",
        output_name="q",
        current_name="h",
        previous_name="ell",
        leading_factor=3,
        output_order=order + 6,
    )
    s_rows = stable_stage(
        r_rows,
        coordinate_name="S",
        output_name="p",
        current_name="q",
        previous_name="h",
        leading_factor=4,
        output_order=order + 4,
    )
    t_rows = stable_stage(
        s_rows,
        coordinate_name="T",
        output_name="r",
        current_name="p",
        previous_name="q",
        leading_factor=5,
        output_order=order + 2,
    )
    u_rows = stable_stage(
        t_rows,
        coordinate_name="U",
        output_name="s",
        current_name="r",
        previous_name="p",
        leading_factor=6,
        output_order=order,
    )
    return {
        "B": b_rows,
        "J": j_rows,
        "R": r_rows,
        "S": s_rows,
        "T": t_rows,
        "U": u_rows,
    }


def localized_dimensionless_bounds(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    order: int = REMAINDER_ORDER,
) -> dict:
    """Hull final normalized coordinates and s remainder on a target block."""
    hierarchy = dimensionless_localized_hierarchy(anchor, h_rows, order)
    central = [
        row
        for row in hierarchy["U"]
        if row["target_t_right"] > anchor and row["target_t_left"] < right
    ]
    if not central:
        raise RuntimeError(f"localized hierarchy misses {anchor}..{right}")
    central.sort(key=lambda row: row["target_t_left"])
    if central[0]["target_t_left"] > anchor or central[-1]["target_t_right"] < right:
        raise RuntimeError(f"localized hierarchy does not cover {anchor}..{right}")
    for previous, current in zip(central, central[1:]):
        if previous["target_t_right"] < current["target_t_left"]:
            raise RuntimeError("localized hierarchy has a central coverage gap")

    coordinates = {
        name: hull([row[name][0] for row in central])
        for name in ("B", "J", "R", "S", "T", "U")
    }
    for name, value in coordinates.items():
        if not bool(value > 0):
            raise RuntimeError(
                f"nonpositive localized common {name} on {anchor}..{right}: "
                f"{value}"
            )
    return {
        "coordinates": coordinates,
        "s_remainder_coefficient": hull([row["s"][order] for row in central]),
        "stage_rows": {name: len(rows) for name, rows in hierarchy.items()},
        "central_rows": len(central),
        "central_t_left": central[0]["target_t_left"],
        "central_t_right": central[-1]["target_t_right"],
    }


def localized_dimensionless_continuation_row(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    *,
    remainder_order: int = REMAINDER_ORDER,
) -> dict:
    """Certify one block using staged localized dimensionless remainders."""
    flint.ctx.prec = PRECISION_BITS
    if not anchor < right:
        raise ValueError("continuation block must have positive width")
    if remainder_order != POINT_ORDER + 1:
        raise ValueError("remainder order must be one above the point-jet order")

    common = localized_dimensionless_bounds(anchor, right, h_rows, remainder_order)
    point = point_hierarchy(anchor, POINT_ORDER)
    coordinate_rows = {}
    for name in ("B", "J", "R", "S", "T", "U"):
        if not bool(point[name][0] > 0):
            raise RuntimeError(f"nonpositive point {name} at t={anchor}")
        coordinate_rows[name] = {
            "point_lower": arb_lower_text(point[name][0]),
            "dimensionless_common_lower": arb_lower_text(common["coordinates"][name]),
            "dimensionless_common_upper": arb_upper_text(common["coordinates"][name]),
        }

    width = arb_rational(right - anchor)
    ratio = width / arb_rational(anchor)
    point_polynomial = 2 * point["s"][2]
    for degree in range(3, POINT_ORDER + 1):
        point_polynomial += (
            upper_absolute(point["s"][degree])
            * degree
            * (degree - 1)
            * width ** (degree - 2)
        )

    normalized_remainder_coefficient = upper_absolute(
        common["s_remainder_coefficient"]
    )
    curvature_remainder = (
        remainder_order
        * (remainder_order - 1)
        * normalized_remainder_coefficient
        / arb_rational(anchor) ** 2
        * ratio ** (remainder_order - 2)
    )
    curvature = point_polynomial + curvature_remainder
    scaled = arb_rational(right) ** 2 * curvature
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    if not bool(margin > 0):
        raise RuntimeError(
            f"localized dimensionless order-eight curvature failed on "
            f"{anchor}..{right}: {scaled}"
        )

    quadrature = point["quadrature"]
    return {
        "anchor": str(anchor),
        "right": str(right),
        "width": str(right - anchor),
        "relative_width": str((right - anchor) / anchor),
        "point_coordinates": {
            name: arb_lower_text(point[name][0])
            for name in ("B", "J", "R", "S", "T", "U")
        },
        "common_coordinates": coordinate_rows,
        "point_scaled_curvature": arb_upper_text(
            arb_rational(anchor) ** 2 * 2 * point["s"][2]
        ),
        "point_polynomial_curvature_upper": arb_upper_text(point_polynomial),
        "normalized_s_remainder_coefficient_upper": arb_upper_text(
            normalized_remainder_coefficient
        ),
        "curvature_remainder_upper": arb_upper_text(curvature_remainder),
        "scaled_curvature_remainder_upper": arb_upper_text(
            arb_rational(right) ** 2 * curvature_remainder
        ),
        "curvature_upper": arb_upper_text(curvature),
        "scaled_curvature_upper": arb_upper_text(scaled),
        "curvature_margin_lower": arb_lower_text(margin),
        "localized_hierarchy": {
            "stage_rows": common["stage_rows"],
            "central_rows": common["central_rows"],
            "central_t_left": str(common["central_t_left"]),
            "central_t_right": str(common["central_t_right"]),
        },
        "quadrature": {
            "shift_count": len(quadrature),
            "maximum_panel_error_upper": max(
                (row["maximum_panel_error_upper"] for row in quadrature),
                key=float,
            ),
            "maximum_tail_moment_upper": max(
                (row["maximum_tail_moment_upper"] for row in quadrature),
                key=float,
            ),
            "mode_brackets": [row["mode_bracket"] for row in quadrature],
        },
        "passed": True,
    }
