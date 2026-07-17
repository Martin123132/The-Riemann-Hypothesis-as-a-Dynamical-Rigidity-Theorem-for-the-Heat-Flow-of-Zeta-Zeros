#!/usr/bin/env python3
"""Localized interval core for the order-eleven eighth stable layer."""

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order10_localized_final_gap_interval_core as order10  # noqa: E402
from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_scale,
    series_sub,
    stable_log_series,
)
from jensen_window_pf_compound_order7_shifted_jet_continuation_core import (  # noqa: E402
    centered_second_difference,
)
from jensen_window_pf_compound_order8_shifted_jet_continuation_core import (  # noqa: E402
    ExactTileHullConvolver,
)
from jensen_window_pf_compound_order9_localized_final_gap_interval_core import (  # noqa: E402
    inherited_coordinate_floor,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    upper_absolute,
)


PRECISION_BITS = order10.PRECISION_BITS
MINIMUM_T = Fraction(1252)
ORDER10_CURVATURE_CONSTANT = 4200
CURVATURE_CONSTANT = 6000
POINT_ORDER = 5
REMAINDER_ORDER = 6


def analytic_x_floor_point(t: Fraction) -> flint.arb:
    """Pointwise floor inherited from B and the order-ten z curvature target."""
    if t < MINIMUM_T:
        raise ValueError("the order-eleven point X floor starts at t=1252")
    return inherited_coordinate_floor(
        t,
        9,
        Fraction(ORDER10_CURVATURE_CONSTANT),
    )


def analytic_x_floor_interval(left: Fraction, right: Fraction) -> flint.arb:
    """Uniform X floor on an exact tile under z_1''(t)<=4200/t^2."""
    if not MINIMUM_T <= left < right:
        raise ValueError("invalid order-eleven X-floor tile")
    left_arb = arb_rational(left)
    right_arb = arb_rational(right)
    floor = (
        flint.arb(9) / (2 * right_arb + 3)
        - flint.arb(ORDER10_CURVATURE_CONSTANT) / (left_arb**2 - 1)
    )
    if not bool(floor > 0):
        raise ValueError(f"nonpositive uniform X floor on {left}..{right}")
    return floor


def point_eighth_hierarchy(
    anchor: Fraction,
    order: int = POINT_ORDER,
    *,
    support_h_rows: list[dict],
    point_h_source: dict[Fraction, tuple[list, dict]],
) -> dict:
    """Build the rigorous eighth-layer Taylor jet at one exact anchor."""
    if anchor < MINIMUM_T:
        raise ValueError("the order-eleven point hierarchy starts at t=1252")
    seventh = {
        shift: order10.point_seventh_hierarchy(
            anchor + shift,
            order,
            support_h_rows=support_h_rows,
            point_h_source=point_h_source,
        )
        for shift in range(-1, 2)
    }
    center = seventh[0]
    z_source = {shift: value["z"] for shift, value in seventh.items()}
    x_coordinate = series_sub(
        series_scale(center["B"][: order + 1], 9),
        centered_second_difference(z_source, 0),
    )
    order10._apply_coordinate_floor(
        x_coordinate,
        analytic_x_floor_point(anchor),
        name="X",
        target=anchor,
    )
    y = series_add(
        series_sub(
            series_scale(center["z"][: order + 1], 2),
            center["w"][: order + 1],
        ),
        stable_log_series(x_coordinate, order),
    )
    coordinates = {
        name: center[name]
        for name in ("B", "J", "R", "S", "T", "U", "V", "W")
    }
    coordinates["X"] = x_coordinate
    for name, coordinate in coordinates.items():
        if not bool(coordinate[0] > 0):
            raise RuntimeError(f"nonpositive point {name} at t={anchor}")
    quadrature = [
        row
        for shift in range(-1, 2)
        for row in seventh[shift]["quadrature"]
    ]
    return {
        **coordinates,
        "w": center["w"],
        "z": center["z"],
        "y": y,
        "quadrature": quadrature,
    }


def dimensionless_localized_eighth_hierarchy(
    scale_t: Fraction,
    h_rows: list[dict],
    output_order: int = REMAINDER_ORDER,
) -> dict:
    """Append an eighth localized stable stage to the order-ten hierarchy."""
    if output_order < 2:
        raise ValueError("eighth-layer output order must be at least two")
    hierarchy = order10.dimensionless_localized_seventh_hierarchy(
        scale_t,
        h_rows,
        output_order + 2,
    )
    scale = arb_rational(scale_t)
    w_rows = hierarchy["W"]
    convolvers = {
        degree: ExactTileHullConvolver(
            w_rows,
            [row["z"][degree] for row in w_rows],
        )
        for degree in range(2, output_order + 3)
    }
    x_rows = []
    for row in w_rows:
        left = row["target_t_left"]
        right = row["target_t_right"]
        averages = [
            convolvers[degree + 2].average(left, right)
            for degree in range(output_order + 1)
        ]
        if any(value is None for value in averages):
            continue
        coordinate = [
            9 * row["B"][degree]
            - (degree + 1) * (degree + 2) * value / scale
            for degree, value in enumerate(averages)
            if value is not None
        ]
        if left < MINIMUM_T:
            continue
        normalized_floor = scale * analytic_x_floor_interval(left, right)
        raw_upper = order10._upper_point(coordinate[0])
        raw_lower = order10._lower_point(coordinate[0])
        improved_lower = normalized_floor if bool(normalized_floor > raw_lower) else raw_lower
        if not bool(raw_upper >= improved_lower):
            raise RuntimeError(
                f"analytic X floor is disjoint from localized X on {left}..{right}"
            )
        coordinate[0] = order10._interval_from_endpoints(improved_lower, raw_upper)
        if not bool(coordinate[0] > 0):
            raise RuntimeError(f"nonpositive localized X on {left}..{right}")
        log_x = stable_log_series(
            [coefficient / scale for coefficient in coordinate],
            output_order,
        )
        y = series_add(
            series_sub(
                series_scale(row["z"][: output_order + 1], 2),
                row["w"][: output_order + 1],
            ),
            log_x,
        )
        x_rows.append({**row, "X": coordinate, "y": y})
    return {**hierarchy, "X": x_rows}


def localized_eighth_component_remainder_bounds(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    remainder_order: int = REMAINDER_ORDER,
    *,
    central_left: Fraction | None = None,
) -> dict:
    """Hull component remainders needed by the explicit y'' formula."""
    left = anchor if central_left is None else central_left
    if not left < right or not left <= anchor <= right:
        raise ValueError("invalid order-eleven component interval")
    hierarchy = dimensionless_localized_eighth_hierarchy(
        anchor,
        h_rows,
        remainder_order,
    )
    central = [
        row
        for row in hierarchy["X"]
        if row["target_t_right"] > left and row["target_t_left"] < right
    ]
    if not central:
        raise RuntimeError(f"localized eighth hierarchy misses {left}..{right}")
    central.sort(key=lambda row: row["target_t_left"])
    if central[0]["target_t_left"] > left or central[-1]["target_t_right"] < right:
        raise RuntimeError("localized eighth hierarchy has incomplete coverage")
    for previous, current in zip(central, central[1:]):
        if previous["target_t_right"] != current["target_t_left"]:
            raise RuntimeError("localized eighth hierarchy has a central gap")
    coordinates = {
        name: order10.hull([row[name][0] for row in central])
        for name in ("B", "J", "R", "S", "T", "U", "V", "W", "X")
    }
    for name, value in coordinates.items():
        if not bool(value > 0):
            raise RuntimeError(f"nonpositive localized common {name} on {left}..{right}")
    return {
        "coordinates": coordinates,
        "w_remainder_coefficient": order10.hull(
            [row["w"][remainder_order] for row in central]
        ),
        "z_remainder_coefficient": order10.hull(
            [row["z"][remainder_order] for row in central]
        ),
        "X_remainder_coefficient": order10.hull(
            [row["X"][remainder_order] for row in central]
        ),
        "y_remainder_coefficient": order10.hull(
            [row["y"][remainder_order] for row in central]
        ),
        "stage_rows": {name: len(rows) for name, rows in hierarchy.items()},
        "central_rows": len(central),
        "central_t_left": central[0]["target_t_left"],
        "central_t_right": central[-1]["target_t_right"],
    }


def localized_eighth_formula_continuation_row(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    *,
    point_order: int = POINT_ORDER,
    remainder_order: int = REMAINDER_ORDER,
    point_h_source: dict[Fraction, tuple[list, dict]],
    require_pass: bool = True,
    block_left: Fraction | None = None,
) -> dict:
    """Certify one block through the explicit stable-log y'' formula."""
    flint.ctx.prec = PRECISION_BITS
    left = anchor if block_left is None else block_left
    if not MINIMUM_T <= left < right or not left <= anchor <= right:
        raise ValueError("invalid order-eleven formula continuation block")
    if remainder_order != point_order + 1:
        raise ValueError("remainder order must be one above point order")
    common = localized_eighth_component_remainder_bounds(
        anchor,
        right,
        h_rows,
        remainder_order,
        central_left=left,
    )
    point = point_eighth_hierarchy(
        anchor,
        point_order,
        support_h_rows=h_rows,
        point_h_source=point_h_source,
    )
    enclosures = {}
    remainders = {}
    for name, coefficients, common_coefficient, degrees in (
        ("w", point["w"], common["w_remainder_coefficient"], (2,)),
        ("z", point["z"], common["z_remainder_coefficient"], (2,)),
        ("X", point["X"], common["X_remainder_coefficient"], (0, 1, 2)),
    ):
        for degree in degrees:
            enclosed, remainder = order10._taylor_derivative_enclosure(
                coefficients,
                common_coefficient,
                derivative=degree,
                anchor=anchor,
                left=left,
                right=right,
                point_order=point_order,
                remainder_order=remainder_order,
            )
            enclosures[f"{name}{degree}"] = enclosed
            remainders[f"{name}{degree}"] = remainder

    floor = analytic_x_floor_interval(left, right)
    raw_x_lower = order10._lower_point(enclosures["X0"])
    x_lower = floor if bool(floor > raw_x_lower) else raw_x_lower
    x_upper = order10._upper_point(enclosures["X0"])
    if not bool(x_upper >= x_lower):
        raise RuntimeError(
            f"analytic X floor is disjoint from formula localization on {left}..{right}"
        )
    phi_upper = order10._upper_point(1 / (x_lower.exp() - 1))
    base_upper = order10._upper_point(2 * enclosures["z2"] - enclosures["w2"])
    x_second_upper = order10._upper_point(enclosures["X2"])
    if not bool(x_second_upper > 0):
        x_second_upper = flint.arb(0)
    positive_term = phi_upper * x_second_upper
    curvature_upper = order10._upper_point(base_upper + positive_term)
    if bool(curvature_upper > 0):
        scaled = order10._upper_point(arb_rational(right) ** 2 * curvature_upper)
    else:
        scaled = flint.arb(0)
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    passed = bool(margin > 0)
    if require_pass and not passed:
        raise RuntimeError(
            f"order-eleven formula continuation failed on {anchor}..{right}: scaled={scaled}"
        )
    quadrature = point["quadrature"]
    return {
        "anchor": str(left),
        "expansion_anchor": str(anchor),
        "right": str(right),
        "width": str(right - left),
        "point_scaled_curvature": arb_upper_text(
            arb_rational(anchor) ** 2 * 2 * point["y"][2]
        ),
        "point_coordinates": {
            name: arb_lower_text(point[name][0])
            for name in ("B", "J", "R", "S", "T", "U", "V", "W", "X")
        },
        "X_lower": arb_lower_text(x_lower),
        "X_upper": arb_upper_text(x_upper),
        "base_second_upper": arb_upper_text(base_upper),
        "positive_phi_X_second_upper": arb_upper_text(positive_term),
        "curvature_upper": arb_upper_text(curvature_upper),
        "scaled_curvature_upper": arb_upper_text(scaled),
        "curvature_margin_lower": arb_lower_text(margin),
        "component_enclosures": {
            name: value.str(50).replace("e", "E") for name, value in enclosures.items()
        },
        "component_remainder_upper": {
            name: arb_upper_text(value) for name, value in remainders.items()
        },
        "common_remainder_coefficient_upper": {
            "w": arb_upper_text(upper_absolute(common["w_remainder_coefficient"])),
            "z": arb_upper_text(upper_absolute(common["z_remainder_coefficient"])),
            "X": arb_upper_text(upper_absolute(common["X_remainder_coefficient"])),
            "y": arb_upper_text(upper_absolute(common["y_remainder_coefficient"])),
        },
        "localized_hierarchy": {
            "stage_rows": common["stage_rows"],
            "central_rows": common["central_rows"],
            "central_t_left": str(common["central_t_left"]),
            "central_t_right": str(common["central_t_right"]),
        },
        "quadrature": {
            "shift_count": len(quadrature),
            "maximum_panel_error_upper": max(
                (row["maximum_panel_error_upper"] for row in quadrature), key=float
            ),
            "maximum_tail_moment_upper": max(
                (row["maximum_tail_moment_upper"] for row in quadrature), key=float
            ),
            "mode_brackets": [row["mode_bracket"] for row in quadrature],
        },
        "proof_formula": "y''=2*z''-w''+phi(X)*X''-chi(X)*(X')^2 <=2*z''-w''+phi(X)*max(X'',0)",
        "passed": passed,
    }
