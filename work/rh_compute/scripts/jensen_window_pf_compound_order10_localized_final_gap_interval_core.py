#!/usr/bin/env python3
"""Localized interval core for the order-ten seventh stable layer.

The completed order-nine hierarchy supplies ``w`` and the global theorem

    w_1''(t) <= 4200/t^2,  t >= 1250.

One exact-tile tent convolution appends

    W(t) = 8 B(t) - Delta^2 w(t),
    z(t) = 2 w(t) - s(t) + log(1-exp(-W(t))).

The routines below combine exact-point Taylor jets with a localized common
remainder.  Their intended output is a rigorous upper bound for ``z_1''``.
"""

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

from jensen_window_pf_compound_order8_shifted_jet_continuation_core import (  # noqa: E402
    ExactTileHullConvolver,
)
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
    centered_second_difference,
)
from jensen_window_pf_compound_order9_localized_final_gap_interval_core import (  # noqa: E402
    PRECISION_BITS,
    _apply_coordinate_floor,
    _interval_from_endpoints,
    _lower_point,
    _symmetric,
    _upper_point,
    analytic_v_floor,
    dimensionless_localized_sixth_hierarchy,
    inherited_coordinate_floor,
    point_sixth_hierarchy,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    upper_absolute,
)


MINIMUM_T = Fraction(1251)
ORDER9_CURVATURE_CONSTANT = 4200
CURVATURE_CONSTANT = 5500
POINT_ORDER = 5
REMAINDER_ORDER = 6


def analytic_w_floor_point(t: Fraction) -> flint.arb:
    """Pointwise floor inherited from B and the order-nine curvature theorem."""
    if t < MINIMUM_T:
        raise ValueError("the order-ten point W floor starts at t=1251")
    return inherited_coordinate_floor(
        t,
        8,
        Fraction(ORDER9_CURVATURE_CONSTANT),
    )


def analytic_w_floor_interval(left: Fraction, right: Fraction) -> flint.arb:
    """Uniform floor for W on an exact tile ``left <= t <= right``.

    The inputs ``B(t)>=1/(2t+3)`` and ``w''(t)<=4200/t^2`` imply

        W(t) >= 8/(2*right+3) - 4200/(left^2-1).

    The second term follows from the unit-tent identity and
    ``-log(1-left^-2) < 1/(left^2-1)``.
    """
    if not MINIMUM_T <= left < right:
        raise ValueError("invalid order-ten W-floor tile")
    left_arb = arb_rational(left)
    right_arb = arb_rational(right)
    floor = (
        flint.arb(8) / (2 * right_arb + 3)
        - flint.arb(ORDER9_CURVATURE_CONSTANT) / (left_arb**2 - 1)
    )
    if not bool(floor > 0):
        raise ValueError(f"nonpositive uniform W floor on {left}..{right}")
    return floor


def point_seventh_hierarchy(
    anchor: Fraction,
    order: int = POINT_ORDER,
    *,
    support_h_rows: list[dict],
    point_h_source: dict[Fraction, tuple[list, dict]],
) -> dict:
    """Build the rigorous seventh-layer Taylor jet at one exact anchor."""
    if anchor < MINIMUM_T:
        raise ValueError("the order-ten point hierarchy starts at t=1251")
    sixth = {
        shift: point_sixth_hierarchy(
            anchor + shift,
            order,
            point_h_source=point_h_source,
        )
        for shift in range(-1, 2)
    }
    center = sixth[0]
    w_source = {shift: value["w"] for shift, value in sixth.items()}
    w_coordinate = series_sub(
        series_scale(center["B"][: order + 1], 8),
        centered_second_difference(w_source, 0),
    )
    _apply_coordinate_floor(
        w_coordinate,
        analytic_w_floor_point(anchor),
        name="W",
        target=anchor,
    )
    z = series_add(
        series_sub(
            series_scale(center["w"][: order + 1], 2),
            center["s"][: order + 1],
        ),
        stable_log_series(w_coordinate, order),
    )
    coordinates = {
        name: center[name]
        for name in ("B", "J", "R", "S", "T", "U", "V")
    }
    coordinates["W"] = w_coordinate
    for name, coordinate in coordinates.items():
        if not bool(coordinate[0] > 0):
            raise RuntimeError(f"nonpositive point {name} at t={anchor}")
    quadrature = [
        row
        for shift in range(-1, 2)
        for row in sixth[shift]["quadrature"]
    ]
    return {
        **coordinates,
        "s": center["s"],
        "w": center["w"],
        "z": z,
        "quadrature": quadrature,
        "V_localizations": {
            str(shift): sixth[shift]["V_localization"]
            for shift in range(-1, 2)
        },
    }


def point_seventh_hierarchy_from_b_source(
    anchor: Fraction,
    order: int,
    *,
    point_b_source: dict[Fraction, tuple[list, dict]],
) -> dict:
    """Build the seventh point jet from direct tent-averaged B jets."""
    b = {}
    quadrature = []
    for shift in range(-7, 8):
        target = anchor + shift
        try:
            series, diagnostics = point_b_source[target]
        except KeyError as exc:
            raise RuntimeError(f"point B source misses target t={target}") from exc
        if len(series) < order + 1:
            raise RuntimeError(f"point B source at t={target} stops before {order}")
        b[shift] = series[: order + 1]
        quadrature.append(diagnostics)

    ell = {
        shift: stable_log_series(b[shift], order)
        for shift in range(-7, 8)
    }
    j_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 2),
            centered_second_difference(ell, shift),
        )
        for shift in range(-6, 7)
    }
    for shift, coordinate in j_coordinate.items():
        target = anchor + shift
        _apply_coordinate_floor(
            coordinate,
            flint.arb(1) / (7 * arb_rational(target)),
            name="J",
            target=target,
        )
    h_stable = {
        shift: series_add(
            series_scale(ell[shift], 2),
            stable_log_series(j_coordinate[shift], order),
        )
        for shift in range(-6, 7)
    }

    r_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 3),
            centered_second_difference(h_stable, shift),
        )
        for shift in range(-5, 6)
    }
    for shift, coordinate in r_coordinate.items():
        target = anchor + shift
        _apply_coordinate_floor(
            coordinate,
            inherited_coordinate_floor(target, 3, Fraction(7, 2)),
            name="R",
            target=target,
        )
    q = {
        shift: series_add(
            series_sub(series_scale(h_stable[shift], 2), ell[shift]),
            stable_log_series(r_coordinate[shift], order),
        )
        for shift in range(-5, 6)
    }

    s_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 4),
            centered_second_difference(q, shift),
        )
        for shift in range(-4, 5)
    }
    for shift, coordinate in s_coordinate.items():
        target = anchor + shift
        _apply_coordinate_floor(
            coordinate,
            inherited_coordinate_floor(target, 4, Fraction(60)),
            name="S",
            target=target,
        )
    p = {
        shift: series_add(
            series_sub(series_scale(q[shift], 2), h_stable[shift]),
            stable_log_series(s_coordinate[shift], order),
        )
        for shift in range(-4, 5)
    }

    t_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 5),
            centered_second_difference(p, shift),
        )
        for shift in range(-3, 4)
    }
    for shift, coordinate in t_coordinate.items():
        target = anchor + shift
        _apply_coordinate_floor(
            coordinate,
            inherited_coordinate_floor(target, 5, Fraction(200)),
            name="T",
            target=target,
        )
    r_final = {
        shift: series_add(
            series_sub(series_scale(p[shift], 2), q[shift]),
            stable_log_series(t_coordinate[shift], order),
        )
        for shift in range(-3, 4)
    }

    u_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 6),
            centered_second_difference(r_final, shift),
        )
        for shift in range(-2, 3)
    }
    for shift, coordinate in u_coordinate.items():
        target = anchor + shift
        _apply_coordinate_floor(
            coordinate,
            inherited_coordinate_floor(target, 6, Fraction(600)),
            name="U",
            target=target,
        )
    s_final = {
        shift: series_add(
            series_sub(series_scale(r_final[shift], 2), p[shift]),
            stable_log_series(u_coordinate[shift], order),
        )
        for shift in range(-2, 3)
    }

    v_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 7),
            centered_second_difference(s_final, shift),
        )
        for shift in range(-1, 2)
    }
    for shift, coordinate in v_coordinate.items():
        target = anchor + shift
        floor = inherited_coordinate_floor(target, 7, Fraction(2500))
        elementary = analytic_v_floor(target)
        if bool(elementary > floor):
            floor = elementary
        _apply_coordinate_floor(
            coordinate,
            floor,
            name="V",
            target=target,
        )
    w_final = {
        shift: series_add(
            series_sub(series_scale(s_final[shift], 2), r_final[shift]),
            stable_log_series(v_coordinate[shift], order),
        )
        for shift in range(-1, 2)
    }

    w_coordinate = series_sub(
        series_scale(b[0], 8),
        centered_second_difference(w_final, 0),
    )
    _apply_coordinate_floor(
        w_coordinate,
        analytic_w_floor_point(anchor),
        name="W",
        target=anchor,
    )
    z = series_add(
        series_sub(series_scale(w_final[0], 2), s_final[0]),
        stable_log_series(w_coordinate, order),
    )
    return {
        "B": b[0],
        "J": j_coordinate[0],
        "R": r_coordinate[0],
        "S": s_coordinate[0],
        "T": t_coordinate[0],
        "U": u_coordinate[0],
        "V": v_coordinate[0],
        "W": w_coordinate,
        "s": s_final[0],
        "w": w_final[0],
        "z": z,
        "quadrature": quadrature,
        "V_localizations": None,
    }


def dimensionless_localized_seventh_hierarchy(
    scale_t: Fraction,
    h_rows: list[dict],
    output_order: int = REMAINDER_ORDER,
) -> dict:
    """Append the localized seventh stable stage to the six-stage hierarchy."""
    if output_order < 2:
        raise ValueError("seventh-layer output order must be at least two")
    hierarchy = dimensionless_localized_sixth_hierarchy(
        scale_t,
        h_rows,
        output_order + 2,
    )
    scale = arb_rational(scale_t)
    v_rows = hierarchy["V"]
    convolvers = {
        degree: ExactTileHullConvolver(
            v_rows,
            [row["w"][degree] for row in v_rows],
        )
        for degree in range(2, output_order + 3)
    }
    w_rows = []
    for row in v_rows:
        left = row["target_t_left"]
        right = row["target_t_right"]
        averages = [
            convolvers[degree + 2].average(left, right)
            for degree in range(output_order + 1)
        ]
        if any(value is None for value in averages):
            continue
        coordinate = [
            8 * row["B"][degree]
            - (degree + 1) * (degree + 2) * value / scale
            for degree, value in enumerate(averages)
            if value is not None
        ]
        if left < MINIMUM_T:
            continue
        normalized_floor = scale * analytic_w_floor_interval(left, right)
        raw_upper = _upper_point(coordinate[0])
        raw_lower = _lower_point(coordinate[0])
        improved_lower = (
            normalized_floor
            if bool(normalized_floor > raw_lower)
            else raw_lower
        )
        if not bool(raw_upper >= improved_lower):
            raise RuntimeError(
                f"analytic W floor is disjoint from localized W on {left}..{right}"
            )
        coordinate[0] = _interval_from_endpoints(improved_lower, raw_upper)
        if not bool(coordinate[0] > 0):
            raise RuntimeError(f"nonpositive localized W on {left}..{right}")
        log_w = stable_log_series(
            [coefficient / scale for coefficient in coordinate],
            output_order,
        )
        z = series_add(
            series_sub(
                series_scale(row["w"][: output_order + 1], 2),
                row["s"][: output_order + 1],
            ),
            log_w,
        )
        w_rows.append({**row, "W": coordinate, "z": z})
    return {**hierarchy, "W": w_rows}


def localized_seventh_remainder_bounds(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    remainder_order: int = REMAINDER_ORDER,
) -> dict:
    """Hull the normalized seventh-layer remainder on a central block."""
    hierarchy = dimensionless_localized_seventh_hierarchy(
        anchor,
        h_rows,
        remainder_order,
    )
    central = [
        row
        for row in hierarchy["W"]
        if row["target_t_right"] > anchor
        and row["target_t_left"] < right
    ]
    if not central:
        raise RuntimeError(f"localized seventh hierarchy misses {anchor}..{right}")
    central.sort(key=lambda row: row["target_t_left"])
    if central[0]["target_t_left"] > anchor or central[-1]["target_t_right"] < right:
        raise RuntimeError("localized seventh hierarchy has incomplete central coverage")
    for previous, current in zip(central, central[1:]):
        if previous["target_t_right"] != current["target_t_left"]:
            raise RuntimeError("localized seventh hierarchy has a central gap")
    coordinates = {
        name: hull([row[name][0] for row in central])
        for name in ("B", "J", "R", "S", "T", "U", "V", "W")
    }
    for name, value in coordinates.items():
        if not bool(value > 0):
            raise RuntimeError(
                f"nonpositive localized common {name} on {anchor}..{right}"
            )
    return {
        "coordinates": coordinates,
        "z_remainder_coefficient": hull(
            [row["z"][remainder_order] for row in central]
        ),
        "stage_rows": {
            name: len(rows)
            for name, rows in hierarchy.items()
        },
        "central_rows": len(central),
        "central_t_left": central[0]["target_t_left"],
        "central_t_right": central[-1]["target_t_right"],
    }


def localized_seventh_continuation_row(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    *,
    point_order: int = POINT_ORDER,
    remainder_order: int = REMAINDER_ORDER,
    point_h_source: dict[Fraction, tuple[list, dict]],
) -> dict:
    """Certify one rightward block by a point jet and localized remainder."""
    flint.ctx.prec = PRECISION_BITS
    if not MINIMUM_T <= anchor < right:
        raise ValueError("invalid order-ten continuation block")
    if remainder_order != point_order + 1:
        raise ValueError("remainder order must be one above point order")
    common = localized_seventh_remainder_bounds(
        anchor,
        right,
        h_rows,
        remainder_order,
    )
    point = point_seventh_hierarchy(
        anchor,
        point_order,
        support_h_rows=h_rows,
        point_h_source=point_h_source,
    )
    coordinate_rows = {}
    for name in ("B", "J", "R", "S", "T", "U", "V", "W"):
        coordinate_rows[name] = {
            "point_lower": arb_lower_text(point[name][0]),
            "dimensionless_common_lower": arb_lower_text(
                common["coordinates"][name]
            ),
            "dimensionless_common_upper": arb_upper_text(
                common["coordinates"][name]
            ),
        }

    width = arb_rational(right - anchor)
    ratio = width / arb_rational(anchor)
    point_polynomial = 2 * point["z"][2]
    for degree in range(3, point_order + 1):
        point_polynomial += (
            upper_absolute(point["z"][degree])
            * degree
            * (degree - 1)
            * width ** (degree - 2)
        )
    normalized_remainder = upper_absolute(
        common["z_remainder_coefficient"]
    )
    curvature_remainder = (
        remainder_order
        * (remainder_order - 1)
        * normalized_remainder
        / arb_rational(anchor) ** 2
        * ratio ** (remainder_order - 2)
    )
    curvature = point_polynomial + curvature_remainder
    curvature_upper = _upper_point(curvature)
    if bool(curvature_upper > 0):
        scaled = _upper_point(arb_rational(right) ** 2 * curvature_upper)
    else:
        scaled = flint.arb(0)
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    if not bool(margin > 0):
        raise RuntimeError(
            f"order-ten localized continuation failed on {anchor}..{right}: "
            f"scaled={scaled}"
        )

    quadrature = point["quadrature"]
    return {
        "anchor": str(anchor),
        "right": str(right),
        "width": str(right - anchor),
        "relative_width": str((right - anchor) / anchor),
        "point_coordinates": {
            name: arb_lower_text(point[name][0])
            for name in ("B", "J", "R", "S", "T", "U", "V", "W")
        },
        "common_coordinates": coordinate_rows,
        "point_scaled_curvature": arb_upper_text(
            arb_rational(anchor) ** 2 * 2 * point["z"][2]
        ),
        "point_polynomial_curvature_upper": arb_upper_text(point_polynomial),
        "normalized_z_remainder_coefficient_upper": arb_upper_text(
            normalized_remainder
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


def localized_component_remainder_bounds(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    remainder_order: int = REMAINDER_ORDER,
    *,
    central_left: Fraction | None = None,
) -> dict:
    """Hull only the component remainders needed by the explicit z'' formula."""
    left = anchor if central_left is None else central_left
    if not left < right or not left <= anchor <= right:
        raise ValueError("invalid component remainder interval")
    hierarchy = dimensionless_localized_seventh_hierarchy(
        anchor,
        h_rows,
        remainder_order,
    )
    central = [
        row
        for row in hierarchy["W"]
        if row["target_t_right"] > left
        and row["target_t_left"] < right
    ]
    if not central:
        raise RuntimeError(f"localized component hierarchy misses {left}..{right}")
    central.sort(key=lambda row: row["target_t_left"])
    if central[0]["target_t_left"] > left or central[-1]["target_t_right"] < right:
        raise RuntimeError("localized component hierarchy has incomplete coverage")
    for previous, current in zip(central, central[1:]):
        if previous["target_t_right"] != current["target_t_left"]:
            raise RuntimeError("localized component hierarchy has a central gap")
    return {
        "s_remainder_coefficient": hull(
            [row["s"][remainder_order] for row in central]
        ),
        "w_remainder_coefficient": hull(
            [row["w"][remainder_order] for row in central]
        ),
        "W_remainder_coefficient": hull(
            [row["W"][remainder_order] for row in central]
        ),
        "stage_rows": {
            name: len(rows)
            for name, rows in hierarchy.items()
        },
        "central_rows": len(central),
        "central_t_left": central[0]["target_t_left"],
        "central_t_right": central[-1]["target_t_right"],
    }


def _taylor_derivative_enclosure(
    point_coefficients: list[flint.arb],
    common_remainder_coefficient: flint.arb,
    *,
    derivative: int,
    anchor: Fraction,
    left: Fraction,
    right: Fraction,
    point_order: int,
    remainder_order: int,
) -> tuple[flint.arb, flint.arb]:
    """Enclose one derivative on ``[anchor,right]`` by Taylor plus remainder."""
    if not 0 <= derivative <= point_order < remainder_order:
        raise ValueError("invalid Taylor derivative orders")
    local_left = arb_rational(left - anchor)
    local_right = arb_rational(right - anchor)
    local = (local_left + local_right) / 2 + flint.arb(
        0,
        _upper_point((local_right - local_left) / 2),
    )
    polynomial = flint.arb(0)
    for degree in range(point_order, derivative - 1, -1):
        falling = math.factorial(degree) // math.factorial(degree - derivative)
        polynomial = polynomial * local + falling * point_coefficients[degree]
    remainder_falling = (
        math.factorial(remainder_order)
        // math.factorial(remainder_order - derivative)
    )
    radius = max(abs(left - anchor), abs(right - anchor))
    ratio = arb_rational(radius) / arb_rational(anchor)
    remainder = (
        remainder_falling
        * upper_absolute(common_remainder_coefficient)
        / arb_rational(anchor) ** derivative
        * ratio ** (remainder_order - derivative)
    )
    return polynomial + _symmetric(remainder), remainder


def localized_seventh_formula_continuation_row(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    *,
    point_order: int = POINT_ORDER,
    remainder_order: int = REMAINDER_ORDER,
    point_h_source: dict[Fraction, tuple[list, dict]] | None = None,
    point_b_source: dict[Fraction, tuple[list, dict]] | None = None,
    require_pass: bool = True,
    block_left: Fraction | None = None,
) -> dict:
    """Certify a block through the explicit stable-log second derivative.

    This keeps the cancellation-preserving point jets for ``s``, ``w``, and
    ``W`` separate.  It avoids the very pessimistic interval for the sixth
    derivative of the composed logarithm.
    """
    flint.ctx.prec = PRECISION_BITS
    left = anchor if block_left is None else block_left
    if not MINIMUM_T <= left < right or not left <= anchor <= right:
        raise ValueError("invalid order-ten formula continuation block")
    if remainder_order != point_order + 1:
        raise ValueError("remainder order must be one above point order")
    common = localized_component_remainder_bounds(
        anchor,
        right,
        h_rows,
        remainder_order,
        central_left=left,
    )
    if point_b_source is not None:
        point = point_seventh_hierarchy_from_b_source(
            anchor,
            point_order,
            point_b_source=point_b_source,
        )
    elif point_h_source is not None:
        point = point_seventh_hierarchy(
            anchor,
            point_order,
            support_h_rows=h_rows,
            point_h_source=point_h_source,
        )
    else:
        raise ValueError("either point_h_source or point_b_source is required")

    enclosures = {}
    remainders = {}
    for name, coefficients, common_coefficient, degrees in (
        (
            "s",
            point["s"],
            common["s_remainder_coefficient"],
            (2,),
        ),
        (
            "w",
            point["w"],
            common["w_remainder_coefficient"],
            (2,),
        ),
        (
            "W",
            point["W"],
            common["W_remainder_coefficient"],
            (0, 1, 2),
        ),
    ):
        for degree in degrees:
            enclosed, remainder = _taylor_derivative_enclosure(
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

    floor = analytic_w_floor_interval(left, right)
    raw_w_lower = _lower_point(enclosures["W0"])
    w_lower = floor if bool(floor > raw_w_lower) else raw_w_lower
    w_upper = _upper_point(enclosures["W0"])
    if not bool(w_upper >= w_lower):
        raise RuntimeError(
            f"analytic W floor is disjoint from formula localization on "
            f"{left}..{right}"
        )
    phi_upper = _upper_point(
        1 / (w_lower.exp() - 1)
    )
    base_upper = _upper_point(2 * enclosures["w2"] - enclosures["s2"])
    w_second_upper = _upper_point(enclosures["W2"])
    if not bool(w_second_upper > 0):
        w_second_upper = flint.arb(0)
    positive_term = phi_upper * w_second_upper
    curvature_upper = _upper_point(base_upper + positive_term)
    if bool(curvature_upper > 0):
        scaled = _upper_point(arb_rational(right) ** 2 * curvature_upper)
    else:
        scaled = flint.arb(0)
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    passed = bool(margin > 0)
    if require_pass and not passed:
        raise RuntimeError(
            f"order-ten formula continuation failed on {anchor}..{right}: "
            f"scaled={scaled}"
        )

    return {
        "anchor": str(left),
        "expansion_anchor": str(anchor),
        "right": str(right),
        "width": str(right - left),
        "point_scaled_curvature": arb_upper_text(
            arb_rational(anchor) ** 2 * 2 * point["z"][2]
        ),
        "point_coordinates": {
            name: arb_lower_text(point[name][0])
            for name in ("B", "J", "R", "S", "T", "U", "V", "W")
        },
        "W_lower": arb_lower_text(w_lower),
        "W_upper": arb_upper_text(w_upper),
        "base_second_upper": arb_upper_text(base_upper),
        "positive_phi_W_second_upper": arb_upper_text(positive_term),
        "curvature_upper": arb_upper_text(curvature_upper),
        "scaled_curvature_upper": arb_upper_text(scaled),
        "curvature_margin_lower": arb_lower_text(margin),
        "component_enclosures": {
            name: value.str(50).replace("e", "E")
            for name, value in enclosures.items()
        },
        "component_remainder_upper": {
            name: arb_upper_text(value)
            for name, value in remainders.items()
        },
        "common_remainder_coefficient_upper": {
            "s": arb_upper_text(
                upper_absolute(common["s_remainder_coefficient"])
            ),
            "w": arb_upper_text(
                upper_absolute(common["w_remainder_coefficient"])
            ),
            "W": arb_upper_text(
                upper_absolute(common["W_remainder_coefficient"])
            ),
        },
        "localized_hierarchy": {
            "stage_rows": common["stage_rows"],
            "central_rows": common["central_rows"],
            "central_t_left": str(common["central_t_left"]),
            "central_t_right": str(common["central_t_right"]),
        },
        "proof_formula": (
            "z''=2*w''-s''+phi(W)*W''-chi(W)*(W')^2 "
            "<=2*w''-s''+phi(W)*max(W'',0)"
        ),
        "passed": passed,
    }
