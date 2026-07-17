#!/usr/bin/env python3
"""Localized interval core for the order-nine sixth stable layer.

The first five stable layers are propagated with exact-tile tent hulls.  The
last centered second difference is then localized by

    |Delta^2 f^(r)(t) - f^(r+2)(t)|
        <= sup_[t-1,t+1] |f^(r+4)| / 12,   r = 0, 1, 2.

This avoids evaluating the sixth logarithm as one highly dependent interval
jet while retaining a rigorous upper bound for its second derivative.
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
    dimensionless_localized_hierarchy,
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
    point_h_series,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    upper_absolute,
)


PRECISION_BITS = 512
LOCALIZED_ORDER = 6
CURVATURE_CONSTANT = 4200
MINIMUM_T = Fraction(1250)
POINT_ORDER = 7
REMAINDER_ORDER = 8


def _upper_point(value: flint.arb) -> flint.arb:
    return flint.arb(value.upper())


def _lower_point(value: flint.arb) -> flint.arb:
    return flint.arb(value.lower())


def _symmetric(radius: flint.arb) -> flint.arb:
    return flint.arb(0, _upper_point(radius))


def _interval_from_endpoints(
    lower: flint.arb,
    upper: flint.arb,
) -> flint.arb:
    lower_point = _lower_point(lower)
    upper_point = _upper_point(upper)
    if not bool(upper_point >= lower_point):
        raise ValueError("reversed interval endpoints")
    return (lower_point + upper_point) / 2 + flint.arb(
        0, _upper_point((upper_point - lower_point) / 2)
    )


def stable_phi(value: flint.arb) -> flint.arb:
    """Return phi(x)=1/(exp(x)-1) on a positive point enclosure."""
    if not bool(value > 0):
        raise ValueError("stable phi requires a positive argument")
    return 1 / (value.exp() - 1)


def stable_chi(value: flint.arb) -> flint.arb:
    """Return chi(x)=exp(x)/(exp(x)-1)^2 on a positive point enclosure."""
    if not bool(value > 0):
        raise ValueError("stable chi requires a positive argument")
    exponential = value.exp()
    return exponential / (exponential - 1) ** 2


def minimum_absolute(value: flint.arb) -> flint.arb:
    """A nonnegative point lower bound for |value|."""
    if bool(value > 0):
        return _lower_point(value)
    if bool(value < 0):
        return -_upper_point(value)
    return flint.arb(0)


def analytic_v_floor(t_right: Fraction) -> flint.arb:
    """Continuous floor V(t)>=4/(3t), minimized at a tile's right edge."""
    if t_right < MINIMUM_T:
        raise ValueError("the analytic sixth-gap floor starts at t=1250")
    return flint.arb(4) / (3 * arb_rational(t_right))


def analytic_v_floor_polynomial(t: Fraction) -> Fraction:
    """Numerator of 7/(2t+3)-2501/t^2-4/(3t)."""
    return 13 * t * t - 15018 * t - 22509


def inherited_coordinate_floor(
    t: Fraction,
    leading_factor: int,
    previous_curvature_constant: Fraction,
) -> flint.arb:
    """Floor m*B-Delta^2 f from B and a prior curvature theorem.

    The inherited inputs are B(t)>=1/(2t+3) and
    f''(x)<=C/x^2.  The unit-tent identity and
    -log(1-t^-2)<1/(t^2-1) then give the displayed floor.
    """
    t_arb = arb_rational(t)
    floor = (
        flint.arb(leading_factor) / (2 * t_arb + 3)
        - arb_rational(previous_curvature_constant) / (t_arb**2 - 1)
    )
    if not bool(floor > 0):
        raise ValueError(
            f"nonpositive inherited coordinate floor at t={t}"
        )
    return floor


def _apply_coordinate_floor(
    coordinate: list[flint.arb],
    floor: flint.arb,
    *,
    name: str,
    target: Fraction,
) -> None:
    raw_lower = _lower_point(coordinate[0])
    raw_upper = _upper_point(coordinate[0])
    improved_lower = floor if bool(floor > raw_lower) else raw_lower
    if not bool(raw_upper >= improved_lower):
        raise RuntimeError(
            f"analytic {name} floor is disjoint at t={target}"
        )
    coordinate[0] = _interval_from_endpoints(improved_lower, raw_upper)


def _physical_gap_derivative(
    coefficients: list[flint.arb],
    degree: int,
    scale: flint.arb,
) -> flint.arb:
    return (
        math.factorial(degree)
        * coefficients[degree]
        / scale ** (degree + 1)
    )


def _physical_log_derivative(
    coefficients: list[flint.arb],
    degree: int,
    scale: flint.arb,
) -> flint.arb:
    return math.factorial(degree) * coefficients[degree] / scale**degree


def localized_final_gap_rows(
    scale_t: Fraction,
    h_rows: list[dict],
    *,
    central_left: Fraction,
    central_right: Fraction,
    order: int = LOCALIZED_ORDER,
) -> dict:
    """Enclose w_1'' tilewise on one exact-t H2-H18 collar.

    Every input row must have exact rational ``target_t_left/right`` fields and
    interval ``H`` derivatives.  The returned rows cover the requested central
    interval whenever the supplied seven-unit collar is complete.
    """
    flint.ctx.prec = PRECISION_BITS
    if order != LOCALIZED_ORDER:
        raise ValueError("the final-gap localization requires order six")
    if not MINIMUM_T <= central_left < central_right:
        raise ValueError("invalid order-nine central interval")
    if analytic_v_floor_polynomial(central_left) <= 0:
        raise ValueError("continuous sixth-gap floor polynomial is not positive")

    hierarchy = dimensionless_localized_hierarchy(scale_t, h_rows, order)
    final_rows = hierarchy["U"]
    scale = arb_rational(scale_t)
    support_convolvers = {
        degree: ExactTileHullConvolver(
            final_rows,
            [row["s"][degree] for row in final_rows],
        )
        for degree in (4, 5, 6)
    }

    rows = []
    for row in final_rows:
        left = row["target_t_left"]
        right = row["target_t_right"]
        if right <= central_left or left >= central_right:
            continue
        if left < central_left or right > central_right:
            raise ValueError(
                "central interval must be aligned with the exact source tiles"
            )

        support = {
            degree: support_convolvers[degree].average(left, right)
            for degree in (4, 5, 6)
        }
        if any(value is None for value in support.values()):
            raise RuntimeError(f"missing final unit collar on {left}..{right}")

        b_derivatives = {
            degree: _physical_gap_derivative(row["B"], degree, scale)
            for degree in range(3)
        }
        r_second = _physical_log_derivative(row["r"], 2, scale)
        s_derivatives = {
            degree: _physical_log_derivative(row["s"], degree, scale)
            for degree in range(2, 5)
        }
        support_s = {
            degree: _physical_log_derivative(
                [flint.arb(0)] * degree + [support[degree]],
                degree,
                scale,
            )
            for degree in (4, 5, 6)
        }

        local_v = {
            degree: 7 * b_derivatives[degree] - s_derivatives[degree + 2]
            for degree in range(3)
        }
        remainders = {
            degree: upper_absolute(support_s[degree + 4]) / 12
            for degree in range(3)
        }
        v_enclosures = {
            degree: local_v[degree] + _symmetric(remainders[degree])
            for degree in range(3)
        }

        floor = analytic_v_floor(right)
        local_v_lower = _lower_point(v_enclosures[0])
        v_lower = floor if bool(floor > local_v_lower) else local_v_lower
        v_upper = _upper_point(v_enclosures[0])
        if not bool(v_upper >= v_lower):
            raise RuntimeError(
                f"analytic V floor is disjoint from localization on {left}..{right}"
            )

        phi_upper = _upper_point(stable_phi(v_lower))
        chi_lower = _lower_point(stable_chi(v_upper))
        v_second_upper = _upper_point(v_enclosures[2])
        if not bool(v_second_upper > 0):
            v_second_upper = flint.arb(0)
        v_first_minimum_absolute = minimum_absolute(v_enclosures[1])

        base = 2 * s_derivatives[2] - r_second
        base_upper = _upper_point(base)
        positive_term = phi_upper * v_second_upper
        negative_term = chi_lower * v_first_minimum_absolute**2
        curvature_upper = _upper_point(
            base_upper + positive_term - negative_term
        )
        if bool(curvature_upper > 0):
            scaled_upper = _upper_point(
                arb_rational(right) ** 2 * curvature_upper
            )
        else:
            scaled_upper = flint.arb(0)
        margin = flint.arb(CURVATURE_CONSTANT) - scaled_upper

        rows.append(
            {
                "target_t_left": str(left),
                "target_t_right": str(right),
                "analytic_V_floor": arb_lower_text(floor),
                "localized_V_lower": arb_lower_text(v_lower),
                "localized_V_upper": arb_upper_text(v_upper),
                "V_prime_ball": v_enclosures[1].str(40).replace("e", "E"),
                "V_second_upper": arb_upper_text(v_enclosures[2]),
                "V_prime_minimum_absolute": arb_lower_text(
                    v_first_minimum_absolute
                ),
                "remainder_E0_upper": arb_upper_text(remainders[0]),
                "remainder_E1_upper": arb_upper_text(remainders[1]),
                "remainder_E2_upper": arb_upper_text(remainders[2]),
                "base_second_upper": arb_upper_text(base),
                "positive_phi_V_second_upper": arb_upper_text(positive_term),
                "negative_chi_V_prime_square_lower": arb_lower_text(
                    negative_term
                ),
                "curvature_upper": arb_upper_text(curvature_upper),
                "scaled_curvature_upper": arb_upper_text(scaled_upper),
                "curvature_margin_lower": arb_lower_text(margin),
                "coordinates": {
                    name: arb_lower_text(row[name][0])
                    for name in ("B", "J", "R", "S", "T", "U")
                },
                "passed": bool(margin > 0),
            }
        )

    if not rows:
        raise RuntimeError("localized hierarchy produced no central rows")
    if Fraction(rows[0]["target_t_left"]) != central_left:
        raise RuntimeError("localized final-gap rows miss the central left edge")
    if Fraction(rows[-1]["target_t_right"]) != central_right:
        raise RuntimeError("localized final-gap rows miss the central right edge")
    for previous, current in zip(rows, rows[1:]):
        if previous["target_t_right"] != current["target_t_left"]:
            raise RuntimeError("localized final-gap rows have a coverage gap")

    failed = [row for row in rows if not row["passed"]]
    return {
        "rows": rows,
        "row_count": len(rows),
        "all_rows_passed": not failed,
        "failed_rows": len(failed),
        "stage_rows": {
            name: len(stage_rows)
            for name, stage_rows in hierarchy.items()
        },
        "proof_formula": (
            "w''=2s''-r''+phi(V)V''-chi(V)(V')^2; "
            "|Delta^2 s^(r)-s^(r+2)|<=sup|s^(r+4)|/12"
        ),
    }


def point_sixth_hierarchy(
    anchor: Fraction,
    order: int = POINT_ORDER,
    *,
    support_h_rows: list[dict] | None = None,
    point_h_source: dict[Fraction, tuple[list, dict]] | None = None,
) -> dict:
    """Build the rigorous sixth-layer Taylor jet at one exact t anchor."""
    internal_order = order + 2 if support_h_rows is not None else order
    h_source: dict[int, list] = {}
    quadrature = []
    for shift in range(-7, 8):
        target = anchor + shift
        if point_h_source is None:
            h_source[shift], diagnostics = point_h_series(
                target,
                internal_order,
            )
        else:
            try:
                cached_series, diagnostics = point_h_source[target]
            except KeyError as exc:
                raise RuntimeError(
                    f"point H cache misses target t={target}"
                ) from exc
            if len(cached_series) < internal_order + 1:
                raise RuntimeError(
                    f"point H cache at t={target} stops before order "
                    f"{internal_order}"
                )
            h_source[shift] = cached_series[: internal_order + 1]
        quadrature.append(diagnostics)

    b = {
        shift: centered_second_difference(h_source, shift)
        for shift in range(-6, 7)
    }
    ell = {
        shift: stable_log_series(b[shift], internal_order)
        for shift in range(-6, 7)
    }
    j_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 2),
            centered_second_difference(ell, shift),
        )
        for shift in range(-5, 6)
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
            stable_log_series(j_coordinate[shift], internal_order),
        )
        for shift in range(-5, 6)
    }
    r_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 3),
            centered_second_difference(h_stable, shift),
        )
        for shift in range(-4, 5)
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
            stable_log_series(r_coordinate[shift], internal_order),
        )
        for shift in range(-4, 5)
    }
    s_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 4),
            centered_second_difference(q, shift),
        )
        for shift in range(-3, 4)
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
            stable_log_series(s_coordinate[shift], internal_order),
        )
        for shift in range(-3, 4)
    }
    t_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 5),
            centered_second_difference(p, shift),
        )
        for shift in range(-2, 3)
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
            stable_log_series(t_coordinate[shift], internal_order),
        )
        for shift in range(-2, 3)
    }
    u_coordinate = {
        shift: series_sub(
            series_scale(b[shift], 6),
            centered_second_difference(r_final, shift),
        )
        for shift in range(-1, 2)
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
            stable_log_series(u_coordinate[shift], internal_order),
        )
        for shift in range(-1, 2)
    }
    if support_h_rows is None:
        v_coordinate = series_sub(
            series_scale(b[0], 7),
            centered_second_difference(s_final, 0),
        )[: order + 1]
        localization_diagnostics = None
    else:
        support_order = order + 4
        support_hierarchy = dimensionless_localized_hierarchy(
            anchor,
            support_h_rows,
            support_order,
        )
        support_rows = [
            row
            for row in support_hierarchy["U"]
            if row["target_t_right"] > anchor - 1
            and row["target_t_left"] < anchor + 1
        ]
        if not support_rows:
            raise RuntimeError(f"localized point V support misses t={anchor}")
        support_rows.sort(key=lambda row: row["target_t_left"])
        if (
            support_rows[0]["target_t_left"] > anchor - 1
            or support_rows[-1]["target_t_right"] < anchor + 1
        ):
            raise RuntimeError(f"localized point V support is incomplete at t={anchor}")
        scale = arb_rational(anchor)
        v_coordinate = []
        remainder_rows = []
        for degree in range(order + 1):
            local = (
                7 * b[0][degree]
                - (degree + 1)
                * (degree + 2)
                * s_final[0][degree + 2]
            )
            support_coefficient = hull(
                [row["s"][degree + 4] for row in support_rows]
            )
            derivative_remainder = (
                math.factorial(degree + 4)
                * upper_absolute(support_coefficient)
                / scale ** (degree + 4)
                / 12
            )
            coefficient_remainder = derivative_remainder / math.factorial(degree)
            v_coordinate.append(local + _symmetric(coefficient_remainder))
            remainder_rows.append(arb_upper_text(coefficient_remainder))
        localization_diagnostics = {
            "support_order": support_order,
            "support_rows": len(support_rows),
            "support_t_left": str(support_rows[0]["target_t_left"]),
            "support_t_right": str(support_rows[-1]["target_t_right"]),
            "coefficient_remainder_upper": remainder_rows,
        }
    point_floor = inherited_coordinate_floor(
        anchor,
        7,
        Fraction(2500),
    )
    elementary_floor = analytic_v_floor(anchor)
    if bool(elementary_floor > point_floor):
        point_floor = elementary_floor
    raw_v_upper = _upper_point(v_coordinate[0])
    raw_v_lower = _lower_point(v_coordinate[0])
    improved_v_lower = (
        point_floor if bool(point_floor > raw_v_lower) else raw_v_lower
    )
    if not bool(raw_v_upper >= improved_v_lower):
        raise RuntimeError(f"analytic point V floor is disjoint at t={anchor}")
    v_coordinate[0] = _interval_from_endpoints(
        improved_v_lower,
        raw_v_upper,
    )
    w_final = series_add(
        series_sub(
            series_scale(s_final[0][: order + 1], 2),
            r_final[0][: order + 1],
        ),
        stable_log_series(v_coordinate, order),
    )
    return {
        "B": b[0],
        "J": j_coordinate[0],
        "R": r_coordinate[0],
        "S": s_coordinate[0],
        "T": t_coordinate[0],
        "U": u_coordinate[0],
        "V": v_coordinate,
        "r": r_final[0],
        "s": s_final[0],
        "w": w_final,
        "quadrature": quadrature,
        "V_localization": localization_diagnostics,
    }


def dimensionless_localized_sixth_hierarchy(
    scale_t: Fraction,
    h_rows: list[dict],
    output_order: int = REMAINDER_ORDER,
) -> dict:
    """Append a localized sixth stable stage to the five-stage hierarchy."""
    if output_order < 2:
        raise ValueError("sixth-layer output order must be at least two")
    first_five_order = output_order + 2
    hierarchy = dimensionless_localized_hierarchy(
        scale_t,
        h_rows,
        first_five_order,
    )
    scale = arb_rational(scale_t)
    u_rows = hierarchy["U"]
    convolvers = {
        degree: ExactTileHullConvolver(
            u_rows,
            [row["s"][degree] for row in u_rows],
        )
        for degree in range(2, output_order + 3)
    }
    v_rows = []
    for row in u_rows:
        left = row["target_t_left"]
        right = row["target_t_right"]
        averages = [
            convolvers[degree + 2].average(left, right)
            for degree in range(output_order + 1)
        ]
        if any(value is None for value in averages):
            continue
        coordinate = [
            7 * row["B"][degree]
            - (degree + 1) * (degree + 2) * value / scale
            for degree, value in enumerate(averages)
            if value is not None
        ]
        if left < MINIMUM_T:
            continue
        normalized_floor = scale * analytic_v_floor(right)
        raw_upper = _upper_point(coordinate[0])
        raw_lower = _lower_point(coordinate[0])
        improved_lower = (
            normalized_floor
            if bool(normalized_floor > raw_lower)
            else raw_lower
        )
        if not bool(raw_upper >= improved_lower):
            raise RuntimeError(
                f"analytic V floor is disjoint from localized V on {left}..{right}"
            )
        coordinate[0] = _interval_from_endpoints(improved_lower, raw_upper)
        if not bool(coordinate[0] > 0):
            raise RuntimeError(f"nonpositive localized V on {left}..{right}")
        log_v = stable_log_series(
            [coefficient / scale for coefficient in coordinate],
            output_order,
        )
        w = series_add(
            series_sub(
                series_scale(row["s"][: output_order + 1], 2),
                row["r"][: output_order + 1],
            ),
            log_v,
        )
        v_rows.append({**row, "V": coordinate, "w": w})
    return {**hierarchy, "V": v_rows}


def localized_sixth_remainder_bounds(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    remainder_order: int = REMAINDER_ORDER,
) -> dict:
    """Hull the normalized sixth-layer remainder on a central block."""
    hierarchy = dimensionless_localized_sixth_hierarchy(
        anchor,
        h_rows,
        remainder_order,
    )
    central = [
        row
        for row in hierarchy["V"]
        if row["target_t_right"] > anchor
        and row["target_t_left"] < right
    ]
    if not central:
        raise RuntimeError(f"localized sixth hierarchy misses {anchor}..{right}")
    central.sort(key=lambda row: row["target_t_left"])
    if central[0]["target_t_left"] > anchor or central[-1]["target_t_right"] < right:
        raise RuntimeError("localized sixth hierarchy has incomplete central coverage")
    for previous, current in zip(central, central[1:]):
        if previous["target_t_right"] != current["target_t_left"]:
            raise RuntimeError("localized sixth hierarchy has a central gap")
    coordinates = {
        name: hull([row[name][0] for row in central])
        for name in ("B", "J", "R", "S", "T", "U", "V")
    }
    for name, value in coordinates.items():
        if not bool(value > 0):
            raise RuntimeError(
                f"nonpositive localized common {name} on {anchor}..{right}"
            )
    return {
        "coordinates": coordinates,
        "w_remainder_coefficient": hull(
            [row["w"][remainder_order] for row in central]
        ),
        "stage_rows": {
            name: len(rows)
            for name, rows in hierarchy.items()
        },
        "central_rows": len(central),
        "central_t_left": central[0]["target_t_left"],
        "central_t_right": central[-1]["target_t_right"],
    }


def localized_sixth_continuation_row(
    anchor: Fraction,
    right: Fraction,
    h_rows: list[dict],
    *,
    point_order: int = POINT_ORDER,
    remainder_order: int = REMAINDER_ORDER,
    point_h_source: dict[Fraction, tuple[list, dict]] | None = None,
) -> dict:
    """Certify one rightward block by a point jet and localized remainder."""
    flint.ctx.prec = PRECISION_BITS
    if not MINIMUM_T <= anchor < right:
        raise ValueError("invalid order-nine continuation block")
    if remainder_order != point_order + 1:
        raise ValueError("remainder order must be one above point order")
    if analytic_v_floor_polynomial(anchor) <= 0:
        raise ValueError("continuous sixth-gap floor polynomial is not positive")

    common = localized_sixth_remainder_bounds(
        anchor,
        right,
        h_rows,
        remainder_order,
    )
    point = point_sixth_hierarchy(
        anchor,
        point_order,
        support_h_rows=h_rows,
        point_h_source=point_h_source,
    )
    coordinate_rows = {}
    for name in ("B", "J", "R", "S", "T", "U", "V"):
        if not bool(point[name][0] > 0):
            raise RuntimeError(f"nonpositive point {name} at t={anchor}")
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
    point_polynomial = 2 * point["w"][2]
    for degree in range(3, point_order + 1):
        point_polynomial += (
            upper_absolute(point["w"][degree])
            * degree
            * (degree - 1)
            * width ** (degree - 2)
        )
    normalized_remainder = upper_absolute(
        common["w_remainder_coefficient"]
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
            f"order-nine localized continuation failed on {anchor}..{right}: "
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
            for name in ("B", "J", "R", "S", "T", "U", "V")
        },
        "common_coordinates": coordinate_rows,
        "point_scaled_curvature": arb_upper_text(
            arb_rational(anchor) ** 2 * 2 * point["w"][2]
        ),
        "point_polynomial_curvature_upper": arb_upper_text(point_polynomial),
        "normalized_w_remainder_coefficient_upper": arb_upper_text(
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
