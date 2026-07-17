#!/usr/bin/env python3
"""Localized tent-convolution core for the order-seven curvature."""

from __future__ import annotations

from bisect import bisect_left
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
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    exact_lower,
    exact_upper,
)


CURVATURE_CONSTANT = 600
TENT_RADIUS = 1


def interval_from_bounds(lower: flint.arb, upper: flint.arb) -> flint.arb:
    """Return an outward ball containing the supplied ordered bounds."""
    lo = exact_lower(lower)
    hi = exact_upper(upper)
    if bool(hi < lo):
        raise ValueError("reversed interval bounds")
    return (lo + hi) / 2 + flint.arb(0, (hi - lo) / 2)


def _tent_cdf_point(value: flint.arb) -> flint.arb:
    """Evaluate the unit-tent CDF at an effectively point-valued ball."""
    if bool(value <= -1):
        return flint.arb(0)
    if bool(value < 0):
        return (value + 1) ** 2 / 2
    if bool(value < 1):
        return 1 - (1 - value) ** 2 / 2
    if bool(value >= 1):
        return flint.arb(1)
    raise ValueError("tent CDF point straddles a branch boundary")


def tent_cdf(value: flint.arb) -> flint.arb:
    """Enclose the monotone unit-tent CDF on an interval."""
    lower = _tent_cdf_point(exact_lower(value))
    upper = _tent_cdf_point(exact_upper(value))
    return interval_from_bounds(lower, upper)


def tent_mass(
    tile_left: flint.arb,
    tile_right: flint.arb,
    center: flint.arb,
) -> flint.arb:
    """Enclose the tent mass of one uncertain t tile for every center."""
    a_lo = exact_lower(tile_left)
    a_hi = exact_upper(tile_left)
    b_lo = exact_lower(tile_right)
    b_hi = exact_upper(tile_right)
    c_lo = exact_lower(center)
    c_hi = exact_upper(center)
    if not bool(a_hi < b_lo):
        raise ValueError("t tile endpoints are not strictly ordered")

    def mass_at(a: flint.arb, b: flint.arb, c: flint.arb) -> flint.arb:
        return tent_cdf(b - c) - tent_cdf(a - c)

    lower_candidates = (
        exact_lower(mass_at(a_hi, b_lo, c_lo)),
        exact_lower(mass_at(a_hi, b_lo, c_hi)),
    )
    lower = min(lower_candidates, key=lambda value: value.lower())
    if bool(lower < 0):
        lower = flint.arb(0)

    midpoint = (a_lo + b_hi) / 2
    if bool(midpoint < c_lo):
        maximizing_center = c_lo
    elif bool(midpoint > c_hi):
        maximizing_center = c_hi
    else:
        maximizing_center = midpoint
    upper = exact_upper(mass_at(a_lo, b_hi, maximizing_center))
    tile_width = exact_upper(b_hi - a_lo)
    if bool(upper > tile_width):
        upper = tile_width
    if bool(upper < lower):
        raise ValueError("tent mass bounds are inconsistent")
    return interval_from_bounds(lower, upper)


def _midpoint_key(value: flint.arb) -> float:
    return float(value.mid())


def support_slice(rows: list[dict], center: flint.arb) -> tuple[int, int] | None:
    """Locate a contiguous tile slice covering center+-1."""
    if not rows:
        return None
    lower = center - TENT_RADIUS
    upper = center + TENT_RADIUS
    right_keys = [_midpoint_key(row["t_right"]) for row in rows]
    left = max(0, bisect_left(right_keys, _midpoint_key(lower)) - 1)
    while left > 0 and not bool(rows[left]["t_left"] < lower):
        left -= 1
    while left + 1 < len(rows) and bool(rows[left]["t_right"] < lower):
        left += 1
    if not bool(rows[left]["t_left"] < lower):
        return None

    right = max(left, bisect_left(right_keys, _midpoint_key(upper)))
    right = min(right, len(rows) - 1)
    while right + 1 < len(rows) and not bool(rows[right]["t_right"] > upper):
        right += 1
    if not bool(rows[right]["t_right"] > upper):
        return None
    return left, right + 1


def tent_masses(rows: list[dict], center: flint.arb) -> list[tuple[dict, flint.arb]] | None:
    """Return rigorous tile masses for a complete unit-tent support."""
    support = support_slice(rows, center)
    if support is None:
        return None
    left, right = support
    selected = rows[left:right]
    masses = [
        (row, tent_mass(row["t_left"], row["t_right"], center))
        for row in selected
    ]
    return masses


def tent_average(
    masses: list[tuple[dict, flint.arb]],
    values,
) -> flint.arb:
    """Enclose a tent average while cancelling its exact unit mass."""
    source_values = [values(row) for row, _ in masses]
    value_hull = hull(source_values)
    reference = flint.arb(value_hull.mid())
    return reference + sum(
        ((value - reference) * mass for value, (_, mass) in zip(source_values, masses)),
        flint.arb(0),
    )


class TentConvolver:
    """Fast rigorous convolution of one interval field with the unit tent."""

    def __init__(
        self,
        rows: list[dict],
        values: list[flint.arb],
        *,
        anchor: int = 320,
    ) -> None:
        if len(rows) != len(values):
            raise ValueError("tent field length mismatch")
        self.rows = rows
        self.left_keys = [_midpoint_key(row["t_left"]) for row in rows]
        self.right_keys = [_midpoint_key(row["t_right"]) for row in rows]
        value_hull = hull(values)
        self.reference = flint.arb(value_hull.mid())
        self.anchor = flint.arb(anchor)
        self.prefix_mass = [flint.arb(0)]
        self.prefix_first = [flint.arb(0)]
        for row, value in zip(rows, values):
            width = row["t_right"] - row["t_left"]
            midpoint = (row["t_left"] + row["t_right"]) / 2
            deviation_mass = (value - self.reference) * width
            self.prefix_mass.append(self.prefix_mass[-1] + deviation_mass)
            self.prefix_first.append(
                self.prefix_first[-1]
                + deviation_mass * (midpoint - self.anchor)
            )

    def _first_left_strictly_after(self, bound: flint.arb) -> int:
        index = bisect_left(self.left_keys, _midpoint_key(bound))
        while index > 0 and bool(self.rows[index - 1]["t_left"] > bound):
            index -= 1
        while index < len(self.rows) and not bool(
            self.rows[index]["t_left"] > bound
        ):
            index += 1
        return index

    def _first_right_not_strictly_before(self, bound: flint.arb) -> int:
        index = bisect_left(self.right_keys, _midpoint_key(bound))
        while index > 0 and not bool(self.rows[index - 1]["t_right"] < bound):
            index -= 1
        while index < len(self.rows) and bool(
            self.rows[index]["t_right"] < bound
        ):
            index += 1
        return index

    def _prefix_range(self, left: int, right: int) -> tuple[flint.arb, flint.arb]:
        return (
            self.prefix_mass[right] - self.prefix_mass[left],
            self.prefix_first[right] - self.prefix_first[left],
        )

    def _support_slice(self, center: flint.arb) -> tuple[int, int] | None:
        lower = center - TENT_RADIUS
        upper = center + TENT_RADIUS
        left = max(0, bisect_left(self.right_keys, _midpoint_key(lower)) - 1)
        while left > 0 and not bool(self.rows[left]["t_left"] < lower):
            left -= 1
        while left + 1 < len(self.rows) and bool(
            self.rows[left]["t_right"] < lower
        ):
            left += 1
        if not bool(self.rows[left]["t_left"] < lower):
            return None
        right = max(left, bisect_left(self.right_keys, _midpoint_key(upper)))
        right = min(right, len(self.rows) - 1)
        while right + 1 < len(self.rows) and not bool(
            self.rows[right]["t_right"] > upper
        ):
            right += 1
        if not bool(self.rows[right]["t_right"] > upper):
            return None
        return left, right + 1

    def average(self, center: flint.arb) -> flint.arb | None:
        """Enclose the tent average for every center in the supplied ball."""
        support = self._support_slice(center)
        if support is None:
            return None
        support_left, support_right = support
        center_lo = exact_lower(center)
        center_hi = exact_upper(center)

        left_start = max(
            support_left,
            self._first_left_strictly_after(center_hi - 1),
        )
        left_end = min(
            support_right,
            self._first_right_not_strictly_before(center_lo),
        )
        right_start = max(
            support_left,
            self._first_left_strictly_after(center_hi),
        )
        right_end = min(
            support_right,
            self._first_right_not_strictly_before(center_lo + 1),
        )
        full_ranges = [
            (left_start, max(left_start, left_end), "left"),
            (right_start, max(right_start, right_end), "right"),
        ]
        full_ranges = [
            item
            for item in full_ranges
            if item[0] < item[1]
        ]
        full_ranges.sort()

        center_shift = center - self.anchor
        total = self.reference
        for left, right, side in full_ranges:
            mass_sum, first_sum = self._prefix_range(left, right)
            if side == "left":
                total += first_sum + (1 - center_shift) * mass_sum
            else:
                total += -first_sum + (1 + center_shift) * mass_sum

        cursor = support_left
        boundary_indices: list[int] = []
        for left, right, _ in full_ranges:
            boundary_indices.extend(range(cursor, left))
            cursor = max(cursor, right)
        boundary_indices.extend(range(cursor, support_right))
        for index in boundary_indices:
            row = self.rows[index]
            mass = tent_mass(row["t_left"], row["t_right"], center)
            value = (
                self.prefix_mass[index + 1] - self.prefix_mass[index]
            ) / (row["t_right"] - row["t_left"])
            total += value * mass
        return total


class RangeHull:
    """Static exact-endpoint sparse table for interval range hulls."""

    def __init__(self, values: list[flint.arb]) -> None:
        if not values:
            raise ValueError("empty range-hull field")
        self.lower_levels = [[value.lower() for value in values]]
        self.upper_levels = [[value.upper() for value in values]]
        span = 2
        while span <= len(values):
            half = span // 2
            previous_lower = self.lower_levels[-1]
            previous_upper = self.upper_levels[-1]
            count = len(values) - span + 1
            self.lower_levels.append(
                [
                    min(previous_lower[index], previous_lower[index + half])
                    for index in range(count)
                ]
            )
            self.upper_levels.append(
                [
                    max(previous_upper[index], previous_upper[index + half])
                    for index in range(count)
                ]
            )
            span *= 2

    def query(self, left: int, right: int) -> flint.arb:
        if not 0 <= left < right:
            raise ValueError("invalid range-hull query")
        length = right - left
        level = length.bit_length() - 1
        span = 1 << level
        lower = min(
            self.lower_levels[level][left],
            self.lower_levels[level][right - span],
        )
        upper = max(
            self.upper_levels[level][left],
            self.upper_levels[level][right - span],
        )
        return interval_from_bounds(flint.arb(lower), flint.arb(upper))


class LocalizedHullConvolver:
    """Unit-tent enclosure by a localized exact endpoint range hull."""

    def __init__(self, rows: list[dict], values: list[flint.arb]) -> None:
        if len(rows) != len(values):
            raise ValueError("localized hull field length mismatch")
        self.rows = rows
        self.right_keys = [_midpoint_key(row["t_right"]) for row in rows]
        self.range_hull = RangeHull(values)

    def _support_slice(self, center: flint.arb) -> tuple[int, int] | None:
        lower = center - TENT_RADIUS
        upper = center + TENT_RADIUS
        left = max(0, bisect_left(self.right_keys, _midpoint_key(lower)) - 1)
        while left > 0 and not bool(self.rows[left]["t_left"] < lower):
            left -= 1
        while left + 1 < len(self.rows) and bool(
            self.rows[left]["t_right"] < lower
        ):
            left += 1
        if not bool(self.rows[left]["t_left"] < lower):
            return None
        right = max(left, bisect_left(self.right_keys, _midpoint_key(upper)))
        right = min(right, len(self.rows) - 1)
        while right + 1 < len(self.rows) and not bool(
            self.rows[right]["t_right"] > upper
        ):
            right += 1
        if not bool(self.rows[right]["t_right"] > upper):
            return None
        return left, right + 1

    def average(self, center: flint.arb) -> flint.arb | None:
        support = self._support_slice(center)
        if support is None:
            return None
        return self.range_hull.query(*support)


def build_b_rows(h_rows: list[dict]) -> list[dict]:
    """Build B and ell jets through order ten from H^(2)..H^(12)."""
    convolvers = {
        order: LocalizedHullConvolver(
            h_rows,
            [row["H"][order] for row in h_rows],
        )
        for order in range(2, 13)
    }
    output: list[dict] = []
    for row in h_rows:
        center = interval_from_bounds(row["t_left"], row["t_right"])
        averages = [
            convolvers[order + 2].average(center)
            for order in range(11)
        ]
        if any(value is None for value in averages):
            continue
        b = [
            value / math.factorial(order)
            for order, value in enumerate(averages)
            if value is not None
        ]
        if not bool(b[0] > 0):
            raise RuntimeError(f"nonpositive B on mode tile {row['mode_left']}")
        output.append({**row, "B": b, "ell": stable_log_series(b, 10)})
    return output


def build_stable_stage(
    rows: list[dict],
    *,
    coordinate_name: str,
    output_name: str,
    current_name: str,
    previous_name: str | None,
    leading_factor: int,
    output_order: int,
) -> list[dict]:
    """Apply one localized centered-difference stable layer."""
    convolvers = {
        order: LocalizedHullConvolver(
            rows,
            [
                row[current_name][order] * math.factorial(order)
                for row in rows
            ],
        )
        for order in range(2, output_order + 3)
    }
    output: list[dict] = []
    for row in rows:
        center = interval_from_bounds(row["t_left"], row["t_right"])
        averages = [
            convolvers[order + 2].average(center)
            for order in range(output_order + 1)
        ]
        if any(value is None for value in averages):
            continue
        coordinate = [
            leading_factor * row["B"][order]
            - value / math.factorial(order)
            for order, value in enumerate(averages)
            if value is not None
        ]
        if not bool(coordinate[0] > 0):
            raise RuntimeError(
                f"nonpositive {coordinate_name} on mode tile {row['mode_left']}: "
                f"{coordinate[0]}"
            )
        stable = stable_log_series(coordinate, output_order)
        base = series_scale(row[current_name][: output_order + 1], 2)
        if previous_name is not None:
            base = series_sub(base, row[previous_name][: output_order + 1])
        new_value = series_add(base, stable)
        output.append(
            {
                **row,
                coordinate_name: coordinate,
                output_name: new_value,
            }
        )
    return output


def build_order7_hierarchy(h_rows: list[dict]) -> dict[str, list[dict]]:
    """Build the four localized stable layers ending in r jets."""
    b_rows = build_b_rows(h_rows)
    j_rows = build_stable_stage(
        b_rows,
        coordinate_name="J",
        output_name="h",
        current_name="ell",
        previous_name=None,
        leading_factor=2,
        output_order=8,
    )
    r_rows = build_stable_stage(
        j_rows,
        coordinate_name="R",
        output_name="q",
        current_name="h",
        previous_name="ell",
        leading_factor=3,
        output_order=6,
    )
    s_rows = build_stable_stage(
        r_rows,
        coordinate_name="S",
        output_name="p",
        current_name="q",
        previous_name="h",
        leading_factor=4,
        output_order=4,
    )
    t_rows = build_stable_stage(
        s_rows,
        coordinate_name="T",
        output_name="r_jet",
        current_name="p",
        previous_name="q",
        leading_factor=5,
        output_order=2,
    )
    return {
        "B": b_rows,
        "J": j_rows,
        "R": r_rows,
        "S": s_rows,
        "T": t_rows,
    }


def build_taylor_b_rows(h_rows: list[dict]) -> list[dict]:
    """Build B/ell through order eighteen with a symmetric Taylor tent bound."""
    correction_hulls = {
        order: LocalizedHullConvolver(
            h_rows,
            [row["H"][order] for row in h_rows],
        )
        for order in range(4, 23)
    }
    output: list[dict] = []
    for row in h_rows:
        center = interval_from_bounds(row["t_left"], row["t_right"])
        corrections = [
            correction_hulls[order + 4].average(center)
            for order in range(19)
        ]
        if any(value is None for value in corrections):
            continue
        b = [
            (
                row["H"][order + 2]
                + correction / 12
            )
            / math.factorial(order)
            for order, correction in enumerate(corrections)
            if correction is not None
        ]
        if not bool(b[0] > 0):
            raise RuntimeError(f"nonpositive Taylor B on {row['mode_left']}")
        output.append({**row, "B": b, "ell": stable_log_series(b, 18)})
    return output


def build_taylor_stable_stage(
    rows: list[dict],
    *,
    coordinate_name: str,
    output_name: str,
    current_name: str,
    previous_name: str | None,
    leading_factor: int,
    output_order: int,
) -> list[dict]:
    """Apply one stable layer using T(f) in f+Hull(f'')/12."""
    correction_hulls = {
        order: LocalizedHullConvolver(
            rows,
            [
                row[current_name][order] * math.factorial(order)
                for row in rows
            ],
        )
        for order in range(4, output_order + 5)
    }
    output: list[dict] = []
    for row in rows:
        center = interval_from_bounds(row["t_left"], row["t_right"])
        corrections = [
            correction_hulls[order + 4].average(center)
            for order in range(output_order + 1)
        ]
        if any(value is None for value in corrections):
            continue
        coordinate = [
            leading_factor * row["B"][order]
            - (
                row[current_name][order + 2] * math.factorial(order + 2)
                + correction / 12
            )
            / math.factorial(order)
            for order, correction in enumerate(corrections)
            if correction is not None
        ]
        if not bool(coordinate[0] > 0):
            raise RuntimeError(
                f"nonpositive Taylor {coordinate_name} on "
                f"{row['mode_left']}: {coordinate[0]}"
            )
        stable = stable_log_series(coordinate, output_order)
        base = series_scale(row[current_name][: output_order + 1], 2)
        if previous_name is not None:
            base = series_sub(base, row[previous_name][: output_order + 1])
        output.append(
            {
                **row,
                coordinate_name: coordinate,
                output_name: series_add(base, stable),
            }
        )
    return output


def build_order7_taylor_hierarchy(h_rows: list[dict]) -> dict[str, list[dict]]:
    """Build the order-seven hierarchy with second-order tent localization."""
    b_rows = build_taylor_b_rows(h_rows)
    j_rows = build_taylor_stable_stage(
        b_rows,
        coordinate_name="J",
        output_name="h",
        current_name="ell",
        previous_name=None,
        leading_factor=2,
        output_order=14,
    )
    r_rows = build_taylor_stable_stage(
        j_rows,
        coordinate_name="R",
        output_name="q",
        current_name="h",
        previous_name="ell",
        leading_factor=3,
        output_order=10,
    )
    s_rows = build_taylor_stable_stage(
        r_rows,
        coordinate_name="S",
        output_name="p",
        current_name="q",
        previous_name="h",
        leading_factor=4,
        output_order=6,
    )
    t_rows = build_taylor_stable_stage(
        s_rows,
        coordinate_name="T",
        output_name="r_jet",
        current_name="p",
        previous_name="q",
        leading_factor=5,
        output_order=2,
    )
    return {
        "B": b_rows,
        "J": j_rows,
        "R": r_rows,
        "S": s_rows,
        "T": t_rows,
    }


def curvature_row(row: dict) -> dict:
    """Evaluate the order-seven curvature ceiling on one final tile."""
    t = interval_from_bounds(row["t_left"], row["t_right"])
    r_second = 2 * row["r_jet"][2]
    scaled = t**2 * r_second
    margin = flint.arb(CURVATURE_CONSTANT) - scaled
    return {
        "passed": bool(margin > 0),
        "t": t,
        "r_second": r_second,
        "scaled": scaled,
        "margin": margin,
        "J": row["J"][0],
        "R": row["R"][0],
        "S": row["S"][0],
        "T": row["T"][0],
    }
