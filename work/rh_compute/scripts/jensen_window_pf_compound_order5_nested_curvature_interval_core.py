#!/usr/bin/env python3
"""Interval core for the nested first-summand order-five curvature bound."""

from __future__ import annotations

import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
)


CURVATURE_CONSTANT = 60
MAX_JET_ORDER = 6


def series_add(left: list[flint.arb], right: list[flint.arb]) -> list[flint.arb]:
    if len(left) != len(right):
        raise ValueError("series lengths differ")
    return [a + b for a, b in zip(left, right)]


def series_sub(left: list[flint.arb], right: list[flint.arb]) -> list[flint.arb]:
    if len(left) != len(right):
        raise ValueError("series lengths differ")
    return [a - b for a, b in zip(left, right)]


def series_scale(values: list[flint.arb], scalar: int | flint.arb) -> list[flint.arb]:
    return [scalar * value for value in values]


def series_mul(
    left: list[flint.arb], right: list[flint.arb], order: int
) -> list[flint.arb]:
    if len(left) <= order or len(right) <= order:
        raise ValueError("series is too short for multiplication")
    return [
        sum(
            (left[index] * right[degree - index] for index in range(degree + 1)),
            flint.arb(0),
        )
        for degree in range(order + 1)
    ]


def series_inverse(values: list[flint.arb], order: int) -> list[flint.arb]:
    if len(values) <= order:
        raise ValueError("series is too short for inversion")
    if not bool(values[0] > 0 or values[0] < 0):
        raise ValueError("series constant term contains zero")
    result = [1 / values[0]]
    for degree in range(1, order + 1):
        result.append(
            -sum(
                (
                    values[index] * result[degree - index]
                    for index in range(1, degree + 1)
                ),
                flint.arb(0),
            )
            / values[0]
        )
    return result


def series_log(values: list[flint.arb], order: int) -> list[flint.arb]:
    if len(values) <= order:
        raise ValueError("series is too short for logarithm")
    if not bool(values[0] > 0):
        raise ValueError("logarithm constant term is not positive")
    if order == 0:
        return [values[0].log()]
    derivative = [
        (degree + 1) * values[degree + 1] for degree in range(order)
    ]
    quotient = series_mul(
        derivative, series_inverse(values, order - 1), order - 1
    )
    return [values[0].log()] + [
        quotient[degree - 1] / degree for degree in range(1, order + 1)
    ]


def series_exp(values: list[flint.arb], order: int) -> list[flint.arb]:
    if len(values) <= order:
        raise ValueError("series is too short for exponential")
    result = [values[0].exp()]
    for degree in range(1, order + 1):
        result.append(
            sum(
                (
                    index * values[index] * result[degree - index]
                    for index in range(1, degree + 1)
                ),
                flint.arb(0),
            )
            / degree
        )
    return result


def stable_log_series(values: list[flint.arb], order: int) -> list[flint.arb]:
    """Return the Taylor jet of log(1-exp(-values))."""
    exponential = series_exp(series_scale(values[: order + 1], -1), order)
    one = [flint.arb(1)] + [flint.arb(0) for _ in range(order)]
    return series_log(series_sub(one, exponential), order)


def interval_text(value: flint.arb) -> list[str]:
    return [arb_lower_text(value), arb_upper_text(value)]


def nested_curvature_from_h_derivatives(
    derivatives: dict[int, flint.arb],
) -> dict:
    """Enclose the nested stable jets from a common H^(2)..H^(8) cover.

    Each B^(r) is a unit-mass tent average of H^(r+2). The same common
    derivative hull covers the successively shifted ell and h averages, so
    replacing each tent average by that hull is inclusion monotone.
    """
    missing = [order for order in range(2, 9) if order not in derivatives]
    if missing:
        raise ValueError(f"missing H derivatives: {missing}")

    B = [
        derivatives[order + 2] / math.factorial(order)
        for order in range(MAX_JET_ORDER + 1)
    ]
    one = [flint.arb(1)] + [flint.arb(0) for _ in range(MAX_JET_ORDER)]
    x = series_exp(series_scale(B, -1), MAX_JET_ORDER)
    d = series_sub(one, x)
    ell = series_log(d, MAX_JET_ORDER)

    J = [
        2 * B[order]
        - ell[order + 2]
        * math.factorial(order + 2)
        / math.factorial(order)
        for order in range(5)
    ]
    stable_J = stable_log_series(J, 4)
    h = series_add(series_scale(ell[:5], 2), stable_J)

    R = [
        3 * B[order]
        - h[order + 2]
        * math.factorial(order + 2)
        / math.factorial(order)
        for order in range(3)
    ]
    stable_R = stable_log_series(R, 2)
    q = series_add(series_sub(series_scale(h[:3], 2), ell[:3]), stable_R)
    q_second = 2 * q[2]

    return {
        "B": B,
        "d": d,
        "ell": ell,
        "J": J,
        "h": h,
        "R": R,
        "q": q,
        "q_second": q_second,
    }


def evaluate_nested_curvature_from_h_cover(
    left,
    right,
    derivatives: dict[int, flint.arb],
    *,
    cover_diagnostics: dict | None = None,
) -> dict:
    """Certify q_1''(t)<=60/t^2 on a central mode interval."""
    mode = (
        (arb_rational(left) + arb_rational(right)) / 2
        + flint.arb(0, (arb_rational(right) - arb_rational(left)) / 2)
    )
    t = potential_jet_arb(mode, 1)[1]
    if not bool(t > 0):
        return {"passed": False, "failure": "nonpositive-t"}
    try:
        jets = nested_curvature_from_h_derivatives(derivatives)
    except Exception as exc:
        return {
            "passed": False,
            "failure": "nested-jet-exception",
            "detail": repr(exc),
        }

    B = jets["B"]
    d = jets["d"]
    J = jets["J"]
    R = jets["R"]
    q_second = jets["q_second"]
    if not bool(B[0] > 0 and d[0] > 0):
        return {"passed": False, "failure": "nonpositive-B-or-d"}
    if not bool(J[0] > 0):
        return {"passed": False, "failure": "nonpositive-J"}
    if not bool(R[0] > 0):
        return {"passed": False, "failure": "nonpositive-R"}

    target = flint.arb(CURVATURE_CONSTANT) / t**2
    margin = target - q_second
    scaled = t**2 * q_second
    if not bool(margin > 0):
        return {
            "passed": False,
            "failure": "curvature-margin",
            "mode": [str(left), str(right)],
            "t_ball": t.str(40).replace("e", "E"),
            "q_second_ball": q_second.str(40).replace("e", "E"),
            "target_ball": target.str(40).replace("e", "E"),
            "scaled_ball": scaled.str(40).replace("e", "E"),
            "J_ball": J[0].str(40).replace("e", "E"),
            "R_ball": R[0].str(40).replace("e", "E"),
        }
    return {
        "passed": True,
        "failure": None,
        "central_mode": [str(left), str(right)],
        "central_t_ball": t.str(40).replace("e", "E"),
        "B_lower": arb_lower_text(B[0]),
        "d_lower": arb_lower_text(d[0]),
        "J_lower": arb_lower_text(J[0]),
        "R_lower": arb_lower_text(R[0]),
        "q_second_upper": arb_upper_text(q_second),
        "target_lower": arb_lower_text(target),
        "margin_lower": arb_lower_text(margin),
        "scaled_curvature_upper": arb_upper_text(scaled),
        "cover_diagnostics": cover_diagnostics,
    }
