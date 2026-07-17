#!/usr/bin/env python3
"""Interval core for the fourth-nested first-summand order-seven curvature."""

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

from jensen_window_pf_compound_order5_nested_curvature_interval_core import (  # noqa: E402
    series_add,
    series_exp,
    series_log,
    series_scale,
    series_sub,
    stable_log_series,
)
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_lower_text,
    arb_rational,
    arb_upper_text,
)


CURVATURE_CONSTANT = 600
MAX_JET_ORDER = 10


def nested_curvature_from_h_derivatives(
    derivatives: dict[int, flint.arb],
) -> dict:
    """Enclose the fourth stable layer from a common H^(2)..H^(12) cover.

    Each centered second difference consumes two Taylor coefficients. The
    four stable logarithms after log(1-exp(-B)) therefore require a common
    t+-5 H cover and derivatives through order twelve.
    """
    missing = [order for order in range(2, 13) if order not in derivatives]
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
        for order in range(9)
    ]
    h = series_add(series_scale(ell[:9], 2), stable_log_series(J, 8))

    R = [
        3 * B[order]
        - h[order + 2]
        * math.factorial(order + 2)
        / math.factorial(order)
        for order in range(7)
    ]
    q = series_add(
        series_sub(series_scale(h[:7], 2), ell[:7]),
        stable_log_series(R, 6),
    )

    S = [
        4 * B[order]
        - q[order + 2]
        * math.factorial(order + 2)
        / math.factorial(order)
        for order in range(5)
    ]
    p = series_add(
        series_sub(series_scale(q[:5], 2), h[:5]),
        stable_log_series(S, 4),
    )

    T = [
        5 * B[order]
        - p[order + 2]
        * math.factorial(order + 2)
        / math.factorial(order)
        for order in range(3)
    ]
    r = series_add(
        series_sub(series_scale(p[:3], 2), q[:3]),
        stable_log_series(T, 2),
    )
    r_second = 2 * r[2]
    return {
        "B": B,
        "d": d,
        "ell": ell,
        "J": J,
        "h": h,
        "R": R,
        "q": q,
        "S": S,
        "p": p,
        "T": T,
        "r": r,
        "r_second": r_second,
    }


def evaluate_nested_curvature_from_h_cover(
    left,
    right,
    derivatives: dict[int, flint.arb],
    *,
    cover_diagnostics: dict | None = None,
) -> dict:
    """Certify r_1''(t)<=600/t^2 on a central mode interval."""
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
    coordinates = {
        "J": jets["J"][0],
        "R": jets["R"][0],
        "S": jets["S"][0],
        "T": jets["T"][0],
    }
    r_second = jets["r_second"]
    if not bool(B[0] > 0 and d[0] > 0):
        return {"passed": False, "failure": "nonpositive-B-or-d"}
    for name, coordinate in coordinates.items():
        if not bool(coordinate > 0):
            return {
                "passed": False,
                "failure": f"nonpositive-{name}",
                "coordinate_ball": coordinate.str(40).replace("e", "E"),
            }

    target = flint.arb(CURVATURE_CONSTANT) / t**2
    margin = target - r_second
    scaled = t**2 * r_second
    if not bool(margin > 0):
        return {
            "passed": False,
            "failure": "curvature-margin",
            "mode": [str(left), str(right)],
            "t_ball": t.str(40).replace("e", "E"),
            "r_second_ball": r_second.str(40).replace("e", "E"),
            "target_ball": target.str(40).replace("e", "E"),
            "scaled_ball": scaled.str(40).replace("e", "E"),
            **{
                f"{name}_ball": value.str(40).replace("e", "E")
                for name, value in coordinates.items()
            },
        }
    return {
        "passed": True,
        "failure": None,
        "central_mode": [str(left), str(right)],
        "central_t_ball": t.str(40).replace("e", "E"),
        "B_lower": arb_lower_text(B[0]),
        "d_lower": arb_lower_text(d[0]),
        **{
            f"{name}_lower": arb_lower_text(value)
            for name, value in coordinates.items()
        },
        "r_second_upper": arb_upper_text(r_second),
        "target_lower": arb_lower_text(target),
        "margin_lower": arb_lower_text(margin),
        "scaled_curvature_upper": arb_upper_text(scaled),
        "cover_diagnostics": cover_diagnostics,
    }
