#!/usr/bin/env python3
"""Interval core for the seventh-nested first-summand order-ten curvature.

The stable coordinates are

    W(t) = 8 B(t) - w''(t),
    z(t) = 2 w(t) - s(t) + log(1-exp(-W(t))).

The dimensionless evaluator consumes one common H^(2)..H^(18) collar and
returns an enclosure for t^2 z''(t).  Its analytic W floor uses only the
already-certified order-nine inequality w''(t) <= 4200/t^2.
"""

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


CURVATURE_CONSTANT = 5500
ORDER_NINE_CURVATURE_CONSTANT = 4200
MAX_JET_ORDER = 16


def _intersect_with_floor(value: flint.arb, floor: flint.arb, *, name: str) -> flint.arb:
    exact_floor = flint.arb(floor.lower())
    if not bool(value.upper() > exact_floor):
        raise ValueError(f"analytic {name} floor is disjoint from the raw enclosure")
    narrowed = value.intersection(exact_floor.union(value.upper()))
    if not bool(narrowed > 0):
        raise ValueError(f"analytic {name} floor did not produce a positive enclosure")
    return narrowed


def nested_seventh_curvature_from_h_derivatives(
    derivatives: dict[int, flint.arb],
    *,
    v_floor: flint.arb | None = None,
    w_floor: flint.arb | None = None,
) -> dict:
    """Enclose the seventh stable layer from one H^(2)..H^(18) cover."""
    missing = [order for order in range(2, 19) if order not in derivatives]
    if missing:
        raise ValueError(f"missing H derivatives: {missing}")

    B = [
        derivatives[order + 2] / math.factorial(order)
        for order in range(MAX_JET_ORDER + 1)
    ]
    one = [flint.arb(1)] + [flint.arb(0) for _ in range(MAX_JET_ORDER)]
    ell = series_log(series_sub(one, series_exp(series_scale(B, -1), 16)), 16)

    J = [
        2 * B[order]
        - ell[order + 2] * math.factorial(order + 2) / math.factorial(order)
        for order in range(15)
    ]
    h = series_add(series_scale(ell[:15], 2), stable_log_series(J, 14))

    R = [
        3 * B[order]
        - h[order + 2] * math.factorial(order + 2) / math.factorial(order)
        for order in range(13)
    ]
    q = series_add(
        series_sub(series_scale(h[:13], 2), ell[:13]),
        stable_log_series(R, 12),
    )

    S = [
        4 * B[order]
        - q[order + 2] * math.factorial(order + 2) / math.factorial(order)
        for order in range(11)
    ]
    p = series_add(
        series_sub(series_scale(q[:11], 2), h[:11]),
        stable_log_series(S, 10),
    )

    T = [
        5 * B[order]
        - p[order + 2] * math.factorial(order + 2) / math.factorial(order)
        for order in range(9)
    ]
    r = series_add(
        series_sub(series_scale(p[:9], 2), q[:9]),
        stable_log_series(T, 8),
    )

    U = [
        6 * B[order]
        - r[order + 2] * math.factorial(order + 2) / math.factorial(order)
        for order in range(7)
    ]
    s = series_add(
        series_sub(series_scale(r[:7], 2), p[:7]),
        stable_log_series(U, 6),
    )

    V = [
        7 * B[order]
        - s[order + 2] * math.factorial(order + 2) / math.factorial(order)
        for order in range(5)
    ]
    if v_floor is not None:
        V[0] = _intersect_with_floor(V[0], v_floor, name="V")
    w = series_add(
        series_sub(series_scale(s[:5], 2), r[:5]),
        stable_log_series(V, 4),
    )

    W = [
        8 * B[order]
        - w[order + 2] * math.factorial(order + 2) / math.factorial(order)
        for order in range(3)
    ]
    raw_w_coordinate = W[0]
    if w_floor is not None:
        W[0] = _intersect_with_floor(W[0], w_floor, name="W")
    z_coordinate = series_add(
        series_sub(series_scale(w[:3], 2), s[:3]),
        stable_log_series(W, 2),
    )
    return {
        "B": B,
        "ell": ell,
        "J": J,
        "h": h,
        "R": R,
        "q": q,
        "S": S,
        "p": p,
        "T": T,
        "r": r,
        "U": U,
        "s": s,
        "V": V,
        "w": w,
        "W": W,
        "raw_W": raw_w_coordinate,
        "z": z_coordinate,
        "z_second": 2 * z_coordinate[2],
    }


def dimensionless_seventh_curvature_from_normalized_b(
    t: flint.arb,
    b: list[flint.arb],
    *,
    v_floor: flint.arb | None = None,
    w_floor: flint.arb | None = None,
) -> dict:
    """Enclose t^2 z'' from a common normalized B-jet collar."""
    if len(b) != 17:
        raise ValueError("dimensionless seventh core requires B orders 0..16")
    if not bool(t > 1):
        raise ValueError("dimensionless seventh core requires t>1")

    inverse_t = 1 / t
    ell = stable_log_series(series_scale(b, inverse_t), 16)

    J = [
        2 * b[degree]
        - inverse_t * (degree + 1) * (degree + 2) * ell[degree + 2]
        for degree in range(15)
    ]
    h = series_add(
        series_scale(ell[:15], 2),
        stable_log_series(series_scale(J, inverse_t), 14),
    )

    R = [
        3 * b[degree]
        - inverse_t * (degree + 1) * (degree + 2) * h[degree + 2]
        for degree in range(13)
    ]
    q = series_add(
        series_sub(series_scale(h[:13], 2), ell[:13]),
        stable_log_series(series_scale(R, inverse_t), 12),
    )

    S = [
        4 * b[degree]
        - inverse_t * (degree + 1) * (degree + 2) * q[degree + 2]
        for degree in range(11)
    ]
    p = series_add(
        series_sub(series_scale(q[:11], 2), h[:11]),
        stable_log_series(series_scale(S, inverse_t), 10),
    )

    T = [
        5 * b[degree]
        - inverse_t * (degree + 1) * (degree + 2) * p[degree + 2]
        for degree in range(9)
    ]
    r = series_add(
        series_sub(series_scale(p[:9], 2), q[:9]),
        stable_log_series(series_scale(T, inverse_t), 8),
    )

    U = [
        6 * b[degree]
        - inverse_t * (degree + 1) * (degree + 2) * r[degree + 2]
        for degree in range(7)
    ]
    s = series_add(
        series_sub(series_scale(r[:7], 2), p[:7]),
        stable_log_series(series_scale(U, inverse_t), 6),
    )

    V = [
        7 * b[degree]
        - inverse_t * (degree + 1) * (degree + 2) * s[degree + 2]
        for degree in range(5)
    ]
    if v_floor is not None:
        V[0] = _intersect_with_floor(V[0], v_floor, name="V")
    w = series_add(
        series_sub(series_scale(s[:5], 2), r[:5]),
        stable_log_series(series_scale(V, inverse_t), 4),
    )

    W = [
        8 * b[degree]
        - inverse_t * (degree + 1) * (degree + 2) * w[degree + 2]
        for degree in range(3)
    ]
    raw_w_coordinate = W[0]
    if w_floor is not None:
        W[0] = _intersect_with_floor(W[0], w_floor, name="W")
    z_coordinate = series_add(
        series_sub(series_scale(w[:3], 2), s[:3]),
        stable_log_series(series_scale(W, inverse_t), 2),
    )
    return {
        "b": b,
        "ell": ell,
        "J": J,
        "h": h,
        "R": R,
        "q": q,
        "S": S,
        "p": p,
        "T": T,
        "r": r,
        "U": U,
        "s": s,
        "V": V,
        "w": w,
        "W": W,
        "raw_W": raw_w_coordinate,
        "z": z_coordinate,
        "scaled_z_second": 2 * z_coordinate[2],
    }


def dimensionless_seventh_curvature_from_h_derivatives(
    t: flint.arb,
    derivatives: dict[int, flint.arb],
    *,
    v_floor: flint.arb | None = None,
    w_floor: flint.arb | None = None,
) -> dict:
    """Enclose t^2 z'' from a dimensionless common H-derivative collar."""
    missing = [order for order in range(2, 19) if order not in derivatives]
    if missing:
        raise ValueError(f"missing H derivatives: {missing}")
    b = [
        t ** (degree + 1)
        * derivatives[degree + 2]
        / math.factorial(degree)
        for degree in range(17)
    ]
    return dimensionless_seventh_curvature_from_normalized_b(
        t,
        b,
        v_floor=v_floor,
        w_floor=w_floor,
    )


def analytic_unscaled_w_floor(t: flint.arb) -> flint.arb:
    """Return a uniform lower bound for W(t) on a positive t interval."""
    lower_t = flint.arb(t.lower())
    upper_t = flint.arb(t.upper())
    return (
        8 / (2 * upper_t + 3)
        - ORDER_NINE_CURVATURE_CONSTANT / (lower_t**2 - 1)
    ).lower()


def analytic_dimensionless_w_floor(t: flint.arb) -> flint.arb:
    """Return a uniform lower bound for t*W(t) on a positive t interval."""
    return (flint.arb(t.lower()) * analytic_unscaled_w_floor(t)).lower()


def dimensionless_formula_upper(
    t: flint.arb,
    jets: dict,
    *,
    unscaled_w_floor: flint.arb,
) -> dict:
    """Apply the cancellation-preserving one-sided formula for t^2 z''."""
    w_second_scaled = 2 * jets["w"][2]
    s_second_scaled = 2 * jets["s"][2]
    w_coordinate_second_scaled = 2 * jets["W"][2] / t
    positive_w_coordinate_second = flint.arb(
        max(flint.arb(0), flint.arb(w_coordinate_second_scaled.upper()))
    )
    phi = 1 / flint.arb(unscaled_w_floor.lower()).expm1()
    base_upper = flint.arb((2 * w_second_scaled - s_second_scaled).upper())
    positive_term_upper = flint.arb(phi.upper()) * positive_w_coordinate_second
    return {
        "w_second_scaled": w_second_scaled,
        "s_second_scaled": s_second_scaled,
        "W_second_scaled": w_coordinate_second_scaled,
        "phi_upper": flint.arb(phi.upper()),
        "positive_term_upper": positive_term_upper,
        "formula_upper": base_upper + positive_term_upper,
    }


def _evaluate_dimensionless_jets(
    left,
    right,
    t: flint.arb,
    jets: dict,
    *,
    v_floor: flint.arb,
    unscaled_w_floor: flint.arb,
    w_floor: flint.arb,
    diagnostics: dict | None,
) -> dict:
    coordinates = {
        "J": jets["J"][0],
        "R": jets["R"][0],
        "S": jets["S"][0],
        "T": jets["T"][0],
        "U": jets["U"][0],
        "V": jets["V"][0],
        "W": jets["W"][0],
    }
    scaled = jets["scaled_z_second"]
    formula = dimensionless_formula_upper(
        t,
        jets,
        unscaled_w_floor=unscaled_w_floor,
    )
    formula_upper = formula["formula_upper"]
    margin = flint.arb(CURVATURE_CONSTANT) - formula_upper
    if not bool(
        jets["b"][0] > 0
        and all(value > 0 for value in coordinates.values())
        and margin > 0
    ):
        return {
            "passed": False,
            "failure": "dimensionless-margin",
            "mode": [str(left), str(right)],
            "t_ball": t.str(40).replace("e", "E"),
            **{
                f"{name}_ball": value.str(40).replace("e", "E")
                for name, value in coordinates.items()
            },
            "raw_W_ball": jets["raw_W"].str(40).replace("e", "E"),
            "scaled_composed_ball": scaled.str(40).replace("e", "E"),
            "scaled_formula_upper": arb_upper_text(formula_upper),
        }
    return {
        "passed": True,
        "failure": None,
        "mode": [str(left), str(right)],
        "t_ball": t.str(40).replace("e", "E"),
        "analytic_V_floor": arb_lower_text(v_floor),
        "analytic_unscaled_W_floor": arb_lower_text(unscaled_w_floor),
        "analytic_dimensionless_W_floor": arb_lower_text(w_floor),
        **{
            f"{name}_lower": arb_lower_text(value)
            for name, value in coordinates.items()
        },
        "raw_W_ball": jets["raw_W"].str(40).replace("e", "E"),
        "scaled_composed_ball": scaled.str(40).replace("e", "E"),
        "scaled_curvature_upper": arb_upper_text(formula_upper),
        "scaled_margin_lower": arb_lower_text(margin),
        "stable_w_second_scaled_upper": arb_upper_text(
            formula["w_second_scaled"]
        ),
        "s_second_scaled_lower": arb_lower_text(formula["s_second_scaled"]),
        "coordinate_W_second_scaled_upper": arb_upper_text(
            formula["W_second_scaled"]
        ),
        "phi_upper": arb_upper_text(formula["phi_upper"]),
        "positive_phi_W_second_upper": arb_upper_text(
            formula["positive_term_upper"]
        ),
        "diagnostics": diagnostics,
    }


def evaluate_dimensionless_seventh_curvature(
    left,
    right,
    derivatives: dict[int, flint.arb],
    *,
    diagnostics: dict | None = None,
) -> dict:
    """Certify z_1''(t)<=5500/t^2 on one central saddle-mode interval."""
    mode = (
        (arb_rational(left) + arb_rational(right)) / 2
        + flint.arb(0, (arb_rational(right) - arb_rational(left)) / 2)
    )
    t = potential_jet_arb(mode, 1)[1]
    if not bool(t > 1):
        return {"passed": False, "failure": "nonpositive-t"}
    v_floor = flint.arb(4) / 3
    unscaled_w_floor = analytic_unscaled_w_floor(t)
    w_floor = (flint.arb(t.lower()) * unscaled_w_floor).lower()
    try:
        jets = dimensionless_seventh_curvature_from_h_derivatives(
            t,
            derivatives,
            v_floor=v_floor,
            w_floor=w_floor,
        )
    except Exception as exc:
        return {
            "passed": False,
            "failure": "nested-jet-exception",
            "detail": repr(exc),
        }

    return _evaluate_dimensionless_jets(
        left,
        right,
        t,
        jets,
        v_floor=v_floor,
        unscaled_w_floor=unscaled_w_floor,
        w_floor=w_floor,
        diagnostics=diagnostics,
    )


def evaluate_dimensionless_seventh_curvature_from_normalized_b(
    left,
    right,
    b: list[flint.arb],
    *,
    diagnostics: dict | None = None,
) -> dict:
    """Certify the seventh curvature from a pre-normalized common B cover."""
    mode = (
        (arb_rational(left) + arb_rational(right)) / 2
        + flint.arb(0, (arb_rational(right) - arb_rational(left)) / 2)
    )
    t = potential_jet_arb(mode, 1)[1]
    if not bool(t > 1):
        return {"passed": False, "failure": "nonpositive-t"}
    v_floor = flint.arb(4) / 3
    unscaled_w_floor = analytic_unscaled_w_floor(t)
    w_floor = (flint.arb(t.lower()) * unscaled_w_floor).lower()
    try:
        jets = dimensionless_seventh_curvature_from_normalized_b(
            t,
            b,
            v_floor=v_floor,
            w_floor=w_floor,
        )
    except Exception as exc:
        return {
            "passed": False,
            "failure": "nested-jet-exception",
            "detail": repr(exc),
        }
    return _evaluate_dimensionless_jets(
        left,
        right,
        t,
        jets,
        v_floor=v_floor,
        unscaled_w_floor=unscaled_w_floor,
        w_floor=w_floor,
        diagnostics=diagnostics,
    )
