#!/usr/bin/env python3
"""Dimensionless interval core for order-eleven first-summand curvature.

The eighth stable layer is

    X(t) = 9 B(t) - z''(t),
    y(t) = 2 z(t) - w(t) + log(1-exp(-X(t))).

The evaluator consumes a common H^(2)..H^(20) collar and encloses
``t^2*y''(t)``. Analytic floors use only the already-certified global
order-nine and order-ten first-summand curvature theorems.
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


CURVATURE_CONSTANT = 6000
ORDER_NINE_CURVATURE_CONSTANT = 4200
ORDER_TEN_CURVATURE_CONSTANT = 4200
MAX_JET_ORDER = 18


def _intersect_with_floor(value: flint.arb, floor: flint.arb, *, name: str) -> flint.arb:
    exact_floor = flint.arb(floor.lower())
    if not bool(value.upper() > exact_floor):
        raise ValueError(f"analytic {name} floor is disjoint from the raw enclosure")
    narrowed = value.intersection(exact_floor.union(value.upper()))
    if not bool(narrowed > 0):
        raise ValueError(f"analytic {name} floor did not produce a positive enclosure")
    return narrowed


def dimensionless_eighth_curvature_from_normalized_b(
    t: flint.arb,
    b: list[flint.arb],
    *,
    v_floor: flint.arb | None = None,
    w_floor: flint.arb | None = None,
    x_floor: flint.arb | None = None,
) -> dict:
    """Enclose ``t^2*y''`` from one normalized B-jet collar."""
    if len(b) != MAX_JET_ORDER + 1:
        raise ValueError("dimensionless eighth core requires B orders 0..18")
    if not bool(t > 1):
        raise ValueError("dimensionless eighth core requires t>1")

    inverse_t = 1 / t
    ell = stable_log_series(series_scale(b, inverse_t), 18)
    J = [
        2 * b[d] - inverse_t * (d + 1) * (d + 2) * ell[d + 2]
        for d in range(17)
    ]
    h = series_add(
        series_scale(ell[:17], 2),
        stable_log_series(series_scale(J, inverse_t), 16),
    )
    R = [
        3 * b[d] - inverse_t * (d + 1) * (d + 2) * h[d + 2]
        for d in range(15)
    ]
    q = series_add(
        series_sub(series_scale(h[:15], 2), ell[:15]),
        stable_log_series(series_scale(R, inverse_t), 14),
    )
    S = [
        4 * b[d] - inverse_t * (d + 1) * (d + 2) * q[d + 2]
        for d in range(13)
    ]
    p = series_add(
        series_sub(series_scale(q[:13], 2), h[:13]),
        stable_log_series(series_scale(S, inverse_t), 12),
    )
    T = [
        5 * b[d] - inverse_t * (d + 1) * (d + 2) * p[d + 2]
        for d in range(11)
    ]
    r = series_add(
        series_sub(series_scale(p[:11], 2), q[:11]),
        stable_log_series(series_scale(T, inverse_t), 10),
    )
    U = [
        6 * b[d] - inverse_t * (d + 1) * (d + 2) * r[d + 2]
        for d in range(9)
    ]
    s = series_add(
        series_sub(series_scale(r[:9], 2), p[:9]),
        stable_log_series(series_scale(U, inverse_t), 8),
    )
    V = [
        7 * b[d] - inverse_t * (d + 1) * (d + 2) * s[d + 2]
        for d in range(7)
    ]
    if v_floor is not None:
        V[0] = _intersect_with_floor(V[0], v_floor, name="V")
    w = series_add(
        series_sub(series_scale(s[:7], 2), r[:7]),
        stable_log_series(series_scale(V, inverse_t), 6),
    )
    W = [
        8 * b[d] - inverse_t * (d + 1) * (d + 2) * w[d + 2]
        for d in range(5)
    ]
    if w_floor is not None:
        W[0] = _intersect_with_floor(W[0], w_floor, name="W")
    z = series_add(
        series_sub(series_scale(w[:5], 2), s[:5]),
        stable_log_series(series_scale(W, inverse_t), 4),
    )
    X = [
        9 * b[d] - inverse_t * (d + 1) * (d + 2) * z[d + 2]
        for d in range(3)
    ]
    raw_x = X[0]
    if x_floor is not None:
        X[0] = _intersect_with_floor(X[0], x_floor, name="X")
    y = series_add(
        series_sub(series_scale(z[:3], 2), w[:3]),
        stable_log_series(series_scale(X, inverse_t), 2),
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
        "z": z,
        "X": X,
        "raw_X": raw_x,
        "y": y,
        "scaled_y_second": 2 * y[2],
    }


def dimensionless_eighth_curvature_from_h_derivatives(
    t: flint.arb,
    derivatives: dict[int, flint.arb],
    **floors,
) -> dict:
    """Convert a common H^(2)..H^(20) collar to the normalized core."""
    missing = [order for order in range(2, 21) if order not in derivatives]
    if missing:
        raise ValueError(f"missing H derivatives: {missing}")
    b = [
        t ** (degree + 1) * derivatives[degree + 2] / math.factorial(degree)
        for degree in range(MAX_JET_ORDER + 1)
    ]
    return dimensionless_eighth_curvature_from_normalized_b(t, b, **floors)


def _analytic_unscaled_floor(
    t: flint.arb,
    multiplier: int,
    inherited_constant: int,
) -> flint.arb:
    lower_t = flint.arb(t.lower())
    upper_t = flint.arb(t.upper())
    return (
        multiplier / (2 * upper_t + 3)
        - inherited_constant / (lower_t**2 - 1)
    ).lower()


def analytic_unscaled_w_floor(t: flint.arb) -> flint.arb:
    return _analytic_unscaled_floor(t, 8, ORDER_NINE_CURVATURE_CONSTANT)


def analytic_unscaled_x_floor(t: flint.arb) -> flint.arb:
    return _analytic_unscaled_floor(t, 9, ORDER_TEN_CURVATURE_CONSTANT)


def dimensionless_formula_upper(
    t: flint.arb,
    jets: dict,
    *,
    unscaled_x_floor: flint.arb,
) -> dict:
    """Apply the cancellation-preserving one-sided formula for ``t^2*y''``."""
    z_second_scaled = 2 * jets["z"][2]
    w_second_scaled = 2 * jets["w"][2]
    x_second_scaled = 2 * jets["X"][2] / t
    positive_x_second = flint.arb(
        max(flint.arb(0), flint.arb(x_second_scaled.upper()))
    )
    phi = 1 / flint.arb(unscaled_x_floor.lower()).expm1()
    base_upper = flint.arb((2 * z_second_scaled - w_second_scaled).upper())
    positive_term = flint.arb(phi.upper()) * positive_x_second
    return {
        "z_second_scaled": z_second_scaled,
        "w_second_scaled": w_second_scaled,
        "X_second_scaled": x_second_scaled,
        "phi_upper": flint.arb(phi.upper()),
        "positive_term_upper": positive_term,
        "formula_upper": base_upper + positive_term,
    }


def _evaluate(
    left,
    right,
    t: flint.arb,
    jets: dict,
    *,
    v_floor: flint.arb,
    unscaled_w_floor: flint.arb,
    w_floor: flint.arb,
    unscaled_x_floor: flint.arb,
    x_floor: flint.arb,
    diagnostics: dict | None,
) -> dict:
    coordinates = {name: jets[name][0] for name in "J R S T U V W X".split()}
    formula = dimensionless_formula_upper(
        t,
        jets,
        unscaled_x_floor=unscaled_x_floor,
    )
    margin = flint.arb(CURVATURE_CONSTANT) - formula["formula_upper"]
    if not bool(jets["b"][0] > 0 and all(v > 0 for v in coordinates.values()) and margin > 0):
        return {
            "passed": False,
            "failure": "dimensionless-margin",
            "mode": [str(left), str(right)],
            "t_ball": t.str(40).replace("e", "E"),
            "scaled_formula_upper": arb_upper_text(formula["formula_upper"]),
        }
    return {
        "passed": True,
        "failure": None,
        "mode": [str(left), str(right)],
        "t_ball": t.str(40).replace("e", "E"),
        "analytic_V_floor": arb_lower_text(v_floor),
        "analytic_unscaled_W_floor": arb_lower_text(unscaled_w_floor),
        "analytic_dimensionless_W_floor": arb_lower_text(w_floor),
        "analytic_unscaled_X_floor": arb_lower_text(unscaled_x_floor),
        "analytic_dimensionless_X_floor": arb_lower_text(x_floor),
        **{f"{name}_lower": arb_lower_text(value) for name, value in coordinates.items()},
        "raw_X_ball": jets["raw_X"].str(40).replace("e", "E"),
        "scaled_composed_ball": jets["scaled_y_second"].str(40).replace("e", "E"),
        "scaled_curvature_upper": arb_upper_text(formula["formula_upper"]),
        "scaled_margin_lower": arb_lower_text(margin),
        "stable_z_second_scaled_upper": arb_upper_text(formula["z_second_scaled"]),
        "w_second_scaled_lower": arb_lower_text(formula["w_second_scaled"]),
        "coordinate_X_second_scaled_upper": arb_upper_text(formula["X_second_scaled"]),
        "phi_upper": arb_upper_text(formula["phi_upper"]),
        "positive_phi_X_second_upper": arb_upper_text(formula["positive_term_upper"]),
        "diagnostics": diagnostics,
    }


def _floors(t: flint.arb) -> tuple[flint.arb, ...]:
    v_floor = flint.arb(4) / 3
    unscaled_w = analytic_unscaled_w_floor(t)
    unscaled_x = analytic_unscaled_x_floor(t)
    lower_t = flint.arb(t.lower())
    return v_floor, unscaled_w, (lower_t * unscaled_w).lower(), unscaled_x, (lower_t * unscaled_x).lower()


def evaluate_dimensionless_eighth_curvature(
    left,
    right,
    derivatives: dict[int, flint.arb],
    *,
    diagnostics: dict | None = None,
) -> dict:
    mode = (
        (arb_rational(left) + arb_rational(right)) / 2
        + flint.arb(0, (arb_rational(right) - arb_rational(left)) / 2)
    )
    t = potential_jet_arb(mode, 1)[1]
    if not bool(t > 1):
        return {"passed": False, "failure": "nonpositive-t"}
    v_floor, unscaled_w, w_floor, unscaled_x, x_floor = _floors(t)
    try:
        jets = dimensionless_eighth_curvature_from_h_derivatives(
            t,
            derivatives,
            v_floor=v_floor,
            w_floor=w_floor,
            x_floor=x_floor,
        )
    except Exception as exc:
        return {"passed": False, "failure": "nested-jet-exception", "detail": repr(exc)}
    return _evaluate(
        left,
        right,
        t,
        jets,
        v_floor=v_floor,
        unscaled_w_floor=unscaled_w,
        w_floor=w_floor,
        unscaled_x_floor=unscaled_x,
        x_floor=x_floor,
        diagnostics=diagnostics,
    )


def evaluate_dimensionless_eighth_curvature_from_normalized_b(
    left,
    right,
    b: list[flint.arb],
    *,
    diagnostics: dict | None = None,
) -> dict:
    mode = (
        (arb_rational(left) + arb_rational(right)) / 2
        + flint.arb(0, (arb_rational(right) - arb_rational(left)) / 2)
    )
    t = potential_jet_arb(mode, 1)[1]
    if not bool(t > 1):
        return {"passed": False, "failure": "nonpositive-t"}
    v_floor, unscaled_w, w_floor, unscaled_x, x_floor = _floors(t)
    try:
        jets = dimensionless_eighth_curvature_from_normalized_b(
            t,
            b,
            v_floor=v_floor,
            w_floor=w_floor,
            x_floor=x_floor,
        )
    except Exception as exc:
        return {"passed": False, "failure": "nested-jet-exception", "detail": repr(exc)}
    return _evaluate(
        left,
        right,
        t,
        jets,
        v_floor=v_floor,
        unscaled_w_floor=unscaled_w,
        w_floor=w_floor,
        unscaled_x_floor=unscaled_x,
        x_floor=x_floor,
        diagnostics=diagnostics,
    )
