#!/usr/bin/env python3
"""Scout the cancellation-preserving height derivative of the transition join.

This diagnostic differentiates the finite transition source before applying
the alternating boundary-jet collapse, then subtracts the exact derivative of
the matching 42-mode Gamma block.  It does not yet bound the complete K_T.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = ROOT / "work" / "rh_compute" / "vendor"
for directory in (SCRIPT_DIR, VENDOR):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_finite_Fresnel_transition_alternating_boundary_jet_packet_gate as transition
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined


HEIGHT = 10_000_000_000
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_transition_height_derivative_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
EVENT_ATLAS = (
    ROOT
    / "work"
    / "rh_compute"
    / "results"
    / "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_local_height_event_atlas_gate.json"
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(55, more=True),
        "imag_ball": value.imag.str(55, more=True),
        "absolute_ball": abs(value).str(55, more=True),
        "real_radius": value.real.rad().str(25, more=True),
        "imag_radius": value.imag.rad().str(25, more=True),
    }


def real_record(value: arb) -> dict[str, str]:
    return {
        "ball": value.str(55, more=True),
        "lower": value.lower().str(40, more=True),
        "upper": value.upper().str(40, more=True),
        "radius": value.rad().str(25, more=True),
    }


def complex_from_record(record: dict[str, Any]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def time_factor_jets(y: arb, order: int, zeroth_shift: acb | None = None) -> list[acb]:
    """Return y-jets of -i log y plus an optional y-constant shift."""

    jets = [-acb(0, 1) * acb(y).log()]
    if zeroth_shift is not None:
        jets[0] += zeroth_shift
    for k in range(1, order + 1):
        real_log_jet = ((-1) ** (k - 1)) * math.factorial(k - 1) / y**k
        jets.append(-acb(0, 1) * real_log_jet)
    return jets


def time_endpoint_derivatives(
    y: arb,
    p: dict[str, Any],
    order: int,
    series_terms: int,
    zeroth_shift: acb | None = None,
) -> list[acb]:
    """Return y-jets of (-i log y + shift)F, 0 <= r < order."""

    F_derivatives = transition.endpoint_derivatives(y, p, order, series_terms)
    factor_jets = time_factor_jets(y, order - 1, zeroth_shift)
    values: list[acb] = []
    for r in range(order):
        value = acb(0)
        for k in range(r + 1):
            value += math.comb(r, k) * factor_jets[k] * F_derivatives[r - k]
        values.append(value)
    return values


def normalized_time_jet(
    y: arb,
    p: dict[str, Any],
    order: int,
    zeroth_shift: acb | None = None,
) -> acb:
    """Return d_y^order[(-i log y + shift)F]/F."""

    normalized_F = transition.bell_jet(transition.log_derivatives(y, p, order))
    factor_jets = time_factor_jets(y, order, zeroth_shift)
    value = acb(0)
    for k in range(order + 1):
        value += math.comb(order, k) * factor_jets[k] * normalized_F[order - k]
    return value


def first_label_time_derivative_integral(
    p: dict[str, Any],
    left: arb,
    right: arb,
    panels: int,
    series_terms: int,
    tolerance: str,
    zeroth_shift: acb | None = None,
) -> tuple[acb, list[dict[str, str]]]:
    total = acb(0)
    rows: list[dict[str, str]] = []
    width = right - left

    def integrand(y: acb, analytic: bool) -> acb:
        F = transition.centered_first_label_value(y, p, series_terms, analytic)
        factor = -acb(0, 1) * y.log()
        if zeroth_shift is not None:
            factor += zeroth_shift
        return factor * F

    for index in range(panels):
        panel_left = left + width * index / panels
        panel_right = left + width * (index + 1) / panels
        value = acb.integral(
            integrand,
            panel_left,
            panel_right,
            abs_tol=arb(tolerance),
            rel_tol=arb(tolerance),
            eval_limit=300_000,
            depth_limit=40,
        )
        total += value
        rows.append(
            {
                "left": panel_left.str(35, more=True),
                "right": panel_right.str(35, more=True),
                "absolute_ball": abs(value).str(35, more=True),
            }
        )
    return total, rows


def time_derivative_remainder_bound(
    p: dict[str, Any],
    left: arb,
    right: arb,
    order: int,
    slabs: int,
    zeroth_shift: acb | None = None,
) -> tuple[arb, arb]:
    width = right - left
    integral_bound = arb(0)
    maximum_normalized_jet = arb(0)
    for index in range(slabs):
        slab_left = left + width * index / slabs
        slab_right = left + width * (index + 1) / slabs
        y = transition.real_interval(slab_left, slab_right)
        normalized_jet = normalized_time_jet(y, p, order, zeroth_shift)
        jet_bound = abs(normalized_jet).upper()
        F_bound = slab_left ** arb("-0.5")
        integral_bound += (slab_right - slab_left) * F_bound * jet_bound
        maximum_normalized_jet = max(maximum_normalized_jet, jet_bound)
    return integral_bound, maximum_normalized_jet


def transition_source_time_derivative_certificate(
    t: arb,
    endpoint: int,
    source_count: int,
    left: arb,
    right: arb,
    order: int,
    derivative_slabs: int,
    integral_panels: int,
    series_terms: int,
    tolerance: str,
    zeroth_shift: acb | None = None,
    operator_identity: str = "d_t F(y,t)=-i*log(y)*F(y,t)",
) -> dict[str, Any]:
    require(source_count >= 3 and source_count % 2 == 1, "source count must be odd")
    later_count = source_count - 1
    require(later_count % 2 == 0, "later-label count parity drift")
    p = transition.centered_parameters(t, endpoint)
    first, panel_rows = first_label_time_derivative_integral(
        p,
        left,
        right,
        integral_panels,
        series_terms,
        tolerance,
        zeroth_shift,
    )
    left_derivatives = time_endpoint_derivatives(left, p, order, series_terms, zeroth_shift)
    right_derivatives = time_endpoint_derivatives(right, p, order, series_terms, zeroth_shift)
    imaginary = acb(0, 1)
    boundary = acb(0)
    alternating_rows: list[dict[str, Any]] = []
    for r in range(order):
        power = r + 1
        alternating = transition.alternating_power_sum(power, later_count)
        derivative_jump = right_derivatives[r] - left_derivatives[r]
        term = ((-1) ** r) * derivative_jump * alternating / (imaginary * 2 * p["pi"]) ** power
        boundary += term
        alternating_rows.append(
            {
                "power": power,
                "alternating_sum_ball": alternating.str(55, more=True),
                "boundary_term": complex_record(term),
            }
        )
    derivative_bound, maximum_jet = time_derivative_remainder_bound(
        p,
        left,
        right,
        order,
        derivative_slabs,
        zeroth_shift,
    )
    harmonic_bound = arb(order).zeta()
    remainder = derivative_bound * harmonic_bound / (2 * p["pi"]) ** order
    complete = transition.add_complex_error(first + boundary, remainder)
    return {
        "height_ball": t.str(55, more=True),
        "source_endpoint": endpoint,
        "source_count": source_count,
        "later_label_count": later_count,
        "cell_interval": [left.str(35, more=True), right.str(35, more=True)],
        "jet_order": order,
        "derivative_slabs": derivative_slabs,
        "integral_panels": integral_panels,
        "logarithm_series_terms": series_terms,
        "exact_operator_weight": operator_identity,
        "first_label_time_derivative_integral_ball": complex_record(first),
        "boundary_time_derivative_jet_sum_ball": complex_record(boundary),
        "time_derivative_integral_bound_ball": derivative_bound.str(55, more=True),
        "maximum_normalized_time_derivative_order_jet_ball": maximum_jet.str(55, more=True),
        "later_label_time_derivative_remainder_bound_ball": remainder.str(55, more=True),
        "complete_transition_source_time_derivative_ball": complex_record(complete),
        "complete_weighted_transition_source_ball": complex_record(complete),
        "alternating_rows": alternating_rows,
        "panel_rows": panel_rows,
    }


def altered_five_label_fixture() -> dict[str, Any]:
    """Compare the derivative collapse with a direct finite source integral."""

    pi = arb.pi()
    t = 2 * pi * arb(16)
    endpoint = 16
    source_count = 5
    left = arb("3.5")
    right = arb("4.5")
    order = 8
    panels = 24
    series_terms = 24
    p = transition.centered_parameters(t, endpoint)
    direct = acb(0)
    imaginary = acb(0, 1)
    width = right - left

    def integrand(y: acb, analytic: bool) -> acb:
        F = transition.centered_first_label_value(y, p, series_terms, analytic)
        finite_roster = acb(0)
        for label in range(source_count):
            finite_roster += (imaginary * 2 * pi * label * y).exp()
        return -imaginary * y.log() * F * finite_roster

    for index in range(panels):
        panel_left = left + width * index / panels
        panel_right = left + width * (index + 1) / panels
        direct += acb.integral(
            integrand,
            panel_left,
            panel_right,
            abs_tol=arb("1e-24"),
            rel_tol=arb("1e-24"),
            eval_limit=100_000,
            depth_limit=35,
        )
    collapsed_certificate = transition_source_time_derivative_certificate(
        t,
        endpoint,
        source_count,
        left,
        right,
        order,
        128,
        panels,
        series_terms,
        "1e-24",
    )
    collapsed = complex_from_record(
        collapsed_certificate["complete_transition_source_time_derivative_ball"]
    )
    require(direct.overlaps(collapsed), "altered derivative source and boundary-jet collapse miss")
    return {
        "height_ball": t.str(55, more=True),
        "source_endpoint": endpoint,
        "source_count": source_count,
        "cell_interval": [left.str(20, more=True), right.str(20, more=True)],
        "direct_finite_source_time_derivative_ball": complex_record(direct),
        "boundary_jet_time_derivative_ball": complex_record(collapsed),
        "overlap": True,
        "purpose": "independent finite-roster check of the differentiated collapse and its sign",
    }


def weighted_finite_dirichlet(
    t: arb,
    start: int,
    end: int,
    theta_prime: arb,
    real_shift: arb,
) -> acb:
    """Sum [i(theta'-log m)+real_shift]m^(-1/2-it) per mode."""

    s = acb(arb("0.5"), t)
    total = acb(0)
    for mode in range(start, end + 1):
        log_mode = arb(mode).log()
        term = (-s * log_mode).exp()
        weight = acb(real_shift, theta_prime - log_mode)
        total += weight * term
    return total


def build(radius_text: str, precision_bits: int, variant: str = "production") -> dict[str, Any]:
    started = time.perf_counter()
    resource_mode = joined.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = precision_bits
    ctx.threads = 1
    radius = arb(radius_text)
    require(radius >= 0 and radius <= arb("0.001"), "radius must lie in [0,0.001]")
    t = arb(arb(HEIGHT), radius)

    event = load_json(EVENT_ATLAS)["certificate"]["primary_open_cell"]
    event_lower = arb(event["lower_ball"]["ball"])
    event_upper = arb(event["upper_ball"]["ball"])
    require(event_lower.upper() < t.lower() < t.upper() < event_upper.lower(), "box left the fixed-roster cell")

    require(variant in ("production", "independent"), "unknown derivative variant")
    if variant == "production":
        order = 14
        derivative_slabs = 256
        integral_panels = 48
        series_terms = 34
        tolerance = "1e-27"
    else:
        order = 16
        derivative_slabs = 384
        integral_panels = 64
        series_terms = 40
        tolerance = "1e-29"
    theta, theta_prime = joined.theta_and_derivative(t)
    H = joined.direct_log_H(t).exp()
    H_log_prime = joined.direct_log_H_derivative(t)
    require(H.lower() > 0, "H lost positivity")
    Hardy_shift = acb(0, theta_prime) - H_log_prime
    source = transition.transition_source_certificate(
        t,
        transition.SOURCE_ENDPOINT,
        transition.SOURCE_COUNT,
        transition.CELL_LEFT,
        transition.CELL_RIGHT,
        order,
        derivative_slabs,
        integral_panels,
        series_terms,
        tolerance,
    )
    source_K_certificate = transition_source_time_derivative_certificate(
        t,
        transition.SOURCE_ENDPOINT,
        transition.SOURCE_COUNT,
        transition.CELL_LEFT,
        transition.CELL_RIGHT,
        order,
        derivative_slabs,
        integral_panels,
        series_terms,
        tolerance,
        Hardy_shift,
        "(d_t+i*theta'-H'/H)F=[i*(theta'-log(y))-H'/H]*F",
    )
    source_prime = transition_source_time_derivative_certificate(
        t,
        transition.SOURCE_ENDPOINT,
        transition.SOURCE_COUNT,
        transition.CELL_LEFT,
        transition.CELL_RIGHT,
        order,
        derivative_slabs,
        integral_panels,
        series_terms,
        tolerance,
    )
    source_value = complex_from_record(source["complete_transition_source_integral_ball"])
    source_derivative = complex_from_record(
        source_prime["complete_transition_source_time_derivative_ball"]
    )
    source_K = complex_from_record(
        source_K_certificate["complete_weighted_transition_source_ball"]
    )

    pi = arb.pi()
    decay = (-2 * pi * t).exp()
    C_G = (1 + decay) ** (-arb(1) / 2)
    C_G_log_prime = pi * decay / (1 + decay)
    C_G_prime = C_G * C_G_log_prime
    D, D_prime = joined.finite_dirichlet(
        t,
        joined.LOWER_TRANSITION_START,
        joined.TARGET_END,
    )
    gamma_value = C_G * D
    gamma_derivative = C_G_prime * D + C_G * D_prime
    gamma_K_assembled = C_G * (
        D_prime + (acb(0, theta_prime) + C_G_log_prime - H_log_prime) * D
    )
    gamma_K_direct = C_G * weighted_finite_dirichlet(
        t,
        joined.LOWER_TRANSITION_START,
        joined.TARGET_END,
        theta_prime,
        C_G_log_prime - H_log_prime,
    )
    require(gamma_K_direct.overlaps(gamma_K_assembled), "direct and assembled Gamma K miss")
    transition_join = source_value - gamma_value
    transition_join_derivative = source_derivative - gamma_derivative
    transition_K_assembled = transition_join_derivative + Hardy_shift * transition_join
    transition_K_direct = source_K - gamma_K_direct
    require(
        transition_K_direct.overlaps(transition_K_assembled),
        "direct and assembled transition K enclosures miss",
    )
    transition_Q_derivative_contribution = joined.hardy_projection(theta, transition_K_direct) / H

    return {
        "kind": "rh_diagnostic_joined_transition_height_derivative_scout",
        "date": "2026-08-28",
        "status": "joined_transition_height_derivative_endpoint_jet_diagnostic_not_complete_K_T",
        "passed": True,
        "resource_mode": resource_mode,
        "workers": 1,
        "precision_bits": precision_bits,
        "variant": variant,
        "center_height": str(HEIGHT),
        "radius": radius_text,
        "height_ball": t.str(55, more=True),
        "exact_identities": {
            "source_time_derivative": "d_t F(y,t)=-i*log(y)*F(y,t)",
            "joined_transition": "J_tr=I_tr-C_G*D_(39853..39894)",
            "joined_transition_derivative": "J_tr'=I_tr'-C_G'*D-C_G*D'",
            "joined_Hardy_derivative_coordinate": "K_tr=J_tr'+(i*theta'-H'/H)*J_tr",
            "source_Hardy_operator": "(d_t+i*theta'-H'/H)F=[i*(theta'-log(y))-H'/H]*F",
            "Gamma_Hardy_operator": (
                "(d_t+i*theta'-H'/H)(C_G*D)=C_G[D'+"
                "(i*theta'+C_G'/C_G-H'/H)D]"
            ),
        },
        "configuration": {
            "jet_order": order,
            "derivative_slabs": derivative_slabs,
            "integral_panels": integral_panels,
            "series_terms": series_terms,
            "tolerance": tolerance,
            "Gamma_mode_range": [joined.LOWER_TRANSITION_START, joined.TARGET_END],
        },
        "source_value": complex_record(source_value),
        "source_derivative": complex_record(source_derivative),
        "source_K": complex_record(source_K),
        "Gamma_value": complex_record(gamma_value),
        "Gamma_derivative": complex_record(gamma_derivative),
        "Gamma_K_direct": complex_record(gamma_K_direct),
        "Gamma_K_assembled": complex_record(gamma_K_assembled),
        "direct_and_assembled_Gamma_K_overlap": True,
        "joined_transition": complex_record(transition_join),
        "joined_transition_derivative": complex_record(transition_join_derivative),
        "phase_and_prefactor": {
            "theta": real_record(theta),
            "theta_prime": real_record(theta_prime),
            "H": real_record(H),
            "H_log_prime": real_record(H_log_prime),
            "C_G": real_record(C_G),
            "C_G_log_prime": real_record(C_G_log_prime),
        },
        "transition_K_direct": complex_record(transition_K_direct),
        "transition_K_assembled": complex_record(transition_K_assembled),
        "direct_and_assembled_transition_K_overlap": True,
        "transition_Q_derivative_contribution": real_record(transition_Q_derivative_contribution),
        "altered_five_label_fixture": altered_five_label_fixture(),
        "source_value_certificate": source,
        "source_derivative_certificate": source_prime,
        "source_K_certificate": source_K_certificate,
        "elapsed_seconds": time.perf_counter() - started,
        "diagnostic_boundary": (
            "This is the exact differentiated transition join on one compact height box. It is not "
            "the complete K_T: lower-cell, upper-arc, ordinary-packet, and tiny target-correction "
            "derivatives remain to be joined before any Q_K-T transport theorem can be promoted."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", default="0.0001")
    parser.add_argument("--precision-bits", type=int, default=384)
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    args = parser.parse_args()
    artifact = build(args.radius, args.precision_bits, args.variant)
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    joined_prime = artifact["joined_transition_derivative"]
    contribution = artifact["transition_Q_derivative_contribution"]
    print(
        "diagnostic joined-transition height derivative "
        f"real={joined_prime['real_ball']} imag={joined_prime['imag_ball']} "
        f"Hardy_contribution={contribution['ball']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
