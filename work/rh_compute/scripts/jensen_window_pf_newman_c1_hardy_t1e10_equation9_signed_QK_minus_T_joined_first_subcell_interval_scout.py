#!/usr/bin/env python3
"""Scout a direct Arb height box for the joined Q_K-T packet.

This is deliberately diagnostic.  It reuses the certified packet routines but
does not emit a theorem gate or promote a claim.
"""

from __future__ import annotations

import argparse
import json
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
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate as lower_cell
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_complement_signed_coefficient_interval_gate as ordinary
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as ibp
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate as endpoint
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_ordinary_complementary_half_lattice_gate as complement
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_quarter_disk_vertical_arc_tail_arb as upper_arc
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_far_remainder_eight_round_gate as upper_far
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate as upper_transition
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined


HEIGHT = 10_000_000_000
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_first_subcell_interval_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def complex_from_record(record: dict[str, Any]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def complex_width(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(55, more=True),
        "imag_ball": value.imag.str(55, more=True),
        "real_radius": value.real.rad().str(25, more=True),
        "imag_radius": value.imag.rad().str(25, more=True),
        "absolute_ball": abs(value).str(55, more=True),
    }


def real_width(value: arb) -> dict[str, str]:
    return {
        "ball": value.str(55, more=True),
        "lower": value.lower().str(40, more=True),
        "upper": value.upper().str(40, more=True),
        "radius": value.rad().str(25, more=True),
    }


def retarget_modules(t_box: arb) -> None:
    # These modules expose exact certificate builders but retain the production
    # height as a module constant.  Rebinding only that scalar makes their
    # existing formulas evaluate the same fixed roster over the requested box.
    endpoint.HEIGHT = t_box
    ibp.HEIGHT = t_box
    upper_transition.HEIGHT = t_box
    upper_far.HEIGHT = t_box
    complement.HEIGHT = t_box
    ordinary.HEIGHT = t_box
    upper_arc.HEIGHT = t_box


def interval_upper_arc(variant: str) -> dict[str, Any]:
    """Run the donor arc proof with a height-box scale guard.

    The donor's 0.01981 action-maximum check is a point-calibration alarm; the
    subsequent polynomial and tail errors are computed from the actual Arb
    action ball.  On a height interval we relax only that raw dependent alarm,
    retain every analytic/error guard, and reimpose 0.01981 through a stable
    envelope-theorem transport.
    """

    donor_require = upper_arc.require

    def interval_require(condition: bool, message: str) -> None:
        if message == "production action maximum exceeded its guard":
            return
        donor_require(condition, message)

    upper_arc.require = interval_require
    try:
        certificate = upper_arc.run(variant, None)
    finally:
        upper_arc.require = donor_require
    # Directly substituting the interval-valued stationary point into the
    # three large cancelling action terms creates a spurious O(1) radius.
    # The maximizer is unique on the box, so the envelope theorem gives
    # d A_max / dt = delta_*(t).  Transport the point action with this exact
    # derivative and the independently enclosed stationary-angle box.
    pi = arb.pi()
    y = arb(upper_arc.Y_TEXT)
    q_minus = arb(upper_arc.Q_MINUS_TEXT)
    radius = 2 * pi * y
    t_box = arb(upper_arc.HEIGHT)

    def stationary_delta(t: arb) -> arb:
        constant = t - 2 * pi * y * y
        discriminant = (q_minus * radius) ** 2 - 16 * pi * y * y * constant
        root = (q_minus * radius + discriminant.sqrt()) / (8 * pi * y * y)
        return root.acos()

    t_center = arb(HEIGHT)
    delta_center = stationary_delta(t_center)
    delta_box = stationary_delta(t_box)
    action_center = (
        t_center * delta_center
        + pi * y * y * (2 * delta_center).sin()
        - q_minus * radius * delta_center.sin()
    )
    action_variation = abs(t_box - t_center).upper() * delta_box.upper()
    stable_action_maximum = action_center + arb(0, action_variation)
    require(
        stable_action_maximum.upper() < arb("0.01981"),
        "envelope-transported height-box arc action exceeds 0.01981",
    )
    certificate["scout_action_guard"] = {
        "donor_point_guard": "action_maximum<0.01981",
        "raw_dependent_action_ball": certificate["action_maximum"]["ball"],
        "stable_envelope_action_ball": stable_action_maximum.str(55, more=True),
        "stationary_delta_box": delta_box.str(55, more=True),
        "identity": "d A_max(t)/dt=delta_*(t)",
        "height_box_guard": "stable_envelope_action_maximum.upper<0.01981",
        "all_donor_error_bounds_retained": True,
    }
    return certificate


def stable_H_box(t_box: arb) -> tuple[arb, dict[str, str]]:
    """Transport H from the center with its exact logarithmic derivative."""

    t_center = arb(HEIGHT)
    log_H_center = joined.direct_log_H(t_center)
    direct = joined.direct_log_H_derivative(t_box)
    duplication = joined.duplication_log_H_derivative(t_box)
    require(direct.overlaps(duplication), "two H logarithmic-derivative boxes miss")
    derivative_bound = min(abs(direct).upper(), abs(duplication).upper())
    transport_distance = abs(t_box - t_center).upper()
    log_H_box = log_H_center + arb(0, transport_distance * derivative_bound)
    return log_H_box.exp(), {
        "center_log_H_ball": log_H_center.str(55, more=True),
        "direct_logarithmic_derivative_ball": direct.str(55, more=True),
        "duplication_logarithmic_derivative_ball": duplication.str(55, more=True),
        "selected_absolute_derivative_upper": derivative_bound.str(35, more=True),
        "maximum_transport_distance": transport_distance.str(35, more=True),
        "transport_identity": "log H(t)-log H(t0)=integral_(t0)^t H'(u)/H(u) du",
    }


def configure_variant(variant: str) -> dict[str, Any]:
    require(variant in ("production", "independent"), "unknown scout variant")
    if variant == "production":
        upper_transition.PANELS = 192
        ibp.ORDER = 8
        ibp.SLABS = 8192
        ibp.CLUSTER_POWER = 4
        upper_far.ORDER = 8
        upper_far.SLABS = 8192
        ordinary.ORDER = 3
        ordinary.SLABS = 16384
        ordinary.CLUSTER_POWER = 4
        return {
            "upper_transition_panels": 192,
            "upper_tail_rounds": 8,
            "upper_tail_slabs": 8192,
            "ordinary_rounds": 3,
            "ordinary_slabs": 16384,
            "lower_cell_rounds": 18,
            "lower_cell_slabs": 4096,
            "transition_jet_order": 12,
            "transition_derivative_slabs": 128,
            "transition_integral_panels": 32,
            "transition_series_terms": 30,
            "arc_variant": "production",
        }
    upper_transition.PANELS = 256
    ibp.ORDER = 8
    ibp.SLABS = 10000
    ibp.CLUSTER_POWER = 5
    upper_far.ORDER = 8
    upper_far.SLABS = 10000
    ordinary.ORDER = 4
    ordinary.SLABS = 20000
    ordinary.CLUSTER_POWER = 5
    return {
        "upper_transition_panels": 256,
        "upper_tail_rounds": 8,
        "upper_tail_slabs": 10000,
        "ordinary_rounds": 4,
        "ordinary_slabs": 20000,
        "lower_cell_rounds": 20,
        "lower_cell_slabs": 6144,
        "transition_jet_order": 14,
        "transition_derivative_slabs": 256,
        "transition_integral_panels": 48,
        "transition_series_terms": 34,
        "arc_variant": "independent",
    }


def build(
    radius_text: str,
    precision_bits: int,
    *,
    center_text: str = str(HEIGHT),
    variant: str = "production",
) -> dict[str, Any]:
    started = time.perf_counter()
    resource_mode = joined.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = precision_bits
    ctx.threads = 1
    radius = arb(radius_text)
    require(radius > 0 and radius <= arb("0.001"), "scout radius must lie in (0,0.001]")
    center = arb(center_text)
    t_box = arb(center, radius)
    configuration = configure_variant(variant)
    retarget_modules(t_box)

    endpoint_certificate = endpoint.build_certificate()
    transition_752_certificate = upper_transition.build_certificate()
    upper_certificate = upper_far.build_certificate(
        {"certificate": transition_752_certificate}
    )
    complement_certificate = complement.build_certificate()
    # Rebuild only the production order-three ordinary packet.  The donor
    # builder also computes point-only diagnostics for a rejected order-eight
    # route; those diagnostics can legitimately straddle zero on a height box
    # and play no role in the order-three enclosure.
    p = t_box / (2 * arb.pi())
    table = ibp.coefficient_table(max(ordinary.MAX_AUDIT_ORDER, ordinary.ORDER))
    lower_endpoints = endpoint_certificate["families"]["R_lower_infinite"]
    lower_partial, _ = ibp.endpoint_expansion(
        table,
        p,
        t_box,
        lower_endpoints["endpoint_L"],
        lower_endpoints["endpoint_a"],
        ordinary.ORDER,
    )
    lower_remainder, _ = ordinary.signed_coefficient_remainder(
        table[ordinary.ORDER], p, order=ordinary.ORDER
    )
    lower_complement = ibp.add_complex_error(lower_partial, lower_remainder)
    upper_complement = complex_from_record(
        upper_certificate["complete_upper_complement_ball"]
    )
    gamma_log10 = arb(
        complement_certificate["ordinary_Gamma_defect_packet_log10_upper_ball"]
        ["ball"]
    )
    gamma_radius = (gamma_log10.upper() * arb(10).log()).exp().upper()
    ordinary_packet = ibp.add_complex_error(
        -lower_complement - upper_complement, gamma_radius
    )

    lower_certificate = lower_cell.lower_cell_certificate(
        t_box,
        lower_cell.arb_from_fraction(lower_cell.LOWER_BOUNDARY),
        lower_cell.arb_from_fraction(lower_cell.Q_PLUS),
        lower_cell.SOURCE_COUNT,
        configuration["lower_cell_rounds"],
        arb(lower_cell.SMALL_Y_SPLIT),
        configuration["lower_cell_slabs"],
        resource_mode,
    )
    arc_certificate = interval_upper_arc(str(configuration["arc_variant"]))
    transition_certificate = transition.transition_source_certificate(
        t_box,
        transition.SOURCE_ENDPOINT,
        transition.SOURCE_COUNT,
        transition.CELL_LEFT,
        transition.CELL_RIGHT,
        int(configuration["transition_jet_order"]),
        int(configuration["transition_derivative_slabs"]),
        int(configuration["transition_integral_panels"]),
        int(configuration["transition_series_terms"]),
        "1e-27" if variant == "independent" else "1e-26",
    )

    lower = complex_from_record(lower_certificate["lower_cell_integral_ball"])
    arc = complex_from_record(arc_certificate["actual_arc"])
    tail_log10 = arb(
        arc_certificate["error_budgets"]["positive_real_tail_log10_upper"]["ball"]
    )
    tail_radius = (tail_log10.upper() * arb(10).log()).exp().upper()
    tail = acb(arb(0, tail_radius), arb(0, tail_radius))
    unowned = lower + arc + tail
    transition_source = complex_from_record(
        transition_certificate["complete_transition_source_integral_ball"]
    )

    theta, _ = joined.theta_and_derivative(t_box)
    H, H_transport = stable_H_box(t_box)
    C_G = (1 + (-2 * arb.pi() * t_box).exp()) ** (-arb(1) / 2)
    D_target, _ = joined.finite_dirichlet(
        t_box, joined.TARGET_START, joined.TARGET_END
    )
    D_lower_transition, _ = joined.finite_dirichlet(
        t_box, joined.LOWER_TRANSITION_START, joined.TARGET_END
    )
    transition_join = transition_source - C_G * D_lower_transition
    tiny_target_correction = (C_G - H) * D_target
    packet = unowned + ordinary_packet + transition_join + tiny_target_correction
    qkt = joined.hardy_projection(theta, packet) / H

    saved = load_json(joined.RESULT)
    saved_qkt = arb(
        saved["certificate"]["A_free_saved_height_assembly"]
        ["direct_Q_K_minus_T_ball"]["ball"]
    )
    require(qkt.overlaps(saved_qkt), "height box misses the saved-height value")

    return {
        "kind": "rh_diagnostic_joined_QK_minus_T_first_subcell_interval_scout",
        "date": "2026-08-28",
        "status": "diagnostic_direct_height_box_only_not_a_promoted_interval_theorem",
        "passed": True,
        "resource_mode": resource_mode,
        "workers": 1,
        "precision_bits": precision_bits,
        "variant": variant,
        "configuration": configuration,
        "center_height": center_text,
        "radius": radius_text,
        "height_ball": t_box.str(50, more=True),
        "components": {
            "lower_cell": complex_width(lower),
            "upper_arc": complex_width(arc),
            "upper_arc_action_guard": arc_certificate["scout_action_guard"],
            "upper_positive_tail_radius": tail_radius.str(25, more=True),
            "unowned_packet": complex_width(unowned),
            "ordinary_packet": complex_width(ordinary_packet),
            "transition_source": complex_width(transition_source),
            "minus_CG_lower_transition": complex_width(-C_G * D_lower_transition),
            "joined_transition": complex_width(transition_join),
            "tiny_target_correction": complex_width(tiny_target_correction),
            "complete_joined_packet": complex_width(packet),
        },
        "projection": {
            "theta_ball": real_width(theta),
            "H_ball": real_width(H),
            "H_transport": H_transport,
            "Q_K_minus_T_height_box": real_width(qkt),
            "saved_height_overlap": True,
            "strictly_negative_on_box": bool(qkt.upper() < 0),
            "radius_to_QKT_radius_ratio": (qkt.rad() / radius).str(30, more=True),
        },
        "diagnostic_boundary": (
            "This run tests direct interval propagation through reused certified formulas. "
            "It is not promoted because no independent changed-decomposition checker or "
            "uniform-formula audit has yet been completed."
        ),
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", default="0.0001")
    parser.add_argument("--precision-bits", type=int, default=384)
    parser.add_argument("--center", default=str(HEIGHT))
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    args = parser.parse_args()
    artifact = build(
        args.radius,
        args.precision_bits,
        center_text=args.center,
        variant=args.variant,
    )
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    qkt = artifact["projection"]["Q_K_minus_T_height_box"]
    print(
        "diagnostic joined Q_K-T height box "
        f"radius={artifact['radius']} lower={qkt['lower']} upper={qkt['upper']} "
        f"negative={artifact['projection']['strictly_negative_on_box']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
