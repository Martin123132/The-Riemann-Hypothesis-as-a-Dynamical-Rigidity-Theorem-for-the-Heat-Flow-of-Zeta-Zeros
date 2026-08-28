#!/usr/bin/env python3
"""Scout the cancellation-preserving lower-cell plus ordinary K packet.

This diagnostic applies L_t=d_t+i*theta'-H'/H inside the lower finite cell,
both ordinary complementary tails, and the grouped 752-label upper collar.
All endpoint partials are joined before their rigorous remainders are added.
It is not yet a theorem gate for K_(L+O) or complete K_T.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
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

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_complement_signed_coefficient_interval_gate as ordinary
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_eight_round_remainder_gate as ibp
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_nonstationary_complement_root_unity_hurwitz_endpoint_gate as endpoint
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_far_remainder_eight_round_gate as upper_far
import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_upper_complement_closed_752_label_exterior_transition_layer_gate as upper_group
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_lower_finite_cell_K_first_subcell_scout as lower_k


HEIGHT = 10_000_000_000
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_lower_ordinary_K_first_subcell_scout"
)
RESULT = ROOT / "work" / "rh_compute" / "results" / f"{STEM}.json"
EVENT_ATLAS = (
    ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_R_after_A_local_height_event_atlas_gate.json"
)

LogMonomial = tuple[int, int, int]
LogLevel = dict[int, dict[LogMonomial, Fraction]]
LogTable = list[LogLevel]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def ordinary_logarithmic_table(order: int) -> LogTable:
    """Exact tables for T^n[log(y)y^-1/2], T(f)=(f/(q-y-p/y))'."""

    require(order >= 1, "positive ordinary logarithmic order required")
    table: LogTable = [{0: {(0, -1, 1): Fraction(1)}}]
    for _ in range(order):
        next_level: LogLevel = {}
        for denominator_power, terms in table[-1].items():
            for (p_power, exponent2, log_power), coefficient in terms.items():
                updates = [
                    (
                        denominator_power + 1,
                        (p_power, exponent2 - 2, log_power),
                        coefficient * Fraction(exponent2, 2),
                    ),
                    (
                        denominator_power + 2,
                        (p_power, exponent2, log_power),
                        coefficient * (denominator_power + 1),
                    ),
                    (
                        denominator_power + 2,
                        (p_power + 1, exponent2 - 4, log_power),
                        -coefficient * (denominator_power + 1),
                    ),
                ]
                if log_power:
                    updates.append(
                        (
                            denominator_power + 1,
                            (p_power, exponent2 - 2, log_power - 1),
                            coefficient * log_power,
                        )
                    )
                for new_power, monomial, value in updates:
                    destination = next_level.setdefault(new_power, {})
                    destination[monomial] = destination.get(monomial, Fraction(0)) + value
        table.append(
            {
                power: {monomial: value for monomial, value in terms.items() if value}
                for power, terms in next_level.items()
            }
        )
    return table


def evaluate_log_coefficient(
    terms: dict[LogMonomial, Fraction] | None, p: arb, y: arb
) -> arb:
    if not terms:
        return arb(0)
    value = arb(0)
    log_y = y.log()
    for (p_power, exponent2, log_power), coefficient in terms.items():
        value += (
            ibp.arb_from_fraction(coefficient)
            * p**p_power
            * ibp.half_power(y, exponent2)
            * log_y**log_power
        )
    return value


def weighted_coefficient(
    base_terms: dict[tuple[int, int], Fraction] | None,
    log_terms: dict[LogMonomial, Fraction] | None,
    p: arb,
    y: arb,
    hardy_shift: acb,
) -> acb:
    base = ibp.evaluate_coefficient(base_terms, p, y) if base_terms else arb(0)
    logarithmic = evaluate_log_coefficient(log_terms, p, y)
    return hardy_shift * base - acb(0, 1) * logarithmic


def endpoint_level_sum(
    base_level: dict[int, dict[tuple[int, int], Fraction]],
    log_level: LogLevel,
    p: arb,
    t: arb,
    endpoint_record: dict[str, Any],
    hardy_shift: acb,
) -> acb:
    y = arb(endpoint_record["endpoint_y"])
    rows = {row["power"]: row for row in endpoint_record["rows"]}
    value = acb(0)
    for denominator_power in set(base_level) | set(log_level):
        required_power = denominator_power + 1
        require(required_power in rows, f"missing weighted endpoint power {required_power}")
        coefficient = weighted_coefficient(
            base_level.get(denominator_power),
            log_level.get(denominator_power),
            p,
            y,
            hardy_shift,
        )
        value += coefficient * ibp.saved_complex(
            rows[required_power]["phase_weighted_sum_ball"]
        )
    return ibp.endpoint_carrier(t, y) * value


def endpoint_expansion(
    base_table: list[dict[int, dict[tuple[int, int], Fraction]]],
    log_table: LogTable,
    p: arb,
    t: arb,
    left_endpoint: dict[str, Any],
    right_endpoint: dict[str, Any],
    order: int,
    hardy_shift: acb,
) -> tuple[acb, list[dict[str, Any]]]:
    partial = acb(0)
    terms: list[dict[str, Any]] = []
    imaginary_scale = acb(0, 1) * 2 * arb.pi()
    for n in range(order):
        left_value = endpoint_level_sum(
            base_table[n], log_table[n], p, t, left_endpoint, hardy_shift
        )
        right_value = endpoint_level_sum(
            base_table[n], log_table[n], p, t, right_endpoint, hardy_shift
        )
        factor = (-1) ** n * imaginary_scale ** (-(n + 1))
        left_contribution = -factor * left_value
        right_contribution = factor * right_value
        term = left_contribution + right_contribution
        partial += term
        terms.append(
            {
                "round": n,
                "left_integral_contribution": lower_k.complex_record(left_contribution),
                "right_integral_contribution": lower_k.complex_record(right_contribution),
                "signed_term_ball": lower_k.complex_record(term),
            }
        )
    return partial, terms


def weighted_remainder(
    base_level: dict[int, dict[tuple[int, int], Fraction]],
    log_level: LogLevel,
    p: arb,
    left: arb,
    right: arb,
    q0: arb,
    orientation: str,
    finite_count: int | None,
    slabs: int,
    cluster_power: int,
    order: int,
    hardy_shift: acb,
) -> tuple[arb, dict[str, Any]]:
    require(right < p.sqrt(), "ordinary weighted interval crossed the Q minimum")
    ibp.CLUSTER_POWER = cluster_power
    cluster_side = "left" if orientation == "upper" else "right"
    integral = arb(0)
    minimum_gap: arb | None = None
    for index in range(slabs):
        slab_left = ibp.clustered_edge(left, right, index, slabs, cluster_side)
        slab_right = ibp.clustered_edge(left, right, index + 1, slabs, cluster_side)
        midpoint = (slab_left + slab_right) / 2
        y_ball = arb(midpoint, (slab_right - slab_left) / 2)
        if orientation == "upper":
            delta = (q0 - ibp.q_curve(slab_left, p)).lower()
        elif orientation == "lower":
            delta = (ibp.q_curve(slab_right, p) - q0).lower()
        else:
            raise RuntimeError(f"unknown weighted orientation: {orientation}")
        require(delta > 0, f"ordinary weighted gap lost in slab {index}")
        if minimum_gap is None or delta < minimum_gap:
            minimum_gap = delta
        slab_majorant = arb(0)
        for power in set(base_level) | set(log_level):
            coefficient = weighted_coefficient(
                base_level.get(power), log_level.get(power), p, y_ball, hardy_shift
            )
            slab_majorant += abs(coefficient).upper() * ibp.denominator_bound(
                delta, power, finite_count
            )
        integral += (slab_right - slab_left) * slab_majorant
    require(minimum_gap is not None and integral.is_finite(), "weighted remainder failed")
    scaled = (integral / (2 * arb.pi()) ** order).upper()
    return scaled, {
        "slabs": slabs,
        "cluster_side": cluster_side,
        "cluster_power": cluster_power,
        "integration_rounds": order,
        "minimum_denominator_gap_ball": lower_k.real_record(minimum_gap),
        "unscaled_integral_majorant_ball": lower_k.real_record(integral.upper()),
        "scaled_remainder_bound_ball": lower_k.real_record(scaled),
    }


def weighted_family(
    *,
    name: str,
    base_table: list[dict[int, dict[tuple[int, int], Fraction]]],
    log_table: LogTable,
    p: arb,
    t: arb,
    endpoints: dict[str, Any],
    left_key: str,
    right_key: str,
    left: arb,
    right: arb,
    q0: arb,
    orientation: str,
    finite_count: int | None,
    slabs: int,
    cluster_power: int,
    order: int,
    hardy_shift: acb,
) -> dict[str, Any]:
    partial, terms = endpoint_expansion(
        base_table,
        log_table,
        p,
        t,
        endpoints[left_key],
        endpoints[right_key],
        order,
        hardy_shift,
    )
    remainder, census = weighted_remainder(
        base_table[order],
        log_table[order],
        p,
        left,
        right,
        q0,
        orientation,
        finite_count,
        slabs,
        cluster_power,
        order,
        hardy_shift,
    )
    ball = ibp.add_complex_error(partial, remainder)
    return {
        "name": name,
        "endpoint_partial_sum_ball": lower_k.complex_record(partial),
        "endpoint_terms": terms,
        "remainder_census": census,
        "weighted_tail_ball": lower_k.complex_record(ball),
    }


def grouped_752_K(
    t: arb, hardy_shift: acb, panels: int
) -> tuple[acb, dict[str, Any]]:
    base_integrand = upper_group.base.grouped_integrand_factory(
        t=t,
        lower=upper_group.LOWER,
        q0=upper_group.Q_PLUS,
        count=upper_group.PACKET_COUNT,
    )

    def integrand(x: acb, analytic: bool) -> acb:
        y = acb(upper_group.LOWER) + x
        weight = hardy_shift - acb(0, 1) * y.log()
        return weight * base_integrand(x, analytic)

    value, rows = upper_group.base.integrate_panels(
        integrand, upper_group.SPLIT - upper_group.LOWER, panels
    )
    require(value.is_finite(), "grouped weighted 752-label integral is nonfinite")
    return value, {
        "panels": panels,
        "weighted_grouped_752_ball": lower_k.complex_record(value),
        "panel_rows": rows,
        "endpoint_zero_guard": (
            "the 752-label geometric amplitude vanishes exactly at L and C; "
            "multiplication by i(theta'-log(y))-H'/H preserves both zeros"
        ),
    }


def gamma_defect_K_radius(t: arb, hardy_shift: acb) -> dict[str, Any]:
    pi = arb.pi()
    decay = (-2 * pi * t).exp()
    root = (1 + decay).sqrt()
    C_G = 1 / root
    one_minus_C_G = decay / (root * (1 + root))
    C_G_prime = C_G * pi * decay / (1 + decay)
    mode_sum_bound = 2 * (
        arb(joined.LOWER_TRANSITION_START - 1).sqrt()
        - arb(joined.TARGET_START - 1).sqrt()
    )
    log_bound = arb(joined.LOWER_TRANSITION_START - 1).log()
    radius = (
        one_minus_C_G * mode_sum_bound * (abs(hardy_shift).upper() + log_bound)
        + C_G_prime * mode_sum_bound
    ).upper()
    return {
        "C_G_ball": lower_k.real_record(C_G),
        "one_minus_C_G_ball": lower_k.real_record(one_minus_C_G),
        "C_G_prime_ball": lower_k.real_record(C_G_prime),
        "finite_mode_absolute_bound": lower_k.real_record(mode_sum_bound),
        "Hardy_operator_component_radius": lower_k.real_record(radius),
        "identity": (
            "L_t[(1-C_G)D_ord]=(1-C_G)(D_ord'+(i*theta'-H'/H)D_ord)"
            "-C_G'D_ord"
        ),
    }


def configure(variant: str) -> dict[str, Any]:
    require(variant in ("production", "independent"), "unknown lower-ordinary variant")
    if variant == "production":
        return {
            "lower_order": 18,
            "lower_split": "300",
            "lower_slabs": 4096,
            "ordinary_lower_order": 3,
            "ordinary_lower_slabs": 16384,
            "upper_order": 8,
            "upper_slabs": 8192,
            "cluster_power": 4,
            "grouped_panels": 192,
            "endpoint_max_power": 16,
        }
    return {
        "lower_order": 20,
        "lower_split": "280",
        "lower_slabs": 6144,
        "ordinary_lower_order": 4,
        "ordinary_lower_slabs": 20000,
        "upper_order": 9,
        "upper_slabs": 10000,
        "cluster_power": 5,
        "grouped_panels": 256,
        "endpoint_max_power": 20,
    }


def build(radius_text: str, precision_bits: int, variant: str) -> dict[str, Any]:
    started = time.perf_counter()
    resource_mode = joined.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = precision_bits
    ctx.threads = 1
    config = configure(variant)
    radius = arb(radius_text)
    require(radius >= 0 and radius <= arb("0.001"), "radius must lie in [0,0.001]")
    t = arb(arb(HEIGHT), radius)
    event = load_json(EVENT_ATLAS)["certificate"]["primary_open_cell"]
    require(
        arb(event["lower_ball"]["ball"]).upper() < t.lower()
        and t.upper() < arb(event["upper_ball"]["ball"]).lower(),
        "height box left the fixed-roster event cell",
    )

    theta, theta_prime = joined.theta_and_derivative(t)
    H_log_prime = joined.direct_log_H_derivative(t)
    require(
        H_log_prime.overlaps(joined.duplication_log_H_derivative(t)),
        "H logarithmic-derivative formulas miss",
    )
    hardy_shift = acb(-H_log_prime, theta_prime)
    p = t / (2 * arb.pi())

    endpoint.HEIGHT = t
    endpoint.MAX_POWER = int(config["endpoint_max_power"])
    endpoint_certificate = endpoint.build_certificate()
    lower_endpoints = endpoint_certificate["families"]["R_lower_infinite"]
    upper_endpoints = upper_far.endpoint_families(p)

    lower_certificate = lower_k.weighted_lower_cell_certificate(
        t,
        lower_k.lower.arb_from_fraction(lower_k.lower.LOWER_BOUNDARY),
        lower_k.lower.arb_from_fraction(lower_k.lower.Q_PLUS),
        lower_k.lower.SOURCE_COUNT,
        int(config["lower_order"]),
        arb(config["lower_split"]),
        int(config["lower_slabs"]),
        hardy_shift,
    )
    lower_ball = lower_k.lower.complex_from_record(
        lower_certificate["weighted_lower_cell_ball"]
    )
    lower_partial = lower_k.lower.complex_from_record(
        lower_certificate["endpoint_partial_sum_ball"]
    )
    lower_remainder = arb(lower_certificate["total_remainder_bound_ball"]["ball"])

    maximum_order = max(
        int(config["ordinary_lower_order"]), int(config["upper_order"])
    )
    base_table = ibp.coefficient_table(maximum_order)
    log_table = ordinary_logarithmic_table(maximum_order)
    lower_tail = weighted_family(
        name="T_lower_K",
        base_table=base_table,
        log_table=log_table,
        p=p,
        t=t,
        endpoints=lower_endpoints,
        left_key="endpoint_L",
        right_key="endpoint_a",
        left=ordinary.L,
        right=ordinary.A,
        q0=ordinary.Q_LOW,
        orientation="lower",
        finite_count=None,
        slabs=int(config["ordinary_lower_slabs"]),
        cluster_power=int(config["cluster_power"]),
        order=int(config["ordinary_lower_order"]),
        hardy_shift=hardy_shift,
    )
    finite_upper = weighted_family(
        name="F_C_K",
        base_table=base_table,
        log_table=log_table,
        p=p,
        t=t,
        endpoints=upper_endpoints["F_C_finite_752"],
        left_key="endpoint_C",
        right_key="endpoint_a",
        left=upper_far.C,
        right=upper_far.A,
        q0=upper_far.Q_PLUS,
        orientation="upper",
        finite_count=upper_far.PACKET_COUNT,
        slabs=int(config["upper_slabs"]),
        cluster_power=int(config["cluster_power"]),
        order=int(config["upper_order"]),
        hardy_shift=hardy_shift,
    )
    remote_upper = weighted_family(
        name="I_L_K",
        base_table=base_table,
        log_table=log_table,
        p=p,
        t=t,
        endpoints=upper_endpoints["I_L_infinite_upper"],
        left_key="endpoint_L",
        right_key="endpoint_a",
        left=upper_far.L,
        right=upper_far.A,
        q0=upper_far.Q_REMOTE,
        orientation="upper",
        finite_count=None,
        slabs=int(config["upper_slabs"]),
        cluster_power=int(config["cluster_power"]),
        order=int(config["upper_order"]),
        hardy_shift=hardy_shift,
    )
    grouped_upper, grouped_certificate = grouped_752_K(
        t, hardy_shift, int(config["grouped_panels"])
    )

    lower_tail_partial = ibp.saved_complex(lower_tail["endpoint_partial_sum_ball"])
    finite_upper_partial = ibp.saved_complex(finite_upper["endpoint_partial_sum_ball"])
    remote_upper_partial = ibp.saved_complex(remote_upper["endpoint_partial_sum_ball"])
    lower_tail_remainder = arb(
        lower_tail["remainder_census"]["scaled_remainder_bound_ball"]["ball"]
    )
    finite_upper_remainder = arb(
        finite_upper["remainder_census"]["scaled_remainder_bound_ball"]["ball"]
    )
    remote_upper_remainder = arb(
        remote_upper["remainder_census"]["scaled_remainder_bound_ball"]["ball"]
    )
    gamma = gamma_defect_K_radius(t, hardy_shift)
    gamma_radius = arb(gamma["Hardy_operator_component_radius"]["ball"])

    ordinary_partial = -(
        lower_tail_partial + grouped_upper + finite_upper_partial + remote_upper_partial
    )
    ordinary_error = (
        lower_tail_remainder
        + finite_upper_remainder
        + remote_upper_remainder
        + gamma_radius
    ).upper()
    ordinary_ball = ibp.add_complex_error(ordinary_partial, ordinary_error)

    joined_partial = lower_partial + ordinary_partial
    joined_error = (lower_remainder + ordinary_error).upper()
    joined_ball = ibp.add_complex_error(joined_partial, joined_error)
    component_sum = lower_ball + ordinary_ball
    require(joined_ball.overlaps(component_sum), "joined and component lower-ordinary K miss")
    H, H_transport = lower_k.stable_H_box(t)
    ordinary_projection = joined.hardy_projection(theta, ordinary_ball) / H
    joined_projection = joined.hardy_projection(theta, joined_ball) / H

    minimum_rounds = min(
        int(config["lower_order"]),
        int(config["ordinary_lower_order"]),
        int(config["upper_order"]),
    )
    L_currents: list[dict[str, Any]] = []
    for n in range(minimum_rounds):
        finite_current = lower_k.lower.complex_from_record(
            lower_certificate["endpoint_terms"][n]["term_ball"]
        )
        lower_complement_current = -ibp.saved_complex(
            lower_tail["endpoint_terms"][n]["left_integral_contribution"]
        )
        upper_remote_current = -ibp.saved_complex(
            remote_upper["endpoint_terms"][n]["left_integral_contribution"]
        )
        joined_current = (
            finite_current + lower_complement_current + upper_remote_current
        )
        triangle = (
            abs(finite_current)
            + abs(lower_complement_current)
            + abs(upper_remote_current)
        )
        L_currents.append(
            {
                "round": n,
                "joined_complete_lattice_current": lower_k.complex_record(joined_current),
                "separate_triangle_upper": triangle.upper().str(35, more=True),
                "joined_absolute_upper": abs(joined_current).upper().str(35, more=True),
            }
        )

    return {
        "kind": "rh_diagnostic_joined_lower_ordinary_K_first_subcell_scout",
        "date": "2026-08-28",
        "status": "diagnostic_lower_plus_ordinary_Hardy_operator_not_complete_K_T",
        "passed": True,
        "resource_mode": resource_mode,
        "workers": 1,
        "precision_bits": precision_bits,
        "variant": variant,
        "radius": radius_text,
        "height_ball": t.str(58, more=True),
        "configuration": config,
        "exact_identities": {
            "ordinary_packet": "O_join=-T_lower-T_upper+Gamma_defect",
            "upper_decomposition": "T_upper=G_752+F_C+I_L",
            "joined_coordinate": (
                "K_(L+O)=L_t[V_L-T_lower-G_752-F_C-I_L+Gamma_defect]"
            ),
            "endpoint_rule": (
                "join all finite, lower-complement, and upper-complement endpoint partials "
                "before adding the separately certified recurrence remainders"
            ),
        },
        "phase_and_prefactor": {
            "theta": lower_k.real_record(theta),
            "theta_prime": lower_k.real_record(theta_prime),
            "H": lower_k.real_record(H),
            "H_transport": H_transport,
            "H_log_prime": lower_k.real_record(H_log_prime),
            "Hardy_shift": lower_k.complex_record(hardy_shift),
        },
        "lower_finite_cell_K": lower_k.complex_record(lower_ball),
        "lower_complement_K": lower_tail,
        "grouped_752_K": grouped_certificate,
        "finite_upper_K": finite_upper,
        "remote_upper_K": remote_upper,
        "Gamma_defect_K": gamma,
        "ordinary_K": lower_k.complex_record(ordinary_ball),
        "ordinary_Q_derivative_contribution": lower_k.real_record(ordinary_projection),
        "joined_lower_ordinary_endpoint_partial": lower_k.complex_record(joined_partial),
        "joined_lower_ordinary_total_error": lower_k.real_record(joined_error),
        "joined_lower_ordinary_K": lower_k.complex_record(joined_ball),
        "component_sum_overlap": True,
        "joined_lower_ordinary_Q_derivative_contribution": lower_k.real_record(
            joined_projection
        ),
        "complete_lattice_L_endpoint_currents": L_currents,
        "endpoint_certificate": endpoint_certificate,
        "elapsed_seconds": time.perf_counter() - started,
        "diagnostic_boundary": (
            "This diagnostic covers only V_L+O_join on I_1. It does not yet constitute an "
            "independently checked theorem gate and omits the already separate transition--arc "
            "pair, positive-real tail derivative, and tiny target correction required for K_T."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", default="0.0001")
    parser.add_argument("--precision", type=int, default=384)
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    args = parser.parse_args()
    artifact = build(args.radius, args.precision, args.variant)
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    contribution = artifact["joined_lower_ordinary_Q_derivative_contribution"]
    print(
        "completed joined lower-plus-ordinary K diagnostic; "
        f"Hardy_contribution={contribution['ball']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
