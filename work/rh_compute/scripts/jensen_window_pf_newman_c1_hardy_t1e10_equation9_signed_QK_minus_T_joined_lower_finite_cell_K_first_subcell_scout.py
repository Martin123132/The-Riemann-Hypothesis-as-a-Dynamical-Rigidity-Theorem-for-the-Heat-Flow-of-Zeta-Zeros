#!/usr/bin/env python3
"""Scout the lower finite-cell contribution to the joined Hardy derivative.

The certified lower-cell endpoint recurrence starts from y^(-1/2).  This
diagnostic extends it to

    [i*(theta'(t)-log(y))-H'(t)/H(t)] y^(-1/2)

before the reversed finite-roster endpoint collapse.  It does not yet join
the lower cell to the ordinary complementary packet and is not a theorem
gate for the complete K_T.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
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

import jensen_window_pf_newman_c1_hardy_t1e10_equation4_joined_packet_lower_cell_reversed_roster_endpoint_gate as lower
import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_joined_height_derivative_identity_gate as joined


HEIGHT = 10_000_000_000
STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_"
    "joined_lower_finite_cell_K_first_subcell_scout"
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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(58, more=True),
        "imag_ball": value.imag.str(58, more=True),
        "absolute_ball": abs(value).str(58, more=True),
        "real_radius": value.real.rad().str(28, more=True),
        "imag_radius": value.imag.rad().str(28, more=True),
    }


def real_record(value: arb) -> dict[str, str]:
    return {
        "ball": value.str(58, more=True),
        "lower": value.lower().str(45, more=True),
        "upper": value.upper().str(45, more=True),
        "radius": value.rad().str(28, more=True),
    }


def logarithmic_coefficient_table(order: int) -> LogTable:
    """Exact tables for T^n[log(y)y^(-1/2)], T(f)=(f/d)'."""

    require(order >= 1, "positive logarithmic recurrence order required")
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
                        -coefficient * (denominator_power + 1),
                    ),
                    (
                        denominator_power + 2,
                        (p_power + 1, exponent2 - 4, log_power),
                        coefficient * (denominator_power + 1),
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
    terms: dict[LogMonomial, Fraction], p: arb, y: arb
) -> arb:
    value = arb(0)
    log_y = y.log()
    for (p_power, exponent2, log_power), coefficient in terms.items():
        value += (
            lower.arb_from_fraction(coefficient)
            * p**p_power
            * lower.half_power(y, exponent2)
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
    base = lower.evaluate_coefficient(base_terms, p, y) if base_terms else arb(0)
    logarithmic = evaluate_log_coefficient(log_terms, p, y) if log_terms else arb(0)
    return hardy_shift * base - acb(0, 1) * logarithmic


def absolute_log_moment(split: arb, denominator2: int, log_power: int) -> arb:
    """Return integral_0^split y^(denominator2/2-1)|log y|^ell dy."""

    require(split > 0, "the explicit absolute-log formula expects split>0")
    require(denominator2 > 0, "nonintegrable small-y logarithmic monomial")
    require(log_power in (0, 1), "unexpected logarithmic power")
    exponent = arb(denominator2) / 2
    if log_power == 0:
        return split**exponent / exponent
    inverse_square = exponent ** (-2)
    if split.upper() <= 1:
        return split**exponent * (-split.log() / exponent + inverse_square)
    require(split.lower() >= 1, "absolute-log split straddles one")
    return (
        split**exponent * (split.log() / exponent - inverse_square)
        + 2 * inverse_square
    )


def logarithmic_small_y_bound(
    level: LogLevel,
    p: arb,
    count: int,
    split: arb,
) -> arb:
    value = arb(0)
    for denominator_power, terms in level.items():
        for (p_power, exponent2, log_power), coefficient in terms.items():
            denominator2 = exponent2 + 2 * denominator_power + 2
            value += (
                abs(lower.arb_from_fraction(coefficient))
                * count
                * arb(2) ** denominator_power
                * abs(p ** (p_power - denominator_power))
                * absolute_log_moment(split, denominator2, log_power)
            )
    return value.upper()


def weighted_upper_remainder_bound(
    base_level: dict[int, dict[tuple[int, int], Fraction]],
    log_level: LogLevel,
    p: arb,
    q_plus: arb,
    count: int,
    split: arb,
    boundary: arb,
    slabs: int,
    hardy_shift: acb,
) -> arb:
    require(slabs >= 32, "too few weighted upper remainder slabs")
    width = (boundary - split) / slabs
    value = arb(0)
    for index in range(slabs):
        left = split + width * index
        right = split + width * (index + 1)
        midpoint = (left + right) / 2
        y_ball = arb(midpoint, (right - left) / 2)
        d_lower = (p / right + right - q_plus + 1).lower()
        require(d_lower > 0, "weighted nonstationary endpoint gap lost")
        for power in set(base_level) | set(log_level):
            coefficient = weighted_coefficient(
                base_level.get(power), log_level.get(power), p, y_ball, hardy_shift
            )
            finite_bound = (arb(count) * d_lower ** (-power)).upper()
            infinite_bound = (
                d_lower ** (-power) + d_lower ** (1 - power) / (power - 1)
            ).upper()
            value += width * abs(coefficient).upper() * min(finite_bound, infinite_bound)
    return value.upper()


def weighted_lower_cell_certificate(
    t: arb,
    boundary: arb,
    q_plus: arb,
    count: int,
    order: int,
    split: arb,
    upper_slabs: int,
    hardy_shift: acb,
) -> dict[str, Any]:
    pi = arb.pi()
    p = t / (2 * pi)
    a = p / boundary + boundary - q_plus
    nearest_gap = a + 1
    require(count > 0 and count % 2 == 1, "odd finite roster required")
    require(boundary < p.sqrt(), "weighted lower branch crossed the central saddle")
    require(nearest_gap.lower() > 0, "weighted source phase became stationary")
    require(split > 0 and split < boundary, "invalid weighted small-y split")
    require((p / (2 * split) + split - q_plus + 1).lower() > 0, "small-y gap guard failed")

    base_table = lower.coefficient_table(order)
    log_table = logarithmic_coefficient_table(order)
    common_phase = acb(
        0,
        -t * boundary.log() - pi * boundary**2 + 2 * pi * q_plus * boundary,
    ).exp()
    partial = acb(0)
    endpoint_terms: list[dict[str, Any]] = []
    for n in range(order):
        endpoint_sum = acb(0)
        for power in set(base_table[n]) | set(log_table[n]):
            coefficient = weighted_coefficient(
                base_table[n].get(power),
                log_table[n].get(power),
                p,
                boundary,
                hardy_shift,
            )
            endpoint_sum += coefficient * lower.finite_alternating_shifted_sum(
                a, count, power + 1
            )
        factor = acb(0, 1) / (2 * pi) * (-acb(0, 1) / (2 * pi)) ** n
        term = common_phase * factor * endpoint_sum
        partial += term
        endpoint_terms.append({"order": n, "term_ball": complex_record(term)})

    base_low_raw = lower.small_y_remainder_bound(
        base_table[order], p, count, split
    )
    log_low_raw = logarithmic_small_y_bound(log_table[order], p, count, split)
    scale = (2 * pi) ** (-order)
    low = (
        scale * (abs(hardy_shift).upper() * base_low_raw + log_low_raw)
    ).upper()
    high = (
        scale
        * weighted_upper_remainder_bound(
            base_table[order],
            log_table[order],
            p,
            q_plus,
            count,
            split,
            boundary,
            upper_slabs,
            hardy_shift,
        )
    ).upper()
    remainder = (low + high).upper()
    enclosure = lower.add_complex_error(partial, remainder)
    return {
        "passed": True,
        "precision_bits": ctx.prec,
        "height_ball": t.str(58, more=True),
        "boundary_L": boundary.str(40, more=True),
        "q_plus": q_plus.str(40, more=True),
        "source_count": count,
        "expansion_order": order,
        "small_y_split": split.str(35, more=True),
        "upper_remainder_slabs": upper_slabs,
        "hardy_shift_ball": complex_record(hardy_shift),
        "nearest_phase_gap_ball": real_record(nearest_gap),
        "endpoint_terms": endpoint_terms,
        "endpoint_partial_sum_ball": complex_record(partial),
        "small_y_remainder_bound_ball": real_record(low),
        "upper_interval_remainder_bound_ball": real_record(high),
        "total_remainder_bound_ball": real_record(remainder),
        "weighted_lower_cell_ball": complex_record(enclosure),
        "operator_identity": (
            "(d_t+i*theta'-H'/H)F="
            "[i*(theta'-log(y))-H'/H]*F"
        ),
        "recurrence_identity": (
            "T^n[(C-i*log(y))*y^(-1/2)]="
            "C*T^n[y^(-1/2)]-i*T^n[log(y)y^(-1/2)]"
        ),
    }


def stable_H_box(t: arb) -> tuple[arb, dict[str, str]]:
    center = arb(HEIGHT)
    log_center = joined.direct_log_H(center)
    direct = joined.direct_log_H_derivative(t)
    duplicate = joined.duplication_log_H_derivative(t)
    require(direct.overlaps(duplicate), "H logarithmic-derivative formulas miss")
    derivative_bound = min(abs(direct).upper(), abs(duplicate).upper())
    distance = abs(t - center).upper()
    log_box = log_center + arb(0, distance * derivative_bound)
    return log_box.exp(), {
        "center_log_H_ball": log_center.str(58, more=True),
        "direct_logarithmic_derivative_ball": direct.str(58, more=True),
        "duplication_logarithmic_derivative_ball": duplicate.str(58, more=True),
        "maximum_transport_distance": distance.str(35, more=True),
    }


def altered_three_label_fixture() -> dict[str, Any]:
    """Changed finite-roster check against direct logarithmic-chart quadrature."""

    import mpmath as mp

    t_text = "100"
    boundary_text = "1.5"
    q_plus_text = "5.5"
    count = 3
    shift = acb(arb("0.3"), arb("0.7"))
    certificate = weighted_lower_cell_certificate(
        arb(t_text),
        arb(boundary_text),
        arb(q_plus_text),
        count,
        6,
        arb("0.5"),
        2048,
        shift,
    )

    def direct(dps: int, cuts: list[int]) -> mp.mpc:
        mp.mp.dps = dps
        t = mp.mpf(t_text)
        boundary = mp.mpf(boundary_text)
        q_plus = mp.mpf(q_plus_text)
        labels = [q_plus - r for r in range(1, count + 1)]
        c = mp.mpc("0.3", "0.7")

        def integrand(v: mp.mpf) -> mp.mpc:
            if not mp.isfinite(v):
                return mp.mpc(0)
            y = boundary * mp.exp(-v)
            roster = mp.fsum(
                mp.exp(-1j * mp.pi * y**2 + 2j * mp.pi * label * y)
                for label in labels
            )
            return (c - 1j * mp.log(y)) * y ** (mp.mpf("0.5") - 1j * t) * roster

        return mp.quad(integrand, [*cuts, mp.inf])

    first = direct(75, [0, 1, 2, 4, 8, 16, 32, 64])
    second = direct(95, [0, 1, 3, 6, 12, 24, 48, 96])
    enclosure = lower.complex_from_record(certificate["weighted_lower_cell_ball"])
    midpoint = complex(float(enclosure.real.mid()), float(enclosure.imag.mid()))
    allowed = float(arb(certificate["total_remainder_bound_ball"]["ball"]).upper())
    discrepancy = abs(complex(second) - midpoint)
    drift = abs(second - first)
    require(discrepancy < allowed, "altered weighted direct integral misses recurrence ball")
    require(discrepancy < 1e-7, "altered weighted midpoint disagrees with direct quadrature")
    require(drift < mp.mpf("1e-7"), "altered weighted direct quadrature is unstable")
    return {
        "height": t_text,
        "boundary": boundary_text,
        "q_plus": q_plus_text,
        "source_count": count,
        "hardy_shift_fixture": "0.3+0.7i",
        "direct_to_recurrence_discrepancy": mp.nstr(discrepancy, 20),
        "direct_precision_drift": mp.nstr(drift, 20),
        "allowed_remainder": mp.nstr(allowed, 20),
        "overlap": True,
    }


def build(radius_text: str, precision_bits: int, variant: str) -> dict[str, Any]:
    started = time.perf_counter()
    resource_mode = joined.set_low_priority()
    require(resource_mode == "below_normal_one_cpu", f"resource cap unavailable: {resource_mode}")
    ctx.prec = precision_bits
    ctx.threads = 1
    require(variant in ("production", "independent"), "unknown lower-cell K variant")
    radius = arb(radius_text)
    require(radius >= 0 and radius <= arb("0.001"), "radius must lie in [0,0.001]")
    t = arb(arb(HEIGHT), radius)
    event = load_json(EVENT_ATLAS)["certificate"]["primary_open_cell"]
    require(
        arb(event["lower_ball"]["ball"]).upper() < t.lower()
        and t.upper() < arb(event["upper_ball"]["ball"]).lower(),
        "height box left the fixed-roster event cell",
    )

    if variant == "production":
        order, split, slabs = 18, arb("300"), 4096
    else:
        order, split, slabs = 20, arb("280"), 6144
    theta, theta_prime = joined.theta_and_derivative(t)
    H_log_prime = joined.direct_log_H_derivative(t)
    require(
        H_log_prime.overlaps(joined.duplication_log_H_derivative(t)),
        "H logarithmic-derivative formulas miss",
    )
    hardy_shift = acb(-H_log_prime, theta_prime)
    boundary = lower.arb_from_fraction(lower.LOWER_BOUNDARY)
    q_plus = lower.arb_from_fraction(lower.Q_PLUS)

    direct = weighted_lower_cell_certificate(
        t,
        boundary,
        q_plus,
        lower.SOURCE_COUNT,
        order,
        split,
        slabs,
        hardy_shift,
    )
    time_derivative = weighted_lower_cell_certificate(
        t,
        boundary,
        q_plus,
        lower.SOURCE_COUNT,
        order,
        split,
        slabs,
        acb(0),
    )
    value_certificate = lower.lower_cell_certificate(
        t,
        boundary,
        q_plus,
        lower.SOURCE_COUNT,
        order,
        split,
        slabs,
        resource_mode,
    )
    direct_ball = lower.complex_from_record(direct["weighted_lower_cell_ball"])
    derivative_ball = lower.complex_from_record(
        time_derivative["weighted_lower_cell_ball"]
    )
    value_ball = lower.complex_from_record(value_certificate["lower_cell_integral_ball"])
    assembled_ball = derivative_ball + hardy_shift * value_ball
    require(direct_ball.overlaps(assembled_ball), "direct and assembled lower-cell K miss")
    H, H_transport = stable_H_box(t)
    contribution = joined.hardy_projection(theta, direct_ball) / H

    return {
        "kind": "rh_diagnostic_joined_lower_finite_cell_K_first_subcell_scout",
        "date": "2026-08-28",
        "status": "lower_finite_cell_Hardy_operator_diagnostic_not_complete_K_T",
        "passed": True,
        "resource_mode": resource_mode,
        "workers": 1,
        "precision_bits": precision_bits,
        "variant": variant,
        "center_height": str(HEIGHT),
        "radius": radius_text,
        "height_ball": t.str(58, more=True),
        "configuration": {
            "expansion_order": order,
            "small_y_split": split.str(30, more=True),
            "upper_remainder_slabs": slabs,
        },
        "exact_identities": {
            "lower_cell": "V_L=integral_0^L F(y,t)dy",
            "time_derivative": "d_t F=-i*log(y)*F",
            "Hardy_operator": (
                "K_L=(d_t+i*theta'-H'/H)V_L="
                "integral_0^L[i*(theta'-log(y))-H'/H]Fdy"
            ),
            "joined_next_coordinate": (
                "V_L+O_join=V_L-T_lower-T_upper+Gamma_defect; "
                "all three label ranges meet at L and must be rejoined before a final norm"
            ),
        },
        "phase_and_prefactor": {
            "theta": real_record(theta),
            "theta_prime": real_record(theta_prime),
            "H": real_record(H),
            "H_transport": H_transport,
            "H_log_prime": real_record(H_log_prime),
            "Hardy_shift": complex_record(hardy_shift),
        },
        "lower_cell_value": complex_record(value_ball),
        "lower_cell_time_derivative": complex_record(derivative_ball),
        "lower_cell_K_direct": complex_record(direct_ball),
        "lower_cell_K_assembled": complex_record(assembled_ball),
        "direct_and_assembled_overlap": True,
        "lower_cell_Q_derivative_contribution": real_record(contribution),
        "direct_certificate": direct,
        "time_derivative_certificate": time_derivative,
        "altered_three_label_fixture": altered_three_label_fixture(),
        "elapsed_seconds": time.perf_counter() - started,
        "diagnostic_boundary": (
            "This certifies no theorem by itself. The lower finite cell is only one ownership "
            "piece of K_T; its complete-lattice endpoint currents must be joined to the lower "
            "and upper ordinary complements, then to the already certified transition--arc pair, "
            "positive-real tail derivative, and tiny target correction."
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
    contribution = artifact["lower_cell_Q_derivative_contribution"]
    print(
        "completed joined lower finite-cell K diagnostic; "
        f"Hardy_contribution={contribution['ball']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
