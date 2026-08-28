#!/usr/bin/env python3
"""Rigorous Arb enclosure for the production upper quarter-disk arc."""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import sys


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import acb, arb  # noqa: E402


HEIGHT = 10_000_000_000
Y_TEXT = "39936.5"
Q_MINUS_TEXT = "79788.5"
COUNT = 2_481_423


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            process.cpu_affinity([process.cpu_affinity()[0]])
            return "below_normal_one_cpu"
        process.nice(10)
        process.cpu_affinity([process.cpu_affinity()[0]])
        return "nice_10_one_cpu"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def configuration(variant: str) -> dict[str, object]:
    if variant == "production":
        return {
            "precision_bits": 288,
            "degree": 24,
            "delta_min": "0.0000000001",
            "delta_cut": "0.003",
            "initial_panels": 128,
            "tolerance": "1e-26",
            "main_points": (
                "0.0000000001",
                "0.000001",
                "0.00001",
                "0.00005",
                "0.0001",
                "0.00015",
                "0.0002",
                "0.0003",
                "0.0005",
                "0.0008",
                "0.0012",
                "0.0018",
                "0.0024",
                "0.003",
            ),
        }
    if variant == "independent":
        return {
            "precision_bits": 336,
            "degree": 26,
            "delta_min": "0.0000000002",
            "delta_cut": "0.0032",
            "initial_panels": 192,
            "tolerance": "1e-27",
            "main_points": (
                "0.0000000002",
                "0.0000015",
                "0.000012",
                "0.00004",
                "0.00008",
                "0.00013",
                "0.00018",
                "0.00026",
                "0.0004",
                "0.00065",
                "0.001",
                "0.00145",
                "0.002",
                "0.0026",
                "0.0032",
            ),
        }
    raise ValueError(f"unknown variant: {variant}")


def polynomial(coefficients: list[acb], z: acb) -> acb:
    value = coefficients[-1]
    for coefficient in reversed(coefficients[1:-1]):
        value = value * z + coefficient
    return value * z


def ball_payload(value: acb) -> dict[str, object]:
    return {
        "real_ball": str(value.real),
        "imag_ball": str(value.imag),
        "real_mid_float": float(value.real.mid()),
        "imag_mid_float": float(value.imag.mid()),
        "real_radius_float": float(value.real.rad()),
        "imag_radius_float": float(value.imag.rad()),
    }


def real_payload(value: arb) -> dict[str, object]:
    return {
        "ball": str(value),
        "mid_float": float(value.mid()),
        "radius_float": float(value.rad()),
        "lower_float": float(value.lower()),
        "upper_float": float(value.upper()),
    }


def run(variant: str, reference: Path | None) -> dict[str, object]:
    cfg = configuration(variant)
    flint.ctx.prec = int(cfg["precision_bits"])
    priority = set_low_priority()
    require(priority == "below_normal_one_cpu", f"resource cap unavailable: {priority}")

    pi = arb.pi()
    imaginary = acb(0, 1)
    t = arb(HEIGHT)
    y = arb(Y_TEXT)
    q_minus = arb(Q_MINUS_TEXT)
    radius = 2 * pi * y
    count = COUNT
    sqrt_y = y.sqrt()
    degree = int(cfg["degree"])
    delta_min = arb(str(cfg["delta_min"]))
    delta_cut = arb(str(cfg["delta_cut"]))

    # The small endpoint coefficient is exposed before interval evaluation.
    slope = t + 2 * pi * y * y - q_minus * radius
    require(slope > 0, "upper endpoint action must initially increase")

    exponent_coefficients = [acb(0) for _ in range(degree + 1)]
    w_coefficients = [acb(0) for _ in range(degree + 1)]
    exponent_coefficients[1] = slope + imaginary / 2
    for order in range(2, degree + 1):
        factorial = arb(math.factorial(order))
        exponent_coefficients[order] = (
            -imaginary * pi * y * y * (2 * imaginary) ** order
            + imaginary * q_minus * radius * imaginary**order
        ) / factorial
    for order in range(1, degree + 1):
        w_coefficients[order] = (
            imaginary * radius * imaginary**order / arb(math.factorial(order))
        )

    def components(z: acb) -> tuple[acb, acb]:
        return polynomial(exponent_coefficients, z), polynomial(w_coefficients, z)

    def full_integrand(z: acb, analytic: bool = True) -> acb:
        exponent, w_exponent = components(z)
        w = w_exponent.exp()
        return sqrt_y * exponent.exp() * (1 + (count * w_exponent).exp()) / (1 + w)

    def reduced_integrand(z: acb, analytic: bool = True) -> acb:
        exponent, w_exponent = components(z)
        w = w_exponent.exp()
        return sqrt_y * exponent.exp() / (1 + w)

    tolerance = arb(str(cfg["tolerance"]))
    initial = acb(0)
    initial_panels = int(cfg["initial_panels"])
    for index in range(initial_panels):
        left = delta_min * index / initial_panels
        right = delta_min * (index + 1) / initial_panels
        initial += acb.integral(
            full_integrand,
            left,
            right,
            abs_tol=tolerance,
            rel_tol=tolerance,
            eval_limit=150_000,
            depth_limit=40,
        )

    main = acb(0)
    points = [arb(text) for text in cfg["main_points"]]
    require(
        (points[0] - delta_min).contains(0) and (points[-1] - delta_cut).contains(0),
        "panel endpoints drifted",
    )
    for left, right in zip(points, points[1:]):
        main += acb.integral(
            reduced_integrand,
            left,
            right,
            abs_tol=tolerance,
            rel_tol=tolerance,
            eval_limit=300_000,
            depth_limit=45,
        )

    # A'(delta) is a quadratic in c=cos(delta).  Its upper root is the
    # unique production action maximum; the lower root is negative.
    constant = t - 2 * pi * y * y
    discriminant = (q_minus * radius) ** 2 - 16 * pi * y * y * constant
    root_low = (q_minus * radius - discriminant.sqrt()) / (8 * pi * y * y)
    root_high = (q_minus * radius + discriminant.sqrt()) / (8 * pi * y * y)
    require(root_low < 0, "lower action-derivative root is not negative")
    require(0 < delta_cut.cos() < root_high < 1, "action monotonic tail guard failed")
    delta_star = root_high.acos()

    def action(delta: arb) -> arb:
        return t * delta + pi * y * y * (2 * delta).sin() - q_minus * radius * delta.sin()

    action_max = action(delta_star)
    action_cut = action(delta_cut)
    require(action_max < arb("0.01981"), "production action maximum exceeded its guard")
    require(action_cut < arb("-80"), "compact arc cutoff is not in the decaying tail")

    phase_deviation = radius * (1 - delta_cut.cos())
    require(phase_deviation < pi / 2, "core geometric denominator left the positive-real sector")

    next_order = degree + 1
    factorial = arb(math.factorial(next_order))
    exponent_remainder = (
        pi * y * y * (2 * delta_cut).exp() * (2 * delta_cut) ** next_order
        + q_minus * radius * delta_cut.exp() * delta_cut**next_order
    ) / factorial
    w_remainder = radius * delta_cut.exp() * delta_cut**next_order / factorial
    polynomial_relative_error = (
        2 * exponent_remainder * exponent_remainder.exp()
        + 8
        * (count + 2)
        * w_remainder
        * ((count + 2) * w_remainder).exp()
    )
    polynomial_integral_error = (
        2
        * delta_cut
        * sqrt_y
        * action_max.exp()
        * polynomial_relative_error
    )

    omitted_numerator_error = (
        (delta_cut - delta_min)
        * sqrt_y
        * action_max.exp()
        * (-count * radius * delta_min.sin()).exp()
    )

    radial_decay = (-radius * delta_cut.sin()).exp()
    quotient_tail = (1 + radial_decay**count) / (1 - radial_decay)
    compact_arc_tail = (
        sqrt_y
        * (pi / 2 - delta_cut)
        * action_cut.exp()
        * quotient_tail
    )
    total_added_error = polynomial_integral_error + omitted_numerator_error + compact_arc_tail
    require(total_added_error < arb("1e-20"), "analytic arc error budget is too large")

    rotated_arc = initial + main + acb(
        arb(0, total_added_error),
        arb(0, total_added_error),
    )
    rotated_arc_absolute = abs(rotated_arc)
    require(rotated_arc_absolute < arb("0.057665"), "upper arc missed its certified scale guard")

    endpoint_phase = pi - t * y.log() - pi * y * y
    actual_arc = (imaginary * endpoint_phase).exp() * rotated_arc

    positive_tail_log10 = (
        pi * t / 2
        - q_minus * radius
        - (2 * pi).log() / 2
        - radius.log() / 2
        - q_minus.log()
        - (1 - (-radius).exp()).log()
    ) / arb(10).log()
    require(positive_tail_log10 < arb("-100000000"), "positive-real tail is not negligible")

    overlaps_reference: bool | None = None
    if reference is not None:
        payload = json.loads(reference.read_text(encoding="utf-8"))
        saved = payload["arb_certificate"]["rotated_arc"]
        reference_ball = acb(arb(saved["real_ball"]), arb(saved["imag_ball"]))
        overlaps_reference = bool(rotated_arc.overlaps(reference_ball))
        require(overlaps_reference, "independent Arb enclosure does not overlap the saved enclosure")

    return {
        "passed": True,
        "variant": variant,
        "resource_mode": priority,
        "precision_bits": int(cfg["precision_bits"]),
        "degree": degree,
        "delta_min": str(delta_min),
        "delta_cut": str(delta_cut),
        "initial_panels": initial_panels,
        "main_panel_count": len(points) - 1,
        "stable_endpoint_exponent": (
            "X(delta)=(t+2*pi*Y^2-q_-R+i/2)delta"
            "-i*pi*Y^2[expm1(2i*delta)-2i*delta]"
            "+i*q_-R[expm1(i*delta)-i*delta]"
        ),
        "stable_geometric_coordinate": (
            "w=exp(i*R*expm1(i*delta)); Delta(u)/Delta(-iR)="
            "exp(-q_-(u+iR))*(1+w^N)/(1+w)"
        ),
        "slope": real_payload(slope),
        "action_stationary_delta": real_payload(delta_star),
        "action_maximum": real_payload(action_max),
        "action_at_cut": real_payload(action_cut),
        "phase_deviation_at_cut": real_payload(phase_deviation),
        "initial_full_numerator_integral": ball_payload(initial),
        "reduced_main_integral": ball_payload(main),
        "error_budgets": {
            "exponent_taylor_remainder": real_payload(exponent_remainder),
            "w_taylor_remainder": real_payload(w_remainder),
            "polynomial_integral_error": real_payload(polynomial_integral_error),
            "omitted_geometric_numerator_error": real_payload(omitted_numerator_error),
            "compact_arc_tail": real_payload(compact_arc_tail),
            "total_added_error": real_payload(total_added_error),
            "positive_real_tail_log10_upper": real_payload(positive_tail_log10),
        },
        "rotated_arc": ball_payload(rotated_arc),
        "rotated_arc_absolute": real_payload(rotated_arc_absolute),
        "actual_arc": ball_payload(actual_arc),
        "endpoint_phase": "phi_U=pi-t*log(Y)-pi*Y^2 for U=39936+1/2",
        "overlaps_reference": overlaps_reference,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", choices=("production", "independent"), default="production")
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.variant, args.reference), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
