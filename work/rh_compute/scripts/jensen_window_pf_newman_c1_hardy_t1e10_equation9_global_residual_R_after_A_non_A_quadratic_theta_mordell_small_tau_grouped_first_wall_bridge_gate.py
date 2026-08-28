#!/usr/bin/env python3
"""Certify the grouped m3=0|1 first-wall bridge by a common finite current."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_small_tau_grouped_first_wall_bridge_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
ADAPTIVE_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_adaptive_width_mirrored_tile_crossover_gate.py"
)
ADAPTIVE_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_adaptive_width_mirrored_tile_crossover_gate.json"
)

A = 159_577
FIRST_P = Fraction(31_916)
FIRST_N = 992_568
CHILD_N = 496_283
S0 = Fraction(1, 3)
TAU0 = Fraction(1, 2 * CHILD_N)
RIGHT_X = Fraction(2 * CHILD_N + 1, 5 * CHILD_N + 2)
LEFT_X = Fraction(2 * CHILD_N + 1, 5 * CHILD_N + 3)
RIGHT_P = Fraction(95_747, 6)
LEFT_P = Fraction(47_873, 3)
RIGHT_A = Fraction(-896_819, 6 * CHILD_N)
LEFT_A = Fraction(544_156, 3 * CHILD_N)
WALL_HALF_WIDTH = Fraction(1, 10**24)
S_HALF_WIDTH = Fraction(1, 10**24)
MOMENT_BLOCK_SIZE = 16
SOURCE_RADIUS_TARGET = "0.06"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def set_low_priority() -> str:
    try:
        import psutil

        process = psutil.Process()
        if os.name == "nt":
            process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            return "below_normal"
        process.nice(10)
        return "nice_10"
    except Exception as exc:  # pragma: no cover
        return f"unavailable:{type(exc).__name__}"


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def max_distance(value: Any, center: Fraction) -> Fraction:
    return max(abs(value.lo - center), abs(value.hi - center))


def blockwise_phase_moments(
    base_a: Fraction,
    base_tau: Fraction,
    n: int,
    block_size: int,
    prior: Any,
) -> tuple[list[acb], dict[str, Any]]:
    """Compute F_j=sum k^j e(a*k+tau*k^2), j=0..3, with exact block resets."""

    require(block_size > 0 and n >= 0, "invalid blockwise phase-moment request")
    started = time.time()
    moments = [acb(0) for _ in range(4)]
    step = prior.acb_cis_pi(4 * base_tau)
    for start in range(0, n + 1, block_size):
        stop = min(n + 1, start + block_size)
        phase = base_a * start + base_tau * start * start
        ratio_phase = base_a + base_tau * (2 * start + 1)
        term = prior.acb_cis_pi(2 * phase)
        ratio = prior.acb_cis_pi(2 * ratio_phase)
        for k in range(start, stop):
            power = 1
            for degree in range(4):
                moments[degree] += power * term
                power *= k
            term *= ratio
            ratio *= step
    require(all(value.is_finite() for value in moments), "nonfinite blockwise phase moment")
    return moments, {
        "method": "exact_phase_reset_plus_quadratic_ratio_recurrence",
        "block_size": block_size,
        "block_count": (n + block_size) // block_size,
        "phase_convention": "e(u)=exp(2*pi*i*u)",
        "elapsed_seconds": round(time.time() - started, 6),
    }


def wall_taylor_current(
    side: str,
    p_box: Any,
    a_box: Any,
    tau_box: Any,
    box: Any,
    small: Any,
    prior: Any,
    block_size: int,
) -> tuple[acb, dict[str, Any]]:
    require(side in ("right", "left"), "unknown wall side")
    base_p = RIGHT_P if side == "right" else LEFT_P
    base_a = RIGHT_A if side == "right" else LEFT_A
    require(p_box.lo <= base_p <= p_box.hi, f"{side} wall p anchor escaped")
    require(a_box.lo <= base_a <= a_box.hi, f"{side} wall a anchor escaped")
    require(tau_box.lo <= TAU0 <= tau_box.hi, f"{side} wall tau anchor escaped")

    f_balls, method = blockwise_phase_moments(base_a, TAU0, CHILD_N, block_size, prior)
    m_balls = [
        acb(prior.arb_rational(base_p)) * f_balls[degree] + f_balls[degree + 1]
        for degree in range(3)
    ]
    rho = box.arb_box(p_box) - prior.arb_rational(base_p)
    alpha = box.arb_box(a_box) - prior.arb_rational(base_a)
    beta = box.arb_box(tau_box) - prior.arb_rational(TAU0)
    b0 = m_balls[0] + acb(rho) * f_balls[0]
    b1 = m_balls[1] + acb(rho) * f_balls[1]
    b2 = m_balls[2] + acb(rho) * f_balls[2]
    first_order = b0 + 2 * acb.pi() * acb(0, 1) * (
        acb(alpha) * b1 + acb(beta) * b2
    )

    delta_p = max_distance(p_box, base_p)
    delta_a = max_distance(a_box, base_a)
    delta_tau = max_distance(tau_box, TAU0)
    masses = {degree: small.weighted_mass(p_box.hi, CHILD_N, degree) for degree in range(2, 5)}
    remainder_mass = (
        delta_a * delta_a * masses[2]
        + 2 * delta_a * delta_tau * masses[3]
        + delta_tau * delta_tau * masses[4]
    )
    remainder_radius = 2 * arb.pi() ** 2 * prior.arb_rational(remainder_mass)
    enclosure = first_order + box.symmetric_complex_error(remainder_radius)

    dp = prior.arb_rational(delta_p)
    da = prior.arb_rational(delta_a)
    dt = prior.arb_rational(delta_tau)
    variation = dp * f_balls[0].abs_upper() + 2 * arb.pi() * (
        da * (m_balls[1].abs_upper() + dp * f_balls[1].abs_upper())
        + dt * (m_balls[2].abs_upper() + dp * f_balls[2].abs_upper())
    )
    anchor = m_balls[0]
    require(enclosure.overlaps(anchor), f"{side} wall Taylor box misses its anchor")
    return enclosure, {
        "side": side,
        "n": CHILD_N,
        "base_p": str(base_p),
        "base_a": str(base_a),
        "base_tau": str(TAU0),
        "minimal_phase_period": 3 * CHILD_N if side == "right" else 6 * CHILD_N,
        "parameter_boxes": {
            "p": box.box_record(p_box),
            "a": box.box_record(a_box),
            "tau": box.box_record(tau_box),
        },
        "delta_p_max": str(delta_p),
        "delta_a_max": str(delta_a),
        "delta_tau_max": str(delta_tau),
        "F_balls": {f"F{degree}": box.acb_record(value) for degree, value in enumerate(f_balls)},
        "M_balls": {f"M{degree}": box.acb_record(value) for degree, value in enumerate(m_balls)},
        "S2": str(masses[2]),
        "S3": str(masses[3]),
        "S4": str(masses[4]),
        "first_order_variation_radius_upper": box.arb_upper_text(variation),
        "second_order_remainder_radius_upper": box.arb_upper_text(remainder_radius),
        "total_deviation_from_anchor_upper": box.arb_upper_text(variation + remainder_radius),
        "anchor_current": box.acb_record(anchor),
        "first_order_enclosure": box.acb_record(first_order),
        "complete_wall_taylor_enclosure": box.acb_record(enclosure),
        "moment_evaluator": method,
        "passed": True,
    }


def grouped_floor_jump_audit(side: str) -> dict[str, Any]:
    """Prove that the m=1 child jump and upper-endpoint jump cancel at 2*n*tau=1."""

    require(side in ("right", "left"), "unknown grouped-jump side")
    p = RIGHT_P if side == "right" else LEFT_P
    z = RIGHT_A if side == "right" else LEFT_A
    t = TAU0
    n = CHILD_N
    child_p = 2 * t * p - z
    require(child_p == (Fraction(1, 3) if side == "right" else Fraction(-1, 3)), "child p drift")
    w = z / (2 * t)
    sigma = -Fraction(1, 4 * t)
    u_one = z + (2 * n + 1) * t - 1 - Fraction(1, 2)
    u_zero = z + (2 * n + 1) * t - Fraction(1, 2)
    require(u_zero == u_one + 1, "upper Mordell arguments are not unit neighbours")

    endpoint_weight = Fraction(2) * p + 2 * n + 1
    recurrence_derivative_factor = -(u_one + Fraction(1, 2)) / t
    simplified_weight = endpoint_weight + recurrence_derivative_factor
    require(simplified_weight == (child_p + 1) / t, "endpoint jump weight did not simplify")

    upper_phase = 2 * (Fraction(n) + Fraction(1, 2)) * (
        z + t * (Fraction(n) + Fraction(1, 2))
    )
    endpoint_jump_phase = Fraction(1, 4) + upper_phase - (u_one + Fraction(1, 2)) ** 2 / (2 * t)
    child_jump_phase = Fraction(5, 4) - z * z / (2 * t) + 2 * (w + sigma)
    phase_difference = endpoint_jump_phase - child_jump_phase
    require(phase_difference == 2 * n, "grouped jump phase difference drift")
    require((phase_difference / 2).denominator == 1, "grouped jump phases are not equal modulo two")
    return {
        "side": side,
        "wall_condition": "2*n*tau=1",
        "m_below": 0,
        "m_at_and_above": 1,
        "child_p_at_wall": str(child_p),
        "child_w_at_wall": str(w),
        "child_sigma_at_wall": str(sigma),
        "upper_argument_m0": str(u_zero),
        "upper_argument_m1": str(u_one),
        "unit_recurrence_argument_identity": "u_plus(0)=u_plus(1)+1",
        "endpoint_weight_before_derivative": str(endpoint_weight),
        "endpoint_weight_after_unit_recurrence": str(simplified_weight),
        "required_child_jump_weight": str((child_p + 1) / t),
        "endpoint_minus_child_phase_difference": str(phase_difference),
        "phase_difference_is_even_integer": True,
        "exact_conclusion": "Delta_upper_endpoint=-Delta_third_main",
        "passed": True,
    }


def parent_complete_box(
    side: str,
    x_box: Any,
    s_box: Any,
    settings: Any,
    box: Any,
    small: Any,
    prior: Any,
    rational: Any,
    adaptive: Any,
    block_size: int,
) -> tuple[acb, dict[str, Any]]:
    require(side in ("right", "left"), "unknown parent side")
    if side == "right":
        a_parent, tau_parent = rational.exact_right_parent_boxes(x_box, s_box, box)
        p_child, a_child, tau_child, cancellation = rational.exact_right_normalized_child_boxes(
            x_box, s_box, box
        )
    else:
        a_parent, tau_parent = adaptive.exact_left_parent_boxes(x_box, s_box, box)
        p_child, a_child, tau_child, cancellation = adaptive.exact_left_normalized_child_boxes(
            x_box, s_box, box, rational
        )
    second_floor = box.RationalBox(2 * FIRST_N * tau_parent.lo, 2 * FIRST_N * tau_parent.hi)
    require(CHILD_N < second_floor.lo <= second_floor.hi < CHILD_N + 1, "second floor wall crossed")
    child_normalized, child_record = wall_taylor_current(
        side, p_child, a_child, tau_child, box, small, prior, block_size
    )
    child_original = child_normalized.conjugate()

    phase_box = box.RationalBox(
        Fraction(1, 4) - a_parent.hi * a_parent.hi / (2 * tau_parent.lo),
        Fraction(1, 4) - a_parent.lo * a_parent.lo / (2 * tau_parent.hi),
    )
    tt = box.arb_box(tau_parent)
    coefficient = box.cis_pi_interval(box.arb_box(phase_box)) / (2 * acb(tt) * acb(2 * tt).sqrt())
    main = coefficient * child_original

    if side == "right":
        c_lo = Fraction(A, 2) + s_box.lo + 1
        c_hi = Fraction(A, 2) + s_box.hi + 1
        lower_argument = box.RationalBox(
            (FIRST_P - Fraction(1, 2)) / x_box.hi - Fraction(A, 2) - s_box.hi + Fraction(1, 2),
            (FIRST_P - Fraction(1, 2)) / x_box.lo - Fraction(A, 2) - s_box.lo + Fraction(1, 2),
        )
        d = FIRST_P + FIRST_N + Fraction(1, 2)
        upper_argument = box.RationalBox(
            d / x_box.hi - c_hi - 2 * FIRST_N - CHILD_N - Fraction(3, 2),
            d / x_box.lo - c_lo - 2 * FIRST_N - CHILD_N - Fraction(3, 2),
        )
        lower_phase = box.RationalBox(
            c_lo - Fraction(1, 2) - (FIRST_P - Fraction(1, 4)) / x_box.lo,
            c_hi - Fraction(1, 2) - (FIRST_P - Fraction(1, 4)) / x_box.hi,
        )
        half_n = Fraction(FIRST_N) + Fraction(1, 2)
        upper_inner = box.RationalBox(
            (FIRST_P + half_n / 2) / x_box.hi - c_hi - half_n,
            (FIRST_P + half_n / 2) / x_box.lo - c_lo - half_n,
        )
    else:
        lower_argument = box.RationalBox(
            Fraction(A, 2) + s_box.lo + Fraction(1, 2) - (FIRST_P - Fraction(1, 2)) / x_box.lo,
            Fraction(A, 2) + s_box.hi + Fraction(1, 2) - (FIRST_P - Fraction(1, 2)) / x_box.hi,
        )
        d = FIRST_P + FIRST_N + Fraction(1, 2)
        upper_argument = box.RationalBox(
            Fraction(A, 2) + s_box.lo + 3 * FIRST_N + Fraction(5, 2) - CHILD_N - d / x_box.lo,
            Fraction(A, 2) + s_box.hi + 3 * FIRST_N + Fraction(5, 2) - CHILD_N - d / x_box.hi,
        )
        lower_phase = box.RationalBox(
            -Fraction(A, 2) - s_box.hi - Fraction(3, 4) + (FIRST_P - Fraction(1, 4)) / x_box.hi,
            -Fraction(A, 2) - s_box.lo - Fraction(3, 4) + (FIRST_P - Fraction(1, 4)) / x_box.lo,
        )
        half_n = Fraction(FIRST_N) + Fraction(1, 2)
        upper_inner = box.RationalBox(
            Fraction(3, 2) + Fraction(A, 2) + s_box.lo + 3 * half_n / 2
            - (FIRST_P + half_n / 2) / x_box.lo,
            Fraction(3, 2) + Fraction(A, 2) + s_box.hi + 3 * half_n / 2
            - (FIRST_P + half_n / 2) / x_box.hi,
        )
    mordell_tau = box.RationalBox(2 * tau_parent.lo, 2 * tau_parent.hi)
    lower_h, lower_hz, lower_record = box.mordell_target_parameter_box(
        lower_argument, mordell_tau, -1, settings, prior
    )
    upper_h, upper_hz, upper_record = box.mordell_target_parameter_box(
        upper_argument, mordell_tau, -1, settings, prior
    )
    upper_phase = box.RationalBox(2 * half_n * upper_inner.lo, 2 * half_n * upper_inner.hi)
    pi_i = acb.pi() * acb(0, 1)
    endpoint = -acb(0, 1) * (
        box.cis_pi_interval(box.arb_box(lower_phase))
        * (acb(prior.arb_rational(2 * FIRST_P - 1)) * lower_h + lower_hz / pi_i)
        + ((-1) ** CHILD_N)
        * box.cis_pi_interval(box.arb_box(upper_phase))
        * (acb(prior.arb_rational(2 * FIRST_P + 2 * FIRST_N + 1)) * upper_h + upper_hz / pi_i)
    ) / 4
    complete = main + endpoint
    third_floor = box.RationalBox(2 * CHILD_N * tau_child.lo, 2 * CHILD_N * tau_child.hi)
    require(third_floor.lo < 1 < third_floor.hi, f"{side} box does not cross the first third floor wall")
    return complete, {
        "side": side,
        "parent_a_box": box.box_record(a_parent),
        "parent_tau_box": box.box_record(tau_parent),
        "second_floor_argument_box": box.box_record(second_floor),
        "normalized_child_p_box": box.box_record(p_child),
        "normalized_child_a_box": box.box_record(a_child),
        "normalized_child_tau_box": box.box_record(tau_child),
        "cancellation_audit": cancellation,
        "common_wall_taylor_child": child_record,
        "original_child_current": box.acb_record(child_original),
        "second_main_phase_box": box.box_record(phase_box),
        "second_main_coefficient": box.acb_record(coefficient),
        "second_recursive_main": box.acb_record(main),
        "second_lower_argument_box": box.box_record(lower_argument),
        "second_upper_argument_box": box.box_record(upper_argument),
        "second_lower_Mordell": lower_record,
        "second_upper_Mordell": upper_record,
        "second_lower_phase_box": box.box_record(lower_phase),
        "second_upper_phase_box": box.box_record(upper_phase),
        "second_joined_endpoint": box.acb_record(endpoint),
        "complete_parent": box.acb_record(complete),
        "third_floor_argument_box": box.box_record(third_floor),
        "third_floor_branches_crossed": [0, 1],
        "passed": True,
    }


def source_wall_box(
    side: str,
    settings: Any,
    box: Any,
    small: Any,
    prior: Any,
    rational: Any,
    adaptive: Any,
    block_size: int,
) -> dict[str, Any]:
    x0 = RIGHT_X if side == "right" else LEFT_X
    x_box = box.RationalBox(x0 - WALL_HALF_WIDTH, x0 + WALL_HALF_WIDTH)
    s_box = box.RationalBox(S0 - S_HALF_WIDTH, S0 + S_HALF_WIDTH)
    parent, parent_record = parent_complete_box(
        side, x_box, s_box, settings, box, small, prior, rational, adaptive, block_size
    )
    original_transformed = parent.conjugate() if side == "right" else parent
    first_phase = box.RationalBox(
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.lo) - FIRST_P * FIRST_P / x_box.lo,
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.hi) - FIRST_P * FIRST_P / x_box.hi,
    )
    xx = box.arb_box(x_box)
    first_multiplier = acb(2 / (xx * xx.sqrt())) * box.cis_pi_interval(box.arb_box(first_phase))
    first_main = first_multiplier * original_transformed
    first_endpoint, endpoint_record = box.physical_endpoint_parameter_box(
        f"{side}_first_positive_third_floor_wall_first_joined_endpoint",
        x_box,
        s_box,
        settings,
        prior,
    )
    complete = first_main + first_endpoint
    radius = complete.rad()
    require(complete.is_finite(), f"{side} complete source box is nonfinite")
    require(
        radius < arb(SOURCE_RADIUS_TARGET),
        f"{side} complete source radius {box.arb_upper_text(radius)} exceeds {SOURCE_RADIUS_TARGET}",
    )
    return {
        "side": side,
        "wall_x": str(x0),
        "x_box": box.box_record(x_box),
        "s_box": box.box_record(s_box),
        "parent": parent_record,
        "original_transformed_parent": box.acb_record(original_transformed),
        "first_phase_box": box.box_record(first_phase),
        "first_multiplier": box.acb_record(first_multiplier),
        "first_recursive_main": box.acb_record(first_main),
        "first_joined_endpoint": endpoint_record,
        "complete_source_box": box.acb_record(complete),
        "complete_source_radius_upper": box.arb_upper_text(radius),
        "source_radius_target": SOURCE_RADIUS_TARGET,
        "source_radius_target_met": True,
        "passed": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    right = artifact["right_first_wall_source_box"]
    left = artifact["left_first_wall_source_box"]
    r_child = right["parent"]["common_wall_taylor_child"]
    l_child = left["parent"]["common_wall_taylor_child"]
    return f"""# Grouped small-t first-floor-wall bridge

Date: 2026-08-25

Status: certified local `m3=0|1` wall crossing; finite recursive cover remains open

## Exact grouped jump

For the affine current `W_n(p;a,t)`, the third Mordell floor changes from
`m=0` to `m=1` at `2*n*t=1`, with `n={CHILD_N}` and

```text
t*=1/{2 * CHILD_N}.
```

The new child term and the change of the upper joined endpoint must be kept
together.  The two upper Mordell arguments satisfy `u_+(0)=u_+(1)+1`.
Applying the exact unit recurrence simultaneously to `h` and `h_z` reduces
the endpoint coefficient to `(p'+1)/t`, exactly the coefficient of the new
child term.  Their phase difference is

```text
2*n={2 * CHILD_N},
```

an even integer in the `exp(i*pi*phase)` convention.  Hence on both physical
sides

```text
Delta upper endpoint = - Delta third main.
```

At the two wall anchors the child affine weights are respectively `1/3` and
`-1/3`.  This proves that the floor switch is removable only after the main
and endpoint currents are grouped; it does not bound either piece separately.

## Common finite-current enclosure

Instead of evaluating the ill-conditioned small-`t` Mordell pieces, both
floor branches are represented by the analytic finite current

```text
W_n(p;a,t)=sum_(k=0)^n (p+k) exp(2*pi*i*(a*k+t*k^2)).
```

Four base moments are computed rigorously at each exact wall anchor.  The
production evaluator resets the exact quadratic phase every
`{MOMENT_BLOCK_SIZE}` terms, so rectangular complex-ball wrapping cannot
accumulate across the full sum.  For `rho=p-p0`, `alpha=a-a0`, and
`beta=t-t0`, the common enclosure is

```text
B0 + 2*pi*i*(alpha*B1+beta*B2) + R2,

|R2| <= 2*pi^2*(delta_a^2*S2
                 +2*delta_a*delta_t*S3
                 +delta_t^2*S4).
```

This is one enclosure of the full current across `m=0|1`, rather than a hull
of separately bounded Mordell components.

## Certified source boxes

Both physical boxes use x half-width `{WALL_HALF_WIDTH}` and s half-width
`{S_HALF_WIDTH}` around their exact first-wall anchors.  All earlier Mordell
mains and joined endpoints are retained through the physical source.

| side | child deviation | source radius | target |
|---|---:|---:|---:|
| right | `{r_child['total_deviation_from_anchor_upper']}` | `{right['complete_source_radius_upper']}` | `<{SOURCE_RADIUS_TARGET}` |
| left | `{l_child['total_deviation_from_anchor_upper']}` | `{left['complete_source_radius_upper']}` | `<{SOURCE_RADIUS_TARGET}` |

The independent checker uses 384 bits, block size 13, `Y=11`, and tolerance
`1e-65`; it recomputes every base moment with different reset boundaries,
rederives the exact jump algebra, and replays both complete source boxes.

## Boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(ADAPTIVE_BUILDER.is_file() and ADAPTIVE_RESULT.is_file(), "missing adaptive dependency")
    require(CHECKER.is_file(), "missing independent checker")
    previous = json.loads(ADAPTIVE_RESULT.read_text(encoding="utf-8"))
    require(previous.get("passed") is True, "adaptive dependency not passed")
    adaptive = load_module("grouped_wall_adaptive_dependency", ADAPTIVE_BUILDER)
    rational = load_module("grouped_wall_rational_dependency", adaptive.RATIONAL_BUILDER)
    two_sided = load_module("grouped_wall_two_sided_dependency", rational.TWO_SIDED_BUILDER)
    right_dependency = load_module("grouped_wall_right_dependency", two_sided.RIGHT_BUILDER)
    small = load_module("grouped_wall_small_dependency", right_dependency.SMALL_BUILDER)
    box = load_module("grouped_wall_box_dependency", right_dependency.BOX_BUILDER)
    prior = load_module("grouped_wall_point_dependency", right_dependency.PRIOR_BUILDER)
    settings = box.PRODUCTION
    flint.ctx.prec = settings.precision_bits

    jump = {side: grouped_floor_jump_audit(side) for side in ("right", "left")}
    right = source_wall_box("right", settings, box, small, prior, rational, adaptive, MOMENT_BLOCK_SIZE)
    left = source_wall_box("left", settings, box, small, prior, rational, adaptive, MOMENT_BLOCK_SIZE)

    proof_boundary = (
        "The exact m3=0|1 third-main/upper-endpoint jump cancellation and one endpoint-complete "
        "physical source box of x half-width 1e-24 and s half-width 1e-24 around each first positive "
        "third-floor wall only. The common finite-current Taylor enclosure avoids separately estimating "
        "small-t Mordell pieces. No cover between the normalization wall and these first floor walls, "
        "continuous adaptive-width theorem, subsequent floor wall, physical non-A quadrature, non-A bound, "
        "joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level "
        "conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "grouped_small_tau_m3_zero_one_first_wall_jump_and_two_physical_source_boxes_certified_finite_cover_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "workers": 1,
            "wall_half_width": str(WALL_HALF_WIDTH),
            "s_half_width": str(S_HALF_WIDTH),
            "source_radius_target": SOURCE_RADIUS_TARGET,
            "production_moment_block_size": MOMENT_BLOCK_SIZE,
            "daytime_cpu_baseline_samples_percent": [25.29, 14.20, 15.34, 18.07, 19.34, 17.72],
            "daytime_cpu_baseline_average_percent": 18.33,
        },
        "production_settings": settings.__dict__,
        "exact_grouped_floor_jump": jump,
        "right_first_wall_source_box": right,
        "left_first_wall_source_box": left,
        "summary": {
            "right_wall_x": str(RIGHT_X),
            "left_wall_x": str(LEFT_X),
            "wall_tau": str(TAU0),
            "right_child_p_after_floor_jump": jump["right"]["child_p_at_wall"],
            "left_child_p_after_floor_jump": jump["left"]["child_p_at_wall"],
            "right_source_radius_upper": right["complete_source_radius_upper"],
            "left_source_radius_upper": left["complete_source_radius_upper"],
        },
        "decision": {
            "m3_zero_one_main_endpoint_jump_grouped_exactly": True,
            "small_tau_Mordell_components_bounded_separately": False,
            "common_finite_current_crosses_first_floor_wall": True,
            "right_first_floor_wall_source_box_certified": True,
            "left_first_floor_wall_source_box_certified": True,
            "finite_recursive_cover_built": False,
            "physical_quadrature_completed": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Bridge the open interval between the local normalization-wall terminal and the first-floor-wall "
            "Taylor boxes. Derive blockwise van der Corput or Abel bounds with explicit parameter derivatives, "
            "or a hierarchy of exact phase-reset Taylor centres, and prove a finite cumulative source-radius "
            "budget before any physical quadrature."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": previous["primary_source"],
        "dependencies": {
            "adaptive_width_result": {"path": relative(ADAPTIVE_RESULT), "sha256": file_hash(ADAPTIVE_RESULT)},
            "adaptive_width_builder": {"path": relative(ADAPTIVE_BUILDER), "sha256": file_hash(ADAPTIVE_BUILDER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.time() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified grouped m3=0|1 first-wall bridge and two source boxes", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
