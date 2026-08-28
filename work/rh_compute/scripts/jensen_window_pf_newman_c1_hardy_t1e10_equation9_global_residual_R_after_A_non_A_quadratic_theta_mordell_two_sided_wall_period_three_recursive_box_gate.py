#!/usr/bin/env python3
"""Certify both sides of the x=2/5 affine Mordell recursion wall."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_two_sided_wall_period_three_recursive_box_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
RIGHT_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_right_wall_period_three_recursive_box_gate.py"
)
RIGHT_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_right_wall_period_three_recursive_box_gate.json"
)

A = 159_577
L = 2_481_422
K = L - 1
FIRST_P = Fraction(31_916)
FIRST_N = 992_568
X0 = Fraction(2, 5)
S0 = Fraction(1, 3)
LEFT_BASE_P = Fraction(47_873, 3)
LEFT_BASE_A = Fraction(2, 3)
LEFT_BASE_TAU = Fraction(-1)
TWO_SIDED_WIDTH = Fraction(1, 10**23)
SCALE_POWERS = (20, 21, 22, 23, 24, 26, 28)


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


def parse_acb(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def conjugate_pair(value: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    """Conjugate u+v*omega using omega^2=-1-omega."""

    return value[0] - value[1], -value[1]


def left_period_three_power_pair(n: int, degree: int, right: Any, small: Any) -> tuple[Fraction, Fraction]:
    return conjugate_pair(right.period_three_power_pair(n, degree, small))


def left_period_three_weighted_pair(
    p: Fraction,
    n: int,
    degree: int,
    right: Any,
    small: Any,
) -> tuple[Fraction, Fraction]:
    return right.pair_add(
        right.pair_scale(left_period_three_power_pair(n, degree, right, small), p),
        left_period_three_power_pair(n, degree + 1, right, small),
    )


def normalized_left_boxes(x_box: Any, s_box: Any, box: Any) -> tuple[Any, Any]:
    a_box = box.RationalBox(
        Fraction(3, 2) + Fraction(A, 2) + s_box.lo - FIRST_P / x_box.lo,
        Fraction(3, 2) + Fraction(A, 2) + s_box.hi - FIRST_P / x_box.hi,
    )
    tau_box = box.RationalBox(
        Fraction(3, 2) - Fraction(1, 2) / x_box.lo,
        Fraction(3, 2) - Fraction(1, 2) / x_box.hi,
    )
    require(Fraction(-1, 2) < a_box.lo <= a_box.hi < Fraction(1, 2) * 2, "left a branch escaped")
    require(Fraction(0) < tau_box.lo <= tau_box.hi <= Fraction(1, 4), "left tau branch escaped")
    require(a_box.lo > 0 and a_box.hi < 1, "left nearest-integer branch crossed")
    return a_box, tau_box


def left_child_taylor_box(
    p_box: Any,
    a_box: Any,
    tau_box: Any,
    n: int,
    box: Any,
    right: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    require(p_box.lo > 0, "left child weight crossed zero")
    require(a_box.lo <= LEFT_BASE_A <= a_box.hi, "left base a missing")
    require(tau_box.lo <= LEFT_BASE_TAU <= tau_box.hi, "left base tau missing")
    f_pairs = [left_period_three_power_pair(n, degree, right, small) for degree in range(4)]
    m_pairs = [
        left_period_three_weighted_pair(LEFT_BASE_P, n, degree, right, small)
        for degree in range(3)
    ]
    f_balls = [right.pair_ball(value, prior) for value in f_pairs]
    m_balls = [right.pair_ball(value, prior) for value in m_pairs]

    rho = box.arb_box(p_box) - prior.arb_rational(LEFT_BASE_P)
    alpha = box.arb_box(a_box) - prior.arb_rational(LEFT_BASE_A)
    beta = box.arb_box(tau_box) - prior.arb_rational(LEFT_BASE_TAU)
    b0 = m_balls[0] + acb(rho) * f_balls[0]
    b1 = m_balls[1] + acb(rho) * f_balls[1]
    b2 = m_balls[2] + acb(rho) * f_balls[2]
    first_order = b0 + 2 * acb.pi() * acb(0, 1) * (acb(alpha) * b1 + acb(beta) * b2)

    delta_p = right.max_distance(p_box, LEFT_BASE_P)
    delta_a = right.max_distance(a_box, LEFT_BASE_A)
    delta_tau = right.max_distance(tau_box, LEFT_BASE_TAU)
    masses = {degree: small.weighted_mass(p_box.hi, n, degree) for degree in range(2, 5)}
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
    linear_radius = dp * f_balls[0].abs_upper() + 2 * arb.pi() * (
        da * (m_balls[1].abs_upper() + dp * f_balls[1].abs_upper())
        + dt * (m_balls[2].abs_upper() + dp * f_balls[2].abs_upper())
    )
    base_current, base_period = affine.weighted_period_ball(
        LEFT_BASE_P, LEFT_BASE_A, LEFT_BASE_TAU, n, prior
    )
    require(base_period == 3, "left base period drift")
    require(enclosure.contains(base_current), "left Taylor box misses its base current")
    return enclosure, {
        "n": n,
        "p_box": box.box_record(p_box),
        "a_box": box.box_record(a_box),
        "tau_box": box.box_record(tau_box),
        "base_p": str(LEFT_BASE_P),
        "base_a": str(LEFT_BASE_A),
        "base_tau": str(LEFT_BASE_TAU),
        "base_period": base_period,
        "delta_p_max": str(delta_p),
        "delta_a_max": str(delta_a),
        "delta_tau_max": str(delta_tau),
        "F_pairs": {f"F{degree}": right.pair_record(value) for degree, value in enumerate(f_pairs)},
        "M_pairs": {f"M{degree}": right.pair_record(value) for degree, value in enumerate(m_pairs)},
        "S2": str(masses[2]),
        "S3": str(masses[3]),
        "S4": str(masses[4]),
        "first_order_variation_radius_upper": box.arb_upper_text(linear_radius),
        "second_order_remainder_radius_upper": box.arb_upper_text(remainder_radius),
        "total_deviation_from_base_upper": box.arb_upper_text(linear_radius + remainder_radius),
        "first_order_enclosure": box.acb_record(first_order),
        "child_enclosure": box.acb_record(enclosure),
        "exact_base_current": box.acb_record(base_current),
        "contains_exact_base_current": True,
    }


def left_parent_branch_box(
    name: str,
    validity_domain: str,
    a_box: Any,
    tau_box: Any,
    child_m: int,
    settings: Any,
    box: Any,
    right: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    flint.ctx.prec = settings.precision_bits
    p_child, a_child, tau_child = right.child_parameter_boxes(a_box, tau_box, box)
    child, child_record = left_child_taylor_box(
        p_child, a_child, tau_child, child_m, box, right, small, affine, prior
    )

    phase_box = box.RationalBox(
        Fraction(1, 4) - a_box.hi * a_box.hi / (2 * tau_box.lo),
        Fraction(1, 4) - a_box.lo * a_box.lo / (2 * tau_box.hi),
    )
    tt = box.arb_box(tau_box)
    coefficient = box.cis_pi_interval(box.arb_box(phase_box)) / (2 * acb(tt) * acb(2 * tt).sqrt())
    main = coefficient * child

    lower_argument = box.RationalBox(
        a_box.lo - tau_box.hi + Fraction(1, 2),
        a_box.hi - tau_box.lo + Fraction(1, 2),
    )
    upper_argument = box.RationalBox(
        a_box.lo + (2 * FIRST_N + 1) * tau_box.lo - child_m - Fraction(1, 2),
        a_box.hi + (2 * FIRST_N + 1) * tau_box.hi - child_m - Fraction(1, 2),
    )
    mordell_tau = box.RationalBox(2 * tau_box.lo, 2 * tau_box.hi)
    lower_h, lower_hz, lower_record = box.mordell_target_parameter_box(
        lower_argument, mordell_tau, -1, settings, prior
    )
    upper_h, upper_hz, upper_record = box.mordell_target_parameter_box(
        upper_argument, mordell_tau, -1, settings, prior
    )
    lower_phase = box.RationalBox(
        -a_box.hi + tau_box.lo / 2,
        -a_box.lo + tau_box.hi / 2,
    )
    half_n = Fraction(FIRST_N) + Fraction(1, 2)
    upper_phase = box.RationalBox(
        2 * half_n * (a_box.lo + tau_box.lo * half_n),
        2 * half_n * (a_box.hi + tau_box.hi * half_n),
    )
    pi_i = acb.pi() * acb(0, 1)
    endpoint = -acb(0, 1) * (
        box.cis_pi_interval(box.arb_box(lower_phase))
        * (acb(prior.arb_rational(2 * FIRST_P - 1)) * lower_h + lower_hz / pi_i)
        + ((-1) ** child_m)
        * box.cis_pi_interval(box.arb_box(upper_phase))
        * (acb(prior.arb_rational(2 * FIRST_P + 2 * FIRST_N + 1)) * upper_h + upper_hz / pi_i)
    ) / 4
    complete = main + endpoint
    return complete, {
        "name": name,
        "validity_domain": validity_domain,
        "original_current_orientation": "direct_after_double_conjugation",
        "parent_p": str(FIRST_P),
        "parent_n": FIRST_N,
        "parent_a_box": box.box_record(a_box),
        "parent_tau_box": box.box_record(tau_box),
        "parent_nearest_integer_shift": 0,
        "child_m": child_m,
        "child_p_box": box.box_record(p_child),
        "child_a_box": box.box_record(a_child),
        "child_tau_box": box.box_record(tau_child),
        "child_Taylor": child_record,
        "main_phase_box": box.box_record(phase_box),
        "main_coefficient": box.acb_record(coefficient),
        "recursive_main": box.acb_record(main),
        "lower_argument_box": box.box_record(lower_argument),
        "upper_argument_box": box.box_record(upper_argument),
        "lower_Mordell": lower_record,
        "upper_Mordell": upper_record,
        "lower_phase_box": box.box_record(lower_phase),
        "upper_phase_box": box.box_record(upper_phase),
        "joined_endpoint": box.acb_record(endpoint),
        "complete_original_parent": box.acb_record(complete),
        "passed": True,
    }


def two_sided_source_box(
    width: Fraction,
    settings: Any,
    box: Any,
    right: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> dict[str, Any]:
    flint.ctx.prec = settings.precision_bits
    s_box = box.RationalBox(S0 - width, S0 + width)
    left_x = box.RationalBox(X0 - width, X0)
    wall_x = box.RationalBox(X0, X0)
    a_left, tau_left = normalized_left_boxes(left_x, s_box, box)
    a_wall, tau_wall = normalized_left_boxes(wall_x, s_box, box)

    lower_floor_argument = 2 * FIRST_N * tau_left.lo
    upper_floor_limit = 2 * FIRST_N * tau_left.hi
    require(496_283 < lower_floor_argument < 496_284, "left open floor lower edge drift")
    require(upper_floor_limit == 496_284, "left wall floor limit drift")
    require(2 * FIRST_N * tau_wall.lo == 496_284, "left wall floor drift")

    left_open, left_open_record = left_parent_branch_box(
        "open_left_branch",
        "x in [2/5-width,2/5), enclosed on the closed parameter hull",
        a_left,
        tau_left,
        496_283,
        settings,
        box,
        right,
        small,
        affine,
        prior,
    )
    left_wall, left_wall_record = left_parent_branch_box(
        "left_wall_orientation_witness",
        "x=2/5 and s in [1/3-width,1/3+width]",
        a_wall,
        tau_wall,
        496_284,
        settings,
        box,
        right,
        small,
        affine,
        prior,
    )
    require(left_open.overlaps(left_wall), "left open closure misses its wall representation")

    right_row = right.right_wall_source_box(width, settings, box, small, affine, prior)
    right_open_normalized = parse_acb(right_row["interior_branch"]["complete_normalized_parent"])
    right_wall_normalized = parse_acb(right_row["wall_branch"]["complete_normalized_parent"])
    right_open_original = right_open_normalized.conjugate()
    right_wall_original = right_wall_normalized.conjugate()
    require(left_wall.overlaps(right_wall_original), "left/right wall orientations do not overlap")
    require(right_open_original.overlaps(right_wall_original), "right open closure misses the wall")

    transformed_hull = left_open.union(left_wall).union(right_wall_original).union(right_open_original)
    center_transformed, transformed_period = affine.weighted_period_ball(
        FIRST_P, Fraction(1, 3), Fraction(1, 4), FIRST_N, prior
    )
    require(left_wall.contains(center_transformed), "left wall misses exact center transformed current")
    require(right_wall_original.contains(center_transformed), "right wall orientation misses exact center transformed current")
    require(transformed_hull.contains(center_transformed), "two-sided transformed hull misses center")

    x_box = box.RationalBox(X0 - width, X0 + width)
    first_phase = box.RationalBox(
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.lo) - FIRST_P * FIRST_P / x_box.lo,
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.hi) - FIRST_P * FIRST_P / x_box.hi,
    )
    xx = box.arb_box(x_box)
    first_multiplier = acb(2 / (xx * xx.sqrt())) * box.cis_pi_interval(box.arb_box(first_phase))
    recursive_first_main = first_multiplier * transformed_hull
    first_endpoint, endpoint_record = box.physical_endpoint_parameter_box(
        "two_sided_wall_first_joined_endpoint", x_box, s_box, settings, prior
    )
    complete_source = recursive_first_main + first_endpoint
    center_source, source_period = prior.source_current_ball(X0, S0)
    require(complete_source.contains(center_source), "two-sided source box misses exact center source")
    return {
        "half_width": str(width),
        "x_box": box.box_record(x_box),
        "s_box": box.box_record(s_box),
        "left_normalization": "a_L=3/2+A/2+s-31916/x, tau_L=3/2-1/(2x)",
        "orientation_identity": "conjugate(conjugate(W_left))=W_left after the parity half-shift",
        "left_open_child_floor": 496_283,
        "left_open_floor_lower_argument": str(lower_floor_argument),
        "left_open_floor_upper_limit_excluded": str(upper_floor_limit),
        "wall_child_floor": 496_284,
        "left_open_branch": left_open_record,
        "left_wall_branch": left_wall_record,
        "right_dependency_projection": {
            "open_branch_child_floor": right_row["open_branch_child_floor"],
            "wall_branch_child_floor": right_row["wall_branch_child_floor"],
            "right_open_original_parent": box.acb_record(right_open_original),
            "right_wall_original_parent": box.acb_record(right_wall_original),
        },
        "left_open_overlaps_left_wall": True,
        "left_wall_overlaps_right_wall_orientation": True,
        "right_open_overlaps_right_wall": True,
        "transformed_two_sided_hull": box.acb_record(transformed_hull),
        "exact_center_transformed_current": box.acb_record(center_transformed),
        "exact_center_transformed_period": transformed_period,
        "transformed_hull_contains_exact_center": True,
        "first_phase_box": box.box_record(first_phase),
        "first_multiplier": box.acb_record(first_multiplier),
        "recursive_first_main": box.acb_record(recursive_first_main),
        "first_joined_endpoint": endpoint_record,
        "complete_source_box": box.acb_record(complete_source),
        "exact_center_source": box.acb_record(center_source),
        "exact_center_source_period": source_period,
        "complete_contains_exact_center_source": True,
        "complete_radius_absolute_upper": box.arb_upper_text(complete_source.rad()),
        "passed": True,
    }


def left_scale_audit(
    settings: Any,
    box: Any,
    right: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> list[dict[str, Any]]:
    flint.ctx.prec = settings.precision_bits
    rows = []
    for power in SCALE_POWERS:
        width = Fraction(1, 10**power)
        x_box = box.RationalBox(X0 - width, X0)
        s_box = box.RationalBox(S0 - width, S0 + width)
        a_box, tau_box = normalized_left_boxes(x_box, s_box, box)
        p_child, a_child, tau_child = right.child_parameter_boxes(a_box, tau_box, box)
        _, record = left_child_taylor_box(
            p_child, a_child, tau_child, 496_283, box, right, small, affine, prior
        )
        parent_tau = box.arb_box(tau_box)
        second_coefficient_upper = (1 / (2 * parent_tau * (2 * parent_tau).sqrt())).abs_upper()
        physical_x = box.arb_box(x_box)
        first_coefficient_upper = (2 / (physical_x * physical_x.sqrt())).abs_upper()
        child_deviation = arb(record["total_deviation_from_base_upper"])
        rows.append(
            {
                "half_width_decimal": f"1e-{power}",
                "half_width": str(width),
                "delta_p_max": record["delta_p_max"],
                "delta_a_max": record["delta_a_max"],
                "delta_tau_max": record["delta_tau_max"],
                "child_first_order_variation_upper": record["first_order_variation_radius_upper"],
                "child_second_order_remainder_upper": record["second_order_remainder_radius_upper"],
                "child_total_deviation_upper": record["total_deviation_from_base_upper"],
                "two_level_child_variation_upper": box.arb_upper_text(
                    first_coefficient_upper * second_coefficient_upper * child_deviation
                ),
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    row = artifact["certified_two_sided_wall_box"]
    left = row["left_open_branch"]["child_Taylor"]
    audit_lines = [
        "| half-width | max delta p | max delta a | max delta tau | left child deviation | two-level propagation |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for item in artifact["left_scale_audit"]:
        audit_lines.append(
            f"| `{item['half_width_decimal']}` | `{item['delta_p_max']}` | `{item['delta_a_max']}` | "
            f"`{item['delta_tau_max']}` | `{item['child_total_deviation_upper']}` | "
            f"`{item['two_level_child_variation_upper']}` |"
        )
    return f"""# Two-sided affine Mordell wall and period-three recursive box

Date: 2026-08-24

Status: certified two-sided local interval lemma; recursive tiling remains open

## Exact left normalization

For `x<2/5`, the right coordinate has `tau_R>1/4`. Apply the exact parity
half-shift and conjugation:

```text
(TW1) a_L=1/2-a_R=3/2+A/2+s-31916/x,
      tau_L=1/2-tau_R=3/2-1/(2x).
```

The new conjugation cancels the first conjugation, so the normalized left
parent is already in the original transformed-current orientation. At the
wall `(a_L,tau_L)=(1/3,1/4)`.

## Left child and orientation witness

The open left child has index `496283` and tends to

```text
(p_0,a_0,tau_0)=(47873/3,2/3,-1).
```

It uses the conjugate period-three root `omega^2`. Exact residue pairs and the
same cancellation-preserving Taylor form give

```text
(TW2) W=B_0(rho)+2*pi*i[alpha B_1(rho)+beta B_2(rho)]+R_2,

(TW3) |R_2|<=2*pi^2[delta_a^2 S_2
                     +2 delta_a delta_tau S_3
                     +delta_tau^2 S_4].
```

At width `10^-23`, the left child deviation is below
`{left['total_deviation_from_base_upper']}`.

A separate `496284` left-wall branch is completed with its own Mordell main
and joined endpoints. It overlaps the conjugate of the independently
certified right-wall branch, proving the orientation match before hulling.

## Two-sided source box

The completed transformed branches satisfy

```text
(TW4) T_two=hull(T_left_open,T_left_wall,
                 conjugate(T_right_wall),conjugate(T_right_open)),

(TW5) SourceBox=P_1*T_two+Endpoint_1.
```

The certified physical box is

```text
x in [2/5-10^-23,2/5+10^-23],
s in [1/3-10^-23,1/3+10^-23].
```

Its complete source radius is below `{row['complete_radius_absolute_upper']}`
and it contains the exact centre source current.

## Scale audit

{chr(10).join(audit_lines)}

## Validation

- Production uses 320 bits, Mordell cutoff 10, and one numerical worker.
- The independent checker derives the double-conjugation map and both floor
  branches, verifies the `omega^2` residue pairs against direct witnesses, and
  replays the four completed parent branches and source assembly at 384 bits.
- Every `pi` comes from `e(u)=exp(2*pi*i*u)` or the inherited exact Mordell
  transform; no circle or polygon constant is inserted.

## Proof boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(RIGHT_BUILDER.is_file() and RIGHT_RESULT.is_file() and CHECKER.is_file(), "right-wall dependency missing")
    right_artifact = json.loads(RIGHT_RESULT.read_text(encoding="utf-8"))
    require(right_artifact.get("passed") is True, "right-wall dependency not passed")
    right = load_module("mordell_two_sided_right_dependency", RIGHT_BUILDER)
    small = load_module("mordell_two_sided_small_dependency", right.SMALL_BUILDER)
    affine = load_module("mordell_two_sided_affine_dependency", right.AFFINE_BUILDER)
    box = load_module("mordell_two_sided_box_dependency", right.BOX_BUILDER)
    prior = load_module("mordell_two_sided_point_dependency", right.PRIOR_BUILDER)
    settings = box.PRODUCTION

    two_sided = two_sided_source_box(
        TWO_SIDED_WIDTH, settings, box, right, small, affine, prior
    )
    audits = left_scale_audit(settings, box, right, small, affine, prior)
    selected = next(row for row in audits if row["half_width_decimal"] == "1e-23")
    require(two_sided["left_open_child_floor"] == 496_283, "left open floor drift")
    require(two_sided["wall_child_floor"] == 496_284, "left wall floor drift")
    require(bool(arb(selected["two_level_child_variation_upper"]) < arb("0.002")), "left child scale target failed")
    require(bool(arb(two_sided["complete_radius_absolute_upper"]) < arb("0.06")), "two-sided source radius target failed")

    proof_boundary = (
        "One two-sided width-1e-23 physical box at the x=2/5 normalization wall, the exact left "
        "parity/double-conjugation map, both 496283 open child branches, two 496284 wall orientation "
        "representations, and their endpoint-complete source hull only. No recursive tiling beyond this "
        "local box, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, "
        "Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "two_sided_width_1e_minus_23_affine_Mordell_wall_period_three_recursive_and_complete_source_box_certified_recursive_tiling_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "L": L,
            "K": K,
            "workers": 1,
            "two_sided_half_width": str(TWO_SIDED_WIDTH),
        },
        "production_settings": settings.__dict__,
        "exact_theorem": {
            "left_normalization": "a_L=3/2+A/2+s-31916/x, tau_L=3/2-1/(2x)",
            "left_orientation": "direct_after_double_conjugation",
            "open_child_indices": {"left": 496_283, "right": 496_283},
            "wall_child_index": 496_284,
            "left_period_three_base": "(p,a,tau)=(47873/3,2/3,-1)",
        },
        "certified_two_sided_wall_box": two_sided,
        "left_scale_audit": audits,
        "summary": {
            "certified_two_sided_half_width": str(TWO_SIDED_WIDTH),
            "left_open_child_index": two_sided["left_open_child_floor"],
            "right_open_child_index": two_sided["right_dependency_projection"]["open_branch_child_floor"],
            "wall_child_index": two_sided["wall_child_floor"],
            "left_wall_overlaps_right_wall_orientation": True,
            "selected_left_two_level_child_variation_upper": selected["two_level_child_variation_upper"],
            "complete_radius_absolute_upper": two_sided["complete_radius_absolute_upper"],
            "complete_contains_exact_center_source": True,
        },
        "decision": {
            "left_parity_half_shift_derived": True,
            "double_conjugation_orientation_proved": True,
            "left_period_three_child_Taylor_boxes_built": True,
            "left_endpoint_complete_recursive_box_built": True,
            "left_and_right_wall_orientations_joined": True,
            "two_sided_interior_source_box_certified": True,
            "recursive_tiling_beyond_local_box_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Use the now-closed local normalization wall as one node of a recursive interval evaluator. "
            "Derive the next branch walls for the two period-three children, measure contraction and "
            "dependency inflation over adjacent boxes, and decide whether adaptive recursive tiling can "
            "reach a useful physical non-A quadrature scale."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": right_artifact["primary_source"],
        "dependencies": {
            "right_wall_result": {"path": relative(RIGHT_RESULT), "sha256": file_hash(RIGHT_RESULT)},
            "right_wall_builder": {"path": relative(RIGHT_BUILDER), "sha256": file_hash(RIGHT_BUILDER)},
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
    print("certified the two-sided Mordell wall and period-three recursive source box", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
