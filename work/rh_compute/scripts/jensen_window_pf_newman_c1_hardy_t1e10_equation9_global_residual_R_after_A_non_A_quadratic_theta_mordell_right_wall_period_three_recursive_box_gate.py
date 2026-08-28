#!/usr/bin/env python3
"""Certify the right side of the x=2/5 affine Mordell recursion wall."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
import math
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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_right_wall_period_three_recursive_box_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
SMALL_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_alternating_taylor_box_gate.py"
)
SMALL_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_alternating_taylor_box_gate.json"
)
AFFINE_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_affine_weighted_second_step_gate.py"
)
AFFINE_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_affine_weighted_second_step_gate.json"
)
BOX_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate.py"
)
BOX_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_parameter_box_phase_scale_gate.json"
)
PRIOR_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_interval_one_step_current_gate.py"
)
PRIOR_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_interval_one_step_current_gate.json"
)

A = 159_577
L = 2_481_422
K = L - 1
FIRST_P = Fraction(31_916)
FIRST_N = 992_568
X0 = Fraction(2, 5)
S0 = Fraction(1, 3)
BASE_CHILD_P = Fraction(95_747, 6)
BASE_CHILD_A = Fraction(1, 3)
BASE_CHILD_TAU = Fraction(-1)
RIGHT_WIDTH = Fraction(1, 10**23)
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


def pair_add(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return left[0] + right[0], left[1] + right[1]


def pair_scale(value: tuple[Fraction, Fraction], scalar: Fraction) -> tuple[Fraction, Fraction]:
    return scalar * value[0], scalar * value[1]


def pair_record(value: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"one_coefficient": str(value[0]), "omega_coefficient": str(value[1])}


def residue_power_sum(n: int, residue: int, modulus: int, degree: int, small: Any) -> int:
    """Return sum k^degree for k=residue mod modulus in 0..n."""

    require(modulus >= 1 and 0 <= residue < modulus, "invalid residue class")
    require(0 <= degree <= 5, "unsupported residue degree")
    if n < residue:
        return 0
    count = (n - residue) // modulus + 1
    return sum(
        math.comb(degree, j)
        * residue ** (degree - j)
        * modulus**j
        * small.power_sum(count, j)
        for j in range(degree + 1)
    )


def period_three_power_pair(n: int, degree: int, small: Any) -> tuple[Fraction, Fraction]:
    """Represent sum_(k=0)^n k^degree exp(2*pi*i*k/3) as u+v*omega."""

    coefficients = [residue_power_sum(n, residue, 3, degree, small) for residue in range(3)]
    return Fraction(coefficients[0] - coefficients[2]), Fraction(coefficients[1] - coefficients[2])


def period_three_weighted_pair(
    p: Fraction,
    n: int,
    degree: int,
    small: Any,
) -> tuple[Fraction, Fraction]:
    return pair_add(
        pair_scale(period_three_power_pair(n, degree, small), p),
        period_three_power_pair(n, degree + 1, small),
    )


def pair_ball(value: tuple[Fraction, Fraction], prior: Any) -> acb:
    omega = prior.acb_cis_pi(Fraction(2, 3))
    return acb(prior.arb_rational(value[0])) + acb(prior.arb_rational(value[1])) * omega


def max_distance(value: Any, center: Fraction) -> Fraction:
    return max(abs(value.lo - center), abs(value.hi - center))


def normalized_right_boxes(x_box: Any, s_box: Any, box: Any) -> tuple[Any, Any]:
    a_box = box.RationalBox(
        FIRST_P / x_box.hi - Fraction(A, 2) - s_box.hi - 1,
        FIRST_P / x_box.lo - Fraction(A, 2) - s_box.lo - 1,
    )
    tau_box = box.RationalBox(
        Fraction(1, 2) / x_box.hi - 1,
        Fraction(1, 2) / x_box.lo - 1,
    )
    require(Fraction(-1, 2) < a_box.lo <= a_box.hi < Fraction(1, 2), "right a branch crossed an integer wall")
    require(Fraction(0) < tau_box.lo <= tau_box.hi <= Fraction(1, 4), "right tau branch escaped")
    return a_box, tau_box


def child_parameter_boxes(a_box: Any, tau_box: Any, box: Any) -> tuple[Any, Any, Any]:
    p_box = box.RationalBox(
        2 * tau_box.lo * FIRST_P - a_box.hi,
        2 * tau_box.hi * FIRST_P - a_box.lo,
    )
    a_child = box.RationalBox(
        a_box.lo / (2 * tau_box.hi),
        a_box.hi / (2 * tau_box.lo),
    )
    tau_child = box.RationalBox(
        -Fraction(1, 4) / tau_box.lo,
        -Fraction(1, 4) / tau_box.hi,
    )
    return p_box, a_child, tau_child


def periodic_child_taylor_box(
    p_box: Any,
    a_box: Any,
    tau_box: Any,
    n: int,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    """Enclose a child near the exact period-three phase (1/3,-1)."""

    require(p_box.lo > 0, "period-three child weight crossed zero")
    require(a_box.lo <= BASE_CHILD_A <= a_box.hi, "base child a missing")
    require(tau_box.lo <= BASE_CHILD_TAU <= tau_box.hi, "base child tau missing")
    f_pairs = [period_three_power_pair(n, degree, small) for degree in range(4)]
    m_pairs = [period_three_weighted_pair(BASE_CHILD_P, n, degree, small) for degree in range(3)]
    f_balls = [pair_ball(value, prior) for value in f_pairs]
    m_balls = [pair_ball(value, prior) for value in m_pairs]

    rho = box.arb_box(p_box) - prior.arb_rational(BASE_CHILD_P)
    alpha = box.arb_box(a_box) - prior.arb_rational(BASE_CHILD_A)
    beta = box.arb_box(tau_box) - prior.arb_rational(BASE_CHILD_TAU)
    b0 = m_balls[0] + acb(rho) * f_balls[0]
    b1 = m_balls[1] + acb(rho) * f_balls[1]
    b2 = m_balls[2] + acb(rho) * f_balls[2]
    first_order = b0 + 2 * acb.pi() * acb(0, 1) * (acb(alpha) * b1 + acb(beta) * b2)

    delta_p = max_distance(p_box, BASE_CHILD_P)
    delta_a = max_distance(a_box, BASE_CHILD_A)
    delta_tau = max_distance(tau_box, BASE_CHILD_TAU)
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
        BASE_CHILD_P, BASE_CHILD_A, BASE_CHILD_TAU, n, prior
    )
    require(base_period == 3, "base child period drift")
    require(enclosure.contains(base_current), "period-three Taylor box misses its base current")
    return enclosure, {
        "n": n,
        "p_box": box.box_record(p_box),
        "a_box": box.box_record(a_box),
        "tau_box": box.box_record(tau_box),
        "base_p": str(BASE_CHILD_P),
        "base_a": str(BASE_CHILD_A),
        "base_tau": str(BASE_CHILD_TAU),
        "base_period": base_period,
        "delta_p_max": str(delta_p),
        "delta_a_max": str(delta_a),
        "delta_tau_max": str(delta_tau),
        "F_pairs": {f"F{degree}": pair_record(value) for degree, value in enumerate(f_pairs)},
        "M_pairs": {f"M{degree}": pair_record(value) for degree, value in enumerate(m_pairs)},
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


def affine_parent_branch_box(
    name: str,
    validity_domain: str,
    a_box: Any,
    tau_box: Any,
    child_m: int,
    settings: Any,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> tuple[acb, dict[str, Any]]:
    """Apply the endpoint-complete affine Mordell step on one fixed-m branch."""

    flint.ctx.prec = settings.precision_bits
    p_child, a_child, tau_child = child_parameter_boxes(a_box, tau_box, box)
    child, child_record = periodic_child_taylor_box(
        p_child, a_child, tau_child, child_m, box, small, affine, prior
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
        "complete_normalized_parent": box.acb_record(complete),
        "passed": True,
    }


def right_wall_source_box(
    width: Fraction,
    settings: Any,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
) -> dict[str, Any]:
    flint.ctx.prec = settings.precision_bits
    x_box = box.RationalBox(X0, X0 + width)
    s_box = box.RationalBox(S0 - width, S0 + width)
    a_right, tau_right = normalized_right_boxes(x_box, s_box, box)
    wall_x = box.RationalBox(X0, X0)
    a_wall, tau_wall = normalized_right_boxes(wall_x, s_box, box)

    open_child_m = 496_283
    wall_child_m = 496_284
    lower_floor_argument = 2 * FIRST_N * tau_right.lo
    upper_floor_limit = 2 * FIRST_N * tau_right.hi
    require(open_child_m < lower_floor_argument < open_child_m + 1, "open child floor lower edge drift")
    require(upper_floor_limit == open_child_m + 1, "open child floor wall limit drift")
    require(2 * FIRST_N * tau_wall.lo == wall_child_m, "wall child floor drift")

    interior, interior_record = affine_parent_branch_box(
        "open_right_branch",
        "x in (2/5,2/5+width], enclosed on the closed parameter hull",
        a_right,
        tau_right,
        open_child_m,
        settings,
        box,
        small,
        affine,
        prior,
    )
    wall, wall_record = affine_parent_branch_box(
        "wall_line_branch",
        "x=2/5 and s in [1/3-width,1/3+width]",
        a_wall,
        tau_wall,
        wall_child_m,
        settings,
        box,
        small,
        affine,
        prior,
    )
    require(interior.overlaps(wall), "open-right analytic limit does not overlap the wall branch")
    normalized_hull = interior.union(wall)
    original_transformed = normalized_hull.conjugate()

    first_phase = box.RationalBox(
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.lo) - FIRST_P * FIRST_P / x_box.lo,
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.hi) - FIRST_P * FIRST_P / x_box.hi,
    )
    xx = box.arb_box(x_box)
    first_multiplier = acb(2 / (xx * xx.sqrt())) * box.cis_pi_interval(box.arb_box(first_phase))
    recursive_first_main = first_multiplier * original_transformed
    first_endpoint, endpoint_record = box.physical_endpoint_parameter_box(
        "right_wall_first_joined_endpoint", x_box, s_box, settings, prior
    )
    complete_source = recursive_first_main + first_endpoint

    center_parent, center_parent_period = affine.weighted_period_ball(
        FIRST_P, Fraction(1, 6), Fraction(1, 4), FIRST_N, prior
    )
    require(wall.contains(center_parent), "wall branch misses exact center normalized current")
    require(interior.contains(center_parent), "open-right closure misses exact center limit")
    center_source, source_period = prior.source_current_ball(X0, S0)
    require(complete_source.contains(center_source), "right-wall source box misses exact center source")
    return {
        "half_width": str(width),
        "x_box": box.box_record(x_box),
        "s_box": box.box_record(s_box),
        "first_integer_shift_r": int(FIRST_P),
        "first_transformed_n": FIRST_N,
        "normalization_identity": "W(w,sigma)=conjugate(W(p/x-(A+2s)/2-1,1/(2x)-1))",
        "open_branch_child_floor": open_child_m,
        "open_branch_floor_lower_argument": str(lower_floor_argument),
        "open_branch_floor_upper_limit_excluded": str(upper_floor_limit),
        "wall_branch_child_floor": wall_child_m,
        "interior_branch": interior_record,
        "wall_branch": wall_record,
        "branch_complete_balls_overlap": True,
        "normalized_recursive_hull": box.acb_record(normalized_hull),
        "normalized_hull_contains_exact_center_parent": True,
        "exact_center_parent_period": center_parent_period,
        "exact_center_normalized_parent": box.acb_record(center_parent),
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


def child_scale_audit(settings: Any, box: Any, small: Any, affine: Any, prior: Any) -> list[dict[str, Any]]:
    flint.ctx.prec = settings.precision_bits
    rows = []
    for power in SCALE_POWERS:
        width = Fraction(1, 10**power)
        x_box = box.RationalBox(X0, X0 + width)
        s_box = box.RationalBox(S0 - width, S0 + width)
        a_box, tau_box = normalized_right_boxes(x_box, s_box, box)
        p_child, a_child, tau_child = child_parameter_boxes(a_box, tau_box, box)
        _, record = periodic_child_taylor_box(
            p_child, a_child, tau_child, 496_283, box, small, affine, prior
        )
        parent_tau = box.arb_box(tau_box)
        second_coefficient_upper = (
            1 / (2 * parent_tau * (2 * parent_tau).sqrt())
        ).abs_upper()
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
    row = artifact["certified_right_wall_box"]
    interior = row["interior_branch"]
    audit_lines = [
        "| half-width | max delta p | max delta a | max delta tau | child deviation | two-level propagation |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for item in artifact["child_scale_audit"]:
        audit_lines.append(
            f"| `{item['half_width_decimal']}` | `{item['delta_p_max']}` | `{item['delta_a_max']}` | "
            f"`{item['delta_tau_max']}` | `{item['child_total_deviation_upper']}` | "
            f"`{item['two_level_child_variation_upper']}` |"
        )
    return f"""# Right-side affine Mordell wall and period-three recursive box

Date: 2026-08-24

Status: certified right-side local interval lemma; the left side remains open

## Status

Certified at `t=10^10` on

```text
x in [2/5,2/5+10^-23],
s in [1/3-10^-23,1/3+10^-23].
```

This is one side of the interior normalization wall. It is not a complete
two-sided neighbourhood, a physical quadrature, a non-A bound, or an RH-level
result.

## Exact wall partition

After the first conjugation and integer shifts,

```text
(RW1) a=31916/x-(A+2s)/2-1,
      tau=1/(2x)-1,
      W_original=conjugate(W_992568(31916;a,tau)).
```

At `x=2/5`, `tau=1/4` and
`floor(2*992568*tau)=496284`. For every `x>2/5` in the box,
`tau<1/4` and the floor is `496283`. The closed source box is therefore the
union of the wall line and the open-right branch, not one silently frozen
floor branch.

## Period-three child enclosure

Both branches have child parameters near

```text
(p_0,a_0,tau_0)=(95747/6,1/3,-1).
```

Since `exp(2*pi*i*(k/3-k^2))=omega^k`, where
`omega=exp(2*pi*i/3)`, every base moment is represented exactly as
`u+v*omega` by splitting `k` modulo three. If `rho=p-p_0`,
`alpha=a-a_0`, and `beta=tau+1`, then

```text
(RW2) W=B_0(rho)+2*pi*i[alpha B_1(rho)+beta B_2(rho)]+R_2,
      B_j(rho)=M_j+rho F_j.
```

Only the quadratic phase remainder is bounded absolutely:

```text
(RW3) |R_2|<=2*pi^2[delta_a^2 S_2
                     +2 delta_a delta_tau S_3
                     +delta_tau^2 S_4].
```

For the open-right branch, the total period-three child deviation is below
`{interior['child_Taylor']['total_deviation_from_base_upper']}`.

## Endpoint-complete recursion

Each fixed-floor branch uses the affine Mordell identity with its own child
length, main coefficient, and both joined `h,h_z` endpoint currents. The
open branch is bounded on its closed analytic parameter hull but is asserted
only for `x>2/5`; the wall line supplies the excluded endpoint. Their complete
normalized-parent balls overlap and are joined by interval hull:

```text
(RW4) T_right = hull(T_m=496283, T_m=496284),
      SourceBox=P_1*conjugate(T_right)+Endpoint_1.
```

The complete physical source radius is below
`{row['complete_radius_absolute_upper']}` and contains the exact source current
at `(x,s)=(2/5,1/3)`.

## Scale audit

{chr(10).join(audit_lines)}

## Validation

- Production uses 320 bits, Mordell cutoff 10, and one numerical worker.
- The checker derives the branch floors independently, checks every residue
  power sum against direct witnesses, verifies both base period-three currents,
  and replays both endpoint-complete branches at 384 bits and cutoff 11.
- Every `pi` comes from `e(u)=exp(2*pi*i*u)` or the inherited exact Mordell
  transform; no circle or polygon constant is inserted.

## Proof boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    for path in (
        SMALL_BUILDER,
        SMALL_RESULT,
        AFFINE_BUILDER,
        AFFINE_RESULT,
        BOX_BUILDER,
        BOX_RESULT,
        PRIOR_BUILDER,
        PRIOR_RESULT,
        CHECKER,
    ):
        require(path.is_file(), f"missing dependency or checker: {path}")
    dependency_artifacts = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in (SMALL_RESULT, AFFINE_RESULT, BOX_RESULT, PRIOR_RESULT)
    ]
    require(all(item.get("passed") is True for item in dependency_artifacts), "Mordell dependency not passed")
    small = load_module("mordell_right_wall_small_dependency", SMALL_BUILDER)
    affine = load_module("mordell_right_wall_affine_dependency", AFFINE_BUILDER)
    box = load_module("mordell_right_wall_box_dependency", BOX_BUILDER)
    prior = load_module("mordell_right_wall_point_dependency", PRIOR_BUILDER)
    settings = box.PRODUCTION

    right = right_wall_source_box(RIGHT_WIDTH, settings, box, small, affine, prior)
    audits = child_scale_audit(settings, box, small, affine, prior)
    selected = next(row for row in audits if row["half_width_decimal"] == "1e-23")
    require(right["open_branch_child_floor"] == 496_283, "open branch child floor drift")
    require(right["wall_branch_child_floor"] == 496_284, "wall branch child floor drift")
    require(bool(arb(selected["two_level_child_variation_upper"]) < arb("0.002")), "selected child scale target failed")
    require(bool(arb(right["complete_radius_absolute_upper"]) < arb("0.06")), "right-wall source radius target failed")

    proof_boundary = (
        "One right-sided width-1e-23 physical box at the x=2/5 normalization wall, its exact split "
        "between child indices 496283 and 496284, two endpoint-complete affine Mordell branch boxes, "
        "and their complete source-current hull only. No left-side branch, full two-sided interior "
        "neighbourhood, recursive tiling, physical quadrature, non-A bound, joined R_after_A, R_Dir, "
        "Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "right_sided_width_1e_minus_23_affine_Mordell_wall_period_three_recursive_and_complete_source_box_certified_left_side_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "L": L,
            "K": K,
            "workers": 1,
            "right_half_width": str(RIGHT_WIDTH),
        },
        "production_settings": settings.__dict__,
        "exact_theorem": {
            "right_normalization": "a=31916/x-(A+2s)/2-1, tau=1/(2x)-1",
            "wall_child_index": 496_284,
            "open_right_child_index": 496_283,
            "period_three_base": "(p,a,tau)=(95747/6,1/3,-1)",
            "child_Taylor": "B0(rho)+2*pi*i*(alpha*B1(rho)+beta*B2(rho))+R2",
        },
        "certified_right_wall_box": right,
        "child_scale_audit": audits,
        "summary": {
            "certified_right_half_width": str(RIGHT_WIDTH),
            "open_branch_child_index": right["open_branch_child_floor"],
            "wall_branch_child_index": right["wall_branch_child_floor"],
            "branch_complete_balls_overlap": True,
            "selected_two_level_child_variation_upper": selected["two_level_child_variation_upper"],
            "complete_radius_absolute_upper": right["complete_radius_absolute_upper"],
            "complete_contains_exact_center_source": True,
        },
        "decision": {
            "right_normalization_branch_partitioned": True,
            "child_floor_wall_partitioned": True,
            "period_three_child_Taylor_boxes_built": True,
            "right_endpoint_complete_recursive_box_built": True,
            "right_complete_source_box_certified": True,
            "left_normalization_branch_built": False,
            "two_sided_interior_neighbourhood_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Construct the x<2/5 parity-half-shift branch, track its double-conjugation orientation, "
            "split its child floor from the shared wall line, and certify the opposite period-three "
            "endpoint-complete box before forming a two-sided hull."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": dependency_artifacts[1]["primary_source"],
        "dependencies": {
            "small_tau_result": {"path": relative(SMALL_RESULT), "sha256": file_hash(SMALL_RESULT)},
            "small_tau_builder": {"path": relative(SMALL_BUILDER), "sha256": file_hash(SMALL_BUILDER)},
            "affine_result": {"path": relative(AFFINE_RESULT), "sha256": file_hash(AFFINE_RESULT)},
            "affine_builder": {"path": relative(AFFINE_BUILDER), "sha256": file_hash(AFFINE_BUILDER)},
            "parameter_box_result": {"path": relative(BOX_RESULT), "sha256": file_hash(BOX_RESULT)},
            "parameter_box_builder": {"path": relative(BOX_BUILDER), "sha256": file_hash(BOX_BUILDER)},
            "point_result": {"path": relative(PRIOR_RESULT), "sha256": file_hash(PRIOR_RESULT)},
            "point_builder": {"path": relative(PRIOR_BUILDER), "sha256": file_hash(PRIOR_BUILDER)},
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
    print("certified right-side Mordell wall, period-three recursion, and complete source box", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
