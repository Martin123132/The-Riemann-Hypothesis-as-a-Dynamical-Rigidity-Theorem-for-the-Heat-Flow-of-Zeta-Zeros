#!/usr/bin/env python3
"""Certify mirrored adaptive-width tiles and locate the local crossover."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_adaptive_width_mirrored_tile_crossover_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
RATIONAL_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_rational_phase_third_level_adjacent_tile_gate.py"
)
RATIONAL_RESULT = REPO_ROOT / (
    "work/rh_compute/results/"
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_rational_phase_third_level_adjacent_tile_gate.json"
)

A = 159_577
FIRST_P = Fraction(31_916)
FIRST_N = 992_568
OPEN_CHILD_N = 496_283
X0 = Fraction(2, 5)
S0 = Fraction(1, 3)
RIGHT_BASE_P = Fraction(95_747, 6)
RIGHT_BASE_A = Fraction(-1, 3)
LEFT_BASE_P = Fraction(47_873, 3)
LEFT_BASE_A = Fraction(1, 3)
WIDTH_UNIT = Fraction(1, 10**23)
SOURCE_RADIUS_TARGET = "0.01"
MAX_SEARCH_UNITS = 64


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


def exact_left_parent_boxes(x_box: Any, s_box: Any, box: Any) -> tuple[Any, Any]:
    a_box = box.RationalBox(
        Fraction(3, 2) + Fraction(A, 2) + s_box.lo - FIRST_P / x_box.lo,
        Fraction(3, 2) + Fraction(A, 2) + s_box.hi - FIRST_P / x_box.hi,
    )
    tau_box = box.RationalBox(
        Fraction(3, 2) - Fraction(1, 2) / x_box.lo,
        Fraction(3, 2) - Fraction(1, 2) / x_box.hi,
    )
    require(Fraction(0) < a_box.lo <= a_box.hi < 1, "left parent crossed an a wall")
    require(Fraction(0) < tau_box.lo <= tau_box.hi < Fraction(1, 4), "left tile escaped open branch")
    return a_box, tau_box


def exact_left_normalized_child_boxes(
    x_box: Any,
    s_box: Any,
    box: Any,
    rational: Any,
) -> tuple[Any, Any, Any, dict[str, Any]]:
    c_lo = Fraction(3, 2) + Fraction(A, 2) + s_box.lo
    c_hi = Fraction(3, 2) + Fraction(A, 2) + s_box.hi
    p_box = box.RationalBox(
        3 * FIRST_P - Fraction(3, 2) - Fraction(A, 2) - s_box.hi,
        3 * FIRST_P - Fraction(3, 2) - Fraction(A, 2) - s_box.lo,
    )
    raw_a = box.RationalBox(
        (c_lo * x_box.lo - FIRST_P) / (3 * x_box.lo - 1),
        (c_hi * x_box.hi - FIRST_P) / (3 * x_box.hi - 1),
    )
    raw_tau = box.RationalBox(
        -x_box.lo / (2 * (3 * x_box.lo - 1)),
        -x_box.hi / (2 * (3 * x_box.hi - 1)),
    )
    a_box = box.RationalBox(1 - raw_a.hi, 1 - raw_a.lo)
    tau_box = box.RationalBox(-raw_tau.hi - 1, -raw_tau.lo - 1)
    require(Fraction(-1, 2) < a_box.lo <= a_box.hi < Fraction(1, 2), "left child a wall crossed")
    require(tau_box.lo > 0, "left child did not enter positive small-tau branch")

    parent_a, parent_tau = exact_left_parent_boxes(x_box, s_box, box)
    old_p, _, _ = rational.box_dependency_child_boxes(parent_a, parent_tau, box)
    require(old_p.lo <= p_box.lo <= p_box.hi <= old_p.hi, "left exact weight escaped legacy hull")
    return p_box, a_box, tau_box, {
        "raw_child_a_box": box.box_record(raw_a),
        "raw_child_tau_box": box.box_record(raw_tau),
        "normalization": [
            "tau -> tau+1",
            "conjugate (a,tau) -> (-a,-tau)",
            "a -> a+1",
        ],
        "exact_weight_identity": "p_(2,L)=3*31916-3/2-A/2-s",
        "legacy_dependency_p_box": box.box_record(old_p),
        "exact_cancellation_p_box": box.box_record(p_box),
        "p_radius_reduction_factor": str(old_p.radius / p_box.radius),
    }


def left_parent_complete_box(
    x_box: Any,
    s_box: Any,
    settings: Any,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
    rational: Any,
) -> tuple[acb, dict[str, Any]]:
    a_parent, tau_parent = exact_left_parent_boxes(x_box, s_box, box)
    floor_box = box.RationalBox(2 * FIRST_N * tau_parent.lo, 2 * FIRST_N * tau_parent.hi)
    require(OPEN_CHILD_N < floor_box.lo <= floor_box.hi < OPEN_CHILD_N + 1, "left second floor wall crossed")
    p_child, a_child, tau_child, cancellation = exact_left_normalized_child_boxes(
        x_box, s_box, box, rational
    )
    child_normalized, child_record = rational.rational_phase_taylor_box(
        p_child,
        a_child,
        tau_child,
        OPEN_CHILD_N,
        LEFT_BASE_P,
        LEFT_BASE_A,
        Fraction(0),
        box,
        small,
        affine,
        prior,
    )
    child_original = child_normalized.conjugate()

    phase_box = box.RationalBox(
        Fraction(1, 4) - a_parent.hi * a_parent.hi / (2 * tau_parent.lo),
        Fraction(1, 4) - a_parent.lo * a_parent.lo / (2 * tau_parent.hi),
    )
    tt = box.arb_box(tau_parent)
    coefficient = box.cis_pi_interval(box.arb_box(phase_box)) / (2 * acb(tt) * acb(2 * tt).sqrt())
    main = coefficient * child_original

    lower_argument = box.RationalBox(
        Fraction(A, 2) + s_box.lo + Fraction(1, 2) - (FIRST_P - Fraction(1, 2)) / x_box.lo,
        Fraction(A, 2) + s_box.hi + Fraction(1, 2) - (FIRST_P - Fraction(1, 2)) / x_box.hi,
    )
    d = FIRST_P + FIRST_N + Fraction(1, 2)
    upper_argument = box.RationalBox(
        Fraction(A, 2) + s_box.lo + 3 * FIRST_N + Fraction(5, 2) - OPEN_CHILD_N - d / x_box.lo,
        Fraction(A, 2) + s_box.hi + 3 * FIRST_N + Fraction(5, 2) - OPEN_CHILD_N - d / x_box.hi,
    )
    mordell_tau = box.RationalBox(2 * tau_parent.lo, 2 * tau_parent.hi)
    lower_h, lower_hz, lower_record = box.mordell_target_parameter_box(
        lower_argument, mordell_tau, -1, settings, prior
    )
    upper_h, upper_hz, upper_record = box.mordell_target_parameter_box(
        upper_argument, mordell_tau, -1, settings, prior
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
    upper_phase = box.RationalBox(2 * half_n * upper_inner.lo, 2 * half_n * upper_inner.hi)
    pi_i = acb.pi() * acb(0, 1)
    endpoint = -acb(0, 1) * (
        box.cis_pi_interval(box.arb_box(lower_phase))
        * (acb(prior.arb_rational(2 * FIRST_P - 1)) * lower_h + lower_hz / pi_i)
        + ((-1) ** OPEN_CHILD_N)
        * box.cis_pi_interval(box.arb_box(upper_phase))
        * (
            acb(prior.arb_rational(2 * FIRST_P + 2 * FIRST_N + 1)) * upper_h
            + upper_hz / pi_i
        )
    ) / 4
    complete = main + endpoint
    third_floor = box.RationalBox(
        2 * OPEN_CHILD_N * tau_child.lo,
        2 * OPEN_CHILD_N * tau_child.hi,
    )
    require(0 < third_floor.lo <= third_floor.hi < 1, "left third floor wall crossed")
    return complete, {
        "parent_a_box": box.box_record(a_parent),
        "parent_tau_box": box.box_record(tau_parent),
        "second_floor_argument_box": box.box_record(floor_box),
        "normalized_child_p_box": box.box_record(p_child),
        "normalized_child_a_box": box.box_record(a_child),
        "normalized_child_tau_box": box.box_record(tau_child),
        "cancellation_audit": cancellation,
        "rational_phase_terminal": child_record,
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
        "complete_original_parent": box.acb_record(complete),
        "third_floor_argument_box": box.box_record(third_floor),
        "passed": True,
    }


def first_source_assembly(
    side: str,
    units: int,
    settings: Any,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
    rational: Any,
    previous_rational: dict[str, Any],
    previous_two_sided: dict[str, Any],
) -> dict[str, Any]:
    require(side in ("right", "left") and units >= 2, "invalid adaptive tile")
    if side == "right":
        x_box = box.RationalBox(X0 + WIDTH_UNIT, X0 + units * WIDTH_UNIT)
    else:
        x_box = box.RationalBox(X0 - units * WIDTH_UNIT, X0 - WIDTH_UNIT)
    s_box = box.RationalBox(S0 - WIDTH_UNIT, S0 + WIDTH_UNIT)

    if side == "right":
        parent, parent_record = rational.right_parent_complete_box(
            x_box, s_box, settings, box, small, affine, prior
        )
        original_transformed = parent.conjugate()
        previous_source = parse_acb(
            previous_rational["certified_adjacent_third_level_tile"]["complete_source_box"]
        )
    else:
        parent, parent_record = left_parent_complete_box(
            x_box, s_box, settings, box, small, affine, prior, rational
        )
        original_transformed = parent
        previous_source = parse_acb(
            previous_two_sided["certified_two_sided_wall_box"]["complete_source_box"]
        )

    first_phase = box.RationalBox(
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.lo) - FIRST_P * FIRST_P / x_box.lo,
        Fraction(1, 4) + int(FIRST_P) * (A + 2 * s_box.hi) - FIRST_P * FIRST_P / x_box.hi,
    )
    xx = box.arb_box(x_box)
    first_multiplier = acb(2 / (xx * xx.sqrt())) * box.cis_pi_interval(box.arb_box(first_phase))
    first_main = first_multiplier * original_transformed
    first_endpoint, endpoint_record = box.physical_endpoint_parameter_box(
        f"adaptive_{side}_{units}_unit_first_joined_endpoint", x_box, s_box, settings, prior
    )
    complete = first_main + first_endpoint
    require(complete.overlaps(previous_source), f"{side} adaptive tile lost shared-boundary overlap")

    terminal = parent_record["rational_phase_terminal"]
    first_upper = first_multiplier.abs_upper()
    second_upper = parse_acb(parent_record["second_main_coefficient"]).abs_upper()
    terminal_deviation = arb(terminal["total_deviation_from_base_upper"])
    propagation = first_upper * second_upper * terminal_deviation
    radius = complete.rad()
    return {
        "side": side,
        "width_units": units,
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
        "terminal_propagation_upper": box.arb_upper_text(propagation),
        "meets_source_radius_target": bool(radius < arb(SOURCE_RADIUS_TARGET)),
        "overlaps_previous_shared_boundary_box": True,
        "passed": True,
    }


def search_side(
    side: str,
    settings: Any,
    box: Any,
    small: Any,
    affine: Any,
    prior: Any,
    rational: Any,
    previous_rational: dict[str, Any],
    previous_two_sided: dict[str, Any],
) -> dict[str, Any]:
    rows = []
    previous_radius: arb | None = None
    first_failure = None
    for units in range(2, MAX_SEARCH_UNITS + 1):
        row = first_source_assembly(
            side,
            units,
            settings,
            box,
            small,
            affine,
            prior,
            rational,
            previous_rational,
            previous_two_sided,
        )
        radius = arb(row["complete_source_radius_upper"])
        if previous_radius is not None:
            require(radius > previous_radius, f"{side} radius search is not strictly increasing")
        previous_radius = radius
        rows.append(row)
        if not row["meets_source_radius_target"]:
            first_failure = row
            break
    require(first_failure is not None and len(rows) >= 2, f"{side} search did not bracket target")
    last_pass = rows[-2]
    require(last_pass["meets_source_radius_target"], f"{side} target has no passing tile")
    return {
        "side": side,
        "target": SOURCE_RADIUS_TARGET,
        "search_lattice_width_unit": str(WIDTH_UNIT),
        "rows": rows,
        "largest_passing_integer_units": last_pass["width_units"],
        "largest_passing_tile": last_pass,
        "first_failing_integer_units": first_failure["width_units"],
        "first_failing_tile": first_failure,
        "maximality_scope": "consecutive integer multiples of the declared width unit under this enclosure",
        "passed": True,
    }


def local_jacobian_and_wall_anchor() -> dict[str, Any]:
    n = OPEN_CHILD_N
    right_wall_x = Fraction(2 * n + 1, 5 * n + 2)
    left_wall_x = Fraction(2 * n + 1, 5 * n + 3)
    right_c = Fraction(A, 2) + S0 + 1
    right_raw_a = (FIRST_P - right_c * right_wall_x) / (1 - 2 * right_wall_x)
    right_a = -right_raw_a
    left_c = Fraction(3, 2) + Fraction(A, 2) + S0
    left_raw_a = (left_c * left_wall_x - FIRST_P) / (3 * left_wall_x - 1)
    left_a = 1 - left_raw_a
    wall_tau = Fraction(1, 2 * n)

    def period_from_tau_multiple(a0: Fraction) -> tuple[int, int]:
        for q in range(1, 13):
            period = n * q
            if (a0 * period + wall_tau * period * period).denominator == 1:
                return period, q
        raise RuntimeError("wall-anchor phase period not found")

    right_period, right_multiple = period_from_tau_multiple(right_a)
    left_period, left_multiple = period_from_tau_multiple(left_a)
    return {
        "wall_jacobian_at_x_2_over_5_s_1_over_3": {
            "right": {
                "d_p2_d_x": "0",
                "d_p2_d_s": "1",
                "d_a3_d_x": str(Fraction(2_393_675, 6)),
                "d_a3_d_s": "2",
                "d_t3_d_x": "25/2",
                "offset_direction": "x-2/5",
            },
            "left": {
                "d_p2_d_x": "0",
                "d_p2_d_s": "-1",
                "d_a3_d_x": str(Fraction(-1_196_825, 3)),
                "d_a3_d_s": "-2",
                "d_t3_d_x": "-25/2",
                "offset_direction": "2/5-x",
            },
        },
        "first_positive_floor_wall_anchor": {
            "child_n": n,
            "tau": str(wall_tau),
            "right_x": str(right_wall_x),
            "right_a": str(right_a),
            "right_minimal_phase_period": right_period,
            "right_period_multiple_of_child_n": right_multiple,
            "left_x": str(left_wall_x),
            "left_a": str(left_a),
            "left_minimal_phase_period": left_period,
            "left_period_multiple_of_child_n": left_multiple,
            "period_exceeds_current_length_on_both_sides": right_period > n and left_period > n,
            "consequence": "root-of-unity compression has no complete cycle before the first m3 wall",
        },
    }


def render_note(artifact: dict[str, Any]) -> str:
    right = artifact["right_width_search"]
    left = artifact["left_width_search"]
    anchor = artifact["local_jacobian_and_wall_anchor"]["first_positive_floor_wall_anchor"]
    rows = [
        "| side | units | source radius | terminal propagation | target met |",
        "|---|---:|---:|---:|---|",
    ]
    for search in (right, left):
        for item in search["rows"]:
            rows.append(
                f"| {item['side']} | `{item['width_units']}` | `{item['complete_source_radius_upper']}` | "
                f"`{item['terminal_propagation_upper']}` | `{item['meets_source_radius_target']}` |"
            )
    return f"""# Mirrored adaptive-width Mordell tile crossover

Date: 2026-08-25

Status: certified local width-lattice crossover; finite cover remains open

## Exact mirrored map

The right child weight from Section 11.470 is independent of `x`.  The left
map has the matching cancellation

```text
(AW1) p_(2,L)=3*31916-3/2-A/2-s.
```

After the integer curvature shift and conjugation, the left normalized child
is

```text
(AW2) a_(3,L)=1-[(3/2+A/2+s)x-31916]/(3x-1),
      t_(3,L)=x/[2(3x-1)]-1.
```

Its parent main, both joined Mordell endpoint currents, the physical main,
and both physical endpoint currents are retained in the same original-current
orientation as Section 11.469.

## Local width search

The fixed `s` half-width is `10^-23`.  For integer `M>=2`, the right search
uses `[2/5+10^-23,2/5+M*10^-23]` and the left search uses
`[2/5-M*10^-23,2/5-10^-23]`.  The declared complete-source radius target is
`{artifact['scope']['source_radius_target']}`.  Every row is a rigorous source
box; a failing row means only that this enclosure does not certify the target.

{chr(10).join(rows)}

The largest passing right and left lattice widths are respectively
`{right['largest_passing_integer_units']}` and
`{left['largest_passing_integer_units']}` units.  Their immediate next integer
candidates are explicitly certified but fail the declared radius target.
This is lattice maximality under the current wall-anchored Taylor enclosure,
not continuous or method-independent maximality.

## First positive floor-wall anchor

At `m_3=1`, both sides have `t=1/(2*496283)`.  The exact normalized anchors
have minimal phase periods

```text
right: {anchor['right_minimal_phase_period']}={anchor['right_period_multiple_of_child_n']}*496283,
left:  {anchor['left_minimal_phase_period']}={anchor['left_period_multiple_of_child_n']}*496283.
```

Both periods exceed the current length `496283`, so residue-period compression
has no complete cycle there.  Conversely, the next Mordell step contracts the
active index to one.  The missing bridge is therefore a cancellation-preserving
small-`t` endpoint expansion or grouped third step, not further fixed-wall
Taylor tiling.

## Validation and boundary

The production search uses 320 bits, `Y=10`, and tolerance `1e-55`.  The
independent checker rederives the left map and Jacobian, verifies both wall
anchor periods, and replays each passing/failing boundary pair at 384 bits,
`Y=11`, and tolerance `1e-65`.

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(RATIONAL_BUILDER.is_file() and RATIONAL_RESULT.is_file() and CHECKER.is_file(), "missing dependency or checker")
    previous_rational = json.loads(RATIONAL_RESULT.read_text(encoding="utf-8"))
    require(previous_rational.get("passed") is True, "rational-phase dependency not passed")
    rational = load_module("adaptive_width_rational_dependency", RATIONAL_BUILDER)
    two_sided = load_module("adaptive_width_two_sided_dependency", rational.TWO_SIDED_BUILDER)
    right_dependency = load_module("adaptive_width_right_dependency", two_sided.RIGHT_BUILDER)
    small = load_module("adaptive_width_small_dependency", right_dependency.SMALL_BUILDER)
    affine = load_module("adaptive_width_affine_dependency", right_dependency.AFFINE_BUILDER)
    box = load_module("adaptive_width_box_dependency", right_dependency.BOX_BUILDER)
    prior = load_module("adaptive_width_point_dependency", right_dependency.PRIOR_BUILDER)
    previous_two_sided = json.loads(rational.TWO_SIDED_RESULT.read_text(encoding="utf-8"))
    settings = box.PRODUCTION
    flint.ctx.prec = settings.precision_bits

    local = local_jacobian_and_wall_anchor()
    right = search_side(
        "right", settings, box, small, affine, prior, rational, previous_rational, previous_two_sided
    )
    left = search_side(
        "left", settings, box, small, affine, prior, rational, previous_rational, previous_two_sided
    )
    require(right["largest_passing_integer_units"] >= 2, "right tile did not widen")
    require(left["largest_passing_integer_units"] >= 2, "left tile did not widen")

    proof_boundary = (
        "Exact mirrored left child and endpoint formulas, exact local Jacobians and first-wall anchor "
        "periods, and discrete maximal passing/failing complete-source tiles on the declared 1e-23 "
        "integer lattice for the diagnostic radius target 0.01 only. No continuous maximal-width theorem, "
        "finite recursive cover, cancellation-preserving third Mordell endpoint expansion, physical "
        "quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, "
        "PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "mirrored_adaptive_width_lattice_crossover_certified_third_Mordell_endpoint_bridge_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "workers": 1,
            "width_unit": str(WIDTH_UNIT),
            "s_half_width": str(WIDTH_UNIT),
            "source_radius_target": SOURCE_RADIUS_TARGET,
            "max_search_units": MAX_SEARCH_UNITS,
        },
        "production_settings": settings.__dict__,
        "local_jacobian_and_wall_anchor": local,
        "right_width_search": right,
        "left_width_search": left,
        "summary": {
            "right_largest_passing_units": right["largest_passing_integer_units"],
            "right_first_failing_units": right["first_failing_integer_units"],
            "left_largest_passing_units": left["largest_passing_integer_units"],
            "left_first_failing_units": left["first_failing_integer_units"],
            "right_wall_anchor_period": local["first_positive_floor_wall_anchor"]["right_minimal_phase_period"],
            "left_wall_anchor_period": local["first_positive_floor_wall_anchor"]["left_minimal_phase_period"],
        },
        "decision": {
            "mirrored_left_exact_map_built": True,
            "child_weight_x_dependency_removed_on_both_sides": True,
            "right_discrete_width_crossover_certified": True,
            "left_discrete_width_crossover_certified": True,
            "fixed_wall_rational_phase_terminal_reaches_first_positive_floor_wall": False,
            "first_wall_root_of_unity_compression_practical": False,
            "cancellation_preserving_third_Mordell_endpoint_bridge_built": False,
            "finite_recursive_cover_built": False,
            "physical_quadrature_completed": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Derive the small-t grouped third Mordell identity in which the m=0 or m=1 main and both "
            "endpoint currents are expanded and cancelled before interval absolute values. Match that "
            "grouped branch to the rational-phase terminal on an overlap annulus, then test whether the "
            "combined evaluator can cross the first positive floor wall with a finite radius budget."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": previous_rational["primary_source"],
        "dependencies": {
            "rational_phase_result": {"path": relative(RATIONAL_RESULT), "sha256": file_hash(RATIONAL_RESULT)},
            "rational_phase_builder": {"path": relative(RATIONAL_BUILDER), "sha256": file_hash(RATIONAL_BUILDER)},
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
    print("certified mirrored adaptive-width tile crossovers and first-wall anchor periods", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
