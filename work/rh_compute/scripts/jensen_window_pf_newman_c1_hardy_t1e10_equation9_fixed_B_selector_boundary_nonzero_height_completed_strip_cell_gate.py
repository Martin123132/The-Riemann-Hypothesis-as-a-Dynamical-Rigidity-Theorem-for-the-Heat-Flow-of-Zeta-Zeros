#!/usr/bin/env python3
"""Certify a nonzero selector-height cell for the completed strip error."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate as point_gate


POINT_GATE = point_gate.RESULT
POINT_BUILDER = point_gate.BUILDER
EVENT_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json"
ENDPOINT_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_upper_endpoint_selector_boundary_stability_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 90
LOWER_ALPHA = 159_577
ADJACENT_ALPHA = LOWER_ALPHA + 2
FIXED_B = 5_122_423
RAY_RADIUS = 12
RAY_PANELS = 24
MAJORANT_PANELS = 960


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def interval_ball(left: arb, right: arb) -> arb:
    return arb((left + right) / 2, (right - left) / 2)


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def derivative_at_center(p: dict[str, arb]) -> tuple[arb, acb, dict[str, arb]]:
    i = acb(0, 1)
    rotation = (i * p["pi"] / 6).exp()

    def integrand(radius: acb, _: bool) -> acb:
        z = rotation * radius
        scaled = z / p["beta"]
        amplitude = scaled.cosh() ** (-arb(3) / 2)
        exact = amplitude * (i * p["beta"] ** 3 * point_gate.argument_minus_tanh_polynomial(scaled)).exp()
        canonical = (-radius**3 / 3).exp()
        return rotation * (-i * z) * (exact - canonical)

    ray = acb(0)
    radius = arb(RAY_RADIUS)
    for index in range(RAY_PANELS):
        left = radius * index / RAY_PANELS
        right = radius * (index + 1) / RAY_PANELS
        ray += acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb("1e-27"),
            rel_tol=arb("1e-27"),
            eval_limit=400_000,
            depth_limit=50,
        )

    rho = radius / p["beta"]
    tanh_tail = 3 * rho**18 / (1 - rho)
    phase_tail = p["beta"] ** 3 * tanh_tail
    replacement = radius**2 * phase_tail * phase_tail.exp()
    exact_tail = 2 * (
        8 * (-radius**3 / 4).exp() / (3 * radius)
        + (-p["beta"] ** 3 / 5).exp()
        * (p["beta"] / (p["beta"] ** 2 / 5) + 1 / (p["beta"] ** 2 / 5) ** 2)
    )
    canonical_tail = 2 * (-radius**3 / 3).exp() / radius
    total_error = replacement + exact_tail + canonical_tail
    derivative = 2 * ray.real + arb(0, total_error)
    return derivative, ray, {
        "polynomial_replacement": replacement,
        "exact_tail": exact_tail,
        "canonical_tail": canonical_tail,
        "total_error": total_error,
    }


def second_derivative_majorant(p: dict[str, arb], lambda_radius: arb) -> dict[str, arb | int]:
    radius = arb(RAY_RADIUS)

    def pointwise_difference_majorant(r: arb) -> arb:
        rho = r / p["beta"]
        cosh_remainder = rho**2 * rho.cosh() / 2
        require(cosh_remainder.upper() < arb("0.001"), "amplitude disk left its certified range")
        amplitude_error = arb(3) / 2 * cosh_remainder / (1 - cosh_remainder) ** (arb(5) / 2)
        phase_error = arb(0)
        for degree, coefficient in point_gate.TANH_COEFFICIENTS.items():
            if degree >= 5:
                phase_error += abs(arb(coefficient.numerator) / coefficient.denominator) * rho**degree
        phase_error = p["beta"] ** 3 * (phase_error + 3 * rho**18 / (1 - rho))
        return (
            (-r**3 / 3 + lambda_radius * r / 2).exp()
            * (amplitude_error + phase_error)
            * phase_error.exp()
        )

    compact_one_ray = arb(0)
    for index in range(MAJORANT_PANELS):
        left = radius * index / MAJORANT_PANELS
        right = radius * (index + 1) / MAJORANT_PANELS
        r = interval_ball(left, right)
        upper = (r**2 * pointwise_difference_majorant(r)).upper()
        compact_one_ray += (right - left) * upper

    exact_stationary_point = (10 * lambda_radius / 3).sqrt()
    exact_growth = lambda_radius * exact_stationary_point / 2 - exact_stationary_point**3 / 20
    canonical_stationary_point = (2 * lambda_radius).sqrt()
    canonical_growth = lambda_radius * canonical_stationary_point / 2 - canonical_stationary_point**3 / 12
    exact_compact_tail = exact_growth.exp() * arb(5) / 3 * (-radius**3 / 5).exp()
    canonical_tail = canonical_growth.exp() * arb(4) / 3 * (-radius**3 / 4).exp()
    linear_rate = p["beta"] ** 2 / 5 - lambda_radius / 2
    require(linear_rate.lower() > 0, "exact linear ray rate lost positivity")
    exact_far_tail = (-linear_rate * p["beta"]).exp() * (
        p["beta"] ** 2 / linear_rate
        + 2 * p["beta"] / linear_rate**2
        + 2 / linear_rate**3
    )
    total = 2 * (compact_one_ray + exact_compact_tail + exact_far_tail + canonical_tail)
    require(total < arb("1.3e-6"), "second-derivative majorant exceeds 1.3e-6")
    return {
        "compact_panel_count": MAJORANT_PANELS,
        "compact_one_ray_bound": compact_one_ray,
        "exact_compact_tail_bound": exact_compact_tail,
        "exact_far_tail_bound": exact_far_tail,
        "canonical_tail_bound": canonical_tail,
        "two_ray_second_derivative_bound": total,
    }


def event_geometry(p: dict[str, arb]) -> dict[str, Any]:
    old_nearest_distance = p["pi"] / 8
    cell_half_width = p["pi"] / 16
    old_margin = old_nearest_distance - cell_half_width

    minimum_numerator: int | None = None
    minimum_mode: int | None = None
    for mode in range(1, ADJACENT_ALPHA // 2 + 1):
        numerator = 8 * mode * (ADJACENT_ALPHA - 2 * mode) - LOWER_ALPHA**2
        if numerator > 0 and (minimum_numerator is None or numerator < minimum_numerator):
            minimum_numerator = numerator
            minimum_mode = mode
    require((minimum_numerator, minimum_mode) == (3103, 40094), "adjacent-selector first event drift")
    new_nearest_distance = p["pi"] * minimum_numerator / 8
    new_margin = new_nearest_distance - cell_half_width
    require(old_margin.lower() > 0 and new_margin.lower() > 0, "selector cell crosses an internal event")
    return {
        "cell_half_width": "pi/16",
        "cell_half_width_ball": cell_half_width.str(PRECISION, more=True),
        "old_selector_nearest_event_distance": "pi/8",
        "old_selector_event_margin_ball": old_margin.str(PRECISION, more=True),
        "adjacent_selector_first_event_mode": minimum_mode,
        "adjacent_selector_first_event_distance": "3103*pi/8",
        "adjacent_selector_event_margin_ball": new_margin.str(PRECISION, more=True),
        "cell_crosses_no_internal_turning_event": True,
    }


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for dependency in (POINT_GATE, POINT_BUILDER, EVENT_GATE, ENDPOINT_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    point = json.loads(POINT_GATE.read_text(encoding="utf-8"))
    require(point.get("passed") is True, "point completed-strip gate did not pass")
    ctx.dps = PRECISION
    p = point_gate.parameters()
    lambda_radius = p["pi"] / (16 * p["beta"])
    derivative, derivative_ray, derivative_errors = derivative_at_center(p)
    formal_derivative = -arb(13) / 60 * 2 * p["pi"] * arb(0).airy_ai() / p["beta"] ** 2
    require(arb("0.999999999") < derivative / formal_derivative < arb("1.000000001"), "center derivative misses formal beta^-2 term")
    majorant = second_derivative_majorant(p, lambda_radius)

    point_difference = arb(point["numerical_certificate"]["reduced_exact_minus_canonical_fold_ball"])
    uniform_reduced = abs(point_difference) + lambda_radius * abs(derivative) + lambda_radius**2 * majorant["two_ray_second_derivative_bound"] / 2
    strip_scale = arb(2).sqrt() / p["pi"] * p["Y"]
    height_lower = p["height"] - p["pi"] / 16
    height_upper = p["height"] + p["pi"] / 16
    maximum_normalizer = (p["pi"] / (32 * height_lower)) ** (arb(1) / 4)
    uniform_physical = strip_scale * uniform_reduced
    uniform_normalized = maximum_normalizer * uniform_physical
    require(uniform_reduced < arb("9.6e-12"), "uniform reduced cell error exceeds 9.6e-12")
    require(uniform_normalized < arb("8.9e-13"), "uniform normalized cell error exceeds 8.9e-13")
    events = event_geometry(p)

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate",
        "status": "event_free_nonzero_selector_height_completed_strip_error_cell_certified",
        "passed": True,
        "height_cell": {
            "center_height_ball": p["height"].str(PRECISION, more=True),
            "height_interval": "t*-pi/16 <= t <= t*+pi/16",
            "lower_height_ball": height_lower.str(PRECISION, more=True),
            "upper_height_ball": height_upper.str(PRECISION, more=True),
            "lambda_law": "lambda=(t*-t)/beta",
            "absolute_lambda_bound_ball": lambda_radius.str(PRECISION, more=True),
            "beta_is_height_independent": "beta^3=t*eta=pi*A^2/8",
            "Fourier_period_is_height_independent": "hY=2pi",
            "fixed_upper_endpoint_B": FIXED_B,
        },
        "center_derivative": {
            "ray_integral": complex_record(derivative_ray),
            "rigorous_derivative_ball": derivative.str(PRECISION, more=True),
            "formal_beta_minus_2_derivative_ball": formal_derivative.str(PRECISION, more=True),
            "exact_over_formal_ratio_ball": (derivative / formal_derivative).str(PRECISION, more=True),
            "polynomial_replacement_error_ball": derivative_errors["polynomial_replacement"].str(PRECISION, more=True),
            "exact_tail_error_ball": derivative_errors["exact_tail"].str(PRECISION, more=True),
            "canonical_tail_error_ball": derivative_errors["canonical_tail"].str(PRECISION, more=True),
            "total_error_radius_ball": derivative_errors["total_error"].str(PRECISION, more=True),
        },
        "second_derivative_majorant": {
            key: value if isinstance(value, int) else value.str(PRECISION, more=True)
            for key, value in majorant.items()
        },
        "uniform_error": {
            "center_reduced_error_ball": point_difference.str(PRECISION, more=True),
            "Taylor_bound": "|delta(lambda)|<=|delta(0)|+L|delta'(0)|+L^2 sup|delta''|/2",
            "uniform_reduced_completed_strip_error_bound_ball": uniform_reduced.str(PRECISION, more=True),
            "uniform_physical_completed_strip_error_bound_ball": uniform_physical.str(PRECISION, more=True),
            "maximum_equation9_normalizer_ball": maximum_normalizer.str(PRECISION, more=True),
            "uniform_equation9_normalized_error_bound_ball": uniform_normalized.str(PRECISION, more=True),
        },
        "event_geometry": events,
        "decision": {
            "nonzero_selector_height_cell_proved": True,
            "cell_crosses_selector_boundary": True,
            "cell_crosses_no_internal_turning_event": True,
            "all_mode_Fourier_completion_preserved_through_cell": True,
            "uniform_equation9_normalized_error_below_8_9e_minus_13": True,
            "ordinary_Morse_join_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "point_completed_strip_gate": {"path": relative(POINT_GATE), "sha256": file_hash(POINT_GATE)},
            "point_gate_builder": {"path": relative(POINT_BUILDER), "sha256": file_hash(POINT_BUILDER)},
            "turning_event_atlas_gate": {"path": relative(EVENT_GATE), "sha256": file_hash(EVENT_GATE)},
            "upper_endpoint_stability_gate": {"path": relative(ENDPOINT_GATE), "sha256": file_hash(ENDPOINT_GATE)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_obligation": (
            "Use the event-free selector buffer as the common endpoint of the two adjacent fixed-selector charts, "
            "then prove the Airy-to-ordinary-Morse remainder match at the first external event buffers without termwise outer-mode division."
        ),
        "proof_boundary": (
            "A rigorous completed-strip exact-minus-canonical error on t*-pi/16<=t<=t*+pi/16 only. "
            "No ordinary-Morse remainder join, complete T_upper theorem, Lambda<=0, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built nonzero selector cell: normalized completed-strip error < 8.9e-13")


def render_note(artifact: dict[str, Any]) -> str:
    h = artifact["height_cell"]
    d = artifact["center_derivative"]
    m = artifact["second_derivative_majorant"]
    u = artifact["uniform_error"]
    e = artifact["event_geometry"]
    return f"""# Fixed-B selector-boundary nonzero completed-strip height cell

Date: 2026-08-12

Status: rigorous event-free nonzero height cell; this is not a proof of the
ordinary-Morse join, `T_upper`, `Lambda<=0`, or RH

For fixed `A={LOWER_ALPHA}`, the fold scale is independent of height:

```text
beta^3=pi*A^2/8=t*,        lambda=(t*-t)/beta.         (HC1)
```

The exact Fourier period `hY=2pi` is therefore also fixed.  The all-mode
midpoint plus endpoint half-current continues to collapse the exact and
canonical completed strips to `y=0` for every real `t` in

```text
{h['height_interval']},
|lambda| <= {h['absolute_lambda_bound_ball']}.          (HC2)
```

Differentiating the already completed exact-minus-canonical fold integral at
the center gives

```text
delta'(0)={d['rigorous_derivative_ball']}.              (HC3)
```

Its formal first term is `-(13/60)2pi Ai(0)/beta^2`; the rigorous-to-formal
ratio is `{d['exact_over_formal_ratio_ball']}`.

A {m['compact_panel_count']}-panel interval majorant retains the exact-minus-
canonical difference before taking moduli and proves

```text
sup_|lambda|<=L |delta''(lambda)|
 <= {m['two_ray_second_derivative_bound']}.             (HC4)
```

Taylor's theorem applied to the completed object now gives

```text
|delta(lambda)|
 <= {u['uniform_reduced_completed_strip_error_bound_ball']},
equation-(9) normalized completed-strip error
 <= {u['uniform_equation9_normalized_error_bound_ball']} < 8.9e-13. (HC5)
```

The cell is geometrically clean.  The nearest old-selector internal event is
`pi/8` below `t*`, leaving margin
`{e['old_selector_event_margin_ball']}`.  The first adjacent-selector event is
mode `{e['adjacent_selector_first_event_mode']}` at `3103pi/8` above `t*`.
The fixed source endpoint remains `B={FIXED_B}` throughout.

This is a genuine nonzero selector-boundary comparison, not a sampled endpoint
check.  It still does not supply the Airy-to-ordinary-Morse remainder match at
the external event buffers.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No complete `T_upper`, `Lambda<=0`, RH, or prize-level conclusion is proved.
"""


if __name__ == "__main__":
    main()
