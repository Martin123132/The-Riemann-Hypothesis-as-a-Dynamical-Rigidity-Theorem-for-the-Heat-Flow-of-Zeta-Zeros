#!/usr/bin/env python3
"""Certify a small-tau alternating Taylor box at the physical corner."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_small_tau_alternating_taylor_box_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
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
CORNER_WIDTH = Fraction(1, 10**22)
SCALE_POWERS = (18, 20, 21, 22, 24, 28)
BASE_A = Fraction(-1, 2)


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


def power_sum(count: int, degree: int) -> int:
    """Return sum_(h=0)^(count-1) h^degree exactly for 0 <= degree <= 5."""

    require(count >= 0, "negative power-sum count")
    require(0 <= degree <= 5, "unsupported power-sum degree")
    n = count
    if degree == 0:
        return n
    if degree == 1:
        return n * (n - 1) // 2
    if degree == 2:
        return n * (n - 1) * (2 * n - 1) // 6
    if degree == 3:
        return (n * (n - 1) // 2) ** 2
    if degree == 4:
        return (n - 1) * n * (2 * n - 1) * (3 * n * n - 3 * n - 1) // 30
    return (n - 1) ** 2 * n * n * (2 * n * n - 2 * n - 1) // 12


def residue_power_sum(n: int, residue: int, degree: int) -> int:
    """Return sum (residue+2h)^degree over that parity in 0..n."""

    require(residue in (0, 1), "parity residue must be zero or one")
    if n < residue:
        return 0
    count = (n - residue) // 2 + 1
    return sum(
        math.comb(degree, j)
        * residue ** (degree - j)
        * 2**j
        * power_sum(count, j)
        for j in range(degree + 1)
    )


def weighted_mass(p: Fraction, n: int, degree: int) -> Fraction:
    """Return S_degree=sum_(k=0)^n (p+k) k^degree without signs."""

    count = n + 1
    return p * power_sum(count, degree) + power_sum(count, degree + 1)


def alternating_weighted_moment(p: Fraction, n: int, degree: int) -> Fraction:
    """Return M_degree=sum_(k=0)^n (-1)^k (p+k) k^degree."""

    even = p * residue_power_sum(n, 0, degree) + residue_power_sum(n, 0, degree + 1)
    odd = p * residue_power_sum(n, 1, degree) + residue_power_sum(n, 1, degree + 1)
    return even - odd


def corner_boxes(width: Fraction, box: Any) -> dict[str, Any]:
    x_box = box.RationalBox(Fraction(1, 2) - width, Fraction(1, 2))
    s_box = box.RationalBox(Fraction(0), width)
    geometry = box.physical_geometry(x_box, s_box)
    p = Fraction(geometry["integer_shift_r"])
    n = geometry["transformed_m"]
    c_lo = Fraction(A) + 2 * s_box.lo
    c_hi = Fraction(A) + 2 * s_box.hi
    w_box = box.RationalBox(c_lo / 2 - p / x_box.lo, c_hi / 2 - p / x_box.hi)
    sigma_box = box.RationalBox(-Fraction(1, 2) / x_box.lo, -Fraction(1, 2) / x_box.hi)
    a_box = box.RationalBox(-w_box.hi, -w_box.lo)
    t_box = box.RationalBox(-sigma_box.hi - 1, -sigma_box.lo - 1)
    require(t_box.lo >= 0, "small-tau corner box crossed below zero")
    require(a_box.lo <= BASE_A <= a_box.hi, "alternating base phase left the normalized box")
    return {
        "width": width,
        "x": x_box,
        "s": s_box,
        "p": p,
        "n": n,
        "geometry": geometry,
        "original_w": w_box,
        "original_sigma": sigma_box,
        "normalized_a": a_box,
        "normalized_t": t_box,
    }


def small_tau_taylor_box(
    p: Fraction,
    n: int,
    a_box: Any,
    t_box: Any,
    box: Any,
) -> tuple[acb, dict[str, Any]]:
    """Enclose W_n(p;a,t) around the alternating point (-1/2,0)."""

    require(p >= 0 and n >= 0, "small-tau branch needs nonnegative weights")
    require(t_box.lo >= 0, "small-tau branch requires t >= 0")
    delta_a = max(abs(a_box.lo - BASE_A), abs(a_box.hi - BASE_A))
    t_max = t_box.hi
    moments = [alternating_weighted_moment(p, n, degree) for degree in range(3)]
    masses = {degree: weighted_mass(p, n, degree) for degree in range(2, 5)}

    alpha = box.arb_box(a_box) - box.arb_rational(BASE_A)
    tt = box.arb_box(t_box)
    phase_linear = alpha * box.arb_rational(moments[1]) + tt * box.arb_rational(moments[2])
    first_order = acb(box.arb_rational(moments[0])) + 2 * acb.pi() * acb(0, 1) * acb(phase_linear)
    remainder_mass = (
        delta_a * delta_a * masses[2]
        + 2 * delta_a * t_max * masses[3]
        + t_max * t_max * masses[4]
    )
    remainder_radius = 2 * arb.pi() ** 2 * box.arb_rational(remainder_mass)
    enclosure = first_order + box.symmetric_complex_error(remainder_radius)

    linear_radius = 2 * arb.pi() * box.arb_rational(
        delta_a * abs(moments[1]) + t_max * abs(moments[2])
    )
    return enclosure, {
        "base_a": str(BASE_A),
        "base_t": "0",
        "p": str(p),
        "n": n,
        "a_box": box.box_record(a_box),
        "t_box": box.box_record(t_box),
        "delta_a_max": str(delta_a),
        "t_max": str(t_max),
        "M0": str(moments[0]),
        "M1": str(moments[1]),
        "M2": str(moments[2]),
        "S2": str(masses[2]),
        "S3": str(masses[3]),
        "S4": str(masses[4]),
        "first_order_linear_radius_upper": box.arb_upper_text(linear_radius),
        "second_order_remainder_radius_upper": box.arb_upper_text(remainder_radius),
        "total_deviation_from_base_upper": box.arb_upper_text(linear_radius + remainder_radius),
        "first_order_enclosure": box.acb_record(first_order),
        "small_tau_enclosure": box.acb_record(enclosure),
    }


def corner_source_box(width: Fraction, settings: Any, box: Any, affine: Any, prior: Any) -> dict[str, Any]:
    flint.ctx.prec = settings.precision_bits
    data = corner_boxes(width, box)
    p = data["p"]
    n = data["n"]
    normalized, taylor_record = small_tau_taylor_box(
        p, n, data["normalized_a"], data["normalized_t"], box
    )
    original_transformed = normalized.conjugate()

    x_box = data["x"]
    s_box = data["s"]
    r = int(p)
    phase_box = box.RationalBox(
        Fraction(1, 4) + r * (A + 2 * s_box.lo) - Fraction(r * r) / x_box.lo,
        Fraction(1, 4) + r * (A + 2 * s_box.hi) - Fraction(r * r) / x_box.hi,
    )
    xx = box.arb_box(x_box)
    q_box = 2 / (xx * xx.sqrt())
    multiplier = acb(q_box) * box.cis_pi_interval(box.arb_box(phase_box))
    main = multiplier * original_transformed
    endpoint, endpoint_record = box.physical_endpoint_parameter_box(
        "small_tau_corner_joined_endpoint", x_box, s_box, settings, prior
    )
    complete = main + endpoint

    center_child, center_period = affine.weighted_period_ball(p, BASE_A, Fraction(0), n, prior)
    center_source, source_period = prior.source_current_ball(Fraction(1, 2), Fraction(0))
    require(normalized.contains(center_child), "small-tau box misses the exact alternating terminal")
    require(complete.contains(center_source), "complete small-tau box misses the exact corner source")
    return {
        "half_width": str(width),
        "x_box": box.box_record(x_box),
        "s_box": box.box_record(s_box),
        "first_integer_shift_r": r,
        "first_transformed_n": n,
        "original_w_box": box.box_record(data["original_w"]),
        "original_sigma_box": box.box_record(data["original_sigma"]),
        "normalization_identity": "W(p;w,sigma)=conjugate(W(p;-w,-sigma-1))",
        "small_tau_normalized_child": taylor_record,
        "normalized_child_contains_exact_terminal": True,
        "exact_terminal_period": center_period,
        "exact_terminal_current": box.acb_record(center_child),
        "first_phase_box": box.box_record(phase_box),
        "first_multiplier_box": box.acb_record(multiplier),
        "recursive_first_main_box": box.acb_record(main),
        "joined_endpoint": endpoint_record,
        "complete_source_box": box.acb_record(complete),
        "exact_corner_source": box.acb_record(center_source),
        "exact_corner_source_period": source_period,
        "complete_contains_exact_corner_source": True,
        "complete_radius_absolute_upper": box.arb_upper_text(complete.rad()),
        "passed": True,
    }


def scale_audit(settings: Any, box: Any) -> list[dict[str, Any]]:
    flint.ctx.prec = settings.precision_bits
    rows = []
    for power in SCALE_POWERS:
        width = Fraction(1, 10**power)
        data = corner_boxes(width, box)
        enclosure, record = small_tau_taylor_box(
            data["p"], data["n"], data["normalized_a"], data["normalized_t"], box
        )
        xx = box.arb_box(data["x"])
        q_upper = (2 / (xx * xx.sqrt())).abs_upper()
        child_deviation = arb(record["total_deviation_from_base_upper"])
        rows.append(
            {
                "half_width_decimal": f"1e-{power}",
                "half_width": str(width),
                "delta_a_max": record["delta_a_max"],
                "t_max": record["t_max"],
                "first_order_linear_radius_upper": record["first_order_linear_radius_upper"],
                "second_order_remainder_radius_upper": record["second_order_remainder_radius_upper"],
                "child_deviation_from_terminal_upper": record["total_deviation_from_base_upper"],
                "child_box_radius_upper": box.arb_upper_text(enclosure.rad()),
                "first_main_child_variation_upper": box.arb_upper_text(q_upper * child_deviation),
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    row = artifact["certified_corner_box"]
    child = row["small_tau_normalized_child"]
    audit_lines = [
        "| half-width | max |a+1/2| | max t | linear radius | quadratic remainder | propagated child variation |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for item in artifact["scale_audit"]:
        audit_lines.append(
            f"| `{item['half_width_decimal']}` | `{item['delta_a_max']}` | `{item['t_max']}` | "
            f"`{item['first_order_linear_radius_upper']}` | `{item['second_order_remainder_radius_upper']}` | "
            f"`{item['first_main_child_variation_upper']}` |"
        )
    return f"""# Small-`tau` alternating Taylor corner box

Date: 2026-08-24

Status: certified local interval lemma; recursive interior boxes remain open

## Status

Certified interval lemma at `t=10^10` for the one-sided physical box

```text
x in [1/2-10^-22, 1/2],   s in [0,10^-22].
```

This closes the local small-`tau` neighbourhood branch at the corner. It is
not a physical quadrature, a non-A bound, or an RH-level result.

## Exact normalization

Write

```text
e(u)=exp(2*pi*i*u),
W_n(p;a,t)=sum_(k=0)^n (p+k) exp(2*pi*i*(a*k+t*k^2)).
```

For this corner, `p={row['first_integer_shift_r']}` and
`n={row['first_transformed_n']}`. Integer curvature periodicity and
conjugation give the exact identity

```text
(ST1)  W_n(p;w,sigma)=conjugate(W_n(p;-w,-sigma-1)).
```

The normalized parameters are centred at `a_0=-1/2`, `t_0=0`, and the raw
`a` coordinate is retained continuously rather than reduced modulo one.

## Cancellation-preserving Taylor enclosure

Put `alpha=a+1/2`. Since
`exp(2*pi*i*(-k/2))=(-1)^k`, define the exact signed moments

```text
M_j=sum_(k=0)^n (-1)^k (p+k) k^j.
```

Then

```text
(ST2)  W_n(p;a,t)
       = M_0 + 2*pi*i*(alpha*M_1+t*M_2) + R_2,

M_0 = {child['M0']},
M_1 = {child['M1']},
M_2 = {child['M2']}.
```

Here `pi` is not an inserted geometric constant: the factor `2*pi` comes
directly from differentiating the defining phase
`exp(2*pi*i*(a*k+t*k^2))`.

For real `d`, the integral Taylor remainder gives
`|exp(i*d)-1-i*d|<=d^2/2`. Therefore, with
`S_j=sum_(k=0)^n (p+k)k^j`,

```text
(ST3)  |R_2| <= 2*pi^2 [delta_a^2*S_2
                         +2*delta_a*t_max*S_3
                         +t_max^2*S_4].
```

All `M_j` and `S_j` are evaluated by exact integer power sums split into the
even and odd residue classes. Only the quadratic Taylor remainder is placed
under an absolute-value bound.

## Certified source reassembly

The normalized child is conjugated back, multiplied by the first exact
Mordell prefactor, and joined to the already certified physical endpoint box:

```text
(ST4)  SourceBox = P_1 * conjugate(TaylorBox) + EndpointBox,
P_1 = 2*x^(-3/2) exp(pi*i*(1/4+r*(A+2s)-r^2/x)).
```

The complete source radius is at most
`{row['complete_radius_absolute_upper']}` and the box contains the exact
periodic source current at `(x,s)=(1/2,0)`. The exact normalized terminal has
period `{row['exact_terminal_period']}`.

## Scale audit

{chr(10).join(audit_lines)}

The certified half-width `10^-22` is one million times the earlier direct
first-difference width `10^-28`. This is a local Taylor gain, not a global
small-`tau` theorem.

## Validation

- Production: `{artifact['production_settings']['precision_bits']}` bits,
  Mordell cutoff `{artifact['production_settings']['cutoff']}`.
- The independent checker uses exact direct witnesses for every power-sum and
  parity formula, verifies the three full moments, and replays the endpoint
  and source box at 384 bits with cutoff 11.
- Builder and dependency hashes are pinned in the JSON artifact.

## Proof boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    for path in (AFFINE_BUILDER, AFFINE_RESULT, BOX_BUILDER, BOX_RESULT, PRIOR_BUILDER, PRIOR_RESULT, CHECKER):
        require(path.is_file(), f"missing dependency or checker: {path}")
    affine_artifact = json.loads(AFFINE_RESULT.read_text(encoding="utf-8"))
    box_artifact = json.loads(BOX_RESULT.read_text(encoding="utf-8"))
    prior_artifact = json.loads(PRIOR_RESULT.read_text(encoding="utf-8"))
    require(
        affine_artifact.get("passed") is True
        and box_artifact.get("passed") is True
        and prior_artifact.get("passed") is True,
        "Mordell dependency not passed",
    )
    affine = load_module("mordell_small_tau_affine_dependency", AFFINE_BUILDER)
    box = load_module("mordell_small_tau_box_dependency", BOX_BUILDER)
    prior = load_module("mordell_small_tau_point_dependency", PRIOR_BUILDER)
    settings = box.PRODUCTION

    corner = corner_source_box(CORNER_WIDTH, settings, box, affine, prior)
    audits = scale_audit(settings, box)
    require(corner["first_integer_shift_r"] == 39_894, "corner integer shift drift")
    require(corner["first_transformed_n"] == 1_240_710, "corner transformed length drift")
    require(corner["small_tau_normalized_child"]["M0"] == "660249", "M0 drift")
    require(corner["small_tau_normalized_child"]["M1"] == "794429714775", "M1 drift")
    require(corner["small_tau_normalized_child"]["M2"] == "985657301007258645", "M2 drift")
    require(
        bool(arb(corner["complete_radius_absolute_upper"]) < arb("0.02")),
        "complete corner radius target failed",
    )
    selected_audit = next(row for row in audits if row["half_width_decimal"] == "1e-22")
    require(
        bool(arb(selected_audit["first_main_child_variation_upper"]) < arb("0.01")),
        "small-tau propagated child target failed",
    )

    proof_boundary = (
        "A one-sided width-1e-22 small-tau affine-current Taylor box at the physical x=1/2, s=0 "
        "corner, its exact period-two terminal join, and one complete source-current box only. No "
        "uniform small-tau branch away from this corner, recursive interior parameter boxes, physical "
        "quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, "
        "PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "one_sided_width_1e_minus_22_small_tau_alternating_Taylor_corner_and_complete_source_box_certified_interior_recursive_boxes_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "L": L,
            "K": K,
            "workers": 1,
            "corner_half_width": str(CORNER_WIDTH),
        },
        "production_settings": settings.__dict__,
        "exact_theorem": {
            "phase_convention": "e(u)=exp(2*pi*i*u)",
            "normalization": "W(p;w,sigma)=conjugate(W(p;-w,-sigma-1))",
            "base_point": "(a_0,t_0)=(-1/2,0)",
            "first_order": "M_0+2*pi*i*((a+1/2)M_1+tM_2)",
            "remainder": "2*pi^2*(delta_a^2*S_2+2*delta_a*t_max*S_3+t_max^2*S_4)",
        },
        "certified_corner_box": corner,
        "scale_audit": audits,
        "summary": {
            "certified_half_width": str(CORNER_WIDTH),
            "width_gain_over_direct_1e_minus_28_box": 1_000_000,
            "complete_radius_absolute_upper": corner["complete_radius_absolute_upper"],
            "selected_first_main_child_variation_upper": selected_audit["first_main_child_variation_upper"],
            "exact_terminal_period": corner["exact_terminal_period"],
            "complete_contains_exact_corner_source": True,
        },
        "decision": {
            "small_tau_corner_neighbourhood_branch_built": True,
            "exact_tau_zero_terminal_joined_continuously": True,
            "complete_corner_source_box_certified": True,
            "uniform_small_tau_branch_away_from_corner_built": False,
            "recursive_interior_parameter_boxes_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Partition the x=2/5 affine normalization wall into exact parity, conjugation, and child-index "
            "branches; certify the one-sided recursive child boxes and hull them only after each joined "
            "endpoint current passes. Then reassess contraction and physical box-width scaling."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": affine_artifact["primary_source"],
        "dependencies": {
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
    print("certified small-tau alternating Taylor corner and complete source box", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
