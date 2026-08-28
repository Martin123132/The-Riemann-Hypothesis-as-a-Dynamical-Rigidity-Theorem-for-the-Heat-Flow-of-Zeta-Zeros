#!/usr/bin/env python3
"""Bound the complete analytic m>=B exterior block by two tangential IBPs."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_phase_coupled_tangential_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "rational_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_three_current_derivative_tail_gate.json",
    "remainder_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_fresnel_remainder_derivative_tail_gate.json",
    "dictionary_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_completed_B_outer_dictionary_derivative_tail_gate.json",
    "tangential_cutoff": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_tangential_exterior_cutoff_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
M = B
DELTA_TEXT = "1e-4"
TRANSITION_PANELS = 16
PHYSICAL_TARGET = "8e-10"

RATIONAL_TERMS = (
    (1, 2, 0, 0, 1), (1, 1, 1, 1, 1),
    (4, 2, 0, 4, 3), (4, 3, -1, 3, 3), (5, 2, 0, 2, 2), (3, 3, -1, 1, 2),
    (48, 3, -1, 7, 5), (48, 4, -2, 6, 5), (84, 3, -1, 5, 4),
    (60, 4, -2, 4, 4), (35, 3, -1, 3, 3), (15, 4, -2, 2, 3),
)
DICTIONARY_TERMS = ((1, 6, 0, 5), (10, 4, 2, 5), (5, 2, 4, 5))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing dependency: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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


def derivative_rows(order: int, p: int, k: int) -> tuple[tuple[int, int, int], ...]:
    candidates = {
        0: ((1, p, k),),
        1: ((p, p - 1, k), (2 * k, p + 1, k + 1)),
        2: ((p * (p - 1), p - 2, k), (2 * k * (2 * p + 1), p, k + 1), (4 * k * (k + 1), p + 2, k + 2)),
    }[order]
    return tuple(row for row in candidates if row[0])


def signed_power(base: arb, exponent: int) -> arb:
    return base**exponent if exponent >= 0 else 1 / base ** (-exponent)


def tail(exponent: int) -> arb:
    require(exponent > 1, "nonsummable outer comparison")
    start = arb(M)
    return start ** (-exponent) + start ** (1 - exponent) / (exponent - 1)


def rational_majorant(x_right: arb, order: int) -> arb:
    endpoint, pi = arb(B), arb.pi()
    cmax = endpoint * x_right / 2
    rho = 1 - (cmax / M) ** 2
    require(rho > arb(15) / 16 - arb("1e-90"), "rational denominator margin lost")
    total = arb(0)
    for numerator, pi_power, endpoint_power, p, k in RATIONAL_TERMS:
        coefficient = arb(numerator) / pi**pi_power * signed_power(endpoint, endpoint_power)
        for factor, cp, kp in derivative_rows(order, p, k):
            total += coefficient * factor * cmax**cp * rho**(-kp) * tail(2 * kp)
    return total * (endpoint / 2) ** order


def remainder_majorant(x_left: arb, x_right: arb, order: int) -> arb:
    pi = arb.pi()
    alpha, beta, gamma = arb(3) / 4, arb(5) / 4, arb(13) / 4
    a0 = arb(15) / (4 * pi**4 * alpha**7)
    a11 = arb(45) * beta / (8 * pi**3 * alpha**6)
    if order == 0:
        one_sign = a0 * x_right**2 * tail(6)
    elif order == 1:
        one_sign = (3 * a0 / 2) * x_right * tail(6) + a11 * tail(4)
    else:
        require(order == 2, "remainder derivative order outside zero through two")
        one_sign = (
            (15 * a0 / 4) * tail(6)
            + 3 * a11 / x_left * tail(4)
            + 45 * beta**2 / (8 * pi**2 * alpha**5 * x_left**2) * tail(2)
            + 15 * beta**2 / (2 * pi**3 * alpha**7 * x_left) * tail(4)
            + 45 * gamma / (16 * pi**3 * alpha**6 * x_left) * tail(4)
        )
    return 2 * one_sign


def dictionary_majorant(x_right: arb, order: int) -> arb:
    endpoint, pi = arb(B), arb.pi()
    cmax = endpoint * x_right / 2
    rho = 1 - (cmax / M) ** 2
    common = arb(3) / (pi**4 * endpoint**2)
    total = arb(0)
    for coefficient, p, mode_power, k in DICTIONARY_TERMS:
        for factor, cp, kp in derivative_rows(order, p, k):
            total += common * coefficient * factor * cmax**cp * rho**(-kp) * tail(2 * kp - mode_power)
    return total * (endpoint / 2) ** order


def weighted_outer_majorants(x_left: arb, x_right: arb) -> tuple[arb, arb, arb]:
    raw = [
        rational_majorant(x_right, order)
        + remainder_majorant(x_left, x_right, order)
        + dictionary_majorant(x_right, order)
        for order in range(3)
    ]
    weight = (x_left * (1 - x_left)) ** (-arb(1) / 4)
    log_first = (1 - 2 * x_left) / (4 * x_left * (1 - x_left))
    second_ratio = (12 * x_left**2 - 12 * x_left + 5) / (16 * x_left**2 * (1 - x_left) ** 2)
    return (
        weight * raw[0],
        weight * (raw[1] + log_first * raw[0]),
        weight * (raw[2] + 2 * log_first * raw[1] + second_ratio * raw[0]),
    )


def symbolic_certificate() -> dict[str, str]:
    x = sp.symbols("x", positive=True, real=True)
    amplitude = sp.Function("A")(x)
    phase_prime = sp.Function("h")(x)
    q = 1 / phase_prime
    second_ibp_remainder = sp.expand(sp.diff(q * sp.diff(q * amplitude, x), x))
    expected = sp.expand(
        q**2 * sp.diff(amplitude, x, 2)
        + 3 * q * sp.diff(q, x) * sp.diff(amplitude, x)
        + (sp.diff(q, x) ** 2 + q * sp.diff(q, x, 2)) * amplitude
    )
    require(sp.simplify(second_ibp_remainder - expected) == 0, "two-IBP operator drift")

    t, endpoint = sp.symbols("t B", positive=True, real=True)
    h = sp.pi * endpoint**2 / 4 - t / (2 * x * (1 - x))
    h1 = sp.factor(sp.diff(h, x))
    h2 = sp.factor(sp.diff(h, x, 2))
    smooth = 6 * x**5 - 15 * x**4 + 10 * x**3
    require(smooth.subs(x, 0) == 0 and sp.diff(smooth, x).subs(x, 0) == 0, "left cutoff jet drift")
    require(smooth.subs(x, 1) == 1 and sp.diff(smooth, x).subs(x, 1) == 0, "right cutoff jet drift")
    return {
        "phase": "Phi_B'(x)=h(x)=pi B^2/4-t/[2x(1-x)]",
        "phase_first_derivative": str(h1),
        "phase_second_derivative": str(h2),
        "two_IBP_boundary_0": "|A/h|",
        "two_IBP_boundary_1": "|A'|/|h|^2+|A||h'|/|h|^3",
        "two_IBP_integrand": "|A''|/|h|^2+3|A'||h'|/|h|^3+|A||h''|/|h|^3+3|A||h'|^2/|h|^4",
        "exterior_domains": "[delta,x_low] union [x_high,1/2]",
        "boundary_ownership": "eta and eta' vanish at delta; x_low, x_high, and 1/2 are all retained",
    }


def phase_geometry(x: arb) -> tuple[arb, arb, arb]:
    t, endpoint, pi = arb(T), arb(B), arb.pi()
    h = pi * endpoint**2 / 4 - t / (2 * x * (1 - x))
    h1 = t * (1 - 2 * x) / (2 * x**2 * (1 - x) ** 2)
    h2_abs = t * (3 * x**2 - 3 * x + 1) / (x**3 * (1 - x) ** 3)
    return h, h1, h2_abs


def eta_majorants(in_transition: bool) -> tuple[arb, arb, arb]:
    if not in_transition:
        return arb(1), arb(0), arb(0)
    delta = arb(DELTA_TEXT)
    return arb(1), arb(15) / (8 * delta), arb(10) * arb(3).sqrt() / (3 * delta**2)


def amplitude_majorants(x_left: arb, x_right: arb, in_transition: bool) -> tuple[arb, arb, arb]:
    u0, u1, u2 = weighted_outer_majorants(x_left, x_right)
    e0, e1, e2 = eta_majorants(in_transition)
    return e0 * u0, e1 * u0 + e0 * u1, e2 * u0 + 2 * e1 * u1 + e0 * u2


def slab_bound(x_left: arb, x_right: arb, left_branch: bool, in_transition: bool) -> arb:
    require(x_right > x_left > 0, "invalid tangential slab")
    a0, a1, a2 = amplitude_majorants(x_left, x_right, in_transition)
    h_edge, _, _ = phase_geometry(x_right if left_branch else x_left)
    h_min = -h_edge if left_branch else h_edge
    _, h1_max, h2_max = phase_geometry(x_left)
    require(h_min > 0 and h1_max >= 0 and h2_max > 0, "phase monotonicity margin lost")
    integrand = (
        a2 / h_min**2
        + 3 * a1 * h1_max / h_min**3
        + a0 * h2_max / h_min**3
        + 3 * a0 * h1_max**2 / h_min**4
    )
    return (x_right - x_left) * integrand


def boundary_bounds(x: arb) -> tuple[arb, arb]:
    a0, a1, _ = amplitude_majorants(x, x, False)
    h, h1, _ = phase_geometry(x)
    h_abs = abs(h)
    require(h_abs > 0, "boundary reached the trace saddle")
    return a0 / h_abs, a1 / h_abs**2 + a0 * h1 / h_abs**3


def partition_geometry() -> tuple[arb, arb, arb, list[arb], list[arb], list[arb]]:
    t, endpoint, pi = arb(T), arb(B), arb.pi()
    delta = arb(DELTA_TEXT)
    x0 = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x0) / (2 * x0**2 * (1 - x0) ** 2)
    root_hessian = hessian.sqrt()
    x_low, x_high = x0 - arb(70) / root_hessian, x0 + arb(70) / root_hessian

    transition = [delta + delta * index / TRANSITION_PANELS for index in range(TRANSITION_PANELS + 1)]

    left_interior = []
    xi = 70
    while True:
        xi *= 2
        candidate = x0 - arb(xi) / root_hessian
        if candidate <= 2 * delta:
            break
        left_interior.append(candidate)
    left = [2 * delta] + list(reversed(left_interior)) + [x_low]

    right = [x_high]
    xi = 70
    while True:
        xi *= 2
        candidate = x0 + arb(xi) / root_hessian
        if candidate >= arb(1) / 2:
            break
        right.append(candidate)
    right.append(arb(1) / 2)
    for points in (transition, left, right):
        require(all(points[index + 1] > points[index] for index in range(len(points) - 1)), "partition ordering failed")
    return x0, root_hessian, x_low, transition, left, right


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    x0, root_hessian, x_low, transition, left, right = partition_geometry()
    x_high = right[0]

    transition_integral = sum(
        (slab_bound(a, b, True, True) for a, b in zip(transition, transition[1:])), arb(0)
    )
    left_integral = sum((slab_bound(a, b, True, False) for a, b in zip(left, left[1:])), arb(0))
    right_integral = sum((slab_bound(a, b, False, False) for a, b in zip(right, right[1:])), arb(0))
    integral_remainder = transition_integral + left_integral + right_integral

    boundary_rows = {}
    boundary0 = arb(0)
    boundary1 = arb(0)
    for name, point in (("x_low", x_low), ("x_high", x_high), ("half", arb(1) / 2)):
        first, second = boundary_bounds(point)
        boundary_rows[name] = {
            "x_ball": point.str(PRECISION, more=True),
            "first_IBP_boundary_ball": first.str(PRECISION, more=True),
            "second_IBP_boundary_ball": second.str(PRECISION, more=True),
        }
        boundary0 += first
        boundary1 += second

    physical_factor = 2 * (arb.pi() / (32 * arb(T))) ** (arb(1) / 4)
    unnormalized = boundary0 + boundary1 + integral_remainder
    physical = physical_factor * unnormalized
    require(physical < arb(PHYSICAL_TARGET), f"outer phase-coupled cost exceeds {PHYSICAL_TARGET}")

    return {
        "height": T,
        "endpoint": B,
        "outer_mode_start": M,
        "delta": DELTA_TEXT,
        "trace_root_ball": x0.str(PRECISION, more=True),
        "sqrt_trace_hessian_ball": root_hessian.str(PRECISION, more=True),
        "window_low_ball": x_low.str(PRECISION, more=True),
        "window_high_ball": x_high.str(PRECISION, more=True),
        "partition": {
            "transition_panels": len(transition) - 1,
            "left_dyadic_panels": len(left) - 1,
            "right_dyadic_panels": len(right) - 1,
            "total_panels": len(transition) + len(left) + len(right) - 3,
        },
        "boundary_rows": boundary_rows,
        "first_IBP_boundary_sum_ball": boundary0.str(PRECISION, more=True),
        "second_IBP_boundary_sum_ball": boundary1.str(PRECISION, more=True),
        "transition_integral_remainder_ball": transition_integral.str(PRECISION, more=True),
        "left_integral_remainder_ball": left_integral.str(PRECISION, more=True),
        "right_integral_remainder_ball": right_integral.str(PRECISION, more=True),
        "two_IBP_integral_remainder_sum_ball": integral_remainder.str(PRECISION, more=True),
        "unnormalized_outer_bound_ball": unnormalized.str(PRECISION, more=True),
        "physical_projection_factor_ball": physical_factor.str(PRECISION, more=True),
        "physical_outer_bound_ball": physical.str(PRECISION, more=True),
        "physical_outer_upper_bound": PHYSICAL_TARGET,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    p = c["partition"]
    return f"""# Phase-coupled tangential bound for the complete analytic outer B block

Date: 2026-08-13

Status: saved-height, cutoff-uniform and Abel-uniform analytic outer-block
bound; not a proof of the complete completed-B estimate

The three preceding gates bound, for every `m>=B={B}`, the paired first
three rational currents, the exact two-sign Fresnel remainder, and the exact
Fresnel-to-boundary dictionary through two derivatives.  Let their joined
Kummer-weighted amplitude be `U(x)`, and put

```text
A(x)=eta_delta(x) U(x),
h(x)=Phi_B'(x)=pi B^2/4-t/[2x(1-x)].                 (OP1)
```

The exterior domain is `[delta,x_low] union [x_high,1/2]`.  Two exact
tangential integrations by parts give boundary majorants

```text
|A/h|,
|A'|/|h|^2+|A||h'|/|h|^3,                           (OP2)
```

and integral remainder

```text
integral {{ |A''|/|h|^2+3|A'||h'|/|h|^3
           +|A||h''|/|h|^3+3|A||h'|^2/|h|^4 }} dx. (OP3)
```

The amplitude and phase are not replaced by unrelated global suprema.
Instead, (OP3) is enclosed on `{p['total_panels']}` deterministic Arb slabs:
`{p['transition_panels']}` across the endpoint cutoff, `{p['left_dyadic_panels']}`
on the left exterior, and `{p['right_dyadic_panels']}` on the right exterior.
The exterior slabs are dyadic in normalized distance from the B saddle, so
the growth of the outer amplitude remains coupled to the much faster growth
of `|h|`.

All actual exterior boundary terms are retained.  The cutoff has
`A(delta)=A'(delta)=0`; the remaining three boundaries contribute

```text
sum |A/h| <= {c['first_IBP_boundary_sum_ball']},
sum second-boundary <= {c['second_IBP_boundary_sum_ball']}.         (OP4)
```

The complete integral remainder is

```text
{c['two_IBP_integral_remainder_sum_ball']}.                         (OP5)
```

After the exact equation-(9) half-domain projection factor,

```text
E_outer,m>=B <= {c['physical_outer_bound_ball']}
              < {c['physical_outer_upper_bound']}.                 (OP6)
```

This proves the complete analytic `m>=B` contribution to the released
completed-B exterior current is negligible at the saved height, uniformly
in every finite cutoff and Abel weight bounded by one.  No process or mode
enumeration is involved.  The same derivative majorants are summable and
independent of the regulator, so dominated convergence also gives this
outer block its own Abel limit without separating any unsummed endpoint
piece.

The finite block `m<B`, including its exact crossing completion, is still
open and must be compressed before norms.  The common-regulator joined
remainder is also open.  Thus this does not prove a complete exterior-B
bound, the `1.4058e-4` joined target, `T_upper`, a height-uniform theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.

Pi provenance: `pi` in (OP1)--(OP6) is inherited from the equation-(9)
Fresnel phase and exact physical projection.  No fitted geometric constant
is introduced.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["dictionary_tail"]["decision"]["complete_outer_mode_derivative_block_bounded"] is True, "outer derivative block drift")
    require(dependencies["tangential_cutoff"]["decision"]["exterior_B_trace_uniformly_tangentially_nonstationary"] is True, "tangential phase drift")

    artifact = {
        "kind": STEM,
        "status": "complete_analytic_m_ge_B_completed_B_exterior_block_below_8e_minus_10_by_phase_coupled_two_IBP",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "complete_outer_mode_derivative_block_used_with_local_phase": True,
            "all_outer_window_and_half_endpoint_IBP_terms_retained": True,
            "complete_analytic_m_ge_B_exterior_physical_contribution_below_8e_minus_10": True,
            "uniform_in_finite_cutoff_and_Abel_weight_bounded_by_one": True,
            "outer_block_separate_Abel_limit_justified_by_summable_derivative_majorants": True,
            "finite_mode_block_below_B_bounded": False,
            "completed_exterior_current_bound_proved": False,
            "joined_remainder_below_1p4058e_minus_4_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Compress the finite completed roster 1<=m<B using the original finite-Poisson or Dirichlet representation before taking first and second amplitude norms. Preserve the local crossing completion and all x_low/x_high/half endpoint terms, then combine that finite contribution with the certified m>=B cost and R_join,comp under the common regulator.",
        "proof_boundary": "Complete phase-coupled two-IBP estimate for the analytic m>=B outer completed-B block at t=1e10 only. No finite m<B exterior block, complete exterior-B or joined-remainder bound, complete T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified phase-coupled analytic outer B block below 8e-10", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
