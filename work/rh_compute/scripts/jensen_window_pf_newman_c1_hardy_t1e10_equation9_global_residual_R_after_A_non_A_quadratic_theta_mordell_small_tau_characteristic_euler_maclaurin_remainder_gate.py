#!/usr/bin/env python3
"""Certify a uniform Euler--Maclaurin remainder for the small-t characteristic currents."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any

import mpmath as mp


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
import flint
from flint import acb, arb, fmpq


STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_characteristic_euler_maclaurin_remainder_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
GROUPED_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_grouped_first_wall_bridge_gate"
)
GROUPED_RESULT = REPO_ROOT / f"work/rh_compute/results/{GROUPED_STEM}.json"
GROUPED_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{GROUPED_STEM}.py"

A = 159_577
FIRST_P = 31_916
N = 496_283
S0 = Fraction(1, 3)
T_MAX = Fraction(1, 2 * N)
RIGHT_P = Fraction(95_747, 6)
LEFT_P = Fraction(47_873, 3)
RIGHT_B = Fraction(-1, 3)
LEFT_B = Fraction(-2, 3)
EM_ORDER = 48
EM_REMAINDER_TARGET = arb("0.001")
PRECISION_BITS = 384
DIRECT_BLOCK_SIZE = 16
SAMPLE_FRACTIONS = (Fraction(0), Fraction(1, 16), Fraction(1, 4), Fraction(1, 2), Fraction(1))
CPU_BASELINE_SAMPLES = [13.71, 13.36, 16.56, 13.70, 21.06, 15.10]


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


def arb_rational(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def acb_cis_2pi(value: Fraction) -> acb:
    reduced = value % 1
    return (acb(0, 1) * 2 * acb.pi() * acb(arb_rational(reduced))).exp()


def arb_upper_text(value: arb, digits: int = 45) -> str:
    midpoint, radius, exponent = value.mid_rad_10exp()
    upper = mp.mpf(abs(int(midpoint)) + int(radius)) * mp.power(10, int(exponent))
    return mp.nstr(upper, digits, min_fixed=-10, max_fixed=10)


def acb_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": str(value.real),
        "imag_ball": str(value.imag),
        "radius_upper": arb_upper_text(value.rad()),
        "absolute_upper": arb_upper_text(value.abs_upper()),
    }


def symmetric_complex_error(radius: arb) -> acb:
    upper = radius.abs_upper()
    return acb(arb(0, upper), arb(0, upper))


def side_data(side: str) -> tuple[Fraction, Fraction]:
    require(side in ("right", "left"), "unknown side")
    return (RIGHT_P, RIGHT_B) if side == "right" else (LEFT_P, LEFT_B)


def phase_derivative_strip(side: str) -> dict[str, Any]:
    p, b = side_data(side)
    lower = b
    upper = b + 2 * T_MAX * (p + N)
    require(lower < 0 < upper, f"{side} phase derivative does not have one crossing")
    d = max(abs(lower), abs(upper))
    require(d < Fraction(7, 10), f"{side} phase derivative escaped the declared strip")
    entry_t = -b / (2 * (p + N))
    root_at_wall = -b / (2 * T_MAX) - p
    require(Fraction(0) < entry_t < T_MAX, f"{side} stationary point does not enter in the interval")
    require(Fraction(0) < root_at_wall < N, f"{side} wall stationary point escaped the current")
    return {
        "p": str(p),
        "integer_normalized_b": str(b),
        "phase": f"q(k)=t*k^2+({b}+2*t*{p})*k",
        "derivative_lower": str(lower),
        "derivative_upper": str(upper),
        "absolute_derivative_ceiling": str(d),
        "absolute_derivative_ceiling_below_0_7": True,
        "stationary_entry_t": str(entry_t),
        "stationary_root_at_first_wall": str(root_at_wall),
        "one_stationary_crossing_only": True,
    }


def hermite_envelope(order: int, d: Fraction, t_ceiling: Fraction) -> arb:
    """Upper bound for the order-th derivative of exp(2*pi*i*q)."""

    require(order >= 0 and d >= 0 and t_ceiling >= 0, "invalid Hermite envelope request")
    dd = arb_rational(d)
    tt = arb_rational(t_ceiling)
    two_pi = 2 * arb.pi()
    total = arb(0)
    factorial = math.factorial(order)
    for j in range(order // 2 + 1):
        coefficient = factorial // (math.factorial(j) * math.factorial(order - 2 * j))
        total += (
            coefficient
            * (two_pi * dd) ** (order - 2 * j)
            * (two_pi * tt) ** j
        )
    return total


def euler_maclaurin_remainder(side: str, em_order: int = EM_ORDER) -> dict[str, Any]:
    """Bound the EM remainder after B_2,...,B_(2M) endpoint corrections."""

    p, _ = side_data(side)
    strip = phase_derivative_strip(side)
    d = Fraction(strip["absolute_derivative_ceiling"])
    derivative_order = 2 * em_order
    h_m = hermite_envelope(derivative_order, d, T_MAX)
    h_prev = hermite_envelope(derivative_order - 1, d, T_MAX)
    first_weight_mass = arb_rational(p * N + Fraction(N * N, 2))
    derivative_integral = first_weight_mass * h_m + derivative_order * N * h_prev
    coefficient = (
        2
        * arb(derivative_order).zeta()
        / (2 * arb.pi()) ** derivative_order
    )
    remainder = coefficient * derivative_integral
    require(remainder < EM_REMAINDER_TARGET, f"{side} EM remainder misses target")
    return {
        **strip,
        "euler_maclaurin_order_M": em_order,
        "differentiation_order": derivative_order,
        "H_2M_upper": arb_upper_text(h_m),
        "H_2M_minus_1_upper": arb_upper_text(h_prev),
        "integrated_derivative_upper": arb_upper_text(derivative_integral),
        "periodic_Bernoulli_fourier_coefficient_upper": arb_upper_text(coefficient),
        "remainder_radius_upper": arb_upper_text(remainder),
        "remainder_target": "0.001",
        "passed": True,
    }


def exponential_derivative_factor(order: int, q_prime: Fraction, t: Fraction) -> acb:
    """Return E_r with d^r exp(2*pi*i*q)=exp(2*pi*i*q) E_r."""

    require(order >= 0, "negative exponential derivative order")
    i_two_pi = 2 * acb.pi() * acb(0, 1)
    q1 = acb(arb_rational(q_prime))
    tt = acb(arb_rational(t))
    total = acb(0)
    factorial = math.factorial(order)
    for j in range(order // 2 + 1):
        coefficient = factorial // (math.factorial(j) * math.factorial(order - 2 * j))
        total += coefficient * (i_two_pi * q1) ** (order - 2 * j) * (i_two_pi * tt) ** j
    return total


def endpoint_derivative(order: int, p: Fraction, b: Fraction, t: Fraction, x: int) -> acb:
    y = p + x
    q_prime = b + 2 * t * y
    phase = t * x * x + (b + 2 * t * p) * x
    e_phase = acb_cis_2pi(phase)
    value = acb(arb_rational(y)) * exponential_derivative_factor(order, q_prime, t)
    if order:
        value += order * exponential_derivative_factor(order - 1, q_prime, t)
    return e_phase * value


def characteristic_integral(p: Fraction, b: Fraction, t: Fraction) -> acb:
    """Evaluate integral_0^N (p+x)e(q(x)) dx, using the exact Fresnel primitive."""

    if t == 0:
        lam = 2 * acb.pi() * acb(0, 1) * acb(arb_rational(b))
        e_end = acb_cis_2pi(b * N)
        p_ball = acb(arb_rational(p))
        y_end = acb(arb_rational(p + N))
        return e_end * (y_end / lam - 1 / lam**2) - (p_ball / lam - 1 / lam**2)

    sqrt_t = arb_rational(t).sqrt()
    z0 = p + b / (2 * t)
    z1 = p + N + b / (2 * t)
    u0 = acb(2 * sqrt_t * arb_rational(z0))
    u1 = acb(2 * sqrt_t * arb_rational(z1))
    fresnel0 = u0.fresnel_c() + acb(0, 1) * u0.fresnel_s()
    fresnel1 = u1.fresnel_c() + acb(0, 1) * u1.fresnel_s()
    j0 = (fresnel1 - fresnel0) / (2 * sqrt_t)
    z_boundary = (
        acb_cis_2pi(t * z1 * z1) - acb_cis_2pi(t * z0 * z0)
    ) / (4 * acb.pi() * acb(0, 1) * acb(arb_rational(t)))
    j1 = z_boundary - acb(arb_rational(b / (2 * t))) * j0
    constant_phase = -t * p * p - b * p - b * b / (4 * t)
    return acb_cis_2pi(constant_phase) * j1


def euler_maclaurin_approximation(side: str, t: Fraction, em_order: int = EM_ORDER) -> acb:
    p, b = side_data(side)
    approximation = characteristic_integral(p, b, t)
    approximation += (endpoint_derivative(0, p, b, t, 0) + endpoint_derivative(0, p, b, t, N)) / 2
    for r in range(1, em_order + 1):
        derivative_order = 2 * r - 1
        bernoulli = fmpq.bernoulli(2 * r)
        coefficient = acb(arb(bernoulli) / math.factorial(2 * r))
        approximation += coefficient * (
            endpoint_derivative(derivative_order, p, b, t, N)
            - endpoint_derivative(derivative_order, p, b, t, 0)
        )
    return approximation


def direct_characteristic_current(side: str, t: Fraction, block_size: int) -> acb:
    p, b = side_data(side)
    a = b + 2 * t * p
    step = acb_cis_2pi(2 * t)
    total = acb(0)
    for start in range(0, N + 1, block_size):
        stop = min(N + 1, start + block_size)
        term = acb_cis_2pi(a * start + t * start * start)
        ratio = acb_cis_2pi(a + t * (2 * start + 1))
        for k in range(start, stop):
            total += acb(arb_rational(p + k)) * term
            term *= ratio
            ratio *= step
    require(total.is_finite(), f"{side} direct current is nonfinite")
    return total


def sample_row(side: str, fraction_of_wall: Fraction, remainder_upper: arb) -> dict[str, Any]:
    t = fraction_of_wall * T_MAX
    direct = direct_characteristic_current(side, t, DIRECT_BLOCK_SIZE)
    approximation = euler_maclaurin_approximation(side, t, EM_ORDER)
    enclosure = approximation + symmetric_complex_error(remainder_upper)
    require(enclosure.overlaps(direct), f"{side} EM enclosure misses direct current at t/T={fraction_of_wall}")
    error = direct - approximation
    return {
        "side": side,
        "t_over_first_wall": str(fraction_of_wall),
        "t": str(t),
        "direct_current": acb_record(direct),
        "euler_maclaurin_partial_sum": acb_record(approximation),
        "observed_difference_ball": acb_record(error),
        "certified_enclosure": acb_record(enclosure),
        "enclosure_overlaps_direct_current": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    right = artifact["uniform_remainder_contract"]["right"]
    left = artifact["uniform_remainder_contract"]["left"]
    return f"""# Small-t characteristic Euler--Maclaurin remainder gate

Date: 2026-08-25

Status: interval certificate; finite parameter transport remains open.

## Scope

This gate treats the two length-`{N + 1}` affine currents between the
`x=2/5` normalization wall and the first positive third-floor walls.  It does
not yet build the physical parameter-cell cover.

## Exact characteristic phases

At `s=1/3`, after removing an integer multiple of `k`, both currents have

```text
W(p,b,t) = sum_(k=0)^n (p+k) e(t k^2 + (b+2pt)k),
e(u)=exp(2 pi i u),   0 <= t <= 1/(2n).
```

The right pair is `(p,b)=({RIGHT_P},{RIGHT_B})`; the left pair is
`(p,b)=({LEFT_P},{LEFT_B})`.  Their phase derivatives stay in the exact
strips

```text
right: [{right['derivative_lower']}, {right['derivative_upper']}],
left:  [{left['derivative_lower']}, {left['derivative_upper']}].
```

Each strip crosses zero once and has absolute ceiling below `0.7`.  Thus the
apparent small-`t` Mordell singularity has disappeared from the finite current;
only one ordinary quadratic stationary point migrates into each interval.

## Explicit remainder

Writing `q'(x)=b+2t(p+x)`, direct differentiation gives

```text
d^r e(q(x))/dx^r = e(q(x)) E_r(q'(x),t),
H_r(D,T) = r! sum_(j=0)^(floor(r/2))
           (2 pi D)^(r-2j) (2 pi T)^j / (j! (r-2j)!),
|E_r| <= H_r(D,T).
```

Since the affine weight is `p+x`,

```text
integral |d^m((p+x)e(q(x)))/dx^m| dx
 <= (pn+n^2/2) H_m + mn H_(m-1).
```

The Fourier series of the periodic Bernoulli polynomial gives the exact
Euler--Maclaurin remainder coefficient

```text
2 zeta(2M)/(2 pi)^(2M).
```

For `M={EM_ORDER}` this proves the uniform bounds

| side | whole-interval EM remainder |
|---|---:|
| right | `< {right['remainder_radius_upper']}` |
| left | `< {left['remainder_radius_upper']}` |

Both are below `0.001`.  The production run also evaluates the exact Fresnel
primitive and all endpoint corrections through `B_96`, then overlaps the
result with direct phase-reset sums at `t/t_* = 0, 1/16, 1/4, 1/2, 1` on
both sides.

## Pi provenance

Every `pi` here has one of two explicit origins.  The factors `2 pi` in the
derivatives and Fresnel primitive come from the declared character
`e(u)=exp(2 pi i u)`.  The denominator `(2 pi)^(2M)` comes from the Fourier
series of the periodic Bernoulli polynomial.  No geometric circle or polygon
constant is inserted into the theta-current model.

## Source audit

Kuznetsov's truncated-theta algorithm identifies Euler--Maclaurin summation as
the small-parameter branch when `|tau|<1/n`.  Our whole interval satisfies
`0<=t<=1/(2n)<1/n`.  The paper supplies the algorithmic regime cue only; all
characteristic identities, derivative envelopes, constants, and interval
checks in this gate are derived here.

## Boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    require(GROUPED_RESULT.is_file() and GROUPED_BUILDER.is_file(), "missing grouped-wall dependency")
    grouped = json.loads(GROUPED_RESULT.read_text(encoding="utf-8"))
    require(grouped.get("passed") is True, "grouped-wall dependency not passed")

    flint.ctx.prec = PRECISION_BITS
    contracts = {side: euler_maclaurin_remainder(side) for side in ("right", "left")}
    sample_matrix = []
    for side in ("right", "left"):
        remainder_upper = arb(contracts[side]["remainder_radius_upper"])
        for sample_fraction in SAMPLE_FRACTIONS:
            sample_matrix.append(sample_row(side, sample_fraction, remainder_upper))

    proof_boundary = (
        "An exact characteristic reduction at s=1/3, one stationary-crossing inventory, and a uniform "
        "order-48 Euler--Maclaurin remainder below 0.001 for each of the two finite child currents on "
        "0<=t<=1/(2n), plus five direct-sum overlaps per side only. The endpoint/Fresnel partial sum is "
        "pointwise; no finite t- or s-cell transport theorem, endpoint-complete physical source cover, "
        "physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, "
        "PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "small_tau_characteristic_euler_maclaurin_uniform_remainder_certified_parameter_transport_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "n": N,
            "s": str(S0),
            "t_interval": ["0", str(T_MAX)],
            "euler_maclaurin_order": EM_ORDER,
            "precision_bits": PRECISION_BITS,
            "workers": 1,
            "direct_phase_reset_block_size": DIRECT_BLOCK_SIZE,
            "daytime_cpu_baseline_samples_percent": CPU_BASELINE_SAMPLES,
            "daytime_cpu_baseline_average_percent": round(sum(CPU_BASELINE_SAMPLES) / len(CPU_BASELINE_SAMPLES), 2),
        },
        "exact_characteristic_identities": {
            "right": {
                "p": str(RIGHT_P),
                "b": str(RIGHT_B),
                "a_identity": "a_R=2*s-1+2*p_R*t",
                "p_identity": "p_R=A/2+s+1-2*31916",
            },
            "left": {
                "p": str(LEFT_P),
                "b": str(LEFT_B),
                "a_identity_mod_Z": "a_L=-2*s+2*p_L*t (mod 1)",
                "p_identity": "p_L=3*31916-3/2-A/2-s",
            },
            "common_square_identity": "t*k^2+(b+2*p*t)*k=t*((p+k)^2-p^2)+b*k",
        },
        "uniform_remainder_contract": contracts,
        "pointwise_feasibility_matrix": sample_matrix,
        "decision": {
            "whole_small_tau_interval_in_one_frequency_strip": True,
            "one_stationary_crossing_per_side": True,
            "uniform_EM_remainder_below_0_001": True,
            "direct_sum_overlap_at_all_ten_samples": True,
            "finite_parameter_cover_built": False,
            "physical_quadrature_completed": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Differentiate the finite Euler--Maclaurin/Fresnel representation with respect to t and s, "
            "certify its parameter transport on geometrically selected cells, and propagate each current "
            "ball through the complete joined-endpoint source assembly."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": {
            "author": "Alexey Kuznetsov",
            "title": "Computing the truncated theta function via Mordell integral",
            "arxiv": "1306.4081v2",
            "url": "https://arxiv.org/abs/1306.4081v2",
            "used_for": "small-parameter Euler--Maclaurin regime cue only",
            "imported_numerical_constants": False,
        },
        "dependencies": {
            "grouped_wall_result": {"path": relative(GROUPED_RESULT), "sha256": file_hash(GROUPED_RESULT)},
            "grouped_wall_builder": {"path": relative(GROUPED_BUILDER), "sha256": file_hash(GROUPED_BUILDER)},
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
    print("certified characteristic Euler-Maclaurin remainder and ten direct overlaps", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
