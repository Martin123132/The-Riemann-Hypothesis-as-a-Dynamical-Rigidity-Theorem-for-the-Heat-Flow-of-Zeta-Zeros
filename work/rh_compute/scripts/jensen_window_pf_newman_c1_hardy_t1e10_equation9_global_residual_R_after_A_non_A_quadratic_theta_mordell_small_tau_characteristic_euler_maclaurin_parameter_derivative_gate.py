#!/usr/bin/env python3
"""Certify t and s derivatives of the small-t characteristic currents."""

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
    "non_A_quadratic_theta_mordell_small_tau_characteristic_euler_maclaurin_parameter_derivative_gate"
)
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
BASE_STEM = (
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_"
    "non_A_quadratic_theta_mordell_small_tau_characteristic_euler_maclaurin_remainder_gate"
)
BASE_RESULT = REPO_ROOT / f"work/rh_compute/results/{BASE_STEM}.json"
BASE_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{BASE_STEM}.py"
BASE_CHECKER = REPO_ROOT / f"work/rh_compute/scripts/check_{BASE_STEM}.py"

N = 496_283
T_MAX = Fraction(1, 2 * N)
RIGHT = (Fraction(95_747, 6), Fraction(-1, 3), 1)
LEFT = (Fraction(47_873, 3), Fraction(-2, 3), -1)
EM_ORDER = 80
PRECISION_BITS = 448
DIRECT_BLOCK_SIZE = 16
SAMPLES = (Fraction(0), Fraction(1, 16), Fraction(1, 4), Fraction(1, 2), Fraction(1))
T_DERIVATIVE_REMAINDER_TARGET = arb("0.05")
S_DERIVATIVE_REMAINDER_TARGET = arb("0.000001")
DIAGNOSTIC_VARIATION_RADIUS = arb("0.001")
CPU_BASELINE_SAMPLES = [39.07, 40.91, 34.61, 52.37, 40.96, 35.30]


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


def rat(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def cis(value: Fraction) -> acb:
    return (2 * acb.pi() * acb(0, 1) * acb(rat(value % 1))).exp()


def upper_text(value: arb, digits: int = 45) -> str:
    midpoint, radius, exponent = value.mid_rad_10exp()
    upper = mp.mpf(abs(int(midpoint)) + int(radius)) * mp.power(10, int(exponent))
    return mp.nstr(upper, digits, min_fixed=-10, max_fixed=10)


def acb_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": str(value.real),
        "imag_ball": str(value.imag),
        "radius_upper": upper_text(value.rad()),
        "absolute_upper": upper_text(value.abs_upper()),
    }


def symmetric_error(radius: arb) -> acb:
    upper = radius.abs_upper()
    return acb(arb(0, upper), arb(0, upper))


def side_data(side: str) -> tuple[Fraction, Fraction, int]:
    require(side in ("right", "left"), "unknown side")
    return RIGHT if side == "right" else LEFT


def derivative_ceiling(side: str) -> Fraction:
    p, b, _ = side_data(side)
    lower = b
    upper = b + 2 * T_MAX * (p + N)
    require(lower < 0 < upper, f"{side} stationary inventory drift")
    d = max(abs(lower), abs(upper))
    require(d < Fraction(7, 10), f"{side} derivative ceiling drift")
    return d


def hermite_envelope(order: int, d: Fraction) -> arb:
    total = arb(0)
    two_pi = 2 * arb.pi()
    factorial = math.factorial(order)
    for j in range(order // 2 + 1):
        coefficient = factorial // (math.factorial(j) * math.factorial(order - 2 * j))
        total += coefficient * (two_pi * rat(d)) ** (order - 2 * j) * (two_pi * rat(T_MAX)) ** j
    return total


def derivative_integrals(kind: str, p: Fraction) -> list[Fraction]:
    y0 = p
    y1 = p + N
    if kind == "one":
        values = [Fraction(N)]
    elif kind == "quadratic":
        values = [
            (y1**3 - y0**3) / 3 - p * (y1**2 - y0**2) / 2,
            (y1**2 - p * y1) - (y0**2 - p * y0),
            Fraction(2 * N),
        ]
    elif kind == "cubic":
        values = [
            (y1**4 - y0**4) / 4 - p * p * (y1**2 - y0**2) / 2,
            (y1**3 - p * p * y1) - (y0**3 - p * p * y0),
            3 * (y1**2 - y0**2),
            Fraction(6 * N),
        ]
    else:
        raise RuntimeError(f"unknown polynomial kind {kind}")
    require(all(value >= 0 for value in values), f"negative derivative integral for {kind}")
    return values


def polynomial_remainder(side: str, kind: str) -> arb:
    p, _, _ = side_data(side)
    d = derivative_ceiling(side)
    m = 2 * EM_ORDER
    integral = arb(0)
    for r, mass in enumerate(derivative_integrals(kind, p)):
        integral += math.comb(m, r) * rat(mass) * hermite_envelope(m - r, d)
    return 2 * arb(m).zeta() * integral / (2 * arb.pi()) ** m


def remainder_contract(side: str) -> dict[str, Any]:
    rem_one = polynomial_remainder(side, "one")
    rem_quadratic = polynomial_remainder(side, "quadratic")
    rem_cubic = polynomial_remainder(side, "cubic")
    rem_t = 2 * arb.pi() * rem_cubic
    rem_s = rem_one + 4 * arb.pi() * (1 + rat(T_MAX)) * rem_quadratic
    require(rem_t < T_DERIVATIVE_REMAINDER_TARGET, f"{side} t-derivative remainder misses target")
    require(rem_s < S_DERIVATIVE_REMAINDER_TARGET, f"{side} s-derivative remainder misses target")
    return {
        "side": side,
        "euler_maclaurin_order_M": EM_ORDER,
        "differentiation_order": 2 * EM_ORDER,
        "phase_derivative_ceiling": str(derivative_ceiling(side)),
        "unit_weight_remainder_upper": upper_text(rem_one),
        "quadratic_weight_remainder_upper": upper_text(rem_quadratic),
        "cubic_weight_remainder_upper": upper_text(rem_cubic),
        "partial_t_W_remainder_upper": upper_text(rem_t),
        "partial_s_W_remainder_upper": upper_text(rem_s),
        "partial_t_target": "0.05",
        "partial_s_target": "0.000001",
        "passed": True,
    }


def exponential_factor(order: int, q_prime: Fraction, t: Fraction) -> acb:
    i_two_pi = 2 * acb.pi() * acb(0, 1)
    total = acb(0)
    factorial = math.factorial(order)
    for j in range(order // 2 + 1):
        coefficient = factorial // (math.factorial(j) * math.factorial(order - 2 * j))
        total += coefficient * (i_two_pi * acb(rat(q_prime))) ** (order - 2 * j) * (
            i_two_pi * acb(rat(t))
        ) ** j
    return total


def polynomial_derivative(coefficients: tuple[Fraction, ...], order: int, y: Fraction) -> Fraction:
    total = Fraction(0)
    for degree in range(order, len(coefficients)):
        total += (
            coefficients[degree]
            * math.factorial(degree)
            / math.factorial(degree - order)
            * y ** (degree - order)
        )
    return total


def endpoint_derivative(
    order: int,
    p: Fraction,
    b: Fraction,
    t: Fraction,
    coefficients: tuple[Fraction, ...],
    x: int,
) -> acb:
    y = p + x
    q_prime = b + 2 * t * y
    phase = t * x * x + (b + 2 * t * p) * x
    total = acb(0)
    for r in range(min(order, len(coefficients) - 1) + 1):
        amplitude = polynomial_derivative(coefficients, r, y)
        total += math.comb(order, r) * acb(rat(amplitude)) * exponential_factor(order - r, q_prime, t)
    return cis(phase) * total


def characteristic_moments(p: Fraction, b: Fraction, t: Fraction, max_degree: int) -> list[acb]:
    require(max_degree <= 3, "moment primitive implemented only through degree three")
    y0 = p
    y1 = p + N
    if t == 0:
        lam = 2 * acb.pi() * acb(0, 1) * acb(rat(b))
        e_end = cis(b * N)
        moments = [(e_end - 1) / lam]
        for degree in range(1, max_degree + 1):
            boundary = acb(rat(y1**degree)) * e_end - acb(rat(y0**degree))
            moments.append(boundary / lam - degree * moments[-1] / lam)
        return moments

    sqrt_t = rat(t).sqrt()
    z0 = p + b / (2 * t)
    z1 = p + N + b / (2 * t)
    u0 = acb(2 * sqrt_t * rat(z0))
    u1 = acb(2 * sqrt_t * rat(z1))
    fresnel0 = u0.fresnel_c() + acb(0, 1) * u0.fresnel_s()
    fresnel1 = u1.fresnel_c() + acb(0, 1) * u1.fresnel_s()
    constant = cis(-t * p * p - b * p - b * b / (4 * t))
    i0 = constant * (fresnel1 - fresnel0) / (2 * sqrt_t)
    i1 = constant * (
        (cis(t * z1 * z1) - cis(t * z0 * z0))
        / (4 * acb.pi() * acb(0, 1) * acb(rat(t)))
        - acb(rat(b / (2 * t))) * (fresnel1 - fresnel0) / (2 * sqrt_t)
    )
    moments = [i0, i1]
    for degree in range(2, max_degree + 1):
        boundary = acb(rat(y1 ** (degree - 1))) * cis(
            t * N * N + (b + 2 * t * p) * N
        ) - acb(rat(y0 ** (degree - 1)))
        numerator = (
            boundary
            - (degree - 1) * moments[degree - 2]
            - 2 * acb.pi() * acb(0, 1) * acb(rat(b)) * moments[degree - 1]
        )
        moments.append(numerator / (4 * acb.pi() * acb(0, 1) * acb(rat(t))))
    return moments[: max_degree + 1]


def polynomial_integral(
    p: Fraction,
    b: Fraction,
    t: Fraction,
    coefficients: tuple[Fraction, ...],
) -> acb:
    moments = characteristic_moments(p, b, t, len(coefficients) - 1)
    return sum((acb(rat(coefficient)) * moments[degree] for degree, coefficient in enumerate(coefficients)), acb(0))


def em_polynomial(
    p: Fraction,
    b: Fraction,
    t: Fraction,
    coefficients: tuple[Fraction, ...],
) -> acb:
    approximation = polynomial_integral(p, b, t, coefficients)
    approximation += (
        endpoint_derivative(0, p, b, t, coefficients, 0)
        + endpoint_derivative(0, p, b, t, coefficients, N)
    ) / 2
    for r in range(1, EM_ORDER + 1):
        order = 2 * r - 1
        coefficient = acb(arb(fmpq.bernoulli(2 * r)) / math.factorial(2 * r))
        approximation += coefficient * (
            endpoint_derivative(order, p, b, t, coefficients, N)
            - endpoint_derivative(order, p, b, t, coefficients, 0)
        )
    return approximation


def direct_polynomials(side: str, t: Fraction) -> tuple[acb, acb, acb]:
    p, b, _ = side_data(side)
    a = b + 2 * t * p
    step = cis(2 * t)
    one = acb(0)
    quadratic = acb(0)
    cubic = acb(0)
    for start in range(0, N + 1, DIRECT_BLOCK_SIZE):
        stop = min(N + 1, start + DIRECT_BLOCK_SIZE)
        term = cis(a * start + t * start * start)
        ratio = cis(a + t * (2 * start + 1))
        for k in range(start, stop):
            y = p + k
            one += term
            quadratic += acb(rat(y * y - p * y)) * term
            cubic += acb(rat(y**3 - p * p * y)) * term
            term *= ratio
            ratio *= step
    return one, quadratic, cubic


def sample_row(side: str, fraction_of_wall: Fraction, contract: dict[str, Any]) -> dict[str, Any]:
    p, b, epsilon = side_data(side)
    t = fraction_of_wall * T_MAX
    one_coefficients = (Fraction(1),)
    quadratic_coefficients = (Fraction(0), -p, Fraction(1))
    cubic_coefficients = (Fraction(0), -p * p, Fraction(0), Fraction(1))

    em_one = em_polynomial(p, b, t, one_coefficients)
    em_quadratic = em_polynomial(p, b, t, quadratic_coefficients)
    em_cubic = em_polynomial(p, b, t, cubic_coefficients)
    em_t = 2 * acb.pi() * acb(0, 1) * em_cubic
    em_s = epsilon * (em_one + 4 * acb.pi() * acb(0, 1) * (1 + acb(rat(t))) * em_quadratic)

    direct_one, direct_quadratic, direct_cubic = direct_polynomials(side, t)
    direct_t = 2 * acb.pi() * acb(0, 1) * direct_cubic
    direct_s = epsilon * (
        direct_one + 4 * acb.pi() * acb(0, 1) * (1 + acb(rat(t))) * direct_quadratic
    )
    remainder_t = arb(contract["partial_t_W_remainder_upper"])
    remainder_s = arb(contract["partial_s_W_remainder_upper"])
    enclosure_t = em_t + symmetric_error(remainder_t)
    enclosure_s = em_s + symmetric_error(remainder_s)
    require(enclosure_t.overlaps(direct_t), f"{side} partial_t overlap failed at {fraction_of_wall}")
    require(enclosure_s.overlaps(direct_s), f"{side} partial_s overlap failed at {fraction_of_wall}")

    dt_diagnostic = DIAGNOSTIC_VARIATION_RADIUS / direct_t.abs_upper()
    ds_diagnostic = DIAGNOSTIC_VARIATION_RADIUS / direct_s.abs_upper()
    return {
        "side": side,
        "t_over_first_wall": str(fraction_of_wall),
        "t": str(t),
        "direct_partial_t_W": acb_record(direct_t),
        "EM_partial_t_W": acb_record(em_t),
        "partial_t_difference": acb_record(direct_t - em_t),
        "partial_t_enclosure": acb_record(enclosure_t),
        "direct_partial_s_W": acb_record(direct_s),
        "EM_partial_s_W": acb_record(em_s),
        "partial_s_difference": acb_record(direct_s - em_s),
        "partial_s_enclosure": acb_record(enclosure_s),
        "diagnostic_half_width_for_current_variation_0_001": {
            "t": upper_text(dt_diagnostic),
            "s": upper_text(ds_diagnostic),
            "certified_transport_width": False,
        },
        "both_derivative_enclosures_overlap_direct_sums": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    right = artifact["uniform_derivative_remainder_contract"]["right"]
    left = artifact["uniform_derivative_remainder_contract"]["left"]
    diagnostics = artifact["transport_scale_diagnostic"]
    return f"""# Small-t characteristic Euler--Maclaurin parameter derivatives

Date: 2026-08-25

Status: interval derivative certificate; complete-source transport remains open.

## Exact derivative identities

For

```text
W(p,b,t)=sum_(k=0)^n (p+k)e(t k^2+(b+2pt)k),
```

the exact characteristic derivatives at fixed `s` and fixed `t` are

```text
partial_t W=2 pi i sum (p+k)(k^2+2pk)e(q)
           =2 pi i sum (y^3-p^2 y)e(q),

partial_s W=epsilon[sum e(q)
           +4 pi i(1+t)sum(y^2-py)e(q)],
y=p+k,
epsilon=+1 right, -1 left.                         (PD1)
```

The second identity uses `p_s=epsilon` and `b_s=2epsilon`.  It retains the
amplitude derivative and phase derivative together.

## Uniform Euler--Maclaurin remainders

Applying the Section 11.473 Hermite envelope to each polynomial derivative
through order `2M=160` proves

| side | `partial_t W` remainder | `partial_s W` remainder |
|---|---:|---:|
| right | `< {right['partial_t_W_remainder_upper']}` | `< {right['partial_s_W_remainder_upper']}` |
| left | `< {left['partial_t_W_remainder_upper']}` | `< {left['partial_s_W_remainder_upper']}` |

The exact Fresnel moments through degree three and the Bernoulli endpoint
corrections through `B_160` overlap direct weighted 496284-term sums at
`t/t_* = 0,1/16,1/4,1/2,1` on both sides.  The independent checker uses
512 bits, reverse Hermite/Bernoulli/block order, and reset blocks of length
13 instead of 16.

## Transport-scale diagnostic

Using only the pointwise derivative magnitudes, a current variation budget
of `0.001` would suggest half-widths no larger than

```text
minimum t diagnostic: {diagnostics['minimum_t_half_width']},
minimum s diagnostic: {diagnostics['minimum_s_half_width']}.              (PD2)
```

These are diagnostics, not certified cell widths: they do not include
derivative variation across a cell.  Their purpose is architectural.  A
small value proves that the characteristic current should not be transported
in isolation.  Its complete grouped source derivative must be audited next.
If that derivative remains microscopic, the source and pole-free kernel must
be integrated together before norms instead of enclosed by near-constant
value boxes.

## Pi provenance

Every `pi` in (PD1) and the endpoint formulas comes from differentiating the
declared character `e(v)=exp(2 pi i v)`.  The Euler--Maclaurin remainder
denominator comes from the Fourier series of the periodic Bernoulli
polynomial.  No geometric or fitted occurrence is introduced.

## Boundary

{artifact['proof_boundary']}
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(CHECKER.is_file(), "missing independent checker")
    require(BASE_RESULT.is_file() and BASE_BUILDER.is_file() and BASE_CHECKER.is_file(), "missing base dependency")
    base = json.loads(BASE_RESULT.read_text(encoding="utf-8"))
    require(base.get("passed") is True, "base characteristic gate did not pass")

    flint.ctx.prec = PRECISION_BITS
    contracts = {side: remainder_contract(side) for side in ("right", "left")}
    rows = [sample_row(side, sample, contracts[side]) for side in ("right", "left") for sample in SAMPLES]
    t_widths = [arb(row["diagnostic_half_width_for_current_variation_0_001"]["t"]) for row in rows]
    s_widths = [arb(row["diagnostic_half_width_for_current_variation_0_001"]["s"]) for row in rows]
    diagnostics = {
        "variation_radius": "0.001",
        "minimum_t_half_width": upper_text(min(t_widths)),
        "maximum_t_half_width": upper_text(max(t_widths)),
        "minimum_s_half_width": upper_text(min(s_widths)),
        "maximum_s_half_width": upper_text(max(s_widths)),
        "certified_transport_width": False,
        "decision": "audit the complete grouped source derivative, then preserve the source-kernel product through quadrature if its transport scale remains microscopic",
    }

    proof_boundary = (
        "Exact partial_t W and partial_s W identities at s=1/3, uniform order-80 Euler--Maclaurin "
        "remainder certificates for both derivatives on 0<=t<=1/(2n), and five direct weighted-sum "
        "overlaps per side only. Reported transport half-widths are pointwise scale diagnostics, not "
        "certified cells. No differentiated joined Mordell endpoint or outer source assembly, finite t- "
        "or s-cover, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, "
        "Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved."
    )
    artifact = {
        "kind": STEM,
        "status": "characteristic_t_s_derivatives_uniformly_certified_complete_source_differentiation_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "n": N,
            "s": "1/3",
            "t_interval": ["0", str(T_MAX)],
            "euler_maclaurin_order": EM_ORDER,
            "precision_bits": PRECISION_BITS,
            "workers": 1,
            "direct_phase_reset_block_size": DIRECT_BLOCK_SIZE,
            "daytime_cpu_baseline_samples_percent": CPU_BASELINE_SAMPLES,
            "daytime_cpu_baseline_average_percent": round(sum(CPU_BASELINE_SAMPLES) / len(CPU_BASELINE_SAMPLES), 2),
        },
        "exact_derivative_identities": {
            "partial_t": "2*pi*i*sum_(k=0)^n (p+k)(k^2+2*p*k)e(q)",
            "partial_t_shifted_weight": "2*pi*i*sum (y^3-p^2*y)e(q), y=p+k",
            "partial_s": "epsilon*(sum e(q)+4*pi*i*(1+t)*sum(y^2-p*y)e(q))",
            "epsilon": {"right": 1, "left": -1},
        },
        "uniform_derivative_remainder_contract": contracts,
        "derivative_overlap_matrix": rows,
        "transport_scale_diagnostic": diagnostics,
        "decision": {
            "partial_t_W_uniformly_certified": True,
            "partial_s_W_uniformly_certified": True,
            "all_ten_direct_derivative_overlaps_pass": True,
            "current_level_transport_promoted_to_source_cover": False,
            "complete_source_transport_scale_audit_required": True,
            "physical_quadrature_completed": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Audit the exact complete-source t and s derivatives, both through the original folded source "
            "current and through the grouped recursive representation. If their largest pieces do not cancel "
            "to a macroscopic transport scale, retire value-ball covering and integrate the pole-free "
            "source-kernel product with oscillation preserved before absolute values."
        ),
        "proof_boundary": proof_boundary,
        "primary_source": base["primary_source"],
        "dependencies": {
            "base_result": {"path": relative(BASE_RESULT), "sha256": file_hash(BASE_RESULT)},
            "base_builder": {"path": relative(BASE_BUILDER), "sha256": file_hash(BASE_BUILDER)},
            "base_checker": {"path": relative(BASE_CHECKER), "sha256": file_hash(BASE_CHECKER)},
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
    print("certified characteristic t and s derivatives with ten direct overlaps", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
