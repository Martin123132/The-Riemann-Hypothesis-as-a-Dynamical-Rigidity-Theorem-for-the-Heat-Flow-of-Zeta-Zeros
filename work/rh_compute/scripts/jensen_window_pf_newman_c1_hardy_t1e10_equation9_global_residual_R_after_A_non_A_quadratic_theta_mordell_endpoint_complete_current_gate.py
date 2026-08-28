#!/usr/bin/env python3
"""Differentiate Kuznetsov's Mordell identity as one tied theta current."""

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
from typing import Any, Callable


for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
from mpmath import mp
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_endpoint_complete_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
REFERENCE_BUILDER = BUILDER.with_name(
    "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.py"
)
DEPENDENCIES = {
    "reference_evaluator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.json",
    "poisson_current_self_duality": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_poisson_current_self_duality_gate.json",
}

A = 159_577
B = 5_122_421
L = 2_481_422
K = L - 1
BLOCK = 256
MORDELL_DPS_LOW = 60
MORDELL_DPS_HIGH = 90
SHORT_ABSOLUTE_ALLOWANCE = mp.mpf("1e-42")
FULL_MASS_NORMALIZED_ALLOWANCE = 2.0e-10
TRANSFORMED_MASS_NORMALIZED_ALLOWANCE = 2.0e-10

SHORT_CASES = (
    {"name": "interior_endpoints", "n": 37, "a": "137/17", "tau": "3/20"},
    {"name": "shifted_upper_endpoint", "n": 73, "a": "19/7", "tau": "2/9"},
)

FULL_CASES = (
    {"x": "1/2", "s": "0"},
    {"x": "2/5", "s": "1/3"},
)


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


def load_reference_module():
    spec = importlib.util.spec_from_file_location("non_a_theta_reference", REFERENCE_BUILDER)
    require(spec is not None and spec.loader is not None, "reference module loader failed")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def nearest_integer_shift(value: Fraction) -> tuple[int, Fraction]:
    shift = floor_fraction(value + Fraction(1, 2))
    reduced = value - shift
    require(Fraction(-1, 2) <= reduced < Fraction(1, 2), "linear phase normalization failed")
    return shift, reduced


def mp_fraction(value: Fraction) -> Any:
    return mp.mpf(value.numerator) / value.denominator


def mp_cis_pi_fraction(value: Fraction) -> Any:
    reduced = value % 2
    angle = mp.pi * mp_fraction(reduced)
    return mp.mpc(mp.cos(angle), mp.sin(angle))


def mp_record(value: Any, digits: int = 70) -> dict[str, str]:
    return {"real": mp.nstr(mp.re(value), digits), "imag": mp.nstr(mp.im(value), digits)}


def complex_record(value: complex) -> dict[str, float]:
    return {"real": value.real, "imag": value.imag}


def symbolic_certificate() -> dict[str, str]:
    a, z, tau, r, k = sp.symbols("a z tau r k", real=True)
    g0, gw = sp.symbols("G G_w")
    pi = sp.pi
    ii = sp.I
    c_factor = sp.exp(ii * pi / 4 - ii * pi * z**2 / (2 * tau)) / sp.sqrt(2 * tau)
    c_derivative = sp.diff(c_factor, z)
    main_current = a * c_factor * g0 / tau + (c_derivative * g0 + c_factor * gw / (2 * tau)) / (pi * ii)
    expected_main = c_factor * ((a - z) * g0 + gw / (2 * pi * ii)) / tau
    require(sp.simplify(main_current - expected_main) == 0, "tied main-current cancellation failed")

    v = r + k
    phase_left = (a**2 - z**2) / (2 * tau) + z * k / tau - k**2 / (2 * tau)
    phase_right = a * v / tau - v**2 / (2 * tau)
    require(sp.simplify((phase_left - phase_right).subs(a, z + r)) == 0, "absolute-dual phase failed")

    n = sp.symbols("n", integer=True, nonnegative=True)
    h0 = sp.Function("h_0")(z)
    e0 = sp.exp(-pi * ii * (z - tau / 2))
    lower = -ii * e0 * h0 / 2
    lower_joined = a * lower / tau + sp.diff(lower, z) / (pi * ii)
    lower_expected = -ii * e0 * ((a / tau - 1) * h0 + sp.diff(h0, z) / (pi * ii)) / 2
    require(sp.simplify(lower_joined - lower_expected) == 0, "lower endpoint derivative failed")

    h1 = sp.Function("h_1")(z)
    e1 = sp.exp(2 * pi * ii * (n + sp.Rational(1, 2)) * (z + tau * (n + sp.Rational(1, 2))))
    upper = -ii * e1 * h1 / 2
    upper_joined = a * upper / tau + sp.diff(upper, z) / (pi * ii)
    upper_expected = -ii * e1 * (
        (a / tau + 2 * n + 1) * h1 + sp.diff(h1, z) / (pi * ii)
    ) / 2
    require(sp.simplify(upper_joined - upper_expected) == 0, "upper endpoint derivative failed")

    return {
        "source_operator": "D_a=(a/tau)+(pi*i)^(-1)*partial_z",
        "normalized_shift": "a=z+r, r in Z, -1/2<=z<1/2",
        "main_cancellation": "D_a[C(z,tau)G(z/(2tau))]=C/tau*[rG+G_w/(2pi*i)]",
        "absolute_dual_moment": "rG+G_w/(2pi*i)=sum_(k=0)^m (r+k)exp(2pi*i*(wk+sigma*k^2))",
        "dual_phase": "(a^2-z^2)/(2tau)+zk/tau-k^2/(2tau)=(a/tau)(r+k)-(r+k)^2/(2tau)",
        "lower_joined_endpoint": "-i*E_-/2*[(a/tau-1)h_-+h_-prime/(pi*i)]",
        "upper_joined_endpoint": "-i*(-1)^m*E_+/2*[(a/tau+2n+1)h_++h_+prime/(pi*i)]",
    }


def central_mordell_integrals(z: Fraction, tau: Fraction, dps: int) -> tuple[Any, Any, float]:
    """Evaluate h(z,tau) and its z derivative on the central strip."""

    require(tau > 0 and Fraction(-1, 2) <= z < Fraction(1, 2), "central Mordell input out of range")
    mp.dps = dps
    zz = mp_fraction(z)
    tt = mp_fraction(tau)
    theta = mp.e ** (mp.pi * 1j / 4)
    cutoff = mp.sqrt(mp.mpf(dps + 25) * mp.log(10) / (mp.pi * tt))
    integer_cutoff = max(2, int(mp.ceil(cutoff)))
    points = [mp.mpf(index) for index in range(integer_cutoff + 1)]

    def common(y: Any) -> Any:
        return mp.e ** (-mp.pi * tt * y * y) / mp.cosh(mp.pi * theta * y)

    def h_integrand(y: Any) -> Any:
        return common(y) * mp.cosh(2 * mp.pi * zz * theta * y)

    def hp_integrand(y: Any) -> Any:
        return y * common(y) * mp.sinh(2 * mp.pi * zz * theta * y)

    h_value = 2 * theta * mp.quad(h_integrand, points)
    hp_value = 4 * mp.pi * 1j * mp.quad(hp_integrand, points)
    return h_value, hp_value, float(cutoff)


def mordell_h_and_derivative(z: Fraction, signed_tau: Fraction, dps: int) -> tuple[Any, Any, dict[str, Any]]:
    """Reduce z by Mordell's unit recurrence, then evaluate h and h'."""

    require(signed_tau != 0, "Mordell parameter must be nonzero")
    tau = abs(signed_tau)
    shift, central = nearest_integer_shift(z)
    h_value, hp_value, cutoff = central_mordell_integrals(central, tau, dps)
    y = central

    if shift > 0:
        for _ in range(shift):
            g = 2 / mp.sqrt(mp_fraction(tau)) * mp_cis_pi_fraction(
                Fraction(1, 4) + (y + Fraction(1, 2)) ** 2 / tau
            )
            gp = g * 2 * mp.pi * 1j * mp_fraction(y + Fraction(1, 2)) / mp_fraction(tau)
            h_value, hp_value = g - h_value, gp - hp_value
            y += 1
    elif shift < 0:
        for _ in range(-shift):
            y_prev = y - 1
            g = 2 / mp.sqrt(mp_fraction(tau)) * mp_cis_pi_fraction(
                Fraction(1, 4) + (y_prev + Fraction(1, 2)) ** 2 / tau
            )
            gp = g * 2 * mp.pi * 1j * mp_fraction(y_prev + Fraction(1, 2)) / mp_fraction(tau)
            h_value, hp_value = g - h_value, gp - hp_value
            y = y_prev

    require(y == z, "Mordell recurrence did not reach requested argument")
    if signed_tau < 0:
        h_value = mp.conj(h_value)
        hp_value = mp.conj(hp_value)
    return h_value, hp_value, {
        "argument": str(z),
        "signed_tau": str(signed_tau),
        "central_argument": str(central),
        "unit_shift_count": shift,
        "dps": dps,
        "finite_cutoff": cutoff,
        "quadrature_is_diagnostic": True,
    }


def stable_mordell_pair(z: Fraction, signed_tau: Fraction) -> tuple[Any, Any, dict[str, Any]]:
    low_h, low_hp, low_meta = mordell_h_and_derivative(z, signed_tau, MORDELL_DPS_LOW)
    high_h, high_hp, high_meta = mordell_h_and_derivative(z, signed_tau, MORDELL_DPS_HIGH)
    mp.dps = MORDELL_DPS_HIGH
    h_delta = abs(high_h - low_h)
    hp_delta = abs(high_hp - low_hp)
    require(h_delta < mp.mpf("1e-45") and hp_delta < mp.mpf("1e-43"), "Mordell cross-precision drift")
    return high_h, high_hp, {
        "low_precision": low_meta,
        "high_precision": high_meta,
        "h_cross_precision_absolute": mp.nstr(h_delta, 20),
        "h_prime_cross_precision_absolute": mp.nstr(hp_delta, 20),
    }


def minimal_theta_period(w: Fraction, sigma: Fraction, limit: int = 131_072) -> int:
    for period in range(1, limit + 1):
        slope = 2 * sigma * period
        constant = w * period + sigma * period * period
        if slope.denominator == 1 and constant.denominator == 1:
            return period
    raise RuntimeError("transformed theta period not found")


def transformed_weighted_recurrence(
    reference: Any,
    w: Fraction,
    sigma: Fraction,
    m: int,
    shift: int,
) -> complex:
    q = reference.cis_pi_fraction(4 * sigma)
    block_real: list[float] = []
    block_imag: list[float] = []
    for start in range(0, m + 1, BLOCK):
        stop = min(m + 1, start + BLOCK)
        term = reference.cis_pi_fraction(2 * (w * start + sigma * start * start))
        ratio = reference.cis_pi_fraction(2 * (w + sigma * (2 * start + 1)))
        local_real: list[float] = []
        local_imag: list[float] = []
        for index in range(start, stop):
            weight = shift + index
            local_real.append(weight * term.real)
            local_imag.append(weight * term.imag)
            term *= ratio
            ratio *= q
        block_real.append(math.fsum(local_real))
        block_imag.append(math.fsum(local_imag))
    return complex(math.fsum(block_real), math.fsum(block_imag))


def transformed_weighted_periodic(
    reference: Any,
    w: Fraction,
    sigma: Fraction,
    m: int,
    shift: int,
    period: int,
) -> Any:
    mp.dps = MORDELL_DPS_HIGH
    cycle = [reference.mp_cis_pi_fraction(2 * (w * index + sigma * index * index)) for index in range(period)]
    z0 = mp.fsum(cycle)
    z1 = mp.fsum([index * cycle[index] for index in range(period)])
    cycles, remainder = divmod(m + 1, period)
    r0 = mp.fsum(cycle[:remainder])
    r1 = mp.fsum([index * cycle[index] for index in range(remainder)])
    s0 = cycles * z0 + r0
    s1 = period * cycles * (cycles - 1) * z0 / 2 + cycles * z1 + cycles * period * r0 + r1
    return shift * s0 + s1


def transformed_weight_mass(m: int, shift: int) -> float:
    require(shift >= 0, "pilot transformed weights must be nonnegative")
    return float((m + 1) * shift + m * (m + 1) // 2)


def direct_short_current(a: Fraction, tau: Fraction, n: int) -> Any:
    mp.dps = MORDELL_DPS_HIGH
    current = mp.mpc(0)
    c = a / tau
    for index in range(n + 1):
        current += mp_fraction(c + 2 * index) * mp_cis_pi_fraction(2 * (a * index + tau * index * index))
    return mp_cis_pi_fraction(a * a / (2 * tau)) * current


def one_step_current(
    a: Fraction,
    tau: Fraction,
    n: int,
    transformed_sum: Callable[[Fraction, Fraction, int, int], Any],
) -> tuple[Any, dict[str, Any]]:
    require(tau > 0, "one-step pilot requires positive tau")
    shift, z = nearest_integer_shift(a)
    require(shift >= 0, "pilot requires a nonnegative absolute-frequency shift")
    m = floor_fraction(2 * n * tau)
    w = z / (2 * tau)
    sigma = -Fraction(1, 1) / (4 * tau)
    weighted = transformed_sum(w, sigma, m, shift)

    lower_argument = z - tau + Fraction(1, 2)
    upper_argument = z + (2 * n + 1) * tau - m - Fraction(1, 2)
    lower_h, lower_hp, lower_meta = stable_mordell_pair(lower_argument, -2 * tau)
    upper_h, upper_hp, upper_meta = stable_mordell_pair(upper_argument, -2 * tau)
    mp.dps = MORDELL_DPS_HIGH

    main_phase = Fraction(1, 4) + (a * a - z * z) / (2 * tau)
    main_prefactor = mp_cis_pi_fraction(main_phase) / (mp_fraction(tau) * mp.sqrt(2 * mp_fraction(tau)))
    main = main_prefactor * weighted

    lower_phase = -(z - tau / 2)
    upper_phase = 2 * (Fraction(n) + Fraction(1, 2)) * (z + tau * (Fraction(n) + Fraction(1, 2)))
    lower_joined = (
        -1j
        * mp_cis_pi_fraction(lower_phase)
        * ((mp_fraction(a / tau) - 1) * lower_h + lower_hp / (mp.pi * 1j))
        / 2
    )
    upper_joined = (
        -1j
        * ((-1) ** m)
        * mp_cis_pi_fraction(upper_phase)
        * ((mp_fraction(a / tau) + 2 * n + 1) * upper_h + upper_hp / (mp.pi * 1j))
        / 2
    )
    outer = mp_cis_pi_fraction(a * a / (2 * tau))
    lower_physical = outer * lower_joined
    upper_physical = outer * upper_joined
    remainder = lower_physical + upper_physical
    complete = main + remainder
    return complete, {
        "n": n,
        "a": str(a),
        "tau": str(tau),
        "integer_shift_r": shift,
        "normalized_z": str(z),
        "transformed_m": m,
        "transformed_w": str(w),
        "transformed_sigma": str(sigma),
        "lower_argument": str(lower_argument),
        "upper_argument": str(upper_argument),
        "main_current": mp_record(main),
        "lower_endpoint_current": mp_record(lower_physical),
        "upper_endpoint_current": mp_record(upper_physical),
        "joined_endpoint_current": mp_record(remainder),
        "complete_current": mp_record(complete),
        "main_plus_endpoint_condition_ratio": float((abs(main) + abs(remainder)) / max(mp.mpf(1), abs(complete))),
        "lower_mordell": lower_meta,
        "upper_mordell": upper_meta,
    }


def run_short_rows(reference: Any) -> list[dict[str, Any]]:
    rows = []
    for case in SHORT_CASES:
        a = Fraction(case["a"])
        tau = Fraction(case["tau"])
        n = int(case["n"])

        def transformed(w: Fraction, sigma: Fraction, m: int, shift: int) -> Any:
            mp.dps = MORDELL_DPS_HIGH
            return mp.fsum(
                [
                    (shift + index) * reference.mp_cis_pi_fraction(2 * (w * index + sigma * index * index))
                    for index in range(m + 1)
                ]
            )

        complete, metadata = one_step_current(a, tau, n, transformed)
        direct = direct_short_current(a, tau, n)
        error = abs(complete - direct)
        require(error <= SHORT_ABSOLUTE_ALLOWANCE, f"short Mordell identity mismatch: {case}")
        rows.append(
            {
                **case,
                **metadata,
                "direct_current": mp_record(direct),
                "absolute_error": mp.nstr(error, 30),
                "allowance": mp.nstr(SHORT_ABSOLUTE_ALLOWANCE, 10),
                "passed": True,
            }
        )
    return rows


def run_full_rows(reference: Any) -> list[dict[str, Any]]:
    rows = []
    for case in FULL_CASES:
        x = Fraction(case["x"])
        s = Fraction(case["s"])
        tau = x / 2
        a = x * (A + 2 * s) / 2
        shift, z = nearest_integer_shift(a)
        m = floor_fraction(2 * K * tau)
        w = z / (2 * tau)
        sigma = -Fraction(1, 1) / (4 * tau)
        period = minimal_theta_period(w, sigma)
        transformed_observed = transformed_weighted_recurrence(reference, w, sigma, m, shift)
        transformed_expected = transformed_weighted_periodic(reference, w, sigma, m, shift, period)
        transformed_mass = transformed_weight_mass(m, shift)
        transformed_error = float(abs(mp.mpc(transformed_observed) - transformed_expected)) / transformed_mass
        require(
            transformed_error <= TRANSFORMED_MASS_NORMALIZED_ALLOWANCE,
            f"transformed weighted recurrence mismatch: {case}",
        )

        def transformed(_w: Fraction, _sigma: Fraction, _m: int, _shift: int) -> complex:
            require((_w, _sigma, _m, _shift) == (w, sigma, m, shift), "transformed callback drift")
            return transformed_observed

        complete, metadata = one_step_current(a, tau, K, transformed)
        reference_current = reference.recurrence_theta_current(x, s, L)[2]
        source_mass = float(Fraction(L) * (A + 2 * s + L - 1))
        mass_error = float(abs(complete - reference_current)) / source_mass
        require(mass_error <= FULL_MASS_NORMALIZED_ALLOWANCE, f"full Mordell one-step mismatch: {case}")
        rows.append(
            {
                **case,
                **metadata,
                "minimal_transformed_period": period,
                "transformed_term_count": m + 1,
                "transformed_compression_ratio": (m + 1) / L,
                "transformed_weighted_recurrence": complex_record(transformed_observed),
                "transformed_weighted_periodic": mp_record(transformed_expected),
                "transformed_mass_normalized_error": transformed_error,
                "reference_current": complex_record(reference_current),
                "source_absolute_term_mass": source_mass,
                "complete_mass_normalized_error": mass_error,
                "diagnostic_only": True,
                "passed": True,
            }
        )
    return rows


def render_note(artifact: dict[str, Any]) -> str:
    summary = artifact["summary"]
    return f"""# Endpoint-complete Mordell one-step theta current

Date: 2026-08-24

Status: exact endpoint-complete one-step current identity certified;
numerical Mordell validation diagnostic; recursive interval evaluator open

## Scope

Let

```text
F_n(a,tau)=sum_(k=0)^n exp(2*pi*i*(a*k+tau*k^2)),
a=tau*c,                         tau=x/2>0,
J_n=exp(i*pi*a^2/(2*tau))*[c F_n+(pi*i)^(-1)partial_a F_n].   (ME1)
```

Kuznetsov's exact Theorem 1 in
[Computing the truncated theta function via Mordell integral](https://arxiv.org/abs/1306.4081v2)
writes `F_n` as one shorter theta sum plus two explicit Mordell-integral
endpoint terms.  This gate differentiates that identity in the one joined
direction required by the non-A source.  No unnamed constants or code from
the paper are imported.

## Exact joined transform

Choose the unique integer `r` and reduced phase `z` with

```text
a=z+r,                         -1/2<=z<1/2,
m=floor(2*n*tau),              w=z/(2*tau),
sigma=-1/(4*tau),
C=exp(i*pi/4-i*pi*z^2/(2*tau))/sqrt(2*tau).          (ME2)
```

Differentiating the exact Mordell identity and keeping the physical current
tied gives

```text
[a/tau+(pi*i)^(-1)partial_z] [C F_m(w,sigma)]
 =C/tau [r F_m(w,sigma)+(2*pi*i)^(-1)partial_w F_m(w,sigma)]
 =C/tau sum_(k=0)^m (r+k)exp(2*pi*i*(w*k+sigma*k^2)). (ME3)
```

Thus the apparent zeroth and first transformed moments are one absolute-
frequency first moment.  With `v=r+k`, the phase identity

```text
(a^2-z^2)/(2*tau)+z*k/tau-k^2/(2*tau)
 =(a/tau)*v-v^2/(2*tau)                              (ME4)
```

recovers the pure dual current of Section 11.462 without evaluating a huge
unreduced Mordell argument.

Write

```text
u_-=z-tau+1/2,
u_+=z+(2n+1)tau-m-1/2,
E_-=exp[-pi*i*(z-tau/2)],
E_+=exp[2*pi*i*(n+1/2)*(z+tau*(n+1/2))].             (ME5)
```

The endpoint-complete joined current is exactly

```text
J_endpoint=exp(i*pi*a^2/(2*tau))*(-i/2)*{{
 E_-[(a/tau-1)h(u_-,-2tau)+h_z(u_-,-2tau)/(pi*i)]
 +(-1)^m E_+[(a/tau+2n+1)h(u_+,-2tau)
              +h_z(u_+,-2tau)/(pi*i)]}}.            (ME6)
```

Equations (ME3) and (ME6) are one exact modular step for the full current;
they do not split `R_0` from `R_1`.

## Validation

Two short rows compare (ME3)+(ME6) with independent 90-digit direct sums.
Their largest absolute discrepancy is
`{summary['maximum_short_absolute_error']}`.  Two full
`L={L}` rows compare the complete one-step expression with the validated
O(`L`) source oracle.  The transformed weighted sums are also checked
against independent exact rational-period compression, with periods
{summary['minimum_transformed_period']} through
{summary['maximum_transformed_period']}.  The largest discrepancies relative
to absolute term mass are

```text
transformed weighted sum: {summary['maximum_transformed_mass_normalized_error']:.6e},
complete source current:  {summary['maximum_complete_mass_normalized_error']:.6e}.  (ME7)
```

The Mordell values and derivatives are evaluated twice, at 60 and 90 digits,
from the rotated defining integral after exact unit-recurrence reduction.
The largest cross-precision changes are

```text
h:   {summary['maximum_h_cross_precision_absolute']},
h_z: {summary['maximum_h_prime_cross_precision_absolute']}.                 (ME8)
```

These numerical rows validate the specialization and implementation only.
They are not interval enclosures.  Kuznetsov explicitly distinguishes his
rigorously analyzable algorithm from the faster practical Gauss-Laguerre
version and notes that he did not obtain rigorous error estimates for that
practical quadrature.  This gate therefore does not promote (ME7)--(ME8) to
a theorem.

## Decision

The previously open endpoint ownership problem for one modular step is now
closed at the exact-identity level: the remainder is the joined pair (ME6),
not separately estimated `R_0,R_1`.  The next obligation is a certified
evaluator for `h,h_z` and a small-`tau` Euler--Maclaurin branch, followed by a
recursive interval implementation.  Only then can this become a uniform
fast evaluator for the physical `(x,s)` quadrature.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (ME1)--(ME8) comes from the inherited quadratic
Fourier phase or Kuznetsov's stated Mordell identity.  No fitted or geometric
constant is inserted.

This gate proves the exact endpoint-complete one-step current identity and
the rational transformed-period formulas only.  The four floating rows are
diagnostics.  It does not prove a certified Mordell evaluator, the small-
`tau` branch, recursive interval control, physical quadrature, the non-A
bound, joined `R_after_A`, `R_Dir`, `Q_K-T`, an all-height theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(CHECKER.is_file() and REFERENCE_BUILDER.is_file(), "checker or reference builder missing")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "Mordell dependency not passed")
    require(
        dependencies["reference_evaluator"]["scope"]["source_cell_count"] == L,
        "reference roster drift",
    )
    require(
        dependencies["poisson_current_self_duality"]["decision"][
            "completed_saddle_current_is_pure_absolute_dual_first_moment"
        ]
        is True,
        "self-duality dependency drift",
    )

    reference = load_reference_module()
    symbolic = symbolic_certificate()
    short_rows = run_short_rows(reference)
    full_rows = run_full_rows(reference)
    all_mordell = []
    for row in short_rows + full_rows:
        all_mordell.extend([row["lower_mordell"], row["upper_mordell"]])
    summary = {
        "short_row_count": len(short_rows),
        "full_row_count": len(full_rows),
        "maximum_short_absolute_error": max(mp.mpf(row["absolute_error"]) for row in short_rows),
        "minimum_transformed_period": min(row["minimal_transformed_period"] for row in full_rows),
        "maximum_transformed_period": max(row["minimal_transformed_period"] for row in full_rows),
        "maximum_transformed_mass_normalized_error": max(
            row["transformed_mass_normalized_error"] for row in full_rows
        ),
        "maximum_complete_mass_normalized_error": max(row["complete_mass_normalized_error"] for row in full_rows),
        "maximum_h_cross_precision_absolute": max(
            mp.mpf(item["h_cross_precision_absolute"]) for item in all_mordell
        ),
        "maximum_h_prime_cross_precision_absolute": max(
            mp.mpf(item["h_prime_cross_precision_absolute"]) for item in all_mordell
        ),
    }
    summary_json = {
        key: (mp.nstr(value, 30) if hasattr(value, "ae") or type(value).__module__.startswith("mpmath") else value)
        for key, value in summary.items()
    }

    artifact = {
        "kind": STEM,
        "status": "exact_endpoint_complete_Mordell_one_step_current_identity_certified_recursive_interval_evaluator_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "B": B,
            "L": L,
            "K": K,
            "physical_x": "0<x<=1/2",
            "physical_s": "0<=s<=1",
            "workers": 1,
        },
        "primary_source": {
            "author": "Alexey Kuznetsov",
            "title": "Computing the truncated theta function via Mordell integral",
            "arxiv": "1306.4081v2",
            "url": "https://arxiv.org/abs/1306.4081v2",
            "used_for": "Theorem 1 exact truncated-theta identity and its stated derivative extension",
            "imported_numerical_constants": False,
            "practical_quadrature_rigorous_error_bound_claimed_by_source": False,
        },
        "symbolic_certificate": symbolic,
        "short_high_precision_rows": short_rows,
        "full_one_step_rows": full_rows,
        "summary": summary_json,
        "decision": {
            "Kuznetsov_Theorem_1_specialized_to_non_A_current": True,
            "normalized_shift_preserves_absolute_dual_first_moment": True,
            "transformed_zeroth_and_first_moments_joined_exactly": True,
            "two_Mordell_endpoint_terms_joined_exactly_with_derivatives": True,
            "endpoint_complete_one_step_identity_built": True,
            "full_length_one_step_cross_checked_against_reference_oracle": True,
            "Mordell_quadrature_rows_are_diagnostic_only": True,
            "certified_Mordell_h_and_h_prime_evaluator_built": False,
            "small_tau_Euler_Maclaurin_branch_built": False,
            "recursive_interval_theta_evaluator_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Replace the diagnostic rotated-integral quadrature for h and h_z by a certified evaluator with explicit tails and rounding enclosures, add the small-tau Euler--Maclaurin branch, and validate a recursive interval current without separating the joined endpoint terms.",
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "blas_threads": 1,
            "process_priority": priority,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "reference_builder": {"path": relative(REFERENCE_BUILDER), "sha256": file_hash(REFERENCE_BUILDER)},
        },
        "proof_boundary": "Exact endpoint-complete one-step Mordell current identity and exact rational transformed-period compression only. The short and full numerical rows use non-interval high-precision quadrature and are diagnostics. No certified Mordell h,h_z evaluator, small-tau Euler--Maclaurin branch, recursive interval theta evaluator, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"certified endpoint-complete Mordell one-step current identity; "
        f"{len(short_rows)} short and {len(full_rows)} full diagnostic rows",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
