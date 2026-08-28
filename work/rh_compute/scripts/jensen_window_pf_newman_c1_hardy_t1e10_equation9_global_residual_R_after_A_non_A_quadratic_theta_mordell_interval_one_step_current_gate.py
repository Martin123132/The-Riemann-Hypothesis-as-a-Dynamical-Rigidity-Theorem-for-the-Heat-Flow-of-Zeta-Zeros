#!/usr/bin/env python3
"""Certify two full Mordell one-step currents with Arb/ACB integration."""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import hashlib
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
from mpmath import mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_interval_one_step_current_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "Mordell_exact_one_step": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_mordell_endpoint_complete_current_gate.json",
    "reference_evaluator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_non_A_quadratic_theta_reference_evaluator_gate.json",
}

A = 159_577
B = 5_122_421
L = 2_481_422
K = L - 1


@dataclass(frozen=True)
class Settings:
    precision_bits: int
    cutoff: int
    abs_tol: str
    deg_limit: int = 64
    eval_limit: int = 500_000
    depth_limit: int = 64


PRODUCTION = Settings(precision_bits=320, cutoff=10, abs_tol="1e-55")
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


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def is_even_integer(value: Fraction) -> bool:
    return value.denominator == 1 and value.numerator % 2 == 0


def nearest_integer_shift(value: Fraction) -> tuple[int, Fraction]:
    shift = floor_fraction(value + Fraction(1, 2))
    reduced = value - shift
    require(Fraction(-1, 2) <= reduced < Fraction(1, 2), "linear phase normalization failed")
    return shift, reduced


def arb_rational(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def acb_cis_pi(value: Fraction) -> acb:
    reduced = value % 2
    return (acb(0, 1) * acb.pi() * acb(arb_rational(reduced))).exp()


def arb_upper_text(value: arb, digits: int = 40) -> str:
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


def central_tail_bounds(tau: Fraction, cutoff: int) -> tuple[arb, arb, dict[str, str]]:
    """Uniform tails for |z|<=1/2 in the rotated Mordell integral."""

    require(tau > 0 and cutoff > 0, "invalid Mordell tail inputs")
    pi = arb.pi()
    tt = arb_rational(tau)
    yy = arb(cutoff)
    u = pi * yy / arb(2).sqrt()
    coth_u = u.coth()
    h_tail = coth_u / tt.sqrt() * ((pi * tt).sqrt() * yy).erfc()
    hp_tail = 2 * coth_u / tt * (-pi * tt * yy * yy).exp()
    return h_tail, hp_tail, {
        "u=pi*Y/sqrt(2)": str(u),
        "coth_u": str(coth_u),
        "h_tail_modulus_upper": arb_upper_text(h_tail),
        "h_prime_tail_modulus_upper": arb_upper_text(hp_tail),
    }


def central_mordell_balls(
    z: Fraction,
    tau: Fraction,
    settings: Settings,
) -> tuple[acb, acb, dict[str, Any]]:
    """Rigorous finite ACB integrals plus explicit analytic tails."""

    require(tau > 0, "central Mordell tau must be positive")
    require(Fraction(-1, 2) <= z < Fraction(1, 2), "central Mordell z out of range")
    flint.ctx.prec = settings.precision_bits
    pi = acb.pi()
    theta = (acb(0, 1) * pi / 4).exp()
    zz = acb(arb_rational(z))
    tt = acb(arb_rational(tau))
    left = acb(0)
    right = acb(settings.cutoff)
    tolerance = arb(settings.abs_tol)

    def common(y: acb) -> acb:
        return (-pi * tt * y * y).exp() / (pi * theta * y).cosh()

    def h_integrand(y: acb, analytic: bool) -> acb:
        return common(y) * (2 * pi * zz * theta * y).cosh()

    def hp_integrand(y: acb, analytic: bool) -> acb:
        return y * common(y) * (2 * pi * zz * theta * y).sinh()

    started = time.time()
    h_finite = 2 * theta * acb.integral(
        h_integrand,
        left,
        right,
        abs_tol=tolerance,
        rel_tol=tolerance,
        deg_limit=settings.deg_limit,
        eval_limit=settings.eval_limit,
        depth_limit=settings.depth_limit,
    )
    hp_finite = 4 * pi * acb(0, 1) * acb.integral(
        hp_integrand,
        left,
        right,
        abs_tol=tolerance,
        rel_tol=tolerance,
        deg_limit=settings.deg_limit,
        eval_limit=settings.eval_limit,
        depth_limit=settings.depth_limit,
    )
    elapsed = time.time() - started
    require(h_finite.is_finite() and hp_finite.is_finite(), "ACB Mordell integration failed")
    h_tail, hp_tail, tail_record = central_tail_bounds(tau, settings.cutoff)
    h_ball = h_finite + symmetric_complex_error(h_tail)
    hp_ball = hp_finite + symmetric_complex_error(hp_tail)
    return h_ball, hp_ball, {
        "central_z": str(z),
        "positive_tau": str(tau),
        "precision_bits": settings.precision_bits,
        "cutoff": settings.cutoff,
        "abs_tol": settings.abs_tol,
        "finite_h": acb_record(h_finite),
        "finite_h_prime": acb_record(hp_finite),
        "tail": tail_record,
        "enclosed_h": acb_record(h_ball),
        "enclosed_h_prime": acb_record(hp_ball),
        "elapsed_seconds": round(elapsed, 6),
    }


def mordell_h_and_derivative_balls(
    z: Fraction,
    signed_tau: Fraction,
    settings: Settings,
) -> tuple[acb, acb, dict[str, Any]]:
    """Apply the exact unit recurrence to a certified central enclosure."""

    require(signed_tau != 0, "Mordell parameter must be nonzero")
    tau = abs(signed_tau)
    shift, central = nearest_integer_shift(z)
    h_ball, hp_ball, central_record = central_mordell_balls(central, tau, settings)
    y = central
    recurrence_rows = []
    if shift > 0:
        for _ in range(shift):
            phase = Fraction(1, 4) + (y + Fraction(1, 2)) ** 2 / tau
            g = 2 * acb_cis_pi(phase) / acb(arb_rational(tau)).sqrt()
            gp = g * 2 * acb.pi() * acb(0, 1) * acb(arb_rational(y + Fraction(1, 2))) / acb(
                arb_rational(tau)
            )
            h_ball, hp_ball = g - h_ball, gp - hp_ball
            recurrence_rows.append({"from": str(y), "to": str(y + 1), "phase": str(phase)})
            y += 1
    elif shift < 0:
        for _ in range(-shift):
            y_previous = y - 1
            phase = Fraction(1, 4) + (y_previous + Fraction(1, 2)) ** 2 / tau
            g = 2 * acb_cis_pi(phase) / acb(arb_rational(tau)).sqrt()
            gp = g * 2 * acb.pi() * acb(0, 1) * acb(
                arb_rational(y_previous + Fraction(1, 2))
            ) / acb(arb_rational(tau))
            h_ball, hp_ball = g - h_ball, gp - hp_ball
            recurrence_rows.append({"from": str(y), "to": str(y_previous), "phase": str(phase)})
            y = y_previous
    require(y == z, "Mordell unit recurrence did not reach target")
    if signed_tau < 0:
        h_ball = h_ball.conjugate()
        hp_ball = hp_ball.conjugate()
    return h_ball, hp_ball, {
        "target_z": str(z),
        "signed_tau": str(signed_tau),
        "central_z": str(central),
        "unit_shift_count": shift,
        "conjugated_for_negative_tau": signed_tau < 0,
        "central_enclosure": central_record,
        "unit_recurrence": recurrence_rows,
        "final_h": acb_record(h_ball),
        "final_h_prime": acb_record(hp_ball),
    }


def minimal_transformed_period(w: Fraction, sigma: Fraction, limit: int = 131_072) -> int:
    for period in range(1, limit + 1):
        if (2 * sigma * period).denominator == 1 and (
            w * period + sigma * period * period
        ).denominator == 1:
            return period
    raise RuntimeError("transformed theta period not found")


def transformed_weighted_ball(
    w: Fraction,
    sigma: Fraction,
    m: int,
    shift: int,
) -> tuple[acb, int]:
    period = minimal_transformed_period(w, sigma)
    cycle = [acb_cis_pi(2 * (w * index + sigma * index * index)) for index in range(period)]
    z0 = sum(cycle, acb(0))
    z1 = sum((index * cycle[index] for index in range(period)), acb(0))
    cycles, remainder = divmod(m + 1, period)
    r0 = sum(cycle[:remainder], acb(0))
    r1 = sum((index * cycle[index] for index in range(remainder)), acb(0))
    s0 = cycles * z0 + r0
    s1 = period * cycles * (cycles - 1) * z0 / 2 + cycles * z1 + cycles * period * r0 + r1
    return shift * s0 + s1, period


def minimal_source_period(x: Fraction, s: Fraction, limit: int = 131_072) -> int:
    if x == 0:
        return 1
    for period in range(1, limit + 1):
        xp = x * period
        if xp.denominator == 1 and is_even_integer(xp * (period + A + 2 * s)):
            return period
    raise RuntimeError("source theta period not found")


def source_current_ball(x: Fraction, s: Fraction) -> tuple[acb, int]:
    period = minimal_source_period(x, s)
    cycle: list[acb] = []
    alpha: list[Fraction] = []
    for index in range(period):
        alpha_index = A + 2 * index + 2 * s
        alpha.append(alpha_index)
        cycle.append(acb_cis_pi(x * alpha_index**2 / 4))
    w0 = sum(cycle, acb(0))
    wa = sum((acb(arb_rational(alpha[index])) * cycle[index] for index in range(period)), acb(0))
    cycles, remainder = divmod(L, period)
    wr0 = sum(cycle[:remainder], acb(0))
    wra = sum(
        (acb(arb_rational(alpha[index])) * cycle[index] for index in range(remainder)),
        acb(0),
    )
    current = cycles * wa + period * cycles * (cycles - 1) * w0
    current += 2 * cycles * period * wr0 + wra
    return current, period


def compute_full_row(x: Fraction, s: Fraction, settings: Settings) -> dict[str, Any]:
    flint.ctx.prec = settings.precision_bits
    tau = x / 2
    c = A + 2 * s
    a = x * c / 2
    shift, z = nearest_integer_shift(a)
    require(shift >= 0, "physical transformed shift must be nonnegative")
    m = floor_fraction(2 * K * tau)
    w = z / (2 * tau)
    sigma = -Fraction(1, 1) / (4 * tau)
    transformed, transformed_period = transformed_weighted_ball(w, sigma, m, shift)

    cache: dict[tuple[Fraction, Fraction], tuple[acb, acb, dict[str, Any]]] = {}

    def get_mordell(argument: Fraction, signed_tau: Fraction) -> tuple[acb, acb, dict[str, Any]]:
        key = (argument, signed_tau)
        if key not in cache:
            cache[key] = mordell_h_and_derivative_balls(argument, signed_tau, settings)
        return cache[key]

    lower_argument = z - tau + Fraction(1, 2)
    upper_argument = z + (2 * K + 1) * tau - m - Fraction(1, 2)
    lower_h, lower_hp, lower_record = get_mordell(lower_argument, -2 * tau)
    upper_h, upper_hp, upper_record = get_mordell(upper_argument, -2 * tau)

    main_phase = Fraction(1, 4) + (a * a - z * z) / (2 * tau)
    main = acb_cis_pi(main_phase) * transformed / (
        acb(arb_rational(tau)) * acb(2 * arb_rational(tau)).sqrt()
    )
    lower_phase = -(z - tau / 2)
    upper_phase = 2 * (Fraction(K) + Fraction(1, 2)) * (
        z + tau * (Fraction(K) + Fraction(1, 2))
    )
    pi_i = acb.pi() * acb(0, 1)
    lower_joined = (
        -acb(0, 1)
        * acb_cis_pi(lower_phase)
        * (acb(arb_rational(a / tau - 1)) * lower_h + lower_hp / pi_i)
        / 2
    )
    upper_joined = (
        -acb(0, 1)
        * ((-1) ** m)
        * acb_cis_pi(upper_phase)
        * (acb(arb_rational(a / tau + 2 * K + 1)) * upper_h + upper_hp / pi_i)
        / 2
    )
    outer = acb_cis_pi(a * a / (2 * tau))
    lower_current = outer * lower_joined
    upper_current = outer * upper_joined
    endpoint = lower_current + upper_current
    complete = main + endpoint
    source, source_period = source_current_ball(x, s)
    difference = complete - source
    require(complete.overlaps(source), "certified Mordell current does not overlap exact-period source")
    require(difference.contains(0), "certified Mordell/source difference excludes zero")

    source_mass = arb(L) * acb(arb_rational(c + L - 1)).real
    normalized_radius = complete.rad() / source_mass
    denominator = complete.abs_lower()
    if not bool(denominator > arb(1)):
        denominator = arb(1)
    condition_ratio = (main.abs_upper() + endpoint.abs_upper()) / denominator
    unique_mordell = [item[2] for item in cache.values()]
    return {
        "x": str(x),
        "s": str(s),
        "a": str(a),
        "tau": str(tau),
        "integer_shift_r": shift,
        "normalized_z": str(z),
        "transformed_m": m,
        "transformed_w": str(w),
        "transformed_sigma": str(sigma),
        "transformed_period": transformed_period,
        "source_period": source_period,
        "lower_argument": str(lower_argument),
        "upper_argument": str(upper_argument),
        "unique_Mordell_evaluations": len(unique_mordell),
        "Mordell_enclosures": unique_mordell,
        "transformed_weighted_current": acb_record(transformed),
        "main_current": acb_record(main),
        "lower_endpoint_current": acb_record(lower_current),
        "upper_endpoint_current": acb_record(upper_current),
        "joined_endpoint_current": acb_record(endpoint),
        "complete_current": acb_record(complete),
        "exact_period_source_current": acb_record(source),
        "complete_minus_source": acb_record(difference),
        "difference_contains_zero": True,
        "complete_overlaps_source": True,
        "source_absolute_term_mass": str(source_mass),
        "complete_radius_over_source_mass_upper": arb_upper_text(normalized_radius),
        "main_plus_endpoint_condition_ratio_upper": arb_upper_text(condition_ratio),
        "settings": {
            "precision_bits": settings.precision_bits,
            "cutoff": settings.cutoff,
            "abs_tol": settings.abs_tol,
            "deg_limit": settings.deg_limit,
            "eval_limit": settings.eval_limit,
            "depth_limit": settings.depth_limit,
        },
        "passed": True,
    }


def render_note(artifact: dict[str, Any]) -> str:
    summary = artifact["summary"]
    return f"""# Certified Mordell interval one-step current

Date: 2026-08-24

Status: uniform central-strip tail bounds and two full-roster pointwise
Mordell current enclosures certified; parameter-box recursion open

## Exact tail bounds

For `theta=exp(i*pi/4)`, Kuznetsov's rotated representation is

```text
h(z,tau)=2 theta integral_0^infinity exp(-pi tau y^2)
          cosh(2pi z theta y)/cosh(pi theta y) dy.     (MI1)
```

For real `|z|<=1/2`, set `u=pi*y/sqrt(2)`.  The elementary identities

```text
|cosh(a+ib)|^2=sinh(a)^2+cos(b)^2,
|sinh(a+ib)|^2=sinh(a)^2+sin(b)^2                   (MI2)
```

give, for `y>=Y>0`,

```text
|cosh(2pi z theta y)/cosh(pi theta y)|<=coth(u),
|sinh(2pi z theta y)/cosh(pi theta y)|<=coth(u).     (MI3)
```

Since `coth(u)` decreases, truncating (MI1) and its `z` derivative at `Y`
has the uniform modulus bounds

```text
T_h(Y,tau)
 <=coth(pi Y/sqrt(2))*erfc(sqrt(pi tau)Y)/sqrt(tau),

T_hz(Y,tau)
 <=2*coth(pi Y/sqrt(2))*exp(-pi tau Y^2)/tau.        (MI4)
```

No asymptotic constant is hidden in (MI4).  The finite segments are enclosed
by Arb/ACB complex ball integration.  Each complex tail disk is added as a
containing real-imaginary rectangle.  Arguments outside the central strip
are transported by the exact Mordell unit recurrence and negative `tau` by
conjugation.

## Full one-step certificates

Production uses {artifact['production_settings']['precision_bits']} bits,
`Y={artifact['production_settings']['cutoff']}`, and absolute integration
tolerance `{artifact['production_settings']['abs_tol']}`.  At
`L={L}`, both `x=1/2,s=0` and `x=2/5,s=1/3` complete-current balls overlap
their independently enclosed exact rational-period source currents, and the
ball for their difference contains zero.

```text
largest complete radius / source mass:
  {summary['maximum_complete_radius_over_source_mass_upper']},

largest h tail bound:
  {summary['maximum_h_tail_modulus_upper']},

largest h_z tail bound:
  {summary['maximum_h_prime_tail_modulus_upper']},

largest main-plus-endpoint condition ratio:
  {summary['maximum_condition_ratio_upper']}.        (MI5)
```

The transformed periods are `{summary['transformed_periods']}` and the
source periods are `{summary['source_periods']}`.  Thus neither million-term
sum is trusted to floating recurrence in this certificate.

## Decision

Equation (MI4) closes the analytic truncation tail uniformly on the central
strip, and the two physical test rows now have end-to-end rigorous one-step
enclosures.  This removes the nonrigorous quadrature qualification from
those two rows.  It does not yet give a parameter-box evaluator over all
`(x,s)`, a recursive interval algorithm, or the small-`tau`
Euler--Maclaurin branch.  Those remain necessary before physical quadrature
or a non-A bound.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Primary-source boundary: (MI1) and the exact unit recurrence are taken from
Alexey Kuznetsov, *Computing the truncated theta function via Mordell
integral*, arXiv:1306.4081v2.  The tail estimates (MI2)--(MI4) are derived
here.  The paper's practical Gauss--Laguerre rule is not used.

Pi provenance: every `pi` comes from the inherited quadratic Fourier phase
or the stated rotated Mordell representation.  No fitted or geometric
constant is inserted.

This gate proves (MI2)--(MI4) and certifies two full-roster pointwise interval
rows.  It does not prove a uniform parameter-box Mordell evaluator, the
small-`tau` branch, recursive interval control, physical quadrature, the
non-A bound, joined `R_after_A`, `R_Dir`, `Q_K-T`, an all-height theorem,
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    require(CHECKER.is_file(), "independent interval checker missing")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "interval dependency not passed")
    require(
        dependencies["Mordell_exact_one_step"]["decision"]["endpoint_complete_one_step_identity_built"] is True,
        "exact Mordell identity dependency drift",
    )

    rows = [compute_full_row(Fraction(case["x"]), Fraction(case["s"]), PRODUCTION) for case in FULL_CASES]
    h_tail_values = []
    hp_tail_values = []
    for row in rows:
        for item in row["Mordell_enclosures"]:
            h_tail_values.append(mp.mpf(item["central_enclosure"]["tail"]["h_tail_modulus_upper"]))
            hp_tail_values.append(mp.mpf(item["central_enclosure"]["tail"]["h_prime_tail_modulus_upper"]))
    summary = {
        "row_count": len(rows),
        "maximum_complete_radius_over_source_mass_upper": mp.nstr(
            max(mp.mpf(row["complete_radius_over_source_mass_upper"]) for row in rows), 30
        ),
        "maximum_h_tail_modulus_upper": mp.nstr(max(h_tail_values), 30),
        "maximum_h_prime_tail_modulus_upper": mp.nstr(max(hp_tail_values), 30),
        "maximum_condition_ratio_upper": mp.nstr(
            max(mp.mpf(row["main_plus_endpoint_condition_ratio_upper"]) for row in rows), 30
        ),
        "transformed_periods": [row["transformed_period"] for row in rows],
        "source_periods": [row["source_period"] for row in rows],
        "all_difference_balls_contain_zero": all(row["difference_contains_zero"] for row in rows),
    }
    artifact = {
        "kind": STEM,
        "status": "uniform_central_strip_Mordell_tail_bounds_and_two_full_pointwise_interval_one_steps_certified_parameter_box_recursion_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "B": B,
            "L": L,
            "K": K,
            "central_tail_domain": "real |z|<=1/2, tau>0, Y>0",
            "certified_full_rows": list(FULL_CASES),
            "workers": 1,
        },
        "primary_source": {
            "author": "Alexey Kuznetsov",
            "title": "Computing the truncated theta function via Mordell integral",
            "arxiv": "1306.4081v2",
            "url": "https://arxiv.org/abs/1306.4081v2",
            "used_for": "Equation (5) rotated Mordell representation and exact unit recurrence",
            "practical_Gauss_Laguerre_used": False,
            "imported_numerical_constants": False,
        },
        "tail_theorem": {
            "ratio_bound": "both central cosh/cosh and sinh/cosh modulus ratios are <=coth(pi*y/sqrt(2))",
            "h_tail": "coth(pi*Y/sqrt(2))*erfc(sqrt(pi*tau)*Y)/sqrt(tau)",
            "h_prime_tail": "2*coth(pi*Y/sqrt(2))*exp(-pi*tau*Y^2)/tau",
            "uniform_in": "real |z|<=1/2",
            "unnamed_constants": False,
        },
        "production_settings": {
            "precision_bits": PRODUCTION.precision_bits,
            "cutoff": PRODUCTION.cutoff,
            "abs_tol": PRODUCTION.abs_tol,
            "deg_limit": PRODUCTION.deg_limit,
            "eval_limit": PRODUCTION.eval_limit,
            "depth_limit": PRODUCTION.depth_limit,
        },
        "full_interval_rows": rows,
        "summary": summary,
        "decision": {
            "uniform_central_strip_tail_bounds_proved": True,
            "pointwise_rational_Mordell_h_and_h_prime_encloser_built": True,
            "exact_unit_recurrence_propagated_with_balls": True,
            "transformed_current_enclosed_by_exact_rational_period": True,
            "source_current_enclosed_by_exact_rational_period": True,
            "two_full_one_step_current_enclosures_certified": True,
            "nonrigorous_Gauss_Laguerre_route_used": False,
            "uniform_parameter_box_Mordell_evaluator_built": False,
            "small_tau_Euler_Maclaurin_branch_built": False,
            "recursive_interval_theta_evaluator_built": False,
            "physical_quadrature_completed": False,
            "non_A_bound_proved": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Promote the pointwise central encloser to parameter boxes or a rigorously scheduled point quadrature, derive the tied small-tau Euler--Maclaurin branch, and validate at least one complete recursive interval current while preserving main/endpoint cancellation.",
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
        },
        "proof_boundary": "Uniform analytic central-strip truncation tails and two full-roster pointwise Arb/ACB one-step current certificates only. No uniform parameter-box Mordell evaluator, small-tau Euler--Maclaurin branch, recursive interval theta evaluator, physical quadrature, non-A bound, joined R_after_A, R_Dir, Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"certified central Mordell tails and {len(rows)} full interval one-step currents",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
