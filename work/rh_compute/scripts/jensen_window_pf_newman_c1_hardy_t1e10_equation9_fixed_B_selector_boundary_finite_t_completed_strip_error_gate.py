#!/usr/bin/env python3
"""Certify the completed finite-t selector strip against its canonical fold."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, Path(__file__).resolve().parent):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx
import sympy as sp


COMPLETION_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate.json"
NORMAL_FORM_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
SELECTOR_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 90
LOWER_ALPHA = 159_577
FIXED_B = 5_122_423
RAY_RADIUS = 12
RAY_PANELS = 24
TOLERANCE = "1e-28"

TANH_COEFFICIENTS = {
    1: Fraction(1),
    3: Fraction(-1, 3),
    5: Fraction(2, 15),
    7: Fraction(-17, 315),
    9: Fraction(62, 2835),
    11: Fraction(-1382, 155925),
    13: Fraction(21844, 6081075),
    15: Fraction(-929569, 638512875),
    17: Fraction(6404582, 10854718875),
}


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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def parameters() -> dict[str, arb]:
    ctx.dps = PRECISION
    pi = arb.pi()
    lower = arb(LOWER_ALPHA)
    height = pi * lower**2 / 8
    beta = height ** (arb(1) / 3)
    sigma = 4 * beta / (pi * lower)
    y_width = 2 / sigma
    h = 4 * beta / lower
    return {
        "pi": pi,
        "lower": lower,
        "height": height,
        "beta": beta,
        "sigma": sigma,
        "Y": y_width,
        "h": h,
    }


def symbolic_certificate() -> dict[str, str]:
    z, epsilon = sp.symbols("z epsilon", real=True)
    a = epsilon * z
    amplitude = sp.cosh(a) ** (-sp.Rational(3, 2))
    phase = (a - sp.tanh(a)) / epsilon**3
    ratio = sp.series(
        amplitude * sp.exp(sp.I * (phase - z**3 / 3)),
        epsilon,
        0,
        6,
    ).removeO().expand()
    q2 = sp.expand(ratio).coeff(epsilon, 2)
    q4 = sp.expand(ratio).coeff(epsilon, 4)
    require(q2 == -sp.Rational(3, 4) * z**2 - sp.Rational(2, 15) * sp.I * z**5, "beta^-2 coefficient drift")
    require(
        q4 == sp.Rational(13, 32) * z**4 + sp.Rational(97, 630) * sp.I * z**7 - sp.Rational(2, 225) * z**10,
        "beta^-4 coefficient drift",
    )

    x = sp.symbols("x", real=True)
    airy, airy_prime = sp.symbols("Ai0 Aip0", real=True)
    p, q = sp.Integer(1), sp.Integer(0)
    moments: list[tuple[sp.Expr, sp.Expr]] = []
    for _ in range(11):
        moments.append((sp.expand(p.subs(x, 0)), sp.expand(q.subs(x, 0))))
        p, q = sp.diff(p, x) + q * x, p + sp.diff(q, x)

    def integrate_polynomial(poly: sp.Expr) -> sp.Expr:
        value = sp.Integer(0)
        for (degree,), coefficient in sp.Poly(poly, z).terms():
            pp, qq = moments[degree]
            value += coefficient * sp.I ** (-degree) * (pp * airy + qq * airy_prime)
        return sp.simplify(value)

    require(integrate_polynomial(q2) == 0, "beta^-2 integrated correction does not cancel")
    require(integrate_polynomial(q4) == -sp.Rational(9, 560) * airy_prime, "beta^-4 integrated coefficient drift")

    return {
        "selector_boundary": "t*=pi*A^2/8=beta^3, eta=1, lambda=0",
        "one_period": "hY=2pi",
        "exact_completed_strip": "K_exact=(sqrt(2)/pi)Y exp(i*t*) J_beta",
        "canonical_completed_strip": "K_0=(sqrt(2)/pi)Y exp(i*t*) 2pi*Ai(0)",
        "exact_fold_integral": "J_beta=int_R cosh(z/beta)^(-3/2) exp(i*beta^3[z/beta-tanh(z/beta)])dz",
        "canonical_fold_integral": "J_0=int_Airy exp(i*z^3/3)dz=2pi*Ai(0)",
        "all_mode_collapse": "The symmetric Fourier midpoint plus endpoint half-current evaluates both completed strips at y=0 before z integration.",
        "source_carrier": "exp(i[t*-pi/8])=1 because A=159577 is 1 modulo 8",
        "formal_beta_minus_2_integral": "0",
        "formal_beta_minus_4_integral": "-(9/560)*2pi*Ai'(0)",
    }


def tanh_polynomial(value: acb) -> acb:
    total = acb(0)
    for degree, coefficient in TANH_COEFFICIENTS.items():
        total += arb(coefficient.numerator) / coefficient.denominator * value**degree
    return total


def argument_minus_tanh_polynomial(value: acb) -> acb:
    total = acb(0)
    for degree, coefficient in TANH_COEFFICIENTS.items():
        if degree >= 3:
            total -= arb(coefficient.numerator) / coefficient.denominator * value**degree
    return total


def sin_pi_over_six_multiple(degree: int) -> arb:
    residues = {
        1: arb(1) / 2,
        3: arb(1),
        5: arb(1) / 2,
        7: -arb(1) / 2,
        9: -arb(1),
        11: -arb(1) / 2,
    }
    return residues[degree % 12]


def interval_ball(left: arb, right: arb) -> arb:
    return arb((left + right) / 2, (right - left) / 2)


def ray_decay_certificate(p: dict[str, arb]) -> dict[str, Any]:
    root_three = arb(3).sqrt()

    def decay(rho: arb) -> arb:
        return rho / 2 - rho.sin() / ((root_three * rho).cosh() + rho.cos())

    def quotient_polynomial(rho: arb) -> arb:
        total = arb(0)
        for degree, coefficient in TANH_COEFFICIENTS.items():
            if degree >= 3:
                total -= (
                    arb(coefficient.numerator)
                    / coefficient.denominator
                    * sin_pi_over_six_multiple(degree)
                    * rho ** (degree - 3)
                )
        return total

    small = arb("0.25", "0.25")
    small_margin = quotient_polynomial(small) - arb(1) / 4 - 3 * small**15 / (1 - small)
    require(small_margin.lower() > 0, "small-ray cubic decay margin failed")

    mid_min: arb | None = None
    for index in range(256):
        left = arb(1) / 2 + arb(index) / 512
        right = arb(1) / 2 + arb(index + 1) / 512
        rho = interval_ball(left, right)
        margin = decay(rho) - rho**3 / 4
        require(margin.lower() > 0, f"mid-ray cubic decay failed at panel {index}")
        if mid_min is None or margin.lower() < mid_min:
            mid_min = margin.lower()

    bridge_min: arb | None = None
    for index in range(128):
        left = arb(1) + arb(index) / 640
        right = arb(1) + arb(index + 1) / 640
        rho = interval_ball(left, right)
        margin = decay(rho) - rho / 5
        require(margin.lower() > 0, f"bridge linear decay failed at panel {index}")
        if bridge_min is None or margin.lower() < bridge_min:
            bridge_min = margin.lower()

    endpoint = arb("1.2")
    analytic_tail_margin = 3 * endpoint / 10 - 1 / ((root_three * endpoint).cosh() - 1)
    require(analytic_tail_margin.lower() > 0, "analytic linear tail margin failed")
    require(arb(RAY_RADIUS) / p["beta"] < arb("0.006"), "compact ray radius leaves polynomial disk")
    return {
        "ray": "z=exp(i*pi/6)r",
        "exact_decay": "D(rho)=rho/2-sin(rho)/(cosh(sqrt(3)rho)+cos(rho))",
        "small_range": "0<=rho<=1/2: D(rho)>=rho^3/4 by degree-17 tanh polynomial plus Cauchy tail",
        "small_range_margin_ball": small_margin.str(PRECISION, more=True),
        "middle_range": "1/2<=rho<=1: D(rho)>=rho^3/4",
        "middle_panel_count": 256,
        "middle_minimum_margin_lower_bound": mid_min.str(PRECISION, more=True),
        "bridge_range": "1<=rho<=6/5: D(rho)>=rho/5",
        "bridge_panel_count": 128,
        "bridge_minimum_margin_lower_bound": bridge_min.str(PRECISION, more=True),
        "tail_range": "rho>=6/5: D(rho)>=rho/5 from |sin rho|<=1 and monotonic cosh(sqrt(3)rho)",
        "analytic_tail_margin_ball": analytic_tail_margin.str(PRECISION, more=True),
        "amplitude_guard": "|cosh(rho*exp(i*pi/6))|^2=[cosh(sqrt(3)rho)+cos(rho)]/2>=1, so the exact amplitude is at most one.",
    }


def integrate_ray_difference(p: dict[str, arb]) -> tuple[acb, list[dict[str, str]]]:
    i = acb(0, 1)
    rotation = (i * p["pi"] / 6).exp()

    def integrand(radius: acb, _: bool) -> acb:
        z = rotation * radius
        scaled = z / p["beta"]
        amplitude = scaled.cosh() ** (-arb(3) / 2)
        exact_surrogate = amplitude * (i * p["beta"] ** 3 * argument_minus_tanh_polynomial(scaled)).exp()
        canonical = (-radius**3 / 3).exp()
        return rotation * (exact_surrogate - canonical)

    total = acb(0)
    rows: list[dict[str, str]] = []
    for index in range(RAY_PANELS):
        left = arb(RAY_RADIUS * index) / RAY_PANELS
        right = arb(RAY_RADIUS * (index + 1)) / RAY_PANELS
        value = acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb(TOLERANCE),
            rel_tol=arb(TOLERANCE),
            eval_limit=400_000,
            depth_limit=50,
        )
        total += value
        rows.append({
            "panel_index": str(index),
            "left": left.str(20, more=True),
            "right": right.str(20, more=True),
            "real_ball": value.real.str(PRECISION, more=True),
            "imag_ball": value.imag.str(PRECISION, more=True),
        })
    return total, rows


def numerical_certificate(p: dict[str, arb], completion: dict[str, Any]) -> dict[str, Any]:
    ray, panels = integrate_ray_difference(p)
    radius = arb(RAY_RADIUS)
    rho = radius / p["beta"]
    tanh_tail = 3 * rho**18 / (1 - rho)
    phase_tail = p["beta"] ** 3 * tanh_tail
    polynomial_replacement_error = 2 * radius * phase_tail * phase_tail.exp()
    exact_ray_tail = 2 * (
        4 * (-radius**3 / 4).exp() / (3 * radius**2)
        + 5 * (-p["beta"] ** 3 / 5).exp() / p["beta"] ** 2
    )
    canonical_ray_tail = 2 * (-radius**3 / 3).exp() / radius**2
    total_error = polynomial_replacement_error + exact_ray_tail + canonical_ray_tail
    reduced_difference = 2 * ray.real + arb(0, total_error)

    airy_zero = arb(0).airy_ai()
    airy_prime_zero = arb(0).airy_ai(derivative=1)
    canonical_fold = 2 * p["pi"] * airy_zero
    leading_difference = -arb(9) / 560 * 2 * p["pi"] * airy_prime_zero / p["beta"] ** 4
    leading_ratio = reduced_difference / leading_difference

    physical_factor = arb(2).sqrt() / p["pi"]
    strip_scale = physical_factor * p["Y"]
    canonical_physical = strip_scale * canonical_fold
    physical_difference = strip_scale * reduced_difference
    exact_physical = canonical_physical + physical_difference
    equation9_normalizer = (p["pi"] / (32 * p["height"])) ** (arb(1) / 4)
    normalized_difference = equation9_normalizer * physical_difference
    normalized_canonical = equation9_normalizer * canonical_physical

    saved_physical = arb(completion["numerical_completion"]["physical_completed_removed_term"]["real_ball"])
    require(canonical_physical.overlaps(saved_physical), "canonical completion scale misses the 399-mode gate")
    require(reduced_difference.lower() > arb("1.2e-15"), "completed finite-t correction lost positivity")
    require(reduced_difference.upper() < arb("1.22e-15"), "completed finite-t correction exceeds scout range")
    require(arb("0.999999") < leading_ratio < arb("1.000001"), "beta^-4 leading term no longer predicts the correction")
    require(abs(physical_difference) < arb("6.4e-14"), "physical completed-strip correction exceeds 6.4e-14")
    require(abs(normalized_difference) < arb("1.2e-16"), "equation-(9) normalized correction exceeds 1.2e-16")

    parity_integer = (LOWER_ALPHA**2 - 1) // 8
    require(LOWER_ALPHA % 8 == 1 and parity_integer % 2 == 0, "source carrier parity failed")
    return {
        "ray_radius": RAY_RADIUS,
        "ray_panels": RAY_PANELS,
        "ray_surrogate_difference": complex_record(ray),
        "degree17_tanh_tail_ball": tanh_tail.str(PRECISION, more=True),
        "phase_tail_ball": phase_tail.str(PRECISION, more=True),
        "two_ray_polynomial_replacement_error_bound_ball": polynomial_replacement_error.str(PRECISION, more=True),
        "two_ray_exact_tail_bound_ball": exact_ray_tail.str(PRECISION, more=True),
        "two_ray_canonical_tail_bound_ball": canonical_ray_tail.str(PRECISION, more=True),
        "total_rigorous_error_radius_ball": total_error.str(PRECISION, more=True),
        "reduced_exact_minus_canonical_fold_ball": reduced_difference.str(PRECISION, more=True),
        "canonical_fold_2pi_Ai0_ball": canonical_fold.str(PRECISION, more=True),
        "formal_beta_minus_4_leading_ball": leading_difference.str(PRECISION, more=True),
        "exact_over_formal_leading_ratio_ball": leading_ratio.str(PRECISION, more=True),
        "strip_scale_Y_sqrt2_over_pi_ball": strip_scale.str(PRECISION, more=True),
        "canonical_completed_strip_physical_ball": canonical_physical.str(PRECISION, more=True),
        "exact_completed_strip_physical_ball": exact_physical.str(PRECISION, more=True),
        "finite_t_minus_canonical_completed_strip_physical_ball": physical_difference.str(PRECISION, more=True),
        "equation9_normalizer_ball": equation9_normalizer.str(PRECISION, more=True),
        "equation9_normalized_canonical_completed_strip_ball": normalized_canonical.str(PRECISION, more=True),
        "equation9_normalized_finite_t_minus_canonical_ball": normalized_difference.str(PRECISION, more=True),
        "source_carrier_exp_i_tstar_minus_pi_over_8": "1 (exact)",
        "source_carrier_parity_integer": parity_integer,
        "panel_enclosures": panels,
    }


def render_note(artifact: dict[str, Any]) -> str:
    p = artifact["parameters"]
    n = artifact["numerical_certificate"]
    r = artifact["ray_decay_certificate"]
    return f"""# Fixed-B selector-boundary finite-t completed-strip error

Date: 2026-08-12

Status: rigorous one-height completed-strip comparison; this is not a proof of
a nonzero height cell, the ordinary-Morse join, `T_upper`, `Lambda<=0`, or RH

At the exact selector boundary

```text
A={LOWER_ALPHA}, B={FIXED_B},
t*=pi*A^2/8={p['selector_boundary_height_ball']},
beta=t*^(1/3)={p['beta_ball']}.                         (FT1)
```

The exact identity `hY=2pi` makes every complete symmetric mode sum a Fourier
midpoint.  Adding the endpoint half-current evaluates the completed strip at
`y=0`.  This applies to the exact finite-`t` strip and to the canonical
Fourier object before either is estimated.  After suppressing their common
carrier, the comparison is therefore only

```text
J_beta=int_R cosh(z/beta)^(-3/2)
       exp(i*beta^3[z/beta-tanh(z/beta)])dz,
J_0=2pi*Ai(0).                                         (FT2)
```

Thus the 399-mode core and the zero/negative/outer-positive complement are
not estimated separately.  Their exact completion from Section 11.333 is
used first.

For the positive half-line deform to `z=exp(i*pi/6)r`.  With
`rho=r/beta`, the exact decay is

```text
D(rho)=rho/2-sin(rho)/(cosh(sqrt(3)rho)+cos(rho)).     (FT3)
```

Interval subdivision plus an analytic tail proves

```text
D(rho)>=rho^3/4,  0<=rho<=1,
D(rho)>=rho/5,    rho>=1.                              (FT4)
```

The compact proof uses {r['middle_panel_count']} panels on `1/2<=rho<=1`
and {r['bridge_panel_count']} panels on `1<=rho<=6/5`; the small range is
handled by the degree-17 tanh polynomial and its Cauchy tail.  Also
`|cosh(rho exp(i*pi/6))|>=1`, so no amplitude growth is discarded.

Twenty-four interval panels on `0<=r<=12`, followed by the explicit
polynomial and contour-tail bounds, give

```text
J_beta-J_0
 ={n['reduced_exact_minus_canonical_fold_ball']}.       (FT5)
```

The total error radius appended to the direct quadrature is

```text
{n['total_rigorous_error_radius_ball']}.                (FT6)
```

The apparent `beta^(-2)` correction cancels exactly after the two Airy rays
are recombined.  The first formal nonzero integrated term is

```text
-(9/560)*2pi*Ai'(0)/beta^4
 ={n['formal_beta_minus_4_leading_ball']},              (FT7)
```

and the rigorous result divided by (FT7) is
`{n['exact_over_formal_leading_ratio_ball']}`.

Restoring the transformed Kummer scale gives

```text
canonical completed strip
 ={n['canonical_completed_strip_physical_ball']},
exact finite-t completed strip
 ={n['exact_completed_strip_physical_ball']},
finite-t minus canonical
 ={n['finite_t_minus_canonical_completed_strip_physical_ball']}. (FT8)
```

After the paper's equation-(9) normalization, the certified difference is

```text
{n['equation9_normalized_finite_t_minus_canonical_ball']} < 1.2e-16. (FT9)
```

There is no hidden phase in (FT9): `A=159577` is `1 mod 8`, so
`exp(i[t*-pi/8])=1` exactly.

This closes the exact finite-`t` versus canonical comparison at the single
selector-boundary height.  It does not yet provide a uniform neighborhood in
height or join the completed strip to the ordinary-Morse modes.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No nonzero selector-height cell, complete `T_upper`, `Lambda<=0`, RH, or
prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for dependency in (COMPLETION_GATE, NORMAL_FORM_GATE, SELECTOR_GATE, CHECKER):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    completion = json.loads(COMPLETION_GATE.read_text(encoding="utf-8"))
    require(completion.get("passed") is True, "399-mode completion dependency did not pass")
    p = parameters()
    require((p["h"] * p["Y"] - 2 * p["pi"]).contains(0), "one-period identity failed")
    require((p["height"] - p["beta"] ** 3).contains(0), "selector fold eta is not one")
    symbolic = symbolic_certificate()
    decay = ray_decay_certificate(p)
    numerical = numerical_certificate(p, completion)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate",
        "status": "exact_one_height_finite_t_completed_strip_error_certified",
        "passed": True,
        "parameters": {
            "lower_selector_A": LOWER_ALPHA,
            "fixed_upper_endpoint_B": FIXED_B,
            "selector_boundary_height_ball": p["height"].str(PRECISION, more=True),
            "beta_ball": p["beta"].str(PRECISION, more=True),
            "sigma_ball": p["sigma"].str(PRECISION, more=True),
            "Y_ball": p["Y"].str(PRECISION, more=True),
            "h_ball": p["h"].str(PRECISION, more=True),
            "hY_minus_2pi_ball": (p["h"] * p["Y"] - 2 * p["pi"]).str(40, more=True),
        },
        "symbolic_certificate": symbolic,
        "ray_decay_certificate": decay,
        "numerical_certificate": numerical,
        "decision": {
            "all_mode_completion_taken_before_absolute_values": True,
            "exact_and_canonical_completed_strips_collapse_to_y0": True,
            "finite_t_completed_strip_error_rigorously_bounded_at_tstar": True,
            "equation9_normalized_error_below_1_2e_minus_16": True,
            "beta_minus_2_integrated_correction_cancels": True,
            "first_formal_correction_is_beta_minus_4": True,
            "source_carrier_is_exactly_one": True,
            "nonzero_height_cell_proved": False,
            "ordinary_Morse_join_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "fixed_B_399_mode_completion_gate": {"path": relative(COMPLETION_GATE), "sha256": file_hash(COMPLETION_GATE)},
            "grouped_characteristic_normal_form_gate": {"path": relative(NORMAL_FORM_GATE), "sha256": file_hash(NORMAL_FORM_GATE)},
            "selector_strip_reassembly_gate": {"path": relative(SELECTOR_GATE), "sha256": file_hash(SELECTOR_GATE)},
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
            "Differentiate or re-center the completed exact-minus-canonical strip in height, preserve the Fourier endpoint completion, "
            "and certify a nonzero selector cell that joins continuously to the ordinary-Morse chart."
        ),
        "proof_boundary": (
            "Exact all-mode collapse and a rigorous completed-strip finite-t error at the single height t*=pi*159577^2/8 only. "
            "No nonzero height cell, ordinary-Morse join, complete T_upper theorem, Lambda<=0, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built completed selector-strip finite-t error: normalized correction < 1.2e-16")


if __name__ == "__main__":
    main()
