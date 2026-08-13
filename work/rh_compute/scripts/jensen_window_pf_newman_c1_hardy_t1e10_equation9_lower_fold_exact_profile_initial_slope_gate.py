#!/usr/bin/env python3
"""Certify the exact-minus-beta4 common-profile slope at y=0."""

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
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_finite_t_completed_strip_error_gate as point_gate


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_profile_initial_slope_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "common_profile_operator": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_exact_finite_t_common_profile_operator_gate.json",
    "beta4_completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "beta6_leading_coefficient": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta6_common_profile_leading_coefficient_gate.json",
}

PRECISION = 95
RAY_RADIUS = 12
RAY_PANELS = 24
TOLERANCE = "1e-29"


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


def symbolic_certificate() -> dict[str, str]:
    q, z, y = sp.symbols("q z y", real=True)
    i = sp.I
    amplitude = sp.cosh(q * z) ** (-sp.Rational(3, 2)) * (1 + q**2 * y / 2)
    phase_remainder = (
        (q * z - sp.tanh(q * z)) / q**3
        - z**3 / 3
        + y * (z - sp.tanh(q * z) / q)
        - q * y**2 * sp.tanh(q * z) / 4
    )
    ratio = sp.series(amplitude * sp.exp(i * phase_remainder), q, 0, 6).removeO().expand()
    c1 = sp.expand(ratio.coeff(q, 2))
    c2 = sp.expand(ratio.coeff(q, 4))
    d4 = sp.expand(
        -i * z
        + q**2 * (sp.diff(c1, y).subs(y, 0) - i * z * c1.subs(y, 0))
        + q**4 * (sp.diff(c2, y).subs(y, 0) - i * z * c2.subs(y, 0))
    )
    expected = (
        -i * z
        + q**2 * (sp.Rational(1, 2) + sp.Rational(13, 12) * i * z**3 - sp.Rational(2, 15) * z**6)
        + q**4 * (
            -sp.Rational(3, 8) * z**2
            - sp.Rational(137, 160) * i * z**5
            + sp.Rational(25, 126) * z**8
            + sp.Rational(2, 225) * i * z**11
        )
    )
    require(sp.expand(d4 - expected) == 0, "beta4 initial-slope multiplier drift")

    return {
        "exact_unphased_profile": (
            "H_ex(lambda,y)=int_R (1+y/(2beta^2)) cosh(z/beta)^(-3/2) "
            "exp(i{beta^3[z/beta-tanh(z/beta)]-lambda*z-beta*y*tanh(z/beta)}) dz"
        ),
        "initial_exact_slope": (
            "d_y H_ex(lambda,0)=int_R cosh(z/beta)^(-3/2) exp(i Phi_0) "
            "[1/(2beta^2)-i beta tanh(z/beta)] dz"
        ),
        "beta4_initial_slope_multiplier": str(expected),
        "shared_carrier": "exp(i[z^3/3-lambda*z])",
        "ray_recombination": "The left ray is minus the conjugate of the right ray, so both real slopes equal twice the right-ray real part.",
    }


def beta4_slope_multiplier(z: acb, epsilon: arb) -> acb:
    i = acb(0, 1)
    return (
        -i * z
        + epsilon * (arb(1) / 2 + arb(13) * i * z**3 / 12 - arb(2) * z**6 / 15)
        + epsilon**2 * (
            -arb(3) * z**2 / 8
            - arb(137) * i * z**5 / 160
            + arb(25) * z**8 / 126
            + arb(2) * i * z**11 / 225
        )
    )


def cubic_tail_power(radius: arb, power: int, rate: arb) -> arb:
    """Upper-bound int_R^infinity r^power exp(-rate*r^3) dr."""
    exponential = (-rate * radius**3).exp()
    if power <= 1:
        return radius ** (power - 2) * exponential / (3 * rate)
    return (
        radius ** (power - 2) * exponential / (3 * rate)
        + arb(power - 2) / (3 * rate) * cubic_tail_power(radius, power - 3, rate)
    )


def formal_beta6_slope(p: dict[str, arb], lam: arb) -> arb:
    pi = p["pi"]
    epsilon = p["beta"] ** -2
    u3 = (-arb(49856) * lam**6 - arb(516535) * lam**3 - 265860) / 9072000
    v3 = -(
        arb(3584) * lam**7 - arb(43880) * lam**4 - arb(127530) * lam
    ) / 9072000
    u3_y = (arb(1344) * lam**5 + arb(5445) * lam**2 + arb(21195) * lam) / 9072000
    v3_y = -(-arb(1792) * lam**6 + arb(18580) * lam**3 + 42570) / 9072000
    airy = (-lam).airy_ai()
    airy_prime = (-lam).airy_ai(derivative=1)
    # V3 multiplies d_X Ai(-X)=-Ai'(-X), not the argument derivative Ai'(-X).
    return 2 * pi * epsilon**3 * ((u3_y - lam * v3) * airy - (u3 + v3_y) * airy_prime)


def numerical_certificate(p: dict[str, arb]) -> dict[str, Any]:
    i = acb(0, 1)
    rotation = (i * p["pi"] / 6).exp()
    epsilon = p["beta"] ** -2
    lambda_max = p["pi"] / (16 * p["beta"])
    lam = arb(lambda_max / 2, lambda_max / 2)

    def integrand(radius: acb, _: bool) -> acb:
        z = rotation * radius
        scaled = z / p["beta"]
        amplitude = scaled.cosh() ** (-arb(3) / 2)
        ratio_phase = p["beta"] ** 3 * point_gate.argument_minus_tanh_polynomial(scaled) - z**3 / 3
        exact_ratio = amplitude * (i * ratio_phase).exp()
        exact_multiplier = arb(1) / (2 * p["beta"] ** 2) - i * p["beta"] * point_gate.tanh_polynomial(scaled)
        carrier = (i * (z**3 / 3 - lam * z)).exp()
        return rotation * carrier * (exact_ratio * exact_multiplier - beta4_slope_multiplier(z, epsilon))

    ray = acb(0)
    rows: list[dict[str, str]] = []
    radius = arb(RAY_RADIUS)
    for index in range(RAY_PANELS):
        left = radius * index / RAY_PANELS
        right = radius * (index + 1) / RAY_PANELS
        value = acb.integral(
            integrand,
            left,
            right,
            abs_tol=arb(TOLERANCE),
            rel_tol=arb(TOLERANCE),
            eval_limit=500_000,
            depth_limit=55,
        )
        ray += value
        rows.append({
            "panel_index": str(index),
            "left": left.str(20, more=True),
            "right": right.str(20, more=True),
            "real_ball": value.real.str(PRECISION, more=True),
            "imag_ball": value.imag.str(PRECISION, more=True),
        })

    rho = radius / p["beta"]
    tanh_tail = 3 * rho**18 / (1 - rho)
    phase_tail = p["beta"] ** 3 * tanh_tail
    tanh_bound = rho.sinh() / rho.cos()
    multiplier_bound = arb(1) / (2 * p["beta"] ** 2) + p["beta"] * tanh_bound
    carrier_growth = (lambda_max * radius / 2).exp()
    polynomial_replacement = (
        2 * radius * carrier_growth * phase_tail.exp()
        * (multiplier_bound * phase_tail * phase_tail.exp() + p["beta"] * tanh_tail)
    )

    near_rate = radius**2 / 5
    near_exponential = (-radius**3 / 5).exp()
    near_i0 = near_exponential / near_rate
    near_i1 = near_exponential * (radius / near_rate + 1 / near_rate**2)
    exact_near_tail = 2 * (near_i0 / (2 * p["beta"] ** 2) + arb("2.2") * near_i1)
    far_rate = p["beta"] ** 2 / 5 - lambda_max / 2
    require(far_rate.lower() > 0, "exact far-ray decay rate lost positivity")
    exact_far_tail = 2 * (arb(1) / (2 * p["beta"] ** 2) + 2 * p["beta"]) * (-far_rate * p["beta"]).exp() / far_rate

    canonical_coefficients = {
        0: epsilon / 2,
        1: arb(1),
        2: arb(3) * epsilon**2 / 8,
        3: arb(13) * epsilon / 12,
        5: arb(137) * epsilon**2 / 160,
        6: arb(2) * epsilon / 15,
        8: arb(25) * epsilon**2 / 126,
        11: arb(2) * epsilon**2 / 225,
    }
    canonical_tail = 2 * sum(
        coefficient * cubic_tail_power(radius, power, arb(1) / 4)
        for power, coefficient in canonical_coefficients.items()
    )
    total_error = polynomial_replacement + exact_near_tail + exact_far_tail + canonical_tail
    slope = 2 * ray.real + arb(0, total_error)
    formal = formal_beta6_slope(p, lam)
    post_beta6 = slope - formal

    require(abs(slope) < arb("5.56e-22"), "uniform initial slope exceeds 5.56e-22")
    require(slope.overlaps(formal), "initial slope does not overlap the beta^-6 coefficient enclosure")
    return {
        "lambda_interval": "0<=lambda<=pi/(16beta)",
        "lambda_max_ball": lambda_max.str(PRECISION, more=True),
        "ray_radius": RAY_RADIUS,
        "ray_panels": RAY_PANELS,
        "right_ray_difference_ball": complex_record(ray),
        "degree17_tanh_tail_ball": tanh_tail.str(PRECISION, more=True),
        "phase_tail_ball": phase_tail.str(PRECISION, more=True),
        "polynomial_replacement_error_bound_ball": polynomial_replacement.str(PRECISION, more=True),
        "two_ray_exact_near_tail_bound_ball": exact_near_tail.str(PRECISION, more=True),
        "two_ray_exact_far_tail_bound_ball": exact_far_tail.str(PRECISION, more=True),
        "two_ray_beta4_tail_bound_ball": canonical_tail.str(PRECISION, more=True),
        "total_appended_error_radius_ball": total_error.str(PRECISION, more=True),
        "uniform_exact_minus_beta4_initial_slope_ball": slope.str(PRECISION, more=True),
        "uniform_formal_beta6_initial_slope_ball": formal.str(PRECISION, more=True),
        "uniform_post_beta6_initial_slope_ball": post_beta6.str(PRECISION, more=True),
        "panel_enclosures": rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    n = artifact["numerical_certificate"]
    return f"""# Exact common-profile initial slope

Date: 2026-08-13

Status: rigorous exact-minus-beta-minus-four initial slope on the top selector
corridor; this is not yet the full profile or a proof of `T_upper` or RH

After removing the common quadratic phase, write the exact transformed fold as

```text
H_ex(lambda,y)=int_R (1+y/(2beta^2)) cosh(z/beta)^(-3/2)
 exp(i{{beta^3[z/beta-tanh(z/beta)]-lambda*z-beta*y*tanh(z/beta)}}) dz. (IS1)
```

At `y=0`, direct differentiation gives

```text
d_y H_ex=int_R cosh(z/beta)^(-3/2) exp(i Phi_0)
 [1/(2beta^2)-i beta tanh(z/beta)] dz.                 (IS2)
```

The exact and beta-minus-four derivatives share the carrier
`exp(i[z^3/3-lambda*z])`.  Subtracting their ratio multipliers before taking
absolute values, then using a degree-17 tanh polynomial on `0<=r<=12`, proves

```text
sup_(0<=lambda<=pi/(16beta)) |d_y(H_ex-H4)(lambda,0)|
 <= abs({n['uniform_exact_minus_beta4_initial_slope_ball']})
 < 5.56e-22.                                             (IS3)
```

The appended compact-replacement and contour-tail radius is
`{n['total_appended_error_radius_ball']}`.  The independently derived
beta-minus-six coefficient gives

```text
d_y H6(lambda,0)={n['uniform_formal_beta6_initial_slope_ball']},
d_y(H_ex-H4-H6)(lambda,0)
 ={n['uniform_post_beta6_initial_slope_ball']}.           (IS4)
```

The left Airy ray is minus the conjugate of the right ray, so (IS3) is a
two-ray real integral, not a one-sided estimate.  Pi enters only through the
Airy-ray angle `pi/6`, `beta^3=pi*C^2/8`, and the selector width
`lambda_max=pi/(16beta)`.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No full exact-minus-beta-minus-four profile bound, weighted finite-height
remainder, complete `T_upper`, height-uniform theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    for name, dependency in dependencies.items():
        require(dependency.get("passed") is True, f"dependency failed: {name}")

    ctx.dps = PRECISION
    p = point_gate.parameters()
    symbolic = symbolic_certificate()
    numerical = numerical_certificate(p)
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "exact_minus_beta4_common_profile_initial_slope_certified",
        "passed": True,
        "parameters": {
            "lower_endpoint_C": point_gate.LOWER_ALPHA,
            "beta_ball": p["beta"].str(PRECISION, more=True),
            "Y_ball": p["Y"].str(PRECISION, more=True),
        },
        "symbolic_certificate": symbolic,
        "numerical_certificate": numerical,
        "decision": {
            "exact_profile_initial_slope_formula_proved": True,
            "exact_and_beta4_integrands_subtracted_before_absolute_values": True,
            "uniform_initial_slope_below_5_56e_minus_22": True,
            "beta6_slope_basis_and_interval_overlap_certified": True,
            "full_exact_minus_beta4_profile_bound_proved": False,
            "complete_T_upper_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in DEPENDENCIES.items()
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
        "next_action": (
            "Combine this initial slope with the completed-point exact-minus-beta4 value, the exact Kummer ODE, "
            "and its exact beta^-6 residual in an Airy-Green/Volterra whole-profile bound."
        ),
        "proof_boundary": (
            "A rigorous exact-minus-beta4 initial slope on one top selector corridor only. No whole-profile weighted "
            "finite-height transfer, complete T_upper theorem, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("certified exact-minus-beta4 common-profile initial slope < 5.56e-22", flush=True)


if __name__ == "__main__":
    main()
