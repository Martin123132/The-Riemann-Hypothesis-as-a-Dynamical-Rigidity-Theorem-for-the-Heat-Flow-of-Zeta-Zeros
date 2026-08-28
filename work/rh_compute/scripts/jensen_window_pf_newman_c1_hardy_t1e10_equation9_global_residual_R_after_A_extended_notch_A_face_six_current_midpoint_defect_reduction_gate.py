#!/usr/bin/env python3
"""Certify a six-current elementary reduction of the A-face defect."""

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

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "midpoint_partition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.json",
    "two_current_expansion": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate.json",
}

A = 159_577
T = 10_000_000_000
FIRST_MODE = 39_895
LAST_MODE = 39_936
LEFT_NUMERATOR = 79_789
RIGHT_NUMERATOR = 79_873
Q_CUT = 16
TERMS = 6
PRECISION = 110


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


def odd_double_factorial(index: int) -> int:
    if index < 0:
        return 1
    result = 1
    for value in range(1, index + 1, 2):
        result *= value
    return result


def symbolic_certificate() -> dict[str, Any]:
    n = sp.symbols("n", integer=True, nonnegative=True)
    q = sp.symbols("q", negative=True, real=True)
    recurrence_coefficient = sp.factorial2(2 * n - 1) / (sp.I * sp.pi) ** (n + 1)
    require(sp.simplify(recurrence_coefficient.subs(n, 0) - 1 / (sp.I * sp.pi)) == 0, "leading current coefficient drift")
    require(sp.simplify(recurrence_coefficient.subs(n, 5) - sp.factorial2(9) / (sp.I * sp.pi) ** 6) == 0, "sixth retained coefficient drift")

    x, y, delta = sp.symbols("x y delta", positive=True, real=True)
    c = sp.sqrt(2 / x)
    q_delta = -c * delta
    h_six = sum(
        sp.factorial2(2 * index - 1)
        / ((sp.I * sp.pi) ** (index + 1) * q_delta ** (2 * index + 1))
        for index in range(TERMS)
    )
    g_from_h = -1 / (sp.I * sp.pi) - y * c * h_six
    g_expected = A * x / (2 * sp.I * sp.pi * delta)
    for index in range(1, TERMS):
        g_expected += (
            y
            * sp.factorial2(2 * index - 1)
            * (x / 2) ** index
            / ((sp.I * sp.pi) ** (index + 1) * delta ** (2 * index + 1))
        )
    require(sp.simplify(g_from_h.subs(y, delta + sp.Rational(A, 2) * x) - g_expected.subs(y, delta + sp.Rational(A, 2) * x)) == 0, "six-current phase-stripped amplitude drift")

    a = sp.symbols("a", positive=True, real=True)
    delta_y = y - a
    for index in range(1, TERMS):
        primitive = delta_y ** (1 - 2 * index) / (1 - 2 * index) - a * delta_y ** (-2 * index) / (2 * index)
        require(sp.simplify(sp.diff(primitive, y) - y / delta_y ** (2 * index + 1)) == 0, f"strip primitive drift at current {index}")

    bound_shape = y * x**TERMS / (y - sp.Rational(A, 2) * x) ** (2 * TERMS + 1)
    derivative_x = sp.factor(sp.diff(bound_shape, x))
    derivative_y = sp.factor(sp.diff(bound_shape, y))
    require(sp.simplify(derivative_x - x ** (TERMS - 1) * y * (TERMS * y + sp.Rational(A, 2) * (TERMS + 1) * x) / (y - sp.Rational(A, 2) * x) ** (2 * TERMS + 2)) == 0, "tail-bound x monotonicity drift")
    require(sp.simplify(derivative_y + x**TERMS * (2 * TERMS * y + sp.Rational(A, 2) * x) / (y - sp.Rational(A, 2) * x) ** (2 * TERMS + 2)) == 0, "tail-bound y monotonicity drift")

    return {
        "tail_hierarchy": "H(q)=sum_(n=0)^5 (2n-1)!!/[(i*pi)^(n+1)q^(2n+1)]+R_H,6(q)",
        "tail_remainder": "|R_H,6(q)|<=2*11!!/[pi^7*|q|^13]",
        "six_current_amplitude": "G_6(y,x)=A*x/[2*i*pi*delta]+sum_(n=1)^5 y(2n-1)!!(x/2)^n/[(i*pi)^(n+1)delta^(2n+1)], delta=y-A*x/2",
        "pointwise_amplitude_error": "|G_A-G_6|<=2*11!!*y*(x/2)^6/[pi^7*delta^13]",
        "error_monotonicity": "The error majorant increases with x and decreases with y on delta>0.",
        "defect": "D_6(x)=sum_(m=39895)^39936 G_6(m,x)-integral_(39894.5)^39936.5 G_6(y,x)dy",
        "oriented_channel": "E_42(x)=-D_42(x), hence E_42(x)=-D_6(x)+R_D,6(x)",
        "strip_current_zero": "I_0=A*x*log(delta_r/delta_l)/(2*i*pi)",
        "strip_current_n": "I_n=b_n{[delta^(1-2n)/(1-2n)-a*delta^(-2n)/(2n)]_(delta_l)^(delta_r)}, n=1,...,5",
        "zero_endpoint_limit": "lim_(x down 0)D_6(x)/x=A[sum_m 1/m-log(r/l)]/(2*i*pi)+[sum_m 1/m^2-(1/l-1/r)]/[2(i*pi)^2]",
    }


def g_six(y: arb, x: arb) -> acb:
    pi = arb.pi()
    i_pi = acb(0, pi)
    delta = y - arb(A) * x / 2
    value = acb(arb(A) * x) / (2 * i_pi * delta)
    for index in range(1, TERMS):
        coefficient = arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / i_pi ** (index + 1)
        value += acb(y) * coefficient / delta ** (2 * index + 1)
    return value


def integral_g_six(x: arb) -> acb:
    pi = arb.pi()
    i_pi = acb(0, pi)
    left = arb(LEFT_NUMERATOR) / 2
    right = arb(RIGHT_NUMERATOR) / 2
    a = arb(A) * x / 2
    delta_left = left - a
    delta_right = right - a
    value = acb(arb(A) * x * (delta_right / delta_left).log()) / (2 * i_pi)
    for index in range(1, TERMS):
        coefficient = arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / i_pi ** (index + 1)
        primitive_difference = (
            (delta_right ** (1 - 2 * index) - delta_left ** (1 - 2 * index)) / (1 - 2 * index)
            - a * (delta_right ** (-2 * index) - delta_left ** (-2 * index)) / (2 * index)
        )
        value += coefficient * primitive_difference
    return value


def defect_six(x: arb) -> acb:
    value = acb(0)
    for mode in range(FIRST_MODE, LAST_MODE + 1):
        value += g_six(arb(mode), x)
    return value - integral_g_six(x)


def numerical_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    a = arb(A)
    t = arb(T)
    left = arb(LEFT_NUMERATOR) / 2
    right = arb(RIGHT_NUMERATOR) / 2
    q_cut = arb(Q_CUT)
    sqrt_x_cut = ((2 * q_cut**2 + 8 * a * left).sqrt() - q_cut * arb(2).sqrt()) / (2 * a)
    x_cut = sqrt_x_cut**2
    x_star = (1 - (1 - 8 * t / (pi * a**2)).sqrt()) / 2

    delta_cut = left - a * x_cut / 2
    pointwise_error = (
        2
        * arb(odd_double_factorial(2 * TERMS - 1))
        * left
        * (x_cut / 2) ** TERMS
        / (pi ** (TERMS + 1) * delta_cut ** (2 * TERMS + 1))
    )
    defect_error = 84 * pointwise_error
    require(pointwise_error.upper() < arb("1.22e-10"), "six-current pointwise error drift")
    require(defect_error.upper() < arb("1.025e-8"), "six-current defect error drift")

    saddle_defect = defect_six(x_star)
    handoff_defect = defect_six(x_cut)

    reciprocal_sum_1 = sum((arb(1) / mode for mode in range(FIRST_MODE, LAST_MODE + 1)), arb(0))
    reciprocal_sum_2 = sum((arb(1) / arb(mode) ** 2 for mode in range(FIRST_MODE, LAST_MODE + 1)), arb(0))
    i_pi = acb(0, pi)
    endpoint_limit = (
        acb(a) * (reciprocal_sum_1 - (right / left).log()) / (2 * i_pi)
        + (reciprocal_sum_2 - (1 / left - 1 / right)) / (2 * i_pi**2)
    )

    require(abs(saddle_defect).lower() > arb("0.2233") and abs(saddle_defect).upper() < arb("0.2234"), "saddle defect scale drift")
    require(abs(handoff_defect).lower() > arb(8) and abs(handoff_defect).upper() < arb("8.04"), "handoff defect scale drift")
    require(abs(endpoint_limit).upper() < arb("1.4e-9"), "regular endpoint limit drift")

    def complex_record(value: acb) -> dict[str, str]:
        return {
            "real_ball": value.real.str(80, more=True),
            "imag_ball": value.imag.str(80, more=True),
            "modulus_ball": abs(value).str(80, more=True),
        }

    return {
        "x_16_ball": x_cut.str(80, more=True),
        "uniform_pointwise_G_minus_G6_ball": pointwise_error.str(80, more=True),
        "uniform_D42_minus_D6_ball": defect_error.str(80, more=True),
        "six_current_saddle_defect": complex_record(saddle_defect),
        "six_current_x16_defect": complex_record(handoff_defect),
        "regularized_x0_D6_over_x_limit": complex_record(endpoint_limit),
    }


def render_note(artifact: dict[str, Any]) -> str:
    numerical = artifact["numerical_certificate"]
    return f"""# Six-current elementary A-face midpoint defect

Date: 2026-08-23

Status: exact certificate for the six-current rational-log reduction and its
uniform tail error on `0<x<=x_16`; tangential integration remains open

Write `delta(y,x)=y-A*x/2`.  Repeated integration by parts in the lower
Fresnel tail gives, for `q<0`,

```text
H(q)=exp(-i*pi*q^2/2)J_-(q)
 =sum_(n=0)^5 (2n-1)!!/[(i*pi)^(n+1)q^(2n+1)]
  +R_(H,6)(q),

|R_(H,6)(q)|<=2*11!!/[pi^7|q|^13].                 (SC1)
```

The final bound exposes the next boundary current and bounds the remaining
`v^-14` integral absolutely; it is a finite remainder theorem, not a formal
asymptotic series.

Substitution into the phase-stripped A current gives the elementary six-term
amplitude

```text
G_6(y,x)=A*x/[2*i*pi*delta]
 +sum_(n=1)^5
  y(2n-1)!!(x/2)^n/[(i*pi)^(n+1)delta^(2n+1)].       (SC2)
```

On the tail region of Section 11.453,

```text
0<x<=x_16,       39894.5<=y<=39936.5,
```

the exact pointwise error obeys

```text
|G_A-G_6|
 <=2*11!!*y*(x/2)^6/[pi^7 delta^13].                (SC3)
```

The right side increases with `x` and decreases with `y`, so its maximum is
at `(x_16,39894.5)`.  Arb certifies

```text
sup |G_A-G_6| <= {numerical['uniform_pointwise_G_minus_G6_ball']}. (SC4)
```

Define the un-oriented midpoint defect

```text
D_6(x)=sum_(m=39895)^39936 G_6(m,x)
       -integral_(39894.5)^39936.5 G_6(y,x)dy.       (SC5)
```

Every strip integral is explicit.  With `a=A*x/2`,
`delta_l=39894.5-a`, `delta_r=39936.5-a`, and
`b_n=(2n-1)!!(x/2)^n/(i*pi)^(n+1)`,

```text
I_0=A*x*log(delta_r/delta_l)/(2*i*pi),

I_n=b_n
 [delta^(1-2n)/(1-2n)-a*delta^(-2n)/(2n)]
       |_(delta_l)^(delta_r),       1<=n<=5.         (SC6)
```

Thus `D_6=sum_m G_6(m,x)-sum_(n=0)^5 I_n` is a finite rational-log
function.  The exact oriented channel from Section 11.453 is
`E_42=-D_42`; consequently

```text
E_42(x)=-D_6(x)+R_(D,6)(x),
|R_(D,6)(x)|
 <={numerical['uniform_D42_minus_D6_ball']}.          (SC7)
```

The factor `84` in (SC7) is the 42 point masses plus the strip length 42.
No sampled cancellation is used in this bound.

At the two route landmarks, direct Arb evaluation of the elementary formula
gives

```text
|D_6(x_*)|  ={numerical['six_current_saddle_defect']['modulus_ball']},
|D_6(x_16)| ={numerical['six_current_x16_defect']['modulus_ball']}. (SC8)
```

The apparent `1/x` in the inherited Kummer mode response is harmless at the
lower endpoint.  The exact regularized limit is

```text
lim_(x down 0)D_6(x)/x
 =A[sum_m 1/m-log(r/l)]/(2*i*pi)
  +[sum_m 1/m^2-(1/l-1/r)]/[2(i*pi)^2]
 ={numerical['regularized_x0_D6_over_x_limit']['real_ball']}
  +i*{numerical['regularized_x0_D6_over_x_limit']['imag_ball']}. (SC9)
```

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (SC1)--(SC9) is inherited from the canonical
Fresnel primitive and exact Kummer/Fourier phase.  No fitted occurrence is
introduced.

Proof boundary: exact six-current expansion, elementary rational-log defect,
uniform tail error below `1.025e-8`, and regular lower-endpoint limit on
`0<x<=x_16` only.  No Morse-core integral, endpoint-layer integration by
parts, full signed A-face bound, quantitative `R_after_A`, `R_Dir`, or
`Q_K-T` estimate, all-height theorem, `Lambda<=0`, PF-infinity, RH, or
prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "an exact dependency gate is not passed")
    require(dependencies["midpoint_partition"]["decision"]["uniform_abs_q_at_least_16_before_x_16_proved"] is True, "tail partition drift")
    require(dependencies["midpoint_partition"]["decision"]["midpoint_Peano_identity_proved"] is True, "midpoint identity drift")
    require(dependencies["two_current_expansion"]["decision"]["first_two_normal_currents_derived_exactly"] is True, "current hierarchy drift")

    artifact = {
        "kind": STEM,
        "status": "exact_six_current_rational_log_A_face_midpoint_defect_and_uniform_tail_remainder_certified_tangential_integrals_open",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "shifted_modes": [FIRST_MODE, LAST_MODE],
            "continuous_strip": ["79789/2", "79873/2"],
            "x_region": "0<x<=x_16",
            "normal_tail_floor": Q_CUT,
            "retained_currents": TERMS,
        },
        "symbolic_certificate": symbolic_certificate(),
        "numerical_certificate": numerical_certificate(),
        "decision": {
            "six_current_Fresnel_tail_with_finite_remainder_proved": True,
            "phase_stripped_G6_rational_formula_proved": True,
            "G6_strip_integral_reduced_to_rational_log_endpoints": True,
            "uniform_G_minus_G6_below_1_point_22e_minus_10_proved": True,
            "uniform_signed_defect_error_below_1_point_025e_minus_8_proved": True,
            "D6_over_x_regular_at_zero_proved": True,
            "Morse_core_integral_bounded": False,
            "endpoint_layer_IBP_remainder_bounded": False,
            "full_signed_A_face_bound_proved": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
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
            "arb_threads": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Insert -D_6/x with its exact Kummer weight into the existing A-face Morse coordinate on 0<x<=x_16. Certify the rational-log amplitude and its six-current remainder on a core/tail interval cover. Separately retain the exact erfcx/Peano amplitude on x_16<=x<=1/2 and close one tangential integration-by-parts remainder using Phi_A'>21000.",
        "proof_boundary": "Exact six-current expansion, elementary rational-log defect, uniform tail error below 1.025e-8, and regular lower-endpoint limit on 0<x<=x_16 only. No Morse-core integral, endpoint-layer integration by parts, full signed A-face bound, quantitative R_after_A, R_Dir, or Q_K-T estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified six-current elementary A-face midpoint defect", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
