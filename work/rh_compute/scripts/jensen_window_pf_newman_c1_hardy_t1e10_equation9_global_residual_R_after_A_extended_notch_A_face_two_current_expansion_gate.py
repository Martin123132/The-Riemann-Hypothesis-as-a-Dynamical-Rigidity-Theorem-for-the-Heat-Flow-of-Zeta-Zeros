#!/usr/bin/env python3
"""Certify the first two A-face normal currents at the translated edge."""

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
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_two_current_expansion_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "normal_tangential_split": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_tangential_split_gate.json",
    "finite_endpoint_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_finite_endpoint_tail_decomposition_gate.json",
    "A_projector_assembly": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate.json",
}

A = 159_577
D_NUMERATOR = 79_873
D = Fraction(D_NUMERATOR, 2)
LAST_SHIFTED_MODE = 39_936
PRECISION = 100


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


def symbolic_certificate() -> dict[str, Any]:
    x, d, delta = sp.symbols("x d delta", positive=True, real=True)
    q = -delta * sp.sqrt(2 / x)
    a = d * sp.sqrt(2 / x)
    first = -1 / (sp.I * sp.pi) - a / (sp.I * sp.pi * q)
    second = -a / ((sp.I * sp.pi) ** 2 * q**3)
    expected_first = A * x / (2 * sp.I * sp.pi * delta)
    expected_second = -d * x / (2 * sp.pi**2 * delta**3)
    require(sp.simplify(first.subs(delta, d - A * x / 2) - expected_first.subs(delta, d - A * x / 2)) == 0, "first current simplification failed")
    require(sp.simplify(second - expected_second) == 0, "second current simplification failed")

    v, p = sp.symbols("v p", real=True)
    exponential = sp.exp(sp.I * sp.pi * v**2 / 2)
    recurrence_lhs = sp.diff(exponential / (sp.I * sp.pi * v ** (p + 1)), v)
    recurrence_rhs = exponential * v ** (-p) - (p + 1) * exponential * v ** (-p - 2) / (sp.I * sp.pi)
    require(sp.simplify(recurrence_lhs - recurrence_rhs) == 0, "Fresnel-tail recurrence failed")

    d_exact = sp.Rational(D_NUMERATOR, 2)
    delta_exact = d_exact - sp.Rational(A, 2) * x
    q_exact = -delta_exact * sp.sqrt(2 / x)
    q2 = sp.simplify(q_exact**2)
    q2_derivative = sp.factor(sp.diff(q2, x))
    require(sp.simplify(q_exact.subs(x, sp.Rational(1, 2)) + sp.Rational(169, 2)) == 0, "corner q drift")
    require(D_NUMERATOR % 4 == 1 and A % 4 == 1, "half-integer phase parity drift")

    c0 = A * x / (2 * sp.pi * delta_exact)
    c1 = d_exact * x / (2 * sp.pi**2 * delta_exact**3)
    remainder = 3 * d_exact * x**2 / (2 * sp.pi**3 * delta_exact**5)
    for expression, label in ((c0, "c0"), (c1, "c1"), (remainder, "remainder")):
        numerator = sp.together(sp.diff(expression, x)).as_numer_denom()[0]
        require(all(coefficient >= 0 for coefficient in sp.Poly(numerator, x).all_coeffs()), f"{label} monotonicity polynomial drift")

    return {
        "translated_edge": "d_A=79873/2=39936.5",
        "normal_gap": "delta_A(x)=d_A-A*x/2",
        "normal_coordinate": "q_A(d_A,x)=-delta_A(x)*sqrt(2/x)",
        "uniform_coordinate_bound": "q_A(d_A,x)<=-169/2 for 0<x<=1/2",
        "Fresnel_tail": "J_-(q)=F(q)+(1+i)/2=integral_(-infinity)^q exp(i*pi*v^2/2)dv, q<0",
        "Fresnel_recurrence": "I_p(q)=exp(i*pi*q^2/2)/(i*pi*q^(p+1))+(p+1)I_(p+2)(q)/(i*pi)",
        "two_term_tail": "J_-(q)=E(q)[1/(i*pi*q)+1/(i*pi)^2*q^3]+r_2(q)",
        "tail_remainder": "|r_2(q)|<=6/(pi^3*|q|^5)",
        "continuous_A_tail_current": "L_A(d,x)=-E(q)/(i*pi)-d*sqrt(2/x)*J_-(q)",
        "first_current": "c_0(x)=A*x/[2*i*pi*delta_A(x)]",
        "second_current": "c_1(x)=-d_A*x/[2*pi^2*delta_A(x)^3]",
        "current_remainder": "L_A(d_A,x)=E(q_A)[c_0(x)+c_1(x)]+R_2(x), |R_2(x)|<=3*d_A*x^2/[2*pi^3*delta_A(x)^5]",
        "quadratic_completion": "exp(-i*pi*d_A^2/x)E(q_A)=exp(i*pi*A^2*x/4-i*pi*d_A*A)",
        "half_integer_phase": "exp(-i*pi*d_A*A)=-i",
        "q_squared_derivative": str(q2_derivative),
    }


def numerical_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    a = arb(A)
    d = arb(D_NUMERATOR) / 2
    x = arb(1) / 2
    delta = d - a * x / 2
    q_abs = delta * (2 / x).sqrt()
    c0 = a * x / (2 * pi * delta)
    c1 = d * x / (2 * pi**2 * delta**3)
    remainder = 3 * d * x**2 / (2 * pi**3 * delta**5)
    third_current = 3 * d * x**2 / (4 * pi**3 * delta**5)
    require(delta.contains(arb(169) / 4), "corner delta drift")
    require(q_abs.contains(arb(169) / 2), "corner q drift")
    require(remainder.upper() < arb("4e-6"), "two-current remainder threshold drift")

    return {
        "monotonicity": "|c_0|, |c_1|, and the stated |R_2| majorant increase on 0<x<=1/2, so their maxima occur at x=1/2.",
        "corner_delta_ball": delta.str(80, more=True),
        "minimum_abs_q_ball": q_abs.str(80, more=True),
        "maximum_first_current_modulus_ball": c0.str(80, more=True),
        "maximum_second_current_modulus_ball": c1.str(80, more=True),
        "uniform_two_current_remainder_ball": remainder.str(80, more=True),
        "first_omitted_third_current_modulus_ball": third_current.str(80, more=True),
        "last_exact_block_mode_corner_q": "q_A(39936,1/2)=-167/2=-83.5",
        "translated_edge_corner_q": "q_A(39936.5,1/2)=-169/2=-84.5",
        "adjacency": "The translated edge begins exactly one q-unit beyond the last exact outside-target A mode.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    numerical = artifact["numerical_certificate"]
    return f"""# Translated A-face two-current expansion

Date: 2026-08-23

Status: exact certificate for the first two normal currents and uniform
Fresnel-tail remainder; signed edge/block Morse assembly remains open

Put

```text
d_A=39936.5,
delta_A(x)=d_A-A*x/2,
q_A(x)=-delta_A(x)sqrt(2/x),
a_A(x)=d_A sqrt(2/x),              0<x<=1/2.          (TC1)
```

Since `d_A-A/4=169/4`, `|q_A|` decreases to its endpoint minimum

```text
q_A(x)<=-169/2=-84.5.                                (TC2)
```

For `q<0`, retain the exact lower Fresnel tail

```text
J_-(q)=F(q)+(1+i)/2
      =integral_(-infinity)^q exp(i*pi*v^2/2)dv.      (TC3)
```

Repeated integration by parts gives

```text
J_-(q)
 =exp(i*pi*q^2/2)
   [1/(i*pi*q)+1/((i*pi)^2*q^3)]+r_2(q),

|r_2(q)|<=6/(pi^3|q|^5).                             (TC4)
```

The bound follows by exposing the next `3/((i*pi)^3q^5)` boundary term and
bounding the remaining `v^-6` integral absolutely.  It is not a formal
asymptotic equality.

Use the continuous lower-endpoint current inherited from the exact finite
endpoint decomposition,

```text
L_A(d,x)=-exp(i*pi*q_A^2/2)/(i*pi)-a_A J_-(q_A).      (TC5)
```

The bare endpoint exponential and the first Fresnel-tail term cancel
algebraically.  With `E_A=exp(i*pi*q_A^2/2)`, the exact two-current form is

```text
L_A(d_A,x)=E_A[c_0(x)+c_1(x)]+R_2(x),

c_0(x)= A*x/[2*i*pi*delta_A(x)],
c_1(x)=-d_A*x/[2*pi^2*delta_A(x)^3],

|R_2(x)|
 <=3*d_A*x^2/[2*pi^3*delta_A(x)^5].                  (TC6)
```

All three moduli in (TC6) increase on the half-domain.  Their certified
endpoint maxima are

```text
|c_0| <= {numerical['maximum_first_current_modulus_ball']}
|c_1| <= {numerical['maximum_second_current_modulus_ball']}
|R_2| <= {numerical['uniform_two_current_remainder_ball']}. (TC7)
```

Quadratic completion retains the exact A-face phase:

```text
exp(-i*pi*d_A^2/x)E_A
 =exp(i*pi*A^2*x/4-i*pi*d_A*A).
```

Here `A` and the numerator `79873` are both `1 mod 4`, so

```text
exp(-i*pi*d_A*A)=-i.                                  (TC8)
```

Thus both displayed currents and the remainder enter the same tangential
phase `-i exp(i Phi_A(x))`; no fitted phase alignment is being used.

The finite-block boundary is exact:

```text
q_A(39936,1/2)=-83.5,
q_A(d_A,1/2)=-84.5.                                  (TC9)
```

The translated tail starts one normal-coordinate unit beyond the final
exact mode `39936`.  Equations (TC6)--(TC9) must therefore be joined to the
already-certified modes `39895..39936` before norms.  The small remainder in
(TC7) does not make the leading current small; its endpoint modulus is about
`{arb(numerical['maximum_first_current_modulus_ball']).mid().str(8)}`.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (TC1)--(TC9) is inherited from the exact
Kummer/Fourier phase and canonical Fresnel primitive.  No fitted or
geometric occurrence is introduced.

Proof boundary: exact continuous-edge Fresnel expansion, first two normal
currents, uniform pointwise remainder, common A-face phase, and roster
adjacency only.  No signed edge/42-mode cancellation bound, tangential Morse
integral, quantitative `R_after_A`, `R_Dir`, or `Q_K-T` estimate, all-height
theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["normal_tangential_split"]["decision"]["normal_then_tangential_analysis_licensed"] is True, "normal/tangential route drift")
    require(dependencies["finite_endpoint_tail"]["decision"]["finite_current_decomposed_exactly"] is True, "finite endpoint-tail definition drift")
    require(dependencies["A_projector_assembly"]["decision"]["projector_completed_exact_A_transition_certified"] is True, "A projector assembly drift")

    symbolic = symbolic_certificate()
    numerical = numerical_certificate()
    artifact = {
        "kind": STEM,
        "status": "exact_translated_A_face_two_normal_currents_and_uniform_Fresnel_tail_remainder_certified_signed_assembly_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "A": A,
            "translated_upper_edge": str(D),
            "last_exact_outside_target_A_mode": LAST_SHIFTED_MODE,
            "half_domain": "0<x<=1/2",
        },
        "symbolic_certificate": symbolic,
        "numerical_certificate": numerical,
        "decision": {
            "continuous_edge_A_tail_defined_by_exact_finite_endpoint_current": True,
            "first_two_normal_currents_derived_exactly": True,
            "uniform_abs_q_at_least_84_point_5_certified": True,
            "two_current_pointwise_remainder_below_4e_minus_6": True,
            "half_integer_phase_equals_minus_i": True,
            "translated_edge_adjacent_to_exact_42_mode_block": True,
            "leading_current_small_enough_for_separate_norm": False,
            "signed_edge_42_mode_cancellation_bound_proved": False,
            "tangential_Morse_integral_bounded": False,
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
        "next_obligation": "Write the fixed-regulator signed reassembly of c_0+c_1+R_2 with the exact outside-target A modes 39895..39936. Use the common -i exp(i Phi_A) phase to cancel the leading current before any modulus. Then map the residual amplitude to the existing exact A-face Morse coordinate and bound or scout only the resulting signed one-dimensional object.",
        "proof_boundary": "Exact continuous-edge Fresnel expansion, first two normal currents, uniform pointwise remainder, common A-face phase, and roster adjacency only. No signed edge/42-mode cancellation bound, tangential Morse integral, quantitative R_after_A, R_Dir, or Q_K-T estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified translated A-face two-current expansion", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
