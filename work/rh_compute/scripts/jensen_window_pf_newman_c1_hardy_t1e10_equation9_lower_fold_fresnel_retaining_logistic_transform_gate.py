#!/usr/bin/env python3
"""Certify the exact Fresnel-retaining logistic-Morse transform."""

from __future__ import annotations

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

from flint import arb, ctx
import sympy as sp


RETENTION_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate.json"
PORTCULLIS_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json"
FRESNEL_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json"
MORSE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
C = 159_577
B = 5_122_423
MODE = 39_894


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


def symbolic_certificate() -> dict[str, str]:
    x, m, endpoint, t, v = sp.symbols("x m D t v", positive=True)
    pi = sp.pi
    q_endpoint = sp.sqrt(x / 2) * (endpoint - 2 * m / x)
    completed_phase = -pi * m**2 / x + pi * q_endpoint**2 / 2
    original_phase = pi * x * endpoint**2 / 4 - pi * m * endpoint
    require(sp.simplify(completed_phase - original_phase) == 0, "Fresnel completion identity failed")

    r = t / (2 * pi * m**2)
    x_v = 1 / (1 + r * v)
    psi_v = -pi * m**2 / x_v + t * sp.log(r * v) / 2
    psi_1 = -pi * m**2 * (1 + r) + t * sp.log(r) / 2
    target = -t * (v - 1 - sp.log(v)) / 2
    require(sp.simplify(sp.expand_log(psi_v - psi_1 - target, force=True)) == 0, "logistic phase identity failed")

    tau = pi * m * (endpoint - 2 * m)
    x_event = sp.simplify(2 * pi * m**2 / (tau + 2 * pi * m**2))
    require(sp.simplify(x_event - 2 * m / endpoint) == 0, "event saddle x identity failed")
    require(sp.simplify(q_endpoint.subs(x, x_event)) == 0, "event endpoint coordinate failed")

    return {
        "quadratic_completion": "E_m(D)=exp(-i*pi*m^2/x)*exp(i*pi*q_D^2/2)",
        "paired_current": (
            "P_m=[exp(i*pi*q_B^2/2)-exp(i*pi*q_C^2/2)]/(i*pi)"
            "+m*sqrt(2/x)[F(q_B)-F(q_C)]"
        ),
        "Poisson_mode": "J_m=(-1)^m*exp(-i*pi*m^2/x)*P_m/x",
        "logistic_coordinate": "r=t/(2*pi*m^2), v=(1-x)/(r*x), x=1/(1+r*v)",
        "Morse_coordinate": "s=sgn(v-1)*sqrt(2*(v-1-log(v)))",
        "exact_phase": "psi_m(x)-psi_m(x_m)=-t*s^2/4",
        "Jacobian": "-dx/ds=r*x^2*dv/ds, dv/ds=v*s/(v-1), (dv/ds)(0)=1",
        "event_regularization": "at t=tau_m, x_m=2m/C and q_C(m,x_m)=0",
    }


def prototype_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    pi = arb.pi()
    m = MODE
    tau = pi * m * (C - 2 * m)
    r = tau / (2 * pi * m**2)
    x_m = 1 / (1 + r)
    q_c = (x_m / 2).sqrt() * (C - 2 * m / x_m)
    q_b = (x_m / 2).sqrt() * (B - 2 * m / x_m)
    jacobian_weight = r * x_m * (x_m * (1 - x_m)) ** (-arb(1) / 4)
    Delta = (2 * pi * m**2 - tau) / (2 * pi * m**2 + tau)

    require(q_c.contains(0), "prototype q_C does not contain zero")
    require((x_m - arb(2 * m) / C).contains(0), "prototype x_m identity failed")
    require((Delta + arb(1) / C).contains(0), "prototype Delta identity failed")
    require(q_b > arb("2480000"), "upper endpoint Fresnel coordinate lost its large margin")
    require(arb("0.70") < jacobian_weight < arb("0.72"), "event Jacobian-amplitude prefactor is not regular")

    return {
        "mode": m,
        "event_height_ball": tau.str(PRECISION, more=True),
        "logistic_r_ball": r.str(PRECISION, more=True),
        "saddle_x_ball": x_m.str(PRECISION, more=True),
        "lower_endpoint_Fresnel_coordinate_ball": q_c.str(PRECISION, more=True),
        "upper_endpoint_Fresnel_coordinate_ball": q_b.str(PRECISION, more=True),
        "logistic_characteristic_ball": Delta.str(PRECISION, more=True),
        "finite_Jacobian_weight_ball": jacobian_weight.str(PRECISION, more=True),
        "lower_endpoint_current_at_event": "exp(i*pi*q_C^2/2)=1",
        "lower_endpoint_Fresnel_primitive_at_event": "F(q_C)=F(0)=0",
        "division_by_characteristic_defect": False,
    }


def render_note(artifact: dict[str, Any]) -> str:
    p = artifact["prototype_certificate"]
    return f"""# Exact Fresnel-retaining logistic-Morse transform

Date: 2026-08-12

Status: exact positive-mode transform certified; this is not a proof of an
Airy/logistic remainder comparison or 399-event propagation theorem

For a positive finite-Poisson mode, define

```text
q_D(m,x)=sqrt(x/2)(D-2m/x),       D in {{C,B}},
F'(q)=exp(i*pi*q^2/2),            F(0)=0,
```

and keep the endpoint current and Fresnel term paired:

```text
P_m(x)= [exp(i*pi*q_B^2/2)-exp(i*pi*q_C^2/2)]/(i*pi)
       +m sqrt(2/x)[F(q_B)-F(q_C)].                    (FL1)
```

Completing the square in the exact mode integral gives

```text
J_m(x)=(-1)^m exp(-i*pi*m^2/x) P_m(x)/x.               (FL2)
```

No approximation has entered (FL1)--(FL2).

Put

```text
r=t/(2*pi*m^2),       v=(1-x)/(r*x),
x=1/(1+r*v),
s=sgn(v-1)sqrt(2[v-1-log(v)]).                         (FL3)
```

For

```text
psi_m(x)=-pi*m^2/x+(t/2)log((1-x)/x),
```

the phase is globally exact:

```text
psi_m(x)-psi_m(x_m)=-t*s^2/4.                          (FL4)
```

If `K_m` denotes the Kummer-weighted `x` integral of `J_m`, reversing the
monotone `x(s)` orientation gives

```text
K_m=(-1)^m exp(i*psi_m(x_m))
    integral_R exp(-i*t*s^2/4) A_m(s) ds,              (FL5)

A_m(s)=r*x(s)*(dv/ds)*[x(s)(1-x(s))]^(-1/4) P_m(x(s)),
dv/ds=v*s/(v-1),             (dv/ds)|_(s=0)=1.         (FL6)
```

Equations (FL5)--(FL6) are the required Fresnel-retaining logistic chart.
The universal Morse phase is Gaussian, while the endpoint/Fresnel current
remains exact inside the amplitude and is summed before absolute values.

## Event regularity

At `tau_m=pi*m(C-2m)`,

```text
x_m=2m/C,       q_C(m,x_m)=0.                          (FL7)
```

For the nearest event, `m={p['mode']}`,

```text
x_m={p['saddle_x_ball']},
q_C={p['lower_endpoint_Fresnel_coordinate_ball']},
q_B={p['upper_endpoint_Fresnel_coordinate_ball']},
r*x_m*[x_m(1-x_m)]^(-1/4)
 ={p['finite_Jacobian_weight_ball']}.
```

Thus at the crossing, the lower endpoint contributes the finite values
`exp(i*pi*q_C^2/2)=1` and `F(q_C)=0`.  The apparent characteristic defect is

```text
Delta={p['logistic_characteristic_ball']},
```

but neither (FL1) nor (FL5)--(FL6) divides by it.

## Next comparison

The valid prototype target is now precise: retain mode `39894` exactly on
its `pi/16` buffer, and on the overlap compare the grouped Airy-Fresnel
integral with (FL5) using the same `P_m` current before taking absolute
values.  The full-interior replacement `F(q_B)-F(q_C) -> 1+i` is excluded by
the preceding Fresnel-retention gate.

## Boundary

This is an exact coordinate and current identity.  It does not bound the
variation of `A_m`, prove an Airy-to-logistic remainder, propagate through
399 events, establish complete `T_upper`, prove `Lambda<=0`, RH, or a
prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    dependencies = {
        "Fresnel_retention_gate": RETENTION_GATE,
        "portcullis_gate": PORTCULLIS_GATE,
        "Fresnel_normal_form_gate": FRESNEL_GATE,
        "logistic_Morse_gate": MORSE_GATE,
    }
    for path in (*dependencies.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    for path in dependencies.values():
        require(json.loads(path.read_text(encoding="utf-8")).get("passed") is True, f"dependency is not passed: {path}")

    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate",
        "status": "exact_Fresnel_retaining_logistic_Morse_transform_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "prototype_certificate": prototype_certificate(),
        "decision": {
            "endpoint_current_and_Fresnel_term_remain_paired": True,
            "global_logistic_Morse_phase_is_exact": True,
            "transform_is_regular_at_turning_event": True,
            "division_by_characteristic_defect_used": False,
            "Airy_to_logistic_remainder_proved": False,
            "propagation_across_399_events_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Derive a cancellation-preserving comparison between the grouped Airy-Fresnel transform and "
            "the exact Fresnel-retaining logistic amplitude A_m(s), beginning with mode 39894 at the two "
            "pi/16 faces. Bound amplitude/phase transport without freezing F(q_B)-F(q_C)."
        ),
        "proof_boundary": (
            "Exact algebraic transform and event regularity only. No amplitude-variation estimate, "
            "Airy-to-logistic remainder, 399-event propagation, complete T_upper, Lambda<=0, RH, or "
            "prize-level conclusion is proved."
        ),
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in dependencies.items()
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
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified exact Fresnel-retaining logistic transform: "
        f"mode={MODE}, q_C(event)=0, no characteristic division; priority={priority}"
    )


if __name__ == "__main__":
    main()
