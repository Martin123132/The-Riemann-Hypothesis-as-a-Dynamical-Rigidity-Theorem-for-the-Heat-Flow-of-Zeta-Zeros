#!/usr/bin/env python3
"""Certify the finite Airy-Fresnel detuning ODE for the lower fold core."""

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

from flint import arb, acb, ctx
import sympy as sp


NORMAL_FORM_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
CORE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate.json"
LOWER_FOLD_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_airy_fresnel_detuning_ode_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_airy_fresnel_detuning_ode_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 70
T = 10_000_000_000
C = 159577
MODE_LO = 39853
MODE_HI = 39936
Y = 64
SAMPLE_MODES = (39853, 39894, 39936)
ABS_TOL = "1e-24"


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


def symbolic_ode() -> dict[str, str]:
    y, d, a, beta, lam = sp.symbols("y d a beta lambda", real=True, positive=True)
    i = sp.I
    w = sp.exp(i * (a * y**2 - d * y))
    w_second_ratio = sp.simplify(sp.diff(w, y, 2) / w)
    expected = 2 * i * a - (2 * a * y - d) ** 2
    require(sp.simplify(w_second_ratio - expected) == 0, "chirp second derivative failed")

    coefficient_second = sp.simplify(4 * a**2).subs(a, 1 / (4 * beta))
    coefficient_first = sp.simplify(i * (1 + 4 * a * d)).subs(a, 1 / (4 * beta))
    coefficient_zero = sp.simplify(lam - d**2 + 2 * i * a).subs(a, 1 / (4 * beta))
    require(coefficient_second == 1 / (4 * beta**2), "second derivative coefficient failed")
    require(coefficient_first == i * (1 + d / beta), "transport coefficient failed")
    require(coefficient_zero == lam - d**2 + i / (2 * beta), "potential coefficient failed")
    return {
        "finite_transform": "G_Y(d)=int_0^Y Ai(-lambda-y) exp(i*y^2/(4beta)-i*d*y)dy",
        "ode": "G_Y''/(4beta^2)+i(1+d/beta)G_Y'+(lambda-d^2+i/(2beta))G_Y=R_0(d)+R_Y(d)",
        "lower_source": "R_0(d)=-Ai'(-lambda)+i*d*Ai(-lambda)",
        "upper_source": "R_Y(d)=exp(i[Y^2/(4beta)-dY]){Ai'(-lambda-Y)+i[Y/(2beta)-d]Ai(-lambda-Y)}",
        "moment_relations": "G_Y'=-i int_0^Y y f_d(y)dy; G_Y''=-int_0^Y y^2 f_d(y)dy",
        "derivation": "Use Ai(-lambda-y)''+(lambda+y)Ai(-lambda-y)=0 and integrate its second derivative twice by parts without separating either finite endpoint source.",
    }


def parameters() -> dict[str, arb]:
    ctx.dps = PRECISION
    pi = arb.pi()
    t = arb(T)
    endpoint_c = arb(C)
    eta = pi * endpoint_c**2 / (8 * t)
    beta = (t * eta) ** (arb(1) / 3)
    lam = (eta - 1) * t / beta
    h = 4 * beta / endpoint_c
    d_min = h * (arb(MODE_LO) - endpoint_c / 4)
    d_max = h * (arb(MODE_HI) - endpoint_c / 4)
    return {
        "pi": pi,
        "beta": beta,
        "lambda": lam,
        "h": h,
        "d_min": d_min,
        "d_max": d_max,
    }


def moment_integral(d: arb, power: int, p: dict[str, arb]) -> acb:
    i = acb(0, 1)
    beta = p["beta"]
    lam = p["lambda"]

    def integrand(y: acb, _: bool) -> acb:
        airy = (-(lam + y)).airy_ai()
        phase = i * (y**2 / (4 * beta) - d * y)
        return y**power * airy * phase.exp()

    return acb.integral(
        integrand,
        arb(0),
        arb(Y),
        abs_tol=arb(ABS_TOL),
        rel_tol=arb(ABS_TOL),
        eval_limit=600_000,
        depth_limit=60,
    )


def ode_residual(mode: int, p: dict[str, arb]) -> dict[str, Any]:
    i = acb(0, 1)
    beta = p["beta"]
    lam = p["lambda"]
    d = p["h"] * (arb(mode) - arb(C) / 4)
    moments = [moment_integral(d, power, p) for power in range(3)]
    g, g_prime, g_second = moments[0], -i * moments[1], -moments[2]
    lhs = (
        g_second / (4 * beta**2)
        + i * (1 + d / beta) * g_prime
        + (lam - d**2 + i / (2 * beta)) * g
    )
    ai0, ai0p, _, _ = (-lam).airy()
    aiy, aiyp, _, _ = (-(lam + Y)).airy()
    upper_phase = (i * (arb(Y) ** 2 / (4 * beta) - d * Y)).exp()
    rhs = -ai0p + i * d * ai0 + upper_phase * (
        aiyp + i * (arb(Y) / (2 * beta) - d) * aiy
    )
    residual = lhs - rhs
    require(residual.real.contains(0) and residual.imag.contains(0), f"ODE residual excludes zero at mode {mode}")
    return {
        "mode": mode,
        "detuning_ball": d.str(PRECISION, more=True),
        "G_Y": complex_record(g),
        "G_Y_prime": complex_record(g_prime),
        "G_Y_second": complex_record(g_second),
        "lhs": complex_record(lhs),
        "rhs": complex_record(rhs),
        "residual": complex_record(residual),
        "residual_abs_ball": abs(residual).str(PRECISION, more=True),
    }


def certified_coefficients() -> dict[str, Any]:
    p = parameters()
    beta = p["beta"]
    lam = p["lambda"]
    d_min = p["d_min"]
    d_max = p["d_max"]
    epsilon = 1 / (4 * beta**2)
    transport_min = 1 + d_min / beta
    transport_max = 1 + d_max / beta
    detuning_abs_max = max(abs(d_min), abs(d_max))
    potential_min = lam - detuning_abs_max**2
    ai0, ai0p, _, _ = (-lam).airy()
    aiy, aiyp, _, _ = (-(lam + Y)).airy()
    require(epsilon < arb("5.4e-8"), "ODE dispersive coefficient exceeds 5.4e-8")
    require(transport_min > arb("0.998"), "ODE transport coefficient approaches zero")
    require(transport_max < arb("1.002"), "ODE transport coefficient exceeds 1.002")
    require(potential_min > arb("0.026"), "ODE real potential loses the positive mode margin")
    residuals = [ode_residual(mode, p) for mode in SAMPLE_MODES]
    return {
        "height": T,
        "Y": Y,
        "mode_range": [MODE_LO, MODE_HI],
        "mode_count": MODE_HI - MODE_LO + 1,
        "beta_ball": beta.str(PRECISION, more=True),
        "lambda_ball": lam.str(PRECISION, more=True),
        "detuning_min_ball": d_min.str(PRECISION, more=True),
        "detuning_max_ball": d_max.str(PRECISION, more=True),
        "second_derivative_coefficient_ball": epsilon.str(PRECISION, more=True),
        "transport_coefficient_min_ball": transport_min.str(PRECISION, more=True),
        "transport_coefficient_max_ball": transport_max.str(PRECISION, more=True),
        "real_potential_min_ball": potential_min.str(PRECISION, more=True),
        "Ai_minus_lambda_ball": ai0.str(PRECISION, more=True),
        "Ai_prime_minus_lambda_ball": ai0p.str(PRECISION, more=True),
        "Ai_minus_lambda_minus_Y_ball": aiy.str(PRECISION, more=True),
        "Ai_prime_minus_lambda_minus_Y_ball": aiyp.str(PRECISION, more=True),
        "rigorous_moment_residuals": residuals,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_coefficients"]
    return f"""# Finite Airy-Fresnel detuning ODE for the lower fold

Date: 2026-08-11

Status: exact finite-core detuning ODE and representative rigorous residuals
validated; not a proof of the fold-uniform source bound

Define the finite lower-fold transform

```text
G_Y(d)=integral_0^Y Ai(-lambda-y)
       exp(i*y^2/(4beta)-i*d*y)dy.                     (DO1)
```

The Airy equation and two integrations by parts give the exact inhomogeneous
ODE

```text
G_Y''/(4beta^2)+i(1+d/beta)G_Y'
 +(lambda-d^2+i/(2beta))G_Y=R_0(d)+R_Y(d),             (DO2)

R_0(d)=-Ai'(-lambda)+i*d*Ai(-lambda),                 (DO3)

R_Y(d)=exp(i[Y^2/(4beta)-dY])
 {{Ai'(-lambda-Y)+i[Y/(2beta)-d]Ai(-lambda-Y)}}.       (DO4)
```

Both endpoint currents remain explicit in (DO3)--(DO4); no endpoint or
Fresnel component is frozen or discarded.  For `Y=64` and all 84 detunings,

```text
1/(4beta^2)={c['second_derivative_coefficient_ball']}<5.4e-8, (DO5)

{c['transport_coefficient_min_ball']}
 <1+d/beta<
{c['transport_coefficient_max_ball']},                 (DO6)

lambda-d^2>={c['real_potential_min_ball']}>0.026.       (DO7)
```

Thus the new equation is a nondegenerate first-order detuning transport with
a small but retained second-derivative term; it does not inherit the
`1/Delta` singularity of the frozen one-dimensional saddle expansion.

The identities

```text
G_Y'=-i integral_0^Y y f_d(y)dy,
G_Y''=-integral_0^Y y^2 f_d(y)dy                       (DO8)
```

were evaluated by rigorous complex ball quadrature at modes `39853`,
`39894`, and `39936`.  At all three representatives, the real and imaginary
balls of the independently assembled left side minus right side contain
zero.

The exact grouped compact-core value is `2pi*sum_m G_64(d_m)`, so (DO2)
provides a detuning-space route to all 84 modes while preserving their finite
sum.  The next obligation is an interval ODE/energy estimate or a sourced
variation-of-constants formula that bounds this complete discrete sum and
matches its outer-Morse limit.

Pi provenance: the chirp in (DO1) comes from the Kummer quadratic and the
detuning Fourier factor comes from finite Poisson summation.  No fitted
geometric normalization is introduced.

Proof boundary: exact canonical finite-core ODE, coefficient margins, and
three rigorous saved-height residual checks only.  No all-detuning ODE
enclosure, outer-chart join, complete `T_upper`, height-uniform source
theorem, `Lambda<=0`, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for dependency in (NORMAL_FORM_GATE, CORE_GATE, LOWER_FOLD_GATE):
        require(dependency.is_file(), f"missing dependency: {dependency}")
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_airy_fresnel_detuning_ode_gate",
        "status": "exact_finite_lower_fold_Airy_Fresnel_detuning_ODE_complete",
        "passed": True,
        "symbolic_ode": symbolic_ode(),
        "certified_coefficients": certified_coefficients(),
        "decision": {
            "exact_finite_detuning_ODE_proved": True,
            "both_finite_endpoint_sources_retained": True,
            "transport_coefficient_nonzero_on_84_mode_interval": True,
            "real_potential_positive_on_84_mode_interval": True,
            "three_rigorous_moment_residuals_overlap_zero": True,
            "all_detuning_ODE_enclosure_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            "normal_form_gate": {"path": relative(NORMAL_FORM_GATE), "sha256": file_hash(NORMAL_FORM_GATE)},
            "core_gate": {"path": relative(CORE_GATE), "sha256": file_hash(CORE_GATE)},
            "lower_fold_gate": {"path": relative(LOWER_FOLD_GATE), "sha256": file_hash(LOWER_FOLD_GATE)},
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
        "next_obligation": "Construct an interval variation-of-constants or energy estimate for the finite ODE across the full detuning interval, sum its values on the 84-point lattice without termwise absolute loss, and match the endpoint data to the ordinary Morse carrier.",
        "proof_boundary": "Exact finite canonical detuning ODE and representative saved-height residual checks only. No full-interval ODE enclosure, outer join, T_upper theorem, Lambda<=0, RH, or prize-level conclusion is proved.",
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("built lower-fold detuning ODE: transport nonzero, potential>0.026, 3 residuals certified")


if __name__ == "__main__":
    main()
