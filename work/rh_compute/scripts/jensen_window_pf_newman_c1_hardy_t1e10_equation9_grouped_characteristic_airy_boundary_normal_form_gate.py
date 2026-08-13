#!/usr/bin/env python3
"""Certify the grouped characteristic Airy-boundary normal form."""

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


AIRY_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.json"
FRESNEL_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 120
T = 10_000_000_000
C = 159577
B = 5122421
MODE_LO = 39853
MODE_HI = 39936
Z_CORE = 4
Y_CORE = 64


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


def symbolic_reduction() -> dict[str, Any]:
    z, y, beta, c, m, t, eta = sp.symbols(
        "z y beta C m t eta", real=True, positive=True
    )
    pi = sp.pi
    tau = sp.tanh(z / beta)
    x = (1 - tau) / 2
    sigma = 4 * beta / (pi * c)
    d = 4 * beta * (m - c / 4) / c
    s = sigma * y

    linear = sp.expand(pi * (x * c / 2 - m) * s)
    expected_linear = -d * y - beta * y * tau
    require(sp.simplify(linear - expected_linear) == 0, "normal linear phase failed")

    quadratic_ratio = sp.simplify((pi * sigma**2 / 4) / (1 / (2 * beta)))
    require(
        sp.simplify(quadratic_ratio - 8 * beta**3 / (pi * c**2)) == 0,
        "normal quadratic phase ratio failed",
    )

    measure = 2 ** (-sp.Rational(3, 2)) * sp.cosh(z / beta) ** (-sp.Rational(3, 2))
    transformed_amplitude = sp.simplify(
        measure * (2 / beta) * sp.Rational(1, 2) * (c + s) * sigma
    )
    expected_amplitude = (
        sp.sqrt(2)
        / pi
        * sp.cosh(z / beta) ** (-sp.Rational(3, 2))
        * (1 + sigma * y / c)
    )
    require(
        sp.simplify(transformed_amplitude - expected_amplitude) == 0,
        "joint Jacobian-amplitude failed",
    )

    airy_lambda = (eta - 1) * t / beta
    endpoint = t * (z / beta - eta * tau)
    canonical = z**3 / 3 - airy_lambda * z
    endpoint_remainder = sp.expand(endpoint - canonical)
    exact_joint = endpoint - d * y - beta * y * tau + x * y**2 / (2 * beta)
    canonical_joint = canonical - (z + d) * y + y**2 / (4 * beta)
    joint_remainder = sp.simplify(exact_joint - canonical_joint)
    expected_remainder = sp.simplify(
        endpoint_remainder
        + y * (z - beta * tau)
        + (x - sp.Rational(1, 2)) * y**2 / (2 * beta)
    )
    require(sp.simplify(joint_remainder - expected_remainder) == 0, "joint remainder failed")

    return {
        "constraints": "eta=pi*C^2/(8t), beta=(t*eta)^(1/3), hence pi*C^2=8*beta^3",
        "coordinates": "u=2z/beta, alpha=C+sigma*y, sigma=4beta/(pi*C)",
        "logit": "x(z)=[1-tanh(z/beta)]/2",
        "mode_detuning": "d_m=4beta(m-C/4)/C",
        "exact_joint_phase": "t[z/beta-eta*tanh(z/beta)]-d_m*y-beta*y*tanh(z/beta)+x(z)y^2/(2beta)",
        "canonical_joint_phase": "z^3/3-lambda*z-(z+d_m)y+y^2/(4beta)",
        "Airy_parameter": "lambda=(eta-1)t/beta",
        "exact_phase_remainder": str(joint_remainder),
        "exact_amplitude": "sqrt(2)/pi*cosh(z/beta)^(-3/2)*(1+sigma*y/C)",
        "mode_kernel": "D_M(y)=sum_(m=M0)^M1 exp(-i*d_m*y)",
        "Dirichlet_form": "D_M(y)=exp[-i(d_0+(L-1)h/2)y] sin(L*h*y/2)/sin(h*y/2), h=4beta/C",
        "removable_kernel_values": "D_M(2*pi*k/h)=L*exp[-i*d_0*2*pi*k/h]",
        "grouping_guard": "The finite mode sum is taken before absolute values. No endpoint current is separated from its compensating Fresnel integral.",
    }


def certified_core() -> dict[str, Any]:
    ctx.dps = PRECISION
    t = arb(T)
    c = arb(C)
    pi = arb.pi()
    eta = pi * c**2 / (8 * t)
    beta = (t * eta) ** (arb(1) / 3)
    sigma = 4 * beta / (pi * c)
    h = 4 * beta / c
    airy_lambda = (eta - 1) * t / beta

    modes = list(range(MODE_LO, MODE_HI + 1))
    detunings = [h * (arb(mode) - c / 4) for mode in modes]
    endpoint_remainder = (
        arb(64) / 15 * arb(Z_CORE) ** 5 / beta**2
    )
    coupling_remainder = (
        arb(Y_CORE) * arb(Z_CORE) ** 3 / (3 * beta**2)
    )
    quadratic_remainder = (
        arb(Z_CORE) * arb(Y_CORE) ** 2 / (4 * beta**2)
    )
    total_remainder = endpoint_remainder + coupling_remainder + quadratic_remainder

    relative_amplitude_upper = sigma * arb(Y_CORE) / c
    relative_amplitude_lower = 1 - (arb(Z_CORE) / beta).cosh() ** (-arb(3) / 2)
    relative_amplitude_error = max(relative_amplitude_upper, relative_amplitude_lower)

    mode_count = len(modes)
    first_kernel_zero = 2 * pi / (arb(mode_count) * h)
    kernel_period = 2 * pi / h
    y_max = arb(B - C) / sigma
    fresnel_y_scale = beta.sqrt()

    require((pi * c**2 / 8 - beta**3).contains(arb(0)), "beta scale identity failed")
    require(arb(5) < airy_lambda < arb(6), "Airy parameter drift")
    require(arb("0.05") < h < arb("0.06"), "mode detuning spacing drift")
    require(arb("-2.3") < detunings[0] < arb("-2.2"), "lower detuning drift")
    require(arb("2.2") < detunings[-1] < arb("2.3"), "upper detuning drift")
    require(total_remainder < arb("0.0022"), "joint phase error exceeds 0.0022")
    require(relative_amplitude_error < arb("0.00001"), "relative amplitude error exceeds 1e-5")
    require(arb(1) < arb(Y_CORE) / fresnel_y_scale < arb(2), "chosen y core lost Fresnel scale")
    require(mode_count == 84 and C % 2 == 1 and (B - C) % 2 == 0, "roster parity drift")

    return {
        "precision_decimal_digits": PRECISION,
        "height": T,
        "lower_endpoint": C,
        "upper_endpoint": B,
        "mode_range": [MODE_LO, MODE_HI],
        "mode_count": mode_count,
        "z_core": [-Z_CORE, Z_CORE],
        "y_core": [0, Y_CORE],
        "eta_ball": eta.str(PRECISION, more=True),
        "beta_ball": beta.str(PRECISION, more=True),
        "sigma_ball": sigma.str(PRECISION, more=True),
        "Airy_lambda_ball": airy_lambda.str(PRECISION, more=True),
        "mode_detuning_spacing_ball": h.str(PRECISION, more=True),
        "mode_detuning_min_ball": detunings[0].str(PRECISION, more=True),
        "mode_detuning_max_ball": detunings[-1].str(PRECISION, more=True),
        "Dirichlet_first_zero_ball": first_kernel_zero.str(PRECISION, more=True),
        "Dirichlet_period_ball": kernel_period.str(PRECISION, more=True),
        "full_y_span_ball": y_max.str(PRECISION, more=True),
        "Fresnel_y_scale_sqrt_beta_ball": fresnel_y_scale.str(PRECISION, more=True),
        "y_core_over_Fresnel_scale_ball": (arb(Y_CORE) / fresnel_y_scale).str(PRECISION, more=True),
        "endpoint_phase_remainder_bound_ball": endpoint_remainder.str(PRECISION, more=True),
        "normal_coupling_remainder_bound_ball": coupling_remainder.str(PRECISION, more=True),
        "quadratic_coefficient_remainder_bound_ball": quadratic_remainder.str(PRECISION, more=True),
        "total_joint_phase_remainder_bound_ball": total_remainder.str(PRECISION, more=True),
        "relative_amplitude_upper_error_bound_ball": relative_amplitude_upper.str(PRECISION, more=True),
        "relative_amplitude_lower_error_bound_ball": relative_amplitude_lower.str(PRECISION, more=True),
        "relative_amplitude_error_bound_ball": relative_amplitude_error.str(PRECISION, more=True),
        "decision": "The legal 84-mode grouping has an exact two-variable characteristic normal form. On |z|<=4 and 0<=y<=64, its cubic-quadratic canonical phase has error below 0.0022 and its positive Jacobian amplitude differs relatively from sqrt(2)/pi by less than 1e-5. The mode dependence is entirely the explicit finite Dirichlet kernel.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_core"]
    return f"""# Grouped characteristic Airy-boundary normal form

Date: 2026-08-10

Status: exact grouped two-variable normal form and compact-core error budget
validated; not a proof of a complete fold-uniform estimate

The previous lower-endpoint calculation found the correct cubic boundary
phase, but it did not justify treating the remaining grouped current as a
slow scalar amplitude.  Keep the endpoint and Fresnel current together and
sum the finite transition roster before taking absolute values.

Put

```text
eta=pi*C^2/(8t),          beta=(t*eta)^(1/3),
u=2z/beta,                x(z)=[1-tanh(z/beta)]/2,
alpha=C+sigma*y,          sigma=4beta/(pi*C),
d_m=4beta(m-C/4)/C.                                      (GB1)
```

After removing the common constant `pi*C^2/8`, oddness of `C` gives the exact
joint phase

```text
Theta_m(z,y)
 =t[z/beta-eta*tanh(z/beta)]
  -d_m*y-beta*y*tanh(z/beta)+x(z)y^2/(2beta).           (GB2)
```

The transformed Kummer measure, alpha current, and both Jacobians combine to

```text
A(z,y)=sqrt(2)/pi cosh(z/beta)^(-3/2)(1+sigma*y/C).     (GB3)
```

Thus there is no separated endpoint series.  For the 84 modes
`{MODE_LO}..{MODE_HI}`, all mode dependence is the exact finite kernel

```text
D_84(y)=sum_m exp(-i*d_m*y)
 =exp[-i(d_0+83h/2)y] sin(84hy/2)/sin(hy/2),
h=4beta/C.                                              (GB4)
```

At removable zeros of the denominator, (GB4) is interpreted by continuity.
This is the cancellation object that the separated endpoint/Fresnel formula
conceals.

The correct local canonical phase is two-dimensional:

```text
Theta_m^0(z,y)
 =z^3/3-lambda*z-(z+d_m)y+y^2/(4beta),
lambda=(eta-1)t/beta.                                   (GB5)
```

The exact remainder is

```text
R=R_Airy(z)+y[z-beta*tanh(z/beta)]
  +[x(z)-1/2]y^2/(2beta).                               (GB6)
```

Using `|a-tanh(a)|<=|a|^3/3`, `|tanh(a)|<=|a|`, and the
certified Airy remainder, the rectangle `|z|<={Z_CORE}`,
`0<=y<={Y_CORE}` satisfies

```text
|R_Airy| <= {c['endpoint_phase_remainder_bound_ball']},
|R-R_Airy| <=
  {c['normal_coupling_remainder_bound_ball']}
 +{c['quadratic_coefficient_remainder_bound_ball']},
|R| <= {c['total_joint_phase_remainder_bound_ball']} < 0.0022. (GB7)
```

On the same rectangle,

```text
|A/(sqrt(2)/pi)-1|
 <= {c['relative_amplitude_error_bound_ball']} < 1e-5. (GB8)
```

The numerical characteristic data are

```text
beta={c['beta_ball']},
sigma={c['sigma_ball']},
lambda={c['Airy_lambda_ball']},
h={c['mode_detuning_spacing_ball']},
{c['mode_detuning_min_ball']} <= d_m <=
{c['mode_detuning_max_ball']}.
```

The normal core reaches
`{c['y_core_over_Fresnel_scale_ball']}` times the natural Fresnel scale
`sqrt(beta)`.  This is enough to validate the local model but not enough to
discard the complement: for some `(z,m)` the normal saddle lies beyond this
rectangle.  The next obligation is therefore to bound the finite-kernel
canonical integral and make an endpoint/interior partition whose outer piece
contains every displaced normal saddle.  A one-dimensional Airy estimate
with a frozen or slowly varying amplitude is not admissible.

Pi provenance: the `pi` in `sigma`, `h`, and the phase is inherited from the
original Kummer factor `exp(i*pi*x*alpha^2/4)` and the integer
Fourier--Poisson character `exp(-i*pi*m*alpha)`.  The factor `sqrt(2)/pi` in
(GB3) is forced by those coordinates and their two Jacobians.  No circle,
polygon, or fitted geometric normalization is introduced here.

Proof boundary: exact coordinate algebra and a compact numerical error budget
at `t=10^10`.  No canonical-integral bound, outer normal tail, complete
`T_upper` assembly, height-uniform source error, `Lambda<=0`, RH, or
prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(AIRY_GATE.is_file() and FRESNEL_GATE.is_file() and PAPER.is_file() and CHECKER.is_file(), "missing dependency or checker")
    symbolic = symbolic_reduction()
    core = certified_core()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate",
        "status": "exact_grouped_characteristic_Airy_boundary_normal_form_and_compact_core_complete",
        "passed": True,
        "scope": {
            "height_center": "1e10",
            "lower_endpoint": C,
            "transition_modes": f"{MODE_LO}..{MODE_HI}",
            "characteristic_core": f"|z|<={Z_CORE}, 0<=y<={Y_CORE}",
        },
        "symbolic_reduction": symbolic,
        "certified_core": core,
        "decision": {
            "exact_joint_coordinate_reduction_proved": True,
            "finite_mode_Dirichlet_grouping_proved": True,
            "compact_joint_phase_error_below_0_0022": True,
            "compact_relative_amplitude_error_below_1e_5": True,
            "scalar_slow_Airy_amplitude_model_admissible": False,
            "two_variable_characteristic_model_required": True,
            "canonical_integral_bound_proved": False,
            "outer_normal_tail_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Bound the finite-Dirichlet-kernel canonical integral generated by (GB4)-(GB5), then construct an endpoint/interior partition whose outer component captures all displaced alpha saddles and joins to the ordinary Morse ranges.",
        "proof_boundary": "Exact grouped coordinate algebra and compact characteristic-core error budget at t=10^10 only. No canonical-integral bound, outer normal tail, complete T_upper assembly, height-uniform source error, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "Airy_core_gate": {"path": relative(AIRY_GATE), "sha256": file_hash(AIRY_GATE)},
            "Fresnel_gate": {"path": relative(FRESNEL_GATE), "sha256": file_hash(FRESNEL_GATE)},
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "workers": 1,
            "sympy_threads": 1,
            "flint_threads": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built grouped characteristic Airy-boundary normal form: modes=84, phase<0.0022, amplitude<1e-5")


if __name__ == "__main__":
    main()
