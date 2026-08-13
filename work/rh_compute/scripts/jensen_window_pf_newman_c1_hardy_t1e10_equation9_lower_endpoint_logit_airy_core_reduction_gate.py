#!/usr/bin/env python3
"""Certify the exact lower-endpoint logit fold and its Airy core scaling."""

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


MORSE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json"
FRESNEL_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json"
PAPER = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/Lewis_Brereton_2607.15310.pdf"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 110
T = 10_000_000_000
A = 159577
CORE_Z = 4


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
    u, t, C, m, z = sp.symbols("u t C m z", real=True, positive=True)
    pi = sp.pi
    x = 1 / (1 + sp.exp(u))
    eta = pi * C**2 / (8 * t)
    mu = eta - 1
    phase = pi * C**2 * x / 4 - pi * m * C + t * u / 2
    constant = pi * C**2 / 8 - pi * m * C
    exact_fold = t * (u / 2 - eta * sp.tanh(u / 2))
    require(sp.simplify(sp.expand_trig(phase - constant - exact_fold).rewrite(sp.exp)) == 0, "exact logit fold phase failed")

    derivative = sp.factor(sp.diff(exact_fold / t, u))
    expected_derivative = (1 - eta * sp.sech(u / 2) ** 2) / 2
    require(sp.simplify(sp.expand_trig(derivative - expected_derivative).rewrite(sp.exp)) == 0, "fold stationary equation failed")
    u_star = 2 * sp.acosh(sp.sqrt(eta))
    require(sp.simplify(expected_derivative.subs(u, u_star)) == 0, "positive fold root failed")

    scale = 2 / (t * eta) ** sp.Rational(1, 3)
    airy_lambda = mu * t ** sp.Rational(2, 3) * eta ** (-sp.Rational(1, 3))
    cubic = -t * mu * u / 2 + t * eta * u**3 / 24
    require(sp.simplify(sp.powsimp(cubic.subs(u, scale * z) - (z**3 / 3 - airy_lambda * z), force=True)) == 0, "Airy scaling failed")

    measure = sp.simplify((-sp.diff(x, u)) * (x * (1 - x)) ** (-sp.Rational(1, 4)))
    expected_measure = 2 ** (-sp.Rational(3, 2)) * sp.cosh(u / 2) ** (-sp.Rational(3, 2))
    require(sp.simplify(sp.expand_trig(measure - expected_measure).rewrite(sp.exp)) == 0, "logit Kummer measure failed")

    lower_mode_factor = sp.exp(sp.I * pi * m) * sp.exp(-sp.I * pi * m * C)
    parity_statement = "For integer m and odd integer C, (-1)^m exp(-i*pi*m*C)=1."

    return {
        "logit_coordinate": "u=log((1-x)/x), x=1/(1+exp(u))",
        "lower_boundary_phase": "Phi_C(u;m,t)=pi*C^2/[4(1+e^u)]-pi*m*C+(t/2)u",
        "phase_constant": "pi*C^2/8-pi*m*C",
        "exact_fold_phase": "Phi_C-constant=t[u/2-eta*tanh(u/2)], eta=pi*C^2/(8t)",
        "stationary_equation": "sech^2(u/2)=1/eta",
        "stationary_points": "u_+-=+/-2 acosh(sqrt(eta)) for eta>=1",
        "fold_value": "eta=1 gives t[u/2-tanh(u/2)]=t[u^3/24-u^5/240+... ]",
        "control_parameter": "mu=eta-1",
        "Airy_scale": "u=2z/(t*eta)^(1/3)",
        "Airy_parameter": "lambda=mu*t^(2/3)*eta^(-1/3)",
        "canonical_cubic": "t[-mu*u/2+eta*u^3/24]=z^3/3-lambda*z",
        "Taylor_guard": "Using |tanh^(5)(y)|<=512 for real y, |tanh(u/2)-u/2+u^3/24|<=(2/15)|u|^5.",
        "scaled_phase_remainder": "|R_phase(z)|<=(64/15)t^(-2/3)eta^(-2/3)|z|^5",
        "Kummer_measure": "[x(1-x)]^(-1/4)|dx|=2^(-3/2)cosh(u/2)^(-3/2)du",
        "lower_endpoint_mode_factor": str(lower_mode_factor),
        "integer_parity_guard": parity_statement,
        "series_grouping_guard": "The separated lower-boundary current is mode-independent after parity and is not summable over modes; the Airy phase must be used only inside the exact grouped boundary/Fresnel current.",
    }


def certified_core() -> dict[str, Any]:
    ctx.dps = PRECISION
    t = arb(T)
    c = arb(A)
    pi = arb.pi()
    eta = pi * c**2 / (8 * t)
    mu = eta - 1
    scale = 2 / (t * eta) ** (arb(1) / 3)
    airy_lambda = mu * t ** (arb(2) / 3) * eta ** (-arb(1) / 3)
    u_star = 2 * eta.sqrt().acosh()
    x_lower = 1 / (1 + u_star.exp())
    x_upper = 1 / (1 + (-u_star).exp())
    n_lower = c * x_lower / 2
    n_upper = c * x_upper / 2
    remainder_bound = arb(64) / 15 * t ** (-arb(2) / 3) * eta ** (-arb(2) / 3) * arb(CORE_Z) ** 5
    require(eta > 1 and mu > 0, "source fold is not in two-saddle regime")
    require(arb(5) < airy_lambda < arb(6), "Airy parameter classification failed")
    require(remainder_bound < arb("0.001"), "Airy core phase remainder exceeds 0.001")
    require(arb(39852) < n_lower < arb(39853), "lower continuous root drift")
    require(arb(39936) < n_upper < arb(39937), "upper continuous root drift")
    return {
        "precision_decimal_digits": PRECISION,
        "height": T,
        "lower_endpoint": A,
        "eta_ball": eta.str(PRECISION, more=True),
        "mu_ball": mu.str(PRECISION, more=True),
        "Airy_scale_ball": scale.str(PRECISION, more=True),
        "Airy_lambda_ball": airy_lambda.str(PRECISION, more=True),
        "positive_stationary_u_ball": u_star.str(PRECISION, more=True),
        "lower_stationary_x_ball": x_lower.str(PRECISION, more=True),
        "upper_stationary_x_ball": x_upper.str(PRECISION, more=True),
        "lower_continuous_mode_ball": n_lower.str(PRECISION, more=True),
        "upper_continuous_mode_ball": n_upper.str(PRECISION, more=True),
        "Airy_core_z_radius": CORE_Z,
        "Airy_core_u_radius_ball": (scale * CORE_Z).str(PRECISION, more=True),
        "conservative_phase_remainder_bound_ball": remainder_bound.str(PRECISION, more=True),
        "decision": "The exact lower-boundary phase is a two-saddle fold with Airy parameter about 5.11. On |z|<=4 its noncubic phase is rigorously below 0.001 under a conservative derivative bound. The remaining work is grouped amplitude control and contour/tail completion, not discovery of the canonical fold phase.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certified_core"]
    return f"""# Exact lower-endpoint logit fold and Airy core

Date: 2026-08-10

Status: exact fold phase and certified Airy core validated; not a proof of the grouped fold remainder

Restrict the joint equation-(9) phase to the lower alpha endpoint `alpha=C`
and introduce the logit coordinate

```text
u=log((1-x)/x),       x=1/(1+exp(u)),
eta=pi*C^2/(8t).
```

After removing the constant `pi*C^2/8-pi*m*C`, the phase is exactly

```text
Phi_C(u)-Phi_C(0)=t[u/2-eta*tanh(u/2)].                 (LA1)
```

This is independent of the Poisson mode except for the removed constant.  Its
stationary equation is

```text
sech^2(u/2)=1/eta.                                     (LA2)
```

At the fold `eta=1`,

```text
t[u/2-tanh(u/2)]=t[u^3/24-u^5/240+...].                (LA3)
```

Put `mu=eta-1` and scale

```text
u=2z/(t*eta)^(1/3),
lambda=mu*t^(2/3)*eta^(-1/3).                          (LA4)
```

The linear and cubic phase becomes the canonical Airy polynomial exactly:

```text
t[-mu*u/2+eta*u^3/24]=z^3/3-lambda*z.                 (LA5)
```

A global real derivative bound `|tanh^(5)|<=512` gives

```text
|R_phase(z)|
 <=(64/15)t^(-2/3)eta^(-2/3)|z|^5.                    (LA6)
```

For `t=10^10`, `C=159577`, interval arithmetic certifies

```text
eta    ={c['eta_ball']},
mu     ={c['mu_ball']},
lambda ={c['Airy_lambda_ball']}.
```

The two exact logit stationary points reproduce the prior portcullis roots:

```text
N_-={c['lower_continuous_mode_ball']},
N_+={c['upper_continuous_mode_ball']}.
```

On the natural core `|z|<={c['Airy_core_z_radius']}`, the conservative phase
remainder satisfies

```text
|R_phase|<={c['conservative_phase_remainder_bound_ball']}<0.001. (LA7)
```

The transformed Kummer measure is also exact and regular:

```text
[x(1-x)]^(-1/4)|dx|
 =2^(-3/2)cosh(u/2)^(-3/2)du.                          (LA8)
```

There is a crucial grouping guard.  Since `C` is odd,

```text
(-1)^m exp(-i*pi*m*C)=1                                (LA9)
```

for every integer mode.  The lower endpoint current is therefore
mode-independent after parity and cannot be summed separately.  Equations
(LA1)--(LA8) must be applied inside the exact grouped boundary/Fresnel current,
where that apparent divergence cancels.

The canonical fold phase is no longer the unknown.  The next obligation is
to transport the full grouped amplitude through the Airy core, bound its
variation and the `|z|>4` contour tails explicitly, and join the result to the
ordinary interior and nonstationary ranges.

Proof boundary: exact lower-boundary fold algebra and a finite certified Airy
phase core at `t=10^10`.  No grouped fold-uniform remainder, complete
`T_upper` assembly, source-aligned height-uniform error, `Lambda<=0`, RH, or
prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(MORSE_GATE.is_file() and FRESNEL_GATE.is_file() and PAPER.is_file() and CHECKER.is_file(), "missing dependency or checker")
    symbolic = symbolic_reduction()
    core = certified_core()
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_endpoint_logit_airy_core_reduction_gate",
        "status": "exact_lower_endpoint_logit_fold_and_certified_Airy_phase_core_complete",
        "passed": True,
        "scope": {"height_center": "1e10", "lower_endpoint": A, "Airy_core": f"|z|<={CORE_Z}", "diagnostic_midpoint_used": False},
        "symbolic_reduction": symbolic,
        "certified_core": core,
        "decision": {
            "exact_logit_fold_phase_proved": True,
            "canonical_Airy_scaling_proved": True,
            "Airy_core_phase_remainder_below_0_001": True,
            "portcullis_roots_recovered": True,
            "separated_endpoint_mode_series_summable": False,
            "grouped_amplitude_remainder_proved": False,
            "rh_implication": False,
        },
        "next_obligation": "Transport the exact grouped boundary/Fresnel amplitude through the |z|<=4 Airy core, derive explicit variation and contour-tail constants uniformly in the fold parameter, and join it to ordinary interior and nonstationary mode estimates.",
        "proof_boundary": "Exact lower-boundary fold algebra and finite Airy phase-core bound at t=10^10 only. No grouped fold-uniform remainder, complete T_upper assembly, source-aligned height-uniform error, Lambda<=0, RH, or prize-level conclusion is proved.",
        "dependencies": {
            "morse_fold_gate": {"path": relative(MORSE_GATE), "sha256": file_hash(MORSE_GATE)},
            "fresnel_gate": {"path": relative(FRESNEL_GATE), "sha256": file_hash(FRESNEL_GATE)},
            "paper": {"path": relative(PAPER), "sha256": file_hash(PAPER)},
        },
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"workers": 1, "sympy_threads": 1, "flint_threads": 1, "process_priority": priority, "elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("built lower-endpoint logit Airy core: lambda=5.1099, roots=2, phase<0.001")


if __name__ == "__main__":
    main()
