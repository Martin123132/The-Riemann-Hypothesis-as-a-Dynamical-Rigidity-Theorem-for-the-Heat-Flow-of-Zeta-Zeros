#!/usr/bin/env python3
"""Conjugate the selected beta^-4 detuning ODE to an exact Airy Green kernel."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_detuning_airy_green_kernel_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "detuning_ode": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_beta4_endpoint_coherent_detuning_ode_gate.json",
    "branch_envelope": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_beta4_scaled_airy_branch_amplitude_envelope_gate.json",
}

C = 159_577
Q = (C - 1) // 4
MODE_LO = 39_696
MODE_HI = 40_094
PRECISION = 90
DLMF_URL = "https://dlmf.nist.gov/9.8"


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
    beta, lam, d = sp.symbols("beta lambda d", positive=True, real=True)
    i = sp.I
    chi = 2 * beta**2 * d + beta * d**2
    r = beta**2 + lam + 2 * beta * d
    h = sp.Function("H")
    g = sp.exp(-i * chi) * h(r)
    operator = sp.diff(g, d, 2) / (4 * beta**2)
    operator += i * (1 + d / beta) * sp.diff(g, d)
    operator += (lam - d**2 + i / (2 * beta)) * g
    expected = sp.exp(-i * chi) * (sp.diff(h(r), d, 2) / (4 * beta**2) + r * h(r))
    require(sp.simplify(operator - expected) == 0, "Airy gauge conjugation failed")

    c = sp.Integer(C)
    q = sp.Integer(Q)
    mode = sp.symbols("m", integer=True)
    d_mode = beta * (4 * mode / c - 1)
    r_mode = sp.simplify(beta**2 + lam + 2 * beta * d_mode)
    step = sp.simplify(r_mode.subs(mode, q + 1) - r_mode.subs(mode, q))
    require(step == 8 * beta**2 / c, "Airy lattice step failed")
    j = sp.symbols("j", integer=True)
    center = sp.simplify(r_mode.subs(mode, q))
    require(sp.simplify(r_mode.subs(mode, q + j) - center - j * step) == 0, "positive Airy lattice symmetry failed")
    require(sp.simplify(r_mode.subs(mode, q - j) - center + j * step) == 0, "negative Airy lattice symmetry failed")
    return {
        "gauge": "G(d)=exp[-i(2beta^2 d+beta d^2)] H(r)",
        "airy_coordinate": "r=beta^2+lambda+2beta d",
        "conjugated_equation": "H_rr+rH=exp[i(2beta^2 d+beta d^2)]S(d)",
        "green_kernel": "K(r,s)=pi[Ai(-r)Bi(-s)-Bi(-r)Ai(-s)]",
        "kernel_initial_data": "K(s,s)=0; partial_r K(s,s)=1",
        "variation_of_constants": (
            "H(r)=-partial_s K(r,r0)H(r0)+K(r,r0)H_r(r0)"
            "+Integral_(r0)^r K(r,s)exp(i chi(d(s)))S(d(s))ds"
        ),
        "mode_coordinate": str(r_mode),
        "lattice_step": str(step),
        "event_pair_symmetry": "r_(q+j)=r_q+j Delta_r; r_(q-j)=r_q-j Delta_r",
    }


def interval_certificate() -> dict[str, Any]:
    pi = arb.pi()
    c = arb(C)
    beta = (pi * c**2 / 8) ** (arb(1) / 3)
    lambda_max = pi / (16 * beta)
    r_min = beta**2 * (arb(8) * MODE_LO / c - 1)
    r_max = beta**2 * (arb(8) * MODE_HI / c - 1) + lambda_max
    step = 8 * beta**2 / c
    center_lower = beta**2 * (1 - 2 / c)
    center_upper = center_lower + lambda_max
    kernel_envelope = 1 / r_min.sqrt()
    require(r_min > arb("4.59e6"), "Airy Green coordinate approaches its turning point")
    require(r_max < arb("4.69e6"), "Airy Green coordinate exceeds certified range")
    require(kernel_envelope < arb("4.67e-4"), "Airy Green kernel envelope exceeds 4.67e-4")
    return {
        "height_interval": "t*-pi/16<=t<=t*",
        "lambda_interval": "0<=lambda<=pi/(16beta)",
        "mode_roster": [MODE_LO, MODE_HI],
        "removed_event_zero_mode": Q,
        "r_min_ball": r_min.str(PRECISION, more=True),
        "r_max_ball": r_max.str(PRECISION, more=True),
        "event_zero_r_interval": [center_lower.str(PRECISION, more=True), center_upper.str(PRECISION, more=True)],
        "exact_lattice_step_ball": step.str(PRECISION, more=True),
        "green_kernel_absolute_envelope_ball": kernel_envelope.str(PRECISION, more=True),
        "modulus_inequality": "M(-r)^2<=1/(pi sqrt(r)) for r>0 from DLMF 9.8.20 with n=0 and its signed first-neglected-term remainder",
        "kernel_bound_derivation": "|K(r,s)|<=pi M(-r)M(-s)<=(rs)^(-1/4)<=r_min^(-1/2)",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Airy Green kernel for the selected-branch detuning equation

Date: 2026-08-13

Status: exact gauge conjugation, variation-of-constants kernel, and uniform
pointwise kernel envelope proved; no grouped Green-operator bound

For the detuning equation (11.387.2), put

```text
chi(d)=2beta^2 d+beta d^2,
r=beta^2+lambda+2beta d,
G(d)=exp[-i chi(d)]H(r).                              (GK1)
```

Direct differentiation gives the exact conjugation

```text
L_beta G=exp[-i chi(d)](H_rr+rH).                     (GK2)
```

Thus if `L_beta G=S`, the transformed equation is simply

```text
H_rr+rH=exp[i chi(d(r))]S(d(r)).                     (GK3)
```

Its causal Green kernel is the Airy determinant

```text
K(r,s)=pi[Ai(-r)Bi(-s)-Bi(-r)Ai(-s)],                (GK4)
K(s,s)=0,       partial_r K(s,s)=1.                  (GK5)
```

Consequently, for any reference point `r_0`,

```text
H(r)=-partial_s K(r,r_0)H(r_0)+K(r,r_0)H_r(r_0)
 +Integral_(r_0)^r K(r,s)e^[i chi(d(s))]S(d(s))ds.   (GK6)
```

For `d_m=beta(4m/C-1)`, the Green coordinate is affine:

```text
r_m=beta^2(8m/C-1)+lambda,
Delta_r=8beta^2/C,
r_(q+j)=r_q+j Delta_r,
r_(q-j)=r_q-j Delta_r.                               (GK7)
```

The event-ordered pairs therefore form an exactly symmetric Airy lattice,
not merely an approximately regular sample.  On the full top corridor,

```text
{c['r_min_ball']} <= r <= {c['r_max_ball']},
Delta_r={c['exact_lattice_step_ball']}.              (GK8)
```

DLMF 9.8.20 and its signed first-neglected-term remainder imply, at order
zero and for positive `r`,

```text
M(-r)^2<=1/(pi sqrt(r)).                              (GK9)
```

Using `M(-r)^2=Ai(-r)^2+Bi(-r)^2` in the determinant gives

```text
|K(r,s)|<=pi M(-r)M(-s)<=(rs)^(-1/4)
          <{c['green_kernel_absolute_envelope_ball']}<4.67e-4. (GK10)
```

The pointwise bound alone is not summed over the `r` interval: doing so by
absolute values would throw away the equally spaced event-pair phases.  The
next estimate must insert the compressed sources (11.387.7)--(11.387.8) and
the grouped forcing projection into (GK6), then use the Airy phase and finite
Abel differences on the symmetric lattice.

Pi provenance: the Airy Wronskian contributes the `pi` in (GK4); the lattice
scale comes from `beta^3=pi C^2/8` and the exact Kummer/Fourier detuning.
The modulus theorem is sourced from `{DLMF_URL}`.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No grouped Green-operator estimate, endpoint-completion cancellation,
finite-integral 398-mode splice, all-corridor continuation, complete
`Q_K-T` or `T_upper`, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion
is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    require(dependencies["detuning_ode"].get("passed") is True, "detuning ODE dependency failed")
    require(dependencies["branch_envelope"].get("passed") is True, "branch envelope dependency failed")

    ctx.dps = PRECISION
    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "selected_detuning_equation_conjugated_to_exact_Airy_Green_lattice",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "external_theorem": {
            "source": "NIST Digital Library of Mathematical Functions, Section 9.8",
            "url": DLMF_URL,
            "used_statements": [
                "Airy modulus definition 9.8.3",
                "Airy modulus expansion 9.8.20",
                "signed first-neglected-term remainder statement following 9.8.23",
            ],
        },
        "decision": {
            "detuning_ODE_exactly_conjugated_to_Airy_equation": True,
            "variation_of_constants_kernel_explicit": True,
            "remaining_event_pairs_form_exact_symmetric_Airy_lattice": True,
            "uniform_pointwise_Green_kernel_bound_below_4_67e_minus_4": True,
            "grouped_Green_operator_bound_proved": False,
            "completed_endpoint_cancellation_proved": False,
            "grouped_398_mode_splice_proved": False,
            "complete_T_upper_proved": False,
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
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "Insert the exact endpoint-source aggregates and grouped forcing from the detuning-ODE gate into (GK6). "
            "Derive finite Abel differences of the Airy Green phase on the symmetric r lattice, while composing the "
            "endpoint terms with g4_lambda(0), the half-current, and the Poisson complement before taking norms."
        ),
        "proof_boundary": (
            "Exact Airy gauge conjugation, Green kernel, symmetric event-pair lattice, and pointwise kernel envelope only. "
            "No grouped Green-operator estimate, completed endpoint cancellation, finite-integral splice, complete Q_K-T or "
            "T_upper, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("conjugated selected detuning ODE to an exact Airy Green lattice", flush=True)


if __name__ == "__main__":
    main()
