#!/usr/bin/env python3
"""Expose the exact bulk/lower-tail/upper-tail decomposition of each mode."""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_finite_endpoint_tail_decomposition_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "portcullis": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_finite_poisson_portcullis_saddle_reduction_gate.json",
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
    "Gamma_bulk": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate.json",
    "frozen_barrier": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_frozen_endpoint_current_aggregate_barrier_gate.json",
}


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


def symbolic_certificate() -> dict[str, str]:
    f_a, f_b, e_a, e_b, a, pi = sp.symbols("F_A F_B E_A E_B a pi", nonzero=True)
    imaginary_unit = sp.I
    full_fresnel = 1 + imaginary_unit
    lower_tail = f_a + full_fresnel / 2
    upper_tail = full_fresnel / 2 - f_b
    finite = (e_b - e_a) / (imaginary_unit * pi) + a * (f_b - f_a)
    bulk = a * full_fresnel
    lower = -e_a / (imaginary_unit * pi) - a * lower_tail
    upper = e_b / (imaginary_unit * pi) - a * upper_tail
    require(sp.simplify(finite - bulk - lower - upper) == 0, "finite-current decomposition failed")

    x, endpoint, mode = sp.symbols("x D m", positive=True)
    q = sp.sqrt(x / 2) * (endpoint - 2 * mode / x)
    completed = sp.simplify(-mode**2 / x + q**2 / 2)
    source_endpoint_phase = x * endpoint**2 / 4 - mode * endpoint
    require(sp.simplify(completed - source_endpoint_phase) == 0, "endpoint phase completion failed")

    t = sp.symbols("t", positive=True)
    phase = sp.pi * endpoint**2 * x / 4 + t * sp.log((1 - x) / x) / 2
    derivative = sp.factor(sp.diff(phase, x))
    expected = sp.pi * endpoint**2 / 4 - t / (2 * x * (1 - x))
    require(sp.simplify(derivative - expected) == 0, "endpoint phase derivative failed")

    return {
        "finite_current": "P_m=(E_B-E_A)/(i*pi)+a_m(F_B-F_A), a_m=m*sqrt(2/x)",
        "full_fresnel_value": "F(+infinity)-F(-infinity)=1+i and F(+/-infinity)=+/-(1+i)/2",
        "bulk": "P_bulk=a_m(1+i)",
        "lower_tail": "P_A=-E_A/(i*pi)-a_m[F_A+(1+i)/2]",
        "upper_tail": "P_B=+E_B/(i*pi)-a_m[(1+i)/2-F_B]",
        "exact_decomposition": "P_m=P_bulk+P_A+P_B",
        "completed_endpoint_phase": "-m^2/x+q_D^2/2=x*D^2/4-m*D",
        "odd_endpoint_character": "for odd integer D and integer m, (-1)^m exp(-i*pi*m*D)=1",
        "endpoint_phase_after_parity": "Phi_D(x)=pi*D^2*x/4+(t/2)log((1-x)/x)",
        "endpoint_phase_derivative": str(derivative),
        "grouping_guard": "The bare E_D endpoint current becomes mode-independent after odd parity and is not an independently summable Poisson object. It must remain paired with its Fresnel tail and the symmetric zero/negative/outer-positive completion.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    return """# Exact finite-endpoint tail decomposition

Date: 2026-08-13

Status: exact reassembly reduction proved; not a proof of the quantitative endpoint-tail bounds

Let

```text
a_m=m sqrt(2/x),
E_D=exp(i*pi*q_D^2/2),
F_D=integral_0^(q_D) exp(i*pi*q^2/2)dq.
```

Since `F(+infinity)=(1+i)/2` and
`F(-infinity)=-(1+i)/2`, every finite current splits exactly as

```text
P_m=(E_B-E_A)/(i*pi)+a_m(F_B-F_A)
   =P_bulk+P_A+P_B,                                    (ET1)

P_bulk=a_m(1+i),
P_A=-E_A/(i*pi)-a_m[F_A+(1+i)/2],
P_B=+E_B/(i*pi)-a_m[(1+i)/2-F_B].                      (ET2)
```

Section 11.349 closes `P_bulk` by an exact Gamma integral.  Equations
(ET1)--(ET2) isolate the only remaining finite-endpoint objects without
freezing either at the outer saddle.

There is a mandatory summation guard.  Quadratic completion gives

```text
-m^2/x+q_D^2/2=xD^2/4-mD.                             (ET3)
```

Both source endpoints are odd, so for integer `m`

```text
(-1)^m exp(-i*pi*mD)=1.                               (ET4)
```

The bare endpoint exponential is therefore independent of `m`.  It is not
an independently summable infinite Poisson series.  The lower and upper
terms in (ET2) must remain paired with their Fresnel tails and then be
reassembled with the symmetric zero, negative, and outer-positive modes plus
the endpoint half-current before absolute values are taken.  This explains
why the frozen endpoint-current model failed while the exact bulk succeeded.

After parity, the endpoint phase is

```text
Phi_D(x)=pi D^2 x/4+(t/2)log((1-x)/x),
Phi_D'(x)=pi D^2/4-t/[2x(1-x)].                        (ET5)
```

The lower `A` phase is the characteristic fold already covered locally by
the Airy/logistic atlas.  The upper `B` phase is the nonstationary endpoint
piece already bounded on the transition chart.  The next quantitative task
is to perform their complete symmetric roster reassembly, not to assign
termwise endpoint errors.

Pi provenance: all occurrences derive from the equation-(9) quadratic
Kummer phase and integer Fourier-Poisson character.  No geometric fit is
used.

Proof boundary: exact per-mode decomposition, endpoint phase, and summation
guard only.  No complete symmetric endpoint-tail bound, ordinary/fold
reassembly, complete `T_upper` theorem, height-uniform theorem, or
`Lambda<=0` theorem is proved.  No claim of PF-infinity, RH, or a prize-level
conclusion is made.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    certificate = symbolic_certificate()
    artifact = {
        "kind": STEM,
        "status": "exact_bulk_lower_tail_upper_tail_decomposition_and_symmetric_grouping_guard_proved",
        "passed": True,
        "symbolic_certificate": certificate,
        "decision": {
            "finite_current_decomposed_exactly": True,
            "bulk_closed_by_Gamma_dependency": True,
            "bare_endpoint_current_may_be_summed_independently": False,
            "symmetric_endpoint_tail_reassembly_required": True,
            "complete_symmetric_endpoint_tail_bound_proved": False,
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
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "process_priority": priority,
        },
        "next_obligation": "Construct the finite symmetric mode partial sums of P_A and P_B together with the endpoint half-current, prove a uniform Dirichlet-kernel/Abel limit, and splice the lower characteristic part to the certified event atlas while extending the upper nonstationary bound to the complete roster.",
        "proof_boundary": "Exact per-mode current decomposition and symmetric summation guard only. No complete endpoint-tail bound, complete T_upper theorem, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified exact finite-current bulk/lower-tail/upper-tail decomposition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
