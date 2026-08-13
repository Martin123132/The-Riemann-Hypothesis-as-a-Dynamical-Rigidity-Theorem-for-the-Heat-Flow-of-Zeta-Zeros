#!/usr/bin/env python3
"""Certify the exact global-Morse Gamma factor for the ordinary bulk current."""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "work/rh_compute/vendor"))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_global_morse_gamma_bulk_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "coverage_ledger": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_ordinary_mode_coverage_ledger_gate.json",
    "universal_morse": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json",
    "fresnel_transform": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.json",
}

PRECISION = 100
T = 10_000_000_000
LOWER_START = 622
UPPER_END = 39_894
TARGET_NORMALIZED = arb("0.000019")
TARGET_PHYSICAL = arb("0.0000086")


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


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def symbolic_reduction() -> dict[str, str]:
    r, v, m = sp.symbols("r v m", positive=True)
    x = 1 / (1 + r * v)
    weight = r * x * (x * (1 - x)) ** (-sp.Rational(1, 4))
    full_current = m * sp.sqrt(2 / x) * (1 + sp.I)
    reduced = sp.powsimp(weight * full_current, force=True)
    target = m * sp.sqrt(2) * (1 + sp.I) * r ** sp.Rational(3, 4) * v ** (-sp.Rational(1, 4))
    require(sp.simplify(reduced / target - 1) == 0, "bulk amplitude reduction failed")
    return {
        "full_line_current": "P_bulk=m*sqrt(2/x)*(1+i)",
        "global_morse_bulk_amplitude": "A_bulk(s)=m*sqrt(2)*(1+i)*r^(3/4)*v(s)^(-1/4)*dv/ds",
        "Abel_integral": "I_bulk(t)=lim_(eps downarrow 0) exp(i*t/2) integral_0^infinity v^(-1/4+i*t/2) exp(-(eps+i*t/2)v)dv",
        "Gamma_value": "I_bulk(t)=exp(i*t/2)*Gamma(3/4+i*t/2)/(i*t/2)^(3/4+i*t/2)",
        "Gaussian_leading_value": "I_G(t)=2*sqrt(pi/t)*exp(-i*pi/4)",
        "common_mode_factor": "B(t)=I_bulk(t)/I_G(t)",
    }


def gamma_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(T)
    pi = arb.pi()
    imaginary_unit = acb(0, 1)
    y = t / 2
    exponent = acb(arb("0.75"), y)

    log_integral = exponent.lgamma() + imaginary_unit * y - exponent * (y.log() + imaginary_unit * pi / 2)
    log_gaussian = arb(2).log() + (pi.log() - t.log()) / 2 - imaginary_unit * pi / 4
    log_factor = log_integral - log_gaussian
    factor = log_factor.exp()

    theta_exact = acb(arb("0.25"), y).lgamma().imag - y * pi.log()
    theta_zero = y * ((t / (2 * pi)).log() - 1) - pi / 8
    theta_correction = theta_exact - theta_zero
    target_phase_factor = (imaginary_unit * theta_correction).exp()
    coefficient_defect = factor - target_phase_factor

    reciprocal_sqrt_sum = arb(0)
    for mode in range(LOWER_START, UPPER_END + 1):
        reciprocal_sqrt_sum += arb(mode).rsqrt()
    aggregate_complex_bound = abs(coefficient_defect) * reciprocal_sqrt_sum
    aggregate_two_real_bound = 2 * aggregate_complex_bound
    physical_bound = arb(2).sqrt() / pi * aggregate_complex_bound

    require(abs(log_factor.imag - theta_correction) < arb("1e-80"), "bulk phase identity lost at working precision")
    require(abs(coefficient_defect) < arb("3.2e-22"), "bulk coefficient defect exceeds 3.2e-22")
    require(aggregate_two_real_bound < arb("3e-19"), "bulk aggregate exceeds 3e-19")
    require(physical_bound < arb("1e-19"), "physical bulk aggregate exceeds 1e-19")
    require(aggregate_complex_bound < TARGET_NORMALIZED / arb("1e14"), "normalized headroom below 1e14")
    require(physical_bound < TARGET_PHYSICAL / arb("1e13"), "physical headroom below 1e13")

    return {
        "height": T,
        "mode_range": [LOWER_START, UPPER_END],
        "mode_count": UPPER_END - LOWER_START + 1,
        "precision_decimal_digits": PRECISION,
        "log_Gamma_bulk_factor_ball": complex_record(log_factor),
        "Gamma_bulk_factor_ball": complex_record(factor),
        "exact_theta_minus_theta_zero_ball": theta_correction.str(PRECISION, more=True),
        "target_phase_factor_ball": complex_record(target_phase_factor),
        "bulk_coefficient_defect_ball": complex_record(coefficient_defect),
        "bulk_coefficient_defect_absolute_ball": abs(coefficient_defect).str(PRECISION, more=True),
        "reciprocal_sqrt_mode_sum_ball": reciprocal_sqrt_sum.str(PRECISION, more=True),
        "aggregate_complex_absolute_bound_ball": aggregate_complex_bound.str(PRECISION, more=True),
        "aggregate_two_real_absolute_bound_ball": aggregate_two_real_bound.str(PRECISION, more=True),
        "aggregate_physical_absolute_bound_ball": physical_bound.str(PRECISION, more=True),
        "normalized_target_headroom_factor_ball": (TARGET_NORMALIZED / aggregate_complex_bound).str(PRECISION, more=True),
        "physical_target_headroom_factor_ball": (TARGET_PHYSICAL / physical_bound).str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# Exact global-Morse Gamma bulk

Date: 2026-08-13

Status: ordinary full-line bulk closed at the saved height; not a proof of the endpoint charts

After the exact logistic substitution, replace only the finite alpha current
by its full-line bulk part

```text
P_bulk=m sqrt(2/x)(1+i).
```

The complete Kummer weight and Jacobian then simplify exactly:

```text
A_bulk(s)=m sqrt(2)(1+i) r^(3/4)
          v(s)^(-1/4) dv/ds.                           (GB1)
```

No stationary-phase truncation is needed.  With Abel damping first and the
limit taken afterward,

```text
I_bulk(t)=exp(i t/2) Gamma(3/4+i t/2)
          /(i t/2)^(3/4+i t/2).                       (GB2)
```

Divide by the leading Gaussian
`I_G=2 sqrt(pi/t) exp(-i*pi/4)` and write `B=I_bulk/I_G`.
At `t=10^10`, rigorous Arb evaluation gives

```text
log B={c['log_Gamma_bulk_factor_ball']['real_ball']}
      + i {c['log_Gamma_bulk_factor_ball']['imag_ball']},
theta-theta_0={c['exact_theta_minus_theta_zero_ball']},
|B-exp(i(theta-theta_0))|
 ={c['bulk_coefficient_defect_absolute_ball']}.        (GB3)
```

Thus the imaginary part of `log B` supplies the exact saved-height
Riemann--Siegel phase correction at the working precision; only a relative
amplitude defect of about `3.125e-22` remains.  Summing it over every source
target mode `622..39894` before taking the real projection gives

```text
complex aggregate <= {c['aggregate_complex_absolute_bound_ball']},
two-real aggregate <= {c['aggregate_two_real_absolute_bound_ball']},
equation-(9) physical <= {c['aggregate_physical_absolute_bound_ball']}. (GB4)
```

This is more than `{c['normalized_target_headroom_factor_ball']}` times below
the normalized target.  The ordinary bulk is therefore not the quantitative
wall at the saved height.

The result does not discard finite-endpoint terms.  Their lower
characteristic current and upper nonstationary current are separate exact
pieces and must be reassembled with the already certified fold and endpoint
charts.  In particular, (GB4) does not turn the frozen endpoint-current
model into a valid approximation.

Pi provenance: `pi` comes from the equation-(9) Kummer/Fourier phase, the
Gaussian Gamma integral, and the Riemann--Siegel theta normalization.  No
geometric fit is introduced.

Proof boundary: exact Abel-regularized full-line bulk at `t=10^10` and its
finite source-mode aggregate only.  No finite-endpoint reassembly, ordinary
endpoint-tail bound, complete `T_upper` theorem, height-uniform theorem, or
`Lambda<=0` theorem is proved.  No claim of PF-infinity, RH, or a prize-level
conclusion is made.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    coverage = dependencies["coverage_ledger"]["mode_coverage"]["source_classical_target"]
    require(coverage == {"range": [LOWER_START, UPPER_END], "count": 39_273}, "source target drift")

    artifact = {
        "kind": STEM,
        "status": "saved_height_ordinary_full_line_bulk_exactly_reduced_to_Gamma_and_closed",
        "passed": True,
        "symbolic_reduction": symbolic_reduction(),
        "certificate": gamma_certificate(),
        "decision": {
            "ordinary_full_line_bulk_closed_at_saved_height": True,
            "stationary_phase_truncation_needed_for_bulk": False,
            "finite_lower_endpoint_current_closed_here": False,
            "finite_upper_endpoint_current_closed_here": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Use the exact identity finite current = full-line bulk plus lower endpoint/Fresnel tail minus upper endpoint/Fresnel tail. Extend the grouped B nonstationary enclosure to the ordinary roster, then match the lower tail to the certified fold charts without freezing it at the saddle.",
        "proof_boundary": "Exact Abel-regularized full-line bulk and finite source-mode aggregate at t=10^10 only. No finite-endpoint reassembly, complete T_upper theorem, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified exact global-Morse Gamma bulk: "
        f"aggregate two-real<{arb(artifact['certificate']['aggregate_two_real_absolute_bound_ball']).upper()}, "
        f"modes={artifact['certificate']['mode_count']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
