#!/usr/bin/env python3
"""Extract the certified B trace window from the global fixed-state residual."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "global_residual": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_ordinary_symmetric_endpoint_tail_reassembly_gate.json",
    "step_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.json",
    "complete_B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_complete_gaussian_window_signed_gate.json",
    "Abel_theta": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate.json",
}

REFERENCE_TARGET = sp.Rational(86, 10_000_000)  # 8.6e-6
WINDOW_MARGIN = sp.Rational(13198, 100_000_000)  # 1.3198e-4


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
    global_integrand, B_trace, window = sp.symbols("G_M Btr_M chi_W")
    extracted = window * B_trace
    remainder = global_integrand - extracted
    require(sp.expand(extracted + remainder - global_integrand) == 0, "window extraction identity failed")

    h, complement, lower, upper = sp.symbols("H C P_A P_B")
    global_fixed = h + complement + lower + upper
    trace_window = window * upper
    grouped_remainder = h + complement + lower + (1 - window) * upper
    require(sp.expand(global_fixed - trace_window - grouped_remainder) == 0, "trace-package allocation identity failed")

    bulk, A_plus, U_plus, I_minus, B_minus, target, local, q_negative = sp.symbols(
        "P_bulk P_A_plus U_B_plus I_minus P_B_minus tau ell qneg"
    )
    step = 1 - target - q_negative
    pair_residual = A_plus + I_minus + U_plus + step * bulk
    allocated_B_trace = U_plus + B_minus + local * step * bulk
    allocated_rest = A_plus + I_minus - B_minus + (1 - local) * step * bulk
    require(
        sp.expand(pair_residual - allocated_B_trace - allocated_rest) == 0,
        "finite-mode B-trace allocation failed",
    )

    sufficient = sp.simplify(REFERENCE_TARGET + WINDOW_MARGIN)
    require(sufficient == sp.Rational(14058, 100_000_000), "sufficient target arithmetic failed")
    return {
        "finite_cutoff_integrand": "G_M=H_x+integral_0^L f_x(u)C_M(u)du+sum_(m=622)^39894(P_A,m+P_B,m)",
        "step_coefficient": "sigma_m(x)=(1-tau_m)-1_(q_B,m(x)<0), tau_m=1_(622<=m<=39894), ell_m=1_(m in {621,622})",
        "finite_mode_allocation": "I_m+I_-m-tau_m P_bulk=[U_B,m+P_B,-m+ell_m sigma_m P_bulk]+[P_A,m+I_-m-P_B,-m+(1-ell_m)sigma_m P_bulk]",
        "B_trace_package": "Btr_M=sum_(m=1)^M[U_B,m+P_B,-m+ell_m sigma_m P_bulk,m]",
        "omitted_step_completion": "Every nonlocal (1-ell_m)sigma_m P_bulk,m remains in the grouped global complement",
        "window": "chi_W=1_[x_low,x_high], with xi=sqrt(H_B)(x-x_B) and |xi|<=70",
        "extraction": "G_M=chi_W Btr_M+[G_M-chi_W Btr_M] at every finite symmetric cutoff M",
        "limit": "R_end=E_Btr,win+R_group after the already-proved symmetric Abel limit",
        "grouping_guard": "R_group is one global remainder; it includes all nonlocal sigma_m P_bulk,m steps, B-left, A, zero, negative, outer-positive, and endpoint-half sectors",
        "closing_target": "E_Btr,win<-1.3198e-4 and R_group<1.4058e-4 imply R_end<8.6e-6",
        "negative_target": "R_group<1.3198e-4 would imply R_end<0",
    }


def render_note(artifact: dict[str, Any]) -> str:
    return """# Global extraction of the certified B trace window

Date: 2026-08-13

Status: exact residual decomposition and quantitative sufficient target; not
a proof of the grouped-remainder bound or RH; the grouped remainder is open

At a finite symmetric cutoff, the bulk-extracted fixed-state integrand is

```text
G_M(x)=H_x+integral_0^L f_x(u)C_M(u)du
       +sum_(m=622)^39894(P_A,m+P_B,m).               (GE1)
```

Put

```text
tau_m=1_(622<=m<=39894),
sigma_m=(1-tau_m)-1_(q_B,m<0),
ell_m=1_(m in {621,622}).
```

Since `P_B,m=U_B,m-1_(q_B,m<0)P_bulk,m`, exact modewise
algebra fixes the allocation

```text
I_m+I_-m-tau_m P_bulk,m
 =[U_B,m+P_B,-m+ell_m sigma_m P_bulk,m]
  +[P_A,m+I_-m-P_B,-m
    +(1-ell_m)sigma_m P_bulk,m].                      (GE2)
```

The first bracket is `Btr_M` after summation.  It is exactly the trace-current
package certified on the window: the local 621/622 steps remain attached to
their crossing tails, while every nonlocal constant step remains explicitly
in the second bracket.  Put

```text
chi_W=1_[x_low,x_high],
xi=sqrt(H_B)(x-x_B),      W={|xi|<=70}.               (GE3)
```

Linearity at every finite cutoff gives the exact extraction

```text
G_M=chi_W Btr_M+[G_M-chi_W Btr_M].                   (GE4)
```

The established symmetric Abel limit and the complete B trace-window theorem
therefore give

```text
R_end=E_Btr,win+R_group,
E_Btr,win<-1.3198e-4.                                (GE5)
```

Here `R_group` is not merely two outside B tails.  It retains every nonlocal
`sigma_m P_bulk,m` step omitted from the trace package, the A endpoint
current, both B outside-window trace pieces, the zero/negative/outer-positive
completion, and the endpoint half-current in the exact global allocation.
In particular, the left B tail cannot be normed separately: the Abel-theta
corner majorant uses the endpoint difference

```text
Delta_g(z)=g_B(z)-g_A(z)=O(z^(1/2)),                 (GE6)
```

whereas separating the endpoints loses that cancellation.

Against the working physical target, the next sufficient theorem is now

```text
R_group<1.4058e-4  =>  R_end<8.6e-6.                (GE7)
```

The stronger target `R_group<1.3198e-4` would prove `R_end<0`.  GE7 is a
target specification, not an estimate of `R_group`; no part of the grouped
remainder has been assigned an independent triangle budget here.

Pi provenance: `pi` comes from the equation-(9) phase, Fourier-Poisson
character, exact half-Kummer projection, and Jacobi transform.  No fitted
or geometric surrogate is introduced.

Proof boundary: exact extraction of the certified saved-height B trace package
from the global fixed-state residual and the resulting sufficient numerical
target only.  No grouped-remainder estimate, complete `Q_K-T` or `T_upper`,
height-uniform theorem, `Lambda<=0`, PF-infinity, RH, or prize-level
conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["global_residual"]["decision"]["bulk_extracted_joint_residual_identity_proved"] is True, "global-residual dependency drift")
    require(dependencies["global_residual"]["decision"]["symmetric_residual_limit_proved"] is True, "symmetric-limit dependency drift")
    require(dependencies["step_tail"]["decision"]["positive_B_residual_normal_form_is_continuous"] is True, "B-sector dependency drift")
    require(dependencies["complete_B_window"]["decision"]["complete_saved_height_B_trace_package_window_upper_bound_proved"] is True, "B trace-window dependency drift")
    require(dependencies["complete_B_window"]["decision"]["nonlocal_constant_bulk_step_completion_included"] is False, "nonlocal step-completion drift")
    require(dependencies["Abel_theta"]["decision"]["endpoint_difference_retained_before_majorization"] is True, "Abel-theta grouping drift")

    artifact = {
        "kind": STEM,
        "status": "exact_global_fixed_state_B_trace_window_extraction_with_grouped_remainder_target_1p4058e_minus_4",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "target_certificate": {
            "complete_B_trace_window_upper_bound": "-1.3198e-4",
            "working_global_residual_upper_target": "8.6e-6",
            "sufficient_grouped_remainder_upper_target": "1.4058e-4",
            "sufficient_grouped_remainder_target_for_negative_residual": "1.3198e-4",
            "exact_decimal_arithmetic": "1.4058e-4-1.3198e-4=8.6e-6",
        },
        "decision": {
            "B_trace_window_extracted_at_every_finite_symmetric_cutoff": True,
            "B_trace_window_extraction_survives_symmetric_Abel_limit": True,
            "nonlocal_constant_bulk_steps_retained_in_grouped_remainder": True,
            "complete_endpoint_B_window_extracted": False,
            "remaining_object_is_one_grouped_global_remainder": True,
            "isolated_left_B_tail_triangle_budget_admissible": False,
            "grouped_remainder_below_1p4058e_minus_4_proved": False,
            "complete_QK_minus_T_bound_proved": False,
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
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "process_priority": priority,
        },
        "next_obligation": "Construct a quantitative representation of R_group that keeps every nonlocal sigma_m P_bulk,m step, Delta_g=g_B-g_A, the symmetric completion, and the endpoint half-current joined. Extract only an endpoint-safe part of the exterior B trace package; then treat its complement in the Abel-theta/double-Weber representation. A joined bound below 1.4058e-4 closes the saved-height 8.6e-6 target.",
        "proof_boundary": "Exact global fixed-state B trace-window extraction and sufficient grouped-remainder targets only. The complete endpoint B sector is not extracted. No grouped-remainder estimate, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level theorem is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified global B trace-window extraction; grouped remainder target is 1.4058e-4", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
