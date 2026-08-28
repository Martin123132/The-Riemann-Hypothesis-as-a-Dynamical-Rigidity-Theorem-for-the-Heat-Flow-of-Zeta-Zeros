#!/usr/bin/env python3
"""Extract the crossing-completed exterior B current away from the Abel corner."""

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

from flint import arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_endpoint_safe_completed_B_extraction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "step_tail": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_positive_B_crossing_step_tail_normal_form_gate.json",
    "scope_guard": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_trace_package_bulk_step_scope_guard_gate.json",
    "global_extraction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate.json",
    "tangential_cutoff": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_tangential_exterior_cutoff_gate.json",
    "Abel_theta": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate.json",
    "double_Weber": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate.json",
}

PRECISION = 90
B = 5_122_421


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
    bulk, p_b_minus, u_b = sp.symbols("P_bulk P_B_minus U_B")
    global_current, trace, window, exterior, cutoff, nonlocal_steps = sp.symbols(
        "G Btr chi_W E_W eta N"
    )
    tau, local, q_negative = sp.symbols("tau ell qneg")
    sigma = 1 - tau - q_negative

    p_b = u_b - q_negative * bulk
    trace_mode = u_b + p_b_minus + local * sigma * bulk
    omitted_mode = (1 - local) * sigma * bulk
    completed_mode = p_b + p_b_minus + (1 - tau) * bulk
    require(
        sp.expand(trace_mode + omitted_mode - completed_mode) == 0,
        "completed B-mode identity failed",
    )

    old_join = global_current - window * trace - cutoff * exterior * trace
    completed_trace = trace + nonlocal_steps
    completed_exterior = cutoff * exterior * completed_trace
    new_join = global_current - window * trace - completed_exterior
    require(
        sp.expand(new_join - (old_join - cutoff * exterior * nonlocal_steps)) == 0,
        "joined-remainder transfer identity failed",
    )
    require(
        sp.expand(window * trace + completed_exterior + new_join - global_current) == 0,
        "completed exterior allocation failed",
    )

    crossing_rows: list[dict[str, int]] = []
    for tau_value, local_value in ((0, 0), (1, 0), (0, 1), (1, 1)):
        sigma_left = 1 - tau_value - 1
        sigma_right = 1 - tau_value
        trace_jump = -1 + local_value * (sigma_right - sigma_left)
        omitted_jump = (1 - local_value) * (sigma_right - sigma_left)
        completed_jump = trace_jump + omitted_jump
        require(completed_jump == 0, "crossing completion failed")
        crossing_rows.append(
            {
                "tau": tau_value,
                "ell": local_value,
                "trace_jump_in_P_bulk_units": trace_jump,
                "omitted_step_jump_in_P_bulk_units": omitted_jump,
                "completed_jump_in_P_bulk_units": completed_jump,
            }
        )

    return {
        "definitions": {
            "tau": "tau_m=1_(622<=m<=39894)",
            "ell": "ell_m=1_(m in {621,622})",
            "sigma": "sigma_m=(1-tau_m)-1_(q_B,m<0)",
            "nonlocal_steps": "N_(M,epsilon)=sum_(m=1)^M w_(m,epsilon)(1-ell_m)sigma_m P_bulk,m",
            "completed_B": "Ctr_(M,epsilon)=Btr_(M,epsilon)+N_(M,epsilon)",
        },
        "mode_identity": "U_B,m+P_B,-m+sigma_m P_bulk,m=P_B,m+P_B,-m+(1-tau_m)P_bulk,m",
        "crossing_rows": crossing_rows,
        "regularity": "P_B,m, P_B,-m, and P_bulk,m are smooth for x>0, so every finite completed mode and Ctr_(M,epsilon) are smooth across q_B,m=0",
        "old_trace_obstruction": "for ell_m=0, Btr has jump -P_bulk,m at q_B,m=0; its classical first and second derivatives do not exist across the full released exterior",
        "completed_allocation": "G=chi_W Btr+eta_delta E_W Ctr+R_join,comp",
        "joined_transfer": "R_join,comp=R_join-eta_delta E_W N",
    }


def geometry_certificate(dependencies: dict[str, dict[str, Any]]) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    scope = dependencies["scope_guard"]["interval_certificate"]
    c_low = arb(scope["c_window_low_ball"])
    c_high = arb(scope["c_window_high_ball"])
    c_transition_high = arb(B) / 10_000
    c_half = arb(B) / 4

    require(c_transition_high < arb(620) < c_low, "mode 620 is not in the fully released left exterior")
    require(c_high < arb(623) < c_half, "mode 623 is not in the fully released right exterior")

    return {
        "c_at_2delta_exact": "512.2421",
        "c_window_low_ball": c_low.str(PRECISION, more=True),
        "c_window_high_ball": c_high.str(PRECISION, more=True),
        "c_at_half_exact": "1280605.25",
        "fully_released_nonlocal_crossing_witnesses": [620, 623],
        "left_positive_eta_crossing_range": "257..620",
        "right_exterior_crossing_range": "623..1280605",
        "left_positive_eta_crossing_count": 364,
        "right_exterior_crossing_count": 1_279_983,
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["geometry_certificate"]
    return f"""# Endpoint-safe extraction of the crossing-completed B current

Date: 2026-08-13

Status: exact finite-cutoff recombination; not a proof artifact

No bound for the completed exterior current, joined remainder, or RH is
asserted.

For a common finite cutoff and Abel weight `w_(m,epsilon)`, retain

```text
tau_m=1_(622<=m<=39894),
ell_m=1_(m in {{621,622}}),
sigma_m=(1-tau_m)-1_(q_B,m<0),

N_(M,epsilon)
 =sum_(m=1)^M w_(m,epsilon)(1-ell_m)sigma_m P_bulk,m. (CE1)
```

The trace package by itself is not differentiable across the full released
exterior.  At every nonlocal crossing `q_B,m=0`, `U_B,m` jumps by
`-P_bulk,m`, while the omitted step in (CE1) jumps by `+P_bulk,m`.  This is
not a remote possibility: modes 620 and 623 lie in the fully-on part of the
released left and right exteriors, respectively.  More broadly, the support
with positive cutoff contains the crossing ranges

```text
left:  {c['left_positive_eta_crossing_range']} ({c['left_positive_eta_crossing_count']} crossings),
right: {c['right_exterior_crossing_range']} ({c['right_exterior_crossing_count']} crossings). (CE2)
```

Consequently a classical first- or second-derivative bound for `Btr` alone
cannot be uniform on those exterior intervals.  Distributional jump terms
would have to be retained and cancelled against `R_join`.

The cancellation can instead be made pointwise before differentiation.  From

```text
P_B,m=U_B,m-1_(q_B,m<0)P_bulk,m
```

one obtains the exact mode identity

```text
U_B,m+P_B,-m+sigma_m P_bulk,m
 =P_B,m+P_B,-m+(1-tau_m)P_bulk,m.                    (CE3)
```

The right side is smooth for `x>0`.  Define the crossing-completed trace

```text
Ctr_(M,epsilon)=Btr_(M,epsilon)+N_(M,epsilon).        (CE4)
```

With the same `chi_W`, `E_W=1-chi_W`, and rational C2 cutoff from the
previous gate, replace the old exterior allocation by

```text
R_B,comp=eta_delta E_W Ctr_(M,epsilon),
R_join,comp=G_(M,epsilon)-chi_W Btr_(M,epsilon)-R_B,comp.

G_(M,epsilon)=chi_W Btr_(M,epsilon)
               +R_B,comp+R_join,comp.                (CE5)
```

Equivalently,

```text
R_join,comp=R_join-eta_delta E_W N_(M,epsilon).       (CE6)
```

Because `eta_delta=0` for `x<=delta`, (CE6) changes nothing at the Abel
corner: the endpoint difference, zero mode, negative modes, outer completion,
and half-current remain in the inherited common prescription.  Because
`E_W=0` on the certified window, the negative estimate
`E_Btr,win<-1.3198e-4` is also unchanged.  Away from both protected regions,
the extracted B current is now crossing-complete and classically smooth at
every finite cutoff.

This corrects the next quantitative obligation.  Tangential integration by
parts must be applied to `Ctr`, not to the discontinuous trace package alone.
One must still derive first and second normalized-amplitude bounds uniform in
`M` and `epsilon`, retain the piecewise window-edge boundary terms, and bound
the common-regulator sum with `R_join,comp` below `1.4058e-4`.

Pi provenance: every `pi` remains inherited from the exact equation-(9)
Fresnel current, its bulk completion, and the B phase.  The transfer (CE1)--
(CE6) is algebraic and introduces no geometric or fitted value of `pi`.

Proof boundary: exact finite-cutoff/common-regulator allocation, a proof that
the uncompleted exterior trace has nonlocal jumps, and smooth modewise
crossing completion only.  This gate does not prove a derivative norm,
exterior-current bound, joined-remainder bound, complete `Q_K-T` or `T_upper`,
or height-uniform theorem.  It makes no claim of `Lambda<=0`, PF-infinity,
RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["step_tail"]["decision"]["B_endpoint_rewritten_as_sign_adapted_tail_plus_bulk_step"] is True, "step-tail identity drift")
    require(dependencies["scope_guard"]["decision"]["nonlocal_bulk_steps_owned_by_grouped_remainder"] is True, "nonlocal-step ownership drift")
    require(dependencies["global_extraction"]["decision"]["nonlocal_constant_bulk_steps_retained_in_grouped_remainder"] is True, "global ownership drift")
    require(dependencies["tangential_cutoff"]["decision"]["C2_corner_cutoff_proved"] is True, "corner cutoff drift")
    require(dependencies["tangential_cutoff"]["decision"]["Abel_endpoint_difference_retained_in_joined_remainder"] is True, "corner grouping drift")
    require(dependencies["Abel_theta"]["decision"]["endpoint_difference_retained_before_majorization"] is True, "Abel-theta grouping drift")
    require(dependencies["double_Weber"]["decision"]["entire_finite_equation9_source_roster_reconstructed_exactly"] is True, "double-Weber grouping drift")

    artifact = {
        "kind": STEM,
        "status": "endpoint_safe_crossing_completed_B_extraction_with_trace_derivative_obstruction_corrected",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "geometry_certificate": geometry_certificate(dependencies),
        "decision": {
            "finite_cutoff_completed_B_mode_identity_proved": True,
            "nonlocal_B_trace_jump_obstruction_proved": True,
            "classical_derivatives_of_uncompleted_exterior_B_trace_available": False,
            "crossing_completed_exterior_B_current_smooth_at_finite_cutoff": True,
            "endpoint_safe_completed_B_allocation_proved": True,
            "Abel_corner_grouping_unchanged": True,
            "negative_B_trace_window_certificate_unchanged": True,
            "completed_exterior_first_second_derivative_bounds_proved": False,
            "joined_remainder_below_1p4058e_minus_4_proved": False,
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
            "flint_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Differentiate the crossing-completed finite-cutoff amplitude Ctr_(M,epsilon)=Btr_(M,epsilon)+N_(M,epsilon), not Btr alone. Derive first and second xi-amplitude bounds uniform in M and epsilon, include both window-edge boundary terms, and combine the resulting tangential estimate with R_join,comp under the inherited common Abel-theta/double-Weber prescription. The joined limsup target remains below 1.4058e-4.",
        "proof_boundary": "Exact finite-cutoff crossing completion, derivative-obstruction correction, and endpoint-safe allocation only. No derivative bounds, completed exterior-current or joined-remainder estimate, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified endpoint-safe crossing-completed B extraction; uncompleted trace derivatives rejected", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
