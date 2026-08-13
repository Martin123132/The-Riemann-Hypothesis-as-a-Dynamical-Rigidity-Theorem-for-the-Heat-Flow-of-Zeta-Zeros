#!/usr/bin/env python3
"""Certify an endpoint-safe tangential extraction of the exterior B trace package."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_tangential_exterior_cutoff_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "global_extraction": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_B_window_extraction_target_gate.json",
    "B_window": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_B_face_gaussian_window_fresnel_remainder_gate.json",
    "Abel_theta": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate.json",
    "double_Weber": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_dual_double_weber_exact_source_reconstruction_gate.json",
}

PRECISION = 100
T = 10_000_000_000
B = 5_122_421
DELTA = sp.Rational(1, 10_000)


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
    u = sp.symbols("u", real=True)
    smoothstep = 6 * u**5 - 15 * u**4 + 10 * u**3
    first = sp.diff(smoothstep, u)
    second = sp.diff(first, u)
    require(smoothstep.subs(u, 0) == 0 and smoothstep.subs(u, 1) == 1, "cutoff values failed")
    require(first.subs(u, 0) == first.subs(u, 1) == 0, "cutoff first-derivative joins failed")
    require(second.subs(u, 0) == second.subs(u, 1) == 0, "cutoff second-derivative joins failed")
    require(sp.simplify(first - 30 * u**2 * (1 - u) ** 2) == 0, "cutoff first derivative failed")
    critical = (sp.Integer(3) - sp.sqrt(3)) / 6
    require(sp.simplify(second.subs(u, critical) - 10 * sp.sqrt(3) / 3) == 0, "cutoff second-derivative maximum failed")

    global_integrand, B_trace, window, cutoff = sp.symbols("G_M Btr_M chi_W eta")
    exterior = 1 - window
    grouped_remainder = global_integrand - window * B_trace
    tangent_B = cutoff * exterior * B_trace
    joined_remainder = grouped_remainder - tangent_B
    require(
        sp.expand(window * B_trace + tangent_B + joined_remainder - global_integrand) == 0,
        "three-part global allocation failed",
    )

    return {
        "cutoff": "eta_delta(x)=0 for x<=delta, S((x-delta)/delta) for delta<x<2delta, and 1 for x>=2delta",
        "smoothstep": "S(u)=6u^5-15u^4+10u^3",
        "cutoff_regular": "S(0)=S'(0)=S''(0)=0 and S(1)=1, S'(1)=S''(1)=0, so eta_delta is C2",
        "derivative_bounds": "|eta_x|<=15/(8delta), |eta_xx|<=10sqrt(3)/(3delta^2)",
        "normalized_derivatives": "xi=sqrt(H_B)(x-x_B), so |eta_xi|<=15/(8delta sqrt(H_B)) and |eta_xixi|<=10sqrt(3)/(3delta^2 H_B)",
        "exterior": "E_W=1-chi_W",
        "tangential_B_trace": "R_Btr,tan(M,epsilon)=eta_delta E_W Btr_(M,epsilon) at a common finite cutoff and Abel regulator",
        "joined_remainder": "R_join(M,epsilon)=G_(M,epsilon)-chi_W Btr_(M,epsilon)-R_Btr,tan(M,epsilon)",
        "global_identity": "R_end=E_Btr,win+lim_(epsilon->0)lim_(M->infinity)[R_Btr,tan(M,epsilon)+R_join(M,epsilon)] in the inherited common prescription",
        "corner_guard": "eta_delta=0 on x<=delta, so the complete extracted B trace package remains joined to the omitted nonlocal steps, A, zero, negative, completion, and half-current sectors there",
        "target": "limsup_(epsilon->0)lim_(M->infinity)[R_Btr,tan(M,epsilon)+R_join(M,epsilon)]<1.4058e-4 is sufficient for R_end<8.6e-6",
    }


def interval_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    t = arb(T)
    endpoint = arb(B)
    pi = arb.pi()
    delta = arb(1) / 10_000
    x_trace = (1 - (1 - 8 * t / (pi * endpoint**2)).sqrt()) / 2
    hessian = t * (1 - 2 * x_trace) / (2 * x_trace**2 * (1 - x_trace) ** 2)
    root_h = hessian.sqrt()
    x_low = x_trace - arb(70) / root_h
    x_high = x_trace + arb(70) / root_h

    def G(x: arb) -> arb:
        return (pi * endpoint**2 / 4 - t / (2 * x * (1 - x))) / root_h

    left_gap = -G(x_low)
    right_gap = G(x_high)
    require(2 * delta < x_low, "cutoff transition reaches the B window")
    require(left_gap > arb("70.06"), "left tangential gap below 70.06")
    require(right_gap > arb("69.92"), "right tangential gap below 69.92")

    d1_xi = (arb(15) / (8 * delta)) / root_h
    d2_xi = (arb(10) * arb(3).sqrt() / (3 * delta**2)) / hessian
    require(d1_xi < arb("6.5e-5"), "normalized first cutoff derivative too large")
    require(d2_xi < arb("7e-9"), "normalized second cutoff derivative too large")
    transition_high_xi = (2 * delta - x_trace) * root_h
    require(transition_high_xi < arb(-12_000), "cutoff transition is not far left of the B window")

    return {
        "height": T,
        "endpoint": B,
        "delta": "1e-4",
        "cutoff_transition": ["1e-4", "2e-4"],
        "trace_root_ball": x_trace.str(PRECISION, more=True),
        "trace_hessian_ball": hessian.str(PRECISION, more=True),
        "sqrt_trace_hessian_ball": root_h.str(PRECISION, more=True),
        "window_low_ball": x_low.str(PRECISION, more=True),
        "window_high_ball": x_high.str(PRECISION, more=True),
        "transition_high_xi_ball": transition_high_xi.str(PRECISION, more=True),
        "left_exterior_abs_G_lower_ball": left_gap.str(PRECISION, more=True),
        "right_exterior_G_lower_ball": right_gap.str(PRECISION, more=True),
        "normalized_cutoff_first_derivative_upper_ball": d1_xi.str(PRECISION, more=True),
        "normalized_cutoff_second_derivative_upper_ball": d2_xi.str(PRECISION, more=True),
        "uniform_exterior_abs_G_lower_bound": "69.92",
    }


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["interval_certificate"]
    return f"""# Endpoint-safe tangential extraction of the exterior B trace package

Date: 2026-08-13

Status: exact global decomposition and saved-height phase/cutoff certificate;
not a proof of either remaining quantitative bound or RH

Section 11.414 gives `R_end=E_Btr,win+R_group`.  A direct norm of the left
B trace tail is inadmissible because the Abel corner uses the joined endpoint
difference.  Put `delta=10^-4` and define

```text
S(u)=6u^5-15u^4+10u^3,
eta_delta(x)=0                         (x<=delta),
             S((x-delta)/delta)       (delta<x<2delta),
             1                         (x>=2delta).       (TE1)
```

The endpoint values and first two derivatives match, so `eta_delta` is C2.
Moreover

```text
|eta_x| <= 15/(8delta),
|eta_xx|<=10sqrt(3)/(3delta^2).                         (TE2)
```

Let `E_W=1-chi_W`.  At every finite symmetric cutoff define

```text
R_Btr,tan,M = eta_delta E_W Btr_M,
R_join,M    = G_M-chi_W Btr_M-eta_delta E_W Btr_M.     (TE3)
```

Then exactly

```text
G_M=chi_W Btr_M+R_Btr,tan,M+R_join,M.                  (TE4)
```

Because `eta_delta=0` for `x<=delta`, `R_join,M` retains the complete
extracted B trace package there with every omitted nonlocal bulk step, the A
endpoint, zero and negative modes, completion, and endpoint half-current at
the singular corner.  Thus the
previous `Delta_g=g_B-g_A=O(z^(1/2))` majorant and the exact double-Weber
zero/half-current collapse are not broken.

The transition ends at `x=2e-4`, while

```text
x_low={c['window_low_ball']},
xi(2delta)={c['transition_high_xi_ball']}<-12000.       (TE5)
```

Outside the removed window the B trace phase is already tangentially
nonstationary.  Monotonicity gives

```text
-G_B(x_low)={c['left_exterior_abs_G_lower_ball']}>70.06,
 G_B(x_high)={c['right_exterior_G_lower_ball']}>69.92,
 |G_B|>69.92 on supp(eta_delta E_W).                   (TE6)
```

In normalized trace coordinate `xi=sqrt(H_B)(x-x_B)`, the cutoff costs only

```text
|eta_xi|  <= {c['normalized_cutoff_first_derivative_upper_ball']} <6.5e-5,
|eta_xixi|<= {c['normalized_cutoff_second_derivative_upper_ball']} <7e-9. (TE7)
```

Keeping one common Abel regulator therefore leaves the exact target

```text
R_end=E_Btr,win+lim_(eps->0)lim_(M->infinity)
                  [R_Btr,tan(M,eps)+R_join(M,eps)],
limsup_(eps->0)lim_(M->infinity)[R_Btr,tan+R_join]
                  <1.4058e-4  ==> R_end<8.6e-6.        (TE8)
```

No separate `eps->0` limit for the two new summands is asserted here.

The next quantitative step is to bound the first two `xi` derivatives of
the complete summed exterior B trace amplitude before absolute values and apply
tangential integration by parts using (TE6).  The residual `R_join` must
remain in its Abel-theta/double-Weber endpoint grouping.

Pi provenance: `pi` enters only through the exact equation-(9) B trace
phase and its Hessian.  The cutoff is a rational polynomial and introduces
no fitted geometric constant.

Proof boundary: exact endpoint-safe exterior decomposition, C2 cutoff, and
saved-height tangential phase margins only.  Neither a bound for
`R_Btr,tan` or `R_join`, a complete `Q_K-T` or `T_upper` result, nor a
height-uniform theorem is obtained.  These results do not prove
`Lambda<=0`, PF-infinity, RH, or a prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["global_extraction"]["decision"]["remaining_object_is_one_grouped_global_remainder"] is True, "global grouping drift")
    require(dependencies["global_extraction"]["decision"]["nonlocal_constant_bulk_steps_retained_in_grouped_remainder"] is True, "global nonlocal-step grouping drift")
    require(dependencies["B_window"]["decision"]["xi_70_window_contains_exactly_B_crossings_621_622"] is True, "B-window geometry drift")
    require(dependencies["Abel_theta"]["decision"]["endpoint_difference_retained_before_majorization"] is True, "Abel endpoint grouping drift")
    require(dependencies["double_Weber"]["decision"]["entire_finite_equation9_source_roster_reconstructed_exactly"] is True, "double-Weber reconstruction drift")

    artifact = {
        "kind": STEM,
        "status": "endpoint_safe_exterior_B_trace_tangential_extraction_with_uniform_abs_G_above_69p92",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(),
        "decision": {
            "C2_corner_cutoff_proved": True,
            "cutoff_transition_disjoint_from_B_window": True,
            "exterior_B_trace_uniformly_tangentially_nonstationary": True,
            "global_grouped_remainder_decomposition_proved_at_common_cutoff_and_regulator": True,
            "Abel_endpoint_difference_retained_in_joined_remainder": True,
            "double_Weber_zero_and_half_current_grouping_retained": True,
            "separate_exterior_and_joined_Abel_limits_proved": False,
            "exterior_B_trace_tangential_bound_proved": False,
            "joined_remainder_bound_proved": False,
            "grouped_remainder_below_1p4058e_minus_4_proved": False,
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
        "next_obligation": "Differentiate the complete finite-cutoff, Abel-regulated exterior B trace package in normalized xi before taking absolute values, prove two bounds uniform in cutoff and regulator through the sign-adapted crossing cancellations, and apply two tangential integrations by parts on eta_delta(1-chi_W)Btr_(M,epsilon). Keep every omitted nonlocal bulk step in R_join under the same Abel-theta/double-Weber prescription. The final target is the joined regulated limsup below 1.4058e-4; separate limits must not be assumed.",
        "proof_boundary": "Exact endpoint-safe exterior decomposition, C2 cutoff, and saved-height tangential phase margins only. No exterior-B or joined-remainder bound, complete Q_K-T or T_upper, height-uniform theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified endpoint-safe exterior B trace extraction with tangential |G|>69.92", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
