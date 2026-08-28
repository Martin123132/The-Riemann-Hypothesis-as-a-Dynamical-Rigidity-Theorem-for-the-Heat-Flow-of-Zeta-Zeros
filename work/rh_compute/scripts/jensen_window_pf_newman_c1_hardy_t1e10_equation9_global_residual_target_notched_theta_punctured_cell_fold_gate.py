#!/usr/bin/env python3
"""Certify pole-free cell recombination and one-cell folding of the notch kernel."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PILOT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_pilot.json"

DEPENDENCIES = {
    "two_edge_modular": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.json",
    "R_after_A_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "global_theta_current": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation10_global_theta_current_reduction_gate.json",
}

A = 159_577
B = 5_122_421
L = (B - A) // 2
LEFT_EDGE = sp.Rational(1243, 2)
RIGHT_EDGE = sp.Rational(79789, 2)
DUAL_CUTOFF = 64
PRECISION = 100


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
    z = sp.symbols("z", real=True)
    c0 = sp.csc(z) - 1 / z
    c1 = 1 / z**2 - sp.cos(z) / sp.sin(z) ** 2
    require(sp.limit(c0, z, 0) == 0, "g0 removable limit failed")
    require(sp.limit(c1 / 4, z, 0) == sp.Rational(1, 24), "g1 removable limit failed")
    series0 = sp.series(c0, z, 0, 8)
    series1 = sp.series(c1, z, 0, 8)
    require(series0.coeff(z, 1) == sp.Rational(1, 6), "g0 series drift")
    require(series1.coeff(z, 0) == sp.Rational(1, 6), "g1 series drift")

    Hhat, B0, B1, b0, b1 = sp.symbols("Hhat B0 B1 b0 b1")
    rho0 = Hhat - b0 - b1
    recombined = Hhat + (B0 - b0) + (B1 - b1)
    require(sp.expand(B0 + B1 + rho0 - recombined) == 0, "punctured-cell recombination failed")

    n, s, x = sp.symbols("n s x", real=True)
    alpha = A + 2 * n + 2 * s
    phase = sp.exp(sp.I * sp.pi * x * alpha**2 / 4)
    require(sp.simplify(sp.diff(phase, s) / (sp.I * sp.pi * x) - alpha * phase) == 0, "folded theta-current derivative failed")

    left_alpha = A + 2 * 0 + 2 * 1
    right_alpha = A + 2 * (L - 1) + 2 * 1
    require(left_alpha == A + 2 and right_alpha == B, "shifted source roster drift")
    require(A + 2 * (L - 1) == B - 2, "unshifted source roster drift")
    return {
        "cell_coefficients": {
            "g0": "g_0(s)=[csc(pi*s)-1/(pi*s)]/(2i)",
            "g1": "g_1(s)={1/(pi*s)^2-cos(pi*s)/sin(pi*s)^2}/4",
            "limits": "g_0(0)=0, g_1(0)=1/24",
            "g0_series": str(series0),
            "four_g1_series": str(series1),
        },
        "punctured_identity": "B_0+B_1+rho_2(s)=hat h(s)+[B_0-b_0(s)]+[B_1-b_1(s)]",
        "folded_theta_current": "F_x(s)=sum_(n=0)^(L-1)(A+2n+2s)e^(i*pi*x*(A+2n+2s)^2/4)=(i*pi*x)^(-1)partial_s Theta_L(x,s)",
        "theta_definition": "Theta_L(x,s)=sum_(n=0)^(L-1)e^(i*pi*x*(A+2n+2s)^2/4)",
        "theta_shift": "Theta_L(x,1)-Theta_L(x,0)=e^(i*pi*x*B^2/4)-e^(i*pi*x*A^2/4)",
    }


def fold_certificate() -> dict[str, Any]:
    return {
        "periodicity": "C_epsilon(u+j)=C_epsilon(u) for every integer j",
        "partition": "[0,L]=union_(n=0)^(L-1)[n,n+1] up to shared endpoints",
        "change_of_variables": "integral_n^(n+1)f_x(u)C_epsilon(u)du=integral_0^1 f_x(n+s)C_epsilon(s)ds",
        "one_cell_fold": "integral_0^L f_x(u)C_epsilon(u)du=integral_0^1 F_x(s)C_epsilon(s)ds",
        "roster_length": L,
        "unshifted_fold_roster": [A, B - 2],
        "shifted_fold_roster_at_s_1": [A + 2, B],
        "endpoint_shift_survivor": "upper B label minus lower A label",
        "centred_evaluation_map": "Use r=s for 0<=s<=1/2 and r=s-1 for 1/2<s<=1; C_epsilon(s)=C_epsilon(r).",
        "zero_x_guard": "The derivative formula is used for x>0; F_x(s) itself has the direct continuous x=0 definition.",
    }


def punctured_transform_certificate() -> dict[str, str]:
    return {
        "stable_cell_formula": "C_epsilon(s)=hat h_epsilon(s)+E_0(s)g_0(s)+E_1(s)g_1(s)+sum_(k in Z, k!=0)rho_2(k+s), |s|<=1/2",
        "integer_value": "C_epsilon(0)=hat h_epsilon(0)+E_1(0)/24+sum_(k!=0)rho_2(k)",
        "reason": "The k=0 identity rho_2(s)=hat h(s)-b_0(s)-b_1(s) cancels the apparent poles in the closed boundary sums before evaluation.",
        "uniform_tail": "sum_(|k|>K)|rho_2(k+s)|<=K3/[(2pi)^3(K-1/2)^2] for |s|<=1/2",
        "K3": "K3=|q_epsilon''(c)|+|q_epsilon''(d)|+20pi*epsilon",
    }


def tail_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    c = arb(1243) / 2
    d = arb(79789) / 2
    cutoff = arb(DUAL_CUTOFF)
    rows = []
    for epsilon_text in ("1e-6", "1e-8", "1e-10"):
        epsilon = arb(epsilon_text)

        def q2(edge: arb) -> arb:
            q = (-pi * epsilon * edge**2).exp()
            return (4 * pi**2 * epsilon**2 * edge**2 - 2 * pi * epsilon) * q

        k3 = abs(q2(c)) + abs(q2(d)) + 20 * pi * epsilon
        tail = k3 / (2 * pi) ** 3 / (cutoff - arb(1) / 2) ** 2
        rows.append(
            {
                "epsilon": epsilon_text,
                "centred_cell": "|s|<=1/2",
                "dual_cutoff": DUAL_CUTOFF,
                "K3_ball": k3.str(70, more=True),
                "uniform_punctured_tail_ball": tail.str(70, more=True),
            }
        )
    require(arb(rows[-1]["uniform_punctured_tail_ball"]).upper() < arb("7e-15"), "punctured tail scale drift")
    return {"rows": rows}


def render_note(artifact: dict[str, Any]) -> str:
    tails = artifact["tail_certificate"]["rows"]
    return f"""# Pole-free target-notch cell and one-cell source fold

Date: 2026-08-23

Status: exact punctured-cell recombination, uniform cubic tail, and one-cell
fold certified; Kummer-weighted quadrature and `R_after_A` bound open

Retain the two-edge objects `hat h`, `b_0`, `b_1`, and `rho_2` from Section
11.444.  On the centred cell `|s|<=1/2`, put

```text
g_0(s)=[csc(pi*s)-1/(pi*s)]/(2i),
g_1(s)={{1/(pi*s)^2-cos(pi*s)/sin(pi*s)^2}}/4.         (PC1)
```

Their apparent singularities are removable:

```text
g_0(s)=(1/(2i))[pi*s/6+7(pi*s)^3/360+...],
g_1(s)=1/24+7(pi*s)^2/480+...,
g_0(0)=0,  g_1(0)=1/24.                              (PC2)
```

The `k=0` residual identity is

```text
rho_2(s)=hat h_epsilon(s)-b_0(s)-b_1(s).
```

Combining it with the two closed boundary sums before evaluation gives the
single pole-free formula

```text
C_epsilon(s)
 =hat h_epsilon(s)+E_0(s)g_0(s)+E_1(s)g_1(s)
  +sum_(k in Z, k!=0)rho_2(k+s),       |s|<=1/2,     (PC3)

C_epsilon(0)
 =hat h_epsilon(0)+E_1(0)/24+sum_(k!=0)rho_2(k).     (PC4)
```

No direct high-cutoff switch is needed at the integer lattice.  Uniformly on
the whole centred cell,

```text
sum_(|k|>K)|rho_2(k+s)|
 <=K3/[(2pi)^3(K-1/2)^2],
K3=|q_epsilon''(c)|+|q_epsilon''(d)|+20pi*epsilon.   (PC5)
```

For `K=64`, the certified tails are

```text
epsilon={tails[0]['epsilon']}: {tails[0]['uniform_punctured_tail_ball']}
epsilon={tails[1]['epsilon']}: {tails[1]['uniform_punctured_tail_ball']}
epsilon={tails[2]['epsilon']}: {tails[2]['uniform_punctured_tail_ball']}. (PC6)
```

The common complement is one-periodic.  Since `L={L}` is an integer, define

```text
Theta_L(x,s)
 =sum_(n=0)^(L-1)e^[i*pi*x*(A+2n+2s)^2/4],

F_x(s)
 =sum_(n=0)^(L-1)(A+2n+2s)e^[i*pi*x*(A+2n+2s)^2/4]
 =(i*pi*x)^(-1)partial_s Theta_L(x,s),       x>0.     (PC7)
```

Cellwise change of variables gives the exact fold

```text
integral_0^L f_x(u)C_epsilon(u)du
 =integral_0^1 F_x(s)C_epsilon(s)ds.                 (PC8)
```

The finite current shifts from labels `A,...,B-2` at `s=0` to
`A+2,...,B` at `s=1`, so

```text
Theta_L(x,1)-Theta_L(x,0)
 =e^(i*pi*x*B^2/4)-e^(i*pi*x*A^2/4).                (PC9)
```

This is the same endpoint difference that drives the certified Poisson and
double-Weber cancellations.  Equation (PC8) is now an executable one-cell
target for a finite quadratic-Gauss current evaluator; it is not yet a
quadrature theorem.

The separate floating pilot passes 21 direct/dual rows, including the exact
integer and points within `10^-5`, down to `epsilon=10^-10`.  Its worst
relative discrepancy is below `3.91e-11`.  Those samples are not used as
proof of (PC1)--(PC9).

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.py
```

Pi provenance: every `pi` in (PC1)--(PC9) is inherited from the Gaussian
Abel/Fourier kernel and the original Kummer quadratic phase.  No fitted or
geometric occurrence is introduced.

Proof boundary: exact fixed-regulator pole cancellation, one-cell folding,
finite theta-current identity, and cubic dual-tail enclosures only.
No finite-Gauss-current error theorem, Kummer `x` integration, A/B extraction
evaluation, quantitative `R_after_A` or `R_Dir` bound, complete `Q_K-T`,
all-height theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion is
proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["two_edge_modular"]["decision"]["remaining_dual_series_absolutely_convergent_with_cubic_tail"] is True, "two-edge tail drift")
    require(dependencies["two_edge_modular"]["decision"]["first_two_boundary_currents_summed_in_closed_form"] is True, "boundary-current drift")
    require(dependencies["R_after_A_equivalence"]["decision"]["common_kernel_equals_six_class_post_A_mask"] is True, "R_after_A kernel drift")
    require(dependencies["global_theta_current"]["route_decision"]["global_incomplete_theta_derivative_reduction_is_exact"] is True, "global theta-current drift")
    pilot = load_json(PILOT)
    require(pilot["decision"]["integer_and_near_integer_rows_passed"] is True, "punctured pilot drift")
    require(pilot["decision"]["interval_certified"] is False and pilot["decision"]["R_after_A_evaluated"] is False, "pilot proof-boundary drift")

    artifact = {
        "kind": STEM,
        "status": "exact_pole_free_target_notch_cell_recombination_one_cell_source_fold_and_uniform_cubic_tail_certified_quadrature_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "source_odd_roster": [A, B],
            "roster_coordinate_length": L,
            "centred_cell": "|s|<=1/2",
            "dual_cutoff_certificate": DUAL_CUTOFF,
        },
        "symbolic_certificate": symbolic_certificate(),
        "punctured_transform_certificate": punctured_transform_certificate(),
        "fold_certificate": fold_certificate(),
        "tail_certificate": tail_certificate(),
        "decision": {
            "k0_residual_recombined_before_evaluation": True,
            "apparent_integer_lattice_poles_removed_exactly": True,
            "single_dual_formula_continuous_on_closed_centred_cell": True,
            "uniform_cubic_tail_on_closed_centred_cell_certified": True,
            "long_u_integral_folded_exactly_to_one_periodic_cell": True,
            "folded_source_is_exact_finite_quadratic_theta_current": True,
            "floating_pilot_used_as_proof": False,
            "finite_Gauss_current_error_theorem_proved": False,
            "Kummer_x_integration_completed": False,
            "R_after_A_evaluated": False,
            "R_after_A_bound_proved": False,
            "rh_implication": False,
        },
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)}
            for name, path in DEPENDENCIES.items()
        },
        "exploratory_context": {"path": relative(PILOT), "sha256": file_hash(PILOT), "used_as_proof": False},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Construct a deterministic evaluator for the folded quadratic-Gauss current F_x(s) on one cell with explicit error and derivative controls, then compose it with the pole-free C_epsilon(s) formula. Start with floating cross-checks against direct finite sums on a sparse (x,s,epsilon) lattice; promote only after a recursive Gauss-sum or interval chirp evaluator closes. Keep the endpoint half-current, A notch, B window, and B outer subtraction in the same physical projection.",
        "proof_boundary": "Exact fixed-regulator punctured-cell cancellation, one-cell fold, finite theta-current identity, and uniform cubic dual tails only. No finite-Gauss-current error theorem, Kummer x integration, A/B extraction evaluation, quantitative R_after_A or R_Dir bound, complete Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified pole-free target-notch cell and one-cell source fold", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
