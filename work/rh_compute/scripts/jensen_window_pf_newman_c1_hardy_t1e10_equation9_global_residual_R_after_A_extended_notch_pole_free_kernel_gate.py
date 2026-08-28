#!/usr/bin/env python3
"""Certify the instantiated pole-free kernel for the post-A extended notch."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_pole_free_kernel_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

DEPENDENCIES = {
    "extended_notch": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_normal_form_gate.json",
    "two_edge_modular": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.json",
    "pole_free_cell": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_punctured_cell_fold_gate.json",
    "A_common_phase": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_common_phase_reduction_gate.json",
}

A = 159_577
T = 10_000_000_000
TARGET_START = 622
TARGET_END = 39_936
TARGET_COUNT = TARGET_END - TARGET_START + 1
LEFT_NUMERATOR = 1_243
RIGHT_NUMERATOR = 79_873
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


def exact_kernel_certificate() -> dict[str, Any]:
    c = sp.Rational(LEFT_NUMERATOR, 2)
    d_a = sp.Rational(RIGHT_NUMERATOR, 2)
    require(c == sp.Rational(1243, 2), "left edge drift")
    require(d_a == sp.Rational(79873, 2), "extended upper edge drift")
    require(d_a - c == TARGET_COUNT, "extended-notch width drift")
    require(LEFT_NUMERATOR % 2 == 1 and RIGHT_NUMERATOR % 2 == 1, "half-integer parity drift")

    z = sp.symbols("z", real=True)
    g0_core = sp.csc(z) - 1 / z
    g1 = (1 / z**2 - sp.cos(z) / sp.sin(z) ** 2) / 4
    require(sp.limit(g0_core, z, 0) == 0, "g0 removable limit failed")
    require(sp.limit(g1, z, 0) == sp.Rational(1, 24), "g1 removable limit failed")
    require(sp.series(g0_core, z, 0, 6).coeff(z, 1) == sp.Rational(1, 6), "g0 series drift")
    require(sp.series(4 * g1, z, 0, 6).coeff(z, 2) == sp.Rational(7, 120), "g1 series drift")

    return {
        "extended_target": [TARGET_START, TARGET_END],
        "target_count": TARGET_COUNT,
        "half_integer_edges": [str(c), str(d_a)],
        "edge_width": str(d_a - c),
        "shared_integer_dual_phase": "exp(-2*pi*i*k*c)=exp(-2*pi*i*k*d_A)=(-1)^k",
        "Gaussian": "q_epsilon(y)=exp(-pi*epsilon*y^2)",
        "interpolant": "h_A,epsilon(y)=q_epsilon(y)[1_(y<c)+1_(y>d_A)]",
        "integer_samples": "h_A,epsilon(m)=q_epsilon(m)[1-chi_U(m)], U={622,...,39936}",
        "Poisson_identity": "C_A,epsilon(s)=sum_(m in Z)[1-chi_U(m)]q_epsilon(m)e^(-2*pi*i*m*s)=sum_(k in Z)hat h_A,epsilon(k+s)",
        "Fourier_transform": "hat h_A,epsilon(xi)=e^(-pi*xi^2/epsilon)/(2sqrt(epsilon)){erfc[-sqrt(pi*epsilon)c-i sqrt(pi/epsilon)xi]+erfc[sqrt(pi*epsilon)d_A+i sqrt(pi/epsilon)xi]}",
        "boundary_currents": "b_A,j(xi)=[q_epsilon^(j)(d_A)e^(-2*pi*i*xi*d_A)-q_epsilon^(j)(c)e^(-2*pi*i*xi*c)]/(2*pi*i*xi)^(j+1), j=0,1",
        "edge_numerator": "E_A,j(s)=q_epsilon^(j)(d_A)e^(-2*pi*i*s*d_A)-q_epsilon^(j)(c)e^(-2*pi*i*s*c)",
        "closed_currents": [
            "sum_k b_A,0(k+s)=E_A,0(s)/(2*i*sin(pi*s))",
            "sum_k b_A,1(k+s)=-E_A,1(s)cos(pi*s)/(4*sin(pi*s)^2)",
        ],
        "pole_free_coefficients": {
            "g0": "g_0(s)=[csc(pi*s)-1/(pi*s)]/(2*i)",
            "g1": "g_1(s)={1/(pi*s)^2-cos(pi*s)/sin(pi*s)^2}/4",
            "limits": "g_0(0)=0, g_1(0)=1/24",
        },
        "pole_free_cell": "C_A,epsilon(s)=hat h_A,epsilon(s)+E_A,0(s)g_0(s)+E_A,1(s)g_1(s)+sum_(k in Z,k!=0)rho_A,2(k+s), |s|<=1/2",
        "integer_value": "C_A,epsilon(0)=hat h_A,epsilon(0)+E_A,1(0)/24+sum_(k!=0)rho_A,2(k)",
        "remainder": "rho_A,2=hat h_A,epsilon-b_A,0-b_A,1",
    }


def tail_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    c = arb(LEFT_NUMERATOR) / 2
    d_a = arb(RIGHT_NUMERATOR) / 2
    cutoff = arb(DUAL_CUTOFF)
    rows = []
    for epsilon_text in ("1e-6", "1e-8", "1e-10"):
        epsilon = arb(epsilon_text)

        def q2(edge: arb) -> arb:
            q = (-pi * epsilon * edge**2).exp()
            return (4 * pi**2 * epsilon**2 * edge**2 - 2 * pi * epsilon) * q

        q2_c = q2(c)
        q2_d = q2(d_a)
        k3 = abs(q2_c) + abs(q2_d) + 20 * pi * epsilon
        tail = k3 / (2 * pi) ** 3 / (cutoff - arb(1) / 2) ** 2
        rows.append(
            {
                "epsilon": epsilon_text,
                "centred_cell": "|s|<=1/2",
                "dual_cutoff": DUAL_CUTOFF,
                "q2_lower_edge_ball": q2_c.str(70, more=True),
                "q2_upper_edge_ball": q2_d.str(70, more=True),
                "K3_ball": k3.str(70, more=True),
                "uniform_punctured_tail_ball": tail.str(70, more=True),
            }
        )
    require(arb(rows[-1]["uniform_punctured_tail_ball"]).upper() < arb("7.1e-15"), "epsilon=1e-10 extended-notch tail drift")
    return {
        "pointwise_bound": "|rho_A,2(xi)|<=[|q_epsilon''(c)|+|q_epsilon''(d_A)|+20*pi*epsilon]/(2*pi*|xi|)^3",
        "uniform_tail": "sum_(|k|>K)|rho_A,2(k+s)|<=K3_A/[(2*pi)^3(K-1/2)^2], |s|<=1/2",
        "rows": rows,
    }


def stationary_guard_certificate() -> dict[str, Any]:
    x = sp.symbols("x", positive=True)
    phase = sp.pi * A**2 * x / 4 + sp.Rational(T, 2) * sp.log((1 - x) / x)
    derivative = sp.simplify(sp.diff(phase, x))
    expected = sp.pi * A**2 / 4 - sp.Rational(T, 2) / (x * (1 - x))
    require(sp.simplify(derivative - expected) == 0, "common A-phase derivative drift")

    root = (1 - sp.sqrt(1 - sp.Rational(8 * T, 1) / (sp.pi * A**2))) / 2
    require(sp.simplify(root * (1 - root) - sp.Rational(2 * T, 1) / (sp.pi * A**2)) == 0, "stationary root identity drift")

    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    discriminant = 1 - arb(8 * T) / (pi * arb(A) ** 2)
    require(discriminant.lower() > arb(0), "A stationary discriminant not positive")
    x_star = (1 - discriminant.sqrt()) / 2
    require(x_star.lower() > arb(0) and x_star.upper() < arb(1) / 2, "A stationary point left half-domain")
    curvature = arb(T) * (1 - 2 * x_star) / (2 * (x_star * (1 - x_star)) ** 2)
    require(curvature.lower() > arb(0), "A stationary curvature not positive")

    return {
        "common_phase": "Phi_A(x)=pi*A^2*x/4+(t/2)log((1-x)/x)",
        "phase_derivative": str(expected),
        "stationary_equation": "x*(1-x)=2*t/(pi*A^2)",
        "lower_half_stationary_point": "x_*=[1-sqrt(1-8*t/(pi*A^2))]/2",
        "x_star_ball": x_star.str(80, more=True),
        "curvature_ball": curvature.str(80, more=True),
        "edge_independence": "Phi_A and x_* depend on A and t, not on the notch edge d_A.",
        "route_guard": "Moving d to d_A removes the A half-boundary collar from the projector edge, but it does not remove the interior stationary point of the surviving common A phase.",
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = artifact["tail_certificate"]["rows"]
    saddle = artifact["stationary_guard_certificate"]
    return f"""# Post-A extended-notch pole-free kernel

Date: 2026-08-23

Status: exact instantiated two-edge transform, pole-free centred-cell formula,
and cubic dual-tail certificate; the joined stationary estimate remains open

After the exact A-window subtraction, put

```text
U={{622,...,39936}},       c=621.5,       d_A=39936.5,
q_epsilon(y)=exp(-pi*epsilon*y^2),
h_A,epsilon(y)=q_epsilon(y)[1_(y<c)+1_(y>d_A)].       (EA1)
```

The edge width is exactly `d_A-c=39315=|U|`.  Both edge numerators are odd,
so for every integer `k`

```text
exp(-2*pi*i*k*c)=exp(-2*pi*i*k*d_A)=(-1)^k.          (EA2)
```

Thus the extended complement has the exact fixed-regulator Poisson form

```text
C_A,epsilon(s)
 =sum_(m in Z)[1-chi_U(m)]q_epsilon(m)e^(-2*pi*i*m*s)
 =sum_(k in Z)hat h_A,epsilon(k+s),                  (EA3)

hat h_A,epsilon(xi)
 =e^(-pi*xi^2/epsilon)/(2sqrt(epsilon))
  {{erfc[-sqrt(pi*epsilon)c-i sqrt(pi/epsilon)xi]
   +erfc[ sqrt(pi*epsilon)d_A+i sqrt(pi/epsilon)xi]}}. (EA4)
```

Define, for `j=0,1`,

```text
b_A,j(xi)
 =[q_epsilon^(j)(d_A)e^(-2*pi*i*xi*d_A)
   -q_epsilon^(j)(c)e^(-2*pi*i*xi*c)]
   /(2*pi*i*xi)^(j+1),

E_A,j(s)=q_epsilon^(j)(d_A)e^(-2*pi*i*s*d_A)
         -q_epsilon^(j)(c)e^(-2*pi*i*s*c).           (EA5)
```

With `rho_A,2=hat h_A,epsilon-b_A,0-b_A,1`, the same `k=0`
recombination as in Section 11.445 gives one formula on the closed centred
cell:

```text
C_A,epsilon(s)
 =hat h_A,epsilon(s)+E_A,0(s)g_0(s)+E_A,1(s)g_1(s)
  +sum_(k in Z,k!=0)rho_A,2(k+s),       |s|<=1/2,    (EA6)

g_0(s)=[csc(pi*s)-1/(pi*s)]/(2i),
g_1(s)={{1/(pi*s)^2-cos(pi*s)/sin(pi*s)^2}}/4,
g_0(0)=0,       g_1(0)=1/24.                        (EA7)
```

Three integrations by parts on the two exterior intervals give

```text
|rho_A,2(xi)|
 <=K3_A/(2*pi*|xi|)^3,
K3_A=|q_epsilon''(c)|+|q_epsilon''(d_A)|+20*pi*epsilon,

sum_(|k|>K)|rho_A,2(k+s)|
 <=K3_A/[(2*pi)^3(K-1/2)^2],       |s|<=1/2.        (EA8)
```

At `K=64`, the fresh Arb enclosures are

```text
epsilon={rows[0]['epsilon']}: {rows[0]['uniform_punctured_tail_ball']}
epsilon={rows[1]['epsilon']}: {rows[1]['uniform_punctured_tail_ball']}
epsilon={rows[2]['epsilon']}: {rows[2]['uniform_punctured_tail_ball']}. (EA9)
```

The 42.25-mode edge displacement is useful, but it is not a nonstationary
certificate.  The surviving exact A common phase is

```text
Phi_A(x)=pi*A^2*x/4+(t/2)log((1-x)/x),
x_*=[1-sqrt(1-8*t/(pi*A^2))]/2
   ={saddle['x_star_ball']}.                          (EA10)
```

Its positive curvature is `{saddle['curvature_ball']}`.  Neither `Phi_A`
nor `x_*` contains `d_A`.  Moving the projector edge out of the A
half-boundary collar therefore does not remove this interior stationary
point.  The upper-edge current and the already-deleted A window must be
estimated in one phase-adapted Morse/Fresnel assembly before norms; blanket
nonstationary integration by parts is invalid.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (EA1)--(EA10) is inherited from the Gaussian
Abel/Fourier kernel and the exact Kummer common phase.  No fitted or geometric
occurrence is introduced.

Proof boundary: exact fixed-regulator extended-notch transform, pole-free
cell recombination, cubic dual tails, and common-phase stationary guard only.
No Kummer-weighted cell integral, joined phase-adapted remainder estimate,
quantitative `R_after_A`, `R_Dir`, or `Q_K-T` bound, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["extended_notch"]["decision"]["extended_notch_one_cell_fold_and_modular_transform_inherited"] is True, "extended-notch inheritance drift")
    require(dependencies["two_edge_modular"]["decision"]["remaining_dual_series_absolutely_convergent_with_cubic_tail"] is True, "two-edge cubic tail drift")
    require(dependencies["pole_free_cell"]["decision"]["apparent_integer_lattice_poles_removed_exactly"] is True, "pole-free cell drift")
    require(dependencies["A_common_phase"]["decision"]["full_exact_exterior_common_phase_reduction_proved"] is True, "A common-phase reduction drift")
    require(dependencies["A_common_phase"]["decision"]["ten_mode_incomplete_stationary_roster_certified"] is True, "A common-phase saddle drift")

    artifact = {
        "kind": STEM,
        "status": "exact_post_A_extended_notch_pole_free_kernel_and_cubic_tail_certified_common_phase_stationary_guard_active",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "extended_target": [TARGET_START, TARGET_END],
            "centred_cell": "|s|<=1/2",
            "dual_cutoff": DUAL_CUTOFF,
        },
        "exact_kernel_certificate": exact_kernel_certificate(),
        "tail_certificate": tail_certificate(),
        "stationary_guard_certificate": stationary_guard_certificate(),
        "decision": {
            "extended_upper_edge_instantiated_exactly": True,
            "shared_half_integer_dual_phase_certified": True,
            "pole_free_closed_centred_cell_formula_certified": True,
            "uniform_cubic_dual_tail_certified": True,
            "A_half_boundary_projector_collar_removed": True,
            "surviving_common_A_phase_has_interior_stationary_point": True,
            "blanket_nonstationary_IBP_allowed": False,
            "Kummer_weighted_cell_integral_evaluated": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
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
            "arb_threads": 1,
            "process_priority": priority,
            "precision_decimal_digits": PRECISION,
        },
        "next_obligation": "Build a bounded floating pilot for the Kummer-weighted post-A upper-edge current using the exact common A phase and its Morse coordinate, retaining the A-window deletion and B-window ownership in the same signed assembly. Use it only to choose the rigorous interval decomposition; then certify the compact stationary core and exterior tails separately.",
        "proof_boundary": "Exact fixed-regulator extended-notch transform, pole-free cell recombination, cubic dual tails, and common-phase stationary guard only. No Kummer-weighted cell integral, joined phase-adapted remainder estimate, quantitative R_after_A, R_Dir, or Q_K-T bound, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified post-A extended-notch pole-free kernel", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
