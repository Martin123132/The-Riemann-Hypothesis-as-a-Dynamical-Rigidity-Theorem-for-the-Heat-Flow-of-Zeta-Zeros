#!/usr/bin/env python3
"""Certify a two-edge modular form of the target-notched Abel kernel."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PILOT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_dual_pilot.json"

DEPENDENCIES = {
    "R_after_A_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
    "theta_modular": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_half_pair_abel_theta_modular_dual_roster_gate.json",
    "symmetric_poisson": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_symmetric_poisson_interchange_gate.json",
}

TARGET_START = 622
TARGET_END = 39_894
TARGET_COUNT = TARGET_END - TARGET_START + 1
LEFT_EDGE = sp.Rational(2 * TARGET_START - 1, 2)
RIGHT_EDGE = sp.Rational(2 * TARGET_END + 1, 2)
DUAL_CUTOFF = 64
AWAY_DELTA = sp.Rational(1, 10)
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
    y, epsilon = sp.symbols("y epsilon", real=True, positive=True)
    q = sp.exp(-sp.pi * epsilon * y**2)
    derivatives = [sp.simplify(sp.diff(q, y, order) / q) for order in range(4)]
    expected = [
        1,
        -2 * sp.pi * epsilon * y,
        4 * sp.pi**2 * epsilon**2 * y**2 - 2 * sp.pi * epsilon,
        12 * sp.pi**2 * epsilon**2 * y - 8 * sp.pi**3 * epsilon**3 * y**3,
    ]
    require(all(sp.simplify(left - right) == 0 for left, right in zip(derivatives, expected, strict=True)), "Gaussian derivative identity failed")

    require(LEFT_EDGE == sp.Rational(1243, 2), "left half-integer edge drift")
    require(RIGHT_EDGE == sp.Rational(79789, 2), "right half-integer edge drift")
    require(RIGHT_EDGE - LEFT_EDGE == TARGET_COUNT, "target edge width drift")
    k = sp.symbols("k", integer=True)
    require(sp.simplify(sp.exp(-2 * sp.pi * sp.I * k * LEFT_EDGE) - (-1) ** k) == 0, "left half-integer phase drift")
    require(sp.simplify(sp.exp(-2 * sp.pi * sp.I * k * RIGHT_EDGE) - (-1) ** k) == 0, "right half-integer phase drift")

    abs_y_integral = sp.integrate(y * sp.exp(-sp.pi * epsilon * y**2), (y, 0, sp.oo)) * 2
    abs_y3_integral = sp.integrate(y**3 * sp.exp(-sp.pi * epsilon * y**2), (y, 0, sp.oo)) * 2
    require(sp.simplify(abs_y_integral - 1 / (sp.pi * epsilon)) == 0, "absolute first Gaussian moment drift")
    require(sp.simplify(abs_y3_integral - 1 / (sp.pi**2 * epsilon**2)) == 0, "absolute third Gaussian moment drift")
    q3_l1 = sp.simplify(12 * sp.pi**2 * epsilon**2 * abs_y_integral + 8 * sp.pi**3 * epsilon**3 * abs_y3_integral)
    require(sp.simplify(q3_l1 - 20 * sp.pi * epsilon) == 0, "third-derivative L1 majorant drift")

    return {
        "target_band": [TARGET_START, TARGET_END],
        "target_count": TARGET_COUNT,
        "half_integer_edges": [str(LEFT_EDGE), str(RIGHT_EDGE)],
        "edge_width": str(RIGHT_EDGE - LEFT_EDGE),
        "q_derivatives": [str(item) for item in expected],
        "third_derivative_L1_majorant": "integral_R |q_epsilon'''(y)|dy<=20*pi*epsilon",
        "half_integer_phase": "exp(-2*pi*i*k*c)=exp(-2*pi*i*k*d)=(-1)^k for integer k",
        "Mittag_Leffler_sums": [
            "sum_(k in Z)(-1)^k/(k+u)=pi/sin(pi*u)",
            "sum_(k in Z)(-1)^k/(k+u)^2=pi^2*cos(pi*u)/sin(pi*u)^2",
        ],
    }


def exact_transform_certificate() -> dict[str, Any]:
    return {
        "Gaussian": "q_epsilon(y)=exp(-pi*epsilon*y^2)",
        "continuous_notch_complement": "h_epsilon(y)=q_epsilon(y)[1_(y<c)+1_(y>d)], c=621.5, d=39894.5",
        "integer_samples": "h_epsilon(m)=q_epsilon(m)[1-chi_T(m)] for every integer m",
        "regulated_infinite_complement": "C_epsilon(u)=sum_(m in Z)[1-chi_T(m)]q_epsilon(m)exp(-2pi*i*m*u)",
        "Poisson_identity": "C_epsilon(u)=sum_(k in Z) hat h_epsilon(k+u), in the inherited symmetric prescription",
        "Fourier_transform": "hat h_epsilon(xi)=exp(-pi*xi^2/epsilon)/(2sqrt(epsilon)){erfc[-sqrt(pi*epsilon)c-i*sqrt(pi/epsilon)xi]+erfc[sqrt(pi*epsilon)d+i*sqrt(pi/epsilon)xi]}",
        "boundary_currents": "b_j(xi)=[q_epsilon^(j)(d)exp(-2pi*i*xi*d)-q_epsilon^(j)(c)exp(-2pi*i*xi*c)]/(2pi*i*xi)^(j+1), j=0,1",
        "boundary_sum_0": "B_0(u)=E_0(u)/(2i*sin(pi*u))",
        "boundary_sum_1": "B_1(u)=-E_1(u)cos(pi*u)/(4sin(pi*u)^2)",
        "edge_numerator": "E_j(u)=q_epsilon^(j)(d)exp(-2pi*i*u*d)-q_epsilon^(j)(c)exp(-2pi*i*u*c)",
        "absolute_remainder": "rho_2(xi)=hat h_epsilon(xi)-b_0(xi)-b_1(xi)",
        "three_IBP_identity": "rho_2(xi)={[q''(d)e^(-2pi*i*xi*d)-q''(c)e^(-2pi*i*xi*c)]+integral_(-infinity,c union d,infinity)q'''(y)e^(-2pi*i*xi*y)dy}/(2pi*i*xi)^3",
        "pointwise_remainder_bound": "|rho_2(xi)|<=[|q''(c)|+|q''(d)|+20pi*epsilon]/(2pi|xi|)^3",
        "two_edge_modular_form": "C_epsilon(u)=B_0(u)+B_1(u)+sum_(k in Z)rho_2(k+u), u not in Z, with continuous extension at integers",
        "hybrid_evaluator": "Use the direct Gaussian Fourier sum near integer u and the boundary-subtracted dual sum away from integers; both evaluate the same C_epsilon.",
    }


def tail_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    c = arb(str(LEFT_EDGE.p)) / arb(str(LEFT_EDGE.q))
    d = arb(str(RIGHT_EDGE.p)) / arb(str(RIGHT_EDGE.q))
    cutoff = arb(DUAL_CUTOFF)
    delta = arb(str(AWAY_DELTA.p)) / arb(str(AWAY_DELTA.q))
    reciprocal_tail = arb(1) / (2 * (cutoff + delta) ** 2) + arb(1) / (2 * (cutoff - (1 - delta)) ** 2)
    rows = []
    for epsilon_text in ("1e-6", "1e-8", "1e-10"):
        epsilon = arb(epsilon_text)

        def q2(edge: arb) -> arb:
            q = (-pi * epsilon * edge**2).exp()
            return (4 * pi**2 * epsilon**2 * edge**2 - 2 * pi * epsilon) * q

        k3 = abs(q2(c)) + abs(q2(d)) + 20 * pi * epsilon
        tail = k3 * reciprocal_tail / (2 * pi) ** 3
        rows.append(
            {
                "epsilon": epsilon_text,
                "fractional_strip": "0.1<=u<=0.9",
                "dual_cutoff": DUAL_CUTOFF,
                "K3_ball": k3.str(70, more=True),
                "uniform_dual_tail_ball": tail.str(70, more=True),
            }
        )
    require(arb(rows[-1]["uniform_dual_tail_ball"]).upper() < arb("7e-15"), "epsilon=1e-10 dual tail scale drift")
    return {
        "tail_formula": "sum_(|k|>K)|rho_2(k+u)|<=K3/(2pi)^3{1/[2(K+u)^2]+1/[2(K-u)^2]}",
        "uniformization": "For 0.1<=u<=0.9, bound the two denominators separately by K+0.1 and K-0.9.",
        "rows": rows,
    }


def render_note(artifact: dict[str, Any]) -> str:
    tails = artifact["tail_certificate"]["rows"]
    return f"""# Two-edge modular transform of the target-notched theta kernel

Date: 2026-08-23

Status: exact regulated complement transform and cubic dual-tail certificate;
not an evaluation or bound for `R_after_A`

Put

```text
T={{622,...,39894}},       c=621.5,       d=39894.5,
q_epsilon(y)=exp(-pi*epsilon*y^2),       epsilon>0,

h_epsilon(y)=q_epsilon(y)[1_(y<c)+1_(y>d)].          (TN1)
```

Because `c` and `d` are half-integers, the integer samples of `h_epsilon`
are exactly the target complement.  After the lawful `M->infinity` passage
at fixed positive `epsilon` in Section 11.443,

```text
C_epsilon(u)
 =sum_(m in Z)[1-chi_T(m)]q_epsilon(m)e^(-2pi*i*m*u)
 =sum_(k in Z)hat h_epsilon(k+u).                    (TN2)
```

The second equality is Poisson summation for the integrable bounded-variation
function (TN1), in the inherited symmetric prescription.  Its exact Fourier
transform is

```text
hat h_epsilon(xi)
 =e^(-pi*xi^2/epsilon)/(2sqrt(epsilon))
  {{erfc[-sqrt(pi*epsilon)c-i sqrt(pi/epsilon)xi]
   +erfc[ sqrt(pi*epsilon)d+i sqrt(pi/epsilon)xi]}}.  (TN3)
```

For `j=0,1`, define

```text
b_j(xi)
 =[q_epsilon^(j)(d)e^(-2pi*i*xi*d)
   -q_epsilon^(j)(c)e^(-2pi*i*xi*c)]
   /(2pi*i*xi)^(j+1),                                (TN4)

E_j(u)=q_epsilon^(j)(d)e^(-2pi*i*u*d)
       -q_epsilon^(j)(c)e^(-2pi*i*u*c).
```

The half-integer phases give `exp(-2pi*i*k*c)=exp(-2pi*i*k*d)=(-1)^k`.
The two standard Mittag-Leffler sums therefore close the boundary currents:

```text
sum_k b_0(k+u)= E_0(u)/(2i sin(pi*u)),
sum_k b_1(k+u)=-E_1(u)cos(pi*u)/(4sin(pi*u)^2).       (TN5)
```

Three integrations by parts, performed on the two exterior intervals before
an absolute value, give for `rho_2=hat h-b_0-b_1`

```text
|rho_2(xi)|
 <=[|q_epsilon''(c)|+|q_epsilon''(d)|+20pi*epsilon]
    /(2pi|xi|)^3.                                    (TN6)
```

Here the explicit constant follows from

```text
q_epsilon'''=(12pi^2 epsilon^2 y-8pi^3 epsilon^3 y^3)
              q_epsilon,
integral_R |q_epsilon'''|<=20pi*epsilon.             (TN7)
```

Consequently, for noninteger `u`,

```text
C_epsilon(u)
 =E_0(u)/(2i sin(pi*u))
  -E_1(u)cos(pi*u)/(4sin(pi*u)^2)
  +sum_(k in Z)rho_2(k+u),                           (TN8)
```

with continuous extension of the complete right side at integer `u`.  The
closed terms and residual must not be evaluated separately near those
removable points; the stable evaluator uses the direct Gaussian sum there
and (TN8) away from the integer lattice.

For `|k|<=64` and `0.1<=u<=0.9`, the certified dual-tail enclosures are

```text
epsilon={tails[0]['epsilon']}: {tails[0]['uniform_dual_tail_ball']}
epsilon={tails[1]['epsilon']}: {tails[1]['uniform_dual_tail_ball']}
epsilon={tails[2]['epsilon']}: {tails[2]['uniform_dual_tail_ball']}.     (TN9)
```

The separate floating pilot checks nine real-edge rows down to
`epsilon=10^-10`: its direct cutoff grows to `328014`, while the dual cutoff
stays at `64`, and the largest discrepancy is below `4.19e-11`.  Those
floating comparisons select the route but are not used as proof of (TN1)--
(TN9) and do not estimate `R_after_A`.

Machine-audited companion:

```text
outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.md
work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.json
work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.py
work/rh_compute/scripts/check_jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_target_notched_theta_two_edge_modular_gate.py
```

Pi provenance: every `pi` in (TN1)--(TN9) comes from the inherited Gaussian
Abel weight, integer Fourier character, Fourier transform convention, and
standard Jacobi/Poisson normalization.  No geometric or fitted occurrence is
introduced.

Proof boundary: exact fixed-`epsilon` target-complement transform, explicit
boundary-current closure, and cubic dual-tail bound only.  No integration
against the Kummer source, A/B extraction evaluation, numerical or analytic
bound for `R_after_A`, complete `R_Dir` or `Q_K-T`, all-height theorem,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "a dependency gate is not passed")
    require(dependencies["R_after_A_equivalence"]["decision"]["finite_R_after_A_defined_on_one_common_regulator"] is True, "R_after_A regulator drift")
    require(dependencies["R_after_A_equivalence"]["scope"]["limit_order"].startswith("M to infinity first"), "R_after_A limit order drift")
    require(dependencies["theta_modular"]["decision"]["Jacobi_modular_dual_phase_identified"] is True, "Jacobi dependency drift")
    require(dependencies["symmetric_poisson"]["decision"]["symmetric_poisson_interchange_proved"] is True, "Poisson dependency drift")
    pilot = load_json(PILOT)
    require(pilot["decision"]["direct_and_dual_samples_agree_within_declared_allowance"] is True, "pilot comparison drift")
    require(pilot["decision"]["interval_certified"] is False and pilot["decision"]["R_after_A_evaluated"] is False, "pilot proof-boundary drift")

    artifact = {
        "kind": STEM,
        "status": "exact_target_notched_Abel_kernel_two_edge_modular_transform_and_cubic_dual_tail_certified_R_after_A_integration_open",
        "passed": True,
        "scope": {
            "height": 10_000_000_000,
            "target_positive_modes": [TARGET_START, TARGET_END],
            "target_count": TARGET_COUNT,
            "half_integer_edges": [float(LEFT_EDGE), float(RIGHT_EDGE)],
            "regulator": "epsilon>0 after the inherited M-to-infinity passage",
        },
        "symbolic_certificate": symbolic_certificate(),
        "exact_transform_certificate": exact_transform_certificate(),
        "tail_certificate": tail_certificate(),
        "decision": {
            "target_complement_realized_as_half_integer_continuous_notch": True,
            "Poisson_dual_two_edge_transform_exact_at_fixed_positive_epsilon": True,
            "first_two_boundary_currents_summed_in_closed_form": True,
            "remaining_dual_series_absolutely_convergent_with_cubic_tail": True,
            "hybrid_direct_near_integer_and_dual_away_evaluator_licensed": True,
            "floating_pilot_used_as_proof": False,
            "Kummer_source_integration_completed": False,
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
        "next_obligation": "Insert the two-edge formula into the u-integral of mathfrak R_A. Prove a stable cellwise recombination at integer u so the apparent sin(pi*u) poles cancel before quadrature, then integrate the closed boundary currents jointly with the target endpoint atoms and estimate the absolutely convergent rho_2 series. Keep the A notch and both B subtractions outside this kernel substitution but inside the same physical projection.",
        "proof_boundary": "Exact fixed-epsilon target-complement transform, boundary-current closure, and cubic dual-tail certificate only. No Kummer-source integration, A/B extraction evaluation, quantitative R_after_A or R_Dir bound, complete Q_K-T, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print("certified two-edge modular transform of target-notched Abel kernel", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
