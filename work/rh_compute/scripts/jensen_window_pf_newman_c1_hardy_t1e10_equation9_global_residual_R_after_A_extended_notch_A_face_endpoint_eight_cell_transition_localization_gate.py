#!/usr/bin/env python3
"""Localize the A-face endpoint transition to exactly eight midpoint cells."""

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

from flint import acb, arb, ctx
import sympy as sp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_eight_cell_transition_localization_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PILOT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_endpoint_layer_pilot.json"

DEPENDENCIES = {
    "midpoint_partition": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_midpoint_Peano_partition_gate.json",
    "six_current_defect": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.json",
    "finite_regulator_equivalence": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_finite_regulator_equivalence_gate.json",
}

A = 159_577
T = 10_000_000_000
FULL_LEFT_NUMERATOR = 79_789
TAIL_LEFT_NUMERATOR = 79_805
RIGHT_NUMERATOR = 79_873
TRANSITION_FIRST_MODE = 39_895
TRANSITION_LAST_MODE = 39_902
TAIL_FIRST_MODE = 39_903
TAIL_LAST_MODE = 39_936
Q_CUT = 16
TERMS = 6
PRECISION = 110


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


def odd_double_factorial(index: int) -> int:
    result = 1
    for value in range(1, index + 1, 2):
        result *= value
    return result


def symbolic_certificate() -> dict[str, Any]:
    full_left = sp.Rational(FULL_LEFT_NUMERATOR, 2)
    tail_left = sp.Rational(TAIL_LEFT_NUMERATOR, 2)
    right = sp.Rational(RIGHT_NUMERATOR, 2)
    transition_modes = tuple(range(TRANSITION_FIRST_MODE, TRANSITION_LAST_MODE + 1))
    tail_modes = tuple(range(TAIL_FIRST_MODE, TAIL_LAST_MODE + 1))
    require(len(transition_modes) == 8 and tail_left - full_left == 8, "transition cardinality drift")
    require(len(tail_modes) == 34 and right - tail_left == 34, "tail cardinality drift")
    require(all(sp.Rational(mode) == full_left + sp.Rational(2 * index + 1, 2) for index, mode in enumerate(transition_modes)), "transition midpoint roster drift")
    require(all(sp.Rational(mode) == tail_left + sp.Rational(2 * index + 1, 2) for index, mode in enumerate(tail_modes)), "tail midpoint roster drift")

    endpoint_gap = tail_left - sp.Rational(A, 4)
    endpoint_q_size = 2 * endpoint_gap
    require(endpoint_gap == sp.Rational(33, 4), "tail endpoint gap drift")
    require(endpoint_q_size == sp.Rational(33, 2), "tail q floor drift")

    x, y = sp.symbols("x y", positive=True, real=True)
    delta = y - sp.Rational(A, 2) * x
    error_shape = y * x**TERMS / delta ** (2 * TERMS + 1)
    dx_positive = x ** (TERMS - 1) * y * (TERMS * y + sp.Rational(A * (TERMS + 1), 2) * x) / delta ** (2 * TERMS + 2)
    dy_negative = -x**TERMS * (2 * TERMS * y + sp.Rational(A, 2) * x) / delta ** (2 * TERMS + 2)
    require(sp.simplify(sp.diff(error_shape, x) - dx_positive) == 0, "endpoint-tail x monotonicity drift")
    require(sp.simplify(sp.diff(error_shape, y) - dy_negative) == 0, "endpoint-tail y monotonicity drift")
    physical_weight = x ** (-sp.Rational(5, 4)) * (1 - x) ** (-sp.Rational(1, 4))
    require(sp.simplify(sp.diff(physical_weight, x) - physical_weight * (6 * x - 5) / (4 * x * (1 - x))) == 0, "physical weight monotonicity drift")

    return {
        "full_strip": "[39894.5,39936.5] with midpoint modes 39895,...,39936",
        "transition_strip": "[39894.5,39902.5] with exactly eight midpoint modes 39895,...,39902",
        "tail_strip": "[39902.5,39936.5] with exactly 34 midpoint modes 39903,...,39936",
        "tail_normal_floor": "min(-q_A)=2*(39902.5-A/4)=33/2=16.5 on x_16<=x<=1/2",
        "tail_six_current_reduction": "D_tail=D_tail,6+R_tail,6 with D_tail,6 a finite rational-log midpoint defect",
        "tail_error_accounting": "|R_tail,6|<=68*sup_tail |G_A-G_6|; 68=34 point masses+strip length 34",
        "physical_weight_monotonicity": "x^(-5/4)(1-x)^(-1/4) decreases for 0<x<5/6, so its endpoint-layer maximum is at x_16.",
        "exact_split": "D_42=D_transition,8+D_tail,34 and E_42=-D_42",
        "remaining_obligation": "Only the compact eight-cell exact Fresnel transition requires special-function enclosure on the endpoint layer.",
    }


def g_six(y: arb, x: arb) -> acb:
    i_pi = acb(0, arb.pi())
    delta = y - arb(A) * x / 2
    value = acb(arb(A) * x) / (2 * i_pi * delta)
    for index in range(1, TERMS):
        value += acb(y) * arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / (i_pi ** (index + 1) * delta ** (2 * index + 1))
    return value


def integral_g_six(left: arb, right: arb, x: arb) -> acb:
    i_pi = acb(0, arb.pi())
    a = arb(A) * x / 2
    dl = left - a
    dr = right - a
    value = acb(arb(A) * x * (dr / dl).log()) / (2 * i_pi)
    for index in range(1, TERMS):
        coefficient = arb(odd_double_factorial(2 * index - 1)) * (x / 2) ** index / i_pi ** (index + 1)
        value += coefficient * (
            (dr ** (1 - 2 * index) - dl ** (1 - 2 * index)) / (1 - 2 * index)
            - a * (dr ** (-2 * index) - dl ** (-2 * index)) / (2 * index)
        )
    return value


def tail_defect_six(x: arb) -> acb:
    left = arb(TAIL_LEFT_NUMERATOR) / 2
    right = arb(RIGHT_NUMERATOR) / 2
    value = acb(0)
    for mode in range(TAIL_FIRST_MODE, TAIL_LAST_MODE + 1):
        value += g_six(arb(mode), x)
    return value - integral_g_six(left, right, x)


def numerical_certificate() -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    pi = arb.pi()
    a = arb(A)
    t = arb(T)
    full_left = arb(FULL_LEFT_NUMERATOR) / 2
    tail_left = arb(TAIL_LEFT_NUMERATOR) / 2
    q_cut = arb(Q_CUT)
    sqrt_x_cut = ((2 * q_cut**2 + 8 * a * full_left).sqrt() - q_cut * arb(2).sqrt()) / (2 * a)
    x_cut = sqrt_x_cut**2
    endpoint = arb(1) / 2
    delta_tail_endpoint = tail_left - a / 4

    pointwise_tail_error = (
        2
        * arb(odd_double_factorial(11))
        * tail_left
        * (endpoint / 2) ** TERMS
        / (pi**7 * delta_tail_endpoint**13)
    )
    defect_tail_error = 68 * pointwise_tail_error
    require(pointwise_tail_error.upper() < arb("8.177e-11"), "endpoint tail pointwise error drift")
    require(defect_tail_error.upper() < arb("5.56e-9"), "endpoint tail defect error drift")

    normalization = 2 * (pi / (32 * t)) ** (arb(1) / 4)
    weight_at_cut = 1 / (x_cut * (x_cut * (1 - x_cut)) ** (arb(1) / 4))
    width = endpoint - x_cut
    physical_replacement_error = normalization * width * weight_at_cut * defect_tail_error
    require(physical_replacement_error.upper() < arb("6e-15"), "physical endpoint-tail replacement error drift")

    tail_at_cut = tail_defect_six(x_cut)
    tail_at_endpoint = tail_defect_six(endpoint)

    def complex_record(value: acb) -> dict[str, str]:
        return {
            "real_ball": value.real.str(80, more=True),
            "imag_ball": value.imag.str(80, more=True),
            "modulus_ball": abs(value).str(80, more=True),
        }

    return {
        "x_16_ball": x_cut.str(80, more=True),
        "endpoint_layer_width_ball": width.str(80, more=True),
        "tail_minimum_abs_q": "33/2",
        "uniform_tail_pointwise_G_minus_G6_ball": pointwise_tail_error.str(80, more=True),
        "uniform_tail_defect_error_ball": defect_tail_error.str(80, more=True),
        "physical_tail_replacement_error_ball": physical_replacement_error.str(80, more=True),
        "tail_D6_at_x16": complex_record(tail_at_cut),
        "tail_D6_at_endpoint": complex_record(tail_at_endpoint),
    }


def render_note(artifact: dict[str, Any]) -> str:
    numerical = artifact["numerical_certificate"]
    telemetry = artifact["floating_route_telemetry"]
    return f"""# Eight-cell localization of the A-face endpoint transition

Date: 2026-08-23

Status: exact certificate localizing the special-function endpoint layer to
eight midpoint cells; the eight-cell interval integral remains open

Split the 42-cell strip exactly as

```text
[39894.5,39936.5]
 =[39894.5,39902.5] union [39902.5,39936.5].          (EL1)
```

The first interval contains precisely the midpoint modes `39895,...,39902`;
the second contains precisely `39903,...,39936`.  Therefore

```text
D_42(x)=D_(tr,8)(x)+D_(tail,34)(x),
E_42(x)=-D_42(x).                                    (EL2)
```

On the complete endpoint layer `x_16<=x<=1/2`, the smallest normal-tail
coordinate in the 34-cell block occurs at `(y,x)=(39902.5,1/2)` and is

```text
-q_A=2*(39902.5-A/4)=33/2=16.5.                     (EL3)
```

Hence the six-current theorem applies uniformly to all 34 tail cells.  Its
pointwise error increases with `x` and decreases with `y`, giving

```text
sup_tail |G_A-G_6|
 <={numerical['uniform_tail_pointwise_G_minus_G6_ball']},

|D_(tail,34)-D_(tail,34),6|
 <={numerical['uniform_tail_defect_error_ball']}.     (EL4)
```

The factor in the second line is `68=34+34`: 34 point masses plus strip
length 34.  The reduced tail defect is again the explicit rational-log
formula of Section 11.454, now restricted to the 34-mode roster.

After the exact equation-(9) physical normalization, replacing the endpoint
tail by its six-current form costs at most

```text
{numerical['physical_tail_replacement_error_ball']}. (EL5)
```

At the layer boundaries the rational-log tail itself has certified moduli

```text
|D_(tail,34),6(x_16)|
 ={numerical['tail_D6_at_x16']['modulus_ball']},

|D_(tail,34),6(1/2)|
 ={numerical['tail_D6_at_endpoint']['modulus_ball']}. (EL6)
```

Only

```text
D_(tr,8)(x)
 =sum_(m=39895)^39902 G_A(m,x)
  -integral_(39894.5)^39902.5 G_A(y,x)dy             (EL7)
```

still requires exact Fresnel enclosure on the endpoint layer.  This is a
compact rectangle in the normal variables: eight unit cells and
`1/4<=39894.5-A*x/2<=7.9993`.  No remote mode or 34-cell tail must enter that
special-function calculation.

Floating route telemetry, not used as proof, gives the complete endpoint
layer candidate `{telemetry['candidate_physical_endpoint_layer']}` with
coarse/fine relative-integral difference
`{telemetry['coarse_fine_difference']:.4e}`.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (EL1)--(EL7) is inherited from the canonical
Fresnel primitive, Gaussian Abel factor, equation-(9) normalization, and
exact Kummer/Fourier phase.  The integer eight-cell split is forced by the
next half-integer edge with endpoint `|q_A|>16`, not fitted to telemetry.

Proof boundary: exact eight-cell/34-cell endpoint decomposition, uniform
six-current tail replacement below `5.56e-9`, and physical replacement cost
below `6e-15` only.  No interval value for the eight-cell transition,
complete endpoint-layer value, full signed A-face bound, quantitative
`R_after_A`, `R_Dir`, or `Q_K-T` estimate, all-height theorem, `Lambda<=0`,
PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependencies = {name: load_json(path) for name, path in DEPENDENCIES.items()}
    require(all(item.get("passed") is True for item in dependencies.values()), "an exact dependency gate is not passed")
    require(dependencies["midpoint_partition"]["decision"]["endpoint_layer_phase_derivative_above_21000_proved"] is True, "endpoint-layer route drift")
    require(dependencies["six_current_defect"]["decision"]["six_current_Fresnel_tail_with_finite_remainder_proved"] is True, "six-current hierarchy drift")
    require(dependencies["finite_regulator_equivalence"]["decision"]["finite_R_after_A_defined_on_one_common_regulator"] is True, "physical normalization ownership drift")

    pilot = load_json(PILOT)
    require(pilot.get("passed") is True and pilot["decision"]["floating_values_used_as_proof"] is False, "endpoint pilot proof-boundary drift")
    artifact = {
        "kind": STEM,
        "status": "exact_A_face_endpoint_transition_localized_to_eight_cells_with_34_cell_six_current_tail_replacement_certified",
        "passed": True,
        "scope": {
            "height": T,
            "A": A,
            "endpoint_x_region": "x_16<=x<=1/2",
            "transition_modes": [TRANSITION_FIRST_MODE, TRANSITION_LAST_MODE],
            "tail_modes": [TAIL_FIRST_MODE, TAIL_LAST_MODE],
            "tail_normal_floor": "33/2",
        },
        "symbolic_certificate": symbolic_certificate(),
        "numerical_certificate": numerical_certificate(),
        "floating_route_telemetry": {
            "path": relative(PILOT),
            "sha256": file_hash(PILOT),
            "used_as_proof": False,
            "candidate_physical_endpoint_layer": pilot["quadrature"]["candidate_physical_endpoint_layer"],
            "coarse_fine_difference": pilot["quadrature"]["coarse_fine_difference"],
        },
        "decision": {
            "endpoint_defect_split_into_exact_eight_cell_transition_and_34_cell_tail": True,
            "tail_uniform_abs_q_at_least_16_point_5_proved": True,
            "tail_six_current_defect_error_below_5_point_56e_minus_9_proved": True,
            "physical_tail_replacement_error_below_6e_minus_15_proved": True,
            "special_function_obligation_reduced_to_eight_cells": True,
            "eight_cell_transition_interval_value_proved": False,
            "complete_endpoint_layer_bound_proved": False,
            "full_signed_A_face_bound_proved": False,
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
        "next_obligation": "Certify only the eight-cell transition D_(tr,8) on x_16<=x<=1/2. Use the stable saddle-relative q phase from the endpoint pilot, cellwise normal quadrature, and either Arb Taylor/Cauchy panels or centre-plus-ODE derivative enclosures for the exact H(q) current. Reattach the 34-cell rational-log tail with the <6e-15 physical replacement error.",
        "proof_boundary": "Exact eight-cell/34-cell endpoint decomposition, uniform six-current tail replacement below 5.56e-9, and physical replacement cost below 6e-15 only. No interval value for the eight-cell transition, complete endpoint-layer value, full signed A-face bound, quantitative R_after_A, R_Dir, or Q_K-T estimate, all-height theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified eight-cell localization of the A-face endpoint transition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
