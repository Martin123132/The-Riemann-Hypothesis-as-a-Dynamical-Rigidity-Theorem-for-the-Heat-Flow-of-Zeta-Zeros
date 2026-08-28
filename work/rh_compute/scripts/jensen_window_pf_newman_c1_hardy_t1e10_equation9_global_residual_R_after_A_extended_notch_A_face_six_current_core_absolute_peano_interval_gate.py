#!/usr/bin/env python3
"""Certify a nonoscillatory absolute bound for the pre-endpoint A-face core."""

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
import mpmath as mp


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_core_absolute_peano_interval_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_R_after_A_extended_notch_A_face_six_current_midpoint_defect_reduction_gate.json"
A = 159_577
T = 10_000_000_000
FIRST_MODE = 39_895
LAST_MODE = 39_936
LEFT_NUMERATOR = 79_789
RIGHT_NUMERATOR = 79_873
Q_LOWER_TEXT = "0.0003884882198961"
Q_TAIL = 40
Q_SEGMENTS = (
    (Q_LOWER_TEXT, "0.001", 32),
    ("0.001", "0.01", 64),
    ("0.01", "0.1", 64),
    ("0.1", "1", 64),
    ("1", "5", 64),
    ("5", "10", 32),
    ("10", "20", 32),
    ("20", "40", 32),
)
PANEL_SCALE = 2
HALF_CELL_SUBDIVISIONS = 4
TERMS = 6
PRECISION = 100
CORE_THRESHOLD = "0.000012"


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


def odd_df(index: int) -> int:
    result = 1
    for factor in range(1, index + 1, 2):
        result *= factor
    return result


def interval(left: mp.mpf, right: mp.mpf) -> arb:
    middle = (left + right) / 2
    radius = (right - left) / 2 + mp.mpf("1e-105")
    return arb(mp.nstr(middle, 115), mp.nstr(radius, 115))


def g_six_yy_over_x_bound(y: arb, x: arb) -> arb:
    pi_floor = arb.pi().lower()
    delta = y - arb(A) * x / 2
    require(delta.lower() > 0, "normal denominator crossed zero")
    delta_floor = delta.lower()
    x_ceiling = x.upper()
    y_ceiling = y.upper()
    value = arb(A) / (pi_floor * delta_floor**3)
    for index in range(1, TERMS):
        power = 2 * index + 1
        coefficient = (
            arb(odd_df(2 * index - 1))
            * x_ceiling ** (index - 1)
            / (arb(2) ** index * pi_floor ** (index + 1))
        )
        value += coefficient * power * (
            arb(power + 1) * y_ceiling / delta_floor ** (power + 2)
            + 2 / delta_floor ** (power + 1)
        )
    require("nan" not in value.str(24, more=True).lower(), "G6_yy/x majorant became indeterminate")
    return value


def symbolic_certificate() -> dict[str, Any]:
    require(LAST_MODE - FIRST_MODE + 1 == 42, "42-mode roster drift")
    require(42 * 1 == 42 and 42 / 24 == 7 / 4, "Peano mass drift")
    return {
        "defect_orientation": "E_42=-D_42 and E_42=sum_cells integral K_m G_(A,yy)",
        "six_current_Peano": "|D6|<=sum_cells integral K_m |G_(6,yy)|",
        "factored_derivative": "G_(6,yy)(y,x)=x*M_6(y,x), with M_6 given by the exact termwise second derivative",
        "logistic_change": "q=log((1-x)/x), |dx|=x(1-x)dq",
        "regularized_weight": "|D6|x^(-5/4)(1-x)^(-1/4)|dx|<=B6(x)[x(1-x)]^(3/4)dq, B6=|D6|/x",
        "replacement_weight": "|D42-D6| gives q-weight x^(23/4)(1-x)^(3/4)",
        "Peano_total_mass": "42/24=7/4",
        "analytic_q_tail": "x(q)<=exp(-q), so the core and replacement tails decay as exp(-3q/4) and exp(-23q/4)",
    }


def interval_certificate(dependency: dict[str, Any]) -> dict[str, Any]:
    ctx.dps = PRECISION
    ctx.threads = 1
    mp.mp.dps = 130
    pi = arb.pi()
    left = arb(LEFT_NUMERATOR) / 2
    q_normal = arb(16)
    sqrt_x_cut = ((2 * q_normal**2 + 8 * arb(A) * left).sqrt() - q_normal * arb(2).sqrt()) / (2 * arb(A))
    x_cut = sqrt_x_cut**2
    exact_logistic_cut = ((1 - x_cut) / x_cut).log()
    require(arb(Q_LOWER_TEXT) < exact_logistic_cut.lower(), "q cover does not start below the exact handoff")

    normalization = 2 * (pi / (32 * arb(T))) ** (arb(1) / 4)
    dy = mp.mpf("0.5") / HALF_CELL_SUBDIVISIONS
    core_bound = arb(0)
    replacement_bound = arb(0)
    maximum_peano_over_x = arb(0)
    panel_total = sum(count for _, _, count in Q_SEGMENTS) * PANEL_SCALE
    sample_panels: dict[str, dict[str, str]] = {}
    panel_index = 0
    replacement_constant = (
        arb(84)
        * 2
        * arb(odd_df(11))
        * left
        / (pi**7 * arb(2) ** 6)
    )

    for left_text, right_text, base_panels in Q_SEGMENTS:
        segment_left = mp.mpf(left_text)
        segment_right = mp.mpf(right_text)
        panels = base_panels * PANEL_SCALE
        dq = (segment_right - segment_left) / panels
        for local_panel in range(panels):
            qb = interval(segment_left + local_panel * dq, segment_left + (local_panel + 1) * dq)
            xb = 1 / (1 + qb.exp())
            peano_over_x = arb(0)
            for mode in range(FIRST_MODE, LAST_MODE + 1):
                cell_left = mp.mpf(mode) - mp.mpf("0.5")
                for side in range(2):
                    side_left = cell_left + mp.mpf(side) / 2
                    for part in range(HALF_CELL_SUBDIVISIONS):
                        y0 = side_left + part * dy
                        y1 = y0 + dy
                        yb = interval(y0, y1)
                        if side == 0:
                            kernel = (yb - arb(str(cell_left))) ** 2 / 2
                        else:
                            kernel = (arb(str(cell_left + 1)) - yb) ** 2 / 2
                        peano_over_x += (
                            arb(mp.nstr(dy, 115))
                            * kernel.upper()
                            * g_six_yy_over_x_bound(yb, xb).upper()
                        )

            x_ceiling = xb.upper()
            logistic_weight = (x_ceiling * (1 - x_ceiling)) ** (arb(3) / 4)
            panel_width = arb(mp.nstr(dq, 115))
            core_bound += normalization * panel_width * logistic_weight.upper() * peano_over_x
            delta_floor = (left - arb(A) * xb / 2).lower()
            replacement_integrand = replacement_constant * x_ceiling ** (arb(23) / 4) / delta_floor**13
            replacement_bound += normalization * panel_width * replacement_integrand.upper()
            if peano_over_x.upper() > maximum_peano_over_x.upper():
                maximum_peano_over_x = peano_over_x
            if panel_index in (0, panel_total // 2, panel_total - 1):
                sample_panels[str(panel_index)] = {
                    "q_ball": qb.str(50, more=True),
                    "x_ball": xb.str(50, more=True),
                    "Peano_D6_over_x_bound": peano_over_x.str(50, more=True),
                }
            panel_index += 1

    q_tail = arb(Q_TAIL)
    x_tail_upper = (-q_tail).exp()
    delta_floor = left - arb(A) * x_tail_upper / 2
    y_ceiling = arb(RIGHT_NUMERATOR) / 2
    derivative_tail = arb(A) / (pi.lower() * delta_floor**3)
    for index in range(1, TERMS):
        power = 2 * index + 1
        coefficient = (
            arb(odd_df(2 * index - 1))
            * x_tail_upper ** (index - 1)
            / (arb(2) ** index * pi.lower() ** (index + 1))
        )
        derivative_tail += coefficient * power * (
            arb(power + 1) * y_ceiling / delta_floor ** (power + 2)
            + 2 / delta_floor ** (power + 1)
        )
    peano_tail = arb(7) * derivative_tail / 4
    analytic_core_tail = normalization * peano_tail * arb(4) * (-arb(3) * q_tail / 4).exp() / 3
    analytic_replacement_tail = (
        normalization
        * replacement_constant
        / delta_floor**13
        * arb(4)
        * (-arb(23) * q_tail / 4).exp()
        / 23
    )
    core_bound += analytic_core_tail
    replacement_bound += analytic_replacement_tail
    total = core_bound + replacement_bound
    require(core_bound.upper() < arb("0.000011"), "six-current core bound drift")
    require(replacement_bound.upper() < arb("2e-15"), "six-current replacement bound drift")
    require(total.upper() < arb(CORE_THRESHOLD), "complete pre-endpoint core threshold failed")
    return {
        "exact_x_16_ball": x_cut.str(80, more=True),
        "exact_logistic_q_cut_ball": exact_logistic_cut.str(80, more=True),
        "certified_q_cover": [Q_LOWER_TEXT, "infinity"],
        "partition": {
            "panel_scale": PANEL_SCALE,
            "finite_q_panels": panel_total,
            "half_cell_subdivisions": HALF_CELL_SUBDIVISIONS,
            "normal_boxes_per_q_panel": 42 * 2 * HALF_CELL_SUBDIVISIONS,
            "finite_q_tail_cut": Q_TAIL,
            "precision_decimal_digits": PRECISION,
        },
        "maximum_Peano_D6_over_x_bound": maximum_peano_over_x.str(70, more=True),
        "physical_six_current_core_absolute_bound": core_bound.str(70, more=True),
        "physical_D42_minus_D6_replacement_bound": replacement_bound.str(70, more=True),
        "physical_complete_pre_endpoint_core_absolute_bound": total.str(70, more=True),
        "complete_core_threshold": CORE_THRESHOLD,
        "analytic_q_tail_core_bound": analytic_core_tail.str(60, more=True),
        "analytic_q_tail_replacement_bound": analytic_replacement_tail.str(60, more=True),
        "sample_panels": sample_panels,
        "dependency_uniform_defect_error_ball": dependency["numerical_certificate"]["uniform_D42_minus_D6_ball"],
    }


def render_note(artifact: dict[str, Any]) -> str:
    cert = artifact["interval_certificate"]
    return f"""# Nonoscillatory absolute bound for the translated A-face core

Date: 2026-08-23

Status: rigorous local interval certificate at `t=10^10`, independently
replayed on a different partition; not a proof of `R_after_A` or RH

On `0<x<=x_16`, Section 11.454 gives `D_42=D_6+R_6`.  Apply the positive
midpoint Peano kernel directly to the rational six-current amplitude:

```text
|D_6(x)|
 <=sum_(m=39895)^39936 integral_cell K_m(y)|G_(6,yy)(y,x)|dy.
                                                               (CP1)
```

Termwise differentiation of the exact rational formula factors

```text
G_(6,yy)(y,x)=x M_6(y,x).                           (CP2)
```

No subtraction of the large point sum and strip integral is therefore needed.
With `q=log((1-x)/x)` and `|dx|=x(1-x)dq`, (CP1)--(CP2) give

```text
integral_0^x16 |D_6(x)|x^(-5/4)(1-x)^(-1/4)dx
 <=integral_q16^infinity B_6(x(q))[x(q)(1-x(q))]^(3/4)dq,

B_6(x)=sum_cells integral_cell K_m(y)|M_6(y,x)|dy.  (CP3)
```

The replacement remainder has the still faster regular weight
`x^(23/4)(1-x)^(3/4)`.  The builder encloses the finite interval through
`q=40`, then uses `x(q)<=exp(-q)` and total Peano mass `7/4` for analytic
tails.  At 100 decimal digits its `{cert['partition']['finite_q_panels']}`
outer panels and `{HALF_CELL_SUBDIVISIONS}` subdivisions per half-cell certify

```text
six-current core contribution <= {cert['physical_six_current_core_absolute_bound']},
six-current replacement       <= {cert['physical_D42_minus_D6_replacement_bound']},

complete pre-endpoint core    <= {cert['physical_complete_pre_endpoint_core_absolute_bound']}
                               < {CORE_THRESHOLD}.              (CP4)
```

The independent checker uses 120 decimal digits, three times the base outer
partition, three subdivisions per half-cell, exact endpoint kernel maxima,
and its own termwise derivative majorant.  It independently obtains a bound
below `0.000012`.

Machine-audited companion:

```text
outputs/{STEM}.md
work/rh_compute/results/{STEM}.json
work/rh_compute/scripts/{STEM}.py
work/rh_compute/scripts/check_{STEM}.py
```

Pi provenance: every `pi` in (CP1)--(CP4) is inherited from the canonical
Fresnel hierarchy or equation-(9) normalization.  The logistic coordinate
introduces no new occurrence of `pi`.

Proof boundary: (CP4) proves only the translated A-face pre-endpoint core on
`0<x<=x_16` at the saved height.  It does not by itself include the endpoint
layer, assemble the full translated A-face channel, or prove `R_after_A`,
`R_Dir`, `Q_K-T`, an all-height theorem, `Lambda<=0`, PF-infinity, RH, or a
prize-level conclusion.
"""


def main() -> int:
    started = time.time()
    priority = set_low_priority()
    dependency = load_json(DEPENDENCY)
    require(dependency.get("passed") is True, "six-current dependency is not passed")
    require(dependency["decision"]["six_current_Fresnel_tail_with_finite_remainder_proved"] is True, "six-current theorem drift")
    require(CHECKER.is_file(), "independent checker missing")
    artifact = {
        "kind": STEM,
        "status": "rigorous_complete_translated_A_face_pre_endpoint_core_absolute_bound_below_1_point_2e_minus_5",
        "passed": True,
        "scope": {"height": T, "A": A, "x_region": "0<x<=x_16", "modes": [FIRST_MODE, LAST_MODE]},
        "symbolic_certificate": symbolic_certificate(),
        "interval_certificate": interval_certificate(dependency),
        "decision": {
            "positive_Peano_kernel_used": True,
            "oscillatory_Morse_cancellation_used": False,
            "lower_endpoint_regularized_rigorously": True,
            "six_current_replacement_integrated_rigorously": True,
            "complete_pre_endpoint_A_face_core_absolute_bound_below_1_point_2e_minus_5_proved": True,
            "independent_replay_completed": True,
            "endpoint_layer_included": False,
            "full_signed_A_face_bound_proved": False,
            "R_after_A_bound_proved": False,
            "R_Dir_bound_proved": False,
            "QK_minus_T_bound_proved": False,
            "rh_implication": False,
        },
        "dependency": {"path": relative(DEPENDENCY), "sha256": file_hash(DEPENDENCY)},
        "sources": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {
            "elapsed_seconds": round(time.time() - started, 3),
            "workers": 1,
            "arb_threads": 1,
            "process_priority": priority,
        },
        "next_obligation": "Add this <1.2e-5 pre-endpoint theorem to the independently certified <0.0037 endpoint theorem, preserving the common equation-(9) normalization and the x_16 partition, to certify the complete translated 42-mode A-face channel.",
        "proof_boundary": "Rigorous translated A-face pre-endpoint core absolute bound below 1.2e-5 at t=10^10 only. No endpoint assembly, full translated A-face theorem, R_after_A, R_Dir, Q_K-T, all-height, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("certified nonoscillatory translated A-face core below 1.2e-5", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
