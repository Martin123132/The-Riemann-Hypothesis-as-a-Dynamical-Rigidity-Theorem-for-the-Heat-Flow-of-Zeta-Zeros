#!/usr/bin/env python3
"""Extend the first fold-event residual through the selector top corridor."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event0_top_corridor_transport_extension_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
DEPENDENCIES = {
    "continuous_event_atlas": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_continuous_height_gate.json",
    "beta4_completed_projection": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_beta4_completed_branch_projection_gate.json",
    "stationary_carrier_bridge": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selected_logistic_stationary_action_bridge_gate.json",
    "turning_event_handoff": REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json",
}

C = 159_577
MODE = 39_894
THETA_MAX = 2
MOMENT_ORDER = 6
PRECISION = 100


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def add_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def polynomial_value(coefficients: list[acb], theta: int) -> acb:
    total = acb(0)
    power = arb(1)
    for coefficient in coefficients:
        total += power * coefficient
        power *= theta
    return total


def extension_certificate(row: dict[str, Any]) -> dict[str, Any]:
    coefficients = [parse_complex(record["transport_coefficient_ball"]) for record in row["moments"]]
    require(len(coefficients) == MOMENT_ORDER + 1, "stored moment order drift")
    mass = arb(row["finite_contour_absolute_mass_bound"])
    expansion_size = arb(row["maximum_transport_argument"])
    event_far = arb(row["event_far_tail_absolute_bound"])
    z_max = (arb(21) * 21 + 1).sqrt()
    kappa = expansion_size / z_max
    scaled_size = THETA_MAX * expansion_size
    finite_remainder = mass * scaled_size.exp() * scaled_size ** (MOMENT_ORDER + 1) / arb(math.factorial(MOMENT_ORDER + 1))
    transported_far = (THETA_MAX * kappa).exp() * event_far
    polynomial_majorant = sum(
        (arb(THETA_MAX) ** degree * abs(coefficient).upper() for degree, coefficient in enumerate(coefficients)),
        arb(0),
    )
    computed_normalized = polynomial_majorant + finite_remainder + transported_far
    computed_physical = arb(2).sqrt() / arb.pi() * computed_normalized
    uniform_normalized = arb("2.141e-6")
    uniform_physical = arb("9.64e-7")
    endpoint = add_error(polynomial_value(coefficients, THETA_MAX), finite_remainder + transported_far)
    require(computed_normalized < uniform_normalized, "event-0 top-corridor residual exceeds certified cap")
    require(computed_physical < uniform_physical, "event-0 top-corridor physical residual exceeds certified cap")
    return {
        "transport_theta_interval": "0<=theta<=2",
        "top_corridor_theta_interval": "1<=theta<=2",
        "moment_order": MOMENT_ORDER,
        "serialized_coefficient_count": len(coefficients),
        "polynomial_absolute_majorant_ball": polynomial_majorant.str(PRECISION, more=True),
        "maximum_transport_argument_ball": scaled_size.str(PRECISION, more=True),
        "finite_exponential_remainder_bound_ball": finite_remainder.str(PRECISION, more=True),
        "transported_far_tail_bound_ball": transported_far.str(PRECISION, more=True),
        "computed_uniform_normalized_bound_ball": computed_normalized.str(PRECISION, more=True),
        "computed_uniform_physical_bound_ball": computed_physical.str(PRECISION, more=True),
        "uniform_normalized_exact_minus_beta4_bound_ball": uniform_normalized.str(PRECISION, more=True),
        "uniform_physical_exact_minus_beta4_bound_ball": uniform_physical.str(PRECISION, more=True),
        "selector_center_theta2_transport_ball": complex_record(endpoint),
    }


def render_note(artifact: dict[str, Any]) -> str:
    g = artifact["geometry"]
    c = artifact["certificate"]
    return f"""# Event-0 transport through the selector top corridor

Date: 2026-08-13

Status: rigorous one-mode fold/Fresnel transport extension; not a complete
ordinary splice, `T_upper` theorem, or proof of RH

For the first event mode `m={MODE}`,

```text
tau_0=t*-pi/8,        h=(pi/16)theta.                  (TC1)
```

The exact-minus-beta-minus-four integrand obeys the already certified identity

```text
D_(0,tau_0+h)(z)=exp(i*h*z/beta)D_(0,tau_0)(z).        (TC2)
```

The existing atlas used moments through degree `{MOMENT_ORDER}` on
`|theta|<=1`.  Reusing those stored interval balls, multiplying the kth
coefficient by `2^k`, and recomputing the exponential and lifted-tail
remainders proves the larger one-sided interval

```text
0<=theta<=2,
tau_0<=t<=t*,
|Delta I_2(0,t)|
 < {c['uniform_normalized_exact_minus_beta4_bound_ball']} < 2.15e-6,
physical < {c['uniform_physical_exact_minus_beta4_bound_ball']}. (TC3)
```

The top corridor is exactly `1<=theta<=2`, because its lower face is
`{g['top_corridor_lower_height']}` and `theta=2` is the selector center
`{g['selector_center_height']}`.  Thus the edge mode remains fold-owned all
the way from its certified event face to the completed selector object; it
need not be compared there with the invalid bare full-saddle carrier.

This does not identify the one-mode transport ball with the complete 399-mode
selector strip.  At `t*`, the mode is only one summand inside the selected
branch projection, whose opposite branch, outer complement, and endpoint
half-current remain grouped.

Machine-audited companion:

```text
{relative(NOTE)}
{relative(RESULT)}
{relative(BUILDER)}
{relative(CHECKER)}
```

No remaining-398-mode grouped estimate, finite-integral ordinary amplitude
theorem, all-corridor continuation, complete `Q_K-T` or `T_upper`,
`Lambda<=0`, PF-infinity, RH, or prize-level conclusion is proved.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    for path in (*DEPENDENCIES.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    dependencies = {name: json.loads(path.read_text(encoding="utf-8")) for name, path in DEPENDENCIES.items()}
    atlas = dependencies["continuous_event_atlas"]
    require(atlas["claims"]["all_399_exact_event_cells_uniformly_certified"] is True, "continuous atlas dependency drift")
    require(dependencies["beta4_completed_projection"].get("passed") is True, "beta4 completion dependency failed")

    ctx.dps = PRECISION
    row = atlas["certificate"]["all_rows"][0]
    require((row["event_index"], row["mode"]) == (0, MODE), "event-0 row drift")
    certificate = extension_certificate(row)
    pi = arb.pi()
    c = arb(C)
    mode = arb(MODE)
    tstar = pi * c * c / 8
    tau = pi * mode * (c - 2 * mode)
    require((tstar - tau - pi / 8).contains(0), "event-0 selector distance drift")
    top_lower = tstar - pi / 16
    require((tau + pi / 16 - top_lower).contains(0), "event face/top-corridor boundary mismatch")
    require((tau + THETA_MAX * pi / 16 - tstar).contains(0), "theta=2 does not reach selector center")

    artifact = {
        "kind": STEM,
        "date": "2026-08-13",
        "status": "event0_beta4_residual_transported_through_complete_selector_top_corridor",
        "passed": True,
        "geometry": {
            "C": C,
            "mode": MODE,
            "event_index": 0,
            "event_center_identity": "tau_0=t*-pi/8",
            "event_center_height": tau.str(PRECISION, more=True),
            "event_upper_face_height": top_lower.str(PRECISION, more=True),
            "top_corridor_lower_height": top_lower.str(PRECISION, more=True),
            "selector_center_height": tstar.str(PRECISION, more=True),
            "theta_height_law": "t=tau_0+(pi/16)theta",
        },
        "certificate": certificate,
        "decision": {
            "stored_event_moments_reused_without_atlas_recomputation": True,
            "event0_exact_minus_beta4_residual_certified_for_0_le_theta_le_2": True,
            "event0_fold_chart_covers_complete_top_corridor": True,
            "bare_full_saddle_comparison_avoided": True,
            "complete_399_mode_selector_object_identified_with_event0": False,
            "remaining_398_mode_grouped_estimate_proved": False,
            "ordinary_Morse_corridor_splice_proved": False,
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
            "workers": 1,
            "process_priority": priority,
            "elapsed_seconds": round(time.perf_counter() - started, 3),
        },
        "next_action": (
            "At t*, remove the certified m=39894 selected summand from the beta^-4 completed branch projection by exact finite algebra. "
            "Then derive a grouped finite-difference bound for the remaining 199 reflected pairs, retaining the completed remainder."
        ),
        "proof_boundary": (
            "One event-mode exact-minus-beta4 transport on tau_0<=t<=t* only. No identification with the complete selector strip, "
            "remaining-mode grouped bound, all-corridor continuation, complete Q_K-T or T_upper theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion."
        ),
    }
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    NOTE.parent.mkdir(parents=True, exist_ok=True)
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    print("extended event 0 through selector top: theta<=2, normalized residual < 2.15e-6", flush=True)


if __name__ == "__main__":
    main()
