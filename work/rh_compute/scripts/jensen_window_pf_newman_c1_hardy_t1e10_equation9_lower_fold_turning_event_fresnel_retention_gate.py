#!/usr/bin/env python3
"""Certify characteristic crowding and mandatory Fresnel retention at turning events."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, acb, ctx


EVENT_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate.json"
FRESNEL_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fresnel_endpoint_normal_form_gate.json"
MORSE_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_universal_logistic_morse_characteristic_fold_reduction_gate.json"
SELECTOR_CELL_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_nonzero_height_completed_strip_cell_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
C = 159_577
Q = 39_894
EVENT_COUNT = 399


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


def mode_for_event(index: int) -> int:
    require(0 <= index < EVENT_COUNT, "event index outside atlas")
    return Q - index // 2 if index % 2 == 0 else Q + (index + 1) // 2


def triangular(index: int) -> int:
    return index * (index + 1) // 2


def event_height(index: int, pi: arb) -> arb:
    return pi * C**2 / 8 - pi * (triangular(index) + arb(1) / 8)


def fresnel_primitive(z: arb, pi: arb) -> acb:
    i = acb(0, 1)
    return (i * pi / 4).exp() / arb(2).sqrt() * ((-i * pi / 4).exp() * (pi / 2).sqrt() * z).erf()


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def endpoint_coordinate(height: arb, mode: int, tau: arb, pi: arb) -> arb:
    """Exact z_C^* from (11.310.4), simplified using the event equation."""
    return (tau - height) / (pi * (height + 2 * pi * mode**2)).sqrt()


def isolation_radius(index: int, pi: arb, selector_lower: arb) -> arb:
    """Largest symmetric radius before another event or the booked top join is met."""
    if index == 0:
        return pi / 16
    upward = pi * index
    if index + 1 < EVENT_COUNT:
        downward = pi * (index + 1)
    else:
        downward = event_height(index, pi) - selector_lower
    return min(upward, downward)


def disjoint_radius(index: int, pi: arb) -> arb:
    """A concrete pairwise-disjoint symmetric event-buffer atlas."""
    return pi / 16 if index == 0 else pi * index / 2


def scan_atlas(pi: arb, selector_lower: arb) -> dict[str, Any]:
    isolation_max = arb(0)
    isolation_record: dict[str, Any] | None = None
    disjoint_max = arb(0)
    disjoint_record: dict[str, Any] | None = None
    minimum_characteristic = arb(1)
    maximum_characteristic = arb(0)

    for index in range(EVENT_COUNT):
        mode = mode_for_event(index)
        tau = event_height(index, pi)
        signed_index = 4 * mode - C
        require(signed_index == (-1 if index % 2 == 0 else 1) * (2 * index + 1), "signed event index drift")

        delta = abs(arb(signed_index) / C)
        minimum_characteristic = min(minimum_characteristic, delta)
        maximum_characteristic = max(maximum_characteristic, delta)

        for family, radius in (
            ("isolation", isolation_radius(index, pi, selector_lower)),
            ("disjoint", disjoint_radius(index, pi)),
        ):
            lower_height = tau - radius
            upper_height = tau + radius
            require(
                lower_height > selector_lower or (lower_height - selector_lower).contains(0),
                f"{family} buffer left selector cell at event {index}",
            )
            if index == 0:
                top_join = pi * C**2 / 8 - pi / 16
                require(
                    upper_height < top_join or (upper_height - top_join).contains(0),
                    "top event buffer crossed booked selector join",
                )
            lower_z = abs(endpoint_coordinate(lower_height, mode, tau, pi))
            upper_z = abs(endpoint_coordinate(upper_height, mode, tau, pi))
            local_max = max(lower_z, upper_z)
            record = {
                "event_index": index,
                "mode": mode,
                "radius_ball": radius.str(PRECISION, more=True),
                "lower_face_abs_z_ball": lower_z.str(PRECISION, more=True),
                "upper_face_abs_z_ball": upper_z.str(PRECISION, more=True),
            }
            if family == "isolation" and local_max > isolation_max:
                isolation_max = local_max
                isolation_record = record
            if family == "disjoint" and local_max > disjoint_max:
                disjoint_max = local_max
                disjoint_record = record

    require(isolation_record is not None and disjoint_record is not None, "atlas scan produced no records")
    require(isolation_max < arb("0.004976"), "isolation envelope exceeds 0.004976 Fresnel units")
    require(disjoint_max < arb("0.002501"), "disjoint envelope exceeds 0.002501 Fresnel units")
    require((minimum_characteristic - arb(1) / C).contains(0), "nearest characteristic defect drift")
    require((maximum_characteristic - arb(797) / C).contains(0), "outer characteristic defect drift")
    return {
        "event_count": EVENT_COUNT,
        "signed_event_index_law": "4*m_n-C=(-1)^(n+1)*(2*n+1)",
        "minimum_abs_logistic_characteristic_ball": minimum_characteristic.str(PRECISION, more=True),
        "maximum_abs_logistic_characteristic_ball": maximum_characteristic.str(PRECISION, more=True),
        "maximum_reciprocal_characteristic_loss": C,
        "isolation_envelope": {
            "definition": "largest symmetric radius before another event, the selector-cell lower edge, or the booked top join",
            "maximum_abs_Fresnel_coordinate_ball": isolation_max.str(PRECISION, more=True),
            "maximizer": isolation_record,
        },
        "pairwise_disjoint_symmetric_buffers": {
            "radius_law": "r_0=pi/16; r_n=pi*n/2 for 1<=n<=398",
            "maximum_abs_Fresnel_coordinate_ball": disjoint_max.str(PRECISION, more=True),
            "maximizer": disjoint_record,
        },
        "full_interior_Fresnel_relative_defect_lower_bound_ball": (
            arb(1) / 2 - isolation_max / arb(2).sqrt()
        ).str(PRECISION, more=True),
        "nonstationary_zero_replacement_absolute_defect_lower_bound_ball": (
            1 / arb(2).sqrt() - isolation_max
        ).str(PRECISION, more=True),
    }


def prototype(pi: arb, beta: arb, selector_top: arb) -> dict[str, Any]:
    mode = Q
    tau = event_height(0, pi)
    radius = pi / 16
    lower_height = tau - radius
    upper_height = tau + radius
    d = 4 * beta * (arb(mode) - arb(C) / 4) / C
    lambda_event = (selector_top - tau) / beta
    require((lambda_event - d**2).contains(0), "prototype threshold lambda=d^2 failed")

    logistic_delta = (2 * pi * mode**2 - tau) / (2 * pi * mode**2 + tau)
    require((logistic_delta - arb(4 * mode - C) / C).contains(0), "prototype logistic characteristic identity failed")
    require((logistic_delta - d / beta).contains(0), "prototype d/beta identity failed")

    lower_lambda = (selector_top - lower_height) / beta
    upper_lambda = (selector_top - upper_height) / beta
    lower_defect = lower_lambda - d**2
    upper_defect = upper_lambda - d**2
    require((lower_defect - radius / beta).contains(0), "lower defect orientation failed")
    require((upper_defect + radius / beta).contains(0), "upper defect orientation failed")

    lower_z = endpoint_coordinate(lower_height, mode, tau, pi)
    upper_z = endpoint_coordinate(upper_height, mode, tau, pi)
    require(lower_z > 0 and upper_z < 0, "prototype Fresnel orientation failed")
    maximum_abs_z = max(abs(lower_z), abs(upper_z))

    half = acb(arb(1) / 2, arb(1) / 2)
    full = acb(1, 1)
    lower_transition = half - fresnel_primitive(lower_z, pi)
    upper_transition = half - fresnel_primitive(upper_z, pi)
    lower_full_defect = abs(lower_transition - full) / arb(2).sqrt()
    upper_full_defect = abs(upper_transition - full) / arb(2).sqrt()
    analytic_relative_lower = arb(1) / 2 - maximum_abs_z / arb(2).sqrt()
    require(analytic_relative_lower > arb("0.4999994"), "prototype full-Morse defect lost 0.4999994 lower bound")
    require(lower_full_defect > analytic_relative_lower and upper_full_defect > analytic_relative_lower, "direct Fresnel values violate analytic bound")

    return {
        "event_index": 0,
        "mode": mode,
        "crossing_height": "tau=t*-pi/8",
        "crossing_height_ball": tau.str(PRECISION, more=True),
        "buffer_radius": "pi/16",
        "buffer_lower_height_ball": lower_height.str(PRECISION, more=True),
        "buffer_upper_height_ball": upper_height.str(PRECISION, more=True),
        "detuning_ball": d.str(PRECISION, more=True),
        "event_lambda_ball": lambda_event.str(PRECISION, more=True),
        "event_threshold_identity": "lambda(tau)=d^2=pi/(8*beta)",
        "signed_fold_defect_identity": "lambda(t)-d^2=(tau-t)/beta",
        "lower_face_fold_defect_ball": lower_defect.str(PRECISION, more=True),
        "upper_face_fold_defect_ball": upper_defect.str(PRECISION, more=True),
        "lower_face_Fresnel_coordinate_ball": lower_z.str(PRECISION, more=True),
        "upper_face_Fresnel_coordinate_ball": upper_z.str(PRECISION, more=True),
        "maximum_face_abs_Fresnel_coordinate_ball": maximum_abs_z.str(PRECISION, more=True),
        "lower_face_transition_factor": complex_record(lower_transition),
        "upper_face_transition_factor": complex_record(upper_transition),
        "lower_face_full_interior_relative_defect_ball": lower_full_defect.str(PRECISION, more=True),
        "upper_face_full_interior_relative_defect_ball": upper_full_defect.str(PRECISION, more=True),
        "analytic_full_interior_relative_defect_lower_bound_ball": analytic_relative_lower.str(PRECISION, more=True),
        "logistic_characteristic_at_event_ball": logistic_delta.str(PRECISION, more=True),
        "logistic_characteristic_identity": "Delta_m(tau)=(4m-C)/C=d/beta=-1/C",
    }


def render_note(artifact: dict[str, Any]) -> str:
    p = artifact["prototype_event"]
    a = artifact["atlas_certificate"]
    return f"""# Turning-event characteristic crowding and Fresnel retention

Date: 2026-08-12

Status: exact event-coordinate conversion and a rigorous obstruction to a
full-interior Fresnel replacement are certified; this is not a proof of the
replacement Airy-to-Fresnel-retaining logistic join

For fixed `C=159577`, put

```text
t*=pi*C^2/8,                 beta^3=t*,
d_m=4*beta*(m-C/4)/C,       lambda(t)=(t*-t)/beta.
```

At the event `tau_m=pi*m(C-2m)`, exact algebra gives

```text
lambda(tau_m)=d_m^2,
lambda(t)-d_m^2=(tau_m-t)/beta,                         (FR1)
Delta_m(tau_m)=(4m-C)/C=d_m/beta.                       (FR2)
```

The reduced-saddle Fresnel endpoint coordinate from Section 11.310 is

```text
z_C^*(m,t)=(tau_m-t)/sqrt(pi*(t+2*pi*m^2)).              (FR3)
```

Thus the Airy fold defect and the logistic/Fresnel transition have the same
sign and vanish at exactly the same height.  No characteristic denominator
is introduced.

## Nearest event

For mode `{p['mode']}`,

```text
tau=t*-pi/8,
d={p['detuning_ball']},
d^2=pi/(8*beta).
```

On the booked exact-mode buffer `|t-tau|<=pi/16`, the two face coordinates
are

```text
z_lower={p['lower_face_Fresnel_coordinate_ball']},
z_upper={p['upper_face_Fresnel_coordinate_ball']}.
```

Let `F'(z)=exp(i*pi*z^2/2)`, `F(0)=0`, and
`Q(z)=(1+i)/2-F(z)`.  For real `z`, `|F(z)|<=|z|`.  Replacing `Q(z)` by the
full-interior value `1+i` therefore has relative defect at least

```text
1/2-|z|/sqrt(2)
 >= {p['analytic_full_interior_relative_defect_lower_bound_ball']}.          (FR4)
```

So the prototype face is still a half-saddle to better than one part in a
million.  It is not an ordinary full Gaussian face.

## All 399 events

The ordered event modes obey

```text
4*m_n-C=(-1)^(n+1)*(2*n+1).
```

Consequently every atlas event remains in the characteristic band

```text
{a['minimum_abs_logistic_characteristic_ball']}
 <= |Delta_m(tau_m)| <=
{a['maximum_abs_logistic_characteristic_ball']}.
```

Even the largest symmetric event-isolation radius anywhere in the selector
cell reaches at most

```text
|z_C^*| <= {a['isolation_envelope']['maximum_abs_Fresnel_coordinate_ball']}.
```

Hence a full-interior Fresnel replacement has atlas-wide relative defect at
least

```text
{a['full_interior_Fresnel_relative_defect_lower_bound_ball']}.               (FR5)
```

This closes the proposed direct `exact buffer -> full ordinary-Morse main`
route as inadmissible.  The next valid chart must keep the exact endpoint
current and Fresnel transition inside the universal logistic Morse phase,
then compare that grouped object with the Airy-Fresnel representation before
taking absolute values.

## Phase orientation

Because `C` is odd,

```text
(-1)^m exp(-i*pi*m*C)=1
```

for every integer mode.  Since `C=1 mod 8`, the source-rotated common lower
boundary carrier is also exact:

```text
exp(i[t*-pi/8])=1.
```

No phase was discarded in obtaining (FR1)--(FR5).

## Boundary

This gate proves exact coordinate identities and rules out one nonuniform
approximation strategy.  It does not prove the Airy-to-Fresnel-retaining
logistic join, propagation through the 399 events, complete `T_upper`,
`Lambda<=0`, RH, or a prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    dependencies = {
        "event_atlas_gate": EVENT_GATE,
        "Fresnel_normal_form_gate": FRESNEL_GATE,
        "logistic_Morse_gate": MORSE_GATE,
        "selector_cell_gate": SELECTOR_CELL_GATE,
    }
    for path in (*dependencies.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    for path in dependencies.values():
        require(json.loads(path.read_text(encoding="utf-8")).get("passed") is True, f"dependency is not passed: {path}")

    ctx.dps = PRECISION
    pi = arb.pi()
    selector_top = pi * C**2 / 8
    selector_lower = pi * (C - 2) ** 2 / 8
    beta = selector_top ** (arb(1) / 3)
    require(C == 4 * Q + 1 and C % 8 == 1, "selector residue drift")
    require(((C * C - 1) // 8) % 2 == 0, "source-rotated common carrier is not one")

    prototype_record = prototype(pi, beta, selector_top)
    atlas_record = scan_atlas(pi, selector_lower)
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate",
        "status": "turning_event_characteristic_crowding_and_Fresnel_retention_certified",
        "passed": True,
        "fixed_selector": {
            "C": C,
            "q": Q,
            "selector_top_ball": selector_top.str(PRECISION, more=True),
            "selector_lower_ball": selector_lower.str(PRECISION, more=True),
            "beta_ball": beta.str(PRECISION, more=True),
            "common_mode_phase_identity": "(-1)^m*exp(-i*pi*m*C)=1 for odd C",
            "source_rotated_boundary_carrier": "exp(i[t*-pi/8])=1 because C=1 mod 8",
        },
        "exact_coordinate_identities": {
            "turning_height": "tau_m=pi*m*(C-2m)",
            "Airy_fold_defect": "lambda(t)-d_m^2=(tau_m-t)/beta",
            "Fresnel_endpoint_coordinate": "z_C^*(m,t)=(tau_m-t)/sqrt(pi*(t+2*pi*m^2))",
            "coordinate_conversion": "z_C^*=beta*(lambda-d_m^2)/sqrt(pi*(t+2*pi*m^2))",
            "event_logistic_characteristic": "Delta_m(tau_m)=(4m-C)/C=d_m/beta",
        },
        "prototype_event": prototype_record,
        "atlas_certificate": atlas_record,
        "decision": {
            "exact_Airy_event_Fresnel_coordinate_join_proved": True,
            "all_399_events_remain_characteristically_crowded": True,
            "naive_full_interior_Fresnel_face_replacement_admissible": False,
            "exact_Fresnel_retention_required": True,
            "Airy_to_Fresnel_retaining_logistic_join_proved": False,
            "propagation_across_399_events_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Construct the one-mode prototype in the exact logistic Morse coordinate while retaining "
            "F(z_B)-F(z_C) and the paired endpoint current; compare its grouped integral with the "
            "Airy-Fresnel transform on the overlap before absolute values. Do not substitute the full "
            "interior Fresnel limit at an event-buffer face."
        ),
        "proof_boundary": (
            "Exact event-coordinate algebra, rigorous 399-event crowding bounds, and a no-go result for "
            "the full-interior Fresnel face replacement only. No Airy-to-Fresnel-retaining logistic join, "
            "399-event propagation theorem, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
        "dependencies": {
            name: {"path": relative(path), "sha256": file_hash(path)} for name, path in dependencies.items()
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
    }
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    print(
        "certified turning-event Fresnel retention: "
        f"events={EVENT_COUNT}, isolation_z<0.004976, full-Morse-defect>0.496; priority={priority}"
    )


if __name__ == "__main__":
    main()
