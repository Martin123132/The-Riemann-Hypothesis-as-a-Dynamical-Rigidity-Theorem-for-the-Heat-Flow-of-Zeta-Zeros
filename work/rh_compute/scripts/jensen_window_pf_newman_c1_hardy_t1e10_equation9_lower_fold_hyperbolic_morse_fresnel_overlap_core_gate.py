#!/usr/bin/env python3
"""Certify an exact hyperbolic fold form and a compact Airy overlap core."""

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

from flint import arb, ctx
import sympy as sp


RETENTION_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_fresnel_retention_gate.json"
LOGISTIC_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_fresnel_retaining_logistic_transform_gate.json"
AIRY_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_grouped_characteristic_airy_boundary_normal_form_gate.json"
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 100
C = 159_577
B = 5_122_423
MODE = 39_894
Y_CUTOFF = 64
DELTA_RADIUS_TEXT = "0.00284"
CORE_U = 200
PANELS = 65_536
SERIES_DEGREE = 18


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


def interval_ball(left: arb, right: arb) -> arb:
    return arb((left + right) / 2, (right - left) / 2)


def upper_abs(value: arb) -> arb:
    return abs(value).upper()


def symbolic_certificate() -> dict[str, str]:
    x, m, t, v, q, Q, Y, beta, Csym = sp.symbols("x m t v q Q Y beta C", positive=True)
    pi = sp.pi
    r = t / (2 * pi * m**2)
    x_v = 1 / (1 + r * v)

    completion = -pi * m**2 / x + pi * q**2 / 2
    endpoint_phase = pi * x * Csym**2 / 4 - pi * m * Csym
    q_c = sp.sqrt(x / 2) * (Csym - 2 * m / x)
    require(sp.simplify(completion.subs(q, q_c) - endpoint_phase) == 0, "quadratic completion failed")

    psi = -pi * m**2 / x_v + t * sp.log(r * v) / 2
    psi_one = sp.simplify(psi.subs(v, 1))
    require(
        sp.simplify(sp.expand_log(psi - psi_one + t * (v - 1 - sp.log(v)) / 2, force=True)) == 0,
        "logistic phase identity failed",
    )

    rho = sp.sqrt(2 * x_v)
    sigma = 4 * beta / (pi * Csym)
    y = sp.sqrt(2 * beta) * Y
    q_increment = sp.sqrt(pi * x_v / 2) * sigma * y
    require(
        sp.simplify(q_increment.subs(beta**3, pi * Csym**2 / 8) - rho * Y) == 0,
        "same-alpha normal scaling failed",
    )

    g = sp.symbols("g", positive=True)
    raw_amplitude = (
        r * x_v**2 * g * sp.sqrt(2 / t)
        * (x_v * (1 - x_v)) ** (-sp.Rational(1, 4))
        * (m * sp.sqrt(2 / pi) * x_v ** (-sp.Rational(3, 2)) + Q / (pi * x_v))
    )
    expected = sp.sqrt(2) / pi * r ** sp.Rational(1, 4) * v ** (-sp.Rational(1, 4)) * g * (
        1 + Q * sp.sqrt(r * x_v / t)
    )
    require(sp.simplify(sp.powsimp(raw_amplitude / expected, force=True) - 1) == 0, "hyperbolic amplitude failed")

    return {
        "scaled_variables": "u=sqrt(t/2)*s, Q=sqrt(pi)*q, Y=y/sqrt(2*beta)",
        "exact_phase": "Phi-Phi_center=(Q^2-Q_0^2-u^2)/2",
        "same_alpha_inner_coordinate": "Q=Q_C(u)+sqrt(2*x(u))*Y",
        "exact_normalized_amplitude": (
            "A=r^(1/4)*v^(-1/4)*(dv/ds)*sqrt(2*x)"
            "*[1+(Q_C+sqrt(2*x)Y)*sqrt(r*x/t)]"
        ),
        "canonical_phase": (
            "Phi_0-Phi_0_center=ell*u+u^2/(2*C)+u^3/(3*C*sqrt(pi))-u*Y+Y^2/2"
        ),
        "canonical_coefficients": "ell=2*(t-tau)/(C*sqrt(pi)), beta^3=pi*C^2/8, d=-beta/C",
        "phase_difference": (
            "R=[(Q_C^2-Q_0^2-u^2)/2-ell*u-u^2/(2*C)-u^3/(3*C*sqrt(pi))]"
            "+[sqrt(2*x)*Q_C+u]Y+(x-1/2)Y^2"
        ),
    }


def k_series(delta: arb) -> arb:
    """Rigorous analytic evaluation of (delta-log(1+delta))/delta^2."""
    value = arb(0)
    for power in range(SERIES_DEGREE, -1, -1):
        coefficient = arb(1 if power % 2 == 0 else -1) / (power + 2)
        value = value * delta + coefficient
    radius = arb(DELTA_RADIUS_TEXT)
    tail = radius ** (SERIES_DEGREE + 1) / ((SERIES_DEGREE + 3) * (1 - radius))
    value += arb(0, tail)
    require(value.lower() > 0, "normalized logarithmic remainder lost positivity")
    return value


def core_row(t: arb, delta: arb, p: dict[str, arb]) -> dict[str, arb]:
    pi = p["pi"]
    m = arb(MODE)
    r = t / (2 * pi * m**2)
    d0 = 1 + r
    denominator = d0 + r * delta
    x = 1 / denominator
    k = k_series(delta)
    u = delta * (t * k).sqrt()

    a0 = (p["tau"] - t) / (pi * m)
    b = t / (pi * m)
    q_c = (a0 - b * delta) * (pi / (2 * denominator)).sqrt()

    # This is (Q_C^2-Q_0^2-u^2)/2 after exact algebraic cancellation.
    inside = (
        -pi * d0 * a0 * b
        - pi * a0**2 * r / 2
        + t * d0 * delta * (r - k * denominator)
    )
    boundary_phase = delta * inside / (2 * d0 * denominator)

    ell = 2 * (t - p["tau"]) / (arb(C) * pi.sqrt())
    boundary_model = ell * u + u**2 / (2 * C) + u**3 / (3 * C * pi.sqrt())
    boundary_remainder = boundary_phase - boundary_model
    rho = (2 * x).sqrt()
    # Rationalize sqrt(2*x)*Q_C+u; its two direct terms are about 200 on the
    # core while their characteristic defect is below 0.096.
    cross_remainder = (
        pi.sqrt() * a0 / denominator
        + delta * t * (k - 2 * r / denominator**2)
        / ((t * k).sqrt() + b * pi.sqrt() / denominator)
    )
    quadratic_remainder = x - arb(1) / 2

    g = (1 + delta) * (2 * k).sqrt()
    base_amplitude = r ** (arb(1) / 4) * (1 + delta) ** (-arb(1) / 4) * g * rho
    epsilon = (r * x / t).sqrt()
    amplitude_zero = base_amplitude * (1 + q_c * epsilon)
    amplitude_top = base_amplitude * (1 + (q_c + rho * p["Ymax"]) * epsilon)
    return {
        "u": u,
        "x": x,
        "q_c": q_c,
        "boundary_remainder": boundary_remainder,
        "cross_remainder": cross_remainder,
        "quadratic_remainder": quadratic_remainder,
        "amplitude_zero": amplitude_zero,
        "amplitude_top": amplitude_top,
    }


def phase_polynomial_bound(row: dict[str, arb], ymax: arb) -> arb:
    r0 = row["boundary_remainder"]
    cross = row["cross_remainder"]
    quad = row["quadratic_remainder"]
    endpoint_zero = upper_abs(r0)
    endpoint_top = upper_abs(r0 + cross * ymax + quad * ymax**2)
    bound = max(endpoint_zero, endpoint_top)

    if not quad.contains(0):
        critical = -cross / (2 * quad)
        if critical.lower() >= 0 and critical.upper() <= ymax.lower():
            bound = max(bound, upper_abs(r0 - cross**2 / (4 * quad)))
        elif not (critical.upper() < 0 or critical.lower() > ymax.upper()):
            bound = max(bound, upper_abs(r0) + ymax * upper_abs(cross) + ymax**2 * upper_abs(quad))
    elif cross.contains(0):
        bound = max(bound, upper_abs(r0) + ymax * upper_abs(cross) + ymax**2 * upper_abs(quad))
    return bound


def face_certificate(label: str, t: arb, p: dict[str, arb]) -> dict[str, Any]:
    radius = arb(DELTA_RADIUS_TEXT)
    maxima = {
        "boundary_phase_remainder": arb(0),
        "linear_inner_remainder": arb(0),
        "quadratic_inner_remainder": arb(0),
        "complete_phase_remainder": arb(0),
        "normalized_amplitude_remainder": arb(0),
    }
    witnesses = {key: 0 for key in maxima}
    for index in range(PANELS):
        left = -radius + 2 * radius * index / PANELS
        right = -radius + 2 * radius * (index + 1) / PANELS
        row = core_row(t, interval_ball(left, right), p)
        values = {
            "boundary_phase_remainder": upper_abs(row["boundary_remainder"]),
            "linear_inner_remainder": upper_abs(row["cross_remainder"]),
            "quadratic_inner_remainder": upper_abs(row["quadratic_remainder"]),
            "complete_phase_remainder": phase_polynomial_bound(row, p["Ymax"]),
            "normalized_amplitude_remainder": max(
                upper_abs(row["amplitude_zero"] - 1),
                upper_abs(row["amplitude_top"] - 1),
            ),
        }
        for key, value in values.items():
            if value > maxima[key]:
                maxima[key] = value
                witnesses[key] = index

    negative = core_row(t, -radius, p)["u"]
    positive = core_row(t, radius, p)["u"]
    require(negative.upper() < -CORE_U and positive.lower() > CORE_U, f"{label} delta cover misses |u|<=200")
    require(maxima["boundary_phase_remainder"] < arb("0.016"), f"{label} boundary phase target failed")
    require(maxima["linear_inner_remainder"] < arb("0.096"), f"{label} inner linear target failed")
    require(maxima["quadratic_inner_remainder"] < arb("0.00073"), f"{label} inner quadratic target failed")
    require(maxima["complete_phase_remainder"] < arb("0.082"), f"{label} complete phase target failed")
    require(maxima["normalized_amplitude_remainder"] < arb("0.001"), f"{label} amplitude target failed")
    return {
        "label": label,
        "height_ball": t.str(PRECISION, more=True),
        "covered_negative_u_ball": negative.str(PRECISION, more=True),
        "covered_positive_u_ball": positive.str(PRECISION, more=True),
        "panel_count": PANELS,
        "bounds": {key: value.str(PRECISION, more=True) for key, value in maxima.items()},
        "witness_panel_indices": witnesses,
    }


def carrier_mismatch(t: arb, p: dict[str, arb]) -> arb:
    r = t / (2 * p["pi"] * MODE**2)
    z0 = p["beta"] * r.log() / 2
    eta = p["tstar"] / t
    exact = t * (z0 / p["beta"] - eta * (z0 / p["beta"]).tanh())
    lam = (p["tstar"] - t) / p["beta"]
    canonical = lam * p["d"] - p["d"] ** 3 / 3
    return exact - canonical


def parameters() -> dict[str, arb]:
    pi = arb.pi()
    tstar = pi * C**2 / 8
    beta = tstar ** (arb(1) / 3)
    tau = pi * MODE * (C - 2 * MODE)
    return {
        "pi": pi,
        "tstar": tstar,
        "beta": beta,
        "tau": tau,
        "d": -beta / C,
        "Ymax": arb(Y_CUTOFF) / (2 * beta).sqrt(),
    }


def render_note(artifact: dict[str, Any]) -> str:
    core = artifact["compact_overlap_core"]
    return f"""# Characteristic hyperbolic Morse-Fresnel overlap core

Date: 2026-08-12

Status: exact hyperbolic fold form and a compact three-height overlap core
certified; this is not a proof of the complete Airy/logistic integral join

Keep the paired finite-Poisson current from Section 11.337 and set

```text
u=sqrt(t/2)s,       Q=sqrt(pi)q,
Y=y/sqrt(2 beta),  beta^3=pi C^2/8.                    (HF1)
```

The universal logistic phase and exact Fresnel square then combine to the
exact hyperbolic phase

```text
Phi-Phi_center=(Q^2-Q_0^2-u^2)/2.                     (HF2)
```

There is no finite-`t` phase remainder in (HF2).  On the same physical alpha
slice as the Airy chart,

```text
Q=Q_C(u)+sqrt(2x(u))Y.                                (HF3)
```

After extracting the common `sqrt(2)/pi` factor, the exact amplitude is

```text
A_t(u,Y)=r^(1/4)v^(-1/4)(dv/ds)sqrt(2x)
 [1+(Q_C+sqrt(2x)Y)sqrt(rx/t)].                       (HF4)
```

For mode `39894`, the canonical Airy-Fresnel phase under
`z=-d+u/sqrt(2 beta)` is

```text
Phi_0-Phi_0,center
 =ell u+u^2/(2C)+u^3/(3C sqrt(pi))-uY+Y^2/2,
ell=2(t-tau)/(C sqrt(pi)).                             (HF5)
```

Thus the exact-minus-canonical phase is

```text
R_t(u,Y)=R_boundary(u)
 +[sqrt(2x)Q_C+u]Y+[x-1/2]Y^2.                        (HF6)
```

The exact mode carrier and (HF5) carrier differ by less than
`{core['maximum_carrier_phase_mismatch_ball']}` at all three certified
heights.

## Certified core

The interval cover `|v-1|<=0.00284` contains `|u|<=200` at

```text
t=tau-pi/16,       t=tau,       t=tau+pi/16,
0<=Y<=64/sqrt(2 beta).
```

Across all three heights, independent interval panels prove

```text
|R_boundary| < {core['uniform_bounds']['boundary_phase_remainder']},
|sqrt(2x)Q_C+u| < {core['uniform_bounds']['linear_inner_remainder']},
|x-1/2| < {core['uniform_bounds']['quadratic_inner_remainder']},
|R_t(u,Y)| < {core['uniform_bounds']['complete_phase_remainder']},
|A_t(u,Y)-1| < {core['uniform_bounds']['normalized_amplitude_remainder']}.   (HF7)
```

The logarithmic quotient
`[v-1-log(v)]/(v-1)^2` was evaluated by a degree-{SERIES_DEGREE} convergent
series with an explicit geometric tail.  This avoids interval cancellation
at `v=1`.

## Interpretation

The Airy fold and ordinary-Morse region are not two unrelated phases.  They
are the same exact quadratic form `(Q^2-u^2)/2` cut by a lower boundary whose
tangent is nearly the characteristic line `Q=-u`.  Straightening that
boundary produces the cubic Airy term in (HF5).  The retained Fresnel current
is therefore the correct transition object on both sides of the event.

## Boundary

This gate proves an exact coordinate identity and compact pointwise phase and
amplitude bounds.  It does not integrate the exact-minus-canonical core,
control `|u|>200`, complete the one-mode Airy/logistic remainder, propagate
through 399 events, establish complete `T_upper`, or prove `Lambda<=0`, RH,
or a prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    dependencies = {
        "Fresnel_retention_gate": RETENTION_GATE,
        "exact_logistic_transform_gate": LOGISTIC_GATE,
        "Airy_normal_form_gate": AIRY_GATE,
    }
    for path in (*dependencies.values(), CHECKER):
        require(path.is_file(), f"missing dependency: {path}")
    for path in dependencies.values():
        require(json.loads(path.read_text(encoding="utf-8")).get("passed") is True, f"dependency is not passed: {path}")

    ctx.dps = PRECISION
    p = parameters()
    require((p["tau"] - (p["tstar"] - p["pi"] / 8)).contains(0), "prototype event height drift")
    require((p["d"] ** 2 - p["pi"] / (8 * p["beta"])).contains(0), "prototype detuning threshold drift")
    heights = {
        "lower_face": p["tau"] - p["pi"] / 16,
        "event": p["tau"],
        "upper_face": p["tau"] + p["pi"] / 16,
    }
    faces = [face_certificate(label, height, p) for label, height in heights.items()]
    mismatches = {label: carrier_mismatch(height, p) for label, height in heights.items()}
    maximum_carrier = max(upper_abs(value) for value in mismatches.values())
    require(maximum_carrier < arb("2e-12"), "carrier mismatch exceeds 2e-12")

    keys = faces[0]["bounds"].keys()
    uniform = {
        key: max(arb(face["bounds"][key]) for face in faces).str(PRECISION, more=True)
        for key in keys
    }
    artifact = {
        "kind": "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_overlap_core_gate",
        "status": "exact_hyperbolic_fold_form_and_three_height_compact_overlap_core_certified",
        "passed": True,
        "symbolic_certificate": symbolic_certificate(),
        "prototype": {
            "C": C,
            "B": B,
            "mode": MODE,
            "event_height_ball": p["tau"].str(PRECISION, more=True),
            "beta_ball": p["beta"].str(PRECISION, more=True),
            "detuning_ball": p["d"].str(PRECISION, more=True),
            "Y64_ball": p["Ymax"].str(PRECISION, more=True),
        },
        "compact_overlap_core": {
            "claimed_u_core": f"|u|<={CORE_U}",
            "covering_v_minus_1_radius": DELTA_RADIUS_TEXT,
            "Y_interval": "0<=Y<=64/sqrt(2*beta)",
            "faces": faces,
            "carrier_phase_mismatches": {
                label: value.str(PRECISION, more=True) for label, value in mismatches.items()
            },
            "maximum_carrier_phase_mismatch_ball": maximum_carrier.str(PRECISION, more=True),
            "uniform_bounds": uniform,
            "logarithmic_quotient_series_degree": SERIES_DEGREE,
        },
        "decision": {
            "exact_hyperbolic_Morse_Fresnel_phase_proved": True,
            "same_physical_alpha_slice_used": True,
            "three_height_compact_overlap_core_proved": True,
            "retained_Fresnel_transition_required": True,
            "core_integral_remainder_proved": False,
            "outer_u_tail_proved": False,
            "propagation_across_399_events_proved": False,
            "rh_implication": False,
        },
        "next_obligation": (
            "Integrate the exact-minus-canonical retained-Fresnel core on |u|<=200 without taking "
            "the separate absolute masses, then construct analytic |u|>200 contours in the exact "
            "hyperbolic variables and compare their tails at both pi/16 faces."
        ),
        "proof_boundary": (
            "Exact coordinate algebra and compact pointwise overlap bounds at three heights only. "
            "No integrated core remainder, outer-u tail, complete one-mode join, 399-event propagation, "
            "complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved."
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
        "certified hyperbolic Morse-Fresnel overlap core: "
        "faces=3, |u|<=200, phase<0.082, amplitude<0.001; "
        f"priority={priority}"
    )


if __name__ == "__main__":
    main()
