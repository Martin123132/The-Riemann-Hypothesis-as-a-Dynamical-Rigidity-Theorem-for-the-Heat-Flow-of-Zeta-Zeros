#!/usr/bin/env python3
"""Certify the beta^-2 corrected normal form at both extreme events."""

from __future__ import annotations

import hashlib
import importlib.util
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

import flint
from flint import acb, acb_series, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_corrected_normal_form_gate"
FULL_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate"
PILOT_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_full_line_pilot"
FULL_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{FULL_STEM}.py"
PILOT = REPO_ROOT / f"work/rh_compute/results/{PILOT_STEM}.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 80
C = 159_577
EXTREMES = ((398, 39_695), (397, 40_093))
TARGET_NORMALIZED = arb("0.000019")
TARGET_PHYSICAL = arb("0.0000086")


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


def load_full_module():
    spec = importlib.util.spec_from_file_location("corrected_normal_form_source", FULL_BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load full-line builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_corrected_face(module):
    class CorrectedFace(module.FullLineFace):
        def normal_moments_three(self, a: acb_series, b: acb_series) -> tuple[acb_series, acb_series, acb_series]:
            if abs(b[0]) < arb("0.6"):
                coefficients = [acb_series([1], prec=module.SERIES_ORDER), self.i * b]
                for degree in range(2, module.INNER_ORDER + 1):
                    coefficients.append((self.i * b * coefficients[-1] + 2 * self.i * a * coefficients[-2]) / degree)
                moments = [acb_series([], prec=module.SERIES_ORDER) for _ in range(3)]
                for degree, coefficient in enumerate(coefficients):
                    for power in range(3):
                        moments[power] += coefficient * self.y ** (degree + power + 1) / (degree + power + 1)
                return moments[0], moments[1], moments[2]

            ib = self.i * b
            edge = (ib * self.y).exp()
            y_moments = [(edge - 1) / ib]
            for degree in range(1, 2 * module.QUADRATIC_ORDER + 3):
                y_moments.append(self.y**degree * edge / ib - degree * y_moments[-1] / ib)
            moments = [acb_series([], prec=module.SERIES_ORDER) for _ in range(3)]
            coefficient = acb_series([1], prec=module.SERIES_ORDER)
            for degree in range(module.QUADRATIC_ORDER + 1):
                if degree:
                    coefficient *= self.i * a / degree
                for power in range(3):
                    moments[power] += coefficient * y_moments[2 * degree + power]
            return moments[0], moments[1], moments[2]

        def correction_coefficients(self, variable: acb_series) -> tuple[acb_series, acb_series, acb_series]:
            p0 = -arb(3) * variable**2 / 4 - 2 * self.i * variable**5 / 15
            p1 = arb(1) / 2 + self.i * variable**3 / 3
            p2 = -self.i * variable / 4
            return p0, p1, p2

        def difference_series(self, center: acb, direction: acb) -> acb_series:
            variable = acb_series([center, direction], prec=module.SERIES_ORDER)
            q = variable / self.beta
            tanh_q = self.tanh_series(q)
            x = (1 - tanh_q) / 2
            exact_phase = self.t * ((1 - self.eta) * q + self.eta * (q - tanh_q))
            exact_h0, exact_h1 = self.normal_moments(x / (2 * self.beta), -(self.detuning + self.beta * tanh_q))
            exact = (
                self.cosh_series(q) ** (-arb(3) / 2)
                * (self.i * exact_phase).exp()
                * (exact_h0 + self.sigma / self.c * exact_h1)
            )

            h0, h1, h2 = self.normal_moments_three(
                acb_series([acb(1) / (4 * self.beta)], prec=module.SERIES_ORDER),
                -(variable + self.detuning),
            )
            p0, p1, p2 = self.correction_coefficients(variable)
            canonical_phase = variable**3 / 3 - self.lam * variable
            corrected = (self.i * canonical_phase).exp() * (h0 + (p0 * h0 + p1 * h1 + p2 * h2) / self.beta**2)
            return direction * (exact - corrected)

        def correction_bound(self, z: acb) -> arb:
            radius = abs(z).upper()
            y = self.y
            return (
                3 * radius**2 / 4
                + y / 2
                + 2 * radius**5 / 15
                + radius**3 * y / 3
                + radius * y**2 / 4
            ) / self.beta**2

        def true_integrand_bound(self, z: acb) -> arb:
            exact_exp, canonical_exp, z_amp, _, _ = self.phase_bounds(z)
            return self.y * (z_amp * exact_exp + canonical_exp * (1 + self.correction_bound(z)))

        def inner_truncation_bound(self, z: acb, use_removable: bool) -> arb:
            q = z / self.beta
            tanh_q, _ = self.tanh_ball(q)
            x = (1 - tanh_q) / 2
            exact_a = x / (2 * self.beta)
            exact_b = -(self.detuning + self.beta * tanh_q)
            canonical_a = acb(1) / (4 * self.beta)
            canonical_b = -(z + self.detuning)

            def one(a: acb, b: acb, removable: bool) -> tuple[arb, arb, arb]:
                if removable:
                    y_radius = arb(128)
                    ratio = self.y / y_radius
                    maximum = (abs(a).upper() * y_radius**2 + abs(b).upper() * y_radius).exp()
                    h0 = maximum * self.y * ratio ** (module.INNER_ORDER + 1) / ((module.INNER_ORDER + 2) * (1 - ratio))
                else:
                    require(b.abs_lower() > 0, "quadratic corrected chart disk contains b=0")
                    size = abs(a).upper() * self.y**2
                    growth = (abs(b.imag).upper() * self.y).exp()
                    h0 = self.y * growth * size.exp() * size ** (module.QUADRATIC_ORDER + 1) / arb(math.factorial(module.QUADRATIC_ORDER + 1))
                return h0, self.y * h0, self.y**2 * h0

            exact_h0, exact_h1, _ = one(exact_a, exact_b, use_removable)
            canonical_h0, canonical_h1, canonical_h2 = one(canonical_a, canonical_b, use_removable)
            radius = abs(z).upper()
            p0 = 3 * radius**2 / 4 + 2 * radius**5 / 15
            p1 = arb(1) / 2 + radius**3 / 3
            p2 = radius / 4
            corrected_error = canonical_h0 + (p0 * canonical_h0 + p1 * canonical_h1 + p2 * canonical_h2) / self.beta**2
            exact_exp, canonical_exp, z_amp, _, _ = self.phase_bounds(z)
            return z_amp * exact_exp * (exact_h0 + self.sigma / self.c * exact_h1) + canonical_exp * corrected_error

    return CorrectedFace


def configure(module, checker: bool = False) -> None:
    if checker:
        module.SERIES_ORDER = 200
        module.INNER_ORDER = 240
        module.QUADRATIC_ORDER = 46
        module.REAL_PANELS = 280
        module.REAL_HALF_WIDTH = arb("0.05")
        module.VERTICAL_PANELS = 20
        module.VERTICAL_HALF_WIDTH = arb("0.025")
        module.HORIZONTAL_PANELS_PER_SIDE = 70
    else:
        module.REAL_PANELS = 224
        module.REAL_HALF_WIDTH = arb("0.0625")
        module.VERTICAL_PANELS = 16
        module.VERTICAL_HALF_WIDTH = arb("0.03125")
        module.HORIZONTAL_PANELS_PER_SIDE = 56
    module.CONTOUR_RADIUS = arb(14)
    module.HORIZONTAL_END = arb(21)
    flint.ctx.cap = module.SERIES_ORDER + 8


def event_face(face_class, module, event_index: int, mode: int):
    face = face_class(f"event_{event_index}_mode_{mode}", 0)
    face.m = arb(mode)
    face.tau = face.pi * face.m * (face.c - 2 * face.m)
    face.t = face.tau
    face.eta = face.tstar / face.t
    face.lam = (face.eta - 1) * face.t / face.beta
    face.detuning = 4 * face.beta * (face.m - face.c / 4) / face.c
    return face


def gaussian_moment_tail(lower: arb, degree: int) -> arb:
    values = [(-lower**2).exp() / (2 * lower), (-lower**2).exp() / 2]
    if degree < 2:
        return values[degree]
    for index in range(2, degree + 1):
        values.append(lower ** (index - 1) * (-lower**2).exp() / 2 + arb(index - 1) * values[index - 2] / 2)
    return values[degree]


def corrected_far_bound(module, face) -> arb:
    base = module.horizontal_far_bound(face)
    lower = module.HORIZONTAL_END
    scale = (1 + 1 / lower**2).sqrt()
    y = face.y
    prefactor = 2 * y * (y + arb(1) / 3 + abs(face.lam).upper()).exp() / face.beta**2
    polynomial = (
        3 * scale**2 * gaussian_moment_tail(lower, 2) / 4
        + y * gaussian_moment_tail(lower, 0) / 2
        + 2 * scale**5 * gaussian_moment_tail(lower, 5) / 15
        + y * scale**3 * gaussian_moment_tail(lower, 3) / 3
        + y**2 * scale * gaussian_moment_tail(lower, 1) / 4
    )
    return base + prefactor * polynomial


def complex_record(value: acb) -> dict[str, str]:
    return {"real_ball": value.real.str(PRECISION, more=True), "imag_ball": value.imag.str(PRECISION, more=True)}


def integrate_face(module, face, event_index: int, mode: int) -> dict[str, Any]:
    core = module.real_core(face)
    right = module.connector(face, 1)
    left = module.connector(face, -1)
    horizontal = module.horizontal_near(face)
    far = corrected_far_bound(module, face)
    full = module.add_error(core + right - left + horizontal, far)
    physical = arb(2).sqrt() / face.pi * full
    normalized_abs = abs(full)
    physical_abs = abs(physical)
    return {
        "event_index": event_index,
        "mode": mode,
        "signed_odd": 4 * mode - C,
        "detuning_ball": face.detuning.str(PRECISION, more=True),
        "normalized_corrected_difference_ball": complex_record(full),
        "normalized_corrected_absolute_ball": normalized_abs.str(PRECISION, more=True),
        "normalized_corrected_absolute_lower": normalized_abs.lower().str(PRECISION, more=True),
        "normalized_corrected_absolute_upper": normalized_abs.upper().str(PRECISION, more=True),
        "physical_corrected_difference_ball": complex_record(physical),
        "physical_corrected_absolute_ball": physical_abs.str(PRECISION, more=True),
        "physical_corrected_absolute_lower": physical_abs.lower().str(PRECISION, more=True),
        "physical_corrected_absolute_upper": physical_abs.upper().str(PRECISION, more=True),
        "corrected_far_tail_absolute_bound": far.str(PRECISION, more=True),
        "prototype_targets_proved": normalized_abs < TARGET_NORMALIZED and physical_abs < TARGET_PHYSICAL,
        "prototype_targets_disproved": normalized_abs.lower() > TARGET_NORMALIZED and physical_abs.lower() > TARGET_PHYSICAL,
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = "\n".join(
        f"event {row['event_index']} (mode {row['mode']}): normalized={row['normalized_corrected_absolute_ball']}, "
        f"physical={row['physical_corrected_absolute_ball']}"
        for row in artifact["certified_extremes"]
    )
    return f"""# Extreme-event corrected normal form

Date: 2026-08-12

Status: rigorous beta^-2 corrected diagnostic at both extreme event
detunings; not a proof of uniform 399-event propagation

At an event `t=tau_m=beta^3-beta*d_m^2`, expansion of the exact common
hyperbolic chart gives

```text
Theta_exact-Theta_0=beta^-2 r_1+O(beta^-4),
A_exact=1+beta^-2 a_1+O(beta^-4),

a_1=y/2-3z^2/4,
r_1=-2z^5/15+z^3y/3-zy^2/4.                         (ECN1)
```

The corrected comparison integrand is therefore

```text
exp(i Theta_0)[1+beta^-2(a_1+i r_1)].                (ECN2)
```

The inner quadratic exponential is integrated with rigorous moments through
`y^2`; real, connector, and lifted-contour Taylor tails are enclosed on the
common `R=14`, `Im z=1`, `|Re z|<=21` contour.  A separate polynomial-
Gaussian estimate controls the corrected far tail.  At the two extreme
event parameters:

```text
{rows}
```

Both recover the prototype targets: `{artifact['claims']['both_extreme_events_below_prototype_targets']}`.
Both rigorously exceed those targets: `{artifact['claims']['both_extreme_events_exceed_prototype_targets']}`.
This tests the first corrected comparison only at the two extreme event
points.  It does not enclose intermediate detunings or certify continuous
height on those cells.  It does not propagate all 399 events.  It does not
establish complete `T_upper` or prove `Lambda<=0` or RH.  It does not
establish a prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(FULL_BUILDER.is_file() and PILOT.is_file() and CHECKER.is_file(), "missing dependency or checker")
    flint.ctx.dps = PRECISION
    module = load_full_module()
    configure(module)
    face_class = make_corrected_face(module)
    certified = [integrate_face(module, event_face(face_class, module, index, mode), index, mode) for index, mode in EXTREMES]
    both_below = all(row["prototype_targets_proved"] for row in certified)
    both_above = all(row["prototype_targets_disproved"] for row in certified)
    artifact: dict[str, Any] = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "extreme_event_beta_minus_2_corrected_normal_form_diagnostic_complete",
        "resource_policy": {"workers": 1, "process_priority": priority},
        "correction": {
            "amplitude_a1": "y/2-3*z^2/4",
            "phase_r1": "-2*z^5/15+z^3*y/3-z*y^2/4",
            "corrected_integrand": "exp(i*Theta0)*(1+(a1+i*r1)/beta^2)",
            "provenance": "Taylor coefficients of cosh(z/beta)^(-3/2), tanh(z/beta), x=(1-tanh(z/beta))/2, and sigma/C=1/(2beta^2)",
        },
        "contour": {"radius": 14, "height": 1, "horizontal_end": 21, "real_panels": 224, "connector_panels_per_side": 16, "horizontal_panels_per_side": 56},
        "certified_extremes": certified,
        "claims": {
            "beta_minus_2_correction_derived": True,
            "both_extreme_events_below_prototype_targets": both_below,
            "both_extreme_events_exceed_prototype_targets": both_above,
            "uniform_intermediate_detuning_proved": False,
            "continuous_height_for_all_events_proved": False,
            "all_399_event_propagation_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": "Two extreme event centers only. No intermediate-detuning enclosure, all-event continuous-height theorem, 399-event propagation, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "next_action": "Add the complete beta^-4 coefficient a2+i*r2+i*a1*r1-r1^2/2, then retest the same extreme-event contour before any interval-atlas promotion.",
        "dependencies": {
            "uncorrected_extreme_pilot": {"path": relative(PILOT), "sha256": file_hash(PILOT)},
            "full_line_builder": {"path": relative(FULL_BUILDER), "sha256": file_hash(FULL_BUILDER)},
        },
        "artifacts": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
        "runtime": {"elapsed_seconds": round(time.perf_counter() - started, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["artifacts"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(
        "certified corrected extreme events: "
        + ", ".join(f"n={row['event_index']} |D_corr|={row['normalized_corrected_absolute_ball']}" for row in certified),
        flush=True,
    )


if __name__ == "__main__":
    main()
