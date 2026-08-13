#!/usr/bin/env python3
"""Certify the complete beta^-4 normal form at both extreme events."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

import flint
from flint import acb, acb_series, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_second_order_normal_form_gate"
FIRST_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_corrected_normal_form_gate"
FIRST_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{FIRST_STEM}.py"
FIRST_RESULT = REPO_ROOT / f"work/rh_compute/results/{FIRST_STEM}.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)
PRECISION = 80


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_first():
    spec = importlib.util.spec_from_file_location("first_corrected_normal_form", FIRST_BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load first-order builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_second_order_face(first, module):
    first_class = first.make_corrected_face(module)

    class SecondOrderFace(first_class):
        def normal_moments_five(self, a: acb_series, b: acb_series) -> list[acb_series]:
            if abs(b[0]) < arb("0.6"):
                coefficients = [acb_series([1], prec=module.SERIES_ORDER), self.i * b]
                for degree in range(2, module.INNER_ORDER + 1):
                    coefficients.append((self.i * b * coefficients[-1] + 2 * self.i * a * coefficients[-2]) / degree)
                moments = [acb_series([], prec=module.SERIES_ORDER) for _ in range(5)]
                for degree, coefficient in enumerate(coefficients):
                    for power in range(5):
                        moments[power] += coefficient * self.y ** (degree + power + 1) / (degree + power + 1)
                return moments

            ib = self.i * b
            edge = (ib * self.y).exp()
            y_moments = [(edge - 1) / ib]
            for degree in range(1, 2 * module.QUADRATIC_ORDER + 5):
                y_moments.append(self.y**degree * edge / ib - degree * y_moments[-1] / ib)
            moments = [acb_series([], prec=module.SERIES_ORDER) for _ in range(5)]
            coefficient = acb_series([1], prec=module.SERIES_ORDER)
            for degree in range(module.QUADRATIC_ORDER + 1):
                if degree:
                    coefficient *= self.i * a / degree
                for power in range(5):
                    moments[power] += coefficient * y_moments[2 * degree + power]
            return moments

        def second_coefficients(self, z: acb_series) -> list[acb_series]:
            return [
                -2 * z**10 / 225 + 97 * self.i * z**7 / 630 + 13 * z**4 / 32,
                2 * z**8 / 45 - 9 * self.i * z**5 / 20 - 3 * z**2 / 8,
                -4 * z**6 / 45 + 7 * self.i * z**3 / 16,
                z**4 / 12 - self.i * z / 8,
                -z**2 / 32,
            ]

        def second_bound(self, z: acb) -> arb:
            r = abs(z).upper()
            y = self.y
            return (
                2 * r**10 / 225
                + 97 * r**7 / 630
                + 13 * r**4 / 32
                + y * (2 * r**8 / 45 + 9 * r**5 / 20 + 3 * r**2 / 8)
                + y**2 * (4 * r**6 / 45 + 7 * r**3 / 16)
                + y**3 * (r**4 / 12 + r / 8)
                + y**4 * r**2 / 32
            )

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
            moments = self.normal_moments_five(
                acb_series([acb(1) / (4 * self.beta)], prec=module.SERIES_ORDER),
                -(variable + self.detuning),
            )
            first_coefficients = list(self.correction_coefficients(variable))
            second_coefficients = self.second_coefficients(variable)
            corrected_inner = moments[0]
            corrected_inner += sum((coefficient * moments[index] for index, coefficient in enumerate(first_coefficients)), acb_series([], prec=module.SERIES_ORDER)) / self.beta**2
            corrected_inner += sum((coefficient * moments[index] for index, coefficient in enumerate(second_coefficients)), acb_series([], prec=module.SERIES_ORDER)) / self.beta**4
            canonical_phase = variable**3 / 3 - self.lam * variable
            corrected = (self.i * canonical_phase).exp() * corrected_inner
            return direction * (exact - corrected)

        def correction_bound(self, z: acb) -> arb:
            return super().correction_bound(z) + self.second_bound(z) / self.beta**4

        def inner_truncation_bound(self, z: acb, use_removable: bool) -> arb:
            q = z / self.beta
            tanh_q, _ = self.tanh_ball(q)
            x = (1 - tanh_q) / 2
            exact_a = x / (2 * self.beta)
            exact_b = -(self.detuning + self.beta * tanh_q)
            canonical_a = acb(1) / (4 * self.beta)
            canonical_b = -(z + self.detuning)

            def h0_error(a: acb, b: acb, removable: bool) -> arb:
                if removable:
                    y_radius = arb(128)
                    ratio = self.y / y_radius
                    maximum = (abs(a).upper() * y_radius**2 + abs(b).upper() * y_radius).exp()
                    return maximum * self.y * ratio ** (module.INNER_ORDER + 1) / ((module.INNER_ORDER + 2) * (1 - ratio))
                require(b.abs_lower() > 0, "second-order quadratic chart disk contains b=0")
                size = abs(a).upper() * self.y**2
                growth = (abs(b.imag).upper() * self.y).exp()
                return self.y * growth * size.exp() * size ** (module.QUADRATIC_ORDER + 1) / arb(math.factorial(module.QUADRATIC_ORDER + 1))

            exact_h0 = h0_error(exact_a, exact_b, use_removable)
            canonical_h0 = h0_error(canonical_a, canonical_b, use_removable)
            exact_exp, canonical_exp, z_amp, _, _ = self.phase_bounds(z)
            exact_error = z_amp * exact_exp * exact_h0 * (1 + self.sigma * self.y / self.c)
            corrected_error = canonical_exp * canonical_h0 * (1 + super().correction_bound(z) + self.second_bound(z) / self.beta**4)
            return exact_error + corrected_error

    return SecondOrderFace


def positive_second_polynomial(y: arb) -> dict[int, arb]:
    return {
        10: arb(2) / 225,
        8: 2 * y / 45,
        7: arb(97) / 630,
        6: 4 * y**2 / 45,
        5: 9 * y / 20,
        4: arb(13) / 32 + y**3 / 12,
        3: 7 * y**2 / 16,
        2: 3 * y / 8 + y**4 / 32,
        1: y**3 / 8,
    }


def second_order_far_bound(first, module, face) -> arb:
    first_bound = first.corrected_far_bound(module, face)
    lower = module.HORIZONTAL_END
    scale = (1 + 1 / lower**2).sqrt()
    prefactor = 2 * face.y * (face.y + arb(1) / 3 + abs(face.lam).upper()).exp() / face.beta**4
    extra = sum(
        (coefficient * scale**degree * first.gaussian_moment_tail(lower, degree) for degree, coefficient in positive_second_polynomial(face.y).items()),
        arb(0),
    )
    return first_bound + prefactor * extra


def integrate(first, module, face, event_index: int, mode: int) -> dict[str, Any]:
    core = module.real_core(face)
    right = module.connector(face, 1)
    left = module.connector(face, -1)
    horizontal = module.horizontal_near(face)
    far = second_order_far_bound(first, module, face)
    full = module.add_error(core + right - left + horizontal, far)
    physical = arb(2).sqrt() / face.pi * full
    require(abs(full) < first.TARGET_NORMALIZED, f"event {event_index} second-order normalized target failed")
    require(abs(physical) < first.TARGET_PHYSICAL, f"event {event_index} second-order physical target failed")
    return {
        "event_index": event_index,
        "mode": mode,
        "signed_odd": 4 * mode - first.C,
        "normalized_second_order_difference_ball": first.complex_record(full),
        "normalized_second_order_absolute_ball": abs(full).str(PRECISION, more=True),
        "physical_second_order_difference_ball": first.complex_record(physical),
        "physical_second_order_absolute_ball": abs(physical).str(PRECISION, more=True),
        "second_order_far_tail_absolute_bound": far.str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = "\n".join(
        f"event {row['event_index']} (mode {row['mode']}): normalized={row['normalized_second_order_absolute_ball']}, "
        f"physical={row['physical_second_order_absolute_ball']}"
        for row in artifact["certified_extremes"]
    )
    return f"""# Extreme-event second-order normal form

Date: 2026-08-12

Status: complete beta^-4 comparison certified at both extreme event
detunings; not a proof of uniform 399-event propagation

Write the exact amplitude and phase relative to the event Airy chart as

```text
A=1+beta^-2 a_1+beta^-4 a_2+O(beta^-6),
Theta=Theta_0+beta^-2 r_1+beta^-4 r_2+O(beta^-6),

a_1=y/2-3z^2/4,
r_1=-2z^5/15+z^3y/3-zy^2/4,
a_2=13z^4/32-3z^2y/8,
r_2=17z^7/315-2z^5y/15+z^3y^2/12.               (ESN1)
```

The complete second-order comparison multiplier is

```text
1+beta^-2(a_1+i r_1)
 +beta^-4(a_2+i r_2+i a_1 r_1-r_1^2/2).             (ESN2)
```

Its degree-four normal polynomial is integrated signed with moments
`H_0,...,H_4`.  The common event contour and a degree-ten polynomial-
Gaussian far-tail bound give

```text
{rows}
```

Both extreme points are below the prototype targets.  This is a rigorous
two-point endpoint result only.  It does not enclose the continuous detuning
interval or intermediate discrete events.  It does not extend exact height
transport to every corrected cell.  It does not establish complete `T_upper`
or prove `Lambda<=0` or RH.  It does not establish a prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    require(FIRST_BUILDER.is_file() and FIRST_RESULT.is_file() and CHECKER.is_file(), "missing dependency or checker")
    first = load_first()
    priority = first.set_low_priority()
    flint.ctx.dps = PRECISION
    module = first.load_full_module()
    first.configure(module)
    face_class = make_second_order_face(first, module)
    certified = [integrate(first, module, first.event_face(face_class, module, index, mode), index, mode) for index, mode in first.EXTREMES]
    artifact: dict[str, Any] = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "extreme_event_beta_minus_4_second_order_normal_form_certified",
        "resource_policy": {"workers": 1, "process_priority": priority},
        "normal_form": {
            "a1": "y/2-3*z^2/4",
            "r1": "-2*z^5/15+z^3*y/3-z*y^2/4",
            "a2": "13*z^4/32-3*z^2*y/8",
            "r2": "17*z^7/315-2*z^5*y/15+z^3*y^2/12",
            "second_multiplier": "a2+i*r2+i*a1*r1-r1^2/2",
        },
        "certified_extremes": certified,
        "claims": {
            "complete_beta_minus_4_multiplier_derived": True,
            "both_extreme_events_below_prototype_targets": True,
            "uniform_intermediate_detuning_proved": False,
            "all_399_event_propagation_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": "Two extreme event centers only. No intermediate-detuning enclosure, all-event continuous-height theorem, 399-event propagation, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "next_action": "Certify the second-order residual on the complete compact detuning interval, then derive uniform corrected height moments for every pi/16 event cell.",
        "dependencies": {"first_order_gate": {"path": relative(FIRST_RESULT), "sha256": file_hash(FIRST_RESULT)}},
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
        "certified second-order extreme events: "
        + ", ".join(f"n={row['event_index']} |D_2|={row['normalized_second_order_absolute_ball']}" for row in certified),
        flush=True,
    )


if __name__ == "__main__":
    main()
