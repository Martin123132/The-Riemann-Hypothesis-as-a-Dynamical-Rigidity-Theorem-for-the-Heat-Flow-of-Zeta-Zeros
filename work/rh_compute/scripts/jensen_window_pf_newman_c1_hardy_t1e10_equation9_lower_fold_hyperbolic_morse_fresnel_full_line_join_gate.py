#!/usr/bin/env python3
"""Certify the full-line prototype Airy/logistic join on a lifted contour."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any
from fractions import Fraction


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import flint
from flint import acb, acb_series, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate"
COMPACT_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_compact_integral_gate.json"
PARTITION_GATE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_exact_normal_saddle_partition_gate.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 80
SERIES_ORDER = 180
INNER_ORDER = 220
QUADRATIC_ORDER = 40
REAL_PANELS = 36
VERTICAL_PANELS = 4
REAL_HALF_WIDTH = arb("0.25")
VERTICAL_HALF_WIDTH = arb("0.125")
CONTOUR_HEIGHT = arb("1")
CONTOUR_RADIUS = arb("9")
C = 159_577
MODE = 39_894
Y_CUTOFF = 64
TARGET_NORMALIZED = arb("0.000019")
TARGET_PHYSICAL = arb("0.0000086")
HORIZONTAL_END = arb("15")
HORIZONTAL_PANELS_PER_SIDE = 12

TANH_COEFFICIENTS = {
    1: Fraction(1),
    3: Fraction(-1, 3),
    5: Fraction(2, 15),
    7: Fraction(-17, 315),
    9: Fraction(62, 2835),
    11: Fraction(-1382, 155925),
    13: Fraction(21844, 6081075),
    15: Fraction(-929569, 638512875),
    17: Fraction(6404582, 10854718875),
}


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


def complex_record(value: acb) -> dict[str, str]:
    return {"real_ball": value.real.str(PRECISION, more=True), "imag_ball": value.imag.str(PRECISION, more=True)}


class FullLineFace:
    def __init__(self, label: str, offset: int) -> None:
        self.label = label
        self.i = acb(0, 1)
        self.pi = arb.pi()
        self.c = arb(C)
        self.m = arb(MODE)
        self.tstar = self.pi * self.c**2 / 8
        self.beta = self.tstar ** (arb(1) / 3)
        self.tau = self.pi * self.m * (self.c - 2 * self.m)
        self.t = self.tau + arb(offset) * self.pi / 16
        self.eta = self.tstar / self.t
        self.lam = (self.eta - 1) * self.t / self.beta
        self.detuning = 4 * self.beta * (self.m - self.c / 4) / self.c
        self.sigma = 4 * self.beta / (self.pi * self.c)
        self.y = arb(Y_CUTOFF)

    def tanh_series(self, value: acb_series) -> acb_series:
        return -self.i * (self.i * value).tan()

    def cosh_series(self, value: acb_series) -> acb_series:
        return (value.exp() + (-value).exp()) / 2

    def normal_moments_y(self, a: acb_series, b: acb_series) -> tuple[acb_series, acb_series]:
        """Use the removable b=0 chart by expanding the full y exponential."""
        coefficients = [acb_series([1], prec=SERIES_ORDER), self.i * b]
        for degree in range(2, INNER_ORDER + 1):
            coefficients.append(
                (self.i * b * coefficients[-1] + 2 * self.i * a * coefficients[-2]) / degree
            )
        h0 = acb_series([], prec=SERIES_ORDER)
        h1 = acb_series([], prec=SERIES_ORDER)
        for degree, coefficient in enumerate(coefficients):
            h0 += coefficient * self.y ** (degree + 1) / (degree + 1)
            h1 += coefficient * self.y ** (degree + 2) / (degree + 2)
        return h0, h1

    def normal_moments_quadratic(self, a: acb_series, b: acb_series) -> tuple[acb_series, acb_series]:
        """Keep exp(i*b*y) exact and expand only exp(i*a*y^2)."""
        ib = self.i * b
        edge = (ib * self.y).exp()
        moments = [(edge - 1) / ib]
        for degree in range(1, 2 * QUADRATIC_ORDER + 2):
            moments.append(self.y**degree * edge / ib - degree * moments[-1] / ib)
        h0 = acb_series([], prec=SERIES_ORDER)
        h1 = acb_series([], prec=SERIES_ORDER)
        coefficient = acb_series([1], prec=SERIES_ORDER)
        for degree in range(QUADRATIC_ORDER + 1):
            if degree:
                coefficient *= self.i * a / degree
            h0 += coefficient * moments[2 * degree]
            h1 += coefficient * moments[2 * degree + 1]
        return h0, h1

    def normal_moments(self, a: acb_series, b: acb_series) -> tuple[acb_series, acb_series]:
        # The two central real panels contain the removable b=0 point.  All
        # other panels have a zero-free b disk, where the quadratic expansion
        # is dramatically better conditioned.
        if abs(b[0]) < arb("0.6"):
            return self.normal_moments_y(a, b)
        return self.normal_moments_quadratic(a, b)

    def difference_series(self, center: acb, direction: acb) -> acb_series:
        variable = acb_series([center, direction], prec=SERIES_ORDER)
        q = variable / self.beta
        tanh_q = self.tanh_series(q)
        x = (1 - tanh_q) / 2
        exact_phase = self.t * ((1 - self.eta) * q + self.eta * (q - tanh_q))
        exact_h0, exact_h1 = self.normal_moments(
            x / (2 * self.beta),
            -(self.detuning + self.beta * tanh_q),
        )
        exact = (
            self.cosh_series(q) ** (-arb(3) / 2)
            * (self.i * exact_phase).exp()
            * (exact_h0 + self.sigma / self.c * exact_h1)
        )
        canonical_h0, _ = self.normal_moments(
            acb_series([acb(1) / (4 * self.beta)], prec=SERIES_ORDER),
            -(variable + self.detuning),
        )
        canonical_phase = variable**3 / 3 - self.lam * variable
        canonical = (self.i * canonical_phase).exp() * canonical_h0
        return direction * (exact - canonical)

    def tanh_ball(self, q: acb) -> tuple[acb, acb]:
        """Return tanh(q) and q-tanh(q) with a degree-17 tail."""
        tanh_value = acb(0)
        difference = acb(0)
        for degree, coefficient in TANH_COEFFICIENTS.items():
            term = arb(coefficient.numerator) / coefficient.denominator * q**degree
            tanh_value += term
            if degree >= 3:
                difference -= term
        radius = abs(q).upper()
        require(radius < arb("0.02"), "tanh Cauchy disk escaped the local chart")
        tail = 3 * radius**18 / (1 - radius)
        error = acb(arb(0, tail), arb(0, tail))
        return tanh_value + error, difference + error

    def phase_bounds(self, z: acb) -> tuple[arb, arb, arb, acb, acb]:
        """Bound both phase moduli and the exact z-amplitude on a z disk."""
        q = z / self.beta
        tanh_q, q_minus_tanh = self.tanh_ball(q)
        y = arb(self.y / 2, self.y / 2)
        x = (1 - tanh_q) / 2
        exact_phase = (
            -self.lam * z
            + self.beta**3 * q_minus_tanh
            - self.detuning * y
            - self.beta * y * tanh_q
            + x * y**2 / (2 * self.beta)
        )
        canonical_phase = z**3 / 3 - self.lam * z - (z + self.detuning) * y + y**2 / (4 * self.beta)
        if z.imag.lower() > arb("0.5") and abs(z.real).lower() > arb("8"):
            # On the lifted horizontal contour every imaginary phase term has
            # a common positive Im(z).  Factoring it avoids a false interval
            # choice of different heights in the cubic gain and -y Im(z).
            real_lower = abs(z.real).lower()
            imag_lower = z.imag.lower()
            imag_upper = z.imag.upper()
            canonical_bracket = real_lower**2 - imag_upper**2 / 3 - abs(self.lam).upper() - self.y
            require(canonical_bracket > 0, "lifted canonical phase lost positivity")
            canonical_phase_lower = imag_lower * canonical_bracket

            scaled_gap = (real_lower**2 - imag_upper**2) / self.beta**2
            require(scaled_gap > 0, "lifted logistic phase lost its hyperbolic gap")
            exact_bracket = (
                self.beta**2 * scaled_gap / (1 + scaled_gap)
                - abs(self.lam).upper()
                - self.y
                - self.y**2 / (4 * self.beta**2)
            )
            require(exact_bracket > 0, "lifted exact phase lost positivity")
            exact_phase_lower = imag_lower * exact_bracket
            exact_exponential = (-exact_phase_lower).exp()
            canonical_exponential = (-canonical_phase_lower).exp()
        else:
            exact_exponential = (-exact_phase.imag.lower()).exp()
            canonical_exponential = (-canonical_phase.imag.lower()).exp()
        cosh_q = (q.exp() + (-q).exp()) / 2
        cosh_lower = cosh_q.abs_lower()
        require(cosh_lower > arb("0.99"), "cosh amplitude disk lost its zero-free margin")
        z_amplitude = cosh_lower ** (-arb(3) / 2) * (1 + self.sigma * self.y / self.c)
        return exact_exponential, canonical_exponential, z_amplitude, tanh_q, x

    def inner_truncation_bound(self, z: acb, use_removable: bool) -> arb:
        """Bound both exact and canonical inner truncations on one z disk."""
        q = z / self.beta
        tanh_q, _ = self.tanh_ball(q)
        x = (1 - tanh_q) / 2
        exact_a = x / (2 * self.beta)
        exact_b = -(self.detuning + self.beta * tanh_q)
        canonical_a = acb(1) / (4 * self.beta)
        canonical_b = -(z + self.detuning)

        def one(a: acb, b: acb, removable: bool) -> tuple[arb, arb]:
            if removable:
                y_radius = arb(128)
                ratio = self.y / y_radius
                maximum = (abs(a).upper() * y_radius**2 + abs(b).upper() * y_radius).exp()
                h0 = (
                    maximum
                    * self.y
                    * ratio ** (INNER_ORDER + 1)
                    / ((INNER_ORDER + 2) * (1 - ratio))
                )
                return h0, self.y * h0
            require(b.abs_lower() > 0, "quadratic inner chart disk contains b=0")
            size = abs(a).upper() * self.y**2
            linear_growth = (abs(b.imag).upper() * self.y).exp()
            h0 = (
                self.y
                * linear_growth
                * size.exp()
                * size ** (QUADRATIC_ORDER + 1)
                / arb(math.factorial(QUADRATIC_ORDER + 1))
            )
            return h0, self.y * h0

        exact_h0, exact_h1 = one(exact_a, exact_b, use_removable)
        canonical_h0, _ = one(canonical_a, canonical_b, use_removable)
        exact_exp, canonical_exp, z_amp, _, _ = self.phase_bounds(z)
        return z_amp * exact_exp * (exact_h0 + self.sigma / self.c * exact_h1) + canonical_exp * canonical_h0

    def true_integrand_bound(self, z: acb) -> arb:
        exact_exp, canonical_exp, z_amp, _, _ = self.phase_bounds(z)
        return self.y * (z_amp * exact_exp + canonical_exp)

    def panel_error(self, center: acb, half_width: arb, direction: acb) -> arb:
        is_lifted_horizontal = direction.imag.contains(0) and center.imag.lower() > arb("0.9")
        disk_radius = arb("1.2" if is_lifted_horizontal else "1.8") * half_width
        z_disk = acb(
            arb(center.real, disk_radius),
            arb(center.imag, disk_radius),
        )
        q0 = center / self.beta
        tanh0 = -self.i * (self.i * q0).tan()
        exact_b0 = -(self.detuning + self.beta * tanh0)
        canonical_b0 = -(center + self.detuning)
        use_removable = min(abs(exact_b0).upper(), abs(canonical_b0).upper()) < arb("0.6")
        inner_error = self.inner_truncation_bound(z_disk, use_removable)
        maximum = self.true_integrand_bound(z_disk) + inner_error
        ratio = half_width / disk_radius
        outer_tail = (
            2
            * maximum
            * half_width
            * ratio**SERIES_ORDER
            / ((SERIES_ORDER + 1) * (1 - ratio))
        )
        return outer_tail + 2 * half_width * inner_error


def integrate_symmetric_series(series: acb_series, half_width: arb) -> acb:
    total = acb(0)
    for degree, coefficient in enumerate(series.coeffs()):
        if degree % 2 == 0:
            total += coefficient * 2 * half_width ** (degree + 1) / (degree + 1)
    return total


def add_error(value: acb, radius: arb) -> acb:
    return acb(arb(value.real, radius), arb(value.imag, radius))


def real_core(face: FullLineFace) -> acb:
    total = acb(0)
    radius = arb(0)
    for index in range(REAL_PANELS):
        center = -CONTOUR_RADIUS + REAL_HALF_WIDTH + 2 * REAL_HALF_WIDTH * index
        series = face.difference_series(acb(center), acb(1))
        total += integrate_symmetric_series(series, REAL_HALF_WIDTH)
        radius += face.panel_error(acb(center), REAL_HALF_WIDTH, acb(1))
    return add_error(total, radius)


def connector(face: FullLineFace, side: int) -> acb:
    total = acb(0)
    radius = arb(0)
    for index in range(VERTICAL_PANELS):
        center_y = VERTICAL_HALF_WIDTH + 2 * VERTICAL_HALF_WIDTH * index
        center = acb(side * CONTOUR_RADIUS, center_y)
        series = face.difference_series(center, face.i)
        total += integrate_symmetric_series(series, VERTICAL_HALF_WIDTH)
        radius += face.panel_error(center, VERTICAL_HALF_WIDTH, face.i)
    return add_error(total, radius)


def horizontal_near(face: FullLineFace) -> acb:
    total = acb(0)
    radius = arb(0)
    for side in (-1, 1):
        start = -HORIZONTAL_END if side < 0 else CONTOUR_RADIUS
        for index in range(HORIZONTAL_PANELS_PER_SIDE):
            center_x = start + REAL_HALF_WIDTH + 2 * REAL_HALF_WIDTH * index
            center = acb(center_x, CONTOUR_HEIGHT)
            series = face.difference_series(center, acb(1))
            total += integrate_symmetric_series(series, REAL_HALF_WIDTH)
            radius += face.panel_error(center, REAL_HALF_WIDTH, acb(1))
    return add_error(total, radius)


def horizontal_far_bound(face: FullLineFace) -> arb:
    """Absolute tails on Im z=1 beyond |Re z|=15."""
    lower = HORIZONTAL_END
    y = face.y
    beta = face.beta
    lam = abs(face.lam).upper()
    amplitude = arb(1).cos()  # overwritten below; keeps an Arb value
    amplitude = (arb(1) / beta).cos() ** (-arb(3) / 2) * (1 + face.sigma * y / face.c)

    exact_shift = y + arb(1) / 2 + y**2 / (4 * beta**2) + lam
    exact_to_beta = amplitude * y * exact_shift.exp() * (-lower**2 / 2).exp() / lower
    canonical = y * (y + arb(1) / 3 + lam).exp() * (-lower**2).exp() / lower

    c = (1 - (-arb(2)).exp()) / 2
    distant_phase = (beta**2 - 1) / (2 - 1 / beta**2) - y - y**2 / (4 * beta**2) - lam
    distant_amplitude_integral = c ** (-arb(3) / 2) * 2 * beta / 3 * (-arb(3) / 2).exp()
    distant = y * (1 + face.sigma * y / face.c) * (-distant_phase).exp() * distant_amplitude_integral
    return 2 * (exact_to_beta + canonical + distant)


def certify_face(face: FullLineFace) -> dict[str, Any]:
    center = real_core(face)
    right = connector(face, 1)
    left = connector(face, -1)
    horizontal = horizontal_near(face)
    far_tail = horizontal_far_bound(face)
    # The left connector is traversed downward in the deformed full contour.
    finite_contour = center + right - left + horizontal
    full_line = add_error(finite_contour, far_tail)
    physical = arb(2).sqrt() / face.pi * full_line
    require(abs(full_line) < TARGET_NORMALIZED, f"{face.label} normalized full-line target failed")
    require(abs(physical) < TARGET_PHYSICAL, f"{face.label} physical full-line target failed")
    return {
        "label": face.label,
        "height_ball": face.t.str(PRECISION, more=True),
        "real_core_ball": complex_record(center),
        "right_connector_ball": complex_record(right),
        "left_connector_ball": complex_record(left),
        "horizontal_near_ball": complex_record(horizontal),
        "finite_contour_ball": complex_record(finite_contour),
        "horizontal_far_tail_absolute_bound": far_tail.str(PRECISION, more=True),
        "normalized_full_line_difference_ball": complex_record(full_line),
        "normalized_full_line_difference_absolute_ball": abs(full_line).str(PRECISION, more=True),
        "physical_full_line_difference_ball": complex_record(physical),
        "physical_full_line_difference_absolute_ball": abs(physical).str(PRECISION, more=True),
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = "\n".join(
        f"{row['label']}: normalized={row['normalized_full_line_difference_absolute_ball']}, "
        f"physical={row['physical_full_line_difference_absolute_ball']}"
        for row in artifact["certified_faces"]
    )
    return f"""# Full-line hyperbolic Morse--Fresnel prototype join

Date: 2026-08-12

Status: full-line one-mode prototype join certified at three heights; this is
not a proof of 399-event propagation or a complete `T_upper` estimate

Return both the exact retained-Fresnel logistic chart and its Airy model to
the same physical variables `(z,y)`.  Their full real domains then coincide:

```text
Delta I(t)=int_R int_0^64 [A_t(z,y) exp(i Theta_t(z,y))
                           -exp(i Theta_0(z,y))] dy dz. (FLJ1)
```

All normal saddles lie inside `|z|<9`.  Deform each real tail to
`Im z=1`, keeping the vertical connectors at `z=+/-9`.  On every finite
segment, Arb Taylor series integrate the signed difference coefficientwise.
The inner quadratic-exponential integral and its first moment are generated
by the exact recurrence

```text
c_0=1, c_1=i b,
c_n=[i b c_(n-1)+2 i a c_(n-2)]/n,                    (FLJ2)
```

followed by exact termwise integration on `0<=y<=64`.

The signed horizontal pieces are Taylor-integrated through `|Re z|=15`.
Beyond that point, the exact phase obeys a logistic-Gaussian lower bound and
the canonical phase a Gaussian lower bound; their combined absolute tail is
recorded in the machine artifact.  Combining core, connectors, and both
tails gives

```text
{rows}
```

and hence uniformly at the event and both `pi/16` faces,

```text
|Delta I(t)| < 0.000019,
(sqrt(2)/pi)|Delta I(t)| < 0.0000086.                 (FLJ3)
```

The `pi` normalization is inherited from the Kummer quadratic phase and
integer Fourier--Poisson character.  No fitted geometric constant is used.

## Boundary

This gate certifies one prototype crossing mode at three heights.  It does
not prove uniformity through all 399 turning events, control the remaining
ordinary modes, establish complete `T_upper`, or prove `Lambda<=0`, RH, or a
prize-level conclusion.
"""


def main() -> None:
    flint.ctx.dps = PRECISION
    flint.ctx.cap = SERIES_ORDER + 8
    priority = set_low_priority()
    require(COMPACT_GATE.is_file(), "missing compact-integral dependency")
    require(PARTITION_GATE.is_file(), "missing saddle-partition dependency")
    require(CHECKER.is_file(), "missing independent checker")

    faces = [FullLineFace("lower_face", -1), FullLineFace("event", 0), FullLineFace("upper_face", 1)]
    certified = [certify_face(face) for face in faces]
    artifact = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "full_line_one_mode_hyperbolic_Morse_Fresnel_Airy_join_certified",
        "precision_decimal_digits": PRECISION,
        "resource_policy": {"worker_cap": 1, "active_workers": 1, "process_priority": priority},
        "method": {
            "common_chart": "physical (z,y) chart with the complete retained endpoint/Fresnel current",
            "real_core": "|z|<=9 split into 36 Arb Taylor panels",
            "connectors": "z=+/-9+i v, 0<=v<=1, four Taylor panels per connector",
            "inner_integral": f"quadratic-exponential coefficient recurrence through degree {INNER_ORDER}",
            "quadratic_expansion_order_off_removable_chart": QUADRATIC_ORDER,
            "outer_series_order": SERIES_ORDER,
            "horizontal_near_end": HORIZONTAL_END.str(PRECISION, more=True),
            "outer_Cauchy_disk_ratio": "panel_half_width/disk_radius=5/9",
        },
        "certified_faces": certified,
        "claims": {
            "common_full_line_chart_exact": True,
            "all_normal_saddles_inside_abs_z_9": True,
            "signed_full_line_difference_integrated_before_absolute_values": True,
            "normalized_full_line_difference_below_0_000019": True,
            "physical_full_line_difference_below_0_0000086": True,
            "prototype_one_mode_join_proved": True,
            "all_399_events_proved": False,
            "complete_T_upper_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": (
            "One-mode full-line prototype join at three heights only. No 399-event propagation, remaining-mode theorem, "
            "complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved."
        ),
        "next_action": (
            "Express the contour constants as functions of the odd event parameter |4m-C|/C, prove monotonic or finite-exception control, "
            "and propagate the one-mode join through all 399 event buffers."
        ),
        "dependencies": {
            "compact_gate": {"path": relative(COMPACT_GATE), "sha256": file_hash(COMPACT_GATE)},
            "partition_gate": {"path": relative(PARTITION_GATE), "sha256": file_hash(PARTITION_GATE)},
        },
        "artifacts": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
        },
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8")
    artifact["artifacts"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    maximum = max(arb(row["normalized_full_line_difference_absolute_ball"]) for row in certified)
    print(f"certified full-line prototype join: faces=3, max |difference|={maximum}; priority={priority}", flush=True)


if __name__ == "__main__":
    main()
