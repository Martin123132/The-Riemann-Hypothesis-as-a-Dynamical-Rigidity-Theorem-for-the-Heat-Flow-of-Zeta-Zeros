#!/usr/bin/env python3
"""Certify the beta^-4 normal-form residual at all 399 event centers."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, wait
from fractions import Fraction
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

import flint
from flint import acb, acb_series, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_second_order_all_event_atlas_gate"
SECOND_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_second_order_normal_form_gate"
SECOND_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{SECOND_STEM}.py"
SECOND_RESULT = REPO_ROOT / f"work/rh_compute/results/{SECOND_STEM}.json"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
CACHE = REPO_ROOT / f"work/rh_compute/results/cache/{STEM}_rows.jsonl"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).resolve()
CHECKER = BUILDER.with_name("check_" + BUILDER.name)

PRECISION = 70
C = 159_577
Q = (C - 1) // 4
EVENT_COUNT = 399
CONFIG_VERSION = "shifted-rational-v2"
CENTER_DENOMINATOR = 2048
MAX_PANEL_WIDTH = Fraction(1, 8)
CENTRAL_HALF_WIDTH = Fraction(1, 16)
REMOVABLE_SWITCH = arb("0.01")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mode_for_event(index: int) -> int:
    return Q - index // 2 if index % 2 == 0 else Q + (index + 1) // 2


def signed_odd_for_event(index: int) -> int:
    return (-1 if index % 2 == 0 else 1) * (2 * index + 1)


def configure(module, checker: bool = False) -> None:
    if checker:
        module.SERIES_ORDER = 70
        module.INNER_ORDER = 96
        module.QUADRATIC_ORDER = 30
    else:
        module.SERIES_ORDER = 47
        module.INNER_ORDER = 76
        module.QUADRATIC_ORDER = 24
    module.CONTOUR_RADIUS = arb(14)
    module.HORIZONTAL_END = arb(21)
    module.VERTICAL_PANELS = 16
    module.VERTICAL_HALF_WIDTH = arb("0.03125")
    module.HORIZONTAL_PANELS_PER_SIDE = 56
    module.REAL_HALF_WIDTH = arb("0.0625")
    flint.ctx.cap = module.SERIES_ORDER + 8


def make_atlas_face(second, first, module):
    second_class = second.make_second_order_face(first, module)

    class AtlasFace(second_class):
        def phase_bounds(self, z: acb):
            # On a disk meeting the real axis, both phases are real at every
            # real base point.  Integrating vertically bounds Im(phase) by the
            # disk height times the derivative.  The factored derivatives
            # retain the event-saddle cancellations that direct interval
            # evaluation of z^3/3-d^2*z would otherwise discard.
            if z.imag.contains(0):
                q = z / self.beta
                tanh_q, _ = self.tanh_ball(q)
                sech_squared = 1 - tanh_q**2
                y = arb(self.y / 2, self.y / 2)
                exact_derivative = (
                    (self.beta * tanh_q - self.detuning)
                    * (self.beta * tanh_q + self.detuning)
                    - sech_squared * (y + y**2 / (4 * self.beta**2))
                )
                canonical_derivative = (z - self.detuning) * (z + self.detuning) - y
                height = abs(z.imag).upper()
                exact_exponential = (height * abs(exact_derivative).upper()).exp()
                canonical_exponential = (height * abs(canonical_derivative).upper()).exp()
                cosh_q = (q.exp() + (-q).exp()) / 2
                cosh_lower = cosh_q.abs_lower()
                require(cosh_lower > arb("0.99"), "cosh amplitude disk lost its zero-free margin")
                z_amplitude = cosh_lower ** (-arb(3) / 2) * (1 + self.sigma * self.y / self.c)
                x = (1 - tanh_q) / 2
                return exact_exponential, canonical_exponential, z_amplitude, tanh_q, x
            if z.imag.lower() > 0 and abs(z.real).lower() > arb(13):
                q = z / self.beta
                tanh_q, _ = self.tanh_ball(q)
                real_lower = abs(z.real).lower()
                imag_lower = z.imag.lower()
                imag_upper = z.imag.upper()
                canonical_bracket = (
                    real_lower**2
                    - imag_upper**2 / 3
                    - abs(self.lam).upper()
                    - self.y
                )
                scaled_gap = (real_lower**2 - imag_upper**2) / self.beta**2
                exact_bracket = (
                    self.beta**2 * scaled_gap / (1 + scaled_gap)
                    - abs(self.lam).upper()
                    - self.y
                    - self.y**2 / (4 * self.beta**2)
                )
                require(canonical_bracket > 0, "connector canonical phase lost positivity")
                require(exact_bracket > 0, "connector exact phase lost positivity")
                exact_exponential = (-imag_lower * exact_bracket).exp()
                canonical_exponential = (-imag_lower * canonical_bracket).exp()
                cosh_q = (q.exp() + (-q).exp()) / 2
                cosh_lower = cosh_q.abs_lower()
                require(cosh_lower > arb("0.99"), "cosh amplitude disk lost its zero-free margin")
                z_amplitude = cosh_lower ** (-arb(3) / 2) * (1 + self.sigma * self.y / self.c)
                x = (1 - tanh_q) / 2
                return exact_exponential, canonical_exponential, z_amplitude, tanh_q, x
            return super().phase_bounds(z)

        def normal_moments(self, a: acb_series, b: acb_series):
            if abs(b[0]) < REMOVABLE_SWITCH:
                return self.normal_moments_y(a, b)
            return self.normal_moments_quadratic(a, b)

        def normal_moments_five(self, a: acb_series, b: acb_series) -> list[acb_series]:
            if abs(b[0]) < REMOVABLE_SWITCH:
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

        def panel_error(self, center: acb, half_width: arb, direction: acb) -> arb:
            is_lifted_horizontal = direction.imag.contains(0) and center.imag.lower() > arb("0.9")
            central = False
            if is_lifted_horizontal:
                disk_multiplier = arb("1.2")
            else:
                q0 = center / self.beta
                tanh0 = -self.i * (self.i * q0).tan()
                exact_b0 = -(self.detuning + self.beta * tanh0)
                canonical_b0 = -(center + self.detuning)
                central = min(abs(exact_b0).upper(), abs(canonical_b0).upper()) < arb("0.002")
                disk_multiplier = arb("1.8" if central else "3.0")
            disk_radius = disk_multiplier * half_width
            z_disk = acb(arb(center.real, disk_radius), arb(center.imag, disk_radius))
            q0 = center / self.beta
            tanh0 = -self.i * (self.i * q0).tan()
            exact_b0 = -(self.detuning + self.beta * tanh0)
            canonical_b0 = -(center + self.detuning)

            q_disk = z_disk / self.beta
            tanh_disk, _ = self.tanh_ball(q_disk)
            exact_b_disk = -(self.detuning + self.beta * tanh_disk)
            canonical_b_disk = -(z_disk + self.detuning)
            use_removable = not (exact_b_disk.abs_lower() > 0 and canonical_b_disk.abs_lower() > 0)
            if use_removable and not central and not is_lifted_horizontal:
                disk_multiplier = arb("1.8")
                disk_radius = disk_multiplier * half_width
                z_disk = acb(arb(center.real, disk_radius), arb(center.imag, disk_radius))
                q_disk = z_disk / self.beta
                tanh_disk, _ = self.tanh_ball(q_disk)
                exact_b_disk = -(self.detuning + self.beta * tanh_disk)
                canonical_b_disk = -(z_disk + self.detuning)
                use_removable = not (exact_b_disk.abs_lower() > 0 and canonical_b_disk.abs_lower() > 0)
            if use_removable:
                require(
                    min(abs(exact_b0).upper(), abs(canonical_b0).upper()) < arb("0.002"),
                    "a noncentral panel entered the removable chart",
                )

            inner_error = self.inner_truncation_bound(z_disk, use_removable)
            maximum = self.true_integrand_bound(z_disk) + inner_error
            ratio = half_width / disk_radius
            outer_tail = (
                2
                * maximum
                * half_width
                * ratio**module.SERIES_ORDER
                / ((module.SERIES_ORDER + 1) * (1 - ratio))
            )
            return outer_tail + 2 * half_width * inner_error

    return AtlasFace


def fraction_arb(value: Fraction) -> arb:
    return arb(value.numerator) / value.denominator


def equal_panels(start: Fraction, end: Fraction, maximum_width: Fraction) -> list[tuple[Fraction, Fraction]]:
    require(end >= start, "reversed real-panel interval")
    if end == start:
        return []
    count = math.ceil((end - start) / maximum_width)
    width = (end - start) / count
    return [(start + (index + Fraction(1, 2)) * width, width / 2) for index in range(count)]


def shifted_real_panels(face) -> tuple[list[tuple[Fraction, Fraction]], Fraction]:
    numerator = round(-float(face.detuning) * CENTER_DENOMINATOR)
    central = Fraction(numerator, CENTER_DENOMINATOR)
    left = central - CENTRAL_HALF_WIDTH
    right = central + CENTRAL_HALF_WIDTH
    radius = Fraction(14, 1)
    require(-radius < left < right < radius, "central event panel escaped the real contour")
    panels = equal_panels(-radius, left, MAX_PANEL_WIDTH)
    panels.append((central, CENTRAL_HALF_WIDTH))
    panels.extend(equal_panels(right, radius, MAX_PANEL_WIDTH))
    offset = abs(face.detuning + fraction_arb(central)).upper()
    require(offset < arb("0.0006"), "rational center did not align with z=-d")
    return panels, central


def shifted_real_core(module, face) -> tuple[acb, Fraction, int]:
    total = acb(0)
    radius = arb(0)
    panels, central = shifted_real_panels(face)
    for center_fraction, half_fraction in panels:
        center = acb(fraction_arb(center_fraction))
        half_width = fraction_arb(half_fraction)
        series = face.difference_series(center, acb(1))
        total += module.integrate_symmetric_series(series, half_width)
        radius += face.panel_error(center, half_width, acb(1))
    return module.add_error(total, radius), central, len(panels)


def complex_record(value: acb) -> dict[str, str]:
    return {
        "real_ball": value.real.str(PRECISION, more=True),
        "imag_ball": value.imag.str(PRECISION, more=True),
    }


def integrate_event(second, first, module, face_class, event_index: int) -> dict[str, Any]:
    mode = mode_for_event(event_index)
    require(4 * mode - C == signed_odd_for_event(event_index), "event roster identity failed")
    face = first.event_face(face_class, module, event_index, mode)
    started = time.perf_counter()
    core, central, real_panel_count = shifted_real_core(module, face)
    right = module.connector(face, 1)
    left = module.connector(face, -1)
    horizontal = module.horizontal_near(face)
    far = second.second_order_far_bound(first, module, face)
    full = module.add_error(core + right - left + horizontal, far)
    physical = arb(2).sqrt() / face.pi * full
    normalized_abs = abs(full)
    physical_abs = abs(physical)
    require(
        normalized_abs < first.TARGET_NORMALIZED,
        f"event {event_index} normalized target failed with {normalized_abs}",
    )
    require(
        physical_abs < first.TARGET_PHYSICAL,
        f"event {event_index} physical target failed with {physical_abs}",
    )
    return {
        "config_version": CONFIG_VERSION,
        "event_index": event_index,
        "mode": mode,
        "signed_odd": 4 * mode - C,
        "detuning_ball": face.detuning.str(PRECISION, more=True),
        "shift_center_rational": f"{central.numerator}/{central.denominator}",
        "real_panel_count": real_panel_count,
        "normalized_second_order_difference_ball": complex_record(full),
        "normalized_second_order_absolute_ball": normalized_abs.str(PRECISION, more=True),
        "normalized_second_order_absolute_upper": normalized_abs.upper().str(PRECISION, more=True),
        "physical_second_order_difference_ball": complex_record(physical),
        "physical_second_order_absolute_ball": physical_abs.str(PRECISION, more=True),
        "physical_second_order_absolute_upper": physical_abs.upper().str(PRECISION, more=True),
        "second_order_far_tail_absolute_bound": far.str(PRECISION, more=True),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }


def load_cache() -> dict[int, dict[str, Any]]:
    if not CACHE.is_file():
        return {}
    rows: dict[int, dict[str, Any]] = {}
    for line_number, line in enumerate(CACHE.read_text(encoding="utf-8").splitlines(), start=1):
        if not line.strip():
            continue
        row = json.loads(line)
        require(row.get("config_version") == CONFIG_VERSION, f"cache config mismatch on line {line_number}")
        index = int(row["event_index"])
        require(index not in rows, f"duplicate cached event {index}")
        rows[index] = row
    return rows


def append_cache(row: dict[str, Any]) -> None:
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    with CACHE.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def parse_events(value: str) -> list[int]:
    events = [int(part.strip()) for part in value.split(",") if part.strip()]
    require(events and all(0 <= index < EVENT_COUNT for index in events), "invalid event selection")
    return events


def upper(row: dict[str, Any], key: str) -> arb:
    return arb(row[key])


def render_note(artifact: dict[str, Any]) -> str:
    c = artifact["certificate"]
    return f"""# All-event second-order normal-form atlas

Date: 2026-08-12

Status: rigorous beta^-4 residual certificate at all 399 exact event centers;
not a proof of continuous-detuning or continuous-height propagation

The real contour is unchanged: it runs from `-14` to `14`, followed by the
same vertical connectors and lifted horizontal pieces through `|Re z|=21`.
For each exact event detuning `d_n`, the real partition is recentered with one
rational panel within `0.0006` of `z=-d_n`.  Thus only the panel containing the
removable normal point uses the full inner exponential recurrence.  Every
neighboring panel keeps `exp(i b y)` exact.  This is the concrete partition
form of the shift `w=z+d_n`; no contour or integral is changed.

All `{c['event_count']}` exact events were enclosed.  The largest normalized
upper bound is `{c['maximum_normalized_absolute_upper']}` at event
`{c['maximum_normalized_event']}`.  The largest physical upper bound is
`{c['maximum_physical_absolute_upper']}`.  The normalized target is
`{c['normalized_target']}` and the physical target is `{c['physical_target']}`.

This proves the beta^-4 comparison at the discrete event centers only.  It
does not prove the residual on a continuous detuning interval or between event
heights.  It does not complete `T_upper`, prove `Lambda<=0`, RH, or a
prize-level conclusion.
"""


def finalize(rows: list[dict[str, Any]], priority: str, workers: int, elapsed: float) -> None:
    require(len(rows) == EVENT_COUNT, "cannot finalize an incomplete all-event atlas")
    require([row["event_index"] for row in rows] == list(range(EVENT_COUNT)), "event cache is not complete and ordered")
    normalized_max = max(rows, key=lambda row: float(upper(row, "normalized_second_order_absolute_upper")))
    physical_max = max(rows, key=lambda row: float(upper(row, "physical_second_order_absolute_upper")))
    artifact: dict[str, Any] = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "beta_minus_4_second_order_normal_form_certified_at_all_399_exact_event_centers",
        "resource_policy": {"workers": workers, "process_priority": priority, "resumable_append_only_cache": True},
        "partition": {
            "coordinate": "w=z+d_n",
            "real_contour_unchanged": True,
            "rational_center_denominator": CENTER_DENOMINATOR,
            "maximum_panel_width": str(MAX_PANEL_WIDTH),
            "central_half_width": str(CENTRAL_HALF_WIDTH),
            "removable_switch": REMOVABLE_SWITCH.str(PRECISION, more=True),
            "outer_series_order": 47,
            "inner_removable_order": 76,
            "quadratic_order_off_removable_chart": 24,
        },
        "certificate": {
            "event_count": EVENT_COUNT,
            "first_event": rows[0],
            "last_event": rows[-1],
            "maximum_normalized_event": normalized_max["event_index"],
            "maximum_normalized_absolute_upper": normalized_max["normalized_second_order_absolute_upper"],
            "maximum_physical_event": physical_max["event_index"],
            "maximum_physical_absolute_upper": physical_max["physical_second_order_absolute_upper"],
            "normalized_target": "0.000019",
            "physical_target": "0.0000086",
            "all_rows": rows,
        },
        "claims": {
            "all_399_exact_event_centers_enclosed": True,
            "all_399_event_center_targets_proved": True,
            "continuous_detuning_interval_proved": False,
            "continuous_height_for_all_events_proved": False,
            "complete_T_upper_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": "All 399 exact event centers only. No continuous-detuning enclosure, all-event continuous-height theorem, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "next_action": "Use the exact exp(i*h*z/beta) factor and corrected contour moments to certify each |h|<=pi/16 event cell, beginning with the event attaining the all-center maximum.",
        "dependencies": {
            "extreme_second_order_gate": {"path": relative(SECOND_RESULT), "sha256": file_hash(SECOND_RESULT)},
        },
        "artifacts": {
            "builder": {"path": relative(BUILDER), "sha256": file_hash(BUILDER)},
            "checker": {"path": relative(CHECKER), "sha256": file_hash(CHECKER)},
            "row_cache": {"path": relative(CACHE), "sha256": file_hash(CACHE)},
        },
        "runtime": {"elapsed_seconds_this_invocation": round(elapsed, 3)},
    }
    NOTE.write_text(render_note(artifact), encoding="utf-8", newline="\n")
    artifact["artifacts"]["note"] = {"path": relative(NOTE), "sha256": file_hash(NOTE)}
    RESULT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")


def run_events(indices: Iterable[int], second, first, module, face_class) -> list[dict[str, Any]]:
    rows = []
    for index in indices:
        row = integrate_event(second, first, module, face_class, index)
        rows.append(row)
        print(
            f"event {index:03d} |D2|<={row['normalized_second_order_absolute_upper']} "
            f"panels={row['real_panel_count']} elapsed={row['elapsed_seconds']}s",
            flush=True,
        )
    return rows


_WORKER_STATE: tuple[Any, Any, Any, Any] | None = None


def initialize_worker() -> None:
    global _WORKER_STATE
    second = load_module(SECOND_BUILDER, f"second_order_worker_{os.getpid()}")
    first = second.load_first()
    first.set_low_priority()
    flint.ctx.dps = PRECISION
    module = first.load_full_module()
    configure(module)
    face_class = make_atlas_face(second, first, module)
    _WORKER_STATE = second, first, module, face_class


def compute_worker(event_index: int) -> dict[str, Any]:
    require(_WORKER_STATE is not None, "worker was not initialized")
    second, first, module, face_class = _WORKER_STATE
    return integrate_event(second, first, module, face_class, event_index)


def run_parallel_atlas(
    pending: list[int],
    cached: dict[int, dict[str, Any]],
    workers: int,
    started: float,
    runtime_limit: float,
) -> bool:
    try:
        import psutil
    except Exception:  # pragma: no cover
        psutil = None
    high_cpu_samples = 0
    parked = False
    with ProcessPoolExecutor(max_workers=workers, initializer=initialize_worker) as executor:
        for offset in range(0, len(pending), workers):
            if time.perf_counter() - started >= runtime_limit:
                parked = True
                break
            indices = pending[offset : offset + workers]
            futures = [executor.submit(compute_worker, index) for index in indices]
            incomplete = set(futures)
            while incomplete:
                done, incomplete = wait(incomplete, timeout=5)
                if incomplete and psutil is not None:
                    total_cpu = psutil.cpu_percent(interval=0.2)
                    high_cpu_samples = high_cpu_samples + 1 if total_cpu > 75 else 0
                    if high_cpu_samples >= 2:
                        parked = True
            rows = sorted((future.result() for future in futures), key=lambda row: row["event_index"])
            for row in rows:
                append_cache(row)
                cached[row["event_index"]] = row
                print(
                    f"event {row['event_index']:03d} |D2|<={row['normalized_second_order_absolute_upper']} "
                    f"panels={row['real_panel_count']} elapsed={row['elapsed_seconds']}s",
                    flush=True,
                )
            if parked:
                break
    return parked


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", help="comma-separated scout events; do not touch the canonical cache")
    parser.add_argument("--reset-cache", action="store_true")
    parser.add_argument("--runtime-limit-seconds", type=float, default=10_800)
    parser.add_argument("--workers", type=int, choices=(1, 2), default=1)
    args = parser.parse_args()

    require(SECOND_BUILDER.is_file() and SECOND_RESULT.is_file() and CHECKER.is_file(), "missing dependency or checker")
    if args.events:
        second = load_module(SECOND_BUILDER, "second_order_source")
        first = second.load_first()
        first.set_low_priority()
        flint.ctx.dps = PRECISION
        module = first.load_full_module()
        configure(module)
        face_class = make_atlas_face(second, first, module)
        run_events(parse_events(args.events), second, first, module, face_class)
        return

    second = load_module(SECOND_BUILDER, "second_order_controller")
    first = second.load_first()
    priority = first.set_low_priority()
    if args.reset_cache and CACHE.exists():
        CACHE.unlink()
    cached = load_cache()
    started = time.perf_counter()
    pending = [index for index in range(EVENT_COUNT) if index not in cached]
    if args.workers == 1:
        flint.ctx.dps = PRECISION
        module = first.load_full_module()
        configure(module)
        face_class = make_atlas_face(second, first, module)
        for index in pending:
            if time.perf_counter() - started >= args.runtime_limit_seconds:
                break
            row = integrate_event(second, first, module, face_class, index)
            append_cache(row)
            cached[index] = row
            print(
                f"event {index:03d} |D2|<={row['normalized_second_order_absolute_upper']} "
                f"panels={row['real_panel_count']} elapsed={row['elapsed_seconds']}s",
                flush=True,
            )
    else:
        run_parallel_atlas(pending, cached, args.workers, started, args.runtime_limit_seconds)

    rows = [cached[index] for index in sorted(cached)]
    if len(rows) == EVENT_COUNT:
        finalize(rows, priority, args.workers, time.perf_counter() - started)
        print(f"certified all {EVENT_COUNT} exact event centers", flush=True)
    else:
        print(
            f"checkpointed {len(rows)}/{EVENT_COUNT} events in {relative(CACHE)}; rerun the same command to resume",
            flush=True,
        )


if __name__ == "__main__":
    main()
