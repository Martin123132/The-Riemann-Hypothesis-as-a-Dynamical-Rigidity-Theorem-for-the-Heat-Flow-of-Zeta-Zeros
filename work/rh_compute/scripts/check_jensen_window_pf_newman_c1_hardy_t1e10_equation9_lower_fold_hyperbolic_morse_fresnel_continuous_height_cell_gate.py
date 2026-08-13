#!/usr/bin/env python3
"""Independently reconstruct the continuous-height prototype cell gate."""

from __future__ import annotations

import importlib.util
import json
import math
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

import flint
from flint import acb, acb_series, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_continuous_height_cell_gate"
FULL_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
FULL_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{FULL_STEM}.py"
MOMENT_ORDER = 14


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


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
    spec = importlib.util.spec_from_file_location("full_line_continuous_checker", FULL_BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load full-line builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def panels(module, face):
    for index in range(module.REAL_PANELS):
        center = -module.CONTOUR_RADIUS + module.REAL_HALF_WIDTH + 2 * module.REAL_HALF_WIDTH * index
        yield acb(center), module.REAL_HALF_WIDTH, acb(1), 1
    for side, orientation in ((1, 1), (-1, -1)):
        for index in range(module.VERTICAL_PANELS):
            center_y = module.VERTICAL_HALF_WIDTH + 2 * module.VERTICAL_HALF_WIDTH * index
            yield acb(side * module.CONTOUR_RADIUS, center_y), module.VERTICAL_HALF_WIDTH, face.i, orientation
    for side in (-1, 1):
        start = -module.HORIZONTAL_END if side < 0 else module.CONTOUR_RADIUS
        for index in range(module.HORIZONTAL_PANELS_PER_SIDE):
            center_x = start + module.REAL_HALF_WIDTH + 2 * module.REAL_HALF_WIDTH * index
            yield acb(center_x, module.CONTOUR_HEIGHT), module.REAL_HALF_WIDTH, acb(1), 1


def segment_ball(center: acb, half_width: arb, direction: acb) -> acb:
    if direction.imag.contains(0):
        return acb(arb(center.real, half_width), center.imag)
    return acb(center.real, arb(center.imag, half_width))


def panel_error(module, face, center: acb, half_width: arb, direction: acb, degree: int) -> arb:
    lifted = direction.imag.contains(0) and center.imag.lower() > arb("0.9")
    disk_radius = arb("1.2" if lifted else "1.8") * half_width
    z_disk = acb(arb(center.real, disk_radius), arb(center.imag, disk_radius))
    q0 = center / face.beta
    tanh0 = -face.i * (face.i * q0).tan()
    exact_b0 = -(face.detuning + face.beta * tanh0)
    canonical_b0 = -(center + face.detuning)
    removable = min(abs(exact_b0).upper(), abs(canonical_b0).upper()) < arb("0.6")
    inner_error = face.inner_truncation_bound(z_disk, removable)
    maximum = (face.true_integrand_bound(z_disk) + inner_error) * abs(z_disk).upper() ** degree
    ratio = half_width / disk_radius
    outer = 2 * maximum * half_width * ratio**module.SERIES_ORDER / ((module.SERIES_ORDER + 1) * (1 - ratio))
    inner = 2 * half_width * abs(segment_ball(center, half_width, direction)).upper() ** degree * inner_error
    return outer + inner


def reconstruct(module, face) -> tuple[list[acb], arb]:
    moments = [acb(0) for _ in range(MOMENT_ORDER + 1)]
    errors = [arb(0) for _ in range(MOMENT_ORDER + 1)]
    mass = arb(0)
    for center, half_width, direction, orientation in panels(module, face):
        base = face.difference_series(center, direction)
        variable = acb_series([center, direction], prec=module.SERIES_ORDER)
        power = acb_series([1], prec=module.SERIES_ORDER)
        for degree in range(MOMENT_ORDER + 1):
            moments[degree] += orientation * module.integrate_symmetric_series(base * power, half_width)
            errors[degree] += panel_error(module, face, center, half_width, direction, degree)
            power *= variable
        mass += 2 * half_width * face.true_integrand_bound(segment_ball(center, half_width, direction))
    return [module.add_error(value, error) for value, error in zip(moments, errors)], mass


def main() -> None:
    flint.ctx.dps = 95
    priority = set_low_priority()
    require(RESULT.is_file(), "missing continuous-height result")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("status") == "one_mode_continuous_height_cell_transport_certified", "bad status")
    require(artifact["claims"]["one_mode_full_height_cell_uniformly_certified"] is True, "continuous claim missing")
    require(artifact["claims"]["all_399_event_cells_proved"] is False, "399-event boundary was promoted")
    require(artifact["claims"]["RH_proved"] is False, "RH boundary was promoted")

    module = load_full_module()
    module.SERIES_ORDER = 200
    module.INNER_ORDER = 240
    module.QUADRATIC_ORDER = 46
    module.REAL_PANELS = 72
    module.REAL_HALF_WIDTH = arb("0.125")
    module.VERTICAL_PANELS = 8
    module.VERTICAL_HALF_WIDTH = arb("0.0625")
    module.HORIZONTAL_PANELS_PER_SIDE = 24
    flint.ctx.cap = module.SERIES_ORDER + 8
    face = module.FullLineFace("event", 0)
    moments, mass = reconstruct(module, face)
    kappa = face.pi / (16 * face.beta)
    coefficients = [(face.i * kappa) ** degree / arb(math.factorial(degree)) * moment for degree, moment in enumerate(moments)]

    stored = artifact["certification"]["moments"]
    for degree, row in enumerate(stored):
        recorded_moment = parse_complex(row["moment_ball"])
        require(moments[degree].real.overlaps(recorded_moment.real), f"moment {degree} real mismatch")
        require(moments[degree].imag.overlaps(recorded_moment.imag), f"moment {degree} imaginary mismatch")

    z_max = (module.HORIZONTAL_END**2 + module.CONTOUR_HEIGHT**2).sqrt()
    size = kappa * z_max
    remainder = mass * size.exp() * size ** (MOMENT_ORDER + 1) / arb(math.factorial(MOMENT_ORDER + 1))
    far = kappa.exp() * module.horizontal_far_bound(face)
    bound = sum((abs(value).upper() for value in coefficients), arb(0)) + remainder + far
    physical = arb(2).sqrt() / face.pi * bound
    require(bound < arb("0.000019"), "independent continuous normalized target failed")
    require(physical < arb("0.0000086"), "independent continuous physical target failed")

    note = (REPO_ROOT / artifact["artifacts"]["note"]["path"]).read_text(encoding="utf-8")
    for token in ("exact transport identity", "every real `t`", "all 399 event cells"):
        require(token in note, f"missing note token: {token}")
    print(
        f"validated continuous height cell independently: |Delta I|<={bound}; "
        f"physical<={physical}; priority={priority}",
        flush=True,
    )


if __name__ == "__main__":
    main()
