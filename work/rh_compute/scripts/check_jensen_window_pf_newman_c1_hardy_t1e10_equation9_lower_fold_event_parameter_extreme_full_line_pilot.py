#!/usr/bin/env python3
"""Independently reconstruct both extreme-detuning full-line pilots."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

import flint
from flint import acb, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_full_line_pilot"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
FULL_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate"
FULL_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{FULL_STEM}.py"
C = 159_577


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
    spec = importlib.util.spec_from_file_location("extreme_full_line_checker_source", FULL_BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load full-line builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure(module) -> None:
    module.SERIES_ORDER = 200
    module.INNER_ORDER = 240
    module.QUADRATIC_ORDER = 46
    module.CONTOUR_RADIUS = arb(14)
    module.REAL_PANELS = 280
    module.REAL_HALF_WIDTH = arb("0.05")
    module.VERTICAL_PANELS = 20
    module.VERTICAL_HALF_WIDTH = arb("0.025")
    module.HORIZONTAL_END = arb(21)
    module.HORIZONTAL_PANELS_PER_SIDE = 70
    flint.ctx.cap = module.SERIES_ORDER + 8


def event_face(module, event_index: int, mode: int):
    face = module.FullLineFace(f"checker_event_{event_index}_mode_{mode}", 0)
    face.m = arb(mode)
    face.tau = face.pi * face.m * (face.c - 2 * face.m)
    face.t = face.tau
    face.eta = face.tstar / face.t
    face.lam = (face.eta - 1) * face.t / face.beta
    face.detuning = 4 * face.beta * (face.m - face.c / 4) / face.c
    return face


def reconstruct(module, face) -> acb:
    core = module.real_core(face)
    right = module.connector(face, 1)
    left = module.connector(face, -1)
    horizontal = module.horizontal_near(face)
    far = module.horizontal_far_bound(face)
    return module.add_error(core + right - left + horizontal, far)


def main() -> None:
    flint.ctx.dps = 90
    priority = set_low_priority()
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("status") == "extreme_event_full_line_signed_pilot_complete", "bad status")
    require([(row["event_index"], row["mode"], row["signed_odd"]) for row in artifact["pilots"]] == [(398, 39695, -797), (397, 40093, 795)], "wrong extreme roster")
    module = load_full_module()
    configure(module)
    maximum = arb(0)
    for row in artifact["pilots"]:
        value = acb(arb(row["normalized_difference_ball"]["real_ball"]), arb(row["normalized_difference_ball"]["imag_ball"]))
        recorded = arb(row["normalized_absolute_ball"])
        require(abs(value).overlaps(recorded), f"absolute ball mismatch at event {row['event_index']}")
        require(arb(row["far_tail_absolute_bound"]) < arb("1e-6"), f"far tail too large at event {row['event_index']}")
        require(bool(row["prototype_normalized_target_proved"]) == (recorded < arb("0.000019")), "normalized decision mismatch")
        require(bool(row["prototype_normalized_target_disproved"]) == (recorded.lower() > arb("0.000019")), "normalized exclusion mismatch")
        independent = reconstruct(module, event_face(module, row["event_index"], row["mode"]))
        require(independent.real.overlaps(value.real), f"independent real mismatch at event {row['event_index']}")
        require(independent.imag.overlaps(value.imag), f"independent imaginary mismatch at event {row['event_index']}")
        require(abs(independent).lower() > arb("0.000019"), f"independent target exclusion failed at event {row['event_index']}")
        maximum = max(maximum, abs(independent).upper())
    require(artifact["decision"]["uniform_399_event_remainder_proved"] is False, "uniform boundary promoted")
    require(artifact["decision"]["RH_proved"] is False, "RH boundary promoted")
    note = (REPO_ROOT / artifact["artifacts"]["note"]["path"]).read_text(encoding="utf-8")
    for token in ("two-endpoint diagnostic", "quantitative propagation", "does not prove"):
        require(token in note, f"missing note token: {token}")
    print(f"validated extreme full-line pilots independently: max |D|={maximum}; priority={priority}")


if __name__ == "__main__":
    main()
