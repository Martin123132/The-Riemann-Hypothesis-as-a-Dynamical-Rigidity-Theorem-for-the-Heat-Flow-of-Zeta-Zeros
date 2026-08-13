#!/usr/bin/env python3
"""Test the signed full-line join at both extreme event detunings."""

from __future__ import annotations

import hashlib
import importlib.util
import json
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
from flint import acb, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_full_line_pilot"
FULL_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate"
GEOMETRY_STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_contour_geometry_gate"
FULL_BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{FULL_STEM}.py"
GEOMETRY = REPO_ROOT / f"work/rh_compute/results/{GEOMETRY_STEM}.json"
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
    spec = importlib.util.spec_from_file_location("extreme_full_line_source", FULL_BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load full-line builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def configure(module) -> None:
    module.CONTOUR_RADIUS = arb(14)
    module.REAL_PANELS = 224
    module.REAL_HALF_WIDTH = arb("0.0625")
    module.VERTICAL_PANELS = 16
    module.VERTICAL_HALF_WIDTH = arb("0.03125")
    module.HORIZONTAL_END = arb(21)
    module.HORIZONTAL_PANELS_PER_SIDE = 56
    flint.ctx.cap = module.SERIES_ORDER + 8


def event_face(module, event_index: int, mode: int):
    face = module.FullLineFace(f"event_{event_index}_mode_{mode}", 0)
    face.m = arb(mode)
    face.tau = face.pi * face.m * (face.c - 2 * face.m)
    face.t = face.tau
    face.eta = face.tstar / face.t
    face.lam = (face.eta - 1) * face.t / face.beta
    face.detuning = 4 * face.beta * (face.m - face.c / 4) / face.c
    return face


def complex_record(value: acb) -> dict[str, str]:
    return {"real_ball": value.real.str(PRECISION, more=True), "imag_ball": value.imag.str(PRECISION, more=True)}


def integrate_face(module, face, event_index: int, mode: int) -> dict[str, Any]:
    core = module.real_core(face)
    right = module.connector(face, 1)
    left = module.connector(face, -1)
    horizontal = module.horizontal_near(face)
    far = module.horizontal_far_bound(face)
    finite = core + right - left + horizontal
    full = module.add_error(finite, far)
    physical = arb(2).sqrt() / face.pi * full
    normalized_abs = abs(full)
    physical_abs = abs(physical)
    return {
        "event_index": event_index,
        "mode": mode,
        "signed_odd": 4 * mode - C,
        "detuning_ball": face.detuning.str(PRECISION, more=True),
        "lambda_event_ball": face.lam.str(PRECISION, more=True),
        "real_core_ball": complex_record(core),
        "right_connector_ball": complex_record(right),
        "left_connector_ball": complex_record(left),
        "horizontal_near_ball": complex_record(horizontal),
        "far_tail_absolute_bound": far.str(PRECISION, more=True),
        "normalized_difference_ball": complex_record(full),
        "normalized_absolute_ball": normalized_abs.str(PRECISION, more=True),
        "normalized_absolute_lower": normalized_abs.lower().str(PRECISION, more=True),
        "normalized_absolute_upper": normalized_abs.upper().str(PRECISION, more=True),
        "physical_difference_ball": complex_record(physical),
        "physical_absolute_ball": physical_abs.str(PRECISION, more=True),
        "physical_absolute_lower": physical_abs.lower().str(PRECISION, more=True),
        "physical_absolute_upper": physical_abs.upper().str(PRECISION, more=True),
        "prototype_normalized_target_proved": normalized_abs < TARGET_NORMALIZED,
        "prototype_normalized_target_disproved": normalized_abs.lower() > TARGET_NORMALIZED,
        "prototype_physical_target_proved": physical_abs < TARGET_PHYSICAL,
        "prototype_physical_target_disproved": physical_abs.lower() > TARGET_PHYSICAL,
    }


def render_note(artifact: dict[str, Any]) -> str:
    rows = "\n".join(
        f"event {row['event_index']} (mode {row['mode']}, 4m-C={row['signed_odd']}): "
        f"normalized={row['normalized_absolute_ball']}, physical={row['physical_absolute_ball']}"
        for row in artifact["pilots"]
    )
    decision = artifact["decision"]
    return f"""# Extreme-detuning full-line pilot

Date: 2026-08-12

Status: rigorous two-endpoint diagnostic on the common event contour; not a
proof of the 399-event propagation theorem

The event-parameter geometry gate proves that `R=14`, `Im z=1`, and
`|Re z|<=21` contain all exact and canonical saddles and retain positive
lifted-contour margins throughout the 399-event range.  This pilot now runs
the signed retained-Fresnel/logistic minus Airy integral at the two extreme
event detunings, before absolute values:

```text
{rows}
```

Prototype target survives both extremes: `{decision['prototype_target_survives_both_extremes']}`.
Prototype target is rigorously excluded at both extremes:
`{decision['prototype_target_disproved_at_both_extremes']}`.

This diagnostic tests quantitative propagation after contour geometry.
It does not prove a uniform bound between the endpoints, a
paired-event cancellation theorem, all 399 event remainders, complete
`T_upper`, `Lambda<=0`, RH, or a prize-level conclusion.
"""


def main() -> None:
    started = time.perf_counter()
    priority = set_low_priority()
    require(FULL_BUILDER.is_file() and GEOMETRY.is_file() and CHECKER.is_file(), "missing dependency or checker")
    flint.ctx.dps = PRECISION
    module = load_full_module()
    configure(module)
    pilots = [integrate_face(module, event_face(module, index, mode), index, mode) for index, mode in EXTREMES]
    survives = all(row["prototype_normalized_target_proved"] and row["prototype_physical_target_proved"] for row in pilots)
    disproved = all(row["prototype_normalized_target_disproved"] and row["prototype_physical_target_disproved"] for row in pilots)
    artifact: dict[str, Any] = {
        "kind": STEM,
        "date": "2026-08-12",
        "status": "extreme_event_full_line_signed_pilot_complete",
        "resource_policy": {"workers": 1, "process_priority": priority},
        "contour": {"radius": 14, "height": 1, "horizontal_end": 21, "real_panels": 224, "connector_panels_per_side": 16, "horizontal_panels_per_side": 56},
        "pilots": pilots,
        "decision": {
            "extreme_event_balls_rigorously_enclosed": True,
            "prototype_target_survives_both_extremes": survives,
            "prototype_target_disproved_at_both_extremes": disproved,
            "uniform_399_event_remainder_proved": False,
            "paired_event_cancellation_proved": False,
            "RH_proved": False,
        },
        "proof_boundary": "Two extreme event points only. No between-event-parameter enclosure, paired-event cancellation theorem, 399-event propagation, complete T_upper, Lambda<=0, RH, or prize-level conclusion is proved.",
        "next_action": "Use the extreme result to choose either a global event-parameter interval enclosure or a symmetry-preserving paired-event/corrected-normal-form route.",
        "dependencies": {
            "geometry_gate": {"path": relative(GEOMETRY), "sha256": file_hash(GEOMETRY)},
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
        "completed extreme full-line pilots: "
        + ", ".join(f"n={row['event_index']} |D|={row['normalized_absolute_ball']}" for row in pilots),
        flush=True,
    )


if __name__ == "__main__":
    main()
