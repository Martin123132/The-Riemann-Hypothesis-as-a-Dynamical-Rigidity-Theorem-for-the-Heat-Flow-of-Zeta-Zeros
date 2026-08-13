#!/usr/bin/env python3
"""Check the full-line one-mode Morse--Fresnel/Airy join independently."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_hyperbolic_morse_fresnel_full_line_join_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"


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


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def load_builder_module():
    spec = importlib.util.spec_from_file_location("full_line_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def independent_face(module, label: str, offset: int) -> acb:
    # Different, narrower partition and larger outer/inner truncations from
    # the builder.  Wider half-unit disks cross the cubic exponential's useful
    # Taylor margin at this order, so the independent stress goes narrower.
    module.SERIES_ORDER = 200
    module.INNER_ORDER = 240
    module.QUADRATIC_ORDER = 46
    module.REAL_PANELS = 72
    module.REAL_HALF_WIDTH = arb("0.125")
    module.VERTICAL_PANELS = 8
    module.VERTICAL_HALF_WIDTH = arb("0.0625")
    module.HORIZONTAL_PANELS_PER_SIDE = 24
    flint.ctx.cap = module.SERIES_ORDER + 8
    face = module.FullLineFace(label, offset)
    center = module.real_core(face)
    right = module.connector(face, 1)
    left = module.connector(face, -1)
    horizontal = module.horizontal_near(face)
    far_tail = module.horizontal_far_bound(face)
    return module.add_error(center + right - left + horizontal, far_tail)


def main() -> None:
    flint.ctx.dps = 90
    priority = set_low_priority()
    require(RESULT.is_file(), "missing full-line result")
    require(BUILDER.is_file(), "missing full-line builder")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("status") == "full_line_one_mode_hyperbolic_Morse_Fresnel_Airy_join_certified", "bad status")
    require(artifact["claims"]["prototype_one_mode_join_proved"] is True, "prototype claim missing")
    require(artifact["claims"]["all_399_events_proved"] is False, "399-event boundary was promoted")
    require(artifact["claims"]["RH_proved"] is False, "RH boundary was promoted")
    stored = {row["label"]: row for row in artifact["certified_faces"]}
    module = load_builder_module()
    maximum = arb(0)
    for label, offset in (("lower_face", -1), ("event", 0), ("upper_face", 1)):
        value = independent_face(module, label, offset)
        recorded = parse_complex(stored[label]["normalized_full_line_difference_ball"])
        require(value.real.overlaps(recorded.real), f"{label} real values do not overlap")
        require(value.imag.overlaps(recorded.imag), f"{label} imaginary values do not overlap")
        require(abs(value) < arb("0.000019"), f"{label} independent normalized target failed")
        require(arb(2).sqrt() / arb.pi() * abs(value) < arb("0.0000086"), f"{label} independent physical target failed")
        maximum = max(maximum, abs(value).upper())
    note = (REPO_ROOT / artifact["artifacts"]["note"]["path"]).read_text(encoding="utf-8")
    for token in (
        "full-line one-mode prototype join certified",
        "|Delta I(t)| < 0.000019",
        "all 399 turning events",
    ):
        require(token in note, f"missing note token: {token}")
    print(f"validated full-line prototype join independently: faces=3, max |difference|={maximum}; priority={priority}", flush=True)


if __name__ == "__main__":
    main()
