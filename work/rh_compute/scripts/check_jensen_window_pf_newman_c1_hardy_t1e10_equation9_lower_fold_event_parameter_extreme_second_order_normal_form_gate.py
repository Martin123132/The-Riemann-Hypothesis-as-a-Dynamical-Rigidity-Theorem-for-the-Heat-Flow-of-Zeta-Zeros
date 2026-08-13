#!/usr/bin/env python3
"""Independently reconstruct the second-order extreme-event gate."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))

import flint
from flint import acb, arb


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_second_order_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_builder():
    spec = importlib.util.spec_from_file_location("second_order_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load second-order builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> None:
    flint.ctx.dps = 90
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("status") == "extreme_event_beta_minus_4_second_order_normal_form_certified", "bad status")
    second = load_builder()
    first = second.load_first()
    module = first.load_full_module()
    first.configure(module, checker=True)
    face_class = second.make_second_order_face(first, module)
    maximum = arb(0)
    for row in artifact["certified_extremes"]:
        face = first.event_face(face_class, module, row["event_index"], row["mode"])
        core = module.real_core(face)
        right = module.connector(face, 1)
        left = module.connector(face, -1)
        horizontal = module.horizontal_near(face)
        far = second.second_order_far_bound(first, module, face)
        independent = module.add_error(core + right - left + horizontal, far)
        recorded = parse_complex(row["normalized_second_order_difference_ball"])
        require(independent.real.overlaps(recorded.real), f"real mismatch at event {row['event_index']}")
        require(independent.imag.overlaps(recorded.imag), f"imag mismatch at event {row['event_index']}")
        require(abs(independent) < first.TARGET_NORMALIZED, f"normalized target failed at event {row['event_index']}")
        require(arb(2).sqrt() / face.pi * abs(independent) < first.TARGET_PHYSICAL, f"physical target failed at event {row['event_index']}")
        maximum = max(maximum, abs(independent).upper())
    require(artifact["claims"]["uniform_intermediate_detuning_proved"] is False, "uniform boundary promoted")
    require(artifact["claims"]["RH_proved"] is False, "RH boundary promoted")
    note = (REPO_ROOT / artifact["artifacts"]["note"]["path"]).read_text(encoding="utf-8")
    for token in ("complete beta^-4 comparison", "degree-four", "two-point endpoint result"):
        require(token in note, f"missing note token: {token}")
    print(f"validated second-order extreme normal form independently: max |D_2|={maximum}")


if __name__ == "__main__":
    main()
