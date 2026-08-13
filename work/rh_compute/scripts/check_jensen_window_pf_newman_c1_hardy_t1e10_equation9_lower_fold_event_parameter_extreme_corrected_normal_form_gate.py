#!/usr/bin/env python3
"""Independently reconstruct the corrected extreme-event certificates."""

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


STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_event_parameter_extreme_corrected_normal_form_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
BUILDER = REPO_ROOT / f"work/rh_compute/scripts/{STEM}.py"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_builder():
    spec = importlib.util.spec_from_file_location("corrected_extreme_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load corrected builder")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_complex(record: dict[str, str]) -> acb:
    return acb(arb(record["real_ball"]), arb(record["imag_ball"]))


def main() -> None:
    flint.ctx.dps = 90
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(artifact.get("status") == "extreme_event_beta_minus_2_corrected_normal_form_diagnostic_complete", "bad status")
    module = load_builder()
    full = module.load_full_module()
    module.configure(full, checker=True)
    face_class = module.make_corrected_face(full)
    maximum = arb(0)
    for row in artifact["certified_extremes"]:
        face = module.event_face(face_class, full, row["event_index"], row["mode"])
        core = full.real_core(face)
        right = full.connector(face, 1)
        left = full.connector(face, -1)
        horizontal = full.horizontal_near(face)
        far = module.corrected_far_bound(full, face)
        independent = full.add_error(core + right - left + horizontal, far)
        recorded = parse_complex(row["normalized_corrected_difference_ball"])
        require(independent.real.overlaps(recorded.real), f"real mismatch at event {row['event_index']}")
        require(independent.imag.overlaps(recorded.imag), f"imag mismatch at event {row['event_index']}")
        require(abs(independent).lower() > arb("0.000019"), f"normalized target exclusion failed at event {row['event_index']}")
        require((arb(2).sqrt() / face.pi * abs(independent)).lower() > arb("0.0000086"), f"physical target exclusion failed at event {row['event_index']}")
        maximum = max(maximum, abs(independent).upper())
    require(artifact["claims"]["uniform_intermediate_detuning_proved"] is False, "uniform boundary promoted")
    require(artifact["claims"]["RH_proved"] is False, "RH boundary promoted")
    note = (REPO_ROOT / artifact["artifacts"]["note"]["path"]).read_text(encoding="utf-8")
    require(artifact["claims"]["both_extreme_events_below_prototype_targets"] is False, "false target recovery")
    require(artifact["claims"]["both_extreme_events_exceed_prototype_targets"] is True, "target exclusion missing")
    for token in ("beta^-2 corrected diagnostic", "Both rigorously exceed", "does not"):
        require(token in note, f"missing note token: {token}")
    print(f"validated corrected extreme normal form independently: max |D_corr|={maximum}")


if __name__ == "__main__":
    main()
