#!/usr/bin/env python3
"""Independently check the full exact A-exterior common-phase reduction."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_common_phase_reduction_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_builder_module():
    spec = importlib.util.spec_from_file_location("A_exact_exterior_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "exact_endpoint_driver_erfc_primitive_proved",
        "phase_stripped_endpoint_tail_ODE_proved",
        "full_exact_exterior_common_phase_reduction_proved",
        "mode_dependent_face_phase_cancelled_exactly",
        "ten_mode_incomplete_stationary_roster_certified",
        "endpoint_IBP_hierarchy_evaluated_at_all_cutoffs",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "diagnostic_exterior_scale_is_rigorous",
        "uniform_endpoint_hierarchy_remainder_bound_proved",
        "outer_IBP_remainder_bound_proved",
        "exact_exterior_value_proved",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    priority = module.set_low_priority()
    require(priority == "below_normal" or os.name != "nt", f"checker priority drift: {priority}")
    require(module.symbolic_certificate() == artifact["symbolic_certificate"], "symbolic certificate drift")
    fresh_roster = module.roster_certificate(precision=120)
    saved_roster = artifact["roster_certificate"]
    require(fresh_roster["positive_stationary_endpoint_modes"] == list(range(39_927, 39_937)), "stationary roster drift")
    for key in ("x_star_ball", "w0_min_ball", "w0_max_ball"):
        require(module.arb(fresh_roster[key]).overlaps(module.arb(saved_roster[key])), f"roster miss: {key}")
    for fresh, saved in zip(
        fresh_roster["max_abs_endpoint_IBP_terms_at_cutoff"],
        saved_roster["max_abs_endpoint_IBP_terms_at_cutoff"],
    ):
        require(module.arb(fresh).overlaps(module.arb(saved)), "cutoff hierarchy maximum miss")

    fresh_pilot = module.diagnostic_pilot(dps=70, gauss_nodes=220, outer_terms=6)
    saved_pilot = artifact["diagnostic_pilot"]
    import mpmath as mp

    mp.mp.dps = 70
    fresh_near = mp.mpc(*fresh_pilot["near_relative_complex"])
    saved_near = mp.mpc(*saved_pilot["near_relative_complex"])
    require(abs(fresh_near - saved_near) < mp.mpf("1e-10"), "stationary pilot drift")
    fresh_scale = mp.mpf(fresh_pilot["final_physical_scale_pilot"])
    saved_scale = mp.mpf(saved_pilot["final_physical_scale_pilot"])
    require(abs(fresh_scale - saved_scale) < mp.mpf("5e-10"), "physical scale pilot drift")
    require(abs(fresh_scale) > mp.mpf("0.002"), "diagnostic scale unexpectedly small")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("single common-phase family" in note, "common-phase reduction missing")
    require("The final displayed scale" in note, "diagnostic boundary missing")
    print("independently checked full exact A-exterior common-phase reduction", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
