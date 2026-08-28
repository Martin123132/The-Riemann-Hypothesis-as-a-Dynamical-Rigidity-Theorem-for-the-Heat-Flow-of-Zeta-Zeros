#!/usr/bin/env python3
"""Independently refine the A-exterior endpoint-tail remainder bound."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_endpoint_tail_remainder_gate"
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
    spec = importlib.util.spec_from_file_location("A_endpoint_tail_builder", BUILDER)
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
        "four_term_endpoint_tail_used",
        "zero_endpoint_power_vanishing_retained",
        "all_84_mode_remainders_certified",
        "complete_physical_replacement_error_below_1p5e_minus_14",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "rational_common_phase_integral_value_proved",
        "outer_IBP_remainder_bound_proved",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    priority = module.set_low_priority()
    require(priority == "below_normal" or os.name != "nt", f"checker priority drift: {priority}")
    fresh = module.certificate(slabs=8192, precision=70)
    saved = artifact["certificate"]
    require(fresh["mode_count"] == saved["mode_count"] == 84, "mode count drift")
    for key in (
        "max_phase_stripped_endpoint_tail_remainder_ball",
        "canonical_complete_exterior_replacement_error_ball",
        "physical_complete_exterior_replacement_error_ball",
    ):
        fresh_ball = module.arb(fresh[key])
        saved_ball = module.arb(saved[key])
        require(fresh_ball.upper() < saved_ball.upper(), f"refinement did not tighten {key}")
    require(
        module.arb(fresh["physical_complete_exterior_replacement_error_ball"]).upper()
        < module.arb("1e-14"),
        "refined physical error misses 1e-14",
    )

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("integral_0^x |b_3'(u)|" in note, "remainder identity missing")
    require("No value of the resulting" in note, "proof boundary missing")
    print("independently refined four-term exact-endpoint tail remainder on A exterior", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
