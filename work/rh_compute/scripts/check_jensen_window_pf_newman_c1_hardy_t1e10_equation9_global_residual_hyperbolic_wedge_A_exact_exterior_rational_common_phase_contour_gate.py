#!/usr/bin/env python3
"""Independently refine the rational common-phase A-exterior contour."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exact_exterior_rational_common_phase_contour_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module():
    spec = importlib.util.spec_from_file_location("a_exact_exterior_rational_contour_gate", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def complex_ball(module, record: dict[str, str]):
    return module.acb(module.arb(record["real_ball"]), module.arb(record["imag_ball"]))


def max_radius(value) -> object:
    return max(value.real.rad(), value.imag.rad())


def main() -> int:
    require(RESULT.is_file(), "missing production result")
    saved_artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(saved_artifact.get("passed") is True, "production artifact is not passed")
    decisions = saved_artifact["decision"]
    for key in (
        "logistic_contour_strip_analytic_for_all_84_modes",
        "grouped_compact_common_phase_integral_certified",
        "finite_vertical_contour_leg_certified",
        "horizontal_ray_below_2e_minus_31_canonical",
        "four_term_rational_exterior_value_certified",
        "full_exact_exterior_value_after_endpoint_replacement_certified",
    ):
        require(decisions[key] is True, f"missing positive decision: {key}")
    for key in (
        "compact_exact_minus_affine_transformed_amplitude_certified",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decisions[key] is False, f"proof boundary promoted: {key}")

    module = load_module()
    module.set_low_priority()
    fresh = module.certificate(
        precision=70,
        series_order=42,
        compact_max_width_text="0.000015",
        vertical_max_width_text="0.0000075",
    )
    saved = saved_artifact["certificate"]

    for key in (
        "compact_relative_integral_ball",
        "relative_rational_integral_ball",
        "canonical_rational_integral_ball",
    ):
        fresh_ball = complex_ball(module, fresh[key])
        saved_ball = complex_ball(module, saved[key])
        require(fresh_ball.overlaps(saved_ball), f"fresh complex enclosure misses saved {key}")
        require(max_radius(fresh_ball) < max_radius(saved_ball), f"fresh complex enclosure is not tighter: {key}")

    for key in (
        "physical_rational_integral_ball",
        "physical_full_exact_exterior_ball",
    ):
        fresh_ball = module.arb(fresh[key])
        saved_ball = module.arb(saved[key])
        require(fresh_ball.overlaps(saved_ball), f"fresh real enclosure misses saved {key}")
        require(fresh_ball.rad() < saved_ball.rad(), f"fresh real enclosure is not tighter: {key}")

    require(fresh["compact_panel_count"] > saved["compact_panel_count"], "checker did not refine compact partition")
    require(fresh["vertical_leg"]["panel_count"] > saved["vertical_leg"]["panel_count"], "checker did not refine vertical partition")
    require(
        module.arb(fresh["horizontal_ray"]["relative_horizontal_integral_bound_ball"]).upper()
        < module.arb("2e-31"),
        "fresh horizontal ray misses cap",
    )
    require(
        module.arb(fresh["physical_full_exact_exterior_ball"]).upper() < module.arb("-0.0023"),
        "fresh full exact exterior sign guard lost",
    )

    for record in saved_artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in saved_artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")

    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("finite vertical leg" in note, "contour identity missing from note")
    require("No compact exact-minus-affine" in note, "proof boundary missing from note")
    print("independently refined rational common-phase contour and full exact A exterior", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
