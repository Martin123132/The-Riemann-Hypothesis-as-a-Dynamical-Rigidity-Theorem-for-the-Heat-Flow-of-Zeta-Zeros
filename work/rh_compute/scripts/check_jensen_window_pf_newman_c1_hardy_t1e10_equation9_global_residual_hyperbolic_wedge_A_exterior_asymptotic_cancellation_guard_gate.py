#!/usr/bin/env python3
"""Independently check the A-exterior asymptotic cancellation guard."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_exterior_asymptotic_cancellation_guard_gate"
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
    spec = importlib.util.spec_from_file_location("A_exterior_cancellation_builder", BUILDER)
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
        "finite_cutoff_affine_amplitude_decomposition_remains_exact",
        "all_84_mode_grouped_affine_leading_coefficient_nonzero",
        "full_exact_endpoint_current_integrable_at_x_zero",
        "common_cutoff_recombination_required_before_exterior_limit",
        "local_compact_affine_amplitude_split_remains_admissible",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in (
        "standalone_affine_exterior_improper_limit_exists",
        "standalone_exact_minus_affine_exterior_improper_limit_exists",
        "exact_exterior_value_proved",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    priority = module.set_low_priority()
    require(priority == "below_normal" or os.name != "nt", f"checker priority drift: {priority}")
    fresh_symbolic = module.symbolic_certificate()
    require(fresh_symbolic == artifact["symbolic_certificate"], "symbolic certificate drift")
    fresh = module.mode_certificate(precision=120)
    saved = artifact["mode_certificate"]
    require(fresh["mode_count"] == saved["mode_count"] == 84, "mode count drift")
    for key in (
        "mu_S_ball",
        "lambda_minus_mu_real_lower_ball",
        "lambda_minus_mu_real_upper_ball",
        "max_abs_lambda_P_imag_ball",
    ):
        require(module.arb(fresh[key]).overlaps(module.arb(saved[key])), f"higher-precision miss: {key}")
    for component in ("real_ball", "imag_ball", "absolute_ball"):
        require(
            module.arb(fresh["grouped_affine_exterior_leading_coefficient_ball"][component]).overlaps(
                module.arb(saved["grouped_affine_exterior_leading_coefficient_ball"][component])
            ),
            f"grouped coefficient miss: {component}",
        )
    fresh_rows = {row["mode"]: row for row in fresh["rows"]}
    saved_rows = {row["mode"]: row for row in saved["rows"]}
    require(set(fresh_rows) == set(saved_rows) == set(range(39_853, 39_937)), "mode roster drift")
    for mode in fresh_rows:
        require(
            module.arb(fresh_rows[mode]["lambda_minus_mu_real_ball"]).overlaps(
                module.arb(saved_rows[mode]["lambda_minus_mu_real_ball"])
            ),
            f"lambda-mu miss at mode {mode}",
        )
        for component in ("real_ball", "imag_ball", "absolute_ball"):
            require(
                module.arb(fresh_rows[mode]["affine_exterior_leading_coefficient_ball"][component]).overlaps(
                    module.arb(saved_rows[mode]["affine_exterior_leading_coefficient_ball"][component])
                ),
                f"leading coefficient miss at mode {mode}: {component}",
            )

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("standalone affine exterior has no improper limit" in note, "affine nonconvergence guard missing")
    require("with the two terms retained under one cutoff" in note, "common-cutoff guard missing")
    print("independently checked A-exterior common-cutoff asymptotic cancellation guard", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
