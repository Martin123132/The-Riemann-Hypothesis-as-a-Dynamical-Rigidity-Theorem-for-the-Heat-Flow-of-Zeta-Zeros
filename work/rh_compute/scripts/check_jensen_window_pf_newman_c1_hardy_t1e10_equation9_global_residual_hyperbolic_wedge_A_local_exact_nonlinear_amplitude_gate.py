#!/usr/bin/env python3
"""Independently replay the compact exact-minus-affine A-amplitude gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_local_exact_nonlinear_amplitude_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER_PRECISION = 70
CHECKER_SERIES_TERMS = 28
CHECKER_TOLERANCE = "3e-16"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_builder_module():
    spec = importlib.util.spec_from_file_location("A_local_nonlinear_builder", BUILDER)
    require(spec is not None and spec.loader is not None, "cannot load builder module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def record_to_acb(module, record: dict[str, str]):
    return module.acb(module.arb(record["real_ball"]), module.arb(record["imag_ball"]))


def main() -> int:
    artifact = load_json(RESULT)
    require(artifact.get("passed") is True and artifact.get("kind") == STEM, "artifact identity drift")
    decisions = artifact["decision"]
    for key in (
        "exact_endpoint_current_normalization_proved",
        "raw_interval_erfc_disk_cancellation_avoided_by_exact_ODE_enclosure",
        "compact_exact_minus_affine_transformed_amplitude_certified",
        "all_84_local_nonlinear_modes_certified",
        "signed_local_nonlinear_sum_certified_before_norms",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    require(decisions.get("exact_exterior_recomputed_here") is False, "exterior scope drift")
    for key in ("complete_A_endpoint_block_proved", "R_Dir_bound_proved", "rh_implication"):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    priority = module.set_low_priority()
    require(priority == "below_normal" or os.name != "nt", f"checker priority drift: {priority}")
    module.ctx.dps = CHECKER_PRECISION
    module.ctx.threads = 1
    module.symbolic_normalization_certificate()
    saved_rows = {row["mode"]: row for row in artifact["certificate"]["rows"]}
    require(set(saved_rows) == set(range(39_853, 39_937)), "saved mode roster drift")
    fresh_rows: list[dict[str, Any]] = []
    for index, mode in enumerate(range(39_853, 39_937), start=1):
        fresh = module.evaluate_mode(mode, CHECKER_PRECISION, CHECKER_SERIES_TERMS, CHECKER_TOLERANCE)
        fresh_rows.append(fresh)
        saved = saved_rows[mode]
        require(
            module.arb(fresh["physical_local_exact_minus_affine_ball"]).overlaps(
                module.arb(saved["physical_local_exact_minus_affine_ball"])
            ),
            f"higher-precision physical local nonlinear value misses at mode {mode}",
        )
        for field in (
            "canonical_local_exact_minus_affine_ball",
            "raw_A_endpoint_local_exact_minus_affine_ball",
        ):
            require(
                record_to_acb(module, fresh[field]).overlaps(record_to_acb(module, saved[field])),
                f"higher-precision {field} misses at mode {mode}",
            )
        if index % 8 == 0:
            print(f"checked local exact nonlinear rows {index}/84", flush=True)

    fresh_sums = module.aggregate_rows(fresh_rows)
    saved_sums = artifact["certificate"]["signed_sums"]
    for field, fresh in fresh_sums.items():
        saved = saved_sums[field]
        if isinstance(fresh, dict):
            require(record_to_acb(module, fresh).overlaps(record_to_acb(module, saved)), f"aggregate misses {field}")
        else:
            require(module.arb(fresh).overlaps(module.arb(saved)), f"aggregate misses {field}")

    cache_path = REPO_ROOT / artifact["cache"]["path"]
    require(cache_path.is_file(), "missing resumable cache")
    require(file_hash(cache_path) == artifact["cache"]["sha256"], "cache hash drift")
    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("elementary ODE derivative" in note, "stable analytic enclosure explanation missing")
    require("complete A endpoint block" in note, "proof boundary missing")
    print("independently checked compact exact-minus-affine transformed amplitude on A", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
