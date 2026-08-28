#!/usr/bin/env python3
"""Independently replay the finite-box local exact-domain affine A carrier."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_affine_A_local_exact_domain_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER_PRECISION = 70
CHECKER_SERIES_TERMS = 28
CHECKER_TOL = "3e-14"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_builder_module():
    spec = importlib.util.spec_from_file_location("A_local_exact_builder", BUILDER)
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
        "lower_Fresnel_tail_implemented_by_exact_erf_primitive",
        "affine_P_moment_integrated_exactly",
        "local_exact_domain_used_instead_of_global_tangent_plus_local_face",
        "all_84_local_exact_affine_modes_certified",
        "signed_local_exact_affine_sum_certified_before_norms",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    require(decisions.get("tangent_exterior_included_in_local_carrier") is False, "tangent exterior guard drift")
    for key in (
        "exact_exterior_affine_current_bound_proved",
        "transformed_amplitude_remainder_bound_proved",
        "complete_A_endpoint_block_proved",
        "R_Dir_bound_proved",
        "rh_implication",
    ):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    priority = module.set_low_priority()
    require(priority == "below_normal" or os.name != "nt", f"checker priority drift: {priority}")
    module.ctx.dps = CHECKER_PRECISION
    module.ctx.threads = 1
    saved_rows = {row["mode"]: row for row in artifact["certificate"]["rows"]}
    require(set(saved_rows) == set(range(39_853, 39_937)), "saved mode roster drift")
    fresh_rows: list[dict[str, Any]] = []
    for index, mode in enumerate(range(39_853, 39_937), start=1):
        fresh = module.evaluate_mode(mode, CHECKER_PRECISION, CHECKER_SERIES_TERMS, CHECKER_TOL)
        fresh_rows.append(fresh)
        require(
            module.arb(fresh["physical_local_exact_affine_ball"]).overlaps(
                module.arb(saved_rows[mode]["physical_local_exact_affine_ball"])
            ),
            f"higher-precision physical local exact value misses at mode {mode}",
        )
        require(
            module.acb(
                module.arb(fresh["canonical_local_exact_affine_ball"]["real_ball"]),
                module.arb(fresh["canonical_local_exact_affine_ball"]["imag_ball"]),
            ).overlaps(
                module.acb(
                    module.arb(saved_rows[mode]["canonical_local_exact_affine_ball"]["real_ball"]),
                    module.arb(saved_rows[mode]["canonical_local_exact_affine_ball"]["imag_ball"]),
                )
            ),
            f"higher-precision canonical local exact value misses at mode {mode}",
        )
        if index % 8 == 0:
            print(f"checked local exact affine rows {index}/84", flush=True)

    fresh_sums = module.aggregate_rows(fresh_rows)
    saved_sums = artifact["certificate"]["signed_sums"]
    for field in fresh_sums:
        require(
            module.arb(fresh_sums[field]).overlaps(module.arb(saved_sums[field])),
            f"higher-precision local exact aggregate misses {field}",
        )

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
    require("contains no tangent exterior" in note, "localized-domain interpretation missing")
    require("No exact exterior affine" in note, "proof boundary missing")
    print("independently checked finite-box local exact-domain affine A carrier", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
