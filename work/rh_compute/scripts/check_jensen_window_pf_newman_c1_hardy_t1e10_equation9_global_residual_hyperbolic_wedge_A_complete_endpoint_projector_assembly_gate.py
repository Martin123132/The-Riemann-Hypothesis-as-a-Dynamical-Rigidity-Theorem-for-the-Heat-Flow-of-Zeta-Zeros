#!/usr/bin/env python3
"""Independently check the complete exact A endpoint/projector assembly."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_c1_hardy_t1e10_equation9_global_residual_hyperbolic_wedge_A_complete_endpoint_projector_assembly_gate"
RESULT = REPO_ROOT / f"work/rh_compute/results/{STEM}.json"
NOTE = REPO_ROOT / f"outputs/{STEM}.md"
BUILDER = Path(__file__).with_name(STEM + ".py")
CHECKER_PRECISION = 110


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_builder_module():
    spec = importlib.util.spec_from_file_location("A_complete_assembly_builder", BUILDER)
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
        "complete_exact_A_endpoint_carrier_certified",
        "positive_full_line_Gamma_jump_certified",
        "projector_completed_exact_A_transition_certified",
        "endpoint_and_bulk_projector_ownership_kept_separate",
    ):
        require(decisions.get(key) is True, f"missing decision: {key}")
    for key in ("complete_R_Dir_proved", "complete_Q_K_minus_T_proved", "rh_implication"):
        require(decisions.get(key) is False, f"proof-boundary drift: {key}")

    module = load_builder_module()
    module.symbolic_certificate()
    dependencies = {name: load_json(module.REPO_ROOT / record["path"]) for name, record in artifact["dependencies"].items()}
    fresh = module.assembly_certificate(dependencies, CHECKER_PRECISION)
    saved = artifact["certificate"]
    require(fresh["outer_Gamma_mode_range"] == saved["outer_Gamma_mode_range"], "outer roster range drift")
    require(fresh["outer_Gamma_mode_count"] == saved["outer_Gamma_mode_count"] == 42, "outer roster count drift")
    for field, value in fresh["components"].items():
        require(module.arb(value).overlaps(module.arb(saved["components"][field])), f"component misses: {field}")
    saved_rows = {row["mode"]: row for row in saved["outer_Gamma_rows"]}
    require(set(saved_rows) == set(range(39_895, 39_937)), "saved outer Gamma roster drift")
    for row in fresh["outer_Gamma_rows"]:
        old = saved_rows[row["mode"]]
        require(old["positive_full_line_bulk_projector_coefficient"] == 1, "projector coefficient drift")
        for field in (
            "exact_Gamma_full_line_physical_ball",
            "exact_classical_full_line_physical_ball",
            "Gamma_minus_classical_physical_ball",
        ):
            require(module.arb(row[field]).overlaps(module.arb(old[field])), f"mode {row['mode']} misses {field}")

    for record in artifact["dependencies"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"dependency hash drift: {path}")
    for record in artifact["sources"].values():
        path = REPO_ROOT / record["path"]
        require(file_hash(path) == record["sha256"], f"source hash drift: {path}")
    require(NOTE.is_file(), "missing note")
    note = NOTE.read_text(encoding="utf-8")
    require("different projector ownership" in note, "projector ownership guard missing")
    require("complete `R_Dir`" in note, "proof boundary missing")
    print("independently checked complete exact A endpoint and projector transition", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
