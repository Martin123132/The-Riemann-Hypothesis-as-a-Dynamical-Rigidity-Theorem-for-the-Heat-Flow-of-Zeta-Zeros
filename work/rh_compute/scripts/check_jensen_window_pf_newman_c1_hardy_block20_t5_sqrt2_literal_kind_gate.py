#!/usr/bin/env python3
"""Validate the t5 default-real sqrt(2.0) correction artifact."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_t5_sqrt2_literal_kind_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(RESULT.is_file() and NOTE.is_file(), "missing sqrt2 gate output")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact["status"] == "rigorous_913_call_default_real_sqrt_two_argument_correction_enclosed",
        "sqrt2 status drift",
    )
    require(artifact["scope"]["recursive_call_count"] == 374, "sqrt2 chain count drift")
    require(artifact["scope"]["retained_erfc_call_count"] == 913, "sqrt2 call count drift")
    require(artifact["scope"]["precision_overlap_count"] == 374, "sqrt2 precision count drift")
    require(artifact["normalization"]["source_sqrt_two"]["bits"] == "3FB504F3", "sqrt2 bits drift")
    require(len(artifact["rows"]) == 374, "sqrt2 row count drift")
    aggregate = artifact["aggregate"]
    require(sum(aggregate["region_histogram"].values()) == 913, "sqrt2 region count drift")
    require(
        aggregate["sqrt_two_improved_call_count"] + aggregate["sqrt_two_worsened_call_count"]
        + aggregate["sqrt_two_indeterminate_call_count"] == 374,
        "sqrt2 classification partition drift",
    )
    require(
        aggregate["stack_improved_call_count"] + aggregate["stack_worsened_call_count"]
        + aggregate["stack_indeterminate_call_count"] == 374,
        "sqrt2 stack classification partition drift",
    )
    require(
        sum(row["residual_after_sqrt_two_excludes_zero"] for row in artifact["rows"])
        == aggregate["after_sqrt_two_excluding_zero_count"],
        "sqrt2 nonzero recount drift",
    )
    require(
        sum(row["residual_after_stack_excludes_zero"] for row in artifact["rows"])
        == aggregate["after_stack_excluding_zero_count"],
        "sqrt2 stack nonzero recount drift",
    )
    require(all(row["precision_overlap"] for row in artifact["rows"]), "sqrt2 precision flag drift")

    source = artifact["source"]
    source_path = REPO_ROOT / source["path"]
    require(source_path.is_file() and file_hash(source_path) == source["sha256"], "sqrt2 source hash drift")
    for name, record in artifact["sources"].items():
        if not isinstance(record, dict) or "path" not in record:
            continue
        path = REPO_ROOT / record["path"]
        require(path.is_file(), f"missing sqrt2 dependency: {name}")
        require(file_hash(path) == record["sha256"], f"sqrt2 dependency hash drift: {name}")

    note = " ".join(NOTE.read_text(encoding="utf-8").split())
    for token in (
        "0x3FB504F3",
        "principal positive mathematical root",
        "delta x = B * (sqrt(2) - sqrt2_binary32)",
        "does not count that residual twice",
        "R_q-delta_q_sqrt2-T_ip",
        "does not prove that the displayed W1 approximation equals",
    ):
        require(token in note, f"sqrt2 note token missing: {token}")
    print(
        "validated Hardy block-20 t5 sqrt2 literal-kind gate: "
        f"913 calls, {aggregate['after_stack_excluding_zero_count']} stacked residuals exclude zero"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
