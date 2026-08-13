#!/usr/bin/env python3
"""Validate the complete Hardy block-20 discrete selector atlas."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing selector-atlas artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_discrete_selector_atlas_gate",
        "selector-atlas kind drift",
    )
    require(
        artifact["status"] == "complete_424_call_discrete_selector_atlas_with_analytic_q_open",
        "selector-atlas status drift",
    )
    roster = artifact["roster"]
    require(roster["sum_count"] == 212 and roster["branch_count"] == 424, "selector roster drift")
    require(roster["first_rae"] == 657064 and roster["last_rae"] == 745684, "selector endpoint drift")

    aggregate = artifact["aggregate"]
    require(aggregate["classified_call_count"] == 424, "classified call count drift")
    require(aggregate["unclassified_call_count"] == 0, "unclassified call promoted")
    require(aggregate["mit1_chain_count"] == 50, "MIT1 count drift")
    require(aggregate["mit2_chain_count"] == 374, "MIT2 count drift")
    require(aggregate["q_cell_count"] == 374, "full q-cell count drift")
    require(sum(aggregate["child_length_histogram"].values()) == 374, "child length histogram drift")
    require(sum(aggregate["q_cell_halving_histogram"].values()) == 374, "halving histogram drift")
    require(aggregate["maximum_q_cell_halvings"] <= 60, "q-cell halving limit drift")
    require(aggregate["distinct_selector_word_count"] > 1, "selector atlas collapsed")
    for value in aggregate["minimum_q_cell_radii"].values():
        require(Decimal(value) > 0, "zero full-atlas q radius")
    for key in (
        "minimum_psi_boundary_gap",
        "minimum_erf_boundary_gap",
        "minimum_denominator_abs_lower",
        "minimum_stationary_discriminant_lower",
    ):
        require(Decimal(aggregate[key]) > 0, f"nonpositive full-atlas margin: {key}")

    rows = artifact["rows"]
    require(len(rows) == 424, "selector-atlas row count drift")
    for chain, row in enumerate(rows, 1):
        sum_index = (chain + 1) // 2
        branch = 1 if chain % 2 else 2
        require(row["chain"] == chain, "selector chain sequence drift")
        require(row["sum_index"] == sum_index and row["branch"] == branch, "selector physical ordering drift")
        require(row["rae"] == 657064 + 420 * (sum_index - 1), "selector RAE drift")
        if row["mit"] == 1:
            require(row["route"] == "direct_kernel_mit1", "MIT1 route drift")
            require(row["child_length"] is None and row["q_cell"] is None, "MIT1 q data drift")
        else:
            require(row["route"] == "recursive_mit2_q_cell", "MIT2 route drift")
            require(row["child_length"] in (1, 2), "MIT2 child length drift")
            require(row["q_cell"]["chain"] == chain, "q-cell chain drift")
            require(row["q_cell"]["q_branches"]["psi_call_count"] == 6, "PSI roster drift")
            require(row["q_cell"]["q_branches"]["erf_call_count"] == 2, "ERF roster drift")

    mit1_total = sum(len(indices) for indices in artifact["mit1_sum_indices_by_branch"].values())
    require(mit1_total == 50, "MIT1 index roster drift")
    for branch in ("1", "2"):
        segments = artifact["selector_segments_by_branch"][branch]
        require(segments[0]["start_sum_index"] == 1, "selector segment origin drift")
        require(segments[-1]["end_sum_index"] == 212, "selector segment endpoint drift")
        for left, right in zip(segments, segments[1:]):
            require(left["end_sum_index"] + 1 == right["start_sum_index"], "selector segment gap")

    handoff = artifact["next_handoff"]
    require("374 occupied q cells" in handoff["recursive_route"], "recursive handoff drift")
    require("50 MIT=1" in handoff["direct_route"], "direct handoff drift")
    require("do not interpolate" in handoff["falsification_rule"], "transition guard drift")
    for key in ("fixture_result", "telemetry", "prefix_q_cells", "rae_gate", "local_cell_builder"):
        dependency = REPO_ROOT / artifact["sources"][key]["path"]
        require(file_hash(dependency) == artifact["sources"][key]["sha256"], f"{key} dependency hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    for token in (
        "Date: 2026-08-06",
        "Status: complete finite selector atlas",
        "all 212 block-20 pivots",
        "50 direct cases",
        "374 recursive chains",
        "424 / 424 classified branch calls",
        "0 unclassified calls",
        "high fragmentation",
        "does not claim continuous coverage",
        "prize-level conclusion",
    ):
        require(token in note, f"selector-atlas note token missing: {token}")

    require("Complete finite discrete selector atlas" in artifact["proof_boundary"], "selector boundary drift")
    print(
        "validated Hardy block-20 discrete selector atlas: "
        "424 classified calls, 374 q cells, 50 direct kernels, 0 unclassified"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
