#!/usr/bin/env python3
"""Validate the Hardy q-cell physical-spacing obstruction gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing physical-spacing artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")

    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_q_cell_physical_spacing_gate",
        "physical-spacing kind drift",
    )
    require(
        artifact["status"] == "rigorous_adjacent_endpoint_q_cell_obstruction_with_discrete_route_open",
        "physical-spacing status drift",
    )
    contract = artifact["test_contract"]
    require(contract["physical_roster"] == "rae_j=657064+420*(j-1), j=1,...,32", "roster contract drift")
    require(contract["branches"] == [1, 2], "branch contract drift")
    require(contract["tested_adjacent_pairs_per_branch"] == 31, "adjacency contract drift")
    require("min_{n in Z}" in contract["phase_metric"], "circle metric drift")

    aggregate = artifact["aggregate"]
    require(aggregate["adjacent_pair_count"] == 62, "adjacent pair count drift")
    require(aggregate["phase_projected_disjoint_count"] == 62, "phase obstruction count drift")
    require(aggregate["endpoint_continuous_cover_count"] == 0, "unexpected endpoint cover")
    require(Decimal(aggregate["minimum_phase_projected_gap"]["decimal"]) > 0, "nonpositive phase gap")
    require(
        Decimal(aggregate["minimum_phase_distance_to_radius_sum_ratio"]["decimal"]) > 1,
        "phase distance/radius ratio is not separated",
    )

    rows = artifact["rows"]
    require(len(rows) == 62, "physical-spacing row count drift")
    expected = [(branch, index) for branch in (1, 2) for index in range(1, 32)]
    require([(row["branch"], row["left_sum_index"]) for row in rows] == expected, "adjacency sequence drift")
    for row in rows:
        require(row["right_sum_index"] == row["left_sum_index"] + 1, "non-adjacent row")
        require(row["right_rae"] - row["left_rae"] == 420, "physical step drift")
        require(row["rae_step"] == 420, "stored physical step drift")
        require(row["endpoint_q_cells_form_continuous_cover"] is False, "false cover promoted")
        phase = row["coordinates"]["a1_mod_one"]
        require(phase["endpoint_projections_disjoint"] is True, "phase arcs not disjoint")
        require(Decimal(phase["signed_projected_gap"]["decimal"]) > 0, "row phase gap nonpositive")
        require(Decimal(phase["distance_to_radius_sum_ratio"]["decimal"]) > 1, "row phase ratio nonpositive")

    decision = artifact["route_decision"]
    require("endpoint q cells" in decision["rejected"], "rejected route drift")
    require("discrete pivot index" in decision["recommended"], "discrete route missing")
    require("adaptive intermediate q cells" in decision["fallback"], "fallback route missing")
    require("exact positive phase gaps" in decision["falsification_rule"], "future falsification rule missing")

    for key in ("q_cells", "rae_gate"):
        dependency = REPO_ROOT / artifact["sources"][key]["path"]
        require(file_hash(dependency) == artifact["sources"][key]["sha256"], f"{key} dependency hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "checker hash drift")

    for token in (
        "Date: 2026-08-06",
        "Status: rigorous adjacent endpoint-cell obstruction",
        "rae_j = 657064 + 420*(j-1)",
        "distance on the circle",
        "every one of the `62`",
        "cannot stay in the union",
        "islands around discrete source calls",
        "all 212 pivots",
        "does not invalidate the q-cell certificates",
        "prize-level conclusion",
    ):
        require(token in note, f"physical-spacing note token missing: {token}")

    require("rejects only adjacent continuous coverage" in artifact["proof_boundary"], "proof boundary drift")
    print(
        "validated Hardy q-cell physical-spacing gate: "
        f"{aggregate['phase_projected_disjoint_count']}/62 phase projections disjoint, "
        "discrete-index route selected"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
