#!/usr/bin/env python3
"""Validate the Hardy q selector-cell gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_q_selector_cell_gate.py"
CHECKER = Path(__file__).resolve()
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/resumable/zeta14cubicmult_resumable.f90"
SOURCE_SHA256 = "0fb64f090c21b185f27edc1c9194254cf471e74bfd6af3b88cb7f0cb40ec0c3d"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER, SOURCE):
        require(path.is_file(), f"missing q-selector artifact: {path}")
    require(file_hash(SOURCE) == SOURCE_SHA256, "accepted source hash drift")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")

    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_q_selector_cell_gate", "q-cell kind drift")
    require(
        artifact["status"] == "rigorous_local_q_selector_cells_with_physical_coverage_and_analytic_values_open",
        "q-cell status drift",
    )
    require(artifact["source"]["sha256"] == SOURCE_SHA256, "q-cell source hash drift")
    require(len(artifact["source"]["locations"]) == 17, "q-cell source inventory drift")

    contract = artifact["cell_contract"]
    require(contract["coordinates"] == ["a1", "xr=2*a2", "a3", "fracL"], "cell coordinates drift")
    require(len(contract["preserved"]) == 8, "preserved selector inventory drift")
    require("No physical rae/t trajectory" in contract["not_preserved_or_proved"], "physical coverage guard missing")
    require("no source PSI/ERF approximation error" in contract["not_preserved_or_proved"], "special-function guard missing")

    aggregate = artifact["aggregate"]
    require(aggregate["cell_count"] == 64, "q-cell count drift")
    require(aggregate["branch_counts"] == {"1": 32, "2": 32}, "q-cell branch count drift")
    require(aggregate["sum_index_range"] == [1, 32], "q-cell sum-index drift")
    require(0 <= aggregate["maximum_halvings"] <= 30, "q-cell shrink count unreasonable")
    require(sum(aggregate["halving_histogram"].values()) == 64, "q-cell halving histogram drift")
    for radius in aggregate["minimum_radii"].values():
        require(Decimal(radius) > 0, "zero aggregate q-cell radius")
    for name in (
        "minimum_linear_selector_boundary_gap",
        "minimum_quadratic_selector_boundary_gap",
        "minimum_psi_boundary_gap",
        "minimum_erf_boundary_gap",
        "minimum_denominator_abs_lower",
        "minimum_stationary_discriminant_lower",
        "minimum_frac_endpoint_margin",
        "minimum_parent_xr_lower",
    ):
        require(Decimal(aggregate[name]) > 0, f"nonpositive aggregate margin: {name}")
    require(Decimal(aggregate["minimum_denominator_abs_lower"]) > Decimal("0.9"), "q-cell denominator margin too small")
    require(Decimal(aggregate["minimum_stationary_discriminant_lower"]) > Decimal("0.99"), "q-cell discriminant margin too small")

    rows = artifact["rows"]
    require(len(rows) == 64, "q-cell row count drift")
    require([row["chain"] for row in rows] == list(range(1, 65)), "q-cell sequence gap")
    require([row["sum_index"] for row in rows[::2]] == list(range(1, 33)), "q-cell paired roster drift")
    for row in rows:
        require(row["child_length"] in (1, 2), f"chain {row['chain']} child length drift")
        require(row["q_branches"]["ip1"] == 1, f"chain {row['chain']} ip1 drift")
        require(row["q_branches"]["newton_used"] is False, f"chain {row['chain']} Newton drift")
        require(row["q_branches"]["psi_call_count"] == 6, f"chain {row['chain']} PSI roster drift")
        require(row["q_branches"]["erf_call_count"] == 2, f"chain {row['chain']} ERF roster drift")
        require(row["child_selectors"]["source_conjugates"] == (row["child_selectors"]["child_xr_sign"] > 0), f"chain {row['chain']} orientation drift")
        require(all(Decimal(radius["decimal"]) > 0 for radius in row["radii"].values()), f"chain {row['chain']} zero radius")
        require(Decimal(row["analytic_margins"]["minimum_denominator_abs_lower"]) > Decimal("0.9"), f"chain {row['chain']} denominator drift")
        require(Decimal(row["analytic_margins"]["minimum_stationary_discriminant_lower"]) > Decimal("0.99"), f"chain {row['chain']} discriminant drift")

    handoff = artifact["next_handoff"]
    require("physical rae parameter" in handoff["first"], "physical map handoff drift")
    require("rigorous special-function balls" in handoff["second"], "special-function handoff drift")
    require("Reject a box immediately" in handoff["falsification_rule"], "falsification rule drift")

    sources = artifact["sources"]
    for key in ("telemetry", "child_adapter"):
        path = REPO_ROOT / sources[key]["path"]
        require(file_hash(path) == sources[key]["sha256"], f"{key} dependency hash drift")
    require(file_hash(BUILDER) == sources["builder"]["sha256"], "builder hash drift")
    require(file_hash(CHECKER) == sources["checker"]["sha256"], "checker hash drift")

    for token in (
        "Date: 2026-08-06",
        "Status: rigorous local selector cells; not a proof and physical coverage/analytic q open",
        "All `64` saved chains admit such a cell.",
        "all six `PSI` branch paths",
        "both `ERF` branch paths",
        "map actual `rae` intervals into these boxes",
        "does not prove that a physical coefficient trajectory enters or remains in any box",
        "prize-level conclusion",
    ):
        require(token in note, f"q-selector note token missing: {token}")

    print(
        "validated Hardy q selector-cell gate: "
        f"{aggregate['cell_count']} cells, maximum {aggregate['maximum_halvings']} halvings, "
        "6 PSI and 2 ERF paths per cell, 0 physical coverage claims"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
