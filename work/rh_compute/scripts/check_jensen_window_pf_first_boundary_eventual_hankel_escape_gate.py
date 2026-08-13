#!/usr/bin/env python3
"""Independently check the first-boundary eventual-Hankel escape gate."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = "jensen_window_pf_first_boundary_eventual_hankel_escape_gate"
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

EXPECTED_COUNTS = {
    "gate_rows": 15,
    "source_artifacts": 9,
    "effective_order_inventory": 11,
    "collision_pairs": 190,
    "propagation_index_checks": 760,
    "global_first_boundary_covered": 190,
    "global_first_boundary_uncovered": 0,
    "local_effective_all": 139,
    "local_effective_tail_only": 6,
    "local_open": 45,
    "order10_low_sensor_crossings": 4,
    "uniform_in_order_thresholds": 0,
    "rh_conclusions": 0,
}

EXPECTED_GATE_ROWS = [
    "fbe_01_sources", "fbe_02_index", "fbe_03_shift_propagation",
    "fbe_04_eventual_tail", "fbe_05_escape", "fbe_06_quantifiers",
    "fbe_07_first_boundary", "fbe_08_bounded_matrix", "fbe_09_local_matrix",
    "fbe_10_order10", "fbe_11_order11", "fbe_12_order12",
    "fbe_13_parallel", "fbe_14_route", "fbe_15_boundary",
]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def independent_local_status(order: int, offset: int) -> tuple[str, int | None]:
    if order < 10:
        return "effective_all_collision_shifts", 0
    if order < 12:
        first_shift = max(0, 4 - offset)
        if first_shift == 0:
            return "effective_all_collision_shifts", 0
        return "effective_tail_collision_shifts_only", first_shift
    return "no_effective_immediate_block_theorem", None


def check_sources(artifact: dict) -> None:
    source = artifact.get("source_audit", {})
    require(len(source) == EXPECTED_COUNTS["source_artifacts"], "source count drifted")
    for key, row in source.items():
        path = REPO_ROOT / row["path"]
        require(path.exists(), f"missing source: {key}")
        require(file_hash(path) == row["sha256"], f"source hash drifted: {key}")


def check_inventory(artifact: dict) -> None:
    inventory = artifact.get("effective_order_inventory", [])
    require([row.get("order") for row in inventory] == list(range(2, 13)), "inventory order drifted")
    for row in inventory:
        order = row["order"]
        require(row.get("first_boundary_status", "").startswith("covered"), f"uncovered eventual order {order}")
        if order <= 9:
            require(row["effective_shift_floor"] == 0, f"low-order floor drifted: {order}")
        elif order <= 11:
            require(row["effective_shift_floor"] == 4, f"delayed floor drifted: {order}")
        else:
            require(row["effective_shift_floor"] is None, "order-twelve floor promoted")


def check_collision_matrix(artifact: dict) -> tuple[int, int, int, int]:
    matrix = artifact.get("collision_matrix", [])
    expected_ids: list[str] = []
    global_covered = 0
    local_all = 0
    local_tail = 0
    local_open = 0
    index = 0
    for degree in range(2, 21):
        for multiplicity in range(2, degree + 1):
            expected_ids.append(f"fb_{degree:02d}_{multiplicity:02d}")
            row = matrix[index]
            index += 1
            p = degree - multiplicity + 1
            order = p + 1
            offset = multiplicity - 1
            status, floor = independent_local_status(order, offset)
            require(row["degree"] == degree, "degree drifted")
            require(row["multiplicity"] == multiplicity, "multiplicity drifted")
            require(row["recurrence_order"] == p, "recurrence order drifted")
            require(row["hankel_order"] == order, "Hankel order drifted")
            require(row["coefficient_offset"] == offset, "coefficient offset drifted")
            require(row["local_effective_status"] == status, "local status drifted")
            require(row["local_collision_shift_floor"] == floor, "local floor drifted")
            require(
                row["first_boundary_status"] == "covered_by_eventual_hankel_shift_escape",
                "global coverage drifted",
            )
            global_covered += 1
            local_all += status == "effective_all_collision_shifts"
            local_tail += status == "effective_tail_collision_shifts_only"
            local_open += status == "no_effective_immediate_block_theorem"
    require([row.get("id") for row in matrix] == expected_ids, "collision row order drifted")
    return global_covered, local_all, local_tail, local_open


def check_gate_rows(artifact: dict) -> None:
    rows = artifact.get("rows", [])
    require([row.get("id") for row in rows] == EXPECTED_GATE_ROWS, "gate row order drifted")
    require(
        [row["id"] for row in rows if row.get("readiness") == "open"]
        == ["fbe_14_route"],
        "open row drifted",
    )
    for row in rows:
        for field in ("role", "readiness", "claim", "certificate", "proof_boundary"):
            require(bool(row.get(field)), f"empty {field}: {row.get('id')}")


def check_note(note: str) -> None:
    markers = [
        "not a proof of RH", "## Index Map", "d-m+2",
        "## Unlimited Shift Propagation", "## Eventual-Hankel Escape",
        "for each fixed d,m", "## First-Boundary Consequence",
        "190 degree/multiplicity pairs", "zero uncovered finite global rows",
        "## Immediate Sensor Matrix", "Four low order-ten determinant sensors",
        "## Independent Route Comparison", "cofinal degree", "## Pi Provenance",
        "## Proof Boundary",
    ]
    for marker in markers:
        require(marker in note, f"note marker missing: {marker}")


def main() -> int:
    require(RESULT_PATH.exists(), "result artifact missing")
    require(NOTE_PATH.exists(), "note artifact missing")
    artifact = json.loads(RESULT_PATH.read_text(encoding="utf-8"))
    note = NOTE_PATH.read_text(encoding="utf-8")
    require(artifact.get("kind") == KIND, "kind drifted")
    require(artifact.get("counts") == EXPECTED_COUNTS, "counts drifted")
    require(artifact.get("audit_range") == {"minimum_degree": 2, "maximum_degree": 20}, "audit range drifted")
    require("cofinal degree limit open" in artifact.get("status", ""), "status drifted")
    require("no threshold uniform in degree" in artifact.get("proof_boundary", ""), "boundary drifted")
    require(
        set(artifact.get("symbolic_certificate", {}))
        == {
            "local_bridge", "unlimited_shift", "eventual_tail", "escape",
            "quantifiers", "first_boundary", "order10_guard", "order12_guard",
            "parallel_route", "open",
        },
        "certificate keys drifted",
    )
    check_sources(artifact)
    check_inventory(artifact)
    global_covered, local_all, local_tail, local_open = check_collision_matrix(artifact)
    require(global_covered == EXPECTED_COUNTS["global_first_boundary_covered"], "global count drifted")
    require(local_all == EXPECTED_COUNTS["local_effective_all"], "local all count drifted")
    require(local_tail == EXPECTED_COUNTS["local_effective_tail_only"], "local tail count drifted")
    require(local_open == EXPECTED_COUNTS["local_open"], "local open count drifted")
    check_gate_rows(artifact)
    check_note(note)
    print(
        "validated first-boundary eventual-Hankel escape gate: "
        f"{EXPECTED_COUNTS['gate_rows']} rows, {EXPECTED_COUNTS['source_artifacts']} sources, "
        f"{EXPECTED_COUNTS['effective_order_inventory']} effective-order records, "
        f"{EXPECTED_COUNTS['collision_pairs']} collision pairs, "
        f"{EXPECTED_COUNTS['propagation_index_checks']} propagation-index checks, "
        f"{global_covered} global rows covered, "
        f"{EXPECTED_COUNTS['global_first_boundary_uncovered']} global rows uncovered, "
        f"{local_all} local all-shift rows, {local_tail} local delayed rows, "
        f"{local_open} local open rows, "
        f"{EXPECTED_COUNTS['uniform_in_order_thresholds']} uniform-order thresholds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
