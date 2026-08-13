#!/usr/bin/env python3
"""Validate the complete later-recursive signed coefficient-cell bridge."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate as gate


SERIALIZATION_TOLERANCE = Fraction(1, 10**52)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing complete later coefficient-cell result")
    require(gate.NOTE.is_file(), "missing complete later coefficient-cell note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_blocks21_28_coefficient_cell_signed_bridge_full_gate", "kind drift")
    require(artifact["status"] == "all_1040_later_recursive_signed_coefficient_cells_and_15_outputs_closed", "status drift")
    require(bool(artifact["passed"]), "full coefficient-cell gate did not pass")
    require(artifact["precisions_decimal_digits"] == [70, 110], "precision ladder drift")
    rows = artifact["rows"]
    outputs = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_CALLS, "full coefficient-cell row count drift")
    require(len({int(row["chain"]) for row in rows}) == gate.EXPECTED_CALLS, "duplicate full coefficient-cell chain")
    require(all(21 <= int(row["block"]) <= 28 and int(row["child_length"]) == 1 for row in rows), "full coefficient-cell roster drift")
    for row in rows:
        require(row["selector"] in ("W3", "W4") and row["phi3_sign"] in ("negative", "positive"), "cell branch drift")
        require(bool(row["radii_nonzero"]), "zero-radius full cell")
        require(bool(row["precision_overlap"]) and bool(row["point_correction_overlap"]), "full cell overlap drift")
        require(bool(row["cell_correction_within_retained_majorant"]), "full cell exceeds retained majorant")
        require(bool(row["signed_bridge_closed_on_analytic_cell"]), "full signed cell bridge drift")
        require(Fraction(row["cell_to_point_inflation_upper"]) <= Fraction(5, 4), "full cell inflation cap drift")
    require(len(outputs) == 15 and [int(row["output_index"]) for row in outputs] == list(range(1, 16)), "output roster drift")
    for row in outputs:
        point = Fraction(row["later_recursive_point_correction_upper"])
        cell = Fraction(row["later_recursive_cell_correction_upper"])
        prior = Fraction(row["prior_all_recursive_exact_complete_majorant_upper"])
        complete = Fraction(row["all_recursive_cell_complete_majorant_upper"])
        require(cell >= point, f"output {row['output_index']} cell column did not enclose point column")
        require(abs(complete - (prior - point + cell)) <= SERIALIZATION_TOLERANCE, f"output {row['output_index']} budget identity drift")
        require(abs(Fraction(row["margin_below_requested_scale"]) - (gate.REQUESTED_SCALE - complete)) <= SERIALIZATION_TOLERANCE, f"output {row['output_index']} margin drift")
        require(bool(row["within_requested_scale"]), f"output {row['output_index']} exceeds requested scale")
        require(set(row["cell_correction_by_block"]) == {str(block) for block in gate.EXPECTED_BY_BLOCK}, "output block ledger drift")
    aggregate = artifact["aggregate"]
    require(aggregate["call_count_by_block"] == {str(block): count for block, count in gate.EXPECTED_BY_BLOCK.items()}, "block histogram drift")
    for key in ("call_count", "closed_cell_count", "precision_overlap_count", "point_overlap_count", "nonzero_radius_cell_count", "within_retained_majorant_count"):
        require(int(aggregate[key]) == gate.EXPECTED_CALLS, f"aggregate {key} drift")
    require(int(aggregate["w3_count"]) + int(aggregate["w4_count"]) == gate.EXPECTED_CALLS, "selector aggregate drift")
    require(int(aggregate["negative_phi3_count"]) + int(aggregate["positive_phi3_count"]) == gate.EXPECTED_CALLS, "cubic-sign aggregate drift")
    require(int(aggregate["outputs_within_requested_scale"]) == 15, "output closure aggregate drift")
    require(not bool(aggregate["uses_signed_cross_call_cancellation"]), "forbidden cross-call cancellation")
    require(int(aggregate["remaining_later_corrected_model_point_only_calls"]) == 0, "remaining point-only call drift")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")
    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("all `1040` recursive", "not a complete recurrence", "Source binary128 rounding", "outer Hardy", "does not prove"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated complete later coefficient-cell signed bridge: "
        f"1040/1040 cells, 15/15 outputs, max={aggregate['maximum_all_recursive_cell_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
