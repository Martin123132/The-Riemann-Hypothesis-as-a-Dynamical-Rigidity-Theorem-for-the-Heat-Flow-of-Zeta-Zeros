#!/usr/bin/env python3
"""Validate the complete later-cell recurrence-tail transport."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate as gate


SERIALIZATION_TOLERANCE = Fraction(1, 10**52)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing complete later recurrence-tail result")
    require(gate.NOTE.is_file(), "missing complete later recurrence-tail note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_cell_full_gate", "kind drift")
    require(artifact["status"] == "all_1040_later_cells_close_complete_tail_and_pinned_source_multiply_add_rounding", "status drift")
    require(bool(artifact["passed"]), "complete later recurrence-tail gate failed")
    rows = artifact["rows"]
    require(len(rows) == gate.EXPECTED_CALLS, "later recurrence-tail row count drift")
    require(len({int(row["chain"]) for row in rows}) == gate.EXPECTED_CALLS, "duplicate later recurrence-tail chain")
    for row in rows:
        require(21 <= int(row["block"]) <= 28 and int(row["child_length"]) == 1, "later recurrence-tail roster drift")
        require(bool(row["precision_overlap"]) and bool(row["point_tail_contained_in_cell_tail"]), "later tail enclosure drift")
        require(Fraction(row["minimum_stationary_discriminant_lower"]) > 0, "later stationary discriminant drift")
        require(Fraction(row["cell_complete_tail_majorant_upper"]) >= Fraction(row["point_complete_tail_majorant_upper"]), "later cell tail lost point")
        for key in (
            "parent_binary128_rounding_preimages_inside_cell",
            "child_binary128_rounding_preimages_inside_reconstructed_cell",
            "multiplier_binary128_rounding_preimages_inside_anchored_cell",
            "child_state_binary128_rounding_preimages_inside_cell_kernel",
            "tpm_rounding_preimage_contains_minus_two_pi",
        ):
            require(bool(row[key]), f"later rounding guard {key} drift")
        require(
            Fraction(row["source_complex_multiply_add_observed_gap"])
            <= Fraction(row["source_complex_multiply_add_point_roundoff_upper"])
            <= Fraction(row["source_complex_multiply_add_cell_roundoff_upper"]),
            "later source multiply-add roundoff drift",
        )

    outputs = artifact["output_rows"]
    require(len(outputs) == 15 and [int(row["output_index"]) for row in outputs] == list(range(1, 16)), "later output roster drift")
    for row in outputs:
        prior = Fraction(row["prior_corrected_model_cell_complete_majorant_upper"])
        tail = Fraction(row["all_later_cell_complete_tail_majorant_upper"])
        roundoff = Fraction(row["all_later_cell_source_multiply_add_roundoff_upper"])
        complete = Fraction(row["later_complete_recurrence_partial_majorant_upper"])
        require(abs(complete - (prior + tail + roundoff)) <= SERIALIZATION_TOLERANCE, "later complete budget identity drift")
        require(abs(Fraction(row["margin_below_requested_scale"]) - (gate.REQUESTED_SCALE - complete)) <= SERIALIZATION_TOLERANCE, "later complete margin drift")
        require(set(row["tail_by_block"]) == {str(block) for block in gate.EXPECTED_BY_BLOCK}, "later tail block ledger drift")
        require(bool(row["within_requested_scale"]), "later complete output exceeds requested scale")

    aggregate = artifact["aggregate"]
    require(aggregate["call_count_by_block"] == {str(block): count for block, count in gate.EXPECTED_BY_BLOCK.items()}, "later block histogram drift")
    for key in (
        "call_count",
        "precision_overlap_count",
        "point_tail_containment_count",
        "parent_rounding_preimage_containment_count",
        "child_rounding_preimage_containment_count",
        "multiplier_rounding_preimage_containment_count",
        "child_state_rounding_preimage_containment_count",
        "tpm_rounding_containment_count",
        "source_multiply_add_roundoff_containment_count",
    ):
        require(int(aggregate[key]) == gate.EXPECTED_CALLS, f"later aggregate {key} drift")
    require(int(aggregate["outputs_within_requested_scale"]) == 15, "later output closure drift")
    require(not bool(aggregate["uses_cross_call_cancellation"]), "forbidden later cross-call cancellation")
    require(int(aggregate["remaining_later_point_only_tail_calls"]) == 0, "remaining later point-tail call drift")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed later dependency {path}")
            require(file_hash(path) == record["sha256"], f"later dependency hash drift {path}")
    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("all `1040`", "round-to-nearest", "source `q`", "No cancellation", "not a proof"):
        require(token in note, f"later note boundary token missing: {token}")
    print(
        "validated complete later recurrence-tail cell transport: "
        f"1040/1040 cells, 15/15 outputs, max={aggregate['maximum_later_complete_recurrence_partial_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
