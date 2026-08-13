#!/usr/bin/env python3
"""Validate the later complete-recurrence four-cell stress pilot."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate as gate


SERIALIZATION_TOLERANCE = Fraction(1, 10**52)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing later stress-cell pilot result")
    require(gate.NOTE.is_file(), "missing later stress-cell pilot note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_later_complete_recurrence_stress_cell_pilot_gate", "kind drift")
    require(artifact["status"] == "four_later_stress_cells_close_complete_tail_and_pinned_binary128_rounding_pilot", "status drift")
    require(bool(artifact["passed"]), "stress-cell pilot failed")
    rows = artifact["rows"]
    require([int(row["chain"]) for row in rows] == list(gate.WITNESS_CHAINS), "stress-cell roster drift")
    require([int(row["block"]) for row in rows] == [21, 21, 28, 28], "stress block roster drift")
    for row in rows:
        require(int(row["child_length"]) == 1, "stress child-length drift")
        require(bool(row["precision_overlap"]) and bool(row["point_tail_contained_in_cell_tail"]), "tail enclosure drift")
        require(Fraction(row["minimum_stationary_discriminant_lower"]) > 0, "stationary discriminant drift")
        require(Fraction(row["cell_complete_tail_majorant_upper"]) >= Fraction(row["point_complete_tail_majorant_upper"]), "cell tail lost point")
        require(bool(row["parent_binary128_rounding"]["all_preimages_inside_cell"]), "parent rounding guard drift")
        for guard in (
            "child_binary128_rounding_preimages_inside_reconstructed_cell",
            "multiplier_binary128_rounding_preimages_inside_anchored_cell",
            "child_state_binary128_rounding_preimages_inside_cell_kernel",
        ):
            require(all(bool(value) for value in row[guard]), f"{guard} drift")
        require(bool(row["tpm_rounding_preimage_contains_minus_two_pi"]), "tpm rounding drift")
        require(
            Fraction(row["source_complex_multiply_add_observed_gap"])
            <= Fraction(row["source_complex_multiply_add_point_roundoff_upper"])
            <= Fraction(row["source_complex_multiply_add_cell_roundoff_upper"]),
            "source multiply-add roundoff drift",
        )

    outputs = artifact["output_rows"]
    require(len(outputs) == 15 and [int(row["output_index"]) for row in outputs] == list(range(1, 16)), "output roster drift")
    for row in outputs:
        prior = Fraction(row["prior_all_recursive_cell_complete_majorant_upper"])
        tail = Fraction(row["four_stress_cell_complete_tail_upper"])
        roundoff = Fraction(row["four_stress_cell_source_roundoff_upper"])
        augmented = Fraction(row["stress_augmented_complete_majorant_upper"])
        require(abs(augmented - (prior + tail + roundoff)) <= SERIALIZATION_TOLERANCE, "stress output budget identity drift")
        require(abs(Fraction(row["margin_below_requested_scale"]) - (gate.REQUESTED_SCALE - augmented)) <= SERIALIZATION_TOLERANCE, "stress output margin drift")
        require(bool(row["within_requested_scale"]), "stress output exceeds requested scale")

    aggregate = artifact["aggregate"]
    require(int(aggregate["stress_cell_count"]) == 4, "stress aggregate count drift")
    for key in (
        "precision_overlap_count",
        "point_tail_containment_count",
        "parent_rounding_preimage_containment_count",
        "child_rounding_preimage_containment_count",
        "multiplier_rounding_preimage_containment_count",
        "child_state_rounding_preimage_containment_count",
        "tpm_rounding_containment_count",
        "source_multiply_add_roundoff_containment_count",
    ):
        require(int(aggregate[key]) == 4, f"aggregate {key} drift")
    require(int(aggregate["outputs_with_stress_augmented_majorant_within_requested_scale"]) == 15, "stress output closure drift")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed dependency {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")
    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("chains 3406 and 3423", "round-to-nearest ties-to-even", "Source `q`", "not a proof of the complete recurrence or RH"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated later complete-recurrence stress-cell pilot: "
        f"4/4 cells, 15/15 outputs, max={aggregate['maximum_stress_augmented_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
