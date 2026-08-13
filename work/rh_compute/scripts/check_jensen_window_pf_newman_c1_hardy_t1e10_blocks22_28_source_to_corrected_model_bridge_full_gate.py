#!/usr/bin/env python3
"""Validate the complete blocks 22--28 signed bridge artifact."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate as gate


SERIALIZATION_TOLERANCE = Fraction(1, 10**52)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), f"missing later-block full result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing later-block full note: {gate.NOTE}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_blocks22_28_source_to_corrected_model_bridge_full_gate", "later-block full kind drift")
    require(artifact["status"] == "all_1414_recursive_exact_point_bridges_and_outer_transport_closed", "later-block full status drift")
    require(bool(artifact["passed"]), "later-block full gate did not pass")
    rows = artifact["rows"]
    outputs = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_CALLS, "later-block full row count drift")
    require(len({int(row["chain"]) for row in rows}) == gate.EXPECTED_CALLS, "later-block full duplicate chain")
    require(all(22 <= int(row["block"]) <= 28 and int(row["child_length"]) == 1 for row in rows), "later-block full roster drift")
    require(all(bool(row["precision_overlap"]) and bool(row["signed_bridge_closed"]) for row in rows), "later-block full closure drift")
    require(all(bool(row["exact_correction_within_retained_majorant"]) for row in rows), "later-block full retained majorant drift")
    require(len(outputs) == 15 and [int(row["output_index"]) for row in outputs] == list(range(1, 16)), "later-block full output roster drift")
    for row in outputs:
        retained = Fraction(row["blocks22_28_retained_majorant_upper"])
        exact = Fraction(row["blocks22_28_exact_correction_majorant_upper"])
        prior = Fraction(row["prior_block21_promoted_complete_upper"])
        final = Fraction(row["all_recursive_exact_complete_majorant_upper"])
        require(exact < retained, f"output {row['output_index']} later exact column did not improve")
        require(abs(final - (prior - retained + exact)) <= SERIALIZATION_TOLERANCE, f"output {row['output_index']} final budget identity drift")
        require(abs(Fraction(row["margin_below_requested_scale"]) - (gate.REQUESTED_ERROR_SCALE - final)) <= SERIALIZATION_TOLERANCE, f"output {row['output_index']} final margin drift")
        require(bool(row["within_requested_scale"]), f"output {row['output_index']} final budget failed")
        require(set(row["exact_correction_by_block"]) == {str(block) for block in range(22, 29)}, f"output {row['output_index']} block ledger drift")
    aggregate = artifact["aggregate"]
    require(int(aggregate["later_block_call_count"]) == gate.EXPECTED_CALLS, "later-block aggregate count drift")
    require(int(aggregate["all_recursive_call_count"]) == 1414, "all-recursive count drift")
    require(int(aggregate["total_call_count"]) == 6784, "complete call count drift")
    require(int(aggregate["signed_bridge_closed_count"]) == gate.EXPECTED_CALLS, "later-block aggregate closure drift")
    require(int(aggregate["exact_corrections_within_retained_count"]) == gate.EXPECTED_CALLS, "later-block aggregate retained drift")
    require(int(aggregate["outputs_within_requested_scale"]) == 15, "later-block aggregate output drift")
    require(int(aggregate["remaining_recursive_exact_point_calls"]) == 0, "remaining recursive exact-point roster drift")
    require(not bool(aggregate["uses_signed_cross_call_cancellation"]), "later-block full used forbidden cancellation")
    provenance = artifact["provenance"]
    for name, expected in (
        ("telemetry", gate.TELEMETRY),
        ("endpoints", gate.ENDPOINTS),
        ("weights", gate.WEIGHTS),
        ("geometry_pilot", gate.PILOT_RESULT),
        ("block21_full", gate.BLOCK21_RESULT),
        ("bridge_evaluator", gate.bridge.BUILDER),
        ("builder", gate.BUILDER),
        ("checker", gate.CHECKER),
    ):
        require(provenance[name]["sha256"] == file_hash(expected), f"later-block full {name} hash drift")
    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("1414", "5,370", "no signed cross-call cancellation", "not a proof"):
        require(token in note, f"later-block full note token missing: {token}")
    print(
        "validated all-recursive exact-point bridge: "
        f"recursive={aggregate['all_recursive_call_count']}, max={aggregate['maximum_all_recursive_exact_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
