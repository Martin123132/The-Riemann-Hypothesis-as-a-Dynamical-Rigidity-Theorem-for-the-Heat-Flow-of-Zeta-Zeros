#!/usr/bin/env python3
"""Validate the complete recursive block-21 signed bridge artifact."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate as gate


SERIALIZATION_TOLERANCE = Fraction(1, 10**52)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), f"missing full block-21 bridge result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing full block-21 bridge note: {gate.NOTE}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_full_gate", "full block-21 bridge kind drift")
    require(artifact["status"] == "complete_block21_exact_point_signed_bridge_and_outer_transport_closed", "full block-21 bridge status drift")
    require(bool(artifact["passed"]), "full block-21 bridge did not pass")
    rows = artifact["rows"]
    outputs = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_CALLS, "full block-21 bridge row count drift")
    require(len({int(row["chain"]) for row in rows}) == gate.EXPECTED_CALLS, "full block-21 bridge duplicate chain")
    require(all(int(row["block"]) == 21 and int(row["parent_length"]) == 119 for row in rows), "full block-21 roster drift")
    require(all(bool(row["precision_overlap"]) and bool(row["signed_bridge_closed"]) for row in rows), "full block-21 closure drift")
    require(all(bool(row["exact_correction_within_retained_majorant"]) for row in rows), "full block-21 retained bound drift")
    require(len(outputs) == 15 and [int(row["output_index"]) for row in outputs] == list(range(1, 16)), "full block-21 output roster drift")
    require(all(bool(row["within_requested_scale"]) for row in outputs), "full block-21 promoted output failed")
    for row in outputs:
        retained = Fraction(row["block21_retained_majorant_upper"])
        exact = Fraction(row["block21_exact_correction_majorant_upper"])
        prior = Fraction(row["prior_hybrid_complete_majorant_upper"])
        promoted = Fraction(row["promoted_hybrid_complete_majorant_upper"])
        require(exact < retained, f"output {row['output_index']} block-21 exact column did not improve")
        require(abs(promoted - (prior - retained + exact)) <= SERIALIZATION_TOLERANCE, f"output {row['output_index']} promoted budget identity drift")
        require(abs(Fraction(row["margin_below_requested_scale"]) - (gate.REQUESTED_ERROR_SCALE - promoted)) <= SERIALIZATION_TOLERANCE, f"output {row['output_index']} margin drift")
    aggregate = artifact["aggregate"]
    require(int(aggregate["call_count"]) == gate.EXPECTED_CALLS, "full block-21 aggregate count drift")
    require(int(aggregate["signed_bridge_closed_count"]) == gate.EXPECTED_CALLS, "full block-21 aggregate closure drift")
    require(int(aggregate["exact_corrections_within_retained_count"]) == gate.EXPECTED_CALLS, "full block-21 aggregate retained drift")
    require(int(aggregate["outputs_within_requested_scale"]) == 15, "full block-21 aggregate output drift")
    require(not bool(aggregate["uses_signed_cross_call_cancellation"]), "full block-21 used forbidden cancellation")
    require(int(aggregate["remaining_recursive_exact_point_calls_blocks22_28"]) == 733, "full block-21 remaining roster drift")
    provenance = artifact["provenance"]
    for name, expected in (
        ("telemetry", gate.TELEMETRY),
        ("endpoints", gate.ENDPOINTS),
        ("weights", gate.WEIGHTS),
        ("pilot_result", gate.PILOT_RESULT),
        ("prior_hybrid", gate.HYBRID_RESULT),
        ("builder", gate.BUILDER),
        ("checker", gate.CHECKER),
    ):
        require(provenance[name]["sha256"] == file_hash(expected), f"full block-21 {name} hash drift")
    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("307", "733", "No correction", "not a proof"):
        require(token in note, f"full block-21 note token missing: {token}")
    print(
        "validated complete block-21 signed bridge: "
        f"calls={len(rows)}, promoted_max={aggregate['maximum_promoted_hybrid_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
