#!/usr/bin/env python3
"""Validate the finite Hardy block-20 source-to-corrected-model bridge."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def rows_by_chain(path: Path) -> dict[int, dict]:
    artifact = json.loads(path.read_text(encoding="utf-8"))
    return {int(row["chain"]): row for row in artifact["rows"]}


def main() -> int:
    require(gate.RESULT.is_file(), f"missing result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing note: {gate.NOTE}")
    result = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(result["kind"] == "jensen_window_pf_newman_c1_hardy_block20_source_to_corrected_model_bridge_gate", "kind drift")
    require(result["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash drift")
    require(result["sources"]["checker"]["sha256"] == gate.file_hash(Path(__file__).resolve()), "checker hash drift")
    for dependency in result["dependencies"].values():
        path = gate.REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(dependency["sha256"] == gate.file_hash(path), f"dependency hash drift: {path}")
    require(gate.algebra_audit() == 12, "bridge algebra audit drift")

    eq_rows = rows_by_chain(gate.EQ69)
    t2_rows = rows_by_chain(gate.T2)
    special_rows = rows_by_chain(gate.SPECIAL)
    w2_rows = rows_by_chain(gate.W2W5)
    endpoint_rows = rows_by_chain(gate.ENDPOINT)
    rows = result["rows"]
    require(len(rows) == 374 and len({int(row["chain"]) for row in rows}) == 374, "bridge roster drift")
    for row in rows:
        chain = int(row["chain"])
        post_w1 = gate.corrected.complex_from_record(eq_rows[chain]["residual_after_exact_w1"])
        t2_prior = gate.corrected.complex_from_record(t2_rows[chain]["prior_exact_w1_residual"])
        post_t2 = gate.corrected.complex_from_record(t2_rows[chain]["residual_after_deletion"])
        w2_prior = gate.corrected.complex_from_record(w2_rows[chain]["prior_post_lplus1_residual"])
        special_delta = gate.corrected.complex_from_record(
            special_rows[chain]["correlated_replacement"]["total_q_delta"]
        )
        nonsaddle = gate.corrected.complex_from_record(w2_rows[chain]["complete_nonsaddle_transport"])
        normalization = gate.corrected.complex_from_record(row["post_special_normalization_remainder"])
        paper_residual = gate.corrected.complex_from_record(w2_rows[chain]["paper_residual"])
        endpoint_correction = gate.corrected.complex_from_record(endpoint_rows[chain]["exact_correction"])
        require(post_w1.overlaps(t2_prior), f"chain {chain} exact-W1/t2 drift")
        require(post_t2.overlaps(w2_prior), f"chain {chain} t2/W2-W5 drift")
        require((nonsaddle - special_delta).overlaps(normalization), f"chain {chain} normalization split drift")
        require(paper_residual.overlaps(endpoint_correction), f"chain {chain} endpoint overlap drift")
        require(bool(row["endpoint_decomposition_overlap"]), f"chain {chain} endpoint flag drift")
        require(bool(row["precision_and_identity_handoffs"]), f"chain {chain} handoff flag drift")

    aggregate = result["aggregate"]
    require(int(aggregate["exact_algebra_checks"]) == 12, "aggregate algebra drift")
    require(int(aggregate["call_count"]) == 374, "aggregate call drift")
    require(int(aggregate["identity_handoff_count"]) == 374, "aggregate handoff drift")
    require(int(aggregate["endpoint_decomposition_overlap_count"]) == 374, "aggregate endpoint drift")
    require(int(aggregate["outputs_within_requested_scale"]) == 15, "aggregate output drift")

    joint = json.loads(gate.JOINT.read_text(encoding="utf-8"))
    joint_outputs = {int(row["output_index"]): row for row in joint["transported_outputs"]}
    for row in result["transported_outputs"]:
        source = joint_outputs[int(row["output_index"])]
        require(
            Fraction(row["corrected_model_majorant_upper"])
            == Fraction(source["joint_corrected_recurrence_majorant_upper"]),
            "corrected output transport drift",
        )
        require(bool(row["within_requested_scale"]), "corrected output closure drift")
    require(
        result["status"] == "finite_source_to_corrected_model_bridge_closes_at_374_saved_points_and_15_outputs",
        "status drift",
    )
    require("not a proof" in gate.NOTE.read_text(encoding="utf-8").lower(), "note proof boundary missing")
    print("validated source-to-corrected-model bridge: 374 calls, 15/15 outputs below 0.005")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
