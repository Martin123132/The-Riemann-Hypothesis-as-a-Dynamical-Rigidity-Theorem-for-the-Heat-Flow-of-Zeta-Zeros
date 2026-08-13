#!/usr/bin/env python3
"""Validate the finite corrected-model source-q dependency cut."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate as gate


SERIALIZATION_TOLERANCE = Fraction(1, 10**52)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing source-q dependency-cut result")
    require(gate.NOTE.is_file(), "missing source-q dependency-cut note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_rigorous_arb_source_q_dependency_cut_gate",
        "kind drift",
    )
    require(
        artifact["status"] == "finite_saved_height_corrected_model_arb_path_has_no_emitted_source_q_payload_dependency",
        "status drift",
    )
    require(bool(artifact["passed"]), "source-q dependency cut failed")
    decision = artifact["decision"]
    require(not bool(decision["source_q_payload_is_corrected_model_antecedent"]), "source q re-entered corrected-model path")
    require(bool(decision["source_q_payload_remains_legacy_reproduction_obligation"]), "legacy reproduction boundary was lost")

    manifest = artifact["evaluator_manifest"]
    require(len(manifest) == 14, "evaluator manifest size drift")
    for row in manifest:
        require(not row["forbidden_payload_token_hits"], f"forbidden token hit in {row['function']}")
        require(set(row["data_keys"]) <= {"levels"}, f"non-level evaluator data access in {row['function']}")

    replay = artifact["stripped_replay_rows"]
    require([int(row["chain"]) for row in replay] == list(gate.STRESS_CHAINS), "stripped replay roster drift")
    for row in replay:
        require(row["supplied_top_level_keys"] == ["levels"], "stripped replay gained source payload")
        require(int(row["precision_decimal_digits"]) == gate.REPLAY_PRECISION, "replay precision drift")
        require(bool(row["saved_correction_overlap"]), "stripped replay overlap drift")

    projection = artifact["projection"]
    require(int(projection["row_count"]) == gate.EXPECTED_CALLS, "projection roster drift")
    require(int(projection["canonical_json_bytes"]) > 0, "empty projection")
    require(len(projection["sha256"]) == 64, "projection hash drift")
    require(int(projection["forbidden_key_hit_count"]) == 0, "projection contains forbidden payload")
    require(set(projection["forbidden_payload_keys"]) == gate.FORBIDDEN_PAYLOAD_KEYS, "forbidden key policy drift")

    outputs = artifact["output_rows"]
    require(len(outputs) == 15 and [int(row["output_index"]) for row in outputs] == list(range(1, 16)), "output roster drift")
    for row in outputs:
        cell = Fraction(row["source_q_free_cell_majorant_upper"])
        complete = Fraction(row["source_q_free_complete_majorant_upper"])
        tail = Fraction(row["later_analytic_complete_tail_upper"])
        margin = Fraction(row["margin_below_requested_scale"])
        require(abs(complete - cell - tail) <= SERIALIZATION_TOLERANCE, f"output {row['output_index']} complete identity drift")
        require(Fraction(row["dominance_over_stored_cell_upper"]) >= 0, "source-q-free cell lost accepted bound")
        require(Fraction(row["dominance_over_stored_partial_upper"]) >= 0, "source-q-free complete lost accepted bound")
        require(abs(margin - (gate.REQUESTED_SCALE - complete)) <= SERIALIZATION_TOLERANCE, "source-q-free margin drift")
        require(bool(row["within_requested_scale"]) and complete < gate.REQUESTED_SCALE, "source-q-free output exceeds scale")

    aggregate = artifact["aggregate"]
    require(int(aggregate["evaluator_function_count"]) == len(manifest), "function aggregate drift")
    require(int(aggregate["forbidden_payload_token_hit_count"]) == 0, "forbidden source token aggregate drift")
    require(int(aggregate["stripped_replay_count"]) == len(gate.STRESS_CHAINS), "replay aggregate drift")
    require(int(aggregate["stripped_replay_overlap_count"]) == len(gate.STRESS_CHAINS), "replay overlap aggregate drift")
    require(int(aggregate["projected_later_cell_count"]) == gate.EXPECTED_CALLS, "cell aggregate drift")
    require(int(aggregate["outputs_within_requested_scale"]) == 15, "output aggregate drift")
    require(not bool(aggregate["uses_source_q_payload"]), "source q usage aggregate drift")
    require(not bool(aggregate["uses_signed_cross_call_cancellation"]), "cross-call cancellation aggregate drift")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("only access", "1040", "source-q-free", "separate task", "outer Hardy", "not RH"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated rigorous Arb source-q dependency cut: "
        f"1040 cells, 4 stripped replays, 15 outputs, max={aggregate['maximum_source_q_free_complete_majorant_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
