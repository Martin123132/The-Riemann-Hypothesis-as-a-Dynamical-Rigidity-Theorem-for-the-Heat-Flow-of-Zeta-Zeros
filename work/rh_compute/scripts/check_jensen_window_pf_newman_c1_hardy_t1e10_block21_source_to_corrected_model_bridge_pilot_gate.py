#!/usr/bin/env python3
"""Validate the block-21 source-to-corrected-model bridge pilot artifact."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), f"missing block-21 bridge result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing block-21 bridge note: {gate.NOTE}")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_block21_source_to_corrected_model_bridge_pilot_gate", "block-21 bridge kind drift")
    require(artifact["status"] == "representative_block21_signed_source_bridge_closed", "block-21 bridge status drift")
    require(bool(artifact["passed"]), "block-21 bridge did not pass")
    rows = artifact["rows"]
    require([int(row["chain"]) for row in rows] == list(gate.WITNESS_CHAINS), "block-21 bridge witness drift")
    require({(row["selector"], row["phi3_sign"]) for row in rows} == {
        ("W3", "negative"), ("W3", "positive"), ("W4", "negative"), ("W4", "positive")
    }, "block-21 bridge branch coverage drift")
    for row in rows:
        require(int(row["block"]) == 21 and int(row["parent_length"]) == 119, "block-21 bridge roster drift")
        require(bool(row["precision_overlap"]), f"chain {row['chain']} precision overlap missing")
        require(bool(row["signed_bridge_closed"]), f"chain {row['chain']} signed bridge open")
        require(bool(row["exact_correction_within_retained_majorant"]), f"chain {row['chain']} retained majorant failed")
        require(Fraction(row["selected_tail_bound_upper"]) < gate.TAIL_TARGET, f"chain {row['chain']} tail target failed")
        require(Fraction(row["formula_target_gap_upper"]) < gate.FORMULA_GAP_REPORTING_CEILING, f"chain {row['chain']} formula target gap failed")
        for name in ("w1_identity_gap_upper", "t2_identity_gap_upper", "corrected_identity_gap_upper", "endpoint_closure_gap_upper"):
            require(Fraction(row[name]) < gate.FORMULA_GAP_REPORTING_CEILING, f"chain {row['chain']} {name} failed")
    aggregate = artifact["aggregate"]
    require(int(aggregate["witness_count"]) == len(gate.WITNESS_CHAINS), "block-21 bridge aggregate count drift")
    require(int(aggregate["closed_witness_count"]) == len(gate.WITNESS_CHAINS), "block-21 bridge aggregate closure drift")
    require(int(aggregate["exact_corrections_within_retained_count"]) == len(gate.WITNESS_CHAINS), "block-21 bridge retained count drift")
    require(Fraction(aggregate["maximum_formula_target_gap_upper"]) == max(Fraction(row["formula_target_gap_upper"]) for row in rows), "block-21 formula maximum drift")
    provenance = artifact["provenance"]
    for name, expected in (
        ("telemetry", gate.TELEMETRY),
        ("endpoints", gate.ENDPOINTS),
        ("weights", gate.WEIGHTS),
        ("builder", gate.BUILDER),
        ("checker", gate.CHECKER),
    ):
        require(provenance[name]["sha256"] == file_hash(expected), f"block-21 bridge {name} hash drift")
    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("Dsrc", "Dcorr", "Qexact", "parent length 119", "not a proof"):
        require(token in note, f"block-21 bridge note token missing: {token}")
    print(
        "validated block-21 signed bridge pilot: "
        f"witnesses={len(rows)}, max_formula_gap={aggregate['maximum_formula_target_gap_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
