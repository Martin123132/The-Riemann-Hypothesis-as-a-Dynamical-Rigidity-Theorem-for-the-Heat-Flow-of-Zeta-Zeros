#!/usr/bin/env python3
"""Validate the block-20 t5 intrinsic real-erfc residual gate."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.md"
BUILDER = REPO_ROOT / "work/rh_compute/scripts/jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate.py"
CHECKER = Path(__file__).resolve()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    for path in (RESULT, NOTE, BUILDER, CHECKER):
        require(path.is_file(), f"missing t5-erfc residual artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    note = NOTE.read_text(encoding="utf-8")
    require(
        artifact["kind"] == "jensen_window_pf_newman_c1_hardy_block20_t5_real_erfc_residual_gate",
        "t5-erfc residual kind drift",
    )
    require(
        artifact["status"]
        == "rigorous_913_exact_argument_intrinsic_erfc_residuals_and_correlated_q_replacement_enclosed",
        "t5-erfc residual status drift",
    )
    scope = artifact["scope"]
    require(scope["recursive_chain_count"] == 374, "t5-erfc chain count drift")
    require(scope["lower_erfc_evaluation_count"] == 539, "t5-erfc lower count drift")
    require(scope["upper_erfc_evaluation_count"] == 374, "t5-erfc upper count drift")
    require(scope["intrinsic_real_erfc_evaluation_count"] == 913, "t5-erfc total count drift")
    require(scope["saved_t5_bit_replay_count"] == 374, "t5-erfc replay count drift")
    require(scope["precision_decimal_digits"] >= 180, "t5-erfc precision too low")

    contract = artifact["comparison_contract"]
    require("exact emitted binary128" in contract["reference"], "t5-erfc exact-argument contract drift")
    require("reported separately" in contract["separation"], "t5-erfc association separation drift")
    require("introduces no independent pi" in artifact["pi_provenance"]["rigorous_reference"], "t5-erfc pi provenance drift")

    aggregate = artifact["aggregate"]
    require(aggregate["region_histogram"] == {"lower": 539, "upper": 374}, "t5-erfc region histogram drift")
    require(sum(aggregate["argument_bucket_histogram"].values()) == 913, "t5-erfc argument histogram drift")
    require(0 <= aggregate["source_output_enclosed_count"] <= 913, "t5-erfc enclosure count drift")
    require(0 <= aggregate["source_output_zero_count"] <= 913, "t5-erfc zero count drift")
    for key in (
        "maximum_erfc_residual_abs_upper",
        "maximum_weight_abs_upper",
        "maximum_weighted_call_delta_abs_upper",
        "maximum_total_q_delta_abs_upper",
        "maximum_source_term_association_gap_abs_upper",
        "maximum_correlated_source_term_association_gap_abs_upper",
    ):
        require(Decimal(aggregate[key]) >= 0, f"negative t5-erfc bound: {key}")
    require(Decimal(aggregate["maximum_weight_abs_upper"]) > 0, "zero t5-erfc weight drift")

    rows = artifact["rows"]
    require(len(rows) == 374, "t5-erfc row count drift")
    require([row["chain"] for row in rows] == sorted(row["chain"] for row in rows), "t5-erfc order drift")
    require(sum(len(row["calls"]) for row in rows) == 913, "t5-erfc nested call count drift")
    for row in rows:
        require(len(row["calls"]) in (2, 3), "t5-erfc per-chain roster drift")
        require(len(row["saved_t5_hex"]) == 2, "t5-erfc saved t5 payload drift")
        replacement = row["correlated_replacement"]
        require("lower_dyadic" in replacement["total_q_delta"]["real"], "t5-erfc real dyadic interval missing")
        require("upper_dyadic" in replacement["total_q_delta"]["imag"], "t5-erfc imag dyadic interval missing")

    for key in ("probe_result", "probe_output", "selector_atlas", "local_cell_builder"):
        item = artifact["sources"][key]
        require(file_hash(REPO_ROOT / item["path"]) == item["sha256"], f"t5-erfc {key} hash drift")
    require(file_hash(BUILDER) == artifact["sources"]["builder"]["sha256"], "t5-erfc builder hash drift")
    require(file_hash(CHECKER) == artifact["sources"]["checker"]["sha256"], "t5-erfc checker hash drift")
    require(file_hash(REPO_ROOT / artifact["source"]["path"]) == artifact["source"]["sha256"], "t5-erfc source drift")

    for token in (
        "Status: rigorous finite exact-argument intrinsic-erfc residuals",
        "913 total intrinsic real-erfc evaluations",
        "inserts no additional pi",
        "source-term association gap",
        "maximum correlated t5/q replacement",
        "does not validate the saddle locations",
        "prize-level theorem",
    ):
        require(token in note, f"t5-erfc note token missing: {token}")
    require("913 exact finite source arguments" in artifact["proof_boundary"], "t5-erfc boundary drift")
    print(
        "validated Hardy block-20 t5 real-erfc residual gate: "
        f"913 calls, max residual {aggregate['maximum_erfc_residual_abs_upper']}, "
        f"max correlated q shift {aggregate['maximum_total_q_delta_abs_upper']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
