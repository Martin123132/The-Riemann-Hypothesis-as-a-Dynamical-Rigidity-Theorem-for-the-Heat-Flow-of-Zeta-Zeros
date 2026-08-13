#!/usr/bin/env python3
"""Independently check the finite source-aligned Hardy hybrid-error ledger."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
sys.path.insert(0, str(VENDOR))
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")

from flint import arb, ctx


RESULT_ROOT = REPO_ROOT / "work/rh_compute/results"
RESULT = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_t1e10_source_aligned_hybrid_error_ledger_gate.md"
CALIBRATION = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_exact_hardy_arb_calibration_gate.json"
COMPONENT_SPLIT = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_lower_rs_hybrid_upper_component_split_gate.json"
CLASSICAL_SPLIT = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_classical_upper_main_remainder_split_gate.json"
BLOCK_PARTITION = RESULT_ROOT / "jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate.json"
SOURCE = REPO_ROOT / "work/rh_compute/external/hardy_fastcode/telemetry/zeta14cubicmult_telemetry.f90"
EXPECTED_OUTPUTS = 15
EXPECTED_BLOCKS = 36
FIRST_OPEN_ID = "source_aligned_upper_block_sum"
TOLERANCE = arb("0.005")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT.resolve()).as_posix()


def minimum(values: list[arb]) -> arb:
    return min(values, key=lambda value: value.lower())


def maximum(values: list[arb]) -> arb:
    return max(values, key=lambda value: value.upper())


def main() -> None:
    ctx.dps = 130
    ctx.threads = 1
    artifact = load(RESULT)
    calibration = load(CALIBRATION)
    component = load(COMPONENT_SPLIT)
    classical = load(CLASSICAL_SPLIT)
    partition = load(BLOCK_PARTITION)

    require(artifact["passed"] is True, "ledger is not passed")
    require(artifact["scope"]["diagnostic_midpoint_used"] is False, "midpoint entered primary ledger")
    require(artifact["decision"]["source_endpoint_is_primary"] is True, "source endpoint decision drift")
    require(artifact["decision"]["first_unbounded_signed_combination_id"] == FIRST_OPEN_ID, "first open target drift")
    require(artifact["decision"]["corrected_internal_model_closes_source_aligned_outer_error"] is False, "internal model promotion drift")
    require(artifact["paper_audit"]["paper_declares_hybrid_not_exact"] is True, "paper exactness drift")
    require(artifact["source_aligned_identity"]["first_tail_start_n_at_center"] == 37946, "endpoint cutoff drift")
    require(len(artifact["output_rows"]) == EXPECTED_OUTPUTS, "ledger output count drift")

    source_lines = SOURCE.read_text(encoding="utf-8", errors="replace").splitlines()
    for audit in artifact["source_audit"]:
        excerpt = "\n".join(source_lines[int(audit["line_start"]) - 1 : int(audit["line_end"])])
        require(all(token in excerpt for token in audit["required_tokens"]), f"source token drift: {audit['id']}")
        require(hashlib.sha256(excerpt.encode("utf-8")).hexdigest() == audit["excerpt_sha256"], f"source hash drift: {audit['id']}")

    totals: list[arb] = []
    block_sums: list[arb] = []
    triangles: list[arb] = []
    retained_values: list[arb] = []
    cancellation_values: list[arb] = []
    subordinate_values: list[arb] = []

    for index in range(EXPECTED_OUTPUTS):
        saved = artifact["output_rows"][index]
        cal_row = calibration["output_rows"][index]
        component_row = component["output_rows"][index]
        classical_row = classical["output_rows"][index]
        partition_row = partition["output_rows"][index]
        require(saved["output_index"] == index + 1, f"ledger row order drift at {index + 1}")
        require(all(arb(saved["target_t"]).overlaps(arb(row["target_t"])) for row in (cal_row, component_row, classical_row, partition_row)), f"target mismatch at {index + 1}")

        total = arb(cal_row["signed_source_minus_hardy_ball"])
        lower = arb(component_row["lower_source_minus_exact_ball"])
        upper = arb(component_row["hybrid_upper_source_minus_exact_ball"])
        addition = arb(component_row["final_addition_rounding_ball"])
        main_gap = arb(classical_row["hybrid_minus_classical_main_ball"])
        correction = arb(classical_row["exact_remaining_correction_ball"])
        increments = [arb(row["increment_residual_ball"]) for row in partition_row["block_rows"]]
        require(len(increments) == EXPECTED_BLOCKS, f"block count drift at {index + 1}")
        block_sum = sum(increments, arb(0))
        triangle = sum((abs(value) for value in increments), arb(0))
        later_net = sum(increments[1:], arb(0))
        retained = abs(block_sum) / triangle
        cancellation = 1 - retained
        subordinate = abs(lower) + abs(correction) + abs(addition)

        require((lower + upper + addition - total).contains(0), f"component identity failed at {index + 1}")
        require((main_gap - correction - upper).contains(0), f"upper identity failed at {index + 1}")
        require((block_sum - main_gap).contains(0), f"block identity failed at {index + 1}")
        require((lower + block_sum - correction + addition - total).contains(0), f"total identity failed at {index + 1}")
        require(total < 0 and block_sum < 0 and increments[0] < 0 and later_net > 0, f"sign pattern drift at {index + 1}")
        require(triangle > TOLERANCE and 0 < retained < 1 and 0 < cancellation < 1, f"cancellation metrics drift at {index + 1}")

        comparisons = {
            "total_source_minus_exact_ball": total,
            "lower_source_minus_exact_ball": lower,
            "hybrid_upper_source_minus_exact_ball": upper,
            "source_aligned_upper_block_sum_ball": block_sum,
            "exact_upper_remaining_correction_ball": correction,
            "final_addition_rounding_ball": addition,
            "block_zero_increment_residual_ball": increments[0],
            "later_blocks_net_increment_residual_ball": later_net,
            "independent_block_triangle_ball": triangle,
            "retained_fraction_ball": retained,
            "cancellation_fraction_ball": cancellation,
            "subordinate_channel_triangle_ball": subordinate,
        }
        for field, recomputed in comparisons.items():
            require(arb(saved[field]).overlaps(recomputed), f"saved field misses recomputation: output {index + 1}, {field}")

        totals.append(abs(total))
        block_sums.append(abs(block_sum))
        triangles.append(triangle)
        retained_values.append(retained)
        cancellation_values.append(cancellation)
        subordinate_values.append(subordinate)

    aggregate_checks = {
        "minimum_total_absolute_error_ball": minimum(totals),
        "maximum_total_absolute_error_ball": maximum(totals),
        "minimum_source_aligned_upper_block_sum_absolute_ball": minimum(block_sums),
        "maximum_source_aligned_upper_block_sum_absolute_ball": maximum(block_sums),
        "minimum_independent_block_triangle_ball": minimum(triangles),
        "maximum_independent_block_triangle_ball": maximum(triangles),
        "minimum_retained_fraction_ball": minimum(retained_values),
        "maximum_retained_fraction_ball": maximum(retained_values),
        "minimum_cancellation_fraction_ball": minimum(cancellation_values),
        "maximum_cancellation_fraction_ball": maximum(cancellation_values),
        "maximum_subordinate_channel_triangle_ball": maximum(subordinate_values),
    }
    for field, recomputed in aggregate_checks.items():
        require(arb(artifact["aggregate"][field]).overlaps(recomputed), f"aggregate drift: {field}")

    require(artifact["aggregate"]["block_row_count"] == EXPECTED_OUTPUTS * EXPECTED_BLOCKS, "block row aggregate drift")
    require(artifact["aggregate"]["total_identity_contains_zero_count"] == EXPECTED_OUTPUTS, "total identity count drift")
    require(artifact["aggregate"]["block_triangle_rejects_0p005_count"] == EXPECTED_OUTPUTS, "triangle count drift")
    channels = {row["id"]: row for row in artifact["channel_ledger"]}
    require(channels[FIRST_OPEN_ID]["state"] == "open_height_uniform_signed_theorem", "open target state drift")
    require(channels["corrected_internal_model_column"]["state"] == "interval_certified_but_not_source_aligned_outer_error", "internal-model boundary drift")

    for dependency in artifact["dependencies"].values():
        path = REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {dependency['path']}")
        require(file_hash(path) == dependency["sha256"], f"dependency hash drift: {dependency['path']}")
    for source in artifact["sources"].values():
        path = REPO_ROOT / source["path"]
        require(file_hash(path) == source["sha256"], f"source hash drift: {source['path']}")

    note = NOTE.read_text(encoding="utf-8")
    for token in (
        "Z_src-Z_exact",
        "E_hyb_main(t)=sum_k e_k(t)=Q_src(t)-T_upper(t)",
        "diagnostic nearest-root midpoint is not used",
        "no RH",
    ):
        require(token in note, f"note token missing: {token}")

    print(
        "validated source-aligned hybrid error ledger independently: "
        f"outputs={EXPECTED_OUTPUTS}, blocks={EXPECTED_OUTPUTS * EXPECTED_BLOCKS}, open-first={FIRST_OPEN_ID}"
    )


if __name__ == "__main__":
    main()
