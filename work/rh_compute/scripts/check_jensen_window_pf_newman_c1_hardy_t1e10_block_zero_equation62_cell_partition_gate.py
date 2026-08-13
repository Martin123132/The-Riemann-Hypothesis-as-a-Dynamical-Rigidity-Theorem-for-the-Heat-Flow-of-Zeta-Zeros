#!/usr/bin/env python3
"""Independently check the equation-(62) cell-partition artifact."""

from __future__ import annotations

from decimal import Decimal
import hashlib
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb

import jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate as gate


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing cell-partition result")
    require(gate.NOTE.is_file(), "missing cell-partition note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_block_zero_equation62_cell_partition_gate", "kind drift")
    require(artifact["status"] == "finite_diagnostic_nearest_root_midpoint_comparison_validated", "status drift")
    require(bool(artifact["passed"]), "cell-partition gate failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")

    midpoint_abs: list[arb] = []
    source_abs: list[arb] = []
    boundary_abs: list[arb] = []
    all_prefix_abs: list[arb] = []
    half_integral_count = 0
    below_count = 0
    rows = artifact["output_rows"]
    require(len(rows) == gate.EXPECTED_OUTPUTS, "output roster drift")
    for output_index, row in enumerate(rows, start=1):
        offset = Decimal(output_index - 8) / Decimal(100)
        fresh = gate.output_atlas(offset, CHECK_PRECISION)
        require(row["output_index"] == output_index, "output index drift")
        require(row["target_t"] == fresh["target_t"], "target height drift")
        require(len(row["root_certificates"]) == len(gate.ROOT_ALPHAS), "root roster length drift")
        require(len(row["prefix_rows"]) == len(gate.ALPHAS), "prefix roster length drift")
        for saved, current in zip(row["root_certificates"], fresh["root_certificates"]):
            for key in ("alpha", "nearest_integer", "floor_integer"):
                require(saved[key] == current[key], f"root discrete drift: {key}")
            for key in (
                "root_ball",
                "distance_from_lower_half_integer_ball",
                "distance_from_upper_half_integer_ball",
                "distance_from_floor_ball",
                "distance_to_next_integer_ball",
            ):
                require(arb(saved[key]).overlaps(arb(current[key])), f"root ball drift: {key}")

        for saved, current in zip(row["prefix_rows"], fresh["prefix_rows"]):
            for key in (
                "prefix_index",
                "upper_alpha",
                "nearest_integer",
                "next_nearest_integer",
                "midpoint_numerator",
                "midpoint_half_integral",
                "midpoint_floor_start",
                "midpoint_ceil_start",
                "source_endpoint_start",
            ):
                require(saved[key] == current[key], f"prefix discrete drift at {saved['prefix_index']}: {key}")
            for key in (
                "equation62_prefix_sum_ball",
                "midpoint_floor_target_ball",
                "midpoint_ceil_target_ball",
                "source_endpoint_target_ball",
                "prefix_minus_midpoint_floor_target_ball",
                "prefix_minus_midpoint_ceil_target_ball",
                "prefix_minus_source_endpoint_target_ball",
                "midpoint_ceil_target_minus_source_endpoint_target_ball",
            ):
                require(arb(saved[key]).overlaps(arb(current[key])), f"prefix ball drift at {saved['prefix_index']}: {key}")
            decomposition = (
                arb(current["prefix_minus_midpoint_ceil_target_ball"])
                + arb(current["midpoint_ceil_target_minus_source_endpoint_target_ball"])
                - arb(current["prefix_minus_source_endpoint_target_ball"])
            )
            require(decomposition.contains(0), f"fresh boundary decomposition failed at prefix={saved['prefix_index']}")
            value_abs = abs(arb(current["prefix_minus_midpoint_ceil_target_ball"]))
            all_prefix_abs.append(value_abs)
            half_integral_count += int(current["midpoint_half_integral"])
            below_count += int(value_abs < gate.REQUESTED_TOLERANCE)

        final = fresh["prefix_rows"][-1]
        require(final["midpoint_floor_start"] == gate.EXPECTED_FINAL_DIAGNOSTIC_START, "final floor midpoint drift")
        require(final["midpoint_ceil_start"] == gate.EXPECTED_FINAL_DIAGNOSTIC_START, "final ceil midpoint drift")
        require(final["source_endpoint_start"] == gate.EXPECTED_FINAL_SOURCE_START, "final source endpoint drift")
        m_abs = abs(arb(final["prefix_minus_midpoint_ceil_target_ball"]))
        s_abs = abs(arb(final["prefix_minus_source_endpoint_target_ball"]))
        b_abs = abs(arb(final["midpoint_ceil_target_minus_source_endpoint_target_ball"]))
        require(m_abs < gate.REQUESTED_TOLERANCE, f"fresh final midpoint tolerance failed at output={output_index}")
        require(s_abs > gate.REQUESTED_TOLERANCE, f"fresh final source rejection failed at output={output_index}")
        midpoint_abs.append(m_abs)
        source_abs.append(s_abs)
        boundary_abs.append(b_abs)

    aggregate = artifact["aggregate"]
    require(aggregate["output_count"] == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(aggregate["prefixes_per_output"] == len(gate.ALPHAS), "prefix aggregate drift")
    require(aggregate["total_prefix_comparisons"] == gate.EXPECTED_OUTPUTS * len(gate.ALPHAS), "total prefix aggregate drift")
    require(aggregate["total_root_certificates"] == gate.EXPECTED_OUTPUTS * len(gate.ROOT_ALPHAS), "root aggregate drift")
    require(aggregate["half_integral_midpoint_count"] == half_integral_count, "half-integral aggregate drift")
    require(aggregate["midpoint_ceil_prefix_residual_below_0p005_count"] == below_count, "prefix tolerance aggregate drift")
    require(min(all_prefix_abs, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_midpoint_ceil_prefix_residual_absolute_ball"])), "minimum prefix residual drift")
    require(max(all_prefix_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_midpoint_ceil_prefix_residual_absolute_ball"])), "maximum prefix residual drift")
    require(min(midpoint_abs, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_final_midpoint_residual_absolute_ball"])), "minimum final midpoint residual drift")
    require(max(midpoint_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_final_midpoint_residual_absolute_ball"])), "maximum final midpoint residual drift")
    require(min(source_abs, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_final_source_endpoint_residual_absolute_ball"])), "minimum final source residual drift")
    require(max(source_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_final_source_endpoint_residual_absolute_ball"])), "maximum final source residual drift")
    require(min(boundary_abs, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_final_boundary_correction_absolute_ball"])), "minimum boundary correction drift")
    require(max(boundary_abs, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_final_boundary_correction_absolute_ball"])), "maximum boundary correction drift")
    require(aggregate["final_boundary_classical_terms"] == list(range(37941, 37946)), "boundary term roster drift")

    decision = artifact["decision"]
    require(bool(decision["diagnostic_midpoint_rule_requires_lattice_completion_at_some_prefixes"]), "midpoint ambiguity hidden")
    require(not bool(decision["final_101_term_midpoint_is_ambiguous"]), "final midpoint falsely ambiguous")
    require(not bool(decision["source_endpoint_and_diagnostic_midpoint_agree_at_final_prefix"]), "endpoint comparison hidden")
    require(bool(decision["final_source_endpoint_residual_rejects_0p005_at_all_outputs"]), "source rejection lost")
    require(bool(decision["final_diagnostic_midpoint_residual_satisfies_0p005_at_all_outputs"]), "diagnostic midpoint result lost")
    require(bool(decision["five_term_cutoff_change_dominates_source_aligned_residual"]), "cutoff comparison lost")
    require(not bool(decision["diagnostic_midpoint_is_a_published_cutoff_rule"]), "diagnostic midpoint promoted to source rule")
    require(not bool(decision["height_uniform_midpoint_remainder_proved"]), "height-uniform overpromotion")
    require(not bool(decision["finite_gate_has_rh_implication"]), "RH overpromotion")

    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file: {path}")
            require(file_hash(path) == record["sha256"], f"hash drift: {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "diagnostic alternative", "39771.5", "37941", "37946", "37941..37945", "no RH implication"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated equation-(62) cell partition independently: "
        f"roots={aggregate['total_root_certificates']}, prefixes={aggregate['total_prefix_comparisons']}, "
        f"final-midpoint<0.005={len(midpoint_abs)}/{gate.EXPECTED_OUTPUTS}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
