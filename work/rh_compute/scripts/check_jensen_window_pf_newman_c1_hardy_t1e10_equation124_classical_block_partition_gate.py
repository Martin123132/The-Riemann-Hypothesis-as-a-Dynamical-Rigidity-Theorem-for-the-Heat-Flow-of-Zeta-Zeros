#!/usr/bin/env python3
"""Independently validate the equation-(124) classical block partition."""

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

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate as gate


CHECK_PRECISION = 130


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    require(gate.RESULT.is_file(), "missing equation-124 partition result")
    require(gate.NOTE.is_file(), "missing equation-124 partition note")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact["kind"] == "jensen_window_pf_newman_c1_hardy_t1e10_equation124_classical_block_partition_gate", "kind drift")
    require(artifact["status"] == "rigorous_fifteen_output_thirty_six_block_equation124_induced_classical_partition_atlas", "status drift")
    require(bool(artifact["passed"]), "partition gate failed")
    require(artifact["scope"]["precisions_decimal_digits"] == list(gate.PRECISIONS), "precision roster drift")
    gate.require_source_cutoff_semantics()
    partition = artifact["partition_definition"]
    require(bool(partition["source_rn1_endpoint_used_without_cell_edge_shift"]), "source endpoint convention drift")
    require("NC=AINT(CE)" in partition["source_code_rule"], "source AINT rule missing")
    require("N_k=NC+1" in partition["source_code_rule"], "source tail-start rule missing")

    source_blocks = gate.load_stage_zero()
    outputs = artifact["output_rows"]
    require(len(outputs) == gate.EXPECTED_OUTPUTS, "output roster drift")
    canonical_cutoffs = artifact["aggregate"]["cutoff_roster"]
    require(len(canonical_cutoffs) == gate.EXPECTED_BLOCKS, "cutoff roster length drift")
    require(canonical_cutoffs[-1] == gate.LOWER_START, "final cutoff drift")
    require(sum(row["exact_block_term_count"] for row in outputs[0]["block_rows"]) == gate.UPPER_END - gate.LOWER_START + 1, "partition term count drift")

    fresh_block_zero: list[arb] = []
    fresh_final: list[arb] = []
    fresh_later: list[arb] = []
    for output_index, output in enumerate(outputs, start=1):
        offset = Decimal(output_index - 8) / Decimal(100)
        target_t = format(Decimal("1e10") + offset, "f")
        require(output["output_index"] == output_index, "output index drift")
        require(output["output_label"] == f"t{offset:+.2f}", "output label drift")
        require(output["target_t"] == target_t, "target height drift")
        block_rows = output["block_rows"]
        require(len(block_rows) == gate.EXPECTED_BLOCKS, "block roster drift")

        fresh_certs = [gate.cutoff_certificate(target_t, row["rn1"], CHECK_PRECISION) for row in source_blocks]
        cutoffs = [int(cert["cutoff_n"]) for cert in fresh_certs]
        require(cutoffs == canonical_cutoffs, f"fresh cutoff roster drift at output {output_index}")
        targets = gate.cumulative_targets(target_t, cutoffs, CHECK_PRECISION)
        ctx.dps = CHECK_PRECISION
        previous_source = arb(0)
        previous_target = arb(0)
        previous_residual = arb(0)
        for block, (saved, source_row, cert, cutoff) in enumerate(zip(block_rows, source_blocks, fresh_certs, cutoffs)):
            require(saved["block"] == block, "block index drift")
            require(saved["cutoff_n"] == cutoff, "saved cutoff drift")
            require(cert["source_nc_aint_positive"] + 1 == cutoff, "source NC/tail boundary drift")
            require(cert["source_tail_start_n"] == cutoff, "source tail start drift")
            require(arb(cert["source_ce_minus_nc_ball"]) > 0, "source CE lower floor margin lost")
            require(arb(cert["source_nc_plus_one_minus_ce_ball"]) > 0, "source CE upper floor margin lost")
            require(arb(cert["rn1_minus_alpha_at_cutoff_ball"]) > 0, "lower cutoff margin lost")
            require(arb(cert["alpha_before_cutoff_minus_rn1_ball"]) > 0, "upper cutoff margin lost")
            target = targets[cutoff]
            require(target.overlaps(arb(saved["high_precision"]["exact_cumulative_target_ball"])), f"fresh target misses saved ball at output {output_index}, block {block}")
            source = arb(source_row["zsum"][output_index - 1])
            residual = source - target
            source_increment = source - previous_source
            target_increment = target - previous_target
            increment_residual = source_increment - target_increment
            identity = increment_residual - (residual - previous_residual)
            require(identity.contains(0), f"fresh increment identity excludes zero at output {output_index}, block {block}")
            require(residual.overlaps(arb(saved["cumulative_residual_ball"])), f"fresh cumulative residual misses saved ball at output {output_index}, block {block}")
            require(increment_residual.overlaps(arb(saved["increment_residual_ball"])), f"fresh increment residual misses saved ball at output {output_index}, block {block}")
            require(gate.sign_name(residual) == saved["cumulative_residual_sign"], "cumulative sign drift")
            require(gate.sign_name(increment_residual) == saved["increment_residual_sign"], "increment sign drift")
            previous_source = source
            previous_target = target
            previous_residual = residual

        block_zero = arb(block_rows[0]["cumulative_residual_ball"])
        final_residual = arb(block_rows[-1]["cumulative_residual_ball"])
        later = final_residual - block_zero
        require(block_zero < -arb("0.005"), f"fresh block-zero tolerance rejection drift at output {output_index}")
        require(later > 0, f"fresh later cancellation drift at output {output_index}")
        require(output["first_nonzero_cumulative_residual_block"] == 0, "earliest cumulative block drift")
        require(output["first_nonzero_increment_residual_block"] == 0, "earliest increment block drift")
        require(later.overlaps(arb(output["later_blocks_net_residual_ball"])), "later net residual drift")
        fresh_block_zero.append(abs(block_zero))
        fresh_final.append(abs(final_residual))
        fresh_later.append(abs(later))

    aggregate = artifact["aggregate"]
    require(aggregate["output_count"] == gate.EXPECTED_OUTPUTS, "output aggregate drift")
    require(aggregate["block_count"] == gate.EXPECTED_BLOCKS, "block aggregate drift")
    require(aggregate["row_count"] == gate.EXPECTED_OUTPUTS * gate.EXPECTED_BLOCKS, "row aggregate drift")
    require(aggregate["precision_overlap_count"] == aggregate["row_count"], "precision aggregate drift")
    require(aggregate["cutoff_inequality_certificate_count"] == aggregate["row_count"], "cutoff aggregate drift")
    require(aggregate["source_ce_floor_certificate_count"] == aggregate["row_count"], "source CE aggregate drift")
    require(aggregate["source_cutoff_semantics_markers_verified"] == len(gate.SOURCE_CUTOFF_MARKERS), "source marker aggregate drift")
    require(aggregate["partition_total_terms"] == gate.UPPER_END - gate.LOWER_START + 1, "partition aggregate drift")
    require(aggregate["block_zero_rejects_0p005_count"] == gate.EXPECTED_OUTPUTS, "block-zero aggregate drift")
    require(aggregate["later_blocks_net_positive_count"] == gate.EXPECTED_OUTPUTS, "later cancellation aggregate drift")
    require(aggregate["final_residual_matches_classical_split_count"] == gate.EXPECTED_OUTPUTS, "final comparison aggregate drift")
    require(min(fresh_block_zero, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_block_zero_residual_absolute_ball"])), "minimum block-zero aggregate drift")
    require(max(fresh_block_zero, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_block_zero_residual_absolute_ball"])), "maximum block-zero aggregate drift")
    require(min(fresh_later, key=lambda x: x.lower()).overlaps(arb(aggregate["minimum_later_net_residual_absolute_ball"])), "minimum later aggregate drift")
    require(max(fresh_final, key=lambda x: x.upper()).overlaps(arb(aggregate["maximum_final_residual_absolute_ball"])), "maximum final aggregate drift")

    decision = artifact["decision"]
    require(bool(decision["source_discrete_boundary_convention_resolved"]), "source boundary convention unresolved")
    require(decision["earliest_block_aligned_discrepancy"] == 0, "earliest block decision drift")
    require(bool(decision["block_zero_alone_rejects_saved_tolerance_on_all_outputs"]), "block-zero decision lost")
    require(bool(decision["blocks_1_through_35_net_partially_cancel_block_zero_on_all_outputs"]), "later cancellation decision lost")
    require(not bool(decision["finite_atlas_is_height_uniform_theorem"]), "height-uniform overpromotion")
    require(not bool(decision["finite_atlas_has_rh_implication"]), "RH overpromotion")
    for section in ("dependencies", "sources"):
        for record in artifact[section].values():
            path = REPO_ROOT / record["path"]
            require(path.is_file(), f"missing hashed file {path}")
            require(file_hash(path) == record["sha256"], f"hash drift {path}")

    note = gate.NOTE.read_text(encoding="utf-8")
    for token in ("not a proof of RH", "NC=AINT(CE)", "No `RN1+1`", "Block zero", "partially cancels", "not exact", "no RH implication"):
        require(token in note, f"note boundary token missing: {token}")
    print(
        "validated equation-124 block partition: "
        f"{aggregate['row_count']} rows, block0-rejects-0.005={aggregate['block_zero_rejects_0p005_count']}, "
        f"later-net-cancels={aggregate['later_blocks_net_positive_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
