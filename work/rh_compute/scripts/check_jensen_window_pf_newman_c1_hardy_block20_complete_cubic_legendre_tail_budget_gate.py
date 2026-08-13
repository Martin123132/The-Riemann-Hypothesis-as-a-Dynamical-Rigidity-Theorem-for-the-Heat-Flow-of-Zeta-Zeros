#!/usr/bin/env python3
"""Validate the complete saved-point cubic Legendre-tail budget."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate as gate
import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as transport
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q
from flint import arb, ctx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_outward_close(exact: Fraction, printed: Fraction, message: str) -> None:
    require(printed >= exact, f"{message}: not outward")
    scale = max(abs(exact), Fraction(1, 10**100))
    require(printed - exact <= scale / 10**70, f"{message}: excessive slack")


def main() -> int:
    require(gate.RESULT.is_file(), f"missing result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing note: {gate.NOTE}")
    result = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(result["kind"] == "jensen_window_pf_newman_c1_hardy_block20_complete_cubic_legendre_tail_budget_gate", "kind drift")
    require(result["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash drift")
    require(result["sources"]["checker"]["sha256"] == gate.file_hash(Path(__file__).resolve()), "checker hash drift")
    for dependency in result["dependencies"].values():
        path = gate.REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(dependency["sha256"] == gate.file_hash(path), f"dependency hash drift: {path}")

    rows = result["rows"]
    require(len(rows) == 374, "recursive row count drift")
    require(sorted(int(row["chain"]) for row in rows) == sorted(set(int(row["chain"]) for row in rows)), "duplicate chain")
    require(sum(int(row["child_length"]) == 1 for row in rows) == 369, "length-one count drift")
    require(sum(int(row["child_length"]) == 2 for row in rows) == 5, "length-two count drift")
    require(sum(len(row["terms"]) for row in rows) == 753, "child-index count drift")
    row_by_chain = {int(row["chain"]): row for row in rows}

    for row in rows:
        require(bool(row["precision_overlap"]), f"chain {row['chain']} precision failure")
        require(Fraction(row["minimum_stationary_discriminant_lower"]) > 0, "nonpositive discriminant")
        require(Fraction(row["minimum_exact_fourth_root_denominator_lower"]) > 0, "nonpositive exact denominator")
        require(Fraction(row["minimum_source_taylor_denominator_lower"]) > 0, "nonpositive source denominator")
        direct_sum = sum(Fraction(term["direct_term_gap_abs_upper"]) for term in row["terms"])
        elementary_sum = sum(Fraction(term["elementary_term_majorant_upper"]) for term in row["terms"])
        stored_direct = Fraction(row["child_termwise_direct_gap_upper"])
        stored_elementary = Fraction(row["child_termwise_elementary_majorant_upper"])
        require(abs(direct_sum - stored_direct) <= direct_sum / 10**70, "direct child serialization drift")
        require(abs(elementary_sum - stored_elementary) <= elementary_sum / 10**70, "elementary child serialization drift")
        for term in row["terms"]:
            require(Fraction(term["direct_term_gap_abs_upper"]) <= Fraction(term["elementary_term_majorant_upper"]), "term escaped elementary bound")
            require(Fraction(term["discriminant_lower"]) > 0, "term discriminant failure")
            require(Fraction(term["exact_denominator_lower"]) > 0, "term exact denominator failure")
            require(Fraction(term["source_denominator_lower"]) > 0, "term source denominator failure")
        multiplier = Fraction(row["multiplier_abs_upper"])
        stored_local = Fraction(row["local_recurrence_tail_majorant_upper"])
        stored_local_elementary = Fraction(row["local_recurrence_tail_elementary_majorant_upper"])
        require(abs(multiplier * direct_sum - stored_local) <= multiplier * direct_sum / 10**70, "local direct serialization drift")
        require(abs(multiplier * elementary_sum - stored_local_elementary) <= multiplier * elementary_sum / 10**70, "local elementary serialization drift")
        actual = native_q.complex_from_record(row["local_recurrence_tail_actual"])
        ctx.dps = 120
        require(gate.upper_abs(actual) <= Fraction(row["local_recurrence_tail_majorant_upper"]), "actual local tail escaped")

    coefficient_result = json.loads(gate.COEFFICIENT_TRANSPORT.read_text(encoding="utf-8"))
    endpoints = {int(row["output_index"]): row for row in coefficient_result["transported_outputs"]}
    ctx.dps = gate.WEIGHT_PRECISION
    ctx.threads = 1
    weights = transport.load_weights()
    outputs = result["transported_outputs"]
    require(len(outputs) == 15, "output count drift")
    for output in outputs:
        output_index = int(output["output_index"])
        recomputed_tail = Fraction(0)
        for row in rows:
            weight = weights[(int(row["sum_index"]), output_index, int(row["branch"]))]
            recomputed_tail += (
                weight["amplitude_fraction"]
                * weight["phase_abs_upper"]
                * Fraction(row["local_recurrence_tail_majorant_upper"])
            )
        endpoint = Fraction(endpoints[output_index]["cell_endpoint_majorant_upper"])
        combined = endpoint + recomputed_tail
        require_outward_close(recomputed_tail, Fraction(output["complete_cubic_legendre_tail_majorant_upper"]), "output tail transport drift")
        require_outward_close(endpoint, Fraction(output["coefficient_cell_endpoint_majorant_upper"]), "endpoint dependency drift")
        require_outward_close(combined, Fraction(output["endpoint_plus_legendre_tail_majorant_upper"]), "combined output drift")
        require(bool(output["combined_within_requested_scale"]) == (combined < gate.REQUESTED_ERROR_SCALE), "output closure flag drift")

    aggregate = result["aggregate"]
    require(int(aggregate["recursive_call_count"]) == 374, "aggregate call count drift")
    require(int(aggregate["child_index_count"]) == 753, "aggregate child-index count drift")
    require(int(aggregate["precision_overlap_count"]) == 374, "aggregate overlap drift")
    closures = sum(bool(row["combined_within_requested_scale"]) for row in outputs)
    require(int(aggregate["outputs_with_combined_majorant_within_requested_scale"]) == closures, "aggregate closure count drift")
    expected_status = (
        "complete_saved_point_cubic_legendre_tail_and_endpoint_budget_close"
        if closures == 15
        else "complete_saved_point_cubic_legendre_tail_breaks_endpoint_budget"
    )
    require(result["status"] == expected_status, "status drift")
    require("not a proof" in gate.NOTE.read_text(encoding="utf-8").lower(), "note proof boundary missing")
    require("RH" in result["proof_boundary"], "JSON proof boundary missing")
    print(
        "validated complete cubic Legendre-tail budget: "
        f"374 calls, 753 child indices, {closures}/15 combined outputs below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
