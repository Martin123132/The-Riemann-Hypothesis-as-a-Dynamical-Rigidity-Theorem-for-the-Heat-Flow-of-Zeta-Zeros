#!/usr/bin/env python3
"""Validate coefficient-cell transport of the complete cubic Legendre tail."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate as gate
import jensen_window_pf_newman_c1_hardy_block20_transported_endpoint_budget_gate as transport
from flint import ctx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def outward_close(exact: Fraction, printed: Fraction, message: str) -> None:
    require(printed >= exact, f"{message}: not outward")
    require(printed - exact <= max(abs(exact), Fraction(1, 10**100)) / 10**70, f"{message}: excessive slack")


def main() -> int:
    require(gate.RESULT.is_file(), f"missing result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing note: {gate.NOTE}")
    result = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(result["kind"] == "jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_complete_cubic_legendre_tail_budget_gate", "kind drift")
    require(result["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash drift")
    require(result["sources"]["checker"]["sha256"] == gate.file_hash(Path(__file__).resolve()), "checker hash drift")
    for dependency in result["dependencies"].values():
        path = gate.REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(dependency["sha256"] == gate.file_hash(path), f"dependency hash drift: {path}")
    require(gate.algebra_audit() == 10, "exact algebra audit drift")

    rows = result["rows"]
    require(len(rows) == 374, "cell count drift")
    require(len({int(row["chain"]) for row in rows}) == 374, "duplicate chain")
    require(sum(len(row["term_gap_abs_uppers"]) for row in rows) == 753, "child-index count drift")
    for row in rows:
        require(bool(row["precision_overlap"]), "precision overlap failure")
        require(Fraction(row["minimum_stationary_discriminant_lower"]) > 0, "nonpositive discriminant")
        require(Fraction(row["minimum_exact_fourth_root_denominator_lower"]) > 0, "nonpositive exact denominator")
        require(Fraction(row["minimum_source_taylor_denominator_lower"]) > 0, "nonpositive source denominator")
        require(Fraction(row["maximum_factored_legendre_phase_tail_upper"]) >= 0, "negative phase tail")
        require(Fraction(row["local_point_tail_majorant_upper"]) <= Fraction(row["local_cell_tail_majorant_upper"]), "point tail escaped cell")
        direct_sum = sum(Fraction(value) for value in row["term_gap_abs_uppers"])
        stored_child = Fraction(row["child_cell_termwise_direct_gap_upper"])
        require(abs(direct_sum - stored_child) <= direct_sum / 10**70, "term serialization drift")
        multiplier = Fraction(row["multiplier_cell_abs_upper"])
        stored_local = Fraction(row["local_cell_tail_majorant_upper"])
        require(abs(multiplier * direct_sum - stored_local) <= multiplier * direct_sum / 10**70, "local transport drift")

    endpoint_result = json.loads(gate.ENDPOINT_CELLS.read_text(encoding="utf-8"))
    endpoints = {int(row["output_index"]): row for row in endpoint_result["transported_outputs"]}
    ctx.dps = gate.WEIGHT_PRECISION
    ctx.threads = 1
    weights = transport.load_weights()
    outputs = result["transported_outputs"]
    require(len(outputs) == 15, "output count drift")
    for output in outputs:
        output_index = int(output["output_index"])
        tail = Fraction(0)
        for row in rows:
            weight = weights[(int(row["sum_index"]), output_index, int(row["branch"]))]
            tail += weight["amplitude_fraction"] * weight["phase_abs_upper"] * Fraction(row["local_cell_tail_majorant_upper"])
        endpoint = Fraction(endpoints[output_index]["cell_endpoint_majorant_upper"])
        combined = endpoint + tail
        outward_close(tail, Fraction(output["coefficient_cell_legendre_tail_majorant_upper"]), "output tail drift")
        outward_close(endpoint, Fraction(output["coefficient_cell_endpoint_majorant_upper"]), "endpoint drift")
        outward_close(combined, Fraction(output["combined_cell_majorant_upper"]), "combined drift")
        require(bool(output["combined_within_requested_scale"]) == (combined < gate.REQUESTED_ERROR_SCALE), "closure flag drift")

    aggregate = result["aggregate"]
    require(int(aggregate["exact_algebra_checks"]) == 10, "aggregate algebra count drift")
    require(int(aggregate["cell_count"]) == 374, "aggregate cell count drift")
    require(int(aggregate["child_index_count"]) == 753, "aggregate child count drift")
    closures = sum(bool(row["combined_within_requested_scale"]) for row in outputs)
    require(int(aggregate["outputs_with_combined_cell_majorant_within_requested_scale"]) == closures, "aggregate closure drift")
    expected_status = (
        "complete_cubic_legendre_tail_transported_over_all_selector_cells_and_combined_budget_closes"
        if closures == 15
        else "complete_cubic_legendre_tail_cell_transport_breaks_combined_budget"
    )
    require(result["status"] == expected_status, "status drift")
    require("not a proof" in gate.NOTE.read_text(encoding="utf-8").lower(), "note proof boundary missing")
    print(
        "validated coefficient-cell complete cubic Legendre-tail budget: "
        f"374 cells, 753 child indices, {closures}/15 combined outputs below 0.005"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
