#!/usr/bin/env python3
"""Validate the block-20 coefficient-cell joint recurrence identity."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


SCRIPT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_ROOT))

import jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate as gate
import jensen_window_pf_newman_c1_hardy_block20_native_q_required_compensation_gate as native_q


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    require(gate.RESULT.is_file(), f"missing result: {gate.RESULT}")
    require(gate.NOTE.is_file(), f"missing note: {gate.NOTE}")
    result = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(result["kind"] == "jensen_window_pf_newman_c1_hardy_block20_coefficient_cell_joint_recurrence_identity_gate", "kind drift")
    require(result["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash drift")
    require(result["sources"]["checker"]["sha256"] == gate.file_hash(Path(__file__).resolve()), "checker hash drift")
    for dependency in result["dependencies"].values():
        path = gate.REPO_ROOT / dependency["path"]
        require(path.is_file(), f"missing dependency: {path}")
        require(dependency["sha256"] == gate.file_hash(path), f"dependency hash drift: {path}")
    require(gate.algebra_audit() == 12, "formal algebra audit drift")

    atlas = json.loads(gate.ATLAS.read_text(encoding="utf-8"))
    atlas_rows = {int(row["chain"]): row for row in atlas["rows"] if int(row["mit"]) == 2}
    recursive = native_q.load_recursive_chains()
    rows = result["rows"]
    require(len(rows) == 374 and len({int(row["chain"]) for row in rows}) == 374, "cell roster drift")
    for row in rows:
        chain = int(row["chain"])
        margins = gate.cell_margins(recursive[chain], atlas_rows[chain])
        comparisons = (
            ("a1_lower_fundamental_margin", "a1_lower_margin"),
            ("a1_upper_fundamental_margin", "a1_upper_margin"),
            ("dual_lower_selector_margin", "dual_lower_margin"),
            ("dual_upper_selector_margin", "dual_upper_margin"),
            ("upper_nonsaddle_gap_lower", "upper_nonsaddle_gap"),
            ("lower_nonsaddle_gap_lower", "lower_nonsaddle_gap"),
            ("curvature_floor_lower", "curvature_floor"),
        )
        for stored, computed in comparisons:
            require(Fraction(row[stored]) == Fraction(cells_value(margins[computed])), f"chain {chain} {stored} drift")
        require(int(row["jbot"]) == margins["jbot"], "jbot drift")
        require(int(row["child_length"]) == margins["child_length"], "child length drift")
        require(int(row["cubic_sign"]) == margins["cubic_sign"], "cubic sign drift")
        require(bool(row["point_equation_69_precision_overlap"]), "equation-(69) point guard drift")
        require(bool(row["point_poisson_identity_contains_zero"]), "Poisson point guard drift")

    endpoint = json.loads(gate.ENDPOINT_CELLS.read_text(encoding="utf-8"))
    endpoint_outputs = {int(row["output_index"]): row for row in endpoint["transported_outputs"]}
    outputs = result["transported_outputs"]
    require(len(outputs) == 15, "output count drift")
    for row in outputs:
        index = int(row["output_index"])
        bound = Fraction(endpoint_outputs[index]["cell_endpoint_majorant_upper"])
        require(Fraction(row["joint_corrected_recurrence_majorant_upper"]) == bound, "endpoint reduction drift")
        require(bool(row["within_requested_scale"]) == (bound < gate.REQUESTED_ERROR_SCALE), "output closure drift")

    aggregate = result["aggregate"]
    require(int(aggregate["exact_algebra_checks"]) == 12, "aggregate algebra drift")
    require(int(aggregate["cell_count"]) == 374, "aggregate cell drift")
    require(int(aggregate["point_equation_69_overlap_count"]) == 374, "aggregate equation-(69) drift")
    require(int(aggregate["point_poisson_identity_count"]) == 374, "aggregate Poisson drift")
    require(int(aggregate["outputs_within_requested_scale"]) == 15, "aggregate closure drift")
    require(
        result["status"]
        == "joint_corrected_recurrence_reduces_exactly_to_cell_endpoint_error_on_all_374_physical_cell_slices",
        "status drift",
    )
    require("F'(N)=L+fracL" in result["scope"], "physical-slice scope guard missing")
    require("not a proof" in gate.NOTE.read_text(encoding="utf-8").lower(), "note proof boundary missing")
    print("validated coefficient-cell joint recurrence identity: 374 cells, 15/15 outputs below 0.005")
    return 0


def cells_value(value: Fraction) -> str:
    return gate.cells.decimal_text(value)


if __name__ == "__main__":
    raise SystemExit(main())
