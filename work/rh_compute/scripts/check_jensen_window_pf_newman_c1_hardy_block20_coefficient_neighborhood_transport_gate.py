#!/usr/bin/env python3
"""Validate the block-20 coefficient-neighborhood transport artifact."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate.json"
NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate.md"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def resolve(path: str) -> Path:
    return REPO_ROOT / Path(path)


def nested(row: dict[str, Any], *path: str) -> Any:
    value: Any = row
    for key in path:
        value = value[key]
    return value


def main() -> int:
    for path in (RESULT, NOTE):
        require(path.is_file(), f"missing coefficient-neighborhood artifact: {path}")
    artifact = json.loads(RESULT.read_text(encoding="utf-8"))
    require(
        artifact.get("kind") == "jensen_window_pf_newman_c1_hardy_block20_coefficient_neighborhood_transport_gate",
        "coefficient-neighborhood kind drift",
    )
    rows = artifact.get("rows", [])
    outputs = artifact.get("transported_outputs", [])
    aggregate = artifact.get("aggregate", {})
    require(len(rows) == aggregate.get("cell_count") == 374, "coefficient-neighborhood full roster drift")
    require(aggregate.get("full_roster") is True, "coefficient-neighborhood artifact is not full")
    require(aggregate.get("precision_overlap_count") == 374, "precision overlap drift")
    require(aggregate.get("w3_count") + aggregate.get("w4_count") == 374, "selector count drift")
    require(len(outputs) == 15, "transported output count drift")
    require([int(row["chain"]) for row in rows] == sorted(int(row["chain"]) for row in rows), "chain ordering drift")
    require(len({int(row["chain"]) for row in rows}) == 374, "duplicate chain")

    for row in rows:
        require(row["selector"] in ("W3", "W4"), f"chain {row['chain']} bad selector")
        require(row["child_column"]["source_kernel_contained"] is True, f"chain {row['chain']} child containment drift")
        require(row["multiplier_column"]["source_multiplier_contained"] is True, f"chain {row['chain']} multiplier containment drift")
        require(row["endpoint_column"]["precision_overlap"] is True, f"chain {row['chain']} endpoint overlap drift")
        require(
            Fraction(row["endpoint_column"]["center_majorant_upper"])
            <= Fraction(row["endpoint_column"]["cell_majorant_upper"]),
            f"chain {row['chain']} center/cell order drift",
        )
        require(
            all(Fraction(value) > 0 for value in row["minimum_positive_parameters"].values()),
            f"chain {row['chain']} positive parameter drift",
        )
        require(row["arithmetic_column"]["status"] == "open_not_allocated_to_endpoint_budget", "arithmetic boundary drift")

    for output in outputs:
        require(
            Fraction(output["center_endpoint_majorant_upper"])
            <= Fraction(output["cell_endpoint_majorant_upper"]),
            f"output {output['output_index']} center/cell order drift",
        )
        require(
            output["within_requested_scale"]
            == (Fraction(output["cell_endpoint_majorant_upper"]) < Fraction(5, 1000)),
            f"output {output['output_index']} scale classification drift",
        )

    maxima = {
        "maximum_center_endpoint_majorant_upper": max(Fraction(nested(row, "endpoint_column", "center_majorant_upper")) for row in rows),
        "maximum_cell_endpoint_majorant_upper": max(Fraction(nested(row, "endpoint_column", "cell_majorant_upper")) for row in rows),
        "maximum_cell_endpoint_inflation_upper": max(Fraction(nested(row, "endpoint_column", "cell_to_center_inflation_upper")) for row in rows),
        "maximum_parent_phase_lipschitz_upper": max(Fraction(nested(row, "phase_column", "parent_phase_lipschitz_upper")) for row in rows),
        "maximum_child_kernel_variation_upper": max(Fraction(nested(row, "child_column", "kernel_variation_upper")) for row in rows),
        "maximum_multiplier_relative_variation_upper": max(Fraction(nested(row, "multiplier_column", "relative_variation_upper")) for row in rows),
        "maximum_saddle_lipschitz_upper": max(Fraction(nested(row, "saddle_column", "equation_69_lipschitz_upper")) for row in rows),
    }
    for name, expected in maxima.items():
        require(Fraction(aggregate[name]) == expected, f"aggregate {name} drift")
    output_count = sum(bool(row["within_requested_scale"]) for row in outputs)
    require(aggregate.get("outputs_within_requested_scale") == output_count, "output closure count drift")
    expected_status = (
        "rigorous_full_cell_endpoint_transport_closes_finite_block20_scale_other_columns_open"
        if output_count == 15
        else "rigorous_full_cell_endpoint_transport_exceeds_finite_block20_scale_other_columns_open"
    )
    require(artifact.get("status") == expected_status, "coefficient-neighborhood status drift")

    for group in ("dependencies", "sources"):
        for name, record in artifact[group].items():
            path = resolve(record["path"])
            require(path.is_file(), f"missing {group} file {name}: {path}")
            require(file_hash(path) == record["sha256"], f"{group} hash drift: {name}")
    note = NOTE.read_text(encoding="utf-8")
    require("separate open budget" in note or "separate open" in note, "note proof boundary drift")
    require("No complete evaluator theorem" in artifact["proof_boundary"], "result proof boundary drift")
    print(
        "validated coefficient-neighborhood transport: "
        f"374 cells, {output_count}/15 endpoint outputs below 0.005, other columns open"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
