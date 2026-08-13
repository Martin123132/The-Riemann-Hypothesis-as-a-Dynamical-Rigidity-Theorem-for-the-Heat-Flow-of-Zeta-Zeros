#!/usr/bin/env python3
"""Independently validate the moving y=64 mode roster."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_moving_roster_dirichlet_partition_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing moving-roster result")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "moving-roster artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["partition_gate"]["sha256"] == gate.file_hash(gate.PARTITION_GATE), "partition hash mismatch")
    require(artifact["dependencies"]["primitive_gate"]["sha256"] == gate.file_hash(gate.PRIMITIVE_GATE), "primitive hash mismatch")

    fresh = gate.certified_moving_roster()
    saved = artifact["certified_moving_roster"]
    for key in (
        "h_ball",
        "delta_ball",
        "real_Fresnel_roster_width_ball",
        "distance_of_width_to_nearest_integer_ball",
        "minimum_event_gap_ball",
    ):
        require(arb(saved[key]).overlaps(arb(fresh[key])), f"saved {key} mismatch")
    require(saved["event_count"] == fresh["event_count"] == 168, "event count mismatch")
    require(saved["open_cell_count"] == fresh["open_cell_count"] == 169, "cell count mismatch")
    require(saved["Fresnel_count_cell_histogram"] == fresh["Fresnel_count_cell_histogram"], "histogram mismatch")
    require(saved["maximum_simultaneous_Fresnel_modes"] == 3, "Fresnel mode cap mismatch")

    for index in (0, 83, 84, 167):
        require(saved["events"][index] == fresh["events"][index], f"event {index} mismatch")
    for index in (0, 1, 84, 167, 168):
        require(saved["cells"][index] == fresh["cells"][index], f"cell {index} mismatch")

    decisions = artifact["decision"]
    require(decisions["chart_memberships_form_contiguous_mode_blocks"] is True, "contiguity decision missing")
    require(decisions["grouped_chart_integral_bound_proved"] is False, "grouped-boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(
        "validated moving roster independently: 168 events, 169 cells, "
        f"Fresnel cap=3; priority={priority}"
    )


if __name__ == "__main__":
    main()
