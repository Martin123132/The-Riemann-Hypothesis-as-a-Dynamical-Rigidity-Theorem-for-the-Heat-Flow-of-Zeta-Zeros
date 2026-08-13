#!/usr/bin/env python3
"""Independently validate the lower-fold turning-event atlas."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_turning_event_atlas_handoff_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.window_gate.cell_gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing turning-event artifact")
    ctx.dps = 120
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "turning-event artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["displayed_window_gate"]["sha256"] == gate.file_hash(gate.WINDOW_GATE), "window hash mismatch")
    require(artifact["dependencies"]["portcullis_gate"]["sha256"] == gate.file_hash(gate.PORTCULLIS_GATE), "portcullis hash mismatch")

    saved = artifact["certified_event_atlas"]
    fresh = gate.certificate()
    for key in (
        "event_count_inside_selector_cell",
        "lower_edge_transition_roster",
        "saved_transition_count",
        "saved_transition_roster",
        "next_upward_event_mode",
        "next_downward_event_mode",
        "adjacent_selector_lower_edge_transition_roster",
    ):
        require(saved[key] == fresh[key], f"{key} mismatch")
    for key in (
        "minimum_event_spacing_ball",
        "next_upward_event_distance_ball",
        "next_downward_event_distance_ball",
        "displayed_window_upper_event_margin_ball",
        "displayed_window_lower_event_margin_ball",
    ):
        require(arb(saved[key]).overlaps(arb(fresh[key])), f"{key} interval mismatch")

    sequence = [gate.mode_for_event(index) for index in range(gate.EVENT_COUNT)]
    require(sequence[:6] == [39894, 39895, 39893, 39896, 39892, 39897], "event alternation mismatch")
    require(sequence[83] == 39936 and sequence[84] == 39852, "nearest event modes mismatch")
    require(set(sequence[:84]) == set(range(39853, 39937)), "84-mode roster mismatch")
    for index in range(gate.EVENT_COUNT - 1):
        spacing = gate.event_height(index) - gate.event_height(index + 1)
        require((spacing - arb.pi() * (index + 1)).contains(0), f"event spacing {index} failed")

    decisions = artifact["decision"]
    require(decisions["all_399_one_mode_handoffs_preserve_exact_finite_sum"] is True, "handoff decision missing")
    require(decisions["selector_jump_reducible_to_one_mode_handoff"] is False, "selector-jump guard corrupted")
    require(decisions["airy_morse_approximation_join_proved"] is False, "approximation boundary corrupted")
    require(decisions["selector_strip_reassembly_proved"] is False, "strip boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated turning-event atlas independently: 399 handoffs and exact spacings pass; priority={priority}")


if __name__ == "__main__":
    main()
