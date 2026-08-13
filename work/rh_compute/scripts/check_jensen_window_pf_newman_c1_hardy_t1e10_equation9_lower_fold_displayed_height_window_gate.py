#!/usr/bin/env python3
"""Independently validate the displayed lower-fold height window."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_displayed_height_window_gate as gate


CHECKER_PANELS = 4099


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.cell_gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing displayed-window artifact")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "displayed-window artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["source_output"]["sha256"] == gate.file_hash(gate.SOURCE_OUTPUT), "source output hash mismatch")
    require(artifact["dependencies"]["portcullis_gate"]["sha256"] == gate.file_hash(gate.PORTCULLIS_GATE), "portcullis gate hash mismatch")
    require(artifact["dependencies"]["first_nonzero_cell_gate"]["sha256"] == gate.file_hash(gate.CELL_GATE), "cell gate hash mismatch")
    require(
        artifact["dependencies"]["height_transport_cache"]["sha256"] == gate.file_hash(gate.cell_gate.HEIGHT_CACHE),
        "height cache hash mismatch",
    )

    saved = artifact["certified_displayed_window"]
    fresh = gate.certificate(CHECKER_PANELS)
    require(saved["source_height_roster"] == fresh["source_height_roster"], "source height roster mismatch")
    require(
        saved["turning_roster_stability"] == fresh["turning_roster_stability"],
        "turning-roster stability mismatch",
    )
    saved_cell = saved["canonical_lower_fold_window"]
    fresh_cell = fresh["canonical_lower_fold_window"]
    for key in ("grouped_center_value", "grouped_center_first_height_derivative", "grouped_center_second_height_derivative"):
        require(
            gate.cell_gate.parse_complex(saved_cell[key]).overlaps(gate.cell_gate.parse_complex(fresh_cell[key])),
            f"{key} mismatch",
        )
    require(arb(fresh_cell["derived_uniform_second_derivative_bound_ball"]) < arb("0.146"), "second bound failed")
    require(arb(fresh_cell["derived_taylor_remainder_bound_ball"]) < arb("0.000358"), "remainder bound failed")
    require(arb(fresh_cell["derived_total_variation_bound_ball"]) < arb("0.009494"), "variation bound failed")

    decisions = artifact["decision"]
    require(decisions["all_15_displayed_heights_inside_one_selector_stable_cell"] is True, "window decision missing")
    require(decisions["all_15_displayed_heights_preserve_84_mode_turning_roster"] is True, "turning-roster decision missing")
    require(decisions["printed_hardy_values_used_as_proof_data"] is False, "source-value boundary corrupted")
    require(decisions["complete_source_formula_enclosed"] is False, "source-formula boundary corrupted")
    require(decisions["selector_transition_join_proved"] is False, "selector boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(
        "validated displayed lower-fold height window independently: "
        f"15 offsets and {CHECKER_PANELS} Airy panels pass; priority={priority}"
    )


if __name__ == "__main__":
    main()
