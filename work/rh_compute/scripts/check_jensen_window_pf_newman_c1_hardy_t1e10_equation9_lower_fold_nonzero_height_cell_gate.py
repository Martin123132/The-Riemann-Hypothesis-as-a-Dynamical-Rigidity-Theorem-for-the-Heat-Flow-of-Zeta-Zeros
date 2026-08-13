#!/usr/bin/env python3
"""Independently validate the first nonzero lower-fold height cell."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_nonzero_height_cell_gate as gate


CHECKER_PANELS = 3073


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing nonzero height-cell artifact")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "height-cell artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(
        artifact["dependencies"]["height_transport_gate"]["sha256"] == gate.file_hash(gate.HEIGHT_GATE),
        "height-transport gate hash mismatch",
    )
    require(
        artifact["dependencies"]["height_transport_cache"]["sha256"] == gate.file_hash(gate.HEIGHT_CACHE),
        "height-transport cache hash mismatch",
    )

    p = gate.ode_gate.parameters()
    rows = gate.height_gate.load_cache()
    fresh = gate.certificate(rows, p, CHECKER_PANELS)
    saved = artifact["certified_height_cell"]
    for key in ("grouped_center_value", "grouped_center_first_height_derivative", "grouped_center_second_height_derivative"):
        require(gate.parse_complex(saved[key]).overlaps(gate.parse_complex(fresh[key])), f"{key} mismatch")

    require(saved["stated_uniform_second_derivative_bound"] == "0.146", "second bound changed")
    require(saved["stated_taylor_remainder_bound"] == "0.0000073", "remainder bound changed")
    require(saved["stated_total_variation_bound"] == "0.001313", "variation bound changed")
    require(fresh["airy_envelope"]["all_panels_strictly_below_envelope"] is True, "fresh Airy envelope failed")
    require(
        arb(fresh["derived_uniform_second_derivative_bound_ball"]) < arb("0.146"),
        "fresh second derivative does not meet the stated bound",
    )
    require(
        arb(fresh["derived_taylor_remainder_bound_ball"]) <= arb("0.0000073"),
        "fresh Taylor remainder does not meet the stated bound",
    )
    require(
        arb(fresh["derived_total_variation_bound_ball"]) < arb("0.001313"),
        "fresh total variation does not meet the stated bound",
    )

    decisions = artifact["decision"]
    require(decisions["nonzero_radius_height_cell_proved"] is True, "height-cell decision missing")
    require(decisions["selector_transition_join_proved"] is False, "selector boundary corrupted")
    require(decisions["source_height_interval_proved"] is False, "source-height boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(
        "validated first nonzero lower-fold height cell independently: "
        f"{CHECKER_PANELS} Airy panels pass; priority={priority}"
    )


if __name__ == "__main__":
    main()
