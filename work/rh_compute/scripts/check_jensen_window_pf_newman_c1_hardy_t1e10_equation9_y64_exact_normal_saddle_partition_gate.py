#!/usr/bin/env python3
"""Independently validate the exact y=64 normal-saddle partition."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_exact_normal_saddle_partition_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing partition result")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "partition artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["normal_form_gate"]["sha256"] == gate.file_hash(gate.NORMAL_FORM_GATE), "normal-form hash mismatch")
    require(artifact["dependencies"]["ray_gate"]["sha256"] == gate.file_hash(gate.RAY_GATE), "ray hash mismatch")

    fresh = gate.certified_partition()
    saved = artifact["certified_partition"]
    scalar_keys = (
        "uniform_g_buffer_delta_ball",
        "minimum_normalized_boundary_cut_ball",
        "left_cut_global_min_ball",
        "crossing_global_max_ball",
        "right_cut_global_max_ball",
        "maximum_exact_minus_canonical_crossing_ball",
        "maximum_normal_saddle_on_R9_ball",
        "minimum_upper_endpoint_normalized_clearance_ball",
    )
    for key in scalar_keys:
        require(arb(saved[key]).overlaps(arb(fresh[key])), f"saved {key} mismatch")

    require(len(saved["roster"]) == 84, "saved roster length mismatch")
    for index in (0, 41, 42, 83):
        old = saved["roster"][index]
        new = fresh["roster"][index]
        require(old["mode"] == new["mode"], "mode identity mismatch")
        for key in (
            "nonstationary_to_Fresnel_z_ball",
            "saddle_crossing_z_ball",
            "Fresnel_to_Morse_z_ball",
        ):
            require(arb(old[key]).overlaps(arb(new[key])), f"mode {old['mode']} {key} mismatch")
        require(
            arb(old["nonstationary_to_Fresnel_z_ball"])
            < arb(old["saddle_crossing_z_ball"])
            < arb(old["Fresnel_to_Morse_z_ball"]),
            f"mode {old['mode']} cut ordering failed",
        )

    decisions = artifact["decision"]
    require(decisions["three_chart_partition_disjoint_and_exhaustive"] is True, "partition decision missing")
    require(decisions["three_chart_integral_estimates_proved"] is False, "estimate boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(
        "validated y=64 partition independently: 84 exact closed-form crossings, "
        f"four-width cuts inside R9; priority={priority}"
    )


if __name__ == "__main__":
    main()
