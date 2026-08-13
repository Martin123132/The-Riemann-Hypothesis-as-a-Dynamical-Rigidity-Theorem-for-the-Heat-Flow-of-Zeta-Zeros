#!/usr/bin/env python3
"""Independently validate the grouped B-endpoint nonstationary enclosure."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_yB_upper_endpoint_nonstationary_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing B-endpoint result")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "B-endpoint artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["reassembly_gate"]["sha256"] == gate.file_hash(gate.REASSEMBLY_GATE), "reassembly hash mismatch")
    require(artifact["dependencies"]["normal_form_gate"]["sha256"] == gate.file_hash(gate.NORMAL_FORM_GATE), "normal-form hash mismatch")

    fresh = gate.certified_bound()
    saved = artifact["certified_bound"]
    for key in (
        "grouped_amplitude_max_ball",
        "grouped_amplitude_derivative_max_ball",
        "phase_derivative_absolute_min_ball",
        "phase_second_derivative_max_ball",
        "grouped_B_endpoint_integral_bound_ball",
    ):
        require(arb(saved[key]).overlaps(arb(fresh[key])), f"saved {key} mismatch")
    require(arb(saved["phase_derivative_absolute_min_ball"]) > arb("4.77e9"), "phase derivative decision failed")
    require(arb(saved["grouped_B_endpoint_integral_bound_ball"]) < arb("7.5e-12"), "IBP decision failed")

    decisions = artifact["decision"]
    require(decisions["exact_transformed_to_source_B_phase_reduction_proved"] is True, "phase reduction missing")
    require(decisions["B_endpoint_nonstationary_on_R9_proved"] is True, "nonstationary decision missing")
    require(decisions["lower_characteristic_fold_proved"] is False, "lower-fold boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated grouped B endpoint independently: |P'|>4.77e9, integral<7.5e-12; priority={priority}")


if __name__ == "__main__":
    main()
