#!/usr/bin/env python3
"""Independently validate the y=64 primitive and three-chart envelopes."""

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
import sympy as sp

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_y64_three_chart_primitive_envelope_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing three-chart result")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "three-chart artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["partition_gate"]["sha256"] == gate.file_hash(gate.PARTITION_GATE), "partition hash mismatch")
    require(artifact["dependencies"]["fresnel_gate"]["sha256"] == gate.file_hash(gate.FRESNEL_GATE), "Fresnel hash mismatch")

    y, k, p = sp.symbols("y k p", positive=True, real=True)
    completed = k * (y - p / k) ** 2 / 2 - p**2 / (2 * k)
    require(sp.simplify(k * y**2 / 2 - p * y - completed) == 0, "independent square completion failed")
    require(
        artifact["symbolic_primitive"]["full_gaussian"]
        == "H_full=exp(-i*p^2/(2*k))*sqrt(2*pi/k)*exp(i*pi/4)",
        "full Gaussian carrier is missing",
    )

    fresh = gate.certified_envelopes()
    saved = artifact["certified_envelopes"]
    keys = (
        "uniform_g_buffer_delta_ball",
        "nonstationary_coefficient_margin_ball",
        "per_mode_nonstationary_outer_bound_ball",
        "transition_normalized_q_max_ball",
        "per_mode_transition_affine_bound_ball",
        "per_mode_grouped_Morse_main_bound_ball",
        "per_mode_Morse_Fresnel_remainder_bound_ball",
        "upper_affine_amplitude_ball",
    )
    for key in keys:
        require(arb(saved[key]).overlaps(arb(fresh[key])), f"saved {key} mismatch")

    require(arb(saved["nonstationary_coefficient_margin_ball"]) > arb("0.99"), "IBP coefficient margin failed")
    require(arb(saved["per_mode_nonstationary_outer_bound_ball"]) < arb(33), "nonstationary decision failed")
    require(arb(saved["per_mode_transition_affine_bound_ball"]) < arb(600), "transition decision failed")
    require(arb(saved["per_mode_Morse_Fresnel_remainder_bound_ball"]) < arb(34), "Morse remainder decision failed")
    require(arb(saved["per_mode_grouped_Morse_main_bound_ball"]) < arb(200), "Morse main decision failed")

    decisions = artifact["decision"]
    require(decisions["artificial_y64_endpoint_cancels_exactly"] is True, "telescope decision missing")
    require(decisions["grouped_84_mode_bound_proved"] is False, "grouped-boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(
        "validated y=64 primitives independently: source pairing exact, "
        f"nonstat<33, transition<600, Morse remainder<34; priority={priority}"
    )


if __name__ == "__main__":
    main()
