#!/usr/bin/env python3
"""Independently check the fixed-B 399-mode canonical Fourier completion."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, acb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_fixed_B_selector_boundary_399_mode_fourier_completion_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file() and gate.CACHE.is_file(), "missing 399-mode artifact or cache")
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "399-mode artifact does not pass")
    for key, path in (
        ("source_transition_scheduling_gate", gate.SCHEDULING_GATE),
        ("selector_strip_gate", gate.SELECTOR_GATE),
        ("turning_event_gate", gate.EVENT_GATE),
    ):
        require(artifact["dependencies"][key]["sha256"] == gate.file_hash(path), f"{key} hash mismatch")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["sources"]["cache"]["sha256"] == gate.file_hash(gate.CACHE), "cache hash mismatch")

    ctx.dps = 80
    p = gate.parameters()
    require(artifact["symbolic_identities"] == gate.symbolic_identities(), "symbolic identity mismatch")
    fresh_geometry = gate.geometry_certificate(p)
    require(artifact["geometry"] == fresh_geometry, "geometry certificate mismatch")

    rows = gate.load_cache()
    require(set(rows) == set(range(gate.MODE_LO, gate.MODE_HI + 1)), "cache roster mismatch")
    coefficients = [gate.parse_complex(rows[mode]["coefficient"]) for mode in range(gate.MODE_LO, gate.MODE_HI + 1)]
    cache_integral = 2 * p["pi"] * sum(coefficients, acb(0))
    saved_integral = gate.parse_complex(artifact["numerical_completion"]["canonical_399_integral_2pi_sum"])
    require(cache_integral.real.overlaps(saved_integral.real) and cache_integral.imag.overlaps(saved_integral.imag), "cache sum misses saved integral")

    for mode in (gate.MODE_LO, 39894, gate.CENTER_MODE, gate.MODE_HI):
        fresh = gate.canonical_coefficient(mode, p, tolerance="1e-22")
        cached = gate.parse_complex(rows[mode]["coefficient"])
        require(fresh.real.overlaps(cached.real) and fresh.imag.overlaps(cached.imag), f"mode {mode} recomputation failed")

    direct = gate.direct_grouped_integral(p, panels=401, tolerance="1e-21")
    require(direct.real.overlaps(cache_integral.real) and direct.imag.overlaps(cache_integral.imag), "401-panel direct kernel misses cache sum")

    ai0 = arb(0).airy_ai()
    ai_y = (-p["Y"]).airy_ai()
    g0 = acb(p["Y"] * ai0)
    g1 = acb(-p["Y"] * ai_y)
    midpoint = (g0 + g1) / 2
    half_current = (g0 - g1) / 2
    coefficient_sum = cache_integral / (2 * p["pi"])
    complement = midpoint - coefficient_sum
    closure = half_current + coefficient_sum + complement - g0
    require(closure.real.contains(0) and closure.imag.contains(0), "independent Fourier closure failed")
    saved_complement = gate.parse_complex(artifact["numerical_completion"]["zero_negative_outer_positive_complement"])
    require(complement.real.overlaps(saved_complement.real) and complement.imag.overlaps(saved_complement.imag), "complement mismatch")

    decisions = artifact["decision"]
    require(decisions["fixed_theorem_endpoint_B"] == gate.FIXED_B, "fixed endpoint decision corrupted")
    require(decisions["all_399_joint_saddles_inside_one_cell_strip"] is True, "stationary roster decision corrupted")
    require(decisions["zero_negative_outer_positive_complement_completed_jointly"] is True, "completion decision corrupted")
    require(decisions["finite_t_canonical_error_bounded"] is False, "finite-t proof boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated fixed-B 399-mode Fourier completion independently; priority={priority}")


if __name__ == "__main__":
    main()
