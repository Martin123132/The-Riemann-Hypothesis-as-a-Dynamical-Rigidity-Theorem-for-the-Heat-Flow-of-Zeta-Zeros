#!/usr/bin/env python3
"""Independently validate the direct R=9 outgoing-ray enclosure."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_ROOT = Path(__file__).resolve().parent
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
for candidate in (VENDOR, SCRIPT_ROOT):
    sys.path.insert(0, str(candidate))

from flint import arb, acb, ctx

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_R9_outgoing_ray_quadrature_gate as gate
from jensen_window_pf_newman_c1_hardy_t1e10_equation9_characteristic_canonical_core_quadrature_gate import (
    CharacteristicCore,
    parse_complex,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    started = time.perf_counter()
    priority = gate.set_low_priority()
    require(gate.RESULT.is_file() and gate.CACHE.is_file(), "missing ray result or cache")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "ray artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["geometry_gate"]["sha256"] == gate.file_hash(gate.GEOMETRY_GATE), "geometry hash mismatch")
    require(artifact["dependencies"]["connector_gate"]["sha256"] == gate.file_hash(gate.CONNECTOR_GATE), "connector hash mismatch")
    require(artifact["dependencies"]["cache"]["sha256"] == gate.file_hash(gate.CACHE), "cache hash mismatch")

    rows = gate.load_cache()
    expected = {(side, index) for side in ("positive", "negative") for index in range(gate.S_PANELS)}
    require(set(rows) == expected, "ray cache key set mismatch")
    for side, index in sorted(expected):
        row = rows[(side, index)]
        left, right = gate.panel_bounds(index)
        require(row["formula_version"] == gate.FORMULA_VERSION, "formula version mismatch")
        require(arb(row["s_left"]) == left and arb(row["s_right"]) == right, "panel boundary mismatch")
        require(row["outer_tolerance"] == gate.OUTER_TOLERANCE, "outer tolerance mismatch")
        require(row["inner_tolerance"] == gate.INNER_TOLERANCE, "inner tolerance mismatch")

    core = CharacteristicCore()
    positive = gate.sum_rows(rows, "positive")
    negative = gate.sum_rows(rows, "negative")
    truncated = positive + negative

    s = arb(gate.S_CUTOFF)
    linear = (arb(gate.RADIUS) ** 2 - core.lam - arb(gate.Y_BOUND)) / 2
    quadratic = arb(3).sqrt() * arb(gate.RADIUS) / 2
    q = linear * s + quadratic * s**2 + s**3 / 3
    q_prime = linear + 2 * quadratic * s + s**2
    tail_bound = 2 * arb(gate.MODE_COUNT) * arb(gate.Y_BOUND) * (-q).exp() / q_prime
    direct = acb(arb(truncated.real, tail_bound), arb(truncated.imag, tail_bound))

    stored = parse_complex(artifact["certified_values"]["combined_ray_complete_enclosure_ball"])
    require(stored.real.overlaps(direct.real) and stored.imag.overlaps(direct.imag), "stored direct enclosure mismatch")
    require(linear > arb("5.9") and tail_bound < arb("1e-14"), "analytic tail decision failed")

    connector = json.loads(gate.CONNECTOR_GATE.read_text(encoding="utf-8"))
    independent = parse_complex(connector["certified_values"]["true_R9_outgoing_tail_ball"])
    require(direct.real.overlaps(independent.real), "direct real component misses independent tail")
    require(direct.imag.overlaps(independent.imag), "direct imaginary component misses independent tail")

    for side, index in (("positive", 0), ("negative", gate.S_PANELS - 1)):
        replay = gate.integrate_ray_panel(core, side, index)
        cached = parse_complex(rows[(side, index)]["value"])
        require(replay.real.overlaps(cached.real), f"{side} spot-check real mismatch")
        require(replay.imag.overlaps(cached.imag), f"{side} spot-check imaginary mismatch")

    print(
        "validated direct R9 rays independently: "
        f"{len(rows)} cache rows, analytic tail <1e-14, both components overlap; "
        f"priority={priority}, elapsed={time.perf_counter()-started:.3f}s"
    )


if __name__ == "__main__":
    main()
