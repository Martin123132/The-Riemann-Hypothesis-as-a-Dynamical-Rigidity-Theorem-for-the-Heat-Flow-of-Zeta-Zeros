#!/usr/bin/env python3
"""Independently validate fixed-selector lower-fold height transport."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_height_transport_identity_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file() and gate.CACHE.is_file(), "missing height-transport result or cache")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "height-transport artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["sources"]["cache"]["sha256"] == gate.file_hash(gate.CACHE), "cache hash mismatch")
    require(artifact["dependencies"]["ode_gate"]["sha256"] == gate.file_hash(gate.ODE_GATE), "ODE hash mismatch")
    require(artifact["dependencies"]["lattice_gate"]["sha256"] == gate.file_hash(gate.LATTICE_GATE), "lattice hash mismatch")
    require(artifact["dependencies"]["lattice_cache"]["sha256"] == gate.file_hash(gate.lattice_gate.CACHE), "lattice cache hash mismatch")

    rows = gate.load_cache()
    require(set(rows) == set(range(gate.MODE_LO, gate.MODE_HI + 1)), "height derivative roster mismatch")
    p = gate.ode_gate.parameters()
    fresh = gate.certificate(rows, p)
    saved = artifact["certified_transport"]
    require(saved["row_count"] == 84, "saved derivative row count mismatch")
    require(
        gate.parse_complex(saved["grouped_2pi_height_derivative_fixed_C"]).overlaps(
            gate.parse_complex(fresh["grouped_2pi_height_derivative_fixed_C"])
        ),
        "grouped height derivative mismatch",
    )

    lattice_rows = gate.lattice_gate.load_cache()
    for mode in (gate.MODE_LO, 39894, gate.MODE_HI):
        recomputed = gate.compute_row(mode, p, lattice_rows)
        require(
            gate.parse_complex(rows[mode]["G64_height_derivative_fixed_C"]).overlaps(
                gate.parse_complex(recomputed["G64_height_derivative_fixed_C"])
            ),
            f"mode {mode} height derivative mismatch",
        )

    decisions = artifact["decision"]
    require(decisions["exact_finite_transform_height_derivative_proved"] is True, "transport identity missing")
    require(decisions["all_84_height_derivatives_rigorously_evaluated"] is True, "complete derivative decision missing")
    require(decisions["nonzero_radius_height_cell_proved"] is False, "height-cell boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated lower-fold height transport independently: 84 rows and 3 recomputations pass; priority={priority}")


if __name__ == "__main__":
    main()
