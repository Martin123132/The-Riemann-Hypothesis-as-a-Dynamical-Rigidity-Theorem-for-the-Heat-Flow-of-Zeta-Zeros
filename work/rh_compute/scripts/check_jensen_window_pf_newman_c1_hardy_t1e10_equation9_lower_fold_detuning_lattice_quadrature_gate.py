#!/usr/bin/env python3
"""Independently validate the complete lower-fold detuning lattice."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_detuning_lattice_quadrature_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file() and gate.CACHE.is_file(), "missing lattice result or cache")
    ctx.dps = gate.PRECISION
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "lattice artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["sources"]["cache"]["sha256"] == gate.file_hash(gate.CACHE), "cache hash mismatch")
    require(artifact["dependencies"]["ode_gate"]["sha256"] == gate.file_hash(gate.ODE_GATE), "ODE hash mismatch")
    require(artifact["dependencies"]["core_gate"]["sha256"] == gate.file_hash(gate.CORE_GATE), "core hash mismatch")

    rows = gate.load_cache()
    require(set(rows) == set(range(gate.MODE_LO, gate.MODE_HI + 1)), "cached mode roster mismatch")
    p = gate.ode_gate.parameters()
    fresh = gate.build_certificate(rows, p)
    saved = artifact["certified_lattice"]
    require(saved["row_count"] == 84, "saved row count mismatch")
    require(gate.parse_complex(saved["grouped_2pi_sum"]).overlaps(gate.parse_complex(fresh["grouped_2pi_sum"])), "grouped sum mismatch")

    for mode in (gate.MODE_LO, 39894, gate.MODE_HI):
        recomputed = gate.compute_row(mode, p)
        require(
            gate.parse_complex(rows[mode]["G64"]).overlaps(gate.parse_complex(recomputed["G64"])),
            f"mode {mode} recomputation mismatch",
        )

    decisions = artifact["decision"]
    require(decisions["all_84_detuning_values_rigorously_integrated"] is True, "complete lattice decision missing")
    require(decisions["grouped_sum_overlaps_independent_Airy_first_value"] is True, "Airy overlap decision missing")
    require(decisions["height_uniform_propagation_proved"] is False, "height boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated complete detuning lattice independently: 84 rows and 3 recomputations pass; priority={priority}")


if __name__ == "__main__":
    main()
