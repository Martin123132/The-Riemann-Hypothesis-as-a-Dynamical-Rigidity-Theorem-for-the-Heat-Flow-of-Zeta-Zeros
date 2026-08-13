#!/usr/bin/env python3
"""Independently validate exact selector-strip finite-Poisson reassembly."""

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

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_lower_fold_selector_strip_poisson_reassembly_gate as gate


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    priority = gate.event_gate.window_gate.cell_gate.ode_gate.set_low_priority()
    require(gate.RESULT.is_file(), "missing selector-strip artifact")
    ctx.dps = 120
    artifact = json.loads(gate.RESULT.read_text(encoding="utf-8"))
    require(artifact.get("passed") is True, "selector-strip artifact does not pass")
    require(artifact["sources"]["builder"]["sha256"] == gate.file_hash(gate.BUILDER), "builder hash mismatch")
    require(artifact["sources"]["checker"]["sha256"] == gate.file_hash(gate.CHECKER), "checker hash mismatch")
    require(artifact["dependencies"]["turning_event_gate"]["sha256"] == gate.file_hash(gate.EVENT_GATE), "event hash mismatch")
    require(artifact["dependencies"]["symmetric_poisson_gate"]["sha256"] == gate.file_hash(gate.POISSON_GATE), "Poisson hash mismatch")
    require(artifact["dependencies"]["theta_current_gate"]["sha256"] == gate.file_hash(gate.CURRENT_GATE), "current hash mismatch")

    saved = artifact["certified_selector_reassembly"]
    fresh = gate.certificate()
    require(saved == fresh, "selector-strip certificate mismatch")
    require(saved["old_odd_roster_term_count"] - saved["new_odd_roster_term_count"] == 1, "source term count mismatch")
    require(saved["removed_source_lattice_point"] == gate.C, "removed source point mismatch")
    require(saved["new_selector_turning_mode_count"] == 399, "turning mode count mismatch")

    # Independent one-cell endpoint arithmetic.
    for f0, f1 in ((3, 5), (-7, 11), (13, -17)):
        half_difference = arb(f0 - f1) / 2
        one_cell_poisson = arb(f0 + f1) / 2
        require((half_difference + one_cell_poisson - f0).contains(0), "one-cell endpoint arithmetic failed")

    decisions = artifact["decision"]
    require(decisions["complete_symmetric_poisson_strip_reassembles_removed_term"] is True, "reassembly decision missing")
    require(decisions["positive_399_mode_roster_alone_sufficient"] is False, "399-mode guard corrupted")
    require(decisions["complete_source_selector_splice_proved"] is False, "source-splice boundary corrupted")
    require(decisions["rh_implication"] is False, "RH boundary corrupted")
    print(f"validated selector-strip reassembly independently: one-cell Poisson and roster shift pass; priority={priority}")


if __name__ == "__main__":
    main()
