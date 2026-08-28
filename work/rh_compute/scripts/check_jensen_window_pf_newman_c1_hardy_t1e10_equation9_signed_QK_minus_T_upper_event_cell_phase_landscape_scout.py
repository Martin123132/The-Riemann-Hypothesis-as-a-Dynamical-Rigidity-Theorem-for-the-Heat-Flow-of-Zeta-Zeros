#!/usr/bin/env python3
"""Validate the cache-backed upper-event-cell phase landscape scout."""

from __future__ import annotations

import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_c1_hardy_t1e10_equation9_signed_QK_minus_T_upper_event_cell_phase_landscape_scout as scout


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    require(scout.RESULT.is_file(), f"missing result: {scout.RESULT}")
    require(scout.NOTE.is_file(), f"missing note: {scout.NOTE}")
    saved = json.loads(scout.RESULT.read_text(encoding="utf-8"))
    rebuilt = scout.build_artifact()
    require(saved == rebuilt, "saved landscape artifact differs from cache-backed rebuild")
    require(scout.NOTE.read_text(encoding="utf-8") == scout.render_note(saved), "landscape note drifted")

    values = saved["value_landscape"]
    derivatives = saved["derivative_landscape"]
    require(len(values["integer_grid"]) == 20, "integer grid count drifted")
    require(len(values["candidate_trough_grid"]) == 13, "trough grid count drifted")
    require(len(values["candidate_peak_grid"]) == 13, "peak grid count drifted")
    require(values["candidate_peak_upper_envelope"].startswith("-"), "peak envelope lost negativity")
    require(
        derivatives["sign_sequence"] == ["positive", "negative", "negative", "positive"],
        "derivative sign sequence drifted",
    )
    require(
        saved["route_decision"]["candidate_full_cell_union"] == ["0", "19.5"],
        "candidate phase-cell union drifted",
    )
    require("not certified extrema" in saved["proof_boundary"], "proof boundary weakened")
    require("imply RH" in saved["proof_boundary"], "RH boundary weakened")

    print(
        "validated upper-event-cell phase landscape scout: "
        "20 integer boxes, 13 trough candidates, 13 peak candidates, "
        "derivative signs +--+, 0 interval-sign promotions, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
