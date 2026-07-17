#!/usr/bin/env python3
"""Build the left thousandth-tile guard for the order-ten pilot."""

from fractions import Fraction
from pathlib import Path

import jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache as base


base.START_T = Fraction(2485, 2)
base.END_T = Fraction(1243)
base.TILE_WIDTH_T = Fraction(1, 1000)
base.DEFAULT_CACHE = (
    base.REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_jet_pilot_h2_h22_thousandth_left_guard.jsonl"
)
base.DEFAULT_MANIFEST = (
    base.REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_jet_pilot_h2_h22_thousandth_left_guard_cache.json"
)
base.CACHE_TILE_LABEL = "order10_pilot_thousandth_left_guard"
base.GENERATOR_PATH = Path(__file__).resolve().relative_to(base.REPO_ROOT).as_posix()


if __name__ == "__main__":
    raise SystemExit(base.main())
