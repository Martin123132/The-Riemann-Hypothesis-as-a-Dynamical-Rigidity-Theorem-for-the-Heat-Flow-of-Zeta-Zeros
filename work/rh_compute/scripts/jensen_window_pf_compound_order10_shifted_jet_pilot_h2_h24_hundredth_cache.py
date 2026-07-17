#!/usr/bin/env python3
"""Build a hundredth-tile H2-H24 pilot collar for order ten."""

from fractions import Fraction
from pathlib import Path

import jensen_window_pf_compound_order7_shifted_jet_compact_h2_h22_cache as base


base.START_T = Fraction(1243)
base.END_T = Fraction(1260)
base.TILE_WIDTH_T = Fraction(1, 100)
base.MAX_MOMENT = 24
base.DEFAULT_CACHE = (
    base.REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_jet_pilot_h2_h24_hundredth_tiles.jsonl"
)
base.DEFAULT_MANIFEST = (
    base.REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order10_shifted_jet_pilot_h2_h24_hundredth_cache.json"
)
base.CACHE_TILE_LABEL = "order10_pilot_h2_h24_hundredth"
base.GENERATOR_PATH = Path(__file__).resolve().relative_to(base.REPO_ROOT).as_posix()


if __name__ == "__main__":
    raise SystemExit(base.main())
