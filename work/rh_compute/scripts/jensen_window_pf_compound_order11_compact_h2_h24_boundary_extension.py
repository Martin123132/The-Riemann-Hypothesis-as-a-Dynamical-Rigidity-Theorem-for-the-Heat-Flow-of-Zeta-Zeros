#!/usr/bin/env python3
"""Build the two H2-H24 boundary tiles needed by the order-eleven compact bridge."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
for candidate in (SCRIPT_DIR, REPO_ROOT / "work/rh_compute/vendor"):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import jensen_window_pf_compound_order10_compact_h2_h24_unit_cache as source  # noqa: E402


DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension.json"
TASKS = ((0, Fraction(5691), Fraction(5692)), (1, Fraction(38028), Fraction(38029)))


def build_artifact() -> dict:
    rows = []
    for task in TASKS:
        row = source.tile_task(task)
        if row.get("passed") is not True:
            raise RuntimeError(f"boundary H tile failed: {row}")
        rows.append({**row, "kind": "order11_compact_h2_h24_boundary_tile"})
    return {
        "kind": "jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension",
        "date": "2026-07-18",
        "status": "rigorous two-tile H2-H24 compact boundary extension",
        "proof_boundary": (
            "This artifact extends only the inherited compact H source by the two "
            "unit tiles required by the eighth localization layer."
        ),
        "source_contract": {
            "inherited_row_contract": source.ROW_CONTRACT,
            "precision_bits": source.PRECISION_BITS,
            "derivative_orders": [2, source.MAX_MOMENT],
        },
        "rows": rows,
        "summary": {"tiles": 2, "failed_tiles": 0, "covered_tiles": [["5691", "5692"], ["38028", "38029"]]},
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_h2_h24_boundary_extension.py",
    }


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote order-eleven compact H boundary extension: 2/2 tiles")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
