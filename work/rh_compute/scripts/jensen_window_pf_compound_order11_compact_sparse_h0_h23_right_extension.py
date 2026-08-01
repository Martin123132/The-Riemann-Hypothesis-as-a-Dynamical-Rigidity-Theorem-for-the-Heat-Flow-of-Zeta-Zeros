#!/usr/bin/env python3
"""Build the exact H0-H23 anchor at t=38028 for the order-eleven compact bridge."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
for candidate in (SCRIPT_DIR, REPO_ROOT / "work/rh_compute/vendor"):
    if str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import jensen_window_pf_compound_order10_compact_sparse_point_h0_h23_cache as source  # noqa: E402


DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order11_compact_sparse_h0_h23_right_extension.json"
TASK = (len(source.deterministic_tasks()), Fraction(38028), "far")


def build_artifact() -> dict:
    row = source.exact_task(TASK)
    if row.get("passed") is not True:
        raise RuntimeError(f"right exact anchor failed: {row}")
    row = {**row, "kind": "order11_compact_sparse_point_h0_h23_right_extension_jet"}
    return {
        "kind": "jensen_window_pf_compound_order11_compact_sparse_h0_h23_right_extension",
        "date": "2026-07-18",
        "status": "rigorous exact H0-H23 compact right-boundary anchor",
        "proof_boundary": "This artifact certifies one exact first-summand H jet only.",
        "source_contract": {
            "inherited_row_contract": source.row_contract("far"),
            "target_t": "38028",
            "derivative_orders": [0, source.MAX_MOMENT],
        },
        "row": row,
        "summary": {"rows": 1, "failed_rows": 0},
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_sparse_h0_h23_right_extension.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_sparse_h0_h23_right_extension.py",
    }


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote order-eleven compact H0-H23 right extension: t=38028")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
