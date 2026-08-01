#!/usr/bin/env python3
"""Build four interleaved exact H0-H23 anchors for the compact step-four pilot."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
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


DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order11_compact_step4_pilot_anchors.json"
TARGETS = (Fraction(20000), Fraction(20008), Fraction(38016), Fraction(38024))


def tasks() -> list[tuple[int, Fraction, str]]:
    return [(index, target, "far") for index, target in enumerate(TARGETS)]


def build_artifact(*, workers: int = 4) -> dict:
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(source.exact_task, tasks()))
    if any(row.get("passed") is not True for row in rows):
        raise RuntimeError(f"a step-four pilot anchor failed: {rows}")
    rows = [
        {**row, "kind": "order11_compact_step4_pilot_h0_h23_jet"}
        for row in rows
    ]
    return {
        "kind": "jensen_window_pf_compound_order11_compact_step4_pilot_anchors",
        "date": "2026-07-18",
        "status": "rigorous exact four-row H0-H23 compact step-four pilot source",
        "proof_boundary": "These four exact jets test source density only and do not cover the compact interval.",
        "source_contract": {"row_contract": source.row_contract("far"), "derivative_orders": [0, 23]},
        "rows": rows,
        "summary": {"rows": 4, "failed_rows": 0, "targets": [str(t) for t in TARGETS]},
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_step4_pilot_anchors.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_step4_pilot_anchors.py",
    }


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote order-eleven compact step-four pilot anchors: 4/4")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
