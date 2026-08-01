#!/usr/bin/env python3
"""Build a high-precision exact H0-H23 stencil around t=20004."""

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


DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_compound_order11_compact_high_precision_stencil_pilot.json"
TARGETS = tuple(Fraction(value) for value in range(19995, 20014))
PROFILE_NAME = "order11_stress"
PROFILE = {
    "precision_bits": 1536,
    "mode_bisections": 320,
    "panels": 128,
    "row_contract": "order11-compact-stress-h0-h23-p1536-b320-n128-w15-t30-v1",
}


def stress_task(task: tuple[int, Fraction]) -> dict:
    index, target = task
    source.profiles.PROFILE_SPECS[PROFILE_NAME] = PROFILE
    row = source.exact_task((index, target, PROFILE_NAME))
    return {**row, "kind": "order11_compact_high_precision_stencil_h0_h23_jet"}


def tasks() -> list[tuple[int, Fraction]]:
    return list(enumerate(TARGETS))


def build_artifact(*, workers: int = 4) -> dict:
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(stress_task, tasks()))
    if any(row.get("passed") is not True for row in rows):
        raise RuntimeError("a high-precision stencil row failed")
    return {
        "kind": "jensen_window_pf_compound_order11_compact_high_precision_stencil_pilot",
        "date": "2026-07-18",
        "status": "rigorous exact nineteen-row high-precision H0-H23 stencil pilot",
        "proof_boundary": "This local stencil diagnoses source precision and does not cover the compact interval.",
        "source_contract": {"profile": PROFILE, "derivative_orders": [0, 23]},
        "rows": rows,
        "summary": {"rows": len(rows), "failed_rows": 0, "targets": [str(t) for t in TARGETS]},
        "generator": "work/rh_compute/scripts/jensen_window_pf_compound_order11_compact_high_precision_stencil_pilot.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_compound_order11_compact_high_precision_stencil_pilot.py",
    }


def main() -> int:
    artifact = build_artifact()
    DEFAULT_OUT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_OUT.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote order-eleven high-precision compact stencil: 19/19")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
