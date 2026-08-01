#!/usr/bin/env python3
"""Run the high-precision refinement of the compact quadratic matrix."""

from __future__ import annotations

from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import jensen_window_pf_newman_theta_modular_tail_arb_quadratic_matrix as matrix


matrix.STEM = (
    "jensen_window_pf_newman_theta_modular_tail_arb_quadratic_refinement"
)
matrix.SCHEMA = "newman_theta_modular_tail_arb_quadratic_refinement_v1"
matrix.PRECISION_BITS = 192
matrix.N_VALUES = (7, 8, 9, 10)
matrix.INTEGRATION_OPTIONS = {
    "deg_limit": 100,
    "eval_limit": 200_000,
    "depth_limit": 35,
    "use_heap": True,
    "abs_tol": "2^-320",
    "rel_tol": "2^-120",
}
matrix.DEFAULT_CACHE = (
    matrix.RESULT_DIR / f"{matrix.STEM}.jsonl"
)
matrix.DEFAULT_OUT = matrix.RESULT_DIR / f"{matrix.STEM}.json"
matrix.DEFAULT_NOTE = matrix.REPO_ROOT / "outputs" / f"{matrix.STEM}.md"


if __name__ == "__main__":
    raise SystemExit(matrix.main())
