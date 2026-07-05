#!/usr/bin/env python3
"""Validate the exact Toeplitz-minor to Jacobi-Trudi reindexing gate.

This is an algebraic theorem-search guard for the coefficient-PF route.  It
does not prove Schur positivity, PF-infinity, RH, or Lambda <= 0.  It checks
that the exact finite reindexing used by the bridge notes is reproducible and
matches the structural-zero support condition used by the Toeplitz/PF probes.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from toeplitz_jacobi_trudi_map import map_minor, validate_grid  # noqa: E402


EXPECTED_GRID = {
    "matrix_size": 10,
    "max_order": 5,
    "tests": 124_129,
    "structurally_nonzero": 39_094,
    "structural_zero": 85_035,
    "mapped_nonzero_minors": 39_094,
}


def require_equal(label: str, actual: object, expected: object) -> None:
    if actual != expected:
        raise AssertionError(f"{label}: {actual!r} != {expected!r}")


def validate_grid_counts(matrix_size: int, max_order: int) -> dict:
    summary = validate_grid(matrix_size, max_order)
    expected = dict(EXPECTED_GRID)
    expected["matrix_size"] = matrix_size
    expected["max_order"] = max_order
    if matrix_size == EXPECTED_GRID["matrix_size"] and max_order == EXPECTED_GRID["max_order"]:
        for key, value in EXPECTED_GRID.items():
            require_equal(key, summary.get(key), value)
    else:
        require_equal("nonzero + zero", summary["structurally_nonzero"] + summary["structural_zero"], summary["tests"])
        require_equal("mapped nonzero", summary["mapped_nonzero_minors"], summary["structurally_nonzero"])
    return summary


def validate_examples() -> None:
    mapped = asdict(map_minor((0, 1), (1, 2)))
    require_equal("example nonzero flag", mapped["structurally_nonzero"], True)
    require_equal("example lambda", mapped["lambda_partition"], (1, 1))
    require_equal("example mu", mapped["mu_partition"], (0, 0))
    require_equal("example indices", mapped["jt_indices"], ((1, 2), (0, 1)))
    require_equal("example scale", mapped["scale_factor"], "c_0^2")

    zero = asdict(map_minor((2, 4), (1, 5)))
    require_equal("zero flag", zero["structurally_nonzero"], False)
    require_equal("zero scale", zero["scale_factor"], "0")
    require_equal("zero indices", zero["jt_indices"], ())


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-size", type=int, default=EXPECTED_GRID["matrix_size"])
    parser.add_argument("--max-order", type=int, default=EXPECTED_GRID["max_order"])
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.matrix_size <= 0 or args.max_order <= 0:
        raise ValueError("matrix size and max order must be positive")
    if args.max_order > args.matrix_size:
        raise ValueError("max order cannot exceed matrix size")

    validate_examples()
    summary = validate_grid_counts(args.matrix_size, args.max_order)
    print(
        "OK Toeplitz/Jacobi-Trudi example maps: "
        "one nonzero skew shape and one structural zero"
    )
    print(
        "validated Toeplitz/Jacobi-Trudi reindexing: "
        f"N={summary['matrix_size']}, orders<={summary['max_order']}, "
        f"{summary['tests']} minors, "
        f"{summary['structurally_nonzero']} nonzero maps, "
        f"{summary['structural_zero']} structural zeros"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
