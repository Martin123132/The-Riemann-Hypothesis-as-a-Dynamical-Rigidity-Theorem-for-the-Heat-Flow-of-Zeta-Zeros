#!/usr/bin/env python
"""
Count finite upper-triangular Toeplitz/PF probe sizes.

This is a planning utility. It counts:

  tests = sum_m binom(N,m)^2
  evaluated_nonstructural_tests = row/column pairs with c_i >= r_i
  structural_zero = tests - evaluated_nonstructural_tests

The nonstructural count is computed by a prefix DP rather than enumerating all
row/column pairs. At each index, we may choose no mark, row only, column only,
or both. The condition c_i >= r_i is equivalent to every prefix having at
least as many selected rows as selected columns.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
import json
from math import comb


def nonstructural_pairs(size: int, order: int) -> int:
    @lru_cache(maxsize=None)
    def dp(pos: int, rows_used: int, cols_used: int, balance: int) -> int:
        if pos == size:
            return int(rows_used == order and cols_used == order and balance == 0)

        total = dp(pos + 1, rows_used, cols_used, balance)
        if rows_used < order:
            total += dp(pos + 1, rows_used + 1, cols_used, balance + 1)
        if cols_used < order and balance > 0:
            total += dp(pos + 1, rows_used, cols_used + 1, balance - 1)
        if rows_used < order and cols_used < order:
            total += dp(pos + 1, rows_used + 1, cols_used + 1, balance)
        return total

    return dp(0, 0, 0, 0)


def frontier_counts(size: int, max_order: int) -> dict[str, object]:
    orders = []
    total_tests = 0
    total_evaluated = 0
    for order in range(1, max_order + 1):
        tests = comb(size, order) ** 2
        evaluated = nonstructural_pairs(size, order)
        structural_zero = tests - evaluated
        total_tests += tests
        total_evaluated += evaluated
        orders.append(
            {
                "order": order,
                "tests": tests,
                "evaluated_nonstructural_tests": evaluated,
                "structural_zero": structural_zero,
            }
        )
    return {
        "kind": "toeplitz_frontier_counts",
        "matrix_size": size,
        "max_order": max_order,
        "tests": total_tests,
        "evaluated_nonstructural_tests": total_evaluated,
        "structural_zero": total_tests - total_evaluated,
        "orders": orders,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-size", type=int, required=True)
    parser.add_argument("--max-order", type=int, required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.matrix_size < 1:
        raise ValueError("--matrix-size must be positive")
    if args.max_order < 1 or args.max_order > args.matrix_size:
        raise ValueError("--max-order must be between 1 and --matrix-size")
    print(json.dumps(frontier_counts(args.matrix_size, args.max_order), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
