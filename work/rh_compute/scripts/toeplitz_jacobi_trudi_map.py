#!/usr/bin/env python3
"""Map upper-triangular Toeplitz minors to Jacobi-Trudi data.

For sorted row and column sets R=(r_i), Q=(q_j), the Toeplitz minor

    det(c_{q_j-r_i})

is a Jacobi-Trudi determinant after normalizing d_k = c_k/c_0.  This
script records the exact finite reindexing.  It is a theorem-search
workbench, not a positivity proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
from dataclasses import asdict, dataclass
from math import comb


@dataclass(frozen=True)
class JacobiTrudiMap:
    rows: tuple[int, ...]
    cols: tuple[int, ...]
    structurally_nonzero: bool
    shift_k: int | None
    lambda_partition: tuple[int, ...]
    mu_partition: tuple[int, ...]
    jt_indices: tuple[tuple[int, ...], ...]
    scale_factor: str
    note: str


def parse_int_tuple(text: str) -> tuple[int, ...]:
    values = tuple(int(part.strip()) for part in text.split(",") if part.strip())
    if tuple(sorted(values)) != values or len(set(values)) != len(values):
        raise ValueError(f"expected a strictly increasing comma list, got {text!r}")
    if values and values[0] < 0:
        raise ValueError("row/column indices must be nonnegative")
    return values


def is_partition(parts: tuple[int, ...]) -> bool:
    return all(part >= 0 for part in parts) and all(parts[i] >= parts[i + 1] for i in range(len(parts) - 1))


def structural_nonzero(rows: tuple[int, ...], cols: tuple[int, ...]) -> bool:
    return len(rows) == len(cols) and all(col >= row for row, col in zip(rows, cols))


def map_minor(rows: tuple[int, ...], cols: tuple[int, ...]) -> JacobiTrudiMap:
    if len(rows) != len(cols):
        raise ValueError("rows and cols must have the same length")
    if not rows:
        raise ValueError("minor order must be positive")

    nonzero = structural_nonzero(rows, cols)
    if not nonzero:
        return JacobiTrudiMap(
            rows=rows,
            cols=cols,
            structurally_nonzero=False,
            shift_k=None,
            lambda_partition=tuple(),
            mu_partition=tuple(),
            jt_indices=tuple(),
            scale_factor="0",
            note="structural zero: at least one sorted column index is smaller than the matching row index",
        )

    # Jacobi-Trudi uses 1-based matrix indices:
    #   h_{lambda_i - mu_j - i + j}.
    # To match q_j-r_i, choose lambda_i=K+i-r_i and mu_i=K+i-q_i.
    # K only shifts the skew diagram horizontally; choose the minimal
    # nonnegative value making both partitions nonnegative.
    candidates = [0]
    for i, row in enumerate(rows, start=1):
        candidates.append(row - i)
    for j, col in enumerate(cols, start=1):
        candidates.append(col - j)
    shift_k = max(candidates)
    lam = tuple(shift_k + i - row for i, row in enumerate(rows, start=1))
    mu = tuple(shift_k + j - col for j, col in enumerate(cols, start=1))
    if not is_partition(lam):
        raise AssertionError(f"lambda is not a partition: {lam}")
    if not is_partition(mu):
        raise AssertionError(f"mu is not a partition: {mu}")
    if any(l_part < m_part for l_part, m_part in zip(lam, mu)):
        raise AssertionError(f"not a skew partition: lambda={lam}, mu={mu}")

    indices: list[tuple[int, ...]] = []
    for i, l_part in enumerate(lam, start=1):
        row_indices = []
        for j, m_part in enumerate(mu, start=1):
            row_indices.append(l_part - m_part - i + j)
        indices.append(tuple(row_indices))
    expected = tuple(tuple(col - row for col in cols) for row in rows)
    if tuple(indices) != expected:
        raise AssertionError(f"Jacobi-Trudi indices {tuple(indices)} != Toeplitz differences {expected}")

    order = len(rows)
    return JacobiTrudiMap(
        rows=rows,
        cols=cols,
        structurally_nonzero=True,
        shift_k=shift_k,
        lambda_partition=lam,
        mu_partition=mu,
        jt_indices=tuple(indices),
        scale_factor=f"c_0^{order}",
        note="det(c_{q_j-r_i}) = c_0^m * s_{lambda/mu}(d), d_k=c_k/c_0, h_k -> d_k",
    )


def validate_grid(matrix_size: int, max_order: int) -> dict:
    total = 0
    nonzero = 0
    structural_zero = 0
    mapped = 0
    examples: list[dict] = []
    for order in range(1, max_order + 1):
        for rows in itertools.combinations(range(matrix_size), order):
            total_row_choices = comb(matrix_size, order)
            if total_row_choices <= 0:
                raise AssertionError("comb returned nonpositive")
            for cols in itertools.combinations(range(matrix_size), order):
                total += 1
                mapped_row = map_minor(rows, cols)
                if mapped_row.structurally_nonzero:
                    nonzero += 1
                    mapped += 1
                    if len(examples) < 8 and rows != cols:
                        examples.append(asdict(mapped_row))
                else:
                    structural_zero += 1
    return {
        "kind": "toeplitz_jacobi_trudi_map_validation",
        "matrix_size": matrix_size,
        "max_order": max_order,
        "tests": total,
        "structurally_nonzero": nonzero,
        "structural_zero": structural_zero,
        "mapped_nonzero_minors": mapped,
        "examples": examples,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", help="comma-separated row indices, e.g. 0,1,3")
    parser.add_argument("--cols", help="comma-separated column indices, e.g. 1,3,4")
    parser.add_argument("--matrix-size", type=int, help="validate all minors up to this matrix size")
    parser.add_argument("--max-order", type=int, default=4, help="maximum order for --matrix-size validation")
    args = parser.parse_args()

    if args.rows is not None or args.cols is not None:
        if args.rows is None or args.cols is None:
            raise SystemExit("--rows and --cols must be supplied together")
        print(json.dumps(asdict(map_minor(parse_int_tuple(args.rows), parse_int_tuple(args.cols))), indent=2))
        return 0

    if args.matrix_size is None:
        raise SystemExit("supply either --rows/--cols or --matrix-size")
    if args.matrix_size <= 0 or args.max_order <= 0:
        raise SystemExit("--matrix-size and --max-order must be positive")
    if args.max_order > args.matrix_size:
        raise SystemExit("--max-order cannot exceed --matrix-size")
    print(json.dumps(validate_grid(args.matrix_size, args.max_order), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
