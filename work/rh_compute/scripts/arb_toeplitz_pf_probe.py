#!/usr/bin/env python
"""
FLINT/Arb propagation probe for finite Toeplitz/PF minors.

Given cached decimal A_k(lambda) values and an assumed absolute radius on
each A_k, build coefficient balls for either:

  beta:   a_k = A_k
  taylor: a_k = A_k / k!

Then propagate those balls through finite upper-triangular Toeplitz minors.

This is not a certificate for the exact analytic coefficients unless the
input balls are known to enclose the exact A_k(lambda). It is a propagation
gate for deciding whether a future interval coefficient generator would be
strong enough for the requested minor range.
"""

from __future__ import annotations

import argparse
from collections import Counter
import itertools
import json
import sys
from dataclasses import asdict, dataclass
from math import comb, factorial
from pathlib import Path


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
import mpmath as mp  # noqa: E402


@dataclass
class ProbeRow:
    kind: str
    lam: str
    sequence: str
    order: int
    rows: str
    cols: str
    abs_radius_on_A: str
    det: str
    classification: str
    ok_for_nonnegative: bool


def permutation_sign(perm: tuple[int, ...]) -> int:
    inversions = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inversions += 1
    return -1 if inversions % 2 else 1


def structural_terms(rows: tuple[int, ...], cols: tuple[int, ...]) -> dict[tuple[tuple[int, int], ...], int]:
    """Return symbolic determinant support for a Toeplitz minor.

    Each non-forced-zero entry is represented by a symbol a_d where
    d = col - row. If all Leibniz monomials cancel, the minor is
    identically zero for every coefficient sequence.
    """
    terms: dict[tuple[tuple[int, int], ...], int] = {}
    order = len(rows)
    for perm in itertools.permutations(range(order)):
        powers: Counter[int] = Counter()
        valid = True
        for i, j_perm in enumerate(perm):
            diff = cols[j_perm] - rows[i]
            if diff < 0:
                valid = False
                break
            powers[diff] += 1
        if not valid:
            continue
        key = tuple(sorted(powers.items()))
        terms[key] = terms.get(key, 0) + permutation_sign(perm)
        if terms[key] == 0:
            del terms[key]
    return terms


def is_structural_zero_exact(rows: tuple[int, ...], cols: tuple[int, ...]) -> bool:
    return not structural_terms(rows, cols)


def is_structural_zero_fast(rows: tuple[int, ...], cols: tuple[int, ...]) -> bool:
    """Fast structural-zero test for upper-triangular Toeplitz minors.

    For sorted row/column sets, the symbolic minor is nonzero exactly when
    c_i >= r_i for every i. This is the Jacobi-Trudi/skew-Schur support
    condition for Toeplitz minors, and matches the exact permutation
    expansion used by is_structural_zero_exact.
    """
    return any(col < row for row, col in zip(rows, cols))


def is_structural_zero(rows: tuple[int, ...], cols: tuple[int, ...], mode: str) -> bool:
    if mode == "fast":
        return is_structural_zero_fast(rows, cols)
    if mode == "exact":
        return is_structural_zero_exact(rows, cols)
    if mode == "validate":
        fast = is_structural_zero_fast(rows, cols)
        exact = is_structural_zero_exact(rows, cols)
        if fast != exact:
            raise RuntimeError(f"Structural classifier mismatch rows={rows} cols={cols}: fast={fast} exact={exact}")
        return fast
    raise ValueError(f"Unknown structural mode: {mode}")


def constrained_nonstructural_cols(rows: tuple[int, ...], size: int):
    """Yield sorted column sets satisfying cols[i] >= rows[i].

    This enumerates exactly the structurally nonzero upper-triangular
    Toeplitz minors under the fast structural criterion. It avoids looping
    over the much larger set of forced-zero column choices.
    """
    order = len(rows)
    current: list[int] = []

    def rec(pos: int, previous: int):
        if pos == order:
            yield tuple(current)
            return
        lower = max(previous + 1, rows[pos])
        upper = size - (order - pos)
        for col in range(lower, upper + 1):
            current.append(col)
            yield from rec(pos + 1, col)
            current.pop()

    yield from rec(0, -1)


def load_coefficients(path: Path, lam: str) -> list[str]:
    target = mp.mpf(lam)
    best: list[str] | None = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") != "coefficients":
                continue
            try:
                cached_lam = mp.mpf(row["lam"])
            except Exception:
                continue
            if cached_lam == target:
                values = list(row["A"])
                if best is None or len(values) > len(best):
                    best = values
    if best is None:
        raise RuntimeError(f"No coefficient row found for lambda={lam} in {path}")
    return best


def ball(value: str, abs_radius: str) -> flint.arb:
    return flint.arb(f"[{value} +/- {abs_radius}]")


def build_sequence(values: list[str], sequence: str, size: int, abs_radius: str) -> list[flint.arb]:
    if size > len(values):
        raise RuntimeError(f"Need {size} coefficients, cache has {len(values)}")
    out: list[flint.arb] = []
    for k, value in enumerate(values[:size]):
        ak = ball(value, abs_radius)
        if sequence == "beta":
            out.append(ak)
        elif sequence == "taylor":
            out.append(ak / factorial(k))
        else:
            raise ValueError(f"Unknown sequence {sequence!r}; expected beta or taylor")
    return out


def toeplitz_matrix(seq: list[flint.arb], size: int) -> flint.arb_mat:
    rows = []
    zero = flint.arb(0)
    for i in range(size):
        row = []
        for j in range(size):
            row.append(seq[j - i] if j >= i else zero)
        rows.append(row)
    return flint.arb_mat(rows)


def extract_matrix(matrix: flint.arb_mat, rows: tuple[int, ...], cols: tuple[int, ...]) -> flint.arb_mat:
    return flint.arb_mat([[matrix[i, j] for j in cols] for i in rows])


def classify(detv: flint.arb) -> str:
    if detv > 0:
        return "positive"
    if detv < 0:
        return "negative"
    if detv == flint.arb(0):
        return "zero"
    if detv.contains(0):
        return "inconclusive_contains_zero"
    return "unknown"


def probe_minor(
    matrix: flint.arb_mat,
    lam: str,
    sequence: str,
    rows: tuple[int, ...],
    cols: tuple[int, ...],
    abs_radius: str,
    structural_mode: str,
) -> ProbeRow:
    if is_structural_zero(rows, cols, structural_mode):
        return ProbeRow(
            kind="arb_toeplitz_minor",
            lam=lam,
            sequence=sequence,
            order=len(rows),
            rows=",".join(map(str, rows)),
            cols=",".join(map(str, cols)),
            abs_radius_on_A=abs_radius,
            det="0 (structural)",
            classification="structural_zero",
            ok_for_nonnegative=True,
        )
    if len(rows) == 1:
        detv = matrix[rows[0], cols[0]]
    else:
        detv = extract_matrix(matrix, rows, cols).det()
    cls = classify(detv)
    return ProbeRow(
        kind="arb_toeplitz_minor",
        lam=lam,
        sequence=sequence,
        order=len(rows),
        rows=",".join(map(str, rows)),
        cols=",".join(map(str, cols)),
        abs_radius_on_A=abs_radius,
        det=str(detv),
        classification=cls,
        ok_for_nonnegative=cls in {"positive", "zero"},
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coeff-cache", type=Path, required=True)
    parser.add_argument("--lambda", dest="lam", default="0")
    parser.add_argument("--sequence", choices=["beta", "taylor"], default="taylor")
    parser.add_argument("--matrix-size", type=int, default=6)
    parser.add_argument("--max-order", type=int, default=3)
    parser.add_argument("--abs-radius", default="1e-80")
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--record-limit", type=int, default=100)
    parser.add_argument("--structural-mode", choices=["fast", "exact", "validate"], default="fast")
    parser.add_argument(
        "--enumeration-mode",
        choices=["all", "nonzero"],
        default="all",
        help=(
            "all checks every row/column pair; nonzero enumerates only structurally "
            "nonzero minors and counts skipped structural zeros via the fast criterion"
        ),
    )
    parser.add_argument("--out-jsonl", type=Path, default=None)
    parser.add_argument("--out-summary", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.max_order < 1 or args.max_order > args.matrix_size:
        raise ValueError("--max-order must be between 1 and --matrix-size")
    if args.enumeration_mode == "nonzero" and args.structural_mode != "fast":
        raise ValueError("--enumeration-mode nonzero requires --structural-mode fast")
    flint.ctx.dps = args.dps

    values = load_coefficients(args.coeff_cache, args.lam)
    seq = build_sequence(values, args.sequence, args.matrix_size, args.abs_radius)
    matrix = toeplitz_matrix(seq, args.matrix_size)

    counts: dict[str, int] = {
        "positive": 0,
        "zero": 0,
        "structural_zero": 0,
        "negative": 0,
        "inconclusive_contains_zero": 0,
        "unknown": 0,
    }
    records: list[ProbeRow] = []
    tests = 0
    for order in range(1, args.max_order + 1):
        total_order_tests = comb(args.matrix_size, order) ** 2
        if args.enumeration_mode == "all":
            for rows in itertools.combinations(range(args.matrix_size), order):
                for cols in itertools.combinations(range(args.matrix_size), order):
                    row = probe_minor(matrix, args.lam, args.sequence, rows, cols, args.abs_radius, args.structural_mode)
                    counts[row.classification] = counts.get(row.classification, 0) + 1
                    tests += 1
                    if row.classification in {"negative", "inconclusive_contains_zero", "unknown"}:
                        if len(records) < args.record_limit:
                            records.append(row)
        else:
            nonstructural_tests = 0
            for rows in itertools.combinations(range(args.matrix_size), order):
                for cols in constrained_nonstructural_cols(rows, args.matrix_size):
                    row = probe_minor(matrix, args.lam, args.sequence, rows, cols, args.abs_radius, args.structural_mode)
                    counts[row.classification] = counts.get(row.classification, 0) + 1
                    nonstructural_tests += 1
                    if row.classification in {"negative", "inconclusive_contains_zero", "unknown"}:
                        if len(records) < args.record_limit:
                            records.append(row)
            counts["structural_zero"] = counts.get("structural_zero", 0) + (total_order_tests - nonstructural_tests)
            tests += total_order_tests

    summary = {
        "kind": "arb_toeplitz_pf_probe_summary",
        "lam": args.lam,
        "sequence": args.sequence,
        "matrix_size": args.matrix_size,
        "max_order": args.max_order,
        "abs_radius_on_A": args.abs_radius,
        "dps": args.dps,
        "structural_mode": args.structural_mode,
        "enumeration_mode": args.enumeration_mode,
        "tests": tests,
        "evaluated_nonstructural_tests": tests - counts.get("structural_zero", 0),
        "counts": counts,
        "ok_for_nonnegative": counts.get("negative", 0) == 0
        and counts.get("inconclusive_contains_zero", 0) == 0
        and counts.get("unknown", 0) == 0,
        "recorded_problem_rows": len(records),
    }

    if args.out_jsonl is not None:
        args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
        with args.out_jsonl.open("w", encoding="utf-8") as f:
            for row in records:
                f.write(json.dumps(asdict(row), sort_keys=True) + "\n")
    if args.out_summary is not None:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        with args.out_summary.open("w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2, sort_keys=True)
            f.write("\n")

    print(json.dumps(summary, sort_keys=True))
    return 0 if counts.get("negative", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
