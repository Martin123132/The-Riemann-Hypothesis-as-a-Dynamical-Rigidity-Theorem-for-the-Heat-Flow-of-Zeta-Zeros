#!/usr/bin/env python
"""
Finite Toeplitz/PF minor audit for coefficient-cache rows.

This script tests the Aissen-Schoenberg-Whitney style upper-triangular
Toeplitz matrix T[i,j] = a[j-i] for j >= i, otherwise 0.

It is a finite exact-rational audit of numerically computed coefficients,
not a rigorous certificate for the underlying analytic coefficients.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import itertools
import json
from dataclasses import asdict, dataclass
from math import factorial
from pathlib import Path
from time import perf_counter

import mpmath as mp
import sympy as sp


@dataclass
class MinorRow:
    kind: str
    sequence: str
    lam: str
    order: int
    rows: str
    cols: str
    sign: int
    log10_abs_det: str
    det_decimal: str
    rational_digits: int
    matrix_size: int
    source_max_k: int


def parse_ints(text: str) -> tuple[int, ...]:
    if not text.strip():
        return tuple()
    return tuple(int(part.strip()) for part in text.split(",") if part.strip())


def parse_sequence_names(text: str) -> list[str]:
    names = [part.strip().lower() for part in text.split(",") if part.strip()]
    allowed = {"beta", "taylor"}
    bad = sorted(set(names) - allowed)
    if bad:
        raise ValueError(f"Unknown sequence(s): {', '.join(bad)}. Allowed: beta,taylor")
    return names


def rationalize_mpf_text(text: str, digits: int) -> sp.Rational:
    value = mp.mpf(text)
    return sp.Rational(mp.nstr(value, digits, min_fixed=-10, max_fixed=10))


def load_coeff_row(path: Path, lam: mp.mpf, needed_max_k: int) -> dict:
    best: dict | None = None
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
            if cached_lam != lam:
                continue
            max_k = int(row["max_k"])
            if max_k < needed_max_k:
                continue
            if best is None or max_k < int(best["max_k"]):
                best = row
    if best is None:
        raise RuntimeError(
            f"No coefficient-cache row for lambda={mp.nstr(lam, 30)} with max_k >= {needed_max_k}"
        )
    return best


def build_sequence(row: dict, name: str, needed_max_k: int, digits: int) -> list[sp.Rational]:
    beta = [rationalize_mpf_text(text, digits) for text in row["A"][: needed_max_k + 1]]
    if name == "beta":
        return beta
    if name == "taylor":
        return [value / factorial(k) for k, value in enumerate(beta)]
    raise ValueError(name)


def toeplitz_matrix(seq: list[sp.Rational], size: int) -> sp.Matrix:
    return sp.Matrix(
        size,
        size,
        lambda i, j: seq[j - i] if j >= i else sp.Rational(0),
    )


def permutation_sign(perm: tuple[int, ...]) -> int:
    inversions = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                inversions += 1
    return -1 if inversions % 2 else 1


def structural_terms(rows: tuple[int, ...], cols: tuple[int, ...]) -> dict[tuple[tuple[int, int], ...], int]:
    """Return the determinant polynomial support for a symbolic Toeplitz minor.

    Entries are independent symbols a_d with forced zeros when d < 0.
    The returned dictionary maps monomial exponent patterns to integer
    coefficients after Leibniz-term cancellation.
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
    """Fast structural-zero test for upper-triangular Toeplitz minors."""
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


def sign(value: sp.Rational) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def log10_abs(value: sp.Rational, digits: int = 24) -> str:
    if value == 0:
        return "-inf"
    mp.mp.dps = max(50, digits + 20)
    q = sp.Rational(value)
    out = mp.log10(mp.mpf(str(abs(q.p)))) - mp.log10(mp.mpf(str(q.q)))
    return mp.nstr(out, digits)


def det_decimal(value: sp.Rational, digits: int = 24) -> str:
    if value == 0:
        return "0"
    return str(sp.N(value, digits))


def minor_det(matrix: sp.Matrix, rows: tuple[int, ...], cols: tuple[int, ...]) -> sp.Rational:
    if len(rows) == 1:
        return sp.Rational(matrix[rows[0], cols[0]])
    sub = matrix.extract(rows, cols)
    return sp.Rational(sub.det(method="bareiss"))


def write_jsonl(path: Path, rows: list[dict], overwrite: bool) -> None:
    mode = "w" if overwrite else "a"
    with path.open(mode, encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict], overwrite: bool) -> None:
    if not rows:
        return
    exists = path.exists() and not overwrite
    fieldnames = sorted({key for row in rows for key in row.keys()})
    mode = "a" if exists else "w"
    with path.open(mode, newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def row_for_minor(
    sequence: str,
    lam: mp.mpf,
    order: int,
    rows: tuple[int, ...],
    cols: tuple[int, ...],
    detv: sp.Rational,
    args: argparse.Namespace,
    source_max_k: int,
) -> MinorRow:
    return MinorRow(
        kind="toeplitz_minor",
        sequence=sequence,
        lam=mp.nstr(lam, 30),
        order=order,
        rows=",".join(map(str, rows)),
        cols=",".join(map(str, cols)),
        sign=sign(detv),
        log10_abs_det=log10_abs(detv),
        det_decimal=det_decimal(detv),
        rational_digits=args.rational_digits,
        matrix_size=args.matrix_size,
        source_max_k=source_max_k,
    )


def audit_sequence(sequence: str, coeff_row: dict, args: argparse.Namespace, lam: mp.mpf) -> dict:
    source_max_k = int(coeff_row["max_k"])
    seq = build_sequence(coeff_row, sequence, args.matrix_size - 1, args.rational_digits)
    matrix = toeplitz_matrix(seq, args.matrix_size)
    start = perf_counter()

    known_rows = parse_ints(args.known_rows)
    known_cols = parse_ints(args.known_cols)
    known: dict | None = None
    if known_rows and known_cols:
        if len(known_rows) != len(known_cols):
            raise ValueError("--known-rows and --known-cols must have the same length")
        if max(known_rows + known_cols) >= args.matrix_size:
            raise ValueError("Known minor indices exceed --matrix-size")
        detv = minor_det(matrix, known_rows, known_cols)
        known = asdict(
            row_for_minor(sequence, lam, len(known_rows), known_rows, known_cols, detv, args, source_max_k)
        )

    sign_counts = {str(order): {"positive": 0, "zero": 0, "negative": 0} for order in range(1, args.max_order + 1)}
    zero_structure_counts = {
        str(order): {"structural_zero": 0, "point_zero_nonstructural": 0, "nonzero": 0}
        for order in range(1, args.max_order + 1)
    }
    tests = 0
    first_negative: dict | None = None
    most_negative: dict | None = None
    smallest_nonzero: dict | None = None
    negative_rows: list[dict] = []
    point_zero_nonstructural_rows: list[dict] = []

    for order in range(1, args.max_order + 1):
        for rows in itertools.combinations(range(args.matrix_size), order):
            for cols in itertools.combinations(range(args.matrix_size), order):
                structural_zero = is_structural_zero(rows, cols, args.structural_mode)
                detv = minor_det(matrix, rows, cols)
                sig = sign(detv)
                tests += 1
                if sig > 0:
                    sign_counts[str(order)]["positive"] += 1
                    zero_structure_counts[str(order)]["nonzero"] += 1
                elif sig == 0:
                    sign_counts[str(order)]["zero"] += 1
                    if structural_zero:
                        zero_structure_counts[str(order)]["structural_zero"] += 1
                    else:
                        zero_structure_counts[str(order)]["point_zero_nonstructural"] += 1
                        if len(point_zero_nonstructural_rows) < args.negative_limit:
                            point_zero_nonstructural_rows.append(
                                asdict(row_for_minor(sequence, lam, order, rows, cols, detv, args, source_max_k))
                            )
                else:
                    sign_counts[str(order)]["negative"] += 1
                    zero_structure_counts[str(order)]["nonzero"] += 1
                    minor = asdict(row_for_minor(sequence, lam, order, rows, cols, detv, args, source_max_k))
                    if first_negative is None:
                        first_negative = minor
                    if most_negative is None or detv < sp.Rational(most_negative["_det_exact"]):
                        most_negative = dict(minor)
                        most_negative["_det_exact"] = str(detv)
                    if len(negative_rows) < args.negative_limit:
                        negative_rows.append(minor)
                if sig != 0:
                    minor_abs = abs(detv)
                    if smallest_nonzero is None or minor_abs < sp.Rational(smallest_nonzero["_abs_exact"]):
                        smallest_nonzero = asdict(
                            row_for_minor(sequence, lam, order, rows, cols, detv, args, source_max_k)
                        )
                        smallest_nonzero["_abs_exact"] = str(minor_abs)

    if most_negative is not None:
        most_negative.pop("_det_exact", None)
    if smallest_nonzero is not None:
        smallest_nonzero.pop("_abs_exact", None)

    elapsed = perf_counter() - start
    return {
        "kind": "toeplitz_pf_audit_summary",
        "sequence": sequence,
        "lam": mp.nstr(lam, 30),
        "matrix_size": args.matrix_size,
        "max_order": args.max_order,
        "rational_digits": args.rational_digits,
        "structural_mode": args.structural_mode,
        "source_cache": str(args.coeff_cache),
        "source_max_k": source_max_k,
        "source_dps": coeff_row.get("dps"),
        "source_n_sum": coeff_row.get("n_sum"),
        "source_cutoff": coeff_row.get("cutoff"),
        "tests": tests,
        "sign_counts": sign_counts,
        "zero_structure_counts": zero_structure_counts,
        "known_minor": known,
        "first_negative": first_negative,
        "most_negative": most_negative,
        "smallest_nonzero": smallest_nonzero,
        "negative_examples": negative_rows,
        "point_zero_nonstructural_examples": point_zero_nonstructural_rows,
        "elapsed_seconds": elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coeff-cache", type=Path, required=True)
    parser.add_argument("--lam", default="0")
    parser.add_argument("--sequences", default="beta,taylor", help="Comma list: beta,taylor")
    parser.add_argument("--matrix-size", type=int, default=10)
    parser.add_argument("--max-order", type=int, default=4)
    parser.add_argument("--rational-digits", type=int, default=70)
    parser.add_argument("--known-rows", default="0,1,2")
    parser.add_argument("--known-cols", default="1,2,3")
    parser.add_argument("--negative-limit", type=int, default=20)
    parser.add_argument("--structural-mode", choices=["fast", "exact", "validate"], default="fast")
    parser.add_argument("--output-dir", type=Path, default=Path("work/rh_compute/results"))
    parser.add_argument("--run-id", default="toeplitz_pf_audit")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.matrix_size < 1:
        raise ValueError("--matrix-size must be positive")
    if args.max_order < 1 or args.max_order > args.matrix_size:
        raise ValueError("--max-order must be between 1 and --matrix-size")
    if args.rational_digits < 10:
        raise ValueError("--rational-digits should be at least 10")

    lam = mp.mpf(args.lam)
    sequences = parse_sequence_names(args.sequences)
    needed_max_k = args.matrix_size - 1
    coeff_row = load_coeff_row(args.coeff_cache, lam, needed_max_k)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = [audit_sequence(sequence, coeff_row, args, lam) for sequence in sequences]

    summary_path = args.output_dir / f"{args.run_id}_summary.json"
    jsonl_path = args.output_dir / f"{args.run_id}_negative_examples.jsonl"
    csv_path = args.output_dir / f"{args.run_id}_negative_examples.csv"

    with summary_path.open("w" if args.overwrite else "x", encoding="utf-8") as f:
        json.dump({"summaries": summaries}, f, indent=2, sort_keys=True)
        f.write("\n")

    negative_rows = [row for summary in summaries for row in summary["negative_examples"]]
    write_jsonl(jsonl_path, negative_rows, overwrite=args.overwrite)
    write_csv(csv_path, negative_rows, overwrite=args.overwrite)

    for summary in summaries:
        negatives = sum(counts["negative"] for counts in summary["sign_counts"].values())
        known = summary["known_minor"]
        known_part = ""
        if known is not None:
            known_part = f", known_sign={known['sign']}, known_log10={known['log10_abs_det']}"
        print(
            f"{summary['sequence']}: tests={summary['tests']}, negatives={negatives}"
            f"{known_part}, elapsed={summary['elapsed_seconds']:.2f}s"
        )
    print(f"summary={summary_path}")
    print(f"negative_examples={jsonl_path}")


if __name__ == "__main__":
    main()
