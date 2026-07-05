#!/usr/bin/env python3
"""Arb reshaped-Hankel sign-consistency certificate probe.

This tests the Grussler-Damm reshaped-Hankel finite condition using rigorous
`A_ball` coefficient enclosures.  For each order k, it checks all k-column
minors of the k-row reshaped Hankel matrix with entries A_{i+j}.

The result is a finite interval certificate only.  It is not an all-order
sign-consistency theorem, not a Jensen hyperbolicity theorem, and not a proof
of RH or Lambda <= 0.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_LAMBDAS = ("0", "1e-6", "1e-4", "1e-2", "1e-1")
DEFAULT_ENCLOSURE_JSONL = (
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
)


@dataclass(frozen=True)
class ProbeSubsummary:
    lam: str
    order_k: int
    n_cols: int
    max_coefficient_index: int
    expected_sign: int
    tests: int
    positive: int
    negative: int
    zero: int
    inconclusive: int
    first_bad_columns: tuple[int, ...] | None
    all_ok: bool


@dataclass(frozen=True)
class ProbeRow:
    kind: str
    lam: str
    order_k: int
    n_cols: int
    columns: tuple[int, ...]
    expected_sign: int
    classification: str
    contains_zero: bool
    ok: bool


def dec(value: object) -> Decimal:
    return Decimal(str(value))


def decimal_equal(left: str, right: str) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


def signed_hankel_sign(order_k: int) -> int:
    return -1 if (order_k * (order_k - 1) // 2) % 2 else 1


def classify(value: flint.arb) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    if value == flint.arb(0):
        return "zero"
    if value.contains(0):
        return "inconclusive_contains_zero"
    return "unknown"


def load_enclosure_balls(paths: list[Path], lam: str, needed_max_k: int) -> dict[int, flint.arb]:
    balls: dict[int, flint.arb] = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                if row.get("kind") != "acb_coefficient_enclosure":
                    continue
                if not decimal_equal(row["lam"], lam):
                    continue
                k = int(row["k"])
                if k <= needed_max_k:
                    balls[k] = flint.arb(row["A_ball"])
    missing = [k for k in range(needed_max_k + 1) if k not in balls]
    if missing:
        raise RuntimeError(f"missing A_k enclosure balls for lambda={lam}: {missing[:10]}")
    return balls


def probe_one(
    balls: dict[int, flint.arb],
    lam: str,
    order_k: int,
    n_cols: int,
    columns: tuple[int, ...],
) -> ProbeRow:
    matrix = flint.arb_mat(
        [[balls[i + column] for column in columns] for i in range(order_k)]
    )
    detv = matrix.det()
    expected = signed_hankel_sign(order_k)
    signed = detv if expected > 0 else -detv
    classification = classify(signed)
    contains_zero = signed.contains(0)
    return ProbeRow(
        kind="arb_hankel_sign_consistency_reduction_row",
        lam=lam,
        order_k=order_k,
        n_cols=n_cols,
        columns=columns,
        expected_sign=expected,
        classification=classification,
        contains_zero=contains_zero,
        ok=(classification == "positive" and not contains_zero),
    )


def probe_block(
    balls: dict[int, flint.arb],
    lam: str,
    order_k: int,
    n_cols: int,
    row_handle,
) -> ProbeSubsummary:
    counts: Counter[str] = Counter()
    first_bad_columns: tuple[int, ...] | None = None
    expected = signed_hankel_sign(order_k)
    for columns in combinations(range(n_cols), order_k):
        row = probe_one(balls, lam, order_k, n_cols, tuple(columns))
        counts[row.classification] += 1
        if not row.ok and first_bad_columns is None:
            first_bad_columns = row.columns
        if row_handle is not None:
            row_handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    tests = comb(n_cols, order_k)
    positive = counts.get("positive", 0)
    negative = counts.get("negative", 0)
    zero = counts.get("zero", 0)
    inconclusive = sum(
        count
        for key, count in counts.items()
        if key not in {"positive", "negative", "zero"}
    )
    return ProbeSubsummary(
        lam=lam,
        order_k=order_k,
        n_cols=n_cols,
        max_coefficient_index=n_cols + order_k - 2,
        expected_sign=expected,
        tests=tests,
        positive=positive,
        negative=negative,
        zero=zero,
        inconclusive=inconclusive,
        first_bad_columns=first_bad_columns,
        all_ok=(positive == tests and negative == 0 and zero == 0 and inconclusive == 0),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=[Path(p) for p in DEFAULT_ENCLOSURE_JSONL])
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--orders", default="2,3,4,5")
    parser.add_argument("--n-cols", type=int, default=18)
    parser.add_argument("--dps", type=int, default=520)
    parser.add_argument("--out-jsonl", type=Path, default=None)
    parser.add_argument("--out-summary", type=Path, default=None)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = args.dps
    lambdas = [part.strip() for part in args.lambdas.split(",") if part.strip()]
    orders = [int(part.strip()) for part in args.orders.split(",") if part.strip()]
    needed_max_k = args.n_cols + max(orders) - 2

    if args.out_jsonl is not None:
        args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
        row_handle = args.out_jsonl.open("w", encoding="utf-8")
    else:
        row_handle = None

    summaries: list[ProbeSubsummary] = []
    try:
        for lam in lambdas:
            balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
            for order_k in orders:
                summaries.append(probe_block(balls, lam, order_k, args.n_cols, row_handle))
    finally:
        if row_handle is not None:
            row_handle.close()

    rows = sum(summary.tests for summary in summaries)
    ok_rows = sum(summary.positive for summary in summaries if summary.all_ok)
    failed = [summary for summary in summaries if not summary.all_ok]
    payload = {
        "kind": "arb_hankel_sign_consistency_reduction_summary",
        "date": "2026-07-05",
        "proof_boundary": (
            "Finite Arb coefficient-enclosure certificate only; not all-order "
            "sign consistency, not a Jensen hyperbolicity theorem, and not a "
            "proof of RH or Lambda <= 0."
        ),
        "enclosure_jsonl": [str(path) for path in args.enclosure_jsonl],
        "lambdas": lambdas,
        "orders": orders,
        "n_cols": args.n_cols,
        "dps": args.dps,
        "needed_max_k": needed_max_k,
        "summary_rows": len(summaries),
        "rows": rows,
        "ok": rows - sum(summary.tests - summary.positive for summary in summaries),
        "failed_or_inconclusive": sum(summary.tests - summary.positive for summary in summaries),
        "all_ok": not failed,
        "summaries": [asdict(summary) for summary in summaries],
    }
    if args.out_summary is not None:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        for summary in summaries:
            status = "OK" if summary.all_ok else "FAIL"
            print(
                f"{status} Arb reshaped-Hankel sign-consistency: "
                f"lambda={summary.lam}, k={summary.order_k}, N={summary.n_cols}, "
                f"{summary.positive}/{summary.tests} positive"
            )
        print(
            "validated "
            f"{sum(1 for summary in summaries if summary.all_ok)}/{len(summaries)} "
            "Arb reshaped-Hankel sign-consistency certificate rows"
        )
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
