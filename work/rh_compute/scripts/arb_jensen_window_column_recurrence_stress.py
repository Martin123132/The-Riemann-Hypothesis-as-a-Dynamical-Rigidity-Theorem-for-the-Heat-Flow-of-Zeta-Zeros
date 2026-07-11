#!/usr/bin/env python3
"""Arb stress test for Jensen-window column recurrence positivity.

For h_j = binom(d,j) A_{n+j}(lambda), the column determinant sequence C_m
satisfies

    C[m] = sum_j (-1)^(j-1) h_0^(j-1) h_j C[m-j].

This script checks that necessary column-recurrence condition on a finite
lambda/shift/degree/size grid using existing Arb coefficient enclosures.  It is
finite evidence only, not an all-order recurrence theorem.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal
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
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl",
)


@dataclass(frozen=True)
class RecurrenceStressRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    recurrence_size_m: int
    max_coefficient_index: int
    classification: str
    contains_zero: bool
    ok: bool


@dataclass(frozen=True)
class RecurrenceStressSubsummary:
    lam: str
    shift_n: int
    degree_d: int
    sizes: tuple[int, ...]
    tests: int
    positive: int
    negative: int
    zero: int
    inconclusive: int
    first_bad_size: int | None
    all_ok: bool


def decimal_equal(left: str, right: str) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


def parse_int_range(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ".." in part:
            left, right = part.split("..", 1)
            start = int(left.strip())
            stop = int(right.strip())
            if stop < start:
                raise ValueError(f"descending range is not supported: {part}")
            out.extend(range(start, stop + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


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


def recurrence_values(
    balls: dict[int, flint.arb],
    shift_n: int,
    degree_d: int,
    max_size_m: int,
) -> list[flint.arb]:
    h = [comb(degree_d, index) * balls[shift_n + index] for index in range(degree_d + 1)]
    values = [flint.arb(1)]
    for size_m in range(1, max_size_m + 1):
        value = flint.arb(0)
        for offset in range(1, min(degree_d, size_m) + 1):
            term = (h[0] ** (offset - 1)) * h[offset] * values[size_m - offset]
            value = value + term if offset % 2 == 1 else value - term
        values.append(value)
    return values


def probe_rows(
    balls: dict[int, flint.arb],
    lam: str,
    shift_n: int,
    degree_d: int,
    sizes: tuple[int, ...],
) -> list[RecurrenceStressRow]:
    values = recurrence_values(balls, shift_n, degree_d, max(sizes))
    rows = []
    for size_m in sizes:
        value = values[size_m]
        classification = classify(value)
        contains_zero = value.contains(0)
        rows.append(
            RecurrenceStressRow(
                kind="arb_jensen_window_column_recurrence_stress_row",
                lam=lam,
                shift_n=shift_n,
                degree_d=degree_d,
                recurrence_size_m=size_m,
                max_coefficient_index=shift_n + degree_d,
                classification=classification,
                contains_zero=contains_zero,
                ok=(classification == "positive" and not contains_zero),
            )
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=[Path(p) for p in DEFAULT_ENCLOSURE_JSONL])
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--shifts", default="0..20")
    parser.add_argument("--degrees", default="3..8")
    parser.add_argument("--sizes", default="1..20")
    parser.add_argument("--dps", type=int, default=520)
    parser.add_argument(
        "--out-jsonl",
        type=Path,
        default=Path("work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520.jsonl"),
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=Path("work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520_summary.json"),
    )
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = args.dps
    lambdas = [part.strip() for part in args.lambdas.split(",") if part.strip()]
    shifts = parse_int_range(args.shifts)
    degrees = parse_int_range(args.degrees)
    sizes = tuple(parse_int_range(args.sizes))
    needed_max_k = max(shifts) + max(degrees)

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)

    rows: list[RecurrenceStressRow] = []
    summaries: list[RecurrenceStressSubsummary] = []
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                probe = probe_rows(balls, lam, shift_n, degree_d, sizes)
                rows.extend(probe)
                counts: Counter[str] = Counter(row.classification for row in probe)
                first_bad_size = next((row.recurrence_size_m for row in probe if not row.ok), None)
                positive = counts.get("positive", 0)
                negative = counts.get("negative", 0)
                zero = counts.get("zero", 0)
                inconclusive = sum(count for key, count in counts.items() if key not in {"positive", "negative", "zero"})
                summaries.append(
                    RecurrenceStressSubsummary(
                        lam=lam,
                        shift_n=shift_n,
                        degree_d=degree_d,
                        sizes=sizes,
                        tests=len(sizes),
                        positive=positive,
                        negative=negative,
                        zero=zero,
                        inconclusive=inconclusive,
                        first_bad_size=first_bad_size,
                        all_ok=(positive == len(sizes) and negative == 0 and zero == 0 and inconclusive == 0),
                    )
                )

    with args.out_jsonl.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    failed = [summary for summary in summaries if not summary.all_ok]
    degree_counts: dict[str, int] = {}
    for row in rows:
        if row.ok:
            degree_counts[str(row.degree_d)] = degree_counts.get(str(row.degree_d), 0) + 1
    payload = {
        "kind": "arb_jensen_window_column_recurrence_stress_summary",
        "date": "2026-07-06",
        "proof_boundary": (
            "Finite Arb stress test for the necessary Jensen-window column "
            "recurrence only; not all degrees, not all shifts, not all lambda "
            "values, not all skew shapes, and not a proof of RH or Lambda <= 0."
        ),
        "enclosure_jsonl": [str(path) for path in args.enclosure_jsonl],
        "lambdas": lambdas,
        "shifts": shifts,
        "degrees": degrees,
        "sizes": list(sizes),
        "dps": args.dps,
        "needed_max_k": needed_max_k,
        "summary_rows": len(summaries),
        "rows": len(rows),
        "ok": sum(1 for row in rows if row.ok),
        "failed_or_inconclusive": sum(1 for row in rows if not row.ok),
        "all_ok": not failed,
        "ok_by_degree": degree_counts,
        "summaries": [asdict(summary) for summary in summaries],
    }
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "validated "
            f"{payload['ok']}/{payload['rows']} Arb Jensen-window column recurrence stress rows"
        )
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
