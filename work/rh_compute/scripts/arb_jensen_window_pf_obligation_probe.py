#!/usr/bin/env python3
"""Arb finite diagnostic for low-degree Jensen-window PF obligations.

For each lambda, shift n, and degree d, form the binomially weighted Jensen
window B_j = binom(d,j) A_{n+j}.  This probe checks contiguous banded Toeplitz
minors det(T[rows 0..m-1, cols 1..m]) for small m.

The result is a finite interval diagnostic only.  It is not Jensen-window
PF-infinity, not all-degree Jensen hyperbolicity, and not a proof of RH or
Lambda <= 0.
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
class JensenWindowPFRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    minor_size_m: int
    max_coefficient_index: int
    classification: str
    contains_zero: bool
    ok: bool


@dataclass(frozen=True)
class JensenWindowPFSubsummary:
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


def parse_degree_sizes(text: str) -> dict[int, tuple[int, ...]]:
    out: dict[int, tuple[int, ...]] = {}
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        degree_text, sizes_text = part.split(":", 1)
        out[int(degree_text.strip())] = tuple(parse_int_range(sizes_text.strip()))
    return out


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


def window_value(balls: dict[int, flint.arb], shift_n: int, degree_d: int, index: int) -> flint.arb:
    if 0 <= index <= degree_d:
        return comb(degree_d, index) * balls[shift_n + index]
    return flint.arb(0)


def contiguous_toeplitz_det(
    balls: dict[int, flint.arb],
    shift_n: int,
    degree_d: int,
    minor_size_m: int,
) -> flint.arb:
    matrix = flint.arb_mat(
        [
            [window_value(balls, shift_n, degree_d, 1 + col - row) for col in range(minor_size_m)]
            for row in range(minor_size_m)
        ]
    )
    return matrix.det()


def probe_one(
    balls: dict[int, flint.arb],
    lam: str,
    shift_n: int,
    degree_d: int,
    minor_size_m: int,
) -> JensenWindowPFRow:
    detv = contiguous_toeplitz_det(balls, shift_n, degree_d, minor_size_m)
    classification = classify(detv)
    contains_zero = detv.contains(0)
    return JensenWindowPFRow(
        kind="arb_jensen_window_pf_obligation_row",
        lam=lam,
        shift_n=shift_n,
        degree_d=degree_d,
        minor_size_m=minor_size_m,
        max_coefficient_index=shift_n + degree_d,
        classification=classification,
        contains_zero=contains_zero,
        ok=(classification == "positive" and not contains_zero),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=[Path(p) for p in DEFAULT_ENCLOSURE_JSONL])
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--shifts", default="0..20")
    parser.add_argument("--degree-sizes", default="3:1..8;4:1..6")
    parser.add_argument("--dps", type=int, default=520)
    parser.add_argument(
        "--out-jsonl",
        type=Path,
        default=Path("work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520.jsonl"),
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=Path("work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520_summary.json"),
    )
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = args.dps
    lambdas = [part.strip() for part in args.lambdas.split(",") if part.strip()]
    shifts = parse_int_range(args.shifts)
    degree_sizes = parse_degree_sizes(args.degree_sizes)
    needed_max_k = max(shifts) + max(degree_sizes)

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)

    rows: list[JensenWindowPFRow] = []
    summaries: list[JensenWindowPFSubsummary] = []
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d, sizes in sorted(degree_sizes.items()):
                counts: Counter[str] = Counter()
                first_bad_size: int | None = None
                for minor_size_m in sizes:
                    row = probe_one(balls, lam, shift_n, degree_d, minor_size_m)
                    rows.append(row)
                    counts[row.classification] += 1
                    if not row.ok and first_bad_size is None:
                        first_bad_size = minor_size_m
                positive = counts.get("positive", 0)
                negative = counts.get("negative", 0)
                zero = counts.get("zero", 0)
                inconclusive = sum(
                    count for key, count in counts.items() if key not in {"positive", "negative", "zero"}
                )
                summaries.append(
                    JensenWindowPFSubsummary(
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
    payload = {
        "kind": "arb_jensen_window_pf_obligation_summary",
        "date": "2026-07-05",
        "proof_boundary": (
            "Finite Arb coefficient-enclosure diagnostic for selected "
            "Jensen-window PF Toeplitz minors only; not all-minor "
            "Jensen-window PF-infinity, not Jensen hyperbolicity, and not a "
            "proof of RH or Lambda <= 0."
        ),
        "enclosure_jsonl": [str(path) for path in args.enclosure_jsonl],
        "lambdas": lambdas,
        "shifts": shifts,
        "degree_sizes": {str(degree): list(sizes) for degree, sizes in sorted(degree_sizes.items())},
        "dps": args.dps,
        "needed_max_k": needed_max_k,
        "summary_rows": len(summaries),
        "rows": len(rows),
        "ok": sum(1 for row in rows if row.ok),
        "failed_or_inconclusive": sum(1 for row in rows if not row.ok),
        "all_ok": not failed,
        "summaries": [asdict(summary) for summary in summaries],
    }
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        for summary in summaries:
            status = "OK" if summary.all_ok else "FAIL"
            print(
                f"{status} Jensen-window PF obligation: lambda={summary.lam}, "
                f"n={summary.shift_n}, d={summary.degree_d}, "
                f"{summary.positive}/{summary.tests} positive"
            )
        print(
            "validated "
            f"{sum(1 for summary in summaries if summary.all_ok)}/{len(summaries)} "
            "Arb Jensen-window PF obligation rows"
        )
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
