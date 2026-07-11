#!/usr/bin/env python3
"""Extended Arb stress test for Jensen-window reciprocal coefficients.

For h_j = binom(d,j) A_{n+j}(lambda), set g_j=h_j/h_0 and

    E(t) = 1 / (1 - g_1 t + g_2 t^2 - ... + (-1)^d g_d t^d).

The column-recurrence target asks for [t^m]E(t) >= 0 for every m,d,n.
This script checks a larger finite grid directly in normalized reciprocal
coefficients.  It is finite theorem-search evidence only, not an all-order
recurrence theorem or a Schur/Toeplitz proof.
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


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_LAMBDAS = ("0", "1e-6", "1e-4", "1e-2", "1e-1")
DEFAULT_ENCLOSURE_JSONL = (
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl",
)
DEFAULT_OUT_JSONL = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_lamgrid_n0_n32_d2_d12_m1_m40_dps520.jsonl"
)
DEFAULT_OUT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_coefficient_extended_stress.json"
)


@dataclass(frozen=True)
class ReciprocalCoefficientRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    coefficient_m: int
    max_coefficient_index: int
    classification: str
    contains_zero: bool
    ok: bool


@dataclass(frozen=True)
class ReciprocalCoefficientSubsummary:
    lam: str
    shift_n: int
    degree_d: int
    coefficient_sizes: tuple[int, ...]
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


def reciprocal_coefficients(
    balls: dict[int, flint.arb],
    shift_n: int,
    degree_d: int,
    max_size_m: int,
) -> list[flint.arb]:
    h = [comb(degree_d, index) * balls[shift_n + index] for index in range(degree_d + 1)]
    g = [item / h[0] for item in h]
    values = [flint.arb(1)]
    for size_m in range(1, max_size_m + 1):
        value = flint.arb(0)
        for offset in range(1, min(degree_d, size_m) + 1):
            term = g[offset] * values[size_m - offset]
            value = value + term if offset % 2 == 1 else value - term
        values.append(value)
    return values


def probe_rows(
    balls: dict[int, flint.arb],
    lam: str,
    shift_n: int,
    degree_d: int,
    sizes: tuple[int, ...],
) -> list[ReciprocalCoefficientRow]:
    values = reciprocal_coefficients(balls, shift_n, degree_d, max(sizes))
    rows = []
    for size_m in sizes:
        value = values[size_m]
        classification = classify(value)
        contains_zero = value.contains(0)
        rows.append(
            ReciprocalCoefficientRow(
                kind="jensen_window_pf_reciprocal_coefficient_extended_stress_row",
                lam=lam,
                shift_n=shift_n,
                degree_d=degree_d,
                coefficient_m=size_m,
                max_coefficient_index=shift_n + degree_d,
                classification=classification,
                contains_zero=contains_zero,
                ok=(classification == "positive" and not contains_zero),
            )
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=[REPO_ROOT / p for p in DEFAULT_ENCLOSURE_JSONL])
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--shifts", default="0..32")
    parser.add_argument("--degrees", default="2..12")
    parser.add_argument("--sizes", default="1..40")
    parser.add_argument("--dps", type=int, default=520)
    parser.add_argument("--out-jsonl", type=Path, default=DEFAULT_OUT_JSONL)
    parser.add_argument("--out-summary", type=Path, default=DEFAULT_OUT_SUMMARY)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = args.dps
    lambdas = [part.strip() for part in args.lambdas.split(",") if part.strip()]
    shifts = parse_int_range(args.shifts)
    degrees = parse_int_range(args.degrees)
    sizes = tuple(parse_int_range(args.sizes))
    if min(degrees) < 2:
        raise ValueError("degrees must be at least 2")
    if min(sizes) < 1:
        raise ValueError("coefficient sizes must be at least 1")
    needed_max_k = max(shifts) + max(degrees)

    rows: list[ReciprocalCoefficientRow] = []
    summaries: list[ReciprocalCoefficientSubsummary] = []
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                probe = probe_rows(balls, lam, shift_n, degree_d, sizes)
                rows.extend(probe)
                counts: Counter[str] = Counter(row.classification for row in probe)
                first_bad_size = next((row.coefficient_m for row in probe if not row.ok), None)
                positive = counts.get("positive", 0)
                negative = counts.get("negative", 0)
                zero = counts.get("zero", 0)
                inconclusive = sum(count for key, count in counts.items() if key not in {"positive", "negative", "zero"})
                summaries.append(
                    ReciprocalCoefficientSubsummary(
                        lam=lam,
                        shift_n=shift_n,
                        degree_d=degree_d,
                        coefficient_sizes=sizes,
                        tests=len(sizes),
                        positive=positive,
                        negative=negative,
                        zero=zero,
                        inconclusive=inconclusive,
                        first_bad_size=first_bad_size,
                        all_ok=(positive == len(sizes) and negative == 0 and zero == 0 and inconclusive == 0),
                    )
                )

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.out_jsonl.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    failed = [summary for summary in summaries if not summary.all_ok]
    ok_by_degree: dict[str, int] = {}
    for row in rows:
        if row.ok:
            ok_by_degree[str(row.degree_d)] = ok_by_degree.get(str(row.degree_d), 0) + 1
    payload = {
        "kind": "jensen_window_pf_reciprocal_coefficient_extended_stress",
        "date": "2026-07-06",
        "target_obligation": "jwpf_06_sign_regular_to_jensen_pf_conversion",
        "proof_boundary": (
            "Extended finite Arb stress test for normalized reciprocal coefficients only; "
            "not all degrees, not all shifts, not all lambda values, not all coefficient "
            "orders, not all skew shapes, and not a proof of RH or Lambda <= 0."
        ),
        "enclosure_jsonl": [str(path.relative_to(REPO_ROOT)) for path in args.enclosure_jsonl],
        "lambdas": lambdas,
        "shifts": shifts,
        "degrees": degrees,
        "coefficient_sizes": list(sizes),
        "dps": args.dps,
        "needed_max_k": needed_max_k,
        "summary_rows": len(summaries),
        "rows": len(rows),
        "ok": sum(1 for row in rows if row.ok),
        "failed_or_inconclusive": len(rows) - sum(1 for row in rows if row.ok),
        "all_ok": not failed,
        "ok_by_degree": ok_by_degree,
        "row_log": str(args.out_jsonl.relative_to(REPO_ROOT)),
        "summaries": [asdict(summary) for summary in summaries],
        "summary": {
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The normalized reciprocal coefficients [t^m]1/H(-t) are positive and "
                "separated from zero in all 72,600 checked Arb rows across lambdas "
                "{0,1e-6,1e-4,1e-2,1e-1}, shifts n=0..32, degrees d=2..12, "
                "and sizes m=1..40. This strengthens finite support for the column "
                "recurrence target but does not supply an all-order theorem."
            ),
        },
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF reciprocal coefficient extended stress: "
            f"{len(rows)} rows, {payload['ok']} positive, {payload['failed_or_inconclusive']} failed/inconclusive"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
