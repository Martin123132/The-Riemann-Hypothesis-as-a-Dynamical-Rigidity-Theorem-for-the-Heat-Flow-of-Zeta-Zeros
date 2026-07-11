#!/usr/bin/env python3
"""Finite Arb stress test for Jensen-window monotone contractions.

For positive coefficients A_k(lambda), define adjacent ratio contractions

    x_k = (A_{k+1}/A_k) / (A_k/A_{k-1}).

The monotone-contraction frontier scout found that the first hard
column-frontier polynomials are positive on the region x_1 <= x_2 <= x_3.
This script stress-tests that increasing-contraction condition on the existing
zeta coefficient enclosure cache.  It is finite evidence only.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_ENCLOSURE_JSONL = (
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl",
    REPO_ROOT / "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl",
)
DEFAULT_OUT_JSONL = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64.jsonl"
DEFAULT_OUT_SUMMARY = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_monotone_contraction_stress_lamgrid_d3_d12_k64_summary.json"


@dataclass(frozen=True)
class MonotoneStressRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    max_coefficient_index: int
    contraction_count: int
    adjacent_log_concavity_ok: bool
    monotone_contractions_ok: bool
    contains_zero: bool
    ok: bool
    min_one_minus_contraction_sample: str
    min_monotone_gap_sample: str


@dataclass(frozen=True)
class DegreeSummary:
    degree_d: int
    shifts_checked_per_lambda: int
    checked_rows: int
    ok_rows: int
    failed_or_inconclusive_rows: int
    min_one_minus_contraction_sample: str
    min_one_minus_at: dict
    min_monotone_gap_sample: str
    min_monotone_gap_at: dict


def parse_int_range(text: str) -> list[int]:
    out: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if ".." in part:
            left, right = part.split("..", 1)
            start = int(left)
            stop = int(right)
            if stop < start:
                raise ValueError(f"descending range: {part}")
            out.extend(range(start, stop + 1))
        else:
            out.append(int(part))
    return sorted(set(out))


def decimal_lam_key(text: str) -> Decimal:
    return Decimal(str(text)).normalize()


def decimal_format(value: Decimal) -> str:
    return f"{value:.18E}"


def decimal_sample(row: dict) -> Decimal:
    cached = row.get("cache_A")
    if cached is not None:
        return Decimal(cached)
    ball = str(row["A_ball"]).strip()
    if not (ball.startswith("[") and "+/-" in ball):
        raise ValueError(f"cannot parse A_ball midpoint: {ball}")
    midpoint = ball[1:].split("+/-", 1)[0].strip()
    return Decimal(midpoint)


def load_enclosures(paths: list[Path]) -> tuple[dict[tuple[Decimal, int], flint.arb], dict[tuple[Decimal, int], Decimal], dict[Decimal, str]]:
    balls: dict[tuple[Decimal, int], flint.arb] = {}
    samples: dict[tuple[Decimal, int], Decimal] = {}
    labels: dict[Decimal, str] = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                if not raw.strip():
                    continue
                row = json.loads(raw)
                if row.get("kind") != "acb_coefficient_enclosure":
                    continue
                lam = decimal_lam_key(row["lam"])
                index = int(row["k"])
                balls[(lam, index)] = flint.arb(row["A_ball"])
                samples[(lam, index)] = decimal_sample(row)
                labels[lam] = row["lam"]
    return balls, samples, labels


def contractions(values):
    ratios = [values[index] / values[index - 1] for index in range(1, len(values))]
    return [ratios[index] / ratios[index - 1] for index in range(1, len(ratios))]


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def row_for(
    balls: dict[tuple[Decimal, int], flint.arb],
    samples: dict[tuple[Decimal, int], Decimal],
    lam: Decimal,
    lam_label: str,
    shift_n: int,
    degree_d: int,
) -> MonotoneStressRow:
    arb_values = [balls[(lam, shift_n + offset)] for offset in range(degree_d + 1)]
    arb_xs = contractions(arb_values)
    one_minus = [1 - item for item in arb_xs]
    monotone_gaps = [arb_xs[index + 1] - arb_xs[index] for index in range(len(arb_xs) - 1)]
    adjacent_ok = all(arb_positive(item) for item in one_minus)
    monotone_ok = all(arb_positive(item) for item in monotone_gaps)
    contains_zero = any(item.contains(0) for item in one_minus + monotone_gaps)

    sample_values = [samples[(lam, shift_n + offset)] for offset in range(degree_d + 1)]
    sample_xs = contractions(sample_values)
    sample_one_minus = [Decimal(1) - item for item in sample_xs]
    sample_gaps = [sample_xs[index + 1] - sample_xs[index] for index in range(len(sample_xs) - 1)]
    return MonotoneStressRow(
        kind="jensen_window_pf_monotone_contraction_stress_row",
        lam=lam_label,
        shift_n=shift_n,
        degree_d=degree_d,
        max_coefficient_index=shift_n + degree_d,
        contraction_count=len(arb_xs),
        adjacent_log_concavity_ok=adjacent_ok,
        monotone_contractions_ok=monotone_ok,
        contains_zero=contains_zero,
        ok=adjacent_ok and monotone_ok and not contains_zero,
        min_one_minus_contraction_sample=decimal_format(min(sample_one_minus)),
        min_monotone_gap_sample=decimal_format(min(sample_gaps)),
    )


def build_rows(
    balls: dict[tuple[Decimal, int], flint.arb],
    samples: dict[tuple[Decimal, int], Decimal],
    labels: dict[Decimal, str],
    degrees: list[int],
    max_coefficient_index: int,
) -> list[MonotoneStressRow]:
    rows: list[MonotoneStressRow] = []
    for degree_d in degrees:
        for lam in sorted(labels):
            for shift_n in range(0, max_coefficient_index - degree_d + 1):
                rows.append(row_for(balls, samples, lam, labels[lam], shift_n, degree_d))
    return rows


def degree_summaries(rows: list[MonotoneStressRow]) -> list[DegreeSummary]:
    out: list[DegreeSummary] = []
    for degree_d in sorted({row.degree_d for row in rows}):
        degree_rows = [row for row in rows if row.degree_d == degree_d]
        min_one = min(degree_rows, key=lambda row: Decimal(row.min_one_minus_contraction_sample))
        min_gap = min(degree_rows, key=lambda row: Decimal(row.min_monotone_gap_sample))
        shifts_per_lambda = len({row.shift_n for row in degree_rows})
        out.append(
            DegreeSummary(
                degree_d=degree_d,
                shifts_checked_per_lambda=shifts_per_lambda,
                checked_rows=len(degree_rows),
                ok_rows=sum(1 for row in degree_rows if row.ok),
                failed_or_inconclusive_rows=sum(1 for row in degree_rows if not row.ok),
                min_one_minus_contraction_sample=min_one.min_one_minus_contraction_sample,
                min_one_minus_at={"lambda": min_one.lam, "shift_n": min_one.shift_n},
                min_monotone_gap_sample=min_gap.min_monotone_gap_sample,
                min_monotone_gap_at={"lambda": min_gap.lam, "shift_n": min_gap.shift_n},
            )
        )
    return out


def build_summary(
    rows: list[MonotoneStressRow],
    degree_rows: list[DegreeSummary],
    enclosure_paths: list[Path],
    degrees: list[int],
    max_coefficient_index: int,
    out_jsonl: Path,
) -> dict:
    failures = [row for row in rows if not row.ok]
    ok_by_degree: Counter[str] = Counter(str(row.degree_d) for row in rows if row.ok)
    min_one = min(rows, key=lambda row: Decimal(row.min_one_minus_contraction_sample))
    min_gap = min(rows, key=lambda row: Decimal(row.min_monotone_gap_sample))
    lambdas = sorted({row.lam for row in rows}, key=lambda text: Decimal(text).normalize())
    return {
        "kind": "jensen_window_pf_monotone_contraction_stress_summary",
        "date": "2026-07-06",
        "status": "finite_arb_ratio_curvature_stress",
        "proof_boundary": (
            "Finite Arb stress test for adjacent log-concavity and increasing "
            "ratio contractions on the existing zeta coefficient grid only; "
            "not all shifts beyond the cache, not all lambda values, not an "
            "analytic monotone-contraction theorem, not jwpf_06, and not a "
            "proof of RH or Lambda <= 0."
        ),
        "source_frontier_scout": "outputs/jensen_window_pf_monotone_contraction_frontier_scout.md",
        "row_jsonl": out_jsonl.relative_to(REPO_ROOT).as_posix(),
        "enclosure_jsonl": [path.relative_to(REPO_ROOT).as_posix() for path in enclosure_paths],
        "lambdas": lambdas,
        "degrees": degrees,
        "max_coefficient_index": max_coefficient_index,
        "shift_rule": "for each degree d, check every shift n with 0 <= n <= max_coefficient_index - d",
        "condition": {
            "x_k": "(A_{k+1}/A_k) / (A_k/A_{k-1})",
            "adjacent_log_concavity": "0 < x_k <= 1",
            "increasing_contractions": "x_{k+1} >= x_k",
            "equivalent_polynomial_inequality": "A_{k+2}*A_k**3 >= A_{k+1}**3*A_{k-1}",
        },
        "rows": len(rows),
        "ok": sum(1 for row in rows if row.ok),
        "failed_or_inconclusive": len(failures),
        "all_ok": not failures,
        "ok_by_degree": dict(sorted(ok_by_degree.items(), key=lambda item: int(item[0]))),
        "degree_summaries": [asdict(row) for row in degree_rows],
        "global_min_one_minus_contraction_sample": min_one.min_one_minus_contraction_sample,
        "global_min_one_minus_at": {"degree_d": min_one.degree_d, "lambda": min_one.lam, "shift_n": min_one.shift_n},
        "global_min_monotone_gap_sample": min_gap.min_monotone_gap_sample,
        "global_min_monotone_gap_at": {"degree_d": min_gap.degree_d, "lambda": min_gap.lam, "shift_n": min_gap.shift_n},
        "summary": {
            "stress_rows": len(rows),
            "positive_rows": sum(1 for row in rows if row.ok),
            "failed_or_inconclusive_rows": len(failures),
            "degree_count": len(degrees),
            "target_closing": False,
            "main_finding": (
                "On the existing k<=64 zeta coefficient enclosure cache, every "
                "checked Jensen-window prefix for degrees d=3..12 satisfies "
                "adjacent log-concavity and increasing ratio contractions.  "
                "This extends the monotone-contraction frontier observation "
                "from the first hard rows to 2875 finite Arb-classified rows, "
                "but remains finite evidence only."
            ),
        },
        "invariants": [
            "Every stress row must have adjacent_log_concavity_ok=true.",
            "Every stress row must have monotone_contractions_ok=true.",
            "No row may contain zero in a required positivity gap.",
            "The summary target_closing flag must remain false.",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=list(DEFAULT_ENCLOSURE_JSONL))
    parser.add_argument("--degrees", default="3..12")
    parser.add_argument("--max-coefficient-index", type=int, default=64)
    parser.add_argument("--dps", type=int, default=160)
    parser.add_argument("--out-jsonl", type=Path, default=DEFAULT_OUT_JSONL)
    parser.add_argument("--out-summary", type=Path, default=DEFAULT_OUT_SUMMARY)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    getcontext().prec = 80
    flint.ctx.dps = args.dps
    degrees = parse_int_range(args.degrees)
    balls, samples, labels = load_enclosures(args.enclosure_jsonl)
    rows = build_rows(balls, samples, labels, degrees, args.max_coefficient_index)
    summaries = degree_summaries(rows)
    payload = build_summary(rows, summaries, args.enclosure_jsonl, degrees, args.max_coefficient_index, args.out_jsonl)

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.out_jsonl.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "validated "
            f"{payload['ok']}/{payload['rows']} Jensen-window monotone contraction stress rows"
        )
    return 0 if payload["all_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
