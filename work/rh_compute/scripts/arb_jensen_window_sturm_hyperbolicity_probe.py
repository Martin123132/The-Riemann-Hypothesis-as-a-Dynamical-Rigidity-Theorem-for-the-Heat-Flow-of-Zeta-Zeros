#!/usr/bin/env python3
"""Arb/Sturm finite diagnostic for low-degree Jensen-window hyperbolicity.

For each lambda, shift n, and degree d, form:

    Q_{d,n,lambda}(y) = P_{d,n,lambda}(-y)
                      = sum_j (-1)^j binom(d,j) A_{n+j}(lambda) y^j.

The selected diagnostic checks that an interval-enclosed Sturm sequence has
exactly d sign-variation drops on (0, infinity), for d=3,4 in the configured
finite grid.  This is a finite diagnostic only; it is not all-degree/all-shift
Jensen hyperbolicity and not a proof of RH or Lambda <= 0.
"""

from __future__ import annotations

import argparse
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
)


@dataclass(frozen=True)
class SturmRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    expected_positive_roots: int
    positive_root_count: int | None
    variations_at_zero: int | None
    variations_at_infinity: int | None
    signs_at_zero: tuple[int | None, ...]
    signs_at_infinity: tuple[int | None, ...]
    sturm_degrees: tuple[int, ...]
    classification: str
    ok: bool


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


def trim(poly: list[flint.arb]) -> list[flint.arb]:
    out = list(poly)
    while len(out) > 1 and out[-1] == flint.arb(0):
        out.pop()
    return out


def degree(poly: list[flint.arb]) -> int:
    return len(trim(poly)) - 1


def derivative(poly: list[flint.arb]) -> list[flint.arb]:
    if len(poly) <= 1:
        return [flint.arb(0)]
    return [i * poly[i] for i in range(1, len(poly))]


def negative(poly: list[flint.arb]) -> list[flint.arb]:
    return [-coeff for coeff in poly]


def remainder(numerator: list[flint.arb], denominator: list[flint.arb]) -> list[flint.arb]:
    """Interval polynomial remainder with algebraic leading cancellation.

    The leading term in each division step is exactly cancelled in exact
    arithmetic.  We pop it explicitly to avoid interval dependency leaving a
    small zero-containing leading ball and stalling the degree reduction.
    """
    rem = trim(numerator)
    den = trim(denominator)
    den_degree = degree(den)
    den_lead = den[-1]
    while degree(rem) >= den_degree and not (len(rem) == 1 and rem[0] == flint.arb(0)):
        rem_degree = degree(rem)
        factor = rem[-1] / den_lead
        shift = rem_degree - den_degree
        for i, coeff in enumerate(den):
            rem[shift + i] = rem[shift + i] - factor * coeff
        rem.pop()
        rem = trim(rem)
    return rem


def sturm_sequence(poly: list[flint.arb]) -> list[list[flint.arb]]:
    seq = [trim(poly), trim(derivative(poly))]
    while degree(seq[-1]) > 0:
        seq.append(trim(negative(remainder(seq[-2], seq[-1]))))
        if len(seq) > 8:
            raise RuntimeError("unexpectedly long Sturm sequence")
    return seq


def sign_arb(value: flint.arb) -> int | None:
    if value > 0:
        return 1
    if value < 0:
        return -1
    if value == flint.arb(0):
        return 0
    return None


def variations(signs: list[int | None]) -> int | None:
    if any(sign is None for sign in signs):
        return None
    nonzero = [sign for sign in signs if sign != 0]
    return sum(1 for left, right in zip(nonzero, nonzero[1:]) if left * right < 0)


def q_coefficients(balls: dict[int, flint.arb], shift_n: int, degree_d: int) -> list[flint.arb]:
    return [((-1) ** j) * comb(degree_d, j) * balls[shift_n + j] for j in range(degree_d + 1)]


def probe_one(balls: dict[int, flint.arb], lam: str, shift_n: int, degree_d: int) -> SturmRow:
    seq = sturm_sequence(q_coefficients(balls, shift_n, degree_d))
    signs_at_zero = [sign_arb(poly[0]) for poly in seq]
    signs_at_infinity = [sign_arb(poly[-1]) for poly in seq]
    v0 = variations(signs_at_zero)
    vinf = variations(signs_at_infinity)
    count = None if v0 is None or vinf is None else v0 - vinf
    ok = count == degree_d
    classification = "positive_root_count_certified" if ok else "inconclusive_or_failed"
    return SturmRow(
        kind="arb_jensen_window_sturm_hyperbolicity_row",
        lam=lam,
        shift_n=shift_n,
        degree_d=degree_d,
        expected_positive_roots=degree_d,
        positive_root_count=count,
        variations_at_zero=v0,
        variations_at_infinity=vinf,
        signs_at_zero=tuple(signs_at_zero),
        signs_at_infinity=tuple(signs_at_infinity),
        sturm_degrees=tuple(degree(poly) for poly in seq),
        classification=classification,
        ok=ok,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=[Path(p) for p in DEFAULT_ENCLOSURE_JSONL])
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--shifts", default="0..20")
    parser.add_argument("--degrees", default="3,4")
    parser.add_argument("--dps", type=int, default=520)
    parser.add_argument(
        "--out-jsonl",
        type=Path,
        default=Path("work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520.jsonl"),
    )
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=Path("work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json"),
    )
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = args.dps
    lambdas = [part.strip() for part in args.lambdas.split(",") if part.strip()]
    shifts = parse_int_range(args.shifts)
    degrees = parse_int_range(args.degrees)
    needed_max_k = max(shifts) + max(degrees)
    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)

    rows: list[SturmRow] = []
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                rows.append(probe_one(balls, lam, shift_n, degree_d))

    with args.out_jsonl.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    ok_rows = [row for row in rows if row.ok]
    failed = [row for row in rows if not row.ok]
    payload = {
        "kind": "arb_jensen_window_sturm_hyperbolicity_summary",
        "date": "2026-07-05",
        "proof_boundary": (
            "Finite Arb/Sturm-style diagnostic for degree-3 and degree-4 "
            "Jensen windows only; not all-degree or all-shift Jensen "
            "hyperbolicity, not Jensen-window PF-infinity, and not a proof of "
            "RH or Lambda <= 0."
        ),
        "enclosure_jsonl": [str(path) for path in args.enclosure_jsonl],
        "lambdas": lambdas,
        "shifts": shifts,
        "degrees": degrees,
        "dps": args.dps,
        "needed_max_k": needed_max_k,
        "rows": len(rows),
        "ok": len(ok_rows),
        "failed_or_inconclusive": len(failed),
        "all_ok": not failed,
        "rows_by_degree": {str(degree): sum(1 for row in rows if row.degree_d == degree) for degree in degrees},
        "ok_by_degree": {str(degree): sum(1 for row in ok_rows if row.degree_d == degree) for degree in degrees},
    }
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        for row in rows:
            status = "OK" if row.ok else "FAIL"
            print(
                f"{status} Jensen-window Sturm: lambda={row.lam}, n={row.shift_n}, "
                f"d={row.degree_d}, roots={row.positive_root_count}/{row.expected_positive_roots}"
            )
        print(f"validated {len(ok_rows)}/{len(rows)} Arb Jensen-window Sturm hyperbolicity rows")
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
