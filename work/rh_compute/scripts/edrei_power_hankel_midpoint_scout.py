#!/usr/bin/env python3
"""Midpoint-only scout for Edrei power-Hankel frontier determinants.

This script deliberately discards Arb radii from the Edrei-log rows and
computes determinants from the interval midpoints only.  It is useful for
frontier triage when the rigorous determinant ball contains zero, but it is
not a certificate and must not be used as PF-infinity evidence.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


@dataclass(frozen=True)
class MidpointScoutRow:
    kind: str
    lam: str
    m: int
    shift: int
    det_midpoint_scout: str
    midpoint_classification: str
    contains_zero: bool
    rigorous: bool


def decimal_equal(left: str, right: str) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


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


def load_midpoint_powers(paths: list[Path], lam: str, needed_max_n: int) -> dict[int, flint.arb]:
    powers: dict[int, flint.arb] = {}
    for path in paths:
        with path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                row = json.loads(line)
                if row.get("kind") != "arb_edrei_log_sign_probe":
                    continue
                if not decimal_equal(row["lam"], lam):
                    continue
                n = int(row["n"])
                if 1 <= n <= needed_max_n:
                    powers[n] = flint.arb(row["signed_q_n"]).mid()
    missing = [n for n in range(1, needed_max_n + 1) if n not in powers]
    if missing:
        raise RuntimeError(f"missing Edrei-log rows for lambda={lam}: {missing[:10]}")
    return powers


def probe_one(powers: dict[int, flint.arb], lam: str, m: int, shift: int) -> MidpointScoutRow:
    mat = flint.arb_mat(
        [[powers[i + j + shift] for j in range(m + 1)] for i in range(m + 1)]
    )
    detv = mat.det()
    return MidpointScoutRow(
        kind="edrei_power_hankel_midpoint_scout",
        lam=lam,
        m=m,
        shift=shift,
        det_midpoint_scout=str(detv),
        midpoint_classification=classify(detv),
        contains_zero=detv.contains(0),
        rigorous=False,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edrei-jsonl", type=Path, nargs="+", required=True)
    parser.add_argument("--lambda", dest="lam", required=True)
    parser.add_argument("--m-min", type=int, required=True)
    parser.add_argument("--m-max", type=int, required=True)
    parser.add_argument("--shift-min", type=int, default=1)
    parser.add_argument("--shift-max", type=int, default=1)
    parser.add_argument("--dps", type=int, default=2000)
    parser.add_argument("--out-jsonl", type=Path, default=None)
    parser.add_argument("--out-summary", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.m_min < 0 or args.m_max < args.m_min:
        raise ValueError("invalid m range")
    if args.shift_min < 1 or args.shift_max < args.shift_min:
        raise ValueError("shift range must start at 1 or higher")
    flint.ctx.dps = args.dps
    needed_max_n = 2 * args.m_max + args.shift_max
    powers = load_midpoint_powers(args.edrei_jsonl, args.lam, needed_max_n)

    rows: list[MidpointScoutRow] = []
    for m in range(args.m_min, args.m_max + 1):
        for shift in range(args.shift_min, args.shift_max + 1):
            rows.append(probe_one(powers, args.lam, m, shift))

    counts = Counter(row.midpoint_classification for row in rows)
    summary = {
        "kind": "edrei_power_hankel_midpoint_scout_summary",
        "rigorous": False,
        "proof_boundary": (
            "midpoint-only warning scout for rigorous frontier determinants; "
            "not a certificate and not PF-infinity evidence"
        ),
        "lam": args.lam,
        "m_min": args.m_min,
        "m_max": args.m_max,
        "shift_min": args.shift_min,
        "shift_max": args.shift_max,
        "needed_max_n": needed_max_n,
        "dps": args.dps,
        "rows": len(rows),
        "counts": dict(counts),
    }

    if args.out_jsonl is not None:
        args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
        with args.out_jsonl.open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")
    if args.out_summary is not None:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        with args.out_summary.open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")

    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
