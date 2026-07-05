#!/usr/bin/env python3
"""Arb sign probe for signed Hankel determinants using A_k enclosure rows.

Unlike arb_hankel_sign_probe.py, this script does not wrap truncated cache
decimals in a uniform radius.  It reads the rigorous A_ball fields produced by
acb_coefficient_enclosures.py and propagates those coefficient balls directly
through the signed Hankel determinants.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal
from pathlib import Path


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


@dataclass(frozen=True)
class ProbeRow:
    kind: str
    lam: str
    m: int
    shift: int
    sigma: int
    det: str
    signed_det: str
    classification: str
    contains_zero: bool
    sign_separated: bool
    expected_positive: bool
    ok: bool


def sigma(m: int) -> int:
    return -1 if ((m * (m + 1) // 2) % 2) else 1


def decimal_equal(left: str, right: str) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


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


def probe_one(balls: dict[int, flint.arb], lam: str, m: int, shift: int) -> ProbeRow:
    mat = flint.arb_mat(
        [[balls[i + j + shift] for j in range(m + 1)] for i in range(m + 1)]
    )
    detv = mat.det()
    sig = sigma(m)
    signed = detv if sig > 0 else -detv
    classification = classify(signed)
    contains_zero = signed.contains(0)
    sign_separated = not contains_zero
    expected_positive = signed > 0
    return ProbeRow(
        kind="arb_hankel_enclosure_sign_probe",
        lam=lam,
        m=m,
        shift=shift,
        sigma=sig,
        det=str(detv),
        signed_det=str(signed),
        classification=classification,
        contains_zero=contains_zero,
        sign_separated=sign_separated,
        expected_positive=expected_positive,
        ok=sign_separated and expected_positive,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", required=True)
    parser.add_argument("--lambda", dest="lam", default="0")
    parser.add_argument("--m-min", type=int, default=0)
    parser.add_argument("--m-max", type=int, default=12)
    parser.add_argument("--shift-min", type=int, default=0)
    parser.add_argument("--shift-max", type=int, default=8)
    parser.add_argument("--dps", type=int, default=180)
    parser.add_argument("--out-jsonl", type=Path, default=None)
    parser.add_argument("--out-summary", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.m_min < 0 or args.m_max < args.m_min:
        raise ValueError("invalid m range")
    if args.shift_min < 0 or args.shift_max < args.shift_min:
        raise ValueError("invalid shift range")

    flint.ctx.dps = args.dps
    needed_max_k = 2 * args.m_max + args.shift_max
    balls = load_enclosure_balls(args.enclosure_jsonl, args.lam, needed_max_k)

    rows: list[ProbeRow] = []
    counts: Counter[str] = Counter()
    for shift in range(args.shift_min, args.shift_max + 1):
        for m in range(args.m_min, args.m_max + 1):
            row = probe_one(balls, args.lam, m, shift)
            rows.append(row)
            counts[row.classification] += 1

    if args.out_jsonl is not None:
        args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
        with args.out_jsonl.open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    failed = [row for row in rows if not row.ok]
    summary = {
        "kind": "arb_hankel_enclosure_sign_probe_summary",
        "lam": args.lam,
        "m_min": args.m_min,
        "m_max": args.m_max,
        "shift_min": args.shift_min,
        "shift_max": args.shift_max,
        "dps": args.dps,
        "needed_max_k": needed_max_k,
        "rows": len(rows),
        "ok": len(rows) - len(failed),
        "failed_or_inconclusive": len(failed),
        "counts": dict(counts),
        "all_ok": not failed,
    }
    if args.out_summary is not None:
        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        args.out_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(json.dumps(summary, sort_keys=True))
    for row in failed[:10]:
        print(json.dumps(asdict(row), sort_keys=True))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
