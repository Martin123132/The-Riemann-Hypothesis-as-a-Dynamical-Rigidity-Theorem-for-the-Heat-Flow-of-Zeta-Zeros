#!/usr/bin/env python3
"""Arb probe for the Edrei/entire-LP logarithmic sign diagnostic.

For the normalized ordinary generating function

    H_lambda(z) = sum_k d_k(lambda) z^k,   d_k = c_k / c_0,

an entire PF-infinity / restricted Laguerre-Polya target with only
nonpositive real zeros would have

    log H(z) = sum_{n>=1} ell_n z^n

with q_n = n ell_n satisfying

    (-1)^(n-1) q_n >= 0.

This script reads rigorous c_k enclosure balls from
acb_coefficient_enclosures.py output and checks that finite necessary
condition with Arb interval arithmetic.  It is not a proof of PF-infinity:
passing finite log-sign checks is only a bridge-theorem diagnostic.
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
class LogSignRow:
    kind: str
    lam: str
    n: int
    log_coeff: str
    q_n: str
    signed_q_n: str
    expected_sign: int
    classification: str
    contains_zero: bool
    ok: bool


def decimal_equal(left: str, right: str) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


def load_c_balls(paths: list[Path], lam: str, max_k: int) -> dict[int, flint.arb]:
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
                if k <= max_k:
                    balls[k] = flint.arb(row["c_ball"])
    missing = [k for k in range(max_k + 1) if k not in balls]
    if missing:
        raise RuntimeError(f"missing c_k enclosure balls for lambda={lam}: {missing[:10]}")
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


def log_coefficients(d: list[flint.arb]) -> list[flint.arb]:
    """Return ell_n for log(sum d_k z^k), with ell_0 unused/zero.

    The recurrence follows from H'(z) = H(z) (log H)'(z):

        n d_n = sum_{k=1}^n k ell_k d_{n-k}.
    """
    ell = [flint.arb(0) for _ in d]
    for n in range(1, len(d)):
        acc = flint.arb(n) * d[n]
        for k in range(1, n):
            acc -= flint.arb(k) * ell[k] * d[n - k]
        ell[n] = acc / flint.arb(n)
    return ell


def run_probe(c_balls: dict[int, flint.arb], lam: str, max_n: int) -> list[LogSignRow]:
    c0 = c_balls[0]
    if c0.contains(0):
        raise RuntimeError(f"c_0 enclosure contains zero for lambda={lam}: {c0}")
    d = [c_balls[k] / c0 for k in range(max_n + 1)]
    ell = log_coefficients(d)

    rows: list[LogSignRow] = []
    for n in range(1, max_n + 1):
        q_n = flint.arb(n) * ell[n]
        expected = 1 if n % 2 else -1
        signed_q = q_n if expected > 0 else -q_n
        classification = classify(signed_q)
        contains_zero = signed_q.contains(0)
        ok = signed_q > 0 and not contains_zero
        rows.append(
            LogSignRow(
                kind="arb_edrei_log_sign_probe",
                lam=lam,
                n=n,
                log_coeff=str(ell[n]),
                q_n=str(q_n),
                signed_q_n=str(signed_q),
                expected_sign=expected,
                classification=classification,
                contains_zero=contains_zero,
                ok=ok,
            )
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", required=True)
    parser.add_argument("--lambda", dest="lam", default="0")
    parser.add_argument("--max-n", type=int, default=32)
    parser.add_argument("--dps", type=int, default=220)
    parser.add_argument("--out-jsonl", type=Path, default=None)
    parser.add_argument("--out-summary", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.max_n < 1:
        raise ValueError("--max-n must be positive")
    flint.ctx.dps = args.dps

    c_balls = load_c_balls(args.enclosure_jsonl, args.lam, args.max_n)
    rows = run_probe(c_balls, args.lam, args.max_n)
    counts = Counter(row.classification for row in rows)
    failed = [row for row in rows if not row.ok]

    summary = {
        "kind": "arb_edrei_log_sign_probe_summary",
        "rigorous": True,
        "proof_boundary": (
            "finite necessary-condition diagnostic for the entire coefficient PF route; "
            "not a proof of PF-infinity or RH"
        ),
        "lam": args.lam,
        "max_n": args.max_n,
        "dps": args.dps,
        "rows": len(rows),
        "counts": dict(counts),
        "failed_or_inconclusive": len(failed),
        "all_ok": len(failed) == 0,
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
    return 0 if summary["all_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
