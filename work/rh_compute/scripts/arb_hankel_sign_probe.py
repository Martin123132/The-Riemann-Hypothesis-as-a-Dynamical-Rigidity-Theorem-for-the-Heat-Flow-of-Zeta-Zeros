#!/usr/bin/env python
"""
FLINT/Arb sign-separation probe for cached signed Hankel determinants.

This is not a certificate for the exact transcendental moments unless the
input coefficient balls are known to enclose the exact A_k(lambda). It is a
propagation gate: given cached decimal coefficients and an assumed absolute
coefficient radius, determine whether the signed determinant ball stays away
from zero.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


@dataclass
class ProbeRow:
    lam: str
    m: int
    shift: int
    sigma: int
    abs_radius: str
    det: str
    signed_det: str
    contains_zero: bool
    sign_separated: bool
    expected_positive: bool
    ok: bool


def sigma(m: int) -> int:
    return -1 if ((m * (m + 1) // 2) % 2) else 1


def load_coefficients(path: Path, lam: str) -> list[str]:
    target = flint.arb(lam)
    best: list[str] | None = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") != "coefficients":
                continue
            if flint.arb(row["lam"]) == target:
                values = list(row["A"])
                if best is None or len(values) > len(best):
                    best = values
    if best is None:
        raise RuntimeError(f"No coefficient row found for lambda={lam} in {path}")
    return best


def ball(value: str, abs_radius: str) -> flint.arb:
    return flint.arb(f"[{value} +/- {abs_radius}]")


def probe_one(values: list[str], lam: str, m: int, shift: int, abs_radius: str) -> ProbeRow:
    needed = 2 * m + shift
    if needed >= len(values):
        raise RuntimeError(f"Need A up to {needed}, cache only has {len(values)-1}")
    mat = flint.arb_mat(
        [[ball(values[i + j + shift], abs_radius) for j in range(m + 1)] for i in range(m + 1)]
    )
    detv = mat.det()
    sig = sigma(m)
    signed = detv if sig > 0 else -detv
    contains_zero = signed.contains(0)
    sign_separated = not contains_zero
    expected_positive = signed > 0
    return ProbeRow(
        lam=lam,
        m=m,
        shift=shift,
        sigma=sig,
        abs_radius=abs_radius,
        det=str(detv),
        signed_det=str(signed),
        contains_zero=contains_zero,
        sign_separated=sign_separated,
        expected_positive=expected_positive,
        ok=sign_separated and expected_positive,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coeff-cache", type=Path, required=True)
    parser.add_argument("--lambda", dest="lam", default="0")
    parser.add_argument("--m-min", type=int, default=0)
    parser.add_argument("--m-max", type=int, default=4)
    parser.add_argument("--shift-min", type=int, default=0)
    parser.add_argument("--shift-max", type=int, default=0)
    parser.add_argument("--abs-radius", default="1e-80")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--out-jsonl", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = args.dps
    values = load_coefficients(args.coeff_cache, args.lam)
    rows: list[ProbeRow] = []
    for shift in range(args.shift_min, args.shift_max + 1):
        for m in range(args.m_min, args.m_max + 1):
            rows.append(probe_one(values, args.lam, m, shift, args.abs_radius))

    if args.out_jsonl is not None:
        args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
        with args.out_jsonl.open("w", encoding="utf-8") as f:
            for row in rows:
                f.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    failed = [row for row in rows if not row.ok]
    print(f"rows={len(rows)} ok={len(rows)-len(failed)} failed_or_inconclusive={len(failed)}")
    for row in failed[:10]:
        print(json.dumps(asdict(row), sort_keys=True))
    return 0 if not failed else 1


if __name__ == "__main__":
    raise SystemExit(main())
