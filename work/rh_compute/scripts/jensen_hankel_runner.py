#!/usr/bin/env python
"""
Clean Jensen/Hankel runner for the RH dynamical-rigidity programme.

This is a reproducibility runner, not a rigorous interval certificate.
It computes high-precision mpmath moments, builds the A_k(lambda)
sequence, tests signed Hankel determinants, and can count positive roots
of Q_{d,n}(y)=P_{d,n}(-y) after rationalizing the computed coefficients.

Outputs JSONL rows so larger runs can be resumed, filtered, and audited.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from math import comb
from pathlib import Path
from time import perf_counter

import mpmath as mp
import sympy as sp


@dataclass
class HankelRow:
    kind: str
    lam: str
    m: int
    shift: int
    sigma: int
    det_sign: int
    signed_sign: int
    log10_abs_det: str
    near_zero: bool
    ok: bool
    dps: int
    n_sum: int
    cutoff: str


@dataclass
class JensenRow:
    kind: str
    lam: str
    d: int
    n: int
    positive_root_count: int
    expected: int
    ok: bool
    dps: int
    rational_digits: int
    n_sum: int
    cutoff: str


def parse_mpf_list(text: str) -> list[mp.mpf]:
    return [mp.mpf(part.strip()) for part in text.split(",") if part.strip()]


def phi_positive(u: mp.mpf, n_sum: int) -> mp.mpf:
    """Phi(u) for u >= 0 using the standard Newman/Riemann kernel series."""
    u = mp.mpf(u)
    total = mp.mpf("0")
    pi = mp.pi
    e4u = mp.exp(4 * u)
    e5u = mp.exp(5 * u)
    e9u = mp.exp(9 * u)
    for r in range(1, n_sum + 1):
        rr = mp.mpf(r)
        term = (2 * pi**2 * rr**4 * e9u - 3 * pi * rr**2 * e5u)
        term *= mp.exp(-pi * rr**2 * e4u)
        total += term
    return total


def phi_lambda(u: mp.mpf, lam: mp.mpf, n_sum: int) -> mp.mpf:
    u = abs(mp.mpf(u))
    return mp.exp(lam * u * u) * phi_positive(u, n_sum)


def integration_intervals(cutoff: mp.mpf) -> list[mp.mpf]:
    base = [mp.mpf("0"), mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf("1"), mp.mpf("2"), mp.mpf("4")]
    out = [x for x in base if x < cutoff]
    if not out or out[-1] != cutoff:
        out.append(cutoff)
    return out


def moment_even(k: int, lam: mp.mpf, n_sum: int, cutoff: mp.mpf) -> mp.mpf:
    """mu_{2k}(lambda) = integral_R u^(2k) Phi_lambda(u) du."""
    k = int(k)

    def integrand(u: mp.mpf) -> mp.mpf:
        return (u ** (2 * k)) * phi_lambda(u, lam, n_sum)

    return 2 * mp.quad(integrand, integration_intervals(cutoff))


def compute_moments(max_k: int, lam: mp.mpf, n_sum: int, cutoff: mp.mpf) -> list[mp.mpf]:
    return [moment_even(k, lam, n_sum, cutoff) for k in range(max_k + 1)]


def compute_a_sequence(moments: list[mp.mpf]) -> list[mp.mpf]:
    return [moment * mp.factorial(k) / mp.factorial(2 * k) for k, moment in enumerate(moments)]


def sigma(m: int) -> int:
    return -1 if ((m * (m + 1) // 2) % 2) else 1


def sign(x: mp.mpf) -> int:
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def log10_abs(x: mp.mpf) -> str:
    if x == 0:
        return "-inf"
    return mp.nstr(mp.log10(abs(x)), 20)


def hankel_det(a_values: list[mp.mpf], m: int, shift: int) -> mp.mpf:
    mat = mp.matrix([[a_values[i + j + shift] for j in range(m + 1)] for i in range(m + 1)])
    return mp.det(mat)


def rationalize(x: mp.mpf, digits: int) -> sp.Rational:
    return sp.Rational(mp.nstr(x, digits, min_fixed=-10, max_fixed=10))


def q_coeffs_ascending(a_values: list[mp.mpf], n: int, d: int) -> list[mp.mpf]:
    return [mp.mpf(comb(d, k)) * a_values[n + k] * ((-1) ** k) for k in range(d + 1)]


def positive_root_count_rationalized(coeffs_ascending: list[mp.mpf], digits: int) -> int:
    """Count positive roots of Q using exact Sturm on rationalized coefficients."""
    coeffs_desc = list(reversed(coeffs_ascending))
    rationals = [rationalize(c, digits) for c in coeffs_desc]
    y = sp.Symbol("y")
    poly = sp.Poly.from_list(rationals, gens=y, domain=sp.QQ)
    return int(poly.count_roots(sp.Rational(0), sp.oo))


def write_jsonl(path: Path, rows: list[dict]) -> None:
    with path.open("a", encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    exists = path.exists()
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def write_coeff_cache(path: Path, lam: mp.mpf, a_values: list[mp.mpf], args: argparse.Namespace, max_k: int) -> None:
    row = {
        "kind": "coefficients",
        "lam": mp.nstr(lam, 30),
        "max_k": max_k,
        "dps": args.dps,
        "n_sum": args.n_sum,
        "cutoff": mp.nstr(args.cutoff, 30),
        "A": [mp.nstr(value, args.coeff_digits) for value in a_values],
    }
    write_jsonl(path, [row])


def load_coeff_cache(path: Path, lam: mp.mpf, needed_max_k: int) -> tuple[list[mp.mpf], int, dict] | None:
    if not path.exists():
        raise FileNotFoundError(f"Coefficient cache not found: {path}")

    best: tuple[list[mp.mpf], int, dict] | None = None
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") != "coefficients":
                continue
            try:
                cached_lam = mp.mpf(row["lam"])
            except Exception:
                continue
            if cached_lam != lam:
                continue
            max_k = int(row["max_k"])
            values = [mp.mpf(x) for x in row["A"]]
            if max_k >= needed_max_k and len(values) >= needed_max_k + 1:
                if best is None or max_k < best[1]:
                    best = (values, max_k, row)
    return best


def compute_needed_max_k(args: argparse.Namespace) -> int:
    max_for_hankel = 2 * args.hankel_m_max + args.shift_max if args.mode in ("hankel", "both") else 0
    max_for_jensen = args.n_max + args.d_max if args.mode in ("jensen", "both") else 0
    return max(max_for_hankel, max_for_jensen)


def get_a_sequence(args: argparse.Namespace, lam: mp.mpf, max_k: int) -> tuple[list[mp.mpf], int, str]:
    if args.coeff_cache is not None:
        cached = load_coeff_cache(args.coeff_cache, lam, max_k)
        if cached is not None:
            values, cached_max_k, _row = cached
            return values, cached_max_k, f"cache:{args.coeff_cache}"
        if args.require_coeff_cache:
            raise RuntimeError(f"No usable coefficient cache row for lambda={mp.nstr(lam, 30)} max_k>={max_k}")

    moments = compute_moments(max_k, lam, args.n_sum, args.cutoff)
    return compute_a_sequence(moments), max_k, "computed"


def run_for_lambda(args: argparse.Namespace, lam: mp.mpf) -> tuple[list[HankelRow], list[JensenRow], list[mp.mpf], int, str]:
    max_k = compute_needed_max_k(args)
    a_values, cache_max_k, coeff_source = get_a_sequence(args, lam, max_k)

    hankel_rows: list[HankelRow] = []
    if args.mode in ("hankel", "both"):
        for shift in range(args.shift_min, args.shift_max + 1):
            for m in range(args.hankel_m_min, args.hankel_m_max + 1):
                detv = hankel_det(a_values, m, shift)
                sig = sigma(m)
                signed = sig * detv
                log_abs = log10_abs(detv)
                near_zero = detv == 0 or mp.log10(abs(detv)) <= args.near_zero_log10
                hankel_rows.append(
                    HankelRow(
                        kind="hankel",
                        lam=mp.nstr(lam, 30),
                        m=m,
                        shift=shift,
                        sigma=sig,
                        det_sign=sign(detv),
                        signed_sign=sign(signed),
                        log10_abs_det=log_abs,
                        near_zero=near_zero,
                        ok=signed > 0,
                        dps=args.dps,
                        n_sum=args.n_sum,
                        cutoff=mp.nstr(args.cutoff, 30),
                    )
                )

    jensen_rows: list[JensenRow] = []
    if args.mode in ("jensen", "both"):
        for d in range(args.d_min, args.d_max + 1):
            for n in range(args.n_min, args.n_max + 1):
                coeffs = q_coeffs_ascending(a_values, n, d)
                count = positive_root_count_rationalized(coeffs, args.rational_digits)
                jensen_rows.append(
                    JensenRow(
                        kind="jensen",
                        lam=mp.nstr(lam, 30),
                        d=d,
                        n=n,
                        positive_root_count=count,
                        expected=d,
                        ok=count == d,
                        dps=args.dps,
                        rational_digits=args.rational_digits,
                        n_sum=args.n_sum,
                        cutoff=mp.nstr(args.cutoff, 30),
                    )
                )

    return hankel_rows, jensen_rows, a_values, cache_max_k, coeff_source


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["hankel", "jensen", "both"], default="hankel")
    parser.add_argument("--lambdas", default="0")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--rational-digits", type=int, default=60)
    parser.add_argument("--n-sum", type=int, default=80)
    parser.add_argument("--cutoff", type=mp.mpf, default=mp.mpf("8"))
    parser.add_argument("--out-dir", type=Path, default=Path("work/rh_compute/results"))
    parser.add_argument("--run-name", default="run")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--write-coeff-cache", action="store_true")
    parser.add_argument("--coeff-cache", type=Path, default=None)
    parser.add_argument("--require-coeff-cache", action="store_true")
    parser.add_argument("--coeff-digits", type=int, default=80)
    parser.add_argument("--near-zero-log10", type=mp.mpf, default=mp.mpf("-250"))
    parser.add_argument("--hankel-m-min", type=int, default=0)
    parser.add_argument("--hankel-m-max", type=int, default=4)
    parser.add_argument("--shift-min", type=int, default=0)
    parser.add_argument("--shift-max", type=int, default=2)
    parser.add_argument("--d-min", type=int, default=2)
    parser.add_argument("--d-max", type=int, default=4)
    parser.add_argument("--n-min", type=int, default=0)
    parser.add_argument("--n-max", type=int, default=5)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    mp.mp.dps = args.dps
    args.out_dir.mkdir(parents=True, exist_ok=True)

    jsonl_path = args.out_dir / f"{args.run_name}.jsonl"
    csv_path = args.out_dir / f"{args.run_name}.csv"
    summary_path = args.out_dir / f"{args.run_name}_summary.json"
    coeff_path = args.out_dir / f"{args.run_name}_coefficients.jsonl"

    if args.overwrite:
        for path in (jsonl_path, csv_path, summary_path, coeff_path):
            if path.exists():
                path.unlink()

    started = perf_counter()
    total_hankel = 0
    total_jensen = 0
    failed_hankel = 0
    failed_jensen = 0
    near_zero_hankel = 0
    coeff_sources: list[str] = []

    for lam in parse_mpf_list(args.lambdas):
        hankel_rows, jensen_rows, a_values, max_k, coeff_source = run_for_lambda(args, lam)
        coeff_sources.append(coeff_source)
        rows = [asdict(row) for row in hankel_rows] + [asdict(row) for row in jensen_rows]
        write_jsonl(jsonl_path, rows)
        write_csv(csv_path, rows)
        if args.write_coeff_cache:
            write_coeff_cache(coeff_path, lam, a_values, args, max_k)
        total_hankel += len(hankel_rows)
        total_jensen += len(jensen_rows)
        failed_hankel += sum(1 for row in hankel_rows if not row.ok)
        failed_jensen += sum(1 for row in jensen_rows if not row.ok)
        near_zero_hankel += sum(1 for row in hankel_rows if row.near_zero)
        print(
            f"lambda={mp.nstr(lam, 12)} "
            f"coeff_source={coeff_source} "
            f"hankel={len(hankel_rows)} fail={sum(1 for row in hankel_rows if not row.ok)} "
            f"near_zero={sum(1 for row in hankel_rows if row.near_zero)} "
            f"jensen={len(jensen_rows)} fail={sum(1 for row in jensen_rows if not row.ok)}"
        )

    summary = {
        "mode": args.mode,
        "lambdas": [mp.nstr(x, 30) for x in parse_mpf_list(args.lambdas)],
        "dps": args.dps,
        "rational_digits": args.rational_digits,
        "n_sum": args.n_sum,
        "cutoff": mp.nstr(args.cutoff, 30),
        "coeff_cache": str(args.coeff_cache) if args.coeff_cache is not None else None,
        "coeff_sources": coeff_sources,
        "total_hankel": total_hankel,
        "failed_hankel": failed_hankel,
        "near_zero_hankel": near_zero_hankel,
        "total_jensen": total_jensen,
        "failed_jensen": failed_jensen,
        "elapsed_seconds": perf_counter() - started,
        "jsonl": str(jsonl_path),
        "csv": str(csv_path),
        "coefficients_jsonl": str(coeff_path) if args.write_coeff_cache else None,
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if failed_hankel == 0 and failed_jensen == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
