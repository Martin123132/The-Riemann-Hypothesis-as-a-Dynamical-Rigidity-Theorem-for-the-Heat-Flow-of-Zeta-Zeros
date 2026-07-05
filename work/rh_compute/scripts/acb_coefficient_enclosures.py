#!/usr/bin/env python
"""
Rigorous retained-integral enclosures for Newman-kernel coefficients.

This script uses python-flint's Arb/ACB rigorous integration to enclose

  mu_{2k}^{retained}(lambda)
    = 2 integral_0^cutoff u^(2k) exp(lambda u^2) Phi_{n_sum}(u) du

where Phi_{n_sum} is the finite n=1..n_sum kernel sum. It then adds the
analytic n-tail and cutoff-tail bounds from coefficient_tail_bounds.py to
enclose the full moment coefficient.

Output balls are produced for:

  A_k = mu_{2k} k! / (2k)!
  c_k = A_k / k! = mu_{2k} / (2k)!

This is intended as the coefficient-enclosure input to the structural-zero
Toeplitz/PF Arb probe.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import asdict, dataclass
from math import factorial
from pathlib import Path
from time import perf_counter

import mpmath as mp


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import acb, arb  # noqa: E402

from coefficient_tail_bounds import bound_for_k, parse_mpf_list  # noqa: E402


@dataclass
class EnclosureRow:
    kind: str
    lam: str
    k: int
    n_sum: int
    cutoff: str
    dps: int
    abs_tol: str
    retained_mu_ball: str
    retained_mu_rad: str
    tail_mu_bound: str
    full_mu_ball: str
    full_mu_rad: str
    A_ball: str
    A_rad: str
    c_ball: str
    c_rad: str
    cache_A: str | None
    cache_A_exact_contained: bool | None
    cache_A_radius: str | None
    cache_A_radius_contains_A_ball: bool | None
    cache_A_radius_overlaps_A_ball: bool | None
    elapsed_seconds: float


def nstr_mpf(value: mp.mpf, digits: int) -> str:
    return mp.nstr(value, digits, min_fixed=-10, max_fixed=10)


def load_coeff_cache(path: Path | None, lam: mp.mpf, needed_max_k: int) -> list[str] | None:
    if path is None:
        return None
    best: list[str] | None = None
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
            values = list(row.get("A", []))
            if len(values) >= needed_max_k + 1:
                if best is None or len(values) < len(best):
                    best = values
    return best


def phi_finite_acb(u: acb, n_sum: int, pi: acb) -> acb:
    total = acb(0)
    e4 = (4 * u).exp()
    e5 = (5 * u).exp()
    e9 = (9 * u).exp()
    for n in range(1, n_sum + 1):
        nn = acb(n)
        total += (2 * pi * pi * nn**4 * e9 - 3 * pi * nn**2 * e5) * (-pi * nn**2 * e4).exp()
    return total


def retained_mu_ball(k: int, lam_text: str, n_sum: int, cutoff_text: str, abs_tol_text: str) -> acb:
    pi = acb.pi()
    lam = acb(lam_text)
    cutoff = acb(cutoff_text)
    zero = acb(0)
    abs_tol = arb(abs_tol_text)

    def integrand(u: acb, analytic: bool) -> acb:
        return (u ** (2 * k)) * (lam * u * u).exp() * phi_finite_acb(u, n_sum, pi)

    return 2 * acb.integral(
        integrand,
        zero,
        cutoff,
        abs_tol=abs_tol,
        rel_tol=abs_tol,
        deg_limit=64,
        eval_limit=200000,
        depth_limit=64,
    )


def add_symmetric_radius(x: arb, radius_text: str) -> arb:
    radius = arb(radius_text)
    return x + arb(f"[0 +/- {radius_text}]")


def arb_upper_bound_text(x: arb, digits: int = 40) -> str:
    mid, rad, exp = x.mid_rad_10exp()
    upper = mp.mpf(abs(int(mid)) + int(rad)) * mp.power(10, int(exp))
    return mp.nstr(upper, digits, min_fixed=-10, max_fixed=10)


def ball_radius_text(x: arb, digits: int = 40) -> str:
    return arb_upper_bound_text(x.rad(), digits)


def contains_real_text(ball: arb, value_text: str | None) -> bool | None:
    if value_text is None:
        return None
    return bool(ball.contains(arb(value_text)))


def cache_radius_ball(value_text: str | None, radius_text: str | None) -> arb | None:
    if value_text is None or radius_text is None:
        return None
    return arb(f"[{value_text} +/- {radius_text}]")


def enclose_one(
    k: int,
    lam: mp.mpf,
    args: argparse.Namespace,
    cache_values: list[str] | None,
) -> EnclosureRow:
    start = perf_counter()
    lam_text = nstr_mpf(lam, args.digits)
    retained = retained_mu_ball(k, lam_text, args.n_sum, args.cutoff, args.abs_tol)
    if not retained.imag.contains(0):
        raise RuntimeError(f"Retained integral has non-real enclosure for lambda={lam_text}, k={k}: {retained}")
    retained_real = retained.real

    tail = bound_for_k(
        k,
        lam,
        args.n_sum,
        mp.mpf(args.cutoff),
        args.constant_terms,
        args.digits,
        args.dps,
    )
    full_mu = add_symmetric_radius(retained_real, tail.mu_total_tail_bound)
    A_ball = full_mu * factorial(k) / factorial(2 * k)
    c_ball = full_mu / factorial(2 * k)

    cache_A = None
    if cache_values is not None and k < len(cache_values):
        cache_A = cache_values[k]
    cache_ball = cache_radius_ball(cache_A, args.cache_radius) if args.cache_radius else None

    elapsed = perf_counter() - start
    return EnclosureRow(
        kind="acb_coefficient_enclosure",
        lam=lam_text,
        k=k,
        n_sum=args.n_sum,
        cutoff=args.cutoff,
        dps=args.dps,
        abs_tol=args.abs_tol,
        retained_mu_ball=str(retained_real),
        retained_mu_rad=ball_radius_text(retained_real),
        tail_mu_bound=tail.mu_total_tail_bound,
        full_mu_ball=str(full_mu),
        full_mu_rad=ball_radius_text(full_mu),
        A_ball=str(A_ball),
        A_rad=ball_radius_text(A_ball),
        c_ball=str(c_ball),
        c_rad=ball_radius_text(c_ball),
        cache_A=cache_A,
        cache_A_exact_contained=contains_real_text(A_ball, cache_A),
        cache_A_radius=args.cache_radius if cache_A is not None else None,
        cache_A_radius_contains_A_ball=bool(cache_ball.contains(A_ball)) if cache_ball is not None else None,
        cache_A_radius_overlaps_A_ball=bool(cache_ball.overlaps(A_ball)) if cache_ball is not None else None,
        elapsed_seconds=elapsed,
    )


def write_jsonl(path: Path, rows: list[dict], overwrite: bool) -> None:
    mode = "w" if overwrite else "a"
    with path.open(mode, encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict], overwrite: bool) -> None:
    if not rows:
        return
    exists = path.exists() and not overwrite
    mode = "a" if exists else "w"
    fieldnames = sorted({key for row in rows for key in row.keys()})
    with path.open(mode, newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        if not exists:
            writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[EnclosureRow], args: argparse.Namespace) -> dict:
    def max_rad(field: str) -> EnclosureRow:
        return max(rows, key=lambda row: mp.mpf(getattr(row, field)))

    max_a = max_rad("A_rad")
    max_c = max_rad("c_rad")
    cache_checks = [
        row.cache_A_radius_contains_A_ball
        for row in rows
        if row.cache_A_radius_contains_A_ball is not None
    ]
    return {
        "kind": "acb_coefficient_enclosure_summary",
        "lambdas": [nstr_mpf(lam, args.digits) for lam in parse_mpf_list(args.lambdas)],
        "k_min": args.k_min,
        "k_max": args.k_max,
        "n_sum": args.n_sum,
        "cutoff": args.cutoff,
        "dps": args.dps,
        "abs_tol": args.abs_tol,
        "rows": len(rows),
        "all_cache_A_radius_contains_A_ball": all(cache_checks) if cache_checks else None,
        "cache_checks": len(cache_checks),
        "max_A_radius": {
            "lam": max_a.lam,
            "k": max_a.k,
            "radius": max_a.A_rad,
            "log10_radius": mp.nstr(mp.log10(mp.mpf(max_a.A_rad)), 30),
        },
        "max_c_radius": {
            "lam": max_c.lam,
            "k": max_c.k,
            "radius": max_c.c_rad,
            "log10_radius": mp.nstr(mp.log10(mp.mpf(max_c.c_rad)), 30),
        },
        "elapsed_seconds_total": sum(row.elapsed_seconds for row in rows),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lambdas", default="0")
    parser.add_argument("--k-min", type=int, default=0)
    parser.add_argument("--k-max", type=int, default=4)
    parser.add_argument("--n-sum", type=int, default=100)
    parser.add_argument("--cutoff", default="8")
    parser.add_argument("--dps", type=int, default=120)
    parser.add_argument("--digits", type=int, default=60)
    parser.add_argument("--abs-tol", default="1e-90")
    parser.add_argument("--constant-terms", type=int, default=40)
    parser.add_argument("--coeff-cache", type=Path, default=None)
    parser.add_argument("--cache-radius", default="1e-80")
    parser.add_argument("--output-dir", type=Path, default=Path("work/rh_compute/results"))
    parser.add_argument("--run-id", default="acb_coefficient_enclosures")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.k_min < 0 or args.k_max < args.k_min:
        raise ValueError("Require 0 <= --k-min <= --k-max")
    if args.n_sum < 1:
        raise ValueError("--n-sum must be positive")

    flint.ctx.dps = args.dps
    mp.mp.dps = max(args.dps, args.digits + 20)

    rows: list[EnclosureRow] = []
    for lam in parse_mpf_list(args.lambdas):
        cache_values = load_coeff_cache(args.coeff_cache, lam, args.k_max) if args.coeff_cache else None
        for k in range(args.k_min, args.k_max + 1):
            row = enclose_one(k, lam, args, cache_values)
            rows.append(row)
            print(
                f"lambda={row.lam} k={k} A_rad={row.A_rad} "
                f"cache_radius_contains={row.cache_A_radius_contains_A_ball} elapsed={row.elapsed_seconds:.3f}s"
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    row_dicts = [asdict(row) for row in rows]
    summary = summarize(rows, args)

    jsonl_path = args.output_dir / f"{args.run_id}.jsonl"
    csv_path = args.output_dir / f"{args.run_id}.csv"
    summary_path = args.output_dir / f"{args.run_id}_summary.json"

    write_jsonl(jsonl_path, row_dicts, overwrite=args.overwrite)
    write_csv(csv_path, row_dicts, overwrite=args.overwrite)
    with summary_path.open("w" if args.overwrite else "x", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")

    print(json.dumps(summary, sort_keys=True))
    print(f"rows={jsonl_path}")
    print(f"csv={csv_path}")
    print(f"summary={summary_path}")


if __name__ == "__main__":
    main()
