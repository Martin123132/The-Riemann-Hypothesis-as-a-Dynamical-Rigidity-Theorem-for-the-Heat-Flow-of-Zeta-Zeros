#!/usr/bin/env python
"""
Analytic truncation/cutoff tail bounds for Newman-kernel moment coefficients.

This bounds two errors for

  mu_{2k}(lambda) = 2 integral_0^infty u^(2k) exp(lambda u^2) Phi(u) du

using the standard kernel series

  Phi(u) = sum_{n>=1} (2*pi^2*n^4*exp(9u) - 3*pi*n^2*exp(5u))
                   * exp(-pi*n^2*exp(4u)).

The bounds cover:
  1. truncating the n-series after n_sum on 0 <= u <= cutoff;
  2. truncating the u-integral after cutoff, with the infinite n-series.

They do not bound quadrature error for the retained finite sum on
[0, cutoff]. This script is a tail/cutoff certification scaffold, not a
complete coefficient enclosure.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from math import factorial
from pathlib import Path

import mpmath as mp


@dataclass
class TailBoundRow:
    kind: str
    lam: str
    k: int
    n_sum: int
    cutoff: str
    dps: int
    mu_series_tail_bound: str
    mu_cutoff_tail_bound: str
    mu_total_tail_bound: str
    A_tail_bound: str
    c_tail_bound: str
    log10_mu_series_tail_bound: str
    log10_mu_cutoff_tail_bound: str
    log10_mu_total_tail_bound: str
    log10_A_tail_bound: str
    log10_c_tail_bound: str
    cutoff_monotone_gate: bool


def parse_mpf_list(text: str) -> list[mp.mpf]:
    return [mp.mpf(part.strip()) for part in text.split(",") if part.strip()]


def log10_text(value: mp.mpf, digits: int = 30) -> str:
    if value == 0:
        return "-inf"
    return mp.nstr(mp.log10(value), digits)


def nstr(value: mp.mpf, digits: int) -> str:
    return mp.nstr(value, digits, min_fixed=-10, max_fixed=10)


def gaussian_power_tail_bound(p: int, lower: mp.mpf) -> mp.mpf:
    """Bound integral_lower^infty x^p exp(-pi*x^2) dx."""
    s = mp.mpf(p + 1) / 2
    return mp.mpf("0.5") * mp.power(mp.pi, -s) * mp.gammainc(s, mp.pi * lower * lower, mp.inf)


def series_tail_constant(p: int, n_sum: int) -> mp.mpf:
    """Bound sum_{n>n_sum} n^p exp(-pi*n^2)."""
    threshold = mp.sqrt(mp.mpf(p) / (2 * mp.pi))
    lower = mp.mpf(n_sum)
    if lower < threshold:
        raise ValueError(f"n_sum={n_sum} too small for monotone integral tail bound with p={p}")
    return gaussian_power_tail_bound(p, lower)


def shifted_total_constant(p: int, terms: int) -> mp.mpf:
    """Bound K_p = sum_{n>=1} n^p exp(-pi*(n^2-1))."""
    total = mp.mpf("0")
    for n in range(1, terms + 1):
        nn = mp.mpf(n)
        total += nn**p * mp.exp(-mp.pi * (nn * nn - 1))
    tail = mp.exp(mp.pi) * gaussian_power_tail_bound(p, mp.mpf(terms))
    return total + tail


def crude_finite_interval_weight(k: int, lam: mp.mpf, a: int, cutoff: mp.mpf) -> mp.mpf:
    """Bound integral_0^cutoff u^(2k) exp(lambda*u^2 + a*u) du."""
    lam_pos = max(lam, mp.mpf("0"))
    return (
        mp.power(cutoff, 2 * k + 1)
        / mp.mpf(2 * k + 1)
        * mp.exp(lam_pos * cutoff * cutoff + mp.mpf(a) * cutoff)
    )


def cutoff_endpoint_integral_bound(k: int, lam: mp.mpf, a: int, cutoff: mp.mpf) -> tuple[mp.mpf, bool]:
    """Bound integral_cutoff^infty u^(2k) exp(lambda*u^2 + a*u - pi*exp(4u)) du.

    If g(u)=pi*exp(4u)-lambda*u^2-a*u-2k*log(u), then
    integral exp(-g(u)) du <= exp(-g(cutoff))/g'(cutoff)
    when g' is positive and increasing. For the programme's cutoff=8
    this gate is very strongly satisfied.
    """
    if cutoff <= 0:
        raise ValueError("cutoff must be positive for endpoint cutoff-tail bound")
    deriv = 4 * mp.pi * mp.exp(4 * cutoff) - 2 * lam * cutoff - mp.mpf(a) - (mp.mpf(2 * k) / cutoff)
    second = 16 * mp.pi * mp.exp(4 * cutoff) - 2 * lam + (mp.mpf(2 * k) / (cutoff * cutoff))
    gate = deriv > 0 and second > 0
    exponent = lam * cutoff * cutoff + mp.mpf(a) * cutoff + mp.mpf(2 * k) * mp.log(cutoff) - mp.pi * mp.exp(
        4 * cutoff
    )
    if not gate:
        return mp.inf, False
    return mp.exp(exponent) / deriv, True


def bound_for_k(
    k: int,
    lam: mp.mpf,
    n_sum: int,
    cutoff: mp.mpf,
    constant_terms: int,
    digits: int,
    dps: int,
) -> TailBoundRow:
    s2_tail = series_tail_constant(2, n_sum)
    s4_tail = series_tail_constant(4, n_sum)
    k2_total = shifted_total_constant(2, constant_terms)
    k4_total = shifted_total_constant(4, constant_terms)

    finite_9 = crude_finite_interval_weight(k, lam, 9, cutoff)
    finite_5 = crude_finite_interval_weight(k, lam, 5, cutoff)
    mu_series_tail = 2 * (2 * mp.pi**2 * s4_tail * finite_9 + 3 * mp.pi * s2_tail * finite_5)

    cutoff_9, gate_9 = cutoff_endpoint_integral_bound(k, lam, 9, cutoff)
    cutoff_5, gate_5 = cutoff_endpoint_integral_bound(k, lam, 5, cutoff)
    mu_cutoff_tail = 2 * (2 * mp.pi**2 * k4_total * cutoff_9 + 3 * mp.pi * k2_total * cutoff_5)

    mu_total_tail = mu_series_tail + mu_cutoff_tail
    A_tail = mu_total_tail * factorial(k) / factorial(2 * k)
    c_tail = mu_total_tail / factorial(2 * k)

    return TailBoundRow(
        kind="coefficient_tail_bound",
        lam=nstr(lam, digits),
        k=k,
        n_sum=n_sum,
        cutoff=nstr(cutoff, digits),
        dps=dps,
        mu_series_tail_bound=nstr(mu_series_tail, digits),
        mu_cutoff_tail_bound=nstr(mu_cutoff_tail, digits),
        mu_total_tail_bound=nstr(mu_total_tail, digits),
        A_tail_bound=nstr(A_tail, digits),
        c_tail_bound=nstr(c_tail, digits),
        log10_mu_series_tail_bound=log10_text(mu_series_tail),
        log10_mu_cutoff_tail_bound=log10_text(mu_cutoff_tail),
        log10_mu_total_tail_bound=log10_text(mu_total_tail),
        log10_A_tail_bound=log10_text(A_tail),
        log10_c_tail_bound=log10_text(c_tail),
        cutoff_monotone_gate=gate_9 and gate_5,
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


def summarize(rows: list[TailBoundRow], args: argparse.Namespace) -> dict:
    def max_row(field: str) -> TailBoundRow:
        return max(rows, key=lambda row: mp.mpf(getattr(row, field)))

    max_mu = max_row("mu_total_tail_bound")
    max_a = max_row("A_tail_bound")
    max_c = max_row("c_tail_bound")
    return {
        "kind": "coefficient_tail_bound_summary",
        "lambdas": [nstr(lam, args.digits) for lam in parse_mpf_list(args.lambdas)],
        "max_k": args.max_k,
        "n_sum": args.n_sum,
        "cutoff": nstr(mp.mpf(args.cutoff), args.digits),
        "dps": args.dps,
        "constant_terms": args.constant_terms,
        "target_abs_radius": args.target_abs_radius,
        "all_cutoff_monotone_gates": all(row.cutoff_monotone_gate for row in rows),
        "max_mu_total_tail": {
            "lam": max_mu.lam,
            "k": max_mu.k,
            "bound": max_mu.mu_total_tail_bound,
            "log10_bound": max_mu.log10_mu_total_tail_bound,
        },
        "max_A_tail": {
            "lam": max_a.lam,
            "k": max_a.k,
            "bound": max_a.A_tail_bound,
            "log10_bound": max_a.log10_A_tail_bound,
        },
        "max_c_tail": {
            "lam": max_c.lam,
            "k": max_c.k,
            "bound": max_c.c_tail_bound,
            "log10_bound": max_c.log10_c_tail_bound,
        },
        "target_abs_radius_log10": log10_text(mp.mpf(args.target_abs_radius)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lambdas", default="0,0.000001,0.0001,0.01,0.1")
    parser.add_argument("--max-k", type=int, default=32)
    parser.add_argument("--n-sum", type=int, default=100)
    parser.add_argument("--cutoff", default="8")
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--digits", type=int, default=50)
    parser.add_argument("--constant-terms", type=int, default=30)
    parser.add_argument("--target-abs-radius", default="1e-80")
    parser.add_argument("--output-dir", type=Path, default=Path("work/rh_compute/results"))
    parser.add_argument("--run-id", default="coefficient_tail_bounds")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.max_k < 0:
        raise ValueError("--max-k must be nonnegative")
    if args.n_sum < 1:
        raise ValueError("--n-sum must be positive")

    mp.mp.dps = args.dps
    cutoff = mp.mpf(args.cutoff)
    rows: list[TailBoundRow] = []
    for lam in parse_mpf_list(args.lambdas):
        for k in range(args.max_k + 1):
            rows.append(bound_for_k(k, lam, args.n_sum, cutoff, args.constant_terms, args.digits, args.dps))

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
