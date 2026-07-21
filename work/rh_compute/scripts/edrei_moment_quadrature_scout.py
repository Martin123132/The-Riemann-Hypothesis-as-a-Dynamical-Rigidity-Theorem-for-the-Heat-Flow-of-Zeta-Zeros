#!/usr/bin/env python3
"""Arb recurrence scout for Edrei log-power moment data.

The Edrei-log diagnostic defines p_n = (-1)^(n-1) n [z^n] log H(z).  If

    H(z) = exp(gamma z) prod_j (1 + beta_j z),   gamma, beta_j >= 0,

then p_n for n >= 2 are moments of the positive measure sum_j delta_{beta_j}
weighted by beta_j^n.  Equivalently, the shifted sequence

    a_n = p_{n+1},  n >= 0,

should look like a Stieltjes moment sequence for x dnu(x).

This script reads existing Arb row logs and checks finite monic recurrence
data for that shifted sequence using Arb interval arithmetic.  It is a
theorem-search scout: positive finite recurrence data is not a proof of the
Edrei representation, PF-infinity, or RH.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import glob
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_LAMBDAS = ("0", "1e-6", "1e-4", "1e-2", "1e-1")


@dataclass(frozen=True)
class RecurrenceRow:
    kind: str
    lam: str
    order: int
    scale: str
    min_norm: str
    min_alpha: str
    min_beta: str
    status: str
    all_positive: bool
    proof_boundary: str


def decimal_equal(left: str, right: str) -> bool:
    return Decimal(str(left)) == Decimal(str(right))


def lam_label(lam: str) -> str:
    if Decimal(str(lam)) == Decimal("0"):
        return "lam0"
    return "lam" + lam.replace("-", "m").replace("+", "").replace(".", "p")


def load_powers(paths: list[Path], lam: str, needed_max_n: int) -> dict[int, flint.arb]:
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
                    powers[n] = flint.arb(row["signed_q_n"])
    missing = [n for n in range(1, needed_max_n + 1) if n not in powers]
    if missing:
        raise RuntimeError(f"missing p_n rows for lambda={lam}: {missing[:10]}")
    return powers


def inner(poly_a: list[flint.arb], poly_b: list[flint.arb], moments: list[flint.arb]) -> flint.arb:
    total = flint.arb(0)
    for i, ai in enumerate(poly_a):
        for j, bj in enumerate(poly_b):
            total += ai * bj * moments[i + j]
    return total


def x_times(poly: list[flint.arb]) -> list[flint.arb]:
    return [flint.arb(0)] + poly


def sub_scaled(left: list[flint.arb], right: list[flint.arb], scale: flint.arb) -> list[flint.arb]:
    n = max(len(left), len(right))
    out = [flint.arb(0) for _ in range(n)]
    for i, value in enumerate(left):
        out[i] += value
    for i, value in enumerate(right):
        out[i] -= scale * value
    return out


def classify_values(values: list[flint.arb]) -> str:
    if all(value > 0 and not value.contains(0) for value in values):
        return "positive"
    if any(value < 0 and not value.contains(0) for value in values):
        return "negative"
    if any(value.contains(0) for value in values):
        return "inconclusive_contains_zero"
    return "unknown"


def min_by_lower(values: list[flint.arb]) -> flint.arb:
    if not values:
        return flint.arb("inf")
    best = values[0]
    for value in values[1:]:
        if value.lower() < best.lower():
            best = value
    return best


def recurrence_data(moments: list[flint.arb], order: int) -> tuple[list[flint.arb], list[flint.arb], list[flint.arb]]:
    polys: list[list[flint.arb]] = [[flint.arb(1)]]
    norms: list[flint.arb] = [inner(polys[0], polys[0], moments)]
    alphas: list[flint.arb] = []
    betas: list[flint.arb] = []

    for n in range(order):
        norm_n = norms[n]
        if norm_n.contains(0) or not norm_n > 0:
            break
        p_n = polys[n]
        alpha_n = inner(x_times(p_n), p_n, moments) / norm_n
        alphas.append(alpha_n)
        if n == order - 1:
            break
        beta_n = flint.arb(0) if n == 0 else norms[n] / norms[n - 1]
        if n > 0:
            betas.append(beta_n)
        next_poly = sub_scaled(x_times(p_n), p_n, alpha_n)
        if n > 0:
            next_poly = sub_scaled(next_poly, polys[n - 1], beta_n)
        polys.append(next_poly)
        norms.append(inner(next_poly, next_poly, moments))
    return norms[:order], alphas, betas


def run_one(lam: str, powers: dict[int, flint.arb], order: int) -> RecurrenceRow:
    # For N-point recurrence data, alpha_{N-1} uses moments through a_{2N-1}.
    shifted = [powers[n + 1] for n in range(0, 2 * order)]
    mass = shifted[0]
    if mass.contains(0) or not mass > 0:
        raise ArithmeticError(f"nonpositive shifted mass for lambda={lam}: {mass}")
    scale = shifted[1] / shifted[0] if shifted[1] > 0 else flint.arb(1)
    moments = [shifted[n] / (mass * (scale ** n)) for n in range(0, 2 * order)]

    norms, alphas, betas = recurrence_data(moments, order)
    values = norms + alphas + betas
    status = classify_values(values)
    all_positive = status == "positive" and len(norms) == order and len(alphas) == order
    return RecurrenceRow(
        kind="edrei_moment_recurrence_scout",
        lam=lam,
        order=order,
        scale=str(scale),
        min_norm=str(min_by_lower(norms)),
        min_alpha=str(min_by_lower(alphas)) if alphas else "nan",
        min_beta=str(min_by_lower(betas)) if betas else "inf",
        status=status if len(norms) == order and len(alphas) == order else "inconclusive_contains_zero",
        all_positive=all_positive,
        proof_boundary=(
            "finite Arb recurrence scout for the Edrei log-power representation; "
            "not a proof of PF-infinity or RH"
        ),
    )


def default_log_paths(results_dir: Path, lambdas: list[str], needed_max_n: int) -> list[Path]:
    log_sign_dir = results_dir / "arb_edrei" / "log_sign"
    paths: list[Path] = []
    for lam in lambdas:
        label = lam_label(lam)
        patterns = []
        if needed_max_n <= 57:
            patterns.append(f"arb_edrei_log_sign_{label}_n1_n57_dps2400_boundary_tol1e-140.jsonl")
        patterns.extend(
            [
                f"arb_edrei_log_sign_{label}_n1_n64_dps340.jsonl",
                f"arb_edrei_log_sign_{label}_n1_n*_dps2400_boundary_tol1e-140.jsonl",
            ]
        )
        matches: list[str] = []
        for pattern in patterns:
            matches = sorted(glob.glob(str(log_sign_dir / pattern)))
            if matches:
                break
        if not matches:
            raise FileNotFoundError(f"no default Edrei log-sign row file for lambda={lam}")
        paths.append(Path(matches[-1]))
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=Path("work/rh_compute/results"))
    parser.add_argument("--edrei-jsonl", type=Path, nargs="*", default=None)
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--order-min", type=int, default=2)
    parser.add_argument("--order-max", type=int, default=12)
    parser.add_argument("--dps", type=int, default=2000)
    parser.add_argument("--out-jsonl", type=Path, default=None)
    parser.add_argument("--out-summary", type=Path, default=None)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.order_min < 1 or args.order_max < args.order_min:
        raise ValueError("invalid order range")
    flint.ctx.dps = args.dps

    lambdas = [part.strip() for part in args.lambdas.split(",") if part.strip()]
    needed_max_n = 2 * args.order_max
    paths = args.edrei_jsonl if args.edrei_jsonl else default_log_paths(args.results_dir, lambdas, needed_max_n)

    rows: list[RecurrenceRow] = []
    for lam in lambdas:
        powers = load_powers(paths, lam, needed_max_n)
        for order in range(args.order_min, args.order_max + 1):
            rows.append(run_one(lam, powers, order))

    nonpositive = [row for row in rows if row.status == "negative"]
    inconclusive = [row for row in rows if row.status != "positive" and row.status != "negative"]
    summary = {
        "kind": "edrei_moment_recurrence_scout_summary",
        "rigorous": True,
        "proof_boundary": (
            "finite recurrence diagnostic for positive Edrei zero-parameter moments; "
            "not a proof of PF-infinity or RH"
        ),
        "lambdas": lambdas,
        "order_min": args.order_min,
        "order_max": args.order_max,
        "needed_max_n": needed_max_n,
        "dps": args.dps,
        "rows": len(rows),
        "positive_rows": sum(1 for row in rows if row.status == "positive"),
        "negative_rows": len(nonpositive),
        "inconclusive_rows": len(inconclusive),
        "all_positive": all(row.status == "positive" for row in rows),
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
    for row in rows:
        if row.status != "positive":
            print(json.dumps(asdict(row), sort_keys=True))
    return 0 if not nonpositive else 1


if __name__ == "__main__":
    raise SystemExit(main())
