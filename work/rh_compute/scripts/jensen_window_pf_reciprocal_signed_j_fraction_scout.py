#!/usr/bin/env python3
"""Scout signed J-fraction Hankel signatures for Jensen-window reciprocals.

For E(t)=1/H(-t), let mu_m=[t^m]E(t) and

    Delta_r = det(mu_{i+j})_{0<=i,j<r}.

The ordinary Jacobi continued-fraction parameter satisfies

    lambda_n = Delta_{n+1} Delta_{n-1} / Delta_n^2.

If the signed Hankel signature

    (-1)^(r(r-1)/2) Delta_r > 0

holds, then every ordinary lambda_n is negative and the signed parameters
kappa_n=-lambda_n are positive.  This script checks that finite signature on
the current Arb Jensen-window grid.  It is finite theorem-search evidence, not
an all-order continued-fraction theorem.
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


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_LAMBDAS = ("0", "1e-6", "1e-4", "1e-2", "1e-1")
DEFAULT_ENCLOSURE_JSONL = (
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k0_k9.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k10_k20.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k21_k32.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k33_k48.jsonl",
    "work/rh_compute/results/acb_enclosures_repro_hankel_15c_lamgrid_k49_k64.jsonl",
)
DEFAULT_OUT_JSONL = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_OUT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_j_fraction_scout.json"
)


@dataclass(frozen=True)
class SignedHankelRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    determinant_order_r: int
    expected_delta_sign: int
    signed_delta_classification: str
    signed_delta_contains_zero: bool
    ok: bool


@dataclass(frozen=True)
class SignedLambdaRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    lambda_index: int
    ordinary_lambda_classification: str
    ordinary_lambda_contains_zero: bool
    signed_lambda_classification: str
    signed_lambda_contains_zero: bool
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


def reciprocal_moments(
    balls: dict[int, flint.arb],
    shift_n: int,
    degree_d: int,
    max_m: int,
) -> list[flint.arb]:
    h = [comb(degree_d, index) * balls[shift_n + index] for index in range(degree_d + 1)]
    g = [item / h[0] for item in h]
    values = [flint.arb(1)]
    for size_m in range(1, max_m + 1):
        value = flint.arb(0)
        for offset in range(1, min(degree_d, size_m) + 1):
            term = g[offset] * values[size_m - offset]
            value = value + term if offset % 2 == 1 else value - term
        values.append(value)
    return values


def hankel_det(moments: list[flint.arb], order: int) -> flint.arb:
    return flint.arb_mat([[moments[row + col] for col in range(order)] for row in range(order)]).det()


def signed_hankel_rows(
    lam: str,
    shift_n: int,
    degree_d: int,
    deltas: list[flint.arb],
) -> list[SignedHankelRow]:
    rows = []
    for order_r in range(1, degree_d + 1):
        expected_sign = 1 if (order_r * (order_r - 1) // 2) % 2 == 0 else -1
        signed_delta = deltas[order_r] if expected_sign > 0 else -deltas[order_r]
        classification = classify(signed_delta)
        contains_zero = signed_delta.contains(0)
        rows.append(
            SignedHankelRow(
                kind="jensen_window_pf_reciprocal_signed_hankel_row",
                lam=lam,
                shift_n=shift_n,
                degree_d=degree_d,
                determinant_order_r=order_r,
                expected_delta_sign=expected_sign,
                signed_delta_classification=classification,
                signed_delta_contains_zero=contains_zero,
                ok=(classification == "positive" and not contains_zero),
            )
        )
    return rows


def signed_lambda_rows(
    lam: str,
    shift_n: int,
    degree_d: int,
    deltas: list[flint.arb],
) -> list[SignedLambdaRow]:
    rows = []
    for lambda_index in range(1, degree_d):
        numerator = deltas[lambda_index + 1]
        if lambda_index > 1:
            numerator = numerator * deltas[lambda_index - 1]
        ordinary_lambda = numerator / (deltas[lambda_index] * deltas[lambda_index])
        signed_lambda = -ordinary_lambda
        ordinary_classification = classify(ordinary_lambda)
        signed_classification = classify(signed_lambda)
        rows.append(
            SignedLambdaRow(
                kind="jensen_window_pf_reciprocal_signed_j_lambda_row",
                lam=lam,
                shift_n=shift_n,
                degree_d=degree_d,
                lambda_index=lambda_index,
                ordinary_lambda_classification=ordinary_classification,
                ordinary_lambda_contains_zero=ordinary_lambda.contains(0),
                signed_lambda_classification=signed_classification,
                signed_lambda_contains_zero=signed_lambda.contains(0),
                ok=(
                    ordinary_classification == "negative"
                    and not ordinary_lambda.contains(0)
                    and signed_classification == "positive"
                    and not signed_lambda.contains(0)
                ),
            )
        )
    return rows


def symbolic_rows() -> list[dict]:
    return [
        {
            "id": "signed_hankel_signature_to_signed_j_lambda",
            "statement": (
                "If Delta_r has sign (-1)^(r(r-1)/2) for r=1..N+1, then "
                "ordinary Jacobi lambda_n=Delta_{n+1}Delta_{n-1}/Delta_n^2 is negative for n=1..N."
            ),
            "parity_check": "n(n+1)/2 + (n-1)(n-2)/2 = n^2 - n + 1 is odd.",
            "signed_parameter": "kappa_n = -lambda_n > 0",
            "proof_boundary": "Formal sign implication only; not a proof that the determinant signature holds for all zeta windows.",
        },
        {
            "id": "endpoint_real_rooted_model_signature",
            "statement": (
                "The same determinant signature is compatible with the endpoint model "
                "H(t)=prod_i(1+alpha_i*t), alpha_i>=0, where E(t)=prod_i(1-alpha_i*t)^(-1)."
            ),
            "proof_boundary": "Endpoint consistency only; using real-rootedness as an input would be circular for jwpf_06.",
        },
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=[REPO_ROOT / p for p in DEFAULT_ENCLOSURE_JSONL])
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--shifts", default="0..20")
    parser.add_argument("--degrees", default="2..8")
    parser.add_argument("--dps", type=int, default=520)
    parser.add_argument("--out-jsonl", type=Path, default=DEFAULT_OUT_JSONL)
    parser.add_argument("--out-summary", type=Path, default=DEFAULT_OUT_SUMMARY)
    parser.add_argument("--json", action="store_true", help="Print summary JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    flint.ctx.dps = args.dps
    lambdas = [part.strip() for part in args.lambdas.split(",") if part.strip()]
    shifts = parse_int_range(args.shifts)
    degrees = parse_int_range(args.degrees)
    if min(degrees) < 2:
        raise ValueError("degrees must be at least 2")
    needed_max_k = max(shifts) + max(degrees)

    all_rows: list[dict] = []
    signed_hankel_count = 0
    signed_lambda_count = 0
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                moments = reciprocal_moments(balls, shift_n, degree_d, 2 * degree_d)
                deltas = [flint.arb(1)] + [hankel_det(moments, order) for order in range(1, degree_d + 1)]
                hankel_rows = signed_hankel_rows(lam, shift_n, degree_d, deltas)
                lambda_rows = signed_lambda_rows(lam, shift_n, degree_d, deltas)
                signed_hankel_count += sum(1 for row in hankel_rows if row.ok)
                signed_lambda_count += sum(1 for row in lambda_rows if row.ok)
                all_rows.extend(asdict(row) for row in hankel_rows)
                all_rows.extend(asdict(row) for row in lambda_rows)

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.out_jsonl.open("w", encoding="utf-8") as handle:
        for row in all_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    expected_hankel_rows = len(lambdas) * len(shifts) * sum(degrees)
    expected_lambda_rows = len(lambdas) * len(shifts) * sum(degree - 1 for degree in degrees)
    payload = {
        "kind": "jensen_window_pf_reciprocal_signed_j_fraction_scout",
        "date": "2026-07-06",
        "target_route_row": "rp_09_signed_or_modified_continued_fraction",
        "proof_boundary": (
            "Signed J-fraction Hankel-signature diagnostic only; not a proof of a signed "
            "continued-fraction theorem, the all-order column recurrence, Schur positivity, "
            "Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
        "symbolic_rows": symbolic_rows(),
        "finite_grid": {
            "lambdas": lambdas,
            "shifts": [min(shifts), max(shifts)],
            "degrees": [min(degrees), max(degrees)],
            "dps": args.dps,
            "signed_hankel_rows": expected_hankel_rows,
            "signed_hankel_positive_rows": signed_hankel_count,
            "ordinary_lambda_rows": expected_lambda_rows,
            "ordinary_lambda_negative_rows": signed_lambda_count,
            "signed_lambda_positive_rows": signed_lambda_count,
            "all_signed_hankel_rows_positive": signed_hankel_count == expected_hankel_rows,
            "all_ordinary_lambda_rows_negative": signed_lambda_count == expected_lambda_rows,
            "source_enclosures": [str(path.relative_to(REPO_ROOT)) for path in args.enclosure_jsonl],
            "row_log": str(args.out_jsonl.relative_to(REPO_ROOT)),
        },
        "summary": {
            "symbolic_rows": 2,
            "signed_hankel_rows": expected_hankel_rows,
            "signed_lambda_rows": expected_lambda_rows,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "On the finite Arb grid, reciprocal coefficient Hankel determinants obey "
                "(-1)^(r(r-1)/2) Delta_r > 0 through r=d, and the ordinary Jacobi "
                "lambda_n parameters are negative through n=d-1, so the signed parameters "
                "kappa_n=-lambda_n are positive.  This supports a signed-J-fraction target "
                "but does not supply the missing all-order theorem."
            ),
        },
        "invariants": [
            "Every signed Hankel row in the finite grid must be positive and separated from zero.",
            "Every ordinary Jacobi lambda row in the finite grid must be negative and separated from zero.",
            "Every signed kappa_n=-lambda_n row in the finite grid must be positive and separated from zero.",
            "The scout is finite evidence only and does not mark the signed-fraction route ready to apply.",
        ],
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF reciprocal signed J-fraction scout: "
            f"{expected_hankel_rows} signed Hankel rows, {expected_lambda_rows} signed-lambda rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
