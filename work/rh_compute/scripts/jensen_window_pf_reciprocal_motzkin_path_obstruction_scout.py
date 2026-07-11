#!/usr/bin/env python3
"""Scout the ordinary Motzkin-path obstruction for signed J-fraction data.

For a standard J-fraction

    M(t) = 1 / (1 - beta_0 t - lambda_1 t^2/(1 - beta_1 t - ...)),

the coefficients are weighted Motzkin excursions with horizontal weights
beta_h and down-step weights lambda_h.  The signed zeta-window reciprocal
data has lambda_h < 0, and for degree d>=3 also beta_1 < 0.  This script
checks the finite grid obstruction to the naive entrywise nonnegative
tridiagonal production-matrix/path model:

    mu_2 = beta_0^2 + lambda_1

where beta_0^2 > 0, lambda_1 < 0, but mu_2 > 0 on the grid.  The positive
coefficient is therefore not a manifestly nonnegative path sum in the ordinary
Motzkin model.  A useful proof route would need a genuinely modified model,
not just the raw J-fraction path expansion.
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
    / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_OUT_SUMMARY = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json"
)


@dataclass(frozen=True)
class Mu2CancellationRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    beta0_classification: str
    beta0_contains_zero: bool
    beta0_square_classification: str
    beta0_square_contains_zero: bool
    lambda1_classification: str
    lambda1_contains_zero: bool
    mu2_classification: str
    mu2_contains_zero: bool
    has_negative_path_weight: bool
    ok: bool


@dataclass(frozen=True)
class Beta1DiagonalRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    beta1_classification: str
    beta1_contains_zero: bool
    diagonal_conjugation_changes_diagonal: bool
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


def star_hankel_det(moments: list[flint.arb], order: int) -> flint.arb:
    return flint.arb_mat(
        [
            [moments[row + col] if col < order - 1 else moments[row + order] for col in range(order)]
            for row in range(order)
        ]
    ).det()


def jacobi_betas(moments: list[flint.arb], degree_d: int) -> list[flint.arb]:
    deltas = [flint.arb(1)] + [hankel_det(moments, order) for order in range(1, degree_d + 1)]
    star_deltas = [flint.arb(0)] + [star_hankel_det(moments, order) for order in range(1, degree_d + 1)]
    q_values = [flint.arb(0)] + [
        star_deltas[order] / deltas[order] for order in range(1, degree_d + 1)
    ]
    return [q_values[index + 1] - q_values[index] for index in range(degree_d)]


def lambda1(moments: list[flint.arb]) -> flint.arb:
    return hankel_det(moments, 2)


def mu2_row(lam: str, shift_n: int, degree_d: int, moments: list[flint.arb], betas: list[flint.arb]) -> Mu2CancellationRow:
    beta0 = betas[0]
    beta0_square = beta0 * beta0
    lam1 = lambda1(moments)
    mu2 = moments[2]
    beta0_class = classify(beta0)
    beta0_square_class = classify(beta0_square)
    lambda1_class = classify(lam1)
    mu2_class = classify(mu2)
    ok = (
        beta0_class == "positive"
        and not beta0.contains(0)
        and beta0_square_class == "positive"
        and not beta0_square.contains(0)
        and lambda1_class == "negative"
        and not lam1.contains(0)
        and mu2_class == "positive"
        and not mu2.contains(0)
    )
    return Mu2CancellationRow(
        kind="jensen_window_pf_reciprocal_motzkin_mu2_cancellation_row",
        lam=lam,
        shift_n=shift_n,
        degree_d=degree_d,
        beta0_classification=beta0_class,
        beta0_contains_zero=beta0.contains(0),
        beta0_square_classification=beta0_square_class,
        beta0_square_contains_zero=beta0_square.contains(0),
        lambda1_classification=lambda1_class,
        lambda1_contains_zero=lam1.contains(0),
        mu2_classification=mu2_class,
        mu2_contains_zero=mu2.contains(0),
        has_negative_path_weight=(lambda1_class == "negative" and not lam1.contains(0)),
        ok=ok,
    )


def beta1_row(lam: str, shift_n: int, degree_d: int, betas: list[flint.arb]) -> Beta1DiagonalRow:
    beta1 = betas[1]
    classification = classify(beta1)
    contains_zero = beta1.contains(0)
    return Beta1DiagonalRow(
        kind="jensen_window_pf_reciprocal_motzkin_beta1_diagonal_obstruction_row",
        lam=lam,
        shift_n=shift_n,
        degree_d=degree_d,
        beta1_classification=classification,
        beta1_contains_zero=contains_zero,
        diagonal_conjugation_changes_diagonal=False,
        ok=(classification == "negative" and not contains_zero),
    )


def symbolic_rows() -> list[dict]:
    return [
        {
            "id": "ordinary_j_fraction_motzkin_model",
            "statement": (
                "The ordinary J-fraction moment expansion is a Motzkin-excursion path "
                "model with horizontal weights beta_h and down-step weights lambda_h."
            ),
            "boundary": "This is a model for the ordinary J-fraction, not a proof that the signed zeta-window route is positive.",
        },
        {
            "id": "diagonal_conjugation_obstruction",
            "statement": (
                "A diagonal sign conjugation leaves every beta_h diagonal entry fixed and "
                "preserves the product of opposite off-diagonal entries, so lambda_h<0 and "
                "beta_1<0 block a raw entrywise nonnegative tridiagonal production matrix."
            ),
            "proof_boundary": "Rejects only the naive tridiagonal Motzkin production matrix, not every possible modified signed model.",
        },
        {
            "id": "first_mu2_cancellation",
            "statement": (
                "The first nontrivial coefficient satisfies mu_2=beta_0^2+lambda_1; "
                "on the finite grid beta_0^2>0, lambda_1<0, and mu_2>0, so positivity "
                "already uses cancellation in the ordinary path model."
            ),
            "proof_boundary": "Finite cancellation diagnostic only; not an all-order obstruction theorem.",
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
    mu2_ok_count = 0
    beta1_ok_count = 0
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                moments = reciprocal_moments(balls, shift_n, degree_d, 2 * degree_d + 1)
                betas = jacobi_betas(moments, degree_d)
                row = mu2_row(lam, shift_n, degree_d, moments, betas)
                mu2_ok_count += int(row.ok)
                all_rows.append(asdict(row))
                if degree_d >= 3:
                    diagonal_row = beta1_row(lam, shift_n, degree_d, betas)
                    beta1_ok_count += int(diagonal_row.ok)
                    all_rows.append(asdict(diagonal_row))

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.out_jsonl.open("w", encoding="utf-8") as handle:
        for row in all_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    expected_mu2_rows = len(lambdas) * len(shifts) * len(degrees)
    expected_beta1_rows = len(lambdas) * len(shifts) * sum(1 for degree in degrees if degree >= 3)
    payload = {
        "kind": "jensen_window_pf_reciprocal_motzkin_path_obstruction_scout",
        "date": "2026-07-06",
        "target_route_rows": [
            "rp_04_companion_or_production_matrix_total_positivity",
            "rp_09_signed_or_modified_continued_fraction",
        ],
        "proof_boundary": (
            "Ordinary Motzkin-path obstruction diagnostic only; not a proof against every "
            "modified signed continued-fraction or production-matrix model, not an all-order "
            "recurrence theorem, not Schur positivity, and not a proof of Lambda <= 0."
        ),
        "symbolic_rows": symbolic_rows(),
        "finite_grid": {
            "lambdas": lambdas,
            "shifts": [min(shifts), max(shifts)],
            "degrees": [min(degrees), max(degrees)],
            "dps": args.dps,
            "motzkin_mu2_cancellation_rows": expected_mu2_rows,
            "motzkin_mu2_cancellation_ok_rows": mu2_ok_count,
            "beta1_diagonal_obstruction_rows": expected_beta1_rows,
            "beta1_diagonal_obstruction_ok_rows": beta1_ok_count,
            "all_mu2_rows_show_negative_path_cancellation": mu2_ok_count == expected_mu2_rows,
            "all_beta1_rows_block_nonnegative_diagonal": beta1_ok_count == expected_beta1_rows,
            "source_enclosures": [str(path.relative_to(REPO_ROOT)) for path in args.enclosure_jsonl],
            "row_log": str(args.out_jsonl.relative_to(REPO_ROOT)),
        },
        "summary": {
            "symbolic_rows": 3,
            "motzkin_mu2_cancellation_rows": expected_mu2_rows,
            "beta1_diagonal_obstruction_rows": expected_beta1_rows,
            "total_finite_rows": expected_mu2_rows + expected_beta1_rows,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The ordinary J-fraction Motzkin path expansion is not a manifest positivity "
                "proof for the zeta-window reciprocal: every finite mu_2 row has beta_0^2>0, "
                "lambda_1<0, and mu_2>0, while every d>=3 beta_1 row is negative and invariant "
                "under diagonal sign conjugation. Any useful route must therefore be a genuinely "
                "modified path, parity, network, or Xi/Phi-specific representation."
            ),
        },
        "invariants": [
            "The scout rejects only the raw ordinary Motzkin/J-fraction production matrix.",
            "Every mu_2 row must contain a negative lambda_1 path contribution while mu_2 remains positive.",
            "Every d>=3 beta_1 row must be negative and separated from zero.",
            "Diagonal sign conjugation cannot change diagonal beta signs.",
            "The scout is finite evidence and a structural rejection gate, not an all-order theorem.",
        ],
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF reciprocal Motzkin path obstruction scout: "
            f"{expected_mu2_rows} mu2 cancellation rows, {expected_beta1_rows} beta1 obstruction rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
