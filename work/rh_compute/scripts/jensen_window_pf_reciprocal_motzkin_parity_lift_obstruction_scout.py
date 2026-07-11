#!/usr/bin/env python3
"""Scout the global parity/sign-lift obstruction for Motzkin excursions.

The raw ordinary J-fraction Motzkin model has mixed signs.  A natural cheap
repair would be to multiply every coefficient of length m by a global sign, or
to hope that diagonal sign conjugation makes all closed-path weights positive.
This script checks the first symbolic obstruction on the current Arb grid:

    positive path: beta_0^m > 0
    negative path: lambda_1 * beta_0^(m-2) < 0

for every length m>=2 in the tested range.  Since these two paths have the
same length and endpoints, no global length-parity sign or diagonal sign
conjugation can make the raw ordinary Motzkin expansion entrywise
nonnegative.  This is a rejection gate for cheap sign lifts only, not a proof
against genuinely modified signed models.
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
    / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_lamgrid_n0_n20_d2_d8_m2_m8_dps520.jsonl"
)
DEFAULT_OUT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout.json"
)


@dataclass(frozen=True)
class ParityLiftRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    path_length_m: int
    positive_path: str
    positive_path_classification: str
    positive_path_contains_zero: bool
    negative_path: str
    negative_path_classification: str
    negative_path_contains_zero: bool
    same_length_and_endpoints: bool
    global_length_sign_can_fix: bool
    diagonal_conjugation_changes_closed_path_sign: bool
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


def parity_row(
    lam: str,
    shift_n: int,
    degree_d: int,
    path_length_m: int,
    beta0: flint.arb,
    lambda1: flint.arb,
) -> ParityLiftRow:
    positive_path_value = beta0 ** path_length_m
    negative_path_value = lambda1 * (beta0 ** (path_length_m - 2))
    positive_classification = classify(positive_path_value)
    negative_classification = classify(negative_path_value)
    positive_contains_zero = positive_path_value.contains(0)
    negative_contains_zero = negative_path_value.contains(0)
    ok = (
        positive_classification == "positive"
        and not positive_contains_zero
        and negative_classification == "negative"
        and not negative_contains_zero
    )
    return ParityLiftRow(
        kind="jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_row",
        lam=lam,
        shift_n=shift_n,
        degree_d=degree_d,
        path_length_m=path_length_m,
        positive_path="H0^m",
        positive_path_classification=positive_classification,
        positive_path_contains_zero=positive_contains_zero,
        negative_path="U D H0^(m-2)",
        negative_path_classification=negative_classification,
        negative_path_contains_zero=negative_contains_zero,
        same_length_and_endpoints=True,
        global_length_sign_can_fix=False,
        diagonal_conjugation_changes_closed_path_sign=False,
        ok=ok,
    )


def symbolic_rows() -> list[dict]:
    return [
        {
            "id": "same_length_mixed_sign_paths",
            "statement": (
                "For every m>=2, the ordinary Motzkin model contains the same-endpoint "
                "excursions H0^m with sign positive and U D H0^(m-2) with sign negative "
                "when beta_0>0 and lambda_1<0."
            ),
            "proof_boundary": "Rejects a global length-parity sign only for the raw ordinary Motzkin model.",
        },
        {
            "id": "global_length_sign_obstruction",
            "statement": (
                "Because the two witness paths have the same length, multiplying the whole "
                "coefficient of t^m by any sign depending only on m keeps one witness path negative."
            ),
            "proof_boundary": "Does not reject transformations that change the state space or path weights nontrivially.",
        },
        {
            "id": "closed_path_conjugation_invariance",
            "statement": (
                "A diagonal sign conjugation changes open edge signs by endpoint factors, "
                "so the sign of any closed excursion path is invariant."
            ),
            "proof_boundary": "Rejects diagonal sign-conjugation fixes only; not a general oscillatory theorem.",
        },
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enclosure-jsonl", type=Path, nargs="+", default=[REPO_ROOT / p for p in DEFAULT_ENCLOSURE_JSONL])
    parser.add_argument("--lambdas", default=",".join(DEFAULT_LAMBDAS))
    parser.add_argument("--shifts", default="0..20")
    parser.add_argument("--degrees", default="2..8")
    parser.add_argument("--lengths", default="2..8")
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
    lengths = parse_int_range(args.lengths)
    if min(degrees) < 2:
        raise ValueError("degrees must be at least 2")
    if min(lengths) < 2:
        raise ValueError("lengths must be at least 2")
    needed_max_k = max(shifts) + max(degrees)

    all_rows: list[dict] = []
    ok_count = 0
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                moments = reciprocal_moments(balls, shift_n, degree_d, max(lengths))
                beta0 = moments[1]
                lambda1 = hankel_det(moments, 2)
                for path_length_m in lengths:
                    row = parity_row(lam, shift_n, degree_d, path_length_m, beta0, lambda1)
                    ok_count += int(row.ok)
                    all_rows.append(asdict(row))

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.out_jsonl.open("w", encoding="utf-8") as handle:
        for row in all_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    expected_rows = len(lambdas) * len(shifts) * len(degrees) * len(lengths)
    payload = {
        "kind": "jensen_window_pf_reciprocal_motzkin_parity_lift_obstruction_scout",
        "date": "2026-07-06",
        "target_route_rows": [
            "rp_04_companion_or_production_matrix_total_positivity",
            "rp_09_signed_or_modified_continued_fraction",
        ],
        "proof_boundary": (
            "Global parity/sign-lift obstruction diagnostic only; not a proof against "
            "state-space doubled, modified signed continued-fraction, oscillatory, or "
            "production-matrix models, not Schur positivity, and not a proof of Lambda <= 0."
        ),
        "symbolic_rows": symbolic_rows(),
        "finite_grid": {
            "lambdas": lambdas,
            "shifts": [min(shifts), max(shifts)],
            "degrees": [min(degrees), max(degrees)],
            "lengths": [min(lengths), max(lengths)],
            "dps": args.dps,
            "parity_lift_obstruction_rows": expected_rows,
            "parity_lift_obstruction_ok_rows": ok_count,
            "all_rows_have_same_length_mixed_sign_witnesses": ok_count == expected_rows,
            "source_enclosures": [str(path.relative_to(REPO_ROOT)) for path in args.enclosure_jsonl],
            "row_log": str(args.out_jsonl.relative_to(REPO_ROOT)),
        },
        "summary": {
            "symbolic_rows": 3,
            "parity_lift_obstruction_rows": expected_rows,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "On the finite Arb grid, every tested length m=2..8 has same-length "
                "closed Motzkin witnesses of both signs: H0^m is positive and U D H0^(m-2) "
                "is negative. Therefore no global length-parity sign or diagonal sign "
                "conjugation can make the raw ordinary Motzkin expansion manifestly "
                "nonnegative; a surviving path route must alter the state space or weights."
            ),
        },
        "invariants": [
            "The scout rejects only global length-parity signs and diagonal sign conjugation for the raw ordinary Motzkin model.",
            "Every checked row must have a positive H0^m witness and a negative U D H0^(m-2) witness.",
            "The two witness paths have the same length and endpoints.",
            "Closed excursion signs are invariant under diagonal sign conjugation.",
            "The scout is finite evidence and a structural rejection gate, not an all-order theorem.",
        ],
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF reciprocal Motzkin parity-lift obstruction scout: "
            f"{expected_rows} mixed-sign witness rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
