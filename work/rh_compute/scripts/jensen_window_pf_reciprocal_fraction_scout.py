#!/usr/bin/env python3
"""Scout standard continued-fraction signs for Jensen-window reciprocals.

For the column recurrence contract we study

    E(t) = 1 / H(-t),
    H(t) = 1 + g1*t + ... + gd*t^d,
    g_j = binom(d,j) A_{n+j} / A_n.

The standard positive S-fraction and J-fraction routes have immediate sign
requirements on their first nontrivial parameters.  This script records the
symbolic obstruction and checks the sign regime on the current finite Arb
zeta-window grid.  It is a theorem-search diagnostic, not a proof of the
all-order recurrence.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
import json
from math import comb
from pathlib import Path
import sys

import sympy as sp


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
    / "work/rh_compute/results/jensen_window_pf_reciprocal_fraction_sign_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_OUT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_fraction_scout.json"
)


@dataclass(frozen=True)
class FractionSignRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    h0_classification: str
    h1_classification: str
    h2_classification: str
    standard_s_a1_classification: str
    standard_s_a2_classification: str
    standard_j_beta0_classification: str
    standard_j_lambda1_classification: str
    standard_positive_fraction_obstructed: bool


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


def symbolic_rows() -> list[dict]:
    t = sp.Symbol("t")
    g1, g2, g3, g4 = sp.symbols("g1 g2 g3 g4", positive=True)
    denominator = 1 - g1 * t + g2 * t**2 - g3 * t**3 + g4 * t**4
    reciprocal = sp.series(1 / denominator, t, 0, 5).removeO()
    mu1 = sp.expand(reciprocal).coeff(t, 1)
    mu2 = sp.expand(reciprocal).coeff(t, 2)
    inverse = denominator
    s_a1 = mu1
    s_tail = sp.series((1 - inverse) / (s_a1 * t), t, 0, 4).removeO()
    s_a2 = sp.expand(s_tail).coeff(t, 1)
    j_beta0 = mu1
    j_lambda1 = sp.simplify(mu2 - mu1**2)
    return [
        {
            "id": "symbolic_standard_s_fraction_first_obstruction",
            "fraction_family": "standard S-fraction",
            "normalization": "F(t)=1/(1-a1*t/(1-a2*t/(1-...)))",
            "reciprocal_series_prefix": sp.sstr(reciprocal),
            "first_parameter": sp.sstr(s_a1),
            "second_parameter": sp.sstr(s_a2),
            "sign_conclusion_under_g1_g2_positive": "a1>0 and a2<0",
            "obstruction": "Standard positive S-fraction coefficients cannot all be nonnegative for positive g1,g2.",
        },
        {
            "id": "symbolic_standard_j_fraction_first_obstruction",
            "fraction_family": "standard J-fraction",
            "normalization": "F(t)=1/(1-beta0*t-lambda1*t^2*F1(t)), F1(0)=1",
            "reciprocal_series_prefix": sp.sstr(reciprocal),
            "first_beta": sp.sstr(j_beta0),
            "first_lambda": sp.sstr(j_lambda1),
            "sign_conclusion_under_g2_positive": "beta0>0 and lambda1<0",
            "obstruction": "Standard positive Jacobi/Stieltjes fraction has the wrong first quadratic sign for positive g2.",
        },
        {
            "id": "signed_or_modified_fraction_remaining_candidate",
            "fraction_family": "signed or modified continued fraction",
            "normalization": "Would need a theorem compatible with negative ordinary lambda1 or alternating signs.",
            "first_signed_lambda": sp.sstr(-j_lambda1),
            "sign_conclusion_under_g2_positive": "-lambda1=g2>0",
            "obstruction": "A modified signed fraction remains possible only with a separate total-positivity or oscillatory theorem.",
        },
    ]


def sign_row(balls: dict[int, flint.arb], lam: str, shift_n: int, degree_d: int) -> FractionSignRow:
    h0 = balls[shift_n]
    h1 = degree_d * balls[shift_n + 1]
    h2 = comb(degree_d, 2) * balls[shift_n + 2]
    standard_s_a1 = h1 / h0
    standard_s_a2 = -h2 / h1
    standard_j_beta0 = standard_s_a1
    standard_j_lambda1 = -h2 / h0
    return FractionSignRow(
        kind="jensen_window_pf_reciprocal_fraction_sign_row",
        lam=lam,
        shift_n=shift_n,
        degree_d=degree_d,
        h0_classification=classify(h0),
        h1_classification=classify(h1),
        h2_classification=classify(h2),
        standard_s_a1_classification=classify(standard_s_a1),
        standard_s_a2_classification=classify(standard_s_a2),
        standard_j_beta0_classification=classify(standard_j_beta0),
        standard_j_lambda1_classification=classify(standard_j_lambda1),
        standard_positive_fraction_obstructed=(
            classify(h0) == "positive"
            and classify(h1) == "positive"
            and classify(h2) == "positive"
            and classify(standard_s_a1) == "positive"
            and classify(standard_s_a2) == "negative"
            and classify(standard_j_beta0) == "positive"
            and classify(standard_j_lambda1) == "negative"
        ),
    )


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
        raise ValueError("degrees must be at least 2 for the first fraction obstruction")
    needed_max_k = max(shifts) + max(degrees)

    rows: list[FractionSignRow] = []
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                rows.append(sign_row(balls, lam, shift_n, degree_d))

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.out_jsonl.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(asdict(row), sort_keys=True) + "\n")

    obstructed = sum(1 for row in rows if row.standard_positive_fraction_obstructed)
    payload = {
        "kind": "jensen_window_pf_reciprocal_fraction_scout",
        "date": "2026-07-06",
        "target_route_rows": [
            "rp_03_positive_stieltjes_or_j_fraction",
            "rp_04_companion_or_production_matrix_total_positivity",
        ],
        "proof_boundary": (
            "Continued-fraction sign diagnostic only; not a proof of the all-order "
            "column recurrence, Schur positivity, Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
        "symbolic_rows": symbolic_rows(),
        "finite_grid": {
            "lambdas": lambdas,
            "shifts": [min(shifts), max(shifts)],
            "degrees": [min(degrees), max(degrees)],
            "dps": args.dps,
            "rows": len(rows),
            "standard_positive_fraction_obstructed_rows": obstructed,
            "all_rows_obstruct_standard_positive_fraction": obstructed == len(rows),
            "source_enclosures": [str(path.relative_to(REPO_ROOT)) for path in args.enclosure_jsonl],
            "row_log": str(args.out_jsonl.relative_to(REPO_ROOT)),
        },
        "summary": {
            "symbolic_rows": 3,
            "finite_rows": len(rows),
            "standard_positive_fraction_obstructed_rows": obstructed,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The ordinary positive S-fraction and J-fraction routes have the wrong "
                "first nontrivial sign for E(t)=1/H(-t): a2=-g2/g1 and lambda1=-g2. "
                "Thus a continued-fraction proof would have to be a signed or modified "
                "fraction with a separate theorem, not the standard positive Stieltjes/Jacobi route."
            ),
        },
        "invariants": [
            "The symbolic S-fraction row must have second_parameter -g2/g1.",
            "The symbolic J-fraction row must have first_lambda -g2.",
            "Every finite row in the checked grid must obstruct the standard positive fraction route.",
            "The scout does not reject signed or modified continued-fraction routes.",
            "No row may set target_closing=true.",
        ],
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF reciprocal fraction scout: "
            f"{len(rows)} finite rows, {obstructed} standard-positive-fraction obstructions"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
