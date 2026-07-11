#!/usr/bin/env python3
"""Scout signed Jacobi diagonal parameters for Jensen-window reciprocals.

For E(t)=1/H(-t), let mu_m=[t^m]E(t).  The ordinary J-fraction

    sum_m mu_m t^m =
    1 / (1 - beta_0 t - lambda_1 t^2/(1 - beta_1 t - ...))

uses diagonal parameters beta_n.  These are not computed by the fully shifted
Hankel determinants det(mu_{i+j+1}).  Instead, if

    Delta_r = det(mu_{i+j})_{0<=i,j<r}
    Delta_r^* = det with the final Hankel column shifted to mu_{i+r},
    Q_r = Delta_r^*/Delta_r, Q_0=0,

then beta_n = Q_{n+1}-Q_n.  This script checks the finite beta signature on
the current Arb Jensen-window grid.  It is theorem-search evidence only.
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
    / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_OUT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_signed_jacobi_beta_scout.json"
)


@dataclass(frozen=True)
class JacobiBetaRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    beta_index: int
    expected_signature: str
    beta_classification: str
    beta_contains_zero: bool
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
    """Hankel determinant with only the final column shifted one step."""
    return flint.arb_mat(
        [
            [moments[row + col] if col < order - 1 else moments[row + order] for col in range(order)]
            for row in range(order)
        ]
    ).det()


def expected_signature(degree_d: int, beta_index: int) -> str:
    if degree_d == 2 and beta_index == 1:
        return "terminal_degree2_zero"
    if beta_index == 1:
        return "negative_second"
    return "positive"


def beta_ok(signature: str, classification: str, contains_zero: bool) -> bool:
    if signature == "terminal_degree2_zero":
        return contains_zero and classification in {"zero", "inconclusive_contains_zero"}
    if signature == "negative_second":
        return classification == "negative" and not contains_zero
    return classification == "positive" and not contains_zero


def jacobi_beta_rows(
    lam: str,
    shift_n: int,
    degree_d: int,
    moments: list[flint.arb],
) -> list[JacobiBetaRow]:
    deltas = [flint.arb(1)] + [hankel_det(moments, order) for order in range(1, degree_d + 1)]
    star_deltas = [flint.arb(0)] + [star_hankel_det(moments, order) for order in range(1, degree_d + 1)]
    q_values = [flint.arb(0)] + [
        star_deltas[order] / deltas[order] for order in range(1, degree_d + 1)
    ]

    rows = []
    for beta_index in range(degree_d):
        beta_value = q_values[beta_index + 1] - q_values[beta_index]
        signature = expected_signature(degree_d, beta_index)
        classification = classify(beta_value)
        contains_zero = beta_value.contains(0)
        rows.append(
            JacobiBetaRow(
                kind="jensen_window_pf_reciprocal_signed_jacobi_beta_row",
                lam=lam,
                shift_n=shift_n,
                degree_d=degree_d,
                beta_index=beta_index,
                expected_signature=signature,
                beta_classification=classification,
                beta_contains_zero=contains_zero,
                ok=beta_ok(signature, classification, contains_zero),
            )
        )
    return rows


def symbolic_rows() -> list[dict]:
    return [
        {
            "id": "last_column_shifted_beta_formula",
            "statement": (
                "For the ordinary J-fraction diagonal parameter, define "
                "Delta_r^* by shifting only the final Hankel column, set "
                "Q_r=Delta_r^*/Delta_r and Q_0=0, then beta_n=Q_{n+1}-Q_n."
            ),
            "formula_boundary": (
                "This is not the fully shifted determinant det(mu_{i+j+1}); using the "
                "fully shifted determinant gives the wrong degree-2 sanity check."
            ),
        },
        {
            "id": "degree_two_formula_sanity_check",
            "statement": (
                "For H(t)=1+g_1*t+g_2*t^2 and E(t)=1/(1-g_1*t+g_2*t^2), "
                "the formula gives beta_0=g_1, lambda_1=-g_2, and terminal beta_1=0."
            ),
            "proof_boundary": "Sanity check only; not an all-degree theorem for zeta windows.",
        },
        {
            "id": "finite_beta_signature",
            "statement": (
                "On the current finite grid the diagonal beta signature is beta_0>0, "
                "beta_1<0 for d>=3, beta_n>0 for n>=2, and the degree-2 terminal beta_1 "
                "contains zero."
            ),
            "proof_boundary": "Finite Arb evidence only; not a signed Jacobi or production-matrix theorem.",
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
    positive_count = 0
    negative_count = 0
    terminal_zero_count = 0
    ok_count = 0
    for lam in lambdas:
        balls = load_enclosure_balls(args.enclosure_jsonl, lam, needed_max_k)
        for shift_n in shifts:
            for degree_d in degrees:
                moments = reciprocal_moments(balls, shift_n, degree_d, 2 * degree_d + 1)
                rows = jacobi_beta_rows(lam, shift_n, degree_d, moments)
                for row in rows:
                    if row.beta_classification == "positive":
                        positive_count += 1
                    if row.beta_classification == "negative":
                        negative_count += 1
                    if row.expected_signature == "terminal_degree2_zero" and row.beta_contains_zero:
                        terminal_zero_count += 1
                    if row.ok:
                        ok_count += 1
                all_rows.extend(asdict(row) for row in rows)

    args.out_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with args.out_jsonl.open("w", encoding="utf-8") as handle:
        for row in all_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    expected_beta_rows = len(lambdas) * len(shifts) * sum(degrees)
    expected_negative_rows = len(lambdas) * len(shifts) * sum(1 for degree in degrees if degree >= 3)
    expected_terminal_zero_rows = len(lambdas) * len(shifts) if 2 in degrees else 0
    expected_positive_rows = expected_beta_rows - expected_negative_rows - expected_terminal_zero_rows
    payload = {
        "kind": "jensen_window_pf_reciprocal_signed_jacobi_beta_scout",
        "date": "2026-07-06",
        "target_route_row": "rp_09_signed_or_modified_continued_fraction",
        "proof_boundary": (
            "Signed Jacobi diagonal-parameter diagnostic only; not a proof of a signed "
            "continued-fraction theorem, a production-matrix model, the all-order column "
            "recurrence, Schur positivity, Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
        "symbolic_rows": symbolic_rows(),
        "finite_grid": {
            "lambdas": lambdas,
            "shifts": [min(shifts), max(shifts)],
            "degrees": [min(degrees), max(degrees)],
            "dps": args.dps,
            "beta_rows": expected_beta_rows,
            "ok_beta_rows": ok_count,
            "positive_beta_rows": positive_count,
            "negative_beta_rows": negative_count,
            "terminal_zero_beta_rows": terminal_zero_count,
            "expected_positive_beta_rows": expected_positive_rows,
            "expected_negative_beta_rows": expected_negative_rows,
            "expected_terminal_zero_beta_rows": expected_terminal_zero_rows,
            "all_beta_rows_match_signature": ok_count == expected_beta_rows,
            "source_enclosures": [str(path.relative_to(REPO_ROOT)) for path in args.enclosure_jsonl],
            "row_log": str(args.out_jsonl.relative_to(REPO_ROOT)),
        },
        "summary": {
            "symbolic_rows": 3,
            "beta_rows": expected_beta_rows,
            "positive_beta_rows": positive_count,
            "negative_beta_rows": negative_count,
            "terminal_zero_beta_rows": terminal_zero_count,
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "On the finite Arb grid, the Jacobi diagonal parameters for E(t)=1/H(-t) "
                "match beta_0>0, beta_1<0 for degrees d>=3, beta_n>0 for n>=2, and a "
                "terminal zero-containing beta_1 in degree d=2. This gives a sharper "
                "signed-Jacobi target but does not supply the missing all-order theorem."
            ),
        },
        "invariants": [
            "The beta formula uses a last-column shifted determinant, not a fully shifted Hankel determinant.",
            "Every nonterminal positive beta row in the finite grid must be positive and separated from zero.",
            "Every d>=3 beta_1 row in the finite grid must be negative and separated from zero.",
            "Every degree-2 terminal beta_1 row must contain zero.",
            "The scout is finite evidence only and does not mark the signed-fraction route ready to apply.",
        ],
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF reciprocal signed Jacobi beta scout: "
            f"{expected_beta_rows} beta rows, {positive_count} positive, "
            f"{negative_count} negative, {terminal_zero_count} terminal-zero rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
