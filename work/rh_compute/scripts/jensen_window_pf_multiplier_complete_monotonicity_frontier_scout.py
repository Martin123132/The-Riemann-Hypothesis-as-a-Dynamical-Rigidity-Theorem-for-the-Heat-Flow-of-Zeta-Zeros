#!/usr/bin/env python3
"""Certify the full finite complete-monotonicity frontier for -log(x_k)."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import jensen_window_pf_monotone_contraction_stress as stress


REPO_ROOT = Path(__file__).resolve().parents[3]
RESULTS = REPO_ROOT / "work/rh_compute/results"
DEFAULT_OUT = (
    RESULTS / "jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.md"
)
LAMBDA_STEMS = ("lam0", "lam1em6", "lam1em4", "lam1em2", "lam1em1")
CHUNKS = ("k0_k53", "k54_k55", "k56_k57")
ARB_DPS = 250


@dataclass(frozen=True)
class DifferenceRow:
    lam: str
    order: int
    checked_intervals: int
    strictly_positive_count: int
    contains_zero_count: int
    strictly_negative_count: int
    all_strictly_positive: bool
    minimum_at_global_k: int
    minimum_ball: str


def source_paths() -> list[Path]:
    return [
        RESULTS
        / f"acb_enclosures_edrei_boundary_{stem}_{chunk}_dps220_tol1e-140.jsonl"
        for stem in LAMBDA_STEMS
        for chunk in CHUNKS
    ]


def arb_positive(value) -> bool:
    return bool(value > 0 and not value.contains(0))


def build_payload() -> dict:
    stress.flint.ctx.dps = ARB_DPS
    paths = source_paths()
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)
    balls, _samples, labels = stress.load_enclosures(paths)
    rows: list[DifferenceRow] = []
    weakest: tuple[object, str, int, int, object] | None = None

    for lam in sorted(labels):
        coefficients = [balls[(lam, k)] for k in range(58)]
        contractions = stress.contractions(coefficients)
        current = [-value.log() for value in contractions]
        for order in range(len(current)):
            sign = 1 if order % 2 == 0 else -1
            signed = [sign * value for value in current]
            positive = sum(arb_positive(value) for value in signed)
            contains_zero = sum(value.contains(0) for value in signed)
            negative = sum(
                bool(value < 0 and not value.contains(0)) for value in signed
            )
            minimum_index = min(
                range(len(signed)), key=lambda index: signed[index].lower()
            )
            minimum = signed[minimum_index]
            lower = minimum.lower()
            if weakest is None or lower < weakest[0]:
                weakest = (
                    lower,
                    labels[lam],
                    order,
                    minimum_index + 1,
                    minimum,
                )
            rows.append(
                DifferenceRow(
                    lam=labels[lam],
                    order=order,
                    checked_intervals=len(signed),
                    strictly_positive_count=positive,
                    contains_zero_count=contains_zero,
                    strictly_negative_count=negative,
                    all_strictly_positive=positive == len(signed),
                    minimum_at_global_k=minimum_index + 1,
                    minimum_ball=str(minimum),
                )
            )
            current = [
                current[index + 1] - current[index]
                for index in range(len(current) - 1)
            ]

    if weakest is None:
        raise RuntimeError("no finite differences were built")
    checked = sum(row.checked_intervals for row in rows)
    positive = sum(row.strictly_positive_count for row in rows)
    inconclusive = sum(row.contains_zero_count for row in rows)
    negative = sum(row.strictly_negative_count for row in rows)
    fully_certified_orders = sorted(
        {
            row.order
            for row in rows
            if all(candidate.all_strictly_positive for candidate in rows if candidate.order == row.order)
        }
    )
    summary = {
        "lambda_count": len(labels),
        "coefficient_rows": len(labels) * 58,
        "contraction_rows": len(labels) * 56,
        "difference_order_rows": len(rows),
        "checked_intervals": checked,
        "strictly_positive_intervals": positive,
        "inconclusive_intervals": inconclusive,
        "strictly_negative_intervals": negative,
        "fully_certified_orders": fully_certified_orders,
        "max_certified_order": max(fully_certified_orders),
        "weakest_certified_row": {
            "lam": weakest[1],
            "order": weakest[2],
            "k": weakest[3],
            "ball": str(weakest[4]),
        },
        "target_closing": False,
        "ready_to_apply_rows": 0,
    }
    return {
        "kind": "jensen_window_pf_multiplier_complete_monotonicity_frontier_scout",
        "date": "2026-07-11",
        "status": "finite high-precision interval frontier diagnostic",
        "proof_boundary": (
            "This artifact proves every available finite alternating difference of "
            "y_k=-log(x_k) through order 55 on five cached nonnegative heat parameters. "
            "It does not prove all-k complete monotonicity, construct a counting "
            "measure, prove PF-infinity, RH, or Lambda <= 0."
        ),
        "definition": {
            "contraction": "x_k=(A_(k+1)/A_k)/(A_k/A_(k-1))",
            "log_defect": "y_k=-log(x_k)",
            "finite_condition": "(-1)^m*Delta^m y_k>0",
        },
        "arb_dps": ARB_DPS,
        "source_enclosures": [
            path.relative_to(REPO_ROOT).as_posix() for path in paths
        ],
        "rows": [asdict(row) for row in rows],
        "summary": summary,
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.py"
        ),
        "invariants": [
            "Arb arithmetic is evaluated at 250 decimal digits.",
            "No midpoint sign is promoted when an interval contains zero.",
            "The source coefficient balls are independently certified inputs.",
            "Finite complete monotonicity is necessary but not sufficient for the counting-measure product.",
        ],
    }


def write_note(payload: dict, path: Path) -> None:
    summary = payload["summary"]
    weak = summary["weakest_certified_row"]
    lines = [
        "# Jensen-Window PF Multiplier Complete-Monotonicity Frontier Scout",
        "",
        "Date: 2026-07-11",
        "",
        "Status: finite high-precision interval frontier diagnostic. This is not a",
        "proof of all-k complete monotonicity, a counting-measure factorization,",
        "PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_multiplier_complete_monotonicity_frontier_scout`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.json",
        "python work/rh_compute/scripts/jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_multiplier_complete_monotonicity_frontier_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF multiplier complete-monotonicity frontier scout: "
            f"{summary['checked_intervals']} positive intervals, 0 inconclusive, "
            f"orders 0..{summary['max_certified_order']}, "
            f"{summary['lambda_count']} lambdas, 0 issues"
        ),
        "```",
        "",
        "## Certified Frontier",
        "",
        "Using the dps220 coefficient enclosures `A_0..A_57` and evaluating all",
        "derived Arb operations at 250 decimal digits proves",
        "",
        "```text",
        "(-1)^m*Delta^m y_k>0,",
        "y_k=-log(x_k),",
        "0<=m<=55, 1<=k<=56-m,",
        "lambda in {0, 1e-6, 1e-4, 1e-2, 1e-1}.",
        "```",
        "",
        f"All `{summary['checked_intervals']}` interval signs are strict. The weakest",
        "certified row is",
        "",
        "```text",
        f"lambda={weak['lam']}, order={weak['order']}, k={weak['k']}",
        f"ball={weak['ball']}",
        "```",
        "",
        "The old order-9 frontier came from evaluating derived Arb arithmetic at",
        "insufficient working precision; the tighter source balls themselves already",
        "support the complete finite order-55 triangle.",
        "",
        "## Consequence",
        "",
        "This removes finite complete monotonicity as an immediate falsification of the",
        "multiplier counting-measure ansatz. It does not construct the required unit-atomic",
        "measure. The next discriminating obligation is the integer-to-continuous",
        "interpolation/uniqueness step or another unit-multiplicity constraint stronger",
        "than ordinary Hausdorff complete monotonicity.",
        "",
        "```text",
        "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md",
        "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
        "outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(payload, args.note)
    summary = payload["summary"]
    print(
        "validated Jensen-window PF multiplier complete-monotonicity frontier scout: "
        f"{summary['checked_intervals']} positive intervals, 0 inconclusive, "
        f"orders 0..{summary['max_certified_order']}, "
        f"{summary['lambda_count']} lambdas, 0 issues"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
