#!/usr/bin/env python3
"""Scout a cheap state-space sign-lift obstruction for the signed model route.

The modified signed-model target leaves a genuine state-space doubled model
open.  This script rejects only a much cheaper candidate: keep the raw
ordinary Motzkin paths, split signs into states, replace each raw signed
weight by its absolute value, and count all lifted paths with nonnegative
readout.

At the first nontrivial coefficient this cannot be exact.  The ordinary
J-fraction path model gives

    mu_2 = beta_0^2 + lambda_1 = beta_0^2 - kappa_1

with kappa_1=-lambda_1>0 on the validated grid.  An absolute-value sign lift
would contribute

    beta_0^2 + kappa_1,

so the gap is 2*kappa_1>0.  A viable state-space model must therefore change
the representation in a deeper way; it cannot be only an absolute-value cover
of the raw Motzkin paths.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_SOURCE_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_OUT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_lamgrid_n0_n20_d2_d8_dps520.jsonl"
)
DEFAULT_OUT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_state_space_sign_lift_obstruction_scout.json"
)


@dataclass(frozen=True)
class SignLiftRow:
    kind: str
    lam: str
    shift_n: int
    degree_d: int
    source_kind: str
    raw_mu2_formula: str
    absolute_sign_lift_formula: str
    gap_formula: str
    lambda1_classification: str
    kappa1_classification: str
    mu2_classification: str
    absolute_lift_gap_classification: str
    absolute_lift_matches_mu2: bool
    ok: bool


def load_mu2_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("kind") == "jensen_window_pf_reciprocal_motzkin_mu2_cancellation_row":
                rows.append(row)
    return rows


def derived_row(source: dict) -> SignLiftRow:
    lambda_negative = (
        source.get("lambda1_classification") == "negative"
        and source.get("lambda1_contains_zero") is False
    )
    mu2_positive = (
        source.get("mu2_classification") == "positive"
        and source.get("mu2_contains_zero") is False
    )
    beta_square_positive = (
        source.get("beta0_square_classification") == "positive"
        and source.get("beta0_square_contains_zero") is False
    )
    ok = bool(lambda_negative and mu2_positive and beta_square_positive and source.get("ok") is True)
    return SignLiftRow(
        kind="jensen_window_pf_state_space_sign_lift_mu2_obstruction_row",
        lam=str(source["lam"]),
        shift_n=int(source["shift_n"]),
        degree_d=int(source["degree_d"]),
        source_kind=str(source["kind"]),
        raw_mu2_formula="mu_2=beta_0^2+lambda_1=beta_0^2-kappa_1",
        absolute_sign_lift_formula="beta_0^2+kappa_1",
        gap_formula="absolute_lift_mu2-mu_2=2*kappa_1=-2*lambda_1",
        lambda1_classification=str(source.get("lambda1_classification")),
        kappa1_classification="positive" if lambda_negative else "unknown",
        mu2_classification=str(source.get("mu2_classification")),
        absolute_lift_gap_classification="positive" if lambda_negative else "unknown",
        absolute_lift_matches_mu2=False if lambda_negative else True,
        ok=ok,
    )


def symbolic_rows() -> list[dict]:
    return [
        {
            "id": "absolute_value_sign_state_cover",
            "statement": (
                "A cheap two-state sign cover keeps the raw Motzkin paths, replaces "
                "negative raw path weights by their absolute values, and counts the "
                "lifted paths with nonnegative transition and readout weights."
            ),
            "proof_boundary": "This is only one naive state-space lift, not every possible state-space model.",
        },
        {
            "id": "mu2_absolute_lift_gap",
            "statement": (
                "For mu_2, the raw signed expression is beta_0^2+lambda_1=beta_0^2-kappa_1, "
                "while the absolute-value sign lift gives beta_0^2+kappa_1, so the gap is "
                "2*kappa_1=-2*lambda_1."
            ),
            "proof_boundary": "Uses the validated finite lambda_1<0 rows as a derived obstruction; not an all-order theorem.",
        },
        {
            "id": "surviving_state_space_requirement",
            "statement": (
                "A surviving state-space doubled model must change the representation, "
                "transition weights, or coefficient extraction so the final formula is "
                "exactly E(t)=1/H(-t) and manifestly nonnegative; it cannot be merely "
                "an absolute-value cover of the raw Motzkin paths."
            ),
            "proof_boundary": "Requirement for future theorem search only; not a construction.",
        },
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-rows", type=Path, default=DEFAULT_SOURCE_ROWS)
    parser.add_argument("--out-rows", type=Path, default=DEFAULT_OUT_ROWS)
    parser.add_argument("--out-summary", type=Path, default=DEFAULT_OUT_SUMMARY)
    parser.add_argument("--json", action="store_true", help="Print summary JSON.")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    source_rows = load_mu2_rows(args.source_rows)
    derived = [asdict(derived_row(row)) for row in source_rows]

    args.out_rows.parent.mkdir(parents=True, exist_ok=True)
    with args.out_rows.open("w", encoding="utf-8") as handle:
        for row in derived:
            handle.write(json.dumps(row, sort_keys=True) + "\n")

    ok_rows = sum(1 for row in derived if row.get("ok") is True)
    payload = {
        "kind": "jensen_window_pf_state_space_sign_lift_obstruction_scout",
        "date": "2026-07-06",
        "target_model_row": "msm_04_state_space_doubled_model",
        "target_route_rows": [
            "rp_04_companion_or_production_matrix_total_positivity",
            "rp_09_signed_or_modified_continued_fraction",
        ],
        "source_artifacts": [
            "outputs/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.md",
            "work/rh_compute/results/jensen_window_pf_reciprocal_motzkin_path_obstruction_scout.json",
            str(args.source_rows.relative_to(REPO_ROOT)),
        ],
        "proof_boundary": (
            "State-space sign-lift obstruction diagnostic only; rejects only the "
            "absolute-value cover of raw Motzkin path weights, not every state-space "
            "doubled model, not a production-matrix theorem, not an all-order recurrence "
            "theorem, not Schur positivity, and not Lambda <= 0."
        ),
        "symbolic_rows": symbolic_rows(),
        "finite_grid": {
            "derived_from": str(args.source_rows.relative_to(REPO_ROOT)),
            "row_log": str(args.out_rows.relative_to(REPO_ROOT)),
            "lambdas": ["0", "1e-6", "1e-4", "1e-2", "1e-1"],
            "shifts": [0, 20],
            "degrees": [2, 8],
            "mu2_sign_lift_obstruction_rows": len(derived),
            "mu2_sign_lift_obstruction_ok_rows": ok_rows,
            "all_absolute_lift_gaps_positive": ok_rows == len(derived),
            "all_absolute_lift_rows_fail_exact_mu2": all(
                row.get("absolute_lift_matches_mu2") is False for row in derived
            ),
        },
        "summary": {
            "symbolic_rows": 3,
            "mu2_sign_lift_obstruction_rows": len(derived),
            "ready_to_apply_rows": 0,
            "target_closing": False,
            "main_finding": (
                "The cheap state-space sign lift that simply takes absolute values of "
                "the raw Motzkin path weights is rejected at mu_2: every derived finite "
                "row has absolute_lift_mu2-mu_2=2*kappa_1>0. A live state-space doubled "
                "route must therefore alter the representation itself, not just split "
                "raw signed paths into positive and negative states."
            ),
        },
        "invariants": [
            "The scout rejects only the absolute-value cover of raw Motzkin paths.",
            "Every derived row depends on a validated lambda_1<0 Motzkin obstruction row.",
            "No row is ready_to_apply.",
            "The state-space doubled route remains live only if a genuinely modified exact positive model is constructed.",
        ],
    }
    args.out_summary.parent.mkdir(parents=True, exist_ok=True)
    args.out_summary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        print(
            "wrote Jensen-window PF state-space sign-lift obstruction scout: "
            f"{len(derived)} mu2 sign-lift obstruction rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
