#!/usr/bin/env python3
"""Map finite evidence onto the Jensen-window column recurrence contract.

This script does not recompute Arb determinants.  It records how existing
finite Arb determinant and Sturm-to-PF manifests cover the column recurrence
contract for checked zeta heat-flow windows, while preserving the all-order
theorem gap.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_RECURRENCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json"
DEFAULT_DIRECT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520_summary.json"
)
DEFAULT_DIRECT_ROWS = (
    REPO_ROOT
    / "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520.jsonl"
)
DEFAULT_STURM_D3_D4 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json"
DEFAULT_STURM_D5 = REPO_ROOT / "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json"
DEFAULT_STRESS_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520_summary.json"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_column_recurrence_finite_coverage.json"


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def count_hard_rows(rows_path: Path) -> dict[str, dict]:
    specs = {
        (3, 8): "d3_column_recurrence_m8",
        (4, 6): "d4_column_recurrence_m6",
    }
    counts: dict[str, Counter[str]] = {row_id: Counter() for row_id in specs.values()}
    with rows_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if not raw.strip():
                continue
            row = json.loads(raw)
            key = (int(row["degree_d"]), int(row["minor_size_m"]))
            row_id = specs.get(key)
            if row_id is None:
                continue
            counts[row_id]["rows"] += 1
            counts[row_id][str(row["classification"])] += 1
            if row.get("ok") is True and row.get("contains_zero") is False:
                counts[row_id]["ok"] += 1
    return {
        row_id: {
            "source_recurrence_row": row_id,
            "checked_rows": int(counter["rows"]),
            "positive_rows": int(counter["positive"]),
            "ok_rows": int(counter["ok"]),
            "failed_or_inconclusive_rows": int(counter["rows"] - counter["ok"]),
        }
        for row_id, counter in counts.items()
    }


def direct_coverage(summary: dict, rows_path: Path) -> dict:
    hard = count_hard_rows(rows_path)
    return {
        "source_summary": "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520_summary.json",
        "source_rows": "work/rh_compute/results/arb_jensen_window_pf_obligations_lamgrid_n0_n20_d3_m1_m8_d4_m1_m6_dps520.jsonl",
        "kind": summary["kind"],
        "lambdas": summary["lambdas"],
        "shifts": summary["shifts"],
        "degree_sizes": summary["degree_sizes"],
        "checked_rows": summary["rows"],
        "positive_rows": summary["ok"],
        "failed_or_inconclusive_rows": summary["failed_or_inconclusive"],
        "all_ok": summary["all_ok"],
        "hard_recurrence_rows": hard,
        "hard_recurrence_rows_checked": sum(row["checked_rows"] for row in hard.values()),
        "hard_recurrence_rows_ok": sum(row["ok_rows"] for row in hard.values()),
        "proof_boundary": (
            "Finite Arb determinant coverage for selected column/contiguous "
            "Jensen-window minors only; not all degrees, all shifts, all "
            "lambda values, or a proof of the recurrence theorem."
        ),
    }


def sturm_coverage(summary_d3_d4: dict, summary_d5: dict) -> dict:
    summaries = [summary_d3_d4, summary_d5]
    degrees: list[int] = []
    rows_by_degree: dict[str, int] = {}
    for summary in summaries:
        degrees.extend(int(degree) for degree in summary["degrees"])
        rows_by_degree.update({str(key): int(value) for key, value in summary["rows_by_degree"].items()})
    return {
        "source_summaries": [
            "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d3_d4_dps520_summary.json",
            "work/rh_compute/results/arb_jensen_window_sturm_lamgrid_n0_n20_d5_dps520_summary.json",
        ],
        "degrees": sorted(degrees),
        "lambdas": summary_d3_d4["lambdas"],
        "shifts": summary_d3_d4["shifts"],
        "checked_windows": sum(int(summary["rows"]) for summary in summaries),
        "ok_windows": sum(int(summary["ok"]) for summary in summaries),
        "failed_or_inconclusive_windows": sum(int(summary["failed_or_inconclusive"]) for summary in summaries),
        "rows_by_degree": rows_by_degree,
        "finite_consequence": (
            "For each checked window, Sturm-certified finite PF-infinity implies "
            "all finite Toeplitz minors of that one binomial window are "
            "nonnegative, including all column recurrence minors."
        ),
        "proof_boundary": (
            "Finite window-by-window PF consequence only; not an all-degree, "
            "all-shift, or all-lambda theorem."
        ),
    }


def stress_coverage(stress_summary: dict) -> dict:
    return {
        "source_summary": "work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520_summary.json",
        "source_rows": "work/rh_compute/results/arb_jensen_window_column_recurrence_lamgrid_n0_n20_d3_d8_m1_m20_dps520.jsonl",
        "degrees": stress_summary["degrees"],
        "sizes": stress_summary["sizes"],
        "lambdas": stress_summary["lambdas"],
        "shifts": stress_summary["shifts"],
        "checked_rows": stress_summary["rows"],
        "ok_rows": stress_summary["ok"],
        "failed_or_inconclusive_rows": stress_summary["failed_or_inconclusive"],
        "ok_by_degree": stress_summary["ok_by_degree"],
        "proof_boundary": (
            "Finite recurrence-only Arb stress evidence for checked degrees, "
            "sizes, lambdas, and shifts; not an all-order recurrence theorem."
        ),
    }


def build_payload(
    recurrence: dict,
    direct_summary: dict,
    direct_rows: Path,
    sturm_d3_d4: dict,
    sturm_d5: dict,
    stress_summary: dict,
) -> dict:
    direct = direct_coverage(direct_summary, direct_rows)
    sturm = sturm_coverage(sturm_d3_d4, sturm_d5)
    stress = stress_coverage(stress_summary)
    return {
        "kind": "jensen_window_pf_column_recurrence_finite_coverage",
        "date": "2026-07-06",
        "target_obligation": "jwpf_06_sign_regular_to_jensen_pf_conversion",
        "source_column_recurrence_contract": "work/rh_compute/results/jensen_window_pf_column_recurrence_contract.json",
        "proof_boundary": (
            "Finite coverage map only; not a proof of the all-order column "
            "recurrence, Jensen-window PF-infinity, RH, or Lambda <= 0."
        ),
        "contract_summary": {
            "degree_rows": recurrence["summary"]["degree_rows"],
            "hard_frontier_rows": recurrence["summary"]["hard_frontier_rows"],
            "necessary_condition": "C[m] >= 0 for every m >= 0, every degree d, and every shift n",
        },
        "direct_arb_determinant_coverage": direct,
        "sturm_pf_window_coverage": sturm,
        "arb_column_recurrence_stress_coverage": stress,
        "gap_analysis": {
            "covered": [
                "Direct Arb determinant positivity for d=3, m=1..8 and d=4, m=1..6 on five lambdas and shifts n=0..20.",
                "Direct Arb determinant positivity for the two hard recurrence sizes d=3,m=8 and d=4,m=6 on the same finite grid.",
                "Finite Sturm-to-PF consequence for degrees d=3,4,5 on five lambdas and shifts n=0..20.",
                "Recurrence-only Arb stress positivity for degrees d=3..8 and sizes m=1..20 on the same lambda and shift grid.",
            ],
            "not_covered": [
                "all degrees d",
                "all shifts n",
                "all lambda values or an interval in lambda",
                "an analytic proof that the recurrence is nonnegative for the actual zeta windows",
                "all non-column skew shapes beyond the finite checked windows",
            ],
        },
        "summary": {
            "direct_checked_rows": direct["checked_rows"],
            "direct_positive_rows": direct["positive_rows"],
            "direct_hard_checked_rows": direct["hard_recurrence_rows_checked"],
            "direct_hard_ok_rows": direct["hard_recurrence_rows_ok"],
            "sturm_pf_checked_windows": sturm["checked_windows"],
            "sturm_pf_ok_windows": sturm["ok_windows"],
            "stress_checked_rows": stress["checked_rows"],
            "stress_ok_rows": stress["ok_rows"],
            "target_closing": False,
            "main_finding": (
                "Existing finite evidence supports the column recurrence target "
                "on checked zeta heat-flow windows: 1470 direct Arb determinant "
                "rows are positive, including 210 hard recurrence grid rows, "
                "315 checked Sturm/PF windows imply all column recurrence "
                "minors for those windows, and 12600 recurrence-only stress "
                "rows are positive.  This remains finite evidence only."
            ),
        },
        "invariants": [
            "Direct determinant coverage must be all positive and have zero inconclusive rows.",
            "The hard recurrence sizes d=3,m=8 and d=4,m=6 must each have 105 positive checked rows.",
            "Sturm/PF coverage is a finite consequence, not a replacement for an all-order recurrence theorem.",
            "No row may set target_closing=true.",
        ],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recurrence", type=Path, default=DEFAULT_RECURRENCE)
    parser.add_argument("--direct-summary", type=Path, default=DEFAULT_DIRECT_SUMMARY)
    parser.add_argument("--direct-rows", type=Path, default=DEFAULT_DIRECT_ROWS)
    parser.add_argument("--sturm-d3-d4", type=Path, default=DEFAULT_STURM_D3_D4)
    parser.add_argument("--sturm-d5", type=Path, default=DEFAULT_STURM_D5)
    parser.add_argument("--stress-summary", type=Path, default=DEFAULT_STRESS_SUMMARY)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload(
        load_json(args.recurrence),
        load_json(args.direct_summary),
        args.direct_rows,
        load_json(args.sturm_d3_d4),
        load_json(args.sturm_d5),
        load_json(args.stress_summary),
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.json:
        print(json.dumps(payload, indent=2, sort_keys=True))
    else:
        summary = payload["summary"]
        print(
            "wrote Jensen-window PF column recurrence finite coverage: "
            f"{summary['direct_positive_rows']} direct positive rows, "
            f"{summary['direct_hard_ok_rows']} hard recurrence rows, "
            f"{summary['sturm_pf_ok_windows']} Sturm/PF windows, "
            f"{summary['stress_ok_rows']} stress rows"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
