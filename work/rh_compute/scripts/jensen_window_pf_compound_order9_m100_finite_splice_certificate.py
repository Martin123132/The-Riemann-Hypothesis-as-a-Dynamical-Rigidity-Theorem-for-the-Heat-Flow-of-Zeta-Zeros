#!/usr/bin/env python3
"""Certify the two lambda=-100 order-nine indices before the analytic bridge."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402

import jensen_window_pf_compound_order9_m100_prefix_certificate as prefix  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_m100_finite_splice_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_compound_order9_m100_finite_splice_certificate.md"
)
BRIDGE_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order9_first_summand_curvature_bridge.json"
)
PREFIX_SOURCE = prefix.DEFAULT_OUT
NEW_COEFFICIENT_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order9_k1257_k1258_dps220.jsonl"
)
NEW_COEFFICIENT_SUMMARY = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "negative_lambda_m100_order9_k1257_k1258_dps220_summary.json"
)
SOURCE_PATHS = (*prefix.SOURCE_PATHS, NEW_COEFFICIENT_SOURCE)
PREFIX_LAST_N = 1242
MAX_COEFFICIENT_INDEX = PREFIX_LAST_N + 16
PRECISION_BITS = prefix.PRECISION_BITS


@dataclass(frozen=True)
class SpliceRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_sources() -> dict:
    bridge = load_json(BRIDGE_SOURCE)
    old_prefix = load_json(PREFIX_SOURCE)
    coefficient_summary = load_json(NEW_COEFFICIENT_SUMMARY)
    if bridge.get("exact", {}).get("finite_splice") != (
        "prove Q_(9,n)(-100)>0 for n=1241,1242"
    ):
        raise RuntimeError("order-nine bridge splice contract changed")
    if old_prefix.get("summary", {}).get("positive_Q9_rows") != 1241:
        raise RuntimeError("order-nine prefix source changed")
    if (
        coefficient_summary.get("rows") != 2
        or coefficient_summary.get("k_min") != 1257
        or coefficient_summary.get("k_max") != 1258
        or coefficient_summary.get("lambdas") != ["-100.0"]
        or coefficient_summary.get("n_sum") != 70
        or coefficient_summary.get("cutoff") != "7"
        or coefficient_summary.get("dps") != 220
    ):
        raise RuntimeError("two-row coefficient source changed")
    return {
        "bridge_splice": bridge["exact"]["finite_splice"],
        "old_prefix": old_prefix["exact"]["prefix"],
        "old_prefix_sha256": prefix.sha256(PREFIX_SOURCE),
        "new_coefficients": "A_1257(-100)>0 and A_1258(-100)>0",
        "new_coefficient_sha256": prefix.sha256(NEW_COEFFICIENT_SOURCE),
        "new_coefficient_summary_sha256": prefix.sha256(
            NEW_COEFFICIENT_SUMMARY
        ),
        "new_coefficient_summary": coefficient_summary,
    }


def build_artifact() -> dict:
    flint.ctx.prec = PRECISION_BITS
    sources = validate_sources()
    values, source_diagnostics = prefix.merged_coefficients(
        SOURCE_PATHS, MAX_COEFFICIENT_INDEX
    )
    finite = prefix.finite_prefix(
        values, PREFIX_LAST_N, MAX_COEFFICIENT_INDEX
    )
    splice_rows = finite["rows"][-2:]
    if [row["n"] for row in splice_rows] != [1241, 1242]:
        raise RuntimeError("finite splice rows are not aligned")
    exact = {
        "signed_condensation": (
            "Q_(9,n)*Q_(7,n+2)=Q_(8,n+1)^2-Q_(8,n)*Q_(8,n+2)"
        ),
        "stable_coordinate": (
            "M_n=Q_(8,n+1)^2/(Q_(8,n)*Q_(8,n+2))-1"
        ),
        "splice": "Q_(9,n)(-100)>0 for n=1241,1242",
        "combined_prefix": "Q_(9,n)(-100)>0 for every 0<=n<=1242",
        "analytic_tail": "Q_(9,n)(-100)>0 for every n>=1243",
    }
    rows = [
        SpliceRow(
            "co9m100fsc_01_bridge_target",
            "theorem_input",
            "ready_to_apply",
            "The first/full bridge isolates exactly two finite signs before its analytic half-line begins.",
            sources["bridge_splice"],
            "Bridge bookkeeping only.",
        ),
        SpliceRow(
            "co9m100fsc_02_new_coefficients",
            "interval_input",
            "ready_to_apply",
            "Retained-integral Arb quadrature rigorously encloses the two additional coefficients.",
            sources["new_coefficients"],
            "Two lambda=-100 coefficients only.",
            sources["new_coefficient_summary"],
        ),
        SpliceRow(
            "co9m100fsc_03_stable_rebuild",
            "interval_certificate",
            "ready_to_apply",
            "The full cancellation-preserving H4/H5/Q6/Q7/Q8 chain remains strictly positive through the enlarged endpoint.",
            "Q_(8,n)(-100)>0 for every 0<=n<=1244",
            "Finite 2048-bit Arb prefix only.",
            finite["positive_coordinate_counts"],
        ),
        SpliceRow(
            "co9m100fsc_04_splice_signs",
            "interval_theorem",
            "ready_to_apply",
            "Both new relative Q8 margins are strictly positive, closing the finite splice.",
            exact["splice"],
            "Finite signed-condensation signs only.",
            {"rows": splice_rows},
        ),
        SpliceRow(
            "co9m100fsc_05_combined_prefix",
            "analytic_composition",
            "ready_to_apply",
            "The established prefix and the two-row splice form one contiguous endpoint theorem.",
            exact["combined_prefix"],
            "Lambda=-100 endpoint only; the analytic tail remains conditional.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order9_m100_finite_splice_certificate",
        "date": "2026-07-13",
        "status": "rigorous two-index order-nine finite splice at lambda=-100",
        "proof_boundary": (
            "This artifact proves the order-nine endpoint signs only through "
            "n=1242. It does not prove the analytic tail, order-nine all-shift "
            "entry, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": sources,
        "source_diagnostics": source_diagnostics,
        "finite": {
            "lambda": finite["lambda"],
            "n_range": finite["n_range"],
            "coefficient_range": finite["coefficient_range"],
            "precision_bits": finite["precision_bits"],
            "positive_coordinate_counts": finite["positive_coordinate_counts"],
            "splice_rows": splice_rows,
            "minimum_relative_n": finite["minimum_relative_n"],
            "minimum_relative_ball": finite["minimum_relative_ball"],
            "minimum_relative_lower": finite["minimum_relative_lower"],
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_rows": len(rows),
            "open_rows": 0,
            "coefficients": len(values),
            "new_coefficient_rows": 2,
            "positive_Q8_rows": finite["positive_coordinate_counts"]["Q8_values"],
            "combined_positive_Q9_rows": len(finite["rows"]),
            "finite_splice_rows": len(splice_rows),
            "finite_splice_theorems": 1,
        },
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order9_m100_finite_splice_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order9_m100_finite_splice_certificate.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    finite = artifact["finite"]
    exact = artifact["exact"]
    lines = [
        "# Order-Nine Lambda=-100 Finite Splice Certificate",
        "",
        "Date: 2026-07-13",
        "",
        "Status: rigorous two-index endpoint splice. This is not a proof of",
        "the analytic order-nine tail, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order9_m100_finite_splice_certificate.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order9_m100_finite_splice_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order9_m100_finite_splice_certificate.py",
        "```",
        "",
        "## Added Coefficients",
        "",
        "Retained-integral quadrature with `n_sum=70`, cutoff 7, and 220",
        "decimal digits gives rigorous balls for `A_1257` and `A_1258`.",
        "Together with the established sources this covers `A_0,...,A_1258`.",
        "",
        "## Stable Signs",
        "",
        "```text",
        exact["signed_condensation"],
        exact["stable_coordinate"],
        exact["splice"],
        exact["combined_prefix"],
        "minimum relative margin at n="
        + str(finite["minimum_relative_n"])
        + ": "
        + finite["minimum_relative_ball"],
        "```",
        "",
        "The finite prefix now meets the analytic bridge with no index gap.",
        "The remaining handoff is the continuous sixth-nested theorem",
        "`w_1''(t)<=4200/t^2` on `t>=1250`.",
        "",
        "```text",
        "outputs/jensen_window_pf_compound_order9_m100_prefix_certificate.md",
        "outputs/jensen_window_pf_compound_order9_first_summand_curvature_bridge.md",
        "outputs/jensen_window_pf_compound_order9_m100_tail_curvature_reduction.md",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "wrote order-nine finite splice: "
        f"{summary['coefficients']} coefficients, "
        f"{summary['finite_splice_rows']} splice rows, "
        f"{summary['combined_positive_Q9_rows']} combined positive signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
