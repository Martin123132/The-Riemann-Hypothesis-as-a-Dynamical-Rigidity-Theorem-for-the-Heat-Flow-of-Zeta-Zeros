#!/usr/bin/env python3
"""Prove inverse-thirteenth-power first-summand dominance."""

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

import jensen_window_pf_negative_lambda_first_summand_power8_rebalanced_dominance_extension as power8  # noqa: E402


DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension.md"
POWER12_SOURCE = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power12_rebalanced_dominance_extension.json"
SOURCE_START_K = 320
START_K = 340
POWER = 13


@dataclass(frozen=True)
class DominanceRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def source_contract() -> dict:
    source = json.loads(POWER12_SOURCE.read_text(encoding="utf-8"))
    if source.get("summary", {}).get("full_tail_power") != 12:
        raise RuntimeError("power-twelve dominance source changed")
    if source.get("summary", {}).get("tail_start_k") != SOURCE_START_K:
        raise RuntimeError("power-twelve dominance start changed")
    if source.get("summary", {}).get("positive_analytic_gates") != 14:
        raise RuntimeError("power-twelve analytic gate count changed")
    bound = source.get("diagnostics", {}).get("full_tail_relative_bound")
    if bound != (
        "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^12 "
        "for every integer k>=320"
    ):
        raise RuntimeError("power-twelve dominance contract changed")
    return {
        "status": source.get("status"),
        "power12_bound": bound,
        "adaptive_split": source["diagnostics"]["parameters"]["adaptive_a"],
    }


def build_artifact() -> dict:
    source = source_contract()
    diagnostics = power8.endpoint_diagnostics(power=POWER, start_k=START_K)
    rows = [
        DominanceRow(
            "fsd13r_01_power12_input",
            "theorem_input",
            "ready_to_apply",
            "The exact summand-ratio geometry and rebalanced split are inherited unchanged.",
            source["power12_bound"] + "; " + source["adaptive_split"],
            "Inherited lambda=-100 kernel geometry only.",
        ),
        DominanceRow(
            "fsd13r_02_high_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The high-region exponential reserve supports a thirteenth inverse power.",
            diagnostics["high_region_bound"],
            "Arb endpoint gates plus exact derivative monotonicity.",
        ),
        DominanceRow(
            "fsd13r_03_saddle_geometry",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The log-integrand derivative remains positive across the comparison collar for every k>=340.",
            "S_k'(a(k))>100 and S_k'(c(k))>0",
            "Arb endpoint gates plus six exact ratio-monotonicity checks.",
        ),
        DominanceRow(
            "fsd13r_04_low_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The tilted first-summand low-region probability supports a thirteenth inverse power.",
            diagnostics["low_region_bound"],
            "Arb endpoint gate and exact half-line derivative comparisons.",
        ),
        DominanceRow(
            "fsd13r_05_full_tail",
            "analytic_theorem",
            "ready_to_apply",
            "The two regions compose to an inverse-thirteenth-power complete-to-first moment defect.",
            diagnostics["full_tail_relative_bound"],
            "Lambda=-100 integer moments only.",
        ),
        DominanceRow(
            "fsd13r_06_wall_error",
            "exact_consequence",
            "ready_to_apply",
            "The adjacent logarithmic wall inherits the inverse-thirteenth-power stencil bound.",
            diagnostics["moment_log_error"] + "; " + diagnostics["adjacent_B_error"],
            "Uses log(1+z)<=z and the exact centered stencil.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension",
        "date": "2026-07-16",
        "status": "rebalanced inverse-thirteenth-power first-summand dominance extension from k=340",
        "proof_boundary": (
            "This strengthens first-summand dominance only for lambda=-100 and "
            "k>=340. It does not prove an order-eleven curvature ceiling, "
            "full-kernel transfer, endpoint entry, PF-infinity, RH, or Lambda<=0."
        ),
        "source_contract": source,
        "diagnostics": diagnostics,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "ready_to_apply_rows": len(rows),
            "positive_analytic_gates": len(diagnostics["positive_gates"]),
            "full_tail_power": POWER,
            "tail_start_k": START_K,
            "dominance_theorems": 1,
            "rebalanced_splits": 1,
        },
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power13_rebalanced_dominance_extension.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    d = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Rebalanced Power-Thirteen First-Summand Dominance",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous inverse-thirteenth-power first-summand dominance at",
        "`lambda=-100` for `k>=340`. This is not a proof of order eleven,",
        "PF-infinity, RH, or `Lambda<=0`.",
        "",
        "```text",
        "a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.",
        d["high_region_bound"],
        d["low_region_bound"],
        d["full_tail_relative_bound"],
        d["moment_log_error"],
        d["adjacent_B_error"],
        "```",
        "",
        "At `k=340`, all fourteen Arb endpoint and derivative gates are strict.",
        "Exact monotonicity then covers every larger integer. This is the",
        "natural perturbation power for an eighth stable-log transfer; that",
        "order-eleven transfer remains a separate theorem.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(f"wrote power-thirteen rebalanced first-summand dominance: {summary['rows']} rows, {summary['positive_analytic_gates']} positive analytic gates")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
