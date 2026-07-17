#!/usr/bin/env python3
"""Prove inverse-eleventh-power first-summand dominance with the rebalanced split."""

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


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/"
    "jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.md"
)
POWER10_SOURCE = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_negative_lambda_first_summand_power10_rebalanced_dominance_extension.json"
)
START_K = 300
POWER = 11


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
    source = json.loads(POWER10_SOURCE.read_text(encoding="utf-8"))
    if source.get("summary", {}).get("full_tail_power") != 10:
        raise RuntimeError("power-ten dominance source changed")
    if source.get("summary", {}).get("tail_start_k") != START_K:
        raise RuntimeError("power-ten dominance start changed")
    if source.get("summary", {}).get("positive_analytic_gates") != 14:
        raise RuntimeError("power-ten analytic gate count changed")
    bound = source.get("diagnostics", {}).get("full_tail_relative_bound")
    if bound != (
        "0<=delta_k=(M_k-M_k^(1))/M_k^(1)<2/k^10 "
        "for every integer k>=300"
    ):
        raise RuntimeError("power-ten dominance contract changed")
    return {
        "status": source.get("status"),
        "power10_bound": bound,
        "adaptive_split": source["diagnostics"]["parameters"]["adaptive_a"],
    }


def build_artifact() -> dict:
    source = source_contract()
    diagnostics = power8.endpoint_diagnostics(power=POWER, start_k=START_K)
    rows = [
        DominanceRow(
            "fsd11r_01_power10_input",
            "theorem_input",
            "ready_to_apply",
            "The exact summand-ratio geometry and rebalanced split are inherited unchanged.",
            source["power10_bound"] + "; " + source["adaptive_split"],
            "Inherited lambda=-100 kernel geometry only.",
        ),
        DominanceRow(
            "fsd11r_02_high_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The higher-summand ratio on the high region has enough exponential reserve for an eleventh inverse power.",
            diagnostics["high_region_bound"],
            "Arb endpoint gates plus exact derivative monotonicity.",
        ),
        DominanceRow(
            "fsd11r_03_saddle_geometry",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The log-integrand derivative remains positive across the comparison collar for every k>=300.",
            "S_k'(a(k))>100 and S_k'(c(k))>0",
            "Arb endpoint gates plus six exact ratio-monotonicity checks.",
        ),
        DominanceRow(
            "fsd11r_04_low_region",
            "interval_analytic_certificate",
            "ready_to_apply",
            "The tilted first-summand probability below the split has enough reserve for an eleventh inverse power.",
            diagnostics["low_region_bound"],
            "Arb endpoint gate and exact half-line derivative comparisons.",
        ),
        DominanceRow(
            "fsd11r_05_full_tail",
            "analytic_theorem",
            "ready_to_apply",
            "The two regions compose to an inverse-eleventh-power complete-to-first moment defect.",
            diagnostics["full_tail_relative_bound"],
            "Lambda=-100 integer moments only.",
        ),
        DominanceRow(
            "fsd11r_06_wall_error",
            "exact_consequence",
            "ready_to_apply",
            "The adjacent logarithmic wall inherits the inverse-eleventh-power stencil bound.",
            diagnostics["moment_log_error"] + "; " + diagnostics["adjacent_B_error"],
            "Uses log(1+z)<=z and the exact centered stencil.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension",
        "date": "2026-07-16",
        "status": "rebalanced all-k inverse-eleventh-power first-summand dominance extension",
        "proof_boundary": (
            "This artifact strengthens first-summand dominance only for lambda=-100 "
            "and k>=300. It does not prove the order-ten curvature ceiling, "
            "full-kernel transfer, PF-infinity, RH, or Lambda<=0."
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
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    d = artifact["diagnostics"]
    lines = [
        "# Jensen-Window PF Rebalanced Power-Eleven First-Summand Dominance",
        "",
        "Date: 2026-07-16",
        "",
        "Status: rigorous inverse-eleventh-power first-summand dominance at",
        "`lambda=-100`. This is not a proof of order ten, PF-infinity, RH,",
        "or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_power11_rebalanced_dominance_extension.py",
        "```",
        "",
        "## Reused Split",
        "",
        "```text",
        "a(k)=log(k)/10, b(k)=a(k)+1/10, c(k)=b(k)+1/100.",
        "```",
        "",
        "At `k=300`, all fourteen Arb endpoint and derivative gates remain",
        "strictly positive, and the six exact ratio-log derivatives cover the",
        "complete half-line.",
        "",
        "```text",
        d["high_region_bound"],
        d["low_region_bound"],
        d["full_tail_relative_bound"],
        d["moment_log_error"],
        d["adjacent_B_error"],
        "```",
        "",
        "The next power requires a slightly later endpoint gate and is certified",
        "separately rather than inferred from this artifact.",
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
        "wrote power-eleven rebalanced first-summand dominance: "
        f"{summary['rows']} rows, "
        f"{summary['positive_analytic_gates']} positive analytic gates"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
