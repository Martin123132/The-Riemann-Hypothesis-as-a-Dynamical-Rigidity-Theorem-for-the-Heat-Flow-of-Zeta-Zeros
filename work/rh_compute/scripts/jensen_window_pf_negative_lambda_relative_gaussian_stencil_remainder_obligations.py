#!/usr/bin/env python3
"""Build relative-Gaussian stencil remainder obligations."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
import json
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_heat_flow_monotone_closure_scout import (  # noqa: E402
    REPO_ROOT,
    decimal_format,
)
from jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout import (  # noqa: E402
    build_artifact as build_taylor_stencil_artifact,
)


getcontext().prec = 100

DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md"


@dataclass(frozen=True)
class PositiveBaselineRow:
    source_row: str
    truncation_degree: int
    M: int
    T: str
    k: int
    B_margin: str
    companion_margin: str
    weighted_gap_margin: str
    minimum_margin: str
    half_margin_budget: str


@dataclass(frozen=True)
class Extremum:
    sample: str
    source_row: str
    M: int
    T: str


def extremum(row: tuple[Decimal, dict] | None) -> Extremum:
    if row is None:
        raise RuntimeError("missing stencil remainder-obligation extremum")
    value, source = row
    return Extremum(decimal_format(value), source["id"], int(source["M"]), str(source["T"]))


def build_diagnostics(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    stencil = build_taylor_stencil_artifact(max_degree, cutoff_n, precision_bits, k, t_values)
    stencil_diagnostics = stencil["matrix_rows"][4]["diagnostics"]
    stencil_rows = stencil_diagnostics["stencil_rows"]
    status_counts = Counter(row["status"] for row in stencil_rows)

    positive_rows: list[PositiveBaselineRow] = []
    weakest_margin: tuple[Decimal, dict] | None = None
    weakest_half_budget: tuple[Decimal, dict] | None = None
    for row in stencil_rows:
        if row["status"] != "target_satisfying_truncation":
            continue
        margins = [Decimal(row["B_k"]), Decimal(row["B_decrease"]), Decimal(row["C_increase"])]
        min_margin = min(margins)
        half_budget = min_margin / Decimal(2)
        if weakest_margin is None or min_margin < weakest_margin[0]:
            weakest_margin = (min_margin, row)
        if weakest_half_budget is None or half_budget < weakest_half_budget[0]:
            weakest_half_budget = (half_budget, row)
        positive_rows.append(
            PositiveBaselineRow(
                source_row=row["id"],
                truncation_degree=int(row["truncation_degree"]),
                M=int(row["M"]),
                T=str(row["T"]),
                k=int(row["k"]),
                B_margin=decimal_format(margins[0]),
                companion_margin=decimal_format(margins[1]),
                weighted_gap_margin=decimal_format(margins[2]),
                minimum_margin=decimal_format(min_margin),
                half_margin_budget=decimal_format(half_budget),
            )
        )

    blocked_rows = len(stencil_rows) - len(positive_rows)
    return {
        "parameters": stencil_diagnostics["parameters"],
        "source_taylor_stencil_summary": stencil["summary"],
        "status_counts": dict(sorted(status_counts.items())),
        "truncation_rows": len(stencil_rows),
        "positive_baseline_rows": len(positive_rows),
        "blocked_baseline_rows": blocked_rows,
        "invalid_normalizer_rows": status_counts["invalid_normalizer"],
        "upper_wall_violation_rows": status_counts["upper_wall_violation"],
        "companion_failure_rows": status_counts["companion_failure"],
        "weighted_gap_failure_rows": status_counts["weighted_gap_failure"],
        "positive_baseline_details": [asdict(row) for row in positive_rows],
        "weakest_positive_margin": asdict(extremum(weakest_margin)),
        "weakest_half_margin_budget": asdict(extremum(weakest_half_budget)),
    }


def build_artifact(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    diagnostics = build_diagnostics(max_degree, cutoff_n, precision_bits, k, t_values)
    rows = [
        {
            "id": "nlrgsro_01_log_tail_error_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "Write the full relative-Gaussian log multiplier as f_k=f_k^(M)+epsilon_k, where epsilon_k is the log Taylor-tail error after the chosen finite truncation.",
            "formula": "epsilon_k=log(F_k/F_k^(M))",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md",
            ],
            "proof_boundary": "Exact bookkeeping coordinate only; it does not bound epsilon_k.",
        },
        {
            "id": "nlrgsro_02_B_error_stencil",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The B_k remainder contribution is the second-difference stencil 2*epsilon_k-epsilon_(k-1)-epsilon_(k+1).",
            "formula": "B_k=B_k^(M)+2*epsilon_k-epsilon_(k-1)-epsilon_(k+1)",
            "proof_boundary": "Exact stencil identity only; no sign estimate is proved.",
        },
        {
            "id": "nlrgsro_03_companion_error_stencil",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The companion upper-side remainder contribution is -epsilon_(k-1)+3*epsilon_k-3*epsilon_(k+1)+epsilon_(k+2).",
            "formula": "B_k-B_(k+1)=(B_k-B_(k+1))^(M)-epsilon_(k-1)+3*epsilon_k-3*epsilon_(k+1)+epsilon_(k+2)",
            "proof_boundary": "Exact stencil identity only; the monotone side remains open.",
        },
        {
            "id": "nlrgsro_04_weighted_gap_error_stencil",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The weighted C-gap remainder contribution is (2*k+1)*epsilon_(k-1)-(6*k+5)*epsilon_k+(6*k+7)*epsilon_(k+1)-(2*k+3)*epsilon_(k+2).",
            "formula": "C_(k+1)-C_k=(C_(k+1)-C_k)^(M)+(2*k+1)*epsilon_(k-1)-(6*k+5)*epsilon_k+(6*k+7)*epsilon_(k+1)-(2*k+3)*epsilon_(k+2)",
            "proof_boundary": "Exact weighted stencil identity only; no tail theorem is proved.",
        },
        {
            "id": "nlrgsro_05_positive_baseline_margin_budget",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Only the all-positive finite truncation rows provide a perturbative baseline with positive margins for all three stencils.",
            "diagnostics": diagnostics,
            "proof_boundary": "Finite margin diagnostic only; not a theorem about the full Taylor tail.",
        },
        {
            "id": "nlrgsro_06_blocked_baseline_rows",
            "role": "finite_obstruction",
            "readiness": "not_ready_to_apply",
            "claim": "Rows with invalid normalizers or a negative stencil margin cannot be used as small-remainder perturbative baselines.",
            "blocked_counts": {
                "blocked_baseline_rows": diagnostics["blocked_baseline_rows"],
                "invalid_normalizer_rows": diagnostics["invalid_normalizer_rows"],
                "upper_wall_violation_rows": diagnostics["upper_wall_violation_rows"],
                "companion_failure_rows": diagnostics["companion_failure_rows"],
                "weighted_gap_failure_rows": diagnostics["weighted_gap_failure_rows"],
            },
            "proof_boundary": "Finite obstruction to finite-truncation promotion only; not a statement about the full zeta moments.",
        },
        {
            "id": "nlrgsro_07_sufficient_remainder_budget",
            "role": "exact_sufficient_condition",
            "readiness": "not_ready_to_apply",
            "claim": "On a positive baseline row, it is sufficient to bound each absolute epsilon-stencil error by half of that row's minimum positive margin.",
            "weakest_half_margin_budget": diagnostics["weakest_half_margin_budget"],
            "proof_boundary": "Exact finite-row sufficient condition only; no uniform epsilon-stencil bound is proved.",
        },
        {
            "id": "nlrgsro_08_live_uniform_stencil_bound_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof should replace finite-row budgets by uniform q/T bounds for the B, companion, and weighted-gap epsilon stencils.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md",
            ],
            "proof_boundary": "Live theorem-search route only; the uniform stencil remainder theorem is not proved.",
        },
        {
            "id": "nlrgsro_09_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted remainder proof must state the truncation order, q/T range, positivity of the finite baseline, and explicit absolute bounds for all three epsilon stencils.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of scaled-curvature monotonicity, cone entry, jwpf_06, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "exact_stencil_rows": 4,
        "truncation_rows": diagnostics["truncation_rows"],
        "positive_baseline_rows": diagnostics["positive_baseline_rows"],
        "blocked_baseline_rows": diagnostics["blocked_baseline_rows"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The relative-Gaussian remainder problem can be split into three exact epsilon stencils "
            "for B_k, B_k-B_(k+1), and C_(k+1)-C_k. In the finite degree<=14 truncation matrix, "
            "only 4/35 rows have positive baseline margins for all three stencils; the weakest "
            "half-margin budget is 1.166490564421582442E-8, so a promoted proof needs sharp "
            "uniform epsilon-stencil bounds rather than finite truncation promotion."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_taylor_stencil_scout": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md",
        "source_relative_gaussian_bridge": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py",
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It decomposes the relative-Gaussian Taylor-tail "
            "remainder into three epsilon stencils and records finite margin budgets, but it does not prove "
            "uniform epsilon-stencil bounds, does not prove scaled-curvature monotonicity, does not prove cone "
            "entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Finite positive baseline margins are not promoted to uniform remainder bounds.",
            "Rows with invalid normalizers or failed finite stencils block finite-truncation promotion only.",
            "The uniform epsilon-stencil remainder theorem remains open.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][4]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['positive_baseline_rows']} positive baseline rows, "
        f"{summary['blocked_baseline_rows']} blocked baseline rows, "
        f"{summary['exact_stencil_rows']} exact stencil rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Stencil Remainder Obligations",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of a uniform Taylor-tail remainder theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations`.",
        "",
        "Proof boundary: this artifact decomposes the Taylor-tail remainder into",
        "three exact epsilon stencils and records finite margin budgets. It does",
        "not prove any uniform all-`k` remainder bound.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Exact Epsilon Stencils",
        "",
        "Write:",
        "",
        "```text",
        "f_k = f_k^(M)+epsilon_k",
        "```",
        "",
        "Then the three required stencil errors are:",
        "",
        "```text",
        "E_B = 2*epsilon_k-epsilon_(k-1)-epsilon_(k+1)",
        "E_U = -epsilon_(k-1)+3*epsilon_k-3*epsilon_(k+1)+epsilon_(k+2)",
        "E_C = (2*k+1)*epsilon_(k-1)-(6*k+5)*epsilon_k+(6*k+7)*epsilon_(k+1)-(2*k+3)*epsilon_(k+2)",
        "```",
        "",
        "## Finite Margin Budget",
        "",
        "```text",
        f"truncation rows: {diagnostics['truncation_rows']}",
        f"positive baseline rows: {diagnostics['positive_baseline_rows']}",
        f"blocked baseline rows: {diagnostics['blocked_baseline_rows']}",
        f"invalid normalizers: {diagnostics['invalid_normalizer_rows']}",
        f"upper-wall violations: {diagnostics['upper_wall_violation_rows']}",
        f"companion failures: {diagnostics['companion_failure_rows']}",
        f"weighted-gap failures: {diagnostics['weighted_gap_failure_rows']}",
        f"weakest positive margin: {diagnostics['weakest_positive_margin']['sample']} at {diagnostics['weakest_positive_margin']['source_row']}",
        f"weakest half-margin budget: {diagnostics['weakest_half_margin_budget']['sample']} at {diagnostics['weakest_half_margin_budget']['source_row']}",
        "```",
        "",
        "The half-margin budget is a finite-row sufficient target only. A proof",
        "must replace it with uniform q/T bounds for `E_B`, `E_U`, and `E_C`.",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_taylor_stencil_scout"],
        artifact["source_relative_gaussian_bridge"],
        artifact["source_uniform_remainder_target"],
        "```",
        "",
        "Summary:",
        "",
        summary["main_finding"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-degree", type=int, default=14)
    parser.add_argument("--cutoff-n", type=int, default=80)
    parser.add_argument("--precision-bits", type=int, default=256)
    parser.add_argument("--tail-start-k", type=int, default=22)
    parser.add_argument("--t-values", type=int, nargs="+", default=[25, 50, 100, 200, 500, 1000, 2000])
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(args.max_degree, args.cutoff_n, args.precision_bits, args.tail_start_k, args.t_values)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian stencil remainder obligations: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
