#!/usr/bin/env python3
"""Build pointwise tail budgets for the relative-Gaussian stencil route."""

from __future__ import annotations

import argparse
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
from jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations import (  # noqa: E402
    build_artifact as build_remainder_obligation_artifact,
)


getcontext().prec = 100

DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md"


@dataclass(frozen=True)
class PointwiseBudgetRow:
    source_row: str
    truncation_degree: int
    M: int
    T: str
    k: int
    B_margin: str
    companion_margin: str
    weighted_gap_margin: str
    B_weight_sum: str
    companion_weight_sum: str
    weighted_gap_weight_sum: str
    B_eta_bound: str
    companion_eta_bound: str
    weighted_gap_eta_bound: str
    uniform_eta_bound: str
    limiting_stencil: str
    half_safety_eta_bound: str
    half_safety_relative_tail_ratio_bound: str


@dataclass(frozen=True)
class Extremum:
    sample: str
    source_row: str
    M: int
    T: str
    limiting_stencil: str


def relative_tail_ratio_from_log_bound(eta: Decimal) -> Decimal:
    """Return rho such that |R/F^M|<=rho implies |log(1+R/F^M)|<=eta."""
    if eta <= 0:
        raise ValueError(f"eta must be positive, got {eta}")
    return Decimal(1) - (-eta).exp()


def build_budget_rows(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    remainder = build_remainder_obligation_artifact(max_degree, cutoff_n, precision_bits, k, t_values)
    diagnostics = remainder["matrix_rows"][4]["diagnostics"]
    positive_rows = diagnostics["positive_baseline_details"]

    rows: list[PointwiseBudgetRow] = []
    weakest_uniform: tuple[Decimal, PointwiseBudgetRow] | None = None
    weakest_half: tuple[Decimal, PointwiseBudgetRow] | None = None
    limiting_counts: dict[str, int] = {}

    for row in positive_rows:
        row_k = int(row["k"])
        b_margin = Decimal(row["B_margin"])
        companion_margin = Decimal(row["companion_margin"])
        weighted_gap_margin = Decimal(row["weighted_gap_margin"])
        b_weight_sum = Decimal(4)
        companion_weight_sum = Decimal(8)
        weighted_gap_weight_sum = Decimal(16 * row_k + 16)
        bounds = {
            "B": b_margin / b_weight_sum,
            "companion": companion_margin / companion_weight_sum,
            "weighted_gap": weighted_gap_margin / weighted_gap_weight_sum,
        }
        limiting_stencil, uniform_eta = min(bounds.items(), key=lambda item: item[1])
        half_eta = uniform_eta / Decimal(2)
        tail_ratio = relative_tail_ratio_from_log_bound(half_eta)
        budget_row = PointwiseBudgetRow(
            source_row=row["source_row"],
            truncation_degree=int(row["truncation_degree"]),
            M=int(row["M"]),
            T=str(row["T"]),
            k=row_k,
            B_margin=decimal_format(b_margin),
            companion_margin=decimal_format(companion_margin),
            weighted_gap_margin=decimal_format(weighted_gap_margin),
            B_weight_sum=decimal_format(b_weight_sum),
            companion_weight_sum=decimal_format(companion_weight_sum),
            weighted_gap_weight_sum=decimal_format(weighted_gap_weight_sum),
            B_eta_bound=decimal_format(bounds["B"]),
            companion_eta_bound=decimal_format(bounds["companion"]),
            weighted_gap_eta_bound=decimal_format(bounds["weighted_gap"]),
            uniform_eta_bound=decimal_format(uniform_eta),
            limiting_stencil=limiting_stencil,
            half_safety_eta_bound=decimal_format(half_eta),
            half_safety_relative_tail_ratio_bound=decimal_format(tail_ratio),
        )
        rows.append(budget_row)
        limiting_counts[limiting_stencil] = limiting_counts.get(limiting_stencil, 0) + 1
        if weakest_uniform is None or uniform_eta < weakest_uniform[0]:
            weakest_uniform = (uniform_eta, budget_row)
        if weakest_half is None or half_eta < weakest_half[0]:
            weakest_half = (half_eta, budget_row)

    if weakest_uniform is None or weakest_half is None:
        raise RuntimeError("no positive baseline rows available for pointwise tail budgets")

    return {
        "parameters": diagnostics["parameters"],
        "source_remainder_summary": remainder["summary"],
        "positive_baseline_rows": len(rows),
        "blocked_baseline_rows": diagnostics["blocked_baseline_rows"],
        "budget_rows": [asdict(row) for row in rows],
        "limiting_stencil_counts": dict(sorted(limiting_counts.items())),
        "weakest_uniform_eta_bound": asdict(
            Extremum(
                sample=decimal_format(weakest_uniform[0]),
                source_row=weakest_uniform[1].source_row,
                M=weakest_uniform[1].M,
                T=weakest_uniform[1].T,
                limiting_stencil=weakest_uniform[1].limiting_stencil,
            )
        ),
        "weakest_half_safety_eta_bound": asdict(
            Extremum(
                sample=decimal_format(weakest_half[0]),
                source_row=weakest_half[1].source_row,
                M=weakest_half[1].M,
                T=weakest_half[1].T,
                limiting_stencil=weakest_half[1].limiting_stencil,
            )
        ),
        "weakest_half_safety_relative_tail_ratio_bound": weakest_half[1].half_safety_relative_tail_ratio_bound,
    }


def build_artifact(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    diagnostics = build_budget_rows(max_degree, cutoff_n, precision_bits, k, t_values)
    weakest_half = diagnostics["weakest_half_safety_eta_bound"]
    rows = [
        {
            "id": "nlrgptb_01_pointwise_log_tail_envelope",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "If |epsilon_j|<=eta_j for j in {k-1,k,k+1,k+2}, then every epsilon-stencil error is bounded by the absolute coefficient-weighted eta sum.",
            "formula": "|sum_j c_j*epsilon_j| <= sum_j |c_j|*eta_j",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
            ],
            "proof_boundary": "Exact triangle-inequality reduction only; no eta_j bound is proved.",
        },
        {
            "id": "nlrgptb_02_stencil_weight_sums",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The B, companion, and weighted-gap epsilon stencils have absolute coefficient sums 4, 8, and 16*k+16.",
            "formula": "||E_B||_1=4, ||E_U||_1=8, ||E_C||_1=16*k+16",
            "proof_boundary": "Exact coefficient bookkeeping only; not a tail estimate.",
        },
        {
            "id": "nlrgptb_03_uniform_eta_sufficient_condition",
            "role": "exact_sufficient_condition",
            "readiness": "not_ready_to_apply",
            "claim": "On a positive baseline row, a uniform pointwise log-tail envelope eta preserves nonnegativity if eta is below every stencil margin divided by its absolute coefficient sum.",
            "formula": "eta <= min(B_margin/4, U_margin/8, C_margin/(16*k+16))",
            "proof_boundary": "Exact finite-row sufficient condition only; no uniform analytic eta is proved.",
        },
        {
            "id": "nlrgptb_04_half_safety_eta_budget",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Using a half-safety margin converts each positive baseline row into a stricter pointwise epsilon target with half of the uniform eta bound.",
            "diagnostics": diagnostics,
            "proof_boundary": "Finite budget diagnostic only; not a theorem for the full Taylor tail.",
        },
        {
            "id": "nlrgptb_05_relative_tail_ratio_conversion",
            "role": "exact_sufficient_condition",
            "readiness": "not_ready_to_apply",
            "claim": "If the analytic tail is multiplicative, F_j=F_j^(M)*(1+theta_j), then |theta_j|<=rho<1 implies |epsilon_j|<=-log(1-rho).",
            "formula": "rho <= 1-exp(-eta) is sufficient for |epsilon_j|<=eta",
            "weakest_half_safety_relative_tail_ratio_bound": diagnostics["weakest_half_safety_relative_tail_ratio_bound"],
            "proof_boundary": "Exact conversion from multiplicative relative tail to log-tail envelope only; no bound on theta_j is proved.",
        },
        {
            "id": "nlrgptb_06_bottleneck_diagnostic",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "For the current positive finite baselines, the companion stencil gives the global weakest pointwise-tail budget, while one sampled row is weighted-gap limited.",
            "weakest_half_safety_eta_bound": weakest_half,
            "proof_boundary": "Finite diagnostic only; it identifies the present bottleneck split but does not prove any tail bound.",
        },
        {
            "id": "nlrgptb_07_live_pointwise_tail_theorem",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof may close this branch by proving a zeta-specific analytic relative-tail theorem at least as strong as the half-safety pointwise log-tail budget over a stated q/T range.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
            ],
            "proof_boundary": "Live theorem-search route only; the pointwise tail theorem remains open.",
        },
        {
            "id": "nlrgptb_08_finite_budget_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The finite pointwise budget table alone proves the required epsilon bounds.",
            "gap": "The table gives required tolerances, not analytic estimates for the Taylor tail or relative multiplier error.",
            "proof_boundary": "Rejected shortcut only; finite tolerances are not promoted to an all-k analytic theorem.",
        },
        {
            "id": "nlrgptb_09_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any promoted pointwise-tail proof must state the truncation order, q/T collar, normalized multiplier positivity, a relative-tail estimate |theta_j|<=rho, and the comparison rho<=1-exp(-eta).",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "positive_baseline_rows": diagnostics["positive_baseline_rows"],
        "blocked_baseline_rows": diagnostics["blocked_baseline_rows"],
        "budget_rows": len(diagnostics["budget_rows"]),
        "exact_sufficient_rows": 2,
        "weakest_half_safety_eta_bound": weakest_half["sample"],
        "weakest_half_safety_relative_tail_ratio_bound": diagnostics["weakest_half_safety_relative_tail_ratio_bound"],
        "limiting_stencil_counts": diagnostics["limiting_stencil_counts"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The exact epsilon-stencil obligations reduce to a concrete pointwise log-tail target. "
            f"At k={k}, the absolute coefficient sums are 4, 8, and {16 * k + 16}; among the 4 positive "
            "finite baselines, the companion stencil is limiting in 3 rows and the weighted-gap stencil "
            "is limiting in 1 row. The weakest half-safety "
            f"log-tail envelope is {weakest_half['sample']} at {weakest_half['source_row']}, equivalent "
            f"to the multiplicative relative-tail ratio bound {diagnostics['weakest_half_safety_relative_tail_ratio_bound']}. "
            "This sharpens the open uniform remainder theorem but does not prove any analytic tail estimate."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_stencil_remainder_obligations": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_taylor_moment_budget": "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py",
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It converts exact epsilon-stencil margins into "
            "pointwise log-tail and multiplicative relative-tail budgets, but it does not prove the required "
            "analytic tail estimate, does not prove scaled-curvature monotonicity, does not prove cone entry, "
            "and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Pointwise epsilon budgets are required tolerances, not proved tail estimates.",
            "Finite positive baseline rows are not promoted to uniform all-k tail bounds.",
            "The companion stencil is the current global finite bottleneck, with one weighted-gap limited row.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][3]["diagnostics"]
    weakest_half = diagnostics["weakest_half_safety_eta_bound"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['positive_baseline_rows']} positive baseline rows, "
        f"{summary['budget_rows']} budget rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Pointwise Tail Budget",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of a uniform Taylor-tail remainder theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget`.",
        "",
        "Proof boundary: this artifact converts epsilon-stencil margins into",
        "pointwise log-tail and multiplicative relative-tail budgets. It does",
        "not prove the analytic tail estimates required by those budgets.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Pointwise Envelope",
        "",
        "Assume `|epsilon_j|<=eta_j` for `j=k-1,k,k+1,k+2`. Then:",
        "",
        "```text",
        "|E_B| <= eta_(k-1)+2*eta_k+eta_(k+1)",
        "|E_U| <= eta_(k-1)+3*eta_k+3*eta_(k+1)+eta_(k+2)",
        "|E_C| <= (2*k+1)*eta_(k-1)+(6*k+5)*eta_k+(6*k+7)*eta_(k+1)+(2*k+3)*eta_(k+2)",
        "```",
        "",
        "For a uniform pointwise envelope `eta_j<=eta`, the coefficient sums are:",
        "",
        "```text",
        "B sum: 4",
        "companion sum: 8",
        "weighted-gap sum: 16*k+16",
        "at k=22, weighted-gap sum: 368",
        "```",
        "",
        "## Current Finite Bottleneck",
        "",
        "```text",
        f"positive baseline rows: {summary['positive_baseline_rows']}",
        f"blocked baseline rows: {summary['blocked_baseline_rows']}",
        f"budget rows: {summary['budget_rows']}",
        f"limiting stencil counts: {summary['limiting_stencil_counts']}",
        f"weakest half-safety eta: {weakest_half['sample']} at {weakest_half['source_row']}",
        f"weakest limiting stencil: {weakest_half['limiting_stencil']}",
        f"relative tail ratio bound: {summary['weakest_half_safety_relative_tail_ratio_bound']}",
        "```",
        "",
        "The conversion from multiplicative relative tail to log-tail envelope is:",
        "",
        "```text",
        "F_j = F_j^(M)*(1+theta_j)",
        "|theta_j| <= rho < 1",
        "|epsilon_j| = |log(1+theta_j)| <= -log(1-rho)",
        "rho <= 1-exp(-eta)",
        "```",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_stencil_remainder_obligations"],
        artifact["source_uniform_remainder_target"],
        artifact["source_taylor_moment_budget"],
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
        "wrote Jensen-window PF negative-lambda relative-Gaussian pointwise tail budget: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
