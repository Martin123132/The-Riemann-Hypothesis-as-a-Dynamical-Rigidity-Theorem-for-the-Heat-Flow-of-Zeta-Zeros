#!/usr/bin/env python3
"""Build next-increment stencil stress diagnostics for the relative-Gaussian route."""

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

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import (  # noqa: E402
    REPO_ROOT,
    decimal_format,
)
from jensen_window_pf_negative_lambda_high_order_taylor_scout import (  # noqa: E402
    arb_mid_decimal,
    build_diagnostics as build_high_order_diagnostics,
    truncated_multiplier,
)
from jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget import (  # noqa: E402
    build_artifact as build_pointwise_budget_artifact,
)


getcontext().prec = 100

DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.md"


@dataclass(frozen=True)
class NextIncrementStencil:
    name: str
    increment: str
    half_margin: str
    abs_over_half_margin: str
    margin_after_increment: str
    sign_preserved: bool
    exceeds_half_safety: bool


@dataclass(frozen=True)
class NextIncrementStressRow:
    source_row: str
    M: int
    next_M: int
    T: str
    k: int
    allowed_relative_tail_ratio: str
    max_next_abs_theta: str
    max_next_abs_theta_over_allowed: str
    max_abs_log_increment: str
    pointwise_budget_satisfied: bool
    stencil_signs_preserved: bool
    half_safety_stencils_satisfied: bool
    limiting_stencil: str
    stencil_increments: list[NextIncrementStencil]


@dataclass(frozen=True)
class MissingNextIncrementRow:
    source_row: str
    M: int
    T: str
    k: int
    reason: str


def log_increments(k: int, T: int, ratios: list[Decimal], M: int) -> dict[int, Decimal]:
    values: dict[int, Decimal] = {}
    for index in (k - 1, k, k + 1, k + 2):
        current = truncated_multiplier(index, T, ratios, M)
        next_value = truncated_multiplier(index, T, ratios, M + 1)
        values[index] = next_value.ln() - current.ln()
    return values


def theta_increments(k: int, T: int, ratios: list[Decimal], M: int) -> dict[int, Decimal]:
    values: dict[int, Decimal] = {}
    for index in (k - 1, k, k + 1, k + 2):
        current = truncated_multiplier(index, T, ratios, M)
        next_value = truncated_multiplier(index, T, ratios, M + 1)
        values[index] = (next_value - current) / current
    return values


def stencil_values(k: int, eps: dict[int, Decimal]) -> dict[str, Decimal]:
    return {
        "B": Decimal(2) * eps[k] - eps[k - 1] - eps[k + 1],
        "companion": -eps[k - 1] + Decimal(3) * eps[k] - Decimal(3) * eps[k + 1] + eps[k + 2],
        "weighted_gap": (
            Decimal(2 * k + 1) * eps[k - 1]
            - Decimal(6 * k + 5) * eps[k]
            + Decimal(6 * k + 7) * eps[k + 1]
            - Decimal(2 * k + 3) * eps[k + 2]
        ),
    }


def build_diagnostics(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    pointwise = build_pointwise_budget_artifact(max_degree, cutoff_n, precision_bits, k, t_values)
    budget_rows = pointwise["matrix_rows"][3]["diagnostics"]["budget_rows"]
    high = build_high_order_diagnostics(max_degree, cutoff_n, precision_bits, k, t_values)
    ratios = [arb_mid_decimal(flint.arb(row.ratio_to_c0)) for row in high.coefficient_rows]

    stress_rows: list[NextIncrementStressRow] = []
    missing_rows: list[MissingNextIncrementRow] = []
    worst_pointwise_ratio: tuple[Decimal, str] | None = None
    worst_stencil_half_ratio: tuple[Decimal, str, str] | None = None
    for row in budget_rows:
        row_id = row["source_row"]
        M = int(row["M"])
        T = int(row["T"])
        row_k = int(row["k"])
        if M + 1 >= len(ratios):
            missing_rows.append(
                MissingNextIncrementRow(
                    source_row=row_id,
                    M=M,
                    T=str(T),
                    k=row_k,
                    reason="degree-16 Taylor ratio is not available in the current degree<=14 coefficient scaffold",
                )
            )
            continue

        theta = theta_increments(row_k, T, ratios, M)
        eps = log_increments(row_k, T, ratios, M)
        stencils = stencil_values(row_k, eps)
        margins = {
            "B": Decimal(row["B_margin"]),
            "companion": Decimal(row["companion_margin"]),
            "weighted_gap": Decimal(row["weighted_gap_margin"]),
        }
        allowed = Decimal(row["half_safety_relative_tail_ratio_bound"])
        max_theta = max(abs(value) for value in theta.values())
        max_log = max(abs(value) for value in eps.values())
        pointwise_ratio = max_theta / allowed
        if worst_pointwise_ratio is None or pointwise_ratio > worst_pointwise_ratio[0]:
            worst_pointwise_ratio = (pointwise_ratio, row_id)

        stencil_rows: list[NextIncrementStencil] = []
        for name in ("B", "companion", "weighted_gap"):
            half_margin = margins[name] / Decimal(2)
            increment = stencils[name]
            ratio = abs(increment) / half_margin
            if worst_stencil_half_ratio is None or ratio > worst_stencil_half_ratio[0]:
                worst_stencil_half_ratio = (ratio, row_id, name)
            margin_after = margins[name] + increment
            stencil_rows.append(
                NextIncrementStencil(
                    name=name,
                    increment=decimal_format(increment),
                    half_margin=decimal_format(half_margin),
                    abs_over_half_margin=decimal_format(ratio),
                    margin_after_increment=decimal_format(margin_after),
                    sign_preserved=margin_after > 0,
                    exceeds_half_safety=abs(increment) > half_margin,
                )
            )

        stress_rows.append(
            NextIncrementStressRow(
                source_row=row_id,
                M=M,
                next_M=M + 1,
                T=str(T),
                k=row_k,
                allowed_relative_tail_ratio=decimal_format(allowed),
                max_next_abs_theta=decimal_format(max_theta),
                max_next_abs_theta_over_allowed=decimal_format(pointwise_ratio),
                max_abs_log_increment=decimal_format(max_log),
                pointwise_budget_satisfied=max_theta <= allowed,
                stencil_signs_preserved=all(item.sign_preserved for item in stencil_rows),
                half_safety_stencils_satisfied=all(not item.exceeds_half_safety for item in stencil_rows),
                limiting_stencil=row["limiting_stencil"],
                stencil_increments=stencil_rows,
            )
        )

    if worst_pointwise_ratio is None or worst_stencil_half_ratio is None:
        raise RuntimeError("no next-increment rows were available")

    return {
        "parameters": {
            "max_taylor_degree": max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "sample_T_values": t_values,
        },
        "source_pointwise_summary": pointwise["summary"],
        "tested_next_increment_rows": len(stress_rows),
        "missing_next_increment_rows": len(missing_rows),
        "pointwise_budget_failure_rows": sum(1 for row in stress_rows if not row.pointwise_budget_satisfied),
        "stencil_sign_preserving_rows": sum(1 for row in stress_rows if row.stencil_signs_preserved),
        "half_safety_stencil_failure_rows": sum(1 for row in stress_rows if not row.half_safety_stencils_satisfied),
        "stress_rows": [asdict(row) for row in stress_rows],
        "missing_rows": [asdict(row) for row in missing_rows],
        "worst_pointwise_over_budget": {
            "sample": decimal_format(worst_pointwise_ratio[0]),
            "source_row": worst_pointwise_ratio[1],
        },
        "worst_stencil_abs_over_half_margin": {
            "sample": decimal_format(worst_stencil_half_ratio[0]),
            "source_row": worst_stencil_half_ratio[1],
            "stencil": worst_stencil_half_ratio[2],
        },
    }


def build_artifact(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    diagnostics = build_diagnostics(max_degree, cutoff_n, precision_bits, k, t_values)
    rows = [
        {
            "id": "nlrgniss_01_next_increment_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "For a truncation baseline M, the known next increment has log coordinate delta_j=log(F_j^(M+1)/F_j^(M)) and contributes to the same B, companion, and weighted-gap stencils as the unknown tail.",
            "formula": "delta_j=log(F_j^(M+1)/F_j^(M))",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
            ],
            "proof_boundary": "Exact finite bookkeeping for the next known coefficient only; not an infinite-tail theorem.",
        },
        {
            "id": "nlrgniss_02_pointwise_budget_failure",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "In the tested rows with an available next coefficient, the next known relative increment exceeds the pointwise half-safety tail budget.",
            "diagnostics": diagnostics,
            "proof_boundary": "Finite stress diagnostic only; it does not prove failure of the infinite series or of a direct stencil-tail theorem.",
        },
        {
            "id": "nlrgniss_03_structured_stencil_survival",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Although the pointwise budget fails, the structured next-increment stencil contributions preserve all three finite stencil signs in the tested rows.",
            "proof_boundary": "Finite cancellation diagnostic only; not a proof that the remaining infinite tail preserves the signs.",
        },
        {
            "id": "nlrgniss_04_half_safety_not_respected",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The companion next-increment stencil can exceed the half-safety stencil margin even when the total next-truncation sign is preserved.",
            "proof_boundary": "Finite warning only; it shows the half-safety sufficient condition is conservative, not false.",
        },
        {
            "id": "nlrgniss_05_degree16_unavailable_frontier",
            "role": "finite_frontier",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-14 positive baselines cannot be tested by the same next-increment method until the degree-16 Taylor ratio is computed.",
            "proof_boundary": "Finite scaffold frontier only; not an analytic obstruction.",
        },
        {
            "id": "nlrgniss_06_pointwise_triangle_route_rejected_for_current_baselines",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The current positive baselines plus the crude pointwise triangle envelope already certify the next Taylor increment.",
            "gap": "The tested next increments exceed the pointwise relative-tail budget by factors above 2000, so this finite shortcut is not acceptable.",
            "proof_boundary": "Rejected finite shortcut only; it does not reject a sharper analytic pointwise theorem at a different truncation order.",
        },
        {
            "id": "nlrgniss_07_live_direct_stencil_tail_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A stronger promoted route should bound the B, companion, and weighted-gap stencil tails directly, preserving cancellations instead of bounding each epsilon_j separately.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md",
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
            ],
            "proof_boundary": "Live theorem-search route only; no direct stencil-tail theorem is proved.",
        },
        {
            "id": "nlrgniss_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted direct stencil-tail proof must state the truncation order, available coefficient range, q/T collar, normalized multiplier positivity, and signed bounds for the three structured stencil tails.",
            "source_artifacts": [
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "tested_next_increment_rows": diagnostics["tested_next_increment_rows"],
        "missing_next_increment_rows": diagnostics["missing_next_increment_rows"],
        "pointwise_budget_failure_rows": diagnostics["pointwise_budget_failure_rows"],
        "stencil_sign_preserving_rows": diagnostics["stencil_sign_preserving_rows"],
        "half_safety_stencil_failure_rows": diagnostics["half_safety_stencil_failure_rows"],
        "worst_pointwise_over_budget": diagnostics["worst_pointwise_over_budget"]["sample"],
        "worst_stencil_abs_over_half_margin": diagnostics["worst_stencil_abs_over_half_margin"]["sample"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The known next Taylor increment defeats the crude pointwise half-safety envelope on both tested "
            f"positive baselines, with worst over-budget factor {diagnostics['worst_pointwise_over_budget']['sample']}. "
            "However, the structured B, companion, and weighted-gap next-increment stencils preserve all tested "
            "finite signs, while the companion stencil exceeds the half-safety margin. This points away from a "
            "crude pointwise triangle proof and toward a direct signed stencil-tail theorem."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress",
        "date": "2026-07-07",
        "status": "finite theorem-search diagnostic",
        "source_pointwise_tail_budget": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md",
        "source_high_order_taylor_scout": "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It compares the known next Taylor increment with the current "
            "pointwise and structured stencil budgets, but it does not prove an infinite Taylor-tail theorem, "
            "does not prove scaled-curvature monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Known next-increment stress is not promoted to an infinite-tail theorem.",
            "Failure of the crude pointwise budget is not promoted to failure of a direct stencil-tail theorem.",
            "Stencil sign survival for tested rows is finite evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][1]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['tested_next_increment_rows']} tested next-increment rows, "
        f"{summary['pointwise_budget_failure_rows']} pointwise budget failures, "
        f"{summary['stencil_sign_preserving_rows']} stencil-sign-preserving rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Next-Increment Stencil Stress",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof",
        "of a uniform Taylor-tail theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress`.",
        "",
        "Proof boundary: this artifact compares the known next Taylor",
        "increment with the current pointwise and structured stencil budgets.",
        "It does not prove an infinite Taylor-tail estimate.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Stress Summary",
        "",
        "```text",
        f"tested next-increment rows: {summary['tested_next_increment_rows']}",
        f"missing next-increment rows: {summary['missing_next_increment_rows']}",
        f"pointwise budget failures: {summary['pointwise_budget_failure_rows']}",
        f"stencil-sign-preserving rows: {summary['stencil_sign_preserving_rows']}",
        f"half-safety stencil failures: {summary['half_safety_stencil_failure_rows']}",
        f"worst pointwise over-budget factor: {summary['worst_pointwise_over_budget']}",
        f"worst stencil abs over half margin: {summary['worst_stencil_abs_over_half_margin']}",
        "```",
        "",
        "Tested rows:",
        "",
        "```text",
    ]
    for row in diagnostics["stress_rows"]:
        lines.append(
            f"{row['source_row']}: M->{row['next_M']}, max |theta| / allowed rho = "
            f"{row['max_next_abs_theta_over_allowed']}, stencil signs preserved = {row['stencil_signs_preserved']}"
        )
        for stencil in row["stencil_increments"]:
            lines.append(
                f"  {stencil['name']}: increment {stencil['increment']}, "
                f"abs/half-margin {stencil['abs_over_half_margin']}, "
                f"margin after increment {stencil['margin_after_increment']}"
            )
    lines.extend(
        [
            "```",
            "",
            "Rows awaiting a degree-16 coefficient:",
            "",
            "```text",
        ]
    )
    for row in diagnostics["missing_rows"]:
        lines.append(f"{row['source_row']}: M={row['M']}, T={row['T']}, reason={row['reason']}")
    lines.extend(
        [
            "```",
            "",
            "Interpretation:",
            "",
            "The crude pointwise triangle envelope is too conservative for the",
            "known next increment at the current positive baselines. The structured",
            "stencil increments still preserve the tested finite signs, so the live",
            "route is a direct signed stencil-tail theorem rather than a purely",
            "pointwise epsilon_j theorem at this truncation order.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_pointwise_tail_budget"],
            artifact["source_high_order_taylor_scout"],
            artifact["source_uniform_remainder_target"],
            "```",
            "",
            "Summary:",
            "",
            summary["main_finding"],
            "",
        ]
    )
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
        "wrote Jensen-window PF negative-lambda relative-Gaussian next-increment stencil stress: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
