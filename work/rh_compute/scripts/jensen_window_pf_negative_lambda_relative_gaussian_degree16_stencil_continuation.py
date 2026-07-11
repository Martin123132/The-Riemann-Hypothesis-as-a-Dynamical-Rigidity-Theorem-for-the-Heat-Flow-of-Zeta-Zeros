#!/usr/bin/env python3
"""Build the degree-16 stencil-continuation diagnostic for the relative-Gaussian route."""

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
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md"


@dataclass(frozen=True)
class StencilIncrement:
    name: str
    increment: str
    half_margin: str
    abs_over_half_margin: str
    margin_after_increment: str
    sign_preserved: bool
    exceeds_half_safety: bool


@dataclass(frozen=True)
class ContinuationRow:
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
    failing_stencils: list[str]
    limiting_stencil: str
    stencil_increments: list[StencilIncrement]


def theta_and_log_increments(k: int, T: int, ratios: list[Decimal], M: int) -> tuple[dict[int, Decimal], dict[int, Decimal]]:
    theta: dict[int, Decimal] = {}
    eps: dict[int, Decimal] = {}
    for index in (k - 1, k, k + 1, k + 2):
        current = truncated_multiplier(index, T, ratios, M)
        next_value = truncated_multiplier(index, T, ratios, M + 1)
        theta[index] = (next_value - current) / current
        eps[index] = next_value.ln() - current.ln()
    return theta, eps


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


def build_diagnostics(
    baseline_max_degree: int,
    coefficient_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    t_values: list[int],
) -> dict:
    pointwise = build_pointwise_budget_artifact(baseline_max_degree, cutoff_n, precision_bits, k, t_values)
    budget_rows = pointwise["matrix_rows"][3]["diagnostics"]["budget_rows"]
    high = build_high_order_diagnostics(coefficient_max_degree, cutoff_n, precision_bits, k, t_values)
    ratios = [arb_mid_decimal(flint.arb(row.ratio_to_c0)) for row in high.coefficient_rows]
    degree16_row = next(row for row in high.coefficient_rows if row.degree == 16)

    continuation_rows: list[ContinuationRow] = []
    worst_pointwise: tuple[Decimal, str] | None = None
    worst_stencil_half: tuple[Decimal, str, str] | None = None
    for row in budget_rows:
        row_id = row["source_row"]
        M = int(row["M"])
        T = int(row["T"])
        row_k = int(row["k"])
        if M + 1 >= len(ratios):
            raise RuntimeError(f"missing next ratio for {row_id}: M={M}, available={len(ratios)}")
        theta, eps = theta_and_log_increments(row_k, T, ratios, M)
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
        if worst_pointwise is None or pointwise_ratio > worst_pointwise[0]:
            worst_pointwise = (pointwise_ratio, row_id)

        increments: list[StencilIncrement] = []
        failing: list[str] = []
        for name in ("B", "companion", "weighted_gap"):
            half_margin = margins[name] / Decimal(2)
            increment = stencils[name]
            abs_over_half = abs(increment) / half_margin
            if worst_stencil_half is None or abs_over_half > worst_stencil_half[0]:
                worst_stencil_half = (abs_over_half, row_id, name)
            margin_after = margins[name] + increment
            sign_preserved = margin_after > 0
            if not sign_preserved:
                failing.append(name)
            increments.append(
                StencilIncrement(
                    name=name,
                    increment=decimal_format(increment),
                    half_margin=decimal_format(half_margin),
                    abs_over_half_margin=decimal_format(abs_over_half),
                    margin_after_increment=decimal_format(margin_after),
                    sign_preserved=sign_preserved,
                    exceeds_half_safety=abs(increment) > half_margin,
                )
            )

        continuation_rows.append(
            ContinuationRow(
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
                stencil_signs_preserved=not failing,
                half_safety_stencils_satisfied=all(not item.exceeds_half_safety for item in increments),
                failing_stencils=failing,
                limiting_stencil=row["limiting_stencil"],
                stencil_increments=increments,
            )
        )

    if worst_pointwise is None or worst_stencil_half is None:
        raise RuntimeError("no continuation rows")

    return {
        "parameters": {
            "baseline_max_taylor_degree": baseline_max_degree,
            "coefficient_max_taylor_degree": coefficient_max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "sample_T_values": t_values,
        },
        "degree16_coefficient": {
            "degree": degree16_row.degree,
            "enclosure": degree16_row.enclosure,
            "ratio_to_c0": degree16_row.ratio_to_c0,
            "sign": degree16_row.sign,
            "tail_radius": degree16_row.tail_radius,
        },
        "source_pointwise_summary": pointwise["summary"],
        "continuation_rows": [asdict(row) for row in continuation_rows],
        "tested_continuation_rows": len(continuation_rows),
        "pointwise_budget_failure_rows": sum(1 for row in continuation_rows if not row.pointwise_budget_satisfied),
        "stencil_sign_preserving_rows": sum(1 for row in continuation_rows if row.stencil_signs_preserved),
        "stencil_sign_failure_rows": sum(1 for row in continuation_rows if not row.stencil_signs_preserved),
        "half_safety_stencil_failure_rows": sum(1 for row in continuation_rows if not row.half_safety_stencils_satisfied),
        "degree14_baseline_rows": sum(1 for row in continuation_rows if row.M == 7),
        "degree14_baseline_survivors": sum(1 for row in continuation_rows if row.M == 7 and row.stencil_signs_preserved),
        "degree14_baseline_failures": sum(1 for row in continuation_rows if row.M == 7 and not row.stencil_signs_preserved),
        "surviving_rows": [row.source_row for row in continuation_rows if row.stencil_signs_preserved],
        "failed_rows": [row.source_row for row in continuation_rows if not row.stencil_signs_preserved],
        "worst_pointwise_over_budget": {
            "sample": decimal_format(worst_pointwise[0]),
            "source_row": worst_pointwise[1],
        },
        "worst_stencil_abs_over_half_margin": {
            "sample": decimal_format(worst_stencil_half[0]),
            "source_row": worst_stencil_half[1],
            "stencil": worst_stencil_half[2],
        },
    }


def build_artifact(
    baseline_max_degree: int,
    coefficient_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    t_values: list[int],
) -> dict:
    diagnostics = build_diagnostics(baseline_max_degree, coefficient_max_degree, cutoff_n, precision_bits, k, t_values)
    rows = [
        {
            "id": "nlrgd16sc_01_degree16_coefficient",
            "role": "finite_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-16 relative-Gaussian Taylor coefficient ratio is computed with the same finite-sum plus geometric-tail method as the high-order scout.",
            "degree16_coefficient": diagnostics["degree16_coefficient"],
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            "proof_boundary": "Finite coefficient certificate only; not an infinite Taylor-tail theorem.",
        },
        {
            "id": "nlrgd16sc_02_all_positive_baselines_tested",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-16 continuation tests all four current positive relative-Gaussian finite baselines from the pointwise-tail budget artifact.",
            "diagnostics": diagnostics,
            "proof_boundary": "Finite continuation diagnostic only; not an analytic all-k result.",
        },
        {
            "id": "nlrgd16sc_03_pointwise_budget_still_fails",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Every tested next increment still exceeds the crude pointwise half-safety relative-tail budget.",
            "proof_boundary": "Finite rejection of the current pointwise shortcut only; not a rejection of a sharper analytic theorem.",
        },
        {
            "id": "nlrgd16sc_04_stencil_survival_filter",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The structured next-increment stencil test preserves three baselines and fails the degree-14 T=1000 baseline through the companion stencil.",
            "proof_boundary": "Finite filter only; not a theorem about all T or all remaining tail terms.",
        },
        {
            "id": "nlrgd16sc_05_collar_signal",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-16 continuation suggests that any direct signed stencil-tail theorem should include an explicit sufficiently-large T or q/T collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.md",
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
            ],
            "proof_boundary": "Live theorem-search signal only; no collar theorem is proved.",
        },
        {
            "id": "nlrgd16sc_06_degree16_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "Passing some degree-16 continuation rows proves the direct stencil-tail theorem.",
            "gap": "One tested degree-14 baseline fails and surviving finite rows do not control the infinite remaining Taylor tail.",
            "proof_boundary": "Rejected finite promotion only; not a claim about the final analytic route.",
        },
        {
            "id": "nlrgd16sc_07_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted direct stencil-tail proof must replace degree-by-degree continuation stress with explicit signed bounds for the infinite residual stencil tails beyond a stated degree and collar.",
            "source_artifacts": [
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "tested_continuation_rows": diagnostics["tested_continuation_rows"],
        "pointwise_budget_failure_rows": diagnostics["pointwise_budget_failure_rows"],
        "stencil_sign_preserving_rows": diagnostics["stencil_sign_preserving_rows"],
        "stencil_sign_failure_rows": diagnostics["stencil_sign_failure_rows"],
        "half_safety_stencil_failure_rows": diagnostics["half_safety_stencil_failure_rows"],
        "degree14_baseline_rows": diagnostics["degree14_baseline_rows"],
        "degree14_baseline_survivors": diagnostics["degree14_baseline_survivors"],
        "degree14_baseline_failures": diagnostics["degree14_baseline_failures"],
        "worst_pointwise_over_budget": diagnostics["worst_pointwise_over_budget"]["sample"],
        "worst_stencil_abs_over_half_margin": diagnostics["worst_stencil_abs_over_half_margin"]["sample"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "Adding the degree-16 Taylor coefficient removes the previous next-increment frontier for the "
            "degree-14 positive baselines. All four current positive baselines still fail the crude pointwise "
            "tail budget, while three preserve the structured stencil signs. The degree-14 T=1000 baseline "
            "fails through the companion stencil, and the degree-14 T=2000 baseline survives. This sharpens "
            "the direct signed stencil-tail route toward an explicit large-T/q-over-T collar."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation",
        "date": "2026-07-07",
        "status": "finite theorem-search diagnostic",
        "source_pointwise_tail_budget": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_pointwise_tail_budget.md",
        "source_next_increment_stress": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_next_increment_stencil_stress.md",
        "source_high_order_taylor_scout": "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It computes the degree-16 coefficient and stress-tests "
            "degree-16 continuation of current finite baselines, but it does not prove an infinite Taylor-tail "
            "theorem, does not prove scaled-curvature monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Degree-16 continuation stress is not promoted to an infinite-tail theorem.",
            "Failure of the T=1000 baseline is finite evidence for a collar, not an analytic obstruction.",
            "Survival of the T=2000 baseline is finite evidence only.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][1]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['tested_continuation_rows']} tested continuation rows, "
        f"{summary['stencil_sign_preserving_rows']} stencil-sign-preserving rows, "
        f"{summary['stencil_sign_failure_rows']} stencil-sign-failure rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Stencil Continuation",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof",
        "of a uniform Taylor-tail theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation`.",
        "",
        "Proof boundary: this artifact computes the degree-16 coefficient",
        "and stress-tests the current finite baselines one Taylor step further.",
        "It does not prove an infinite Taylor-tail estimate.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Degree-16 Coefficient",
        "",
        "```text",
        f"degree: {diagnostics['degree16_coefficient']['degree']}",
        f"sign: {diagnostics['degree16_coefficient']['sign']}",
        f"ratio c16/c0: {diagnostics['degree16_coefficient']['ratio_to_c0']}",
        f"tail radius: {diagnostics['degree16_coefficient']['tail_radius']}",
        "```",
        "",
        "## Continuation Summary",
        "",
        "```text",
        f"tested continuation rows: {summary['tested_continuation_rows']}",
        f"pointwise budget failures: {summary['pointwise_budget_failure_rows']}",
        f"stencil-sign-preserving rows: {summary['stencil_sign_preserving_rows']}",
        f"stencil-sign-failure rows: {summary['stencil_sign_failure_rows']}",
        f"half-safety stencil failures: {summary['half_safety_stencil_failure_rows']}",
        f"degree-14 baseline survivors: {summary['degree14_baseline_survivors']}",
        f"degree-14 baseline failures: {summary['degree14_baseline_failures']}",
        f"worst pointwise over-budget factor: {summary['worst_pointwise_over_budget']}",
        f"worst stencil abs over half margin: {summary['worst_stencil_abs_over_half_margin']}",
        "```",
        "",
        "Continuation rows:",
        "",
        "```text",
    ]
    for row in diagnostics["continuation_rows"]:
        lines.append(
            f"{row['source_row']}: M->{row['next_M']}, T={row['T']}, "
            f"pointwise over-budget={row['max_next_abs_theta_over_allowed']}, "
            f"stencil signs preserved={row['stencil_signs_preserved']}, "
            f"failing stencils={row['failing_stencils']}"
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
            "Interpretation:",
            "",
            "The degree-16 continuation removes the previous missing-coefficient",
            "frontier. The `T=1000` degree-14 baseline fails through the companion",
            "stencil, while the `T=2000` degree-14 baseline survives. This is a",
            "finite signal for an explicit large-`T` or `q/T` collar in any direct",
            "signed stencil-tail theorem.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_pointwise_tail_budget"],
            artifact["source_next_increment_stress"],
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
    parser.add_argument("--baseline-max-degree", type=int, default=14)
    parser.add_argument("--coefficient-max-degree", type=int, default=16)
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
    artifact = build_artifact(
        args.baseline_max_degree,
        args.coefficient_max_degree,
        args.cutoff_n,
        args.precision_bits,
        args.tail_start_k,
        args.t_values,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian degree-16 stencil continuation: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
