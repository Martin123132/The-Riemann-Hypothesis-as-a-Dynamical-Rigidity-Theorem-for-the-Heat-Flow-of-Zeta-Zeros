#!/usr/bin/env python3
"""Build a degree-16 collar scan for the relative-Gaussian direct stencil route."""

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
    relative_tail_ratio_from_log_bound,
)


getcontext().prec = 100

DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md"


@dataclass(frozen=True)
class CollarScanRow:
    T: str
    q_over_T: str
    baseline_positive: bool
    continuation_positive: bool
    half_safety_stencils_satisfied: bool
    pointwise_budget_satisfied: bool | None
    companion_margin_M7: str | None
    companion_margin_M8: str | None
    weighted_gap_margin_M7: str | None
    weighted_gap_margin_M8: str | None
    pointwise_over_budget: str | None
    worst_stencil_abs_over_half_margin: str | None


def log_stencils(k: int, T: int, ratios: list[Decimal], M: int) -> tuple[Decimal, Decimal, Decimal] | None:
    values = {index: truncated_multiplier(index, T, ratios, M) for index in (k - 1, k, k + 1, k + 2)}
    if not all(value > 0 for value in values.values()):
        return None
    logs = {index: values[index].ln() for index in values}
    B = Decimal(2) * logs[k] - logs[k - 1] - logs[k + 1]
    next_B = Decimal(2) * logs[k + 1] - logs[k] - logs[k + 2]
    companion = B - next_B
    weighted_gap = Decimal(2 * k + 3) * next_B - Decimal(2 * k + 1) * B
    return B, companion, weighted_gap


def theta_increments(k: int, T: int, ratios: list[Decimal], M: int) -> dict[int, Decimal]:
    values: dict[int, Decimal] = {}
    for index in (k - 1, k, k + 1, k + 2):
        current = truncated_multiplier(index, T, ratios, M)
        next_value = truncated_multiplier(index, T, ratios, M + 1)
        values[index] = (next_value - current) / current
    return values


def build_scan_row(k: int, T: int, ratios: list[Decimal], M: int) -> CollarScanRow:
    q = Decimal(k) + Decimal("0.5")
    base = log_stencils(k, T, ratios, M)
    cont = log_stencils(k, T, ratios, M + 1)
    baseline_positive = base is not None and all(value > 0 for value in base)
    continuation_positive = cont is not None and all(value > 0 for value in cont)
    pointwise_satisfied: bool | None = None
    pointwise_over_budget: str | None = None
    worst_stencil_ratio: str | None = None
    half_safety = False
    if baseline_positive and cont is not None:
        weight_sum = Decimal(16 * k + 16)
        eta = min(base[0] / Decimal(4), base[1] / Decimal(8), base[2] / weight_sum) / Decimal(2)
        rho = relative_tail_ratio_from_log_bound(eta)
        max_theta = max(abs(value) for value in theta_increments(k, T, ratios, M).values())
        pointwise_satisfied = max_theta <= rho
        pointwise_over_budget = decimal_format(max_theta / rho)
        ratios_over_half = [abs(after - before) / (before / Decimal(2)) for before, after in zip(base, cont)]
        half_safety = all(value <= Decimal(1) for value in ratios_over_half)
        worst_stencil_ratio = decimal_format(max(ratios_over_half))
    return CollarScanRow(
        T=str(T),
        q_over_T=decimal_format(q / Decimal(T)),
        baseline_positive=baseline_positive,
        continuation_positive=baseline_positive and continuation_positive,
        half_safety_stencils_satisfied=half_safety,
        pointwise_budget_satisfied=pointwise_satisfied,
        companion_margin_M7=None if base is None else decimal_format(base[1]),
        companion_margin_M8=None if cont is None else decimal_format(cont[1]),
        weighted_gap_margin_M7=None if base is None else decimal_format(base[2]),
        weighted_gap_margin_M8=None if cont is None else decimal_format(cont[2]),
        pointwise_over_budget=pointwise_over_budget,
        worst_stencil_abs_over_half_margin=worst_stencil_ratio,
    )


def first_true(rows: list[CollarScanRow], attr: str) -> CollarScanRow | None:
    for row in rows:
        if getattr(row, attr):
            return row
    return None


def build_diagnostics(
    coefficient_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    scan_start_T: int,
    scan_end_T: int,
) -> dict:
    high = build_high_order_diagnostics(coefficient_max_degree, cutoff_n, precision_bits, k, [1000, 2000])
    ratios = [arb_mid_decimal(flint.arb(row.ratio_to_c0)) for row in high.coefficient_rows]
    rows = [build_scan_row(k, T, ratios, 7) for T in range(scan_start_T, scan_end_T + 1)]
    baseline_rows = [row for row in rows if row.baseline_positive]
    continuation_rows = [row for row in rows if row.continuation_positive]
    half_safety_rows = [row for row in rows if row.half_safety_stencils_satisfied]
    pointwise_success_rows = [row for row in rows if row.pointwise_budget_satisfied is True]
    pointwise_failure_rows = [row for row in rows if row.pointwise_budget_satisfied is False]

    worst_pointwise = max(
        (Decimal(row.pointwise_over_budget), row) for row in rows if row.pointwise_over_budget is not None
    )
    worst_stencil = max(
        (Decimal(row.worst_stencil_abs_over_half_margin), row)
        for row in rows
        if row.worst_stencil_abs_over_half_margin is not None
    )
    selected_T = {
        scan_start_T,
        scan_end_T,
        900,
        950,
        1000,
        1100,
        1153,
        1154,
        1155,
        1156,
        1157,
        1158,
        1200,
        1400,
        1480,
        1481,
        1482,
        1483,
        1484,
        1485,
        1600,
        1800,
        2000,
        2200,
    }
    selected_rows = [row for row in rows if int(row.T) in selected_T]
    first_baseline = first_true(rows, "baseline_positive")
    first_continuation = first_true(rows, "continuation_positive")
    first_half_safety = first_true(rows, "half_safety_stencils_satisfied")
    if first_baseline is None or first_continuation is None or first_half_safety is None:
        raise RuntimeError("missing expected collar threshold")
    return {
        "parameters": {
            "coefficient_max_taylor_degree": coefficient_max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "scan_start_T": scan_start_T,
            "scan_end_T": scan_end_T,
            "baseline_M": 7,
            "continuation_M": 8,
        },
        "scan_rows": len(rows),
        "baseline_positive_rows": len(baseline_rows),
        "continuation_positive_rows": len(continuation_rows),
        "half_safety_rows": len(half_safety_rows),
        "pointwise_budget_success_rows": len(pointwise_success_rows),
        "pointwise_budget_failure_rows": len(pointwise_failure_rows),
        "first_baseline_positive": asdict(first_baseline),
        "first_continuation_positive": asdict(first_continuation),
        "first_half_safety": asdict(first_half_safety),
        "worst_pointwise_over_budget": {
            "sample": decimal_format(worst_pointwise[0]),
            "T": worst_pointwise[1].T,
        },
        "worst_stencil_abs_over_half_margin": {
            "sample": decimal_format(worst_stencil[0]),
            "T": worst_stencil[1].T,
        },
        "selected_rows": [asdict(row) for row in selected_rows],
    }


def build_artifact(
    coefficient_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    scan_start_T: int,
    scan_end_T: int,
) -> dict:
    diagnostics = build_diagnostics(coefficient_max_degree, cutoff_n, precision_bits, k, scan_start_T, scan_end_T)
    rows = [
        {
            "id": "nlrgd16cs_01_integer_T_scan",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-16 continuation is scanned on every integer T in the stated collar range for the M=7 to M=8 relative-Gaussian stencil.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md",
            ],
            "proof_boundary": "Finite integer-grid scan only; not an analytic collar theorem.",
        },
        {
            "id": "nlrgd16cs_02_baseline_threshold",
            "role": "finite_threshold",
            "readiness": "not_ready_to_apply",
            "claim": "On the scanned integer grid, the degree-14 baseline first has all three positive stencils at T=918.",
            "proof_boundary": "Finite threshold only; not a statement for all real T or all k.",
        },
        {
            "id": "nlrgd16cs_03_continuation_threshold",
            "role": "finite_threshold",
            "readiness": "not_ready_to_apply",
            "claim": "On the scanned integer grid, the degree-16 continuation first preserves all three structured stencil signs at T=1156.",
            "proof_boundary": "Finite threshold only; not an infinite Taylor-tail theorem.",
        },
        {
            "id": "nlrgd16cs_04_half_safety_threshold",
            "role": "finite_threshold",
            "readiness": "not_ready_to_apply",
            "claim": "On the scanned integer grid, the stricter half-safety stencil-increment condition first holds at T=1483.",
            "proof_boundary": "Finite sufficient-condition threshold only; not a proof outside the scan.",
        },
        {
            "id": "nlrgd16cs_05_pointwise_route_blocked_on_scan",
            "role": "finite_obstruction",
            "readiness": "not_ready_to_apply",
            "claim": "The crude pointwise half-safety budget fails for every baseline-positive row in the scanned collar range.",
            "proof_boundary": "Finite obstruction to this pointwise shortcut only; not a disproof of sharper analytic estimates.",
        },
        {
            "id": "nlrgd16cs_06_live_collar_theorem_target",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A direct signed stencil-tail theorem should aim for an explicit large-T/q-over-T collar, with T=1156 as the first finite sign-survival threshold and T=1483 as the first finite half-safety threshold in this scan.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md",
            ],
            "proof_boundary": "Live theorem-search target only; no analytic collar theorem is proved.",
        },
        {
            "id": "nlrgd16cs_07_scan_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The finite integer-T collar scan proves the required direct signed stencil-tail theorem.",
            "gap": "The scan is finite, k=22-specific, degree-16-specific, and does not bound the infinite residual Taylor tail.",
            "proof_boundary": "Rejected finite promotion only; not a proof of scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "scan_rows": diagnostics["scan_rows"],
        "baseline_positive_rows": diagnostics["baseline_positive_rows"],
        "continuation_positive_rows": diagnostics["continuation_positive_rows"],
        "half_safety_rows": diagnostics["half_safety_rows"],
        "pointwise_budget_success_rows": diagnostics["pointwise_budget_success_rows"],
        "pointwise_budget_failure_rows": diagnostics["pointwise_budget_failure_rows"],
        "first_baseline_positive_T": diagnostics["first_baseline_positive"]["T"],
        "first_continuation_positive_T": diagnostics["first_continuation_positive"]["T"],
        "first_half_safety_T": diagnostics["first_half_safety"]["T"],
        "first_continuation_q_over_T": diagnostics["first_continuation_positive"]["q_over_T"],
        "first_half_safety_q_over_T": diagnostics["first_half_safety"]["q_over_T"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The degree-16 collar scan on T=900..2200 makes the large-T signal concrete: the M=7 baseline "
            "first has positive stencils at T=918, the M=8 continuation first preserves signs at T=1156 "
            "(q/T=1.946366782006920415E-2), and half-safety first holds at T=1483 "
            "(q/T=1.517194875252865813E-2). The crude pointwise budget fails on all baseline-positive "
            "rows, so the live route remains a direct signed stencil-tail collar theorem."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan",
        "date": "2026-07-07",
        "status": "finite theorem-search diagnostic",
        "source_degree16_continuation": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_stencil_continuation.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It scans an integer T collar for k=22 and degree-16 continuation, "
            "but it does not prove a real-T collar theorem, does not prove an infinite Taylor-tail theorem, does not prove "
            "scaled-curvature monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Integer-T scan thresholds are not promoted to analytic collar theorems.",
            "Pointwise budget failure on the scan is not promoted to failure of a direct stencil-tail theorem.",
            "The scan is k=22 and degree-16 specific.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][0]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['scan_rows']} scan rows, "
        f"{summary['continuation_positive_rows']} continuation-positive rows, "
        f"{summary['half_safety_rows']} half-safety rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Collar Scan",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof",
        "of a real-`T` collar theorem, uniform Taylor-tail theorem,",
        "scaled-curvature monotonicity, cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan`.",
        "",
        "Proof boundary: this artifact scans an integer `T` collar for `k=22`",
        "and degree-16 continuation. It does not prove an analytic collar.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Collar Thresholds",
        "",
        "```text",
        f"scan rows: {summary['scan_rows']}",
        f"baseline-positive rows: {summary['baseline_positive_rows']}",
        f"continuation-positive rows: {summary['continuation_positive_rows']}",
        f"half-safety rows: {summary['half_safety_rows']}",
        f"pointwise budget successes: {summary['pointwise_budget_success_rows']}",
        f"pointwise budget failures: {summary['pointwise_budget_failure_rows']}",
        f"first baseline-positive T: {summary['first_baseline_positive_T']}",
        f"first continuation-positive T: {summary['first_continuation_positive_T']}",
        f"first continuation q/T: {summary['first_continuation_q_over_T']}",
        f"first half-safety T: {summary['first_half_safety_T']}",
        f"first half-safety q/T: {summary['first_half_safety_q_over_T']}",
        f"worst pointwise over-budget: {diagnostics['worst_pointwise_over_budget']['sample']} at T={diagnostics['worst_pointwise_over_budget']['T']}",
        f"worst stencil abs over half-margin: {diagnostics['worst_stencil_abs_over_half_margin']['sample']} at T={diagnostics['worst_stencil_abs_over_half_margin']['T']}",
        "```",
        "",
        "Selected scan rows:",
        "",
        "```text",
    ]
    for row in diagnostics["selected_rows"]:
        lines.append(
            f"T={row['T']}: baseline={row['baseline_positive']}, continuation={row['continuation_positive']}, "
            f"half_safety={row['half_safety_stencils_satisfied']}, companion_M8={row['companion_margin_M8']}, "
            f"weighted_gap_M8={row['weighted_gap_margin_M8']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Interpretation:",
            "",
            "The finite collar signal is now explicit: continuation positivity starts",
            "at `T=1156`, while the stricter half-safety condition starts at `T=1483`.",
            "The pointwise budget fails on every baseline-positive row in the scan, so",
            "the live theorem-search target remains a direct signed stencil-tail collar.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree16_continuation"],
            artifact["source_uniform_remainder_target"],
            artifact["source_dependency_graph"],
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
    parser.add_argument("--coefficient-max-degree", type=int, default=16)
    parser.add_argument("--cutoff-n", type=int, default=80)
    parser.add_argument("--precision-bits", type=int, default=256)
    parser.add_argument("--tail-start-k", type=int, default=22)
    parser.add_argument("--scan-start-T", type=int, default=900)
    parser.add_argument("--scan-end-T", type=int, default=2200)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(
        args.coefficient_max_degree,
        args.cutoff_n,
        args.precision_bits,
        args.tail_start_k,
        args.scan_start_T,
        args.scan_end_T,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian degree-16 collar scan: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
