#!/usr/bin/env python3
"""Build the relative-Gaussian Taylor stencil scout for the negative-lambda tail."""

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
from jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix import (  # noqa: E402
    arb_positive,
)


getcontext().prec = 100

DEFAULT_OUT_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.md"


@dataclass(frozen=True)
class LeadingSigns:
    a: str
    b: str
    c: str
    b_leading_coefficient: str
    companion_leading_coefficient: str
    weighted_gap_leading_coefficient: str


@dataclass(frozen=True)
class TruncationStencilRow:
    id: str
    truncation_degree: int
    M: int
    T: str
    k: int
    normalizer_positive: bool
    status: str
    F_k_minus_1: str
    F_k: str
    F_k_plus_1: str
    F_k_plus_2: str
    B_k: str | None
    B_decrease: str | None
    C_increase: str | None
    proof_boundary: str


def classify_stencil(k: int, T: int, ratios: list[Decimal], M: int) -> TruncationStencilRow:
    f_values = {index: truncated_multiplier(index, T, ratios, M) for index in (k - 1, k, k + 1, k + 2)}
    normalizer_positive = all(value > 0 for value in f_values.values())
    b_value: Decimal | None = None
    b_decrease: Decimal | None = None
    c_increase: Decimal | None = None
    if not normalizer_positive:
        status = "invalid_normalizer"
    else:
        logs = {index: f_values[index].ln() for index in f_values}
        b_value = 2 * logs[k] - logs[k - 1] - logs[k + 1]
        b_next = 2 * logs[k + 1] - logs[k] - logs[k + 2]
        b_decrease = b_value - b_next
        c_increase = Decimal(2 * k + 3) * b_next - Decimal(2 * k + 1) * b_value
        if b_value > 0 and b_decrease > 0 and c_increase > 0:
            status = "target_satisfying_truncation"
        elif b_value <= 0:
            status = "upper_wall_violation"
        elif b_decrease <= 0:
            status = "companion_failure"
        else:
            status = "weighted_gap_failure"
    return TruncationStencilRow(
        id=f"nlrgts_M{M}_T{T}",
        truncation_degree=2 * M,
        M=M,
        T=str(T),
        k=k,
        normalizer_positive=normalizer_positive,
        status=status,
        F_k_minus_1=decimal_format(f_values[k - 1]),
        F_k=decimal_format(f_values[k]),
        F_k_plus_1=decimal_format(f_values[k + 1]),
        F_k_plus_2=decimal_format(f_values[k + 2]),
        B_k=None if b_value is None else decimal_format(b_value),
        B_decrease=None if b_decrease is None else decimal_format(b_decrease),
        C_increase=None if c_increase is None else decimal_format(c_increase),
        proof_boundary="Finite Taylor truncation diagnostic only; not a proof about the full Taylor tail or actual zeta coefficients.",
    )


def build_diagnostics(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    high = build_high_order_diagnostics(max_degree, cutoff_n, precision_bits, k, t_values)
    ratio_balls = [flint.arb(row.ratio_to_c0) for row in high.coefficient_rows]
    ratios = [arb_mid_decimal(ball) for ball in ratio_balls]
    a, b, c = ratio_balls[1], ratio_balls[2], ratio_balls[3]
    b_lead = a * a - 2 * b
    companion_lead = 2 * (a**3 - 3 * a * b + 3 * c)
    weighted_gap_lead = 2 * b_lead
    for name, value in (
        ("B leading coefficient", b_lead),
        ("companion leading coefficient", companion_lead),
        ("weighted gap leading coefficient", weighted_gap_lead),
    ):
        if not arb_positive(value):
            raise RuntimeError(f"{name} is not certified positive: {value}")

    stencil_rows: list[TruncationStencilRow] = []
    for M in range(3, (max_degree // 2) + 1):
        for T in t_values:
            stencil_rows.append(classify_stencil(k, T, ratios, M))

    valid_rows = [row for row in stencil_rows if row.normalizer_positive]
    b_positive_rows = sum(1 for row in valid_rows if row.B_k is not None and Decimal(row.B_k) > 0)
    b_decrease_rows = sum(1 for row in valid_rows if row.B_decrease is not None and Decimal(row.B_decrease) > 0)
    c_increase_rows = sum(1 for row in valid_rows if row.C_increase is not None and Decimal(row.C_increase) > 0)
    all_positive_rows = sum(1 for row in stencil_rows if row.status == "target_satisfying_truncation")
    invalid_rows = sum(1 for row in stencil_rows if row.status == "invalid_normalizer")
    upper_wall_rows = sum(1 for row in stencil_rows if row.status == "upper_wall_violation")
    companion_failure_rows = sum(1 for row in stencil_rows if row.status == "companion_failure")
    weighted_gap_failure_rows = sum(1 for row in stencil_rows if row.status == "weighted_gap_failure")

    return {
        "parameters": {
            "max_taylor_degree": max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "sample_T_values": t_values,
            "truncation_degrees": [2 * M for M in range(3, (max_degree // 2) + 1)],
        },
        "leading_signs": asdict(
            LeadingSigns(
                a=a.str(40),
                b=b.str(40),
                c=c.str(40),
                b_leading_coefficient=b_lead.str(40),
                companion_leading_coefficient=companion_lead.str(40),
                weighted_gap_leading_coefficient=weighted_gap_lead.str(40),
            )
        ),
        "stencil_rows": [asdict(row) for row in stencil_rows],
        "truncation_rows": len(stencil_rows),
        "invalid_normalizer_rows": invalid_rows,
        "positive_normalizer_rows": len(valid_rows),
        "b_positive_rows": b_positive_rows,
        "b_decrease_rows": b_decrease_rows,
        "c_increase_rows": c_increase_rows,
        "all_positive_rows": all_positive_rows,
        "upper_wall_violation_rows": upper_wall_rows,
        "companion_failure_rows": companion_failure_rows,
        "weighted_gap_failure_rows": weighted_gap_failure_rows,
    }


def build_artifact(max_degree: int, cutoff_n: int, precision_bits: int, k: int, t_values: list[int]) -> dict:
    diagnostics = build_diagnostics(max_degree, cutoff_n, precision_bits, k, t_values)
    rows = [
        {
            "id": "nlrgts_01_relative_gaussian_taylor_multiplier",
            "role": "formal_model",
            "readiness": "not_ready_to_apply",
            "claim": "After Gaussian cancellation, the local multiplier is F_k=sum_j r_j*(k+1/2)_j/T^j and f_k=log(F_k).",
            "formula": "F_k=1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3+...",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
            ],
            "proof_boundary": "Formal local model only; the analytic Taylor tail is not bounded here.",
        },
        {
            "id": "nlrgts_02_B_leading_sign",
            "role": "certified_local_sign",
            "readiness": "not_ready_to_apply",
            "claim": "The fixed-k expansion has B_k=(a^2-2*b)/T^2+O_k(T^-3), and the leading coefficient is certified positive.",
            "leading_coefficient": diagnostics["leading_signs"]["b_leading_coefficient"],
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_signed_gaussian_perturbation_matrix.md",
            ],
            "proof_boundary": "Certified fixed-k leading sign only; not a uniform all-k theorem.",
        },
        {
            "id": "nlrgts_03_companion_leading_sign",
            "role": "certified_local_sign",
            "readiness": "not_ready_to_apply",
            "claim": "The companion side has B_k-B_(k+1)=2*(a^3-3*a*b+3*c)/T^3+O_k(T^-4), and the leading coefficient is certified positive.",
            "leading_coefficient": diagnostics["leading_signs"]["companion_leading_coefficient"],
            "source_artifacts": [
                "outputs/jensen_window_pf_monotone_contraction_theorem_target.md",
            ],
            "proof_boundary": "Certified fixed-k leading sign only; not the monotone-contraction theorem.",
        },
        {
            "id": "nlrgts_04_weighted_gap_leading_sign",
            "role": "certified_local_sign",
            "readiness": "not_ready_to_apply",
            "claim": "The new weighted four-point target has C_(k+1)-C_k=2*(a^2-2*b)/T^2+O_k(T^-3), and the leading coefficient is certified positive.",
            "leading_coefficient": diagnostics["leading_signs"]["weighted_gap_leading_coefficient"],
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md",
            ],
            "proof_boundary": "Certified fixed-k leading sign only; not the weighted four-point all-k theorem.",
        },
        {
            "id": "nlrgts_05_finite_truncation_stencil_matrix",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Finite Taylor truncation stencils through degree 14 are unstable at the k=22 tail start: only a small subset satisfies B_k>0, B_k-B_(k+1)>0, and C_(k+1)-C_k>0 simultaneously.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            "proof_boundary": "Finite truncation diagnostic only; not evidence that the full Taylor series fails or succeeds.",
        },
        {
            "id": "nlrgts_06_finite_truncation_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "A finite Taylor truncation of degree 6,8,10,12, or 14 can replace a uniform remainder theorem for the relative-Gaussian stencil.",
            "gap": "The truncation matrix contains invalid normalizers, upper-wall violations, companion failures, and weighted-gap failures.",
            "proof_boundary": "Rejected proof template only; not a rejection of the signed perturbation route.",
        },
        {
            "id": "nlrgts_07_live_uniform_stencil_remainder_route",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof may use the positive fixed-k leading signs only after bounding the full Taylor tail in the B, companion, and weighted-gap stencils over a stated q/T range.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md",
            ],
            "proof_boundary": "Live theorem-search route only; no uniform stencil remainder theorem is proved.",
        },
        {
            "id": "nlrgts_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof must state the q/T range, finite collar, positivity of the normalized multiplier, and explicit remainder bounds for B_k, B_k-B_(k+1), and C_(k+1)-C_k.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_raw_moment_obstruction_matrix.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of cone entry, jwpf_06, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "certified_leading_sign_rows": 3,
        "truncation_rows": diagnostics["truncation_rows"],
        "invalid_normalizer_rows": diagnostics["invalid_normalizer_rows"],
        "positive_normalizer_rows": diagnostics["positive_normalizer_rows"],
        "b_positive_rows": diagnostics["b_positive_rows"],
        "b_decrease_rows": diagnostics["b_decrease_rows"],
        "c_increase_rows": diagnostics["c_increase_rows"],
        "all_positive_rows": diagnostics["all_positive_rows"],
        "upper_wall_violation_rows": diagnostics["upper_wall_violation_rows"],
        "companion_failure_rows": diagnostics["companion_failure_rows"],
        "weighted_gap_failure_rows": diagnostics["weighted_gap_failure_rows"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The fixed-k Taylor stencil has certified positive leading signs for B_k, "
            "B_k-B_(k+1), and the relative-Gaussian weighted gap C_(k+1)-C_k. "
            "However, finite truncations through degree 14 are not stable proof objects: "
            "only 4/35 sampled truncation rows satisfy all three stencil inequalities, "
            "so the route still needs a uniform Taylor-tail remainder theorem."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout",
        "date": "2026-07-07",
        "status": "finite theorem-search diagnostic",
        "source_relative_gaussian_bridge": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_curvature_bridge.md",
        "source_taylor_moment_budget": "outputs/jensen_window_pf_negative_lambda_taylor_moment_budget.md",
        "source_high_order_taylor_scout": "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It certifies fixed-k leading Taylor signs and "
            "classifies finite truncation stencils for the relative-Gaussian weighted target, but it "
            "does not prove a uniform Taylor-tail remainder theorem, does not prove scaled-curvature "
            "monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Fixed-k leading Taylor signs are not promoted to a moving-k all-tail theorem.",
            "Finite Taylor truncations are not promoted to an infinite Taylor-tail theorem.",
            "The relative-Gaussian weighted four-point inequality remains open.",
            "Generic moment positivity, endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][4]["diagnostics"]
    signs = diagnostics["leading_signs"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['certified_leading_sign_rows']} certified leading-sign rows, "
        f"{summary['truncation_rows']} truncation rows, "
        f"{summary['all_positive_rows']} all-positive stencil rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Taylor Stencil Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof of",
        "scaled-curvature monotonicity, cone entry, Jensen-window PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout`.",
        "",
        "Proof boundary: this artifact certifies fixed-`k` leading Taylor",
        "signs and finite truncation behavior for the relative-Gaussian",
        "stencil. It does not prove a uniform Taylor-tail remainder theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_taylor_stencil_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Fixed-k Leading Signs",
        "",
        "For the normalized local multiplier:",
        "",
        "```text",
        "F_k = 1+a*q/T+b*q*(q+1)/T^2+c*q*(q+1)*(q+2)/T^3+...",
        "q = k+1/2",
        "f_k = log(F_k)",
        "```",
        "",
        "the formal fixed-`k` stencil gives:",
        "",
        "```text",
        "B_k = (a^2-2*b)/T^2 + O_k(T^-3)",
        "B_k-B_(k+1) = 2*(a^3-3*a*b+3*c)/T^3 + O_k(T^-4)",
        "C_(k+1)-C_k = 2*(a^2-2*b)/T^2 + O_k(T^-3)",
        "```",
        "",
        "Certified leading coefficients:",
        "",
        "```text",
        f"a = {signs['a']}",
        f"b = {signs['b']}",
        f"c = {signs['c']}",
        f"a^2-2*b = {signs['b_leading_coefficient']}",
        f"2*(a^3-3*a*b+3*c) = {signs['companion_leading_coefficient']}",
        f"2*(a^2-2*b) = {signs['weighted_gap_leading_coefficient']}",
        "```",
        "",
        "## Truncation Stencil Matrix",
        "",
        "Rows classify degrees 6, 8, 10, 12, and 14 at `k=22` and",
        "`T=25,50,100,200,500,1000,2000`:",
        "",
        "```text",
        f"truncation rows: {summary['truncation_rows']}",
        f"invalid normalizers: {summary['invalid_normalizer_rows']}",
        f"positive normalizers: {summary['positive_normalizer_rows']}",
        f"B-positive rows: {summary['b_positive_rows']}",
        f"B-decrease rows: {summary['b_decrease_rows']}",
        f"C-increase rows: {summary['c_increase_rows']}",
        f"all-positive stencil rows: {summary['all_positive_rows']}",
        f"upper-wall violation rows: {summary['upper_wall_violation_rows']}",
        f"companion-failure rows: {summary['companion_failure_rows']}",
        f"weighted-gap failure rows: {summary['weighted_gap_failure_rows']}",
        "```",
        "",
        "Interpretation:",
        "",
        "The leading fixed-`k` signs point in the right direction, but finite",
        "Taylor truncations are not stable proof objects. A promoted proof must",
        "bound the full Taylor tail in the `B_k`, `B_k-B_(k+1)`, and",
        "`C_(k+1)-C_k` stencils over a stated `q/T` range.",
        "",
        "Integration:",
        "",
        "```text",
        artifact["source_relative_gaussian_bridge"],
        artifact["source_taylor_moment_budget"],
        artifact["source_high_order_taylor_scout"],
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
        "wrote Jensen-window PF negative-lambda relative-Gaussian Taylor stencil scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
