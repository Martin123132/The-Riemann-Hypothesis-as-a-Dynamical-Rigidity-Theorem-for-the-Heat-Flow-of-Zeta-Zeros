#!/usr/bin/env python3
"""Build a formal-tail obstruction scout for the relative-Gaussian collar."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
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

from jensen_window_pf_heat_flow_monotone_closure_scout import REPO_ROOT  # noqa: E402
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    arb_rising,
    build_ratio_rows,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"

VALUE_BUDGET_A = "5.382819486765314521E-01"
DERIVATIVE_BUDGET_B = "9.315354075509573936E-03"


@dataclass(frozen=True)
class FormalTailProfileRow:
    index: int
    term_j_range: str
    value_term_j21: str
    value_term_j40: str
    value_term_j80: str
    value_term_j103: str
    value_term_j104: str
    value_term_j120: str
    derivative_term_j21: str
    derivative_term_j40: str
    derivative_term_j80: str
    derivative_term_j103: str
    derivative_term_j104: str
    derivative_term_j120: str
    value_sum_j21_to_j40: str
    value_sum_j41_to_j80: str
    value_sum_j81_to_j120: str
    derivative_sum_j21_to_j40: str
    derivative_sum_j41_to_j80: str
    derivative_sum_j81_to_j120: str
    value_budget_fraction_j21_to_j40: str
    value_budget_fraction_j41_to_j80: str
    value_budget_fraction_j81_to_j120: str
    derivative_budget_fraction_j21_to_j40: str
    derivative_budget_fraction_j41_to_j80: str
    derivative_budget_fraction_j81_to_j120: str
    least_value_term_j: int
    least_value_term: str
    least_derivative_term_j: int
    least_derivative_term: str
    first_value_growth_j_from_21: int
    first_value_growth_ratio_from_21: str
    first_value_growth_j_from_80: int
    first_value_growth_ratio_from_80: str
    first_derivative_growth_j_from_21: int
    first_derivative_growth_ratio_from_21: str
    first_derivative_growth_j_from_80: int
    first_derivative_growth_ratio_from_80: str
    formal_monotone_tail_from_j80_rejected: bool
    proof_boundary: str


def arb_text(value: flint.arb, digits: int = 18) -> str:
    return value.str(digits)


def sci(value: float) -> str:
    return f"{value:.18E}"


def first_growth(seq: list[tuple[int, flint.arb]], start_j: int) -> tuple[int, flint.arb]:
    for (j, left), (_next_j, right) in zip(seq, seq[1:]):
        if j >= start_j and right > left:
            return j, right / left
    raise RuntimeError(f"no growth found from j>={start_j}")


def sum_range(seq: list[tuple[int, flint.arb]], start_j: int, end_j: int) -> flint.arb:
    total = flint.arb(0)
    for j, value in seq:
        if start_j <= j <= end_j:
            total += value
    return total


def build_profile_rows(
    ratios: list[flint.arb],
    k: int,
    collar_start_T: int,
    first_j: int,
    last_j: int,
) -> list[FormalTailProfileRow]:
    u = flint.arb(1) / flint.arb(collar_start_T)
    value_budget = flint.arb(VALUE_BUDGET_A)
    derivative_budget = flint.arb(DERIVATIVE_BUDGET_B)
    rows: list[FormalTailProfileRow] = []

    for index in (k - 1, k, k + 1, k + 2):
        value_terms: list[tuple[int, flint.arb]] = []
        derivative_terms: list[tuple[int, flint.arb]] = []
        for j in range(first_j, last_j + 1):
            coeff = ratios[j] * arb_rising(index, j)
            value_terms.append((j, (coeff * (u**j)).abs_upper() / (u**3)))
            derivative_terms.append((j, (flint.arb(j) * coeff * (u ** (j - 1))).abs_upper() / u))

        values = dict(value_terms)
        derivatives = dict(derivative_terms)
        least_value_j, least_value = min(value_terms, key=lambda item: item[1])
        least_derivative_j, least_derivative = min(derivative_terms, key=lambda item: item[1])
        value_growth_21_j, value_growth_21 = first_growth(value_terms, first_j)
        value_growth_80_j, value_growth_80 = first_growth(value_terms, 80)
        derivative_growth_21_j, derivative_growth_21 = first_growth(derivative_terms, first_j)
        derivative_growth_80_j, derivative_growth_80 = first_growth(derivative_terms, 80)

        value_sum_21_40 = sum_range(value_terms, 21, 40)
        value_sum_41_80 = sum_range(value_terms, 41, 80)
        value_sum_81_120 = sum_range(value_terms, 81, 120)
        derivative_sum_21_40 = sum_range(derivative_terms, 21, 40)
        derivative_sum_41_80 = sum_range(derivative_terms, 41, 80)
        derivative_sum_81_120 = sum_range(derivative_terms, 81, 120)

        rows.append(
            FormalTailProfileRow(
                index=index,
                term_j_range=f"{first_j}..{last_j}",
                value_term_j21=arb_text(values[21]),
                value_term_j40=arb_text(values[40]),
                value_term_j80=arb_text(values[80]),
                value_term_j103=arb_text(values[103]),
                value_term_j104=arb_text(values[104]),
                value_term_j120=arb_text(values[120]),
                derivative_term_j21=arb_text(derivatives[21]),
                derivative_term_j40=arb_text(derivatives[40]),
                derivative_term_j80=arb_text(derivatives[80]),
                derivative_term_j103=arb_text(derivatives[103]),
                derivative_term_j104=arb_text(derivatives[104]),
                derivative_term_j120=arb_text(derivatives[120]),
                value_sum_j21_to_j40=arb_text(value_sum_21_40),
                value_sum_j41_to_j80=arb_text(value_sum_41_80),
                value_sum_j81_to_j120=arb_text(value_sum_81_120),
                derivative_sum_j21_to_j40=arb_text(derivative_sum_21_40),
                derivative_sum_j41_to_j80=arb_text(derivative_sum_41_80),
                derivative_sum_j81_to_j120=arb_text(derivative_sum_81_120),
                value_budget_fraction_j21_to_j40=arb_text(value_sum_21_40 / value_budget),
                value_budget_fraction_j41_to_j80=arb_text(value_sum_41_80 / value_budget),
                value_budget_fraction_j81_to_j120=arb_text(value_sum_81_120 / value_budget),
                derivative_budget_fraction_j21_to_j40=arb_text(derivative_sum_21_40 / derivative_budget),
                derivative_budget_fraction_j41_to_j80=arb_text(derivative_sum_41_80 / derivative_budget),
                derivative_budget_fraction_j81_to_j120=arb_text(derivative_sum_81_120 / derivative_budget),
                least_value_term_j=least_value_j,
                least_value_term=arb_text(least_value),
                least_derivative_term_j=least_derivative_j,
                least_derivative_term=arb_text(least_derivative),
                first_value_growth_j_from_21=value_growth_21_j,
                first_value_growth_ratio_from_21=arb_text(value_growth_21),
                first_value_growth_j_from_80=value_growth_80_j,
                first_value_growth_ratio_from_80=arb_text(value_growth_80),
                first_derivative_growth_j_from_21=derivative_growth_21_j,
                first_derivative_growth_ratio_from_21=arb_text(derivative_growth_21),
                first_derivative_growth_j_from_80=derivative_growth_80_j,
                first_derivative_growth_ratio_from_80=arb_text(derivative_growth_80),
                formal_monotone_tail_from_j80_rejected=True,
                proof_boundary=(
                    "Finite formal-term diagnostic only; it rejects a monotone/geometric formal-tail shortcut "
                    "but does not bound the actual asymptotic remainder."
                ),
            )
        )
    return rows


def build_diagnostics(max_degree: int, cutoff_n: int, precision_bits: int, k: int, collar_start_T: int) -> dict:
    if max_degree < 240 or max_degree % 2:
        raise ValueError("max_degree must be an even integer at least 240")
    flint.ctx.prec = precision_bits
    ratio_rows, ratios = build_ratio_rows(max_degree, cutoff_n)
    first_j = 21
    last_j = max_degree // 2
    profile_rows = build_profile_rows(ratios, k, collar_start_T, first_j, last_j)
    value_growth_ratios = [flint.arb(row.first_value_growth_ratio_from_80) for row in profile_rows]
    derivative_growth_ratios = [flint.arb(row.first_derivative_growth_ratio_from_80) for row in profile_rows]
    max_value_growth = max(value_growth_ratios)
    max_derivative_growth = max(derivative_growth_ratios)
    return {
        "parameters": {
            "max_taylor_degree": max_degree,
            "last_tested_j": last_j,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "collar_start_T": collar_start_T,
            "real_interval_u": f"[0, 1/{collar_start_T}]",
            "profile_indices": [k - 1, k, k + 1, k + 2],
            "formal_residual_first_j": first_j,
            "value_budget_A": VALUE_BUDGET_A,
            "derivative_budget_B": DERIVATIVE_BUDGET_B,
        },
        "ratio_ball_rows": len(ratio_rows),
        "highest_ratio_degree": int(ratio_rows[-1].degree),
        "highest_ratio_to_c0": ratio_rows[-1].ratio_to_c0,
        "highest_ratio_sign": ratio_rows[-1].sign,
        "profile_rows": [asdict(row) for row in profile_rows],
        "profile_row_count": len(profile_rows),
        "formal_tail_turnaround_rows": sum(1 for row in profile_rows if row.formal_monotone_tail_from_j80_rejected),
        "all_first_value_growth_from_j80_at_j103": all(row.first_value_growth_j_from_80 == 103 for row in profile_rows),
        "all_first_derivative_growth_from_j80_at_j103": all(
            row.first_derivative_growth_j_from_80 == 103 for row in profile_rows
        ),
        "max_first_value_growth_ratio_from_j80": arb_text(max_value_growth),
        "max_first_derivative_growth_ratio_from_j80": arb_text(max_derivative_growth),
        "fixed_radius_cauchy_obstruction": (
            "A fixed-radius Cauchy coefficient estimate |c_(2j)|<=C*R^(-2j) gives moment terms bounded by "
            "C*(a)_j*(u/R^2)^j. For any u/R^2>0, the consecutive ratio (a+j)*(u/R^2) tends to infinity, "
            "so this termwise majorant cannot be summed as an infinite absolute tail."
        ),
        "proof_boundary_note": (
            "This artifact rejects naive formal-tail summation templates. It does not reject an actual "
            "contour, saddle, integral-remainder, or optimally truncated asymptotic remainder theorem."
        ),
    }


def build_artifact(max_degree: int, cutoff_n: int, precision_bits: int, k: int, collar_start_T: int) -> dict:
    diagnostics = build_diagnostics(max_degree, cutoff_n, precision_bits, k, collar_start_T)
    rows = [
        {
            "id": "nlrgftos_01_formal_residual_terms",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The formal post-degree-40 terms for F_i use |(c_(2j)/c_0)*(i+1/2)_j*u^j|, with derivative terms multiplied by j/u.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
            ],
            "proof_boundary": "Exact formal-term bookkeeping only; it is not an actual remainder theorem.",
        },
        {
            "id": "nlrgftos_02_degree240_coefficient_profile",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "Arb coefficient-ratio balls through degree 240 extend the formal residual-term profile from j=21 through j=120.",
            "diagnostics": diagnostics,
            "proof_boundary": "Finite coefficient profile only; not a bound for all later formal terms or the actual remainder.",
        },
        {
            "id": "nlrgftos_03_finite_window_decay_until_least_term",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "In the tested window, the scaled formal terms keep shrinking to a least term near j=103 before rising.",
            "proof_boundary": "Finite least-term diagnostic only; not an infinite-tail estimate.",
        },
        {
            "id": "nlrgftos_04_turnaround_after_j103",
            "role": "finite_obstruction",
            "readiness": "not_ready_to_apply",
            "claim": "For each i=21..24, both value and derivative formal terms first grow again from j=103 to j=104 after j>=80.",
            "proof_boundary": "Finite obstruction to monotone formal-tail summation only; not a theorem about the exact heat-flow remainder.",
        },
        {
            "id": "nlrgftos_05_monotone_geometric_tail_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "A monotone decreasing or fixed-ratio geometric majorant beginning at degree 80 proves the degree-40 residual budget.",
            "gap": "The tested formal terms grow from j=103 to j=104, with growth ratios above 4 in all four index rows.",
            "proof_boundary": "Rejected shortcut only; it does not rule out a different non-termwise remainder proof.",
        },
        {
            "id": "nlrgftos_06_fixed_radius_cauchy_termwise_rejected",
            "role": "exact_obstruction",
            "readiness": "available_exact",
            "claim": "A fixed-radius Cauchy coefficient bound combined term-by-term with the moment rising factors cannot give a convergent infinite absolute tail.",
            "formula": "((a)_(j+1)*x^(j+1))/((a)_j*x^j) = (a+j)*x -> infinity for x>0",
            "proof_boundary": "Exact method obstruction only; not an obstruction to contour, saddle, or integral-remainder proofs.",
        },
        {
            "id": "nlrgftos_07_required_remainder_theorem",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "The residual-budget branch now requires an actual asymptotic remainder theorem, likely using contour/saddle/integral remainder control or an optimally truncated least-term estimate.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md",
            ],
            "proof_boundary": "Live theorem-search route only; the required remainder theorem is not proved.",
        },
        {
            "id": "nlrgftos_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any promoted residual proof must distinguish the divergent formal asymptotic series from the actual analytic remainder and must compare interval-safe bounds to the degree-40 budgets.",
            "source_artifacts": [
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    profile_rows = diagnostics["profile_rows"]
    summary = {
        "matrix_rows": len(rows),
        "profile_rows": diagnostics["profile_row_count"],
        "formal_tail_turnaround_rows": diagnostics["formal_tail_turnaround_rows"],
        "last_tested_j": diagnostics["parameters"]["last_tested_j"],
        "max_taylor_degree": max_degree,
        "max_first_value_growth_ratio_from_j80": diagnostics["max_first_value_growth_ratio_from_j80"],
        "max_first_derivative_growth_ratio_from_j80": diagnostics["max_first_derivative_growth_ratio_from_j80"],
        "least_value_term_js": {f"F_{row['index']}": row["least_value_term_j"] for row in profile_rows},
        "least_derivative_term_js": {f"F_{row['index']}": row["least_derivative_term_j"] for row in profile_rows},
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The degree-240 formal-term scout rejects the naive infinite formal-tail summation route. "
            "After the degree-40 residual budget, formal terms continue decreasing through a least-term "
            "region near j=103, but every index F_21..F_24 then grows from j=103 to j=104 after j>=80. "
            "A fixed-radius Cauchy coefficient bound is also structurally insufficient once multiplied by "
            "the moment rising factors. The live route is therefore an actual asymptotic remainder theorem, "
            "not an infinite absolute sum of formal Taylor terms."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout",
        "date": "2026-07-07",
        "status": "finite theorem-search obstruction",
        "source_degree40_residual_tail_budget": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md"
        ),
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py"
        ),
        "proof_boundary": (
            "Finite theorem-search obstruction only. It rejects naive formal-tail summation templates, but it "
            "does not bound the actual residual remainder, does not prove scaled-curvature monotonicity, does "
            "not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The formal asymptotic term profile is not promoted to an actual remainder theorem.",
            "A monotone or fixed-ratio geometric formal-tail proof from degree 80 is rejected.",
            "The degree-40 residual budget remains a target, not a proved estimate.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][1]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['profile_rows']} profile rows, "
        f"{summary['formal_tail_turnaround_rows']} formal-tail turnaround rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Formal-Tail Obstruction Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search obstruction. This is not a proof",
        "of an analytic residual estimate, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout`.",
        "",
        "Proof boundary: this artifact rejects naive formal-tail summation",
        "templates. It does not reject or prove an actual contour, saddle,",
        "integral-remainder, or optimally truncated asymptotic remainder theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Formal-Term Profile",
        "",
        "```text",
        f"max Taylor degree: {summary['max_taylor_degree']}",
        f"tested j range: 21..{summary['last_tested_j']}",
        f"u interval: {diagnostics['parameters']['real_interval_u']}",
        f"value budget A: {diagnostics['parameters']['value_budget_A']}",
        f"derivative budget B: {diagnostics['parameters']['derivative_budget_B']}",
        f"highest coefficient ratio: c{diagnostics['highest_ratio_degree']}/c0 = {diagnostics['highest_ratio_to_c0']}",
        "```",
        "",
        "Per-index profile:",
        "",
        "```text",
    ]
    for row in diagnostics["profile_rows"]:
        lines.append(
            f"F_{row['index']}: least value j={row['least_value_term_j']} {row['least_value_term']}, "
            f"least derivative j={row['least_derivative_term_j']} {row['least_derivative_term']}, "
            f"first value growth after j>=80: j={row['first_value_growth_j_from_80']} ratio={row['first_value_growth_ratio_from_80']}, "
            f"first derivative growth after j>=80: j={row['first_derivative_growth_j_from_80']} ratio={row['first_derivative_growth_ratio_from_80']}"
        )
        lines.append(
            f"  value sums j21..40={row['value_sum_j21_to_j40']} "
            f"(budget fraction {row['value_budget_fraction_j21_to_j40']}), "
            f"j41..80={row['value_sum_j41_to_j80']} "
            f"(budget fraction {row['value_budget_fraction_j41_to_j80']}), "
            f"j81..120={row['value_sum_j81_to_j120']} "
            f"(budget fraction {row['value_budget_fraction_j81_to_j120']})"
        )
        lines.append(
            f"  derivative sums j21..40={row['derivative_sum_j21_to_j40']} "
            f"(budget fraction {row['derivative_budget_fraction_j21_to_j40']}), "
            f"j41..80={row['derivative_sum_j41_to_j80']} "
            f"(budget fraction {row['derivative_budget_fraction_j41_to_j80']}), "
            f"j81..120={row['derivative_sum_j81_to_j120']} "
            f"(budget fraction {row['derivative_budget_fraction_j81_to_j120']})"
        )
    lines.extend(
        [
            "```",
            "",
            "## Rejected Shortcut",
            "",
            "A monotone decreasing or fixed-ratio geometric formal tail from",
            "degree 80 is rejected in the tested window, because all four",
            "index rows grow again from `j=103` to `j=104`.",
            "",
            "The fixed-radius Cauchy termwise route is also structurally",
            "insufficient:",
            "",
            "```text",
            diagnostics["fixed_radius_cauchy_obstruction"],
            "```",
            "",
            "This does not close the degree-40 residual budget. It tells us",
            "what the next proof must look like: a genuine asymptotic remainder",
            "estimate for the actual analytic object, not an infinite absolute",
            "sum of the formal Taylor terms.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree40_residual_tail_budget"],
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
    parser.add_argument("--max-degree", type=int, default=240)
    parser.add_argument("--cutoff-n", type=int, default=80)
    parser.add_argument("--precision-bits", type=int, default=256)
    parser.add_argument("--tail-start-k", type=int, default=22)
    parser.add_argument("--collar-start-T", type=int, default=1156)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(args.max_degree, args.cutoff_n, args.precision_bits, args.tail_start_k, args.collar_start_T)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian formal-tail obstruction scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
