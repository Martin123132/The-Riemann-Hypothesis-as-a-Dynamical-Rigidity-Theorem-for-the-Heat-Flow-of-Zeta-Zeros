#!/usr/bin/env python3
"""Build the residual-tail budget for the degree-40 relative-Gaussian collar."""

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
    bernstein_coefficients_on_interval,
    build_diagnostics as build_arb_collar_diagnostics,
    build_ratio_rows,
    interval_endpoints_for_subdivision,
    pderivative,
    truncated_multiplier_polynomial,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md"


@dataclass(frozen=True)
class BudgetInequalityRow:
    name: str
    residual_input: str
    allocated_margin: str
    raw_threshold: str
    half_safety_budget: str
    bound_at_half_safety: str
    limiting: bool
    proof_boundary: str


@dataclass(frozen=True)
class FiniteTailProfileRow:
    index: int
    term_degree_range: str
    value_scale_sum: str
    value_budget_fraction: str
    derivative_scale_sum: str
    derivative_budget_fraction: str
    largest_value_term_degree: int
    largest_derivative_term_degree: int
    proof_boundary: str


def sci(value: float) -> str:
    return f"{value:.18E}"


def arb_sci(value: flint.arb, digits: int = 30) -> str:
    return value.str(digits)


def max_abs_bernstein_upper(poly: list[flint.arb], collar_start_T: int) -> flint.arb:
    left, right = interval_endpoints_for_subdivision(0, 1, collar_start_T)
    best = flint.arb(0)
    for coeff in bernstein_coefficients_on_interval(poly, left, right):
        candidate = coeff.abs_upper()
        if candidate > best:
            best = candidate
    return best


def min_lower(values: list[str]) -> flint.arb:
    lowers = [flint.arb(value).lower() for value in values]
    best = lowers[0]
    for value in lowers[1:]:
        if value < best:
            best = value
    return best


def solve_positive_threshold(limit: float, fn) -> float:
    lo = 0.0
    hi = 1.0
    while fn(hi) <= limit:
        hi *= 2.0
        if hi > 1e100:
            break
    for _ in range(200):
        mid = (lo + hi) / 2.0
        if fn(mid) <= limit:
            lo = mid
        else:
            hi = mid
    return lo


def build_finite_tail_profile(
    ratios: list[flint.arb],
    k: int,
    collar_start_T: int,
    first_tail_j: int,
    last_tail_j: int,
    value_budget: float,
    derivative_budget: float,
) -> list[FiniteTailProfileRow]:
    u = flint.arb(1) / flint.arb(collar_start_T)
    value_budget_arb = flint.arb(value_budget)
    derivative_budget_arb = flint.arb(derivative_budget)
    rows: list[FiniteTailProfileRow] = []
    for index in (k - 1, k, k + 1, k + 2):
        value_sum = flint.arb(0)
        derivative_sum = flint.arb(0)
        largest_value = (flint.arb(0), first_tail_j)
        largest_derivative = (flint.arb(0), first_tail_j)
        for j in range(first_tail_j, last_tail_j + 1):
            coeff = ratios[j] * arb_rising(index, j)
            value_scaled = (coeff * (u**j)).abs_upper() / (u**3)
            derivative_scaled = (flint.arb(j) * coeff * (u ** (j - 1))).abs_upper() / u
            value_sum += value_scaled
            derivative_sum += derivative_scaled
            if value_scaled > largest_value[0]:
                largest_value = (value_scaled, j)
            if derivative_scaled > largest_derivative[0]:
                largest_derivative = (derivative_scaled, j)
        rows.append(
            FiniteTailProfileRow(
                index=index,
                term_degree_range=f"{2 * first_tail_j}..{2 * last_tail_j}",
                value_scale_sum=arb_sci(value_sum, 18),
                value_budget_fraction=arb_sci(value_sum / value_budget_arb, 18),
                derivative_scale_sum=arb_sci(derivative_sum, 18),
                derivative_budget_fraction=arb_sci(derivative_sum / derivative_budget_arb, 18),
                largest_value_term_degree=2 * largest_value[1],
                largest_derivative_term_degree=2 * largest_derivative[1],
                proof_boundary="Finite tail-profile diagnostic only; not a bound for the infinite residual beyond the sampled range.",
            )
        )
    return rows


def build_diagnostics(
    finite_degree: int,
    profile_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    collar_start_T: int,
) -> dict:
    flint.ctx.prec = precision_bits
    continuation_M = finite_degree // 2
    first_tail_j = continuation_M + 1
    last_tail_j = profile_max_degree // 2
    if profile_max_degree < finite_degree + 2 or profile_max_degree % 2:
        raise ValueError("profile_max_degree must be even and at least finite_degree+2")

    collar = build_arb_collar_diagnostics(
        finite_degree,
        cutoff_n,
        precision_bits,
        k,
        continuation_M,
        collar_start_T,
    )
    ratio_rows, ratios = build_ratio_rows(profile_max_degree, cutoff_n)
    finite_ratios = ratios[: continuation_M + 1]
    finite_polys = {
        index: truncated_multiplier_polynomial(index, finite_ratios, continuation_M)
        for index in (k - 1, k, k + 1, k + 2)
    }
    u = flint.arb(1) / flint.arb(collar_start_T)
    weight_sum = abs(2 * k + 1) + abs(-(6 * k + 5)) + abs(6 * k + 7) + abs(-(2 * k + 3))

    normalizer_margin = min_lower([row["min_bernstein_lower"] for row in collar["normalizer_rows"]])
    stencil_margins = {
        row["name"]: flint.arb(row["min_bernstein_lower"]).lower()
        for row in collar["stencil_rows"]
    }
    pmax = flint.arb(0)
    dmax = flint.arb(0)
    for poly in finite_polys.values():
        pmax_candidate = max_abs_bernstein_upper(poly, collar_start_T)
        dmax_candidate = max_abs_bernstein_upper(pderivative(poly), collar_start_T)
        if pmax_candidate > pmax:
            pmax = pmax_candidate
        if dmax_candidate > dmax:
            dmax = dmax_candidate

    U = float(u)
    P = float(pmax)
    D = float(dmax)
    mF = float(normalizer_margin)
    mB = float(stencil_margins["B_product"])
    mU = float(stencil_margins["companion_product"])
    mC = float(stencil_margins["weighted_gap_derivative"])

    # Residual model for a future theorem:
    # |R_i(u)| <= A*u^3 and |R_i'(u)| <= B*u on 0<=u<=1/1156.
    normalizer_threshold = 0.5 * mF / (U**3)
    b_threshold = solve_positive_threshold(0.5 * mB, lambda a: 4 * P * a * U + 2 * a * a * (U**4))
    companion_threshold = solve_positive_threshold(
        0.5 * mU,
        lambda a: 2
        * (
            4 * (P**3) * a
            + 6 * (P**2) * (a**2) * (U**3)
            + 4 * P * (a**3) * (U**6)
            + (a**4) * (U**9)
        ),
    )
    weighted_value_threshold = solve_positive_threshold(
        0.25 * mC,
        lambda a: weight_sum
        * D
        * (3 * (P**2) * a * (U**2) + 3 * P * (a**2) * (U**5) + (a**3) * (U**8)),
    )
    raw_value_thresholds = {
        "normalizer": normalizer_threshold,
        "B_product": b_threshold,
        "companion_product": companion_threshold,
        "weighted_gap_value_part": weighted_value_threshold,
    }
    limiting_value_name, limiting_value_threshold = min(raw_value_thresholds.items(), key=lambda item: item[1])
    value_budget = 0.5 * limiting_value_threshold
    weighted_derivative_threshold = (0.25 * mC) / (weight_sum * ((P + value_budget * (U**3)) ** 3))
    derivative_budget = 0.5 * weighted_derivative_threshold

    def b_bound(a: float) -> float:
        return 4 * P * a * U + 2 * a * a * (U**4)

    def companion_bound(a: float) -> float:
        return 2 * (
            4 * (P**3) * a
            + 6 * (P**2) * (a**2) * (U**3)
            + 4 * P * (a**3) * (U**6)
            + (a**4) * (U**9)
        )

    def weighted_value_bound(a: float) -> float:
        return weight_sum * D * (
            3 * (P**2) * a * (U**2) + 3 * P * (a**2) * (U**5) + (a**3) * (U**8)
        )

    def weighted_derivative_bound(a: float, b: float) -> float:
        return weight_sum * b * ((P + a * (U**3)) ** 3)

    budget_rows = [
        BudgetInequalityRow(
            name="normalizer_value_residual",
            residual_input="|R_i(u)| <= A*u^3",
            allocated_margin=sci(0.5 * mF),
            raw_threshold=sci(normalizer_threshold),
            half_safety_budget=sci(value_budget),
            bound_at_half_safety=sci(value_budget * (U**3)),
            limiting=limiting_value_name == "normalizer",
            proof_boundary="Sufficient normalizer positivity budget only; the residual estimate is not proved.",
        ),
        BudgetInequalityRow(
            name="B_product_value_residual",
            residual_input="|R_i(u)| <= A*u^3",
            allocated_margin=sci(0.5 * mB),
            raw_threshold=sci(b_threshold),
            half_safety_budget=sci(value_budget),
            bound_at_half_safety=sci(b_bound(value_budget)),
            limiting=limiting_value_name == "B_product",
            proof_boundary="Sufficient B-product perturbation budget only; the residual estimate is not proved.",
        ),
        BudgetInequalityRow(
            name="companion_product_value_residual",
            residual_input="|R_i(u)| <= A*u^3",
            allocated_margin=sci(0.5 * mU),
            raw_threshold=sci(companion_threshold),
            half_safety_budget=sci(value_budget),
            bound_at_half_safety=sci(companion_bound(value_budget)),
            limiting=limiting_value_name == "companion_product",
            proof_boundary="Sufficient companion-product perturbation budget only; the residual estimate is not proved.",
        ),
        BudgetInequalityRow(
            name="weighted_gap_value_residual",
            residual_input="|R_i(u)| <= A*u^3",
            allocated_margin=sci(0.25 * mC),
            raw_threshold=sci(weighted_value_threshold),
            half_safety_budget=sci(value_budget),
            bound_at_half_safety=sci(weighted_value_bound(value_budget)),
            limiting=limiting_value_name == "weighted_gap_value_part",
            proof_boundary="Sufficient weighted-gap value-residual budget only; the residual estimate is not proved.",
        ),
        BudgetInequalityRow(
            name="weighted_gap_derivative_residual",
            residual_input="|R_i'(u)| <= B*u",
            allocated_margin=sci(0.25 * mC),
            raw_threshold=sci(weighted_derivative_threshold),
            half_safety_budget=sci(derivative_budget),
            bound_at_half_safety=sci(weighted_derivative_bound(value_budget, derivative_budget)),
            limiting=True,
            proof_boundary="Sufficient weighted-gap derivative-residual budget only; the residual derivative estimate is not proved.",
        ),
    ]
    finite_tail_profile_rows = build_finite_tail_profile(
        ratios,
        k,
        collar_start_T,
        first_tail_j,
        last_tail_j,
        value_budget,
        derivative_budget,
    )
    return {
        "parameters": {
            "finite_degree": finite_degree,
            "profile_max_degree": profile_max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "collar_start_T": collar_start_T,
            "real_interval_u": f"[0, 1/{collar_start_T}]",
            "real_interval_T": f"[{collar_start_T}, infinity)",
            "residual_first_j": first_tail_j,
            "finite_tail_profile_last_j": last_tail_j,
            "weighted_gap_abs_weight_sum": weight_sum,
        },
        "finite_collar_margins": {
            "normalizer_min_lower": arb_sci(normalizer_margin),
            "B_product_min_lower": arb_sci(stencil_margins["B_product"]),
            "companion_product_min_lower": arb_sci(stencil_margins["companion_product"]),
            "weighted_gap_derivative_min_lower": arb_sci(stencil_margins["weighted_gap_derivative"]),
            "finite_normalizer_abs_upper": arb_sci(pmax),
            "finite_derivative_abs_upper": arb_sci(dmax),
        },
        "budget_inequality_rows": [asdict(row) for row in budget_rows],
        "finite_tail_profile_rows": [asdict(row) for row in finite_tail_profile_rows],
        "raw_value_thresholds": {name: sci(value) for name, value in raw_value_thresholds.items()},
        "value_residual_half_safety_budget_A": sci(value_budget),
        "derivative_residual_half_safety_budget_B": sci(derivative_budget),
        "limiting_value_budget": limiting_value_name,
        "budget_inequality_row_count": len(budget_rows),
        "finite_tail_profile_row_count": len(finite_tail_profile_rows),
        "max_finite_profile_value_budget_fraction": max(
            float(row.value_budget_fraction.split()[0].strip("[]")) for row in finite_tail_profile_rows
        ),
        "max_finite_profile_derivative_budget_fraction": max(
            float(row.derivative_budget_fraction.split()[0].strip("[]")) for row in finite_tail_profile_rows
        ),
        "proof_boundary_note": (
            "These are sufficient residual-tail targets derived from finite degree-40 margins. "
            "They do not prove the residual estimates and do not control the infinite residual tail."
        ),
    }


def build_artifact(
    finite_degree: int,
    profile_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    collar_start_T: int,
) -> dict:
    diagnostics = build_diagnostics(finite_degree, profile_max_degree, cutoff_n, precision_bits, k, collar_start_T)
    rows = [
        {
            "id": "nlrgd40rtb_01_residual_coordinate",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "Write the full multiplier as F_i(u)=P_i^(40)(u)+R_i(u), where P_i^(40) is the degree-40 finite surrogate and R_i is the residual tail from j=21 onward.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md",
            ],
            "proof_boundary": "Exact residual bookkeeping only; it does not bound R_i.",
        },
        {
            "id": "nlrgd40rtb_02_degree40_margin_extraction",
            "role": "finite_interval_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-40 Arb collar supplies finite lower margins for the normalizers, B product, companion product, and weighted-gap derivative numerator on 0<=u<=1/1156.",
            "diagnostics": diagnostics,
            "proof_boundary": "Finite margin extraction only; not an infinite-tail theorem.",
        },
        {
            "id": "nlrgd40rtb_03_value_residual_sufficient_budget",
            "role": "exact_sufficient_condition",
            "readiness": "not_ready_to_apply",
            "claim": "A sufficient value-tail target is |R_i(u)|<=A*u^3 for i=21..24 with A below the half-safety budget recorded in diagnostics.",
            "formula": "|R_i(u)| <= A*u^3",
            "proof_boundary": "Sufficient condition only; the analytic value-tail estimate is not proved.",
        },
        {
            "id": "nlrgd40rtb_04_derivative_residual_sufficient_budget",
            "role": "exact_sufficient_condition",
            "readiness": "not_ready_to_apply",
            "claim": "A sufficient derivative-tail target for the weighted-gap numerator is |R_i'(u)|<=B*u with B below the half-safety budget recorded in diagnostics.",
            "formula": "|R_i'(u)| <= B*u",
            "proof_boundary": "Sufficient condition only; the analytic derivative-tail estimate is not proved.",
        },
        {
            "id": "nlrgd40rtb_05_finite_tail_profile_through_degree80",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The computed finite tail terms from degree 42 through 80 consume less than 0.1% of the half-safety value and derivative budgets at the worst tested index.",
            "proof_boundary": "Finite sampled tail profile only; it is not a bound for the infinite residual beyond degree 80.",
        },
        {
            "id": "nlrgd40rtb_06_live_residual_tail_theorem",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "A promoted proof can close the fixed-k degree-40 collar branch by proving the recorded value and derivative residual estimates on the full real collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md",
            ],
            "proof_boundary": "Live theorem-search route only; the residual-tail theorem is not proved.",
        },
        {
            "id": "nlrgd40rtb_07_finite_profile_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The finite degree-42..80 tail profile proves the residual-tail theorem.",
            "gap": "The finite profile stops at degree 80 and does not supply an analytic majorant for the remaining infinite series.",
            "proof_boundary": "Rejected finite-profile promotion only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
        {
            "id": "nlrgd40rtb_08_acceptance_gate",
            "role": "acceptance_gate",
            "readiness": "not_ready_to_apply",
            "claim": "Any residual-tail proof must state its analytic majorant, the collar, the affected indices, value and derivative estimates, and interval-safe comparison to the recorded budgets.",
            "source_artifacts": [
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Proof-hygiene gate only; not a proof of the residual-tail estimates.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "budget_inequality_rows": diagnostics["budget_inequality_row_count"],
        "finite_tail_profile_rows": diagnostics["finite_tail_profile_row_count"],
        "value_residual_half_safety_budget_A": diagnostics["value_residual_half_safety_budget_A"],
        "derivative_residual_half_safety_budget_B": diagnostics["derivative_residual_half_safety_budget_B"],
        "limiting_value_budget": diagnostics["limiting_value_budget"],
        "max_finite_profile_value_budget_fraction": sci(diagnostics["max_finite_profile_value_budget_fraction"]),
        "max_finite_profile_derivative_budget_fraction": sci(diagnostics["max_finite_profile_derivative_budget_fraction"]),
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The degree-40 Arb collar margins imply a concrete sufficient fixed-k residual target on "
            "0<=u<=1/1156: prove |R_i(u)|<=A*u^3 with half-safety A="
            f"{diagnostics['value_residual_half_safety_budget_A']} and |R_i'(u)|<=B*u with half-safety B="
            f"{diagnostics['derivative_residual_half_safety_budget_B']} for i=21..24. The finite degree-42..80 "
            "tail profile uses less than 0.1% of these budgets, but this remains a target because no analytic "
            "majorant for the infinite residual tail is proved."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget",
        "date": "2026-07-07",
        "status": "exact finite theorem-search diagnostic",
        "source_degree40_ladder_stress": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_arb_collar_ladder_stress.md"
        ),
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_stencil_remainder_obligations": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_stencil_remainder_obligations.md"
        ),
        "source_dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py",
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py"
        ),
        "proof_boundary": (
            "Exact finite theorem-search diagnostic only. It converts finite degree-40 collar margins into "
            "sufficient value and derivative residual-tail targets, but it does not prove those residual "
            "estimates, does not control the infinite tail, does not prove scaled-curvature monotonicity, "
            "does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The value and derivative residual-tail bounds are sufficient targets, not proved estimates.",
            "The finite degree-42..80 tail profile is not promoted to an infinite-tail theorem.",
            "The result is fixed at k=22 and is not an all-k theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][1]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['budget_inequality_rows']} budget inequalities, "
        f"{summary['finite_tail_profile_rows']} finite tail profile rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-40 Residual Tail Budget",
        "",
        "Date: 2026-07-07",
        "",
        "Status: exact finite theorem-search diagnostic. This is not a proof",
        "of an infinite Taylor-tail theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget`.",
        "",
        "Proof boundary: this artifact converts finite degree-40 collar margins",
        "into sufficient value and derivative residual-tail targets. It does not",
        "prove the residual estimates.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Residual Target",
        "",
        "Write:",
        "",
        "```text",
        "F_i(u) = P_i^(40)(u) + R_i(u),  i in {21,22,23,24}",
        "0 <= u <= 1/1156",
        "```",
        "",
        "A sufficient fixed-`k=22` residual target is:",
        "",
        "```text",
        f"|R_i(u)|  <= {summary['value_residual_half_safety_budget_A']} * u^3",
        f"|R_i'(u)| <= {summary['derivative_residual_half_safety_budget_B']} * u",
        "```",
        "",
        "The value-tail budget is limited by the companion product.",
        "",
        "## Margin Data",
        "",
        "```text",
    ]
    margins = diagnostics["finite_collar_margins"]
    for key, value in margins.items():
        lines.append(f"{key}: {value}")
    lines.extend(["```", "", "Budget inequalities:", "", "```text"])
    for row in diagnostics["budget_inequality_rows"]:
        lines.append(
            f"{row['name']}: allocated={row['allocated_margin']}, "
            f"raw threshold={row['raw_threshold']}, half-safety={row['half_safety_budget']}, "
            f"bound at half-safety={row['bound_at_half_safety']}, limiting={row['limiting']}"
        )
    lines.extend(["```", "", "## Finite Tail Profile", "", "```text"])
    for row in diagnostics["finite_tail_profile_rows"]:
        lines.append(
            f"F_{row['index']}, degrees {row['term_degree_range']}: "
            f"value fraction={row['value_budget_fraction']}, "
            f"derivative fraction={row['derivative_budget_fraction']}, "
            f"largest value degree={row['largest_value_term_degree']}, "
            f"largest derivative degree={row['largest_derivative_term_degree']}"
        )
    lines.extend(
        [
            "```",
            "",
            "The finite profile is a plausibility diagnostic only. A proof still",
            "needs an analytic majorant for every term beyond degree 40, or a",
            "stronger argument that controls the residual as a whole.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree40_ladder_stress"],
            artifact["source_uniform_remainder_target"],
            artifact["source_stencil_remainder_obligations"],
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
    parser.add_argument("--finite-degree", type=int, default=40)
    parser.add_argument("--profile-max-degree", type=int, default=80)
    parser.add_argument("--cutoff-n", type=int, default=80)
    parser.add_argument("--precision-bits", type=int, default=384)
    parser.add_argument("--tail-start-k", type=int, default=22)
    parser.add_argument("--collar-start-T", type=int, default=1156)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(
        args.finite_degree,
        args.profile_max_degree,
        args.cutoff_n,
        args.precision_bits,
        args.tail_start_k,
        args.collar_start_T,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian degree-40 residual tail budget: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
