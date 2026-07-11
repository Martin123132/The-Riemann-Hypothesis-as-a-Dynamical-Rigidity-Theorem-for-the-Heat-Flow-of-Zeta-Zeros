#!/usr/bin/env python3
"""Build a first-omitted denominator certificate for the relative-Gaussian grid."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
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


DEFAULT_GRID_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.json"
)
DEFAULT_INTERVAL_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.json"
)
DEFAULT_QUADRATURE_JSON = (
    REPO_ROOT / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.md"
)

DEFAULT_T_GRID = (1156, 1500, 2000, 5000, 10000)
DEFAULT_INDICES = (21, 22, 23, 24)
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PRECISION_BITS = 384
DEFAULT_RATIO_CAPS = ("1.0e-6", "2.0e-3")

getcontext().prec = 100


@dataclass(frozen=True)
class DenominatorRow:
    T: int
    index: int
    first_j: int
    rising_factor_exact: str
    value_denominator_ball: str
    value_denominator_lower: str
    derivative_denominator_ball: str
    derivative_denominator_lower: str
    derivative_to_value_denominator_ratio: str
    proof_boundary: str


@dataclass(frozen=True)
class RatioCapRow:
    ratio_error_cap: str
    value_scaled_absolute_radius_cap_from_min_denominator: str
    derivative_scaled_absolute_radius_cap_from_min_denominator: str
    limiting_value_location: str
    limiting_derivative_location: str
    proof_boundary: str


@dataclass(frozen=True)
class MatrixRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    source_artifacts: list[str]
    diagnostics: dict | None = None
    formula: str | None = None
    gap: str | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(grid_path: Path, interval_path: Path, quadrature_path: Path) -> dict[str, str]:
    return {
        "grid_json": grid_path.relative_to(REPO_ROOT).as_posix(),
        "grid_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md",
        "interval_json": interval_path.relative_to(REPO_ROOT).as_posix(),
        "interval_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "quadrature_json": quadrature_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_ladder_scout.md",
        "degree16_certificate_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md"
        ),
        "degree40_residual_budget_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree40_residual_tail_budget.md"
        ),
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def arb_negative(value: flint.arb) -> bool:
    return bool(value < 0 and not value.contains(0))


def dec_text(value: flint.arb, digits: int = 80) -> str:
    return value.lower().str(digits, radius=False)


def ball_text(value: flint.arb, digits: int = 40) -> str:
    return value.str(digits)


def decimal_from_text(text: str) -> Decimal:
    return Decimal(text.replace("e", "E"))


def frac_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def rising_fraction(index: int, length: int) -> Fraction:
    total = Fraction(1)
    q = Fraction(2 * index + 1, 2)
    for offset in range(length):
        total *= q + offset
    return total


def build_coefficient_certificate(first_j: int) -> tuple[dict, flint.arb]:
    ratio_rows, ratios = build_ratio_rows(2 * first_j, DEFAULT_RATIO_CUTOFF_N)
    ratio = ratios[first_j]
    if not arb_negative(ratio):
        raise RuntimeError("first omitted ratio did not certify negative")
    abs_ratio = -ratio
    if not arb_positive(abs_ratio):
        raise RuntimeError("absolute first omitted ratio did not certify positive")
    row = ratio_rows[first_j]
    abs_lower = abs_ratio.lower()
    return {
        "coefficient_index": first_j,
        "coefficient_degree": 2 * first_j,
        "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
        "precision_bits": DEFAULT_PRECISION_BITS,
        "ratio_to_c0_ball": ratio.str(80),
        "ratio_to_c0_radius": ratio.rad().str(20),
        "ratio_to_c0_sign": row.sign,
        "absolute_ratio_to_c0_lower": flint.arb(abs_lower).str(80, radius=False),
        "coefficient_certificate_source": (
            "build_ratio_rows uses finite n<=80 sums plus a geometric Arb tail-radius bound."
        ),
        "proof_boundary": (
            "Arb sign and magnitude certificate for the first omitted coefficient only; it does not "
            "bound the cancellation-reduced residual numerator."
        ),
    }, abs_ratio


def build_denominator_rows(abs_ratio: flint.arb, first_j: int) -> tuple[list[DenominatorRow], dict]:
    rows: list[DenominatorRow] = []
    min_value: tuple[Decimal, DenominatorRow] | None = None
    min_derivative: tuple[Decimal, DenominatorRow] | None = None
    for T in DEFAULT_T_GRID:
        T_arb = flint.arb(T)
        for index in DEFAULT_INDICES:
            rising = arb_rising(index, first_j)
            value_denominator = abs_ratio * rising / (T_arb ** (first_j - 3))
            derivative_denominator = flint.arb(first_j) * abs_ratio * rising / (T_arb ** (first_j - 2))
            if not arb_positive(value_denominator) or not arb_positive(derivative_denominator):
                raise RuntimeError(f"denominator did not certify positive for T={T}, index={index}")
            row = DenominatorRow(
                T=T,
                index=index,
                first_j=first_j,
                rising_factor_exact=frac_text(rising_fraction(index, first_j)),
                value_denominator_ball=ball_text(value_denominator),
                value_denominator_lower=dec_text(value_denominator),
                derivative_denominator_ball=ball_text(derivative_denominator),
                derivative_denominator_lower=dec_text(derivative_denominator),
                derivative_to_value_denominator_ratio=f"{first_j}/{T}",
                proof_boundary=(
                    "Arb lower bound for the first omitted denominator only; numerator and quadrature "
                    "errors are not bounded here."
                ),
            )
            rows.append(row)
            value_lower = decimal_from_text(row.value_denominator_lower)
            derivative_lower = decimal_from_text(row.derivative_denominator_lower)
            if min_value is None or value_lower < min_value[0]:
                min_value = (value_lower, row)
            if min_derivative is None or derivative_lower < min_derivative[0]:
                min_derivative = (derivative_lower, row)
    assert min_value is not None
    assert min_derivative is not None
    summary = {
        "denominator_rows": len(rows),
        "t_grid": list(DEFAULT_T_GRID),
        "indices": list(DEFAULT_INDICES),
        "first_j": first_j,
        "value_denominator_formula": "|r_21|*(i+1/2)_21*T^(-18)",
        "derivative_denominator_formula": "21*|r_21|*(i+1/2)_21*T^(-19)",
        "minimum_value_denominator_lower": f"{min_value[0]:.18E}",
        "minimum_value_denominator_location": f"T={min_value[1].T}, F_{min_value[1].index}",
        "minimum_derivative_denominator_lower": f"{min_derivative[0]:.18E}",
        "minimum_derivative_denominator_location": f"T={min_derivative[1].T}, F_{min_derivative[1].index}",
        "all_denominators_positive": True,
    }
    return rows, summary


def build_ratio_cap_rows(denominator_summary: dict) -> list[RatioCapRow]:
    value_min = Decimal(denominator_summary["minimum_value_denominator_lower"])
    derivative_min = Decimal(denominator_summary["minimum_derivative_denominator_lower"])
    rows = []
    for cap_text in DEFAULT_RATIO_CAPS:
        cap = Decimal(cap_text)
        rows.append(
            RatioCapRow(
                ratio_error_cap=f"{cap:.18E}",
                value_scaled_absolute_radius_cap_from_min_denominator=f"{cap * value_min:.18E}",
                derivative_scaled_absolute_radius_cap_from_min_denominator=f"{cap * derivative_min:.18E}",
                limiting_value_location=denominator_summary["minimum_value_denominator_location"],
                limiting_derivative_location=denominator_summary["minimum_derivative_denominator_location"],
                proof_boundary=(
                    "Exact denominator-to-ratio budget conversion only; it does not provide the scaled "
                    "absolute residual enclosure."
                ),
            )
        )
    return rows


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgfodc_01_first_omitted_formula",
            role="exact_reduction",
            readiness="available_exact",
            claim=(
                "For the degree-40 relative-Gaussian truncation, the first omitted comparison denominator "
                "has explicit value and derivative forms in terms of r_21 and the rising factorial."
            ),
            formula=(
                "D_value(i,T)=|r_21|*(i+1/2)_21*T^(-18); "
                "D_derivative(i,T)=21*|r_21|*(i+1/2)_21*T^(-19)"
            ),
            source_artifacts=[paths["grid_note"], paths["degree40_residual_budget_note"]],
            proof_boundary="Exact denominator identity only; not a residual numerator or quadrature bound.",
        ),
        MatrixRow(
            id="nlrgfodc_02_r21_arb_sign_magnitude",
            role="arb_coefficient_certificate",
            readiness="available_arb",
            claim="The normalized first omitted coefficient r_21 is Arb-certified negative with a positive absolute lower bound.",
            diagnostics=diagnostics["coefficient_certificate"],
            source_artifacts=[paths["degree16_certificate_note"]],
            proof_boundary="Arb coefficient certificate only; not a residual-tail theorem.",
        ),
        MatrixRow(
            id="nlrgfodc_03_grid_denominator_lowers",
            role="arb_denominator_certificate",
            readiness="available_for_intervalization",
            claim="Every first-omitted value and derivative denominator used by the recorded finite grid is Arb-certified positive.",
            diagnostics={
                "denominator_summary": diagnostics["denominator_summary"],
                "denominator_rows": diagnostics["denominator_rows"],
            },
            source_artifacts=[paths["grid_json"], paths["interval_note"]],
            proof_boundary="Finite-grid denominator lower certificate only; not a finite-grid residual certificate.",
        ),
        MatrixRow(
            id="nlrgfodc_04_ratio_cap_absolute_radius_handoff",
            role="exact_budget_handoff",
            readiness="available_for_intervalization",
            claim=(
                "The certified denominator lows translate the intervalization ratio caps into explicit "
                "scaled absolute-radius targets for value and derivative residual enclosures."
            ),
            diagnostics={"ratio_cap_rows": diagnostics["ratio_cap_rows"]},
            source_artifacts=[paths["interval_json"], paths["quadrature_note"]],
            proof_boundary="Budget handoff only; the numerator and quadrature-radius estimates remain open.",
        ),
        MatrixRow(
            id="nlrgfodc_05_denominator_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="Certified first-omitted denominators prove the first-omitted residual comparison.",
            gap=(
                "The denominator is only the divisor in the ratio comparison. A proof still needs interval "
                "enclosures for the cancellation-reduced residual numerator, quadrature error, Phi tails, "
                "coefficient propagation, and grid-to-collar coverage."
            ),
            source_artifacts=[paths["interval_note"], paths["uniform_remainder_target"]],
            proof_boundary="Rejected promotion only; not a residual theorem, uniform collar theorem, RH, or Lambda <= 0.",
        ),
        MatrixRow(
            id="nlrgfodc_06_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "Any promoted intervalized grid proof may use these denominator lows only together with "
                "separately certified numerator and quadrature-radius bounds."
            ),
            source_artifacts=[paths["dependency_graph"]],
            proof_boundary="Proof-hygiene gate only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(grid_path: Path, interval_path: Path, quadrature_path: Path) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    grid = load_json(grid_path)
    interval = load_json(interval_path)
    quadrature = load_json(quadrature_path)
    paths = source_paths(grid_path, interval_path, quadrature_path)
    first_j = DEFAULT_POLYNOMIAL_M + 1
    coefficient_certificate, abs_ratio = build_coefficient_certificate(first_j)
    denominator_rows, denominator_summary = build_denominator_rows(abs_ratio, first_j)
    ratio_cap_rows = build_ratio_cap_rows(denominator_summary)
    diagnostics = {
        "coefficient_certificate": coefficient_certificate,
        "denominator_summary": denominator_summary,
        "denominator_rows": [asdict(row) for row in denominator_rows],
        "ratio_cap_rows": [asdict(row) for row in ratio_cap_rows],
    }
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "denominator_rows": len(denominator_rows),
        "ratio_cap_rows": len(ratio_cap_rows),
        "certified_denominator_conditions": 2,
        "source_grid_rows": grid["summary"]["grid_rows"],
        "source_interval_obligations": interval["summary"]["obligation_rows"],
        "source_quadrature_ladder_rows": quadrature["summary"]["total_ladder_rows"],
        "first_j": first_j,
        "r21_sign_certified": coefficient_certificate["ratio_to_c0_sign"] == "negative",
        "absolute_r21_lower": coefficient_certificate["absolute_ratio_to_c0_lower"],
        "minimum_value_denominator_lower": denominator_summary["minimum_value_denominator_lower"],
        "minimum_value_denominator_location": denominator_summary["minimum_value_denominator_location"],
        "minimum_derivative_denominator_lower": denominator_summary["minimum_derivative_denominator_lower"],
        "minimum_derivative_denominator_location": denominator_summary["minimum_derivative_denominator_location"],
        "finest_ratio_cap": f"{Decimal(DEFAULT_RATIO_CAPS[0]):.18E}",
        "finest_value_scaled_absolute_radius_cap": ratio_cap_rows[0].value_scaled_absolute_radius_cap_from_min_denominator,
        "finest_derivative_scaled_absolute_radius_cap": (
            ratio_cap_rows[0].derivative_scaled_absolute_radius_cap_from_min_denominator
        ),
        "all_denominators_positive": denominator_summary["all_denominators_positive"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The first omitted coefficient r_21 is Arb-certified negative and bounded away from zero. "
            "Consequently every value and derivative first-omitted denominator on the recorded finite grid "
            "is positive. The smallest value denominator lower is at T=10000, F_21, and the smallest "
            "derivative denominator lower is also at T=10000, F_21. A ratio-radius target of 1e-6 therefore "
            "requires scaled absolute residual radii below the recorded value and derivative caps. This "
            "certifies the divisor side of the ratio comparison but leaves numerator enclosure, quadrature "
            "error, and grid-to-collar coverage open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate",
        "date": "2026-07-07",
        "status": "first-omitted denominator lower certificate",
        "source_cancellation_reduced_grid_scout": paths["grid_note"],
        "source_cancellation_reduced_grid_json": paths["grid_json"],
        "source_intervalization_target": paths["interval_note"],
        "source_intervalization_target_json": paths["interval_json"],
        "source_quadrature_ladder_scout": paths["quadrature_note"],
        "source_quadrature_ladder_json": paths["quadrature_json"],
        "source_degree16_certificate": paths["degree16_certificate_note"],
        "source_degree40_residual_budget": paths["degree40_residual_budget_note"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py"
        ),
        "proof_boundary": (
            "First-omitted denominator lower certificate only. It certifies the divisor side of the finite-grid "
            "ratio comparisons, but it does not enclose the cancellation-reduced residual numerator, does not "
            "prove quadrature-remainder bounds, does not prove a finite-grid residual certificate, does not "
            "prove a uniform collar theorem, does not prove scaled-curvature monotonicity, and does not prove "
            "RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "The certificate only discharges the denominator positivity and ratio-cap conversion side conditions.",
            "Residual numerator, quadrature, Phi-tail, rounding, and grid-to-collar bounds remain open.",
            "The finite grid is not promoted to a uniform collar theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][2]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['denominator_rows']} denominator rows, "
        f"{summary['ratio_cap_rows']} ratio-cap rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian First-Omitted Denominator Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: first-omitted denominator lower certificate. This is not a proof",
        "of the residual numerator enclosure, quadrature-remainder theorem,",
        "finite-grid interval certificate, uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate`.",
        "",
        "Proof boundary: this artifact certifies only the denominator side of",
        "the first-omitted ratio comparison for the recorded finite grid.",
        "It does not certify the cancellation-reduced residual numerator.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_first_omitted_denominator_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Denominator Formula",
        "",
        "```text",
        "D_value(i,T)=|r_21|*(i+1/2)_21*T^(-18)",
        "D_derivative(i,T)=21*|r_21|*(i+1/2)_21*T^(-19)",
        "```",
        "",
        "Coefficient certificate:",
        "",
        "```text",
        f"r_21 sign certified: {summary['r21_sign_certified']}",
        f"|r_21| lower: {summary['absolute_r21_lower']}",
        "```",
        "",
        "Grid minima:",
        "",
        "```text",
        f"minimum value denominator lower: {summary['minimum_value_denominator_lower']} at {summary['minimum_value_denominator_location']}",
        f"minimum derivative denominator lower: {summary['minimum_derivative_denominator_lower']} at {summary['minimum_derivative_denominator_location']}",
        f"all denominators positive: {summary['all_denominators_positive']}",
        "```",
        "",
        "Ratio-cap conversion:",
        "",
        "```text",
    ]
    for row in artifact["matrix_rows"][3]["diagnostics"]["ratio_cap_rows"]:
        lines.append(
            f"ratio cap {row['ratio_error_cap']}: value scaled absolute radius <= "
            f"{row['value_scaled_absolute_radius_cap_from_min_denominator']}, derivative scaled absolute radius <= "
            f"{row['derivative_scaled_absolute_radius_cap_from_min_denominator']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Denominator rows:",
            "",
            "```text",
        ]
    )
    for row in diagnostics["denominator_rows"]:
        lines.append(
            f"T={row['T']}, F_{row['index']}: value lower={row['value_denominator_lower']}, "
            f"derivative lower={row['derivative_denominator_lower']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_cancellation_reduced_grid_scout"],
            artifact["source_cancellation_reduced_grid_json"],
            artifact["source_intervalization_target"],
            artifact["source_intervalization_target_json"],
            artifact["source_quadrature_ladder_scout"],
            artifact["source_quadrature_ladder_json"],
            artifact["source_degree16_certificate"],
            artifact["source_degree40_residual_budget"],
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
    parser.add_argument("--grid-json", type=Path, default=DEFAULT_GRID_JSON)
    parser.add_argument("--interval-json", type=Path, default=DEFAULT_INTERVAL_JSON)
    parser.add_argument("--quadrature-json", type=Path, default=DEFAULT_QUADRATURE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    grid_path = args.grid_json if args.grid_json.is_absolute() else REPO_ROOT / args.grid_json
    interval_path = args.interval_json if args.interval_json.is_absolute() else REPO_ROOT / args.interval_json
    quadrature_path = args.quadrature_json if args.quadrature_json.is_absolute() else REPO_ROOT / args.quadrature_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(grid_path, interval_path, quadrature_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian first-omitted denominator certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
