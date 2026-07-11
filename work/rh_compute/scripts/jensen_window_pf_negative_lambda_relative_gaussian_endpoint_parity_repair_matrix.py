#!/usr/bin/env python3
"""Build an endpoint parity-repair matrix for the relative-Gaussian compact route."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_INTERPOLATION_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.json"
)
DEFAULT_PHI_TAIL_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.md"

DEFAULT_T = 10000
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_SERIES_ORDER = 15
DEFAULT_PRECISION_BITS = 8192
DEFAULT_FIRST_PANEL_Y_RIGHT = 1
DEFAULT_COMPACT_Y_RIGHT = 200


@dataclass(frozen=True)
class MatrixRow:
    id: str
    role: str
    readiness: str
    claim: str
    proof_boundary: str
    source_artifacts: list[str]
    diagnostics: dict | None = None
    gap: str | None = None


def arb_text(value: flint.arb, digits: int = 50) -> str:
    return value.str(digits).replace("e", "E")


def arb_mid_text(value: flint.arb, digits: int = 50) -> str:
    return value.str(digits, radius=False).replace("e", "E")


def arb_upper_text(value: flint.arb, digits: int = 50) -> str:
    return value.upper().str(digits, radius=False).replace("e", "E")


def abs_upper(value: flint.arb) -> flint.arb:
    lower_abs = abs(flint.arb(value.lower()))
    upper_abs = abs(flint.arb(value.upper()))
    return lower_abs if lower_abs > upper_abs else upper_abs


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def source_paths(interpolation_route_path: Path, phi_tail_path: Path) -> dict[str, str]:
    return {
        "interpolation_route_json": interpolation_route_path.relative_to(REPO_ROOT).as_posix(),
        "interpolation_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_interpolation_remainder_route_matrix.md"
        ),
        "arb_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_arb_chebyshev_interpolant_moment_scout.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
        ),
        "phi_tail_json": phi_tail_path.relative_to(REPO_ROOT).as_posix(),
        "phi_tail_note": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_phi_tail_grid_certificate.md",
        "compact_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md"
        ),
        "intervalization_target": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "formal_tail_obstruction": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_formal_tail_obstruction_scout.md"
        ),
    }


def exp_linear_coeffs(a: int, order: int) -> list[flint.arb]:
    coeffs: list[flint.arb] = []
    factorial = 1
    for degree in range(order + 1):
        if degree:
            factorial *= degree
        coeffs.append((flint.arb(a) ** degree) / flint.arb(factorial))
    return coeffs


def exp_series(coefficients: list[flint.arb], order: int) -> list[flint.arb]:
    result = [flint.arb(0)] * (order + 1)
    result[0] = coefficients[0].exp()
    for degree in range(1, order + 1):
        total = flint.arb(0)
        for k in range(1, degree + 1):
            total += flint.arb(k) * coefficients[k] * result[degree - k]
        result[degree] = total / flint.arb(degree)
    return result


def add_series(left: list[flint.arb], right: list[flint.arb], order: int) -> list[flint.arb]:
    return [left[index] + right[index] for index in range(order + 1)]


def scale_series(series: list[flint.arb], scalar: flint.arb) -> list[flint.arb]:
    return [scalar * item for item in series]


def multiply_series(left: list[flint.arb], right: list[flint.arb], order: int) -> list[flint.arb]:
    result = [flint.arb(0)] * (order + 1)
    for i in range(order + 1):
        for j in range(order + 1 - i):
            result[i + j] += left[i] * right[j]
    return result


def finite_phi_series(term_count: int, order: int) -> list[flint.arb]:
    pi = flint.arb.pi()
    exp9 = exp_linear_coeffs(9, order)
    exp5 = exp_linear_coeffs(5, order)
    total = [flint.arb(0)] * (order + 1)
    for n in range(1, term_count + 1):
        n2 = flint.arb(n * n)
        n4 = n2 * n2
        a = flint.arb(2) * pi * pi * n4
        b = flint.arb(3) * pi * n2
        q = pi * n2
        prefactor = add_series(scale_series(exp9, a), scale_series(exp5, -b), order)
        exponent = []
        factorial = 1
        for degree in range(order + 1):
            if degree:
                factorial *= degree
            exponent.append(-q * (flint.arb(4) ** degree) / flint.arb(factorial))
        damping = exp_series(exponent, order)
        total = add_series(total, multiply_series(prefactor, damping, order), order)
    return total


def partial_odd_sum(odd_rows: list[dict], x_max: flint.arb, c0_lower: flint.arb) -> flint.arb:
    total = flint.arb(0)
    for row in odd_rows:
        total += flint.arb(row["odd_coefficient_abs_upper"]) * (x_max ** int(row["degree"])) / c0_lower
    return total


def build_diagnostics(interpolation_route: dict, phi_tail: dict) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    coefficients = finite_phi_series(DEFAULT_PHI_TERM_COUNT, DEFAULT_SERIES_ORDER)
    c0 = coefficients[0]
    if c0.contains(0) or not bool(c0 > 0):
        raise RuntimeError("finite Phi(0) did not certify positive")
    c0_lower = c0.lower()
    odd_rows = []
    max_row = None
    for degree in range(1, DEFAULT_SERIES_ORDER + 1, 2):
        coefficient = coefficients[degree]
        upper = abs_upper(coefficient)
        normalized = upper / c0_lower
        row = {
            "degree": degree,
            "odd_coefficient_ball": arb_text(coefficient, 80),
            "odd_coefficient_abs_upper": arb_mid_text(upper, 50),
            "normalized_abs_upper": arb_mid_text(normalized, 50),
            "proof_boundary": "Finite-truncation low-order odd Taylor coefficient only; not an all-order parity theorem.",
        }
        odd_rows.append(row)
        if max_row is None or upper > flint.arb(max_row["odd_coefficient_abs_upper"]):
            max_row = row
    assert max_row is not None
    first_panel_x_max = flint.arb(DEFAULT_FIRST_PANEL_Y_RIGHT) / flint.arb(DEFAULT_T)
    first_panel_x_max = first_panel_x_max.sqrt()
    compact_x_max = flint.arb(DEFAULT_COMPACT_Y_RIGHT) / flint.arb(DEFAULT_T)
    compact_x_max = compact_x_max.sqrt()
    first_panel_partial = partial_odd_sum(odd_rows, first_panel_x_max, c0_lower)
    compact_partial = partial_odd_sum(odd_rows, compact_x_max, c0_lower)
    phi_tail_summary = phi_tail["summary"]
    return {
        "T": DEFAULT_T,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "series_order": DEFAULT_SERIES_ORDER,
        "precision_bits": DEFAULT_PRECISION_BITS,
        "compact_interval": interpolation_route["summary"]["compact_interval"],
        "heaviest_panel": interpolation_route["summary"]["heaviest_panel"],
        "finite_c0_ball": arb_text(c0, 80),
        "finite_c0_lower": c0_lower.str(50, radius=False).replace("e", "E"),
        "odd_taylor_rows": odd_rows,
        "odd_taylor_row_count": len(odd_rows),
        "first_odd_coefficient_abs_upper": odd_rows[0]["odd_coefficient_abs_upper"],
        "max_low_order_odd_degree": max_row["degree"],
        "max_low_order_odd_abs_upper": max_row["odd_coefficient_abs_upper"],
        "first_panel_x_max": arb_mid_text(first_panel_x_max, 50),
        "compact_x_max": arb_mid_text(compact_x_max, 50),
        "first_panel_low_order_normalized_odd_partial_sum_upper": arb_mid_text(first_panel_partial, 50),
        "compact_low_order_normalized_odd_partial_sum_upper": arb_mid_text(compact_partial, 50),
        "phi_tail_certificate_rows": phi_tail_summary["matrix_rows"],
        "phi_tail_grid_tail_sources": phi_tail_summary["tail_source_rows"],
        "endpoint_problem": (
            "A y-panel Bernstein proof for finite Phi(sqrt(y/T)) needs exact evenness in x at y=0; "
            "finite-truncation odd Taylor coefficients are tail-sized but not a theorem of exact evenness."
        ),
        "infinite_even_core_route": (
            "Use the theta functional equation to split the infinite Phi kernel into an exact even analytic "
            "core in x, hence an analytic core in y=x^2, and charge the finite truncation/tail parity "
            "defect to the existing Phi-tail source."
        ),
        "x_variable_endpoint_route": (
            "Alternatively handle the first panel in x with y=T*x^2 and Gamma density "
            "2*T^(alpha+1)*x^(2*alpha+1)*exp(-T*x^2)/Gamma(alpha+1), avoiding a y-plane branch theorem."
        ),
        "near_evenness_is_proof": False,
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgeprm_01_endpoint_branch_obligation",
            role="route_constraint",
            readiness="open_requirement",
            claim=(
                "The first compact panel touches y=0, so a y-plane interpolation-remainder theorem for "
                "Phi(sqrt(y/T)) must supply exact evenness/analyticity or avoid the branch by changing variables."
            ),
            diagnostics={"endpoint_problem": diagnostics["endpoint_problem"]},
            source_artifacts=[paths["interpolation_route_note"], paths["finite_part_note"]],
            proof_boundary="Route constraint only; no endpoint analytic-continuation theorem is proved.",
        ),
        MatrixRow(
            id="nlrgeprm_02_finite_truncation_odd_taylor_audit",
            role="arb_series_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "Arb power-series expansion of the finite n<=30 Phi truncation at x=0 finds tiny but "
                "nonzero low-order odd Taylor coefficients, consistent with parity defect living in the omitted tail."
            ),
            diagnostics={
                "phi_term_count": diagnostics["phi_term_count"],
                "series_order": diagnostics["series_order"],
                "finite_c0_ball": diagnostics["finite_c0_ball"],
                "odd_taylor_rows": diagnostics["odd_taylor_rows"],
                "first_odd_coefficient_abs_upper": diagnostics["first_odd_coefficient_abs_upper"],
                "max_low_order_odd_degree": diagnostics["max_low_order_odd_degree"],
                "max_low_order_odd_abs_upper": diagnostics["max_low_order_odd_abs_upper"],
            },
            source_artifacts=[paths["finite_part_note"], paths["phi_tail_note"]],
            proof_boundary="Low-order finite-truncation parity diagnostic only; not an all-order evenness theorem.",
        ),
        MatrixRow(
            id="nlrgeprm_03_low_order_branch_size_scout",
            role="arb_series_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "The normalized low-order odd Taylor contribution through degree 15 is far below numerical "
                "budgets on the first panel and on the whole compact interval, but this finite-order sum is not a proof."
            ),
            diagnostics={
                "first_panel_x_max": diagnostics["first_panel_x_max"],
                "compact_x_max": diagnostics["compact_x_max"],
                "first_panel_low_order_normalized_odd_partial_sum_upper": diagnostics[
                    "first_panel_low_order_normalized_odd_partial_sum_upper"
                ],
                "compact_low_order_normalized_odd_partial_sum_upper": diagnostics[
                    "compact_low_order_normalized_odd_partial_sum_upper"
                ],
            },
            source_artifacts=[paths["compact_scout_note"], paths["interpolation_route_note"]],
            proof_boundary="Finite-order size scout only; it does not bound all odd coefficients or the analytic remainder.",
        ),
        MatrixRow(
            id="nlrgeprm_04_infinite_even_core_tail_charge_route",
            role="live_route",
            readiness="not_ready_to_apply",
            claim=(
                "One admissible repair is to prove/use exact evenness of the infinite theta kernel, work with "
                "the resulting analytic y-core, and charge finite-truncation parity defect to the Phi-tail certificate."
            ),
            diagnostics={
                "infinite_even_core_route": diagnostics["infinite_even_core_route"],
                "phi_tail_certificate_rows": diagnostics["phi_tail_certificate_rows"],
                "phi_tail_grid_tail_sources": diagnostics["phi_tail_grid_tail_sources"],
            },
            source_artifacts=[paths["phi_tail_json"], paths["phi_tail_note"], paths["uniform_remainder_target"]],
            proof_boundary="Live repair route only; this artifact does not prove the exact even-kernel split or tail charge.",
        ),
        MatrixRow(
            id="nlrgeprm_05_x_variable_endpoint_route",
            role="live_route",
            readiness="not_ready_to_apply",
            claim=(
                "A second admissible repair is to certify the first panel directly in x, using y=T*x^2 and "
                "the transformed Gamma density, then hand off to y-panel Bernstein/Taylor models away from zero."
            ),
            diagnostics={"x_variable_endpoint_route": diagnostics["x_variable_endpoint_route"]},
            source_artifacts=[paths["interpolation_route_json"], paths["intervalization_target"]],
            proof_boundary="Live endpoint route only; no x-panel quadrature or Taylor-model certificate is produced.",
        ),
        MatrixRow(
            id="nlrgeprm_06_near_evenness_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="Near-zero low-order odd coefficients prove finite Phi(sqrt(y/T)) is analytic in y at the endpoint.",
            diagnostics={"near_evenness_is_proof": diagnostics["near_evenness_is_proof"]},
            gap=(
                "Near-evenness and low-order cancellation are calibration evidence. A referee-usable proof needs "
                "an exact even-kernel theorem with tail charges or an x-variable endpoint certificate."
            ),
            source_artifacts=[paths["interpolation_route_note"], paths["formal_tail_obstruction"]],
            proof_boundary="Rejected promotion only; near-evenness is not an endpoint analyticity proof.",
        ),
        MatrixRow(
            id="nlrgeprm_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted compact-interval interpolation theorem must choose one endpoint repair and validate it "
                "before applying y-plane Bernstein/Taylor remainders to the panel touching y=0."
            ),
            diagnostics={"target_closing": diagnostics["target_closing"], "ready_to_apply_rows": 0},
            source_artifacts=[paths["dependency_graph"], paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Proof-hygiene gate only; not a compact certificate, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(interpolation_route_path: Path, phi_tail_path: Path) -> dict:
    interpolation_route = load_json(interpolation_route_path)
    phi_tail = load_json(phi_tail_path)
    paths = source_paths(interpolation_route_path, phi_tail_path)
    diagnostics = build_diagnostics(interpolation_route, phi_tail)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "compact_interval": diagnostics["compact_interval"],
        "heaviest_panel": diagnostics["heaviest_panel"],
        "phi_term_count": diagnostics["phi_term_count"],
        "series_order": diagnostics["series_order"],
        "precision_bits": diagnostics["precision_bits"],
        "finite_c0_lower": diagnostics["finite_c0_lower"],
        "odd_taylor_row_count": diagnostics["odd_taylor_row_count"],
        "first_odd_coefficient_abs_upper": diagnostics["first_odd_coefficient_abs_upper"],
        "max_low_order_odd_degree": diagnostics["max_low_order_odd_degree"],
        "max_low_order_odd_abs_upper": diagnostics["max_low_order_odd_abs_upper"],
        "first_panel_x_max": diagnostics["first_panel_x_max"],
        "compact_x_max": diagnostics["compact_x_max"],
        "first_panel_low_order_normalized_odd_partial_sum_upper": diagnostics[
            "first_panel_low_order_normalized_odd_partial_sum_upper"
        ],
        "compact_low_order_normalized_odd_partial_sum_upper": diagnostics[
            "compact_low_order_normalized_odd_partial_sum_upper"
        ],
        "live_route_rows": 2,
        "rejected_route_rows": 1,
        "ready_to_apply_rows": 0,
        "target_closing": diagnostics["target_closing"],
        "near_evenness_is_proof": diagnostics["near_evenness_is_proof"],
        "main_finding": (
            "The endpoint parity-repair matrix makes the first-panel branch obligation explicit. "
            "For the finite n<=30 Phi truncation, Arb series coefficients through odd degree 15 are "
            "tail-sized rather than exactly zero: the first odd coefficient has absolute upper "
            f"{diagnostics['first_odd_coefficient_abs_upper']}, and the largest recorded low-order odd "
            f"coefficient is degree {diagnostics['max_low_order_odd_degree']} with upper "
            f"{diagnostics['max_low_order_odd_abs_upper']}. This supports, but does not prove, the "
            "correct repair: either use exact evenness of the infinite theta kernel plus a certified tail "
            "charge, or handle the endpoint panel in the x variable."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix",
        "date": "2026-07-08",
        "status": "endpoint parity-repair route matrix",
        "source_interpolation_remainder_route_matrix": paths["interpolation_route_note"],
        "source_interpolation_remainder_route_json": paths["interpolation_route_json"],
        "source_arb_chebyshev_interpolant_moment_scout": paths["arb_scout_note"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_phi_tail_grid_certificate": paths["phi_tail_note"],
        "source_phi_tail_grid_json": paths["phi_tail_json"],
        "source_compact_interval_integration_scout": paths["compact_scout_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "source_formal_tail_obstruction": paths["formal_tail_obstruction"],
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py",
        "proof_boundary": (
            "Endpoint parity-repair route matrix only. It audits low-order odd Taylor coefficients of the "
            "finite Phi truncation and records two admissible endpoint repair routes, but it does not prove "
            "exact evenness of the infinite theta kernel, does not certify the finite-truncation tail charge, "
            "does not produce an x-variable endpoint-panel certificate, does not prove an interpolation "
            "remainder theorem, does not prove a finite-grid interval certificate, does not prove a uniform "
            "collar theorem, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "Low-order near-evenness is not promoted to endpoint analyticity.",
            "The finite Phi truncation is not treated as exactly even.",
            "A compact y-plane interpolation theorem must repair the y=0 endpoint branch.",
            "The accepted repairs are exact infinite-even-core plus tail charge, or an x-variable endpoint-panel certificate.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['odd_taylor_row_count']} odd Taylor rows, "
        f"order {summary['series_order']}, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    odd_rows = artifact["matrix_rows"][1]["diagnostics"]["odd_taylor_rows"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Endpoint Parity-Repair Matrix",
        "",
        "Date: 2026-07-08",
        "",
        "Status: endpoint parity-repair route matrix. This is not a proof",
        "of endpoint analyticity, an interpolation-remainder theorem, a finite-grid",
        "interval certificate, a uniform collar theorem, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix`.",
        "",
        "Proof boundary: this artifact audits low-order odd Taylor coefficients",
        "of the finite Phi truncation and records admissible endpoint repair",
        "routes. It does not prove exact evenness or an endpoint certificate.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_endpoint_parity_repair_matrix.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Endpoint Audit",
        "",
        "```text",
        f"T: {summary['T']}",
        f"compact interval: {summary['compact_interval']}",
        f"finite Phi terms: {summary['phi_term_count']}",
        f"series order: {summary['series_order']}",
        f"finite c0 lower: {summary['finite_c0_lower']}",
        f"first odd coefficient abs upper: {summary['first_odd_coefficient_abs_upper']}",
        f"max low-order odd degree: {summary['max_low_order_odd_degree']}",
        f"max low-order odd abs upper: {summary['max_low_order_odd_abs_upper']}",
        f"first-panel x max: {summary['first_panel_x_max']}",
        f"compact x max: {summary['compact_x_max']}",
        f"first-panel normalized odd partial sum upper: {summary['first_panel_low_order_normalized_odd_partial_sum_upper']}",
        f"compact normalized odd partial sum upper: {summary['compact_low_order_normalized_odd_partial_sum_upper']}",
        "```",
        "",
        "Odd Taylor rows:",
        "",
        "```text",
    ]
    for row in odd_rows:
        lines.append(
            f"x^{row['degree']}: abs upper {row['odd_coefficient_abs_upper']}; "
            f"normalized {row['normalized_abs_upper']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Accepted repairs:",
            "",
            "```text",
            artifact["matrix_rows"][3]["diagnostics"]["infinite_even_core_route"],
            artifact["matrix_rows"][4]["diagnostics"]["x_variable_endpoint_route"],
            "```",
            "",
            "Rejected shortcut:",
            "",
            "```text",
            "Near-evenness and low-order cancellation are not an endpoint analyticity proof.",
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_interpolation_remainder_route_matrix"],
            artifact["source_interpolation_remainder_route_json"],
            artifact["source_arb_chebyshev_interpolant_moment_scout"],
            artifact["source_finite_part_weighted_sum_interval_certificate"],
            artifact["source_phi_tail_grid_certificate"],
            artifact["source_phi_tail_grid_json"],
            artifact["source_compact_interval_integration_scout"],
            artifact["source_intervalization_target"],
            artifact["source_uniform_remainder_target"],
            artifact["source_dependency_graph"],
            artifact["source_formal_tail_obstruction"],
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
    parser.add_argument("--interpolation-route-json", type=Path, default=DEFAULT_INTERPOLATION_ROUTE_JSON)
    parser.add_argument("--phi-tail-json", type=Path, default=DEFAULT_PHI_TAIL_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    interpolation_route_path = (
        args.interpolation_route_json
        if args.interpolation_route_json.is_absolute()
        else REPO_ROOT / args.interpolation_route_json
    )
    phi_tail_path = args.phi_tail_json if args.phi_tail_json.is_absolute() else REPO_ROOT / args.phi_tail_json
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(interpolation_route_path, phi_tail_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian endpoint parity-repair matrix: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
