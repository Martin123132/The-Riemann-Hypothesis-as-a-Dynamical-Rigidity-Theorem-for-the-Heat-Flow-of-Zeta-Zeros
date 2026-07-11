#!/usr/bin/env python3
"""Build a worst-row Chebyshev panel-moment integration scout."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
import mpmath as mp  # noqa: E402

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate import (  # noqa: E402
    build_ratio_rows,
)


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_COMPACT_SCOUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.json"
)
DEFAULT_QUADRATURE_ROUTE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.md"
)

DEFAULT_T = 10000
DEFAULT_INDEX = 21
DEFAULT_ALPHA = "20.5"
DEFAULT_MPMATH_DPS = 100
DEFAULT_FLINT_PRECISION_BITS = 384
DEFAULT_POLYNOMIAL_M = 20
DEFAULT_RATIO_CUTOFF_N = 80
DEFAULT_PHI_TERM_COUNT = 30
DEFAULT_PANELS = ((0, 1), (1, 5), (5, 20), (20, 50), (50, 100), (100, 200))
DEFAULT_DEGREES = (12, 16, 20, 24, 32)


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


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def mp_text(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(value, digits, min_fixed=-8, max_fixed=8).replace("e", "E")


def source_paths(compact_scout_path: Path, quadrature_route_path: Path) -> dict[str, str]:
    return {
        "compact_scout_json": compact_scout_path.relative_to(REPO_ROOT).as_posix(),
        "compact_scout_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_compact_interval_integration_scout.md"
        ),
        "quadrature_route_json": quadrature_route_path.relative_to(REPO_ROOT).as_posix(),
        "quadrature_route_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_quadrature_remainder_route_matrix.md"
        ),
        "far_tail_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_far_tail_split_certificate.md"
        ),
        "finite_part_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_finite_part_weighted_sum_interval_certificate.md"
        ),
        "cancellation_grid_note": (
            "outputs/jensen_window_pf_negative_lambda_relative_gaussian_cancellation_reduced_remainder_grid_scout.md"
        ),
        "intervalization_target": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_intervalization_target.md",
        "uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
    }


def configure_precision() -> None:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    flint.ctx.prec = DEFAULT_FLINT_PRECISION_BITS


def ratio_midpoints() -> list[mp.mpf]:
    _ratio_rows, ratios_arb = build_ratio_rows(2 * (DEFAULT_POLYNOMIAL_M + 1), DEFAULT_RATIO_CUTOFF_N)
    return [
        mp.mpf(ratio.str(DEFAULT_MPMATH_DPS + 20, radius=False).replace("e", "E"))
        for ratio in ratios_arb[: DEFAULT_POLYNOMIAL_M + 1]
    ]


def finite_phi_mp(x_value: mp.mpf) -> mp.mpf:
    x = mp.mpf(x_value)
    pi = mp.pi
    total = mp.mpf("0")
    for n in range(1, DEFAULT_PHI_TERM_COUNT + 1):
        n2 = mp.mpf(n * n)
        q = pi * n2
        total += (
            2 * pi * pi * n2 * n2 * mp.exp(9 * x) - 3 * pi * n2 * mp.exp(5 * x)
        ) * mp.exp(-q * mp.exp(4 * x))
    return total


def finite_phip_mp(x_value: mp.mpf) -> mp.mpf:
    x = mp.mpf(x_value)
    pi = mp.pi
    total = mp.mpf("0")
    for n in range(1, DEFAULT_PHI_TERM_COUNT + 1):
        n2 = mp.mpf(n * n)
        q = pi * n2
        exp9 = mp.exp(9 * x)
        exp5 = mp.exp(5 * x)
        exp4 = mp.exp(4 * x)
        prefactor = 2 * pi * pi * n2 * n2 * exp9 - 3 * pi * n2 * exp5
        prefactor_derivative = 18 * pi * pi * n2 * n2 * exp9 - 15 * pi * n2 * exp5
        total += (prefactor_derivative - prefactor * 4 * q * exp4) * mp.exp(-q * exp4)
    return total


def make_core_functions(ratios: list[mp.mpf]):
    t = mp.mpf(DEFAULT_T)
    c0 = finite_phi_mp(mp.mpf("0"))

    def value_core(y_value: mp.mpf) -> mp.mpf:
        y = mp.mpf(y_value)
        v = y / t
        x = mp.sqrt(v)
        polynomial = mp.fsum(ratios[j] * (v**j) for j in range(DEFAULT_POLYNOMIAL_M + 1))
        return finite_phi_mp(x) / c0 - polynomial

    def derivative_core(y_value: mp.mpf) -> mp.mpf:
        y = mp.mpf(y_value)
        v = y / t
        x = mp.sqrt(v)
        polynomial = mp.fsum(
            mp.mpf(j) * ratios[j] * (v**j) for j in range(1, DEFAULT_POLYNOMIAL_M + 1)
        )
        return x * finite_phip_mp(x) / (2 * c0) - polynomial

    return value_core, derivative_core


CHEB_POWER_CACHE: dict[int, list[mp.mpf]] = {}


def chebyshev_power_coefficients(degree: int) -> list[mp.mpf]:
    if degree in CHEB_POWER_CACHE:
        return CHEB_POWER_CACHE[degree]
    if degree == 0:
        coefficients = [mp.mpf(1)]
    elif degree == 1:
        coefficients = [mp.mpf(0), mp.mpf(1)]
    else:
        previous = [mp.mpf(1)]
        current = [mp.mpf(0), mp.mpf(1)]
        for _ in range(1, degree):
            twice_t_current = [mp.mpf(0)] + [2 * item for item in current]
            if len(previous) < len(twice_t_current):
                previous = previous + [mp.mpf(0)] * (len(twice_t_current) - len(previous))
            next_coefficients = [twice_t_current[i] - previous[i] for i in range(len(twice_t_current))]
            previous, current = current, next_coefficients
        coefficients = current
    CHEB_POWER_CACHE[degree] = coefficients
    return coefficients


def panel_chebyshev_moments(left: int, right: int, degree: int) -> list[mp.mpf]:
    a = mp.mpf(left)
    b = mp.mpf(right)
    center = (a + b) / 2
    half_width = (b - a) / 2
    alpha = mp.mpf(DEFAULT_ALPHA)
    gamma_alpha = mp.gamma(alpha + 1)
    y_moments = [
        mp.gammainc(alpha + p + 1, a, b) / gamma_alpha for p in range(degree + 1)
    ]
    t_moments: list[mp.mpf] = []
    for m in range(degree + 1):
        t_moments.append(
            mp.fsum(
                mp.binomial(m, p) * ((-center) ** (m - p)) * y_moments[p]
                for p in range(m + 1)
            )
            / (half_width**m)
        )
    chebyshev_moments = []
    for k in range(degree + 1):
        coefficients = chebyshev_power_coefficients(k)
        chebyshev_moments.append(
            mp.fsum(coefficients[m] * t_moments[m] for m in range(len(coefficients)))
        )
    return chebyshev_moments


def panel_integral(func, left: int, right: int, degree: int) -> mp.mpf:
    n_nodes = degree + 1
    a = mp.mpf(left)
    b = mp.mpf(right)
    theta_values: list[mp.mpf] = []
    function_values: list[mp.mpf] = []
    for j in range(n_nodes):
        theta = mp.pi * (mp.mpf(j) + mp.mpf("0.5")) / n_nodes
        y = (a + b) / 2 + (b - a) * mp.cos(theta) / 2
        theta_values.append(theta)
        function_values.append(func(y))
    coefficients = [
        (mp.mpf(2) / n_nodes)
        * mp.fsum(function_values[j] * mp.cos(k * theta_values[j]) for j in range(n_nodes))
        for k in range(n_nodes)
    ]
    moments = panel_chebyshev_moments(left, right, degree)
    return coefficients[0] * moments[0] / 2 + mp.fsum(
        coefficients[k] * moments[k] for k in range(1, n_nodes)
    )


def build_ladder_rows(value_core, derivative_core) -> tuple[list[dict], dict[int, tuple[mp.mpf, mp.mpf]]]:
    estimates: dict[int, tuple[mp.mpf, mp.mpf]] = {}
    ladder_rows: list[dict] = []
    for degree in DEFAULT_DEGREES:
        value_panel_contributions = [
            panel_integral(value_core, left, right, degree) for left, right in DEFAULT_PANELS
        ]
        derivative_panel_contributions = [
            panel_integral(derivative_core, left, right, degree) for left, right in DEFAULT_PANELS
        ]
        value_total = mp.fsum(value_panel_contributions)
        derivative_total = mp.fsum(derivative_panel_contributions)
        estimates[degree] = (value_total, derivative_total)
        ladder_rows.append(
            {
                "degree": degree,
                "panel_count": len(DEFAULT_PANELS),
                "value_integral_estimate": mp_text(value_total, 60),
                "derivative_integral_estimate": mp_text(derivative_total, 60),
                "proof_boundary": (
                    "High-precision floating Chebyshev panel-moment estimate only; not an Arb "
                    "coefficient-ball integration certificate."
                ),
            }
        )
    return ladder_rows, estimates


def build_cauchy_rows(
    estimates: dict[int, tuple[mp.mpf, mp.mpf]],
    value_cap: mp.mpf,
    derivative_cap: mp.mpf,
) -> list[dict]:
    rows: list[dict] = []
    for lower, upper in zip(DEFAULT_DEGREES, DEFAULT_DEGREES[1:]):
        value_delta = abs(estimates[upper][0] - estimates[lower][0])
        derivative_delta = abs(estimates[upper][1] - estimates[lower][1])
        value_ratio = value_delta / value_cap
        derivative_ratio = derivative_delta / derivative_cap
        rows.append(
            {
                "pair": f"{lower}->{upper}",
                "lower_degree": lower,
                "upper_degree": upper,
                "value_delta": mp_text(value_delta, 50),
                "derivative_delta": mp_text(derivative_delta, 50),
                "value_delta_to_cap_ratio": mp_text(value_ratio, 40),
                "derivative_delta_to_cap_ratio": mp_text(derivative_ratio, 40),
                "both_deltas_below_caps": bool(value_ratio < 1 and derivative_ratio < 1),
                "proof_boundary": (
                    "Cauchy-ladder calibration only; a small consecutive-degree delta is not an "
                    "interpolation-remainder proof."
                ),
            }
        )
    return rows


def build_panel_rows(value_core, derivative_core, degree: int) -> list[dict]:
    rows: list[dict] = []
    for left, right in DEFAULT_PANELS:
        rows.append(
            {
                "panel": f"{left}<=y<={right}",
                "degree": degree,
                "value_contribution_estimate": mp_text(panel_integral(value_core, left, right, degree), 50),
                "derivative_contribution_estimate": mp_text(
                    panel_integral(derivative_core, left, right, degree), 50
                ),
                "proof_boundary": "Panel contribution estimate only; not a rigorous panel enclosure.",
            }
        )
    return rows


def build_diagnostics(compact_scout: dict, quadrature_route: dict) -> dict:
    configure_precision()
    ratios = ratio_midpoints()
    value_core, derivative_core = make_core_functions(ratios)
    value_cap = mp.mpf(quadrature_route["summary"]["value_unscaled_expectation_error_cap"])
    derivative_cap = mp.mpf(quadrature_route["summary"]["derivative_unscaled_expectation_error_cap"])
    ladder_rows, estimates = build_ladder_rows(value_core, derivative_core)
    cauchy_rows = build_cauchy_rows(estimates, value_cap, derivative_cap)
    reference_degree = DEFAULT_DEGREES[-1]
    panel_rows = build_panel_rows(value_core, derivative_core, reference_degree)
    cap_safe_pairs = [row for row in cauchy_rows if row["both_deltas_below_caps"]]
    first_cap_safe_pair = cap_safe_pairs[0]["pair"] if cap_safe_pairs else None
    raw_summary = compact_scout["summary"]
    return {
        "T": DEFAULT_T,
        "index": DEFAULT_INDEX,
        "u": "1/10000",
        "alpha": DEFAULT_ALPHA,
        "mpmath_dps": DEFAULT_MPMATH_DPS,
        "flint_precision_bits_for_coefficients": DEFAULT_FLINT_PRECISION_BITS,
        "polynomial_M": DEFAULT_POLYNOMIAL_M,
        "ratio_cutoff_n": DEFAULT_RATIO_CUTOFF_N,
        "phi_term_count": DEFAULT_PHI_TERM_COUNT,
        "panel_count": len(DEFAULT_PANELS),
        "degree_count": len(DEFAULT_DEGREES),
        "degrees": list(DEFAULT_DEGREES),
        "cauchy_pair_count": len(cauchy_rows),
        "cap_safe_pair_count": len(cap_safe_pairs),
        "first_cap_safe_pair": first_cap_safe_pair,
        "reference_degree": reference_degree,
        "compact_interval": raw_summary["compact_interval"],
        "raw_value_width_to_cap_ratio_upper": raw_summary["value_raw_width_to_cap_ratio_upper"],
        "raw_derivative_width_to_cap_ratio_upper": raw_summary[
            "derivative_raw_width_to_cap_ratio_upper"
        ],
        "value_unscaled_expectation_error_cap": quadrature_route["summary"][
            "value_unscaled_expectation_error_cap"
        ],
        "derivative_unscaled_expectation_error_cap": quadrature_route["summary"][
            "derivative_unscaled_expectation_error_cap"
        ],
        "reference_value_integral_estimate": mp_text(estimates[reference_degree][0], 60),
        "reference_derivative_integral_estimate": mp_text(estimates[reference_degree][1], 60),
        "ladder_rows": ladder_rows,
        "cauchy_rows": cauchy_rows,
        "panel_rows": panel_rows,
        "recommended_upgrade": (
            "Promote this scout only by replacing floating Chebyshev coefficients with Arb coefficient "
            "balls and by adding a rigorous interpolation or Taylor remainder bound on each panel."
        ),
        "target_closing": False,
    }


def build_rows(paths: dict[str, str], diagnostics: dict) -> list[dict]:
    rows = [
        MatrixRow(
            id="nlrgwrcpms_01_raw_panel_rejection_import",
            role="scope_reduction",
            readiness="diagnostic_only",
            claim=(
                "The prior compact-interval scout rejects raw Arb interval-Riemann hulls on 0<=y<=200, "
                "so the next route must use local panel models rather than direct hull evaluation."
            ),
            diagnostics={
                "compact_interval": diagnostics["compact_interval"],
                "raw_value_width_to_cap_ratio_upper": diagnostics[
                    "raw_value_width_to_cap_ratio_upper"
                ],
                "raw_derivative_width_to_cap_ratio_upper": diagnostics[
                    "raw_derivative_width_to_cap_ratio_upper"
                ],
            },
            source_artifacts=[paths["compact_scout_json"], paths["compact_scout_note"]],
            proof_boundary="Scope import only; not a compact interval-integration proof.",
        ),
        MatrixRow(
            id="nlrgwrcpms_02_panel_moment_formula",
            role="floating_panel_moment_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "Each Chebyshev interpolant is integrated termwise against the normalized Gamma weight by "
                "expanding T_k((2y-a-b)/(b-a)) into monomials and using incomplete-Gamma panel moments."
            ),
            diagnostics={
                "panel_count": diagnostics["panel_count"],
                "degrees": diagnostics["degrees"],
                "formula": (
                    "Integral T_k(t(y))*y^alpha*exp(-y)/Gamma(alpha+1) dy is reduced to a finite "
                    "linear combination of incomplete-Gamma moments."
                ),
            },
            source_artifacts=[paths["quadrature_route_note"], paths["finite_part_note"]],
            proof_boundary=(
                "Moment-formula diagnostic only; the Chebyshev coefficients and interpolation remainder "
                "are not Arb-certified."
            ),
        ),
        MatrixRow(
            id="nlrgwrcpms_03_degree_ladder",
            role="floating_panel_moment_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "The degree ladder 12,16,20,24,32 gives stable high-precision floating panel-moment "
                "estimates for the compact value and derivative integrals."
            ),
            diagnostics={"ladder_rows": diagnostics["ladder_rows"]},
            source_artifacts=[paths["cancellation_grid_note"], paths["quadrature_route_json"]],
            proof_boundary="Floating degree-ladder diagnostic only; not an interval enclosure.",
        ),
        MatrixRow(
            id="nlrgwrcpms_04_cauchy_delta_vs_caps",
            role="floating_panel_moment_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "Consecutive Chebyshev panel-moment deltas from degree 16->20 onward are below the "
                "unscaled quadrature route caps in both value and derivative channels."
            ),
            diagnostics={
                "cauchy_rows": diagnostics["cauchy_rows"],
                "cap_safe_pair_count": diagnostics["cap_safe_pair_count"],
                "first_cap_safe_pair": diagnostics["first_cap_safe_pair"],
            },
            source_artifacts=[paths["quadrature_route_json"], paths["quadrature_route_note"]],
            proof_boundary=(
                "Cauchy convergence calibration only; small floating deltas are not a proof of quadrature "
                "or interpolation error."
            ),
        ),
        MatrixRow(
            id="nlrgwrcpms_05_reference_panel_contributions",
            role="floating_panel_moment_diagnostic",
            readiness="diagnostic_only",
            claim=(
                "At reference degree 32, the compact integral estimate is concentrated in the central "
                "Gamma-mass panels and is stable enough to guide a formal Arb panel partition."
            ),
            diagnostics={
                "reference_degree": diagnostics["reference_degree"],
                "reference_value_integral_estimate": diagnostics["reference_value_integral_estimate"],
                "reference_derivative_integral_estimate": diagnostics[
                    "reference_derivative_integral_estimate"
                ],
                "panel_rows": diagnostics["panel_rows"],
            },
            source_artifacts=[paths["far_tail_note"], paths["compact_scout_note"]],
            proof_boundary="Reference floating estimate only; not a certified integral value.",
        ),
        MatrixRow(
            id="nlrgwrcpms_06_promotion_rejected",
            role="rejected_route",
            readiness="not_ready_to_apply",
            claim="The floating Chebyshev panel-moment ladder proves the compact-interval integration source.",
            gap=(
                "The ladder is numerical calibration. A proof still needs Arb coefficient balls plus a "
                "rigorous interpolation or Taylor remainder bound on every panel."
            ),
            source_artifacts=[paths["intervalization_target"], paths["uniform_remainder_target"]],
            proof_boundary="Rejected promotion only; not a compact interval certificate, RH, or Lambda <= 0.",
        ),
        MatrixRow(
            id="nlrgwrcpms_07_acceptance_gate",
            role="acceptance_gate",
            readiness="not_ready_to_apply",
            claim=(
                "A promoted compact-interval certificate may reuse this partition and degree ladder only "
                "after replacing all floating components with interval-certified local models."
            ),
            diagnostics={"recommended_upgrade": diagnostics["recommended_upgrade"]},
            source_artifacts=[paths["dependency_graph"], paths["intervalization_target"]],
            proof_boundary="Proof-hygiene gate only; not a finite-grid certificate, RH, or Lambda <= 0.",
        ),
    ]
    return [asdict(row) for row in rows]


def build_artifact(compact_scout_path: Path, quadrature_route_path: Path) -> dict:
    compact_scout = load_json(compact_scout_path)
    quadrature_route = load_json(quadrature_route_path)
    paths = source_paths(compact_scout_path, quadrature_route_path)
    diagnostics = build_diagnostics(compact_scout, quadrature_route)
    rows = build_rows(paths, diagnostics)
    summary = {
        "matrix_rows": len(rows),
        "T": diagnostics["T"],
        "index": diagnostics["index"],
        "u": diagnostics["u"],
        "alpha": diagnostics["alpha"],
        "mpmath_dps": diagnostics["mpmath_dps"],
        "flint_precision_bits_for_coefficients": diagnostics["flint_precision_bits_for_coefficients"],
        "polynomial_M": diagnostics["polynomial_M"],
        "ratio_cutoff_n": diagnostics["ratio_cutoff_n"],
        "phi_term_count": diagnostics["phi_term_count"],
        "panel_count": diagnostics["panel_count"],
        "degree_count": diagnostics["degree_count"],
        "degrees": diagnostics["degrees"],
        "cauchy_pair_count": diagnostics["cauchy_pair_count"],
        "cap_safe_pair_count": diagnostics["cap_safe_pair_count"],
        "first_cap_safe_pair": diagnostics["first_cap_safe_pair"],
        "reference_degree": diagnostics["reference_degree"],
        "compact_interval": diagnostics["compact_interval"],
        "value_unscaled_expectation_error_cap": diagnostics["value_unscaled_expectation_error_cap"],
        "derivative_unscaled_expectation_error_cap": diagnostics[
            "derivative_unscaled_expectation_error_cap"
        ],
        "reference_value_integral_estimate": diagnostics["reference_value_integral_estimate"],
        "reference_derivative_integral_estimate": diagnostics[
            "reference_derivative_integral_estimate"
        ],
        "ready_to_apply_rows": 0,
        "target_closing": diagnostics["target_closing"],
        "recommended_upgrade": diagnostics["recommended_upgrade"],
        "main_finding": (
            "The floating Chebyshev panel-moment route is numerically plausible: on the six-panel "
            "compact interval 0<=y<=200, consecutive degree deltas from 16->20 onward are below the "
            "unscaled quadrature caps in both channels, with reference degree-32 estimates "
            f"{diagnostics['reference_value_integral_estimate']} and "
            f"{diagnostics['reference_derivative_integral_estimate']}. This is calibration only; "
            "a referee-usable proof still needs Arb coefficient balls and panel remainder bounds."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout",
        "date": "2026-07-07",
        "status": "worst-row floating Chebyshev panel-moment scout",
        "source_compact_interval_integration_scout": paths["compact_scout_note"],
        "source_compact_interval_integration_scout_json": paths["compact_scout_json"],
        "source_quadrature_remainder_route_matrix": paths["quadrature_route_note"],
        "source_quadrature_remainder_route_json": paths["quadrature_route_json"],
        "source_far_tail_split_certificate": paths["far_tail_note"],
        "source_finite_part_weighted_sum_interval_certificate": paths["finite_part_note"],
        "source_cancellation_reduced_grid_scout": paths["cancellation_grid_note"],
        "source_intervalization_target": paths["intervalization_target"],
        "source_uniform_remainder_target": paths["uniform_remainder_target"],
        "source_dependency_graph": paths["dependency_graph"],
        "generator": (
            "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.py"
        ),
        "checker": (
            "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.py"
        ),
        "proof_boundary": (
            "Worst-row floating Chebyshev panel-moment scout only. It calibrates a local-model route on "
            "0<=y<=200 using high-precision floating Chebyshev coefficients and incomplete-Gamma moment "
            "formulas, but it does not provide Arb coefficient balls, does not prove interpolation "
            "remainder bounds, does not prove a compact interval-integration certificate, does not prove "
            "a quadrature-remainder theorem, does not cover all rows or quadrature orders, does not prove "
            "a finite-grid interval certificate, does not prove a uniform collar theorem, and does not "
            "prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply for the full intervalization target.",
            "The Chebyshev panel-moment ladder is floating calibration only.",
            "Small consecutive-degree deltas are not promoted to interpolation-remainder bounds.",
            "The compact interval 0<=y<=200 remains open until Arb local models and remainders are certified.",
            "All-row coverage, rounding aggregation, and grid-to-collar coverage remain separate.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def result_line(artifact: dict) -> str:
    summary = artifact["summary"]
    return (
        "validated Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment "
        f"scout: {summary['matrix_rows']} rows, 0 issues, {summary['degree_count']} degrees, "
        f"{summary['cauchy_pair_count']} Cauchy pairs, {summary['cap_safe_pair_count']} cap-safe pairs, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Worst-Row Chebyshev Panel-Moment Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: worst-row floating Chebyshev panel-moment scout. This is not a proof",
        "of a compact interval-integration certificate, quadrature-remainder",
        "theorem, finite-grid interval certificate, uniform collar theorem, RH,",
        "or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout`.",
        "",
        "Proof boundary: this artifact is high-precision floating calibration",
        "for a local panel-model route on `0<=y<=200`. It does not provide",
        "Arb Chebyshev coefficient balls or interpolation-remainder bounds.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_worst_row_chebyshev_panel_moment_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(artifact),
        "```",
        "",
        "## Degree Ladder",
        "",
        "```text",
        f"T: {summary['T']}",
        f"index: F_{summary['index']}",
        f"compact interval: {summary['compact_interval']}",
        f"panel count: {summary['panel_count']}",
        f"degrees: {summary['degrees']}",
        f"value unscaled cap: {summary['value_unscaled_expectation_error_cap']}",
        f"derivative unscaled cap: {summary['derivative_unscaled_expectation_error_cap']}",
        f"first cap-safe pair: {summary['first_cap_safe_pair']}",
        f"cap-safe pair count: {summary['cap_safe_pair_count']}",
        f"reference degree: {summary['reference_degree']}",
        f"reference value integral estimate: {summary['reference_value_integral_estimate']}",
        f"reference derivative integral estimate: {summary['reference_derivative_integral_estimate']}",
        f"target closing: {summary['target_closing']}",
        "```",
        "",
        "Cauchy rows:",
        "",
        "```text",
    ]
    for row in artifact["matrix_rows"][3]["diagnostics"]["cauchy_rows"]:
        lines.extend(
            [
                f"{row['pair']}: value delta/cap {row['value_delta_to_cap_ratio']}; "
                f"derivative delta/cap {row['derivative_delta_to_cap_ratio']}; "
                f"below caps {row['both_deltas_below_caps']}",
            ]
        )
    lines.extend(
        [
            "```",
            "",
            "Required upgrade:",
            "",
            "```text",
            summary["recommended_upgrade"],
            "```",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_compact_interval_integration_scout"],
            artifact["source_compact_interval_integration_scout_json"],
            artifact["source_quadrature_remainder_route_matrix"],
            artifact["source_quadrature_remainder_route_json"],
            artifact["source_far_tail_split_certificate"],
            artifact["source_finite_part_weighted_sum_interval_certificate"],
            artifact["source_cancellation_reduced_grid_scout"],
            artifact["source_intervalization_target"],
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
    parser.add_argument("--compact-scout-json", type=Path, default=DEFAULT_COMPACT_SCOUT_JSON)
    parser.add_argument("--quadrature-route-json", type=Path, default=DEFAULT_QUADRATURE_ROUTE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    compact_scout_path = (
        args.compact_scout_json if args.compact_scout_json.is_absolute() else REPO_ROOT / args.compact_scout_json
    )
    quadrature_route_path = (
        args.quadrature_route_json
        if args.quadrature_route_json.is_absolute()
        else REPO_ROOT / args.quadrature_route_json
    )
    out_json = args.out_json if args.out_json.is_absolute() else REPO_ROOT / args.out_json
    note = args.note if args.note.is_absolute() else REPO_ROOT / args.note
    artifact = build_artifact(compact_scout_path, quadrature_route_path)
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian worst-row Chebyshev panel-moment scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
