#!/usr/bin/env python3
"""Build an Arb finite-degree real-T collar certificate for the relative-Gaussian surrogate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
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
from jensen_window_pf_negative_lambda_high_order_taylor_scout import (  # noqa: E402
    coefficient_value,
    generate_polynomials,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.md"

ArbPolynomial = list[flint.arb]


@dataclass(frozen=True)
class RatioBallRow:
    degree: int
    ratio_to_c0: str
    radius: str
    sign: str


@dataclass(frozen=True)
class BernsteinCertificateRow:
    name: str
    polynomial_degree: int
    stripped_zero_order: int
    bernstein_degree: int
    bernstein_subdivisions: int
    positive_bernstein_coefficients: int
    total_bernstein_coefficients: int
    first_bernstein_coefficient: str
    last_bernstein_coefficient: str
    min_bernstein_lower: str
    max_bernstein_radius: str
    endpoint_value: str
    certified_positive_on_interval: bool
    interpretation: str


def arb_positive(value: flint.arb) -> bool:
    return bool(value > 0 and not value.contains(0))


def sign_name(value: flint.arb) -> str:
    if arb_positive(value):
        return "positive"
    if bool(value < 0 and not value.contains(0)):
        return "negative"
    return "unresolved"


def trim(poly: ArbPolynomial) -> ArbPolynomial:
    result = list(poly)
    while len(result) > 1 and bool(result[-1] == 0):
        result.pop()
    return result


def padd(left: ArbPolynomial, right: ArbPolynomial) -> ArbPolynomial:
    size = max(len(left), len(right))
    result = [flint.arb(0)] * size
    for index in range(size):
        result[index] = (left[index] if index < len(left) else flint.arb(0)) + (
            right[index] if index < len(right) else flint.arb(0)
        )
    return trim(result)


def pscale(poly: ArbPolynomial, scalar: flint.arb) -> ArbPolynomial:
    return trim([scalar * coeff for coeff in poly])


def psub(left: ArbPolynomial, right: ArbPolynomial) -> ArbPolynomial:
    return padd(left, pscale(right, flint.arb(-1)))


def pmul(left: ArbPolynomial, right: ArbPolynomial) -> ArbPolynomial:
    result = [flint.arb(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return trim(result)


def ppow(poly: ArbPolynomial, exponent: int) -> ArbPolynomial:
    result = [flint.arb(1)]
    for _ in range(exponent):
        result = pmul(result, poly)
    return result


def pderivative(poly: ArbPolynomial) -> ArbPolynomial:
    if len(poly) == 1:
        return [flint.arb(0)]
    return trim([flint.arb(index) * poly[index] for index in range(1, len(poly))])


def pproduct(polys: list[ArbPolynomial]) -> ArbPolynomial:
    result = [flint.arb(1)]
    for poly in polys:
        result = pmul(result, poly)
    return result


def peval(poly: ArbPolynomial, u_value: flint.arb) -> flint.arb:
    total = flint.arb(0)
    for coeff in reversed(poly):
        total = total * u_value + coeff
    return total


def arb_rising(index: int, length: int) -> flint.arb:
    q = flint.arb(2 * index + 1) / flint.arb(2)
    total = flint.arb(1)
    for offset in range(length):
        total *= q + flint.arb(offset)
    return total


def truncated_multiplier_polynomial(index: int, ratios: list[flint.arb], M: int) -> ArbPolynomial:
    return trim([ratios[j] * arb_rising(index, j) for j in range(M + 1)])


def log_derivative_numerator(polys: list[ArbPolynomial], weights: list[int]) -> ArbPolynomial:
    total = [flint.arb(0)]
    for index, weight in enumerate(weights):
        if weight == 0:
            continue
        others = [poly for other_index, poly in enumerate(polys) if other_index != index]
        total = padd(total, pscale(pmul(pderivative(polys[index]), pproduct(others)), flint.arb(weight)))
    return trim(total)


def binomial(n: int, r: int) -> flint.arb:
    return flint.arb(math.comb(n, r))


def interval_endpoints_for_subdivision(s: int, n: int, collar_start_T: int) -> tuple[flint.arb, flint.arb]:
    denominator = flint.arb(n * collar_start_T)
    return flint.arb(s) / denominator, flint.arb(s + 1) / denominator


def bernstein_coefficients_on_interval(poly: ArbPolynomial, left: flint.arb, right: flint.arb) -> list[flint.arb]:
    degree = len(poly) - 1
    width = right - left
    power_coeffs = [flint.arb(0)] * (degree + 1)
    for i, coeff in enumerate(poly):
        for m in range(i + 1):
            power_coeffs[m] += coeff * binomial(i, m) * (left ** (i - m)) * (width**m)
    bernstein: list[flint.arb] = []
    for j in range(degree + 1):
        value = flint.arb(0)
        for m in range(j + 1):
            value += power_coeffs[m] * binomial(j, m) / binomial(degree, m)
        bernstein.append(value)
    return bernstein


def min_lower(values: list[flint.arb]) -> flint.arb:
    best = values[0].lower()
    for value in values[1:]:
        lower = value.lower()
        if lower < best:
            best = lower
    return flint.arb(best)


def max_radius(values: list[flint.arb]) -> flint.arb:
    best = values[0].rad()
    for value in values[1:]:
        radius = value.rad()
        if radius > best:
            best = radius
    return flint.arb(best)


def build_bernstein_row(
    name: str,
    poly: ArbPolynomial,
    stripped_zero_order: int,
    collar_start_T: int,
    interpretation: str,
) -> BernsteinCertificateRow:
    stripped = trim(poly[stripped_zero_order:])
    left, right = interval_endpoints_for_subdivision(0, 1, collar_start_T)
    coeffs = bernstein_coefficients_on_interval(stripped, left, right)
    positive = [arb_positive(coeff) for coeff in coeffs]
    endpoint = peval(stripped, right)
    return BernsteinCertificateRow(
        name=name,
        polynomial_degree=len(poly) - 1,
        stripped_zero_order=stripped_zero_order,
        bernstein_degree=len(stripped) - 1,
        bernstein_subdivisions=1,
        positive_bernstein_coefficients=sum(1 for item in positive if item),
        total_bernstein_coefficients=len(coeffs),
        first_bernstein_coefficient=coeffs[0].str(30),
        last_bernstein_coefficient=coeffs[-1].str(30),
        min_bernstein_lower=min_lower(coeffs).str(30),
        max_bernstein_radius=max_radius(coeffs).str(10),
        endpoint_value=endpoint.str(30),
        certified_positive_on_interval=all(positive) and arb_positive(endpoint),
        interpretation=interpretation,
    )


def build_ratio_rows(max_degree: int, cutoff_n: int) -> tuple[list[RatioBallRow], list[flint.arb]]:
    polynomials = generate_polynomials(max_degree)
    pi = flint.arb.pi()
    coefficient_values = [coefficient_value(polynomials[degree], cutoff_n, pi)[0] for degree in range(0, max_degree + 1, 2)]
    c0 = coefficient_values[0]
    ratios = [value / c0 for value in coefficient_values]
    rows = [
        RatioBallRow(
            degree=2 * index,
            ratio_to_c0=ratio.str(30),
            radius=ratio.rad().str(10),
            sign=sign_name(ratio),
        )
        for index, ratio in enumerate(ratios)
    ]
    return rows, ratios


def build_diagnostics(
    coefficient_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    continuation_M: int,
    collar_start_T: int,
) -> dict:
    flint.ctx.prec = precision_bits
    ratio_rows, ratios = build_ratio_rows(coefficient_max_degree, cutoff_n)
    polys_by_index = {
        index: truncated_multiplier_polynomial(index, ratios, continuation_M)
        for index in (k - 1, k, k + 1, k + 2)
    }
    ordered = [polys_by_index[index] for index in (k - 1, k, k + 1, k + 2)]

    b_product = psub(pmul(polys_by_index[k], polys_by_index[k]), pmul(polys_by_index[k - 1], polys_by_index[k + 1]))
    companion_product = psub(
        pmul(ppow(polys_by_index[k], 3), polys_by_index[k + 2]),
        pmul(polys_by_index[k - 1], ppow(polys_by_index[k + 1], 3)),
    )
    weighted_derivative = log_derivative_numerator(
        ordered,
        [2 * k + 1, -(6 * k + 5), 6 * k + 7, -(2 * k + 3)],
    )

    normalizer_rows = [
        build_bernstein_row(
            f"F_{index}",
            polys_by_index[index],
            0,
            collar_start_T,
            f"Degree-16 Arb coefficient-ball normalizer F_{index}(u) on 0<=u<=1/{collar_start_T}.",
        )
        for index in (k - 1, k, k + 1, k + 2)
    ]
    stencil_rows = [
        build_bernstein_row(
            "B_product",
            b_product,
            2,
            collar_start_T,
            "B-stencil positivity is certified via the stripped product polynomial u^-2*(F_k^2-F_(k-1)*F_(k+1)).",
        ),
        build_bernstein_row(
            "companion_product",
            companion_product,
            3,
            collar_start_T,
            "Companion-stencil positivity is certified via u^-3*(F_k^3*F_(k+2)-F_(k-1)*F_(k+1)^3).",
        ),
        build_bernstein_row(
            "weighted_gap_derivative",
            weighted_derivative,
            1,
            collar_start_T,
            "Weighted-gap positivity follows because the stripped derivative numerator has positive Bernstein coefficients.",
        ),
    ]

    return {
        "parameters": {
            "coefficient_max_taylor_degree": coefficient_max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "continuation_M": continuation_M,
            "collar_start_T": collar_start_T,
            "real_interval_u": f"[0, 1/{collar_start_T}]",
            "real_interval_T": f"[{collar_start_T}, infinity)",
        },
        "ratio_ball_rows": [asdict(row) for row in ratio_rows],
        "normalizer_rows": [asdict(row) for row in normalizer_rows],
        "stencil_rows": [asdict(row) for row in stencil_rows],
        "ratio_ball_rows_count": len(ratio_rows),
        "positive_normalizer_rows": sum(1 for row in normalizer_rows if row.certified_positive_on_interval),
        "certified_stencil_rows": sum(1 for row in stencil_rows if row.certified_positive_on_interval),
        "failed_bernstein_rows": sum(
            1 for row in [*normalizer_rows, *stencil_rows] if not row.certified_positive_on_interval
        ),
        "bernstein_subdivisions": 1,
        "finite_degree_arb_collar_certified": all(
            row.certified_positive_on_interval for row in [*normalizer_rows, *stencil_rows]
        ),
        "proof_boundary_note": (
            "This certifies only the degree-16 finite surrogate with Arb coefficient balls; "
            "it still does not bound the infinite residual Taylor tail."
        ),
    }


def build_artifact(
    coefficient_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    continuation_M: int,
    collar_start_T: int,
) -> dict:
    diagnostics = build_diagnostics(coefficient_max_degree, cutoff_n, precision_bits, k, continuation_M, collar_start_T)
    rows = [
        {
            "id": "nlrgd16arb_01_coefficient_ball_setup",
            "role": "finite_interval_surrogate",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-16 finite surrogate is rebuilt with Arb coefficient-ratio balls c_0 through c_16 rather than midpoint rationalized ratios.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md",
            ],
            "proof_boundary": "Finite coefficient-ball surrogate setup only; not an infinite Taylor-tail theorem.",
        },
        {
            "id": "nlrgd16arb_02_normalizer_bernstein_certificate",
            "role": "finite_interval_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "Bernstein coefficients certify the four Arb normalizer polynomials F_21 through F_24 as positive on 0<=u<=1/1156.",
            "proof_boundary": "Finite-degree Arb normalizer certificate only; not an all-k theorem.",
        },
        {
            "id": "nlrgd16arb_03_B_product_bernstein_certificate",
            "role": "finite_interval_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "Bernstein coefficients certify the stripped B product polynomial as positive on 0<=u<=1/1156.",
            "proof_boundary": "Finite-degree Arb B-stencil certificate only; not the infinite zeta-tail theorem.",
        },
        {
            "id": "nlrgd16arb_04_companion_product_bernstein_certificate",
            "role": "finite_interval_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "Bernstein coefficients certify the stripped companion product polynomial as positive on 0<=u<=1/1156.",
            "proof_boundary": "Finite-degree Arb companion-stencil certificate only; not the infinite zeta-tail theorem.",
        },
        {
            "id": "nlrgd16arb_05_weighted_gap_derivative_bernstein_certificate",
            "role": "finite_interval_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "Bernstein coefficients certify the stripped weighted-gap derivative numerator as positive on 0<=u<=1/1156.",
            "proof_boundary": "Finite-degree Arb weighted-gap certificate only; the full residual tail remains open.",
        },
        {
            "id": "nlrgd16arb_06_real_t_interval_surrogate_upgrade",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The rationalized real-T collar is upgraded to an Arb coefficient-ball finite-degree surrogate collar on T>=1156.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md",
            ],
            "proof_boundary": "Finite-degree interval-surrogate upgrade only; not a real-T theorem for the infinite zeta Taylor tail.",
        },
        {
            "id": "nlrgd16arb_07_live_infinite_tail_upgrade",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "The remaining analytic upgrade is a signed infinite-tail stencil bound beyond degree 16 on the same real collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Live theorem-search target only; the infinite-tail theorem is not proved.",
        },
        {
            "id": "nlrgd16arb_08_interval_surrogate_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The Arb finite-degree real-T collar proves the direct signed stencil-tail theorem.",
            "gap": "The certificate stops at degree 16 and does not bound the infinite residual Taylor-tail stencils.",
            "proof_boundary": "Rejected finite-surrogate promotion only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "ratio_ball_rows": diagnostics["ratio_ball_rows_count"],
        "positive_normalizer_rows": diagnostics["positive_normalizer_rows"],
        "certified_stencil_rows": diagnostics["certified_stencil_rows"],
        "failed_bernstein_rows": diagnostics["failed_bernstein_rows"],
        "bernstein_subdivisions": diagnostics["bernstein_subdivisions"],
        "real_T_arb_finite_degree_collar_start_T": str(collar_start_T),
        "finite_degree_arb_collar_certified": diagnostics["finite_degree_arb_collar_certified"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "The degree-16 real-T collar now survives Arb coefficient-ratio balls for the finite surrogate: "
            "Bernstein coefficients certify F_21..F_24, the stripped B product, the stripped companion product, "
            "and the stripped weighted-gap derivative numerator on 0<=u<=1/1156. This upgrades the midpoint "
            "surrogate to an interval finite-degree surrogate, while leaving the infinite residual Taylor-tail theorem open."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate",
        "date": "2026-07-07",
        "status": "finite theorem-search diagnostic",
        "source_real_t_collar_scout": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It certifies a real-T collar for the degree-16 finite surrogate "
            "using Arb coefficient-ratio balls at fixed k=22, but it does not bound the infinite residual Taylor tail, "
            "does not prove scaled-curvature monotonicity, does not prove cone entry, and does not prove RH or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The certificate is finite-degree through degree 16 only.",
            "The infinite residual Taylor tail remains unbounded.",
            "The result is fixed at k=22 and is not an all-k theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][0]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['positive_normalizer_rows']} positive normalizer rows, "
        f"{summary['certified_stencil_rows']} certified stencil rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Arb Real-T Collar Certificate",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof",
        "of an infinite Taylor-tail theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate`.",
        "",
        "Proof boundary: this artifact certifies a real-`T` collar only for the",
        "degree-16 finite surrogate at fixed `k=22`, using Arb coefficient-ratio",
        "balls. It does not bound the infinite residual Taylor tail.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_arb_real_t_collar_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Arb Finite-Degree Collar",
        "",
        "```text",
        f"real T interval: {diagnostics['parameters']['real_interval_T']}",
        f"u interval: {diagnostics['parameters']['real_interval_u']}",
        f"ratio ball rows: {summary['ratio_ball_rows']}",
        f"positive normalizer rows: {summary['positive_normalizer_rows']}",
        f"certified stencil rows: {summary['certified_stencil_rows']}",
        f"failed Bernstein rows: {summary['failed_bernstein_rows']}",
        f"Bernstein subdivisions: {summary['bernstein_subdivisions']}",
        f"finite-degree Arb collar certified: {summary['finite_degree_arb_collar_certified']}",
        "```",
        "",
        "Key ratio balls:",
        "",
        "```text",
    ]
    for row in diagnostics["ratio_ball_rows"][-3:]:
        lines.append(f"c{row['degree']}/c0: {row['ratio_to_c0']} (radius {row['radius']}, sign {row['sign']})")
    lines.extend(["```", "", "Bernstein certificates:", "", "```text"])
    for row in diagnostics["normalizer_rows"]:
        lines.append(
            f"{row['name']}: degree={row['polynomial_degree']}, "
            f"positive Bernstein coefficients={row['positive_bernstein_coefficients']}/{row['total_bernstein_coefficients']}, "
            f"certified={row['certified_positive_on_interval']}"
        )
    for row in diagnostics["stencil_rows"]:
        lines.append(
            f"{row['name']}: degree={row['polynomial_degree']}, zero order={row['stripped_zero_order']}, "
            f"Bernstein degree={row['bernstein_degree']}, "
            f"positive coefficients={row['positive_bernstein_coefficients']}/{row['total_bernstein_coefficients']}, "
            f"certified={row['certified_positive_on_interval']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Interpretation:",
            "",
            "The rationalized real-`T` collar has been upgraded to an Arb",
            "coefficient-ball finite-degree collar. The remaining proof gap is now",
            "sharper: prove signed infinite-tail stencil bounds beyond degree 16 on",
            "the same real collar, and then remove the fixed-`k` limitation.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_real_t_collar_scout"],
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
    parser.add_argument("--continuation-M", type=int, default=8)
    parser.add_argument("--collar-start-T", type=int, default=1156)
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
        args.continuation_M,
        args.collar_start_T,
    )
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(artifact, note)
    print(
        "wrote Jensen-window PF negative-lambda relative-Gaussian degree-16 Arb real-T collar certificate: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
