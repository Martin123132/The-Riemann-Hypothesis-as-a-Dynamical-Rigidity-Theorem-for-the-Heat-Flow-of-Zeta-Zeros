#!/usr/bin/env python3
"""Build a real-T collar scout for the degree-16 relative-Gaussian surrogate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
import json
from pathlib import Path
import sys
from time import perf_counter


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
import sympy as sp  # noqa: E402

from jensen_window_pf_heat_flow_monotone_closure_scout import (  # noqa: E402
    REPO_ROOT,
    decimal_format,
)
from jensen_window_pf_negative_lambda_high_order_taylor_scout import (  # noqa: E402
    arb_mid_decimal,
    build_diagnostics as build_high_order_diagnostics,
)
from jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan import (  # noqa: E402
    log_stencils,
)


getcontext().prec = 100

DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.md"

Polynomial = list[Fraction]


@dataclass(frozen=True)
class NormalizerRow:
    index: int
    polynomial_degree: int
    zero_order_at_u0: int
    stripped_degree: int
    roots_in_open_interval: int
    sign_at_u0_after_stripping: str
    sign_at_endpoint: str
    value_at_endpoint: str
    certified_positive_on_interval: bool


@dataclass(frozen=True)
class ProductStencilRow:
    name: str
    polynomial_degree: int
    zero_order_at_u0: int
    stripped_degree: int
    roots_in_open_interval: int
    sign_at_u0_after_stripping: str
    sign_at_endpoint: str
    value_at_endpoint: str
    certified_positive_on_interval: bool
    interpretation: str


@dataclass(frozen=True)
class DerivativeStencilRow:
    name: str
    numerator_degree: int
    zero_order_at_u0: int
    stripped_degree: int
    roots_in_open_interval: int
    sign_at_u0_after_stripping: str
    sign_at_endpoint: str
    numerator_value_at_endpoint: str
    certified_positive_derivative_on_interval: bool
    certified_positive_stencil_on_interval: bool
    interpretation: str


@dataclass(frozen=True)
class RootCountResult:
    roots_in_open_interval: int
    elapsed_seconds: str


def rationalize_decimal(value: Decimal) -> Fraction:
    return Fraction(str(value))


def rising(q: Fraction, length: int) -> Fraction:
    total = Fraction(1)
    for index in range(length):
        total *= q + index
    return total


def trim(poly: Polynomial) -> Polynomial:
    result = list(poly)
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def padd(left: Polynomial, right: Polynomial) -> Polynomial:
    size = max(len(left), len(right))
    result = [Fraction(0)] * size
    for index in range(size):
        result[index] = (left[index] if index < len(left) else Fraction(0)) + (
            right[index] if index < len(right) else Fraction(0)
        )
    return trim(result)


def pscale(poly: Polynomial, scalar: Fraction) -> Polynomial:
    return trim([scalar * value for value in poly])


def psub(left: Polynomial, right: Polynomial) -> Polynomial:
    return padd(left, pscale(right, Fraction(-1)))


def pmul(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return trim(result)


def ppow(poly: Polynomial, exponent: int) -> Polynomial:
    result = [Fraction(1)]
    for _ in range(exponent):
        result = pmul(result, poly)
    return result


def pderivative(poly: Polynomial) -> Polynomial:
    if len(poly) == 1:
        return [Fraction(0)]
    return trim([index * poly[index] for index in range(1, len(poly))])


def pproduct(polys: list[Polynomial]) -> Polynomial:
    result = [Fraction(1)]
    for poly in polys:
        result = pmul(result, poly)
    return result


def valuation(poly: Polynomial) -> int:
    for index, value in enumerate(poly):
        if value:
            return index
    return len(poly)


def strip_zero_at_origin(poly: Polynomial) -> tuple[Polynomial, int]:
    order = valuation(poly)
    if order >= len(poly):
        return [Fraction(0)], order
    return trim(poly[order:]), order


def sign_name(value: Fraction) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def peval(poly: Polynomial, u_value: Fraction) -> Fraction:
    total = Fraction(0)
    for coeff in reversed(poly):
        total = total * u_value + coeff
    return total


def fraction_decimal(value: Fraction) -> str:
    return decimal_format(Decimal(value.numerator) / Decimal(value.denominator))


def sympy_poly(poly: Polynomial) -> sp.Poly:
    u_symbol = sp.symbols("u")
    coeffs = [sp.Rational(value.numerator, value.denominator) for value in reversed(poly)]
    return sp.Poly.from_list(coeffs, gens=u_symbol, domain=sp.QQ)


def count_roots(poly: Polynomial, endpoint_T: int) -> RootCountResult:
    stripped, _order = strip_zero_at_origin(poly)
    start = perf_counter()
    count = sympy_poly(stripped).count_roots(sp.Rational(0), sp.Rational(1, endpoint_T))
    elapsed = perf_counter() - start
    return RootCountResult(roots_in_open_interval=int(count), elapsed_seconds=f"{elapsed:.3f}")


def truncated_multiplier_polynomial(k: int, ratios: list[Fraction], M: int) -> Polynomial:
    q = Fraction(2 * k + 1, 2)
    return trim([ratios[j] * rising(q, j) for j in range(M + 1)])


def log_derivative_numerator(polys: list[Polynomial], weights: list[int]) -> Polynomial:
    total = [Fraction(0)]
    for index, weight in enumerate(weights):
        if weight == 0:
            continue
        others = [poly for other_index, poly in enumerate(polys) if other_index != index]
        numerator_piece = pmul(pderivative(polys[index]), pproduct(others))
        total = padd(total, pscale(numerator_piece, Fraction(weight)))
    return trim(total)


def build_normalizer_row(index: int, poly: Polynomial, endpoint_T: int) -> NormalizerRow:
    endpoint = Fraction(1, endpoint_T)
    stripped, order = strip_zero_at_origin(poly)
    roots = count_roots(poly, endpoint_T)
    sign0 = sign_name(stripped[0])
    sign_endpoint = sign_name(peval(stripped, endpoint))
    return NormalizerRow(
        index=index,
        polynomial_degree=len(poly) - 1,
        zero_order_at_u0=order,
        stripped_degree=len(stripped) - 1,
        roots_in_open_interval=roots.roots_in_open_interval,
        sign_at_u0_after_stripping=sign0,
        sign_at_endpoint=sign_endpoint,
        value_at_endpoint=fraction_decimal(peval(poly, endpoint)),
        certified_positive_on_interval=roots.roots_in_open_interval == 0 and sign0 == "positive" and sign_endpoint == "positive",
    )


def build_product_stencil_row(name: str, poly: Polynomial, endpoint_T: int, interpretation: str) -> ProductStencilRow:
    endpoint = Fraction(1, endpoint_T)
    stripped, order = strip_zero_at_origin(poly)
    roots = count_roots(poly, endpoint_T)
    sign0 = sign_name(stripped[0])
    sign_endpoint = sign_name(peval(stripped, endpoint))
    return ProductStencilRow(
        name=name,
        polynomial_degree=len(poly) - 1,
        zero_order_at_u0=order,
        stripped_degree=len(stripped) - 1,
        roots_in_open_interval=roots.roots_in_open_interval,
        sign_at_u0_after_stripping=sign0,
        sign_at_endpoint=sign_endpoint,
        value_at_endpoint=fraction_decimal(peval(poly, endpoint)),
        certified_positive_on_interval=roots.roots_in_open_interval == 0 and sign0 == "positive" and sign_endpoint == "positive",
        interpretation=interpretation,
    )


def build_derivative_stencil_row(name: str, numerator: Polynomial, endpoint_T: int, interpretation: str) -> DerivativeStencilRow:
    endpoint = Fraction(1, endpoint_T)
    stripped, order = strip_zero_at_origin(numerator)
    roots = count_roots(numerator, endpoint_T)
    sign0 = sign_name(stripped[0])
    sign_endpoint = sign_name(peval(stripped, endpoint))
    positive_derivative = roots.roots_in_open_interval == 0 and sign0 == "positive" and sign_endpoint == "positive"
    return DerivativeStencilRow(
        name=name,
        numerator_degree=len(numerator) - 1,
        zero_order_at_u0=order,
        stripped_degree=len(stripped) - 1,
        roots_in_open_interval=roots.roots_in_open_interval,
        sign_at_u0_after_stripping=sign0,
        sign_at_endpoint=sign_endpoint,
        numerator_value_at_endpoint=fraction_decimal(peval(numerator, endpoint)),
        certified_positive_derivative_on_interval=positive_derivative,
        certified_positive_stencil_on_interval=positive_derivative,
        interpretation=interpretation,
    )


def build_diagnostics(
    coefficient_max_degree: int,
    cutoff_n: int,
    precision_bits: int,
    k: int,
    continuation_M: int,
    collar_start_T: int,
) -> dict:
    high = build_high_order_diagnostics(coefficient_max_degree, cutoff_n, precision_bits, k, [1000, 2000])
    decimal_ratios = [arb_mid_decimal(flint.arb(row.ratio_to_c0)) for row in high.coefficient_rows]
    rational_ratios = [rationalize_decimal(value) for value in decimal_ratios]
    polys_by_index = {
        index: truncated_multiplier_polynomial(index, rational_ratios, continuation_M)
        for index in (k - 1, k, k + 1, k + 2)
    }
    ordered_polys = [polys_by_index[index] for index in (k - 1, k, k + 1, k + 2)]

    b_product = psub(pmul(polys_by_index[k], polys_by_index[k]), pmul(polys_by_index[k - 1], polys_by_index[k + 1]))
    companion_product = psub(
        pmul(ppow(polys_by_index[k], 3), polys_by_index[k + 2]),
        pmul(polys_by_index[k - 1], ppow(polys_by_index[k + 1], 3)),
    )
    weighted_derivative = log_derivative_numerator(
        ordered_polys,
        [2 * k + 1, -(6 * k + 5), 6 * k + 7, -(2 * k + 3)],
    )

    normalizer_rows = [
        build_normalizer_row(index, polys_by_index[index], collar_start_T)
        for index in (k - 1, k, k + 1, k + 2)
    ]
    product_rows = [
        build_product_stencil_row(
            "B",
            b_product,
            collar_start_T,
            "B=2*log(F_k)-log(F_(k-1))-log(F_(k+1)) is positive iff F_k^2-F_(k-1)*F_(k+1)>0.",
        ),
        build_product_stencil_row(
            "companion",
            companion_product,
            collar_start_T,
            "The companion stencil is positive iff F_k^3*F_(k+2)-F_(k-1)*F_(k+1)^3>0.",
        ),
    ]
    derivative_rows = [
        build_derivative_stencil_row(
            "weighted_gap",
            weighted_derivative,
            collar_start_T,
            "The weighted-gap log stencil starts at 0 as u->0; positive derivative on the interval implies positivity.",
        )
    ]
    endpoint_stencils = log_stencils(k, collar_start_T, decimal_ratios, continuation_M)
    if endpoint_stencils is None:
        raise RuntimeError("endpoint stencils are undefined")

    certified_normalizers = sum(1 for row in normalizer_rows if row.certified_positive_on_interval)
    certified_product_stencils = sum(1 for row in product_rows if row.certified_positive_on_interval)
    certified_derivative_stencils = sum(1 for row in derivative_rows if row.certified_positive_stencil_on_interval)
    root_count_failures = sum(
        row.roots_in_open_interval
        for row in [*normalizer_rows, *product_rows, *derivative_rows]
    )

    return {
        "parameters": {
            "coefficient_max_taylor_degree": coefficient_max_degree,
            "tail_cutoff_n": cutoff_n,
            "precision_bits": precision_bits,
            "tail_start_k": k,
            "continuation_M": continuation_M,
            "collar_start_T": collar_start_T,
            "real_interval_u": f"(0, 1/{collar_start_T}]",
            "real_interval_T": f"[{collar_start_T}, infinity)",
        },
        "surrogate_note": (
            "Coefficient ratios are rationalized from the degree-16 high-order Taylor midpoint decimals; "
            "this is an exact root-count analysis of that finite surrogate, not an Arb coefficient-ball theorem."
        ),
        "normalizer_rows": [asdict(row) for row in normalizer_rows],
        "product_stencil_rows": [asdict(row) for row in product_rows],
        "derivative_stencil_rows": [asdict(row) for row in derivative_rows],
        "endpoint_log_stencils": {
            "T": str(collar_start_T),
            "B": decimal_format(endpoint_stencils[0]),
            "companion": decimal_format(endpoint_stencils[1]),
            "weighted_gap": decimal_format(endpoint_stencils[2]),
        },
        "certified_positive_normalizer_rows": certified_normalizers,
        "certified_product_stencil_rows": certified_product_stencils,
        "certified_derivative_stencil_rows": certified_derivative_stencils,
        "certified_surrogate_stencil_rows": certified_product_stencils + certified_derivative_stencils,
        "root_count_failures": root_count_failures,
        "real_T_surrogate_collar_certified": (
            certified_normalizers == 4
            and certified_product_stencils == 2
            and certified_derivative_stencils == 1
            and root_count_failures == 0
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
            "id": "nlrgd16rt_01_surrogate_setup",
            "role": "finite_surrogate",
            "readiness": "not_ready_to_apply",
            "claim": "The degree-16 real-T collar scout converts the midpoint coefficient-ratio finite truncation into exact rational polynomial data in u=1/T.",
            "diagnostics": diagnostics,
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md",
                "outputs/jensen_window_pf_negative_lambda_high_order_taylor_scout.md",
            ],
            "proof_boundary": "Exact finite-surrogate setup only; not an Arb coefficient-ball theorem and not an infinite Taylor-tail theorem.",
        },
        {
            "id": "nlrgd16rt_02_normalizer_positivity",
            "role": "finite_surrogate_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "The four degree-16 surrogate normalizer polynomials F_(k-1), F_k, F_(k+1), F_(k+2) have no roots on 0<u<=1/1156 and are positive there.",
            "proof_boundary": "Finite surrogate real-T certificate only; not a statement about the full zeta tail.",
        },
        {
            "id": "nlrgd16rt_03_B_product_collar",
            "role": "finite_surrogate_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "The B product polynomial F_k^2-F_(k-1)F_(k+1), after removing its u=0 zero, has no roots on 0<u<=1/1156 and is positive there.",
            "proof_boundary": "Finite surrogate real-T B-stencil certificate only; not an all-k theorem.",
        },
        {
            "id": "nlrgd16rt_04_companion_product_collar",
            "role": "finite_surrogate_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "The companion product polynomial F_k^3F_(k+2)-F_(k-1)F_(k+1)^3, after removing its u=0 zero, has no roots on 0<u<=1/1156 and is positive there.",
            "proof_boundary": "Finite surrogate real-T companion-stencil certificate only; not an infinite-tail theorem.",
        },
        {
            "id": "nlrgd16rt_05_weighted_gap_derivative_collar",
            "role": "finite_surrogate_certificate",
            "readiness": "not_ready_to_apply",
            "claim": "The weighted-gap log-stencil derivative numerator has no roots on 0<u<=1/1156 and is positive there, so the weighted-gap surrogate stencil is positive on the real collar.",
            "proof_boundary": "Finite surrogate derivative certificate only; not a proof of the analytic weighted-gap tail.",
        },
        {
            "id": "nlrgd16rt_06_integer_scan_to_real_surrogate",
            "role": "finite_diagnostic",
            "readiness": "not_ready_to_apply",
            "claim": "The integer threshold T=1156 from the collar scan is upgraded to a real-T sign collar for the degree-16 rationalized finite surrogate.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md",
            ],
            "proof_boundary": "Finite surrogate upgrade only; not a real-T theorem for the infinite zeta Taylor tail.",
        },
        {
            "id": "nlrgd16rt_07_live_tail_upgrade",
            "role": "live_route",
            "readiness": "not_ready_to_apply",
            "claim": "The next proof task is to replace rationalized finite-surrogate signs with interval coefficient control and signed infinite-tail stencil bounds on the same real-T collar.",
            "source_artifacts": [
                "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
                "outputs/signed_hankel_jensen_dependency_graph.md",
            ],
            "proof_boundary": "Live theorem-search target only; the interval coefficient and infinite-tail upgrades remain open.",
        },
        {
            "id": "nlrgd16rt_08_surrogate_promotion_rejected",
            "role": "rejected_route",
            "readiness": "not_ready_to_apply",
            "claim": "The rationalized degree-16 finite-surrogate real-T collar proves the direct signed stencil-tail theorem.",
            "gap": "The result uses midpoint coefficient ratios, finite degree 16 truncation, fixed k=22, and no infinite residual Taylor-tail bound.",
            "proof_boundary": "Rejected finite-surrogate promotion only; not scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0.",
        },
    ]
    summary = {
        "matrix_rows": len(rows),
        "positive_normalizer_rows": diagnostics["certified_positive_normalizer_rows"],
        "certified_product_stencil_rows": diagnostics["certified_product_stencil_rows"],
        "certified_derivative_stencil_rows": diagnostics["certified_derivative_stencil_rows"],
        "certified_surrogate_stencil_rows": diagnostics["certified_surrogate_stencil_rows"],
        "root_count_failures": diagnostics["root_count_failures"],
        "real_T_surrogate_collar_start_T": str(collar_start_T),
        "real_T_surrogate_collar_certified": diagnostics["real_T_surrogate_collar_certified"],
        "ready_to_apply_rows": 0,
        "target_closing": False,
        "main_finding": (
            "For the rationalized degree-16 finite surrogate at k=22, the integer collar threshold T=1156 "
            "extends to the whole real half-line T>=1156: the four normalizers are positive, the B and "
            "companion product stencils have no stripped roots on 0<u<=1/1156, and the weighted-gap derivative "
            "numerator is positive there. This is a finite-surrogate real-T collar, not an infinite Taylor-tail theorem."
        ),
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout",
        "date": "2026-07-07",
        "status": "finite theorem-search diagnostic",
        "source_degree16_collar_scan": "outputs/jensen_window_pf_negative_lambda_relative_gaussian_degree16_collar_scan.md",
        "source_uniform_remainder_target": "outputs/jensen_window_pf_negative_lambda_uniform_remainder_target.md",
        "source_dependency_graph": "outputs/signed_hankel_jensen_dependency_graph.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py",
        "proof_boundary": (
            "Finite theorem-search diagnostic only. It certifies a real-T collar for a rationalized degree-16 finite "
            "surrogate at fixed k=22, but it does not prove an Arb coefficient-ball theorem, an infinite Taylor-tail "
            "theorem, scaled-curvature monotonicity, cone entry, RH, or Lambda <= 0."
        ),
        "matrix_rows": rows,
        "summary": summary,
        "invariants": [
            "No row is ready_to_apply.",
            "The real-T collar is only for the rationalized degree-16 finite surrogate.",
            "The result does not control the infinite residual Taylor tail.",
            "The result is fixed at k=22 and is not an all-k theorem.",
            "Endpoint PF, RH, Laguerre-Polya membership, and Lambda <= 0 are forbidden as inputs.",
        ],
    }


def write_note(artifact: dict, path: Path) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["matrix_rows"][0]["diagnostics"]
    result_line = (
        "validated Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: "
        f"{summary['matrix_rows']} rows, 0 issues, "
        f"{summary['positive_normalizer_rows']} positive normalizer rows, "
        f"{summary['certified_surrogate_stencil_rows']} certified surrogate stencil rows, "
        f"{summary['ready_to_apply_rows']} ready-to-apply rows"
    )
    lines = [
        "# Jensen-Window PF Negative-Lambda Relative-Gaussian Degree-16 Real-T Collar Scout",
        "",
        "Date: 2026-07-07",
        "",
        "Status: finite theorem-search diagnostic. This is not a proof",
        "of an infinite Taylor-tail theorem, scaled-curvature monotonicity,",
        "cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout`.",
        "",
        "Proof boundary: this artifact certifies a real-`T` collar only for the",
        "rationalized degree-16 finite surrogate at fixed `k=22`. It does not",
        "prove the corresponding analytic zeta-tail theorem.",
        "",
        "Machine-readable artifact:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.json",
        "```",
        "",
        "Generator:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py",
        "```",
        "",
        "Checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_relative_gaussian_degree16_real_t_collar_scout.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line,
        "```",
        "",
        "## Real-T Surrogate Collar",
        "",
        "```text",
        f"real T interval: {diagnostics['parameters']['real_interval_T']}",
        f"u interval: {diagnostics['parameters']['real_interval_u']}",
        f"positive normalizer rows: {summary['positive_normalizer_rows']}",
        f"certified product stencil rows: {summary['certified_product_stencil_rows']}",
        f"certified derivative stencil rows: {summary['certified_derivative_stencil_rows']}",
        f"certified surrogate stencil rows: {summary['certified_surrogate_stencil_rows']}",
        f"root-count failures: {summary['root_count_failures']}",
        f"real-T surrogate collar certified: {summary['real_T_surrogate_collar_certified']}",
        "```",
        "",
        "Endpoint log-stencils at `T=1156`:",
        "",
        "```text",
        f"B: {diagnostics['endpoint_log_stencils']['B']}",
        f"companion: {diagnostics['endpoint_log_stencils']['companion']}",
        f"weighted gap: {diagnostics['endpoint_log_stencils']['weighted_gap']}",
        "```",
        "",
        "Normalizer root counts:",
        "",
        "```text",
    ]
    for row in diagnostics["normalizer_rows"]:
        lines.append(
            f"F_{row['index']}: degree={row['polynomial_degree']}, roots={row['roots_in_open_interval']}, "
            f"endpoint sign={row['sign_at_endpoint']}, certified={row['certified_positive_on_interval']}"
        )
    lines.extend(["```", "", "Stencil root-count certificates:", "", "```text"])
    for row in diagnostics["product_stencil_rows"]:
        lines.append(
            f"{row['name']}: degree={row['polynomial_degree']}, zero order={row['zero_order_at_u0']}, "
            f"stripped degree={row['stripped_degree']}, roots={row['roots_in_open_interval']}, "
            f"certified={row['certified_positive_on_interval']}"
        )
    for row in diagnostics["derivative_stencil_rows"]:
        lines.append(
            f"{row['name']} derivative: degree={row['numerator_degree']}, zero order={row['zero_order_at_u0']}, "
            f"stripped degree={row['stripped_degree']}, roots={row['roots_in_open_interval']}, "
            f"certified={row['certified_positive_stencil_on_interval']}"
        )
    lines.extend(
        [
            "```",
            "",
            "Interpretation:",
            "",
            "The integer scan threshold `T=1156` is no longer just a sampled",
            "integer threshold for the degree-16 surrogate: within the rationalized",
            "finite model, all three structured stencil signs persist on the full",
            "real half-line `T>=1156`. The live proof route is now to replace this",
            "surrogate statement with interval coefficient control and signed",
            "infinite-tail stencil bounds.",
            "",
            "Integration:",
            "",
            "```text",
            artifact["source_degree16_collar_scan"],
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
        "wrote Jensen-window PF negative-lambda relative-Gaussian degree-16 real-T collar scout: "
        f"{out_json.relative_to(REPO_ROOT).as_posix()} and {note.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
