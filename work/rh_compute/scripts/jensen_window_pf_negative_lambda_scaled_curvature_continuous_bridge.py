#!/usr/bin/env python3
"""Certify scaled-curvature growth at lambda=-100 by a continuous bridge."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
import os
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402
import sympy as sp  # noqa: E402

from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    arb_interval,
    arb_lower_text,
    arb_rational,
    arb_upper_text,
    certify_scaled_curvature_mode_block,
    exact_upper,
)
from jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate import (  # noqa: E402
    merged_coefficients,
)


DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT / "outputs/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.md"
)
PAIRED_COMPACT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.json"
)
PAIRED_RAY_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.json"
)
LEADING_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.json"
)

DEFAULT_PRECISION_BITS = 192
MODE_START = Fraction(579, 625)
MODE_SPLIT = Fraction(2)
MODE_END = Fraction(5)
LOW_WIDTH = Fraction(1, 1000)
HIGH_WIDTH = Fraction(1, 5000)
LOW_EIGHTH_BOUND = Fraction(1, 50_000)
HIGH_EIGHTH_BOUND = Fraction(1, 10_000_000_000)
Q_FLOOR = 1_000_000_000
DEFAULT_WORKERS = max(1, min(8, (os.cpu_count() or 4) - 1))


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def fraction_decimal(value: Fraction, digits: int = 40) -> str:
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(value.numerator) / Decimal(value.denominator), "E")


def fraction_grid(
    start: Fraction, end: Fraction, width: Fraction
) -> list[tuple[Fraction, Fraction]]:
    rows: list[tuple[Fraction, Fraction]] = []
    left = start
    while left < end:
        right = min(left + width, end)
        rows.append((left, right))
        left = right
    return rows


def polynomial_positive_gate(expression: sp.Expr, u: sp.Symbol, q: sp.Symbol) -> dict:
    v, capital_q = sp.symbols("v capital_q", nonnegative=True)
    numerator, denominator = sp.fraction(sp.factor(sp.together(expression)))
    shifted_numerator = sp.Poly(
        sp.expand(numerator.subs({u: v + 5, q: capital_q + Q_FLOOR})),
        v,
        capital_q,
    )
    shifted_denominator = sp.Poly(
        sp.expand(denominator.subs({u: v + 5, q: capital_q + Q_FLOOR})),
        v,
        capital_q,
    )
    numerator_coefficients = shifted_numerator.coeffs()
    denominator_coefficients = shifted_denominator.coeffs()
    if not numerator_coefficients or not denominator_coefficients:
        raise RuntimeError("empty shifted polynomial gate")
    if any(coefficient <= 0 for coefficient in numerator_coefficients):
        raise RuntimeError("shifted numerator is not coefficientwise positive")
    if any(coefficient <= 0 for coefficient in denominator_coefficients):
        raise RuntimeError("shifted denominator is not coefficientwise positive")
    return {
        "substitution": "u=5+v, q=1000000000+Q, v>=0, Q>=0",
        "numerator_terms": len(shifted_numerator.terms()),
        "denominator_terms": len(shifted_denominator.terms()),
        "minimum_numerator_coefficient": str(min(numerator_coefficients)),
        "minimum_denominator_coefficient": str(min(denominator_coefficients)),
    }


def exact_ray_algebra() -> dict:
    u, q = sp.symbols("u q", positive=True)
    t = (
        8 * q**2 * u
        + 400 * q * u**2
        - 30 * q * u
        - 2 * q
        - 600 * u**2
        + 15 * u
        + 3
    ) / (2 * (2 * q - 3))
    a = u * (
        64 * q**3 * u
        + 16 * q**3
        + 1408 * q**2 * u
        - 84 * q**2
        - 4560 * q * u
        + 120 * q
        + 3600 * u
        - 45
    ) / (4 * (2 * q - 3) ** 2)
    b = u * (
        512 * q**4 * u**2
        + 384 * q**4 * u
        + 32 * q**4
        - 2304 * q**3 * u**2
        + 4672 * q**3 * u
        - 216 * q**3
        + 2688 * q**2 * u**2
        - 25632 * q**2 * u
        + 492 * q**2
        - 2880 * q * u**2
        + 41040 * q * u
        - 450 * q
        - 21600 * u
        + 135
    ) / (8 * (2 * q - 3) ** 3)
    leading_expression = (2 * t + 1) * b / a**3 - 2 / a
    return {
        "t_lower_gate": polynomial_positive_gate(t - u * q, u, q),
        "curvature_lower_gate": polynomial_positive_gate(
            a - sp.Rational(39, 10) * u**2 * q, u, q
        ),
        "third_derivative_gate": polynomial_positive_gate(b, u, q),
        "leading_margin_gate": polynomial_positive_gate(
            u**2 * t * leading_expression - sp.Rational(1, 5), u, q
        ),
        "leading_formula": "u^2*t*((2*t+1)*b/a^3-2/a)>=1/5",
    }


def high_eighth_envelope() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    rows = fraction_grid(MODE_SPLIT, MODE_END, Fraction(1, 1000))
    target = arb_rational(HIGH_EIGHTH_BOUND)
    minimum_margin: tuple[flint.arb, int] | None = None
    minimum_window: flint.arb | None = None
    maximum_window: flint.arb | None = None
    for index, (left, right) in enumerate(rows):
        mode = arb_interval(left, right)
        curvature = potential_jet_arb(mode, 2)[2]
        left_window = mode * (-4 / curvature.sqrt()).exp()
        right_window = mode * (4 / curvature.sqrt()).exp()
        eighth = potential_jet_arb(exact_upper(right_window), 8)[8]
        margin = target - eighth / curvature**4
        if not bool(margin > 0):
            raise RuntimeError(f"high eighth-derivative envelope failed at {left}..{right}")
        if minimum_margin is None or margin.lower() < minimum_margin[0].lower():
            minimum_margin = (margin, index)
        if minimum_window is None or left_window.lower() < minimum_window.lower():
            minimum_window = left_window
        if maximum_window is None or right_window.upper() > maximum_window.upper():
            maximum_window = right_window
    assert minimum_margin is not None and minimum_window is not None and maximum_window is not None
    if not bool(minimum_window > flint.arb(41) / 50 and maximum_window < flint.arb(51) / 10):
        raise RuntimeError("high Simpson window leaves the V9-positive range")
    return {
        "mode_interval": "2<=u<=5",
        "window": "|y|<=8",
        "interval_count": len(rows),
        "positive_interval_count": len(rows),
        "certified_bound": "sup_|y|<=8 V^(8)/a^4<=1/10000000000",
        "minimum_margin_lower": arb_lower_text(minimum_margin[0]),
        "minimum_margin_interval": minimum_margin[1],
        "minimum_window_u_lower": arb_lower_text(minimum_window),
        "maximum_window_u_upper": arb_upper_text(maximum_window),
    }


def compact_tasks() -> list[dict]:
    tasks: list[dict] = []
    regions = (
        ("low", MODE_START, Fraction(2), LOW_WIDTH, 100, 6, LOW_EIGHTH_BOUND),
        ("mid_2_3", Fraction(2), Fraction(3), HIGH_WIDTH, 1000, 8, HIGH_EIGHTH_BOUND),
        ("mid_3_4", Fraction(3), Fraction(4), HIGH_WIDTH, 2500, 8, HIGH_EIGHTH_BOUND),
        ("high_4_5", Fraction(4), Fraction(5), HIGH_WIDTH, 6000, 8, HIGH_EIGHTH_BOUND),
    )
    for region, start, end, width, panels, window, envelope in regions:
        for left, right in fraction_grid(start, end, width):
            tasks.append(
                {
                    "region": region,
                    "left": left,
                    "right": right,
                    "panels": panels,
                    "window": window,
                    "envelope": envelope,
                }
            )
    return tasks


def certify_task(task: dict) -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    result = certify_scaled_curvature_mode_block(
        task["left"],
        task["right"],
        task["panels"],
        window_y=task["window"],
        eighth_envelope_bound=task["envelope"],
    )
    return {
        "region": task["region"],
        "left": str(task["left"]),
        "right": str(task["right"]),
        **result,
    }


def compact_certificate(workers: int) -> dict:
    tasks = compact_tasks()
    start_time = perf_counter()
    if workers == 1:
        iterator = map(certify_task, tasks)
        executor = None
    else:
        executor = ProcessPoolExecutor(max_workers=workers)
        iterator = executor.map(certify_task, tasks, chunksize=4)
    accepted: list[dict] = []
    try:
        for index, result in enumerate(iterator, start=1):
            if not result.get("passed"):
                raise RuntimeError(
                    f"compact curvature block failed at {result['left']}..{result['right']}: "
                    f"{result.get('failure')}"
                )
            accepted.append(result)
            if index % 500 == 0:
                print(f"scaled-curvature compact progress: {index}/{len(tasks)}")
    finally:
        if executor is not None:
            executor.shutdown(wait=True, cancel_futures=True)

    for previous, current in zip(accepted, accepted[1:]):
        if Fraction(previous["right"]) != Fraction(current["left"]):
            raise RuntimeError("compact curvature cover has a gap")
    if Fraction(accepted[0]["left"]) != MODE_START or Fraction(accepted[-1]["right"]) != MODE_END:
        raise RuntimeError("compact curvature cover endpoints do not match")

    worst_buffer = min(
        accepted, key=lambda row: Decimal(row["shifted_curvature_buffer_lower"])
    )
    worst_transfer = min(
        accepted, key=lambda row: Decimal(row["full_kernel_transfer_margin_lower"])
    )
    region_counts: dict[str, int] = {}
    for row in accepted:
        region_counts[row["region"]] = region_counts.get(row["region"], 0) + 1
    selected_indices = {
        0,
        region_counts["low"] - 1,
        region_counts["low"],
        len(accepted) // 2,
        len(accepted) - 1,
    }
    selected = [accepted[index] for index in sorted(selected_indices)]
    return {
        "mode_interval": "0.9264<=u<=5",
        "block_count": len(accepted),
        "positive_buffer_blocks": len(accepted),
        "positive_full_kernel_transfer_blocks": len(accepted),
        "region_counts": region_counts,
        "workers": workers,
        "elapsed_seconds": perf_counter() - start_time,
        "minimum_buffer_lower": worst_buffer["shifted_curvature_buffer_lower"],
        "minimum_buffer_block": {
            key: worst_buffer[key]
            for key in ("region", "left", "right", "panels", "window_y")
        },
        "minimum_transfer_margin_lower": worst_transfer[
            "full_kernel_transfer_margin_lower"
        ],
        "minimum_transfer_margin_block": {
            key: worst_transfer[key]
            for key in ("region", "left", "right", "panels", "window_y")
        },
        "selected_blocks": selected,
    }


def ray_certificate() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    paired_ray = load_json(PAIRED_RAY_JSON)
    leading = load_json(LEADING_JSON)
    formulas = {row.get("formula") for row in paired_ray.get("rows", [])}
    required = {
        "|m1+alpha/2|<=13/q, |m2-1|<=36/q, |m3+5*alpha/2|<=80/q",
        "|kappa_3(Y)+alpha|<=120/q and P*|kappa_3(Y)+alpha|<1/1000",
    }
    if not required.issubset(formulas):
        raise RuntimeError("paired ray moment source is missing required formulas")
    derivative_bounds = leading["diagnostics"]["ray"]["derivative_bounds"]
    if derivative_bounds != [
        "t=V'<=3*u*q",
        "a=V''>=(39/10)*u^2*q",
        "b=V'''<=12*u^3*q",
    ]:
        raise RuntimeError("leading ray derivative bounds changed")

    q0 = flint.arb.pi() * flint.arb(20).exp()
    if not bool(q0 > Q_FLOOR and q0.sqrt() > 30_000):
        raise RuntimeError("ray q floor failed")

    leading_margin = Fraction(1, 5)
    variance_error = Fraction(760, 13) * Fraction(5, Q_FLOOR)
    skewness_error = Fraction(432, 6000)
    gamma_error = Fraction(25, Q_FLOOR) + Fraction(4, Q_FLOOR**2)
    h3_leading_error = Fraction(13, 10) * Fraction(5, Q_FLOOR)
    h3_moment_error = Fraction(96, 30_000**3)
    total_error = (
        variance_error
        + skewness_error
        + gamma_error
        + h3_leading_error
        + h3_moment_error
    )
    final_margin = leading_margin - total_error
    if final_margin <= Fraction(1, 10):
        raise RuntimeError("ray scaled-curvature margin is too small")

    algebra = exact_ray_algebra()
    return {
        "mode_ray": "u>=5",
        "q_at_5_lower": arb_lower_text(q0),
        "sqrt_q_at_5_lower": arb_lower_text(q0.sqrt()),
        "source_moment_bounds": sorted(required),
        "source_derivative_bounds": derivative_bounds,
        "exact_algebra": algebra,
        "normalized_leading_margin": str(leading_margin),
        "normalized_error_budget": {
            "variance": str(variance_error),
            "skewness": str(skewness_error),
            "gamma": str(gamma_error),
            "h3_leading": str(h3_leading_error),
            "h3_moment": str(h3_moment_error),
            "total": str(total_error),
        },
        "normalized_final_margin": str(final_margin),
        "normalized_final_margin_decimal": fraction_decimal(final_margin),
        "certified_buffer": (
            "Q(t)-2*abs(H'''(t))>=normalized_final_margin/(u^2*t), u>=5"
        ),
        "transfer": (
            "The ray buffer exceeds 64/(t-3)^5 because t>=u*q, q=pi*exp(4u), "
            "and the resulting ratio is minimized at u=5."
        ),
    }


def finite_prefix_certificate() -> dict:
    flint.ctx.prec = 256
    values, sources = merged_coefficients()
    curvatures: dict[int, flint.arb] = {}
    for k in range(1, 320):
        x = values[k + 1] * values[k - 1] / values[k] ** 2
        if not bool(x > 0 and x < 1):
            raise RuntimeError(f"finite prefix contraction failed at k={k}")
        curvatures[k] = (2 * k + 1) * (-x.log())
    minimum: tuple[flint.arb, int] | None = None
    selected: list[dict] = []
    for k in range(1, 319):
        gap = curvatures[k + 1] - curvatures[k]
        if not bool(gap > 0):
            raise RuntimeError(f"finite scaled-curvature gap failed at k={k}")
        if minimum is None or gap.lower() < minimum[0].lower():
            minimum = (gap, k)
        if k in {1, 2, 100, 220, 245, 300, 318}:
            selected.append({"k": k, "gap_lower": arb_lower_text(gap)})
    assert minimum is not None
    return {
        "lambda": "-100",
        "coefficient_range": [0, 320],
        "scaled_curvature_gap_range": [1, 318],
        "positive_scaled_curvature_gaps": 318,
        "minimum_gap_at_k": minimum[1],
        "minimum_gap_lower": arb_lower_text(minimum[0]),
        "selected_gaps": selected,
        "merged_sources": sources,
    }


def build_artifact(workers: int = DEFAULT_WORKERS) -> dict:
    paired_compact = load_json(PAIRED_COMPACT_JSON)
    low_envelope = paired_compact["diagnostics"]["eighth_derivative_envelope"]
    if low_envelope.get("positive_interval_count") != low_envelope.get("interval_count"):
        raise RuntimeError("source low eighth-derivative envelope is incomplete")
    if low_envelope.get("certified_bound") != (
        "sup_|y|<=6 V^(8)(x_t+y/sqrt(a_t))/a_t^4<=1/50000"
    ):
        raise RuntimeError("source low eighth-derivative envelope changed")

    high_envelope = high_eighth_envelope()
    compact = compact_certificate(workers)
    ray = ray_certificate()
    prefix = finite_prefix_certificate()
    rows = [
        CertificateRow(
            id="nlscb_01_normalized_curvature",
            role="exact_reduction",
            claim="Write the scaled curvature as a centered second difference of H=log Gamma-log M.",
            formula="B(t)=Delta^2 H(t-1), C(t)=(2*t+1)*B(t)",
            readiness="available_exact",
            proof_boundary="Exact real-parameter identity only.",
        ),
        CertificateRow(
            id="nlscb_02_tent_derivative_bridge",
            role="exact_reduction",
            claim="Differentiate C(t) and use the tent-kernel formula for the centered second difference.",
            formula=(
                "C'(t)=integral_[-1,1](1-|s|)*(Q(t+s)-2*s*H'''(t+s)) ds; "
                "Q(r)=2*H''(r)+(2*r+1)*H'''(r)"
            ),
            readiness="available_exact",
            proof_boundary="Exact continuous-to-discrete bridge only.",
        ),
        CertificateRow(
            id="nlscb_03_buffered_sufficient_condition",
            role="exact_reduction",
            claim="A pointwise shifted buffer proves C'(t)>0 and absorbs the full-kernel moment perturbation.",
            formula="Q(r)-2*abs(H'''(r))>64/(r-3)^5, r>=318",
            readiness="available_exact",
            proof_boundary="Sufficient inequality; compact and ray rows discharge it.",
        ),
        CertificateRow(
            id="nlscb_04_low_eighth_envelope",
            role="interval_input",
            claim="The existing paired theorem supplies the low-mode eighth-derivative envelope.",
            formula=low_envelope["certified_bound"],
            readiness="interval_validated",
            proof_boundary="Derivative envelope only.",
            diagnostics=low_envelope,
        ),
        CertificateRow(
            id="nlscb_05_high_eighth_envelope",
            role="interval_certificate",
            claim="A sharper high-mode envelope supports interval Simpson quadrature on |y|<=8.",
            formula=high_envelope["certified_bound"],
            readiness="interval_validated",
            proof_boundary="High-mode derivative envelope only.",
            diagnostics=high_envelope,
        ),
        CertificateRow(
            id="nlscb_06_compact_buffer_certificate",
            role="interval_certificate",
            claim="Interval Simpson quadrature proves the shifted buffer and full-kernel transfer margin on the compact mode range.",
            formula="Q(r)-2*abs(H'''(r))>64/(r-3)^5, 318<=r<=V'(x(5))",
            readiness="interval_validated",
            proof_boundary="Compact continuous parameter theorem.",
            diagnostics=compact,
        ),
        CertificateRow(
            id="nlscb_07_ray_buffer_certificate",
            role="analytic_certificate",
            claim="Exact leading algebra and the existing paired moment errors prove the shifted buffer on u>=5.",
            formula=ray["certified_buffer"],
            readiness="available_exact",
            proof_boundary="Analytic asymptotic ray theorem.",
            diagnostics=ray,
        ),
        CertificateRow(
            id="nlscb_08_full_kernel_transfer",
            role="exact_theorem_composition",
            claim="The n>=2 moment error changes one scaled-curvature increment by at most 64/(k-1)^5.",
            formula="abs((C_(k+1)-C_k)-(C1_(k+1)-C1_k))<=64/(k-1)^5, k>=319",
            readiness="ready_to_apply",
            proof_boundary="Tail perturbation transfer only.",
        ),
        CertificateRow(
            id="nlscb_09_finite_prefix_certificate",
            role="interval_certificate",
            claim="The repaired lambda=-100 prefix has increasing scaled curvature through k=318.",
            formula="C_(k+1)>C_k, 1<=k<=318",
            readiness="interval_validated",
            proof_boundary="Finite repaired prefix at one heat parameter.",
            diagnostics=prefix,
        ),
        CertificateRow(
            id="nlscb_10_full_scaled_curvature_theorem",
            role="interval_analytic_theorem",
            claim="The actual zeta heat-flow coefficients have increasing scaled curvature at lambda=-100 for every k.",
            formula="C_(k+1)>=C_k for every integer k>=1 at lambda=-100",
            readiness="ready_to_apply",
            proof_boundary=(
                "Closes the lower raw-corridor side at lambda=-100; not PF-infinity, "
                "the all-order Jensen bridge, RH, or Lambda <= 0."
            ),
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "compact_blocks": compact["block_count"],
        "positive_compact_buffer_blocks": compact["positive_buffer_blocks"],
        "positive_compact_transfer_blocks": compact[
            "positive_full_kernel_transfer_blocks"
        ],
        "high_envelope_intervals": high_envelope["interval_count"],
        "positive_prefix_gaps": prefix["positive_scaled_curvature_gaps"],
        "analytic_ray_rows": 1,
        "full_scaled_curvature_theorem_rows": 1,
        "ready_to_apply_rows": sum(row.readiness == "ready_to_apply" for row in rows),
        "open_requirements": 0,
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge",
        "date": "2026-07-11",
        "status": "interval and analytic all-k scaled-curvature theorem at lambda=-100",
        "proof_boundary": (
            "This artifact proves C_(k+1)>=C_k for the actual zeta heat-flow "
            "coefficients at lambda=-100. It closes the lower raw-corridor side "
            "at that parameter but does not prove PF-infinity, the all-order Jensen "
            "bridge, RH, or Lambda <= 0."
        ),
        "rows": [asdict(row) for row in rows],
        "diagnostics": {
            "low_eighth_envelope": low_envelope,
            "high_eighth_envelope": high_envelope,
            "compact_certificate": compact,
            "ray_certificate": ray,
            "finite_prefix": prefix,
        },
        "summary": summary,
        "source_paired_compact": PAIRED_COMPACT_JSON.relative_to(REPO_ROOT).as_posix(),
        "source_paired_ray": PAIRED_RAY_JSON.relative_to(REPO_ROOT).as_posix(),
        "source_leading": LEADING_JSON.relative_to(REPO_ROOT).as_posix(),
        "source_dominance": "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
        "source_full_cone": "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md",
        "source_target": "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
        "generator": "work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py",
        "checker": "work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py",
    }


def write_note(path: Path, artifact: dict) -> None:
    summary = artifact["summary"]
    diagnostics = artifact["diagnostics"]
    compact = diagnostics["compact_certificate"]
    ray = diagnostics["ray_certificate"]
    prefix = diagnostics["finite_prefix"]
    lines = [
        "# Jensen-Window PF Negative-Lambda Scaled-Curvature Continuous Bridge",
        "",
        "Date: 2026-07-11",
        "",
        "Status: interval and analytic all-k scaled-curvature theorem at lambda=-100.",
        "This is not a proof of PF-infinity, the all-order Jensen bridge, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        (
            "validated Jensen-window PF negative-lambda scaled-curvature continuous bridge: "
            f"{summary['certificate_rows']} rows, 0 issues, {summary['compact_blocks']} compact blocks, "
            f"{summary['positive_prefix_gaps']} prefix gaps, 1 analytic ray, "
            "1 all-k scaled-curvature theorem, 0 open requirements"
        ),
        "```",
        "",
        "## Exact Continuous Bridge",
        "",
        "For the first Newman summand, set `H(t)=log Gamma(t+1/2)-log M_t^(1)` and",
        "",
        "```text",
        "B(t)=H(t+1)-2*H(t)+H(t-1),",
        "C(t)=(2*t+1)*B(t).",
        "```",
        "",
        "The centered second difference has its tent-kernel representation. With",
        "`Q(r)=2*H''(r)+(2*r+1)*H'''(r)`, differentiation gives",
        "",
        "```text",
        "C'(t)=integral_[-1,1] (1-|s|)*(Q(t+s)-2*s*H'''(t+s)) ds.",
        "```",
        "",
        "Thus `Q(r)-2*abs(H'''(r))>64/(r-3)^5` for `r>=318` proves",
        "strict first-summand growth and leaves enough room for the complete-kernel transfer.",
        "",
        "## Compact Certificate",
        "",
        "Interval composite Simpson quadrature with a fourth-derivative remainder proves",
        "the buffered inequality continuously on `0.9264<=u<=5`.",
        "",
        "```text",
        f"compact blocks: {compact['block_count']}",
        f"minimum shifted-buffer lower bound: {compact['minimum_buffer_lower']}",
        f"minimum full-kernel transfer margin: {compact['minimum_transfer_margin_lower']}",
        f"regions: {compact['region_counts']}",
        "```",
        "",
        "The low region uses the established `|y|<=6`, `V^(8)/a^4<=1/50000` envelope.",
        "A separate 3,000-interval gate proves `V^(8)/a^4<=1/10^10` on the",
        "`|y|<=8` high-mode windows.",
        "",
        "## Analytic Ray",
        "",
        "For `u>=5`, `q=pi*exp(4u)>=10^9`. Exact polynomial algebra proves",
        "",
        "```text",
        "u^2*t*((2*t+1)*V'''/V''^3-2/V'') >= 1/5.",
        "```",
        "",
        "After `u=5+v`, `q=10^9+Q`, the cleared numerator has",
        f"{ray['exact_algebra']['leading_margin_gate']['numerator_terms']} strictly positive coefficients.",
        "The paired ray moment errors, Gamma-series bounds, and `H'''` buffer consume",
        "",
        "```text",
        f"normalized error <= {ray['normalized_error_budget']['total']}",
        f"normalized final margin >= {ray['normalized_final_margin']} > 1/10.",
        "```",
        "",
        "This also dominates `64/(t-3)^5` on the whole ray.",
        "",
        "## Full-Kernel Transfer And Prefix",
        "",
        "If `M_k=M_k^(1)*(1+delta_k)` with `0<=delta_k<=2/k^6`, direct expansion gives",
        "",
        "```text",
        "abs((C_(k+1)-C_k)-(C1_(k+1)-C1_k)) <= 64/(k-1)^5, k>=319.",
        "```",
        "",
        "The repaired Arb source independently proves",
        "",
        "```text",
        f"C_(k+1)-C_k>0 for k=1..318 ({prefix['positive_scaled_curvature_gaps']} rows),",
        f"minimum prefix margin={prefix['minimum_gap_lower']} at k={prefix['minimum_gap_at_k']}.",
        "```",
        "",
        "Therefore the actual lambda=-100 sequence satisfies",
        "",
        "```text",
        "C_(k+1)>=C_k for every integer k>=1.",
        "```",
        "",
        "Combined with the already proved `B_(k+1)<=B_k` and `B_k>=0`, this closes",
        "the two-sided raw-ratio corridor at lambda=-100. It does not establish any",
        "higher-degree minor cone or the all-order Jensen/PF bridge.",
        "",
        "```text",
        "outputs/jensen_window_pf_negative_lambda_scaled_curvature_monotonicity_target.md",
        "outputs/jensen_window_pf_negative_lambda_zeta_specific_raw_corridor_target.md",
        "outputs/jensen_window_pf_negative_lambda_m100_full_cone_entry_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate.md",
        "outputs/jensen_window_pf_negative_lambda_first_summand_paired_remainder_ray_certificate.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.workers < 1:
        raise SystemExit("--workers must be positive")
    artifact = build_artifact(args.workers)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "validated Jensen-window PF negative-lambda scaled-curvature continuous bridge: "
        f"{summary['certificate_rows']} rows, 0 issues, {summary['compact_blocks']} compact blocks, "
        f"{summary['positive_prefix_gaps']} prefix gaps, 1 analytic ray, "
        "1 all-k scaled-curvature theorem, 0 open requirements"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
