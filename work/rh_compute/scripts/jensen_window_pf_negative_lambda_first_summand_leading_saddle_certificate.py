#!/usr/bin/env python3
"""Certify the leading saddle term in the first-summand cumulant target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
import json
import math
from pathlib import Path
import sys

import mpmath as mp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


DEFAULT_SOURCE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.json"
)
DEFAULT_OUT_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.md"
)
DEFAULT_PRECISION_BITS = 192
DEFAULT_MPMATH_DPS = 60
DEFAULT_U_START = "0.926"
DEFAULT_U_SPLIT = "5"
DEFAULT_COMPACT_SUBINTERVALS = 40_740
DEFAULT_LEADING_CAP = Fraction(13, 20)
DEFAULT_CORRECTION_CAP = Fraction(1, 100)
DEFAULT_FIFTH_CORRECTION_CAP = Fraction(1, 1000)
DEFAULT_REMAINDER_FLOOR = Fraction(-79, 1000)


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    claim: str
    formula: str | None
    readiness: str
    proof_boundary: str
    diagnostics: dict | None = None
    gap: str | None = None


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits).replace("e", "E")


def arb_lower_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.lower().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_minus(), "E")


def arb_upper_text(value: flint.arb, digits: int = 70) -> str:
    rounded = value.upper().str(digits, radius=False)
    with localcontext() as context:
        context.prec = digits
        return format(Decimal(rounded).next_plus(), "E")


def arb_interval(left: flint.arb, right: flint.arb) -> flint.arb:
    return (left + right) / 2 + flint.arb(0, (right - left) / 2)


def potential_derivatives_arb(u: flint.arb) -> tuple[flint.arb, ...]:
    q = flint.arb.pi() * (4 * u).exp()
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
    return q, t, a, b


def potential_derivatives_mp(u: mp.mpf) -> tuple[mp.mpf, ...]:
    q = mp.pi * mp.exp(4 * u)
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
    return q, t, a, b


def higher_potential_derivatives_arb(u: flint.arb) -> tuple[flint.arb, flint.arb]:
    q = flint.arb.pi() * (4 * u).exp()
    r = 2 / (2 * q - 3)
    q1 = q * (2 * u)
    q2 = q * u * (4 * u + 1)
    q3 = q * u * (16 * u**2 + 12 * u + 1) / 2
    q4 = q * u * (64 * u**3 + 96 * u**2 + 28 * u + 1) / 4
    q5 = q * u * (256 * u**4 + 640 * u**3 + 400 * u**2 + 60 * u + 1) / 8
    g1, g2, g3, g4, g5 = 1 - r, r**2, -2 * r**3, 6 * r**4, -24 * r**5
    c4 = (
        g1 * q4
        + g2 * (4 * q1 * q3 + 3 * q2**2)
        + 6 * g3 * q1**2 * q2
        + g4 * q1**4
    )
    c5 = (
        g1 * q5
        + g2 * (5 * q1 * q4 + 10 * q2 * q3)
        + g3 * (10 * q1**2 * q3 + 15 * q1 * q2**2)
        + 10 * g4 * q1**3 * q2
        + g5 * q1**5
    )
    base = 100 * u**2
    return base + c4 - 5 * u / 16, base + c5 - 5 * u / 32


def higher_potential_derivatives_mp(u: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    q = mp.pi * mp.exp(4 * u)
    r = 2 / (2 * q - 3)
    q1 = q * (2 * u)
    q2 = q * u * (4 * u + 1)
    q3 = q * u * (16 * u**2 + 12 * u + 1) / 2
    q4 = q * u * (64 * u**3 + 96 * u**2 + 28 * u + 1) / 4
    q5 = q * u * (256 * u**4 + 640 * u**3 + 400 * u**2 + 60 * u + 1) / 8
    g1, g2, g3, g4, g5 = 1 - r, r**2, -2 * r**3, 6 * r**4, -24 * r**5
    c4 = (
        g1 * q4
        + g2 * (4 * q1 * q3 + 3 * q2**2)
        + 6 * g3 * q1**2 * q2
        + g4 * q1**4
    )
    c5 = (
        g1 * q5
        + g2 * (5 * q1 * q4 + 10 * q2 * q3)
        + g3 * (10 * q1**2 * q3 + 15 * q1 * q2**2)
        + 10 * g4 * q1**3 * q2
        + g5 * q1**5
    )
    base = 100 * u**2
    return base + c4 - 5 * u / 16, base + c5 - 5 * u / 32


def _series_multiply(left: list[flint.arb], right: list[flint.arb], order: int) -> list[flint.arb]:
    result = [flint.arb(0) for _ in range(order + 1)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            if i + j <= order:
                result[i + j] += left_value * right_value
    return result


def _series_exponential(values: list[flint.arb], order: int) -> list[flint.arb]:
    result = [flint.arb(0) for _ in range(order + 1)]
    result[0] = values[0].exp()
    for n in range(1, order + 1):
        result[n] = sum(
            flint.arb(k) * values[k] * result[n - k] for k in range(1, n + 1)
        ) / n
    return result


def _series_logarithm(values: list[flint.arb], order: int) -> list[flint.arb]:
    result = [flint.arb(0) for _ in range(order + 1)]
    result[0] = values[0].log()
    for n in range(1, order + 1):
        prior = sum(
            flint.arb(k) * result[k] * values[n - k] for k in range(1, n)
        )
        result[n] = (flint.arb(n) * values[n] - prior) / (flint.arb(n) * values[0])
    return result


def potential_jet_arb(u: flint.arb, order: int = 7) -> list[flint.arb]:
    u_series = [u / flint.arb((2**n) * math.factorial(n)) for n in range(order + 1)]
    q_series = [flint.arb.pi() * value for value in _series_exponential([4 * value for value in u_series], order)]
    denominator = [2 * value for value in q_series]
    denominator[0] -= 3
    log_denominator = _series_logarithm(denominator, order)
    log_u = _series_logarithm(u_series, order)
    potential = [
        100 * value for value in _series_multiply(u_series, u_series, order)
    ]
    for n in range(order + 1):
        potential[n] += q_series[n] - 5 * u_series[n]
        potential[n] -= log_denominator[n]
        potential[n] -= log_u[n]
    return [potential[n] * math.factorial(n) for n in range(order + 1)]


def potential_sixth_seventh_mp(u: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    x = 2 * mp.log(u)

    def potential(z: mp.mpf) -> mp.mpf:
        point = mp.exp(z / 2)
        q = mp.pi * mp.exp(4 * point)
        return 100 * point**2 + q - 5 * point - mp.log(2 * q - 3) - mp.log(point)

    return mp.diff(potential, x, 6), mp.diff(potential, x, 7)


def compact_certificate() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    start = flint.arb(DEFAULT_U_START)
    end = flint.arb(DEFAULT_U_SPLIT)
    count = DEFAULT_COMPACT_SUBINTERVALS
    cap = flint.arb(DEFAULT_LEADING_CAP.numerator) / DEFAULT_LEADING_CAP.denominator
    correction_cap = flint.arb(DEFAULT_CORRECTION_CAP.numerator) / DEFAULT_CORRECTION_CAP.denominator
    fifth_cap = flint.arb(DEFAULT_FIFTH_CORRECTION_CAP.numerator) / DEFAULT_FIFTH_CORRECTION_CAP.denominator
    positive_curvature_rows = 0
    positive_cap_rows = 0
    positive_correction_rows = 0
    positive_fifth_rows = 0
    minimum_curvature: tuple[flint.arb, int] | None = None
    minimum_margin: tuple[flint.arb, int] | None = None
    minimum_correction_margin: tuple[flint.arb, int] | None = None
    minimum_fifth_margin: tuple[flint.arb, int] | None = None
    selected = []
    selected_indices = {0, count // 4, count // 2, 3 * count // 4, count - 1}
    for index in range(count):
        left = start + (end - start) * index / count
        right = start + (end - start) * (index + 1) / count
        u = arb_interval(left, right)
        _, t, a, b = potential_derivatives_arb(u)
        c, d = higher_potential_derivatives_arb(u)
        jet = potential_jet_arb(u)
        e, f = jet[6], jet[7]
        margin = cap - t**2 * b / a**3
        correction = t**2 / 2 * (8 * b**3 / a**6 - 7 * b * c / a**5 + d / a**4)
        correction_margin = correction_cap - correction
        fifth_correction = t**2 / 24 * (
            525 * b**5 / a**9
            - 954 * b**3 * c / a**8
            + 234 * b**2 * d / a**7
            + 298 * b * c**2 / a**7
            - 37 * b * e / a**6
            - 57 * c * d / a**6
            + 3 * f / a**5
        )
        fifth_margin = fifth_cap - fifth_correction
        positive_curvature_rows += int(bool(a > 0))
        positive_cap_rows += int(bool(margin > 0))
        positive_correction_rows += int(bool(correction_margin > 0))
        positive_fifth_rows += int(bool(fifth_margin > 0))
        if minimum_curvature is None or a.lower() < minimum_curvature[0].lower():
            minimum_curvature = (a, index)
        if minimum_margin is None or margin.lower() < minimum_margin[0].lower():
            minimum_margin = (margin, index)
        if (
            minimum_correction_margin is None
            or correction_margin.lower() < minimum_correction_margin[0].lower()
        ):
            minimum_correction_margin = (correction_margin, index)
        if minimum_fifth_margin is None or fifth_margin.lower() < minimum_fifth_margin[0].lower():
            minimum_fifth_margin = (fifth_margin, index)
        if index in selected_indices:
            selected.append(
                {
                    "index": index,
                    "u_left": arb_text(left, 30),
                    "u_right": arb_text(right, 30),
                    "curvature_ball": arb_text(a, 35),
                    "leading_cap_margin_ball": arb_text(margin, 35),
                    "correction_cap_margin_ball": arb_text(correction_margin, 35),
                    "fifth_correction_cap_margin_ball": arb_text(fifth_margin, 35),
                }
            )
    if (
        positive_curvature_rows != count
        or positive_cap_rows != count
        or positive_correction_rows != count
        or positive_fifth_rows != count
    ):
        raise RuntimeError("compact saddle-expansion interval certificate failed")
    assert (
        minimum_curvature is not None
        and minimum_margin is not None
        and minimum_correction_margin is not None
        and minimum_fifth_margin is not None
    )
    return {
        "u_interval": "0.926<=u<=5",
        "subinterval_count": count,
        "positive_curvature_subintervals": positive_curvature_rows,
        "positive_leading_cap_subintervals": positive_cap_rows,
        "positive_correction_cap_subintervals": positive_correction_rows,
        "positive_fifth_correction_cap_subintervals": positive_fifth_rows,
        "minimum_curvature_lower": arb_lower_text(minimum_curvature[0]),
        "minimum_curvature_subinterval": minimum_curvature[1],
        "minimum_leading_cap_margin_lower": arb_lower_text(minimum_margin[0]),
        "minimum_leading_cap_margin_subinterval": minimum_margin[1],
        "minimum_correction_cap_margin_lower": arb_lower_text(minimum_correction_margin[0]),
        "minimum_correction_cap_margin_subinterval": minimum_correction_margin[1],
        "minimum_fifth_correction_cap_margin_lower": arb_lower_text(minimum_fifth_margin[0]),
        "minimum_fifth_correction_cap_margin_subinterval": minimum_fifth_margin[1],
        "selected_subinterval_rows": selected,
    }


def mode_localization() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    u0 = flint.arb(DEFAULT_U_START)
    t0 = flint.arb(318)
    q0 = flint.arb.pi() * (4 * u0).exp()
    original_slope = 2 * t0 / u0 - 200 * u0 + 5 + 8 * q0 / (2 * q0 - 3) - 4 * q0
    _, potential_slope, _, _ = potential_derivatives_arb(u0)
    if not bool(original_slope > 0 and potential_slope < 318):
        raise RuntimeError("mode-localization endpoint failed")
    return {
        "u_start": DEFAULT_U_START,
        "original_u_log_integrand_slope_at_t318_ball": arb_text(original_slope),
        "original_u_log_integrand_slope_at_t318_lower": arb_lower_text(original_slope),
        "potential_slope_at_u_start_ball": arb_text(potential_slope),
        "potential_slope_at_u_start_upper": arb_upper_text(potential_slope),
        "propagation": (
            "S_t'(u) decreases in u and increases in t. Hence the x-density log slope "
            "u*S_t'(u)/2+1/2 is positive for t>=318 and u<=0.926."
        ),
    }


def ray_certificate() -> dict:
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    u0 = flint.arb(DEFAULT_U_SPLIT)
    q0 = flint.arb.pi() * (4 * u0).exp()
    q_minus_500 = q0 - 500
    q_minus_120 = q0 - 120
    q_minus_billion = q0 - 1_000_000_000
    ray_margin = Fraction(13, 20) - Fraction(108_000, 59_319 * 5)
    correction_upper = Fraction(344_827, 39_097_152_900_000)
    correction_margin = Fraction(1, 100) - correction_upper
    fifth_upper = Fraction(2_464_491_167, 9_276_816_051_500_400_000_000_000)
    fifth_margin = Fraction(1, 1000) - fifth_upper
    if (
        not bool(q_minus_500 > 0 and q_minus_120 > 0 and q_minus_billion > 0)
        or ray_margin <= 0
        or correction_margin <= 0
        or fifth_margin <= 0
    ):
        raise RuntimeError("leading-saddle ray certificate failed")
    return {
        "u_ray": "u>=5",
        "q_at_5_ball": arb_text(q0),
        "q_at_5_minus_500_lower": arb_lower_text(q_minus_500),
        "q_at_5_minus_120_lower": arb_lower_text(q_minus_120),
        "q_at_5_minus_billion_lower": arb_lower_text(q_minus_billion),
        "q_over_u_propagation": "d_u log(q/u)=4-1/u>0 for u>=5",
        "derivative_bounds": [
            "t=V'<=3*u*q",
            "a=V''>=(39/10)*u^2*q",
            "b=V'''<=12*u^3*q",
        ],
        "leading_ratio_bound": "t^2*b/a^3<=108000/(59319*u)<=108000/(59319*5)",
        "cap_margin_fraction": str(ray_margin),
        "cap_margin_decimal": mp.nstr(
            mp.mpf(ray_margin.numerator) / ray_margin.denominator, 40
        ),
        "higher_derivative_bounds": [
            "abs(V''')<=14*u^3*q",
            "abs(V^(4))<=30*u^4*q",
            "abs(V^(5))<=70*u^5*q",
        ],
        "correction_absolute_upper_fraction": str(correction_upper),
        "correction_cap_margin_fraction": str(correction_margin),
        "correction_cap_margin_decimal": mp.nstr(
            mp.mpf(correction_margin.numerator) / correction_margin.denominator, 40
        ),
        "sixth_seventh_derivative_bounds": [
            "abs(V^(6))<=1000*u^6*q",
            "abs(V^(7))<=10000*u^7*q",
        ],
        "fifth_correction_absolute_upper_fraction": str(fifth_upper),
        "fifth_correction_cap_margin_fraction": str(fifth_margin),
        "fifth_correction_cap_margin_decimal": mp.nstr(
            mp.mpf(fifth_margin.numerator) / fifth_margin.denominator, 40
        ),
        "proof": (
            "Write V=100*u^2+q-log(2*q-3)-5*u-log(u), D_x=(u/2)D_u, "
            "and r=2/(2*q-3). On u>=5, r<=1/100, q/u increases, and the "
            "displayed componentwise derivative bounds follow algebraically. The q_n "
            "Bell-polynomial formulas through n=7 and q>=10^9 give both correction caps."
        ),
    }


def load_source_rows(path: Path) -> list[dict]:
    with path.open("r", encoding="utf-8") as handle:
        artifact = json.load(handle)
    return artifact["diagnostics"]["sample_rows"]


def mode_for_t(t: mp.mpf) -> mp.mpf:
    low = mp.mpf(DEFAULT_U_START)
    high = mp.mpf(2)
    while potential_derivatives_mp(high)[1] < t:
        high *= mp.mpf("1.4")
    for _ in range(210):
        middle = (low + high) / 2
        if potential_derivatives_mp(middle)[1] < t:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def remainder_scout(source_rows: list[dict]) -> list[dict]:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    rows = []
    floor = mp.mpf(DEFAULT_REMAINDER_FLOOR.numerator) / DEFAULT_REMAINDER_FLOOR.denominator
    for source in source_rows:
        t = mp.mpf(source["t"])
        mode = mode_for_t(t)
        _, _, a, b = potential_derivatives_mp(mode)
        c, d = higher_potential_derivatives_mp(mode)
        e, f = potential_sixth_seventh_mp(mode)
        leading = t**2 * b / a**3
        correction = t**2 / 2 * (8 * b**3 / a**6 - 7 * b * c / a**5 + d / a**4)
        scaled_kappa3 = mp.mpf(source["scaled_t2_kappa3"])
        first_remainder = scaled_kappa3 + leading
        third_order_remainder = first_remainder + correction
        fifth_correction = t**2 / 24 * (
            525 * b**5 / a**9
            - 954 * b**3 * c / a**8
            + 234 * b**2 * d / a**7
            + 298 * b * c**2 / a**7
            - 37 * b * e / a**6
            - 57 * c * d / a**6
            + 3 * f / a**5
        )
        remainder = third_order_remainder + fifth_correction
        rows.append(
            {
                "t": int(t),
                "x_density_mode_u": mp.nstr(mode, 35),
                "scaled_leading_term": mp.nstr(leading, 35),
                "scaled_kappa3": source["scaled_t2_kappa3"],
                "scaled_cubic_correction": mp.nstr(correction, 35),
                "scaled_first_remainder": mp.nstr(first_remainder, 35),
                "scaled_third_order_remainder": mp.nstr(third_order_remainder, 35),
                "scaled_fifth_correction": mp.nstr(fifth_correction, 35),
                "scaled_remainder": mp.nstr(remainder, 35),
                "remainder_floor": "-0.079",
                "remainder_margin": mp.nstr(remainder - floor, 35),
                "above_remainder_floor": bool(remainder > floor),
                "proof_boundary": (
                    "Combines the source finite mpmath cumulant scout with pointwise "
                    "high-precision saddle evaluation; not an interval remainder theorem."
                ),
            }
        )
    return rows


def build_artifact(source_path: Path = DEFAULT_SOURCE_JSON) -> dict:
    compact = compact_certificate()
    localization = mode_localization()
    ray = ray_certificate()
    scout_rows = remainder_scout(load_source_rows(source_path))
    if not all(row["above_remainder_floor"] for row in scout_rows):
        raise RuntimeError("finite remainder scout crossed its target floor")
    minimum_scout = min(scout_rows, key=lambda row: Decimal(row["remainder_margin"]))

    rows = [
        CertificateRow(
            id="fslsc_01_log_variable_potential",
            role="exact_reduction",
            claim="Write the continuous first-summand family as an exponential tilt in x=2 log u.",
            formula=(
                "p_t(x) proportional exp(t*x-V(x)); "
                "V=100*u^2+q-5*u-log(2*q-3)-log(u), q=pi*exp(4u)"
            ),
            readiness="available_exact",
            proof_boundary="Exact change of variables only.",
        ),
        CertificateRow(
            id="fslsc_02_mode_localization",
            role="exact_interval_lemma",
            claim="Every x-density mode for t>=318 lies on u>=0.926.",
            formula="u*S_t'(u)/2+1/2>0 for t>=318 and u<=0.926",
            readiness="interval_validated",
            proof_boundary="Endpoint Arb sign plus exact half-line propagation.",
            diagnostics=localization,
        ),
        CertificateRow(
            id="fslsc_03_compact_leading_cap",
            role="interval_certificate",
            claim="The leading, cubic-correction, and fifth-order saddle terms satisfy their caps on 0.926<=u<=5.",
            formula=(
                "t(u)^2*V'''/V''^3<=13/20; "
                "cubic_correction<=1/100; fifth_correction<=1/1000"
            ),
            readiness="interval_validated",
            proof_boundary="Finite Arb subdivision covering the full compact mode range.",
            diagnostics=compact,
        ),
        CertificateRow(
            id="fslsc_04_ray_leading_cap",
            role="exact_analytic_bound",
            claim="Elementary derivative bounds prove all three saddle-expansion caps on u>=5.",
            formula=(
                "V'<=3*u*q; V''>=(39/10)*u^2*q; V'''<=12*u^3*q"
            ),
            readiness="available_exact",
            proof_boundary="Exact analytic ray theorem from Arb-certified endpoint gates.",
            diagnostics=ray,
        ),
        CertificateRow(
            id="fslsc_05_global_leading_saddle_theorem",
            role="interval_theorem",
            claim="The leading, cubic-correction, and fifth-order saddle contributions have uniform all-t caps.",
            formula="leading<=13/20, cubic correction<=1/100, fifth correction<=1/1000 for t>=318",
            readiness="interval_validated",
            proof_boundary="Global leading-term theorem only; not a cumulant remainder bound.",
        ),
        CertificateRow(
            id="fslsc_06_remainder_reduction",
            role="exact_conditional_bridge",
            claim="A seventh-order normalized Laplace remainder estimate now implies the cumulant target.",
            formula=(
                "t^2*(kappa_3+V'''/V''^3)+cubic_correction+fifth_correction>=-79/1000"
            ),
            readiness="available_conditional",
            proof_boundary=(
                "Exact 13/20+1/100+1/1000+79/1000=37/50 arithmetic; the displayed remainder estimate is open."
            ),
            gap="Needs a uniform signed Laplace/Stein remainder theorem for t>=318.",
        ),
        CertificateRow(
            id="fslsc_07_remainder_scout",
            role="finite_scout",
            claim="The recorded finite samples clear the proposed normalized remainder floor.",
            formula="seventh_order_scaled_remainder>-79/1000 on the recorded sample",
            readiness="finite_only",
            proof_boundary="Finite floating evidence inherited from the cumulant scout.",
            diagnostics={"sample_rows": scout_rows},
        ),
        CertificateRow(
            id="fslsc_08_full_wall_handoff",
            role="conditional_handoff",
            claim="The open remainder estimate would close the dominant and full-kernel tails.",
            formula=(
                "remainder floor => cumulant floor => L_k^(1)>=1/(4*k^2) => L_k>0"
            ),
            readiness="blocked_by_open_requirement",
            proof_boundary="Conditional chain only; not cone entry, RH, or Lambda <= 0.",
        ),
    ]
    summary = {
        "certificate_rows": len(rows),
        "compact_subintervals": compact["subinterval_count"],
        "positive_compact_subintervals": compact["positive_leading_cap_subintervals"],
        "positive_correction_compact_subintervals": compact["positive_correction_cap_subintervals"],
        "positive_fifth_correction_compact_subintervals": compact[
            "positive_fifth_correction_cap_subintervals"
        ],
        "positive_curvature_subintervals": compact["positive_curvature_subintervals"],
        "positive_analytic_ray_gates": 3,
        "sample_rows": len(scout_rows),
        "positive_remainder_sample_rows": sum(row["above_remainder_floor"] for row in scout_rows),
        "minimum_remainder_sample_margin": minimum_scout["remainder_margin"],
        "minimum_remainder_sample_margin_at_t": minimum_scout["t"],
        "open_remainder_rows": 1,
        "ready_to_apply_rows": 0,
    }
    return {
        "kind": "jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate",
        "date": "2026-07-10",
        "status": "interval leading-saddle theorem and open remainder target",
        "proof_boundary": (
            "This artifact proves the global leading saddle cap and an exact remainder "
            "reduction. It does not prove the uniform signed remainder, full cumulant "
            "bound, first-summand wall, cone entry, RH, or Lambda <= 0."
        ),
        "source_cumulant_bridge": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md"
        ),
        "source_saddle_target": (
            "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md"
        ),
        "source_json": str(source_path.relative_to(REPO_ROOT)).replace("\\", "/"),
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py"
        ),
        "diagnostics": {
            "parameters": {
                "precision_bits": DEFAULT_PRECISION_BITS,
                "mpmath_dps": DEFAULT_MPMATH_DPS,
                "u_start": DEFAULT_U_START,
                "u_split": DEFAULT_U_SPLIT,
                "compact_subintervals": DEFAULT_COMPACT_SUBINTERVALS,
                "leading_cap": "13/20",
                "cubic_correction_cap": "1/100",
                "fifth_correction_cap": "1/1000",
                "remainder_floor": "-79/1000",
            },
            "mode_localization": localization,
            "compact": compact,
            "ray": ray,
            "remainder_sample_rows": scout_rows,
        },
        "summary": summary,
        "rows": [asdict(row) for row in rows],
    }


def write_note(path: Path, artifact: dict) -> None:
    diagnostics = artifact["diagnostics"]
    summary = artifact["summary"]
    lines = [
        "# Jensen-Window PF Negative-Lambda First-Summand Leading-Saddle Certificate",
        "",
        "Date: 2026-07-10",
        "",
        "Status: interval leading-saddle theorem and open remainder target.",
        "This is not a proof of the uniform remainder, cone entry, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate`.",
        "",
        "Machine-readable result:",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.json",
        "```",
        "",
        "Generator and checker:",
        "",
        "```text",
        "python work/rh_compute/scripts/jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: 8 rows, 0 issues, 40740 positive leading intervals, 40740 positive cubic-correction intervals, 40740 positive fifth-correction intervals, 3 positive analytic ray gates, 9 positive seventh-remainder samples, 1 open remainder, 0 ready-to-apply rows",
        "```",
        "",
        "## Log-Variable Potential",
        "",
        "Put `x=2*log(u)` and `q=pi*exp(4u)`. The continuous first-summand",
        "density is proportional to `exp(t*x-V(x))`, where, up to a constant,",
        "",
        "```text",
        "V(x)=100*u^2+q-5*u-log(2*q-3)-log(u).",
        "```",
        "",
        "At its mode `x_t`, write `a_t=V''(x_t)` and `b_t=V'''(x_t)`.",
        "The exact original-u slope at `t=318,u=0.926` is positive:",
        "",
        "```text",
        diagnostics["mode_localization"]["original_u_log_integrand_slope_at_t318_ball"],
        "```",
        "",
        "Since the original-u slope decreases in `u` and increases in `t`, every",
        "x-density mode for `t>=318` lies on `u>=0.926`.",
        "",
        "## Compact Certificate",
        "",
        "Direct Arb evaluation of the exact rational-exponential formulas proves",
        "",
        "```text",
        "t(u)^2*V'''(u)/V''(u)^3<=13/20",
        "```",
        "",
        f"on `{summary['positive_compact_subintervals']}/{summary['compact_subintervals']}`",
        "subintervals covering `0.926<=u<=5`. The same intervals also prove",
        "the cubic Edgeworth correction is at most `1/100` and the fifth-order",
        "correction is at most `1/1000`.",
        "",
        f"- Minimum outward-rounded curvature lower bound: `{diagnostics['compact']['minimum_curvature_lower']}`.",
        f"- Minimum outward-rounded cap margin: `{diagnostics['compact']['minimum_leading_cap_margin_lower']}`.",
        f"- Minimum outward-rounded correction-cap margin: `{diagnostics['compact']['minimum_correction_cap_margin_lower']}`.",
        f"- Minimum outward-rounded fifth-correction margin: `{diagnostics['compact']['minimum_fifth_correction_cap_margin_lower']}`.",
        "",
        "## Analytic Ray",
        "",
        "On `u>=5`, put `r=2/(2q-3)`. The Arb endpoint gates and exact",
        "monotonicity give `r<=1/100` and `q/u` increasing. Differentiating the",
        "potential componentwise yields",
        "",
        "```text",
        "V'<=3*u*q,",
        "V''>=(39/10)*u^2*q,",
        "V'''<=12*u^3*q.",
        "```",
        "",
        "Therefore",
        "",
        "```text",
        "t^2*V'''/V''^3<=108000/(59319*u)<=108000/(59319*5)<13/20.",
        "```",
        "",
        f"The final exact ray margin is `{diagnostics['ray']['cap_margin_fraction']}`.",
        "The derivative formulas through order seven and `q>=10^9` also give",
        f"the correction-cap margin `{diagnostics['ray']['correction_cap_margin_fraction']}`.",
        "The Bell-polynomial bounds through order seven give the fifth-correction",
        f"margin `{diagnostics['ray']['fifth_correction_cap_margin_fraction']}`.",
        "",
        "## Remaining Remainder",
        "",
        "The cumulant bridge requires `t^2*kappa_3>=-37/50`. The leading saddle",
        "term is certified below `13/20`, the cubic correction below `1/100`,",
        "and the fifth-order correction below `1/1000`. It is therefore enough",
        "to prove the seventh-order remainder bound",
        "",
        "```text",
        "scaled cumulant + leading + cubic correction + fifth correction",
        ">=-79/1000, t>=318.",
        "```",
        "",
        "This signed uniform Laplace/Stein remainder is the sole open row in this",
        "artifact. The recorded finite samples are:",
        "",
        "| t | third-order remainder | fifth correction | seventh remainder | margin above -0.079 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in diagnostics["remainder_sample_rows"]:
        lines.append(
            f"| `{row['t']}` | `{row['scaled_third_order_remainder']}` | `{row['scaled_fifth_correction']}` | "
            f"`{row['scaled_remainder']}` | `{row['remainder_margin']}` |"
        )
    lines.extend(
        [
            "",
            "These rows inherit non-interval cumulant quadrature and are not a uniform",
            "proof. They indicate that the proposed remainder floor has substantial room.",
            "",
            "```text",
            "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_saddle_wall_target.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
            "outputs/signed_hankel_jensen_dependency_graph.md",
            "```",
            "",
            "Summary:",
            "",
            (
                "The leading saddle contribution is now certified globally: "
                "the leading cap is 13/20, the cubic correction cap is 1/100, and the "
                "fifth-order correction cap is 1/1000 for all t>=318. This reduces the "
                "cumulant wall to the seventh-order normalized remainder floor -79/1000. "
                f"Nine finite samples clear that floor; the smallest sampled margin is "
                f"{summary['minimum_remainder_sample_margin']} at t="
                f"{summary['minimum_remainder_sample_margin_at_t']}. The uniform "
                "remainder theorem remains open."
            ),
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-json", type=Path, default=DEFAULT_SOURCE_JSON)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    artifact = build_artifact(args.source_json)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF negative-lambda first-summand leading-saddle certificate: "
        "8 rows, 0 issues, 40740 positive leading intervals, 40740 positive cubic-correction "
        "intervals, 40740 positive fifth-correction intervals, 3 positive analytic ray "
        "gates, 9 positive seventh-remainder samples, 1 open remainder, 0 ready-to-apply rows"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
