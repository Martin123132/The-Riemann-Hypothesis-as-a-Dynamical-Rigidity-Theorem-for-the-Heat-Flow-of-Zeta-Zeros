#!/usr/bin/env python3
"""Build the theta compact-transversality theorem/scout handoff."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_newman_theta_compact_transversality_scout.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_theta_compact_transversality_scout.md"
)

TIME_STEP = 0.005
X_STEP = 0.025
X_MIN = 0.25
X_MAX = 38.0
COARSE_ORDER = 600
FINE_ORDER = 900
CUTOFF = 2.0


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def exp_upper_rational(x: Fraction, degree: int) -> Fraction:
    """Upper-bound exp(x) by a Taylor sum plus a geometric tail."""
    if x < 0 or x >= degree + 2:
        raise ValueError("Taylor-tail ratio must lie in [0,1)")
    term = Fraction(1)
    partial = term
    for k in range(1, degree + 1):
        term *= x / k
        partial += term
    next_term = term * x / (degree + 1)
    ratio = x / (degree + 2)
    return partial + next_term / (1 - ratio)


def exp_lower_rational(x: Fraction, degree: int) -> Fraction:
    term = Fraction(1)
    partial = term
    for k in range(1, degree + 1):
        term *= x / k
        partial += term
    return partial


def build_exact() -> dict:
    exp_one_25_upper = exp_upper_rational(Fraction(1, 25), 8)
    if not exp_one_25_upper < Fraction(25, 24):
        raise RuntimeError("exp(1/25) rational upper bound failed")
    endpoint_exponent_upper = Fraction(22, 7) * Fraction(25, 24)
    if not endpoint_exponent_upper < Fraction(10, 3):
        raise RuntimeError("endpoint exponent upper bound failed")
    exp_ten_thirds_upper = exp_upper_rational(Fraction(10, 3), 48)
    if not exp_ten_thirds_upper < 29:
        raise RuntimeError("exp(10/3) rational upper bound failed")
    exp_three_lower = exp_lower_rational(Fraction(3), 12)
    if not exp_three_lower > 20:
        raise RuntimeError("exp(3) rational lower bound failed")
    if not Fraction(22, 7) ** 2 < 10:
        raise RuntimeError("pi-squared rational upper bound failed")

    m0_lower = Fraction(9, 2900)
    m2_upper = Fraction(40, 513)
    origin_correction_upper = Fraction(1, 32) * m2_upper
    origin_margin = m0_lower - origin_correction_upper
    if origin_margin <= 0:
        raise RuntimeError("near-origin positivity margin failed")

    return {
        "kernel": {
            "definition": (
                "phi_t(u)=exp(t*u^2)*Phi(u), "
                "H_t(x)=integral_0^infinity phi_t(u)*cos(xu)du"
            ),
            "summands": (
                "phi_n(u)=pi*n^2*exp(5u)"
                "*(2*pi*n^2*exp(4u)-3)*exp(-pi*n^2*exp(4u))"
            ),
            "positivity": (
                "phi_n(u)>0 for n>=1 and u>=0, hence phi_t(u)>0 "
                "for 0<=t<=1/5."
            ),
        },
        "zeroth_moment_floor": {
            "interval": "0<=u<=1/100",
            "rational_chain": (
                "pi>3, exp(4u)<=exp(1/25)<25/24, "
                "pi*exp(4u)<10/3, exp(10/3)<29"
            ),
            "first_summand_floor": "phi_1(u)>9/29",
            "conclusion": (
                "M_0(t)=integral phi_t(u)du>9/2900 "
                "for 0<=t<=1/5"
            ),
            "value": str(m0_lower),
        },
        "second_moment_ceiling": {
            "exponential_bound": (
                "exp(4u)>=1+4u+8u^2 and "
                "t*u^2+9u-pi*n^2*exp(4u)"
                "<=-pi*n^2-3*n^2*u"
            ),
            "integral_bound": (
                "M_2(t)<=(4*pi^2/27)"
                "*sum_(n>=1)exp(-pi*n^2)/n^2"
            ),
            "series_bound": (
                "sum exp(-pi*n^2)/n^2"
                "<sum exp(-3n)=1/(exp(3)-1)<1/19"
            ),
            "conclusion": "M_2(t)<40/513 for 0<=t<=1/5",
            "value": str(m2_upper),
        },
        "near_origin_theorem": {
            "cosine_bound": "cos(y)>=1-y^2/2 for every real y",
            "lower_bound": (
                "H_t(x)>=M_0(t)-x^2*M_2(t)/2"
                ">9/2900-5/2052"
            ),
            "interval": "0<=t<=1/5 and |x|<=1/4",
            "margin": str(origin_margin),
            "conclusion": (
                "H_t(x)>0 on the stated rectangle, so H_t and H_t' "
                "cannot vanish together there."
            ),
        },
        "rational_audit": {
            "exp_one_25_upper": str(exp_one_25_upper),
            "endpoint_exponent_upper": str(endpoint_exponent_upper),
            "exp_ten_thirds_upper": str(exp_ten_thirds_upper),
            "exp_three_lower": str(exp_three_lower),
            "origin_correction_upper": str(origin_correction_upper),
        },
    }


def component_ratio_grid(
    order: int,
    times: np.ndarray,
    xs: np.ndarray,
    *,
    cutoff: float = CUTOFF,
) -> dict:
    nodes, weights = roots_legendre(order)
    u = (cutoff / 2.0) * (nodes + 1.0)
    quadrature_weights = (cutoff / 2.0) * weights
    pi = math.pi
    phi_one = (
        2 * pi * pi * np.exp(9 * u) - 3 * pi * np.exp(5 * u)
    ) * np.exp(-pi * np.exp(4 * u))

    delta_0 = 1 / 2800
    delta_1 = 1 / 137200
    delta_2 = 1 / 3361400
    delta_3 = 1 / 54900000

    minimum = {
        "ratio": math.inf,
        "t": None,
        "x": None,
        "value_ratio": None,
        "derivative_ratio": None,
        "J_1": None,
        "J_1_prime": None,
        "B_0": None,
        "B_1": None,
    }
    failures = 0
    time_minima: list[dict] = []
    ratios = np.empty((len(times), len(xs)), dtype=np.float64)

    for time_index, t in enumerate(times):
        weighted_kernel = (
            np.exp(t * u * u) * phi_one * quadrature_weights
        )
        row_parts: list[np.ndarray] = []
        metadata_parts: list[tuple[np.ndarray, ...]] = []
        for start in range(0, len(xs), 500):
            x = xs[start : start + 500]
            xu = np.outer(x, u)
            cosine = np.cos(xu)
            sine = np.sin(xu)
            h = cosine @ weighted_kernel
            h_prime = -(sine @ (u * weighted_kernel))
            j_one = 16 * x**4 * h
            j_one_prime = 64 * x**3 * h + 16 * x**4 * h_prime

            p_poly = x**4 + (6 * t + 1) * x**2 + 24 * t**2
            r_poly = (6 * t + 1) * x**2 + 24 * t**2
            b_zero = (
                4 * t**2 * x**2 * delta_2
                + 4 * t * x * (x**2 + 4 * t) * delta_1
                + (x**4 + 2 * r_poly) * delta_0
            )
            b_one = (
                4 * t**2 * x**2 * delta_3
                + 4 * t * x * (x**2 + 2 * t) * delta_2
                + np.abs(p_poly - 4 * t * (3 * x**2 + 4 * t))
                * delta_1
                + (4 * x**3 + 4 * (6 * t + 1) * x) * delta_0
            )
            value_ratio = np.abs(j_one) / b_zero
            derivative_ratio = np.abs(j_one_prime) / b_one
            row_parts.append(np.maximum(value_ratio, derivative_ratio))
            metadata_parts.append(
                (
                    value_ratio,
                    derivative_ratio,
                    j_one,
                    j_one_prime,
                    b_zero,
                    b_one,
                )
            )

        row = np.concatenate(row_parts)
        ratios[time_index] = row
        failures += int(np.count_nonzero(row <= 1.0))
        row_index = int(np.argmin(row))
        time_minima.append(
            {
                "t": round(float(t), 12),
                "minimum_ratio": float(row[row_index]),
                "x": round(float(xs[row_index]), 12),
            }
        )
        if row[row_index] < minimum["ratio"]:
            chunk_index = row_index // 500
            within_chunk = row_index % 500
            metadata = metadata_parts[chunk_index]
            minimum = {
                "ratio": float(row[row_index]),
                "t": round(float(t), 12),
                "x": round(float(xs[row_index]), 12),
                "value_ratio": float(metadata[0][within_chunk]),
                "derivative_ratio": float(metadata[1][within_chunk]),
                "J_1": float(metadata[2][within_chunk]),
                "J_1_prime": float(metadata[3][within_chunk]),
                "B_0": float(metadata[4][within_chunk]),
                "B_1": float(metadata[5][within_chunk]),
            }

    return {
        "minimum": minimum,
        "failures": failures,
        "time_minima": time_minima,
        "ratios": ratios,
    }


def build_scout() -> dict:
    times = np.arange(0.0, 0.2 + TIME_STEP / 2, TIME_STEP)
    xs = np.arange(X_MIN, X_MAX + X_STEP / 2, X_STEP)
    coarse = component_ratio_grid(COARSE_ORDER, times, xs)
    fine = component_ratio_grid(FINE_ORDER, times, xs)
    maximum_ratio_difference = float(
        np.max(np.abs(fine["ratios"] - coarse["ratios"]))
    )

    outer_times = np.array([0.0, 0.05, 0.1, 0.15, 0.2])
    outer_xs = np.arange(38.0, 60.0 + 0.025, 0.05)
    outer = component_ratio_grid(FINE_ORDER, outer_times, outer_xs)
    failure_indices = np.argwhere(outer["ratios"] <= 1.0)
    first_failures: list[dict] = []
    for time_index, t in enumerate(outer_times):
        row_failures = np.flatnonzero(outer["ratios"][time_index] <= 1.0)
        first_failures.append(
            {
                "t": round(float(t), 12),
                "first_x": (
                    round(float(outer_xs[row_failures[0]]), 12)
                    if len(row_failures)
                    else None
                ),
                "failure_count": int(len(row_failures)),
            }
        )

    selected_time_minima = [
        fine["time_minima"][index]
        for index in (0, 10, 20, 30, 40)
    ]
    return {
        "method": {
            "description": (
                "Gauss-Legendre quadrature of the n=1 Phi summand on "
                "[0,2], followed by the exact identities "
                "J_1=16*x^4*H_1 and "
                "J_1'=64*x^3*H_1+16*x^4*H_1'."
            ),
            "coarse_order": COARSE_ORDER,
            "fine_order": FINE_ORDER,
            "cutoff": CUTOFF,
            "status": (
                "high-precision diagnostic only; no interval quadrature "
                "or between-grid derivative enclosure"
            ),
        },
        "compact_grid": {
            "time_min": 0.0,
            "time_max": 0.2,
            "time_step": TIME_STEP,
            "time_rows": int(len(times)),
            "x_min": X_MIN,
            "x_max": X_MAX,
            "x_step": X_STEP,
            "x_rows": int(len(xs)),
            "total_points": int(len(times) * len(xs)),
            "fine_minimum": fine["minimum"],
            "coarse_minimum": coarse["minimum"],
            "fine_failures": fine["failures"],
            "coarse_failures": coarse["failures"],
            "maximum_coarse_fine_ratio_difference": (
                maximum_ratio_difference
            ),
            "selected_time_minima": selected_time_minima,
        },
        "outer_scope_guard": {
            "time_rows": int(len(outer_times)),
            "x_min": 38.0,
            "x_max": 60.0,
            "x_step": 0.05,
            "total_points": int(len(outer_times) * len(outer_xs)),
            "failure_points": int(len(failure_indices)),
            "first_failures": first_failures,
            "interpretation": (
                "The raw first-block disjunction loses its margin near "
                "x=39.6 and cannot be promoted to an all-frequency theorem."
            ),
        },
    }


def build_rows(exact: dict, scout: dict) -> list[GateRow]:
    compact = scout["compact_grid"]
    outer = scout["outer_scope_guard"]
    return [
        GateRow(
            "ntcts_01_positive_phi_kernel",
            "exact_inequality",
            "ready_to_apply",
            "Every Xi kernel summand is positive on the half-line.",
            exact["kernel"]["positivity"],
            "Elementary from 2*pi*n^2*exp(4u)-3>0.",
        ),
        GateRow(
            "ntcts_02_zeroth_moment_floor",
            "exact_inequality",
            "ready_to_apply",
            "The full deformed Xi kernel has a uniform positive mass floor.",
            exact["zeroth_moment_floor"]["conclusion"],
            "Uses only the n=1 summand on 0<=u<=1/100.",
        ),
        GateRow(
            "ntcts_03_second_moment_ceiling",
            "exact_inequality",
            "ready_to_apply",
            "The full deformed Xi kernel has a uniform second-moment ceiling.",
            exact["second_moment_ceiling"]["conclusion"],
            "Uses a rational exponential envelope summed over every n.",
        ),
        GateRow(
            "ntcts_04_near_origin_no_contact",
            "exact_theorem",
            "ready_to_apply",
            "The Newman heat-flow transform is strictly positive near the origin.",
            exact["near_origin_theorem"]["lower_bound"],
            exact["near_origin_theorem"]["interval"],
        ),
        GateRow(
            "ntcts_05_component_contact_test",
            "exact_implication",
            "ready_to_apply",
            "The first theta block supplies a rigorous pointwise contact exclusion.",
            (
                "|J_1|>B_0 or |J_1'|>B_1 implies "
                "(J,J')!=(0,0)"
            ),
            "Imported from the independently checked 20-row theta/operator gate.",
        ),
        GateRow(
            "ntcts_06_compact_grid_scout",
            "numerical_evidence",
            "diagnostic_only",
            "The pointwise disjunction has a robust sampled compact margin.",
            (
                f"{compact['total_points']} points, "
                f"minimum ratio={compact['fine_minimum']['ratio']:.16g}, "
                f"failures={compact['fine_failures']}"
            ),
            "Point quadrature and a finite grid are not an interval certificate.",
        ),
        GateRow(
            "ntcts_07_quadrature_crosscheck",
            "numerical_crosscheck",
            "diagnostic_only",
            "Independent quadrature orders agree on the compact grid.",
            (
                "maximum coarse/fine normalized-ratio difference="
                f"{compact['maximum_coarse_fine_ratio_difference']:.6e}"
            ),
            "Agreement does not enclose truncation or between-grid variation.",
        ),
        GateRow(
            "ntcts_08_outer_scope_guard",
            "nonpromotion_gate",
            "guard_validated",
            "The raw first-block bars visibly lose their margin at larger frequency.",
            (
                f"{outer['failure_points']} failures on "
                f"{outer['total_points']} outer diagnostic points"
            ),
            outer["interpretation"],
        ),
        GateRow(
            "ntcts_09_dominant_saddle_composition",
            "proved_external_composition",
            "ready_to_apply",
            "The existing dominant-saddle certificate closes its complete ray.",
            "L>=50 and t*L>=25 imply L_t(x)>0",
            (
                "Uses the published Polymath-15 effective approximation and "
                "the checked cutoff-transition repair."
            ),
        ),
        GateRow(
            "ntcts_10_oscillatory_zeta_composition",
            "proved_asymptotic_composition",
            "ready_to_apply",
            "The existing oscillatory-zeta theorem closes every fixed scaled ray above c_*.",
            (
                "For every epsilon>0, sufficiently large L and "
                "t*L>=c_*+epsilon imply L_t(x)>0, "
                "c_*=4911678521/1933561194."
            ),
            "The threshold depends on epsilon and is not a bounded-L certificate.",
        ),
        GateRow(
            "ntcts_11_critical_phase_target",
            "open_theorem_target",
            "open",
            "The residual high-frequency layer requires corrected phase-critical-value avoidance.",
            (
                "L>=50, 0<t*L<=c_*+o(1): exclude "
                "theta=pi/2 mod pi together with theta'=0"
            ),
            "This arithmetic transversality theorem remains RH-level and open.",
        ),
        GateRow(
            "ntcts_12_compact_interval_upgrade",
            "open_certification_target",
            "open",
            "Upgrade the sampled first-block margin to a rigorous two-parameter interval cover.",
            (
                "Certify |J_1|>B_0 or |J_1'|>B_1 on "
                "1/4<=x<=38, 0<=t<=1/5"
            ),
            (
                "Requires interval quadrature or a Taylor/Cauchy enclosure "
                "between every sampled point."
            ),
        ),
    ]


def build_payload() -> dict:
    exact = build_exact()
    scout = build_scout()
    rows = build_rows(exact, scout)
    return {
        "kind": "jensen_window_pf_newman_theta_compact_transversality_scout",
        "date": "2026-07-24",
        "status": (
            "exact uniform near-origin no-contact theorem plus a "
            "non-rigorous compact first-block scout and high-frequency "
            "proof-partition composition; not a proof of Lambda<=0 or RH"
        ),
        "proof_boundary": (
            "Only |x|<=1/4 is proved here. The 61951-point compact grid is "
            "numerical evidence, not an interval certificate. Existing "
            "Polymath-15 gates close specified high-frequency regions, while "
            "the corrected critical phase layer and the intervening bounded "
            "region remain open."
        ),
        "exact": exact,
        "scout": scout,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    scout = payload["scout"]
    compact = scout["compact_grid"]
    outer = scout["outer_scope_guard"]
    success = (
        "validated Newman theta compact-transversality scout: 12 rows, "
        "0 issues, 3 exact moment/kernel inequalities, "
        "1 uniform near-origin no-contact theorem, "
        f"{compact['total_points']} compact diagnostic points, "
        f"{compact['fine_failures']} compact grid failures, "
        "1 independent quadrature crosscheck, 1 outer nonpromotion guard, "
        "2 high-frequency composition theorems, 2 open certification targets"
    )
    return "\n".join(
        [
            "# Newman Theta Compact-Transversality Scout",
            "",
            "Date: 2026-07-24",
            "",
            "Status: one exact near-origin theorem, one compact numerical",
            "scout, and an exact proof-partition handoff. This is not a proof",
            "of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_theta_compact_transversality_scout.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_theta_compact_transversality_scout.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_theta_compact_transversality_scout.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            success,
            "```",
            "",
            "## Exact Near-Origin Theorem",
            "",
            exact["kernel"]["definition"],
            "",
            exact["kernel"]["positivity"],
            "",
            "Put `M_k(t)=integral u^k*phi_t(u)du`. Rational exponential",
            "bounds give",
            "",
            "```text",
            exact["zeroth_moment_floor"]["conclusion"],
            exact["second_moment_ceiling"]["conclusion"],
            exact["near_origin_theorem"]["lower_bound"],
            "```",
            "",
            f"Since the final rational margin is `{exact['near_origin_theorem']['margin']}>0`,",
            "",
            "```text",
            exact["near_origin_theorem"]["interval"],
            exact["near_origin_theorem"]["conclusion"],
            "```",
            "",
            "This removes the artificial degeneracy caused by multiplying",
            "the contact functional by `x^4`.",
            "",
            "## Compact First-Block Scout",
            "",
            "The independently checked pointwise implication is",
            "",
            "```text",
            "|J_1|>B_0 or |J_1'|>B_1 implies (J,J')!=(0,0).",
            "```",
            "",
            "On the sampled rectangle",
            "",
            "```text",
            (
                f"0<=t<=1/5 in steps of {compact['time_step']}, "
                f"1/4<=x<=38 in steps of {compact['x_step']}"
            ),
            f"points={compact['total_points']}",
            (
                "minimum max(|J_1|/B_0,|J_1'|/B_1)="
                f"{compact['fine_minimum']['ratio']:.16g}"
            ),
            (
                f"minimum location: t={compact['fine_minimum']['t']}, "
                f"x={compact['fine_minimum']['x']}"
            ),
            f"grid failures={compact['fine_failures']}",
            (
                "maximum 600/900-node ratio difference="
                f"{compact['maximum_coarse_fine_ratio_difference']:.6e}"
            ),
            "```",
            "",
            "These are point values. They do not control the continuum between",
            "points and therefore are not promoted to a compact theorem.",
            "",
            "## Scope Guard",
            "",
            outer["interpretation"],
            "",
            "```text",
            f"outer diagnostic failures={outer['failure_points']}/{outer['total_points']}",
            *[
                (
                    f"t={row['t']}: first failure x={row['first_x']}, "
                    f"failure points={row['failure_count']}"
                )
                for row in outer["first_failures"]
            ],
            "```",
            "",
            "This confirms that raw moment dominance is a compact tool only.",
            "",
            "## Global Proof Partition",
            "",
            "The existing Polymath-15 certificates already supply the correct",
            "cancellation-preserving high-frequency framework:",
            "",
            "```text",
            "L=log(x/(4*pi)), c=t*L",
            "L>=50 and c>=25: exact dominant-saddle ray closed",
            "for every epsilon>0, c>=c_*+epsilon: asymptotically closed",
            "c_*=4911678521/1933561194",
            "L>=50 and 0<c<=c_*+o(1): corrected phase transversality open",
            "bounded/intervening L: compact certification still open",
            "```",
            "",
            "Endpoint subtraction should not replace this phase-aware machinery:",
            "absolute decay of both `H_t` and `H_t'` cannot itself exclude their",
            "simultaneous vanishing.",
            "",
            "## Live Handoff",
            "",
            "First promote the robust sampled rectangle to a rigorous",
            "two-parameter interval certificate using Arb quadrature or a",
            "Taylor/Cauchy enclosure. In parallel, retain the corrected",
            "Riemann-Siegel phase target on the residual scaled layer. Neither",
            "task may assume zero simplicity, `Lambda<=0`, or RH.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Newman theta compact-transversality scout: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
