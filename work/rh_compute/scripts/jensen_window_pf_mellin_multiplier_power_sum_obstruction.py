#!/usr/bin/env python3
"""Certify a power-sum obstruction for the natural Mellin interpolation.

This rejects the stronger continuous-index identity

    B(s) = product_j M_s^(alpha_j)

near s=0.  It does not reject an identity asserted only at integer indices;
that promotion would require a separate interpolation-uniqueness theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import comb, factorial
import json
from pathlib import Path
import sys


VENDOR = Path(__file__).resolve().parents[1] / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402
from flint import acb, arb, arb_mat  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_mellin_multiplier_power_sum_obstruction.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_mellin_multiplier_power_sum_obstruction.md"
DEFAULT_DPS = 100
DEFAULT_N_SUM = 10
DEFAULT_MAX_M = 10
DEFAULT_X_CUTOFF = 200
DEFAULT_U_CUTOFF = 8
DEFAULT_TAIL_RADIUS = "1e-45"
DEFAULT_ABS_TOL = "1e-65"
NEGATIVE_X_PANELS = ((0, 1), (1, 2), (2, 4), (4, 8), (8, 16), (16, 32), (32, 64), (64, 128), (128, 200))


@dataclass(frozen=True)
class Configuration:
    dps: int = DEFAULT_DPS
    n_sum: int = DEFAULT_N_SUM
    max_m: int = DEFAULT_MAX_M
    x_cutoff: int = DEFAULT_X_CUTOFF
    u_cutoff: int = DEFAULT_U_CUTOFF
    tail_radius: str = DEFAULT_TAIL_RADIUS
    abs_tol: str = DEFAULT_ABS_TOL


def upper_text(value: arb, digits: int = 30) -> str:
    return value.upper().str(digits, radius=False)


def lower_text(value: arb, digits: int = 30) -> str:
    return value.lower().str(digits, radius=False)


def phi_finite(u: acb, n_sum: int, pi: acb) -> acb:
    e4 = (4 * u).exp()
    e5 = (5 * u).exp()
    e9 = (9 * u).exp()
    total = acb(0)
    for n in range(1, n_sum + 1):
        nn = acb(n)
        total += (2 * pi * pi * nn**4 * e9 - 3 * pi * nn**2 * e5) * (
            -pi * nn**2 * e4
        ).exp()
    return total


def transformed_integral(a: object, b: object, m: int, negative: bool, cfg: Configuration) -> acb:
    pi = acb.pi()
    abs_tol = arb(cfg.abs_tol)

    def integrand(x: acb, analytic: bool) -> acb:
        del analytic
        u = (-x).exp() if negative else x.exp()
        jacobian = (-x).exp() if negative else x.exp()
        log_power = (-2 * x) ** m if negative else (2 * x) ** m
        return log_power * jacobian * phi_finite(u, cfg.n_sum, pi)

    return acb.integral(
        integrand,
        acb(str(a)),
        acb(str(b)),
        abs_tol=abs_tol,
        rel_tol=abs_tol,
        deg_limit=64,
        eval_limit=200000,
        depth_limit=64,
    )


def raw_log_moment_balls(cfg: Configuration) -> list[arb]:
    log_u_cutoff = arb(cfg.u_cutoff).log().mid()
    inflation = arb(f"[0 +/- {cfg.tail_radius}]")
    raw: list[arb] = []
    for m in range(cfg.max_m + 1):
        value = acb(0)
        for a, b in NEGATIVE_X_PANELS:
            value += transformed_integral(a, b, m, True, cfg)
        value += transformed_integral(0, log_u_cutoff, m, False, cfg)
        if not value.imag.contains(0):
            raise RuntimeError(f"log-moment integral has nonreal enclosure at m={m}: {value}")
        raw.append(value.real + inflation)
    return raw


def geometric_series_bound(p: int, first_n: int, pi: arb) -> tuple[arb, arb]:
    """Bound sum_{n>=first_n} n^p exp(-pi*n^2) by its first term and a ratio."""
    ratio = arb(2) ** p * (-pi * (2 * first_n + 1)).exp()
    first = arb(first_n) ** p * (-pi * first_n**2).exp()
    return first / (1 - ratio), ratio


def tail_diagnostics(cfg: Configuration) -> dict:
    pi = arb.pi()
    n_first = cfg.n_sum + 1
    s2_tail, r2_tail = geometric_series_bound(2, n_first, pi)
    s4_tail, r4_tail = geometric_series_bound(4, n_first, pi)
    u_cutoff = arb(cfg.u_cutoff)
    series_sup = (
        2 * pi**2 * (9 * u_cutoff).exp() * s4_tail
        + 3 * pi * (5 * u_cutoff).exp() * s2_tail
    )

    series_errors: list[arb] = []
    for m in range(cfg.max_m + 1):
        log_integral_bound = arb(2) ** m * (
            factorial(m) + (u_cutoff - 1) * u_cutoff.log() ** m
        )
        series_errors.append(series_sup * log_integral_bound)

    # On 0<=u<=1, use the same ratio argument from n=1 and a deliberately
    # generous round bound for the full kernel.
    s2_total, r2_total = geometric_series_bound(2, 1, pi)
    s4_total, r4_total = geometric_series_bound(4, 1, pi)
    phi_compact_bound = 2 * pi**2 * arb(9).exp() * s4_total + 3 * pi * arb(5).exp() * s2_total
    phi_round_bound = arb("1e10")
    if not phi_compact_bound < phi_round_bound:
        raise RuntimeError("compact Phi round bound failed")

    x_cutoff = arb(cfg.x_cutoff)
    small_u_errors = [
        phi_round_bound * arb(2) ** m * x_cutoff.gamma_upper(m + 1)
        for m in range(cfg.max_m + 1)
    ]

    # For u>=U, (log u)^m<=u^m and convex endpoint decay gives an
    # elementary one-sided integral bound.
    k2 = arb(1) / (1 - r2_total)
    k4 = arb(1) / (1 - r4_total)
    large_u_errors: list[arb] = []
    for m in range(cfg.max_m + 1):
        pieces = []
        for coefficient, kp, exponent in ((2 * pi**2, k4, 9), (3 * pi, k2, 5)):
            derivative = 4 * pi * (4 * u_cutoff).exp() - exponent - arb(m) / u_cutoff
            if not derivative > 0:
                raise RuntimeError(f"large-u endpoint derivative failed at m={m}, a={exponent}")
            endpoint = u_cutoff**m * (exponent * u_cutoff - pi * (4 * u_cutoff).exp()).exp()
            pieces.append(coefficient * kp * endpoint / derivative)
        large_u_errors.append(arb(2) ** m * sum(pieces, arb(0)))

    totals = [series_errors[m] + small_u_errors[m] + large_u_errors[m] for m in range(cfg.max_m + 1)]
    radius = arb(cfg.tail_radius)
    if not all(value < radius for value in totals):
        raise RuntimeError("analytic tail budget exceeds the declared inflation radius")
    return {
        "series_tail_ratio_p2_upper": upper_text(r2_tail),
        "series_tail_ratio_p4_upper": upper_text(r4_tail),
        "series_tail_sup_upper": upper_text(series_sup),
        "compact_phi_bound_upper": upper_text(phi_compact_bound),
        "compact_phi_round_bound": str(phi_round_bound),
        "max_series_error_upper": upper_text(max(series_errors, key=lambda x: float(x.mid()))),
        "max_small_u_error_upper": upper_text(max(small_u_errors, key=lambda x: float(x.mid()))),
        "max_large_u_error_upper": upper_text(max(large_u_errors, key=lambda x: float(x.mid()))),
        "max_total_error_upper": upper_text(max(totals, key=lambda x: float(x.mid()))),
        "inflation_radius": cfg.tail_radius,
        "all_tail_errors_below_inflation": True,
        "series_tail_formula": (
            "sum_(n>N)n^p*exp(-pi*n^2)<=a_(N+1)/(1-2^p*exp(-pi*(2*N+3)))"
        ),
        "small_u_tail_formula": "1e10*2^m*Gamma(m+1,X)",
        "large_u_tail_formula": (
            "2^m*C_p*U^m*exp(a*U-pi*exp(4*U))/(4*pi*exp(4*U)-a-m/U)"
        ),
    }


def power_sums_and_hankel(raw: list[arb], cfg: Configuration) -> tuple[dict[int, arb], list[dict]]:
    normalized = [value / raw[0] for value in raw]
    cumulants = [arb(0) for _ in raw]
    for n in range(1, len(raw)):
        correction = sum(
            (arb(comb(n - 1, j - 1)) * cumulants[j] * normalized[n - j] for j in range(1, n)),
            arb(0),
        )
        cumulants[n] = normalized[n] - correction

    powers: dict[int, arb] = {}
    for m in range(2, cfg.max_m + 1):
        polygamma_half = (-1) ** m * factorial(m - 1) * (2**m - 1) * arb(m).zeta()
        log_a_derivative = cumulants[m] - polygamma_half
        powers[m] = (-1) ** (m - 1) * log_a_derivative / factorial(m - 1)

    hankel_specs = ((2, 2), (2, 3), (2, 4), (3, 3), (4, 3), (7, 2))
    rows = []
    for shift, size in hankel_specs:
        matrix = arb_mat(
            [[powers[shift + i + j] for j in range(size)] for i in range(size)]
        )
        determinant = matrix.det()
        classification = "positive" if determinant > 0 else "negative" if determinant < 0 else "inconclusive"
        rows.append(
            {
                "shift": shift,
                "size": size,
                "determinant": str(determinant),
                "classification": classification,
                "contains_zero": bool(determinant.contains(0)),
            }
        )
    return powers, rows


def build_payload(cfg: Configuration = Configuration()) -> dict:
    flint.ctx.dps = cfg.dps
    tails = tail_diagnostics(cfg)
    raw = raw_log_moment_balls(cfg)
    powers, hankel = power_sums_and_hankel(raw, cfg)
    negative = [row for row in hankel if row["classification"] == "negative"]
    if len(negative) != 3:
        raise RuntimeError(f"expected three separated negative Hankel determinants, found {len(negative)}")
    rows = [
        {
            "id": "mmps_01_mellin_interpolation",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "The Gamma-normalized Mellin interpolation A(s) agrees with A_k at every nonnegative integer k.",
            "formula": "A(s)=2*sqrt(pi)*integral_0^infinity u^(2s)*Phi(u)du/(4^s*Gamma(s+1/2))",
        },
        {
            "id": "mmps_02_log_derivative_reduction",
            "role": "exact_reduction",
            "readiness": "available_exact",
            "claim": "For m>=2, geometric normalization drops out and the log derivative is a log-U cumulant minus the half-Gamma polygamma term.",
            "formula": "d^m log B(0)=kappa_m(2*log U)-psi_(m-1)(1/2)",
        },
        {
            "id": "mmps_03_continuous_product_necessity",
            "role": "exact_necessary_condition",
            "readiness": "available_exact",
            "claim": "A continuous-index elementary multiplier product forces p_m to be positive power sums and every shifted Hankel matrix of p_m to be positive semidefinite.",
            "formula": "p_m=(-1)^(m-1)*d^m log B(0)/(m-1)!=sum_j alpha_j^(-m)",
        },
        {
            "id": "mmps_04_arb_log_moments",
            "role": "interval_certificate",
            "readiness": "ready_to_apply",
            "claim": "Arb quadrature plus explicit series, small-u, and large-u tails encloses raw log moments through order 10.",
        },
        {
            "id": "mmps_05_negative_hankel_obstruction",
            "role": "interval_counterexample_gate",
            "readiness": "guard_validated",
            "claim": "Three shifted power-sum Hankel determinants are strictly negative, including the shift-2 size-4 determinant.",
        },
        {
            "id": "mmps_06_integer_uniqueness_boundary",
            "role": "proof_boundary",
            "readiness": "not_ready_to_apply",
            "claim": "Integer coefficient equality alone does not identify the multiplier product with this Mellin interpolation; a separate uniqueness theorem would be required.",
        },
    ]
    return {
        "kind": "jensen_window_pf_mellin_multiplier_power_sum_obstruction",
        "date": "2026-07-10",
        "status": "interval-certified continuous-interpolation obstruction",
        "proof_boundary": (
            "This rigorously rejects the natural continuous-index Mellin multiplier-product identity. "
            "It does not reject a product identity asserted only for integer coefficients, does not prove or disprove PF-infinity, Jensen hyperbolicity, RH, or Lambda <= 0."
        ),
        "configuration": {
            "dps": cfg.dps,
            "n_sum": cfg.n_sum,
            "max_log_moment_order": cfg.max_m,
            "negative_log_x_cutoff": cfg.x_cutoff,
            "positive_u_cutoff": cfg.u_cutoff,
            "quadrature_abs_tol": cfg.abs_tol,
            "analytic_tail_inflation": cfg.tail_radius,
            "negative_x_panels": [list(panel) for panel in NEGATIVE_X_PANELS],
        },
        "definitions": {
            "mellin_moment": "M(s)=integral_0^infinity u^(2s)*Phi(u)du",
            "coefficient_interpolation": "A(s)=2*sqrt(pi)*M(s)/(4^s*Gamma(s+1/2))",
            "geometric_normalization": "B(s)=A(s)*A(0)^(s-1)/A(1)^s",
            "raw_log_moment": "R_m=integral_0^infinity (2*log u)^m*Phi(u)du",
            "power_sum_candidate": "p_m=(-1)^(m-1)*(kappa_m-psi_(m-1)(1/2))/(m-1)!",
        },
        "tail_diagnostics": tails,
        "raw_log_moment_balls": [str(value) for value in raw],
        "power_sum_candidate_balls": {str(m): str(value) for m, value in powers.items()},
        "hankel_rows": hankel,
        "rows": rows,
        "source_artifacts": [
            "outputs/jensen_window_pf_multiplier_counting_measure_target.md",
            "outputs/jensen_window_pf_rank_two_boundary_family_lemma.md",
            "outputs/jensen_window_pf_defect_complete_monotonicity_scout.md",
        ],
        "summary": {
            "rows": len(rows),
            "raw_log_moments": len(raw),
            "power_sum_candidates": len(powers),
            "hankel_determinants": len(hankel),
            "positive_hankel_determinants": sum(row["classification"] == "positive" for row in hankel),
            "negative_hankel_determinants": len(negative),
            "inconclusive_hankel_determinants": sum(row["classification"] == "inconclusive" for row in hankel),
            "continuous_product_ruled_out": True,
            "integer_only_product_ruled_out": False,
            "ready_to_apply_rows": 1,
            "target_closing": False,
            "main_finding": (
                "The natural Gamma-normalized Mellin interpolation cannot itself be the canonical continuous elementary-multiplier product: "
                "its candidate power sums have three strictly negative shifted Hankel determinants. The integer-only counting-measure target remains open because no interpolation-uniqueness theorem is supplied."
            ),
        },
        "invariants": [
            "All promoted signs use Arb balls with analytic tail inflation.",
            "The continuous-index product is rejected by a necessary positive-semidefinite Hankel condition.",
            "No claim promotes equality at integers to equality of analytic interpolations.",
            "The discrete counting-measure target remains open.",
            "PF-infinity, Jensen hyperbolicity, RH, and Lambda <= 0 are not proved.",
        ],
    }


def result_line(payload: dict) -> str:
    summary = payload["summary"]
    return (
        "validated Jensen-window PF Mellin multiplier power-sum obstruction: "
        f"{summary['raw_log_moments']} log moments, {summary['power_sum_candidates']} power sums, "
        f"{summary['hankel_determinants']} Hankel determinants, {summary['negative_hankel_determinants']} negative, "
        f"{summary['inconclusive_hankel_determinants']} inconclusive, 1 continuous route ruled out, "
        "0 discrete routes ruled out, 0 issues"
    )


def write_note(payload: dict, path: Path) -> None:
    powers = payload["power_sum_candidate_balls"]
    hankel = payload["hankel_rows"]
    lines = [
        "# Jensen-Window PF Mellin Multiplier Power-Sum Obstruction",
        "",
        "Date: 2026-07-10",
        "",
        "Status: interval-certified continuous-interpolation obstruction. This is not a",
        "proof of PF-infinity, Jensen hyperbolicity, RH, or `Lambda <= 0`, and it does",
        "not rule out a multiplier product asserted only at integer coefficient indices.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_mellin_multiplier_power_sum_obstruction.json",
        "python work/rh_compute/scripts/jensen_window_pf_mellin_multiplier_power_sum_obstruction.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_mellin_multiplier_power_sum_obstruction.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        result_line(payload),
        "```",
        "",
        "## Exact Reduction",
        "",
        "For `Re(s)>-1/2`, define",
        "",
        "```text",
        "M(s)=integral_0^infinity u^(2s)*Phi(u)du,",
        "A(s)=2*sqrt(pi)*M(s)/(4^s*Gamma(s+1/2)),",
        "B(s)=A(s)*A(0)^(s-1)/A(1)^s.",
        "```",
        "",
        "The duplication formula gives `A(k)=A_k` at every nonnegative integer. For",
        "`m>=2`, geometric normalization is linear in `s`, so",
        "",
        "```text",
        "d^m log B(0)=kappa_m(2*log U)-psi_(m-1)(1/2).",
        "```",
        "",
        "If this interpolation were a continuous elementary multiplier product, then",
        "",
        "```text",
        "p_m=(-1)^(m-1)*d^m log B(0)/(m-1)!=sum_j alpha_j^(-m).",
        "```",
        "",
        "Every shifted Hankel matrix `[p_(r+i+j)]` would therefore be positive semidefinite.",
        "",
        "## Arb Certificate",
        "",
        "Arb integrates the finite `n<=10` kernel after `u=exp(-x)` on nine",
        "panels through `x=200`, and after `u=exp(t)` through `u=8`. A common",
        "`1e-45` radius covers the geometric `n`-tail, omitted small-`u` tail, and",
        "double-exponential large-`u` tail for log moments through order 10.",
        "",
        "Selected candidate power sums:",
        "",
        "```text",
        f"p_2 = {powers['2']}",
        f"p_4 = {powers['4']}",
        f"p_8 = {powers['8']}",
        "```",
        "",
        "Separated Hankel failures:",
        "",
        "```text",
    ]
    for row in hankel:
        if row["classification"] == "negative":
            lines.append(f"shift={row['shift']}, size={row['size']}: {row['determinant']}")
    lines.extend(
        [
            "```",
            "",
            "The shift-2 size-4 determinant is strictly negative, so the natural Mellin",
            "interpolation is not such a continuous multiplier product.",
            "",
            "## Proof Boundary",
            "",
            "The original target asks for equality of coefficient sequences at integer",
            "indices. Equality there does not by itself force equality with this Mellin",
            "interpolation between the integers. A Carlson-type or other uniqueness theorem",
            "with verified growth hypotheses would be needed for that promotion, and none is",
            "claimed here. The discrete counting-measure route therefore remains open, while",
            "its most direct continuous-index strengthening is rigorously closed.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_note(payload, args.note)
    print(result_line(payload))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
