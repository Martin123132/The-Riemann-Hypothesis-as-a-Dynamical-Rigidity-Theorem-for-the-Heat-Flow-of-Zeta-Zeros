#!/usr/bin/env python3
"""Reduce the lambda=-100 order-four tail to one continuous curvature bound."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
VENDOR = Path(__file__).resolve().parents[1] / "vendor"
for candidate in (SCRIPT_DIR, VENDOR):
    if candidate.exists() and str(candidate) not in sys.path:
        sys.path.insert(0, str(candidate))

import flint  # noqa: E402
import jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate as order3  # noqa: E402
from jensen_window_pf_negative_lambda_first_summand_leading_saddle_certificate import (  # noqa: E402
    potential_jet_arb,
)
from jensen_window_pf_negative_lambda_first_summand_paired_remainder_certificate import (  # noqa: E402
    certify_scaled_curvature_mode_block,
)


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_first_summand_curvature_bridge.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_compound_order4_first_summand_curvature_bridge.md"
)
COMPACT_CURVATURE_JSON = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_compound_order4_localized_curvature_compact_certificate.json"
)
DEFAULT_PRECISION_BITS = 1024
DEFAULT_MPMATH_DPS = 70
DEFAULT_SAMPLE_T = (319, 400, 700, 1000, 2000, 5000, 20000, 100000)
T_HEAT = 100


@dataclass(frozen=True)
class BridgeRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def load_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def compact_curvature_diagnostics() -> dict:
    artifact = load_json(COMPACT_CURVATURE_JSON)
    if artifact.get("kind") != (
        "jensen_window_pf_compound_order4_localized_curvature_compact_certificate"
    ):
        raise RuntimeError("localized-curvature compact source has the wrong kind")
    summary = artifact["summary"]
    handoff = artifact["mode_handoff"]
    if summary.get("all_blocks_passed") is not True:
        raise RuntimeError("localized-curvature compact source is not fully positive")
    return {
        "source": COMPACT_CURVATURE_JSON.relative_to(REPO_ROOT).as_posix(),
        "certified_t_range": handoff["certified_t_range"],
        "central_start_t_ball": handoff["central_start_t_ball"],
        "central_end_t_ball": handoff["central_end_t_ball"],
        "h_tile_count": artifact["parameters"]["tile_count"],
        "accepted_central_blocks": summary["accepted_central_blocks"],
        "worst_margin_lower": summary["worst_margin_lower"],
        "largest_scaled_localized_upper": summary[
            "largest_scaled_localized_upper"
        ],
        "remaining_ray": handoff["open_ray"],
    }


def sci(value: mp.mpf, digits: int = 36) -> str:
    return mp.nstr(value, n=digits, min_fixed=-6, max_fixed=6)


def arb_lower_text(value: flint.arb, digits: int = 60) -> str:
    return order3.lower_text(value, digits)


def moment_error(index: int) -> Fraction:
    """Bound |B_index-B_index^(1)| from 0<=log(1+delta_j)<=2/j^6."""
    return 2 * (
        Fraction(1, (index - 1) ** 6)
        + Fraction(2, index**6)
        + Fraction(1, (index + 1) ** 6)
    )


def log_defect_error(index: int) -> Fraction:
    return 4 * index * moment_error(index)


def gap_log_error(index: int) -> Fraction:
    return (
        2 * moment_error(index)
        + 2 * log_defect_error(index)
        + log_defect_error(index - 1)
        + log_defect_error(index + 1)
    )


def finite_j319_diagnostics() -> dict:
    """Retain an independent integer cross-check for the continuous J floor."""
    flint.ctx.prec = DEFAULT_PRECISION_BITS
    values, sources = order3.merged_order3_coefficients()
    ratios = {index: values[index + 1] / values[index] for index in range(321)}
    contractions = {
        index: ratios[index] / ratios[index - 1] for index in range(1, 321)
    }
    defects = {index: 1 - contractions[index] for index in contractions}
    center = 319
    full_j = -(
        contractions[center] ** 2
        * defects[center - 1]
        * defects[center + 1]
        / defects[center] ** 2
    ).log()

    path_floor_margins = {}
    for index in (318, 319, 320):
        full_b = -contractions[index].log()
        margin = full_b - flint.arb(moment_error(index).numerator) / moment_error(index).denominator
        margin -= flint.arb(1) / (4 * index)
        if not bool(margin > 0):
            raise RuntimeError(f"first/full B path floor failed at {index}")
        path_floor_margins[str(index)] = arb_lower_text(margin)

    q_bound = gap_log_error(center)
    first_j_lower = full_j - flint.arb(q_bound.numerator) / q_bound.denominator
    target = flint.arb(1) / (7 * center)
    target_margin = first_j_lower - target
    if not bool(target_margin > 0):
        raise RuntimeError("first-summand J_319 floor failed")
    return {
        "lambda": "-100",
        "center": center,
        "full_J_319_ball": full_j.str(60).replace("e", "E"),
        "moment_gap_log_error_bound": str(q_bound),
        "first_J_319_lower": arb_lower_text(first_j_lower),
        "target": "1/(7*319)",
        "target_margin_lower": arb_lower_text(target_margin),
        "path_B_floor": "B_theta,j>1/(4*j), j=318,319,320",
        "path_B_floor_margins": path_floor_margins,
        "merged_sources": sources,
    }


def lower_mode_collar_diagnostics() -> dict:
    """Extend H'''<0 and the scaled-curvature buffer below the old mode start."""
    flint.ctx.prec = 192
    blocks = (
        (Fraction(37, 40), Fraction(463, 500)),
        (Fraction(463, 500), Fraction(579, 625)),
    )
    accepted = []
    for left, right in blocks:
        result = certify_scaled_curvature_mode_block(
            left,
            right,
            200,
            window_y=6,
            eighth_envelope_bound=Fraction(1, 50_000),
        )
        if not result.get("passed") or not result.get("h_third_negative"):
            raise RuntimeError(f"lower-mode curvature collar failed at {left}..{right}")
        accepted.append({"left": str(left), "right": str(right), **result})
    start_t = potential_jet_arb(flint.arb("0.925"), 1)[1]
    end_t = potential_jet_arb(flint.arb("0.9264"), 1)[1]
    if not bool(start_t < 316 and end_t < 318):
        raise RuntimeError("lower-mode collar does not cover the required t range")
    minimum = min(
        accepted,
        key=lambda row: mp.mpf(row["shifted_curvature_buffer_lower"]),
    )
    return {
        "mode_interval": "0.925<=u<=0.9264",
        "block_count": len(accepted),
        "positive_blocks": len(accepted),
        "panels_per_block": 200,
        "window_y": 6,
        "eighth_envelope": "1/50000",
        "mode_start_t_ball": start_t.str(50).replace("e", "E"),
        "mode_end_t_ball": end_t.str(50).replace("e", "E"),
        "minimum_shifted_buffer_lower": minimum[
            "shifted_curvature_buffer_lower"
        ],
        "blocks": accepted,
        "consequence": (
            "Combined with the existing u>=0.9264 theorem, H'''<0 and the "
            "shifted scaled-curvature buffer hold for real parameter r>=316; "
            "hence B_1'<0 and ((2t+1)B_1)'>0 for t>=317."
        ),
    }


def exact_diagnostics() -> dict:
    n = sp.symbols("n", nonnegative=True)
    t = n + 319
    j_floor_polynomial = sp.Poly(
        sp.expand(4 * t**3 - 416 * t**2 - 87 * t + 63), n
    )
    k = n + 320
    transfer_polynomial = sp.Poly(
        sp.expand((k - 3) ** 6 - 25280 * k**2 * (k + 1) ** 2), n
    )
    if any(value <= 0 for value in j_floor_polynomial.all_coeffs()):
        raise RuntimeError("shifted J-floor polynomial is not coefficient-positive")
    if any(value <= 0 for value in transfer_polynomial.all_coeffs()):
        raise RuntimeError("shifted perturbation polynomial is not coefficient-positive")
    return {
        "center_coordinate": (
            "g(t)=d(t)^2-x(t)^2*d(t-1)*d(t+1), "
            "x(t)=exp(-B(t)), d(t)=1-x(t)"
        ),
        "stable_factorization": (
            "ell(t)=log(d(t)); J(t)=2*B(t)-ell(t-1)+2*ell(t)-ell(t+1); "
            "g(t)=d(t)^2*(1-exp(-J(t)))"
        ),
        "derivative_functions": (
            "phi(z)=1/(exp(z)-1), chi(z)=exp(z)/(exp(z)-1)^2"
        ),
        "log_defect_curvature": (
            "ell''=phi(B)*B''-chi(B)*(B')^2"
        ),
        "stable_log_gap_curvature": (
            "K(t)=(log g)''=2*ell''+phi(J)*J''-chi(J)*(J')^2"
        ),
        "j_derivatives": (
            "J^(r)=2*B^(r)-ell^(r)(t-1)+2*ell^(r)(t)-ell^(r)(t+1), r=0,1,2"
        ),
        "local_jets": (
            "j_0=2*B-ell''; j_1=2*B'-ell'''; j_2=2*B''-ell''''"
        ),
        "local_remainders": (
            "E_r=(1/12)*sup_|s|<=1 |ell^(r+4)(t+s)|; "
            "|J^(r)-j_r|<=E_r, r=0,1,2"
        ),
        "localized_upper": (
            "U=2*ell''+phi(j_0-E_0)*max(j_2+E_2,0)-"
            "chi(j_0+E_0)*max(abs(j_1)-E_1,0)^2"
        ),
        "localized_sufficient_condition": (
            "j_0>E_0 and U<=7/(2*t^2) imply J>0 and K<=7/(2*t^2)"
        ),
        "tent_identity": (
            "P_n=Delta^2 log(g)(k)=integral_[-1,1](1-|s|)*K(k+s)ds, k=n+3"
        ),
        "continuous_target": "K(t)<=7/(2*t^2) for every real t>=319",
        "first_summand_transfer": (
            "K(t)<=7/(2*t^2) => P_n^(1)<=7/(2*(k^2-1))<=18/(5*k^2), k>=320"
        ),
        "log_bound": "-log(1-z)<=z/(1-z), 0<=z<1",
        "first_summand_rational_margin": (
            "18/(5*k^2)-7/(2*(k^2-1))=(k^2-36)/(10*k^2*(k^2-1))"
        ),
        "integer_B_bounds": (
            "1/(2*m+1)<=B_1(m)<=3/(2*m-1), integer m>=319"
        ),
        "real_B_bounds": (
            "1/(2*t+3)<=B_1(t)<=3/(2*t-3), real t>=319"
        ),
        "j_floor_inputs": (
            "B_1'<0, ((2*t+1)*B_1)'>0, and "
            "B_1(t)-B_1(t+1)>=1/(4*t^2)"
        ),
        "j_floor_lower": (
            "J_1(t)>=2/(2*t+3)-2/(2*t-1)+(t-3)/(6*t^2)>=1/(7*t), t>=319"
        ),
        "j_floor_shifted_polynomial": str(j_floor_polynomial.as_expr()),
        "j_floor_shifted_coefficients_descending": [
            int(value) for value in j_floor_polynomial.all_coeffs()
        ],
        "moment_error": (
            "a_j=2*((j-1)^(-6)+2*j^(-6)+(j+1)^(-6)); |B_j-B_j^(1)|<=a_j"
        ),
        "coarse_error_bounds": (
            "a_j<=8/(j-1)^6; |ell_j-ell_j^(1)|<=32*j/(j-1)^6; "
            "|J_j-J_j^(1)|<=176*j/(j-2)^6<=1/(14*j)"
        ),
        "single_log_gap_error": (
            "|log(g_j)-log(g_j^(1))|<=2528*j^2/(j-2)^6, j>=319"
        ),
        "penalty_error": (
            "|P_n-P_n^(1)|<=10112*(k+1)^2/(k-3)^6<=2/(5*k^2), k>=320"
        ),
        "perturbation_shifted_polynomial": str(transfer_polynomial.as_expr()),
        "perturbation_shifted_coefficients_descending": [
            int(value) for value in transfer_polynomial.all_coeffs()
        ],
        "full_tail_transfer": (
            "P_n<=18/(5*k^2)+2/(5*k^2)=4/k^2, k=n+3>=320"
        ),
    }


def log_phi1(u: mp.mpf) -> mp.mpf:
    q = mp.pi * mp.exp(4 * u)
    return mp.log(mp.pi) + 5 * u + mp.log(2 * q - 3) - q


def saddle_derivative(u: mp.mpf, t: mp.mpf) -> mp.mpf:
    q = mp.pi * mp.exp(4 * u)
    return 2 * t / u - 2 * T_HEAT * u + 5 + 8 * q / (2 * q - 3) - 4 * q


def saddle_point(t: mp.mpf) -> mp.mpf:
    low = mp.mpf("0.01")
    high = mp.mpf(1)
    while saddle_derivative(high, t) > 0:
        high *= mp.mpf("1.4")
    for _ in range(220):
        middle = (low + high) / 2
        if saddle_derivative(middle, t) > 0:
            low = middle
        else:
            high = middle
    return (low + high) / 2


def h_jet(t: mp.mpf) -> tuple[mp.mpf, mp.mpf, mp.mpf]:
    saddle = saddle_point(t)
    maximum = 2 * t * mp.log(saddle) - T_HEAT * saddle**2 + log_phi1(saddle)
    center = 2 * mp.log(saddle)

    def weight(u: mp.mpf) -> mp.mpf:
        if not u:
            return mp.mpf(0)
        return mp.exp(
            2 * t * mp.log(u) - T_HEAT * u**2 + log_phi1(u) - maximum
        )

    upper = saddle + mp.mpf("1.3")
    moments = [
        mp.quad(
            lambda u, power=power: weight(u) * (2 * mp.log(u) - center) ** power,
            [0, saddle, upper],
        )
        for power in range(3)
    ]
    shifted_mean = moments[1] / moments[0]
    shifted_variance = moments[2] / moments[0] - shifted_mean**2
    log_moment = mp.log(2) + maximum + mp.log(moments[0])
    return (
        mp.loggamma(t + mp.mpf("0.5")) - log_moment,
        mp.digamma(t + mp.mpf("0.5")) - center - shifted_mean,
        mp.polygamma(1, t + mp.mpf("0.5")) - shifted_variance,
    )


def jet_add(left: tuple, right: tuple) -> tuple:
    return tuple(left[index] + right[index] for index in range(3))


def jet_scale(value: tuple, scalar: mp.mpf | int) -> tuple:
    return tuple(scalar * value[index] for index in range(3))


def jet_sub(left: tuple, right: tuple) -> tuple:
    return tuple(left[index] - right[index] for index in range(3))


def jet_mul(left: tuple, right: tuple) -> tuple:
    return (
        left[0] * right[0],
        left[1] * right[0] + left[0] * right[1],
        left[2] * right[0] + 2 * left[1] * right[1] + left[0] * right[2],
    )


def jet_exp(value: tuple) -> tuple:
    exponential = mp.exp(value[0])
    return (
        exponential,
        exponential * value[1],
        exponential * (value[2] + value[1] ** 2),
    )


def jet_log(value: tuple) -> tuple:
    ratio = value[1] / value[0]
    return (
        mp.log(value[0]),
        ratio,
        value[2] / value[0] - ratio**2,
    )


def scout_row(t_value: int) -> dict:
    t = mp.mpf(t_value)
    h = {shift: h_jet(t + shift) for shift in range(-2, 3)}
    b = {
        shift: jet_add(
            jet_add(h[shift + 1], jet_scale(h[shift], -2)), h[shift - 1]
        )
        for shift in (-1, 0, 1)
    }
    x = {shift: jet_exp(jet_scale(b[shift], -1)) for shift in (-1, 0, 1)}
    one = (mp.mpf(1), mp.mpf(0), mp.mpf(0))
    d = {shift: jet_sub(one, x[shift]) for shift in (-1, 0, 1)}
    ell = {shift: jet_log(d[shift]) for shift in (-1, 0, 1)}
    j_jet = jet_sub(
        jet_scale(b[0], 2),
        jet_add(jet_add(ell[-1], ell[1]), jet_scale(ell[0], -2)),
    )
    stable_tail = jet_log(jet_sub(one, jet_exp(jet_scale(j_jet, -1))))
    log_g = jet_add(jet_scale(ell[0], 2), stable_tail)
    g = jet_sub(
        jet_mul(d[0], d[0]),
        jet_mul(jet_mul(x[0], x[0]), jet_mul(d[-1], d[1])),
    )[0]
    scaled_curvature = t**2 * log_g[2]
    return {
        "t": t_value,
        "saddle": sci(saddle_point(t)),
        "B": sci(b[0][0]),
        "t_times_J": sci(t * j_jet[0]),
        "g": sci(g),
        "t2_log_g_second": sci(scaled_curvature),
        "margin_below_7_over_2": sci(mp.mpf("3.5") - scaled_curvature),
        "positive_gap": bool(g > 0 and j_jet[0] > 0),
        "below_target": bool(scaled_curvature < mp.mpf("3.5")),
        "proof_boundary": (
            "Finite-upper high-precision mpmath quadrature; not an interval enclosure "
            "or a uniform curvature theorem."
        ),
    }


def finite_scout() -> dict:
    mp.mp.dps = DEFAULT_MPMATH_DPS
    rows = [scout_row(value) for value in DEFAULT_SAMPLE_T]
    if not all(row["positive_gap"] and row["below_target"] for row in rows):
        raise RuntimeError("first-summand curvature scout crossed the target")
    worst = min(rows, key=lambda row: mp.mpf(row["margin_below_7_over_2"]))
    return {
        "mpmath_dps": DEFAULT_MPMATH_DPS,
        "sample_t": list(DEFAULT_SAMPLE_T),
        "rows": rows,
        "minimum_target_margin": worst["margin_below_7_over_2"],
        "minimum_target_margin_at_t": worst["t"],
        "observed_scaled_curvature_range": [
            min(row["t2_log_g_second"] for row in rows),
            max(row["t2_log_g_second"] for row in rows),
        ],
        "formal_asymptotic_scout": "t^2*(log g)'' appears to increase to 3 from below",
    }


def build_artifact() -> dict:
    exact = exact_diagnostics()
    finite_j319 = finite_j319_diagnostics()
    lower_mode_collar = lower_mode_collar_diagnostics()
    compact_curvature = compact_curvature_diagnostics()
    scout = finite_scout()
    rows = [
        BridgeRow(
            id="co4fcb_01_stable_gap_coordinate",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Factor the cancellation-prone order-three gap through a positive logarithmic gap coordinate.",
            formula=exact["stable_factorization"],
            proof_boundary="Exact first-summand coordinate only.",
        ),
        BridgeRow(
            id="co4fcb_02_curvature_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The logarithmic curvature of the stable gap uses only two derivatives of B and J.",
            formula=exact["stable_log_gap_curvature"],
            proof_boundary="Exact differentiation identity only.",
        ),
        BridgeRow(
            id="co4fcb_03_tent_bridge",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The discrete order-four penalty is a tent average of the continuous stable-gap curvature.",
            formula=exact["tent_identity"],
            proof_boundary="Exact continuous-to-discrete reduction.",
        ),
        BridgeRow(
            id="co4fcb_04_gap_floor",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Existing adjacent-wall and scaled-curvature theorems give a quantitative first-summand gap floor.",
            formula="J_1(t)>=1/(7*t), every real t>=319",
            proof_boundary="Uses the certified first-summand theorems and a two-block lower-mode interval collar.",
            diagnostics={
                "analytic": exact["j_floor_lower"],
                "lower_mode_collar": lower_mode_collar,
                "finite_crosscheck": finite_j319,
            },
        ),
        BridgeRow(
            id="co4fcb_04b_localized_curvature",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="A centered-difference Taylor bound localizes the remaining curvature ceiling to same-point derivatives and explicit envelopes.",
            formula=exact["localized_sufficient_condition"],
            proof_boundary="Exact sufficient localization only; its derivative envelopes remain to be certified.",
            diagnostics={
                "local_jets": exact["local_jets"],
                "remainders": exact["local_remainders"],
                "upper": exact["localized_upper"],
            },
        ),
        BridgeRow(
            id="co4fcb_04c_compact_curvature_certificate",
            role="interval_certificate",
            readiness="ready_to_apply",
            claim="Arb quadrature and the localized inequality prove the curvature ceiling from the required start through the u=2 handoff.",
            formula="K_1(t)<=7/(2*t^2), 319<=t<=V'(2)",
            proof_boundary="Rigorous compact real-parameter theorem only; the analytic mode ray u>=2 remains open.",
            diagnostics=compact_curvature,
        ),
        BridgeRow(
            id="co4fcb_05_continuous_curvature_target",
            role="open_analytic_handoff",
            readiness="not_ready_to_apply",
            claim="Extend the compact first-summand stable-gap curvature ceiling across the remaining asymptotic mode ray.",
            formula="K_1(t)<=7/(2*t^2) on the mode ray u>=2",
            proof_boundary="Sole open analytic entry condition in this bridge; the compact range is certified.",
        ),
        BridgeRow(
            id="co4fcb_06_first_summand_transfer",
            role="exact_conditional_bridge",
            readiness="conditional",
            claim="The continuous ceiling leaves a two-fifths inverse-square budget for the full-kernel perturbation.",
            formula=exact["first_summand_transfer"],
            proof_boundary="Conditional only on co4fcb_05.",
        ),
        BridgeRow(
            id="co4fcb_07_full_kernel_transfer",
            role="exact_perturbation_theorem",
            readiness="ready_to_apply",
            claim="The certified first-summand dominance changes the order-four penalty by at most two-fifths inverse square.",
            formula=exact["penalty_error"],
            proof_boundary="Exact perturbation theorem, using the proved gap floor; no curvature claim.",
        ),
        BridgeRow(
            id="co4fcb_08_all_tail_implication",
            role="exact_conditional_bridge",
            readiness="conditional",
            claim="The single continuous target implies the complete lambda=-100 order-four tail estimate.",
            formula=exact["full_tail_transfer"],
            proof_boundary="Conditional on co4fcb_05; contiguous order four only.",
        ),
        BridgeRow(
            id="co4fcb_09_finite_scout",
            role="finite_scout",
            readiness="finite_only",
            claim="High-precision samples clear the proposed continuous curvature ceiling with substantial room.",
            formula="t^2*(log g_1(t))''<7/2 on the recorded sample",
            proof_boundary="Non-interval finite evidence only.",
            diagnostics=scout,
        ),
    ]
    return {
        "kind": "jensen_window_pf_compound_order4_first_summand_curvature_bridge",
        "date": "2026-07-12",
        "status": (
            "exact order-four first-summand curvature reduction with proved gap floor, "
            "a compact curvature certificate, closed full-kernel perturbation budget, "
            "and one open analytic ray"
        ),
        "proof_boundary": (
            "This artifact proves the stable coordinate, the first-summand gap floor, "
            "the compact curvature ceiling through u=2, and the complete-kernel "
            "perturbation budget. It does not prove the remaining analytic curvature "
            "ray, order-four entry, forward "
            "invariance, arbitrary-column order four, PF-infinity, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_compound_order4_condensation_gate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_dominance_certificate.md",
            "outputs/jensen_window_pf_negative_lambda_m100_compound_order3_entry_certificate.md",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
        ],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_compound_order4_first_summand_curvature_bridge.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_compound_order4_first_summand_curvature_bridge.py"
        ),
        "exact": exact,
        "finite_j319": finite_j319,
        "lower_mode_collar": lower_mode_collar,
        "compact_curvature": compact_curvature,
        "scout": scout,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "bridge_rows": len(rows),
            "exact_identity_rows": 2,
            "exact_reduction_rows": 2,
            "localized_curvature_rows": 1,
            "compact_interval_certificate_rows": 1,
            "proved_gap_floor_rows": 1,
            "lower_mode_collar_blocks": lower_mode_collar["block_count"],
            "proved_perturbation_rows": 1,
            "conditional_bridge_rows": 2,
            "open_analytic_handoffs": 1,
            "finite_scout_rows": len(scout["rows"]),
            "ready_to_apply_rows": sum(
                row.readiness == "ready_to_apply" for row in rows
            ),
        },
    }


def write_note(path: Path, artifact: dict) -> None:
    exact = artifact["exact"]
    scout = artifact["scout"]
    finite = artifact["finite_j319"]
    lower_mode_collar = artifact["lower_mode_collar"]
    compact_curvature = artifact["compact_curvature"]
    lines = [
        "# Jensen-Window PF Compound Order-Four First-Summand Curvature Bridge",
        "",
        "Date: 2026-07-12",
        "",
        "Status: exact reduction with a proved gap floor, a compact curvature",
        "certificate, a closed full-kernel perturbation budget, and one open analytic",
        "ray. This is",
        "not a proof of order-four entry, PF-infinity, RH, or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_compound_order4_first_summand_curvature_bridge`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_compound_order4_first_summand_curvature_bridge.json",
        "python work/rh_compute/scripts/jensen_window_pf_compound_order4_first_summand_curvature_bridge.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_compound_order4_first_summand_curvature_bridge.py",
        "```",
        "",
        "## Stable Coordinate",
        "",
        "For the continuous first-summand family put",
        "",
        "```text",
        exact["center_coordinate"],
        "ell(t)=log(d(t)),",
        "J(t)=2*B(t)-ell(t-1)+2*ell(t)-ell(t+1).",
        "```",
        "",
        "Then the cancellation in the order-three gap factors exactly:",
        "",
        "```text",
        "g(t)=d(t)^2*(1-exp(-J(t))).",
        "```",
        "",
        "With `phi(z)=1/(exp(z)-1)` and",
        "`chi(z)=exp(z)/(exp(z)-1)^2`, differentiation gives",
        "",
        "```text",
        exact["log_defect_curvature"],
        exact["stable_log_gap_curvature"],
        exact["j_derivatives"],
        "```",
        "",
        "This replaces a difference of two terms of size `t^-2` by positive",
        "coordinates of sizes `d~t^-1` and `J~t^-1`.",
        "",
        "## Proved Gap Floor",
        "",
        "The existing first-summand cumulant wall, negative `H'''`, and continuous",
        "scaled-curvature growth imply the following. Two extra interval-Simpson",
        "blocks cover `0.925<=u<=0.9264`, with minimum shifted-buffer margin",
        f"`{lower_mode_collar['minimum_shifted_buffer_lower']}`. Combined with the",
        "existing compact and ray theorem, they give `B_1'<0` and increasing",
        "`(2t+1)B_1(t)` early enough for the complete real tail.",
        "",
        "```text",
        exact["integer_B_bounds"],
        exact["real_B_bounds"],
        exact["j_floor_lower"],
        "```",
        "",
        "After `t=319+n`, the final numerator is",
        "",
        "```text",
        exact["j_floor_shifted_polynomial"] + ".",
        "```",
        "",
        "Every coefficient is positive for `n>=0`. Independently, the repaired Arb",
        "coefficients and certified moment perturbation give the cross-check",
        "",
        "```text",
        f"J_1(319)-1/(7*319) > {finite['target_margin_lower']}.",
        "```",
        "",
        "Thus `J_1(t)>=1/(7*t)` for every real `t>=319`. This is already a",
        "theorem and is not part of the remaining curvature ceiling.",
        "",
        "## Localized Curvature Form",
        "",
        "For `r=0,1,2`, the tent representation and Taylor's theorem give",
        "",
        "```text",
        exact["local_jets"],
        exact["local_remainders"],
        "```",
        "",
        "Because both `phi` and `chi` decrease on the positive half-line,",
        "",
        "```text",
        exact["localized_upper"],
        exact["localized_sufficient_condition"],
        "```",
        "",
        "This is the intervalization coordinate: all central terms share the same",
        "real parameter, while shifted dependence is confined to three explicit",
        "derivative envelopes. It avoids independent enclosure of nearly equal",
        "adjacent moments on the large-parameter tail.",
        "",
        "## Compact Curvature Theorem",
        "",
        "A deterministic cache of outward-rounded Arb quadrature tiles now proves",
        "",
        "```text",
        "K_1(t)<=7/(2*t^2), 319<=t<=V'(2).",
        "```",
        "",
        f"The cache contains `{compact_curvature['h_tile_count']}` derivative tiles and",
        f"assembles `{compact_curvature['accepted_central_blocks']}` positive central blocks.",
        f"Its weakest margin is `{compact_curvature['worst_margin_lower']}`, and the",
        f"largest certified `t^2 U(t)` upper bound is `{compact_curvature['largest_scaled_localized_upper']}`.",
        "This is a real-parameter interval theorem, not finite sampling.",
        "",
        "## Remaining Analytic Ray",
        "",
        "For `k=n+3`, the centered second difference has the exact tent form",
        "",
        "```text",
        exact["tent_identity"],
        "```",
        "",
        "so the complete first-summand theorem would follow after proving only",
        "",
        "```text",
        "K_1(t)<=7/(2*t^2) on the mode ray u>=2.",
        "```",
        "",
        "The checked Gaussian-cumulant target expands the standardized tilted",
        "partition exactly through `epsilon^6`, records the alternating factorial",
        "cumulant signature through order eight, and gives explicit candidate",
        "corridors that clear seven conditional `t+-2` ray collars. A 1.8-million-block",
        "finite certificate and coefficient-positive asymptotic theorem now prove the",
        "formal corridors for every `u>=2`. Only the exact-minus-formal central",
        "remainder and two-tail estimates remain open.",
        "",
        "It would give",
        "",
        "```text",
        exact["first_summand_transfer"],
        "```",
        "",
        "using `-log(1-z)<=z/(1-z)` and the exact margin",
        "",
        "```text",
        exact["first_summand_rational_margin"] + ".",
        "```",
        "",
        "## Full-Kernel Budget",
        "",
        "Write `M_j=M_j^(1)*(1+delta_j)` with `0<=delta_j<=2/j^6`.",
        "The exact Lipschitz transfer through `B`, `ell`, and `J` gives",
        "",
        "```text",
        exact["moment_error"],
        exact["coarse_error_bounds"],
        exact["single_log_gap_error"],
        exact["penalty_error"],
        "```",
        "",
        "After `k=320+n`, the last comparison is certified by the",
        "coefficient-positive polynomial",
        "",
        "```text",
        exact["perturbation_shifted_polynomial"] + ".",
        "```",
        "",
        "Consequently the sole continuous target would imply",
        "",
        "```text",
        exact["full_tail_transfer"] + ".",
        "```",
        "",
        "## High-Precision Scout",
        "",
        "These saddle-centered rows use finite-upper mpmath quadrature. They are",
        "not interval enclosures and are not promoted to the continuous theorem.",
        "",
        "| t | saddle | t J(t) | t^2 (log g)'' | margin below 7/2 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in scout["rows"]:
        lines.append(
            f"| `{row['t']}` | `{row['saddle']}` | `{row['t_times_J']}` | "
            f"`{row['t2_log_g_second']}` | `{row['margin_below_7_over_2']}` |"
        )
    lines.extend(
        [
            "",
            "The sampled scaled curvature rises from about `1.795` toward the",
            "formal limit `3`, leaving more than `0.58` at the largest sample.",
            "The compact range is now rigorously interval-certified. The formal cumulant",
            "model and its first omitted coefficient layer are closed globally. The next",
            "proof task is the cancellation-preserving exact-density central and two-tail",
            "theorem at the explicit finite and asymptotic remainder budgets.",
            "",
            "```text",
            "outputs/jensen_window_pf_compound_order4_condensation_gate.md",
            "outputs/jensen_window_pf_negative_lambda_first_summand_cumulant_bridge.md",
            "outputs/jensen_window_pf_negative_lambda_scaled_curvature_continuous_bridge.md",
            "outputs/jensen_window_pf_compound_order4_localized_curvature_compact_certificate.md",
            "outputs/jensen_window_pf_compound_order4_gaussian_cumulant_ray_target.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_corridor_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_formal_cumulant_next_parity_asymptotic_ray_certificate.md",
            "outputs/jensen_window_pf_compound_order4_exact_cumulant_remainder_budget.md",
            "outputs/signed_hankel_jensen_dependency_graph.md",
            "```",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_note(args.note, artifact)
    print(
        "validated Jensen-window PF compound order-four first-summand curvature bridge: "
        "11 rows, 0 issues, 2 exact identities, 2 exact reductions, 1 proved "
        "gap floor, 1 compact interval theorem, 1 proved full-kernel perturbation "
        f"theorem, 1 open analytic ray, {len(DEFAULT_SAMPLE_T)} positive finite scouts"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
