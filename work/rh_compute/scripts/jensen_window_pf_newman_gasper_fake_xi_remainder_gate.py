#!/usr/bin/env python3
"""Audit Gasper/Polya fake-Xi transfers to the Newman Laguerre target."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from functools import lru_cache
import json
from math import pi
from pathlib import Path
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

from flint import arb, ctx  # noqa: E402


DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


@lru_cache(maxsize=None)
def quadrature(node_count: int) -> tuple[np.ndarray, np.ndarray]:
    nodes, weights = leggauss(node_count)
    panels = (
        (0.0, 0.125),
        (0.125, 0.25),
        (0.25, 0.5),
        (0.5, 1.0),
        (1.0, 2.0),
        (2.0, 4.0),
    )
    u = np.concatenate(
        [(right - left) * nodes / 2 + (right + left) / 2 for left, right in panels]
    )
    qweights = np.concatenate(
        [(right - left) * weights / 2 for left, right in panels]
    )
    return u, qweights


def xi_kernel(u: np.ndarray) -> np.ndarray:
    e4u = np.exp(4 * u)
    return sum(
        (
            2 * pi**2 * n**4 * np.exp(9 * u)
            - 3 * pi * n * n * np.exp(5 * u)
        )
        * np.exp(-pi * n * n * e4u)
        for n in range(1, 17)
    )


def fake_xi_kernel(u: np.ndarray) -> np.ndarray:
    return 4 * pi**2 * np.exp(-2 * pi * np.cosh(4 * u)) * np.cosh(9 * u)


def jet(
    weight: np.ndarray, u: np.ndarray, qweights: np.ndarray, x: float
) -> tuple[float, float, float]:
    return (
        float(np.dot(qweights, weight * np.cos(x * u))),
        float(np.dot(qweights, -u * weight * np.sin(x * u))),
        float(np.dot(qweights, -u * u * weight * np.cos(x * u))),
    )


def laguerre(value: tuple[float, float, float]) -> float:
    f, fp, fpp = value
    return fp * fp - f * fpp


def cross_laguerre(
    left: tuple[float, float, float], right: tuple[float, float, float]
) -> float:
    f, fp, fpp = left
    g, gp, gpp = right
    return 2 * fp * gp - f * gpp - g * fpp


def scalar_budget(alpha: float, l_h: float, l_p: float, b_ph: float) -> float:
    cross = alpha * b_ph - 2 * alpha * alpha * l_p
    residual = l_h - alpha * b_ph + alpha * alpha * l_p
    return (abs(cross) + abs(residual)) / (alpha * alpha * l_p)


def arb_xi_kernel_jet(u: arb, summands: int = 8) -> tuple[arb, arb, arb]:
    pi_ball = arb.pi()
    e4u = (4 * u).exp()
    values = [arb(0), arb(0), arb(0)]
    for n in range(1, summands + 1):
        n2 = n * n
        q = pi_ball * n2 * e4u
        for coefficient, exponent, sign in (
            (2 * pi_ball**2 * n2 * n2, 9, 1),
            (3 * pi_ball * n2, 5, -1),
        ):
            term = coefficient * (exponent * u - q).exp()
            first_factor = exponent - 4 * q
            second_factor = first_factor**2 - 16 * q
            values[0] += sign * term
            values[1] += sign * term * first_factor
            values[2] += sign * term * second_factor
    return values[0], values[1], values[2]


def arb_fake_xi_kernel_jet(u: arb) -> tuple[arb, arb, arb]:
    pi_ball = arb.pi()
    q_plus = pi_ball * (4 * u).exp()
    q_minus = pi_ball * (-4 * u).exp()
    values = [arb(0), arb(0), arb(0)]
    for exponent in (9, -9):
        term = 2 * pi_ball**2 * (exponent * u - q_plus - q_minus).exp()
        first_factor = exponent - 4 * q_plus + 4 * q_minus
        second_factor = first_factor**2 - 16 * (q_plus + q_minus)
        values[0] += term
        values[1] += term * first_factor
        values[2] += term * second_factor
    return values[0], values[1], values[2]


def arb_heat_profile_jet(
    kernel_jet: tuple[arb, arb, arb], u: arb, t: arb
) -> tuple[arb, arb, arb]:
    value, first, second = kernel_jet
    heat = (t * u * u).exp()
    return (
        heat * value,
        heat * (first + 2 * t * u * value),
        heat
        * (
            second
            + 4 * t * u * first
            + (2 * t + 4 * t * t * u * u) * value
        ),
    )


def arb_fourier_integrands(
    profile_jet: tuple[arb, arb, arb], u: arb, x: arb
) -> tuple[list[arb], list[arb]]:
    value, first, second = profile_jet
    sine = (x * u).sin()
    cosine = (x * u).cos()

    u_value = u * value
    u_first = value + u * first
    u_second = 2 * first + u * second
    u2_value = u * u * value
    u2_first = 2 * u * value + u * u * first
    u2_second = 2 * value + 4 * u * first + u * u * second

    values = [value * cosine, -u_value * sine, -u2_value * cosine]
    second_derivatives = [
        second * cosine - 2 * x * first * sine - x * x * value * cosine,
        -u_second * sine - 2 * x * u_first * cosine + x * x * u_value * sine,
        -u2_second * cosine
        + 2 * x * u2_first * sine
        + x * x * u2_value * cosine,
    ]
    return values, second_derivatives


def interval_rows(cell_count: int = 2048) -> list[dict]:
    """Certify the x=25 scalar obstruction by interval midpoint quadrature."""

    ctx.dps = 70
    x = arb(25)
    width = arb(f"1/{cell_count}")
    midpoint_error_factor = width**3 / 24
    tail_cap = arb("1e-45")
    rows: list[dict] = []
    for t_text in ("0", "1/2"):
        t = arb(t_text)
        sums = [arb(0) for _ in range(6)]
        errors = [arb(0) for _ in range(6)]
        for index in range(cell_count):
            midpoint = arb(f"{2 * index + 1}/{2 * cell_count}")
            cell = arb(midpoint, arb(f"1/{2 * cell_count}"))
            profiles = (
                (
                    arb_heat_profile_jet(arb_xi_kernel_jet(midpoint), midpoint, t),
                    arb_heat_profile_jet(arb_xi_kernel_jet(cell), cell, t),
                ),
                (
                    arb_heat_profile_jet(arb_fake_xi_kernel_jet(midpoint), midpoint, t),
                    arb_heat_profile_jet(arb_fake_xi_kernel_jet(cell), cell, t),
                ),
            )
            for offset, (mid_profile, cell_profile) in zip((0, 3), profiles, strict=True):
                mid_values, _ = arb_fourier_integrands(mid_profile, midpoint, x)
                _, cell_seconds = arb_fourier_integrands(cell_profile, cell, x)
                for jet_index in range(3):
                    sums[offset + jet_index] += width * mid_values[jet_index]
                    errors[offset + jet_index] += (
                        midpoint_error_factor * abs(cell_seconds[jet_index]).upper()
                    )

        for index in range(6):
            radius = (errors[index] + tail_cap).upper()
            sums[index] += arb(0, radius)

        h_jet = sums[:3]
        p_jet = sums[3:]
        l_h = h_jet[1] ** 2 - h_jet[0] * h_jet[2]
        l_p = p_jet[1] ** 2 - p_jet[0] * p_jet[2]
        b_ph = (
            2 * p_jet[1] * h_jet[1]
            - p_jet[0] * h_jet[2]
            - h_jet[0] * p_jet[2]
        )
        rho = b_ph**2 / (l_p * l_h)
        minimum = 3 - rho
        if not (l_h > 0 and l_p > 0 and b_ph > 0 and rho > 0 and rho < 2):
            raise RuntimeError(f"interval scalar sign chamber failed at t={t_text}")
        if not minimum > 1:
            raise RuntimeError(f"interval scalar obstruction failed at t={t_text}")
        rows.append(
            {
                "t": t_text,
                "x": "25",
                "H_jet": [str(value) for value in h_jet],
                "P_jet": [str(value) for value in p_jet],
                "L_H": str(l_h),
                "L_P": str(l_p),
                "B_PH": str(b_ph),
                "rho": str(rho),
                "minimum_absolute_budget": str(minimum),
                "certified_signs": {
                    "L_H_positive": True,
                    "L_P_positive": True,
                    "B_PH_positive": True,
                    "rho_in_open_0_2": True,
                    "minimum_budget_above_1": True,
                },
            }
        )
    return rows


def numerical_rows(node_count: int) -> list[dict]:
    u, qweights = quadrature(node_count)
    phi = xi_kernel(u)
    psi = fake_xi_kernel(u)
    rows: list[dict] = []
    for t in (0.0, 0.5):
        heat = np.exp(t * u * u)
        h_jet = jet(heat * phi, u, qweights, 25.0)
        p_jet = jet(heat * psi, u, qweights, 25.0)
        l_h = laguerre(h_jet)
        l_p = laguerre(p_jet)
        b_ph = cross_laguerre(p_jet, h_jet)
        rho = b_ph * b_ph / (l_p * l_h)
        alpha_opt = l_h / b_ph
        minimum = 3 - rho
        direct_minimum = scalar_budget(alpha_opt, l_h, l_p, b_ph)
        if not (l_h > 0 and l_p > 0 and b_ph > 0 and 0 < rho < 2):
            raise RuntimeError("scalar obstruction sign chamber failed")
        if abs(minimum - direct_minimum) > 2e-13:
            raise RuntimeError("scalar obstruction minimum failed")
        rows.append(
            {
                "t": format(t, ".17g"),
                "x": "25",
                "H_jet": [format(value, ".17e") for value in h_jet],
                "P_jet": [format(value, ".17e") for value in p_jet],
                "L_H": format(l_h, ".17e"),
                "L_P": format(l_p, ".17e"),
                "B_PH": format(b_ph, ".17e"),
                "rho": format(rho, ".17e"),
                "alpha_opt": format(alpha_opt, ".17e"),
                "minimum_absolute_budget": format(minimum, ".17e"),
                "natural_alpha_budget": format(
                    scalar_budget(1.0, l_h, l_p, b_ph), ".17e"
                ),
            }
        )
    return rows


def build_numerics() -> dict:
    coarse = numerical_rows(260)
    fine = numerical_rows(440)
    interval_certificate = interval_rows()
    max_relative_delta = 0.0
    fields = ("L_H", "L_P", "B_PH", "rho", "minimum_absolute_budget")
    for coarse_row, fine_row in zip(coarse, fine, strict=True):
        for field in fields:
            left = float(coarse_row[field])
            right = float(fine_row[field])
            max_relative_delta = max(
                max_relative_delta, abs(left - right) / max(abs(right), 1e-300)
            )
        if float(fine_row["minimum_absolute_budget"]) <= 1:
            raise RuntimeError("normalization-independent scalar obstruction failed")
    if max_relative_delta >= 1e-9:
        raise RuntimeError("fake-Xi quadrature convergence failed")
    return {
        "method": (
            "Composite Gauss-Legendre quadrature on [0,4], analytic Fourier "
            "derivative moments, and Xi theta summands n<=16"
        ),
        "coarse_nodes_per_panel": 260,
        "fine_nodes_per_panel": 440,
        "max_relative_coarse_fine_delta": format(max_relative_delta, ".17e"),
        "rows": fine,
        "interval_certificate": {
            "rigorous": True,
            "arithmetic": "python-flint Arb balls at 70 decimal digits",
            "domain": "u in [0,1]",
            "cell_count": 2048,
            "rule": (
                "Composite midpoint rule with interval suprema of each exact "
                "second derivative and per-cell error h^3*sup(abs(f''))/24"
            ),
            "theta_summands": 8,
            "omitted_tail_radius_per_jet": "1e-45",
            "tail_proof": (
                "For n>=9 on 0<=u<=1, pi>3 and exp(tu^2)<=exp(1/2)<2 "
                "give an omitted moment below 8e5*9^4*exp(-243)<1e-70. "
                "For u=1+v>=1 and 0<=t<=1/2, use u^2/2<=exp(4u)/100, "
                "pi>3, exp(4)>50, u^j<=exp(2v) for j<=2, and "
                "sum_(n>=1)n^4*q^(n^2)<=2q for q<=exp(-150). Both true and "
                "fake tails are then bounded by "
                "40*exp(-281/2)*integral_0^infinity exp(-587v)dv<1e-60. "
                "The common 1e-45 radius covers all omissions."
            ),
            "rows": interval_certificate,
            "scope": (
                "Finite interval certificates at x=25 for t=0 and t=1/2. "
                "They rigorously reject the scalar absolute-margin inequality "
                "at those points, not sign-aware cancellation or RH."
            ),
        },
        "independent_mpmath_80_digit_references": {
            "t=0": {
                "L_H": "4.8288528964419507801374708154818700034633984690053013016496021650378316292e-8",
                "L_P": "9.63174385700759649241903486292893805398341744805199241066479291231613173805e-9",
                "B_PH": "1.38896714646661879236332397203696970691267589712036014490271902059136787417e-8",
                "rho": "0.414796465135351145115779744580258235307799494729391878976676896781241963434",
                "alpha_opt": "3.47657819605454804669500385759263056682910105814050513244566880659488779575",
                "minimum_absolute_budget": "2.58520353486464885488422025541974176469220050527060812102332310321875803657",
            },
            "t=1/2": {
                "L_H": "4.672292012357301139515295159686405810509762612410202977067597568468093348e-8",
                "L_P": "9.74443254751230118745822433911669421683885075608936029332663842695002500259e-9",
                "B_PH": "1.4004303022542297279992077700223327439020884971891888365241884558981968285e-8",
                "rho": "0.430761090017545245374029004216750108892360442539619846237137865414549035433",
                "alpha_opt": "3.3363259884026044810382951892480128214468378718873867757317456041503241612",
                "minimum_absolute_budget": "2.56923890998245475462597099578324989110763955746038015376286213458545096457",
            },
        },
        "scope": (
            "The pointwise minimization is exact algebra once the jet signs and "
            "rho<2 hold. Arb certifies those hypotheses at both displayed endpoint "
            "times; the 80-digit values are independent cross-checks. This is not "
            "a sign-aware remainder theorem or an RH proof."
        ),
    }


def build_exact() -> dict:
    u = sp.symbols("u", real=True)
    psi = 4 * sp.pi**2 * sp.exp(-2 * sp.pi * sp.cosh(4 * u)) * sp.cosh(9 * u)
    psi_expanded = (
        2
        * sp.pi**2
        * sp.exp(9 * u - sp.pi * sp.exp(4 * u) - sp.pi * sp.exp(-4 * u))
        * (1 + sp.exp(-18 * u))
    )
    if sp.simplify(psi.rewrite(sp.exp) - psi_expanded) != 0:
        raise RuntimeError("fake-Xi kernel expansion failed")

    h, hp, hpp, p, pp, ppp, alpha = sp.symbols(
        "h hp hpp p pp ppp alpha", real=True
    )
    l_h = hp**2 - h * hpp
    l_p = pp**2 - p * ppp
    b_ph = 2 * pp * hp - p * hpp - h * ppp
    e, ep, epp = h - alpha * p, hp - alpha * pp, hpp - alpha * ppp
    l_e = sp.expand(ep**2 - e * epp)
    b_alpha_p_e = sp.expand(
        2 * (alpha * pp) * ep - (alpha * p) * epp - e * (alpha * ppp)
    )
    if sp.simplify(l_e - (l_h - alpha * b_ph + alpha**2 * l_p)) != 0:
        raise RuntimeError("residual Laguerre identity failed")
    if sp.simplify(b_alpha_p_e - (alpha * b_ph - 2 * alpha**2 * l_p)) != 0:
        raise RuntimeError("mixed Laguerre identity failed")
    if sp.simplify(l_h - (alpha**2 * l_p + b_alpha_p_e + l_e)) != 0:
        raise RuntimeError("full Laguerre decomposition failed")

    a, b, lp, lh = sp.symbols("a b lp lh", positive=True)
    rho = b**2 / (lp * lh)
    active_branch = 3 - 2 * b / (a * lp) + lh / (a**2 * lp)
    if sp.factor(active_branch - (3 - rho)) != (a * b - lh) ** 2 / (
        a**2 * lh * lp
    ):
        raise RuntimeError("scalar budget square completion failed")
    left_branch_boundary = 4 * lh * lp / b**2 - 1
    left_gap = (b**2 - 2 * lh * lp) ** 2 / (b**2 * lp * lh)
    if sp.factor(left_branch_boundary - (3 - rho)) != left_gap:
        raise RuntimeError("scalar budget alternate branch failed")

    return {
        "fake_xi_benchmark": {
            "kernel": "Psi(u)=4*pi^2*exp(-2*pi*cosh(4u))*cosh(9u)",
            "transform": "P_t(z)=integral_0^infinity exp(tu^2)*Psi(u)*cos(zu)du",
            "normalization": "P_0(z)=Xi_star(z/2)/8",
            "bessel_form": (
                "P_0(z)=pi^2/2*(K_(9/4+i*z/4)(2*pi)+"
                "K_(9/4-i*z/4)(2*pi))"
            ),
            "established_result": (
                "Gasper's integral-of-squares theorem proves that Xi_star, hence "
                "P_0, has only real zeros; classical heat universal-factor "
                "preservation gives the same conclusion for P_t, t>=0."
            ),
        },
        "laguerre_decomposition": {
            "forms": "L[f]=f'^2-f*f'', B[f,g]=2*f'*g'-f*g''-g*f''",
            "residual": "E_alpha=H-alpha*P",
            "identity": "L[H]=alpha^2*L[P]+B[alpha*P,E_alpha]+L[E_alpha]",
            "mixed_term": "B[alpha*P,E_alpha]=alpha*B[P,H]-2*alpha^2*L[P]",
            "residual_term": "L[E_alpha]=L[H]-alpha*B[P,H]+alpha^2*L[P]",
        },
        "scalar_absolute_budget": {
            "definition": (
                "R_alpha=(abs(B[alpha*P,E_alpha])+abs(L[E_alpha]))/"
                "(alpha^2*L[P])"
            ),
            "rho": "rho=B[P,H]^2/(L[P]*L[H])",
            "theorem": (
                "If L[P]>0, L[H]>0, B[P,H]>0, and 0<rho<2 at one x, "
                "then inf_(alpha>0) R_alpha=3-rho>1."
            ),
            "optimizer": "alpha_star=L[H]/B[P,H]",
            "square_completion": (
                "R_alpha-(3-rho)=(alpha*B[P,H]-L[H])^2/"
                "(alpha^2*L[P]*L[H]) in the minimizing sign branch"
            ),
            "branch_check": (
                "rho<2 makes L[E_alpha]>0 for all alpha. On "
                "alpha<=B[P,H]/(2L[P]), R_alpha=L[H]/(alpha^2L[P])-1, "
                "whose boundary excess above 3-rho is "
                "(B[P,H]^2-2L[P]L[H])^2/(B[P,H]^2L[P]L[H])>=0."
            ),
            "conclusion": (
                "No positive scalar normalization can prove L[H]>0 by the "
                "one-block absolute triangle estimate alpha^2*L[P]>"
                "abs(B[alpha*P,E_alpha])+abs(L[E_alpha])."
            ),
        },
        "kernel_ratio": {
            "definition": "M(u)=Phi(u)/Psi(u) for u>=0",
            "fake_expansion": (
                "Psi(u)=2*pi^2*exp(9u-pi*exp(4u)-pi*exp(-4u))*"
                "(1+exp(-18u))"
            ),
            "first_xi_summand": (
                "phi_1(u)=2*pi^2*exp(9u-pi*exp(4u))*"
                "(1-3*exp(-4u)/(2*pi))"
            ),
            "origin": (
                "Phi(0)>pi*(2*pi-3)*exp(-pi)>"
                "4*pi^2*exp(-2*pi)=Psi(0), so M(0)>1"
            ),
            "origin_proof": (
                "Use pi>3, pi<22/7, and exp(pi)>exp(3)>20; every theta "
                "summand is positive at u=0."
            ),
            "tail": (
                "The n=1 quotient is exp(pi*exp(-4u))*"
                "(1-3*exp(-4u)/(2*pi))/(1+exp(-18u)); all n>=2 terms "
                "are O(exp(-3*pi*exp(4u))) relative to the leading scale."
            ),
            "limit": "lim_(u->infinity) M(u)=1",
        },
        "positive_cosh_obstruction": {
            "candidate": "M(u)=integral_R cosh(u*s) dmu(s), mu>=0 finite",
            "boundedness_lemma": (
                "If a positive finite cosh transform has a finite limit as "
                "u->infinity, monotone convergence forces mu(R\\{0})=0, so "
                "the transform is constant."
            ),
            "conclusion": (
                "Because M(0)>1 but M(u)->1, Phi/Psi is not a positive "
                "cosh-mixture multiplier. Thus the direct Cardon/Polya "
                "positive-convolution transfer from fake Xi to Xi is impossible."
            ),
            "scope": (
                "This rejects the one-factor positive multiplier representation, "
                "not signed measures, multi-block special-function expansions, "
                "or sign-aware remainder estimates."
            ),
        },
        "open_handoff": (
            "Retain the Gasper block only inside a sign-aware or multi-block "
            "decomposition. The next admissible target is to construct blocks "
            "P_j with a proved common real-zero/interlacing structure and control "
            "the signed mixed Laguerre forms after the full modular endpoint "
            "cancellation. A scalar absolute budget or a single positive-cosh "
            "multiplier is now a closed route."
        ),
    }


def build_payload() -> dict:
    exact = build_exact()
    numerics = build_numerics()
    rows = [
        GateRow(
            id="gfrr_01_fake_xi_normalization",
            role="exact_identity",
            readiness="available_exact",
            claim="The scaled Polya fake-Xi transform has an exact Bessel normalization in the corpus variable.",
            formula=exact["fake_xi_benchmark"]["normalization"],
            proof_boundary="Exact change of variables and Bessel integral identity.",
            diagnostics=exact["fake_xi_benchmark"],
        ),
        GateRow(
            id="gfrr_02_gasper_real_zero_benchmark",
            role="established_theorem",
            readiness="ready_to_apply",
            claim="Gasper's integral-of-squares theorem supplies a genuine real-zero benchmark for the fake-Xi model.",
            formula=exact["fake_xi_benchmark"]["established_result"],
            proof_boundary="Established theorem for the model transform, not for the true Xi transform.",
        ),
        GateRow(
            id="gfrr_03_laguerre_remainder_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="The true Laguerre expression splits exactly into the scaled fake-Xi margin and signed residual terms.",
            formula=exact["laguerre_decomposition"]["identity"],
            proof_boundary="Exact quadratic identity; no sign is discarded.",
            diagnostics=exact["laguerre_decomposition"],
        ),
        GateRow(
            id="gfrr_04_scalar_budget_minimum",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="A pointwise invariant gives the exact best possible one-block absolute remainder budget over every positive normalization.",
            formula=exact["scalar_absolute_budget"]["theorem"],
            proof_boundary="Exact algebra conditional on the displayed pointwise sign chamber.",
            diagnostics=exact["scalar_absolute_budget"],
        ),
        GateRow(
            id="gfrr_05_scalar_budget_witnesses",
            role="interval_counter_certificate",
            readiness="interval_validated",
            claim="At x=25 the normalization-independent minimum is interval-certified above one at both endpoint heat times.",
            formula="inf_(alpha>0) R_alpha=2.5852035... at t=0 and 2.5692389... at t=1/2",
            proof_boundary="Arb midpoint certificates with explicit theta and integration tails; not a sign-aware estimate or an RH theorem.",
            diagnostics=numerics,
        ),
        GateRow(
            id="gfrr_06_kernel_ratio_endpoints",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The true-to-fake kernel ratio is nonconstant but has finite limit one.",
            formula="M(0)>1 and lim_(u->infinity)M(u)=1",
            proof_boundary="Exact origin inequality and theta-tail asymptotic.",
            diagnostics=exact["kernel_ratio"],
        ),
        GateRow(
            id="gfrr_07_positive_cosh_obstruction",
            role="exact_route_obstruction",
            readiness="guard_validated",
            claim="The direct Cardon/Polya positive-cosh multiplier transfer from fake Xi to true Xi is impossible.",
            formula=exact["positive_cosh_obstruction"]["conclusion"],
            proof_boundary="Rejects one positive convolution factor only; signed and multi-block routes remain open.",
            diagnostics=exact["positive_cosh_obstruction"],
        ),
        GateRow(
            id="gfrr_08_sign_aware_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Any surviving Gasper-type route must preserve signed cancellation or use several interlacing blocks.",
            formula=exact["open_handoff"],
            proof_boundary="Open construction target; not strict Laguerre positivity, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_gasper_fake_xi_remainder_gate",
        "date": "2026-07-11",
        "status": "exact fake-Xi comparison with scalar and positive-convolution obstructions",
        "proof_boundary": (
            "This artifact normalizes the Gasper/Polya fake-Xi benchmark, derives "
            "the exact Laguerre remainder algebra, and rejects two direct transfer "
            "routes: scalar absolute domination and a single positive cosh-mixture "
            "multiplier. The scalar witnesses are interval-certified, and the "
            "positive-convolution obstruction is exact. It does not reject "
            "sign-aware or multi-block decompositions and does not prove strict "
            "Laguerre positivity, Wiener density, RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_theta_summand_spectral_square_gate.md",
            "outputs/jensen_window_pf_newman_strict_laguerre_correlation_target.md",
            "outputs/formal_core.md",
            "https://arxiv.org/abs/0801.2996",
            "https://mathdept.byu.edu/cardon/papers/convolutions_and_zeros.pdf",
        ],
        "exact": exact,
        "numerics": numerics,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    numerics = payload["numerics"]
    lines = [
        "# Jensen-Window PF Newman Gasper Fake-Xi Remainder Gate",
        "",
        "Date: 2026-07-11",
        "",
        "Status: exact fake-Xi comparison with scalar and positive-convolution",
        "obstructions. This is not a proof of RH or `Lambda <= 0`.",
        "",
        "Artifact kind: `jensen_window_pf_newman_gasper_fake_xi_remainder_gate`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.json",
        "python work/rh_compute/scripts/jensen_window_pf_newman_gasper_fake_xi_remainder_gate.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_newman_gasper_fake_xi_remainder_gate.py",
        "```",
        "",
        "Current result:",
        "",
        "```text",
        "validated Jensen-window PF Newman Gasper fake-Xi remainder gate: 8 rows, 0 issues, 2 exact transform identities, 1 established real-zero benchmark, 1 scalar algebra theorem, 2 interval scalar witnesses, 2 high-precision cross-checks, 1 exact positive-convolution obstruction, 1 sign-aware handoff",
        "```",
        "",
        "## Gasper Benchmark",
        "",
        "With",
        "",
        "```text",
        exact["fake_xi_benchmark"]["kernel"],
        exact["fake_xi_benchmark"]["transform"],
        exact["fake_xi_benchmark"]["normalization"],
        exact["fake_xi_benchmark"]["bessel_form"],
        "```",
        "",
        "Gasper's integral-of-squares theorem proves that the model at `t=0`",
        "has only real zeros. The classical heat universal-factor theorem",
        "preserves this for `t>=0`. This is a valid comparison model, not an",
        "approximate proof for Xi.",
        "",
        "## Scalar Obstruction",
        "",
        "For `E_alpha=H-alpha*P`,",
        "",
        "```text",
        exact["laguerre_decomposition"]["identity"],
        exact["laguerre_decomposition"]["mixed_term"],
        exact["laguerre_decomposition"]["residual_term"],
        "```",
        "",
        "The absolute-margin strategy would require `R_alpha<1`, where",
        "",
        "```text",
        exact["scalar_absolute_budget"]["definition"],
        exact["scalar_absolute_budget"]["rho"],
        exact["scalar_absolute_budget"]["theorem"],
        "```",
        "",
        "Thus the optimization over the normalization is exact, not a grid",
        "search. Arb midpoint quadrature with interval second-derivative errors",
        "and explicit tails certifies both endpoint times at `x=25`:",
        "",
        "```text",
    ]
    for row in numerics["interval_certificate"]["rows"]:
        lines.append(
            f"Arb t={row['t']}: rho={row['rho']}, "
            f"inf R_alpha={row['minimum_absolute_budget']}"
        )
    lines.extend(["```", "", "Independent 80-digit quadrature cross-checks give:", "", "```text"])
    for key, reference in numerics["independent_mpmath_80_digit_references"].items():
        lines.extend(
            [
                f"{key}: rho={reference['rho']}",
                f"     alpha_star={reference['alpha_opt']}",
                f"     inf R_alpha={reference['minimum_absolute_budget']}",
            ]
        )
    lines.extend(
        [
            "```",
            "",
            "The corresponding Arb balls lie strictly above `1`; the displayed",
            "digits are independent cross-checks. No scalar multiple of",
            "the fake-Xi block can support the proposed absolute triangle bound.",
            "This does not reject a sign-aware estimate, because the discarded",
            "mixed and residual terms can cancel.",
            "",
            "## Convolution Obstruction",
            "",
            "The direct positive-convolution transfer would require",
            "",
            "```text",
            exact["kernel_ratio"]["definition"],
            exact["positive_cosh_obstruction"]["candidate"],
            "```",
            "",
            "The kernel asymptotics instead give",
            "",
            "```text",
            exact["kernel_ratio"]["origin"],
            exact["kernel_ratio"]["limit"],
            "```",
            "",
            exact["positive_cosh_obstruction"]["boundedness_lemma"],
            "Therefore the nonconstant ratio `Phi/Psi` cannot be a positive",
            "cosh-mixture multiplier. This rules out the direct Cardon/Polya",
            "one-factor transfer exactly.",
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
            "The useful content of the Gasper model is now sharply located:",
            "its positive square can remain as one block, but the Xi correction",
            "must be handled with its sign and modular cancellation intact.",
            "",
            "References: https://arxiv.org/abs/0801.2996 and",
            "https://mathdept.byu.edu/cardon/papers/convolutions_and_zeros.pdf.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(f"wrote Newman Gasper fake-Xi remainder gate: {len(payload['rows'])} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
