#!/usr/bin/env python3
"""Build an exact normal form for the complete outer quartic contact chart."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_quartic_outer_contact_normal_form_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_quartic_outer_contact_normal_form_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def lower_process_priority() -> None:
    try:
        import psutil

        psutil.Process(os.getpid()).nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass


def require_zero(expression: sp.Expr, message: str) -> None:
    if sp.factor(expression) != 0:
        raise RuntimeError(message)


def build_exact() -> dict[str, object]:
    delta, q, s = sp.symbols("delta q s")
    a = 1 - delta
    curvature = 4 * q * delta**2
    p = sp.factor(4 * a - 3 * a**2 + curvature)
    simple_discriminant = sp.factor((4 - 2 * a) ** 2 - 4 * p)

    root_A = sp.factor(-3 * a**2 + 8 * a + p)
    root_B = sp.factor(-a**2 + 2 * a + p)
    K = sp.factor(root_A / 2)
    L = sp.factor(root_B / 2)
    x2 = sp.factor(root_A / 6)
    x3 = sp.factor(18 * a * root_B / root_A**2)
    x4 = sp.factor(2 * p * root_A / (3 * root_B**2))
    x5 = 1 - delta**2 * s
    threshold = sp.factor(
        root_B * (3 * a**2 - 5 * a + 5 * p) / (6 * p**2)
    )

    expected = {
        "p": 1 + 2 * delta + (4 * q - 3) * delta**2,
        "simple_discriminant": 16 * delta**2 * (1 - q),
        "K": 3 - (3 - 2 * q) * delta**2,
        "L": 1 + delta - 2 * (1 - q) * delta**2,
    }
    require_zero(p - expected["p"], "p normal form failed")
    require_zero(
        simple_discriminant - expected["simple_discriminant"],
        "simple-root discriminant normal form failed",
    )
    require_zero(K - expected["K"], "K normal form failed")
    require_zero(L - expected["L"], "L normal form failed")

    e2 = sp.factor((1 - x2) / delta**2)
    e3 = sp.factor((1 - x3) / delta**2)
    e4 = sp.factor((1 - x4) / delta**2)
    threshold_defect = sp.factor((1 - threshold) / delta**2)
    for index, contraction, defect in (
        (2, x2, e2),
        (3, x3, e3),
        (4, x4, e4),
    ):
        require_zero(
            1 - contraction - delta**2 * defect,
            f"x_{index} defect scaling failed",
        )
    require_zero(
        x5 - threshold - delta**2 * (threshold_defect - s),
        "outer-threshold scaling failed",
    )

    d2 = delta**2 * e2
    d3 = delta**2 * e3
    d4 = delta**2 * e4
    d5 = delta**2 * s
    gap1 = sp.factor(d3**2 - x3**2 * d2 * d4)
    gap2 = sp.factor(d4**2 - x4**2 * d3 * d5)
    h1 = sp.factor(gap1 / delta**6)
    h2 = sp.factor(gap2 / delta**4)
    expected_h1 = sp.factor(8 * q**3 / K**3)
    require_zero(h1 - expected_h1, "first gap scaling failed")

    h2_at_threshold = sp.factor(h2.subs(s, threshold_defect))
    expected_h2_at_threshold = sp.factor(
        8 * delta**2 * q**3 / (27 * L**3)
    )
    require_zero(
        h2_at_threshold - expected_h2_at_threshold,
        "threshold gap remainder failed",
    )
    require_zero(
        h2
        - (
            x4**2 * e3 * (threshold_defect - s)
            + expected_h2_at_threshold
        ),
        "second gap decomposition failed",
    )

    scaled_wall_gap = sp.factor(
        d5**2 - sp.Rational(11, 13) * d5 * x5**2 * d4
    )
    w6 = sp.factor(scaled_wall_gap / delta**4)
    cap6 = sp.factor(gap2**2 / (x4**3 * gap1))
    c6 = sp.factor(cap6 / delta**2)
    y6_upper = sp.factor(scaled_wall_gap / cap6)
    require_zero(
        y6_upper - delta**2 * w6 / c6,
        "first extension cap comparison failed",
    )

    rho, e_prev2, e_prev1, h_prev2, h_prev1 = sp.symbols(
        "rho e_prev2 e_prev1 h_prev2 h_prev1"
    )
    x_prev2 = 1 - delta**2 * e_prev2
    x_prev1 = 1 - delta**2 * e_prev1
    normalized_cap = sp.factor(
        h_prev1**2 / (x_prev2**3 * h_prev2)
    )
    next_gap = rho * normalized_cap
    next_defect = sp.factor(
        (e_prev1**2 - next_gap) / (x_prev1**2 * e_prev2)
    )
    physical_cap = sp.factor(
        (delta**4 * h_prev1) ** 2
        / (x_prev2**3 * delta**4 * h_prev2)
    )
    physical_defect = sp.factor(
        (
            (delta**2 * e_prev1) ** 2
            - rho * physical_cap
        )
        / (x_prev1**2 * delta**2 * e_prev2)
    )
    require_zero(
        physical_defect - delta**2 * next_defect,
        "normalized tail recurrence failed",
    )

    e12, e13, h9, h10, x12, x13 = sp.symbols(
        "e12 e13 h9 h10 x12 x13"
    )
    normalized_delta14 = sp.factor(
        h10**2 / (x12**3 * h9)
        - (e13**2 - x13**2 * e12 * e13)
    )
    physical_delta14 = sp.factor(
        (delta**4 * h10) ** 2
        / (x12**3 * delta**4 * h9)
        - (
            (delta**2 * e13) ** 2
            - x13**2 * (delta**2 * e12) * (delta**2 * e13)
        )
    )
    require_zero(
        physical_delta14 - delta**4 * normalized_delta14,
        "Delta_14 scaling failed",
    )

    r, qr = sp.symbols("r qr")
    right_polynomial = (
        4 * r**2 * qr**2
        - 12 * r**2 * qr
        + 9 * r**2
        + 6 * r * qr
        - 9 * r
        + 18 * qr
        - 18
    )
    right_endpoint_0 = sp.factor(right_polynomial.subs(qr, 0))
    right_endpoint_1 = sp.factor(right_polynomial.subs(qr, 1))
    right_second_q = sp.factor(sp.diff(right_polynomial, qr, 2))
    right_delta = sp.factor((x3 - x2).subs({delta: -r, q: qr}))
    right_factor = sp.factor(
        r**3
        * (3 + r * (3 - 2 * qr))
        * right_polynomial
        / (
            3
            * (
                3 - (3 - 2 * qr) * r**2
            )
            ** 2
        )
    )
    require_zero(
        right_delta - right_factor,
        "right-branch contraction factorization failed",
    )

    w, aa, pp, uu, epsilon = sp.symbols("w aa pp uu epsilon")
    contact_A = -3 * aa**2 + 8 * aa + pp
    contact_B = -aa**2 + 2 * aa + pp
    contact_x2 = contact_A / 6
    contact_x3 = 18 * aa * contact_B / contact_A**2
    contact_x4 = 2 * pp * contact_A / (3 * contact_B**2)
    contact_U = sp.factor(
        contact_B * (3 * aa**2 - 5 * aa + 5 * pp) / (6 * pp**2)
    )
    contact_P5 = (
        1
        + 5 * w
        + 10 * contact_x2 * w**2
        + 10 * contact_x2**2 * contact_x3 * w**3
        + 5
        * contact_x2**3
        * contact_x3**2
        * contact_x4
        * w**4
        + contact_x2**4
        * contact_x3**3
        * contact_x4**2
        * uu
        * w**5
    )
    contact_discriminant = sp.factor(
        sp.discriminant(contact_P5, w)
    )
    normalized_quintic_discriminant = sp.factor(
        contact_discriminant.subs(
            {
                aa: a,
                pp: p,
                uu: threshold + delta**2 * epsilon,
            }
        )
    )
    discriminant_prefactor = sp.factor(
        3125
        * delta**8
        * a**6
        * p**4
        / (9 * L**4)
    )
    discriminant_quadratic = sp.factor(
        normalized_quintic_discriminant
        / (discriminant_prefactor * epsilon**2)
    )
    q_factor = 15 * q**2 - 40 * q + 24
    R0 = sp.factor(
        256
        * delta**6
        * q**4
        * (16 * q - 15)
        * L**2
    )
    R1 = sp.factor(
        96 * delta**3 * a**3 * q_factor * L * p**2
    )
    R2 = sp.factor(9 * a**6 * p**4)
    discriminant_polynomial = sp.Poly(discriminant_quadratic, epsilon)
    for degree, expected_coefficient in enumerate((R0, R1, R2)):
        require_zero(
            discriminant_polynomial.coeff_monomial(epsilon**degree)
            - expected_coefficient,
            f"adjacent quintic discriminant R_{degree} failed",
        )
    slack_discriminant = sp.factor(R1**2 - 4 * R2 * R0)
    expected_slack_discriminant = sp.factor(
        147456
        * delta**6
        * a**6
        * (6 - q) ** 2
        * (1 - q) ** 3
        * L**2
        * p**4
    )
    require_zero(
        slack_discriminant - expected_slack_discriminant,
        "adjacent quintic slack discriminant failed",
    )
    positive_slack_root = (
        "16*delta^3*L/(3*(1-delta)^3*p^2)"
        "*(-(15*q^2-40*q+24)+4*(6-q)*(1-q)^(3/2))"
    )

    return {
        "chart": {
            "delta": "1-a",
            "q": "C/(4*delta^2)",
            "C": str(curvature),
            "a": str(a),
            "p": str(p),
            "simple_root_discriminant": str(simple_discriminant),
            "open_left_outer_domain": "0<delta<1 and 0<q<1",
        },
        "right_branch_exclusion": {
            "substitution": "delta=-r, 0<r<1, 0<q<1",
            "x3_minus_x2_factor": str(right_factor),
            "convex_polynomial": str(right_polynomial),
            "second_q_derivative": str(right_second_q),
            "q0_endpoint": str(right_endpoint_0),
            "q1_endpoint": str(right_endpoint_1),
            "conclusion": "x_3-x_2<0 on the right outer branch",
        },
        "denominators": {
            "K": str(K),
            "L": str(L),
            "p": str(p),
            "root_A": "2*K",
            "root_B": "2*L",
        },
        "contractions": {
            "x_2": str(x2),
            "x_3": str(x3),
            "x_4": str(x4),
            "x_5": str(x5),
        },
        "normalized_defects": {
            "e_2": str(e2),
            "e_3": str(e3),
            "e_4": str(e4),
            "e_5": str(s),
            "identity": "1-x_j=delta^2*e_j",
        },
        "outer_threshold": {
            "U": str(threshold),
            "T": str(threshold_defect),
            "identity": "x_5-U=delta^2*(T-s)",
            "outward_condition": "s<T",
        },
        "initial_gaps": {
            "G_1": "delta^6*h_1",
            "h_1": str(h1),
            "G_2": "delta^4*h_2",
            "h_2": str(h2),
            "h_2_decomposition": (
                "h_2=x_4^2*e_3*(T-s)+"
                + str(expected_h2_at_threshold)
            ),
        },
        "first_extension": {
            "scaled_wall_gap": "W_6=delta^4*w_6",
            "w_6": str(w6),
            "order4_cap": "C_6=delta^2*c_6",
            "c_6": str(c6),
            "Y_6": str(y6_upper),
            "identity": "Y_6=delta^2*w_6/c_6",
            "admissible_gap": "0<G_3<min(W_6,C_6)",
            "required_split": "Y_6<=1 versus Y_6>=1",
        },
        "normalized_tail": {
            "range": "k>=7",
            "gap_recurrence": (
                "h_(k-3)=y_k*h_(k-4)^2/"
                "(x_(k-2)^3*h_(k-5))"
            ),
            "defect_recurrence": (
                "e_k=(e_(k-1)^2-h_(k-3))/"
                "(x_(k-1)^2*e_(k-2))"
            ),
            "contraction": "x_k=1-delta^2*e_k",
        },
        "normalized_delta14": {
            "identity": "Delta_14=delta^4*delta14_hat",
            "delta14_hat": str(normalized_delta14),
        },
        "adjacent_quintic_discriminant": {
            "outward_slack": "epsilon=(x_5-U)/delta^2=T-s",
            "identity": "Disc(P_5)=positive_prefactor*epsilon^2*R(epsilon)",
            "positive_prefactor": str(discriminant_prefactor),
            "quadratic": "R(epsilon)=R_0+R_1*epsilon+R_2*epsilon^2",
            "R_0": str(R0),
            "R_1": str(R1),
            "R_2": str(R2),
            "quadratic_discriminant": str(slack_discriminant),
            "positive_root_for_q_below_15_over_16": positive_slack_root,
            "low_q_phase": (
                "0<q<15/16 => R_0<0<R_2, exactly one positive "
                "root epsilon_+, and Disc(P_5)<0 for 0<epsilon<epsilon_+"
            ),
            "contact_phase_boundary": "q=15/16",
        },
    }


def rows(exact: dict[str, object]) -> list[GateRow]:
    return [
        GateRow(
            "qocn_01_outer_chart",
            "exact_reduction",
            "proved_exact",
            "Every nondegenerate left outer double-root contact has a compact root chart delta=1-a and q=C/(4 delta^2).",
            "0<delta<1, 0<q<1",
            "Root-contact chart only; no Xi contact exclusion.",
            exact["chart"],
        ),
        GateRow(
            "qocn_02_right_exclusion",
            "exact_sign",
            "proved_exact",
            "The opposite outer-root branch is incompatible with the required increase x_2<x_3.",
            "delta<0 => x_3-x_2<0",
            "Uses contraction ordering as an explicit hypothesis, not as an automatic property of every quartic.",
            exact["right_branch_exclusion"],
        ),
        GateRow(
            "qocn_03_contact_map",
            "exact_identity",
            "proved_exact",
            "The complete contact prefix x_2 through x_5 is rational in delta,q,s.",
            "x_j=x_j(delta,q,s)",
            "The remaining scalar, cubic, and signed-minor inequalities still define the admissible subcell.",
            {
                "denominators": exact["denominators"],
                "contractions": exact["contractions"],
            },
        ),
        GateRow(
            "qocn_04_defect_scale",
            "exact_identity",
            "proved_exact",
            "All four contact defects share the exact scale delta^2.",
            "1-x_j=delta^2 e_j, j=2,3,4,5",
            "Normalized defects must still satisfy every imported strict inequality.",
            exact["normalized_defects"],
        ),
        GateRow(
            "qocn_05_outer_threshold",
            "exact_identity",
            "proved_exact",
            "The outward threshold becomes one normalized contact inequality.",
            "x_5-U=delta^2(T-s)",
            "The outward countercontact region is s<T; this does not prove such an Xi contact exists.",
            exact["outer_threshold"],
        ),
        GateRow(
            "qocn_06_initial_gaps",
            "exact_identity",
            "proved_exact",
            "The double-root contact forces an exceptional delta^6 first gap and a delta^4 second gap with an exact positive threshold remainder.",
            "G_1=delta^6 h_1; G_2=delta^4 h_2",
            "Strict positivity is required on the open admissible cell.",
            exact["initial_gaps"],
        ),
        GateRow(
            "qocn_07_first_extension",
            "exact_reduction",
            "proved_exact",
            "The first extension is controlled by the smaller of the order-four cap and the scaled-defect wall.",
            "0<G_3<min(W_6,C_6)",
            "A complete certificate must split the semialgebraic regimes Y_6<=1 and Y_6>=1.",
            exact["first_extension"],
        ),
        GateRow(
            "qocn_08_tail_homogeneity",
            "exact_reduction",
            "proved_exact",
            "After x_6, the complete signed tail recurrence closes in normalized defects and gaps.",
            "G_j=delta^4 h_j; d_j=delta^2 e_j",
            "Conditional on positive denominators and the strict corridor inequalities.",
            exact["normalized_tail"],
        ),
        GateRow(
            "qocn_09_delta14_scale",
            "exact_identity",
            "proved_exact",
            "The terminal compatibility has a removable contact scale.",
            "Delta_14=delta^4*delta14_hat",
            "The sign of delta14_hat over the full admissible contact-tail cell remains unproved.",
            exact["normalized_delta14"],
        ),
        GateRow(
            "qocn_10_quintic_discriminant",
            "exact_identity",
            "proved_exact",
            "The full adjacent-quintic discriminant factors into a positive contact prefactor, the square of the normalized outward slack, and one explicit quadratic in that slack.",
            "Disc(P_5)=positive_prefactor*epsilon^2*R(epsilon)",
            "Adjacent degree five on the quartic double-root chart only.",
            exact["adjacent_quintic_discriminant"],
        ),
        GateRow(
            "qocn_11_slack_quadratic",
            "exact_identity",
            "proved_exact",
            "The slack quadratic has three factored coefficients and a strictly positive discriminant throughout the open contact chart.",
            "Disc_epsilon(R)>0 for 0<delta,q<1",
            "This locates discriminant sign changes but does not prove adjacent-quintic hyperbolicity.",
            {
                key: exact["adjacent_quintic_discriminant"][key]
                for key in ("R_0", "R_1", "R_2", "quadratic_discriminant")
            },
        ),
        GateRow(
            "qocn_12_low_q_complex_collar",
            "exact_theorem",
            "proved_exact",
            "For 0<q<15/16, every sufficiently small strict outward slack lies in a negative-discriminant quintic collar.",
            "0<epsilon<epsilon_+(delta,q) => Disc(P_5)<0",
            "Contact-level adjacent-degree exclusion only; no Xi contact or uniform all-q collar is asserted.",
            {
                key: exact["adjacent_quintic_discriminant"][key]
                for key in (
                    "positive_root_for_q_below_15_over_16",
                    "low_q_phase",
                    "contact_phase_boundary",
                )
            },
        ),
        GateRow(
            "qocn_13_handoff",
            "scope_gate",
            "proved_exact",
            "The full-contact tail problem is reduced to bounded semialgebraic cells, while the adjacent quintic has an exact discriminant phase diagram in normalized outward slack.",
            "(delta,q,s), the split x_6 corridor, y_7,...,y_13, and R(T-s)",
            "Reduction and low-q degree-five exclusion only. No all-contact tail obstruction, Xi threshold, PF-infinity, Lambda<=0, RH, or Clay-prize conclusion.",
        ),
    ]


def build_payload() -> dict[str, object]:
    lower_process_priority()
    exact = build_exact()
    gate_rows = rows(exact)
    return {
        "kind": "jensen_window_pf_quartic_outer_contact_normal_form_gate",
        "date": "2026-07-25",
        "status": "exact outer-contact normal form with open contact-cell sign",
        "proof_boundary": (
            "This exact reduction removes the vanishing contact scale and "
            "identifies the two first-extension regimes and the exact "
            "adjacent-quintic discriminant phase. It does not certify "
            "Delta_14 on the complete contact family, prove the Xi threshold, "
            "PF-infinity, Lambda<=0, RH, or a Clay-prize result."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in gate_rows],
        "summary": {
            "rows": len(gate_rows),
            "contact_variables": 3,
            "excluded_outer_branches": 1,
            "normalized_defects": 4,
            "initial_gap_scalings": 2,
            "first_extension_regimes": 2,
            "normalized_terminal_compatibilities": 1,
            "adjacent_quintic_discriminant_factorizations": 1,
            "slack_quadratic_coefficients": 3,
            "low_q_negative_discriminant_collars": 1,
            "uniform_all_contact_theorems": 0,
        },
    }


def render_note(payload: dict[str, object]) -> str:
    exact = payload["exact"]
    return "\n".join(
        [
            "# Quartic Outer-Contact Normal-Form Gate",
            "",
            "Date: 2026-07-25",
            "",
            "Status: exact contact-level reduction; not a proof artifact for",
            "the complete contact-cell sign of `Delta_14`, which remains open.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_quartic_outer_contact_normal_form_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_quartic_outer_contact_normal_form_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_quartic_outer_contact_normal_form_gate.py",
            "```",
            "",
            "## Root Chart",
            "",
            "Put `delta=1-a` and `q=C/(4 delta^2)`. On the left outer branch,",
            "",
            "```text",
            "0<delta<1, 0<q<1,",
            f"p={exact['chart']['p']},",
            f"Disc_(b,c)={exact['chart']['simple_root_discriminant']}.",
            "```",
            "",
            "The right outer branch has `delta=-r`. The exact factorization of",
            "`x_3-x_2` contains a quadratic polynomial convex in `q`; both",
            "endpoint values are negative for `0<r<1`. Hence that branch has",
            "`x_3<x_2` and is excluded whenever the required contraction",
            "ordering is imposed.",
            "",
            "## Scale Removal",
            "",
            "The contact contractions have the exact form",
            "",
            "```text",
            "1-x_j=delta^2 e_j,  j=2,3,4,5,",
            "x_5-U=delta^2(T-s).",
            "```",
            "",
            "Thus the outward region is `s<T`. The first two order-three gaps",
            "have different contact scales:",
            "",
            "```text",
            "G_1=delta^6 h_1,",
            "G_2=delta^4 h_2,",
            f"h_1={exact['initial_gaps']['h_1']}.",
            "```",
            "",
            "The second normalized gap splits exactly into an outer-threshold",
            "slack plus a positive remainder:",
            "",
            "```text",
            exact["initial_gaps"]["h_2_decomposition"],
            "```",
            "",
            "## Adjacent Quintic Phase",
            "",
            "Put `epsilon=T-s=(x_5-U)/delta^2`. The adjacent normalized",
            "quintic has the exact discriminant factorization",
            "",
            "```text",
            "Disc(P_5)=positive_prefactor*epsilon^2*R(epsilon),",
            "R(epsilon)=R_0+R_1*epsilon+R_2*epsilon^2,",
            f"R_0={exact['adjacent_quintic_discriminant']['R_0']},",
            f"R_1={exact['adjacent_quintic_discriminant']['R_1']},",
            f"R_2={exact['adjacent_quintic_discriminant']['R_2']}.",
            "```",
            "",
            "The prefactor is strictly positive on the open contact chart,",
            "and",
            "",
            "```text",
            f"Disc_epsilon(R)={exact['adjacent_quintic_discriminant']['quadratic_discriminant']}>0.",
            "```",
            "",
            "For `0<q<15/16`, `R_0<0<R_2`. Hence `R` has exactly one",
            "positive root `epsilon_+`, and",
            "",
            "```text",
            "0<epsilon<epsilon_+ => Disc(P_5)<0.",
            "```",
            "",
            "This is an exact negative-discriminant collar of strict outward",
            "contacts. It proves adjacent degree-five nonhyperbolicity there,",
            "not that an Xi contact enters the collar.",
            "",
            "## First Extension",
            "",
            "At `x_6`, the scaled-defect wall and order-four cap compete:",
            "",
            "```text",
            "W_6=delta^4 w_6,",
            "C_6=delta^2 c_6,",
            "Y_6=W_6/C_6=delta^2 w_6/c_6,",
            "0<G_3<min(W_6,C_6).",
            "```",
            "",
            "A complete contact certificate must therefore split the exact",
            "regimes `Y_6<=1` and `Y_6>=1`; the fixed-contact theorem lies in",
            "the former regime.",
            "",
            "For every later step, normalized gaps and defects obey",
            "",
            "```text",
            exact["normalized_tail"]["gap_recurrence"],
            exact["normalized_tail"]["defect_recurrence"],
            "x_k=1-delta^2 e_k.",
            "```",
            "",
            "Finally,",
            "",
            "```text",
            "Delta_14=delta^4 delta14_hat.",
            "```",
            "",
            "This removes the degenerating contact scale. It does not prove",
            "that `delta14_hat<0` over every admissible contact-tail cell. The",
            "complete outer-contact obstruction, Xi threshold, PF-infinity,",
            "`Lambda<=0`, RH, and Clay-prize problem remain open.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote quartic outer-contact normal-form gate: "
        "13 rows, 3 contact variables, 1 excluded outer branch, "
        "4 normalized defects, 2 initial gap scalings, "
        "2 first-extension regimes, 1 normalized Delta_14, "
        "1 adjacent-quintic discriminant factorization, "
        "1 low-q negative-discriminant collar, "
        "0 uniform all-contact theorems"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
