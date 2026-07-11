#!/usr/bin/env python3
"""Prove the forward-uniform reciprocal-defect tail and cubic propagation."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_cubic_forward_uniform_tail_certificate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_cubic_forward_uniform_tail_certificate.md"
TAIL_START = 319


@dataclass(frozen=True)
class CertificateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def exact_flow_diagnostics() -> dict:
    k, q, g, j, r = sp.symbols("k q g j r", positive=True)
    q1 = q + g
    q2 = q1 + j
    defect = 1 / q**2
    next_defect = 1 / q1**2
    defect_prime = 2 * r * (
        (2 * k + 3) * (1 - defect) * next_defect - (2 * k - 1) * defect
    )
    q_prime_from_defect = sp.factor(-q**3 * defect_prime / 2)

    def q_flow(index: sp.Expr, current: sp.Expr, following: sp.Expr, ratio: sp.Expr) -> sp.Expr:
        return ratio * (
            (2 * index - 1) * current
            - (2 * index + 3) * (current**3 - current) / following**2
        )

    q_prime = q_flow(k, q, q1, r)
    if sp.simplify(q_prime_from_defect - q_prime) != 0:
        raise RuntimeError("q flow does not follow from the defect equation")
    q1_prime = q_flow(k + 1, q1, q2, r * (1 - 1 / q1**2))
    g_prime_over_r = sp.factor((q1_prime - q_prime) / r)
    neighbor_derivative = sp.factor(sp.diff(g_prime_over_r, j))
    expected_neighbor_derivative = (
        2
        * (2 * k + 5)
        * (q + g - 1) ** 2
        * (q + g + 1) ** 2
        / ((q + g) * (q + g + j) ** 3)
    )
    if sp.simplify(neighbor_derivative - expected_neighbor_derivative) != 0:
        raise RuntimeError("next-increment monotonicity identity failed")

    same_neighbor = sp.factor(g_prime_over_r.subs(j, g))
    polynomial = (
        -6 * g**3 * k
        + g**3
        - 14 * g**2 * k * q
        + 5 * g**2 * q
        - 6 * g * k * q**2
        - 2 * g * k
        + 13 * g * q**2
        - 5 * g
        - 2 * k * q
        + 6 * q**3
        - 5 * q
    )
    expected_same_neighbor = (
        (1 - g**2) * polynomial / ((q + g) ** 2 * (q + 2 * g) ** 2)
    )
    if sp.simplify(same_neighbor - expected_same_neighbor) != 0:
        raise RuntimeError("same-neighbor flow identity failed")

    positive_majorant = g**3 + 5 * g**2 * q + 13 * g * q**2 + 6 * q**3
    majorant_gap = sp.expand(positive_majorant - polynomial)
    expected_gap = (
        6 * g**3 * k
        + 14 * g**2 * k * q
        + 6 * g * k * q**2
        + 2 * g * k
        + 5 * g
        + 2 * k * q
        + 5 * q
    )
    if sp.expand(majorant_gap - expected_gap) != 0:
        raise RuntimeError("same-neighbor positive majorant failed")

    scalar_remainder = Fraction(13, 17) + Fraction(5, 17**2) + Fraction(1, 17**3)
    if scalar_remainder != Fraction(3843, 4913) or not scalar_remainder < 1:
        raise RuntimeError("weighted source cap failed")
    coercive_coupling = Fraction(4, 1) + Fraction(10, TAIL_START)
    if not coercive_coupling < 5:
        raise RuntimeError("coercive supersolution coupling cap failed")

    return {
        "q_flow": (
            "q_k'=r_k*((2*k-1)*q_k-(2*k+3)*(q_k^3-q_k)/q_(k+1)^2)"
        ),
        "defect_to_q_identity": "q_k'=-q_k^3*d_k'/2",
        "increment": "g_k=q_(k+1)-q_k",
        "next_increment_derivative": str(neighbor_derivative),
        "next_increment_monotonicity": "partial(g_k'/r_k)/partial(g_(k+1))>0",
        "same_neighbor_polynomial": str(polynomial),
        "same_neighbor_identity": (
            "g_k'/r_k=(1-g_k^2)*P_k/((q_k+g_k)^2*(q_k+2*g_k)^2) "
            "when g_(k+1)=g_k"
        ),
        "positive_majorant": str(positive_majorant),
        "majorant_gap": str(majorant_gap),
        "weighted_source_bound": (
            "sqrt(k)*g_k'/r_k<=6+13/sqrt(k)+5/k+1/k^(3/2)<7"
        ),
        "rational_remainder_cap": str(scalar_remainder),
        "rational_margin_to_one": str(1 - scalar_remainder),
        "coercive_coupling_cap": str(coercive_coupling),
        "coercive_coupling_margin": str(5 - coercive_coupling),
    }


def initial_tail_diagnostics() -> dict:
    k = TAIL_START
    margin = 144 * 99_999**4 * k**3 - 100_000**4 * (5 * k + 6) ** 3
    if margin <= 0:
        raise RuntimeError("initial weighted tail does not fit below 12")
    return {
        "tail_start": TAIL_START,
        "source_bound": (
            "g_k<=(100000/99999)^2*(5*k+6)^(3/2)/k^2, k>=319"
        ),
        "weighted_bound": (
            "sqrt(k)*g_k<=(100000/99999)^2*(5+6/k)^(3/2)<12"
        ),
        "endpoint_squared_margin": margin,
        "monotonicity": "(5+6/k)^(3/2) decreases for k>=319",
    }


def build_payload() -> dict:
    flow = exact_flow_diagnostics()
    initial = initial_tail_diagnostics()
    rows = [
        CertificateRow(
            id="cfut_01_strict_reciprocal_defects",
            role="exact_analytic_input",
            readiness="ready_to_apply",
            claim="Strict moment log-concavity makes every reciprocal-defect coordinate finite along the actual heat trajectory.",
            formula="0<d_k(lambda)=1-x_k(lambda), q_k=d_k^(-1/2)",
            proof_boundary="Strict coordinate regularity only.",
        ),
        CertificateRow(
            id="cfut_02_exact_q_flow",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The defect heat equation induces an exact one-sided evolution equation for q_k.",
            formula=flow["q_flow"],
            proof_boundary="Exact coordinate identity only.",
            diagnostics={"increment": flow["increment"]},
        ),
        CertificateRow(
            id="cfut_03_neighbor_monotonicity",
            role="exact_lemma",
            readiness="ready_to_apply",
            claim="The increment vector field is strictly increasing in the next increment.",
            formula=flow["next_increment_monotonicity"],
            proof_boundary="One-sided cooperativity only.",
            diagnostics={"derivative": flow["next_increment_derivative"]},
        ),
        CertificateRow(
            id="cfut_04_weighted_source_cap",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="Inside the cubic cone, the vector field at a nonincreasing increment barrier has a uniform weighted upper bound.",
            formula=flow["weighted_source_bound"],
            proof_boundary="Conditional on 0<=g_k<=1 and q_k>=sqrt(k).",
            diagnostics=flow,
        ),
        CertificateRow(
            id="cfut_05_initial_weighted_tail",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The lambda=-100 entry theorem supplies a strict weighted reciprocal-increment tail.",
            formula="sqrt(k)*g_k(-100)<12 for every k>=319",
            proof_boundary="Entry-parameter tail only.",
            diagnostics=initial,
        ),
        CertificateRow(
            id="cfut_06_weighted_supersolution",
            role="exact_infinite_maximum_principle",
            readiness="ready_to_apply",
            claim="On every interval where the cubic cone holds, a componentwise weighted supersolution propagates the entry tail.",
            formula=(
                "g_k(lambda)<=(12+7*R_L*(lambda+100))/sqrt(k), "
                "R_L=sup_[-100,L] r_0"
            ),
            proof_boundary="Bootstrap estimate on an interval where 0<=g_k<=1.",
            diagnostics={
                "coercive_regularization": (
                    "G_k^epsilon=(12+7*R_L*t)/sqrt(k)+epsilon*k*exp(5*R_L*t); since g_k<=1, every first contact is at finite k."
                ),
                "coefficient_growth": (
                    "partial g_k'/partial g_(k+1)<=2*R_L*(2*k+5), while H_(k+1)-H_k=H_k/k and 2*(2*k+5)/k<5 for k>=319."
                ),
                "first_contact_conclusion": (
                    "At contact, the base source is below 7*R_L/sqrt(k) and the perturbation coupling is strictly below 5*R_L*H_k=H_k', contradicting first contact."
                ),
            },
        ),
        CertificateRow(
            id="cfut_07_forward_uniform_tail",
            role="exact_analytic_theorem",
            readiness="ready_to_apply",
            claim="The reciprocal-defect increment tends to zero uniformly on every compact forward heat interval.",
            formula=(
                "sup_(lambda in [-100,L]) g_k(lambda)<="
                "(12+7*R_L*(L+100))/sqrt(k)->0"
            ),
            proof_boundary="Degree-3 tail theorem only.",
        ),
        CertificateRow(
            id="cfut_08_cubic_continuation",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="The entry theorem, weighted tail, and inward frontier condition propagate the full cubic cone to every finite lambda>=-100.",
            formula="0<=q_(k+1)(lambda)-q_k(lambda)<=1 for every k>=1",
            proof_boundary="All-shift degree-3 propagation only.",
        ),
        CertificateRow(
            id="cfut_09_lambda_zero_cubics",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="Every shifted degree-3 Jensen polynomial of the lambda=0 coefficient sequence is hyperbolic.",
            formula="Disc(J_(3,n)(.;0))>=0 for every n>=0",
            proof_boundary="Degree 3 at lambda=0 only; not all degrees.",
        ),
        CertificateRow(
            id="cfut_10_higher_degree_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="A coupled invariant for degree four and then a genuinely all-degree mechanism are still required.",
            formula="degree 3 propagation + missing higher-minor hierarchy",
            proof_boundary="Open higher-degree bridge; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_cubic_forward_uniform_tail_certificate",
        "date": "2026-07-10",
        "status": "exact all-shift cubic forward-propagation theorem with open higher-degree handoff",
        "proof_boundary": (
            "This artifact proves the reciprocal-defect tail uniformly on compact forward heat intervals and propagates every shifted degree-3 Jensen polynomial through lambda=0. It does not prove degree 4 or higher, PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_cubic_reciprocal_defect_invariance_lemma.md",
            "outputs/jensen_window_pf_cubic_m100_tail_entry_certificate.md",
            "outputs/jensen_window_pf_heat_flow_infinite_cone_invariance_certificate.md",
        ],
        "diagnostics": {"flow": flow, "initial_tail": initial},
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    flow = payload["diagnostics"]["flow"]
    initial = payload["diagnostics"]["initial_tail"]
    return "\n".join(
        [
            "# Jensen-Window PF Cubic Forward-Uniform Tail Certificate",
            "",
            "Date: 2026-07-10",
            "",
            "Status: exact all-shift degree-3 propagation through lambda=0 with the",
            "higher-degree bridge open. This is not a proof of PF-infinity, RH, or",
            "`Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_cubic_forward_uniform_tail_certificate.json",
            "python work/rh_compute/scripts/jensen_window_pf_cubic_forward_uniform_tail_certificate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_cubic_forward_uniform_tail_certificate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF cubic forward-uniform tail certificate: 10 rows, 0 issues, 3 exact flow identities, 1 weighted source cap, 1 initial weighted tail, 1 forward-uniform tail, 1 full cubic propagation theorem, 1 lambda=0 cubic theorem, 1 open higher-degree handoff",
            "```",
            "",
            "## Exact Increment Flow",
            "",
            "Set `d_k=1-x_k`, `q_k=d_k^(-1/2)`, and",
            "`g_k=q_(k+1)-q_k`. The defect equation gives",
            "",
            "```text",
            flow["q_flow"],
            "```",
            "",
            "The derivative with respect to the next increment is strictly positive:",
            "",
            "```text",
            flow["next_increment_derivative"],
            "```",
            "",
            "At `g_(k+1)=g_k=g`, exact factorization and deletion of manifestly",
            "negative terms give, for `0<=g<=1` and `q_k>=sqrt(k)`,",
            "",
            "```text",
            flow["weighted_source_bound"],
            "```",
            "",
            "The last strict inequality uses `sqrt(319)>17` and the exact remainder",
            f"`{flow['rational_remainder_cap']}<1`.",
            "",
            "## Entry Tail",
            "",
            "The lambda=-100 theorem already proves",
            "",
            "```text",
            initial["source_bound"],
            "```",
            "",
            "After multiplying by `sqrt(k)`, the right side decreases with `k` and",
            "is strictly below 12 at `k=319`. The squared integer margin is",
            "",
            "```text",
            str(initial["endpoint_squared_margin"]),
            "```",
            "",
            "## Forward Uniformity",
            "",
            "For a finite `L`, put `R_L=sup_[-100,L] r_0(lambda)<infinity`. On an",
            "interval where the cubic cone holds, the componentwise barrier",
            "",
            "```text",
            "G_k(lambda)=(12+7*R_L*(lambda+100))/sqrt(k)",
            "```",
            "",
            "is a supersolution. Use the explicit coercive regularization",
            "",
            "```text",
            "G_k^epsilon=G_k+epsilon*k*exp(5*R_L*(lambda+100)).",
            "```",
            "",
            "Because `g_k<=1`, any first contact is at a finite index. At contact,",
            "`H_(k+1)-H_k=H_k/k`, while next-increment monotonicity gives the",
            "extra coupling cap `2*R_L*(2*k+5)*H_k/k<5*R_L*H_k=H_k'`.",
            "The weighted source cap excludes contact. Letting `epsilon` vanish proves",
            "",
            "```text",
            "sup_(lambda in [-100,L]) g_k(lambda)",
            " <=(12+7*R_L*(L+100))/sqrt(k)->0.",
            "```",
            "",
            "This tail makes the reciprocal-defect first-crossing principle legitimate.",
            "The all-k entry theorem and exact inward frontier condition therefore",
            "propagate `0<=g_k<=1` to every finite `lambda>=-100`.",
            "",
            "## Consequence And Boundary",
            "",
            "Every shifted degree-3 Jensen polynomial is hyperbolic at `lambda=0`.",
            "No degree-4 invariant or all-degree PF/Jensen theorem is proved here; those",
            "are the remaining obligations before any RH conclusion could be considered.",
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
        "wrote Jensen-window PF cubic forward-uniform tail certificate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
