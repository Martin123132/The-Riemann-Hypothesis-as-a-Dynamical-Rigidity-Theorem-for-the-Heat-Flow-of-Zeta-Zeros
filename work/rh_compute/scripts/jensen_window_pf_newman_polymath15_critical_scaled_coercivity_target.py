#!/usr/bin/env python3
"""Build the corrected finite coercivity target in the critical scaled layer."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.md"
)
SOURCE_URL = "https://arxiv.org/abs/1904.12438"
OSCILLATORY_HANDOFF_NOTE = (
    "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md"
)
C_STAR_EXACT = "4911678521/1933561194"
C_STAR_DECIMAL = "2.540223984760008"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    x, y, u, v, w, phase_first, phase_second, potential = sp.symbols(
        "x y u v w phase_first phase_second potential", real=True
    )
    value = 2 * x
    first = 2 * (u - phase_first * y)
    second = 2 * (
        w
        - 2 * phase_first * v
        - phase_second * y
        - phase_first**2 * x
    )
    curvature = sp.expand(first**2 - value * second + potential * value**2)
    coercivity = 4 * (
        phase_first**2 * (x**2 + y**2)
        + 2 * phase_first * (x * v - u * y)
        + u**2
        - x * w
        + phase_second * x * y
        + potential * x**2
    )
    if sp.simplify(curvature - coercivity) != 0:
        raise RuntimeError("critical finite coercivity identity failed")

    p0, p1, p2, q0, q1, q2 = sp.symbols(
        "p0 p1 p2 q0 q1 q2", real=True
    )
    c = lambda f0, f1, f2: f1**2 - f0 * f2 + potential * f0**2
    corrected_delta = sp.expand(c(p0 - q0, p1 - q1, p2 - q2) - c(p0, p1, p2))
    expected_delta = (
        -2 * p1 * q1
        + q1**2
        + p0 * q2
        + q0 * p2
        - q0 * q2
        + potential * (-2 * p0 * q0 + q0**2)
    )
    if sp.simplify(corrected_delta - expected_delta) != 0:
        raise RuntimeError("Riemann-Siegel correction curvature identity failed")

    return {
        "critical_region": (
            f"L=log(x/(4*pi))->infinity, c=t*L in [0,c_*+o(1)], "
            f"c_*={C_STAR_EXACT}={C_STAR_DECIMAL}..., "
            "N=floor(sqrt(x/(4*pi)+t/16))"
        ),
        "dirichlet_main": (
            "D_(N,t)(x)=sum_(n<=N)exp((t/4)log(n)^2)*n^(-s_*(x)), "
            "P_(N,t)=2Re(exp(i*beta_t)D_(N,t))"
        ),
        "components": (
            "exp(i*beta)D=X+iY, exp(i*beta)D'=U+iV1, "
            "exp(i*beta)D''=W+iV2, B=beta'"
        ),
        "coercivity": (
            "C[P]/4=B^2(X^2+Y^2)+2B(X*V1-U*Y)+U^2-X*W+"
            "beta''*X*Y+V_t*X^2"
        ),
        "invariant_cross": "X*V1-U*Y=Im(conj(D)*D')",
        "correction_parameters": (
            "T'=x/2+pi*t/8, a=sqrt(T'/(2*pi)), N=floor(a), "
            "p=1-2(a-N), U_RS=exp(-i((T'/2)log(T'/(2*pi))-T'/2-pi/8))"
        ),
        "C0": (
            "C0(p)=(exp(pi*i*(p^2/2+3/8))-i*sqrt(2)*cos(pi*p/2))/"
            "(2*cos(pi*p)), with removable values at p=+-1/2"
        ),
        "raw_correction": (
            "C_t(x)=2*(-1)^N*exp(t*pi^2/64)*"
            "Re(M_0(i*T')*C0(p)*U_RS*exp(pi*i/8))"
        ),
        "corrected_main": (
            "Q_(N,t)=C_t/A_t, J_(N,t)=P_(N,t)-Q_(N,t), A_t=|B_t|"
        ),
        "corrected_curvature": (
            "C[J]-C[P]=-2P'Q'+Q'^2+P*Q''+Q*P''-Q*Q''+"
            "V_t*(-2P*Q+Q^2)"
        ),
        "refined_remainder": (
            "H_t/A_t=J_(N,t)+r_ref; the published e_C replaces the e_C0 "
            "factor 1 by tilde_epsilon(s_-)+tilde_epsilon(s_+)=O(1/N)"
        ),
        "shrinking_collar": (
            "On radius rho=1/L collars, the refined normalized C2 remainder is "
            "o(L^2) after including one adjacent-N correction"
        ),
        "live_target": (
            "Prove C_t[J_(N,t)] greater than its explicit refined C2 remainder "
            "budget uniformly for 0<=t*L<=c_*+o(1)"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15csct_01_critical_region",
            role="exact_reduction",
            readiness="ready_to_apply",
            claim="The remaining asymptotic obstruction is a bounded scaled-time finite Riemann-Siegel problem.",
            formula=exact["critical_region"],
            proof_boundary=(
                "The c>c_*+epsilon region is closed by the oscillatory zeta "
                "handoff theorem."
            ),
        ),
        GateRow(
            id="np15csct_02_uncorrected_coercivity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The finite-main Laguerre curvature has an exact component formula valid even when D vanishes.",
            formula=exact["coercivity"],
            proof_boundary="No division by D and no zero-free assumption.",
        ),
        GateRow(
            id="np15csct_03_phase_free_cross",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The leading mixed phase term is invariant under the normalizer rotation.",
            formula=exact["invariant_cross"],
            proof_boundary="Exact complex multiplication identity.",
        ),
        GateRow(
            id="np15csct_04_rs_correction",
            role="published_explicit_input",
            readiness="available_published",
            claim="Polymath 15 supplies an explicit leading Riemann-Siegel endpoint correction.",
            formula=exact["raw_correction"],
            proof_boundary="Imported from the published A+B-C approximation.",
        ),
        GateRow(
            id="np15csct_05_corrected_main",
            role="exact_definition",
            readiness="ready_to_apply",
            claim="The correct critical finite main object includes the normalized endpoint correction.",
            formula=exact["corrected_main"],
            proof_boundary="Omitting Q loses one asymptotic order in the critical layer.",
        ),
        GateRow(
            id="np15csct_06_correction_curvature",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The endpoint correction enters the finite curvature through an explicit six-term perturbation.",
            formula=exact["corrected_curvature"],
            proof_boundary="Exact symbolic identity.",
        ),
        GateRow(
            id="np15csct_07_refined_remainder",
            role="published_theorem_transfer",
            readiness="ready_to_apply",
            claim="After extracting Q, the published contour remainder gains a full factor of the saddle index.",
            formula=exact["refined_remainder"],
            proof_boundary="Uses the published tilde-epsilon bound; derivative constants remain to be made explicit.",
        ),
        GateRow(
            id="np15csct_08_shrinking_collar",
            role="analytic_handoff",
            readiness="conditional_ready",
            claim="A scale-adapted Cauchy collar keeps normalizer growth bounded while differentiating the refined remainder.",
            formula=exact["shrinking_collar"],
            proof_boundary="Asymptotic transfer schema; an explicit finite threshold is not supplied here.",
        ),
        GateRow(
            id="np15csct_09_live_coercivity",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="The global tail is reduced to one corrected finite oscillatory coercivity inequality.",
            formula=exact["live_target"],
            proof_boundary="Open and RH-level; no sign is asserted.",
        ),
        GateRow(
            id="np15csct_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim=(
                "Absolute coefficient bounds stop at c=4; the exponent-pair "
                "handoff lowers the open boundary to c_*, but neither argument "
                "certifies the corrected coercivity at or below c_*."
            ),
            formula="worst coefficient exponent=1/2+c/8<=1",
            proof_boundary="Cancellation, not additional truncation depth, is required.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target",
        "date": "2026-07-17",
        "status": (
            "exact corrected finite coercivity target for the critical scaled "
            "Newman layer; not a proof of its sign, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the uncorrected component identity and the exact "
            "curvature contribution of the published Riemann-Siegel correction. "
            "It identifies the refined remainder order and the live finite target. "
            "It does not prove the corrected coercivity inequality, positivity "
            "throughout the critical layer, positivity for all x and t, Lambda<=0, "
            "or RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            SOURCE_URL,
            OSCILLATORY_HANDOFF_NOTE,
            "outputs/jensen_window_pf_newman_polymath15_scaled_boundary_layer_asymptotic_theorem.md",
            "outputs/jensen_window_pf_newman_polymath15_dominant_saddle_global_ray_certificate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical Scaled Coercivity Target",
            "",
            "Date: 2026-07-17",
            "",
            "Status: exact corrected finite target for the remaining critical scaled",
            "layer. This is not a proof of its sign, `Lambda <= 0`, or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.py",
            "```",
            "",
            "The finite approximation and endpoint correction are imported from",
            f"D. H. J. Polymath, [Effective approximation of heat flow evolution of the Riemann xi function]({SOURCE_URL}),",
            "especially the A+B-C approximation in Corollary 6.4.",
            "",
            "## Critical Region",
            "",
            "```text",
            exact["critical_region"],
            "```",
            "",
            f"The oscillatory zeta handoff closes every fixed `c > c_*`, where",
            f"`c_* = {C_STAR_EXACT} = {C_STAR_DECIMAL}...`. Thus only",
            "`0 <= c <= c_* + o(1)` remains asymptotically open. The earlier",
            "absolute zeta-moment argument stops at `c=4`; published exponent-pair",
            "cancellation is what lowers the boundary from `4` to `c_*`.",
            "",
            "## Finite Coercivity",
            "",
            "Set",
            "",
            "```text",
            exact["dirichlet_main"],
            exact["components"],
            "```",
            "",
            "Direct symbolic differentiation, without dividing by `D`, gives",
            "",
            "```text",
            exact["coercivity"],
            exact["invariant_cross"],
            "```",
            "",
            "This remains valid at zeros of the complex finite sum. It exposes the",
            "precise competition between the `beta'^2*|D|^2` phase floor, the",
            "symplectic cross term, and the first two Dirichlet jets.",
            "",
            "## Endpoint Correction",
            "",
            "The uncorrected A+B theorem treats the leading Riemann-Siegel endpoint",
            "term as error. That is too coarse in the critical layer. Define",
            "",
            "```text",
            exact["correction_parameters"],
            exact["C0"],
            exact["raw_correction"],
            exact["corrected_main"],
            "```",
            "",
            "Its exact curvature contribution is",
            "",
            "```text",
            exact["corrected_curvature"],
            "```",
            "",
            "## Refined Transfer",
            "",
            "After extracting the correction, the paper replaces the `e_C0` factor",
            "`1+tilde_epsilon(s_-)+tilde_epsilon(s_+)` by",
            "`tilde_epsilon(s_-)+tilde_epsilon(s_+)`. Since the latter is `O(1/N)`,",
            "the refined remainder gains one full saddle-index factor:",
            "",
            "```text",
            exact["refined_remainder"],
            exact["shrinking_collar"],
            "```",
            "",
            "A radius `1/L` keeps the normalizer ratio bounded, while Cauchy costs",
            "only the powers of `L` naturally present in the curvature. The refined",
            "remainder is then lower order even under crude absolute main-jet caps.",
            "",
            "## Live Target",
            "",
            "The remaining theorem obligation is now concrete:",
            "",
            "```text",
            exact["live_target"],
            "```",
            "",
            "This is an oscillatory arithmetic inequality for a finite sum of",
            "length about `sqrt(x/(4*pi))`, with an explicit endpoint correction",
            "and explicit error. Absolute tails cannot settle it; proving its sign",
            "would be the genuinely new RH-level step.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman Polymath-15 critical scaled coercivity target: "
        f"{len(artifact['rows'])} rows, 2 exact curvature identities, 1 live target"
    )


if __name__ == "__main__":
    main()
