#!/usr/bin/env python3
"""Build the first signed correction for the critical Dirichlet blocks."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.md"
)
POLYMATH_SOURCE = "https://arxiv.org/abs/1904.12438"
DLMF_SOURCE = "https://dlmf.nist.gov/5.11"
ENDPOINT_NOTE = (
    "outputs/"
    "jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.md"
)
CONTACT_NOTE = (
    "outputs/"
    "jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate.md"
)


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def build_exact() -> dict:
    time = sp.symbols("t", positive=True, real=True)
    v = sp.symbols("v", real=True)
    alpha_n, alpha_first, s = sp.symbols(
        "alpha_n alpha_first s", complex=True
    )
    w = sp.sqrt(time) * v + time * alpha_n / 2

    leading_exponent = sp.simplify(
        -time * alpha_n**2 / 4
        - sp.sqrt(time) * v * alpha_n
        + alpha_n * w
    )
    if sp.simplify(leading_exponent - time * alpha_n**2 / 4) != 0:
        raise RuntimeError("shifted Dirichlet leading exponent failed")

    w_squared = sp.Poly(sp.expand(w**2), v)
    gaussian_w2 = sp.simplify(
        w_squared.coeff_monomial(v**0)
        + sp.Rational(1, 2) * w_squared.coeff_monomial(v**2)
    )
    expected_w2 = time / 2 + time**2 * alpha_n**2 / 4
    if sp.simplify(gaussian_w2 - expected_w2) != 0:
        raise RuntimeError("shifted Gaussian second moment failed")

    first_correction = sp.simplify(
        1 / (6 * s) + alpha_first * gaussian_w2 / 2
    )
    expected_correction = (
        1 / (6 * s)
        + alpha_first * (time / 4 + time**2 * alpha_n**2 / 8)
    )
    if sp.simplify(first_correction - expected_correction) != 0:
        raise RuntimeError("first Dirichlet correction failed")

    return {
        "heat_shift": (
            "r_(t,n)(s)=exp(-t*alpha_n^2/4)*integral "
            "exp(-sqrt(t)*v*alpha_n)r_(0,n)(s+w(v))dmu(v), "
            "w(v)=sqrt(t)*v+t*alpha_n/2"
        ),
        "definitions": (
            "alpha_n=alpha(s)-log(n), "
            "dmu(v)=pi^(-1/2)exp(-v^2)dv"
        ),
        "scaled_gamma_factorization": (
            "r_(0,n)(s)=M_0(s)*n^(-s)*GammaStar(s/2), "
            "GammaStar(z)=Gamma(z)/(exp(-z)z^z*sqrt(2*pi/z))"
        ),
        "leading_main": (
            "m_(t,n)(s)=M_0(s)*n^(-s)*exp(t*alpha_n^2/4)"
        ),
        "exact_relative_ratio": (
            "r_(t,n)(s)/m_(t,n)(s)=integral GammaStar((s+w)/2)*"
            "exp(log(M_0(s+w))-log(M_0(s))-alpha(s)w)dmu(v)"
        ),
        "leading_cancellation": (
            "-t*alpha_n^2/4-sqrt(t)*v*alpha_n+alpha_n*w"
            "=t*alpha_n^2/4"
        ),
        "gaussian_moment": (
            "integral w(v)^2 dmu(v)=t/2+t^2*alpha_n^2/4"
        ),
        "first_correction": (
            "d_(t,n)(s)=1/(6s)+alpha'(s)"
            "*(t/4+t^2*alpha_n^2/8)"
        ),
        "refined_term": (
            "m_[1](t,n;s)=m_(t,n)(s)*(1+d_(t,n)(s))"
        ),
        "alpha_derivatives": (
            "alpha'(s)=-1/(2s^2)-1/(s-1)^2+1/(2s); "
            "alpha''(s)=1/s^3+2/(s-1)^3-1/(2s^2)"
        ),
        "log_gamma_remainder": (
            "log GammaStar(z)=1/(12z)+R_G(z), "
            "|R_G(z)|<=sec(arg(z)/2)^4/(360*|z|^3)"
        ),
        "upper_half_plane_gamma_bound": (
            "For z=(s+w)/2 with Im(s+w)>0, "
            "|R_G(z)|<=16*|s+w|/(45*Im(s+w)^4)"
        ),
        "m0_taylor_remainder": (
            "R_M=log(M_0(s+w))-log(M_0(s))-alpha(s)w"
            "-alpha'(s)w^2/2, "
            "|R_M|<=|w|^3*sup_(0<=u<=1)|alpha''(s+u*w)|/6"
        ),
        "log_residual": (
            "rho(v)=R_M+R_G+1/(6(s+w))-1/(6s)"
        ),
        "central_tail_split": (
            "Split |v|<=V and |v|>V with V=T^(1/3), T=Im(s): "
            "use the second-order expansions on the central part and the "
            "published Proposition 6.1 quadratic Gaussian envelope on the tail"
        ),
        "per_term_target": (
            "|r_(t,n)/m_(t,n)-(1+d_(t,n))|<=C_D/T^2 "
            "uniformly for n<=N on the critical radius-1/L collar"
        ),
        "finite_sum_scale": (
            "sum of the two absolute finite-main coefficient masses "
            "<=50*exp(L/4), T>=2*pi*exp(L-1)"
        ),
        "dirichlet_residual_scale": (
            "C_D/T^2 times 50*exp(L/4)=O(C_D*exp(-7L/4))"
        ),
        "combined_signed_update": (
            "DeltaJ=(DeltaA+DeltaB-DeltaC)/A_t, "
            "J_[1]=J_[0]+DeltaJ, r_[0]=DeltaJ+r_[1]"
        ),
        "global_scale_handoff": (
            "The peeled Dirichlet residual is O(exp(-7L/4)); the RS-C_1 "
            "endpoint and cutoff residuals are O(exp(-5L/4)) and therefore "
            "set the next global C^1 scale"
        ),
        "uniform_target": (
            "Prove explicit constants eta_0,eta_1 with "
            "|r_[1]|<=eta_0*exp(-5L/4), "
            "|r_[1],x|<=eta_1*L*exp(-5L/4) on L>=50, 0<tL<=25, "
            "including adjacent cutoffs, then use r_0=DeltaJ in the signed "
            "contact-normal inequality"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15cdfc_01_heat_shift",
            role="published_primary_input",
            readiness="available_published",
            claim=(
                "Polymath 15 gives an exact contour-shifted Gaussian formula "
                "for every evolved Dirichlet block."
            ),
            formula=exact["heat_shift"] + "; " + exact["definitions"],
            proof_boundary="Equation (40), with its stated contour condition.",
        ),
        GateRow(
            id="np15cdfc_02_scaled_gamma",
            role="published_primary_input",
            readiness="available_published",
            claim=(
                "The exact gamma prefactor is the Polymath normalizer times "
                "the standard scaled gamma function."
            ),
            formula=exact["scaled_gamma_factorization"],
            proof_boundary=(
                "Algebraic comparison of the exact xi gamma factor with M_0; "
                "DLMF supplies the signed scaled-gamma expansion."
            ),
        ),
        GateRow(
            id="np15cdfc_03_exact_relative_ratio",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "After the saddle shift, every discarded finite-sum factor is "
                "contained in one exact centered Gaussian expectation."
            ),
            formula=(
                exact["leading_main"]
                + "; "
                + exact["exact_relative_ratio"]
                + "; "
                + exact["leading_cancellation"]
            ),
            proof_boundary="No asymptotic truncation in the displayed ratio.",
        ),
        GateRow(
            id="np15cdfc_04_gaussian_moment",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The shifted Gaussian second moment is explicit.",
            formula=exact["gaussian_moment"],
            proof_boundary="Uses integral v dmu=0 and integral v^2 dmu=1/2.",
        ),
        GateRow(
            id="np15cdfc_05_first_correction",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The complete first signed relative correction combines the "
                "scaled-gamma term with the quadratic log-M_0 term."
            ),
            formula=(
                exact["first_correction"]
                + "; "
                + exact["refined_term"]
                + "; "
                + exact["alpha_derivatives"]
            ),
            proof_boundary=(
                "Exact coefficient extraction; the size of the remainder "
                "after truncation is a separate obligation."
            ),
        ),
        GateRow(
            id="np15cdfc_06_remainder_decomposition",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The unretained logarithmic error splits into an explicit "
                "M_0 Taylor remainder, scaled-gamma remainder, and rational "
                "shift difference."
            ),
            formula=(
                exact["m0_taylor_remainder"]
                + "; "
                + exact["log_residual"]
            ),
            proof_boundary="Exact Taylor remainder identity on the chosen branch.",
        ),
        GateRow(
            id="np15cdfc_07_gamma_guard",
            role="analytic_guard",
            readiness="conditional_ready",
            claim=(
                "A published complex-sector log-gamma remainder reduces to an "
                "integrable upper-half-plane bound."
            ),
            formula=(
                exact["log_gamma_remainder"]
                + "; "
                + exact["upper_half_plane_gamma_bound"]
            ),
            proof_boundary=(
                "The DLMF bound is published; its uniform insertion into the "
                "full heat collar still needs explicit branch and denominator "
                "bookkeeping."
            ),
        ),
        GateRow(
            id="np15cdfc_08_central_tail_guard",
            role="analytic_guard",
            readiness="conditional_ready",
            claim=(
                "A central/tail split prevents the cubic Taylor remainder from "
                "being exponentiated on the entire Gaussian line."
            ),
            formula=(
                exact["central_tail_split"]
                + "; "
                + exact["per_term_target"]
                + "; "
                + exact["finite_sum_scale"]
                + "; "
                + exact["dirichlet_residual_scale"]
            ),
            proof_boundary=(
                "Architecture only. The finite constant C_D and collar "
                "monotonicity inequalities are not certified here."
            ),
        ),
        GateRow(
            id="np15cdfc_09_signed_main_update",
            role="exact_definition",
            readiness="ready_to_apply",
            claim=(
                "The two Dirichlet corrections and the endpoint correction "
                "combine into the signed remainder model required at contact."
            ),
            formula=(
                exact["combined_signed_update"]
                + "; "
                + exact["global_scale_handoff"]
            ),
            proof_boundary=(
                "Exact repartition once DeltaA, DeltaB, and DeltaC are defined; "
                "the claimed scale remains conditional on the open bounds."
            ),
        ),
        GateRow(
            id="np15cdfc_10_uniform_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim=(
                "The first-order asymptotic programme reduces to one explicit "
                "global second-order value/first-jet estimate."
            ),
            formula=exact["uniform_target"],
            proof_boundary=(
                "Open. No strict signed contact inequality, contact exclusion, "
                "Lambda<=0, or RH is asserted."
            ),
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_newman_polymath15_critical_"
            "dirichlet_first_correction_gate"
        ),
        "date": "2026-07-26",
        "status": (
            "exact first signed correction for the critical Dirichlet blocks "
            "with a guarded second-order remainder programme"
        ),
        "proof_boundary": (
            "This artifact proves the shifted-ratio identity, Gaussian moment, "
            "first correction formula, and logarithmic remainder decomposition. "
            "It does not certify the uniform C_D/T^2 remainder, the resulting "
            "global exp(-5L/4) C^1 constants, the cutoff collar, the strict "
            "signed contact inequality, contact exclusion, Lambda<=0, or RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_identities": 4,
            "primary_source_inputs": 2,
            "analytic_guards": 2,
            "signed_main_definitions": 1,
            "open_uniform_targets": 1,
        },
        "sources": [
            POLYMATH_SOURCE,
            DLMF_SOURCE,
            ENDPOINT_NOTE,
            CONTACT_NOTE,
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_cell_remainder_certificate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical Dirichlet First Correction Gate",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact first signed finite-sum correction and guarded",
            "second-order remainder programme. This is not a proof of",
            "`Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_dirichlet_first_correction_gate.py",
            "```",
            "",
            "## Exact Centered Ratio",
            "",
            f"[Polymath 15]({POLYMATH_SOURCE}), equation (40), gives",
            "",
            "```text",
            exact["definitions"],
            exact["heat_shift"],
            exact["scaled_gamma_factorization"],
            "```",
            "",
            "Substitution and exact cancellation yield",
            "",
            "```text",
            exact["leading_main"],
            exact["leading_cancellation"],
            exact["exact_relative_ratio"],
            "```",
            "",
            "This ratio contains every term discarded by Proposition 6.1.",
            "",
            "## First Signed Term",
            "",
            "The centered Gaussian moment is",
            "",
            "```text",
            exact["gaussian_moment"],
            "```",
            "",
            "The scaled-gamma coefficient `1/(6s)` and the quadratic",
            "`log M_0` term therefore give",
            "",
            "```text",
            exact["first_correction"],
            exact["refined_term"],
            "```",
            "",
            "This is signed and componentwise. It is the finite-sum analogue",
            "of retaining the first omitted Riemann-Siegel endpoint coefficient.",
            "",
            "## Remainder Anatomy",
            "",
            f"The [DLMF log-gamma expansion]({DLMF_SOURCE}) gives",
            "",
            "```text",
            exact["log_gamma_remainder"],
            exact["upper_half_plane_gamma_bound"],
            exact["alpha_derivatives"],
            exact["m0_taylor_remainder"],
            exact["log_residual"],
            "```",
            "",
            "A global cubic Taylor bound cannot simply be exponentiated under",
            "the whole Gaussian. The correct proof split is",
            "",
            "```text",
            exact["central_tail_split"],
            "```",
            "",
            "The central range uses the signed second-order expansion; the tail",
            "reuses the published quadratic Gaussian envelope.",
            "",
            "## Scale Forecast",
            "",
            "The theorem to certify per component is",
            "",
            "```text",
            exact["per_term_target"],
            exact["finite_sum_scale"],
            exact["dirichlet_residual_scale"],
            "```",
            "",
            "Thus the peeled Dirichlet sums should be smaller than the new",
            "endpoint remainder:",
            "",
            "```text",
            exact["global_scale_handoff"],
            "```",
            "",
            "The exponent comparison is conditional until an explicit `C_D`",
            "is proved.",
            "",
            "## Signed Contact Handoff",
            "",
            "```text",
            exact["combined_signed_update"],
            exact["uniform_target"],
            "```",
            "",
            "This is the next quantitative job: certify the central Gaussian",
            "moments, tail, first derivative, and cutoff collars with explicit",
            "constants. Only the residual then pays an absolute-value budget.",
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
        "built Newman Polymath-15 critical Dirichlet first correction gate: "
        f"{len(artifact['rows'])} rows, 4 exact identities, "
        "2 primary-source inputs, 2 analytic guards, "
        "1 signed-main definition, 1 open uniform target"
    )


if __name__ == "__main__":
    main()
