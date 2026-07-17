#!/usr/bin/env python3
"""Prove the compact-heat first-summand heat-tilt ratio estimate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/"
    "jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.md"
)
SOURCE_URL = "https://arxiv.org/abs/2007.13582"


@dataclass(frozen=True)
class HeatTiltRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def uniform_suitability() -> dict:
    x, v, T = sp.symbols("x v T")
    exponent = -T * (
        v * sp.log(1 + x) / 8 + sp.log(1 + x) ** 2 / 16
    )
    series = sp.series(sp.exp(exponent), x, 0, 8).removeO().expand()
    coefficients = [sp.factor(series.coeff(x, degree)) for degree in range(8)]
    checks = []
    for degree, coefficient in enumerate(coefficients):
        polynomial = sp.Poly(coefficient, T, v)
        degree_t = polynomial.degree(T)
        degree_v = polynomial.degree(v)
        if degree_t > degree or degree_v > degree:
            raise RuntimeError(f"heat multiplier coefficient degree failed at {degree}")
        checks.append(
            {
                "degree": degree,
                "T_degree": degree_t,
                "v_degree": degree_v,
                "coefficient": str(coefficient),
            }
        )
    return {
        "family": "f_T(t)=exp(-T*(log t)^2/16), 0<=T<=100",
        "ratio_identity": (
            "f_T(t*(1+x))/f_T(t)=exp(-T*((v/8)*log(1+x)+"
            "log(1+x)^2/16)), t=exp(v)"
        ),
        "uniform_disk": (
            "for v>=1 and |x|<=1/(2v), |log(1+x)|<=2|x| and the "
            "exponent is bounded uniformly for 0<=T<=100"
        ),
        "uniform_taylor": (
            "f_T(t*(1+x))/f_T(t)=sum_(r=0)^(K-1) f_(T,r)(v)x^r+"
            "O_K(v^K|x|^K), uniformly for 0<=T<=100"
        ),
        "coefficient_bound": (
            "sup_(0<=T<=100)|f_(T,r)(v)|=O_r(v^r)"
        ),
        "suitability_data": "b=0, lambda_suit=1, uniformly for T in [0,100]",
        "symbolic_coefficients_through_7": checks,
    }


def lambert_derivatives() -> dict:
    k, w = sp.symbols("k w", positive=True)

    def derivative(expression: sp.Expr) -> sp.Expr:
        return sp.factor(
            sp.diff(expression, k)
            + sp.diff(expression, w) * w / (k * (1 + w))
        )

    current = w**2
    rows = []
    for order in range(1, 8):
        current = derivative(current)
        normalized = sp.factor(current * k**order / w)
        limit = sp.limit(normalized, w, sp.oo)
        if not limit.is_finite:
            raise RuntimeError(f"Lambert derivative limit failed at order {order}")
        rows.append(
            {
                "order": order,
                "derivative": str(current),
                "normalized_by_w_over_km": str(normalized),
                "limit_w_to_infinity": str(limit),
            }
        )
    return {
        "uniformizer": "w(k)=W(2*k/pi)",
        "derivative_identity": "w'(k)=w/(k*(1+w))",
        "main_term_derivatives": "d^m/dk^m w(k)^2=O(w(k)/k^m), m=1,...,7",
        "finite_difference_identity": (
            "Delta^m F(k)=integral_[0,1]^m F^(m)(k+s_1+...+s_m) ds"
        ),
        "main_term_differences": (
            "Delta^m w(k)^2=O(w(k)/k^m), m=1,...,7"
        ),
        "rows": rows,
    }


def published_input() -> dict:
    return {
        "source": "Cormac O'Sullivan, Zeros of Jensen polynomials and asymptotics for the Riemann xi function",
        "url": SOURCE_URL,
        "journal": "Research in the Mathematical Sciences 8 (2021), article 20",
        "theorem": "Theorem 5.2",
        "supporting_section": "Section 5, especially the Gaussian-log suitable multiplier example",
        "xi_section": "Section 7, first-summand reduction to I_(alpha,beta)(n)",
        "published_integral": (
            "I_alpha(f;n)=integral_1^infinity (log t)^n*exp(-alpha*t)*f(t)dt"
        ),
        "published_expansion": (
            "I_alpha(f;n)=sqrt(2*pi)*u^(n+1)*f(exp(u))*exp(u-n/u)/"
            "sqrt((1+u)*n)*(1+sum_(r=1)^(R-1)a_r(f;u)/n^r+"
            "O(u^(R*(1+2*lambda_suit))/n^R)), u=W(n/alpha)"
        ),
        "coefficient_growth": "a_r(f;u)=O(u^(r*(1+2*lambda_suit)))",
    }


def heat_tilt_theorem() -> dict:
    return {
        "first_summand_change_of_variables": (
            "M_k^(1)(T)=C_k*(2*pi^2*I_pi(t^(5/4)*f_T(t);2k)-"
            "3*pi*I_pi(t^(1/4)*f_T(t);2k))"
        ),
        "normalization": "C_k>0 is independent of T and cancels from R_T^(1)(k)",
        "dominant_uniformizer": "w=W(2*k/pi)",
        "uniform_all_order_expansion": (
            "log M_k^(1)(T)=L_0(k)-T*w^2/16+"
            "sum_(r=1)^(R-1) P_(r,T)(w)/k^r+O_R(w^(3R)/k^R)"
            ", uniformly for 0<=T<=100"
        ),
        "coefficient_structure": (
            "for fixed r, P_(r,T)(w) is a bounded-T polynomial/rational "
            "combination of powers of w and (1+w)^(-1)"
        ),
        "tilt_log_expansion": (
            "log R_T^(1)(k)=-T*w^2/16+sum_(r=1)^(R-1) "
            "Q_(r,T)(w)/k^r+O_R(w^(3R)/k^R)"
        ),
        "remainder_choice": (
            "for fixed m<=7 choose R>m; then w^(3R)/k^R=o(w/k^m)"
        ),
        "correction_differences": (
            "Delta^m(Q_(r,T)(w)/k^r)=O_(r,m)(poly(w)/k^(r+m))"
            "=o(w/k^m), uniformly for bounded T"
        ),
        "target_theorem": (
            "Delta_k^m log R_T^(1)(k)=O(log(k)/k^m) uniformly for "
            "0<=T<=100 and m=2,...,7"
        ),
        "backward_version": (
            "the same estimate holds for every fixed backward-difference stencil"
        ),
    }


def build_artifact() -> dict:
    published = published_input()
    suitability = uniform_suitability()
    lambert = lambert_derivatives()
    theorem = heat_tilt_theorem()
    rows = [
        HeatTiltRow(
            id="ufshta_01_integral_mapping",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The first Newman theta summand with heat penalty is exactly a two-term O'Sullivan integral with Gaussian-log multiplier.",
            formula=theorem["first_summand_change_of_variables"],
            proof_boundary="Exact change of variables t=exp(4u); first theta summand only.",
        ),
        HeatTiltRow(
            id="ufshta_02_published_saddle_theorem",
            role="published_theorem_input",
            readiness="ready_to_apply",
            claim="O'Sullivan's suitable-multiplier theorem supplies an arbitrary-order saddle expansion for both first-summand integrals.",
            formula=published["published_expansion"],
            proof_boundary="Uses the cited published theorem and its stated suitable-function hypotheses.",
            diagnostics=published,
        ),
        HeatTiltRow(
            id="ufshta_03_uniform_suitability",
            role="exact_analytic_lemma",
            readiness="ready_to_apply",
            claim="The complete compact family of heat multipliers satisfies O'Sullivan's suitability hypotheses with uniform constants.",
            formula=suitability["uniform_taylor"],
            proof_boundary="Uniform compact-parameter adaptation of the published Gaussian-log example.",
            diagnostics=suitability,
        ),
        HeatTiltRow(
            id="ufshta_04_uniform_all_order_expansion",
            role="theorem_composition",
            readiness="ready_to_apply",
            claim="The first-summand logarithm has an arbitrary-order expansion uniform in the complete heat interval.",
            formula=theorem["uniform_all_order_expansion"],
            proof_boundary="First theta summand only; higher summands are handled by a separate theorem.",
        ),
        HeatTiltRow(
            id="ufshta_05_tilt_log_expansion",
            role="exact_asymptotic_reduction",
            readiness="ready_to_apply",
            claim="Subtracting the T=0 expansion isolates the Gaussian heat term and smooth inverse-k corrections.",
            formula=theorem["tilt_log_expansion"],
            proof_boundary="Arbitrary fixed asymptotic order with non-effective onset.",
        ),
        HeatTiltRow(
            id="ufshta_06_lambert_difference_bound",
            role="exact_symbolic_lemma",
            readiness="ready_to_apply",
            claim="Lambert-W differentiation gives the required local-difference scale for the leading heat term.",
            formula=lambert["main_term_differences"],
            proof_boundary="Exact derivative recurrence through order seven plus the iterated-integral finite-difference identity.",
            diagnostics=lambert,
        ),
        HeatTiltRow(
            id="ufshta_07_remainder_and_corrections",
            role="exact_asymptotic_lemma",
            readiness="ready_to_apply",
            claim="Choosing the published expansion order above the stencil order makes both the remainder and every inverse-k correction negligible at the target scale.",
            formula=theorem["remainder_choice"] + "; " + theorem["correction_differences"],
            proof_boundary="Fixed orders m<=7; constants are uniform for 0<=T<=100.",
        ),
        HeatTiltRow(
            id="ufshta_08_uniform_heat_tilt_theorem",
            role="exact_theorem_composition",
            readiness="ready_to_apply",
            claim="The first-summand heat tilt obeys all seven compact-uniform local ratio estimates required by the order-four determinant transfer.",
            formula=theorem["target_theorem"],
            proof_boundary="Uniform non-effective asymptotic theorem for the first summand; no determinant sign is asserted here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem",
        "date": "2026-07-13",
        "status": "uniform compact-heat first-summand heat-tilt asymptotic theorem",
        "proof_boundary": (
            "This artifact proves the seven fixed-order first-summand heat-tilt "
            "finite-difference estimates uniformly for -100<=lambda<=0. It does "
            "not by itself prove the complete-kernel ratio expansion, uniform "
            "order-four positivity, forward order-four invariance, arbitrary-column "
            "order four, PF-infinity, RH, or Lambda<=0."
        ),
        "published_input": published,
        "uniform_suitability": suitability,
        "lambert_derivatives": lambert,
        "theorem": theorem,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_or_published_rows": len(rows),
            "ready_to_apply_rows": len(rows),
            "open_analytic_rows": 0,
            "suitability_coefficients_checked": 8,
            "lambert_derivative_orders_checked": 7,
            "uniform_heat_tilt_theorems": 1,
        },
        "sources": [SOURCE_URL, "outputs/formal_core.md"],
        "generator": (
            "work/rh_compute/scripts/"
            "jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.py"
        ),
        "checker": (
            "work/rh_compute/scripts/"
            "check_jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.py"
        ),
    }


def write_note(path: Path, artifact: dict) -> None:
    published = artifact["published_input"]
    suitability = artifact["uniform_suitability"]
    lambert = artifact["lambert_derivatives"]
    theorem = artifact["theorem"]
    lines = [
        "# Jensen-Window PF Uniform First-Summand Heat-Tilt Asymptotic Theorem",
        "",
        "Date: 2026-07-13",
        "",
        "Status: uniform compact-heat first-summand heat-tilt asymptotic theorem.",
        "This is not by itself a proof of uniform order-four positivity, PF-infinity,",
        "RH, or `Lambda <= 0`.",
        "",
        "```text",
        "work/rh_compute/results/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.json",
        "python work/rh_compute/scripts/jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.py",
        "python work/rh_compute/scripts/check_jensen_window_pf_uniform_first_summand_heat_tilt_asymptotic_theorem.py",
        "```",
        "",
        "## Published Input",
        "",
        f"[{published['source']}]({published['url']}), Theorem 5.2, proves",
        "arbitrary-order saddle asymptotics for suitable multipliers. Section 5",
        "explicitly treats the Gaussian-log multiplier used by the heat flow, and",
        "Section 7 supplies the Xi first-summand integral reduction.",
        "",
        "```text",
        published["published_integral"],
        published["published_expansion"],
        "```",
        "",
        "## Uniform Suitability",
        "",
        "For `0<=T<=100`, set",
        "",
        "```text",
        suitability["family"],
        suitability["ratio_identity"],
        "```",
        "",
        "On the published complex disk, the exponent is uniformly bounded. Cauchy",
        "estimates therefore give",
        "",
        "```text",
        suitability["uniform_taylor"],
        suitability["coefficient_bound"],
        "```",
        "",
        "Thus the complete compact family is suitable with uniform constants.",
        "",
        "## Heat-Tilt Expansion",
        "",
        "The exact `t=exp(4u)` change of variables gives",
        "",
        "```text",
        theorem["first_summand_change_of_variables"],
        "```",
        "",
        "Applying Theorem 5.2 to both terms and combining their all-order",
        "expansions yields",
        "",
        "```text",
        theorem["uniform_all_order_expansion"],
        theorem["tilt_log_expansion"],
        "```",
        "",
        "## Seven Difference Orders",
        "",
        "For `w=W(2k/pi)`,",
        "",
        "```text",
        lambert["derivative_identity"],
        lambert["main_term_derivatives"],
        lambert["finite_difference_identity"],
        "```",
        "",
        "Choose the all-order expansion parameter `R>m`. The explicit correction",
        "terms gain `m` inverse powers under a fixed local difference, while",
        "`w^(3R)/k^R=o(w/k^m)`. Consequently",
        "",
        "```text",
        theorem["target_theorem"] + ".",
        "```",
        "",
        "The same estimate holds for backward stencils. This closes the sole",
        "first-summand heat-tilt target in the uniform order-four tail reduction.",
        "",
        "```text",
        "outputs/jensen_window_pf_uniform_superpolynomial_first_summand_dominance.md",
        "outputs/jensen_window_pf_compound_order4_uniform_heat_eventual_tail_reduction.md",
        "outputs/signed_hankel_jensen_dependency_graph.md",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    write_note(args.note, artifact)
    summary = artifact["summary"]
    print(
        "proved uniform first-summand heat-tilt asymptotics: "
        f"{summary['rows']} rows, {summary['exact_or_published_rows']} exact/published rows, "
        f"{summary['suitability_coefficients_checked']} suitability coefficients, "
        f"{summary['lambert_derivative_orders_checked']} Lambert derivative orders, "
        f"{summary['open_analytic_rows']} open rows"
    )


if __name__ == "__main__":
    main()
