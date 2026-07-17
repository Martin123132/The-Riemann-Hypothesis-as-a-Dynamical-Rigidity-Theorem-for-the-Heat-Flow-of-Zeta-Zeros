#!/usr/bin/env python3
"""Build a C1 corrected-main target excluding positive Newman collisions."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_transversality_target.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_transversality_target.md"
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


def build_exact() -> dict:
    amplitude, ell, z0, z1 = sp.symbols(
        "amplitude ell z0 z1", real=True, nonzero=True
    )
    h0 = amplitude * z0
    h1 = amplitude * (ell * z0 + z1)
    if sp.simplify(h0.subs(z0, 0)) != 0:
        raise RuntimeError("normalized value transfer failed")
    if sp.simplify(h1.subs({z0: 0, z1: 0})) != 0:
        raise RuntimeError("normalized first-jet transfer failed")

    x, y, u, phase_first, q0, q1, scale = sp.symbols(
        "x y u phase_first q0 q1 scale", real=True, nonzero=True
    )
    p0 = 2 * x
    p1 = 2 * (u - phase_first * y)
    j0 = p0 - q0
    j1 = p1 - q1
    transversality = sp.expand(j0**2 + (j1 / scale) ** 2)
    expected = (2 * x - q0) ** 2 + (
        (2 * (u - phase_first * y) - q1) / scale
    ) ** 2
    if sp.simplify(transversality - expected) != 0:
        raise RuntimeError("corrected transversality identity failed")

    r0, r1, eps0, eps1 = sp.symbols(
        "r0 r1 eps0 eps1", real=True, nonnegative=True
    )
    error_norm = r0**2 + (r1 / scale) ** 2
    error_cap = eps0**2 + eps1**2
    if sp.simplify(
        error_norm.subs({r0: eps0, r1: scale * eps1}) - error_cap
    ) != 0:
        raise RuntimeError("C1 transversality error budget failed")

    return {
        "normalization": (
            "H_t=A_t*Z_t with A_t>0; H_t'=A_t*((log A_t)'*Z_t+Z_t')"
        ),
        "double_zero": "H_t=H_t'=0 if and only if Z_t=Z_t'=0",
        "refined_split": (
            "Z_t=J_(N,t)+r_ref, J_(N,t)=P_(N,t)-Q_(N,t), "
            "P_(N,t)=2Re(exp(i*beta_t)D_(N,t))"
        ),
        "components": (
            "exp(i*beta)D=X+iY, exp(i*beta)D'=U+iV1, B=beta'"
        ),
        "main_first_jet": "J=2X-Q, J'=2*(U-BY)-Q'",
        "transversality": (
            "T_L[J]=(2X-Q)^2+((2*(U-BY)-Q')/L)^2"
        ),
        "c1_caps": "|r_ref|<=epsilon_0, |r_ref'|<=L*epsilon_1",
        "collision_implication": (
            "At a double zero, T_L[J]=r_ref^2+(r_ref'/L)^2 "
            "<=epsilon_0^2+epsilon_1^2"
        ),
        "strict_target": (
            "T_L[J]>epsilon_0^2+epsilon_1^2 uniformly on the critical region"
        ),
        "cauchy_gain": (
            "A radius-1/L collar transfers the refined scalar remainder through "
            "one derivative with one factor L; no second-derivative remainder is needed"
        ),
        "global_composition": (
            "Combine the oscillatory-zeta theorem for every fixed "
            "tL>=c_*+epsilon at sufficiently large L, the critical "
            "transversality target on the residual high-frequency layer "
            "0<tL<=c_*+o(1), and compact no-double-zero certificates for every "
            "bounded-L remainder, where c_*=4911678521/1933561194"
        ),
        "conditional_endgame": (
            "If all three regions are closed, positive-boundary attainment gives "
            "Lambda<=0, hence RH"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15ctt_01_normalized_double_zero",
            role="exact_equivalence",
            readiness="ready_to_apply",
            claim="The positive normalizer preserves the value-and-first-derivative double-zero condition.",
            formula=f"{exact['normalization']}; {exact['double_zero']}",
            proof_boundary="Exact on the real axis because A_t is strictly positive.",
        ),
        GateRow(
            id="np15ctt_02_corrected_split",
            role="published_theorem_transfer",
            readiness="ready_to_apply",
            claim="The refined Polymath-15 approximation splits the normalized heat flow into a corrected main and a smaller remainder.",
            formula=exact["refined_split"],
            proof_boundary="Uses the published A+B-C approximation with fixed cutoff on a collar.",
        ),
        GateRow(
            id="np15ctt_03_first_jet_identity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The corrected finite value and first derivative have an exact two-component formula.",
            formula=f"{exact['components']}; {exact['main_first_jet']}",
            proof_boundary="No division by D and valid when the complex finite sum vanishes.",
        ),
        GateRow(
            id="np15ctt_04_sum_of_squares",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="Distance of the corrected main from a double zero is a positive sum of two squares.",
            formula=exact["transversality"],
            proof_boundary="L=log(x/(4*pi))>0 in the Polymath high-frequency region.",
        ),
        GateRow(
            id="np15ctt_05_c1_collision_budget",
            role="exact_exclusion_lemma",
            readiness="ready_to_apply",
            claim="A strict corrected-main first-jet lower bound excludes an exact double zero.",
            formula=f"{exact['c1_caps']}; {exact['collision_implication']}",
            proof_boundary="Requires rigorous normalized C1 remainder caps, including cutoff collars.",
        ),
        GateRow(
            id="np15ctt_06_derivative_saving",
            role="strategy_improvement",
            readiness="conditional_ready",
            claim="The collision target saves one Cauchy derivative compared with global Laguerre-curvature transfer.",
            formula=exact["cauchy_gain"],
            proof_boundary="The scalar refined remainder still needs explicit uniform constants.",
        ),
        GateRow(
            id="np15ctt_07_live_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim="The critical asymptotic problem is reduced to corrected finite first-jet transversality.",
            formula=exact["strict_target"],
            proof_boundary="Open and RH-level; no uniform lower bound is asserted.",
        ),
        GateRow(
            id="np15ctt_08_global_partition",
            role="proof_architecture",
            readiness="conditional_ready",
            claim="The no-collision endgame separates into the proved outer asymptotic region, residual critical high frequency, and bounded-L certificate regions.",
            formula=exact["global_composition"],
            proof_boundary=(
                "The residual critical and bounded-L regions remain open; the "
                "outer theorem has an existential, not numerical, L_epsilon."
            ),
        ),
        GateRow(
            id="np15ctt_09_conditional_endgame",
            role="conditional_bridge",
            readiness="not_ready_to_apply",
            claim="Closing the three-region no-double-zero partition would prove the missing Newman direction.",
            formula=exact["conditional_endgame"],
            proof_boundary="Conditional on the open critical and compact certificates.",
        ),
        GateRow(
            id="np15ctt_10_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Positive point scouts and the exact sum-of-squares identity do not prove the strict lower bound.",
            formula="finite diagnostics != uniform transversality",
            proof_boundary="Not Lambda<=0 or RH.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_polymath15_critical_transversality_target",
        "date": "2026-07-17",
        "status": (
            "exact corrected C1 transversality reduction for positive Newman "
            "collisions; not a proof of the lower bound, Lambda<=0, or RH"
        ),
        "proof_boundary": (
            "This artifact proves the normalized double-zero equivalence, the "
            "corrected finite first-jet identity, and the exact C1 exclusion lemma. "
            "It does not prove uniform first-jet transversality, supply complete "
            "refined remainder constants, close the compact region, exclude a "
            "positive Newman boundary, prove Lambda<=0, or prove RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "sources": [
            "https://arxiv.org/abs/1904.12438",
            "outputs/jensen_window_pf_newman_positive_boundary_attainment_lemma.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_scaled_coercivity_target.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_lehmer_margin_gate.md",
            "outputs/jensen_window_pf_newman_polymath15_oscillatory_zeta_handoff_theorem.md",
            "outputs/jensen_window_pf_newman_polymath15_cancellation_zero_free_wall_gate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical Transversality Target",
            "",
            "Date: 2026-07-17",
            "",
            "Status: exact corrected `C1` reduction of the positive-boundary",
            "collision problem. This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_transversality_target.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_transversality_target.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_transversality_target.py",
            "```",
            "",
            "## Double-Zero Reduction",
            "",
            "On the real axis write",
            "",
            "```text",
            exact["normalization"],
            exact["double_zero"],
            "```",
            "",
            "Thus the Newman endgame does not require a global lower bound for the",
            "full Laguerre curvature. It is enough to keep the normalized value and",
            "first derivative from vanishing simultaneously.",
            "",
            "## Corrected First Jet",
            "",
            "Use the refined Polymath-15 split",
            "",
            "```text",
            exact["refined_split"],
            exact["components"],
            exact["main_first_jet"],
            "```",
            "",
            "The scale-adapted distance from a double zero is exactly",
            "",
            "```text",
            exact["transversality"],
            "```",
            "",
            "This is positive by construction and avoids the indefinite second-jet",
            "terms in the corrected curvature formula.",
            "",
            "## C1 Exclusion Lemma",
            "",
            "Suppose the refined remainder satisfies",
            "",
            "```text",
            exact["c1_caps"],
            "```",
            "",
            "At an exact double zero `Z=Z'=0`, so `J=-r_ref` and `J'=-r_ref'`.",
            "Consequently",
            "",
            "```text",
            exact["collision_implication"],
            "```",
            "",
            "Therefore the concrete live target is",
            "",
            "```text",
            exact["strict_target"],
            "```",
            "",
            "A radius-`1/L` collar now needs only one Cauchy derivative:",
            "",
            "```text",
            exact["cauchy_gain"],
            "```",
            "",
            "## Global Composition",
            "",
            "```text",
            exact["global_composition"],
            exact["conditional_endgame"],
            "```",
            "",
            "The Lehmer stress gate explains why this first-jet target is preferable",
            "to a fixed `kappa*L^2` curvature floor: close pairs make the exact",
            "curvature margin genuinely small. The lower bound above remains open",
            "and is the new proof-facing obligation, not an established sign theorem.",
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
        "built Newman Polymath-15 critical transversality target: "
        f"{len(artifact['rows'])} rows, 3 exact identities, 1 C1 exclusion lemma, "
        "1 open first-jet target"
    )


if __name__ == "__main__":
    main()
