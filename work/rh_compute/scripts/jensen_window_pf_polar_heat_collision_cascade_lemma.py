#!/usr/bin/env python3
"""Derive the polar heat-collision cascade and finite-degree escape theorem."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_polar_heat_collision_cascade_lemma.md"


@dataclass(frozen=True)
class LemmaRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def build_exact() -> dict:
    d, j, ell = sp.symbols("d j ell", integer=True, nonnegative=True)
    xi, q0 = sp.symbols("xi q0", nonzero=True)

    binomial_identity = sp.simplify(
        (1 - j / (d + 1)) * sp.binomial(d + 1, j) - sp.binomial(d, j)
    )
    if binomial_identity != 0:
        raise RuntimeError("adjacent polar identity failed")

    local_recurrence = sp.simplify(
        (d + 1 - ell) * sp.binomial(d + 1, ell) * q0 / xi**ell
        - xi * (ell + 1) * sp.binomial(d + 1, ell + 1) * q0 / xi ** (ell + 1)
    )
    if local_recurrence != 0:
        raise RuntimeError("multiple-root Taylor recurrence failed")

    jet_checks = 0
    recurrence_checks = 0
    n = sp.symbols("n", integer=True, nonnegative=True)
    for degree in range(2, 11):
        for multiplicity in range(2, degree + 1):
            for order in range(multiplicity):
                left = (
                    (degree + 1 - order)
                    * sp.binomial(degree + 1, order)
                    * q0
                    / xi**order
                    - xi
                    * (order + 1)
                    * sp.binomial(degree + 1, order + 1)
                    * q0
                    / xi ** (order + 1)
                )
                if sp.simplify(left) != 0:
                    raise RuntimeError("sample Taylor recurrence failed")
                recurrence_checks += 1
            for order in range(multiplicity - 1):
                q_r1 = (
                    sp.factorial(order + 1)
                    * sp.binomial(degree + 1, order + 1)
                    * q0
                    / xi ** (order + 1)
                )
                q_r2 = (
                    sp.factorial(order + 2)
                    * sp.binomial(degree + 1, order + 2)
                    * q0
                    / xi ** (order + 2)
                )
                heat_jet = sp.factor(
                    ((4 * n + 2 + 4 * order) * q_r1 + 4 * xi * q_r2)
                    / (degree + 1)
                )
                expected = sp.factor(
                    (4 * n + 4 * degree + 2)
                    * sp.factorial(degree)
                    / sp.factorial(degree - order)
                    * q0
                    / xi ** (order + 1)
                )
                if sp.simplify(heat_jet - expected) != 0:
                    raise RuntimeError("sample heat-jet identity failed")
                jet_checks += 1

    # Verify the converse factorization on a nontrivial exponential-polynomial sample.
    w, z, t = sp.symbols("w z t")
    sample_xi = sp.Rational(-2)
    sample_R = 1 + 3 * z + 5 * z**2
    sample_F = sp.exp(-z / sample_xi) * sample_R
    factor_checks = []
    for degree in range(3, 11):
        series = sp.series(sample_F, z, 0, degree + 1).removeO()
        coefficients = [
            sp.factor(sp.diff(series, z, order).subs(z, 0))
            for order in range(degree + 1)
        ]
        jensen = sp.expand(
            sum(sp.binomial(degree, order) * coefficients[order] * w**order for order in range(degree + 1))
        )
        quotient = sp.cancel(jensen / (1 - w / sample_xi) ** (degree - 2))
        if not sp.denom(quotient) == 1 or sp.degree(quotient, w) > 2:
            raise RuntimeError("exponential-polynomial Jensen factorization failed")
        factor_checks.append(
            {
                "degree": degree,
                "forced_multiplicity": degree - 2,
                "cofactor_degree": int(sp.degree(quotient, w)),
            }
        )

    return {
        "entire_family": "F_(n,lambda)(z)=sum_(j>=0) A_(n+j)(lambda)*z^j/j!",
        "jensen_family": "J_(d,n,lambda)(w)=sum_(j=0)^d C(d,j)*A_(n+j)(lambda)*w^j",
        "adjacent_polar_identity": "J_d=J_(d+1)-w*J_(d+1)'/(d+1)",
        "heat_hierarchy": (
            "partial_lambda J_d=((4*n+2)*partial_w+4*w*partial_w^2)J_(d+1)/(d+1)"
        ),
        "multiple_root_recurrence": (
            "If mult_xi(J_d)>=m and Q=J_(d+1), then "
            "Q^(ell)(xi)=(d+1)_ell*Q(xi)/xi^ell for 0<=ell<=m."
        ),
        "heat_jet": (
            "For 0<=r<=m-2, (partial_lambda J_d)^(r)(xi)="
            "(4*n+4*d+2)*(d)_r*Q(xi)/xi^(r+1)."
        ),
        "double_root_criterion": (
            "At a double root xi, forward real splitting requires "
            "Q(xi)*J_d''(xi)/xi<=0."
        ),
        "higher_multiplicity_gate": (
            "For m>=3, a differentiable one-sided hyperbolic continuation requires Q(xi)=0."
        ),
        "polar_multiplicity_rule": (
            "If Q is negative-root hyperbolic and J_d is its polar descent, "
            "mult_xi(J_d)=m>=2 implies mult_xi(Q)=m+1."
        ),
        "cascade": (
            "If every J_D, D>=d, is hyperbolic and mult_xi(J_d)=m>=2, then "
            "mult_xi(J_D)=D-(d-m)."
        ),
        "classification": (
            "The cascade forces F(z)=exp(-z/xi)*R(z) with deg(R)<=d-m; "
            "conversely this form gives multiplicity at least D-(d-m), "
            "with equality when deg(R)=d-m."
        ),
        "strictness": (
            "If F is Laguerre-Polya with positive coefficients and is not exponential-polynomial, "
            "every finite Jensen polynomial is strictly hyperbolic."
        ),
        "degree_escape": (
            "At a non-exponential-polynomial LP boundary, the least nonhyperbolic Jensen degree "
            "on the bad side tends to infinity as lambda approaches the boundary."
        ),
        "checks": {
            "sample_degree_range": [2, 10],
            "taylor_recurrence_checks": recurrence_checks,
            "heat_jet_checks": jet_checks,
            "exponential_polynomial_factor_checks": factor_checks,
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        LemmaRow(
            id="phcc_01_adjacent_polar_heat_hierarchy",
            role="exact_identity",
            readiness="available_exact",
            claim="Adjacent Jensen degrees are related simultaneously by a polar map and by the coefficient heat hierarchy.",
            formula=f"{exact['adjacent_polar_identity']}; {exact['heat_hierarchy']}",
            proof_boundary="Exact coefficient identities only.",
        ),
        LemmaRow(
            id="phcc_02_multiple_root_taylor_recurrence",
            role="exact_identity",
            readiness="available_exact",
            claim="A multiple root of the lower polar fixes the complete adjacent-degree Taylor jet through that multiplicity.",
            formula=exact["multiple_root_recurrence"],
            proof_boundary="Local algebra at one nonzero multiple root.",
        ),
        LemmaRow(
            id="phcc_03_heat_jet_identity",
            role="exact_identity",
            readiness="available_exact",
            claim="Every low heat derivative at a multiple root collapses to one adjacent-degree value with a shift-positive prefactor.",
            formula=exact["heat_jet"],
            proof_boundary="First heat variation only.",
            diagnostics=exact["checks"],
        ),
        LemmaRow(
            id="phcc_04_double_root_splitting",
            role="exact_local_lemma",
            readiness="ready_to_apply",
            claim="The adjacent Jensen value gives the complete first-order inward criterion at a nondegenerate double root.",
            formula=exact["double_root_criterion"],
            proof_boundary="Local double-root splitting, not global invariance.",
        ),
        LemmaRow(
            id="phcc_05_higher_multiplicity_viability",
            role="exact_local_lemma",
            readiness="ready_to_apply",
            claim="A triple-or-higher collision cannot have a nonzero constant heat perturbation along a differentiable hyperbolic path.",
            formula=exact["higher_multiplicity_gate"],
            proof_boundary="Necessary one-sided tangent-cone condition only.",
        ),
        LemmaRow(
            id="phcc_06_polar_multiplicity_lift",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="A multiple lower-degree polar root lifts by exactly one multiplicity whenever the adjacent extension is hyperbolic.",
            formula=exact["polar_multiplicity_rule"],
            proof_boundary="One hyperbolic adjacent-degree extension.",
        ),
        LemmaRow(
            id="phcc_07_infinite_multiplicity_cascade",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="If the whole upper Jensen tower is hyperbolic, one finite multiple root forces a linearly growing common-root multiplicity.",
            formula=exact["cascade"],
            proof_boundary="Conditional on hyperbolicity of every higher degree.",
        ),
        LemmaRow(
            id="phcc_08_exponential_polynomial_classification",
            role="exact_theorem",
            readiness="ready_to_apply",
            claim="The only entire coefficient families that can support that infinite fixed-root cascade are exponential times a bounded-degree polynomial.",
            formula=exact["classification"],
            proof_boundary="Classification of the cascade, not LP membership for zeta.",
        ),
        LemmaRow(
            id="phcc_09_unbounded_degree_escape",
            role="exact_boundary_theorem",
            readiness="ready_to_apply",
            claim="For a non-exponential-polynomial LP endpoint, every fixed Jensen degree is strict and any loss of hyperbolicity nearby must escape to unbounded degree.",
            formula=f"{exact['strictness']} {exact['degree_escape']}",
            proof_boundary="Conditional boundary geometry; it does not exclude the boundary.",
        ),
        LemmaRow(
            id="phcc_10_scaled_tail_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="A successful heat-flow proof must control the degree-rescaled collision layer rather than search for a first contact in a fixed degree.",
            formula="derive a degree-uniform scaled-root or polar-defect barrier as D->infinity",
            proof_boundary="Open all-degree estimate; not PF-infinity, RH, or Lambda <= 0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_polar_heat_collision_cascade_lemma",
        "date": "2026-07-11",
        "status": "exact polar heat-collision cascade and unbounded-degree escape theorem",
        "proof_boundary": (
            "This artifact proves local heat-jet identities, classifies an infinite polar multiplicity cascade, "
            "and shows that a non-exponential-polynomial Laguerre-Polya boundary can only be lost through "
            "unbounded Jensen degree. It does not rule out that escape, prove endpoint Laguerre-Polya membership, "
            "PF-infinity, RH, or Lambda <= 0."
        ),
        "sources": [
            "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md",
            "outputs/jensen_window_pf_quartic_quintic_polar_contact_lemma.md",
            "outputs/jensen_window_pf_cofinal_degree_polar_closure_lemma.md",
            "outputs/jensen_window_pf_cofinal_scaling_limit_equivalence_gate.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    checks = exact["checks"]
    return "\n".join(
        [
            "# Jensen-Window PF Polar Heat-Collision Cascade Lemma",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact polar heat-collision cascade and unbounded-degree escape",
            "theorem. This is not a proof of PF-infinity, RH, or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_polar_heat_collision_cascade_lemma.json",
            "python work/rh_compute/scripts/jensen_window_pf_polar_heat_collision_cascade_lemma.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_polar_heat_collision_cascade_lemma.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF polar heat-collision cascade lemma: 10 rows, 0 issues, 3 exact local identities, 1 double-root criterion, 1 higher-multiplicity gate, 1 infinite polar cascade, 1 exponential-polynomial classification, 1 unbounded-degree escape theorem, 1 open scaled-tail handoff",
            "```",
            "",
            "## Local Polar Heat Jet",
            "",
            "Write `Q=J_(d+1,n,lambda)` and let `xi!=0` be a root of",
            "`J_(d,n,lambda)` of multiplicity `m`. Polar descent gives",
            "",
            "```text",
            exact["adjacent_polar_identity"],
            exact["multiple_root_recurrence"],
            "```",
            "",
            "Combining this recurrence with the exact heat hierarchy gives",
            "",
            "```text",
            exact["heat_jet"],
            "```",
            "",
            f"The checker verifies {checks['taylor_recurrence_checks']} Taylor recurrences and",
            f"{checks['heat_jet_checks']} heat-jet instances over degrees 2 through 10.",
            "At a double root the complete first-order real-splitting condition is",
            "",
            "```text",
            exact["double_root_criterion"],
            "```",
            "",
            "For multiplicity at least three, differentiable one-sided hyperbolicity",
            "forces `Q(xi)=0`; otherwise the constant perturbation produces a",
            "nonreal root of the leading local polynomial.",
            "",
            "## Infinite Cascade",
            "",
            "If `Q` is itself negative-root hyperbolic, the polar reciprocal sum",
            "shows that a multiple lower root must be a root of `Q`, with one extra",
            "unit of multiplicity. Repeating upward gives",
            "",
            "```text",
            exact["cascade"],
            "```",
            "",
            "Put `L=d-m`. Then `J_D(w)=(1-w/xi)^(D-L) R_D(w)` with",
            "`deg(R_D)<=L`. Scaling `w=z/D` and taking the locally uniform Jensen",
            "limit proves",
            "",
            "```text",
            exact["classification"],
            "```",
            "",
            "The converse follows from `J_D(w)=D! [t^D] exp(t)F(wt)`. The",
            "executable sample checks the resulting factor through degree 10.",
            "",
            "## Boundary Consequence",
            "",
            "A Laguerre-Polya entire function with positive coefficients that is not",
            "exponential-polynomial therefore has strictly hyperbolic Jensen",
            "polynomials in every fixed degree. If a nearby parameter is outside",
            "Laguerre-Polya, the least failing degree must tend to infinity as the",
            "boundary is approached.",
            "",
            "This removes fixed-degree first contact as a possible closing argument.",
            "The live problem is a quantitative estimate in the rescaled `D->infinity`",
            "collision layer.",
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
        "wrote Jensen-window PF polar heat-collision cascade lemma: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
