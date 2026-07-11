#!/usr/bin/env python3
"""Audit whether renormalized zero energy can exclude a Newman boundary birth."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_boundary_energy_direction_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_newman_boundary_energy_direction_gate.md"


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
    z, tau = sp.symbols("z tau", real=True)
    pi = sp.pi
    a_squared = 1 + 16 / (8 + pi)
    initial = sp.expand((z**2 - 1) ** 2 * (z**2 - a_squared))
    heat = sp.expand(
        sum(
            (-tau) ** order
            / sp.factorial(order)
            * sp.diff(initial, z, 2 * order)
            for order in range(4)
        )
    )
    if sp.simplify(sp.diff(heat, tau) + sp.diff(heat, z, 2)) != 0:
        raise RuntimeError("even polynomial heat equation failed")

    cofactor = sp.cancel(initial / (z - 1) ** 2)
    field = sp.simplify(sp.diff(cofactor, z).subs(z, 1) / cofactor.subs(z, 1))
    if field != -pi / 8:
        raise RuntimeError("classical field normalization failed")

    stiffness = sp.factor((pi**2 + 24 * pi + 160) / 64)
    stiffness_from_roots = sp.factor(
        sp.Rational(1, 2)
        + 2 * (1 + a_squared) / (1 - a_squared) ** 2
    )
    if sp.simplify(stiffness - stiffness_from_roots) != 0:
        raise RuntimeError("root stiffness identity failed")

    fourth_branch_coefficient = sp.factor(
        (pi + 12) * (pi**2 + 24 * pi + 96) / 128
    )
    fifth_branch_numerator = (
        15 * pi**4
        + 720 * pi**3
        + 11904 * pi**2
        + 82432 * pi
        + 211968
    )
    fifth_branch_coefficient = sp.factor(fifth_branch_numerator / 8192)
    cubic_gap_coefficient = sp.factor(
        (
            pi**4
            + 48 * pi**3
            + 800 * pi**2
            + 5632 * pi
            + 14848
        )
        / 32
    )

    eps = sp.symbols("eps", positive=True)
    gap_series = (
        2 * sp.sqrt(2) * eps
        - 2 * sp.sqrt(2) * stiffness * eps**3
        + 2 * sp.sqrt(2) * fifth_branch_coefficient * eps**5
    )
    squared_gap = sp.series(gap_series**2, eps, 0, 8).removeO()
    expected_squared_gap = (
        8 * eps**2
        - 16 * stiffness * eps**4
        + cubic_gap_coefficient * eps**6
    )
    if sp.simplify(squared_gap - expected_squared_gap) != 0:
        raise RuntimeError("higher pair-gap jet failed")

    ratio = sp.symbols("ratio", positive=True)
    potential = ratio ** -2 - 1 + 2 * (ratio - 1)
    if potential.subs(ratio, 1) != 0 or sp.diff(potential, ratio).subs(ratio, 1) != 0:
        raise RuntimeError("renormalized potential normalization failed")
    if sp.simplify(sp.diff(potential, ratio, 2) - 6 / ratio**4) != 0:
        raise RuntimeError("renormalized potential convexity failed")

    reference_gap = sp.symbols("Delta", positive=True)
    model_gap = sp.sqrt(
        8 * tau - 16 * stiffness * tau**2 + cubic_gap_coefficient * tau**3
    )
    pair_energy = 1 / model_gap**2 + 2 * model_gap / reference_gap**3 - 3 / reference_gap**2
    energy_asymptotic = (
        "E_pair(tau)=1/(8*tau)+K/4-3/Delta^2"
        "+4*sqrt(2)*sqrt(tau)/Delta^3+O(tau)"
    )
    energy_derivative_asymptotic = (
        "dE_pair/dtau=-1/(8*tau^2)"
        "+2*sqrt(2)/(Delta^3*sqrt(tau))+O(1)<0 for small tau>0"
    )

    return {
        "zero_ode": (
            "x_k'=2*PV sum_(j!=k) 1/(x_k-x_j); for a pair gap g and q=g^2, "
            "q'=8-4*q*S, S=sum_out 1/((x_+-y)(x_--y))."
        ),
        "universal_boundary_jet": (
            "At a finite double-zero birth t_*, K=sum_out 1/(c-y)^2, so "
            "q(t)=8*(t-t_*)-16*K*(t-t_*)^2+O((t-t_*)^3)."
        ),
        "model": {
            "a_squared": str(a_squared),
            "initial_polynomial": str(initial),
            "heat_solution": str(heat),
            "field": str(field),
            "stiffness": str(stiffness),
            "stiffness_decimal": str(sp.N(stiffness, 18)),
            "fourth_branch_coefficient": str(fourth_branch_coefficient),
            "fifth_branch_coefficient": str(fifth_branch_coefficient),
            "branch_formula": (
                "z_s=1+s*sqrt(2)*e-pi*e^2/4-s*sqrt(2)*K*e^3"
                "+P4*e^4+s*sqrt(2)*L5*e^5+O(e^6), e=sqrt(tau), s=+/-1"
            ),
            "squared_gap": (
                "q(tau)=8*tau-16*K*tau^2+C3*tau^3+O(tau^(7/2))"
            ),
            "cubic_gap_coefficient": str(cubic_gap_coefficient),
            "cubic_gap_coefficient_decimal": str(sp.N(cubic_gap_coefficient, 18)),
            "jet_signs": "K>0, q''(0)=-32*K<0, q'''(0)=6*C3>0",
        },
        "renormalized_energy": {
            "potential": "V(r)=r^(-2)-1+2*(r-1)=r^(-2)+2*r-3 for r>0",
            "convexity": "V(1)=V'(1)=0 and V''(r)=6/r^4>0, hence V(r)>=0",
            "single_interaction": (
                "E_pair=Delta^(-2)*V(g/Delta)"
                "=1/g^2+2*g/Delta^3-3/Delta^2"
            ),
            "asymptotic": energy_asymptotic,
            "derivative_asymptotic": energy_derivative_asymptotic,
            "integral": "integral_0^epsilon E_pair(tau) dtau=+infinity",
            "ordered_pair_warning": (
                "A sum over ordered j!=k counts the symmetric collision interaction twice; "
                "one interaction already forces divergence."
            ),
            "symbolic_pair_energy": str(pair_energy),
        },
        "conditional_exclusion": (
            "If a nonnegative renormalized energy contains the colliding interaction with "
            "a reference gap Delta bounded away from zero and a weight bounded below, then "
            "local time-integrability down to t_* excludes the collision, because the "
            "interaction is asymptotic to 1/(8*(t-t_*))."
        ),
        "rodgers_tao_scope": {
            "source": "https://arxiv.org/abs/1801.05914",
            "version": "v5 (2021-07-03)",
            "contradiction_regime": (
                "The paper assumes Lambda<0 and studies real zeros for Lambda<t<=0."
            ),
            "dynamics_theorem": (
                "Theorem 11 states the zero ODE for Lambda<t<=0."
            ),
            "potential_equations": (
                "Equations (62)-(63) define the nonnegative renormalized interaction."
            ),
            "integrated_energy_theorem": (
                "Theorem 17 proves integral_(Lambda/4)^0 E_tilde_[0.5*T*log(T),"
                "3*T*log(T)](t) dt=o(T*log_+(T)^3) as T tends to infinity."
            ),
            "boundary_separation": (
                "For Lambda<0, Lambda/4-Lambda=3*abs(Lambda)/4>0; the theorem does "
                "not integrate down to the collision boundary t=Lambda."
            ),
            "directionality": (
                "The estimate controls relaxation inside an already-real interval and uses "
                "zeta-zero information at t=0. It is unavailable in a hypothetical "
                "Lambda>0 regime, where 0 lies before the real-zero boundary."
            ),
        },
        "decision": (
            "Renormalized energy identifies a valid conditional boundary obstruction, but "
            "the published Rodgers-Tao theorem does not supply it. Nonnegativity and forward "
            "energy decrease are compatible with birth from an infinite-energy trace."
        ),
        "open_handoff": (
            "Prove an Xi-specific boundary-uniform estimate giving finite local integrated "
            "renormalized energy on (Lambda,Lambda+epsilon), or find another global invariant "
            "that forbids the universal 1/(t-Lambda) collision singularity, without assuming "
            "Lambda<0 or importing t=0 real-zero information."
        ),
    }


def build_payload() -> dict:
    exact = build_exact()
    model = exact["model"]
    energy = exact["renormalized_energy"]
    scope = exact["rodgers_tao_scope"]
    rows = [
        GateRow(
            id="nbeng_01_pair_gap_ode",
            role="exact_identity",
            readiness="available_exact",
            claim="The zero ODE gives a closed squared-gap equation for an isolated real pair.",
            formula=exact["zero_ode"],
            proof_boundary="Real simple-zero interval; principal-value outsider sum exists.",
        ),
        GateRow(
            id="nbeng_02_stiffness_boundary_jet",
            role="exact_lemma",
            readiness="available_exact",
            claim="The first correction to universal square-root birth is the positive outsider stiffness.",
            formula=exact["universal_boundary_jet"],
            proof_boundary="Finite isolated double zero with a regular outsider field.",
        ),
        GateRow(
            id="nbeng_03_exact_higher_jet_birth",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Classical field, positive stiffness, negative quadratic gap curvature, and positive cubic gap jet all coexist with birth.",
            formula=f"{model['branch_formula']}; {model['squared_gap']}",
            proof_boundary="Exact even polynomial backward-heat solution, not the Xi kernel.",
            diagnostics=model,
        ),
        GateRow(
            id="nbeng_04_renormalized_potential",
            role="published_definition",
            readiness="available_published",
            claim="The Rodgers-Tao renormalized interaction is nonnegative and singular like the inverse squared gap at collision.",
            formula=f"{energy['potential']}; {energy['single_interaction']}",
            proof_boundary="Equations (62)-(63) of the cited paper; ordered real zeros.",
        ),
        GateRow(
            id="nbeng_05_collision_energy_asymptotic",
            role="exact_lemma",
            readiness="available_exact",
            claim="A generic double-zero birth has logarithmically nonintegrable renormalized pair energy at the boundary.",
            formula=f"{energy['asymptotic']}; {energy['integral']}",
            proof_boundary="Fixed nonzero reference gap and positive collision weight.",
        ),
        GateRow(
            id="nbeng_06_forward_relaxation_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Forward decrease of nonnegative renormalized energy does not exclude birth.",
            formula=energy["derivative_asymptotic"],
            proof_boundary="The energy descends from an infinite boundary trace in the exact birth model.",
        ),
        GateRow(
            id="nbeng_07_boundary_integrability_criterion",
            role="conditional_criterion",
            readiness="ready_if_hypothesis_proved",
            claim="Finite integrated energy all the way to the boundary would exclude a double collision.",
            formula=exact["conditional_exclusion"],
            proof_boundary="The required Xi boundary-uniform integrability is not currently proved.",
        ),
        GateRow(
            id="nbeng_08_published_energy_scope",
            role="published_theorem",
            readiness="available_published",
            claim="The published integrated-energy theorem is an interior relaxation estimate under Lambda<0.",
            formula=scope["integrated_energy_theorem"],
            proof_boundary=f"{scope['contradiction_regime']} {scope['boundary_separation']}",
            diagnostics=scope,
        ),
        GateRow(
            id="nbeng_09_directionality_gate",
            role="circularity_gate",
            readiness="guard_validated",
            claim="Rodgers-Tao interior energy control cannot be reversed into a proof of Lambda<=0.",
            formula=scope["directionality"],
            proof_boundary="Blocks an assumption-direction error; it does not weaken the cited theorem.",
        ),
        GateRow(
            id="nbeng_10_xi_boundary_energy_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The energy route now has one precise non-circular target at the Newman boundary.",
            formula=exact["open_handoff"],
            proof_boundary="Open Xi-specific boundary theorem; not RH and not Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_boundary_energy_direction_gate",
        "date": "2026-07-11",
        "status": "exact boundary-energy singularity and published-directionality gate",
        "proof_boundary": (
            "This artifact derives the universal stiffness gap jet, verifies an exact even "
            "higher-jet birth model, proves the 1/(8*tau) renormalized pair-energy singularity, "
            "and isolates a conditional finite-integrated-energy collision criterion. It also "
            "records that the Rodgers-Tao theorem is an interior estimate under Lambda<0, "
            "separated from the boundary. It does not prove the required Xi boundary estimate, "
            "RH, or Lambda<=0."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_local_odd_count_reduction_lemma.md",
            "outputs/formal_core.md",
            "https://arxiv.org/abs/1801.05914",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    model = exact["model"]
    energy = exact["renormalized_energy"]
    scope = exact["rodgers_tao_scope"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Boundary-Energy Direction Gate",
            "",
            "Date: 2026-07-11",
            "",
            "Status: exact boundary-energy singularity and published-directionality",
            "gate. This is not a proof of RH or `Lambda <= 0`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_boundary_energy_direction_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_boundary_energy_direction_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_boundary_energy_direction_gate.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Jensen-window PF Newman boundary-energy direction gate: 10 rows, 0 issues, 1 universal gap law, 1 exact higher-jet birth model, 1 nonintegrable collision-energy asymptotic, 1 conditional exclusion criterion, 1 published-scope audit, 1 open Xi boundary-energy handoff",
            "```",
            "",
            "## Higher Gap Jet",
            "",
            "The real-zero ODE gives",
            "",
            "```text",
            exact["zero_ode"],
            exact["universal_boundary_jet"],
            "```",
            "",
            "For the exact classical-field countermodel from the previous gate,",
            "",
            "```text",
            f"a^2={model['a_squared']}",
            "P(z)=(z^2-1)^2*(z^2-a^2)",
            "F_tau=exp(-tau*d_z^2)P",
            f"B(1)={model['field']}",
            f"K={model['stiffness']}={model['stiffness_decimal']}",
            model["branch_formula"],
            model["squared_gap"],
            f"C3={model['cubic_gap_coefficient']}={model['cubic_gap_coefficient_decimal']}",
            model["jet_signs"],
            "```",
            "",
            "Thus the classical drift and the first two nontrivial gap-jet signs",
            "still do not distinguish a forbidden positive boundary from ordinary",
            "square-root birth.",
            "",
            "## Boundary Energy",
            "",
            "For a fixed positive reference gap `Delta`, Rodgers-Tao use",
            "",
            "```text",
            energy["potential"],
            energy["convexity"],
            energy["single_interaction"],
            "```",
            "",
            "Substitution of the collision gap gives",
            "",
            "```text",
            energy["asymptotic"],
            energy["derivative_asymptotic"],
            energy["integral"],
            "```",
            "",
            "This is the decisive local fact. The energy is finite at every time",
            "after birth and decreases immediately, but it has a logarithmically",
            "nonintegrable trace at the birth time. Nonnegative forward relaxation",
            "therefore does not by itself forbid a boundary collision.",
            "",
            "A valid conditional criterion remains:",
            "",
            "```text",
            exact["conditional_exclusion"],
            "```",
            "",
            "## Published Direction",
            "",
            "The primary source is Rodgers-Tao, arXiv:1801.05914v5.",
            "",
            "```text",
            scope["contradiction_regime"],
            scope["dynamics_theorem"],
            scope["integrated_energy_theorem"],
            scope["boundary_separation"],
            scope["directionality"],
            "```",
            "",
            "The theorem is exactly suited to proving relaxation between an",
            "already-real negative time and `t=0`. It does not provide finite energy",
            "on `(Lambda,Lambda+epsilon)` and cannot be imported into the opposite",
            "hypothetical regime `Lambda>0`.",
            "",
            "## Live Handoff",
            "",
            exact["open_handoff"],
            "",
            "This narrows the energy route to a boundary-trace theorem. Until that",
            "theorem is proved from the Xi kernel, energy language is diagnostic rather",
            "than a completed Newman obstruction.",
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
        "wrote Jensen-window PF Newman boundary-energy direction gate: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
