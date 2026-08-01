#!/usr/bin/env python3
"""Build the exact Edrei-moment heat-flow boundary and countermodel gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json"
)
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_edrei_heat_flow_boundary_gate.md"


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def symbolic_audit() -> dict:
    z, lam, beta = sp.symbols("z lambda beta", positive=True)
    m = sp.symbols("m", integer=True, positive=True)
    n, s = sp.symbols("n s", integer=True, nonnegative=True)

    def radial_heat(poly: sp.Expr) -> sp.Expr:
        return sp.expand(4 * z * sp.diff(poly, z, 2) + 2 * sp.diff(poly, z))

    base = (1 + beta * z) ** 2
    heated = sp.expand(
        base + lam * radial_heat(base) + lam**2 * radial_heat(radial_heat(base)) / 2
    )
    if radial_heat(radial_heat(radial_heat(base))) != 0:
        raise RuntimeError("quadratic heat series did not terminate")
    expected_heated = (
        beta**2 * z**2
        + (2 * beta + 12 * lam * beta**2) * z
        + 1
        + 4 * lam * beta
        + 12 * lam**2 * beta**2
    )
    if sp.expand(heated - expected_heated) != 0:
        raise RuntimeError("quadratic heat polynomial identity failed")

    discriminant = sp.factor(sp.discriminant(heated, z))
    expected_discriminant = 32 * lam * beta**3 * (1 + 3 * lam * beta)
    if sp.expand(discriminant - expected_discriminant) != 0:
        raise RuntimeError("quadratic discriminant identity failed")

    def rank_one_flow(index: sp.Expr) -> sp.Expr:
        return (
            2
            * (index + 1)
            * m
            * beta ** (index + 2)
            * (2 * m * (index + 1) - (2 * index + 3))
        )

    determinant_derivative = sp.factor(
        rank_one_flow(s) * (m * beta ** (s + 3))
        + (m * beta ** (s + 1)) * rank_one_flow(s + 2)
        - 2 * (m * beta ** (s + 2)) * rank_one_flow(s + 1)
    )
    expected_derivative = 8 * m**2 * (m - 1) * beta ** (2 * s + 5)
    if sp.expand(determinant_derivative - expected_derivative) != 0:
        raise RuntimeError("rank-one determinant orientation identity failed")

    constant = 1 + 4 * lam * beta + 12 * lam**2 * beta**2
    exact_shifted_minor = sp.factor(
        beta ** (2 * s + 2) * discriminant / constant ** (s + 3)
    )

    sample = {}
    for name, value in (
        ("forward_lambda_1_over_10", sp.Rational(1, 10)),
        ("backward_lambda_minus_1_over_10", sp.Rational(-1, 10)),
    ):
        sample[name] = {
            "lambda": str(value),
            "constant": str(sp.factor(constant.subs({lam: value, beta: 1}))),
            "discriminant": str(
                sp.factor(discriminant.subs({lam: value, beta: 1}))
            ),
            "shift0_minor": str(
                sp.factor(exact_shifted_minor.subs({lam: value, beta: 1, s: 0}))
            ),
        }

    return {
        "radial_heat_operator": "L=4*z*partial_z^2+2*partial_z",
        "quadratic_heat_polynomial": str(heated),
        "quadratic_discriminant": str(discriminant),
        "rank_one_minor_derivative": str(determinant_derivative),
        "exact_shifted_minor": str(exact_shifted_minor),
        "sample_witnesses": sample,
    }


def exact_statements() -> dict:
    return {
        "radial_heat": "partial_lambda F_lambda=(4*z*partial_z^2+2*partial_z)F_lambda",
        "log_derivative_pde": "For R=partial_z log(F): partial_lambda R=4*z*R_zz+8*z*R*R_z+6*R_z+4*R^2",
        "moment_ode": "For R(z)=sum_(n>=0)(-1)^n*a_n*z^n: a_n'= -2*(n+1)*(2*n+3)*a_(n+1)+4*(n+1)*sum_(k=0)^n a_k*a_(n-k)",
        "generating_pde": "For A(t)=sum_(n>=0)a_n*t^n=R(-t): partial_lambda A=-4*t*A_tt+(8*t*A-6)*A_t+4*A^2",
        "rank_one_boundary": "At H(z)=(1+beta*z)^m, a_n=m*beta^(n+1) and every shifted 2x2 Stieltjes minor Delta_s=a_s*a_(s+2)-a_(s+1)^2 vanishes",
        "rank_one_orientation": "At that boundary, partial_lambda Delta_s=8*m^2*(m-1)*beta^(2*s+5)",
        "quadratic_heat": "exp(lambda*L)(1+beta*z)^2=beta^2*z^2+(2*beta+12*lambda*beta^2)*z+1+4*lambda*beta+12*lambda^2*beta^2",
        "quadratic_discriminant": "Disc_z=32*lambda*beta^3*(1+3*lambda*beta)",
        "quadratic_shifted_minor": "Delta_s(lambda)=32*lambda*beta^(2*s+5)*(1+3*lambda*beta)/(1+4*lambda*beta+12*lambda^2*beta^2)^(s+3)",
        "backward_failure": "For -1/(3*beta)<lambda<0 every Delta_s(lambda)<0 and the quadratic zeros are nonreal; generic backward Stieltjes-cone invariance is false",
        "open_handoff": "A backward argument from a de Bruijn real-zero time must use Xi/Phi-specific all-order rigidity or a uniform no-escape theorem, not generic invariance of the Stieltjes moment cone",
    }


def build_payload() -> dict:
    exact = exact_statements()
    rows = [
        GateRow(
            "ehfb_01_radial_heat_operator",
            "exact_identity",
            "available_exact",
            "The squared-variable Phi transform obeys the radial heat equation.",
            exact["radial_heat"],
            "This is an evolution identity, not a zero-location statement.",
        ),
        GateRow(
            "ehfb_02_log_derivative_pde",
            "exact_identity",
            "available_exact",
            "Its logarithmic derivative obeys a closed radial Burgers equation.",
            exact["log_derivative_pde"],
            "The nonlinear PDE does not preserve the desired cone backward by itself.",
        ),
        GateRow(
            "ehfb_03_edrei_moment_ode",
            "exact_identity",
            "available_exact",
            "Coefficient extraction gives a closed hierarchy for the exact Edrei moments a_n=p_(n+1).",
            exact["moment_ode"],
            "All shifts are coupled through a_(n+1) and the convolution term.",
        ),
        GateRow(
            "ehfb_04_stieltjes_generating_pde",
            "exact_identity",
            "available_exact",
            "The ordinary moment generating series obeys an equivalent nonlinear PDE.",
            exact["generating_pde"],
            "Formal equivalence only; it does not provide a positive representing measure.",
        ),
        GateRow(
            "ehfb_05_rank_one_boundary",
            "exact_boundary_model",
            "available_exact",
            "A repeated type-I factor is a rank-one Stieltjes boundary point at every shift.",
            exact["rank_one_boundary"],
            "The multiplicity m is an integer for an entire LP+ factor.",
        ),
        GateRow(
            "ehfb_06_rank_one_orientation",
            "exact_boundary_orientation",
            "available_exact",
            "Every shifted 2x2 Stieltjes wall points strictly into the cone under forward heat when m>=2.",
            exact["rank_one_orientation"],
            "The same sign means immediate exit under backward heat.",
        ),
        GateRow(
            "ehfb_07_double_zero_heat_polynomial",
            "exact_countermodel",
            "guard_validated",
            "The radial heat orbit of a double negative zero is an explicit quadratic.",
            exact["quadratic_heat"],
            "A finite LP+ model, not the Xi function.",
        ),
        GateRow(
            "ehfb_08_double_zero_discriminant",
            "exact_countermodel",
            "guard_validated",
            "The quadratic discriminant changes sign exactly at the repeated-zero boundary.",
            exact["quadratic_discriminant"],
            "Shows the local Newman-type collision mechanism exactly.",
        ),
        GateRow(
            "ehfb_09_all_shift_minor_crossing",
            "exact_countermodel",
            "guard_validated",
            "Every shifted 2x2 Edrei/Stieltjes minor has the same crossing sign.",
            exact["quadratic_shifted_minor"],
            "Rejects rescue by moving to a different finite shift.",
        ),
        GateRow(
            "ehfb_10_backward_invariance_guard",
            "rejected_shortcut",
            "rejected_by_countermodel",
            "The Stieltjes moment cone is not generically invariant under backward radial heat.",
            exact["backward_failure"],
            "Does not exclude a stronger Xi/Phi-specific invariant.",
        ),
        GateRow(
            "ehfb_11_xi_specific_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The surviving dynamical route needs a uniform Xi/Phi rigidity mechanism preventing an all-order boundary escape.",
            exact["open_handoff"],
            "No such all-order Xi/Phi theorem is proved here.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_edrei_heat_flow_boundary_gate",
        "date": "2026-07-23",
        "status": "exact heat-flow lemma, exact backward-invariance countermodel, and open Xi/Phi handoff",
        "exact": exact,
        "symbolic_audit": symbolic_audit(),
        "rows": [asdict(row) for row in rows],
        "sources": [
            "outputs/formal_core.md",
            "outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md",
            "outputs/jensen_window_pf_heat_flow_jensen_hierarchy_lemma.md",
        ],
        "proof_boundary": (
            "This artifact proves the radial log-derivative and Edrei-moment heat "
            "hierarchies and gives an exact LP+ polynomial countermodel to generic "
            "backward Stieltjes-cone invariance. It does not prove an Xi/Phi-specific "
            "invariant, all-order Stieltjes positivity at lambda zero, PF-infinity, "
            "Jensen hyperbolicity for zeta, RH, or Lambda<=0."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    audit = payload["symbolic_audit"]
    return "\n".join(
        [
            "# Jensen-Window PF Edrei Heat-Flow Boundary Gate",
            "",
            "Date: 2026-07-23",
            "",
            "Status: exact heat-flow lemma, exact backward-invariance countermodel,",
            "and open Xi/Phi handoff. This is not a proof of PF-infinity, Jensen",
            "hyperbolicity for zeta, RH, or `Lambda <= 0`.",
            "",
            "Artifact kind: `jensen_window_pf_edrei_heat_flow_boundary_gate`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_edrei_heat_flow_boundary_gate.json",
            "python work/rh_compute/scripts/jensen_window_pf_edrei_heat_flow_boundary_gate.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_edrei_heat_flow_boundary_gate.py",
            "```",
            "",
            "## Exact Heat Hierarchy",
            "",
            "For the unnormalized squared-variable Phi transform,",
            "",
            "```text",
            exact["radial_heat"],
            "```",
            "",
            "and normalization does not change `R=F'/F`. Direct differentiation gives",
            "",
            "```text",
            exact["log_derivative_pde"],
            exact["moment_ode"],
            exact["generating_pde"],
            "```",
            "",
            "This is the heat hierarchy for the exact endpoint moments",
            "`a_n=p_(n+1)`, not for the original linear signed-Hankel sequence.",
            "",
            "## Rank-One Boundary Orientation",
            "",
            "At a repeated type-I factor,",
            "",
            "```text",
            exact["rank_one_boundary"],
            exact["rank_one_orientation"],
            "```",
            "",
            "Thus a multiplicity `m>=2` points strictly into every shifted `2x2`",
            "Stieltjes wall under forward heat and strictly out under backward heat.",
            "",
            "## Exact Double-Zero Countermodel",
            "",
            "For `L=4z partial_z^2+2 partial_z`, the heat series terminates:",
            "",
            "```text",
            exact["quadratic_heat"],
            exact["quadratic_discriminant"],
            exact["quadratic_shifted_minor"],
            "```",
            "",
            "The denominator in the shifted-minor formula is positive for real",
            "`lambda`. Hence, for `-1/(3*beta)<lambda<0`, every shifted minor is",
            "negative and the two zeros are nonreal. At `lambda=0` they collide at",
            "`z=-1/beta`; for `lambda>0` they are distinct and negative.",
            "",
            "Exact beta=1 witnesses:",
            "",
            "```text",
            f"lambda=1/10:  Disc={audit['sample_witnesses']['forward_lambda_1_over_10']['discriminant']}, Delta_0={audit['sample_witnesses']['forward_lambda_1_over_10']['shift0_minor']}",
            f"lambda=-1/10: Disc={audit['sample_witnesses']['backward_lambda_minus_1_over_10']['discriminant']}, Delta_0={audit['sample_witnesses']['backward_lambda_minus_1_over_10']['shift0_minor']}",
            "```",
            "",
            "## Consequence",
            "",
            exact["backward_failure"] + ".",
            "",
            "Therefore the known real-zero regime at a positive de Bruijn time cannot",
            "be propagated to lambda zero by a generic backward-invariance assertion",
            "for the Stieltjes cone. The surviving target is:",
            "",
            "```text",
            exact["open_handoff"],
            "```",
            "",
            "The countermodel is a finite LP+ polynomial and not the Xi function. It",
            "rejects only the generic shortcut; it leaves genuinely Xi/Phi-specific",
            "rigidity open.",
            "",
            "## Proof Boundary",
            "",
            payload["proof_boundary"],
            "",
        ]
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "wrote Jensen-window PF Edrei heat-flow boundary gate: "
        "11 rows, 4 exact flow identities, 1 all-shift boundary formula, "
        "2 exact heat witnesses, 1 rejected shortcut, 1 open Xi/Phi handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
