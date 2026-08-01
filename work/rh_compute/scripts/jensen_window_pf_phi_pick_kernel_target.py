#!/usr/bin/env python3
"""Build the exact Phi Pick-kernel endpoint target and mixture guard."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results/jensen_window_pf_phi_pick_kernel_target.json"
DEFAULT_NOTE = REPO_ROOT / "outputs/jensen_window_pf_phi_pick_kernel_target.md"


@dataclass(frozen=True)
class TargetRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def mixture_audit() -> dict:
    weights = (Fraction(9, 10), Fraction(1, 10))
    scales = (Fraction(1, 4), Fraction(5, 2))
    moment1 = sum(weight * scale for weight, scale in zip(weights, scales))
    moment2 = sum(weight * scale**2 for weight, scale in zip(weights, scales))
    wall = moment2 - 3 * moment1**2
    derivative = wall / 12
    if moment1 != Fraction(19, 40):
        raise RuntimeError("unexpected first scale moment")
    if moment2 != Fraction(109, 160):
        raise RuntimeError("unexpected second scale moment")
    if wall != Fraction(7, 1600) or derivative != Fraction(7, 19200):
        raise RuntimeError("two-scale Pick-wall witness failed")
    return {
        "measure": "(9/10)*delta_(1/4)+(1/10)*delta_(5/2) in y=u^2",
        "E_y": fraction_text(moment1),
        "E_y2": fraction_text(moment2),
        "E_y2_minus_3_E_y_squared": fraction_text(wall),
        "R_prime_0": fraction_text(derivative),
        "conclusion": "R'(0)>0, so R is not Stieltjes and the Pick numerator is negative for z=i*epsilon with sufficiently small epsilon>0",
    }


def exact_statements() -> dict:
    return {
        "phi_transform": "F(z)=integral_R Phi(u)*cosh(u*sqrt(z))*du, F(x)>0 for x>=0",
        "kernel_atoms": "C_u(z)=cosh(u*sqrt(z)), D_u(z)=partial_z C_u(z)=u*sinh(u*sqrt(z))/(2*sqrt(z))",
        "log_derivative": "R(z)=F'(z)/F(z)=integral Phi(u)D_u(z)du / integral Phi(u)C_u(z)du",
        "pick_numerator": "P_Phi(z):=-Im(F'(z)*conj(F(z)))=|F(z)|^2*(-Im R(z))",
        "polarized_kernel": "K_z(u,v):=-(1/2)*Im(D_u(z)*conj(C_v(z))+D_v(z)*conj(C_u(z)))",
        "polarization": "P_Phi(z)=double_integral Phi(u)*Phi(v)*K_z(u,v)du dv",
        "pick_equivalence": "For the nonexponential Phi transform: F in LP+ <=> P_Phi(z)>0 for every Im(z)>0",
        "stieltjes_characterization": "R is Stieltjes <=> R(x)>=0 for x>0, R is holomorphic on C\\(-infinity,0], and Im R(z)<=0 for Im(z)>0",
        "first_quadrant": "For z=w^2 and w=a+i*b with a,b>0: 2*|w|^2*P_Phi(w^2)=-Im(conj(w)*M'(w)*conj(M(w))), M(w)=integral Phi(u)cosh(u*w)du",
        "tilted_identity": "For y=u^2, C_x(y)=cosh(sqrt(x*y)), r_x(y)=partial_x log C_x(y), and dPi_x=C_x*dmu/F(x): R'(x)=E_Pi[r_x']+Var_Pi(r_x)",
        "first_wall": "At x=0: R'(0)=(E[y^2]-3*E[y]^2)/12, so the first real-axis Stieltjes wall is E[y^2]<=3*E[y]^2",
        "operator_route": "A direct identity R(z)=gamma+<v,(I+z*T)^(-1)v> with T positive self-adjoint would give the required Stieltjes representation by the spectral theorem",
        "open_target": "Prove P_Phi(z)>0 on the whole upper half-plane, or construct the positive resolvent representation, using verified Xi/Phi structure rather than zero reality",
    }


def build_payload() -> dict:
    exact = exact_statements()
    rows = [
        TargetRow(
            "ppkt_01_phi_transform",
            "exact_identity",
            "available_exact",
            "The squared-variable endpoint is a positive Phi mixture of fixed-scale cosh kernels.",
            exact["phi_transform"],
            "Positive mixing alone does not preserve LP+.",
        ),
        TargetRow(
            "ppkt_02_kernel_atoms",
            "exact_identity",
            "available_exact",
            "The transform and its derivative have explicit entire fixed-scale atoms.",
            exact["kernel_atoms"],
            "The removable value at z=0 is understood by its power series.",
        ),
        TargetRow(
            "ppkt_03_log_derivative_ratio",
            "exact_identity",
            "available_exact",
            "The Edrei endpoint is the exact Phi logarithmic-derivative ratio.",
            exact["log_derivative"],
            "A ratio of positive integrals on the positive axis is not yet a Stieltjes function.",
        ),
        TargetRow(
            "ppkt_04_pick_numerator",
            "exact_identity",
            "available_exact",
            "The anti-Pick sign is equivalent to positivity of one real numerator.",
            exact["pick_numerator"],
            "Strict positivity also excludes off-axis zeros of F.",
        ),
        TargetRow(
            "ppkt_05_polarized_double_kernel",
            "exact_identity",
            "available_exact",
            "Polarization writes the Pick numerator as a quadratic form in the known Phi density.",
            exact["polarization"],
            "The cross-scale kernel K_z(u,v) is not pointwise nonnegative in general.",
        ),
        TargetRow(
            "ppkt_06_pick_endpoint_equivalence",
            "exact_equivalence",
            "ready_to_apply",
            "Krein's Stieltjes/Pick characterization and Sokal's criterion make the global Phi-kernel sign endpoint-equivalent to LP+.",
            exact["pick_equivalence"],
            "This identifies the exact target but proves none of its upper-half-plane signs.",
        ),
        TargetRow(
            "ppkt_07_first_quadrant_form",
            "exact_identity",
            "available_exact",
            "The square-root map converts the target to a real two-parameter first-quadrant inequality.",
            exact["first_quadrant"],
            "All a,b>0 are required; finite grids are nonpromotable.",
        ),
        TargetRow(
            "ppkt_08_tilted_variance_wall",
            "exact_identity",
            "available_exact",
            "On the positive boundary the first Pick wall is a tilted variance-versus-curvature inequality.",
            exact["tilted_identity"],
            "This is necessary but not sufficient for the global Pick condition.",
        ),
        TargetRow(
            "ppkt_09_origin_kurtosis_wall",
            "exact_necessary_condition",
            "available_exact",
            "At the origin the first wall reduces to a sharp scale-kurtosis inequality.",
            exact["first_wall"],
            "Only the first real-axis wall, not an all-order theorem.",
        ),
        TargetRow(
            "ppkt_10_two_scale_pick_failure",
            "exact_countermodel",
            "guard_validated",
            "A positive two-scale mixture violates the first Pick wall exactly.",
            "For (9/10)delta_(1/4)+(1/10)delta_(5/2): R'(0)=7/19200>0",
            "Abstract positive scale measure, not the Xi/Phi density.",
        ),
        TargetRow(
            "ppkt_11_positive_resolvent_route",
            "sufficient_theorem_target",
            "not_ready_to_apply",
            "A Phi-derived positive self-adjoint resolvent realization would close the Stieltjes endpoint.",
            exact["operator_route"],
            "No noncircular operator T and cyclic vector v are constructed here.",
        ),
        TargetRow(
            "ppkt_12_xi_phi_pick_handoff",
            "open_theorem_target",
            "not_ready_to_apply",
            "The surviving direct target is global positivity of the polarized Phi Pick kernel or an equivalent positive resolvent.",
            exact["open_target"],
            "The target must be proved without assuming real zeros, LP+, RH, or Lambda<=0.",
        ),
    ]
    return {
        "kind": "jensen_window_pf_phi_pick_kernel_target",
        "date": "2026-07-23",
        "status": "exact endpoint-equivalent Pick-kernel target, exact first-wall identity, and positive-mixture guard",
        "exact": exact,
        "mixture_audit": mixture_audit(),
        "rows": [asdict(row) for row in rows],
        "sources": [
            {
                "title": "Christian Berg, Stieltjes-Pick-Bernstein-Schoenberg and their connection to complete monotonicity, Theorem 3.2",
                "url": "https://web.math.ku.dk/~berg/manus/castellon.pdf",
            },
            {
                "title": "Kalmykov and Karp, When does a hypergeometric function belong to LP+, Proposition 6",
                "url": "https://doi.org/10.1016/j.jmaa.2022.126432",
            },
            {
                "title": "Local Edrei-Stieltjes endpoint equivalence gate",
                "path": "outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md",
            },
        ],
        "proof_boundary": (
            "This artifact proves an exact Pick-kernel reformulation of the common "
            "LP+/coefficient-PF/Jensen/Stieltjes endpoint, a tilted first-wall "
            "identity, and an exact positive-mixture countermodel. It does not prove "
            "the global Phi-kernel sign, construct a positive resolvent, prove "
            "PF-infinity, Jensen hyperbolicity for zeta, RH, or Lambda<=0."
        ),
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    audit = payload["mixture_audit"]
    return "\n".join(
        [
            "# Jensen-Window PF Phi Pick-Kernel Target",
            "",
            "Date: 2026-07-23",
            "",
            "Status: exact endpoint-equivalent Pick-kernel target, exact first-wall",
            "identity, and positive-mixture guard. This is not a proof of PF-infinity,",
            "Jensen hyperbolicity for zeta, RH, or `Lambda <= 0`.",
            "",
            "Artifact kind: `jensen_window_pf_phi_pick_kernel_target`.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_phi_pick_kernel_target.json",
            "python work/rh_compute/scripts/jensen_window_pf_phi_pick_kernel_target.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_phi_pick_kernel_target.py",
            "```",
            "",
            "## Exact Pick Reformulation",
            "",
            "Set",
            "",
            "```text",
            exact["phi_transform"],
            exact["kernel_atoms"],
            exact["log_derivative"],
            "```",
            "",
            "For `Im(z)>0`, define",
            "",
            "```text",
            exact["pick_numerator"],
            exact["polarized_kernel"],
            exact["polarization"],
            "```",
            "",
            "Christian Berg's Theorem 3.2 gives the Stieltjes/Pick characterization",
            "",
            "```text",
            exact["stieltjes_characterization"],
            "```",
            "",
            "while Sokal's logarithmic-derivative criterion identifies `R` being",
            "Stieltjes with `F` being LP+. Since the nonexponential Phi transform is",
            "positive on the nonnegative real axis, these facts give the exact target",
            "",
            "```text",
            exact["pick_equivalence"],
            "```",
            "",
            "Strict positivity excludes an upper-half-plane zero because such a zero",
            "would make `P_Phi=0`; reflection excludes lower-half-plane zeros, and",
            "`F(x)>0` excludes nonnegative real zeros.",
            "",
            "Under `z=w^2`, `w=a+i*b`, this becomes",
            "",
            "```text",
            exact["first_quadrant"],
            "```",
            "",
            "This is a real two-parameter inequality in the known Phi kernel. It is",
            "endpoint-equivalent, not a completed estimate.",
            "",
            "## First Real-Axis Wall",
            "",
            "Push the normalized Phi scale measure to `y=u^2` and tilt it by",
            "`C_x(y)=cosh(sqrt(x*y))`. Then",
            "",
            "```text",
            exact["tilted_identity"],
            exact["first_wall"],
            "```",
            "",
            "Thus even the first boundary wall is a concentration inequality: mixing",
            "variance must not exceed the average fixed-scale negative curvature.",
            "Passing this wall for every x is still only necessary for the global",
            "upper-half-plane Pick condition.",
            "",
            "## Exact Positive-Mixture Guard",
            "",
            "For the positive two-scale law",
            "",
            "```text",
            audit["measure"],
            f"E[y]={audit['E_y']}",
            f"E[y^2]={audit['E_y2']}",
            f"E[y^2]-3E[y]^2={audit['E_y2_minus_3_E_y_squared']}",
            f"R'(0)={audit['R_prime_0']}>0",
            "```",
            "",
            "so the Stieltjes derivative sign and the Pick numerator fail near",
            "`z=i*epsilon`. This proves that the polarized cross-scale kernel cannot",
            "be pointwise nonnegative in general. The witness is not the Xi density.",
            "",
            "## Surviving Structural Routes",
            "",
            "A direct global route is",
            "",
            "```text",
            exact["open_target"],
            "```",
            "",
            "An operator route is",
            "",
            "```text",
            exact["operator_route"],
            "```",
            "",
            "The spectral theorem would then supply the positive Stieltjes measure;",
            "entireness and Sokal's criterion would recover the discrete Edrei form.",
            "Constructing such an operator from Phi is open and cannot use the zero",
            "reality it is meant to prove.",
            "",
            "## Sources",
            "",
            "- Christian Berg, `Stieltjes-Pick-Bernstein-Schoenberg and their connection to complete monotonicity`, Theorem 3.2: https://web.math.ku.dk/~berg/manus/castellon.pdf",
            "- Kalmykov and Karp, `When does a hypergeometric function belong to LP+?`, Proposition 6: https://doi.org/10.1016/j.jmaa.2022.126432",
            "- `outputs/jensen_window_pf_edrei_stieltjes_equivalence_gate.md`",
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
        "wrote Jensen-window PF Phi Pick-kernel target: "
        "12 rows, 7 exact identities, 1 endpoint equivalence, "
        "1 exact mixture guard, 2 open structural routes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
