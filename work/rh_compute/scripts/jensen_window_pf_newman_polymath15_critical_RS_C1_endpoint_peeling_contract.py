#!/usr/bin/env python3
"""Build the first omitted Riemann-Siegel endpoint peeling contract."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUT = (
    REPO_ROOT
    / "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.md"
)
POLYMATH_SOURCE = "https://arxiv.org/abs/1904.12438"
FLINT_DOC_SOURCE = "https://flintlib.org/doc/acb_dirichlet.html"
FLINT_D_SOURCE = (
    "https://github.com/flintlib/flint/blob/master/"
    "src/acb_dirichlet/zeta_rs_d_coeffs.c"
)
FLINT_F_SOURCE = (
    "https://github.com/flintlib/flint/blob/master/"
    "src/acb_dirichlet/zeta_rs_f_coeffs.c"
)
FLINT_R_SOURCE = (
    "https://github.com/flintlib/flint/blob/master/"
    "src/acb_dirichlet/zeta_rs_r.c"
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
    sigma, v, heat, pi, p = sp.symbols(
        "sigma v heat pi p", real=True, positive=True
    )
    f1, f3 = sp.symbols("F1 F3", complex=True)

    # The FLINT recurrence starts k=1 from d_0^(0)=1, d_1^(0)=0.
    recurrence_u = 1 - 2 * sigma
    previous_d0 = sp.Integer(1)
    previous_d1 = sp.Integer(0)
    d1 = sp.simplify((previous_d1 / 2 + previous_d0 * recurrence_u) / 2)
    d0 = sp.simplify((previous_d0 / 2) / 6)
    if sp.simplify(d0 - sp.Rational(1, 12)) != 0:
        raise RuntimeError("d_0^(1) recurrence failed")
    if sp.simplify(d1 - (1 - 2 * sigma) / 2) != 0:
        raise RuntimeError("d_1^(1) recurrence failed")
    assembled = sp.expand(
        (d0 * f3 + pi * d1 * f1 / (2 * sp.I)) / pi**2
    )
    expected = f3 / (12 * pi**2) + sp.I * (2 * sigma - 1) * f1 / (4 * pi)
    if sp.simplify(assembled - expected) != 0:
        raise RuntimeError("first Riemann-Siegel coefficient assembly failed")

    critical = sp.simplify(expected.subs(sigma, sp.Rational(1, 2)))
    if sp.simplify(critical - f3 / (12 * pi**2)) != 0:
        raise RuntimeError("critical-line C1 specialization failed")

    fluctuating = sp.expand(
        expected.subs(sigma, sp.Rational(1, 2) + sp.sqrt(heat) * v)
    )
    gaussian_mean = sp.simplify(fluctuating.subs(v, 0))
    if sp.simplify(gaussian_mean - critical) != 0:
        raise RuntimeError("centered Gaussian C1 cancellation failed")

    f = (
        sp.exp(sp.I * pi * (p**2 / 2 + sp.Rational(3, 8)))
        - sp.I * sp.sqrt(2) * sp.cos(pi * p / 2)
    ) / (2 * sp.cos(pi * p))
    if sp.simplify(f.subs(p, -p) - f) != 0:
        raise RuntimeError("C0 parity identity failed")
    critical_c1 = sp.diff(f, p, 3) / (12 * pi**2)
    if sp.simplify(critical_c1.subs(p, -p) + critical_c1) != 0:
        raise RuntimeError("critical C1 oddness failed")

    temperature, saddle = sp.symbols("T_0 a", positive=True)
    saddle_relation = sp.simplify(
        (saddle**-2).subs(saddle, sp.sqrt(temperature / (2 * pi)))
    )
    if sp.simplify(saddle_relation - 2 * pi / temperature) != 0:
        raise RuntimeError("Riemann-Siegel scale conversion failed")

    return {
        "published_expansion": (
            "R_0,N(s)=Dirichlet_N(s)+(-1)^(N-1)U*a^(-sigma)"
            "*(sum_(k=0)^K C_k(p,sigma)/a^k+RS_K(s))"
        ),
        "coefficient_assembly": (
            "C_k(p,sigma)=pi^(-2k)*sum_(j=0)^floor(3k/2)"
            "(pi/(2i))^j*d_j^(k)(sigma)*F^(3k-2j)(p)"
        ),
        "first_recurrence": (
            "d_0^(0)=1, d_1^(0)=0 => "
            "d_0^(1)=1/12, d_1^(1)=(1-2sigma)/2"
        ),
        "first_coefficient": (
            "C_1(p,sigma)=F'''(p)/(12*pi^2)"
            "+i*(2sigma-1)*F'(p)/(4*pi), with F=C_0"
        ),
        "critical_specialization": (
            "C_1(p,1/2)=F'''(p)/(12*pi^2)"
        ),
        "heat_fluctuation": (
            "C_1(p,1/2+sqrt(t)*v)=F'''(p)/(12*pi^2)"
            "+i*sqrt(t)*v*F'(p)/(2*pi)"
        ),
        "centered_gaussian": (
            "For dmu(v)=pi^(-1/2)exp(-v^2)dv, "
            "integral C_1(p,1/2+sqrt(t)*v)dmu(v)="
            "F'''(p)/(12*pi^2)"
        ),
        "scale_conversion": (
            "a=sqrt(T_0/(2*pi)); a^(-1)=sqrt(2*pi/T_0), "
            "a^(-2)=2*pi/T_0"
        ),
        "critical_endpoint_correction": (
            "DeltaC_t(x)=2*(-1)^N*exp(t*pi^2/64)*"
            "Re(M_0(i*T_0)*U*exp(pi*i/8)*"
            "F'''(p)/(12*pi^2*a))"
        ),
        "normalized_update": (
            "DeltaQ=DeltaC_t/A_t; Q_[1]=Q_[0]+DeltaQ; "
            "J_[1]=P-Q_[1]=J_[0]-DeltaQ"
        ),
        "remainder_peel": (
            "Z=J_[0]+r_[0]=J_[1]+r_[1], hence "
            "r_[0]=-DeltaQ+r_[1] and "
            "r_[0],x=-DeltaQ_x+r_[1],x"
        ),
        "critical_cutoff_parity": (
            "F(-p)=F(p) => C_1(-p,1/2)=-C_1(p,1/2)"
        ),
        "critical_cutoff_value_match": (
            "At a=m, p_left=-1, p_right=1, and "
            "(-1)^(m-1)C_1(-1,1/2)=(-1)^m C_1(1,1/2)"
        ),
        "critical_cutoff_derivative_jump": (
            "For E_N(T)=(-1)^N*a^(-1)*C_1(p_N(T),1/2), "
            "p_N'(T)=-a/T and "
            "[E_N']_(right-left)=-(-1)^m*F''''(1)/(6*pi^2*T)"
        ),
        "order_handoff": (
            "After retaining C_1, the fixed-sigma Riemann-Siegel tail begins "
            "at a^(-2)=2*pi/T_0; the heat-integrated proof must also retain "
            "the O(T_0^-1) gamma/Stirling and log-M_0 Taylor terms"
        ),
        "uniform_target": (
            "Prove explicit value and x-derivative bounds for r_[1] on "
            "L>=50 and 0<t*L<=25, uniformly across fixed-N cells and "
            "adjacent-cutoff splices, then insert r_0=-DeltaQ into the "
            "signed contact-normal peeling inequality"
        ),
    }


def build_artifact() -> dict:
    exact = build_exact()
    rows = [
        GateRow(
            id="np15c1ep_01_published_expansion",
            role="published_primary_input",
            readiness="available_published",
            claim=(
                "Polymath 15 records the Arias de Reyna endpoint expansion "
                "with coefficients C_k(p,sigma) and an explicit remainder."
            ),
            formula=exact["published_expansion"],
            proof_boundary=(
                "Proposition 6.2; the paper uses only C_0 in its heat-flow "
                "endpoint approximation."
            ),
        ),
        GateRow(
            id="np15c1ep_02_flint_coefficient_assembly",
            role="primary_source_extraction",
            readiness="available_exact",
            claim=(
                "FLINT's certified Riemann-Siegel implementation exposes the "
                "normalization relating d_j^(k), F derivatives, and C_k."
            ),
            formula=exact["coefficient_assembly"],
            proof_boundary=(
                "Read directly from zeta_rs_r.c together with the documented "
                "Riemann-Siegel normalization."
            ),
        ),
        GateRow(
            id="np15c1ep_03_first_recurrence",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The in-place d-coefficient recurrence gives the complete k=1 row.",
            formula=exact["first_recurrence"],
            proof_boundary="Exact rational recurrence; no numerical inference.",
        ),
        GateRow(
            id="np15c1ep_04_first_coefficient",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="The first omitted endpoint coefficient has an explicit two-jet formula.",
            formula=exact["first_coefficient"],
            proof_boundary=(
                "The removable values of F=C_0 and its derivatives are "
                "understood by analytic continuation at p=+-1/2."
            ),
        ),
        GateRow(
            id="np15c1ep_05_critical_specialization",
            role="exact_identity",
            readiness="ready_to_apply",
            claim="On the critical line the F' part of C_1 vanishes identically.",
            formula=exact["critical_specialization"],
            proof_boundary="Exact substitution sigma=1/2.",
        ),
        GateRow(
            id="np15c1ep_06_centered_heat_cancellation",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "The remaining sigma fluctuation in C_1 is odd and has zero "
                "centered Gaussian mean."
            ),
            formula=(
                exact["heat_fluctuation"] + "; " + exact["centered_gaussian"]
            ),
            proof_boundary=(
                "Exact for the frozen-prefactor coefficient layer. The full "
                "endpoint integrand has additional O(T_0^-1) factors."
            ),
        ),
        GateRow(
            id="np15c1ep_07_corrected_endpoint_definition",
            role="exact_definition",
            readiness="ready_to_apply",
            claim=(
                "The critical corrected main can absorb the complete first "
                "Riemann-Siegel endpoint coefficient explicitly."
            ),
            formula=(
                exact["critical_endpoint_correction"]
                + "; "
                + exact["normalized_update"]
                + "; "
                + exact["remainder_peel"]
            ),
            proof_boundary=(
                "This is an exact repartition of the old remainder once "
                "DeltaC_t is defined; it does not yet bound the new remainder."
            ),
        ),
        GateRow(
            id="np15c1ep_08_second_order_handoff",
            role="published_theorem_transfer",
            readiness="conditional_ready",
            claim=(
                "Extracting C_1 removes the a^-1 endpoint obstruction and "
                "moves every remaining formal endpoint contribution to order "
                "a^-2 or smaller."
            ),
            formula=exact["scale_conversion"] + "; " + exact["order_handoff"],
            proof_boundary=(
                "The fixed-sigma order follows from the published expansion. "
                "Uniform heat integration, constants, and first derivatives "
                "still require a new proof."
            ),
        ),
        GateRow(
            id="np15c1ep_09_signed_contact_handoff",
            role="analytic_handoff",
            readiness="conditional_ready",
            claim=(
                "The extracted endpoint block supplies the explicit signed "
                "remainder model required by the contact-normal equation."
            ),
            formula=exact["remainder_peel"],
            proof_boundary=(
                "Useful only after explicit C1 bounds for r_[1] and r_[1],x "
                "are proved."
            ),
        ),
        GateRow(
            id="np15c1ep_10_critical_cutoff_parity",
            role="exact_identity",
            readiness="ready_to_apply",
            claim=(
                "Critical-line parity makes the added C_1 endpoint value "
                "continuous across every Riemann-Siegel cutoff."
            ),
            formula=(
                exact["critical_cutoff_parity"]
                + "; "
                + exact["critical_cutoff_value_match"]
                + "; "
                + exact["critical_cutoff_derivative_jump"]
            ),
            proof_boundary=(
                "The one-sided derivative mismatch is explicit and of order "
                "T^-1; a holomorphic collar splice must still absorb it."
            ),
        ),
        GateRow(
            id="np15c1ep_11_cutoff_guard",
            role="route_guard",
            readiness="guard_validated",
            claim=(
                "Exact value matching does not by itself provide a global "
                "holomorphic C_1 collar."
            ),
            formula=(
                "critical value parity + fixed-N analyticity + adjacent-N "
                "derivative control are jointly required"
            ),
            proof_boundary=(
                "The remaining mismatch is already at the desired T^-1 "
                "scale, but its explicit collar constant is not certified here."
            ),
        ),
        GateRow(
            id="np15c1ep_12_uniform_quantitative_target",
            role="open_theorem_target",
            readiness="not_ready_to_apply",
            claim=(
                "The next proof obligation is a uniform second-order endpoint "
                "value/first-jet theorem compatible with signed contact peeling."
            ),
            formula=exact["uniform_target"],
            proof_boundary=(
                "Open. No global contact exclusion, Lambda<=0, or RH is "
                "asserted."
            ),
        ),
    ]
    return {
        "kind": (
            "jensen_window_pf_newman_polymath15_critical_RS_C1_"
            "endpoint_peeling_contract"
        ),
        "date": "2026-07-26",
        "status": (
            "exact extraction of the first omitted Riemann-Siegel endpoint "
            "coefficient and its critical Gaussian cancellation; the uniform "
            "second-order heat remainder remains open"
        ),
        "proof_boundary": (
            "This artifact derives C_1 exactly from the published/implemented "
            "Riemann-Siegel coefficient recurrence, specializes it to the "
            "critical line, and identifies the signed remainder peel. It does "
            "not prove the required uniform second-order value/derivative "
            "bounds, the adjacent-cutoff splice, contact exclusion, Lambda<=0, "
            "or RH."
        ),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_coefficient_identities": 5,
            "primary_source_inputs": 2,
            "corrected_main_definitions": 1,
            "conditional_transfers": 2,
            "route_guards": 1,
            "open_quantitative_targets": 1,
        },
        "sources": [
            POLYMATH_SOURCE,
            FLINT_DOC_SOURCE,
            FLINT_D_SOURCE,
            FLINT_F_SOURCE,
            FLINT_R_SOURCE,
            "outputs/jensen_window_pf_newman_polymath15_endpoint_holomorphic_lift.md",
            "outputs/jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate.md",
        ],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 Critical C1 Endpoint Peeling Contract",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact first-coefficient extraction and critical-line",
            "cancellation. The uniform second-order heat remainder is open.",
            "This is not a proof of `Lambda <= 0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_polymath15_critical_RS_C1_endpoint_peeling_contract.py",
            "```",
            "",
            "## Source Normalization",
            "",
            f"[Polymath 15]({POLYMATH_SOURCE}), Proposition 6.2, writes",
            "",
            "```text",
            exact["published_expansion"],
            "```",
            "",
            "The current [FLINT Riemann-Siegel documentation]",
            f"({FLINT_DOC_SOURCE}) and implementation expose the coefficient",
            "assembly suppressed in the Polymath notation:",
            "",
            "```text",
            exact["coefficient_assembly"],
            exact["first_recurrence"],
            "```",
            "",
            "Substitution gives the exact first omitted coefficient",
            "",
            "```text",
            exact["first_coefficient"],
            "```",
            "",
            "The identity is source-derived, not fitted from numerical data.",
            "",
            "## Critical Cancellation",
            "",
            "At the real Newman boundary the zeta argument has `sigma=1/2`, so",
            "",
            "```text",
            exact["critical_specialization"],
            exact["heat_fluctuation"],
            exact["centered_gaussian"],
            "```",
            "",
            "Thus the first coefficient's fluctuating `F'` term is odd under",
            "the centered heat Gaussian. This cancellation applies to the",
            "frozen-prefactor coefficient layer. The exact gamma/Stirling and",
            "`log M_0` defects must be retained at the next order.",
            "",
            "## Explicit Peel",
            "",
            "Define the added critical endpoint block by",
            "",
            "```text",
            exact["critical_endpoint_correction"],
            exact["normalized_update"],
            exact["remainder_peel"],
            "```",
            "",
            "The last line is the key contact-normal handoff: relative to the",
            "old corrected main, the explicit signed remainder model is",
            "`r_0=-DeltaQ`, with derivative `r_0,x=-DeltaQ_x`.",
            "",
            "## Order Gain",
            "",
            "```text",
            exact["scale_conversion"],
            exact["order_handoff"],
            "```",
            "",
            "The old endpoint error starts at `a^-1`. Keeping `C_1` removes",
            "that term. The published fixed-sigma expansion then starts at",
            "`a^-2`, the same order as the first gamma/Stirling and",
            "`log M_0` Taylor defects. This is an asymptotic order audit, not",
            "yet the uniform heat-integrated bound needed by the proof.",
            "",
            "## Cutoff Parity",
            "",
            "The critical endpoint coefficient has the parity needed at a",
            "Riemann-Siegel cutoff:",
            "",
            "```text",
            exact["critical_cutoff_parity"],
            exact["critical_cutoff_value_match"],
            exact["critical_cutoff_derivative_jump"],
            "```",
            "",
            "Thus switching from the left cutoff to the right cutoff creates",
            "no `a^-1` value jump. The first one-sided derivative mismatch is",
            "already proportional to `T^-1`, at the desired second-order",
            "scale. A holomorphic collar still needs an explicit bound for it.",
            "",
            "## Remaining Theorem",
            "",
            "```text",
            exact["uniform_target"],
            "```",
            "",
            "Three pieces must be proved together: explicit Gaussian-integrated",
            "constants, the `x` derivative of the new remainder, and a",
            "second-order adjacent-cutoff derivative splice. Only then can the smaller",
            "remainder be inserted into the signed contact equation.",
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
        "built Newman Polymath-15 critical C1 endpoint peeling contract: "
        f"{len(artifact['rows'])} rows, 5 exact coefficient identities, "
        "2 primary-source inputs, 1 corrected-main definition, "
        "2 conditional transfers, 1 route guard, "
        "1 open quantitative target"
    )


if __name__ == "__main__":
    main()
