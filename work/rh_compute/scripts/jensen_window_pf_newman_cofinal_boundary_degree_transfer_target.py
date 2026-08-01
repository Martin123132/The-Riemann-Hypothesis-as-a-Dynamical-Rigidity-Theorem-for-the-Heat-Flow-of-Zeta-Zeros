#!/usr/bin/env python3
"""Build the cofinal Newman boundary-degree transfer target."""

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
    "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json"
)
DEFAULT_NOTE = (
    REPO_ROOT
    / "outputs/jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.md"
)


@dataclass(frozen=True)
class TargetRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | list | None = None


def build_exact() -> dict:
    amplitude, log_slope, scale = sp.symbols(
        "amplitude log_slope scale", positive=True
    )
    normalizer_matrix = amplitude * sp.Matrix(
        [[1, 0], [log_slope / scale, 1]]
    )
    normalizer_determinant = sp.factor(normalizer_matrix.det())
    if normalizer_determinant != amplitude**2:
        raise RuntimeError("normalizer orientation determinant failed")

    f_x, f_xx, f_xt = sp.symbols("f_x f_xx f_xt", real=True)
    scale_t, scale_x = sp.symbols("scale_t scale_x", real=True)
    f_t = -f_xx
    scaled_contact_jacobian_x_t = sp.Matrix(
        [
            [f_x, f_t],
            [
                f_xx / scale - f_x * scale_x / scale**2,
                f_xt / scale - f_x * scale_t / scale**2,
            ],
        ]
    )
    contact_determinant_x_t = sp.factor(
        scaled_contact_jacobian_x_t.det().subs(f_x, 0)
    )
    expected_contact_determinant_x_t = f_xx**2 / scale
    if sp.simplify(
        contact_determinant_x_t - expected_contact_determinant_x_t
    ) != 0:
        raise RuntimeError("standard-orientation scaled contact determinant failed")
    contact_determinant_t_x = -contact_determinant_x_t

    t, x = sp.symbols("t x", real=True)
    epsilon, eta = sp.symbols("epsilon eta", positive=True)
    endpoint_model = x**2 + epsilon - 2 * t
    if sp.simplify(
        sp.diff(endpoint_model, t) + sp.diff(endpoint_model, x, 2)
    ) != 0:
        raise RuntimeError("endpoint-safe heat model failed")
    lower_time = epsilon / 2 + eta
    endpoint_jet_norm = sp.factor(
        endpoint_model.subs({t: lower_time, x: 0}) ** 2
        + sp.diff(endpoint_model, x).subs({t: lower_time, x: 0}) ** 2
    )
    if endpoint_jet_norm != 4 * eta**2:
        raise RuntimeError("endpoint-safe first-jet margin failed")
    endpoint_discriminant = sp.factor(
        sp.discriminant(endpoint_model, x)
    )
    if sp.simplify(endpoint_discriminant - (8 * t - 4 * epsilon)) != 0:
        raise RuntimeError("endpoint-safe root discriminant failed")

    return {
        "normalizer_orientation": {
            "normalization": "H=A*Z with A>0 and a=partial_x log(A)",
            "scaled_jets": (
                "V_H=(H,H_x/ell)=A*[[1,0],[a/ell,1]]"
                "*(Z,Z_x/ell)=A*S_a,ell*V_Z"
            ),
            "matrix": str(normalizer_matrix),
            "determinant": str(normalizer_determinant),
            "conclusion": (
                "The positive normalizer, derivative shear, and every positive "
                "jet scale preserve boundary nonvanishing, local degree, and winding."
            ),
        },
        "scaled_contact_index": {
            "heat_equation": "partial_t F=-partial_x^2 F",
            "scaled_map": "V_F=(F,F_x/ell), where ell(t,x)>0",
            "orientation_convention": (
                "Use the standard (x,t) orientation, matching the "
                "counterclockwise bottom-edge traversal x=0 to x=R."
            ),
            "double_contact_determinant_x_t": str(contact_determinant_x_t),
            "double_contact_determinant_t_x": str(contact_determinant_t_x),
            "all_multiplicity_index": (
                "At a spatial multiplicity-m contact, "
                "ind_(x,t)(V_F)=+floor(m/2)>0; reversing to (t,x) "
                "reverses both the index and boundary orientation."
            ),
            "global_degree": (
                "deg(V_F,Omega,0)=wind((F+i*F_x/ell)(partial Omega),0)"
                "=+sum_p floor(m_p/2) in the standard (x,t) orientation."
            ),
        },
        "boundary_rouche": {
            "split": "V_Z=V_P+V_R on partial Omega",
            "proxy": "V_P=(P,P_x/ell)",
            "error": "V_R=(Z-P,(Z_x-P_x)/ell)",
            "strict_domination": (
                "||V_R||_2<||V_P||_2 pointwise on partial Omega"
            ),
            "homotopy": "V_s=V_P+s*V_R for 0<=s<=1",
            "reverse_triangle": (
                "||V_s||_2>=||V_P||_2-s*||V_R||_2>0"
            ),
            "conclusion": (
                "V_Z and V_P are boundary-nonzero and have the same winding."
            ),
        },
        "zero_degree_exclusion": {
            "criterion": (
                "If the boundary domination holds and wind(V_P)=0, "
                "then Omega contains no real multiple-zero contact of Z or H."
            ),
            "reason": (
                "The transferred full winding is zero, while every possible "
                "interior contact has strictly positive local index in the "
                "standard (x,t) orientation."
            ),
            "scope": (
                "The proxy may have interior zeros or mixed local indices; only "
                "its boundary winding enters the argument."
            ),
        },
        "cofinal_contract": {
            "rectangles": "Q_j=[1/(5*j),1/4]x[0,38+j]",
            "known_equivalence": (
                "Lambda<=0 iff the full Xi first jet is boundary-nonzero "
                "with zero winding on every Q_j."
            ),
            "current_base": (
                "Q_207=[1/1035,1/4]x[0,245] is already contact-free "
                "with zero first-jet winding."
            ),
            "open_boundary_theorem": (
                "For every j>=208, construct a continuous boundary proxy P_j, "
                "a positive scale ell_j, and certified errors satisfying strict "
                "boundary domination and wind(V_(P_j))=0."
            ),
            "conditional_endgame": (
                "The open boundary theorem for every j>=208, together with "
                "Q_207 and positive-boundary attainment, implies Lambda<=0 and RH."
            ),
        },
        "ordinary_proxy": {
            "domain": "x>=245, 0<t<=1/5, L=log(x/(4*pi))",
            "main": "O_(h,t)=J_(N_(F,h),t)^F/(16*x^4*A_t)",
            "value_error": "E_F0=exp(-3h)/(25*x^(23/4))",
            "derivative_error": (
                "E_F1=E_F0*(1/2+(1+4/x)/L)"
            ),
            "boundary_target": (
                "O_h^2+(O_h'/L)^2>E_F0^2+E_F1^2 "
                "only on the selected high-frequency boundary arcs."
            ),
            "loose_sufficient_target": (
                "O_h^2+(O_h'/L)^2>"
                "2*exp(-6h)/(625*x^(23/2))"
            ),
        },
        "corrected_proxy": {
            "domain": "L>=50, 0<t<=1/5, t*L<=25",
            "main": "J_hat_(N,t)",
            "error_budget": (
                "r_RS^2+(r_RS'/L)^2<32000000*exp(-3L/2)"
            ),
            "boundary_target": (
                "J_hat^2+(J_hat'/L)^2>"
                "32000000*exp(-3L/2) only on the selected boundary arcs."
            ),
            "guard": (
                "The ordinary and corrected certified envelopes differ by five "
                "exponents; their quantitative targets must not be interchanged."
            ),
        },
        "continuity_contract": {
            "requirement": (
                "P_j must be continuous around the complete closed boundary."
            ),
            "fixed_cutoff_option": (
                "Use a fixed retained cutoff on each finite boundary arc whenever "
                "the tail theorem remains valid after retaining extra terms."
            ),
            "chart_option": (
                "Alternatively join adjacent cutoff charts through overlap "
                "homotopies that retain strict first-jet domination."
            ),
            "warning": (
                "A pointwise family with unbridged cutoff jumps has no defined "
                "boundary winding and cannot discharge the theorem."
            ),
        },
        "endpoint_safe_model": {
            "flow": "G_(epsilon,t)(x)=x^2+epsilon-2t",
            "heat_equation": "partial_t G=-partial_x^2 G",
            "collision_time": "t=epsilon/2",
            "simple_region": (
                "For t>epsilon/2 the zeros are "
                "+/-sqrt(2t-epsilon), real and simple."
            ),
            "lower_boundary": "t=epsilon/2+eta",
            "jet_norm_at_center": str(endpoint_jet_norm),
            "limit": "4*eta^2 tends to zero as eta tends to zero",
            "conclusion": (
                "Contact-free rectangles need not have a time-independent "
                "interior first-jet floor; boundary margins may deteriorate "
                "with the rectangle."
            ),
        },
        "route_decision": {
            "primary": (
                "Use the cofinal boundary-degree contract as the minimal direct "
                "collision route. It keeps the sign-definite heat-contact charge "
                "and asks arithmetic separation only on one-dimensional edges."
            ),
            "secondary": (
                "Retain the global corrected phase-critical small-ball theorem "
                "as a stronger alternative and the all-degree Jensen/PF route "
                "as an independent structural alternative."
            ),
            "next_bounded_input": (
                "A permanent low-time compact base through x=245 and a continuous "
                "proxy/winding pilot for the first boundary beyond Q_207 would "
                "test the contract without pretending that one finite stage is global."
            ),
        },
    }


def build_payload() -> dict:
    exact = build_exact()
    rows = [
        TargetRow(
            id="ncbdt_01_normalizer_orientation",
            role="exact_lemma",
            readiness="available_exact",
            claim="Positive Xi normalization and positive derivative scaling preserve first-jet degree and winding.",
            formula=exact["normalizer_orientation"]["scaled_jets"],
            proof_boundary="Linear first-jet identity only; it proves no boundary margin.",
            diagnostics=exact["normalizer_orientation"],
        ),
        TargetRow(
            id="ncbdt_02_scaled_contact_charge",
            role="exact_composition",
            readiness="available_exact",
            claim="Every multiple heat contact retains a strictly positive index in the standard (x,t) orientation under a positive scaled first jet.",
            formula=exact["scaled_contact_index"]["global_degree"],
            proof_boundary="Composes the independently checked all-multiplicity heat-contact theorem.",
            diagnostics=exact["scaled_contact_index"],
        ),
        TargetRow(
            id="ncbdt_03_boundary_rouche",
            role="exact_lemma",
            readiness="available_exact",
            claim="Strict first-jet error domination on the boundary transfers winding from a continuous proxy to the full function.",
            formula=exact["boundary_rouche"]["reverse_triangle"],
            proof_boundary="Finite-dimensional boundary homotopy; strict domination must still be proved for Xi.",
            diagnostics=exact["boundary_rouche"],
        ),
        TargetRow(
            id="ncbdt_04_zero_degree_exclusion",
            role="exact_composition",
            readiness="available_exact",
            claim="Zero proxy winding plus boundary domination excludes every hidden full heat contact.",
            formula=exact["zero_degree_exclusion"]["criterion"],
            proof_boundary="The conclusion is conditional on both proxy obligations.",
            diagnostics=exact["zero_degree_exclusion"],
        ),
        TargetRow(
            id="ncbdt_05_cofinal_contract",
            role="exact_reduction",
            readiness="available_exact",
            claim="A cofinal family of boundary proxy certificates would prove Lambda<=0 without an interior small-ball theorem.",
            formula=exact["cofinal_contract"]["open_boundary_theorem"],
            proof_boundary="The cofinal proxy family is open and is not constructed by this artifact.",
            diagnostics=exact["cofinal_contract"],
        ),
        TargetRow(
            id="ncbdt_06_q207_base",
            role="exact_composition",
            readiness="available_exact",
            claim="The existing Q207 theorem supplies the current finite base of the cofinal contract.",
            formula=exact["cofinal_contract"]["current_base"],
            proof_boundary="Finite Q207 only; no stage j>=208 follows.",
        ),
        TargetRow(
            id="ncbdt_07_boundary_only_improvement",
            role="route_decision",
            readiness="available_exact",
            claim="The proxy lower bound is required only on one-dimensional boundary arcs, not throughout each two-dimensional rectangle.",
            formula=exact["zero_degree_exclusion"]["criterion"],
            proof_boundary="A weaker sufficient architecture, not a completed Xi estimate.",
        ),
        TargetRow(
            id="ncbdt_08_ordinary_proxy_contract",
            role="exact_composition",
            readiness="available_exact",
            claim="The adaptive ordinary-theta tail theorem supplies an explicit normalized boundary error budget.",
            formula=exact["ordinary_proxy"]["boundary_target"],
            proof_boundary="Retained ordinary first-jet separation and winding remain open.",
            diagnostics=exact["ordinary_proxy"],
        ),
        TargetRow(
            id="ncbdt_09_corrected_proxy_contract",
            role="exact_composition",
            readiness="available_exact",
            claim="The corrected Riemann-Siegel C1 certificate supplies a second explicit boundary error budget.",
            formula=exact["corrected_proxy"]["boundary_target"],
            proof_boundary="Corrected phase-critical boundary separation and winding remain open.",
            diagnostics=exact["corrected_proxy"],
        ),
        TargetRow(
            id="ncbdt_10_proxy_continuity_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="A cutoff-dependent arithmetic proxy must be joined continuously before its winding is meaningful.",
            formula=exact["continuity_contract"]["requirement"],
            proof_boundary="Architecture guard; it does not reject continuous fixed-cutoff or overlap-chart constructions.",
            diagnostics=exact["continuity_contract"],
        ),
        TargetRow(
            id="ncbdt_11_endpoint_safe_model",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="Positive-time simplicity does not require a time-independent first-jet floor down to the zero-time endpoint.",
            formula=(
                f"{exact['endpoint_safe_model']['flow']}; "
                f"jet norm={exact['endpoint_safe_model']['jet_norm_at_center']}"
            ),
            proof_boundary="Universal quadratic heat model, not the Xi heat flow.",
            diagnostics=exact["endpoint_safe_model"],
        ),
        TargetRow(
            id="ncbdt_12_open_xi_boundary_theorem",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="Construct the cofinal continuous Xi boundary proxies and prove their strict domination and zero winding.",
            formula=exact["cofinal_contract"]["open_boundary_theorem"],
            proof_boundary="Open RH-level arithmetic theorem; neither Lambda<=0 nor RH is proved here.",
            diagnostics=exact["route_decision"],
        ),
    ]
    return {
        "kind": "jensen_window_pf_newman_cofinal_boundary_degree_transfer_target",
        "date": "2026-07-26",
        "status": (
            "exact cofinal boundary-degree transfer target and "
            "endpoint-uniformity guard"
        ),
        "proof_boundary": (
            "This artifact proves an exact vector first-jet boundary homotopy "
            "and composes it with the sign-definite heat-contact index and "
            "diagonal exhaustion. It reduces the missing global input to "
            "continuous arithmetic boundary proxies with strict C1 domination "
            "and zero winding. It does not construct that cofinal proxy family, "
            "prove retained first-jet separation, prove Lambda<=0, or prove RH."
        ),
        "sources": [
            "outputs/jensen_window_pf_newman_first_jet_winding_gate.md",
            "outputs/jensen_window_pf_newman_positive_boundary_diagonal_exhaustion_gate.md",
            "outputs/jensen_window_pf_newman_theta_forward_adaptive_sqrt_tail_gate.md",
            "outputs/jensen_window_pf_newman_theta_forward_sqrt_to_corrected_rs_C1_transfer_gate.md",
            "outputs/jensen_window_pf_newman_polymath15_critical_C1_global_remainder_certificate.md",
            "outputs/jensen_window_pf_newman_theta_forward_six_term_finite_bridge_interval_certificate.md",
            "outputs/jensen_window_pf_newman_counterfactual_birth_signature_atlas.md",
        ],
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    normalizer = exact["normalizer_orientation"]
    contact = exact["scaled_contact_index"]
    rouche = exact["boundary_rouche"]
    exclusion = exact["zero_degree_exclusion"]
    cofinal = exact["cofinal_contract"]
    ordinary = exact["ordinary_proxy"]
    corrected = exact["corrected_proxy"]
    continuity = exact["continuity_contract"]
    endpoint = exact["endpoint_safe_model"]
    route = exact["route_decision"]
    return "\n".join(
        [
            "# Newman Cofinal Boundary-Degree Transfer Target",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact cofinal boundary-degree transfer target and",
            "endpoint-uniformity guard. This is not a proof of `Lambda<=0` or RH.",
            "",
            "```text",
            "work/rh_compute/results/jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.json",
            "python work/rh_compute/scripts/jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.py",
            "python work/rh_compute/scripts/check_jensen_window_pf_newman_cofinal_boundary_degree_transfer_target.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Newman cofinal boundary-degree transfer target: 12 rows, 0 issues, 2 orientation identities, 1 vector boundary-Rouche theorem, 1 sign-definite degree composition, 1 Q207 base, 2 explicit C1 proxy budgets, 1 endpoint-uniformity countermodel, 1 open cofinal Xi boundary theorem",
            "```",
            "",
            "## Orientation Invariance",
            "",
            "Write `H=A*Z`, where `A>0`, and let `ell(t,x)>0`. Then",
            "",
            "```text",
            normalizer["scaled_jets"],
            f"det(A*S_a,ell)={normalizer['determinant']}>0",
            "```",
            "",
            normalizer["conclusion"],
            "",
            "At a resolved double contact of a backward-heat solution,",
            "",
            "```text",
            contact["scaled_map"],
            contact["orientation_convention"],
            f"det D_(x,t)V_F={contact['double_contact_determinant_x_t']}>0",
            f"det D_(t,x)V_F={contact['double_contact_determinant_t_x']}<0",
            contact["all_multiplicity_index"],
            contact["global_degree"],
            "```",
            "",
            "Positive normalization and positive scaling therefore do not alter",
            "the sign-definite contact charge.",
            "",
            "## Boundary Rouché Transfer",
            "",
            "On a compact boundary split the normalized scaled jet as",
            "",
            "```text",
            rouche["split"],
            rouche["proxy"],
            rouche["error"],
            rouche["strict_domination"],
            rouche["homotopy"],
            rouche["reverse_triangle"],
            "```",
            "",
            rouche["conclusion"],
            "",
            "Combining that homotopy with the same-sign contact charge gives",
            "",
            "```text",
            exclusion["criterion"],
            exclusion["reason"],
            "```",
            "",
            "The proxy is allowed to have interior zeros of its own. Only its",
            "boundary winding is used.",
            "",
            "## Cofinal Contract",
            "",
            "```text",
            cofinal["rectangles"],
            cofinal["known_equivalence"],
            cofinal["current_base"],
            cofinal["open_boundary_theorem"],
            "```",
            "",
            cofinal["conditional_endgame"],
            "",
            "This is strictly less demanding than a positive first-jet margin",
            "throughout every two-dimensional rectangle. The arithmetic work is",
            "moved to the bottom and right edges plus one integer winding.",
            "",
            "## Ordinary Proxy Budget",
            "",
            "On the existing adaptive ordinary-theta domain,",
            "",
            "```text",
            ordinary["main"],
            ordinary["value_error"],
            ordinary["derivative_error"],
            ordinary["boundary_target"],
            ordinary["loose_sufficient_target"],
            "```",
            "",
            "Only selected high-frequency boundary arcs need this inequality.",
            "",
            "## Corrected Proxy Budget",
            "",
            "On the corrected Riemann--Siegel overlap,",
            "",
            "```text",
            corrected["main"],
            corrected["error_budget"],
            corrected["boundary_target"],
            "```",
            "",
            corrected["guard"],
            "",
            "## Continuity Gate",
            "",
            "```text",
            continuity["requirement"],
            continuity["fixed_cutoff_option"],
            continuity["chart_option"],
            continuity["warning"],
            "```",
            "",
            "A collection of good pointwise estimates is not yet a winding",
            "argument if the retained cutoff jumps without a certified join.",
            "",
            "## Endpoint Guard",
            "",
            "The exact model",
            "",
            "```text",
            endpoint["flow"],
            endpoint["heat_equation"],
            endpoint["simple_region"],
            endpoint["lower_boundary"],
            f"first-jet norm at x=0: {endpoint['jet_norm_at_center']}",
            endpoint["limit"],
            "```",
            "",
            endpoint["conclusion"],
            "",
            "Thus boundary constants may deteriorate with the cofinal rectangle.",
            "No endpoint-simplicity burden is inserted.",
            "",
            "## Route Decision",
            "",
            route["primary"],
            "",
            route["secondary"],
            "",
            "The next bounded input is:",
            "",
            "```text",
            route["next_bounded_input"],
            "```",
            "",
            "The open cofinal proxy construction is the proof-level obligation.",
            "No finite pilot is to be promoted into that theorem.",
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
        "wrote Newman cofinal boundary-degree transfer target: "
        f"{args.out.relative_to(REPO_ROOT).as_posix()}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
