#!/usr/bin/env python3
"""Build the exact parabolic-frequency energy/current audit."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = "jensen_window_pf_newman_parabolic_frequency_energy_current_gate"
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "ray_aligned": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_ray_aligned_parabolic_frequency_reduction.json"
    ),
    "contact_normal": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_parabolic_frequency_contact_normal_hierarchy_gate.json"
    ),
    "q209_shell": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_q208_base_q209_shell_coverage_gate.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }


def audit_sources(payloads: dict[str, dict]) -> dict:
    ray = payloads["ray_aligned"]
    if ray.get("exact", {}).get("scale") != (
        "s_pf(t,x)=(L(x)^2+(2t)^(-1))^(-1/2)"
    ):
        raise RuntimeError("ray-aligned parabolic-frequency scale drifted")
    if ray.get("exact", {}).get("dimensionless_parameter") != "q=2tL^2":
        raise RuntimeError("ray-aligned q coordinate drifted")
    ray_summary = ray.get("summary", {})
    if ray_summary.get("new_strip_open_antecedents") != 0:
        raise RuntimeError("ray-aligned new-strip closure drifted")
    if ray_summary.get("old_collar_open_antecedents") != 1:
        raise RuntimeError("ray-aligned outer-collar handoff drifted")

    contact = payloads["contact_normal"]
    contact_exact = contact.get("exact", {})
    if "H_tau(y)=y^2-2tau" not in contact_exact.get(
        "quadratic_countermodel", ""
    ):
        raise RuntimeError("contact-normal quadratic countermodel drifted")
    if "scale cannot regularize a genuine collision" not in (
        contact_exact.get("scale_derivative_guard", "")
    ):
        raise RuntimeError("contact-normal scale guard drifted")
    if contact.get("summary", {}).get("open_xi_targets") != 1:
        raise RuntimeError("contact-normal open Xi target drifted")

    shell = payloads["q209_shell"]
    shell_summary = shell.get("summary", {})
    if shell_summary.get("q209_certified") is not False:
        raise RuntimeError("Q209 shell source was unexpectedly promoted")
    if shell_summary.get("ray_aligned_old_collar_antecedents") != 1:
        raise RuntimeError("Q209 ray-aligned handoff drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "source_status": {
            key: payload.get("status", "")
            for key, payload in payloads.items()
        },
        "inherited_boundary": (
            "Q208 remains the proved finite base, Q209 remains uncertified, "
            "and the ray-aligned successor route still has exactly one open "
            "Xi outer-descendant theorem."
        ),
    }


def symbolic_audit() -> dict:
    H, H_x, H_xx, H_xxx = sp.symbols(
        "H H_x H_xx H_xxx", real=True
    )
    s, s_t, s_x = sp.symbols("s s_t s_x", real=True)
    energy_t = (
        -2 * H * H_xx
        + 2 * s * s_t * H_x**2
        - 2 * s**2 * H_x * H_xxx
    )
    flux_x = (
        H_x**2
        + H * H_xx
        + 2 * s * s_x * H_x * H_xx
        + s**2 * H_xx**2
        + s**2 * H_x * H_xxx
    )
    bulk = (
        2 * (s * H_xx + s_x * H_x) ** 2
        + 2 * (1 + s * s_t - s_x**2) * H_x**2
    )
    if sp.expand(energy_t + 2 * flux_x - bulk) != 0:
        raise RuntimeError("general energy/current identity failed")

    t, x = sp.symbols("t x", positive=True)
    L = sp.log(x / (4 * sp.pi))
    q = 2 * t * L**2
    s_pf = sp.sqrt(2 * t) / sp.sqrt(1 + q)
    identities = {
        "time_log": sp.simplify(sp.diff(sp.log(s_pf), t)),
        "space_log": sp.simplify(sp.diff(sp.log(s_pf), x)),
        "s_s_t": sp.simplify(s_pf * sp.diff(s_pf, t)),
        "s_x_squared": sp.simplify(sp.diff(s_pf, x) ** 2),
    }
    targets = {
        "time_log": 1 / (2 * t * (1 + q)),
        "space_log": -2 * t * L / (x * (1 + q)),
        "s_s_t": 1 / (1 + q) ** 2,
        "s_x_squared": 4 * q * t**2 / (x**2 * (1 + q) ** 3),
    }
    for key, value in identities.items():
        if sp.simplify(value - targets[key]) != 0:
            raise RuntimeError(f"parabolic-frequency derivative failed: {key}")

    y, tau, tau_0, s_0 = sp.symbols(
        "y tau tau_0 s_0", real=True, positive=True
    )
    model = y**2 - 2 * (tau - tau_0)
    if sp.simplify(sp.diff(model, tau) + sp.diff(model, y, 2)) != 0:
        raise RuntimeError("shifted quadratic backward-heat model failed")
    model_x = sp.diff(model, y)
    model_xx = sp.diff(model, y, 2)
    model_energy = model**2 + s_0**2 * model_x**2
    model_flux = model * model_x + s_0**2 * model_x * model_xx
    model_rhs = 2 * s_0**2 * model_xx**2 + 2 * model_x**2
    if sp.expand(
        sp.diff(model_energy, tau)
        + 2 * sp.diff(model_flux, y)
        - model_rhs
    ) != 0:
        raise RuntimeError("shifted quadratic energy identity failed")
    contact = {y: 0, tau: tau_0}
    contact_values = {
        "E": sp.simplify(model_energy.subs(contact)),
        "E_t": sp.simplify(sp.diff(model_energy, tau).subs(contact)),
        "J": sp.simplify(model_flux.subs(contact)),
        "J_x": sp.simplify(sp.diff(model_flux, y).subs(contact)),
        "rhs": sp.simplify(model_rhs.subs(contact)),
    }
    expected_contact = {
        "E": 0,
        "E_t": 0,
        "J": 0,
        "J_x": 4 * s_0**2,
        "rhs": 8 * s_0**2,
    }
    for key, value in contact_values.items():
        if sp.simplify(value - expected_contact[key]) != 0:
            raise RuntimeError(f"quadratic contact audit failed: {key}")

    return {
        "general_balance_residual": "0",
        "scale_derivative_residuals": {
            key: "0" for key in identities
        },
        "quadratic_heat_residual": "0",
        "quadratic_balance_residual": "0",
        "quadratic_contact_values": {
            "E": "0",
            "E_t": "0",
            "J": "0",
            "J_x": "4*s_0^2",
            "rhs": "8*s_0^2",
        },
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    ratio_cap = Fraction(1, 4 * 38**2)
    if ratio_cap != Fraction(1, 5776):
        raise RuntimeError("outer-domain ratio cap failed")
    parabolic_surplus = (
        Fraction(1, 1) - ratio_cap
    ) * Fraction(1, 4)
    parabolic_coefficient = Fraction(1, 1) + parabolic_surplus
    if parabolic_coefficient != Fraction(28879, 23104):
        raise RuntimeError("q<=1 coefficient bound failed")

    return {
        "assumptions": (
            "H is real and sufficiently smooth, H_t=-H_xx, and "
            "s(t,x)>0 is C1 on a fixed spatial interval."
        ),
        "energy_density": "E_s=H^2+s^2 H_x^2",
        "direct_time_derivative": (
            "partial_t E_s=-2H H_xx+2s s_t H_x^2"
            "-2s^2 H_x H_xxx"
        ),
        "flux": "J_s=H H_x+s^2 H_x H_xx",
        "flux_derivative": (
            "partial_x J_s=H_x^2+H H_xx+2s s_x H_x H_xx"
            "+s^2 H_xx^2+s^2 H_x H_xxx"
        ),
        "local_balance": (
            "partial_t E_s+2 partial_x J_s="
            "2(s H_xx+s_x H_x)^2"
            "+2(1+s s_t-s_x^2)H_x^2"
        ),
        "bulk_density": (
            "B_s=(s H_xx+s_x H_x)^2"
            "+(1+s s_t-s_x^2)H_x^2"
        ),
        "scale": (
            "L=log(x/(4*pi)), q=2tL^2, "
            "s_pf=sqrt(2t)/sqrt(1+q)"
        ),
        "pi_provenance": (
            "pi is the ordinary circle constant already present in the "
            "established Xi coordinate L=log(x/(4*pi)); it is not "
            "introduced by a polygon, curvature image, or chosen circle."
        ),
        "scale_derivatives": {
            "partial_t_log_s_pf": "1/[2t(1+q)]",
            "partial_x_log_s_pf": "-2tL/[x(1+q)]",
            "s_pf_times_s_pf_t": "1/(1+q)^2",
            "s_pf_x_squared": "4q t^2/[x^2(1+q)^3]",
        },
        "bulk_coefficient": {
            "definition": "A_pf=1+s_pf s_pf,t-s_pf,x^2",
            "factorization": (
                "A_pf=1+[1-r(t,x)]/(1+q)^2, "
                "r(t,x)=4q t^2/[x^2(1+q)]"
            ),
            "outer_domain": "x>=38 and 0<t<=1/4",
            "ratio_chain": (
                "0<=r<=(4t^2/x^2)<="
                f"{ratio_cap.numerator}/{ratio_cap.denominator}"
            ),
            "uniform_bound": (
                "A_pf>=1+5775/[5776(1+q)^2]>1"
            ),
            "q_le_1": (
                "q<=1 implies A_pf>="
                f"{parabolic_coefficient.numerator}/"
                f"{parabolic_coefficient.denominator}>1"
            ),
            "q_ge_1": (
                "q>=1 implies A_pf>1, but the certified surplus "
                "above 1 tends to 0 as q tends to infinity"
            ),
        },
        "fixed_slice_balance": (
            "d/dt integral_a^b E_pf dx"
            "+2[J_pf(t,b)-J_pf(t,a)]"
            "=2 integral_a^b B_pf dx"
        ),
        "spacetime_balance": (
            "integral_a^b(E_pf(t_+,x)-E_pf(t_-,x))dx"
            "+2 integral_t_-^t_+[J_pf(t,b)-J_pf(t,a)]dt"
            "=2 integral_t_-^t_+ integral_a^b B_pf dxdt"
        ),
        "zero_flux_consequence": (
            "If J_pf(t,b)=J_pf(t,a), then the spatially integrated "
            "energy is nondecreasing in Newman time."
        ),
        "spatial_gradient_bridge": {
            "jet": "V_pf=(H,s_pf H_x), E_pf=|V_pf|^2",
            "gradient": (
                "partial_x V_pf=(H_x,s_pf H_xx+s_pf,x H_x)"
            ),
            "domination": (
                "B_pf=|partial_x V_pf|^2+(A_pf-1)H_x^2"
                ">=|partial_x V_pf|^2"
            ),
        },
        "contact_observability": (
            "If V_pf(t,x_0)=0 for some x_0 in [a,b], then "
            "integral_a^b B_pf dx>="
            "[sqrt(E_pf(t,a))+sqrt(E_pf(t,b))]^2/(b-a)"
        ),
        "conditional_no_contact": (
            "If both endpoint jets are certified and "
            "d/dt integral_a^b E_pf dx+2[J_pf(t,b)-J_pf(t,a)]"
            "<2[sqrt(E_pf(t,a))+sqrt(E_pf(t,b))]^2/(b-a), "
            "then V_pf(t,x) has no zero on [a,b]."
        ),
        "contact_observability_proof": (
            "At a contact x_0, Cauchy-Schwarz on [a,x_0] and "
            "[x_0,b] gives integral |V_x|^2>=E(a)/(x_0-a)"
            "+E(b)/(b-x_0)>="
            "[sqrt(E(a))+sqrt(E(b))]^2/(b-a), with endpoint "
            "cases read by continuity."
        ),
        "contact_cancellation": {
            "contact": "H=H_x=0",
            "values": (
                "E_pf=0, partial_t E_pf=0, J_pf=0, "
                "partial_x J_pf=s_pf^2 H_xx^2"
            ),
            "balance": (
                "2 partial_x J_pf=2s_pf^2 H_xx^2; the positive bulk "
                "is exactly canceled by the flux divergence"
            ),
        },
        "quadratic_countermodel": {
            "field": "H(t,x)=(x-x_0)^2-2(t-t_0)",
            "heat_equation": "H_t=-H_xx",
            "contact": "H(t_0,x_0)=H_x(t_0,x_0)=0 and H_xx=2",
            "constant_scale_values": (
                "E=E_t=J=0, J_x=4s_0^2, "
                "2J_x=8s_0^2=positive bulk"
            ),
            "meaning": (
                "A genuine backward-heat double contact obeys the "
                "coercive balance, so that balance alone cannot exclude "
                "a first-jet collision."
            ),
        },
        "nonpromotion": {
            "density_tautology": (
                "E_pf(t,x)>0 is exactly the desired pointwise "
                "noncontact statement, so it cannot be assumed."
            ),
            "integral_guard": (
                "A positive or monotone spatial integral of E_pf does "
                "not rule out E_pf=0 at an isolated point."
            ),
            "maximum_principle_guard": (
                "The flux contains H_xx and the local balance is not a "
                "closed scalar parabolic inequality for E_pf."
            ),
            "relative_quotient_guard": (
                "No kappa is defined by dividing through E_pf or "
                "||(H,s_pf H_x)||."
            ),
        },
        "conditional_utility": (
            "The exact contact-observability floor converts the multiplier "
            "into a pointwise criterion once Xi-specific endpoint margins "
            "and a strict upper bound for the integrated bulk or equivalent "
            "energy-current numerator are supplied."
        ),
        "route_decision": (
            "Retain the balance as an exact multiplier component, reject "
            "energy positivity alone as the descendant theorem, and "
            "return the pointwise step to the Xi contact-normal arithmetic "
            "hierarchy unless an independent trace/observability theorem "
            "is proved."
        ),
        "open_handoff": (
            "On C_j^out=[t_(j+1),t_j]x[38,R_j], derive source-level "
            "endpoint lower margins and prove the strict bulk-budget "
            "inequality D_j(t)<2[sqrt(E_pf(t,38))+"
            "sqrt(E_pf(t,R_j))]^2/(R_j-38), where "
            "D_j=d/dt integral_38^R_j E_pf dx"
            "+2[J_pf(t,R_j)-J_pf(t,38)]. Do not define D_j by a "
            "relative quotient or assume the desired contact exclusion."
        ),
        "symbolic_audit": symbolic,
    }


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    exact = build_exact()
    coefficient = exact["bulk_coefficient"]
    cancellation = exact["contact_cancellation"]
    countermodel = exact["quadratic_countermodel"]
    nonpromotion = exact["nonpromotion"]
    rows = [
        GateRow(
            id="pfec_01_heat_setup",
            role="exact_identity",
            readiness="available_exact",
            claim="The audit starts only from the real backward-heat equation and a positive C1 scale.",
            formula=exact["assumptions"],
            proof_boundary="No Xi-specific lower bound is imported.",
        ),
        GateRow(
            id="pfec_02_direct_derivative",
            role="exact_identity",
            readiness="available_exact",
            claim="Differentiating the scaled first-jet energy retains every scale-time term.",
            formula=exact["direct_time_derivative"],
            proof_boundary="Pointwise differential identity.",
        ),
        GateRow(
            id="pfec_03_flux_derivative",
            role="exact_identity",
            readiness="available_exact",
            claim="The third derivative is placed in one exact spatial flux.",
            formula=f"{exact['flux']}; {exact['flux_derivative']}",
            proof_boundary="No boundary sign is asserted.",
        ),
        GateRow(
            id="pfec_04_completed_balance",
            role="exact_identity",
            readiness="available_exact",
            claim="The remaining bulk completes to a square plus one first-derivative coefficient.",
            formula=exact["local_balance"],
            proof_boundary="Valid for every sufficiently smooth backward-heat field.",
        ),
        GateRow(
            id="pfec_05_scale_derivatives",
            role="exact_identity",
            readiness="available_exact",
            claim="The parabolic-frequency scale derivatives are exact in the common q coordinate.",
            formula="; ".join(exact["scale_derivatives"].values()),
            proof_boundary=exact["pi_provenance"],
        ),
        GateRow(
            id="pfec_06_coefficient_factorization",
            role="exact_identity",
            readiness="available_exact",
            claim="The only unsigned bulk coefficient has an exact positive-factor form.",
            formula=coefficient["factorization"],
            proof_boundary=coefficient["outer_domain"],
        ),
        GateRow(
            id="pfec_07_outer_ratio_bound",
            role="exact_inequality",
            readiness="proved_exact",
            claim="The spatial-scale correction is uniformly tiny on every ray-aligned outer collar.",
            formula=f"{coefficient['ratio_chain']}; {coefficient['uniform_bound']}",
            proof_boundary="Uses only x>=38, 0<t<=1/4, and q>=0.",
        ),
        GateRow(
            id="pfec_08_parabolic_chart",
            role="exact_inequality",
            readiness="proved_exact",
            claim="The q<=1 chart has a uniform coefficient surplus above one.",
            formula=coefficient["q_le_1"],
            proof_boundary="Coarse exact lower bound, not an Xi amplitude estimate.",
        ),
        GateRow(
            id="pfec_09_frequency_chart",
            role="exact_inequality",
            readiness="proved_exact",
            claim="The q>=1 chart remains strictly positive through the chart interface and at every finite q.",
            formula=coefficient["q_ge_1"],
            proof_boundary="No q-independent surplus above one is claimed as q tends to infinity.",
        ),
        GateRow(
            id="pfec_10_fixed_slice_balance",
            role="exact_identity",
            readiness="available_exact",
            claim="Spatial integration exposes only the two endpoint fluxes.",
            formula=exact["fixed_slice_balance"],
            proof_boundary="Fixed finite interval with sufficient trace regularity.",
        ),
        GateRow(
            id="pfec_11_spacetime_balance",
            role="exact_identity",
            readiness="available_exact",
            claim="Time integration gives an exact outer-collar multiplier identity.",
            formula=exact["spacetime_balance"],
            proof_boundary="The endpoint flux has no proved Xi sign or bound.",
        ),
        GateRow(
            id="pfec_12_spatial_gradient_bridge",
            role="exact_contact_lemma",
            readiness="available_exact",
            claim="The coercive bulk dominates the full spatial derivative of the scaled first jet.",
            formula=exact["spatial_gradient_bridge"]["domination"],
            proof_boundary="Exact one-dimensional identity; no endpoint margin is asserted.",
        ),
        GateRow(
            id="pfec_13_conditional_no_contact",
            role="conditional_criterion",
            readiness="ready_if_hypothesis_proved",
            claim="A contact forces an explicit bulk floor determined by the two endpoint jet energies.",
            formula=f"{exact['contact_observability']}; {exact['conditional_no_contact']}",
            proof_boundary="The required Xi endpoint lower margins and strict bulk upper bound are open.",
        ),
        GateRow(
            id="pfec_14_contact_cancellation",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="At a first-jet contact the positive curvature bulk is exactly balanced by flux divergence.",
            formula=f"{cancellation['values']}; {cancellation['balance']}",
            proof_boundary="Blocks a false pointwise inference from bulk positivity.",
        ),
        GateRow(
            id="pfec_15_quadratic_countermodel",
            role="exact_countermodel",
            readiness="guard_validated",
            claim="A shifted quadratic backward-heat solution realizes the contact cancellation exactly.",
            formula=f"{countermodel['field']}; {countermodel['constant_scale_values']}",
            proof_boundary="Generic heat countermodel, not an Xi counterexample.",
        ),
        GateRow(
            id="pfec_16_integral_nonpromotion",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Integrated energy control cannot be promoted to pointwise first-jet noncontact.",
            formula=(
                f"{nonpromotion['integral_guard']} "
                f"{nonpromotion['maximum_principle_guard']}"
            ),
            proof_boundary=nonpromotion["relative_quotient_guard"],
        ),
        GateRow(
            id="pfec_17_conditional_utility",
            role="route_component",
            readiness="available_exact",
            claim="The balance can contribute derivative control if Xi-specific flux and trace estimates are supplied.",
            formula=exact["conditional_utility"],
            proof_boundary="Still requires a separate pointwise conversion theorem.",
        ),
        GateRow(
            id="pfec_18_route_handoff",
            role="open_handoff",
            readiness="not_ready_to_apply",
            claim="The remaining useful target is a source-derived Xi trace/observability theorem coupled to the contact-normal margin.",
            formula=exact["open_handoff"],
            proof_boundary="Open Xi theorem; Q209 and the cofinal successor theorem remain unproved.",
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact parabolic-frequency energy-current balance with "
            "positive bulk and pointwise nonpromotion guard"
        ),
        "proof_boundary": (
            "This artifact proves the local and integrated multiplier "
            "identities, exact parabolic-frequency derivatives, and "
            "positive outer-domain bulk coefficient. It also proves that "
            "bulk positivity or integrated energy control alone does not "
            "exclude a first-jet contact. It does not prove Xi endpoint "
            "lower margins or the strict bulk-budget inequality, certify Q209, "
            "prove the ray-aligned descendant theorem, prove Lambda<=0, "
            "prove PF-infinity, prove RH, or establish a Clay-prize "
            "conclusion."
        ),
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "exact_identity_rows": sum(
                row.role == "exact_identity" for row in rows
            ),
            "coefficient_bound_rows": sum(
                row.role == "exact_inequality" for row in rows
            ),
            "pointwise_bridge_rows": sum(
                row.role in {"exact_contact_lemma", "conditional_criterion"}
                for row in rows
            ),
            "nonpromotion_guard_rows": sum(
                row.role in {"nonpromotion_gate", "exact_countermodel"}
                for row in rows
            ),
            "open_xi_targets": sum(
                row.role == "open_handoff" for row in rows
            ),
            "pointwise_contact_exclusion": False,
            "q209_certified": False,
            "cofinal_descendant_theorem": False,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    coefficient = exact["bulk_coefficient"]
    cancellation = exact["contact_cancellation"]
    countermodel = exact["quadratic_countermodel"]
    nonpromotion = exact["nonpromotion"]
    return "\n".join(
        [
            "# Parabolic-Frequency Energy/Current Gate",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact multiplier identity with positive bulk and an",
            "explicit pointwise nonpromotion guard. This is not a proof of",
            "Q209, the descendant theorem, Lambda<=0, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "Current result:",
            "",
            "```text",
            "validated Newman parabolic-frequency energy/current gate: "
            "18 rows, 8 exact identities, 3 coefficient bounds, "
            "2 pointwise bridge rows, 3 nonpromotion guards, "
            "1 open Xi bulk-margin target, "
            "0 pointwise contact exclusions",
            "```",
            "",
            "## Exact Local Balance",
            "",
            exact["assumptions"],
            "",
            "Set",
            "",
            "```text",
            exact["energy_density"],
            exact["flux"],
            "```",
            "",
            "Direct differentiation gives",
            "",
            "```text",
            exact["direct_time_derivative"],
            exact["flux_derivative"],
            "```",
            "",
            "and exact square completion gives",
            "",
            "```text",
            exact["local_balance"],
            "```",
            "",
            "No term involving `s_t` or `s_x` has been discarded.",
            "",
            "## Parabolic-Frequency Coefficient",
            "",
            "```text",
            exact["scale"],
            *[
                f"{key}={value}"
                for key, value in exact["scale_derivatives"].items()
            ],
            coefficient["factorization"],
            "```",
            "",
            "On every ray-aligned outer collar, `x>=38` and",
            "`0<t<=1/4`. Therefore",
            "",
            "```text",
            coefficient["ratio_chain"],
            coefficient["uniform_bound"],
            "```",
            "",
            "The requested two-chart audit is",
            "",
            "```text",
            coefficient["q_le_1"],
            coefficient["q_ge_1"],
            "```",
            "",
            "Thus the bulk is genuinely coercive in both charts. The",
            "frequency chart has no claimed uniform surplus above one as",
            "`q` tends to infinity.",
            "",
            "Pi provenance:",
            "",
            exact["pi_provenance"],
            "",
            "## Integrated Identities",
            "",
            "For a fixed finite interval `[a,b]`,",
            "",
            "```text",
            exact["fixed_slice_balance"],
            "```",
            "",
            "and on a time slab,",
            "",
            "```text",
            exact["spacetime_balance"],
            "```",
            "",
            exact["zero_flux_consequence"],
            "",
            "This is integrated derivative control, not a pointwise lower",
            "bound for `E_pf`.",
            "",
            "## Exact Pointwise Bridge",
            "",
            "The bulk has the additional exact interpretation",
            "",
            "```text",
            exact["spatial_gradient_bridge"]["jet"],
            exact["spatial_gradient_bridge"]["gradient"],
            exact["spatial_gradient_bridge"]["domination"],
            "```",
            "",
            "Consequently, a contact at any unknown `x_0` forces",
            "",
            "```text",
            exact["contact_observability"],
            "```",
            "",
            exact["contact_observability_proof"],
            "",
            "Combining this floor with the fixed-slice balance gives the",
            "strict conditional criterion",
            "",
            "```text",
            exact["conditional_no_contact"],
            "```",
            "",
            "This is a real pointwise bridge, but its Xi-specific strict",
            "upper bound and endpoint margins have not been proved.",
            "",
            "## Contact Cancellation",
            "",
            "At `H=H_x=0`,",
            "",
            "```text",
            cancellation["values"],
            cancellation["balance"],
            "```",
            "",
            "There is no contradiction: the positive curvature term is",
            "carried entirely by the local flux divergence.",
            "",
            "The exact shifted quadratic model makes this obstruction",
            "concrete:",
            "",
            "```text",
            countermodel["field"],
            countermodel["heat_equation"],
            countermodel["contact"],
            countermodel["constant_scale_values"],
            "```",
            "",
            countermodel["meaning"],
            "",
            "## Nonpromotion Decision",
            "",
            nonpromotion["density_tautology"],
            "",
            nonpromotion["integral_guard"],
            "",
            nonpromotion["maximum_principle_guard"],
            "",
            nonpromotion["relative_quotient_guard"],
            "",
            "The identity is still useful:",
            "",
            exact["conditional_utility"],
            "",
            "Route decision:",
            "",
            exact["route_decision"],
            "",
            "Open handoff:",
            "",
            exact["open_handoff"],
            "",
            "## Boundary",
            "",
            payload["proof_boundary"],
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
    args.out.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(payload), encoding="utf-8")
    print(
        "built Newman parabolic-frequency energy/current gate: "
        "18 rows, 8 exact identities, 3 coefficient bounds, "
        "2 pointwise bridge rows, 3 nonpromotion guards, "
        "1 open Xi bulk-margin target, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
