#!/usr/bin/env python3
"""Build the first-order corrected-main target reconciliation gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_"
    "first_order_target_reconciliation_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "zeroth_order_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_C1_global_remainder_certificate.json"
    ),
    "first_order_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_global_remainder_certificate.json"
    ),
    "first_order_signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_signed_contact_reduction.json"
    ),
    "wronskian_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_wronskian_phase_reduction.json"
    ),
    "component_wronskian": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_component_wronskian_gate.json"
    ),
    "absolute_phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_absolute_phase_anchor_reduction.json"
    ),
    "real_residual": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
    "abel_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_carrier_kernel_abel_prefix_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_abel_scalar_shear_flux_reduction.json"
    ),
    "oscillatory_splice": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_oscillatory_spliced_outer_collar_reduction.json"
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
    old = payloads["zeroth_order_remainder"]
    if old.get("exact", {}).get("remaining_target") != (
        "Prove T_L[J]>32000000*exp(-3L/2) "
        "for the corrected finite main"
    ):
        raise RuntimeError("zeroth-order target drifted")

    first = payloads["first_order_remainder"]
    first_exact = first.get("exact", {})
    if first_exact.get("global_remainder") != (
        "For every critical disk, including cutoff crossings, "
        "|r_[1](x)|<100000*exp(-5L/4) and "
        "|partial_x r_[1](x)|<200000*L*exp(-5L/4)"
    ):
        raise RuntimeError("first-order global remainder drifted")
    first_handoff = first_exact.get("contact_handoff", "")
    if (
        "exact signed contact-normal inequality" not in first_handoff
        or "strict Xi arithmetic reversal" not in first_handoff
    ):
        raise RuntimeError("first-order contact handoff drifted")

    signed = payloads["first_order_signed_contact"]
    signed_exact = signed.get("exact", {})
    if signed_exact.get("corrected_complex_main") != (
        "E_[1]=g_0+sum_(n=1)^N f_n, J_[1]=2*Re(E_[1])"
    ):
        raise RuntimeError("first-order corrected complex main drifted")
    if signed_exact.get("certified_contact_box") != (
        "|X_[1]|<50000*exp(-5L/4), "
        "|U_[1]|<100000*L*exp(-5L/4)"
    ):
        raise RuntimeError("first-order contact box drifted")
    if signed_exact.get("signed_band_target") != (
        "|X_[1]|<=50000*exp(-5L/4) implies "
        "|U_[1]|>100000*L*exp(-5L/4)"
    ):
        raise RuntimeError("first-order signed band drifted")
    if "adjacent-lift difference is retained in r_[1]" not in (
        signed_exact.get("local_lift", "")
    ):
        raise RuntimeError("first-order local-lift boundary drifted")

    phase = payloads["wronskian_phase"].get("exact", {})
    if phase.get("cartesian") != (
        "For E=X+iY and E'=U+iV, "
        "T_L[J]=4*(X^2+(U/L)^2)"
    ):
        raise RuntimeError("Cartesian phase energy drifted")

    component = payloads["component_wronskian"].get("exact", {})
    if "genuine double zero" not in (
        component.get("ordered_speed_countermodel", "")
    ):
        raise RuntimeError("ordered-speed nonpromotion guard drifted")

    anchor = payloads["absolute_phase_anchor"].get("exact", {})
    if "X/|f_1|=Re(eta*Z_0)" not in (
        anchor.get("real_projection", "")
    ):
        raise RuntimeError("absolute phase normalization drifted")

    residual = payloads["real_residual"].get("exact", {})
    contact_and_count = residual.get("contact_and_count", {})
    if contact_and_count.get("contact_box") != (
        "A full contact forces |X|<50000*exp(-5L/4) and "
        "|U|<100000L*exp(-5L/4)."
    ):
        raise RuntimeError("real-residual contact box drifted")
    if residual.get("core_scalar", {}).get("band_approximation") != (
        "For |X|<=50000*exp(-5L/4), "
        "|U-A_a|<1e-6*exp(-5L/4)."
    ):
        raise RuntimeError("real-residual band approximation drifted")

    prefix = payloads["abel_prefix"].get("exact", {})
    if "delta_L=50000*exp(-5L/4)/|f_1|" not in (
        prefix.get("live_target", "")
    ):
        raise RuntimeError("Abel-prefix delta drifted")
    if "A_L=(100000L+1)*exp(-5L/4)/|f_1|" not in (
        prefix.get("live_target", "")
    ):
        raise RuntimeError("Abel-prefix derivative margin drifted")
    if prefix.get("contact_scalar", "").find(
        "mathsf_A=Re(W_A)=mathcal_C_N+c*u_N*mathsf_X"
    ) < 0:
        raise RuntimeError("Abel terminal shear identity drifted")

    shear = payloads["abel_shear"].get("exact", {})
    if shear.get("q_ge_1_pointwise_target") != (
        "Pointwise Xi target: prove |mathsf_X|<=delta_L implies "
        "|mathcal_C_N|>A_L+epsilon_term in each prescribed canonical "
        "q>=1 chart, retaining W_0=0 and the endpoint."
    ):
        raise RuntimeError("Abel shear pointwise target drifted")

    splice = payloads["oscillatory_splice"].get("exact", {})
    if "0<tL<=c_*+epsilon" not in splice.get("open_outer_target", ""):
        raise RuntimeError("oscillatory-splice low-c target drifted")
    if "q=2tL^2>=1" not in splice.get("open_outer_target", ""):
        raise RuntimeError("oscillatory-splice q-domain drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCES.items()
        },
        "chronology": {
            "zeroth_order_date": old.get("date"),
            "first_order_date": first.get("date"),
            "signed_contact_date": signed.get("date"),
        },
        "zeroth_order_target": old["exact"]["remaining_target"],
        "first_order_remainder": first_exact["global_remainder"],
        "first_order_main": signed_exact["corrected_complex_main"],
        "first_order_contact_box": signed_exact["certified_contact_box"],
        "cofinal_cutoff_coverage": (
            "The first-order remainder theorem covers every critical disk, "
            "including adjacent-cutoff crossings."
        ),
        "inherited_boundary": (
            "The first-order theorem supersedes the older remainder scale. "
            "It does not prove the signed main-jet separation."
        ),
    }


def build_exact() -> dict:
    old_constant = 32_000_000
    new_constant = 50_000_000_000
    ratio = Fraction(new_constant, old_constant)
    if ratio != Fraction(3125, 2):
        raise RuntimeError("old/new threshold ratio failed")

    delta_constant = 50_000
    gamma_constant = 100_000
    radial_constant = 4 * (
        delta_constant**2 + gamma_constant**2
    )
    if radial_constant != new_constant:
        raise RuntimeError("first-order radial threshold failed")

    return {
        "source_reconciliation": {
            "historical_main": "J_[0]=P_[0]-Q_[0]",
            "historical_target": (
                "T_L[J_[0]]>32000000*exp(-3L/2)"
            ),
            "current_main": (
                "E_[1]=g_0+sum_(n=1)^N f_n, "
                "J_[1]=2*Re(E_[1])"
            ),
            "current_remainder": (
                "|r_[1]|<100000*exp(-5L/4), "
                "|r_[1],x|<200000*L*exp(-5L/4)"
            ),
            "supersession": (
                "The first-order corrected main and exp(-5L/4) C1 "
                "remainder supersede the older exp(-3L/4) handoff on "
                "L>=50, 0<tL<25."
            ),
        },
        "local_lift": {
            "rule": (
                "On each critical radius-1/L disk choose one prescribed-N "
                "analytic lift E_[1],N; at a cutoff, the adjacent-lift "
                "difference is retained in r_[1]."
            ),
            "physical_split": (
                "Z_t=J_[1],N+r_[1],N on the real axis"
            ),
            "pointwise_quantifier": (
                "At each point, one certified local lift with its matched "
                "remainder is enough for pointwise contact exclusion."
            ),
            "degree_quantifier": (
                "Boundary-degree composition must include equality "
                "boundaries and every prescribed adjacent chart/homotopy."
            ),
        },
        "division_free_components": {
            "components": (
                "z_0=g_0; z_n=f_n for 1<=n<=N; "
                "c_j=Re(z_j), d_j=Re(z_(j,x))"
            ),
            "sums": (
                "X=sum_(j=0)^N c_j=Re(E_[1]), "
                "U=sum_(j=0)^N d_j=Re(E_[1],x)"
            ),
            "first_jet": (
                "J_[1]=2X, J_[1],x=2U"
            ),
            "energy": (
                "T_L[J_[1]]=4[X^2+(U/L)^2]"
            ),
        },
        "phase_amplitude": {
            "component": (
                "For z_j=a_j*exp(i*theta_j)!=0 and "
                "z_(j,x)=(u_j+i*v_j)z_j, "
                "c_j=a_j*cos(theta_j), "
                "d_j=a_j[u_j*cos(theta_j)-v_j*sin(theta_j)]."
            ),
            "sum": (
                "T_L[J_[1]]/4="
                "[sum_j a_j cos(theta_j)]^2"
                "+L^-2[sum_j a_j"
                "{u_j cos(theta_j)-v_j sin(theta_j)}]^2"
            ),
            "zero_component_guard": (
                "The Cartesian c_j,d_j formula is primary when an endpoint "
                "component vanishes; no phase or rate quotient is taken."
            ),
        },
        "cross_term_expansion": {
            "vectors": "w_j=(c_j,d_j/L) in R^2",
            "gram": (
                "T_L[J_[1]]/4=|sum_j w_j|^2"
                "=sum_j|w_j|^2+2sum_(j<k) w_j dot w_k"
            ),
            "rank": (
                "The component Gram form has rank at most two, independent "
                "of the number of components."
            ),
            "cancellation_guard": (
                "w_1=(A,0), w_2=(-A,0) gives "
                "sum_j|w_j|^2=2A^2 but T_L[J_[1]]=0."
            ),
            "boundary": (
                "This is an algebraic nonpromotion guard, not an actual Xi "
                "configuration."
            ),
        },
        "contact_box": {
            "value_halfwidth": (
                "delta_0(L)=50000*exp(-5L/4)"
            ),
            "derivative_halfwidth": (
                "gamma_0(L)=100000*L*exp(-5L/4)"
            ),
            "implication": (
                "H=H_x=0 => |X|<delta_0(L) and |U|<gamma_0(L)"
            ),
            "rectangle": (
                "R_L={|X|<=delta_0(L), |U|<=gamma_0(L)}"
            ),
        },
        "direct_radial_target": {
            "contact_upper": (
                "At contact, T_L[J_[1]]"
                "<50000000000*exp(-5L/2)."
            ),
            "sufficient": (
                "T_L[J_[1]]>=50000000000*exp(-5L/2) "
                "excludes contact."
            ),
            "constant_identity": (
                "4*(50000^2+100000^2)=50000000000"
            ),
            "comparison": (
                "new_threshold/old_threshold"
                "=(3125/2)*exp(-L)<1 for L>=50"
            ),
            "meaning": (
                "The first-order target is exponentially smaller and hence "
                "strictly easier as a lower-bound obligation."
            ),
        },
        "box_optimal_target": {
            "gauge": (
                "G_L(X,U)=max(|X|/delta_0(L),|U|/gamma_0(L))"
            ),
            "contact": "A contact forces G_L(X,U)<1.",
            "exact_target": (
                "G_L(X,U)>1, equivalently "
                "|X|<=delta_0(L) => |U|>gamma_0(L)"
            ),
            "minimality": (
                "This excludes exactly the certified axis-aligned contact "
                "box; it does not replace it by a larger ellipse."
            ),
            "radial_overstrength": (
                "In normalized units delta=1, gamma=2, "
                "(X,U)=(0,21/10) has G=21/20>1 but "
                "X^2+U^2=441/100<5=delta^2+gamma^2."
            ),
        },
        "absolute_phase_normalization": {
            "anchor": (
                "W_0=E_[1]/|f_1|=mathsf_X+i*mathsf_Y"
            ),
            "value": "mathsf_X=X/|f_1|",
            "centered_derivative": "mathsf_A=A_a/|f_1|",
            "bands": (
                "delta_L=50000*exp(-5L/4)/|f_1|, "
                "A_L=(100000L+1)*exp(-5L/4)/|f_1|"
            ),
            "real_residual": (
                "|U-A_a|<1e-6*exp(-5L/4) "
                "on |X|<=50000*exp(-5L/4)"
            ),
        },
        "abel_target": {
            "shear": (
                "mathsf_A=mathcal_C_N+c*u_N*mathsf_X"
            ),
            "terminal_error": (
                "epsilon_term=25000*exp(-11L/4)/|f_1|"
            ),
            "pointwise": (
                "|mathsf_X|<=delta_L => "
                "|mathcal_C_N|>A_L+epsilon_term"
            ),
            "composition": (
                "The Abel target gives |mathsf_A|>A_L, hence "
                "|A_a|>(100000L+1)exp(-5L/4); the real-residual bound then "
                "gives |U|>100000L exp(-5L/4)."
            ),
            "conclusion": (
                "The division-free Abel gap implies the box-optimal signed "
                "band and therefore excludes a first-order contact."
            ),
        },
        "outer_domain_reduction": {
            "threshold": "c_*=4911678521/1933561194",
            "fixed_epsilon": "0<epsilon<25-c_*",
            "height": "B_epsilon=max(50,L_epsilon)",
            "live_outer_wedge": (
                "L>=B_epsilon, q=2tL^2>=1, "
                "0<tL<=c_*+epsilon"
            ),
            "closed_middle": (
                "c_*+epsilon<tL<=25 is closed by the oscillatory-zeta "
                "theorem beyond L_epsilon."
            ),
            "closed_dominant": (
                "tL>=25 is closed by the dominant-saddle theorem."
            ),
            "inner_open": (
                "q<1, the bounded-L shoulder, finite joins, and the "
                "multiplicity-compatible inner degree theorem remain "
                "separate."
            ),
        },
        "nonpromotion": {
            "diagonal": (
                "Positive diagonal component energies do not lower-bound "
                "the two-dimensional sum because the cross terms can "
                "cancel them exactly."
            ),
            "ordered_phase": (
                "Positive component amplitudes and ordered phase speeds do "
                "not exclude a double real crossing; the checked "
                "four-component countermodel already realizes one."
            ),
            "finite": (
                "A finite phase plot or sampled minimum cannot prove the "
                "uniform low-c, cutoff-complete Abel gap."
            ),
        },
        "route_decision": (
            "Retire the old exp(-3L/2) target and do not attack the stronger "
            "radial first-order target by diagonal positivity. The exact "
            "main-only obligation is rectangle avoidance; the current "
            "proof-facing sufficient theorem is its robust division-free "
            "Abel form. Restrict that arithmetic work to the oscillatory-"
            "spliced low-c q>=1 wedge, and keep q<1/bounded-L degree closure "
            "as a separate theorem."
        ),
        "open_handoff": (
            "On one fixed epsilon collar, derive a source-specific lower "
            "bound for the endpoint-complete Abel scalar mathcal_C_N in the "
            "band |mathsf_X|<=delta_L for "
            "L>=B_epsilon, q>=1, and 0<tL<=c_*+epsilon. Preserve W_0=0, "
            "the recurrent endpoint, q=1, equality boundaries, and "
            "adjacent-cutoff homotopies. Start from the exact prefix sum "
            "U=u_N*S+sum_k log((k+1)/k)F_k and test signed cancellation "
            "before absolute values."
        ),
    }


def build_payload() -> dict:
    payloads = load_sources()
    source_audit = audit_sources(payloads)
    exact = build_exact()
    rows = [
        GateRow(
            id="fotr_01_chronology",
            role="source_reconciliation",
            readiness="proved_exact",
            claim="The old C1 target predates and is superseded by the first-order global remainder theorem.",
            formula=exact["source_reconciliation"]["supersession"],
            proof_boundary="This changes the active handoff, not the validity of the older conditional lemma.",
        ),
        GateRow(
            id="fotr_02_current_main",
            role="exact_definition",
            readiness="proved_exact",
            claim="The live finite main retains the first Dirichlet and endpoint corrections.",
            formula=exact["source_reconciliation"]["current_main"],
            proof_boundary="The endpoint component g_0 is retained.",
        ),
        GateRow(
            id="fotr_03_current_remainder",
            role="source_audit",
            readiness="certified",
            claim="The first-order C1 remainder is globally controlled across cutoff crossings.",
            formula=exact["source_reconciliation"]["current_remainder"],
            proof_boundary="Valid on L>=50 and 0<tL<25.",
        ),
        GateRow(
            id="fotr_04_local_lift",
            role="exact_chart_contract",
            readiness="proved_exact",
            claim="The corrected main is a cutoff-local analytic lift with its chart difference retained in the remainder.",
            formula=exact["local_lift"]["rule"],
            proof_boundary=exact["local_lift"]["physical_split"],
        ),
        GateRow(
            id="fotr_05_cutoff_quantifiers",
            role="exact_chart_contract",
            readiness="proved_exact",
            claim="Pointwise exclusion and boundary-degree composition have different chart quantifiers.",
            formula=(
                f"{exact['local_lift']['pointwise_quantifier']} "
                f"{exact['local_lift']['degree_quantifier']}"
            ),
            proof_boundary="No global complex-main chart invariance is asserted.",
        ),
        GateRow(
            id="fotr_06_cartesian_energy",
            role="exact_identity",
            readiness="proved_exact",
            claim="The first-order main jet has a division-free Cartesian energy.",
            formula=(
                f"{exact['division_free_components']['sums']}; "
                f"{exact['division_free_components']['energy']}"
            ),
            proof_boundary=exact["division_free_components"]["components"],
        ),
        GateRow(
            id="fotr_07_phase_amplitude",
            role="exact_identity",
            readiness="proved_exact",
            claim="Every nonzero component has an exact phase/amplitude contribution to the two real jet coordinates.",
            formula=exact["phase_amplitude"]["component"],
            proof_boundary=exact["phase_amplitude"]["zero_component_guard"],
        ),
        GateRow(
            id="fotr_08_phase_energy",
            role="exact_identity",
            readiness="proved_exact",
            claim="The full phase energy is the square of two linked signed component sums.",
            formula=exact["phase_amplitude"]["sum"],
            proof_boundary="All endpoint and first-correction components remain in the sums.",
        ),
        GateRow(
            id="fotr_09_cross_terms",
            role="exact_identity",
            readiness="proved_exact",
            claim="The phase energy expands into diagonal masses plus signed cross terms.",
            formula=exact["cross_term_expansion"]["gram"],
            proof_boundary=exact["cross_term_expansion"]["vectors"],
        ),
        GateRow(
            id="fotr_10_rank_cancellation_guard",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="Componentwise positive energy cannot provide the required aggregate lower bound.",
            formula=exact["cross_term_expansion"]["cancellation_guard"],
            proof_boundary=exact["cross_term_expansion"]["boundary"],
        ),
        GateRow(
            id="fotr_11_contact_box",
            role="exact_composition",
            readiness="available_exact",
            claim="An exact Xi contact places the first-order main jet inside one explicit rectangle.",
            formula=exact["contact_box"]["implication"],
            proof_boundary=exact["contact_box"]["rectangle"],
        ),
        GateRow(
            id="fotr_12_radial_target",
            role="conditional_criterion",
            readiness="not_ready_to_apply",
            claim="One radial first-order energy floor is sufficient for contact exclusion.",
            formula=exact["direct_radial_target"]["sufficient"],
            proof_boundary=exact["direct_radial_target"]["constant_identity"],
        ),
        GateRow(
            id="fotr_13_threshold_improvement",
            role="exact_asymptotic",
            readiness="proved_exact",
            claim="The first-order direct threshold is exponentially smaller than the historical threshold.",
            formula=exact["direct_radial_target"]["comparison"],
            proof_boundary=exact["direct_radial_target"]["meaning"],
        ),
        GateRow(
            id="fotr_14_box_gauge",
            role="exact_reduction",
            readiness="proved_exact",
            claim="The contact rectangle has an exact anisotropic maximum-norm gauge.",
            formula=(
                f"{exact['box_optimal_target']['gauge']}; "
                f"{exact['box_optimal_target']['contact']}"
            ),
            proof_boundary=exact["box_optimal_target"]["minimality"],
        ),
        GateRow(
            id="fotr_15_signed_band",
            role="conditional_criterion",
            readiness="not_ready_to_apply",
            claim="Rectangle avoidance is exactly the signed first-order band target.",
            formula=exact["box_optimal_target"]["exact_target"],
            proof_boundary="This is weaker than the radial sufficient target.",
        ),
        GateRow(
            id="fotr_16_radial_overstrength",
            role="nonpromotion_gate",
            readiness="guard_validated",
            claim="The radial target strictly overcovers the certified contact rectangle.",
            formula=exact["box_optimal_target"]["radial_overstrength"],
            proof_boundary="An exact normalized counterexample to the reverse implication.",
        ),
        GateRow(
            id="fotr_17_absolute_anchor",
            role="exact_normalization",
            readiness="proved_exact",
            claim="The first carrier supplies a branch-free normalization of the value and centered derivative.",
            formula=(
                f"{exact['absolute_phase_normalization']['anchor']}; "
                f"{exact['absolute_phase_normalization']['value']}; "
                f"{exact['absolute_phase_normalization']['centered_derivative']}"
            ),
            proof_boundary="No division by E_[1] or its first jet is used.",
        ),
        GateRow(
            id="fotr_18_abel_shear",
            role="exact_identity",
            readiness="proved_exact",
            claim="The endpoint-complete Abel scalar differs from the centered derivative by one controlled terminal shear.",
            formula=(
                f"{exact['abel_target']['shear']}; "
                f"{exact['abel_target']['terminal_error']}"
            ),
            proof_boundary="The identity includes W_0=0 and the recurrent endpoint.",
        ),
        GateRow(
            id="fotr_19_abel_implies_band",
            role="exact_composition",
            readiness="available_exact",
            claim="The robust division-free Abel gap implies the box-optimal signed band.",
            formula=exact["abel_target"]["composition"],
            proof_boundary=exact["abel_target"]["conclusion"],
        ),
        GateRow(
            id="fotr_20_low_c_splice",
            role="exact_domain_reduction",
            readiness="proved_exact",
            claim="Existing theorems reduce the asymptotic outer arithmetic burden to the low-c q>=1 wedge.",
            formula=exact["outer_domain_reduction"]["live_outer_wedge"],
            proof_boundary=(
                f"{exact['outer_domain_reduction']['closed_middle']} "
                f"{exact['outer_domain_reduction']['closed_dominant']}"
            ),
        ),
        GateRow(
            id="fotr_21_inner_guard",
            role="route_guard",
            readiness="guard_validated",
            claim="The q<1 and bounded-L inner theorem remains a separate multiplicity-compatible obligation.",
            formula=exact["outer_domain_reduction"]["inner_open"],
            proof_boundary="The outer Abel target does not close the inner successor degree.",
        ),
        GateRow(
            id="fotr_22_route_decision",
            role="route_decision",
            readiness="available_exact",
            claim="Use the robust Abel rectangle-avoidance target only on the oscillatory-spliced low-c wedge.",
            formula=exact["open_handoff"],
            proof_boundary=exact["route_decision"],
        ),
    ]
    return {
        "kind": STEM,
        "date": "2026-07-29",
        "status": (
            "exact first-order corrected-main target reconciliation, "
            "box-optimal contact reduction, and low-c Abel handoff"
        ),
        "proof_boundary": (
            "This artifact reconciles the zeroth- and first-order targets, "
            "proves the cutoff-local Cartesian and phase/amplitude energy "
            "identities, derives the first-order contact box, radial "
            "sufficient threshold, box-optimal signed criterion, and exact "
            "Abel-to-band implication, and restricts the outer arithmetic "
            "target to the oscillatory-spliced low-c wedge. It does not "
            "prove the Xi Abel gap, a one-turn phase budget, the q<1 or "
            "bounded-L inner theorem, contact exclusion, Q209, the cofinal "
            "descendant theorem, Lambda<=0, PF-infinity, RH, or a Clay-prize "
            "conclusion."
        ),
        "source_audit": source_audit,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "rows": len(rows),
            "superseded_handoffs": 1,
            "cutoff_uniform_first_order_remainders": 1,
            "exact_energy_identities": 4,
            "contact_box_coordinates": 2,
            "direct_radial_targets": 1,
            "box_optimal_targets": 1,
            "normalized_abel_targets": 1,
            "nonpromotion_guards": 2,
            "open_outer_abel_gaps": 1,
            "pointwise_contact_exclusions": 0,
            "cofinal_descendant_theorem": False,
        },
    }


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    return "\n".join(
        [
            "# First-Order Target Reconciliation Gate",
            "",
            "Date: 2026-07-29",
            "",
            "Status: exact target reconciliation and route reduction. This",
            "is not a proof of contact exclusion, Lambda<=0, or RH.",
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
            "validated Newman first-order target reconciliation gate: "
            "22 rows, 1 superseded handoff, "
            "1 cutoff-uniform first-order remainder, "
            "4 exact energy identities, 2 contact-box coordinates, "
            "1 radial target, 1 box-optimal target, "
            "1 normalized Abel target, 2 nonpromotion guards, "
            "0 pointwise contact exclusions",
            "```",
            "",
            "## Supersession",
            "",
            "The older handoff used",
            "",
            "```text",
            exact["source_reconciliation"]["historical_target"],
            "```",
            "",
            "The later cutoff-uniform theorem instead gives",
            "",
            "```text",
            exact["source_reconciliation"]["current_main"],
            exact["source_reconciliation"]["current_remainder"],
            "```",
            "",
            exact["source_reconciliation"]["supersession"],
            "",
            "## Cutoff-Local Main",
            "",
            exact["local_lift"]["rule"],
            "",
            "```text",
            exact["local_lift"]["physical_split"],
            "```",
            "",
            exact["local_lift"]["pointwise_quantifier"],
            "",
            exact["local_lift"]["degree_quantifier"],
            "",
            "## Phase Energy",
            "",
            "With the endpoint and every corrected Dirichlet component",
            "retained,",
            "",
            "```text",
            exact["division_free_components"]["components"],
            exact["division_free_components"]["sums"],
            exact["division_free_components"]["energy"],
            "```",
            "",
            "For each nonzero component,",
            "",
            "```text",
            exact["phase_amplitude"]["component"],
            exact["phase_amplitude"]["sum"],
            "```",
            "",
            "The division-free Cartesian formula remains primary when a",
            "component vanishes.",
            "",
            "The cross-term structure is",
            "",
            "```text",
            exact["cross_term_expansion"]["gram"],
            "```",
            "",
            exact["cross_term_expansion"]["cancellation_guard"],
            "",
            exact["cross_term_expansion"]["boundary"],
            "",
            "## Contact Box",
            "",
            "A contact forces",
            "",
            "```text",
            exact["contact_box"]["implication"],
            "```",
            "",
            "One radial sufficient target is",
            "",
            "```text",
            exact["direct_radial_target"]["sufficient"],
            exact["direct_radial_target"]["comparison"],
            "```",
            "",
            "But the exact remainder information is rectangular. Define",
            "",
            "```text",
            exact["box_optimal_target"]["gauge"],
            exact["box_optimal_target"]["exact_target"],
            "```",
            "",
            exact["box_optimal_target"]["minimality"],
            "",
            "The radial condition is strictly stronger:",
            "",
            "```text",
            exact["box_optimal_target"]["radial_overstrength"],
            "```",
            "",
            "## Abel Target",
            "",
            "The branch-free normalization and terminal shear are",
            "",
            "```text",
            exact["absolute_phase_normalization"]["anchor"],
            exact["absolute_phase_normalization"]["value"],
            exact["absolute_phase_normalization"]["bands"],
            exact["abel_target"]["shear"],
            exact["abel_target"]["terminal_error"],
            "```",
            "",
            "The proof-facing target is",
            "",
            "```text",
            exact["abel_target"]["pointwise"],
            "```",
            "",
            exact["abel_target"]["composition"],
            "",
            exact["abel_target"]["conclusion"],
            "",
            "## Reduced Domain",
            "",
            "For one fixed epsilon, the remaining outer wedge is",
            "",
            "```text",
            exact["outer_domain_reduction"]["fixed_epsilon"],
            exact["outer_domain_reduction"]["height"],
            exact["outer_domain_reduction"]["live_outer_wedge"],
            "```",
            "",
            exact["outer_domain_reduction"]["closed_middle"],
            "",
            exact["outer_domain_reduction"]["closed_dominant"],
            "",
            exact["outer_domain_reduction"]["inner_open"],
            "",
            "## Route Decision",
            "",
            exact["route_decision"],
            "",
            "Open handoff:",
            "",
            exact["open_handoff"],
            "",
            "Pi provenance is unchanged: `pi` is inherited from the",
            "completed-zeta normalization and the Riemann-Siegel saddle.",
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
        "built Newman first-order target reconciliation gate: "
        "22 rows, 1 superseded handoff, "
        "1 cutoff-uniform first-order remainder, "
        "4 exact energy identities, 2 contact-box coordinates, "
        "1 radial target, 1 box-optimal target, "
        "1 normalized Abel target, 2 nonpromotion guards, "
        "0 pointwise contact exclusions"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
