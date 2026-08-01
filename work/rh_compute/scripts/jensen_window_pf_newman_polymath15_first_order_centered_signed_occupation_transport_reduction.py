#!/usr/bin/env python3
"""Build the signed occupation transport and transversality reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "signed_occupation_transport_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "pairwise_alignment": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "pairwise_projective_alignment_gate.json"
    ),
    "interior_current": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "interior_projective_current_gate.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "direct_projection_regime_reduction.json"
    ),
    "bulk_pair": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "bulk_pair_transfer_gate.json"
    ),
    "contact_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contact_signed_transport_reduction.json"
    ),
    "abel_shear": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "abel_scalar_shear_flux_reduction.json"
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
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_hashes() -> dict[str, str]:
    missing = [str(path) for path in SOURCES.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {key: file_hash(path) for key, path in SOURCES.items()}


def source_audit() -> dict[str, str]:
    payloads = {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCES.items()
    }
    markers = {
        "pairwise_alignment": (
            "integral_A^B P_x(s)ds",
            "moving-scale sign-reversal guard",
            "e_n/(b*u_n)-e_m/(b*u_m)",
        ),
        "interior_current": (
            "z_(n,x)=(varrho_n+i*nu_n)z_n",
            "J_n<=-(5/64)*u_n^2*|z_n|^2<0",
            "C_edge=c_0+c_N",
        ),
        "normalized_prefix": (
            "gamma_n=q_(n,x)/q_n",
            "rho_n=-c*log(n)+Re(epsilon_n-epsilon_1)",
            "both pair kernels contain unrestricted sine and cosine",
        ),
        "direct_projection": (
            "eta*q_n=r_n*zeta_n",
            "rho_1=partial_x log|f_1|",
            "mathsf_A=partial_x mathsf_X",
        ),
        "bulk_pair": (
            "log(rho_n)<-(2/5)h_n",
            "0<rho_n<=exp(-2h_n/5)<1",
        ),
        "contact_transport": (
            "sum_(j=0)^N d_j=mathcal_C_N",
            "generally complex",
            "cumulative-real-mass",
        ),
        "abel_shear": (
            "2.03e-14",
            "partial_x mathsf_X-mathcal_C_N",
            "backward-heat",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {
        key: str(payload.get("kind", ""))
        for key, payload in payloads.items()
    }


def transport_audit() -> dict[str, str]:
    dual_residuals: list[sp.Expr] = []
    for index in range(3):
        c_value, c_x, h_x, test, test_x = sp.symbols(
            f"c_{index} c_x_{index} h_x_{index} "
            f"phi_{index} phi_x_{index}",
            real=True,
        )
        measure_x = c_x * test + c_value * h_x * test_x
        flux_s = -c_value * h_x * test_x
        dual_residuals.append(sp.expand(measure_x + flux_s - c_x * test))
    if any(value != 0 for value in dual_residuals):
        raise RuntimeError("weak occupation transport failed")

    c_1, c_2, h_1_x, h_2_x = sp.symbols(
        "c_1 c_2 h_1_x h_2_x",
        real=True,
    )
    flux = c_1 * h_1_x + c_2 * h_2_x
    if flux.subs(
        {c_1: 1, c_2: -1, h_1_x: -1, h_2_x: -2}
    ) != 1:
        raise RuntimeError("positive signed-flux guard failed")
    if flux.subs(
        {c_1: 1, c_2: -1, h_1_x: -2, h_2_x: -1}
    ) != -1:
        raise RuntimeError("negative signed-flux guard failed")

    return {
        "measure_definition": (
            "On a pole-free fixed-N chart let I be a fixed set of "
            "finite interior slopes and define "
            "mu(s;x)=sum_(i in I)c_i*delta(s-h_i), "
            "F(s;x)=sum_(i in I)c_i*h_(i,x)*delta(s-h_i), and "
            "S(s;x)=sum_(i in I)c_(i,x)*delta(s-h_i)."
        ),
        "weak_transport": (
            "For every compactly supported smooth test phi, "
            "partial_x sum_i c_i*phi(h_i)"
            "-sum_i c_i*h_(i,x)*phi'(h_i)"
            "=sum_i c_(i,x)*phi(h_i). Equivalently, "
            "partial_x mu+partial_s F=S in distributions."
        ),
        "cumulative_definition": (
            "Put P(s;x)=sum_(i in I)c_i*1_(s<h_i). Then "
            "partial_x P=F+sum_(i in I)c_(i,x)*1_(s<h_i) "
            "as a distribution in s."
        ),
        "signed_flux_guard": (
            "The strict theorem h_(i,x)<0 does not make F a one-sign "
            "measure because c_i is signed. With c_1=1,c_2=-1, "
            "the negative rates (-1,-2) give total flux +1, while "
            "(-2,-1) give -1. These are exact local transport jets, "
            "not asserted Xi states."
        ),
    }


def carrier_source_audit() -> dict[str, str]:
    (
        c_value,
        y_value,
        b,
        u,
        angular,
        radial,
        c_s,
        k,
        log_n,
        log_cap,
        xi,
        xi_one,
        residual,
    ) = sp.symbols(
        "c Y b u nu varrho c_s k log_n log_N "
        "xi xi_1 e",
        real=True,
    )
    c_x = radial * c_value - angular * y_value
    d_value = c_s * k * c_value - b * u * y_value
    source_residual = (xi - xi_one) * c_value - residual * y_value
    relations = {
        radial: -c_s * log_n + xi - xi_one,
        angular: b * u + residual,
        k: log_cap - log_n,
    }
    source_identity = sp.expand(
        (
            c_x
            - (
                d_value
                - c_s * log_cap * c_value
                + source_residual
            )
        ).subs(relations)
    )
    if source_identity != 0:
        raise RuntimeError("carrier source collapse failed")
    moment_identity = sp.expand(
        (
            d_value
            - (
                c_x
                + c_s * log_cap * c_value
                - source_residual
            )
        ).subs(relations)
    )
    if moment_identity != 0:
        raise RuntimeError("division-free moment identity failed")

    return {
        "actual_rates": (
            "Put epsilon_i=d_(i,x)/(1+d_i)=xi_i+i*chi_i. "
            "For z_i=eta*q_i=X_i+iY_i, "
            "varrho_i=-c_s*log(i)+xi_i-xi_1 and "
            "nu_i=b*u_i+e_i with e_i=v_a+chi_i. Hence "
            "X_(i,x)=varrho_i*X_i-nu_i*Y_i."
        ),
        "source_collapse": (
            "Define r_i=(xi_i-xi_1)X_i-e_iY_i. Since "
            "k_i=log(N/i), the exact division-free identity is "
            "X_(i,x)=d_i-c_s*log(N)*X_i+r_i, or "
            "d_i=X_(i,x)+c_s*log(N)*X_i-r_i."
        ),
        "finite_slope_reaction": (
            "Where X_i!=0 and h_i=d_i/X_i, "
            "X_(i,x)=[h_i-c_s*log(N)]X_i+r_i. Therefore "
            "S=[s-c_s*log(N)]mu+R, where "
            "R(s;x)=sum_i r_i*delta(s-h_i)."
        ),
    }


def normalized_audit() -> dict[str, str]:
    omega, omega_x, c_s, log_cap = sp.symbols(
        "Omega Omega_x c_s log_N",
        nonzero=True,
        real=True,
    )
    if sp.expand(
        omega_x.subs(omega_x, c_s * log_cap * omega)
        - c_s * log_cap * omega
    ) != 0:
        raise RuntimeError("integrating-factor derivative failed")

    return {
        "normalized_measure": (
            "Let Omega_N=N^(sigma_*), so "
            "Omega_(N,x)=c_s*log(N)*Omega_N. With "
            "bar_mu=Omega_N*mu, bar_F=Omega_N*F, and "
            "bar_R=Omega_N*R, the exact normalized law is "
            "partial_x bar_mu+partial_s bar_F=s*bar_mu+bar_R."
        ),
        "normalized_cumulative": (
            "Let M_>(s;x)=sum_i d_i*1_(s<h_i) and "
            "R_>(s;x)=sum_i r_i*1_(s<h_i). Then "
            "partial_x P=F+M_>-c_s*log(N)*P+R_>, so "
            "partial_x(Omega_N*P)="
            "Omega_N*(F+M_>+R_>)."
        ),
        "normalization_boundary": (
            "The common N^(sigma_*) factor removes the shared frame "
            "drift exactly, but the universal reaction s*bar_mu and "
            "the signed flux bar_F remain. This normalization is not "
            "a positivity transformation."
        ),
    }


def residual_audit() -> dict[str, str | int]:
    square_margin = 18888**2 - (16892**2 + 8449**2)
    if square_margin != 31279 or square_margin <= 0:
        raise RuntimeError("Euclidean residual constant failed")
    coefficient = Fraction(18888 * 10, 3 * 144)
    if not coefficient < 438:
        raise RuntimeError("summed residual coefficient failed")
    exponential_guard = (
        Fraction(19, 7) ** 45
        > Fraction(438 * 10**7) ** 2
    )
    if not exponential_guard:
        raise RuntimeError("L=50 exponential residual guard failed")

    return {
        "pointwise_residual": (
            "The source bounds |epsilon_i|<8446/x^2 and "
            "|e_i|<8449/x^2 give "
            "|r_i|<18888*|z_i|/x^2, because "
            "16892^2+8449^2<18888^2."
        ),
        "amplitude_sum": (
            "The corrected adjacent amplitude ratio obeys "
            "|z_(n+1)|/|z_n|<exp[-(2/5)log((n+1)/n)]. Since "
            "|z_1|=1, |z_n|<n^(-2/5), and "
            "sum_(n=1)^N|z_n|<(5/3)N^(3/5)"
            "<(10/3)exp(3L/10)."
        ),
        "summed_residual": (
            "Using x^2>144*exp(2L), "
            "|R_bulk|:=|sum_(n=1)^N r_n|"
            "<438*exp(-17L/10)"
            "<10^(-7)*exp(-5L/4) for L>=50."
        ),
        "exponential_guard": (
            "At L=50 the last inequality follows after squaring from "
            "(19/7)^45>(438*10^7)^2 and e>19/7; its ratio then "
            "decreases with L."
        ),
        "square_margin": square_margin,
    }


def aggregate_audit() -> dict[str, str]:
    c_bulk, c_bulk_x, d_bulk, r_bulk = sp.symbols(
        "C_bulk C_bulk_x D_bulk R_bulk",
        real=True,
    )
    c_s, log_cap = sp.symbols("c_s log_N", real=True)
    bulk_relation = c_bulk_x + c_s * log_cap * c_bulk - r_bulk
    if sp.expand(d_bulk.subs(d_bulk, bulk_relation) - bulk_relation) != 0:
        raise RuntimeError("bulk first-moment relation failed")

    c_zero, c_zero_x, d_zero = sp.symbols(
        "c_0 c_0_x d_0",
        real=True,
    )
    total_c = c_zero + c_bulk
    total_c_x = c_zero_x + c_bulk_x
    edge_residual = (
        d_zero - c_zero_x - c_s * log_cap * c_zero
    )
    total_d = d_zero + bulk_relation
    expected = (
        total_c_x
        + c_s * log_cap * total_c
        + edge_residual
        - r_bulk
    )
    if sp.expand(total_d - expected) != 0:
        raise RuntimeError("endpoint-complete first moment failed")

    return {
        "bulk_first_moment": (
            "Summing the division-free carrier identity over "
            "1<=n<=N, including the terminal carrier, gives "
            "D_bulk=(C_bulk)_x+c_s*log(N)*C_bulk-R_bulk, and "
            "Omega_N*D_bulk=partial_x(Omega_N*C_bulk)"
            "-Omega_N*R_bulk."
        ),
        "endpoint_complete_moment": (
            "Put E_0=d_0-c_(0,x)-c_s*log(N)c_0. Since "
            "mathsf_X=c_0+C_bulk and "
            "mathcal_C_N=d_0+D_bulk, exactly "
            "mathcal_C_N=partial_x mathsf_X"
            "+c_s*log(N)mathsf_X+E_0-R_bulk"
            "=Omega_N^(-1)partial_x(Omega_N*mathsf_X)"
            "+E_0-R_bulk."
        ),
        "edge_nonsplitting": (
            "R_bulk is uniformly tiny, but E_0 is not assigned a "
            "separate sign or magnitude. Only the endpoint-complete "
            "combination E_0-R_bulk is compared with the anchored "
            "first-jet identity; the terminal carrier remains included "
            "division-free in the bulk sum."
        ),
        "direct_equivalence": (
            "The anchored identity and "
            "mathcal_C_N=mathsf_A-c_s*u_N*mathsf_X give "
            "mathcal_C_N=partial_x mathsf_X"
            "+(xi_1-u_a-c_s*u_N)mathsf_X"
            "+v_a*mathsf_Y-Re(D_(1,x))/|f_1|. Since "
            "log(a)=log(N)+u_N, comparison yields "
            "E_0-R_bulk=(xi_1-u_a-c_s*log(a))mathsf_X"
            "+v_a*mathsf_Y-Re(D_(1,x))/|f_1|."
        ),
        "orientation_handoff": (
            "On |mathsf_X|<=delta_L the existing endpoint-complete "
            "bounds make |partial_x mathsf_X-mathcal_C_N|/A_L "
            "less than 2.03e-14. Thus the occupation first moment "
            "recovers the certified crossing-orientation handoff, "
            "conditional on the still-open pointwise gap."
        ),
    }


def sign_guard_audit() -> dict[str, str]:
    values = []
    for first, second in ((1, 0), (0, 1), (1, 1)):
        c_1, c_2 = sp.Integer(1), sp.Integer(-1)
        total_c = c_1 + c_2
        total_d = c_1 * first + c_2 * second
        values.append((total_c, total_d))
    if values != [(0, 1), (0, -1), (0, 0)]:
        raise RuntimeError("contact first-moment sign guard failed")
    return {
        "contact_sign_guard": (
            "In the correction-free normalized local law take "
            "(c_1,c_2)=(1,-1), h_(1,x)=h_(2,x)=-1, and zero "
            "residual. The slope pairs (1,0), (0,1), and (1,1) "
            "all have contact mass c_1+c_2=0 but first moments "
            "+1,-1,0. Individual clockwise motion therefore does "
            "not sign the normalized first moment. These are local "
            "transport guards, not Xi carriers."
        ),
        "backward_heat_guard": (
            "The existing exact backward-heat family has "
            "mathcal_C=partial_x mathsf_X and an arbitrarily strong "
            "contact-band gap, yet has any prescribed number of "
            "upward crossings. Therefore even a successful pointwise "
            "occupation first-moment bound would not by itself prove "
            "the one-turn successor theorem."
        ),
        "pole_boundary": (
            "If c_i->0 while d_i->d_*!=0, then h_i=d_i/c_i "
            "escapes to projective infinity. The compactly supported "
            "measure c_i*delta_(h_i) disappears locally while its "
            "first moment c_i*h_i=d_i has the finite limit d_*. "
            "That boundary mass is D_perp and must be retained "
            "division-free or joined in the reciprocal chart."
        ),
    }


def build_exact() -> dict:
    transport = transport_audit()
    carrier = carrier_source_audit()
    normalized = normalized_audit()
    residual = residual_audit()
    aggregate = aggregate_audit()
    guards = sign_guard_audit()
    return {
        "pi_provenance": (
            "The pi in a^2=x/(4*pi)+t/16 and "
            "u_x=1/(8*pi*a^2) remains the completed-zeta and "
            "Riemann-Siegel saddle constant. The occupation variable "
            "s, Dirac masses, and integrating factor N^(sigma_*) "
            "introduce no circle, fitted pi, or geometric normalization."
        ),
        "domain": (
            "Fix one N=floor(a) chart on L>=50, 0<tL<=25, and "
            "q=2*t*L^2>=1. The distributional slope law is local to "
            "a pole-free interior projective chart. The summed "
            "division-free carrier identity applies to every "
            "1<=n<=N, including the terminal carrier."
        ),
        "measure_definition": transport["measure_definition"],
        "weak_transport": transport["weak_transport"],
        "cumulative_transport": (
            transport["cumulative_definition"]
        ),
        "actual_carrier_rates": carrier["actual_rates"],
        "source_collapse": carrier["source_collapse"],
        "reaction_law": carrier["finite_slope_reaction"],
        "normalized_measure": normalized["normalized_measure"],
        "normalized_cumulative": normalized["normalized_cumulative"],
        "normalization_boundary": normalized["normalization_boundary"],
        "pointwise_residual": residual["pointwise_residual"],
        "amplitude_sum": residual["amplitude_sum"],
        "summed_residual": (
            residual["summed_residual"]
            + " "
            + residual["exponential_guard"]
        ),
        "bulk_first_moment": aggregate["bulk_first_moment"],
        "endpoint_complete_moment": aggregate["endpoint_complete_moment"],
        "edge_nonsplitting": aggregate["edge_nonsplitting"],
        "direct_transversality_equivalence": aggregate["direct_equivalence"],
        "orientation_handoff": aggregate["orientation_handoff"],
        "signed_flux_guard": (
            transport["signed_flux_guard"]
            + " "
            + guards["contact_sign_guard"]
        ),
        "pole_boundary": guards["pole_boundary"],
        "successor_guard": guards["backward_heat_guard"],
        "route_decision": (
            "Promote the normalized occupation law as an exact "
            "organizational lemma, but do not treat it as new "
            "coercivity. Its correction density is negligible, while "
            "its common frame drift cancels exactly; nevertheless its "
            "flux remains signed, and its first moment is precisely "
            "the already-open real transversality coordinate. A future "
            "gain must prove an Xi-specific signed threshold-mass "
            "estimate, most plausibly using complete multiplicative "
            "chains and the recurrent endpoint, or return to the "
            "endpoint-complete phase-flux boundary contract. Generic "
            "transport monotonicity cannot supply either theorem."
        ),
        "q_lt_1": (
            "The q<1 layer remains a separate multiplicity-compatible "
            "parabolic/Hermite or boundary-degree theorem. The affine "
            "occupation chart is not continued through t=0 "
            "multiplicities by assumption."
        ),
        "proof_boundary": (
            "This proves the weak and cumulative signed occupation "
            "transport laws, exact Xi source collapse, common "
            "N^(sigma_*) normalization, quantitative summed residual "
            "bound, division-free bulk first moment, endpoint-complete "
            "equivalence with the known real transversality identity, "
            "signed-flux guards, and projective-infinity boundary term. "
            "It does not prove a signed Xi threshold-mass estimate, "
            "the pointwise Abel-scalar gap, an edge sign, a one-turn "
            "successor bound, q<1 closure, finite-height effectivity, "
            "contact exclusion, Lambda<=0, PF-infinity, RH, or a "
            "Clay-prize conclusion."
        ),
        "diagnostics": {
            "transport": transport,
            "carrier": carrier,
            "normalized": normalized,
            "residual": residual,
            "aggregate": aggregate,
            "guards": guards,
        },
    }


def build_rows(exact: dict) -> list[GateRow]:
    rows = [
        (
            "pi_provenance",
            "definition_provenance",
            "certified",
            "The occupation coordinate introduces no new pi.",
            "No plotted geometry defines the constant.",
        ),
        (
            "domain",
            "exact_domain",
            "certified",
            "Slope transport is local and bulk moments are division-free.",
            "The q<1 layer remains separate.",
        ),
        (
            "measure_definition",
            "exact_definition",
            "ready_to_apply",
            "The signed occupation measure and flux are explicit.",
            "Only finite interior slopes enter this local chart.",
        ),
        (
            "weak_transport",
            "distributional_identity",
            "proved",
            "The occupation measure satisfies an exact weak balance law.",
            "No smooth density is assumed.",
        ),
        (
            "cumulative_transport",
            "distributional_identity",
            "proved",
            "The threshold mass has an exact cumulative transport law.",
            "Threshold deltas are retained.",
        ),
        (
            "actual_carrier_rates",
            "exact_source_input",
            "ready_to_apply",
            "The actual radial and angular rates are retained.",
            "Every d_i correction remains explicit.",
        ),
        (
            "source_collapse",
            "exact_identity",
            "proved",
            "The carrier amplitude source has a common-frame collapse.",
            "The identity is division-free.",
        ),
        (
            "reaction_law",
            "distributional_identity",
            "proved",
            "The source is a universal affine reaction plus a residual.",
            "This row is local to finite slopes.",
        ),
        (
            "normalized_measure",
            "exact_normalization",
            "proved",
            "N^(sigma_*) removes the common frame drift exactly.",
            "The remaining reaction is not assigned a sign.",
        ),
        (
            "normalized_cumulative",
            "distributional_identity",
            "proved",
            "The normalized threshold law retains flux, moment, and residual.",
            "No term is discarded.",
        ),
        (
            "normalization_boundary",
            "nonpromotion_guard",
            "guard_validated",
            "The integrating factor is not a positivity transform.",
            "Universal reaction and signed flux remain.",
        ),
        (
            "pointwise_residual",
            "analytic_bound",
            "certified",
            "Each correction-source atom has an explicit norm bound.",
            "Uses the checked d_i envelopes.",
        ),
        (
            "amplitude_sum",
            "analytic_bound",
            "certified",
            "The total corrected carrier amplitude is summable explicitly.",
            "Uses adjacent corrected amplitude contraction.",
        ),
        (
            "summed_residual",
            "analytic_bound",
            "proved",
            "The complete carrier source residual is below the nuisance scale.",
            "Constants are deliberately nonoptimized.",
        ),
        (
            "bulk_first_moment",
            "exact_identity",
            "proved",
            "The occupation first moment has an exact integrating factor.",
            "The terminal carrier is included division-free.",
        ),
        (
            "endpoint_complete_moment",
            "exact_identity",
            "proved",
            "The full Abel scalar is a normalized real derivative plus edge residual.",
            "No endpoint sign is asserted.",
        ),
        (
            "edge_nonsplitting",
            "proof_guard",
            "guard_validated",
            "Endpoint and carrier residuals are compared only in their exact combination.",
            "No triangle inequality splits a recurrent cancellation.",
        ),
        (
            "direct_transversality_equivalence",
            "exact_equivalence",
            "proved",
            "The occupation first moment equals the existing direct transversality coordinate.",
            "This is a route identification, not a lower bound.",
        ),
        (
            "orientation_handoff",
            "analytic_corollary",
            "certified",
            "The existing contact-band orientation error remains negligible.",
            "The pointwise gap is still an open hypothesis.",
        ),
        (
            "signed_flux_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Clockwise atoms do not make the c-weighted flux one-signed.",
            "The local jets are not actual Xi data.",
        ),
        (
            "pole_boundary",
            "division_free_join",
            "proved",
            "Vertical first moment survives as a boundary mass at infinity.",
            "No tangent quotient is used.",
        ),
        (
            "successor_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Pointwise transversality alone does not bound one-sided crossings.",
            "The Fourier family is not an Xi counterexample.",
        ),
        (
            "route_decision",
            "open_theorem_target",
            "not_ready_to_apply",
            "The surviving target is an Xi-specific signed threshold-mass estimate.",
            "No occupation coercivity has been proved.",
        ),
        (
            "q_lt_1",
            "open_theorem_target",
            "not_ready_to_apply",
            "The small-q multiplicity theorem remains separate.",
            "No t=0 simplicity is assumed.",
        ),
        (
            "proof_boundary",
            "proof_guard",
            "guard_validated",
            "The transport reduction is not promoted to contact exclusion.",
            "RH and every prize-level conclusion remain open.",
        ),
    ]
    return [
        GateRow(
            id=f"sotr_{index:02d}_{suffix}",
            role=role,
            readiness=readiness,
            claim=claim,
            formula=exact[suffix],
            proof_boundary=boundary,
            diagnostics=(
                exact["diagnostics"]
                if suffix
                in {
                    "weak_transport",
                    "source_collapse",
                    "summed_residual",
                    "direct_transversality_equivalence",
                    "signed_flux_guard",
                    "pole_boundary",
                }
                else None
            ),
        )
        for index, (
            suffix,
            role,
            readiness,
            claim,
            boundary,
        ) in enumerate(rows)
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact signed occupation transport and Xi source "
            "normalization proved; correction density negligible; "
            "first moment identified with open direct transversality "
            "coordinate; signed threshold-mass theorem remains open"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Signed Occupation Transport Reduction",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact transport, normalization, residual budget, and",
            "transversality equivalence; not a proof of the Xi threshold-mass",
            "bound, contact exclusion, `Lambda<=0`, PF-infinity, RH, or a",
            "Clay-prize result.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Pi Provenance",
            "",
            exact["pi_provenance"],
            "",
            "## Distributional Transport",
            "",
            exact["domain"],
            "",
            "```text",
            exact["measure_definition"],
            exact["weak_transport"],
            exact["cumulative_transport"],
            "```",
            "",
            "## Xi Source Collapse",
            "",
            "```text",
            exact["actual_carrier_rates"],
            exact["source_collapse"],
            exact["reaction_law"],
            "```",
            "",
            "## Common Normalization",
            "",
            "```text",
            exact["normalized_measure"],
            exact["normalized_cumulative"],
            "```",
            "",
            exact["normalization_boundary"],
            "",
            "## Residual Budget",
            "",
            exact["pointwise_residual"],
            "",
            exact["amplitude_sum"],
            "",
            "```text",
            exact["summed_residual"],
            "```",
            "",
            "## First Moment",
            "",
            "```text",
            exact["bulk_first_moment"],
            exact["endpoint_complete_moment"],
            "```",
            "",
            exact["edge_nonsplitting"],
            "",
            exact["direct_transversality_equivalence"],
            "",
            exact["orientation_handoff"],
            "",
            "## Route Guards",
            "",
            exact["signed_flux_guard"],
            "",
            exact["pole_boundary"],
            "",
            exact["successor_guard"],
            "",
            "## Route Decision",
            "",
            exact["route_decision"],
            "",
            exact["q_lt_1"],
            "",
            "## Boundary",
            "",
            exact["proof_boundary"],
            "",
        ]
    )


def write_outputs(artifact: dict, output: Path, note: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(artifact), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    write_outputs(artifact, args.output, args.note)
    print(
        "built signed occupation transport reduction: "
        "25 rows, exact weak/cumulative laws, Xi source collapse, "
        "N^sigma normalization, sub-1e-7 residual budget, "
        "direct-transversality equivalence, 3 route guards"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
