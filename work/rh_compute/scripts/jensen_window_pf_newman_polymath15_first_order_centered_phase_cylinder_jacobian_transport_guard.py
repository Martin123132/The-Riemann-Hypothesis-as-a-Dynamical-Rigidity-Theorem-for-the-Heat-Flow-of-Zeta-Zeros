#!/usr/bin/env python3
"""Build the phase-cylinder Jacobian and winding-transport guard."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "phase_cylinder_jacobian_transport_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "joined_odd_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "prime_power_guard": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "prime_power_heat_block_composition_guard.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
}

L_MIN = 50


@dataclass(frozen=True)
class TransportRow:
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
        "joined_odd_prefix": (
            "P_N(z_*,omega_*,1)=H_0",
            "H_0,H_1,H_2,D_0,D_1",
            "B_K*(R_0-1)^K",
        ),
        "prime_power_guard": (
            "p^(-12/25)",
            '"joined_winding": 3',
        ),
        "normalized_prefix": (
            "Gamma_ell=mathsf_X+i*mathsf_A/ell",
            "partial_x arg(Gamma_ell)",
        ),
        "adjacent_recurrence": (
            "Q_N=f_(N+1)+kappa_N*J_a",
            "Delta A_a=Re[-s_*'*ell*f_(N+1)",
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


def pair_kernel_audit() -> dict[str, str]:
    k_1, k_2, l_1, l_2, a_1, a_2 = sp.symbols(
        "k_1 k_2 l_1 l_2 a_1 a_2", real=True
    )
    phi_1, phi_2 = sp.symbols("phi_1 phi_2", real=True)
    term_1 = a_1 * sp.exp(sp.I * phi_1)
    term_2 = a_2 * sp.exp(sp.I * phi_2)
    p_theta = sp.I * (k_1 * term_1 + k_2 * term_2)
    p_xi = sp.I * (l_1 * term_1 + l_2 * term_2)
    jacobian = sp.expand_complex(sp.conjugate(p_theta) * p_xi).as_real_imag()[1]
    expected = (
        a_1
        * a_2
        * (k_1 * l_2 - k_2 * l_1)
        * sp.sin(phi_2 - phi_1)
    )
    if sp.trigsimp(jacobian - expected) != 0:
        raise RuntimeError("two-term phase-cylinder kernel failed")
    return {
        "partials": (
            "P_theta=i*sum_alpha k_alpha w_alpha exp(i phi_alpha), "
            "P_xi=i*sum_alpha l_alpha w_alpha exp(i phi_alpha)"
        ),
        "jacobian": "J_(theta,xi)=Im(conj(P_theta)*P_xi)",
        "pair_kernel": (
            "sum_(alpha<beta) w_alpha*w_beta"
            "*(k_alpha*l_beta-k_beta*l_alpha)"
            "*sin(phi_beta-phi_alpha)"
        ),
        "phase": (
            "phi_alpha=k_alpha*theta+l_alpha*xi, "
            "l_alpha=log(m_alpha)"
        ),
    }


def xi_support_sign_guard() -> dict[str, str]:
    h = sp.log(2)
    l_3 = sp.log(3)
    determinant = sp.simplify(1 * l_3 - 0 * 0)
    linked_phase = sp.simplify((0 - 1) * h + (l_3 - 0))
    if determinant != sp.log(3):
        raise RuntimeError("n=2,3 determinant failed")
    if linked_phase != sp.log(sp.Rational(3, 2)):
        raise RuntimeError("n=2,3 linked phase failed")
    return {
        "support": "n=2=(k,m)=(1,1), n=3=(k,m)=(0,3)",
        "positive_amplitudes": "w_2>0, w_3>0",
        "linked_kernel": (
            "J_23=w_2*w_3*log(3)"
            "*sin(xi*log(3/2))"
        ),
        "positive_point": "xi*log(3/2)=pi/2 modulo 2*pi",
        "negative_point": "xi*log(3/2)=3*pi/2 modulo 2*pi",
        "interpretation": (
            "Even genuine dyadic/odd support with positive amplitudes "
            "has an oscillating pair current. The full Jacobian may "
            "still have structure, but no termwise-positive proof is "
            "available."
        ),
    }


def linked_monotonicity_guard() -> dict[str, str]:
    r = sp.sqrt(2) / 2
    minimum_current = sp.simplify(r**2 - r)
    if not minimum_current < 0:
        raise RuntimeError("linked angular-current guard failed")
    return {
        "prefix": "P(phi)=1+r*exp(i*phi), r=2^(-1/2)",
        "nonvanishing": "|P(phi)|>=1-r>0",
        "first_moment": "H_1=log(2)*r*exp(i*phi)",
        "linked_argument_numerator": (
            "Re(conj(P)*H_1)="
            "log(2)*(r^2+r*cos(phi))"
        ),
        "negative_value": (
            "at phi=pi it equals "
            "log(2)*(1/2-1/sqrt(2))<0"
        ),
        "interpretation": (
            "A zero-free complete dyadic block need not have monotone "
            "argument along the linked physical phase."
        ),
    }


def endpoint_augmentation_guard() -> dict[str, str]:
    z = sp.symbols("z")
    r = sp.Rational(2, 3)
    prefix = 1 + r * z
    root = sp.solve(prefix, z)[0]
    endpoint = -(1 + r)
    if root != -sp.Rational(3, 2):
        raise RuntimeError("endpoint guard prefix root failed")
    if sp.simplify(prefix.subs(z, 1) + endpoint) != 0:
        raise RuntimeError("endpoint augmentation zero failed")
    return {
        "prefix": "P(z)=1+(2/3)z",
        "prefix_zero": "z=-3/2, outside the closed unit disk",
        "endpoint": "e=-5/3",
        "augmented_zero": "P(1)+e=0",
        "interpretation": (
            "Unit-disk zero-freeness of the finite prefix does not "
            "survive an arbitrary endpoint addition. The actual "
            "Riemann-Siegel endpoint needs its own transport current."
        ),
    }


def orientation_audit() -> dict[str, str]:
    # F(theta,xi)=exp(i theta)-r(xi), with r increasing through one.
    # At the crossing theta=0, J=Im(conj(i)*(-r'))=r'>0, while the
    # theta winding drops from one to zero.
    theta = sp.symbols("theta", real=True)
    r_prime = sp.symbols("r_prime", positive=True, real=True)
    p_theta = sp.I * sp.exp(sp.I * theta)
    p_xi_at_crossing = -r_prime
    jacobian = sp.expand_complex(
        sp.conjugate(p_theta.subs(theta, 0)) * p_xi_at_crossing
    ).as_real_imag()[1]
    if jacobian != r_prime:
        raise RuntimeError("winding-transport orientation failed")
    return {
        "orientation_check": (
            "For P=exp(i theta)-r(xi) with r increasing through one, "
            "J>0 at the crossing while wind_theta drops from one to "
            "zero."
        ),
        "transport_sign": (
            "wind(theta,xi_1)-wind(theta,xi_0)="
            "-sum_(P=0) sign J_(theta,xi)"
        ),
    }


def symbolic_audit() -> dict[str, dict | str]:
    pair_kernel = pair_kernel_audit()
    z, xi, h = sp.symbols("z xi h")
    k, log_m = sp.symbols("k log_m", real=True)
    monomial = z**k * sp.exp(sp.I * xi * log_m)
    p_theta = sp.I * z * sp.diff(monomial, z)
    p_xi = sp.diff(monomial, xi)
    linked = sp.simplify(h * p_theta + p_xi)
    expected = sp.I * (h * k + log_m) * monomial
    if sp.simplify(linked - expected) != 0:
        raise RuntimeError("linked physical derivative failed")

    p, h_1 = sp.symbols("P H_1")
    # The displayed identity follows from Im(conj(P)*iH_1)=Re(conj(P)H_1).
    real_part, imag_part, h_real, h_imag = sp.symbols(
        "p_r p_i h_r h_i", real=True
    )
    p_complex = real_part + sp.I * imag_part
    h_complex = h_real + sp.I * h_imag
    lhs = sp.expand_complex(
        sp.conjugate(p_complex) * sp.I * h_complex
    ).as_real_imag()[1]
    rhs = sp.expand_complex(
        sp.conjugate(p_complex) * h_complex
    ).as_real_imag()[0]
    if sp.expand(lhs - rhs) != 0:
        raise RuntimeError("linked argument-current identity failed")

    mu, d_n, d_1 = sp.symbols("mu d_n d_1")
    correction = (1 + mu * d_n) / (1 + mu * d_1)
    correction_mu = sp.factor(sp.diff(correction, mu))
    expected_mu = (d_n - d_1) / (1 + mu * d_1) ** 2
    if sp.simplify(correction_mu - expected_mu) != 0:
        raise RuntimeError("correction-cylinder derivative failed")

    return {
        "pair_kernel": pair_kernel,
        "linked_derivative": (
            "On theta=h*xi, dP/dxi=h*P_theta+P_xi="
            "i*mathcal D_hP=i*H_1."
        ),
        "linked_argument_current": (
            "Where H_0=P is nonzero, "
            "partial_xi arg P_link="
            "Re(conj(H_0)H_1)/|H_0|^2."
        ),
        "correction_derivative": (
            "P_mu=sum_(k,m) z^k a_(k,m)exp(i xi log m)"
            "*(d_(2^k m)-d_1)/(1+mu d_1)^2."
        ),
        "orientation": orientation_audit(),
        "xi_support_guard": xi_support_sign_guard(),
        "linked_monotonicity_guard": linked_monotonicity_guard(),
        "endpoint_guard": endpoint_augmentation_guard(),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a), with the phase "
            "transport performed at fixed physical x,t,N before the "
            "physical x derivative is restored"
        ),
        "pi_provenance": (
            "The pi in a^2=x/(4*pi)+t/16 comes from the completed "
            "zeta normalization and Riemann-Siegel saddle; 4*pi(2N+1) "
            "is the cutoff-cell difference and 2*pi is the period of "
            "exp(i theta). The phase cylinder and its Jacobians do not "
            "define pi."
        ),
        "phase_lift": (
            "Retain P_N(z,xi,mu)=sum_(k=0)^K z^k"
            "*sum_(m<=M_k,m odd)a_(k,m)exp(i xi log m)"
            "*(1+mu d_(2^k m))/(1+mu d_1), with z=exp(i theta). "
            "The synchronized base is (xi,mu)=(0,0); the actual finite "
            "prefix is (theta,xi,mu)=(omega_*log 2,omega_*,1)."
        ),
        "phase_partials": symbolic["pair_kernel"]["partials"],
        "phase_jacobian": (
            symbolic["pair_kernel"]["jacobian"]
            + ". It is the oriented real Jacobian "
            "det partial_(theta,xi)(Re P_N,Im P_N)."
        ),
        "phase_pair_kernel": symbolic["pair_kernel"]["pair_kernel"]
        + ", with "
        + symbolic["pair_kernel"]["phase"]
        + ". For complex correction coefficients, replace "
        "w_alpha*w_beta*sin(phi_beta-phi_alpha) by "
        "Im(conj(c_alpha)c_beta exp(i(phi_beta-phi_alpha))).",
        "linked_derivative": symbolic["linked_derivative"],
        "linked_argument_current": symbolic["linked_argument_current"],
        "physical_x_current": (
            "At the actual point, H_(0,x)=-s_*'H_1+D_0. Thus the "
            "linked phase current iH_1 is only the phase component; "
            "the sigma_* amplitude drift and full d_n current D_0 "
            "must be restored before an x-transversality claim."
        ),
        "winding_transport": (
            "Assume P_N is nonzero on the two boundary circles "
            "xi=xi_0,xi_1 and has only isolated regular zeros in "
            "S^1_theta times (xi_0,xi_1). With "
            "J_(theta,xi)=Im(conj(P_theta)P_xi), the degree theorem "
            "gives wind_theta P(.,xi_1)-wind_theta P(.,xi_0)="
            "-sum_(P=0)sign J_(theta,xi). "
            + symbolic["orientation"]["orientation_check"]
        ),
        "xi_support_sign_guard": symbolic["xi_support_guard"],
        "linked_monotonicity_guard": symbolic[
            "linked_monotonicity_guard"
        ],
        "correction_cylinder": (
            symbolic["correction_derivative"]
            + " Define J_(theta,mu)=Im(conj(P_theta)P_mu). "
            "If the boundary circles are nonzero and all interior "
            "zeros are regular, the same orientation gives "
            "wind_theta P(.,mu_1)-wind_theta P(.,mu_0)="
            "-sum_(P=0)sign J_(theta,mu)."
        ),
        "endpoint_cylinder": (
            "At fixed physical x, after the finite prefix is restored, "
            "define F_N(z,tau)=P_N(z,omega_*,1)+tau*r_0 for "
            "0<=tau<=1. Then F_theta=P_theta, F_tau=r_0, and "
            "J_(theta,tau)=Im(conj(P_theta)r_0). The same cylinder "
            "degree formula transports the finite-prefix winding to "
            "the endpoint-augmented value Z_0. This endpoint cylinder "
            "is indispensable."
        ),
        "endpoint_augmentation_guard": symbolic["endpoint_guard"],
        "three_cylinder_programme": (
            "A sufficient strong route starts from the proved "
            "P_N(z,0,0), transports odd phase xi:0->omega_*, then "
            "correction mu:0->1, then endpoint tau:0->1. Each stage "
            "has an exact local Jacobian and signed winding-change "
            "formula. The final physical theorem may use a weaker "
            "linked-point argument, but no stage may be skipped merely "
            "because the preceding circle was zero-free."
        ),
        "endpoint_first_jet": (
            "The final physical derivatives remain "
            "Z_(0,x)=r_(0,x)-s_*'H_1+D_0 and "
            "Z_(A,x)=r_(A,x)-s_*''P_1-s_*'P_(1,x), where "
            "P_(1,x)=-s_*'(H_2-log(a)H_1)+D_1-log(a)D_0"
            "-H_0/(8*pi*a^2). After eta projection, these enter "
            "Gamma_ell=mathsf_X+i mathsf_A/ell and its exact argument "
            "flux. The adjacent recurrence is retained at cutoff "
            "boundaries."
        ),
        "route_decision": (
            "Reject a global or termwise sign theorem for "
            "J_(theta,xi): its genuine n=2,3 pair kernel oscillates. "
            "Reject monotone linked argument: even the zero-free "
            "1+2^(-1/2)exp(i phi) block has negative angular current "
            "on part of the phase cycle. Retain the Jacobians as local "
            "degree currents. The live theorem is a signed crossing "
            "budget, no-zero cylinder theorem, or Xi-specific "
            "cancellation theorem across the odd-phase, correction, "
            "and endpoint cylinders."
        ),
        "live_q_ge_1_target": (
            "On q=2tL^2>=1, prove enough boundary nonvanishing and "
            "signed-degree control for the three cylinders to transport "
            "the synchronized base to the endpoint-augmented physical "
            "first jet. A sufficient result is zero interior degree "
            "together with nonvanishing final boundary; a more local "
            "result may control only the linked point and complete "
            "successor path. It must still imply "
            "|mathsf_X|<=delta_L => "
            "|mathcal C_N|>A_L+epsilon_term and 0<=kappa_j<1."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 layer still needs the separate "
            "multiplicity-compatible parabolic/Hermite first-jet chart. "
            "The cylinder identities remain exact, but no endpoint "
            "simplicity or global Jacobian sign is imposed there."
        ),
        "proof_boundary": (
            "The phase partials, pair-kernel Jacobian, linked derivative "
            "and argument current, cylinder winding-transport theorem, "
            "correction and endpoint cylinders, and three route guards "
            "are exact. No sign or nonzero lower bound for the full Xi "
            "Jacobian, no three-cylinder crossing budget, no joined Xi "
            "lower bound, no strict successor flux upper bound, no q<1 "
            "closure, no contact exclusion, no Lambda<=0 or PF-infinity "
            "theorem, no RH proof, and no Clay-prize conclusion is "
            "asserted."
        ),
        "diagnostics": {
            "pair_kernel_terms": "one unordered term per distinct lifted pair",
            "xi_support_pair": "n=2 and n=3",
            "linked_zero_free_ratio": format(2 ** -0.5, ".17e"),
            "winding_transport_orientation": symbolic["orientation"][
                "transport_sign"
            ],
        },
    }


def build_rows(exact: dict) -> list[TransportRow]:
    return [
        TransportRow(
            "pcjtg_00_pi_provenance",
            "definition_provenance",
            "certified",
            "Every pi in the cylinder reduction has an explicit completed-zeta, saddle, or phase-period source.",
            exact["pi_provenance"],
            "The phase cylinder does not define pi.",
        ),
        TransportRow(
            "pcjtg_01_phase_lift",
            "exact_recall",
            "ready_to_apply",
            "The joined finite prefix has a canonical dyadic/odd/correction phase lift.",
            exact["phase_lift"],
            "The endpoint is added in a separate exact cylinder.",
        ),
        TransportRow(
            "pcjtg_02_phase_partials",
            "exact_differential_identity",
            "ready_to_apply",
            "Both phase-cylinder tangent vectors are explicit finite sums.",
            exact["phase_partials"],
            "No sign is assigned termwise.",
        ),
        TransportRow(
            "pcjtg_03_pair_jacobian",
            "exact_pair_kernel",
            "ready_to_apply",
            "The oriented phase-cylinder Jacobian has an exact unordered-pair kernel.",
            exact["phase_jacobian"] + " " + exact["phase_pair_kernel"],
            "The sine kernel is retained.",
        ),
        TransportRow(
            "pcjtg_04_linked_current",
            "exact_linked_direction",
            "ready_to_apply",
            "The physical dyadic/odd phase linkage is exactly the total logarithmic moment direction.",
            exact["linked_derivative"]
            + " "
            + exact["linked_argument_current"]
            + " "
            + exact["physical_x_current"],
            "Amplitude and d_n currents must be restored for x motion.",
        ),
        TransportRow(
            "pcjtg_05_winding_transport",
            "exact_topological_theorem",
            "ready_to_apply",
            "Regular cylinder zeros transport the unit-circle winding by their signed Jacobian degree.",
            exact["winding_transport"],
            "Boundary nonvanishing and regularity are hypotheses, not conclusions.",
        ),
        TransportRow(
            "pcjtg_06_xi_pair_sign_guard",
            "countermodel",
            "guard_validated",
            "A genuine dyadic/odd pair kernel changes sign on the linked phase.",
            json.dumps(exact["xi_support_sign_guard"], sort_keys=True),
            "This rejects termwise positivity, not every possible full-sum identity.",
            exact["xi_support_sign_guard"],
        ),
        TransportRow(
            "pcjtg_07_linked_monotonicity_guard",
            "countermodel",
            "guard_validated",
            "A zero-free complete dyadic block need not have monotone linked argument.",
            json.dumps(exact["linked_monotonicity_guard"], sort_keys=True),
            "Zero-freeness and angular-current sign are different statements.",
            exact["linked_monotonicity_guard"],
        ),
        TransportRow(
            "pcjtg_08_correction_cylinder",
            "exact_transport_identity",
            "ready_to_apply",
            "The d_n restoration has its own exact cylinder Jacobian and degree formula.",
            exact["correction_cylinder"],
            "A quantitative correction margin remains open.",
        ),
        TransportRow(
            "pcjtg_09_endpoint_cylinder",
            "exact_transport_identity",
            "ready_to_apply",
            "Endpoint augmentation has an exact local degree current.",
            exact["endpoint_cylinder"],
            "The finite-prefix theorem cannot bypass this cylinder.",
        ),
        TransportRow(
            "pcjtg_10_endpoint_guard",
            "countermodel",
            "guard_validated",
            "Adding an endpoint can destroy unit-disk zero-freeness.",
            json.dumps(exact["endpoint_augmentation_guard"], sort_keys=True),
            "The countermodel endpoint is generic, not the Xi endpoint.",
            exact["endpoint_augmentation_guard"],
        ),
        TransportRow(
            "pcjtg_11_three_cylinders",
            "route_reduction",
            "ready_to_apply",
            "Odd phase, correction, and endpoint restoration form three exact winding-transport stages.",
            exact["three_cylinder_programme"],
            "A full-circle theorem is sufficient and may be stronger than the linked-point target.",
        ),
        TransportRow(
            "pcjtg_12_endpoint_first_jet",
            "exact_endpoint_handoff",
            "ready_to_apply",
            "The transported value must feed the complete five-current endpoint first jet.",
            exact["endpoint_first_jet"],
            "H_2, D_0, D_1, and adjacent recurrence remain present.",
        ),
        TransportRow(
            "pcjtg_13_route_decision",
            "route_decision",
            "guard_validated",
            "The Jacobian survives as a degree current, not a globally signed density.",
            exact["route_decision"],
            "No signed Xi crossing budget is yet proved.",
        ),
        TransportRow(
            "pcjtg_14_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 obligation is a three-cylinder signed-degree or no-zero theorem with endpoint first jet.",
            exact["live_q_ge_1_target"],
            "The Xi lower bound and strict successor bound remain open.",
        ),
        TransportRow(
            "pcjtg_15_q_lt_1_nonpromotion",
            "open_theorem_target",
            "not_ready_to_apply",
            "The small-q chart and every global conclusion remain open.",
            exact["q_lt_1_target"] + " " + exact["proof_boundary"],
            "No Lambda<=0, PF-infinity, RH, or Clay-prize conclusion is asserted.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-27",
        "status": (
            "exact phase-cylinder pair Jacobian, linked-direction "
            "current, winding-transport theorem, correction and endpoint "
            "cylinders, and three route guards; the Xi signed-degree "
            "budget and successor theorem remain open"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "constants": {"L_min": L_MIN},
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Phase-Cylinder Jacobian Transport Guard",
            "",
            "Date: 2026-07-27",
            "",
            "Status: exact Jacobian and winding-transport reduction with",
            "route guards; not a proof of `Lambda<=0`, RH, PF-infinity,",
            "or a Clay-prize result.",
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
            "## Phase Cylinder",
            "",
            "```text",
            exact["phase_lift"],
            exact["phase_partials"],
            "```",
            "",
            "## Pair-Kernel Jacobian",
            "",
            "```text",
            exact["phase_jacobian"],
            exact["phase_pair_kernel"],
            "```",
            "",
            "This Jacobian is the local degree density for unit-circle",
            "zeros crossing as the odd Mellin phase changes.",
            "",
            "## Physical Linked Direction",
            "",
            "```text",
            exact["linked_derivative"],
            exact["linked_argument_current"],
            exact["physical_x_current"],
            "```",
            "",
            "The phase-only linked direction is exact, but physical x",
            "motion also changes amplitudes and correction factors.",
            "",
            "## Winding Transport",
            "",
            "```text",
            exact["winding_transport"],
            "```",
            "",
            "This is the surviving role of the Jacobian. It counts signed",
            "boundary crossings; it is not globally positive.",
            "",
            "## Sign Guards",
            "",
            "### Genuine dyadic/odd pair",
            "",
            "```json",
            json.dumps(exact["xi_support_sign_guard"], indent=2, sort_keys=True),
            "```",
            "",
            "### Zero-free but nonmonotone linked block",
            "",
            "```json",
            json.dumps(
                exact["linked_monotonicity_guard"],
                indent=2,
                sort_keys=True,
            ),
            "```",
            "",
            "## Correction Cylinder",
            "",
            "```text",
            exact["correction_cylinder"],
            "```",
            "",
            "## Endpoint Cylinder",
            "",
            "```text",
            exact["endpoint_cylinder"],
            "```",
            "",
            "The need for this final transport is exact:",
            "",
            "```json",
            json.dumps(
                exact["endpoint_augmentation_guard"],
                indent=2,
                sort_keys=True,
            ),
            "```",
            "",
            "## Three-Cylinder Programme",
            "",
            "```text",
            exact["three_cylinder_programme"],
            exact["endpoint_first_jet"],
            "```",
            "",
            "## Route Decision",
            "",
            "```text",
            exact["route_decision"],
            "```",
            "",
            "## Live Theorem",
            "",
            "```text",
            exact["live_q_ge_1_target"],
            "```",
            "",
            "The separate small-`q` obligation is",
            "",
            "```text",
            exact["q_lt_1_target"],
            "```",
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
        "built Newman phase-cylinder Jacobian transport guard: "
        "16 rows, 1 exact pair-kernel Jacobian, "
        "1 cylinder winding-transport theorem, "
        "3 exact route guards, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
