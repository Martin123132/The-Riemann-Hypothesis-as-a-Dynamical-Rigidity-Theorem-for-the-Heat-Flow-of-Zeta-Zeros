#!/usr/bin/env python3
"""Build the endpoint-complete contact signed-transport reduction."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "contact_signed_transport_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "carrier_abel": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_first_order_signed_contact_reduction.json"
    ),
    "joined_first_jet": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "normalized_flux": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
    "legacy_symmetry": (
        REPO_ROOT
        / "work/rh_compute/results/legacy_prime_curvature_symmetry_audit.json"
    ),
}


@dataclass(frozen=True)
class ReductionRow:
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
        "carrier_abel": (
            "mathcal_C_N",
            "W_A=s_*'*u_N*W_0",
        ),
        "signed_contact": (
            "M*(h_+-h_-)",
            "rate-free near-crossing identity",
        ),
        "joined_first_jet": (
            "H_0,H_1,H_2,D_0,D_1",
            "Delta Z_A",
        ),
        "normalized_flux": (
            "B_(N,x)",
            "mathsf_A_x",
        ),
        "legacy_symmetry": (
            "transpose",
            "indefinite",
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


def centered_component_audit() -> dict[str, str]:
    c, b, u_n, u_terminal = sp.symbols(
        "c b u_n u_N", real=True
    )
    e_r, g_r, x_n, y_n = sp.symbols(
        "e_r g_r x_n y_n", real=True
    )
    alpha = c * u_terminal
    slope = g_r + u_n * (c * x_n - b * y_n)
    value = e_r + x_n
    centered = (
        g_r
        - alpha * e_r
        + c * (u_n - u_terminal) * x_n
        - b * u_n * y_n
    )
    if sp.simplify(centered - (slope - alpha * value)) != 0:
        raise RuntimeError("terminal-centered component identity failed")
    return {
        "centered_endpoint": (
            "c_0=Re(e), d_0=Re(g)-alpha*Re(e)"
        ),
        "centered_carrier": (
            "c_n=Re(z_n), "
            "d_n=c*(u_n-u_N)*Re(z_n)-b*u_n*Im(z_n)"
        ),
        "sum": (
            "sum_(j=0)^N c_j=mathsf_X and "
            "sum_(j=0)^N d_j=mathcal_C_N, alpha=c*u_N."
        ),
    }


def mass_transport_audit() -> dict[str, str]:
    m_plus, m_minus, h_plus, h_minus, d_perp = sp.symbols(
        "M_plus M_minus h_plus h_minus D_perp", real=True
    )
    x_value = m_plus - m_minus
    mass = (m_plus + m_minus) / 2
    direct = m_plus * h_plus - m_minus * h_minus + d_perp
    transported = (
        mass * (h_plus - h_minus)
        + x_value * (h_plus + h_minus) / 2
        + d_perp
    )
    if sp.simplify(direct - transported) != 0:
        raise RuntimeError("signed mass-transport identity failed")
    return {
        "partition": (
            "P={j:c_j>0}, N_-={j:c_j<0}, Z={j:c_j=0}; "
            "M_+=sum_P c_j, M_-=sum_N_(-c_j), "
            "h_+=sum_P d_j/M_+, h_-=sum_N_(-d_j)/M_-, "
            "D_perp=sum_Z d_j, and M=(M_++M_-)/2."
        ),
        "band_identity": (
            "mathcal_C_N=M*(h_+-h_-)"
            "+(mathsf_X/2)*(h_++h_-)+D_perp."
        ),
        "contact_identity": (
            "On mathsf_X=0, M_+=M_-=M and "
            "mathcal_C_N=M*(h_+-h_-)+D_perp."
        ),
    }


def ordered_transport_audit() -> dict[str, str]:
    h = sp.symbols("h1:5", real=True)
    weights = sp.symbols("w1:5", real=True)
    total = sum(weights)
    prefixes = [sum(weights[:index]) for index in range(1, 4)]
    rhs = h[-1] * total - sum(
        (h[index + 1] - h[index]) * prefixes[index]
        for index in range(3)
    )
    lhs = sum(h_j * w_j for h_j, w_j in zip(h, weights))
    if sp.simplify(lhs - rhs) != 0:
        raise RuntimeError("slope-ordered Abel identity failed")
    return {
        "ordering": (
            "Order the nonzero-real components so "
            "h_(1)<=...<=h_(m), where h_(j)=d_(j)/c_(j), "
            "and put P_k=sum_(j<=k)c_(j)."
        ),
        "band_identity": (
            "mathcal_C_N=D_perp+h_(m)*mathsf_X"
            "-sum_(k=1)^(m-1)(h_(k+1)-h_k)*P_k."
        ),
        "contact_identity": (
            "On mathsf_X=0, mathcal_C_N=D_perp"
            "-sum_(k=1)^(m-1)(h_(k+1)-h_k)*P_k."
        ),
    }


def endpoint_audit() -> dict[str, str]:
    kappa, h_a, j_re, j_im, time_0, c, u_terminal = sp.symbols(
        "kappa H_a J_re J_im T_0 c u_N", real=True
    )
    j_a = j_re + sp.I * j_im
    e = kappa * h_a * (time_0 + sp.I)
    g = kappa * j_a * (time_0 + sp.I)
    alpha = c * u_terminal
    c_0 = sp.re(e)
    d_0 = sp.re(g) - alpha * sp.re(e)
    expected_c = kappa * h_a * time_0
    expected_d = kappa * (
        time_0 * j_re - j_im - alpha * time_0 * h_a
    )
    effective_slope = (
        j_re / h_a - j_im / (time_0 * h_a) - alpha
    )
    if sp.simplify(c_0 - expected_c) != 0:
        raise RuntimeError("endpoint real-mass identity failed")
    if sp.simplify(d_0 - expected_d) != 0:
        raise RuntimeError("endpoint centered-slope identity failed")
    if sp.simplify(d_0 - c_0 * effective_slope) != 0:
        raise RuntimeError("endpoint effective-slope identity failed")
    return {
        "shape": (
            "e=kappa*H_a*(T_0+i), g=kappa*J_a*(T_0+i), "
            "with real kappa,H_a and generally complex "
            "J_a=H_(a,x)+mu_a*H_a."
        ),
        "components": (
            "c_0=kappa*T_0*H_a and "
            "d_0=kappa*[T_0*Re(J_a)-Im(J_a)"
            "-c*u_N*T_0*H_a]."
        ),
        "label": (
            "If H_a!=0 and c_0!=0, the endpoint effective slope is "
            "h_0=d_0/c_0=Re(J_a)/H_a"
            "-Im(J_a)/(T_0*H_a)-c*u_N."
        ),
        "zero_branch": (
            "If H_a=0, the endpoint value vanishes but "
            "d_0=kappa*[T_0*Re(J_a)-Im(J_a)] remains in D_perp; "
            "no division by H_a is allowed."
        ),
    }


def endpoint_null_guard_audit() -> dict[str, str]:
    delta = sp.Rational(1, 100)
    root = sp.sqrt(65)
    omega = (-8 + sp.I) / root
    s_prime = sp.Rational(1, 16) - sp.I / 2
    endpoint = delta * (8 + sp.I)
    endpoint_slope = endpoint / 8
    carriers = [
        omega,
        -sp.Rational(3, 5) * omega,
        -sp.Rational(2, 5) * omega,
        -endpoint,
    ]
    distances = [4, 3, 2, 1]
    value = sp.simplify(endpoint + sum(carriers))
    moment = sp.simplify(
        sum(u * z for u, z in zip(distances, carriers))
    )
    slope = sp.simplify(endpoint_slope + s_prime * moment)
    expected = sp.I * (
        sp.Rational(7, 80) * root + sp.Rational(13, 320)
    )
    amplitudes = [sp.simplify(sp.Abs(z)) for z in carriers]
    if value != 0:
        raise RuntimeError("endpoint-augmented zero-fibre value failed")
    if sp.simplify(slope - expected) != 0:
        raise RuntimeError("endpoint-augmented zero-fibre slope failed")
    if sp.re(slope) != 0 or sp.im(slope) <= 0:
        raise RuntimeError("endpoint-augmented contact null guard failed")
    numeric_amplitudes = [float(sp.N(value, 30)) for value in amplitudes]
    if not all(
        numeric_amplitudes[index] > numeric_amplitudes[index + 1]
        for index in range(3)
    ):
        raise RuntimeError("endpoint-augmented amplitude ordering failed")
    return {
        "s_prime": "1/16-i/2",
        "endpoint": "e=(8+i)/100, g=e/8",
        "carriers": (
            "z_1=omega, z_2=-(3/5)omega, z_3=-(2/5)omega, "
            "z_4=-e, omega=(-8+i)/sqrt(65)"
        ),
        "distances": "u_1,u_2,u_3,u_4=4,3,2,1",
        "amplitudes": (
            "1,3/5,2/5,sqrt(65)/100 are strictly decreasing."
        ),
        "value": "W_0=0",
        "slope": (
            "W_A=i*(7*sqrt(65)/80+13/320)!=0, "
            "so mathcal_C_N=Re(W_A)=0."
        ),
        "boundary": (
            "The endpoint has the exact factored (T_0+i) shape with a "
            "real synthetic J_a/H_a, but "
            "the chosen H_a,J_a and carriers are synthetic and do not "
            "satisfy the full Xi recurrence or d_n phase chain."
        ),
    }


def cutoff_audit() -> dict[str, str]:
    a_old, a_new, x_value, d_x, d_a, a_value = sp.symbols(
        "alpha_N alpha_next X Delta_X Delta_A A", real=True
    )
    old = a_value - a_old * x_value
    new = a_value + d_a - a_new * (x_value + d_x)
    jump = d_a - a_new * d_x + (a_old - a_new) * x_value
    if sp.simplify((new - old) - jump) != 0:
        raise RuntimeError("adjacent centered-scalar jump failed")
    return {
        "jump": (
            "For alpha_N=c*u_N and "
            "mathcal_C_N=mathsf_A_N-alpha_N*mathsf_X_N, "
            "mathcal_C_(N+1)-mathcal_C_N="
            "Delta mathsf_A-alpha_(N+1)*Delta mathsf_X"
            "+(alpha_N-alpha_(N+1))*mathsf_X_N."
        ),
        "contact_jump": (
            "On a contact of the N chart this becomes "
            "Delta mathcal_C=Delta mathsf_A"
            "-alpha_(N+1)*Delta mathsf_X."
        ),
        "source_law": (
            "Delta Z_0=q_(N+1)+j_0 and "
            "Delta Z_A=-s_*'*log((N+1)/a)q_(N+1)+j_A; "
            "only their certified real projections are transferred."
        ),
    }


def build_exact() -> dict:
    components = centered_component_audit()
    transport = mass_transport_audit()
    ordered = ordered_transport_audit()
    endpoint = endpoint_audit()
    null_guard = endpoint_null_guard_audit()
    cutoff = cutoff_audit()
    return {
        "pi_provenance": (
            "The pi in a^2=x/(4*pi)+t/16, T_0, u_(N,x), and the "
            "cutoff recurrence is inherited from the completed-zeta "
            "normalization, Gaussian/theta Fourier normalization, and "
            "Riemann-Siegel saddle. The signed-mass partition, slope "
            "ordering, transport identity, and Gram audit introduce no "
            "new pi, fitted circle, or polygonal definition."
        ),
        "domain": (
            "The identities are fixed-chart algebra and remain valid at "
            "q=1, W_0=0, H_a=0, and endpoint/cutoff equality. The open "
            "outer theorem is required on L>=B_epsilon, "
            "q=2*t*L^2>=1, and 0<tL<=c_*+epsilon."
        ),
        "coordinates": (
            "W_0=e+sum_n z_n=mathsf_X+i*mathsf_Y, "
            "W_A=g+s_*'*sum_n u_n*z_n, s_*'=c+i*b, "
            "mathsf_A=Re(W_A), alpha=c*u_N, and "
            "mathcal_C_N=mathsf_A-alpha*mathsf_X."
        ),
        "centered_components": "; ".join(components.values()) + ".",
        "endpoint_atom": "; ".join(endpoint.values()) + ".",
        "signed_partition": transport["partition"],
        "band_transport": transport["band_identity"],
        "contact_covariance": transport["contact_identity"],
        "ordered_transport": (
            ordered["ordering"]
            + " "
            + ordered["band_identity"]
            + " "
            + ordered["contact_identity"]
        ),
        "conditional_margin": (
            "The exact band identity gives "
            "|mathcal_C_N|>=M*|h_+-h_-|"
            "-(|mathsf_X|/2)*|h_++h_-|-|D_perp|. "
            "Therefore on |mathsf_X|<=delta_L it is sufficient to prove "
            "M*|h_+-h_-|>"
            "A_L+epsilon_term+(delta_L/2)*|h_++h_-|+|D_perp|. "
            "Equivalently, the slope-ordered formula may be bounded "
            "directly before taking absolute values."
        ),
        "gram_energy": (
            "For component amplitude variables a_j and real feature "
            "vectors g_j=(p_j,q_j/ell), "
            "mathsf_X^2+(mathcal_C_N/ell)^2="
            "a^T*(p*p^T+q*q^T/ell^2)*a. "
            "The Gram matrix is positive semidefinite with rank at most "
            "two. On p.a=0 it reduces to (q.a)^2/ell^2, rank at most "
            "one; for m>=3 its contact kernel has dimension at least "
            "m-2. Squaring the scalar proves nonnegativity but supplies "
            "no positive lower bound."
        ),
        "endpoint_null_guard": (
            f"{null_guard['s_prime']}; {null_guard['endpoint']}; "
            f"{null_guard['carriers']}; {null_guard['distances']}; "
            f"{null_guard['amplitudes']} {null_guard['value']}; "
            f"{null_guard['slope']} {null_guard['boundary']}"
        ),
        "first_jet_closure": (
            "Let A=log(a), u_x=1/(8*pi*a^2), "
            "P_1=H_1-A*H_0, Z_0=r_0+H_0, "
            "Z_A=r_A-s_*'*P_1. Then "
            "Z_(0,x)=r_(0,x)-s_*'*H_1+D_0, "
            "P_(1,x)=-s_*'*(H_2-A*H_1)+D_1-A*D_0-u_x*H_0, and "
            "Z_(A,x)=r_(A,x)-s_*''*P_1-s_*'*P_(1,x). "
            "With eta_x=i*omega_eta*eta, "
            "mathsf_X_x=Pi_eta(Z_(0,x))-omega_eta*Im(eta*Z_0), "
            "mathsf_A_x=Pi_eta(Z_(A,x))-omega_eta*Im(eta*Z_A), and "
            "(mathcal_C_N)_x=mathsf_A_x-alpha_x*mathsf_X"
            "-alpha*mathsf_X_x, "
            "alpha_x=Re(s_*'')*u_N+c/(4*T_0). "
            "This retains H_0,H_1,H_2,D_0,D_1 together with "
            "r_0,r_A and all displayed endpoint derivatives."
        ),
        "cutoff_law": " ".join(cutoff.values()),
        "zero_projection_guard": (
            "The labels h_j=d_j/c_j are used only when c_j!=0. "
            "All c_j=0 contributions remain in the division-free "
            "D_perp. At a carrier cosine zero or H_a=0 the mass formula "
            "continues without a pole, while any h-order chart must end "
            "and restart."
        ),
        "legacy_symmetry_verdict": (
            "The legacy radial/transpose symmetry and graph-Laplacian "
            "positivity correspond only to a geometry-side norm or "
            "positive Gram construction. The contact restriction leaves "
            "one signed feature and a large nullspace. The missing input "
            "is arithmetic separation of the positive and negative "
            "effective-slope mass distributions, not another symmetric "
            "smoothing or skeletonization."
        ),
        "source_inequality": (
            "A proof may use either of two equivalent signed targets: "
            "(i) separate the contact-weighted means h_+ and h_- at the "
            "required scale while controlling D_perp, or (ii) order by "
            "effective slope and prove a one-sided cumulative-real-mass "
            "estimate for P_k strong enough that "
            "D_perp+h_m*mathsf_X-sum_k(h_(k+1)-h_k)P_k "
            "cannot enter [-(A_L+epsilon_term),A_L+epsilon_term]. "
            "Both targets include the endpoint atom and use the contact "
            "condition before absolute values."
        ),
        "route_decision": (
            "Retain the contact signed-transport formulation as a "
            "proof-facing diagnostic, not as a solved inequality. "
            "Partition a fixed-N physical chart only at real-component "
            "zeros and effective-slope order changes, keep the endpoint "
            "atom, and seek Xi-specific oscillatory control of the "
            "cumulative masses. Test every candidate first at q=1, "
            "W_0=0, H_a=0, the recurrent endpoint, and both sides of "
            "every cutoff recurrence. Do not infer coercivity from the "
            "positive Gram square."
        ),
        "q_lt_1": (
            "The q<1 layer remains a separate multiplicity-compatible "
            "parabolic/Hermite or boundary-degree theorem. A uniform "
            "positive slope floor down to t=0 would impose endpoint "
            "simplicity and is not required by RH."
        ),
        "proof_boundary": (
            "The terminal-centered component identity, endpoint atom, "
            "rate-free contact covariance, slope-ordered Abel transport, "
            "conditional margin, Gram-rank collapse, endpoint-shaped "
            "generic null guard, full five-current x jet, and adjacent "
            "centered-scalar jump are exact or reproducibly finite. The "
            "null guard is not an actual Xi counterexample. No Xi "
            "cumulative-mass estimate, Abel-scalar gap, successor inner "
            "degree theorem, q<1 closure, finite-height effectivity, "
            "contact exclusion, Lambda<=0, PF-infinity, RH proof, or "
            "Clay-prize conclusion is obtained."
        ),
        "diagnostics": {
            "endpoint_null_guard": null_guard,
            "rank_three_null_vector": (
                "For p=(1,1,1), q=(0,1,2), "
                "a=(1,-2,1) has p.a=q.a=0."
            ),
        },
    }


def build_rows(exact: dict) -> list[ReductionRow]:
    return [
        ReductionRow(
            "cstr_00_pi_provenance",
            "definition_provenance",
            "certified",
            "No geometric visualization introduces a new pi.",
            exact["pi_provenance"],
            "All pi factors remain source-traced.",
        ),
        ReductionRow(
            "cstr_01_domain",
            "exact_domain",
            "certified",
            "The algebra includes every closed-edge exceptional fibre.",
            exact["domain"],
            "The arithmetic inequality is still open.",
        ),
        ReductionRow(
            "cstr_02_coordinates",
            "exact_normal_form",
            "ready_to_apply",
            "The Abel scalar is the terminal-centered physical slope.",
            exact["coordinates"],
            "No ratio by W_0 is used.",
        ),
        ReductionRow(
            "cstr_03_centered_components",
            "exact_decomposition",
            "ready_to_apply",
            "Each endpoint/carrier contributes one real mass and one centered slope.",
            exact["centered_components"],
            "All endpoint and imaginary carrier terms remain.",
        ),
        ReductionRow(
            "cstr_04_endpoint_atom",
            "exact_endpoint",
            "ready_to_apply",
            "The recurrent endpoint is one effective-slope atom with a division-free zero branch.",
            exact["endpoint_atom"],
            "No sign for the corrected real endpoint quotient is asserted.",
        ),
        ReductionRow(
            "cstr_05_signed_partition",
            "exact_partition",
            "ready_to_apply",
            "The contact condition balances positive and negative real mass.",
            exact["signed_partition"],
            "Zero-real components are retained separately.",
        ),
        ReductionRow(
            "cstr_06_band_transport",
            "exact_identity",
            "ready_to_apply",
            "The Abel scalar is a near-contact signed covariance.",
            exact["band_transport"],
            "The identity alone has no sign.",
        ),
        ReductionRow(
            "cstr_07_contact_covariance",
            "exact_identity",
            "ready_to_apply",
            "At contact only the effective-slope mean gap and perpendicular terms remain.",
            exact["contact_covariance"],
            "This includes W_0=0.",
        ),
        ReductionRow(
            "cstr_08_ordered_transport",
            "exact_reindexing",
            "ready_to_apply",
            "Slope ordering turns the covariance into positive gaps times signed cumulative masses.",
            exact["ordered_transport"],
            "The cumulative masses remain the arithmetic unknown.",
        ),
        ReductionRow(
            "cstr_09_conditional_margin",
            "conditional_theorem",
            "not_ready_to_apply",
            "A quantitative mean separation is sufficient for the Abel gap.",
            exact["conditional_margin"],
            "No Xi separation estimate is proved.",
        ),
        ReductionRow(
            "cstr_10_gram_rank",
            "nonpromotion_guard",
            "guard_validated",
            "Contact conditioning collapses the positive Gram energy to rank one.",
            exact["gram_energy"],
            "Symmetry and PSD alone cannot prove coercivity.",
            exact["diagnostics"],
        ),
        ReductionRow(
            "cstr_11_endpoint_null_guard",
            "countermodel",
            "guard_validated",
            "Even an endpoint-shaped anchored model can hit W_0=0 and mathcal_C_N=0.",
            exact["endpoint_null_guard"],
            "The model is synthetic, not actual Xi.",
            exact["diagnostics"]["endpoint_null_guard"],
        ),
        ReductionRow(
            "cstr_12_first_jet",
            "exact_first_jet",
            "ready_to_apply",
            "The transport coordinate keeps H_2,D_0,D_1 and every endpoint derivative.",
            exact["first_jet_closure"],
            "No derivative term is deleted.",
        ),
        ReductionRow(
            "cstr_13_cutoff",
            "exact_boundary_law",
            "ready_to_apply",
            "The centered scalar has an explicit adjacent-cutoff jump.",
            exact["cutoff_law"],
            "Real projection bounds are not promoted to a complex jump.",
        ),
        ReductionRow(
            "cstr_14_zero_projection",
            "proof_guard",
            "guard_validated",
            "The rate-free formulation remains valid when an effective slope label diverges.",
            exact["zero_projection_guard"],
            "Order cells must be joined division-free.",
        ),
        ReductionRow(
            "cstr_15_legacy_symmetry",
            "route_decision",
            "guard_validated",
            "The old symmetry identifies a norm, not the missing arithmetic coercivity.",
            exact["legacy_symmetry_verdict"],
            "The legacy field is not substituted for Xi coefficients.",
        ),
        ReductionRow(
            "cstr_16_source_inequality",
            "open_theorem_target",
            "not_ready_to_apply",
            "The missing theorem is now a signed mass-transport estimate.",
            exact["source_inequality"],
            "Either form must be uniform on the physical closed collar.",
        ),
        ReductionRow(
            "cstr_17_route_decision",
            "route_decision",
            "ready_to_apply",
            "Future work should attack physical order cells and cumulative masses.",
            exact["route_decision"],
            "No broad geometry-only or brute-force sweep is promoted.",
        ),
        ReductionRow(
            "cstr_18_q_lt_1",
            "open_theorem_target",
            "not_ready_to_apply",
            "The small-q multiplicity layer remains separate.",
            exact["q_lt_1"],
            "No endpoint simplicity assumption is introduced.",
        ),
        ReductionRow(
            "cstr_19_boundary",
            "proof_guard",
            "guard_validated",
            "The exact reduction is not promoted beyond its evidence.",
            exact["proof_boundary"],
            "RH and every prize-level conclusion remain open.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-28",
        "status": (
            "exact endpoint-complete contact signed-transport reduction, "
            "conditional margin, Gram-rank nonpromotion guard, and "
            "cutoff/first-jet closure; the Xi cumulative-mass theorem "
            "remains open"
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
            "# Newman Contact Signed-Transport Reduction",
            "",
            "Date: 2026-07-28",
            "",
            "Status: exact endpoint-complete reduction and route guard;",
            "not a proof of the Xi contact gap, `Lambda<=0`, RH,",
            "PF-infinity, or a Clay-prize result.",
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
            "## Centered Components",
            "",
            "```text",
            exact["coordinates"],
            exact["centered_components"],
            exact["endpoint_atom"],
            "```",
            "",
            "## Contact Transport",
            "",
            "```text",
            exact["signed_partition"],
            exact["band_transport"],
            exact["contact_covariance"],
            "```",
            "",
            "The contact equation is therefore a signed covariance, not",
            "a geometry-only positive energy.",
            "",
            "## Slope-Ordered Form",
            "",
            "```text",
            exact["ordered_transport"],
            "```",
            "",
            "The positive gaps are exact. Their cumulative real masses",
            "are the open arithmetic input.",
            "",
            "## Conditional Margin",
            "",
            "```text",
            exact["conditional_margin"],
            "```",
            "",
            "## Gram Rank Guard",
            "",
            exact["gram_energy"],
            "",
            "## Endpoint-Shaped Null Guard",
            "",
            "```text",
            exact["endpoint_null_guard"],
            "```",
            "",
            "## First Jet",
            "",
            "```text",
            exact["first_jet_closure"],
            "```",
            "",
            "## Cutoff Law",
            "",
            "```text",
            exact["cutoff_law"],
            "```",
            "",
            exact["zero_projection_guard"],
            "",
            "## Legacy Symmetry Verdict",
            "",
            exact["legacy_symmetry_verdict"],
            "",
            "## Open Source Inequality",
            "",
            exact["source_inequality"],
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
        "built Newman contact signed-transport reduction: "
        "20 rows, 4 exact transport/endpoint identities, "
        "1 conditional margin, 1 Gram-rank guard, "
        "1 endpoint-shaped null guard, 2 open theorem targets"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
