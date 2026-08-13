#!/usr/bin/env python3
"""Build the physical P_lin joined-coefficient envelope gate."""

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
    "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
    "physical_plin_joined_coefficient_envelope_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "critical_ray": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "critical_ray_finite_height_real_edge_gate.json"
    ),
    "prefix_flux": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
    "two_carrier": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_two_carrier_"
        "kernel_reduction.json"
    ),
    "physical_plin": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_amplitude_gate.json"
    ),
    "correction_envelope": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_morse_fresnel_"
        "physical_plin_correction_coefficient_envelope_gate.json"
    ),
}

H_DENOMINATOR = 72_000_000_000


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.cancel(sp.expand(expression)) != 0:
        raise RuntimeError(label)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    critical = payloads["critical_ray"]
    domain = critical.get("exact", {}).get("domain", {})
    if "0<t<=1/2" not in domain.get("effective_box", ""):
        raise RuntimeError("critical-ray t cap drifted")
    if "h<exp(-25)<1/72000000000" not in domain.get("effective_box", ""):
        raise RuntimeError("critical-ray h cap drifted")
    saddle = critical.get("exact", {}).get("terminal", {}).get("alpha", "")
    for token in (
        "s_*'=-i/2-it alpha'/4",
        "s_*''=-t alpha''/8",
    ):
        if token not in saddle:
            raise RuntimeError(f"critical-ray saddle identity drifted: {token}")
    critical_bounds = critical.get("majorant_certificate", {}).get(
        "source_bounds", {}
    )
    if critical_bounds.get("alpha_prime") != "|alpha'|<h^2":
        raise RuntimeError("critical-ray alpha-prime source drifted")
    if critical_bounds.get("alpha_second") != "|alpha''|<h^4":
        raise RuntimeError("critical-ray alpha-second source drifted")

    flux = payloads["prefix_flux"].get("exact", {}).get(
        "normalized_flux_derivative", ""
    )
    if "u_(N,x)=1/(8*pi*a^2)=1/(4*T_0)" not in flux:
        raise RuntimeError("prefix-flux u_(N,x) identity drifted")

    two_carrier = payloads["two_carrier"].get("finite_sum_certificate", {})
    radial = two_carrier.get("physical_radial_defect_bound", "")
    for token in (
        "u_N=-log(1-h*theta)",
        "|s_*'+i/2|<h^2/6000",
        "0<=u_N-h*theta<h^2",
        "|i*b*u_N-chi_N|<4h^2",
    ):
        if token not in radial:
            raise RuntimeError(f"two-carrier physical source drifted: {token}")

    plin_counts = payloads["physical_plin"].get("counts", {})
    if plin_counts.get("exact_centered_complex_coefficients") != 6:
        raise RuntimeError("physical P_lin coefficient count drifted")
    exact = payloads["physical_plin"].get("coefficient_certificate", {})
    if "p_5=i*rho_2*(X_T+V_p)*(s_*')^2" not in exact.get(
        "leading_coefficient", ""
    ):
        raise RuntimeError("physical P_lin top factor drifted")

    correction = payloads["correction_envelope"]
    expected = {
        "c_0": "|c_0-1|<4379h^2<1/3; hence 2/3<|c_0|<4/3",
        "c_1": "|c_1|<h^2/4",
        "c_2": "|c_2|<h^2/16",
    }
    value_bounds = correction.get("bound_certificate", {}).get(
        "centered_value_bounds", {}
    )
    if value_bounds != expected:
        raise RuntimeError("centered C envelope drifted")
    expected_d = {
        "d_0": "|d_0|<16893h^4",
        "d_1": "|d_1|<5h^4/16",
        "d_2": "|d_2|<h^4/16",
    }
    derivative_bounds = correction.get("bound_certificate", {}).get(
        "centered_derivative_bounds", {}
    )
    if derivative_bounds != expected_d:
        raise RuntimeError("centered D envelope drifted")

    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def symbolic_certificate() -> dict:
    y = sp.symbols("y")
    r_0, r_1, t_0, t_1, delta, u_x = sp.symbols(
        "r_0 r_1 t_0 t_1 delta u_x"
    )
    f_v, f_n, f_a, f_q = sp.symbols("f_v f_n f_a f_q", real=True)
    R = r_0 + r_1 * y
    R_x = t_0 + t_1 * y
    G = sp.expand(R * (R + delta) + R_x)
    U = sp.expand(f_n * G + f_q * (R + delta) + f_a * R + f_v)
    V = sp.expand(f_n * R + f_q)
    U_poly = sp.Poly(U, y)
    V_poly = sp.Poly(V, y)

    u = [U_poly.coeff_monomial(y**j) for j in range(3)]
    v = [V_poly.coeff_monomial(y**j) for j in range(2)]
    expected_u = [
        f_n * (r_0 * (r_0 + delta) + t_0)
        + f_q * (r_0 + delta)
        + f_a * r_0
        + f_v,
        f_n * (r_1 * (2 * r_0 + delta) + t_1)
        + (f_q + f_a) * r_1,
        f_n * r_1**2,
    ]
    expected_v = [f_n * r_0 + f_q, f_n * r_1]
    for index, (actual, expected) in enumerate(zip(u, expected_u)):
        require_zero(actual - expected, f"U row {index} drifted")
    for index, (actual, expected) in enumerate(zip(v, expected_v)):
        require_zero(actual - expected, f"V row {index} drifted")

    ideal_subs = {
        r_0: 0,
        r_1: sp.I / 2,
        t_0: -sp.I * u_x / 2,
        t_1: 0,
        delta: 0,
    }
    ideal_u = [sp.expand(item.subs(ideal_subs)) for item in u]
    ideal_v = [sp.expand(item.subs(ideal_subs)) for item in v]
    required_ideal_u = [
        f_v - sp.I * u_x * f_n / 2,
        sp.I * (f_q + f_a) / 2,
        -f_n / 4,
    ]
    required_ideal_v = [f_q, sp.I * f_n / 2]
    for index, (actual, expected) in enumerate(
        zip(ideal_u, required_ideal_u)
    ):
        require_zero(actual - expected, f"ideal U row {index} drifted")
    for index, (actual, expected) in enumerate(
        zip(ideal_v, required_ideal_v)
    ):
        require_zero(actual - expected, f"ideal V row {index} drifted")

    alpha = sp.symbols("alpha_V alpha_N alpha_A alpha_Q", real=True)
    beta = sp.symbols("beta_V beta_N beta_A beta_Q", real=True)

    def row_coefficients(row: tuple[sp.Symbol, ...]) -> tuple[list, list]:
        substitutions = dict(zip((f_v, f_n, f_a, f_q), row))
        return (
            [sp.expand(item.subs(substitutions)) for item in u],
            [sp.expand(item.subs(substitutions)) for item in v],
        )

    u_alpha, v_alpha = row_coefficients(alpha)
    u_beta, v_beta = row_coefficients(beta)
    gamma = [
        u_alpha[0],
        u_alpha[1] + sp.I * u_beta[0],
        u_alpha[2] + sp.I * u_beta[1],
        sp.I * u_beta[2],
    ]
    eta = [
        v_alpha[0],
        v_alpha[1] + sp.I * v_beta[0],
        sp.I * v_beta[1],
    ]
    ideal_gamma = [sp.expand(item.subs(ideal_subs)) for item in gamma]
    ideal_eta = [sp.expand(item.subs(ideal_subs)) for item in eta]
    required_gamma = [
        alpha[0] - sp.I * u_x * alpha[1] / 2,
        sp.I * (alpha[3] + alpha[2]) / 2
        + sp.I * beta[0]
        + u_x * beta[1] / 2,
        -alpha[1] / 4 - (beta[3] + beta[2]) / 2,
        -sp.I * beta[1] / 4,
    ]
    required_eta = [
        alpha[3],
        sp.I * (alpha[1] / 2 + beta[3]),
        -beta[1] / 2,
    ]
    for index, (actual, expected) in enumerate(
        zip(ideal_gamma, required_gamma)
    ):
        require_zero(actual - expected, f"ideal gamma {index} drifted")
    for index, (actual, expected) in enumerate(zip(ideal_eta, required_eta)):
        require_zero(actual - expected, f"ideal eta {index} drifted")

    c = sp.symbols("c_0:3")
    d = sp.symbols("d_0:3")
    C = sum(c[j] * y**j for j in range(3))
    D = sum(d[j] * y**j for j in range(3))
    P = sp.expand(
        C * sum(gamma[j] * y**j for j in range(4))
        + D * sum(eta[j] * y**j for j in range(3))
    )
    p = [sp.Poly(P, y).coeff_monomial(y**j) for j in range(6)]
    expected_p = [
        c[0] * gamma[0] + d[0] * eta[0],
        c[0] * gamma[1]
        + c[1] * gamma[0]
        + d[0] * eta[1]
        + d[1] * eta[0],
        c[0] * gamma[2]
        + c[1] * gamma[1]
        + c[2] * gamma[0]
        + d[0] * eta[2]
        + d[1] * eta[1]
        + d[2] * eta[0],
        c[0] * gamma[3]
        + c[1] * gamma[2]
        + c[2] * gamma[1]
        + d[1] * eta[2]
        + d[2] * eta[1],
        c[1] * gamma[3] + c[2] * gamma[2] + d[2] * eta[2],
        c[2] * gamma[3],
    ]
    for index, (actual, expected) in enumerate(zip(p, expected_p)):
        require_zero(actual - expected, f"p_{index} convolution drifted")

    rho_2, s_prime = sp.symbols("rho_2 s_prime")
    p_5 = p[5].subs({c[2]: rho_2, r_1: -s_prime})
    require_zero(
        p_5 - sp.I * rho_2 * beta[1] * s_prime**2,
        "physical p_5 factorization drifted",
    )

    return {
        "carrier_rates": (
            "r_0=-c*u_N, r_1=-s_*', t_0=-c_x*u_N+i*b*u_(N,x), "
            "t_1=-s_*'', and delta=chi_N-i*b*u_N."
        ),
        "ideal_rate_point": (
            "(r_0,r_1,t_0,t_1,delta)="
            "(0,i/2,-i*u_(N,x)/2,0,0)."
        ),
        "ideal_row_coefficients": {
            "u_f": (
                "bar(u)_(f,0)=f_V-i*u_(N,x)*f_N/2, "
                "bar(u)_(f,1)=i(f_Q+f_A)/2, bar(u)_(f,2)=-f_N/4."
            ),
            "v_f": "bar(v)_(f,0)=f_Q, bar(v)_(f,1)=i*f_N/2.",
        },
        "ideal_gamma": [
            "bar(gamma)_0=alpha_V-i*u_(N,x)*alpha_N/2",
            "bar(gamma)_1=i(alpha_Q+alpha_A)/2+i*beta_V+u_(N,x)*beta_N/2",
            "bar(gamma)_2=-alpha_N/4-(beta_Q+beta_A)/2",
            "bar(gamma)_3=-i*beta_N/4",
        ],
        "ideal_eta": [
            "bar(eta)_0=alpha_Q",
            "bar(eta)_1=i(alpha_N/2+beta_Q)",
            "bar(eta)_2=-beta_N/2",
        ],
        "physical_rows": (
            "alpha=(N_(p,xi),V_(p,xi),-Q_(p,xi),-A_(p,xi)) and "
            "beta=(N_p+A_(T,x),X_T+V_p,-Q_p-X_(T,x),-A_T-A_p)."
        ),
        "p_convolution": [str(sp.expand(item)) for item in expected_p],
        "p_5": "p_5=i*rho_2*beta_N*(s_*')^2=i*rho_2*(X_T+V_p)*(s_*')^2",
        "symbolic_zero_checks": {
            "row_recurrences": 5,
            "ideal_row_coefficients": 5,
            "ideal_gamma_eta": 7,
            "p_convolution": 6,
            "physical_p5": 1,
        },
    }


def bound_certificate() -> dict:
    h_cap = Fraction(1, H_DENOMINATOR)
    r_1_abs = Fraction(3001, 6000)

    t_0_at_half = Fraction(1, 64) + Fraction(3001, 144000)
    if not t_0_at_half < Fraction(1, 16):
        raise RuntimeError("t_0 rounding failed")

    g_0_at_half = (
        Fraction(1, 16)
        + Fraction(1, 144000)
        + Fraction(1, 1500)
        + Fraction(1, 36_000_000)
    )
    if not g_0_at_half < Fraction(1, 15):
        raise RuntimeError("G_0 deviation rounding failed")

    g_1_raw = r_1_abs * Fraction(6001, 1500) + Fraction(1, 16)
    if not g_1_raw < Fraction(17, 8):
        raise RuntimeError("G_1 rounding failed")

    g_2_raw = Fraction(1, 6000) + Fraction(1, 36_000_000)
    if not g_2_raw < Fraction(1, 3000):
        raise RuntimeError("G_2 deviation rounding failed")

    p_5_raw = Fraction(1, 16) * r_1_abs**2
    if not p_5_raw < Fraction(1, 63):
        raise RuntimeError("p_5 rounding failed")

    return {
        "domain": (
            "L>=50, 0<tL<=25, 0<t<=1/2, 0<=theta<=1, "
            "h=1/a<1/72000000000, and 3<pi<22/7."
        ),
        "rate_box": {
            "r_0": "|r_0|<h^3/3000",
            "r_1": (
                "|r_1-i/2|<h^2/6000 and |r_1|<3001/6000"
            ),
            "t_0": "|t_0|<h^2/16",
            "t_1": "|t_1|<h^4/16",
            "delta": "|delta|<4h^2",
        },
        "quadratic_rate_box": {
            "G_0": "|G_0+i*u_(N,x)/2|<h^4/15",
            "G_1": "|G_1|<17h^2/8",
            "G_2": "|G_2+1/4|<h^2/3000",
        },
        "row_error_envelopes": {
            "e_u0": (
                "e_(u,0)(f)=h^4|f_N|/15+"
                "(4h^2+h^3/3000)|f_Q|+h^3|f_A|/3000"
            ),
            "e_u1": (
                "e_(u,1)(f)=17h^2|f_N|/8+"
                "h^2(|f_Q|+|f_A|)/6000"
            ),
            "e_u2": "e_(u,2)(f)=h^2|f_N|/3000",
            "e_v0": "e_(v,0)(f)=h^3|f_N|/3000",
            "e_v1": "e_(v,1)(f)=h^2|f_N|/6000",
        },
        "joined_errors": {
            "gamma": (
                "E_gamma=(e_(u,0)(alpha), "
                "e_(u,1)(alpha)+e_(u,0)(beta), "
                "e_(u,2)(alpha)+e_(u,1)(beta), e_(u,2)(beta))."
            ),
            "eta": (
                "E_eta=(e_(v,0)(alpha), "
                "e_(v,1)(alpha)+e_(v,0)(beta), e_(v,1)(beta))."
            ),
            "majorants": (
                "Gamma_j=|bar(gamma)_j|+E_(gamma,j), "
                "H_j=|bar(eta)_j|+E_(eta,j). Then "
                "|gamma_j|<Gamma_j and |eta_j|<H_j."
            ),
        },
        "p_deviation_envelopes": {
            "Pi_0": (
                "Pi_0=E_(gamma,0)+4379h^2Gamma_0+16893h^4H_0"
            ),
            "Pi_1": (
                "Pi_1=E_(gamma,1)+4379h^2Gamma_1+h^2Gamma_0/4+"
                "16893h^4H_1+5h^4H_0/16"
            ),
            "Pi_2": (
                "Pi_2=E_(gamma,2)+4379h^2Gamma_2+h^2Gamma_1/4+"
                "h^2Gamma_0/16+16893h^4H_2+5h^4H_1/16+h^4H_0/16"
            ),
            "Pi_3": (
                "Pi_3=E_(gamma,3)+4379h^2Gamma_3+h^2Gamma_2/4+"
                "h^2Gamma_1/16+5h^4H_2/16+h^4H_1/16"
            ),
            "Pi_4": "Pi_4=h^2Gamma_3/4+h^2Gamma_2/16+h^4H_2/16",
            "Pi_5": (
                "p_5=i*rho_2*beta_N*(s_*')^2 and "
                "|p_5|<h^2|beta_N|/63"
            ),
        },
        "p_conclusions": (
            "|p_j-bar(gamma)_j|<Pi_j for 0<=j<=3, |p_4|<Pi_4, "
            "and p_5 has the exact displayed fibre factor."
        ),
        "rational_diagnostics": {
            "h_upper": fraction_text(h_cap),
            "r_1_absolute_cap": fraction_text(r_1_abs),
            "t_0_raw_at_h_half": fraction_text(t_0_at_half),
            "G_0_raw_at_h_half": fraction_text(g_0_at_half),
            "G_1_raw": fraction_text(g_1_raw),
            "G_2_raw": fraction_text(g_2_raw),
            "p_5_raw_over_h2_betaN": fraction_text(p_5_raw),
        },
    }


def build_rows(symbolic: dict, bounds: dict) -> list[GateRow]:
    return [
        GateRow(
            "pjce_01_domain",
            "effective_domain",
            "ready_to_apply",
            "The joined envelope uses the full physical critical chart.",
            bounds["domain"],
            "No extension outside this source chart is asserted.",
        ),
        GateRow(
            "pjce_02_rates",
            "exact_identity",
            "certified",
            "All five carrier rates are terminal-centered exact quantities.",
            symbolic["carrier_rates"],
            "The physical height derivative is taken inside one fixed-N chart.",
        ),
        GateRow(
            "pjce_03_rate_box",
            "analytic_envelope",
            "certified",
            "The five physical carrier rates have explicit h-scale bounds.",
            "; ".join(bounds["rate_box"].values()),
            "These are coefficient bounds, not observation-row bounds.",
        ),
        GateRow(
            "pjce_04_quadratic_rates",
            "analytic_envelope",
            "certified",
            "The three quadratic rate coefficients stay near the ideal cubic chart.",
            "; ".join(bounds["quadratic_rate_box"].values()),
            "No moment or terminal amplitude is estimated here.",
        ),
        GateRow(
            "pjce_05_ideal_rows",
            "exact_specialization",
            "certified",
            "Every real row has an exact ideal-rate U/V coefficient vector.",
            (
                symbolic["ideal_row_coefficients"]["u_f"]
                + " "
                + symbolic["ideal_row_coefficients"]["v_f"]
            ),
            "The row entries remain the actual joined physical observations.",
        ),
        GateRow(
            "pjce_06_row_errors",
            "row_sensitive_envelope",
            "certified",
            "Each U/V departure is bounded only after the row is formed.",
            "; ".join(bounds["row_error_envelopes"].values()),
            "No separate numerical caps for alpha or beta components are inserted.",
        ),
        GateRow(
            "pjce_07_physical_rows",
            "joined_observation_definition",
            "certified",
            "The propagation uses the two physical rows from the flow identity.",
            symbolic["physical_rows"],
            "The terminal coordinates remain inside beta.",
        ),
        GateRow(
            "pjce_08_ideal_gamma",
            "joined_exact_identity",
            "certified",
            "The four ideal cubic coefficients preserve physical row cancellation.",
            "; ".join(symbolic["ideal_gamma"]),
            "The correction-free observation identities are not imposed on the actual rows.",
        ),
        GateRow(
            "pjce_09_ideal_eta",
            "joined_exact_identity",
            "certified",
            "The derivative channel has three exact ideal coefficients.",
            "; ".join(symbolic["ideal_eta"]),
            "D remains present in the physical polynomial.",
        ),
        GateRow(
            "pjce_10_joined_errors",
            "joined_envelope",
            "certified",
            "The alpha/beta row errors compose into gamma and eta before propagation.",
            (
                bounds["joined_errors"]["gamma"]
                + " "
                + bounds["joined_errors"]["eta"]
            ),
            "These are symbolic row-sensitive envelopes, not uniform observation bounds.",
        ),
        GateRow(
            "pjce_11_majorants",
            "joined_envelope",
            "certified",
            "Gamma and H majorize the exact gamma and eta coefficients.",
            bounds["joined_errors"]["majorants"],
            "Cancellation inside each bar(gamma)_j is retained until its modulus.",
        ),
        GateRow(
            "pjce_12_p0_p3",
            "coefficient_propagation",
            "certified",
            "The four leading coefficients remain within explicit errors of the ideal cubic.",
            "; ".join(
                bounds["p_deviation_envelopes"][f"Pi_{j}"]
                for j in range(4)
            ),
            "No claim is made that the ideal cubic coefficients are themselves small.",
        ),
        GateRow(
            "pjce_13_p4",
            "coefficient_propagation",
            "certified",
            "The degree-four channel is correction-suppressed.",
            bounds["p_deviation_envelopes"]["Pi_4"],
            "Its physical row scale remains explicit.",
        ),
        GateRow(
            "pjce_14_p5",
            "fibre_factorization",
            "certified",
            "The degree-five channel retains the exact total-value factor.",
            bounds["p_deviation_envelopes"]["Pi_5"],
            "The vanishing fibre is a degree drop, not a contact theorem.",
        ),
        GateRow(
            "pjce_15_basis_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Centered y coefficients are not silently inserted into the old lambda-moment tail template.",
            (
                "The p_j multiply y=lambda-log(a). Any use of the Section 11.172 "
                "K_j template for lambda^j must first perform the exact binomial "
                "basis translation or rederive centered tail weights."
            ),
            "This prevents an untracked log(a) loss.",
        ),
        GateRow(
            "pjce_16_observation_handoff",
            "route_decision",
            "open_quantitative_handoff",
            "The remaining coefficient input is one retained-carrier observation envelope.",
            (
                "Prove a source-specific joined bound for bar(gamma), Gamma, and H, "
                "or estimate the grouped Morse functional directly with these rows "
                "left inside its phase sum."
            ),
            "Phi_B or one scalar current value does not supply this bound automatically.",
        ),
        GateRow(
            "pjce_17_boundary",
            "proof_boundary",
            "guard_validated",
            "The joined coefficient theorem is not promoted to a flow estimate.",
            (
                "There are zero numerical observation-row bounds, grouped Morse "
                "bounds, quadratic residual bounds, Phi_B bounds, or Xi signs here."
            ),
            "No contact exclusion, Lambda<=0, PF-infinity, RH, or prize theorem is claimed.",
        ),
    ]


def render_note(artifact: dict) -> str:
    symbolic = artifact["symbolic_certificate"]
    bounds = artifact["bound_certificate"]
    p_bounds = bounds["p_deviation_envelopes"]
    return "\n".join(
        [
            "# Physical P_lin Joined-Coefficient Envelope Gate",
            "",
            "Date: 2026-08-02",
            "",
            "Status: five physical carrier-rate bounds, three ideal-rate deviation bounds, and joined alpha/beta propagation through all six centered P_lin coefficients. This is not a proof of a grouped Morse estimate, signed flow, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Carrier Rates",
            "",
            bounds["domain"],
            "",
            symbolic["carrier_rates"],
            "",
            "The source identities u_N=-log(1-h theta), u_(N,x)=h^2/(8pi), s_*'=-i/2-it alpha'/4, and s_*''=-t alpha''/8 give",
            "",
            "```text",
            *bounds["rate_box"].values(),
            "```",
            "",
            "Here pi is the inherited completed-zeta/Riemann-Siegel constant in the saddle coordinate. Only pi>3 is used to obtain u_(N,x)<h^2/24; no new circle or polygon is introduced.",
            "",
            "For G_0=r_0(r_0+delta)+t_0, G_1=r_1(2r_0+delta)+t_1, and G_2=r_1^2, the ideal point is (-i u_(N,x)/2,0,-1/4), and",
            "",
            "```text",
            *bounds["quadratic_rate_box"].values(),
            "```",
            "",
            "## Row-Sensitive Errors",
            "",
            "For a real row f=(f_V,f_N,f_A,f_Q), define the ideal coefficients",
            "",
            "```text",
            symbolic["ideal_row_coefficients"]["u_f"],
            symbolic["ideal_row_coefficients"]["v_f"],
            "```",
            "",
            "and the nonnegative error envelopes",
            "",
            "```text",
            *bounds["row_error_envelopes"].values(),
            "```",
            "",
            "Then |u_(f,j)-bar(u)_(f,j)|<e_(u,j)(f) and |v_(f,j)-bar(v)_(f,j)|<e_(v,j)(f). The row is formed before these moduli are taken.",
            "",
            "## Joined Rows",
            "",
            symbolic["physical_rows"],
            "",
            "At the ideal rate point their composed coefficients are",
            "",
            "```text",
            *symbolic["ideal_gamma"],
            *symbolic["ideal_eta"],
            "```",
            "",
            "The actual errors compose as",
            "",
            "```text",
            bounds["joined_errors"]["gamma"],
            bounds["joined_errors"]["eta"],
            bounds["joined_errors"]["majorants"],
            "```",
            "",
            "No independent numerical bound for an alpha or beta component has been assumed.",
            "",
            "## Six Coefficients",
            "",
            "Combining these joined rows with the certified C/D box from Section 11.174 gives",
            "",
            "```text",
            p_bounds["Pi_0"],
            p_bounds["Pi_1"],
            p_bounds["Pi_2"],
            p_bounds["Pi_3"],
            p_bounds["Pi_4"],
            p_bounds["Pi_5"],
            "```",
            "",
            bounds["p_conclusions"],
            "",
            "Thus the actual polynomial is a row-sensitive O(h^2) perturbation of its ideal cubic in coefficient space, while the top channel keeps its exact total-value factor. This does not say the retained rows or the cubic coefficients are small.",
            "",
            "## Basis Guard",
            "",
            "The p_j above multiply y=lambda-log(a). The K_j outside-tail template in Section 11.172 was stated for lambda powers. It may be used only after exact binomial translation or after deriving centered tail weights; substituting the p_j directly would hide powers of log(a).",
            "",
            "## Next Target",
            "",
            "Either prove a source-specific joined envelope for bar(gamma), Gamma, and H from the endpoint-composed retained carrier, or leave those rows inside the grouped Morse sum and prove a coefficient-aware cancellation estimate directly. A scalar Phi_B estimate, even once available, cannot by itself be treated as a gradient-row norm.",
            "",
            "## Proof Boundary",
            "",
            "This gate proves no numerical retained-observation bound, centered-to-lambda tail bound, grouped c-prime sum, reciprocal cancellation estimate, quadratic residual bound, h^2 flow error, Phi_B upper bound, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
            "",
            "## Reproduce",
            "",
            "```powershell",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
        ]
    )


def build_artifact() -> dict:
    payloads = load_sources()
    audit = source_audit(payloads)
    symbolic = symbolic_certificate()
    bounds = bound_certificate()
    rows = build_rows(symbolic, bounds)
    return {
        "kind": STEM,
        "date": "2026-08-02",
        "status": (
            "certified physical rate box and joined alpha/beta propagation "
            "through all six centered P_lin coefficients"
        ),
        "source_audit": audit,
        "symbolic_certificate": symbolic,
        "bound_certificate": bounds,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "sources": len(audit),
            "carrier_rate_bounds": 5,
            "quadratic_rate_bounds": 3,
            "row_error_envelopes": 5,
            "ideal_joined_coefficients": 7,
            "centered_plin_coefficients_propagated": 6,
            "exact_top_fibre_factors": 1,
            "numerical_observation_row_bounds": 0,
            "grouped_interior_bounds": 0,
            "signed_flow_bounds": 0,
        },
        "proof_boundary": (
            "Rate and row-sensitive centered-coefficient envelopes only; no "
            "numerical retained-observation bound, basis-translated tail bound, "
            "grouped Morse estimate, quadratic residual estimate, signed flow "
            "theorem, Phi_B bound, contact exclusion, Xi theorem, Lambda<=0, "
            "PF-infinity, RH, or prize-level conclusion."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()

    artifact = build_artifact()
    atomic_write(args.output, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(artifact))
    counts = artifact["counts"]
    print(
        "built physical P_lin joined-coefficient envelope gate: "
        f"{counts['rows']} rows, {counts['carrier_rate_bounds']} rate bounds, "
        f"{counts['quadratic_rate_bounds']} quadratic-rate bounds, "
        f"{counts['centered_plin_coefficients_propagated']} propagated p_j, "
        f"{counts['numerical_observation_row_bounds']} numerical row bounds, "
        f"{counts['signed_flow_bounds']} signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
