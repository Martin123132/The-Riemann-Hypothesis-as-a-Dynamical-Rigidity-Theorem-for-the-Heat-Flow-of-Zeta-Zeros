#!/usr/bin/env python3
"""Build the normalized-prefix and first-jet phase-flux reduction."""

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
    "normalized_prefix_phase_flux_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "carrier_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "carrier_kernel_abel_prefix_reduction.json"
    ),
    "phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "direct_projection_regime_reduction.json"
    ),
    "wronskian_crossing": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_wronskian_crossing_reduction.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
}

L_MIN = 50
D_X_CONSTANT = 4_223
D_CONSTANT = 2_189


@dataclass(frozen=True)
class FluxRow:
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
        "carrier_prefix": (
            "pi is the universal constant",
            "mathcal_C_N=Re(V_N)-b*u_N*mathsf_Y",
            "Delta x=4*pi*(2N+1)",
        ),
        "phase_anchor": (
            "q_n=f_n/f_1",
            "q_1=1",
            "Z_0=P_0+r_0=E_[1]/f_1",
        ),
        "direct_projection": (
            "eta*q_n=r_n*zeta_n",
            "W_0=eta*Z_0",
            "W_A=eta*Z_A",
        ),
        "wronskian_crossing": (
            "N_up=N_(X=0,E_[1]!=0",
            "W_[1]=V*X-U*Y",
        ),
        "oriented_successor": (
            "kappa_j=N_up(t_j;R_j)-N_up(t_(j+1);R_j)",
            "wind(proxy_j)<1",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(value.get("kind", "")) for key, value in payloads.items()}


def symbolic_audit() -> dict[str, str]:
    alpha = sp.symbols("alpha", real=True)
    eta = sp.cos(alpha) + sp.I * sp.sin(alpha)
    eta_bar = sp.conjugate(eta)
    if sp.trigsimp(eta * eta_bar - 1) != 0:
        raise RuntimeError("unit-anchor identity failed")

    q1r, q1i, q2r, q2i = sp.symbols(
        "q1r q1i q2r q2i", real=True
    )
    q1 = q1r + sp.I * q1i
    q2 = q2r + sp.I * q2i
    normalized_prefix = q1 + q2
    physical_prefix = eta * q1 + eta * q2
    if sp.trigsimp(eta_bar * physical_prefix - normalized_prefix) != 0:
        raise RuntimeError("normalized-prefix identity failed")

    c, b, u = sp.symbols("c b u", real=True)
    r0r, r0i, rar, rai = sp.symbols(
        "r0r r0i rar rai", real=True
    )
    r0 = r0r + sp.I * r0i
    ra = rar + sp.I * rai
    s_prime = c + sp.I * b
    g = normalized_prefix
    z0 = r0 + g
    v_tilde = ra - s_prime * u * r0
    b_shape = v_tilde + s_prime * g
    za = s_prime * u * z0 + b_shape
    projected_value = sp.re(sp.expand_complex(eta * z0))
    projected_slope = sp.re(sp.expand_complex(eta * za))
    projected_imag = sp.im(sp.expand_complex(eta * z0))
    projected_b = sp.re(sp.expand_complex(eta * b_shape))
    contact = projected_b - b * u * projected_imag
    difference = sp.trigsimp(
        sp.expand_complex(
            projected_slope - contact - c * u * projected_value
        )
    )
    if difference != 0:
        raise RuntimeError("normalized projected contact identity failed")

    rho1, rho2, nu1, nu2, theta = sp.symbols(
        "rho1 rho2 nu1 nu2 theta", real=True
    )
    amp1, amp2 = sp.symbols("amp1 amp2", positive=True, real=True)
    gamma1 = rho1 + sp.I * nu1
    gamma2 = rho2 + sp.I * nu2
    carrier1 = amp1 * (sp.cos(theta) + sp.I * sp.sin(theta))
    carrier2 = amp2
    prefix = carrier1 + carrier2
    prefix_x = gamma1 * carrier1 + gamma2 * carrier2
    current = sp.expand_complex(prefix_x * sp.conjugate(prefix))
    expected_radial = (
        rho1 * amp1**2
        + rho2 * amp2**2
        + amp1
        * amp2
        * (
            (rho1 + rho2) * sp.cos(theta)
            - (nu1 - nu2) * sp.sin(theta)
        )
    )
    expected_angular = (
        nu1 * amp1**2
        + nu2 * amp2**2
        + amp1
        * amp2
        * (
            (rho1 - rho2) * sp.sin(theta)
            + (nu1 + nu2) * sp.cos(theta)
        )
    )
    if sp.trigsimp(sp.re(current) - expected_radial) != 0:
        raise RuntimeError("normalized prefix radial current failed")
    if sp.trigsimp(sp.im(current) - expected_angular) != 0:
        raise RuntimeError("normalized prefix angular current failed")

    countermodel = recrossing_countermodel()
    if countermodel["zero_count"] != 6:
        raise RuntimeError("recrossing zero count failed")
    if countermodel["upward_count"] != 3:
        raise RuntimeError("recrossing upward count failed")
    if countermodel["signed_count"] != 0:
        raise RuntimeError("recrossing signed count failed")
    completion = dyadic_completion_scout()
    if completion["completed_upward_count"] != 1:
        raise RuntimeError("dyadic completion scout failed")

    return {
        "normalized_prefix": (
            "Because z_n=eta*q_n, e=eta*r_0, and g=eta*r_A, "
            "G_k=conj(eta)*F_k=sum_(n=1)^k q_n with G_0=0 and G_1=1. "
            "Thus Z_0=r_0+G_N=conj(eta)*W_0. The prefix polygon is the "
            "piecewise-linear path through G_0,G_1,...,G_N; it has no "
            "role in defining pi."
        ),
        "normalized_slope": (
            "Let V_tilde_N=r_A-s_*'*u_N*r_0 and "
            "B_N=V_tilde_N+s_*'*sum_(k=1)^(N-1)h_k*G_k. Then "
            "Z_A=s_*'*u_N*Z_0+B_N=conj(eta)*W_A. For "
            "Pi_eta(w)=Re(eta*w), mathsf_X=Pi_eta(Z_0), "
            "mathsf_A=Pi_eta(Z_A), and "
            "mathsf_A=Pi_eta(B_N)-b*u_N*Im(eta*Z_0)"
            "+c*u_N*mathsf_X. Hence at mathsf_X=0 the normalized "
            "formula is exactly mathcal_C_N=mathsf_A, including Z_0=0."
        ),
        "coefficient_current": (
            "On a fixed-N cell put epsilon_n=d_(n,x)/(1+d_n). Then "
            "gamma_n=q_(n,x)/q_n=-s_*'*log(n)+epsilon_n-epsilon_1, "
            "gamma_1=0, rho_n=Re(gamma_n), "
            "vartheta_n=Im(gamma_n), and "
            "G_(k,x)=sum_(n=1)^k gamma_n*q_n."
        ),
        "prefix_currents": (
            "Writing q_n=R_n*exp(i*theta_n), "
            "delta_nm=theta_n-theta_m, gamma_n=rho_n+i*vartheta_n, "
            "Re(G_(k,x)*conj(G_k))="
            "sum_n rho_n*R_n^2+sum_(n<m)R_n*R_m*"
            "[(rho_n+rho_m)*cos(delta_nm)"
            "-(vartheta_n-vartheta_m)*sin(delta_nm)]; "
            "Im(G_(k,x)*conj(G_k))="
            "sum_n vartheta_n*R_n^2+sum_(n<m)R_n*R_m*"
            "[(rho_n-rho_m)*sin(delta_nm)"
            "+(vartheta_n+vartheta_m)*cos(delta_nm)]."
        ),
    }


def recrossing_countermodel() -> dict:
    theta = sp.symbols("theta", real=True)
    core = sp.cos(theta) + sp.cos(3 * theta) / 2
    factored = sp.cos(theta) * (
        2 * sp.cos(theta) ** 2 - sp.Rational(1, 2)
    )
    if sp.trigsimp(sp.expand_trig(core) - factored) != 0:
        raise RuntimeError("recrossing factorization failed")
    derivative = sp.diff(core, theta)
    zeros = [
        sp.pi / 3,
        sp.pi / 2,
        2 * sp.pi / 3,
        4 * sp.pi / 3,
        3 * sp.pi / 2,
        5 * sp.pi / 3,
    ]
    derivative_values = [sp.simplify(derivative.subs(theta, z)) for z in zeros]
    upward = [z for z, value in zip(zeros, derivative_values) if value > 0]
    downward = [z for z, value in zip(zeros, derivative_values) if value < 0]
    if len(upward) + len(downward) != len(zeros):
        raise RuntimeError("recrossing model has a multiple zero")
    return {
        "coefficients": (
            "q_2=exp(-epsilon*theta)*exp(i*theta), "
            "q_8=(1/2)*exp(-epsilon*theta)*exp(3*i*theta), epsilon>0"
        ),
        "projection": (
            "X_epsilon(theta)=exp(-epsilon*theta)*"
            "[cos(theta)+(1/2)cos(3theta)]"
        ),
        "factorization": (
            "cos(theta)+(1/2)cos(3theta)="
            "cos(theta)*(2cos(theta)^2-1/2)"
        ),
        "zeros": "pi/3,pi/2,2pi/3,4pi/3,3pi/2,5pi/3",
        "upward_zeros": "pi/2,4pi/3,5pi/3",
        "downward_zeros": "pi/3,2pi/3,3pi/2",
        "zero_count": len(zeros),
        "upward_count": len(upward),
        "downward_count": len(downward),
        "signed_count": len(upward) - len(downward),
        "first_jet_winding": -len(upward),
        "endpoint_relation": (
            "Gamma_epsilon(2pi)=exp(-2pi*epsilon)*Gamma_epsilon(0) "
            "for Gamma_epsilon=X_epsilon+i*partial_theta X_epsilon"
        ),
        "interpretation": (
            "Both carriers move strictly inward, their angular currents "
            "are 1 and 3, the second amplitude is half the first, and the "
            "relative phase makes two turns. Nevertheless one phase cell "
            "has six simple crossings, three upward crossings, zero net "
            "signed crossing count, and first-jet winding -3. The simple "
            "zeros persist when the faster angular current is increased "
            "slightly, so a strict-more-than-two-turn version also exists."
        ),
    }


def dyadic_completion_scout() -> dict:
    c = sp.symbols("c", real=True)
    r = sp.sqrt(2) / 2
    completed = (
        sp.chebyshevt(1, c)
        + r * sp.chebyshevt(2, c)
        + r**2 * sp.chebyshevt(3, c)
    )
    polynomial = sp.expand(2 * sp.sqrt(2) * completed)
    expected = (
        4 * sp.sqrt(2) * c**3
        + 4 * c**2
        - sp.sqrt(2) * c
        - 2
    )
    if sp.simplify(polynomial - expected) != 0:
        raise RuntimeError("dyadic completion polynomial failed")
    discriminant = sp.discriminant(polynomial, c)
    left = sp.simplify(polynomial.subs(c, -1))
    right = sp.simplify(polynomial.subs(c, 1))
    if discriminant != -1696 or not left < 0 or not right > 0:
        raise RuntimeError("dyadic completion root certificate failed")
    return {
        "sparse_block": (
            "cos(theta)+(1/2)cos(3theta), from the leading "
            "n=2 and n=8 coefficients after scaling by q_2"
        ),
        "required_intermediate": (
            "The actual t=0 dyadic chain also contains n=4 with "
            "q_4/q_2=2^(-1/2) and phase 2theta."
        ),
        "completed_block": (
            "C_3(theta)=cos(theta)+2^(-1/2)cos(2theta)"
            "+(1/2)cos(3theta)"
        ),
        "chebyshev_polynomial": (
            "For c=cos(theta), 2sqrt(2)*C_3="
            "4sqrt(2)c^3+4c^2-sqrt(2)c-2=:P(c)."
        ),
        "discriminant": str(discriminant),
        "endpoint_signs": "P(-1)=2-3sqrt(2)<0<P(1)=2+3sqrt(2)",
        "real_root_count": 1,
        "unit_interval_root_count": 1,
        "completed_zero_count": 2,
        "completed_upward_count": 1,
        "interpretation": (
            "Restoring the arithmetically required intermediate power "
            "removes the sparse model's extra two upward crossings. This "
            "does not prove the full Xi chain, but it shows that exact "
            "multiplicative completion can supply rigidity that generic "
            "amplitude and phase ordering miss."
        ),
        "next_hypothesis": (
            "Test complete prime-power chains, then odd-part times dyadic "
            "chain decompositions, with the heat-quadratic and d_n "
            "perturbations retained. A theorem must also control sums of "
            "different chains and the endpoint; blockwise one-turn "
            "behavior alone is not promotable."
        ),
    }


def numeric_audit() -> dict[str, str | int]:
    l_value = float(L_MIN)
    x_min = 4 * math.pi * math.exp(l_value)
    t_max = 25 / l_value
    a_max = math.sqrt(math.exp(l_value) + t_max / 16)
    n_value = math.floor(a_max)
    cell_width = 4 * math.pi * (2 * n_value + 1)

    radial_floor = math.log(2) / (16 * l_value**2 * x_min)
    coefficient_error = 4 * D_X_CONSTANT / x_min**2
    radial_error_ratio = coefficient_error / radial_floor
    if not radial_error_ratio < 1:
        raise RuntimeError("q>=1 inward-current margin failed")

    adjacent_surplus = 2 * math.pi / n_value
    adjacent_error = coefficient_error * cell_width
    adjacent_error_ratio = adjacent_error / adjacent_surplus
    if not adjacent_error_ratio < 1:
        raise RuntimeError("adjacent full-cell rotation margin failed")
    adjacent_rotation_lower = (
        cell_width * (1 / (2 * n_value) - coefficient_error)
    )
    if not adjacent_rotation_lower > 4 * math.pi:
        raise RuntimeError("adjacent rotation did not exceed two turns")

    return {
        "x_min": format(x_min, ".17e"),
        "N_at_L50_tL25": n_value,
        "cell_width_at_boundary": format(cell_width, ".17e"),
        "q_ge_1_inward_floor": format(radial_floor, ".17e"),
        "coefficient_current_error": format(coefficient_error, ".17e"),
        "radial_error_to_floor_ratio": format(radial_error_ratio, ".17e"),
        "adjacent_rotation_lower": format(
            adjacent_rotation_lower, ".17e"
        ),
        "adjacent_surplus_over_4pi": format(adjacent_surplus, ".17e"),
        "adjacent_integrated_error": format(adjacent_error, ".17e"),
        "adjacent_error_to_surplus_ratio": format(
            adjacent_error_ratio, ".17e"
        ),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    numeric = numeric_audit()
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a); inward radial current uses "
            "q=2tL^2>=1 and all fixed-cell derivatives hold N constant"
        ),
        "pi_provenance": (
            "The pi in a^2=T_0/(2*pi)=x/(4*pi)+t/16 is inherited from "
            "xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2 and the "
            "Riemann-Siegel saddle. The cutoff width "
            "4*pi*(2N+1) is its algebraic cell difference. A phase cell "
            "uses the same constant because exp(i*(theta+2*pi))="
            "exp(i*theta). No circle or prefix polygon is selected."
        ),
        "normalized_prefix": symbolic["normalized_prefix"],
        "normalized_slope": symbolic["normalized_slope"],
        "coefficient_current": symbolic["coefficient_current"],
        "q_ge_1_inward": (
            "Since |epsilon_n|<8446/x^2, "
            "rho_n=-c*log(n)+Re(epsilon_n-epsilon_1). On q>=1, "
            "c=t*D_x/4>1/(8L^2*x). Therefore for n>=2, "
            "rho_n<-log(2)/(16L^2*x)<0; at L=50 the discarded "
            "coefficient-current error is less than 1.5e-14 of this "
            "certified inward floor."
        ),
        "ordered_rotation": (
            "For 1<=n<m<=N, "
            "vartheta_m-vartheta_n=(-b)*log(m/n)"
            "+Im(epsilon_m-epsilon_n)"
            ">=log(m/n)/2-16892/x^2>0. In a complete fixed-t cutoff "
            "cell, h_n=log((n+1)/n)>1/N and "
            "Delta x=4*pi*(2N+1), so every adjacent relative phase "
            "advances by more than 4*pi. The worst L=50 integrated error "
            "is below 8.3e-20 of the surplus over 4*pi."
        ),
        "prefix_currents": symbolic["prefix_currents"],
        "current_sign_guard": (
            "The diagonal angular terms are positive, but both pair "
            "kernels contain unrestricted sine and cosine factors. "
            "Strict inward motion and ordered angular currents therefore "
            "do not assign a sign to a prefix radial or angular current."
        ),
        "recrossing_countermodel": recrossing_countermodel(),
        "dyadic_completion_scout": dyadic_completion_scout(),
        "signed_flux_identity": (
            "For a regular real scalar X on [A,B], "
            "N_up-N_down=[sgn(X(B))-sgn(X(A))]/2. This endpoint flux "
            "controls only the net signed count. In the exact recrossing "
            "model it is zero while N_up=3."
        ),
        "first_jet_flux": (
            "For ell>0 define Gamma_ell=mathsf_X+i*mathsf_A/ell. "
            "Where Gamma_ell!=0, "
            "partial_x arg(Gamma_ell)="
            "{mathsf_X*partial_x(mathsf_A/ell)"
            "-(mathsf_A/ell)*partial_x mathsf_X}/"
            "{mathsf_X^2+(mathsf_A/ell)^2}. "
            "At mathsf_X=0, mathsf_A>0, and "
            "partial_x mathsf_X>0, the positive-imaginary-ray "
            "intersection has sign -1. Thus the first-jet argument flux, "
            "with the prescribed connector and half-open conventions, "
            "is the correct object for the upward count."
        ),
        "normalized_flux_derivative": (
            "Let eta_x=i*omega_eta*eta, "
            "Z_(0,x)=r_(0,x)+sum_n gamma_n*q_n, "
            "V_tilde_(N,x)=r_(A,x)-s_*''*u_N*r_0"
            "-s_*'*u_(N,x)*r_0-s_*'*u_N*r_(0,x), "
            "u_(N,x)=1/(8*pi*a^2)=1/(4*T_0), and "
            "B_(N,x)=V_tilde_(N,x)+s_*''*sum_k h_k*G_k"
            "+s_*'*sum_k h_k*sum_(n<=k)gamma_n*q_n. Then "
            "Z_(A,x)=(s_*''*u_N+s_*'*u_(N,x))*Z_0"
            "+s_*'*u_N*Z_(0,x)+B_(N,x), "
            "mathsf_X_x=Pi_eta(Z_(0,x))"
            "-omega_eta*Im(eta*Z_0), and "
            "mathsf_A_x=Pi_eta(Z_(A,x))"
            "-omega_eta*Im(eta*Z_A). This supplies an exact O(N) "
            "integrand for the first-jet phase flux."
        ),
        "live_q_ge_1_target": (
            "Retain the pointwise all-fiber lower bound "
            "|mathcal_C_N|>A_L+epsilon_term on "
            "|mathsf_X|<=delta_L. For the integer, do not count turns of "
            "an individual carrier. Insert the normalized O(N) flux into "
            "the complete successor boundary composition and prove "
            "0<=kappa_j=(2*pi)^(-1)*Delta_arg(Gamma_j)<1 after top, "
            "bottom, vertical, chart-join, and finite-shoulder terms are "
            "combined. The lower inequality is topological; the strict "
            "upper inequality is the open arithmetic theorem."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 layer still requires a separate "
            "multiplicity-compatible parabolic/Hermite first-jet flux "
            "chart. No uniform positive slope at t=0 is introduced."
        ),
        "numeric_audit": numeric,
    }


def build_rows(exact: dict) -> list[FluxRow]:
    return [
        FluxRow(
            "nfnpfr_00_pi_provenance",
            "definition_provenance",
            "certified",
            "The saddle, cell-width, and phase-period appearances of pi have an explicit common provenance.",
            exact["pi_provenance"],
            "The prefix polygon neither defines nor approximates pi.",
        ),
        FluxRow(
            "nfnpfr_01_normalized_prefix",
            "exact_normal_form",
            "ready_to_apply",
            "Removing the first coefficient phase turns every carrier prefix into a literal relative-coefficient partial sum.",
            exact["normalized_prefix"],
            "The nonzero unit anchor is retained as the external real projector.",
        ),
        FluxRow(
            "nfnpfr_02_normalized_slope",
            "exact_reduction",
            "ready_to_apply",
            "The endpoint-complete Abel scalar has an exact normalized-prefix form and still includes the zero fiber.",
            exact["normalized_slope"],
            "No argument division or slope sign is assumed.",
        ),
        FluxRow(
            "nfnpfr_03_coefficient_current",
            "exact_differential_identity",
            "ready_to_apply",
            "Every normalized coefficient and every prefix derivative has an exact branch-free logarithmic current.",
            exact["coefficient_current"],
            "All derivatives are inside one fixed cutoff chart.",
        ),
        FluxRow(
            "nfnpfr_04_inward_current",
            "asymptotic_certificate",
            "certified",
            "On q>=1 every nonfirst normalized Xi coefficient moves strictly inward as x increases.",
            exact["q_ge_1_inward"],
            "This is a coefficientwise radial statement, not a prefix-current sign.",
            exact["numeric_audit"],
        ),
        FluxRow(
            "nfnpfr_05_ordered_rotation",
            "asymptotic_certificate",
            "certified",
            "All normalized Xi angular currents are strictly ordered and every adjacent pair makes over two relative turns per complete cutoff cell.",
            exact["ordered_rotation"],
            "A complete cell must lie inside the stated L>=50 domain.",
            exact["numeric_audit"],
        ),
        FluxRow(
            "nfnpfr_06_prefix_currents",
            "exact_kernel",
            "ready_to_apply",
            "The normalized prefix radial and angular currents have explicit diagonal and pairwise kernels.",
            exact["prefix_currents"],
            exact["current_sign_guard"],
        ),
        FluxRow(
            "nfnpfr_07_recrossing_countermodel",
            "countermodel",
            "guard_validated",
            "Inward ordered rotating carriers can produce three upward crossings and three first-jet turns in one phase cell.",
            json.dumps(exact["recrossing_countermodel"], sort_keys=True),
            "This is an exact generic logarithmic-frequency route guard, not an Xi counterexample.",
            exact["recrossing_countermodel"],
        ),
        FluxRow(
            "nfnpfr_07a_dyadic_completion",
            "diagnostic",
            "diagnostic_validated",
            "Restoring the missing leading dyadic carrier exactly repairs the smallest sparse recrossing obstruction.",
            json.dumps(exact["dyadic_completion_scout"], sort_keys=True),
            "This revives multiplicative block completion as a route; it does not control the sum of blocks or the Xi endpoint.",
            exact["dyadic_completion_scout"],
        ),
        FluxRow(
            "nfnpfr_08_signed_flux_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Endpoint signed flux controls upward minus downward crossings, not the one-sided upward count.",
            exact["signed_flux_identity"],
            "A no-recrossing theorem or a first-jet argument bound is still required.",
        ),
        FluxRow(
            "nfnpfr_09_first_jet_flux",
            "exact_topological_reduction",
            "ready_to_apply",
            "The scaled first-jet argument, not an individual carrier phase, is the correct crossing flux.",
            exact["first_jet_flux"],
            "Endpoint, connector, and half-open conventions remain part of the count.",
        ),
        FluxRow(
            "nfnpfr_10_flux_integrand",
            "exact_differential_identity",
            "ready_to_apply",
            "The complete first-jet phase-flux density has an endpoint-complete O(N) normalized-prefix representation.",
            exact["normalized_flux_derivative"],
            "No sign or integral bound for this density is asserted.",
        ),
        FluxRow(
            "nfnpfr_11_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 integer problem is a composed successor first-jet flux bound, not a carrier-cell crossing bound.",
            exact["live_q_ge_1_target"],
            "The strict upper flux bound and pointwise Xi prefix lower bound remain open.",
        ),
        FluxRow(
            "nfnpfr_12_q_lt_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1 multiplicity-compatible first-jet flux chart remains separate.",
            exact["q_lt_1_target"],
            "No simplicity assumption is promoted to t=0.",
        ),
        FluxRow(
            "nfnpfr_13_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "Normalized spiralling and exact phase flux are kept separate from the missing Xi inequalities and RH.",
            exact["current_sign_guard"],
            (
                "No prefix lower bound, successor flux upper bound, q<1 "
                "theorem, finite shoulder closure, contact exclusion, "
                "Lambda<=0, PF-infinity, RH, or Clay-prize conclusion."
            ),
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-27",
        "status": (
            "exact normalized-prefix and first-jet flux reduction, "
            "certified Xi coefficient spiralling, and a recrossing "
            "nonpromotion guard, with an exact dyadic-completion diagnostic; "
            "the composed successor bound remains open"
        ),
        "proof_boundary": (
            "This artifact proves the normalized-prefix identities, exact "
            "coefficient and prefix currents, q>=1 inward coefficient "
            "motion, strict ordered relative rotation, the O(N) first-jet "
            "flux integrand, an exact recrossing countermodel, and the "
            "smallest exact dyadic-completion repair diagnostic. It does "
            "not prove the Xi prefix lower bound, the strict composed "
            "successor flux bound, the q<1 multiplicity-compatible theorem, "
            "finite shoulder closure, contact exclusion, Lambda<=0, "
            "PF-infinity, RH, or a Clay-prize conclusion."
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "constants": {
            "L_min": L_MIN,
            "d_x": D_X_CONSTANT,
            "d": D_CONSTANT,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    numeric = exact["numeric_audit"]
    countermodel = exact["recrossing_countermodel"]
    completion = exact["dyadic_completion_scout"]
    return "\n".join(
        [
            "# Newman Normalized-Prefix Phase-Flux Reduction",
            "",
            "Date: 2026-07-27",
            "",
            "Status: exact normalized-prefix and first-jet flux formulas,",
            "with a certified Xi spiralling law and a decisive recrossing",
            "guard. The successor flux theorem remains open. This is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Pi Provenance",
            "",
            "```text",
            exact["pi_provenance"],
            "```",
            "",
            "The same universal constant is obtained from every Euclidean",
            "circle, but no circle is selected in this calculation. The",
            "prefix polygon is only a complex partial-sum path.",
            "",
            "## Normalized Polygon",
            "",
            "```text",
            exact["normalized_prefix"],
            exact["normalized_slope"],
            "```",
            "",
            "The normalization removes a nonzero unit phase; it does not",
            "discard the absolute projector or the endpoint.",
            "",
            "## Carrier Currents",
            "",
            "```text",
            exact["coefficient_current"],
            exact["q_ge_1_inward"],
            exact["ordered_rotation"],
            "```",
            "",
            "At the L=50 boundary the inward error-to-floor ratio is",
            f"`{numeric['radial_error_to_floor_ratio']}`, and the",
            "adjacent integrated error-to-4pi-surplus ratio is",
            f"`{numeric['adjacent_error_to_surplus_ratio']}`.",
            "",
            "The prefix currents themselves are",
            "",
            "```text",
            exact["prefix_currents"],
            exact["current_sign_guard"],
            "```",
            "",
            "## Recrossing Guard",
            "",
            "```text",
            json.dumps(countermodel, indent=2, sort_keys=True),
            "```",
            "",
            "This exact model has the same leading logarithmic frequency",
            "ratio `log(8)/log(2)=3` and amplitude ratio",
            "`8^(-1/2)/2^(-1/2)=1/2`. It is a route guard, not an Xi",
            "counterexample: the actual endpoint and full coefficient chain",
            "are deliberately absent.",
            "",
            "The first arithmetic completion test gives",
            "",
            "```text",
            json.dumps(completion, indent=2, sort_keys=True),
            "```",
            "",
            "Thus the generic shortcut is false, but the actual",
            "multiplicative chain repairs its smallest sparse obstruction.",
            "The full prime-power, cross-chain, heat, and endpoint theorem",
            "is still open.",
            "",
            "```text",
            exact["signed_flux_identity"],
            "```",
            "",
            "So a signed endpoint flux cannot replace the one-sided count.",
            "",
            "## First-Jet Flux",
            "",
            "```text",
            exact["first_jet_flux"],
            exact["normalized_flux_derivative"],
            "```",
            "",
            "This is the surviving exact O(N) phase-flux representation.",
            "",
            "## Live Theorem",
            "",
            "```text",
            exact["live_q_ge_1_target"],
            exact["q_lt_1_target"],
            "```",
            "",
            "No Xi prefix lower bound, successor flux upper bound, q<1",
            "chart, finite shoulder closure, contact exclusion,",
            "`Lambda<=0`, PF-infinity, RH, or Clay-prize conclusion is",
            "claimed.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman normalized-prefix phase-flux reduction: "
        "15 rows, 2 certified spiral bounds, 1 exact recrossing guard, "
        "1 dyadic completion diagnostic, "
        "1 O(N) first-jet flux, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
