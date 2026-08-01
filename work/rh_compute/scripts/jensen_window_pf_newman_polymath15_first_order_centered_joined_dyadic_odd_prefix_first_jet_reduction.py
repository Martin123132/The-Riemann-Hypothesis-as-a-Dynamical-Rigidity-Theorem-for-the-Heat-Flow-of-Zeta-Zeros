#!/usr/bin/env python3
"""Build the joined dyadic odd-prefix first-jet reduction."""

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
    "joined_dyadic_odd_prefix_first_jet_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
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
    "adjacent_chart": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_chart_stability_certificate.json"
    ),
}

L_MIN = 50
UNCORRECTED_DYADIC_RATIO = "2^(-49/100)"
CORRECTED_DYADIC_RATIO = "2^(-12/25)"


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
        "prime_power_guard": (
            "p^(-12/25)",
            "O_M^(p)(s_Z)",
            '"joined_winding": 3',
        ),
        "normalized_prefix": (
            "Gamma_ell=mathsf_X+i*mathsf_A/ell",
            "u_(N,x)=1/(8*pi*a^2)=1/(4*T_0)",
            "Z_(0,x)=r_(0,x)+sum_n gamma_n*q_n",
        ),
        "adjacent_recurrence": (
            "Q_N=f_(N+1)+kappa_N*J_a",
            "Delta A_a=Re[-s_*'*ell*f_(N+1)",
        ),
        "adjacent_chart": (
            "|Delta X|<1100*exp(-5L/4)",
            "|A_(a,N+1)-A_(a,N)|<5000*exp(-7L/4)",
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


def valuation_audit(n_max: int = 96) -> dict[str, int]:
    reconstructed: list[int] = []
    layer_counts: dict[int, int] = {}
    for k in range(int(math.log2(n_max)) + 1):
        bound = n_max // (2**k)
        for odd in range(1, bound + 1, 2):
            reconstructed.append((2**k) * odd)
            layer_counts[k] = layer_counts.get(k, 0) + 1
    if sorted(reconstructed) != list(range(1, n_max + 1)):
        raise RuntimeError("dyadic odd-prefix valuation audit failed")

    for n in range(1, n_max + 1):
        k = 0
        odd = n
        while odd % 2 == 0:
            odd //= 2
            k += 1
        matches = [
            (layer, base)
            for layer in range(int(math.log2(n)) + 1)
            for base in range(1, n // (2**layer) + 1, 2)
            if (2**layer) * base == n
        ]
        if matches != [(k, odd)]:
            raise RuntimeError(f"nonunique dyadic valuation at n={n}")
    return {
        "n_max": n_max,
        "reconstructed_terms": len(reconstructed),
        "layers": len(layer_counts),
    }


def two_fibre_winding_guard() -> dict[str, str | int]:
    z = sp.symbols("z")
    r = sp.sqrt(2) / 2
    first = z * (1 + r * z + r**2 * z**2)
    second = -sp.Rational(4, 5) * z * (1 + r * z)
    mean_inner = sp.expand(10 * (first + second) / (2 * z))
    expected = 1 + r * z + 5 * r**2 * z**2
    if sp.expand(mean_inner - expected) != 0:
        raise RuntimeError("two-fibre mean factorization failed")
    if sp.discriminant(mean_inner, z) != -sp.Rational(19, 2):
        raise RuntimeError("two-fibre discriminant failed")
    coefficients = sp.Poly(mean_inner, z).all_coeffs()
    root_product = sp.simplify(coefficients[-1] / coefficients[0])
    if root_product != sp.Rational(2, 5):
        raise RuntimeError("two-fibre root product failed")
    return {
        "first_fibre": "A(z)=z*(1+r*z+r^2*z^2)",
        "second_fibre": "B(z)=-(4/5)*z*(1+r*z)",
        "fibre_windings": "wind(A)=wind(B)=1 on |z|=1",
        "equal_weight_mean": "(A+B)/2",
        "mean_inner_discriminant": "-19/2",
        "mean_inner_root_modulus": "sqrt(2/5)<1",
        "mean_winding": 3,
        "interpretation": (
            "Averaging boundary-nonzero winding-one fibres can produce "
            "winding three. A per-Gaussian-fibre theorem cannot be "
            "averaged without additional Gaussian/Xi structure."
        ),
    }


def second_moment_guard() -> dict[str, str]:
    h = sp.symbols("h", nonzero=True, real=True)
    weights = (sp.Integer(1), sp.Integer(-2), sp.Integer(1))
    logs = (sp.Integer(0), h, 2 * h)
    moments = [
        sp.expand(sum(weight * log_value**order for weight, log_value in zip(
            weights, logs, strict=True
        )))
        for order in range(3)
    ]
    if moments != [0, 0, 2 * h**2]:
        raise RuntimeError("second-moment guard failed")
    return {
        "support": "logarithmic nodes 0,h,2h",
        "weights": "1,-2,1",
        "moments": "H_0=0, H_1=0, H_2=2h^2",
        "interpretation": (
            "For unrestricted joined coefficients, value and first "
            "logarithmic moment do not determine the second moment "
            "entering the first-jet derivative. Xi-specific structure "
            "must control H_2 rather than delete it."
        ),
    }


def symbolic_audit() -> dict[str, str | dict]:
    time, h, k, log_m, spectral = sp.symbols(
        "time h k log_m spectral"
    )
    full_log = k * h + log_m
    full_exponent = time * full_log**2 / 4 - spectral * full_log
    shifted_exponent = (
        time * (k * h) ** 2 / 4
        - spectral * k * h
        + time * log_m**2 / 4
        - (spectral - time * k * h / 2) * log_m
    )
    if sp.expand(full_exponent - shifted_exponent) != 0:
        raise RuntimeError("deterministic heat-shift factorization failed")

    z, xi = sp.symbols("z xi")
    integer_k = sp.symbols("integer_k", integer=True, nonnegative=True)
    monomial = z**integer_k * sp.exp(sp.I * xi * log_m)

    def total_log_operator(value: sp.Expr) -> sp.Expr:
        return sp.expand(
            h * z * sp.diff(value, z) - sp.I * sp.diff(value, xi)
        )

    first = sp.simplify(total_log_operator(monomial) / monomial)
    second = sp.simplify(
        total_log_operator(total_log_operator(monomial)) / monomial
    )
    expected_first = integer_k * h + log_m
    if sp.simplify(first - expected_first) != 0:
        raise RuntimeError("first phase-lift moment operator failed")
    if sp.simplify(second - expected_first**2) != 0:
        raise RuntimeError("second phase-lift moment operator failed")

    s_1, s_2, log_a, u_x = sp.symbols("s_1 s_2 log_a u_x")
    h_0, h_1, h_2, d_0, d_1 = sp.symbols(
        "H_0 H_1 H_2 D_0 D_1"
    )
    h_0_x = -s_1 * h_1 + d_0
    h_1_x = -s_1 * h_2 + d_1
    centered = h_1 - log_a * h_0
    centered_x = h_1_x - u_x * h_0 - log_a * h_0_x
    expected_centered_x = (
        -s_1 * (h_2 - log_a * h_1)
        + d_1
        - log_a * d_0
        - u_x * h_0
    )
    if sp.expand(centered_x - expected_centered_x) != 0:
        raise RuntimeError("centered moment derivative failed")

    mu, d_n, d_anchor = sp.symbols("mu d_n d_anchor")
    correction = (1 + mu * d_n) / (1 + mu * d_anchor)
    correction_mu = sp.factor(sp.diff(correction, mu))
    expected_mu = (d_n - d_anchor) / (1 + mu * d_anchor) ** 2
    if sp.simplify(correction_mu - expected_mu) != 0:
        raise RuntimeError("correction homotopy derivative failed")

    return {
        "heat_shift": (
            "For A_t(n;s)=exp[t*log(n)^2/4-s*log(n)], h=log(2), "
            "A_t(2^k*m;s)=A_t(2^k;s)"
            "*A_t(m;s-t*k*h/2)."
        ),
        "moment_operator": (
            "For D_h=h*z*partial_z-i*partial_xi, "
            "D_h[z^k*exp(i*xi*log(m))]="
            "log(2^k*m)*z^k*exp(i*xi*log(m)); "
            "D_h^r supplies the r-th total logarithmic moment."
        ),
        "centered_derivative": (
            "If H_(r,x)=-s_*'*H_(r+1)+D_r and "
            "P_1=H_1-log(a)*H_0, then "
            "P_(1,x)=-s_*'*(H_2-log(a)*H_1)"
            "+D_1-log(a)*D_0-u_x*H_0."
        ),
        "correction_homotopy": (
            "partial_mu[(1+mu*d_n)/(1+mu*d_1)]="
            "(d_n-d_1)/(1+mu*d_1)^2."
        ),
        "valuation_audit": valuation_audit(),
        "two_fibre_guard": two_fibre_winding_guard(),
        "second_moment_guard": second_moment_guard(),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    root_radius = 2 ** (49 / 100)
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a), and x differentiation "
            "is performed inside one fixed-N cutoff chart"
        ),
        "pi_provenance": (
            "The pi in a^2=x/(4*pi)+t/16 is inherited from "
            "xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2 and "
            "the Riemann-Siegel saddle. The cutoff-cell difference is "
            "4*pi*(2N+1), while 2*pi is the period of exp(i*theta). "
            "The dyadic base 2, odd prefixes, phase lift, and "
            "phase-synchronized hypothetical do not define pi."
        ),
        "coefficient_chain": (
            "Put A_t(n;s)=exp[t*log(n)^2/4-s*log(n)], "
            "c_n=(1+d_n)/(1+d_1), epsilon_n=d_(n,x)/(1+d_n), "
            "delta_n=epsilon_n-epsilon_1, and "
            "q_n=A_t(n;s_*)c_n. Then "
            "q_(n,x)/q_n=-s_*'*log(n)+delta_n."
        ),
        "heat_shift": symbolic["heat_shift"],
        "joined_odd_prefix": (
            "Let h=log(2), K=floor(log_2 N), M_k=floor(N/2^k), "
            "and for j>=0 define "
            "O_(k,j)=sum_(m<=M_k,m odd)(log m)^j"
            "*A_t(m;s_*-t*k*h/2)c_(2^k*m). "
            "Then, for H_r=sum_(n<=N)(log n)^r q_n, "
            "H_r=sum_(k=0)^K A_t(2^k;s_*)"
            "*sum_(j=0)^r binom(r,j)(k*h)^(r-j)O_(k,j). "
            "This is an exact deterministic heat-shifted odd-prefix "
            "identity, not a conditional or averaged approximation."
        ),
        "joined_correction_current": (
            "Define E_(k,j) by replacing c_(2^k*m) in O_(k,j) "
            "with delta_(2^k*m)c_(2^k*m), and put "
            "D_r=sum_k A_t(2^k;s_*)"
            "*sum_(j=0)^r binom(r,j)(k*h)^(r-j)E_(k,j). "
            "Then D_r=sum_(n<=N)(log n)^r delta_n q_n exactly."
        ),
        "five_current_ladder": (
            "The fixed-chart logarithmic moments obey "
            "H_(r,x)=-s_*'*H_(r+1)+D_r. Hence the endpoint-complete "
            "first jet closes on H_0,H_1,H_2,D_0,D_1 together with "
            "r_0,r_A and their x derivatives; no individual carrier "
            "phase or Gaussian fibre has to be discarded."
        ),
        "phase_lift": (
            "Write sigma_*=Re(s_*), omega_*=-Im(s_*), "
            "a_(k,m)=exp[t*log(2^k*m)^2/4"
            "-sigma_*log(2^k*m)]>0, and define "
            "P_N(z,xi,mu)=sum_(k=0)^K z^k"
            "*sum_(m<=M_k,m odd)a_(k,m)exp(i*xi*log m)"
            "*(1+mu*d_(2^k*m))/(1+mu*d_1). "
            "For z_*=exp(i*omega_*h), "
            "P_N(z_*,omega_*,1)=H_0. With "
            "D_h=h*z*partial_z-i*partial_xi, "
            "D_h^r P_N(z_*,omega_*,1)=H_r."
        ),
        "phase_lift_derivatives": (
            symbolic["moment_operator"]
            + " Also partial_mu P_N is obtained termwise from "
            + symbolic["correction_homotopy"]
            + " The physical fixed-chart derivative is "
            "H_(0,x)=-s_*'*D_h P_N+D_0 at the actual point."
        ),
        "centered_endpoint_jet": (
            "Let A=log(a), u_x=A_x=1/(8*pi*a^2), "
            "P_1=H_1-AH_0, Z_0=r_0+H_0, and "
            "Z_A=r_A-s_*'P_1. Then "
            "Z_(0,x)=r_(0,x)-s_*'H_1+D_0, "
            "P_(1,x)=-s_*'(H_2-AH_1)+D_1-AD_0-u_xH_0, and "
            "Z_(A,x)=r_(A,x)-s_*''P_1-s_*'P_(1,x). "
            "These formulas include the endpoint before any absolute "
            "value or winding estimate is taken."
        ),
        "projected_flux": (
            "For eta_x=i*omega_eta*eta and Pi_eta(w)=Re(eta*w), "
            "mathsf_X=Pi_eta(Z_0), mathsf_A=Pi_eta(Z_A), "
            "mathsf_X_x=Pi_eta(Z_(0,x))-omega_eta*Im(eta*Z_0), "
            "and mathsf_A_x=Pi_eta(Z_(A,x))"
            "-omega_eta*Im(eta*Z_A). For ell>0, "
            "Gamma_ell=mathsf_X+i*mathsf_A/ell and "
            "partial_x arg Gamma_ell="
            "[mathsf_X*partial_x(mathsf_A/ell)"
            "-(mathsf_A/ell)*mathsf_X_x]"
            "/[mathsf_X^2+(mathsf_A/ell)^2]."
        ),
        "cutoff_jump": (
            "At N->N+1 write n=N+1=2^v*m with m odd. Exactly one "
            "odd-prefix layer, v, gains one term, so "
            "Delta H_r=(log n)^r q_n and "
            "Delta D_r=(log n)^r delta_n q_n. With "
            "j_0=kappa_N*J_a/f_1 and "
            "j_A=kappa_N*(J_(a,x)+mu_a*J_a)/f_1, "
            "Delta Z_0=q_n+j_0=Q_N/f_1 and "
            "Delta Z_A=-s_*'*log(n/a)q_n+j_A. "
            "Differentiating gives "
            "Delta Z_(0,x)=gamma_nq_n+(j_0)_x and "
            "Delta Z_(A,x)=(j_A)_x-s_*''log(n/a)q_n"
            "-s_*'[-u_xq_n+log(n/a)gamma_nq_n]."
        ),
        "adjacent_control": (
            "The imported exact adjacent recurrence is therefore the "
            "cutoff boundary law of the joined representation, not a "
            "separate approximation. Its real projections retain "
            "|Delta X|<1100exp(-5L/4) and "
            "|Delta A_a|<5000exp(-7L/4); no complex jump bound is "
            "invented."
        ),
        "phase_frozen_hypothetical": (
            "Set xi=0 and mu=0 while keeping z free. Then "
            "P_fr(z)=sum_(k=0)^K B_kz^k with "
            "B_k=sum_(m<=M_k,m odd)a_(k,m)>0. "
            "For every m retained in layer k+1, the certified "
            "correction-free dyadic ratio is below 2^(-49/100), "
            "and the odd index set shrinks. Thus "
            "0<B_(k+1)<2^(-49/100)B_k. Enestrom-Kakeya puts every "
            "zero outside |z|=2^(49/100)>1. Therefore P_fr is "
            "unit-disk zero-free with winding zero, and zP_fr has "
            "winding one. This is a proved theorem for the stipulated "
            "phase-synchronized hypothetical, not for the actual Xi "
            "odd phases."
        ),
        "phase_frozen_margin": (
            "If rho_1,...,rho_K are the roots of P_fr, then "
            "|rho_j|>R_0=2^(49/100). Hence on |z|=1, "
            "|P_fr(z)|=B_K*product_j|z-rho_j|"
            ">B_K*(R_0-1)^K. This explicit margin is positive but may "
            "be too small for the full odd-phase homotopy."
        ),
        "homotopy_target": (
            "The exact deformation from the proved base to the actual "
            "finite prefix is (xi,mu):(0,0)->(omega_*,1) in P_N, "
            "with z=z_*. A viable new theorem must exclude "
            "P_N(z_*,xi,mu)=0 along a specified path or control every "
            "boundary crossing using the D_h first and second moments. "
            "The endpoint r_0 and centered jet r_A-s_*'P_1 must then be "
            "inserted before the successor argument is counted. Crude "
            "L1 phase perturbation is not assumed to be smaller than "
            "the phase-frozen margin."
        ),
        "two_fibre_guard": symbolic["two_fibre_guard"],
        "second_moment_guard": symbolic["second_moment_guard"],
        "route_decision": (
            "Retain the dyadic heat-shift, five-current closure, and "
            "phase-lift homotopy. Retire any route that proves winding "
            "one separately for Gaussian fibres and then averages, or "
            "that estimates the first-jet flux from H_0,H_1 while "
            "deleting H_2. The next arithmetic target is a linked-curve "
            "boundary-zero or signed-current theorem for P_N with the "
            "actual odd Mellin phases, followed by the exact endpoint "
            "augmentation and adjacent recurrence."
        ),
        "live_q_ge_1_target": (
            "On q=2tL^2>=1, prove an Xi-specific transport theorem from "
            "the phase-frozen base to (z_*,omega_*,1), strong enough "
            "to yield |mathsf_X|<=delta_L => "
            "|mathcal_C_N|>A_L+epsilon_term and the complete composed "
            "successor bound 0<=kappa_j<1. Candidate mechanisms include "
            "a Hermite-Biehler/de Branges exponential-polynomial "
            "criterion, a signed two-parameter Jacobian, or a "
            "Rouche/argument-principle estimate that uses the linked "
            "dyadic and odd phases rather than independent block "
            "windings."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 layer remains a separate "
            "multiplicity-compatible parabolic/Hermite first-jet chart. "
            "The joined dyadic identities remain valid there, but no "
            "uniform endpoint simplicity or positive slope is inferred."
        ),
        "proof_boundary": (
            "The heat-shift factorization, joined odd-prefix moments, "
            "five-current derivative closure, phase lift, cutoff jump, "
            "endpoint first jet, phase-frozen zero-free theorem and "
            "margin, and two generic route guards are exact or imported "
            "from certified prerequisites. The actual odd-phase/correction "
            "homotopy, joined Xi lower bound, strict successor flux "
            "upper bound, q<1 chart, finite shoulders and connectors, "
            "contact exclusion, Lambda<=0, PF-infinity, RH, and a "
            "Clay-prize conclusion remain open."
        ),
        "diagnostics": {
            "valuation": symbolic["valuation_audit"],
            "phase_frozen_root_radius_lower": format(root_radius, ".17e"),
            "phase_frozen_root_gap_lower": format(root_radius - 1, ".17e"),
            "uncorrected_layer_ratio": UNCORRECTED_DYADIC_RATIO,
            "corrected_layer_ratio": CORRECTED_DYADIC_RATIO,
        },
    }


def build_rows(exact: dict) -> list[ReductionRow]:
    return [
        ReductionRow(
            "jdopfj_00_pi_provenance",
            "definition_provenance",
            "certified",
            "Every pi in the joined reduction has an explicit completed-zeta, saddle, or phase-period source.",
            exact["pi_provenance"],
            "The dyadic and odd-prefix constructions do not define pi.",
        ),
        ReductionRow(
            "jdopfj_01_coefficient_chain",
            "exact_identity",
            "ready_to_apply",
            "The corrected normalized coefficient and its x-current are fixed before reindexing.",
            exact["coefficient_chain"],
            "No d_n correction is dropped.",
        ),
        ReductionRow(
            "jdopfj_02_heat_shift",
            "exact_factorization",
            "ready_to_apply",
            "Each dyadic layer is a deterministic heat-shifted odd prefix.",
            exact["heat_shift"],
            "This is not a probabilistic approximation.",
        ),
        ReductionRow(
            "jdopfj_03_joined_odd_prefix",
            "exact_reindexing",
            "ready_to_apply",
            "All coefficient moments rejoin through nested odd-prefix layers.",
            exact["joined_odd_prefix"],
            "The odd prefixes remain complex at the physical phase.",
            exact["diagnostics"],
        ),
        ReductionRow(
            "jdopfj_04_correction_current",
            "exact_reindexing",
            "ready_to_apply",
            "The complete d_n logarithmic derivative has the same joined odd-prefix form.",
            exact["joined_correction_current"],
            "Correction currents are retained in D_0 and D_1.",
        ),
        ReductionRow(
            "jdopfj_05_phase_lift",
            "exact_phase_lift",
            "ready_to_apply",
            "The physical prefix is one linked point of a polynomial-Mellin phase lift.",
            exact["phase_lift"] + " " + exact["phase_lift_derivatives"],
            "The xi variable is a Mellin phase, not an independent Xi zero coordinate.",
        ),
        ReductionRow(
            "jdopfj_06_five_current_closure",
            "exact_differential_reduction",
            "ready_to_apply",
            "The joined value and first jet close on five finite currents plus endpoint data.",
            exact["five_current_ladder"],
            "H_2, D_0, and D_1 may not be deleted.",
        ),
        ReductionRow(
            "jdopfj_07_endpoint_jet",
            "exact_endpoint_reduction",
            "ready_to_apply",
            "The centered endpoint-complete first jet is explicit in the five currents.",
            exact["centered_endpoint_jet"] + " " + exact["projected_flux"],
            "The formula is division-free on Z_0=0.",
        ),
        ReductionRow(
            "jdopfj_08_cutoff_jump",
            "exact_chart_transition",
            "ready_to_apply",
            "A cutoff change enters exactly one valuation layer and the recurrent endpoint simultaneously.",
            exact["cutoff_jump"] + " " + exact["adjacent_control"],
            "No complex adjacent-jump estimate is asserted.",
        ),
        ReductionRow(
            "jdopfj_09_phase_frozen_theorem",
            "hypothetical_exact_theorem",
            "certified",
            "The stipulated phase-synchronized joined model is unit-disk zero-free.",
            exact["phase_frozen_hypothetical"],
            "This hypothetical is a homotopy base, not the physical odd-phase theorem.",
        ),
        ReductionRow(
            "jdopfj_10_phase_frozen_margin",
            "exact_quantitative_corollary",
            "certified",
            "The phase-frozen theorem has an explicit boundary-modulus margin.",
            exact["phase_frozen_margin"],
            "No claim is made that crude odd-phase perturbation fits inside this margin.",
            exact["diagnostics"],
        ),
        ReductionRow(
            "jdopfj_11_homotopy_target",
            "new_theorem_target",
            "not_ready_to_apply",
            "The physical arithmetic problem is a linked odd-phase and correction homotopy with endpoint augmentation.",
            exact["homotopy_target"],
            "Boundary-zero exclusion along this homotopy is open.",
        ),
        ReductionRow(
            "jdopfj_12_averaging_guard",
            "countermodel",
            "guard_validated",
            "Winding-one fibres cannot be averaged into a winding-one conclusion.",
            json.dumps(exact["two_fibre_guard"], sort_keys=True),
            "This is a generic exact guard, not Xi phase data.",
            exact["two_fibre_guard"],
        ),
        ReductionRow(
            "jdopfj_13_second_moment_guard",
            "countermodel",
            "guard_validated",
            "Value and first logarithmic moment do not determine the second moment.",
            json.dumps(exact["second_moment_guard"], sort_keys=True),
            "Xi-specific coefficient constraints may add information, but it must be proved.",
            exact["second_moment_guard"],
        ),
        ReductionRow(
            "jdopfj_14_route_decision",
            "route_decision",
            "guard_validated",
            "The next route is a linked phase-lift transport theorem, not independent block winding.",
            exact["route_decision"],
            "The endpoint and second moment remain joined.",
        ),
        ReductionRow(
            "jdopfj_15_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 Xi obligation is now a precise phase-lift transport and endpoint flux theorem.",
            exact["live_q_ge_1_target"],
            "No Xi lower bound or strict successor upper bound is yet proved.",
        ),
        ReductionRow(
            "jdopfj_16_q_lt_1_nonpromotion",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1 multiplicity-compatible chart and every global proof stage remain open.",
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
            "exact joined dyadic odd-prefix heat-shift and five-current "
            "first-jet reduction with an endpoint-complete cutoff law, "
            "a proved phase-synchronized zero-free homotopy base, and "
            "two exact route guards; the physical Xi phase transport "
            "and successor theorem remain open"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "constants": {
            "L_min": L_MIN,
            "uncorrected_dyadic_ratio": UNCORRECTED_DYADIC_RATIO,
            "corrected_dyadic_ratio": CORRECTED_DYADIC_RATIO,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Joined Dyadic Odd-Prefix First-Jet Reduction",
            "",
            "Date: 2026-07-27",
            "",
            "Status: exact joined reduction and hypothetical homotopy-base",
            "theorem; not a proof of `Lambda<=0`, RH, PF-infinity, or a",
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
            "## Deterministic Heat Shift",
            "",
            "```text",
            exact["heat_shift"],
            "```",
            "",
            "This factorization is stronger than merely writing the heat",
            "quadratic as a Gaussian expectation: each dyadic layer is an",
            "exactly shifted odd prefix before any averaging.",
            "",
            "## Joined Odd Prefixes",
            "",
            "```text",
            exact["coefficient_chain"],
            exact["joined_odd_prefix"],
            exact["joined_correction_current"],
            "```",
            "",
            "Unique dyadic valuation is doing the bookkeeping. At a cutoff",
            "jump, exactly one odd-prefix layer gains the entering integer.",
            "",
            "## Phase Lift",
            "",
            "```text",
            exact["phase_lift"],
            exact["phase_lift_derivatives"],
            "```",
            "",
            "The variables `z` and `xi` separate the dyadic phase from the",
            "odd Mellin phase. The physical Xi prefix sits on the linked",
            "curve `z=exp(i*omega_*log(2))`, `xi=omega_*`.",
            "",
            "## Five-Current First Jet",
            "",
            "```text",
            exact["five_current_ladder"],
            exact["centered_endpoint_jet"],
            exact["projected_flux"],
            "```",
            "",
            "Thus the exact first jet needs only `H_0,H_1,H_2,D_0,D_1`",
            "plus endpoint value and derivative data. This is a reduction",
            "of organization, not a sign theorem.",
            "",
            "## Cutoff And Endpoint",
            "",
            "```text",
            exact["cutoff_jump"],
            exact["adjacent_control"],
            "```",
            "",
            "The terminal recurrence is part of the joined representation.",
            "It is not appended after a bulk winding calculation.",
            "",
            "## Phase-Synchronized Hypothetical",
            "",
            "The deliberately deformed case `xi=mu=0` can be solved exactly:",
            "",
            "```text",
            exact["phase_frozen_hypothetical"],
            exact["phase_frozen_margin"],
            "```",
            "",
            "This is the useful thought experiment. Synchronizing all odd",
            "phases produces a zero-free joined polynomial. The physical",
            "problem is therefore an explicit transport question from that",
            "proved base to the actual Mellin phases, rather than an",
            "undefined appeal to multiplicativity.",
            "",
            "## Homotopy Target",
            "",
            "```text",
            exact["homotopy_target"],
            "```",
            "",
            "A successful theorem may require exponential-polynomial",
            "Hermite-Biehler/de Branges machinery, a signed two-parameter",
            "Jacobian, or a linked-curve argument-principle estimate.",
            "",
            "## Route Guards",
            "",
            "### Averaging does not preserve winding",
            "",
            "```json",
            json.dumps(exact["two_fibre_guard"], indent=2, sort_keys=True),
            "```",
            "",
            "### The second logarithmic moment is indispensable",
            "",
            "```json",
            json.dumps(exact["second_moment_guard"], indent=2, sort_keys=True),
            "```",
            "",
            "These are generic exact guards, not Xi counterexamples.",
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
        "built Newman joined dyadic odd-prefix first-jet reduction: "
        "17 rows, 1 exact heat-shift factorization, "
        "1 five-current closure, "
        "1 phase-frozen joined zero-free theorem, "
        "2 exact route guards, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
