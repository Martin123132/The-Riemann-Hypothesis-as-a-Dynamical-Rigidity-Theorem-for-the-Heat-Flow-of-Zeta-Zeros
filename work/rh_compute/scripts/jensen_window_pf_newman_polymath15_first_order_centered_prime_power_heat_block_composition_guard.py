#!/usr/bin/env python3
"""Build the prime-power heat-block and cross-block composition guard."""

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
    "prime_power_heat_block_composition_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
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
    "adjacent_chart": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_chart_stability_certificate.json"
    ),
    "phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
}

L_MIN = 50
D_ABSOLUTE_CONSTANT = 2_189
D_X_CONSTANT = 4_223
TERMINAL_RATE_CONSTANT = 16_900
TERMINAL_WIDTH_CONSTANT = 50_701


@dataclass(frozen=True)
class BlockRow:
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
        "normalized_prefix": (
            "every adjacent relative phase advances by more than 4*pi",
            "dyadic_completion_scout",
            "u_(N,x)=1/(8*pi*a^2)=1/(4*T_0)",
        ),
        "direct_projection": (
            "u_n=log(a/n)=-ell_n",
            "eta*q_n=r_n*zeta_n",
            "|v_a|<3/x^2",
        ),
        "adjacent_chart": (
            "|A_(a,N+1)-A_(a,N)|<5000*exp(-7L/4)",
            "5000*exp(-7L/4)<10^-7*exp(-5L/4)",
        ),
        "phase_anchor": (
            "eta=f_1/|f_1|",
            "q_n=f_n/f_1=",
            "q_1=1",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    direct_constants = payloads["direct_projection"].get("constants", {})
    if direct_constants.get("d_absolute") != D_ABSOLUTE_CONSTANT:
        raise RuntimeError("direct-projection d absolute constant drifted")
    if direct_constants.get("d_x") != D_X_CONSTANT:
        raise RuntimeError("direct-projection d derivative constant drifted")
    return {
        key: str(value.get("kind", ""))
        for key, value in payloads.items()
    }


def two_block_countermodel() -> dict[str, str | int]:
    z = sp.symbols("z")
    r = sp.sqrt(2) / 2
    scale = sp.Rational(4, 5)
    long_block = sp.expand(z * (1 + r * z + r**2 * z**2))
    short_block = sp.expand(-scale * z * (1 + r * z))
    joined = sp.factor(long_block + short_block)
    joined_inner = sp.expand(5 * joined / z)
    expected_inner = 1 + r * z + 5 * r**2 * z**2
    if sp.simplify(joined_inner - expected_inner) != 0:
        raise RuntimeError("two-block joined polynomial failed")

    joined_discriminant = sp.discriminant(joined_inner, z)
    joined_product = sp.simplify(
        sp.Poly(joined_inner, z).all_coeffs()[-1]
        / sp.Poly(joined_inner, z).all_coeffs()[0]
    )
    if joined_discriminant != -sp.Rational(19, 2):
        raise RuntimeError("two-block discriminant failed")
    if joined_product != sp.Rational(2, 5):
        raise RuntimeError("two-block root product failed")

    long_inner = 1 + r * z + r**2 * z**2
    short_inner = 1 + r * z
    long_product = sp.simplify(
        sp.Poly(long_inner, z).all_coeffs()[-1]
        / sp.Poly(long_inner, z).all_coeffs()[0]
    )
    short_root = sp.solve(short_inner, z)[0]
    if long_product != 2 or short_root != -sp.sqrt(2):
        raise RuntimeError("individual block root radii failed")

    return {
        "dyadic_ratio": "r=2^(-1/2)",
        "long_block": "A(z)=z*(1+r*z+r^2*z^2)",
        "short_block": "B(z)=-(4/5)*z*(1+r*z)",
        "individual_root_locations": (
            "A/z has roots of modulus sqrt(2)>1; "
            "B/z has its root at -sqrt(2)"
        ),
        "individual_windings": "wind(A)=wind(B)=1 on |z|=1",
        "joined_factorization": (
            "A+B=(z/5)*(1+r*z+5*r^2*z^2)"
        ),
        "joined_quadratic_discriminant": "-19/2",
        "joined_quadratic_root_product": "2/5",
        "joined_quadratic_root_modulus": "sqrt(2/5)<1",
        "joined_winding": 3,
        "robust_scale_interval": (
            "For A-c*z*(1+r*z), every 1/2<c<1 gives a "
            "conjugate root pair of modulus sqrt(2*(1-c))<1."
        ),
        "interpretation": (
            "Two complete geometric dyadic blocks can each be "
            "origin-free with winding one while their sum has winding "
            "three. The opposite external block phase and scale are a "
            "generic composition guard, not asserted Xi data."
        ),
    }


def symbolic_audit() -> dict[str, str | dict]:
    omega_phi, b, log_a, log_n, correction = sp.symbols(
        "omega_phi b log_a log_n correction", real=True
    )
    u_n = log_a - log_n
    v_a = omega_phi - b * log_a
    absolute_rate = omega_phi - b * log_n + correction
    centered_rate = v_a + b * u_n + correction
    if sp.expand(absolute_rate - centered_rate) != 0:
        raise RuntimeError("absolute saddle-frame rate failed")

    b_amp, time, u_value, h = sp.symbols(
        "B_a time u_value h", real=True
    )
    log_amp = b_amp * u_value + time * u_value**2 / 4
    next_log_amp = (
        b_amp * (u_value - h)
        + time * (u_value - h) ** 2 / 4
    )
    ratio = sp.factor(next_log_amp - log_amp)
    expected_ratio = (
        -b_amp * h
        - time * u_value * h / 2
        + time * h**2 / 4
    )
    if sp.expand(ratio - expected_ratio) != 0:
        raise RuntimeError("prime-power amplitude ratio failed")

    sigma = sp.symbols("sigma", positive=True, real=True)
    if sp.simplify((sp.sqrt(time / 2)) ** 2 / 2 - time / 4) != 0:
        raise RuntimeError("Gaussian heat variance failed")

    joined_pairs: list[tuple[int, int, int]] = []
    n_max = 24
    prime = 2
    for base in range(1, n_max + 1):
        if base % prime == 0:
            continue
        power = 0
        value = base
        while value <= n_max:
            joined_pairs.append((value, base, power))
            power += 1
            value *= prime
    if sorted(value for value, _, _ in joined_pairs) != list(
        range(1, n_max + 1)
    ):
        raise RuntimeError("joined p-free valuation decomposition failed")

    return {
        "absolute_rate": (
            "For z_n=eta*q_n=r_n*zeta_n and "
            "zeta_n=phi*exp(i*omega_*log(n))"
            "*(1+d_n)/|1+d_n|, put "
            "omega_phi=Im(phi'/phi), b=Im(s_*'), "
            "u_n=log(a/n), and "
            "v_a=omega_phi-b*log(a). Since omega_*'=-b, "
            "nu_n=partial_x arg(z_n)="
            "v_a+b*u_n+Im[d_(n,x)/(1+d_n)]."
        ),
        "prime_power_ratio": (
            "For p prime, p does not divide m, n_k=m*p^k, "
            "h=log(p), and u_k=log(a/n_k), the correction-free "
            "absolute amplitudes satisfy "
            "log(a_(k+1)/a_k)="
            "-B_a*h-(t/2)*u_k*h+(t/4)*h^2. "
            "If n_(k+1)<=N<=a, then u_k>=h."
        ),
        "gaussian_heat": (
            "For Z standard normal and sigma=sqrt(t/2), "
            "exp(t*log(n)^2/4)=E exp(sigma*Z*log(n))."
        ),
        "joined_p_free": (
            "With s_Z=s_*-sigma*Z and "
            "O_M^(p)(s_Z)=sum_(m<=M,p does not divide m)m^(-s_Z), "
            "sum_(n<=N)n^(-s_Z)="
            "sum_(k=0)^floor(log_p N)p^(-k*s_Z)"
            "*O_floor(N/p^k)^(p)(s_Z)."
        ),
        "two_block_countermodel": two_block_countermodel(),
    }


def numeric_audit() -> dict[str, str]:
    x_min = 4 * math.pi * math.exp(L_MIN)
    d_upper = D_ABSOLUTE_CONSTANT / x_min
    correction_quotient = (
        D_X_CONSTANT / x_min**2 / (1 - d_upper)
    )
    rate_error = 3 / x_min**2 + correction_quotient
    terminal_negative_margin = (
        TERMINAL_RATE_CONSTANT / 2
        - rate_error * x_min**2
    )
    if not terminal_negative_margin > 0:
        raise RuntimeError("terminal rate threshold failed")
    if not x_min ** 1.5 > TERMINAL_RATE_CONSTANT:
        raise RuntimeError("nonterminal rate separation failed")

    worst_amplitude_ratio = (
        2 ** (-49 / 100) * (1 + d_upper) / (1 - d_upper)
    )
    amplitude_target = 2 ** (-12 / 25)
    if not worst_amplitude_ratio < amplitude_target:
        raise RuntimeError("corrected amplitude contraction failed")

    phase_defect = 2 * d_upper / (1 - d_upper)
    if not phase_defect < 8_756 / x_min:
        raise RuntimeError("relative correction phase failed")

    collar_scaled_upper = 50_700 * (1 + math.pi / (8 * x_min))
    if not collar_scaled_upper < TERMINAL_WIDTH_CONSTANT:
        raise RuntimeError("terminal collar width failed")

    gaussian_tail_at_one = math.erfc(1 / math.sqrt(2)) / 2
    if not 0.158 < gaussian_tail_at_one < 0.159:
        raise RuntimeError("Gaussian expanding-ratio calibration failed")

    return {
        "x_min": format(x_min, ".17e"),
        "d_n_upper_at_L50": format(d_upper, ".17e"),
        "d_log_derivative_upper_at_L50": format(
            correction_quotient, ".17e"
        ),
        "terminal_rate_threshold": "16900/x^2",
        "terminal_negative_margin_scaled_by_x2": format(
            terminal_negative_margin, ".17e"
        ),
        "terminal_x_collar_width": "<50701/x",
        "worst_corrected_p_power_ratio": format(
            worst_amplitude_ratio, ".17e"
        ),
        "certified_ratio_target": format(amplitude_target, ".17e"),
        "relative_correction_phase_at_L50": format(
            phase_defect, ".17e"
        ),
        "gaussian_tail_at_terminal_t_half": format(
            gaussian_tail_at_one, ".17e"
        ),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    numeric = numeric_audit()
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a); q>=1 means "
            "q=2tL^2>=1"
        ),
        "pi_provenance": (
            "The pi in a^2=x/(4*pi)+t/16 comes from the completed "
            "zeta normalization "
            "xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2 "
            "and its Riemann-Siegel saddle. The cutoff-cell difference "
            "is 4*pi*(2N+1). The 2*pi in a winding is the period of "
            "exp(i*theta). No circle, prime-power block, or prefix "
            "polygon defines or approximates pi."
        ),
        "absolute_rate": symbolic["absolute_rate"],
        "terminal_rate_collar": (
            "The proved bounds |v_a|<3/x^2, "
            "|d_(n,x)|<4223/x^2, and |d_n|<2189/x<1/2 give "
            "nu_n<=8449/x^2-u_n/2. Therefore "
            "u_n>16900/x^2 implies nu_n<-1/x^2. For n<=N-1, "
            "u_n>=log(N/(N-1))>1/N>16900/x^2 on L>=50, so only "
            "the terminal n=N carrier can be near stationary. If "
            "u_N<=16900/x^2, then at fixed t its distance from the "
            "left cutoff boundary obeys "
            "0<=x-x_N=4*pi*N^2*(exp(2u_N)-1)<50701/x."
        ),
        "prime_power_ratio": symbolic["prime_power_ratio"],
        "prime_power_contraction": (
            "Because B_a>49/100 and u_k>=h whenever the next "
            "p-power remains in the cutoff, the correction-free ratio "
            "is at most exp[-49h/100-th^2/4]<p^(-49/100). "
            "After |1+d_(k+1)|/|1+d_k| is restored, the ratio is "
            "still strictly below p^(-12/25) for L>=50. Thus every "
            "actual complete p-power block has strictly decreasing "
            "positive amplitudes."
        ),
        "correction_phase": (
            "After the common absolute saddle phase and affine "
            "p-power phase are removed, "
            "|arg(1+d_(mp^k))-arg(1+d_m)|<8756/x. "
            "This is an exponentially small phase defect, but a "
            "uniform boundary-modulus margin is still required before "
            "it can be removed by homotopy."
        ),
        "uncorrected_block_winding": (
            "Let Q_R(z)=sum_(k=0)^R a_k z^k with "
            "a_0>a_1>...>a_R>0. Enestrom-Kakeya puts every zero of "
            "Q_R outside |z|=1 because min_k a_k/a_(k+1)>1. Hence "
            "Q_R is nonzero on the unit circle and wind(Q_R,0)=0. "
            "The shifted pure-power block P_R(z)=z*Q_R(z) is also "
            "nonzero there and has wind(P_R,0)=1. For a general "
            "p-free base m, Q_R is the normalized internal block after "
            "its external m-phase is factored out; that external phase "
            "still moves with x and is not covered by the internal "
            "winding statement."
        ),
        "gaussian_heat": symbolic["gaussian_heat"],
        "conditional_geometric_block": (
            "With sigma=sqrt(t/2), the correction-free block is "
            "m^(-s_*) E[m^(sigma Z)"
            "*sum_(k=0)^R(p^(-s_*+sigma Z))^k]. "
            "The geometric quotient is interpreted by continuity as "
            "R+1 when its ratio is one."
        ),
        "tilted_expanding_mass": (
            "Conditionally, |p^(-s_*+sigma Z)|>=1 exactly when "
            "Z>=Re(s_*)/sigma. After the factor m^(sigma Z) is "
            "included, exponential tilting makes the relative bad mass "
            "barPhi((1/2+(t/2)*(u_m-delta_a))/sigma), where "
            "delta_a=log(a)-Re(alpha). For a terminal block at "
            "t=1/2 this threshold approaches 1 and the mass approaches "
            "barPhi(1)=0.158655..., so conditional contraction cannot "
            "be asserted almost surely or discarded uniformly from the "
            "Gaussian identity alone."
        ),
        "joined_p_free": symbolic["joined_p_free"],
        "two_block_countermodel": symbolic["two_block_countermodel"],
        "composition_guard": (
            "The one-block theorem is real, but winding is not additive "
            "under vector addition. The exact length-three/length-two "
            "dyadic model has two complete origin-free winding-one "
            "blocks and a winding-three sum. Therefore a successful Xi "
            "proof must retain the p-free outer prefixes, their relative "
            "phases and amplitudes, the terminal recurrence, and the "
            "endpoint inside one joined first-jet argument theorem."
        ),
        "live_q_ge_1_target": (
            "Use the joined p-free decomposition, preferably first for "
            "p=2, to derive an endpoint-complete common argument-principle "
            "or Schur-Cohn/Rouche inequality for the full Xi prefix. It "
            "must prove both |mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term and the composed successor "
            "bound 0<=kappa_j<1. The terminal u_N<=16900/x^2 collar "
            "must be paired with the adjacent endpoint recurrence. "
            "Independent block windings and the raw Gaussian mixture "
            "are insufficient."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 layer still requires its separate "
            "multiplicity-compatible parabolic/Hermite first-jet chart. "
            "No terminal simplicity or uniform positive endpoint slope "
            "is assumed."
        ),
        "numeric_audit": numeric,
    }


def build_rows(exact: dict) -> list[BlockRow]:
    return [
        BlockRow(
            "pphbcg_00_pi_provenance",
            "definition_provenance",
            "certified",
            "Every appearance of pi in this block reduction has an explicit completed-zeta, saddle, or phase-period source.",
            exact["pi_provenance"],
            "Neither a block nor a polygon supplies pi.",
        ),
        BlockRow(
            "pphbcg_01_absolute_rate",
            "exact_differential_identity",
            "ready_to_apply",
            "The absolute carrier phase rate has an exact saddle-centered form.",
            exact["absolute_rate"],
            "The d_n logarithmic derivative is retained.",
        ),
        BlockRow(
            "pphbcg_02_terminal_collar",
            "asymptotic_certificate",
            "certified",
            "Every nonterminal carrier rotates strictly in the same absolute direction; only one O(1/x) cutoff collar can be near stationary.",
            exact["terminal_rate_collar"],
            "The terminal carrier still has to be paired with the endpoint recurrence.",
            exact["numeric_audit"],
        ),
        BlockRow(
            "pphbcg_03_prime_power_amplitude",
            "exact_reduction",
            "ready_to_apply",
            "Every complete p-power chain has an exact heat-quadratic ratio law.",
            exact["prime_power_ratio"],
            "The law is inside one prescribed cutoff chart.",
        ),
        BlockRow(
            "pphbcg_04_corrected_contraction",
            "asymptotic_certificate",
            "certified",
            "The actual corrected amplitudes in every complete p-power chain decrease geometrically.",
            exact["prime_power_contraction"],
            exact["correction_phase"],
            exact["numeric_audit"],
        ),
        BlockRow(
            "pphbcg_05_one_block_winding",
            "exact_theorem",
            "ready_to_apply",
            "A complete correction-free normalized prime-power block is zero-free in its internal phase, and its one-step shifted pure-power block has winding one.",
            exact["uncorrected_block_winding"],
            "The moving p-free external phase and d_n phase defects remain outside this internal theorem.",
        ),
        BlockRow(
            "pphbcg_06_gaussian_mixture",
            "exact_identity",
            "ready_to_apply",
            "The heat quadratic is an exact Gaussian mixture of geometric Dirichlet blocks.",
            exact["gaussian_heat"] + " " + exact["conditional_geometric_block"],
            "Expectation does not preserve zero-free winding under addition.",
        ),
        BlockRow(
            "pphbcg_07_tilted_mass_guard",
            "nonpromotion_guard",
            "guard_validated",
            "The Gaussian representation has a non-negligible terminal conditional expanding-ratio sector.",
            exact["tilted_expanding_mass"],
            "The exact mixture identity alone is not a conditional multiplicativity theorem.",
            exact["numeric_audit"],
        ),
        BlockRow(
            "pphbcg_08_joined_p_free",
            "exact_reindexing",
            "ready_to_apply",
            "All complete p-power blocks rejoin exactly through nested p-free outer prefixes.",
            exact["joined_p_free"],
            "The p-free prefixes remain complex and need an Xi-specific estimate.",
        ),
        BlockRow(
            "pphbcg_09_two_block_countermodel",
            "countermodel",
            "guard_validated",
            "Two complete winding-one dyadic blocks can have a winding-three sum.",
            json.dumps(exact["two_block_countermodel"], sort_keys=True),
            "This is an exact generic composition guard, not an Xi counterexample.",
            exact["two_block_countermodel"],
        ),
        BlockRow(
            "pphbcg_10_composition_guard",
            "route_decision",
            "guard_validated",
            "Independent prime-power windings cannot be promoted to the full successor flux.",
            exact["composition_guard"],
            "Cross-block, terminal, endpoint, and connector terms must remain joined.",
        ),
        BlockRow(
            "pphbcg_11_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 problem is now a joined p-free-prefix and endpoint first-jet theorem.",
            exact["live_q_ge_1_target"],
            "The Xi lower bound and strict successor flux upper bound remain open.",
        ),
        BlockRow(
            "pphbcg_12_q_lt_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1 multiplicity-compatible chart remains a separate theorem obligation.",
            exact["q_lt_1_target"],
            "No endpoint simplicity is promoted.",
        ),
        BlockRow(
            "pphbcg_13_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "The block identities and guards do not prove contact exclusion or RH.",
            exact["composition_guard"],
            (
                "No full Xi prefix bound, successor theorem, Lambda<=0, "
                "PF-infinity, RH, or Clay-prize conclusion is asserted."
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
            "exact absolute-rate, prime-power heat-block, Gaussian-mixture, "
            "and joined p-free reductions with an uncorrected normalized-"
            "block zero-free and shifted-winding theorem plus an exact "
            "two-block composition guard; "
            "the joined Xi successor theorem remains open"
        ),
        "proof_boundary": (
            "This artifact proves the absolute saddle-frame carrier rate, "
            "the nonterminal rotation and terminal-collar theorem, exact "
            "prime-power heat ratios, corrected amplitude contraction, "
            "the correction-free normalized-block zero-free and shifted-"
            "winding theorem, the Gaussian "
            "mixture and joined p-free identities, and an exact generic "
            "two-block winding-three countermodel. It does not prove a "
            "uniform corrected-block homotopy, the joined Xi prefix lower "
            "bound, the strict composed successor flux bound, the q<1 "
            "multiplicity-compatible theorem, finite shoulder closure, "
            "contact exclusion, Lambda<=0, PF-infinity, RH, or a "
            "Clay-prize conclusion."
        ),
        "sources": {key: str(path.relative_to(REPO_ROOT)) for key, path in SOURCES.items()},
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "constants": {
            "L_min": L_MIN,
            "d_absolute": D_ABSOLUTE_CONSTANT,
            "d_x": D_X_CONSTANT,
            "terminal_rate": TERMINAL_RATE_CONSTANT,
            "terminal_width": TERMINAL_WIDTH_CONSTANT,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    guard = exact["two_block_countermodel"]
    lines = [
        "# Newman Prime-Power Heat-Block Composition Guard",
        "",
        "Date: 2026-07-27",
        "",
        "Status: Exact reduction and countermodel-guard artifact; not a proof of",
        "`Lambda<=0`, RH, PF-infinity, or a Clay-prize result.",
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
        "## Absolute Carrier Rate",
        "",
        "Write the absolute carrier as `z_n=eta q_n=r_n zeta_n`. The exact",
        "phase-rate identity is",
        "",
        "```text",
        "nu_n=partial_x arg(z_n)",
        "    =v_a+b u_n+Im[d_(n,x)/(1+d_n)],",
        "u_n=log(a/n).                                      (1)",
        "```",
        "",
        "The certified bounds give",
        "",
        "```text",
        "nu_n<=8449/x^2-u_n/2.                              (2)",
        "```",
        "",
        "Thus every `n<=N-1` carrier rotates strictly in the same absolute",
        "direction. Only `n=N` can be near stationary, and that can happen",
        "only in",
        "",
        "```text",
        "u_N<=16900/x^2,",
        "0<=x-x_N<50701/x.                                  (3)",
        "```",
        "",
        "This is the terminal collar that must be combined with the already",
        "proved adjacent endpoint recurrence.",
        "",
        "## Prime-Power Blocks",
        "",
        "Fix a prime `p`, a `p`-free integer `m`, and",
        "`n_k=m p^k<=N`. Put `h=log p` and `u_k=log(a/n_k)`. Before the",
        "`d_n` correction, consecutive amplitudes satisfy",
        "",
        "```text",
        "log(a_(k+1)/a_k)",
        " =-B_a h-(t/2)u_k h+(t/4)h^2.                     (4)",
        "```",
        "",
        "Whenever the next power remains below the cutoff, `u_k>=h`. Since",
        "`B_a>49/100`, restoring the absolute `1+d_n` amplitudes still gives",
        "",
        "```text",
        "a_(k+1)/a_k<p^(-12/25)<1.                          (5)",
        "```",
        "",
        "The residual phase defect relative to one affine block phase is",
        "less than `8756/x`. It is tiny, but a quantitative boundary-modulus",
        "margin is needed before a homotopy may discard it.",
        "",
        "## One-Block Theorem",
        "",
        "For positive strictly decreasing coefficients, let",
        "",
        "```text",
        "Q_R(z)=sum_(k=0)^R a_k z^k,",
        "P_R(z)=z Q_R(z).                                   (6)",
        "```",
        "",
        "Enestrom-Kakeya puts every zero of `Q_R` outside the closed unit",
        "disk because `min a_k/a_(k+1)>1`. Thus `Q_R` is nonzero with",
        "internal winding zero, while the shifted pure-power block `P_R`",
        "has winding one. For a general p-free base, its external phase was",
        "factored out and still moves with `x`; this is not an actual",
        "blockwise `x`-winding theorem for the joined Xi sum.",
        "",
        "## Gaussian Heat Identity",
        "",
        "For a standard normal `Z` and `sigma=sqrt(t/2)`,",
        "",
        "```text",
        "exp(t log(n)^2/4)=E exp(sigma Z log n).             (7)",
        "```",
        "",
        "Hence a correction-free block is exactly an expectation of finite",
        "geometric blocks. The tempting pointwise argument nevertheless",
        "fails: after exponential tilting by the block base, the conditional",
        "expanding-ratio mass is",
        "",
        "```text",
        "barPhi((1/2+(t/2)(u_m-delta_a))/sigma).             (8)",
        "```",
        "",
        "For a terminal block at `t=1/2`, this approaches",
        "`barPhi(1)=0.158655...`. The Gaussian identity therefore does not",
        "make every conditional block contract, and expectation does not",
        "preserve zero-free winding.",
        "",
        "## Two-Block Guard",
        "",
        "Let `r=2^(-1/2)` and",
        "",
        "```text",
        f"{guard['long_block']},",
        f"{guard['short_block']}.",
        "```",
        "",
        "Each block is nonzero on the unit circle and has winding one. Their",
        "sum is",
        "",
        "```text",
        "A+B=(z/5)(1+r z+5r^2 z^2).                         (9)",
        "```",
        "",
        "The inner quadratic has discriminant `-19/2`; its two roots have",
        "modulus `sqrt(2/5)<1`. Consequently the sum has winding three.",
        "More generally this obstruction persists for every opposite block",
        "scale `1/2<c<1`. This is an exact generic composition guard, not an Xi counterexample.",
        "",
        "## Joined Structure",
        "",
        "The blocks must be rejoined before winding is estimated. For",
        "`s_Z=s_*-sigma Z`, define the `p`-free prefix",
        "",
        "```text",
        "O_M^(p)(s_Z)=sum_(m<=M, p does not divide m)m^(-s_Z).",
        "```",
        "",
        "Then",
        "",
        "```text",
        "sum_(n<=N)n^(-s_Z)",
        " =sum_k p^(-k s_Z) O_floor(N/p^k)^(p)(s_Z).        (10)",
        "```",
        "",
        "This identity retains the cross-block phases that the rejected",
        "blockwise count loses.",
        "",
        "## Live Theorem",
        "",
        "The next `q>=1` theorem must combine the joined `p`-free prefixes,",
        "the terminal recurrence, the exact endpoint, and the full first-jet",
        "argument. It must establish",
        "",
        "```text",
        "|mathsf_X|<=delta_L => |mathcal_C_N|>A_L+epsilon_term,",
        "0<=kappa_j<1,                                      (11)",
        "```",
        "",
        "with all connector, chart-join, and shoulder terms retained. The",
        "`q<1` multiplicity-compatible parabolic/Hermite chart remains a",
        "separate obligation.",
        "",
        "## Boundary",
        "",
        artifact["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman prime-power heat-block guard: "
        "14 rows, 1 absolute-rate collar, "
        "1 exact normalized-block zero-free/shifted-winding theorem, "
        "1 Gaussian-mixture audit, "
        "1 joined p-free decomposition, "
        "1 exact two-block winding-3 countermodel, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
