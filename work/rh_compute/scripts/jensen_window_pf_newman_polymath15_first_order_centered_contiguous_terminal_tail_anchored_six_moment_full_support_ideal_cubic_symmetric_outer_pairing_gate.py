#!/usr/bin/env python3
"""Build the ideal-cubic symmetric outer-frequency pairing gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_ideal_cubic_"
    "symmetric_outer_pairing_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "finite_cell_inversion": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_finite_cell_inversion_gate.json"
    ),
    "outer_endpoint_kernel": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_outer_endpoint_"
        "kernel_gate.json"
    ),
    "joined_recombination": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_joined_lower_"
        "interior_abel_recombination_gate.json"
    ),
    "ridge_carrier": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_ideal_cubic_double_"
        "morse_ridge_carrier_inversion_gate.json"
    ),
}


@dataclass(frozen=True)
class GateRow:
    id: str
    role: str
    readiness: str
    claim: str
    certificate: str
    proof_boundary: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(sp.factor(sp.together(expression)))
    if value != 0 and value.equals(0) is not True:
        raise RuntimeError(f"{label}: {value}")


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        require(path.exists(), f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    require(
        payloads["finite_cell_inversion"]["counts"]["carrier_inversion_identities"]
        == 3,
        "finite-cell inversion source drift",
    )
    require(
        payloads["outer_endpoint_kernel"]["counts"]
        ["complete_outer_complement_bounds"]
        == 0,
        "outer-endpoint source drift",
    )
    require(
        payloads["joined_recombination"]["counts"]["global_recombinations"]
        == 4,
        "joined-recombination source drift",
    )
    require(
        payloads["ridge_carrier"]["counts"]["inversion_identities"] == 5,
        "ridge/carrier source drift",
    )
    return payloads


def source_audit() -> dict[str, dict[str, str]]:
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def recombination_certificate() -> dict[str, str | int]:
    carrier, exterior, tie, terminal = sp.symbols("M E tau T")
    band = tie + carrier - exterior
    outer = exterior - tie
    defect = band - carrier
    checks = [
        band + outer - carrier,
        outer - (carrier - band),
        defect + outer,
        (carrier - tie + exterior) - (carrier + outer),
        (carrier - tie + exterior + terminal) - (carrier + outer + terminal),
    ]
    for index, expression in enumerate(checks):
        require_zero(expression, f"outer recombination {index}")
    return {
        "poisson": "Full-support Poisson summation gives mathscr F_N[P]+mathscr O_N[P]=mathscr M_N[P], where mathscr O_N is the common symmetric sum over all integer modes outside mathcal R_N.",
        "exterior": "Because mathscr F_N=tau_band+mathscr M_N-E_band, one has mathscr O_N=E_band-tau_band=-mathscr D_N exactly.",
        "hermitian": "The ideal Hermitian join is mathcal J_H^0=mathscr M_N[P_H^0]+mathscr O_N[P_H^0].",
        "transpose": "The ideal transpose join is mathcal J_T^0=mathscr M_N[P_T^0]+mathscr O_N[P_T^0]+2iu_N mathcal T_N^[N][B_(T,N)^0].",
        "guard": "These are recombinations of the same Fourier involution. They do not imply that either complete join is small or has a sign.",
        "identities": len(checks),
    }


def cutoff_certificate() -> dict[str, str | int]:
    exact_roster_checks = 0
    for m_value in range(1, 5):
        for n_value in range(m_value, 8):
            for radius in range(n_value, n_value + 4):
                direct = [
                    r
                    for r in range(-radius, radius + 1)
                    if not m_value <= r <= n_value
                ]
                split = list(range(1, m_value)) + [0]
                split += [-k for k in range(1, n_value + 1)]
                for k in range(n_value + 1, radius + 1):
                    split.extend([-k, k])
                require(sorted(direct) == sorted(split), "finite outer roster split")
                require(len(split) == len(set(split)), "duplicate outer roster mode")
                exact_roster_checks += 1

    mp.mp.dps = 70
    m_value, n_value, radius = 2, 5, 9
    character_checks = 0
    for u_value in map(
        mp.mpf,
        ["0.071", "0.193", "0.347", "0.511", "0.733", "0.889", "1.137"],
    ):
        direct = mp.fsum(
            mp.e ** (-2j * mp.pi * r * u_value)
            for r in range(-radius, radius + 1)
            if not m_value <= r <= n_value
        )
        near = mp.fsum(
            mp.e ** (2j * mp.pi * j * u_value)
            for j in range(-(m_value - 1), n_value + 1)
        )
        remote = 2 * mp.fsum(
            mp.cos(2 * mp.pi * k * u_value)
            for k in range(n_value + 1, radius + 1)
        )
        require(abs(direct - near - remote) < mp.mpf("1e-60"), "character split")
        character_checks += 1

    mean_checks = 0
    for cell in range(-3, 5):
        integral = mp.quad(
            lambda u: mp.fsum(
                mp.e ** (2j * mp.pi * j * u)
                for j in range(-(m_value - 1), n_value + 1)
            ),
            [cell, cell + 1],
        )
        require(abs(integral - 1) < mp.mpf("1e-58"), "near-kernel mean")
        mean_checks += 1

    return {
        "indices": "Put m_N=floor(alpha_P/(N+1/2))+1 and n_N=floor(2alpha_P), with mathcal R_N={m_N,...,n_N}. For every R>=n_N, the finite symmetric outer sum has one canonical disjoint split.",
        "finite_split": "sum_(-R<=r<=R,r notin mathcal R_N)I_P(r)=sum_(r=1)^(m_N-1)I_P(r)+I_P(0)+sum_(k=1)^n_N I_P(-k)+sum_(k=n_N+1)^R{I_P(-k)+I_P(k)}.",
        "near": "Define mathscr N_N[P]=sum_(r=1)^(m_N-1)I_P(r)+I_P(0)+sum_(k=1)^n_N I_P(-k). It is finite and equals integral_1^N f_P(u)K_near,N(u)du.",
        "kernel": "K_near,N(u)=sum_(j=-(m_N-1))^n_N e(ju). Its mean on every integer unit cell is exactly 1 because its only zero Fourier mode has coefficient 1.",
        "remote": "The remaining remote object is the paired series mathscr C_N[P]=sum_(k=n_N+1)^infinity{I_P(-k)+I_P(k)}. The common symmetric cutoff, not a rearrangement of two one-sided series, defines the passage to this limit.",
        "exact_roster_checks": exact_roster_checks,
        "character_checks": character_checks,
        "mean_checks": mean_checks,
    }


def paired_tail_certificate() -> dict[str, str | int]:
    mp.mp.dps = 70
    terminal = mp.mpf(4)
    pair_checks = 0
    for k in [6, 7, 9, 13, 17]:
        direct = mp.quad(
            lambda u: u * u * (
                mp.e ** (2j * mp.pi * k * u)
                + mp.e ** (-2j * mp.pi * k * u)
            ),
            [1, 2, 3, 4],
        )
        expected = (terminal - 1) / (mp.pi**2 * k * k)
        require(abs(direct - expected) < mp.mpf("1e-58"), "paired tail formula")
        pair_checks += 1

    tail_bound_checks = 0
    for cutoff in [2, 4, 7, 12, 20]:
        zeta_tail = mp.zeta(2) - mp.fsum(
            mp.mpf(1) / (k * k) for k in range(1, cutoff + 1)
        )
        exact_tail = (terminal - 1) * zeta_tail / (mp.pi**2)
        derivative_bound = 2 * (terminal - 1) * zeta_tail / (mp.pi**2)
        coarse_bound = 2 * (terminal - 1) / (mp.pi**2 * cutoff)
        require(exact_tail <= derivative_bound, "absolute derivative tail bound")
        require(derivative_bound <= coarse_bound, "coarse reciprocal tail bound")
        tail_bound_checks += 2

    return {
        "pair": "For every positive integer k, I_P(-k)+I_P(k)=2integral_1^N f_P(u)cos(2pi ku)du.",
        "twice_ibp": "If f_P is C^2 on [1,N], integer endpoints make both sine boundary values zero, and I_P(-k)+I_P(k)={2/(2pi k)^2}{f_P'(N)-f_P'(1)-integral_1^N f_P''(u)cos(2pi ku)du}.",
        "bound": "Hence |mathscr C_N[P]|<={|f_P'(N)-f_P'(1)|+||f_P''||_1}/{2pi^2} sum_(k>n_N)k^(-2)<={|f_P'(N)-f_P'(1)|+||f_P''||_1}/{2pi^2 n_N}.",
        "scope": "This proves absolute convergence of the paired remote tail. It is not an h^2 estimate: alpha_P occurs in f_P' and f_P'', so the displayed derivative numerator must still be analyzed in the physical scaling.",
        "pair_formula_checks": pair_checks,
        "tail_bound_checks": tail_bound_checks,
    }


def endpoint_certificate() -> dict[str, str | int]:
    x, log_a, log_n, u_x = sp.symbols("x A L ux", real=True)
    u_n = log_a - log_n
    p_h = -sp.I * (x + u_n) * (x - u_n) ** 2 / 4
    p_t = (
        -sp.I * (x - u_n) * (x + u_n) ** 2 / 4
        + u_x * (x - u_n)
    )
    c_h = log_n * (2 * log_a - log_n) ** 2 / 4
    delta_h = -sp.I * c_h
    kappa = 1 / (2 * sp.pi * sp.I)
    checks = [
        p_h.subs(x, -u_n),
        p_h.subs(x, -log_a) - sp.I * c_h,
        delta_h + p_h.subs(x, -log_a),
        -kappa * delta_h - c_h / (2 * sp.pi),
        kappa * delta_h + c_h / (2 * sp.pi),
        p_t + p_h.subs(x, -x) - u_x * (x - u_n),
    ]
    for index, expression in enumerate(checks):
        require_zero(expression, f"ideal endpoint identity {index}")

    return {
        "ibp": "For Delta_P=f_P(N)-f_P(1), one integration by parts at integer endpoints gives I_P(k)=-kappa Delta_P/k+O(k^(-2)) and I_P(-k)=kappa Delta_P/k+O(k^(-2)), with kappa=1/(2pi i).",
        "criterion": "If Delta_P is nonzero, the positive and negative one-sided tails have opposite nonzero harmonic leading terms. Neither tail converges by itself; only their common symmetric pairing cancels the 1/k term.",
        "hermitian": "For x=log(u/a), u_N=log(a/N), S(0)=0, and P_H^0=-i(x+u_N)(x-u_N)^2/4, one has P_H^0(log N)=0 and P_H^0(0)=i log N{log(a^2/N)}^2/4.",
        "delta": "Thus Delta_H=-i c_H, where c_H=log N{log(a^2/N)}^2/4>0 in the physical regime a>=N>=2.",
        "asymptotic": "Consequently I_H(k)=c_H/(2pi k)+O(k^(-2)) and I_H(-k)=-c_H/(2pi k)+O(k^(-2)). The two separate outer infinities diverge with opposite real signs even though their paired sum is O(k^(-2)).",
        "linear_witness": "The exact nonphysical source f(u)=u on integer [1,N] gives I_f(k)=-kappa(N-1)/k and I_f(-k)=kappa(N-1)/k for every k>=1. This is an exact one-sided-divergence witness, not merely an asymptotic example.",
        "identity_checks": len(checks),
        "one_sided_divergence_witnesses": 1,
    }


def symmetry_certificate() -> dict[str, str | int]:
    mp.mp.dps = 65
    alpha = mp.mpf("2.73")

    def s_value(log_u: mp.mpf) -> mp.mpf:
        return -mp.mpf("0.11") * log_u * log_u

    def p_value(log_u: mp.mpf) -> mp.mpc:
        return 1 + (mp.mpf("1.7") - mp.mpf("0.4") * 1j) * log_u

    conjugation_checks = 0
    for r_value in [-5, -1, 2, 7]:
        original = mp.quad(
            lambda u: (
                mp.e ** (2j * mp.pi * alpha * mp.log(u))
                * mp.e ** s_value(mp.log(u))
                * p_value(mp.log(u))
                * mp.e ** (-2j * mp.pi * r_value * u)
            ),
            [1, 2, 3, 4],
        )
        transformed = mp.quad(
            lambda u: (
                mp.e ** (-2j * mp.pi * alpha * mp.log(u))
                * mp.e ** s_value(mp.log(u))
                * mp.conj(p_value(mp.log(u)))
                * mp.e ** (2j * mp.pi * r_value * u)
            ),
            [1, 2, 3, 4],
        )
        require(abs(mp.conj(original) - transformed) < mp.mpf("1e-50"), "conjugation")
        conjugation_checks += 1

    return {
        "cubic": "The algebraic identity P_T^0(x)=-P_H^0(-x)+u_(N,x)(x-u_N) is exact, but it is only a polynomial identity.",
        "coordinate": "The reflection x->-x sends u to a^2/u because x=log(u/a). It sends [1,N] to [a^2/N,a^2], which is not [1,N]; for a>=N>=2 the image lies at or above N.",
        "phase": "Under u=a^2/v, alpha_Plog u-ru becomes 2alpha_Plog a-alpha_Plog v-r a^2/v. Its v-derivative cannot equal alpha_P/v-s for a constant Fourier mode s on an interval when alpha_P>0.",
        "conjugation": "For real S, conjugate(I_P^(alpha_P)(r))=I_(conjugate P)^(-alpha_P)(-r). Conjugation reverses alpha_P as well as r, so it does not pair two modes inside the fixed physical alpha_P family.",
        "abel": "The Hermitian terminal Abel coefficient still vanishes at q=N, but that finite identity neither changes Delta_H nor regularizes either one-sided outer series.",
        "reflection_identities": 2,
        "conjugation_checks": conjugation_checks,
    }


def fourier_guard_certificate() -> dict[str, str | int]:
    mp.mp.dps = 65
    terminal = 4
    alpha = 6
    m_value = int(mp.floor(alpha / (terminal + mp.mpf("0.5")))) + 1
    n_value = int(mp.floor(2 * alpha))
    require((m_value, n_value) == (2, 12), "Fourier guard roster")
    mode_checks = 0
    for source_mode in [5, 15]:
        for r_value in [0, 1, 5, 12, 13, 15, 18]:
            integral = mp.quad(
                lambda u: mp.e ** (2j * mp.pi * (source_mode - r_value) * u),
                [1, 2, 3, 4],
            )
            expected = terminal - 1 if source_mode == r_value else 0
            require(abs(integral - expected) < mp.mpf("1e-54"), "Fourier mode guard")
            mode_checks += 1

    return {
        "source": "For the nonphysical Fourier source f_K(u)=e(Ku) on integer [1,N], I_K(r)=N-1 when r=K and is zero for every other integer r; the starred carrier is also M=N-1.",
        "inside": "If K lies in mathcal R_N, then mathscr O_N=0 and M+mathscr O_N=N-1.",
        "outside": "If K lies outside mathcal R_N, then mathscr O_N=N-1 and M+mathscr O_N=2(N-1). Negating f_K reverses both values.",
        "scope": "Therefore the outer pullback, finite Abel identity, and Fourier algebra alone provide neither universal cancellation nor a sign. Any useful result must exploit the physical ideal-cubic coefficient and its joined terminal structure.",
        "mode_checks": mode_checks,
        "witnesses": 2,
    }


def handoff_certificate() -> dict[str, str]:
    return {
        "decision": "Retire the requested three-way split into independent negative, low-positive, and high-positive infinite tails. Keep the finite low modes and the absolutely convergent remote pairs as the canonical decomposition of the symmetric outer sum.",
        "exact_target": "The complete ideal targets are mathcal J_H^0=mathscr M_N[P_H^0]+mathscr N_N[P_H^0]+mathscr C_N[P_H^0] and mathcal J_T^0=mathscr M_N[P_T^0]+mathscr N_N[P_T^0]+mathscr C_N[P_T^0]+2iu_Nmathcal T_N^[N][B_(T,N)^0].",
        "next": "Keep K_near,N intact and test a cellwise Euler/Abel composition of the starred carrier with its mean-one finite kernel. The remote paired tail may be attached afterward through its absolute C^2 bound; the transpose terminal survivor remains explicit.",
        "obligation": "Prove or rigorously obstruct a source-specific signed estimate for the carrier-plus-near package, including endpoint half-weights. Do not assign budgets to the divergent one-sided outer tails.",
        "reserve": "No signed complete ideal join, physical derivative-size bound, completed-current theorem, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
    }


def build_rows(certificate: dict) -> list[GateRow]:
    rec = certificate["recombination"]
    cutoff = certificate["cutoff"]
    pair = certificate["paired_tail"]
    endpoint = certificate["endpoint"]
    symmetry = certificate["symmetry"]
    guard = certificate["fourier_guard"]
    handoff = certificate["handoff"]
    rows = [
        GateRow("sop_01_sources", "source chain", "proved", "All four parent gates replay and are hash-pinned.", "Source hashes and parent counts are recorded.", "No parent estimate is strengthened."),
        GateRow("sop_02_poisson", "full-support identity", "proved", "The exterior is the missing Poisson complement.", rec["poisson"], rec["guard"]),
        GateRow("sop_03_exterior", "exterior pullback", "proved", "E_band-tau_band equals the symmetric outer sum.", rec["exterior"], rec["guard"]),
        GateRow("sop_04_hermitian_join", "ideal Hermitian join", "proved", "The Hermitian join is carrier plus outer complement.", rec["hermitian"], rec["guard"]),
        GateRow("sop_05_transpose_join", "ideal transpose join", "proved", "The transpose terminal survivor remains separate.", rec["transpose"], rec["guard"]),
        GateRow("sop_06_roster", "reciprocal roster", "proved", "The complement uses the certified half-open positive band.", cutoff["indices"], "Assume the nonempty physical roster m_N<=n_N."),
        GateRow("sop_07_finite_split", "symmetric cutoff", "proved", "Every finite symmetric cutoff splits without rearrangement.", cutoff["finite_split"], "The common cutoff is part of the theorem."),
        GateRow("sop_08_near", "finite near modes", "proved", "All unpaired modes form one finite object.", cutoff["near"], "No sign is assigned to this finite sum."),
        GateRow("sop_09_remote", "remote pairs", "proved", "Only the remote opposite modes are paired.", cutoff["remote"], "The two one-sided infinities are not defined separately."),
        GateRow("sop_10_near_kernel", "finite kernel", "proved", "The near object has a contiguous Dirichlet kernel.", cutoff["kernel"], "The kernel is generally complex and nonsymmetric."),
        GateRow("sop_11_kernel_mean", "kernel mean", "proved", "The near kernel has unit cell mean one.", cutoff["kernel"], "Mean one is not positivity or pointwise control."),
        GateRow("sop_12_pair_cosine", "paired phase", "proved", "Each remote pair is a cosine integral.", pair["pair"], "No physical derivative bound follows yet."),
        GateRow("sop_13_pair_ibp", "paired integration", "proved", "Two integrations expose k^(-2) decay.", pair["twice_ibp"], pair["scope"]),
        GateRow("sop_14_pair_bound", "absolute paired tail", "proved", "The paired remote series converges absolutely.", pair["bound"], pair["scope"]),
        GateRow("sop_15_scaling_guard", "physical scaling", "guard_validated", "The C^2 numerator remains a physical obligation.", pair["scope"], "No h^2 reserve is claimed."),
        GateRow("sop_16_endpoint_ibp", "one-sided asymptotic", "proved", "Endpoint mismatch creates opposite harmonic terms.", endpoint["ibp"], endpoint["criterion"]),
        GateRow("sop_17_split_guard", "summation guard", "guard_validated", "The separate outer infinities need not converge.", endpoint["criterion"], "They must not receive independent error budgets."),
        GateRow("sop_18_h_terminal", "Hermitian endpoint", "proved", "The upper ideal Hermitian value vanishes.", endpoint["hermitian"], "This does not remove the lower endpoint."),
        GateRow("sop_19_h_lower", "Hermitian lower endpoint", "proved", "The lower ideal Hermitian value is nonzero physically.", endpoint["hermitian"], endpoint["delta"]),
        GateRow("sop_20_h_delta", "Hermitian jump", "proved", "The physical endpoint jump is explicit and nonzero.", endpoint["delta"], "Smooth interior changes do not alter the leading endpoint term."),
        GateRow("sop_21_h_positive", "positive tail", "proved", "The positive Hermitian coefficients have a positive harmonic lead.", endpoint["asymptotic"], "The one-sided series diverges."),
        GateRow("sop_22_h_negative", "negative tail", "proved", "The negative Hermitian coefficients have the opposite harmonic lead.", endpoint["asymptotic"], "The one-sided series diverges."),
        GateRow("sop_23_h_pair", "harmonic cancellation", "proved", "The leading endpoint terms cancel only pairwise.", endpoint["asymptotic"], "No sign follows for the O(k^(-2)) pair."),
        GateRow("sop_24_linear_witness", "exact divergence witness", "guard_validated", "A linear source realizes the harmonic obstruction exactly.", endpoint["linear_witness"], "The witness is nonphysical and only guards the summation method."),
        GateRow("sop_25_cubic_reflection", "polynomial symmetry", "proved", "The ideal cubics obey an algebraic reflection identity.", symmetry["cubic"], "Polynomial symmetry alone is not operator symmetry."),
        GateRow("sop_26_coordinate_reflection", "coordinate map", "guard_validated", "The reflection moves the physical interval.", symmetry["coordinate"], "The transformed domain is not the original domain."),
        GateRow("sop_27_phase_reflection", "phase map", "guard_validated", "The reflection does not preserve a Fourier mode.", symmetry["phase"], "No Hermitian/transpose outer cancellation follows."),
        GateRow("sop_28_conjugation", "conjugation map", "proved", "Conjugation reverses both frequency parameters.", symmetry["conjugation"], "The fixed-alpha physical family is not closed under this map."),
        GateRow("sop_29_abel", "finite Abel guard", "guard_validated", "The terminal Hermitian Abel null remains finite only.", symmetry["abel"], "It does not regularize the outer half-series."),
        GateRow("sop_30_fourier_inside", "inside-band witness", "guard_validated", "A pure inside-band mode gives a nonzero ideal package.", guard["inside"], guard["scope"]),
        GateRow("sop_31_fourier_outside", "outside-band witness", "guard_validated", "A pure outside-band mode changes the package coefficient.", guard["outside"], guard["scope"]),
        GateRow("sop_32_near_open", "near-package estimate", "open", "Control the carrier plus mean-one near kernel without blockwise moduli.", handoff["obligation"], handoff["reserve"]),
        GateRow("sop_33_terminal_open", "transpose completion", "open", "Join the transpose endpoint atom and terminal conditional block.", handoff["next"], handoff["reserve"]),
        GateRow("sop_34_boundary", "proof boundary", "guard_validated", "This gate is not a proof of RH.", handoff["reserve"], "No prize-level conclusion is claimed."),
    ]
    require(len(rows) == 34, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    rec = cert["recombination"]
    cutoff = cert["cutoff"]
    pair = cert["paired_tail"]
    endpoint = cert["endpoint"]
    symmetry = cert["symmetry"]
    guard = cert["fourier_guard"]
    handoff = cert["handoff"]
    return f"""# Ideal-Cubic Symmetric Outer-Pairing Gate

Date: 2026-08-03

Status: exact symmetric outer decomposition and one-sided-divergence obstruction proved; signed carrier-plus-near estimate open; not a proof of RH.

## Exterior Pullback

{rec['poisson']}

{rec['exterior']}

{rec['hermitian']}

{rec['transpose']}

{rec['guard']}

## Canonical Symmetric Split

{cutoff['indices']}

{cutoff['finite_split']}

{cutoff['near']}

{cutoff['kernel']}

{cutoff['remote']}

## Absolutely Convergent Remote Pairs

{pair['pair']}

{pair['twice_ibp']}

{pair['bound']}

{pair['scope']}

## One-Sided Obstruction

{endpoint['ibp']}

{endpoint['criterion']}

{endpoint['hermitian']}

{endpoint['delta']}

{endpoint['asymptotic']}

{endpoint['linear_witness']}

## Symmetry Audit

{symmetry['cubic']}

{symmetry['coordinate']}

{symmetry['phase']}

{symmetry['conjugation']}

{symmetry['abel']}

## Coefficient-Blind Guard

{guard['source']}

{guard['inside']}

{guard['outside']}

{guard['scope']}

## Route Decision

{handoff['decision']}

{handoff['exact_target']}

{handoff['next']}

{handoff['obligation']}

## Pi Provenance

Every pi in this gate comes from the fixed Fourier character e(x)=exp(2pi i x) and kappa=1/(2pi i). The factor 2pi k is the derivative of that character with respect to u. No circle, polygon, fitted constant, or plotted symmetry is inserted.

## Proof Boundary

{handoff['reserve']} This gate is not a proof of RH.
"""


def main() -> int:
    load_sources()
    certificate = {
        "recombination": recombination_certificate(),
        "cutoff": cutoff_certificate(),
        "paired_tail": paired_tail_certificate(),
        "endpoint": endpoint_certificate(),
        "symmetry": symmetry_certificate(),
        "fourier_guard": fourier_guard_certificate(),
        "handoff": handoff_certificate(),
    }
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "recombination_identities": certificate["recombination"]["identities"],
        "finite_roster_checks": certificate["cutoff"]["exact_roster_checks"],
        "character_split_checks": certificate["cutoff"]["character_checks"],
        "near_kernel_mean_checks": certificate["cutoff"]["mean_checks"],
        "remote_pair_formula_checks": certificate["paired_tail"]["pair_formula_checks"],
        "remote_tail_bound_checks": certificate["paired_tail"]["tail_bound_checks"],
        "ideal_endpoint_identities": certificate["endpoint"]["identity_checks"],
        "one_sided_divergence_witnesses": certificate["endpoint"]["one_sided_divergence_witnesses"],
        "reflection_identities": certificate["symmetry"]["reflection_identities"],
        "conjugation_checks": certificate["symmetry"]["conjugation_checks"],
        "fourier_mode_checks": certificate["fourier_guard"]["mode_checks"],
        "fourier_mode_witnesses": certificate["fourier_guard"]["witnesses"],
        "signed_complete_ideal_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-03",
        "status": "symmetric outer decomposition and Hermitian one-sided-divergence obstruction proved; signed carrier-plus-near estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact exterior pullback, canonical common-cutoff decomposition, mean-one finite near kernel, absolute convergence of remote pairs, the physical ideal-Hermitian endpoint obstruction to separate one-sided sums, and reflection/conjugation guards. It proves no signed complete ideal join, physical derivative-size bound, completed current, Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built ideal-cubic symmetric outer-pairing gate: "
        f"{counts['rows']} rows, {counts['recombination_identities']} recombinations, "
        f"{counts['finite_roster_checks']} roster checks, "
        f"{counts['character_split_checks']} character checks, "
        f"{counts['remote_pair_formula_checks']} remote-pair checks, "
        f"{counts['remote_tail_bound_checks']} tail-bound checks, "
        f"{counts['ideal_endpoint_identities']} endpoint identities, "
        f"{counts['one_sided_divergence_witnesses']} one-sided divergence witness, "
        f"{counts['signed_complete_ideal_bounds']} signed complete ideal bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
