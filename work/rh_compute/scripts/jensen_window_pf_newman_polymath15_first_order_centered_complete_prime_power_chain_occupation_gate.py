#!/usr/bin/env python3
"""Build the complete-prime-power-chain occupation gate."""

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
    "complete_prime_power_chain_occupation_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "occupation": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "signed_occupation_transport_reduction.json"
    ),
    "prime_power": (
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
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "direct_projection_regime_reduction.json"
    ),
    "contact_transport": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contact_signed_transport_reduction.json"
    ),
    "adjacent_chart": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_chart_stability_certificate.json"
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
        "occupation": (
            "partial_x bar_mu+partial_s bar_F=s*bar_mu+bar_R",
            "existing direct transversality coordinate",
            "projective infinity",
        ),
        "prime_power": (
            "log(a_(k+1)/a_k)",
            "joined_p_free",
            "p^(-12/25)",
        ),
        "normalized_prefix": (
            "dyadic_completion_scout",
            "G_k=conj(eta)*F_k",
            "first-jet phase flux",
        ),
        "direct_projection": (
            "eta*q_n=r_n*zeta_n",
            "B_a=1/2-t*delta_a/2",
            "The amplitudes are explicit and positive",
        ),
        "contact_transport": (
            "sum_(j=0)^N d_j=mathcal_C_N",
            "cumulative-real-mass",
            "endpoint atom",
        ),
        "adjacent_chart": (
            "rho=exp(Omega)",
            "adjacent_jet",
            "scalar_stability",
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


def chain_factorization_audit() -> dict[str, str]:
    t, log_m, ell, sigma, tau = sp.symbols(
        "t log_m ell sigma tau", real=True
    )
    k = sp.symbols("k", integer=True, nonnegative=True)
    s_star = sigma + sp.I * tau
    direct = (
        t * (log_m + k * ell) ** 2 / 4
        - s_star * (log_m + k * ell)
    )
    factored = (
        t * log_m**2 / 4
        - s_star * log_m
        + k * ell * (t * log_m / 2 - s_star)
        + t * ell**2 * k**2 / 4
    )
    if sp.expand(direct - factored) != 0:
        raise RuntimeError("prime-power chain factorization failed")

    c_s, b, kappa_0, u_0 = sp.symbols(
        "c_s b kappa_0 u_0", real=True
    )
    z_re, z_im, K_re, K_im = sp.symbols(
        "z_re z_im K_re K_im", real=True
    )
    direct_moment = (
        c_s * (kappa_0 * z_re - ell * K_re)
        - b * (u_0 * z_im - ell * K_im)
    )
    complex_moment = sp.re(
        (c_s * kappa_0 + sp.I * b * u_0)
        * (z_re + sp.I * z_im)
        - ell
        * (c_s + sp.I * b)
        * (K_re + sp.I * K_im)
    ).expand(complex=True)
    if sp.simplify(direct_moment - complex_moment) != 0:
        raise RuntimeError("chain Euler-moment identity failed")

    return {
        "coefficient_family": (
            "The active coefficients are "
            "q_n=exp[t*log(n)^2/4-s_*log(n)]"
            "*(1+d_n)/(1+d_1), q_1=1. They carry coefficient +1; "
            "there is no Mobius factor in this Xi prefix."
        ),
        "exact_factorization": (
            "Fix p prime, p does not divide m, ell=log(p), "
            "n_k=m*p^k, and R_m=floor(log_p(N/m)). With "
            "B_m=phi*exp[t*log(m)^2/4-s_*log(m)]/"
            "|1+d_1|, w_(m,p)=exp[ell*(t*log(m)/2-s_*)], and "
            "a_(m,k)=exp[t*ell^2*k^2/4]*(1+d_(m*p^k)), one has "
            "z_(m*p^k)=eta*q_(m*p^k)=B_m*a_(m,k)*w_(m,p)^k "
            "exactly."
        ),
        "chain_polynomial": (
            "Q_(m,p)(w)=sum_(k=0)^R_m a_(m,k)w^k, "
            "Z_(m,p)=B_m*Q_(m,p)(w_(m,p)), and "
            "K_(m,p)=B_m*w_(m,p)*Q_(m,p)'(w_(m,p))"
            "=sum_k k*z_(m*p^k)."
        ),
        "chain_moment": (
            "Put kappa_0=log(N/m), u_0=log(a/m), "
            "C_(m,p)=Re Z_(m,p), and s_*'=c_s+i*b. Then "
            "D_(m,p)=Re[(c_s*kappa_0+i*b*u_0)Z_(m,p)"
            "-ell*s_*'*K_(m,p)]."
        ),
        "chain_transport": (
            "On every fixed-N chart and fixed chain, "
            "D_(m,p)=partial_x C_(m,p)+c_s*log(N)*C_(m,p)"
            "-R_(m,p), where R_(m,p)=sum_k r_(m*p^k). "
            "Thus a chain contact still reduces to transversality, "
            "not positivity."
        ),
    }


def threshold_guard_audit() -> dict[str, str | dict]:
    root_two = sp.sqrt(2)
    ell = sp.symbols("ell", positive=True)
    theta = sp.pi / 3
    ratio = 1 / root_two
    masses = [
        sp.cos(theta),
        ratio * sp.cos(2 * theta),
        ratio**2 * sp.cos(3 * theta),
    ]
    slopes = [
        sp.Rational(5, 2) * ell * sp.tan(theta),
        sp.Rational(3, 2) * ell * sp.tan(2 * theta),
        sp.Rational(1, 2) * ell * sp.tan(3 * theta),
    ]
    expected_masses = [
        sp.Rational(1, 2),
        -root_two / 4,
        -sp.Rational(1, 2),
    ]
    expected_slopes = [
        5 * sp.sqrt(3) * ell / 2,
        -3 * sp.sqrt(3) * ell / 2,
        0,
    ]
    if any(
        sp.simplify(left - right) != 0
        for left, right in zip(masses, expected_masses)
    ):
        raise RuntimeError("dyadic threshold masses failed")
    if any(
        sp.simplify(left - right) != 0
        for left, right in zip(slopes, expected_slopes)
    ):
        raise RuntimeError("dyadic threshold slopes failed")
    profile = [
        sp.simplify(sum(masses)),
        sp.simplify(masses[0] + masses[2]),
        sp.simplify(masses[0]),
    ]
    if profile != [-root_two / 4, 0, sp.Rational(1, 2)]:
        raise RuntimeError("dyadic threshold profile failed")

    first_moment = sp.simplify(
        sum(mass * slope for mass, slope in zip(masses, slopes))
    )
    if first_moment != ell * (10 * sp.sqrt(3) + 3 * sp.sqrt(6)) / 8:
        raise RuntimeError("dyadic threshold first moment failed")

    return {
        "general_construction": (
            "Let A_0,A_1,A_2>0 with A_2<A_0, put "
            "theta_0=arccos(A_2/A_0), "
            "varphi=(pi-theta_0)/2, and "
            "theta_k=theta_0+k*varphi. Then theta_2=pi, "
            "c_0=A_2, c_2=-A_2, and c_1<0."
        ),
        "slope_order": (
            "For u_0>2ell, u_k=u_0-k*ell, "
            "kappa_k=kappa_0-k*ell, beta=-b>=1/2, and "
            "0<=c_s<1/4, h_k=c_s*kappa_k+beta*u_k*tan(theta_k) "
            "satisfies h_1<h_2<h_0. Indeed "
            "h_1-h_2<c_s*ell-beta*u_1<ell/4-ell/2<0, "
            "while h_0-h_2=2c_s*ell+beta*u_0*tan(theta_0)>0."
        ),
        "threshold_profile": (
            "For P_(m,p)(s)=sum_(k=0)^2 c_k*1_(s<h_k), "
            "P=c_1<0 below h_1, P=c_0+c_2=0 on (h_1,h_2), "
            "and P=c_0=A_2>0 on (h_2,h_0). Thus neither a fixed "
            "sign nor a strict positive lower modulus follows from "
            "complete-chain amplitudes and affine phases."
        ),
        "external_phase_flip": (
            "Adding pi to every theta_k leaves every h_k unchanged "
            "but sends every c_k, P_(m,p), C_(m,p), and D_(m,p) "
            "to its negative. Any composable sign theorem must "
            "control the p-free external phase B_m; internal chain "
            "contraction cannot do so."
        ),
        "correction_guard": (
            "For fixed correction phases delta_k, the two affine "
            "parameters alpha,varphi can still set "
            "alpha+delta_0=theta_0 and "
            "alpha+2varphi+delta_2=pi. The middle phase is "
            "(theta_0+pi)/2+delta_1-(delta_0+delta_2)/2. "
            "Small d_n corrections therefore do not create a "
            "coefficient-only sign; this does not assert that the "
            "actual Xi trajectory realizes the guard phases."
        ),
        "exact_dyadic_control": {
            "parameters": (
                "A=(1,2^(-1/2),1/2), "
                "theta=(pi/3,2pi/3,pi), "
                "u=(5ell/2,3ell/2,ell/2), c_s=0, b=-1"
            ),
            "slope_order": "h_1=-3sqrt(3)ell/2<h_2=0<h_0=5sqrt(3)ell/2",
            "ordered_profile": (
                "-sqrt(2)/4 below h_1; 0 on (h_1,h_2); "
                "1/2 on (h_2,h_0)"
            ),
            "first_moment": "ell*(10sqrt(3)+3sqrt(6))/8>0",
            "status": (
                "Exact correction-free complete three-level chain "
                "guard, not an actual Xi point."
            ),
        },
        "pi_provenance": (
            "The pi in the Xi saddle and a^2=x/(4*pi)+t/16 comes "
            "from the completed-zeta Gamma normalization. The pi "
            "in theta_2=pi is separately the half-period of "
            "exp(i*theta); no polygon or prime chain defines pi."
        ),
    }


def outer_decomposition_audit() -> dict[str, str | dict]:
    primes = (2, 3, 5, 7)
    for prime in primes:
        for n_max in range(1, 2049):
            middle = n_max // prime
            formula = (
                n_max
                - 2 * middle
                + middle // prime
            )
            actual = sum(
                1
                for base in range(middle + 1, n_max + 1)
                if base % prime != 0
            )
            if formula != actual:
                raise RuntimeError("singleton-chain count failed")
    for n_max in range(1, 2049):
        middle = n_max // 2
        count = n_max - 2 * middle + middle // 2
        if count < n_max // 4:
            raise RuntimeError("dyadic singleton lower bound failed")

    q = sp.symbols("q")
    prefix = sp.symbols("S_0:7")
    telescoping = sp.expand(
        sum(
            q**k * (prefix[k] - q * prefix[k + 1])
            for k in range(6)
        )
    )
    if telescoping != prefix[0] - q**6 * prefix[6]:
        raise RuntimeError("p-free prefix telescoping failed")

    return {
        "exact_partition": (
            "For every fixed prime p, each 1<=n<=N has the unique "
            "form n=m*p^k with p not dividing m. Therefore "
            "mu=sum_(p not dividing m)mu_(m,p), "
            "P=sum_m P_(m,p), C_bulk=sum_m C_(m,p), and "
            "D_bulk=sum_m D_(m,p) exactly."
        ),
        "singleton_count": (
            "Let M=floor(N/p). The number of p-free complete chains "
            "with R_m=0 is S_p(N)=N-2M+floor(M/p). For p=2, "
            "S_2(N)>=floor(N/4). Thus a fixed-prime chain split "
            "leaves order-N singleton chains with no internal "
            "completion at all."
        ),
        "p_free_prefix": (
            "For the correction-free Dirichlet prefix "
            "S_M(s)=sum_(n<=M)n^(-s), "
            "O_M^(p)(s)=S_M(s)-p^(-s)S_floor(M/p)(s). If "
            "R=floor(log_p N), then "
            "sum_(k=0)^R p^(-ks)O_floor(N/p^k)^(p)(s)=S_N(s) "
            "by exact telescoping."
        ),
        "no_mobius_reinterpretation": (
            "The subtraction in O_M^(p)=S_M-p^(-s)S_floor(M/p) "
            "is a one-prime p-free inclusion identity. It does not "
            "turn the original carrier coefficients into Mobius "
            "coefficients and supplies no positivity."
        ),
        "heat_rejoin": (
            "The Gaussian identity writes the heat factor as an "
            "expectation of tilted Dirichlet prefixes. The p-free "
            "decomposition telescopes for every tilt before taking "
            "the expectation. This exactly reconstructs the original "
            "prefix; expectation cannot be passed through the "
            "phase-dependent indicators 1_(s<h_n)."
        ),
        "singleton_diagnostics": {
            "checked_primes": list(primes),
            "checked_N_range": "1<=N<=2048",
            "dyadic_lower_bound": "S_2(N)>=floor(N/4)",
        },
    }


def build_rows(
    chain: dict[str, str],
    guard: dict[str, str | dict],
    outer: dict[str, str | dict],
) -> list[GateRow]:
    return [
        GateRow(
            "cppco_01_coefficient_family",
            "source contract",
            "proved",
            "The active Xi prefix has coefficient +1 and no Mobius factor.",
            chain["coefficient_family"],
            "A p-free inclusion subtraction must not be relabelled as a Mobius coefficient.",
        ),
        GateRow(
            "cppco_02_chain_factorization",
            "exact algebra",
            "proved",
            "Every fixed prime-power chain has an exact heat-quadratic polynomial factorization.",
            chain["exact_factorization"],
            "The complex correction factors remain in a_(m,k).",
        ),
        GateRow(
            "cppco_03_chain_polynomial",
            "exact algebra",
            "proved",
            "The chain sum and logarithmic index moment are Q and its Euler derivative.",
            chain["chain_polynomial"],
            "This is a fixed-x identity, not a zero-free theorem.",
        ),
        GateRow(
            "cppco_04_chain_moment",
            "exact algebra",
            "proved",
            "The physical centered first moment has an exact chain-Euler representation.",
            chain["chain_moment"],
            "It is a real projection and has no automatic sign.",
        ),
        GateRow(
            "cppco_05_chain_transport",
            "exact transport",
            "proved",
            "Each fixed chain inherits the normalized occupation transversality identity.",
            chain["chain_transport"],
            "At chain contact this remains a derivative statement.",
        ),
        GateRow(
            "cppco_06_general_guard",
            "nonpromotion guard",
            "proved",
            "Arbitrary contracted endpoint amplitudes admit a complete three-level affine-phase threshold guard.",
            str(guard["general_construction"]),
            "The phase parameters are not asserted to be attained by the Xi trajectory.",
        ),
        GateRow(
            "cppco_07_slope_order",
            "nonpromotion guard",
            "proved",
            "The three guard slopes are strictly ordered in the physical c_s,b range.",
            str(guard["slope_order"]),
            "The guard is pole-free and uses u_0>2log(p).",
        ),
        GateRow(
            "cppco_08_threshold_profile",
            "route obstruction",
            "proved",
            "A complete chain threshold profile can be negative, zero on an open interval, and positive.",
            str(guard["threshold_profile"]),
            "This disproves chain-completion positivity, not an Xi-specific joined-prefix theorem.",
            guard["exact_dyadic_control"],
        ),
        GateRow(
            "cppco_09_external_phase",
            "route obstruction",
            "proved",
            "A common half-turn reverses every projected chain sign without changing its slopes.",
            str(guard["external_phase_flip"]),
            "The p-free external phase must be controlled jointly.",
        ),
        GateRow(
            "cppco_10_correction_guard",
            "nonpromotion guard",
            "proved",
            "Fixed small coefficient phase corrections do not create coefficient-only threshold positivity.",
            str(guard["correction_guard"]),
            "No actual Xi phase-attainment statement is made.",
        ),
        GateRow(
            "cppco_11_pi_provenance",
            "normalization provenance",
            "proved",
            "Every occurrence of pi has an explicit source.",
            str(guard["pi_provenance"]),
            "The chain does not define or approximate pi.",
        ),
        GateRow(
            "cppco_12_exact_partition",
            "exact arithmetic",
            "proved",
            "The full occupation measure is the sum of its p-free prime-power chains.",
            str(outer["exact_partition"]),
            "Linearity of the measure does not imply positivity of its summands.",
        ),
        GateRow(
            "cppco_13_singletons",
            "route obstruction",
            "proved",
            "Every fixed-prime split retains order-N complete singleton chains.",
            str(outer["singleton_count"]),
            "Internal prime-power rigidity has no content on these chains.",
            outer["singleton_diagnostics"],
        ),
        GateRow(
            "cppco_14_p_free_telescoping",
            "exact arithmetic",
            "proved",
            "The correction-free p-free rejoin telescopes to the original prefix.",
            str(outer["p_free_prefix"]),
            "The decomposition alone creates no new coercive quantity.",
        ),
        GateRow(
            "cppco_15_no_mobius",
            "source contract",
            "proved",
            "The p-free subtraction is not a Mobius-weighted Xi expansion.",
            str(outer["no_mobius_reinterpretation"]),
            "No alternating-sign theorem may be imported from a different coefficient family.",
        ),
        GateRow(
            "cppco_16_heat_rejoin",
            "exact arithmetic",
            "proved",
            "Gaussian heat mixing reconstructs the same joined prefix and does not linearize threshold indicators.",
            str(outer["heat_rejoin"]),
            "Expectation preserves the sum identity, not zero-free winding or occupation sign.",
        ),
        GateRow(
            "cppco_17_global_contact",
            "endpoint contract",
            "proved",
            "The contact condition is global rather than chainwise.",
            (
                "With the recurrent endpoint retained, "
                "mathsf_X=c_0+sum_m C_(m,p) and "
                "mathcal_C_N=d_0+sum_m D_(m,p). "
                "mathsf_X=0 does not imply C_(m,p)=0 for any m."
            ),
            "A chain-contact estimate cannot be summed under only the global contact equation.",
        ),
        GateRow(
            "cppco_18_cutoff_join",
            "chart contract",
            "proved",
            "Prime-power chains do not evolve independently across cutoff changes.",
            (
                "When N changes, chain lengths and every "
                "kappa_k=log(N/(m*p^k)) change. Recompute terminal "
                "centering and use the exact global adjacent laws "
                "F_(N+1)-F_N=j_0+u_nz^v and "
                "Q_N=f_(N+1)+kappa_N*J_a^(adj)."
            ),
            "No blockwise cutoff homotopy replaces the recurrent endpoint join.",
        ),
        GateRow(
            "cppco_19_route_decision",
            "route decision",
            "closed for chainwise positivity",
            "Retire complete-chain threshold positivity as a standalone route.",
            (
                "Complete chains remain an exact organization of the "
                "prefix, but their projected threshold masses have no "
                "chainwise sign and a fixed-p split leaves order-N "
                "singletons. Return to the endpoint-complete joined "
                "first-jet boundary flux, where p-free external phases "
                "and the recurrent endpoint are controlled together."
            ),
            "A stronger source-specific joined-prefix theorem remains logically possible.",
        ),
        GateRow(
            "cppco_20_next_target",
            "open theorem",
            "open",
            "The next target is a whole-boundary joined-prefix first-jet estimate.",
            (
                "Insert the exact O(N) normalized first-jet flux into "
                "the cofinal boundary-degree contract. Seek strict "
                "boundary domination or a winding cap for the full "
                "joined prefix plus endpoint, beginning on the first "
                "unresolved compact boundary beyond Q_207."
            ),
            "One finite pilot is diagnostic and cannot establish the cofinal theorem.",
        ),
        GateRow(
            "cppco_21_proof_boundary",
            "proof boundary",
            "audited",
            "The chain gate proves an exact obstruction, not RH.",
            (
                "Proved: exact physical chain factorization, Euler "
                "moment, chain transport, three-level threshold guard, "
                "external-phase reversal, p-free partition, singleton "
                "count, and telescoping rejoin. Open: a joined Xi "
                "boundary margin, winding cap, q<1 chart, finite-height "
                "cofinal closure, contact exclusion, Lambda<=0, "
                "PF-infinity, RH, and any Clay-prize conclusion."
            ),
            "No generic guard is labelled an actual Xi counterexample.",
        ),
    ]


def build_payload() -> dict:
    source_kinds = source_audit()
    chain = chain_factorization_audit()
    guard = threshold_guard_audit()
    outer = outer_decomposition_audit()
    rows = build_rows(chain, guard, outer)
    return {
        "kind": "complete_prime_power_chain_occupation_gate",
        "date": "2026-07-28",
        "status": "exact_reduction_and_route_obstruction_not_an_rh_proof",
        "source_hashes": source_hashes(),
        "source_kinds": source_kinds,
        "chain": chain,
        "threshold_guard": guard,
        "outer_decomposition": outer,
        "rows": [asdict(row) for row in rows],
        "summary": {
            "row_count": len(rows),
            "exact_chain_factorizations": 1,
            "exact_euler_moments": 1,
            "threshold_sign_obstructions": 1,
            "external_phase_guards": 1,
            "singleton_count_theorems": 1,
            "p_free_telescoping_identities": 1,
            "open_joined_boundary_theorems": 1,
        },
        "proof_boundary": rows[-1].formula,
    }


def render_note(payload: dict) -> str:
    chain = payload["chain"]
    guard = payload["threshold_guard"]
    outer = payload["outer_decomposition"]
    control = guard["exact_dyadic_control"]
    lines = [
        "# Complete Prime-Power-Chain Occupation Gate",
        "",
        "Date: 2026-07-28",
        "",
        "Status: exact chain reduction and route obstruction. This is not a proof of RH.",
        "",
        "```text",
        f"work/rh_compute/results/{STEM}.json",
        f"python work/rh_compute/scripts/{STEM}.py",
        f"python work/rh_compute/scripts/check_{STEM}.py",
        "```",
        "",
        "## Coefficient Contract",
        "",
        chain["coefficient_family"],
        "",
        "## Exact Chain Algebra",
        "",
        "```text",
        chain["exact_factorization"],
        chain["chain_polynomial"],
        chain["chain_moment"],
        chain["chain_transport"],
        "```",
        "",
        "## Complete-Chain Threshold Guard",
        "",
        guard["general_construction"],
        "",
        guard["slope_order"],
        "",
        guard["threshold_profile"],
        "",
        "```text",
        control["parameters"],
        control["slope_order"],
        control["ordered_profile"],
        control["first_moment"],
        "```",
        "",
        guard["external_phase_flip"],
        "",
        guard["correction_guard"],
        "",
        "## Pi Provenance",
        "",
        guard["pi_provenance"],
        "",
        "## P-Free Rejoining",
        "",
        "```text",
        outer["exact_partition"],
        outer["singleton_count"],
        outer["p_free_prefix"],
        "```",
        "",
        outer["no_mobius_reinterpretation"],
        "",
        outer["heat_rejoin"],
        "",
        "## Route Decision",
        "",
        (
            "Complete prime-power chains remain a useful exact grouping, "
            "but chainwise threshold positivity is false. The p-free "
            "external phases, order-N singleton chains, sharp cutoff, "
            "and recurrent endpoint must be controlled together. The "
            "next target is the endpoint-complete joined-prefix first-jet "
            "boundary flux, beginning with the first unresolved compact "
            "boundary beyond Q_207."
        ),
        "",
        "## Boundary",
        "",
        payload["proof_boundary"],
        "",
    ]
    return "\n".join(lines)


def write_payload(payload: dict, output: Path, note: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    note.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    note.write_text(render_note(payload), encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    write_payload(payload, args.output, args.note)
    summary = payload["summary"]
    print(
        "built complete prime-power-chain occupation gate: "
        f"{summary['row_count']} rows, exact chain factorization "
        "and Euler moment, negative-zero-positive threshold guard, "
        "external-phase reversal, order-N singleton theorem, "
        "p-free telescoping, joined-boundary handoff"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
