#!/usr/bin/env python3
"""Build the six-moment finite-Poisson and stationary-transport reduction."""

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
    "contiguous_terminal_tail_anchored_six_moment_finite_poisson_"
    "transport_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
VANDEHEY_URL = "https://arxiv.org/abs/1205.0090"
SOURCE_PATHS = {
    "flow_matrix_phase": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_flow_matrix_"
        "phase_reduction.json"
    ),
    "saddle_flow": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_saddle_flow_"
        "reduction.json"
    ),
    "hermitian_pairing": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_six_moment_hermitian_"
        "reciprocal_pairing_reduction.json"
    ),
    "signed_reciprocal_symbol": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "signed_handoff_reciprocal_symbol_guard.json"
    ),
}


@dataclass(frozen=True)
class PoissonRow:
    id: str
    role: str
    readiness: str
    claim: str
    formula: str
    proof_boundary: str
    diagnostics: dict | None = None


def file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def require_zero(expression: sp.Expr, label: str) -> None:
    if sp.simplify(sp.expand_complex(expression)) != 0:
        raise RuntimeError(label)


def load_sources() -> dict[str, dict]:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        key: json.loads(path.read_text(encoding="utf-8"))
        for key, path in SOURCE_PATHS.items()
    }


def source_audit(payloads: dict[str, dict]) -> dict:
    flow = payloads["flow_matrix_phase"]
    counts = flow.get("counts", {})
    if counts.get("correction_free_moments") != 6:
        raise RuntimeError("six-moment source drifted")
    if counts.get("real_flow_observations") != 8:
        raise RuntimeError("eight-observation source drifted")
    if "Ghat_j'(xi)=iGhat_(j+1)(xi)" not in flow.get(
        "common_phase_certificate", {}
    ).get("normalized_moments", ""):
        raise RuntimeError("normalized-moment source drifted")

    saddle = payloads["saddle_flow"].get("counts", {})
    if saddle.get("max_moment_order") != 5:
        raise RuntimeError("order-five source drifted")

    pairing = payloads["hermitian_pairing"].get("counts", {})
    if pairing.get("hermitian_swap_pair_identities") != 1:
        raise RuntimeError("Hermitian pairing source drifted")
    if pairing.get("poisson_remainder_bounds") != 0:
        raise RuntimeError("pairing proof boundary drifted")

    reciprocal = payloads["signed_reciprocal_symbol"].get("exact", {})
    if "pi/4" not in reciprocal.get("pi_provenance", ""):
        raise RuntimeError("negative-curvature phase source drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "primary_analytic_source": {
            "author": "Joseph Vandehey",
            "title": "Error term improvements for van der Corput transforms",
            "url": VANDEHEY_URL,
            "specialized_items": [
                "full starred Poisson formula (paper equation (10))",
                "weighted B-process transform (paper Theorem 1.1)",
                "endpoint-sensitive transform discussion and Corollary 1.4",
            ],
        },
        "imported_explicit_numeric_remainder_constants": 0,
        "imported_signed_flow_bounds": 0,
    }


def _amplitude(order: int, logarithm: sp.Expr, t: sp.Expr, sigma: sp.Expr) -> sp.Expr:
    return logarithm**order * sp.exp(t * logarithm**2 / 4 - sigma * logarithm)


def moment_amplitude_certificate() -> dict:
    logarithm, t, sigma = sp.symbols("lambda t sigma", real=True)
    differences: dict[str, str] = {}
    for order in range(6):
        amplitude = _amplitude(order, logarithm, t, sigma)
        previous = 0 if order == 0 else _amplitude(order - 1, logarithm, t, sigma)
        expected = (
            order * previous
            + (t * logarithm / 2 - sigma) * amplitude
        )
        difference = sp.simplify(sp.diff(amplitude, logarithm) - expected)
        require_zero(difference, f"amplitude recurrence failed at order {order}")
        differences[str(order)] = str(difference)

    lower_values = [sp.simplify(_amplitude(order, 0, t, sigma)) for order in range(6)]
    if lower_values != [1, 0, 0, 0, 0, 0]:
        raise RuntimeError("lower endpoint moments drifted")

    return {
        "bare_moments": (
            "For alpha=xi/(2*pi), define "
            "A_j(u)=(log u)^j exp[t(log u)^2/4-sigma log u] and "
            "F_(j,xi)(u)=A_j(u)e(alpha log u). Then the six bare sums "
            "H_j(xi)=sum_(1<=n<=B)F_(j,xi)(n), 0<=j<=5, differ from "
            "Ghat_j only by the fixed normalized source factor."
        ),
        "logarithmic_recurrence": (
            "Writing lambda=log u, partial_lambda A_j="
            "jA_(j-1)+(t*lambda/2-sigma)A_j, with A_(-1)=0."
        ),
        "frequency_shift": (
            "Since e(alpha log u)=exp(i*xi*log u), "
            "partial_xi F_(j,xi)=iF_(j+1,xi) and therefore "
            "H_j'(xi)=iH_(j+1)(xi) for 0<=j<=4."
        ),
        "lower_endpoint": (
            "At u=1, A_0(1)=1 and A_j(1)=0 for 1<=j<=5. Thus the "
            "lower hard endpoint occurs only in the zeroth moment."
        ),
        "symbolic_differences": differences,
        "lower_endpoint_values": [str(value) for value in lower_values],
    }


def finite_poisson_certificate() -> dict:
    alpha, t, sigma, b = sp.symbols(
        "alpha t sigma B", positive=True, real=True
    )
    log_b = sp.log(b)
    phase_b = sp.exp(2 * sp.pi * sp.I * alpha * log_b)
    endpoint_differences: dict[str, str] = {}
    for order in range(5):
        lower = sp.Integer(1) if order == 0 else sp.Integer(0)
        lower_next = sp.Integer(0)
        endpoint = (lower + _amplitude(order, log_b, t, sigma) * phase_b) / 2
        endpoint_next = (
            lower_next + _amplitude(order + 1, log_b, t, sigma) * phase_b
        ) / 2
        difference = sp.simplify(
            sp.diff(endpoint, alpha) / (2 * sp.pi) - sp.I * endpoint_next
        )
        require_zero(difference, f"endpoint shift failed at order {order}")
        endpoint_differences[str(order)] = str(difference)

    return {
        "starred_source_formula": (
            "For a C2 amplitude on [1,B], full Poisson summation in the "
            "symmetric limit gives sum^*_(1<=n<=B)F(n)="
            "lim_(R->infinity)sum_(r=-R)^R integral_1^B "
            "F(u)e(-ru)du. The star gives half weight to an integral "
            "primal endpoint. This is equation (10) in Vandehey's paper."
        ),
        "full_endpoint_formula": (
            "The physical bulk sum uses both endpoints with full weight. "
            "Therefore H_j=E_j+lim_(R->infinity)sum_(r=-R)^R I_(j,r), "
            "where I_(j,r)=integral_1^B A_j(u)e(alpha log u-ru)du and "
            "E_j=[F_(j,xi)(1)+F_(j,xi)(B)]/2."
        ),
        "explicit_endpoint": (
            "E_j=(1/2){delta_(j,0)+A_j(B)e(alpha log B)}. In particular, "
            "E_j'(xi)=iE_(j+1)(xi) for 0<=j<=4."
        ),
        "modewise_shift": (
            "For every integer r, differentiation under the finite integral "
            "gives I_(j,r)'(xi)=iI_(j+1,r)(xi). Applying the sourced Poisson "
            "identity separately at orders j and j+1 avoids an unsupported "
            "termwise differentiation of a conditionally convergent series."
        ),
        "endpoint_symbolic_differences": endpoint_differences,
        "convergence_convention": (
            "The mode sum is a symmetric R->infinity limit. It is not "
            "silently replaced by an absolutely convergent sum."
        ),
    }


def stationary_phase_certificate() -> dict:
    alpha, u, mode = sp.symbols("alpha u r", positive=True, real=True)
    phase = alpha * sp.log(u) - mode * u
    saddle = alpha / mode
    require_zero(
        sp.diff(phase, u).subs(u, saddle),
        "stationary point failed",
    )
    curvature = sp.simplify(sp.diff(phase, u, 2).subs(u, saddle))
    third = sp.simplify(sp.diff(phase, u, 3).subs(u, saddle))
    fourth = sp.simplify(sp.diff(phase, u, 4).subs(u, saddle))
    saddle_phase = sp.expand_log(phase.subs(u, saddle), force=True)
    expected_phase = sp.expand_log(
        alpha * (sp.log(alpha / mode) - 1), force=True
    )
    require_zero(saddle_phase - expected_phase, "stationary phase failed")
    prefactor = sp.sqrt(alpha) / mode
    require_zero(
        prefactor**2 * (-curvature) - 1,
        "stationary prefactor failed",
    )

    return {
        "poisson_phase": (
            "For positive mode r, phi_r(u)=alpha log u-ru has the unique "
            "stationary point u_r=alpha/r. Its curvature is negative. "
            "Modes r<=0 have no stationary point on [1,B]."
        ),
        "active_window": (
            "The saddle lies in [1,B] exactly when alpha/B<=r<=alpha."
        ),
        "stationary_data": (
            "At u_r=alpha/r, phi_r=alpha[log(alpha/r)-1], "
            "phi_r''=-r^2/alpha, and 1/sqrt(|phi_r''|)="
            "sqrt(alpha)/r=u_r/sqrt(alpha)."
        ),
        "signature": (
            "Negative curvature contributes e(-1/8)=exp(-i*pi/4). This "
            "restores the stationary signature suppressed when Section "
            "11.165 only located the saddle."
        ),
        "stationary_main": (
            "The order-j interior main is M_(j,r)="
            "[sqrt(alpha)/r]A_j(u_r)e(alpha[log(u_r)-1]-1/8)."
        ),
        "symbolic": {
            "saddle": str(saddle),
            "curvature": str(curvature),
            "third_derivative": str(third),
            "fourth_derivative": str(fourth),
            "saddle_phase": str(saddle_phase),
            "stationary_prefactor": str(prefactor),
        },
    }


def stationary_transport_certificate() -> dict:
    alpha, mode, t, sigma = sp.symbols(
        "alpha r t sigma", positive=True, real=True
    )
    logarithm = sp.log(alpha / mode)
    prefactor = sp.sqrt(alpha) / mode
    phase = sp.exp(
        2
        * sp.pi
        * sp.I
        * (alpha * (logarithm - 1) - sp.Rational(1, 8))
    )
    differences: dict[str, str] = {}
    transports: dict[str, str] = {}
    for order in range(5):
        amplitude = _amplitude(order, logarithm, t, sigma)
        next_amplitude = _amplitude(order + 1, logarithm, t, sigma)
        previous = 0 if order == 0 else _amplitude(order - 1, logarithm, t, sigma)
        main = prefactor * amplitude * phase
        next_main = prefactor * next_amplitude * phase
        transport = (
            prefactor
            * phase
            / (2 * sp.pi * alpha)
            * (
                order * previous
                + (sp.Rational(1, 2) - sigma + t * logarithm / 2)
                * amplitude
            )
        )
        difference = sp.simplify(
            sp.diff(main, alpha) / (2 * sp.pi)
            - sp.I * next_main
            - transport
        )
        require_zero(difference, f"stationary transport failed at order {order}")
        differences[str(order)] = str(difference)
        transports[str(order)] = str(sp.simplify(transport))

    i_j, i_next, m_j, m_next, correction = sp.symbols(
        "I_j I_next M_j M_next C_j"
    )
    residual_derivative = sp.I * i_next - (sp.I * m_next + correction)
    expected_residual = sp.I * (i_next - m_next) - correction
    require_zero(
        residual_derivative - expected_residual,
        "stationary residual recurrence failed",
    )

    return {
        "main_transport": (
            "On a fixed dual chart, M_(j,r)'(xi)=iM_(j+1,r)+C_(j,r), "
            "where C_(j,r)=[sqrt(alpha)/r]e(alpha[log(u_r)-1]-1/8)"
            "/(2*pi*alpha) times "
            "{jA_(j-1)+(1/2-sigma+(t/2)log(u_r))A_j}."
        ),
        "transport_interpretation": (
            "The iM_(j+1,r) term differentiates the phase. C_(j,r) is the "
            "saddle-motion and stationary-prefactor transport. Dropping it "
            "is not a valid first-frequency derivative."
        ),
        "mode_residual_recurrence": (
            "For R_(j,r)=I_(j,r)-M_(j,r), the exact fixed-chart recurrence "
            "is R_(j,r)'=iR_(j+1,r)-C_(j,r), 0<=j<=4."
        ),
        "six_moment_closure": (
            "The first-frequency remainder through order four requires only "
            "value remainders R_0,...,R_5 and transport terms involving "
            "A_0,...,A_4. No seventh logarithmic moment is introduced."
        ),
        "symbolic_differences": differences,
        "symbolic_transport_terms": transports,
        "residual_recurrence_difference": "0",
    }


def cutoff_chart_certificate() -> dict:
    return {
        "dual_star_sum": (
            "Define P_j=E_j+sum^*_(alpha/B<=r<=alpha)M_(j,r). The dual "
            "star gives half weight when alpha/B or alpha is an integer, "
            "matching the standard B-process endpoint convention."
        ),
        "fixed_chart_remainder": (
            "With Delta_j=H_j-P_j, on every open xi interval where neither "
            "alpha nor alpha/B crosses an integer, "
            "Delta_j'=iDelta_(j+1)-sum^* C_(j,r), 0<=j<=4."
        ),
        "upper_crossing": (
            "As alpha increases through an integer q, the upper condition "
            "r<=alpha admits mode q. The one-sided dual-main jump is "
            "+M_(j,q); the starred value at alpha=q is its midpoint, and "
            "Delta_j has the opposite jump."
        ),
        "lower_crossing": (
            "As alpha increases through Bq, the lower condition r>=alpha/B "
            "removes mode q. The one-sided dual-main jump is -M_(j,q); the "
            "starred value is the midpoint, and Delta_j has the opposite "
            "jump."
        ),
        "c1_guard": (
            "The sharp starred main is piecewise smooth, not globally C1. "
            "A proof across a crossing must either retain these exact jumps "
            "with the hard endpoint/adjacent recurrence or replace the sharp "
            "main by a sourced uniform endpoint transition."
        ),
        "crossing_frequencies": (
            "Since xi=2*pi*alpha, upper crossings occur at xi=2*pi*q and "
            "lower crossings at xi=2*pi*Bq. Their possible absence on the "
            "short auxiliary segment may not be assumed from interval length."
        ),
    }


def dyadic_b_process_certificate() -> dict:
    alpha, u, block, b = sp.symbols(
        "alpha u U B", positive=True, real=True
    )
    phase = alpha * sp.log(u)
    second = sp.diff(phase, u, 2)
    third = sp.diff(phase, u, 3)
    fourth = sp.diff(phase, u, 4)
    global_ratio = sp.simplify(
        sp.Abs(second.subs(u, 1)) / sp.Abs(second.subs(u, b))
    )
    if global_ratio != b**2:
        raise RuntimeError("global curvature ratio failed")
    dyadic_ratio = sp.simplify(
        sp.Abs(second.subs(u, block))
        / sp.Abs(second.subs(u, 2 * block))
    )
    if dyadic_ratio != 4:
        raise RuntimeError("dyadic curvature ratio failed")

    return {
        "sourced_weighted_transform": (
            "After conjugating the negative-curvature phase, Vandehey "
            "Theorem 1.1 gives on a regular block [A,D] the starred "
            "stationary main plus an error of size "
            "O((V+|g(A)|){M/sqrt(T)+log(f'(A)-f'(D)+2)}), with the "
            "implicit constant depending on the derivative-comparability "
            "constants. This is an imported qualitative bound, not an "
            "explicit finite-height constant."
        ),
        "dyadic_hypotheses": (
            "On [U,2U], |f''| lies between alpha/(4U^2) and alpha/U^2, "
            "|f'''|<=2alpha/U^3, and |f''''|<=6alpha/U^4. Thus the sourced "
            "weighted theorem applies with M=U and T=alpha, after "
            "conjugation, uniformly in the block location."
        ),
        "global_guard": (
            "On [1,B], max|f''|/min|f''|=B^2. Therefore a theorem requiring "
            "uniform curvature comparability cannot be applied globally "
            "with B-independent implicit constants. A dyadic or smooth "
            "localization is compulsory."
        ),
        "internal_boundary_guard": (
            "Hard dyadic blocks introduce artificial half-endpoints. They "
            "must be rejoined algebraically so internal endpoint terms "
            "cancel, or replaced by a smooth partition of unity. Bounding "
            "each artificial endpoint separately changes the signed problem."
        ),
        "coarse_error_guard": (
            "Theorem 1.1 retains a logarithmic derivative-range loss and "
            "inexplicit constants. It validates the transform architecture "
            "but is not yet the h^2/400 remainder theorem. Vandehey's sharper "
            "endpoint analysis is a candidate for removing that loss."
        ),
        "symbolic": {
            "second_derivative": str(second),
            "third_derivative": str(third),
            "fourth_derivative": str(fourth),
            "global_curvature_ratio": str(global_ratio),
            "dyadic_curvature_ratio": str(dyadic_ratio),
        },
    }


def exact_payload() -> dict:
    return {
        "domain": (
            "Use the physical q=1, L>=50, fixed-N terminal/bulk chart, "
            "B=N-M-1, and T_0-epsilon<=xi<=T_0. Hold sigma, t, the "
            "normalizer, all terminal jets, and the Section 11.165 "
            "coefficient rows fixed along the auxiliary xi segment."
        ),
        "observation_handoff": (
            "The four physical observations V,N,A,Q and their four xi "
            "derivatives are fixed real linear forms of G_0,...,G_5. The "
            "physical x differentiation needed to define N and Q has already "
            "been performed exactly in their coefficient rows. Therefore one "
            "needs six endpoint-complete value transforms, not eight separate "
            "value/physical-derivative transforms. Uniform dependence on the "
            "fixed physical chart parameters is still required."
        ),
        "normalizer_guard": (
            "The fixed factor eta/S_a and the external common phase c_xi "
            "remain outside every Poisson integral. Endpoint carries c_xi, "
            "Hermitian products cancel it, and transpose products carry "
            "c_xi^2 exactly as in Section 11.165."
        ),
        "pi_provenance": (
            "The phase e(alpha log u-ru) uses e(y)=exp(2*pi*i*y) and "
            "alpha=xi/(2*pi). Negative stationary curvature contributes "
            "e(-1/8)=exp(-i*pi/4). These are the standard Poisson and "
            "stationary-signature normalizations recorded by Vandehey and "
            "the existing reciprocal-symbol guard; no new pi is inserted."
        ),
        "next_target": (
            "Choose a smooth dyadic localization or Vandehey's uniform "
            "endpoint transition, derive explicit value-remainder constants "
            "for H_0,...,H_5 uniformly on the full physical chart, rejoin all "
            "internal boundaries, and compose the two real hard endpoints "
            "with the terminal recurrence. Insert those six transforms into "
            "the eight observation rows, then apply the Section 11.166 "
            "Hermitian pairing before any absolute value."
        ),
        "exact_target": (
            "Prove Psi(T_0)-epsilon*integral_0^1 "
            "Psi'(T_0-theta*epsilon)dtheta<=h^2/400, preferably h^2/800."
        ),
        "proof_boundary": (
            "This proves the specialization of full starred Poisson summation "
            "to the six smooth bare moments, restoration of both full primal "
            "endpoints, modewise frequency closure, the negative-curvature "
            "stationary main, five exact stationary-transport and residual "
            "recurrences, two active-window jump laws, the reduction from "
            "eight observations to six value transforms, and dyadic "
            "admissibility for a sourced weighted B-process theorem. It "
            "proves no explicit uniform Poisson/B-process remainder constant, "
            "smooth endpoint-transition estimate, imaginary Hermitian "
            "coefficient bound, signed offset-sum gain, signed flow estimate, "
            "Phi_B bound, contact exclusion, retained aggregate or Xi theorem, "
            "Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or "
            "prize-level conclusion."
        ),
    }


def build_rows(
    exact: dict,
    moments: dict,
    poisson: dict,
    stationary: dict,
    transport: dict,
    cutoff: dict,
    dyadic: dict,
) -> list[PoissonRow]:
    return [
        PoissonRow("fptr_01_domain", "fixed-chart domain", "proved", "The transform is specialized on the exact auxiliary-frequency chart.", exact["domain"], "No physical parameter is varied inside this reduction."),
        PoissonRow("fptr_02_source", "primary analytic source", "imported_theorem", "The starred Poisson and weighted B-process conventions are sourced.", poisson["starred_source_formula"], VANDEHEY_URL),
        PoissonRow("fptr_03_moments", "six bare moments", "proved", "All required arithmetic values reduce to six smooth logarithmic sums.", moments["bare_moments"], exact["observation_handoff"]),
        PoissonRow("fptr_04_amplitude", "amplitude recurrence", "proved", "The six amplitudes form one exact logarithmic differential ladder.", moments["logarithmic_recurrence"], "Orders zero through five are audited."),
        PoissonRow("fptr_05_endpoint", "full primal endpoints", "proved", "The physical full sum is the starred Poisson sum plus two explicit half-endpoints.", poisson["full_endpoint_formula"] + " " + poisson["explicit_endpoint"], "The lower endpoint occurs only at order zero."),
        PoissonRow("fptr_06_poisson", "finite Poisson identity", "imported_theorem", "Every bare moment has one exact symmetric mode-integral representation.", poisson["starred_source_formula"] + " " + poisson["convergence_convention"], "The symmetric limit is retained."),
        PoissonRow("fptr_07_shift", "modewise frequency shift", "proved", "Both endpoints and every Poisson integral preserve the six-moment shift.", moments["frequency_shift"] + " " + poisson["modewise_shift"], "No conditionally convergent series is differentiated termwise."),
        PoissonRow("fptr_08_phase", "Poisson phase", "proved", "Only positive modes have reciprocal stationary points.", stationary["poisson_phase"], stationary["active_window"]),
        PoissonRow("fptr_09_saddle", "stationary data", "proved", "The reciprocal saddle, phase, curvature, and prefactor are explicit.", stationary["stationary_data"], "Continuous stationary algebra only."),
        PoissonRow("fptr_10_signature", "negative-curvature signature", "proved", "The one-variable main carries the required e(-1/8) factor.", stationary["signature"] + " " + stationary["stationary_main"], exact["pi_provenance"]),
        PoissonRow("fptr_11_transport", "stationary transport", "proved", "Differentiating the main produces phase shift plus one explicit transport correction.", transport["main_transport"], transport["transport_interpretation"]),
        PoissonRow("fptr_12_residual", "residual recurrence", "proved", "The fixed-chart first-frequency residual closes on six value remainders.", transport["mode_residual_recurrence"], transport["six_moment_closure"]),
        PoissonRow("fptr_13_observations", "eight-observation handoff", "proved", "No separate physical-derivative Poisson transform is needed after exact coefficient closure.", exact["observation_handoff"], "Uniform physical-chart dependence remains open."),
        PoissonRow("fptr_14_external", "external normalizer phases", "proved", "The fixed source factor and common phase remain outside all transforms.", exact["normalizer_guard"], "No normalizer phase is absorbed into an arithmetic sum."),
        PoissonRow("fptr_15_dual", "dual starred main", "proved", "The sharp active main has the standard two dual half-weights.", cutoff["dual_star_sum"], cutoff["fixed_chart_remainder"]),
        PoissonRow("fptr_16_upper", "upper cutoff crossing", "proved", "An upper dual mode enters with one exact jump law.", cutoff["upper_crossing"], cutoff["crossing_frequencies"]),
        PoissonRow("fptr_17_lower", "lower cutoff crossing", "proved", "A lower dual mode exits with the opposite exact jump law.", cutoff["lower_crossing"], cutoff["crossing_frequencies"]),
        PoissonRow("fptr_18_c1", "cutoff C1 guard", "proved", "The sharp-star main cannot be differentiated globally across a mode crossing.", cutoff["c1_guard"], "Use exact jumps or a uniform endpoint transition."),
        PoissonRow("fptr_19_dyadic", "dyadic B-process admissibility", "proved", "Every multiplicative block satisfies the sourced derivative hypotheses.", dyadic["dyadic_hypotheses"] + " " + dyadic["sourced_weighted_transform"], "The imported constant is not explicit."),
        PoissonRow("fptr_20_global", "global curvature guard", "proved", "A single global bounded-comparability invocation is invalid.", dyadic["global_guard"], dyadic["internal_boundary_guard"]),
        PoissonRow("fptr_21_coarse", "coarse remainder guard", "proved", "The sourced coarse theorem validates architecture but does not close the required scale.", dyadic["coarse_error_guard"], "No h^2/400 remainder bound is claimed."),
        PoissonRow("fptr_22_handoff", "analytic handoff", "open", "The six value remainders and hard endpoints are ready for explicit localization.", exact["next_target"], exact["exact_target"]),
        PoissonRow("fptr_23_boundary", "proof boundary", "proved", "Exact and imported statements are separated from the missing quantitative theorem.", exact["proof_boundary"], "This is not a proof of RH."),
    ]


def render_note(payload: dict) -> str:
    exact = payload["exact"]
    moments = payload["moment_amplitude_certificate"]
    poisson = payload["finite_poisson_certificate"]
    stationary = payload["stationary_phase_certificate"]
    transport = payload["stationary_transport_certificate"]
    cutoff = payload["cutoff_chart_certificate"]
    dyadic = payload["dyadic_b_process_certificate"]
    return f"""# Six-Moment Finite-Poisson And Stationary-Transport Reduction

Date: 2026-08-01

Status: exact reduction note with one imported weighted B-process theorem,
`0 explicit uniform remainder constants`, `0 signed flow bounds`, and this is not a proof of RH.

Primary analytic source:

```text
{VANDEHEY_URL}
```

## Domain And Moments

{exact['domain']}

{moments['bare_moments']}

{moments['logarithmic_recurrence']}

{moments['frequency_shift']}

{moments['lower_endpoint']}

## Full Poisson Formula

{poisson['starred_source_formula']}

{poisson['full_endpoint_formula']}

{poisson['explicit_endpoint']}

{poisson['modewise_shift']}

{poisson['convergence_convention']}

## Reciprocal Main

{stationary['poisson_phase']}

{stationary['active_window']}

{stationary['stationary_data']}

{stationary['signature']}

{stationary['stationary_main']}

## Stationary Transport

{transport['main_transport']}

{transport['transport_interpretation']}

{transport['mode_residual_recurrence']}

{transport['six_moment_closure']}

## Observation Handoff

{exact['observation_handoff']}

{exact['normalizer_guard']}

## Cutoff Charts

{cutoff['dual_star_sum']}

{cutoff['fixed_chart_remainder']}

{cutoff['upper_crossing']}

{cutoff['lower_crossing']}

{cutoff['c1_guard']}

{cutoff['crossing_frequencies']}

## Sourced Remainder Architecture

{dyadic['sourced_weighted_transform']}

{dyadic['dyadic_hypotheses']}

{dyadic['global_guard']}

{dyadic['internal_boundary_guard']}

{dyadic['coarse_error_guard']}

## Pi Provenance

{exact['pi_provenance']}

## Handoff

{exact['next_target']}

{exact['exact_target']}

## Proof Boundary

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    sources = load_sources()
    moments = moment_amplitude_certificate()
    poisson = finite_poisson_certificate()
    stationary = stationary_phase_certificate()
    transport = stationary_transport_certificate()
    cutoff = cutoff_chart_certificate()
    dyadic = dyadic_b_process_certificate()
    exact = exact_payload()
    rows = build_rows(
        exact,
        moments,
        poisson,
        stationary,
        transport,
        cutoff,
        dyadic,
    )
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "six-moment full-Poisson, stationary-main, transport, cutoff, "
            "and dyadic-admissibility reduction complete; explicit uniform "
            "remainder and signed flow estimates open"
        ),
        "source_audit": source_audit(sources),
        "moment_amplitude_certificate": moments,
        "finite_poisson_certificate": poisson,
        "stationary_phase_certificate": stationary,
        "stationary_transport_certificate": transport,
        "cutoff_chart_certificate": cutoff,
        "dyadic_b_process_certificate": dyadic,
        "exact": exact,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "bare_logarithmic_moments": 6,
            "full_poisson_identities": 6,
            "explicit_primal_half_endpoints": 2,
            "stationary_main_families": 6,
            "stationary_transport_identities": 5,
            "residual_transport_recurrences": 5,
            "physical_flow_observations": 8,
            "required_value_transforms": 6,
            "dual_cutoff_jump_laws": 2,
            "dyadic_b_process_hypothesis_audits": 1,
            "imported_inexplicit_local_b_process_bounds": 1,
            "explicit_uniform_remainder_constants": 0,
            "smooth_endpoint_transition_bounds": 0,
            "signed_flow_bounds": 0,
            "phi_b_bounds": 0,
        },
        "proof_boundary": exact["proof_boundary"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "built six-moment finite-Poisson/transport reduction: 23 rows, "
        "6 moments, 5 transport recurrences, 2 cutoff jump laws, "
        "0 explicit uniform remainder constants, 0 signed flow bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
