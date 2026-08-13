#!/usr/bin/env python3
"""Build the endpoint-resolved terminal/nonterminal relative-lift gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_terminal_nonterminal_"
    "relative_lift_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "terminal_bulk_phase": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_bulk_phase_"
        "completion_gate.json"
    ),
    "terminal_conditional": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_conditional_"
        "quadrature_gate.json"
    ),
    "terminal_obstruction": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_terminal_endpoint_"
        "obstruction_gate.json"
    ),
    "outer_endpoint": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_full_support_outer_endpoint_"
        "kernel_gate.json"
    ),
    "observation_image": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_observation_image_"
        "compression_gate.json"
    ),
    "physical_plin": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "amplitude_gate.json"
    ),
    "full_support": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_terminal_full_support_homotopy_gate.json"
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


def file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text, encoding="utf-8")
    os.replace(temporary, path)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def require_zero(expression: sp.Expr, label: str) -> None:
    value = sp.simplify(expression)
    if value != 0:
        raise RuntimeError(f"{label}: {value}")


def complex_symbol(name: str) -> tuple[sp.Expr, tuple[sp.Symbol, sp.Symbol]]:
    real, imag = sp.symbols(f"{name}r {name}i", real=True)
    return real + sp.I * imag, (real, imag)


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        if not path.exists():
            raise RuntimeError(f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["terminal_bulk_phase"]["counts"][
            "observation_coefficient_completions"
        ]
        == 1,
        "terminal-bulk phase source drifted",
    )
    require(
        payloads["terminal_conditional"]["counts"][
            "bulk_terminal_kernel_decompositions"
        ]
        == 1,
        "terminal conditional source drifted",
    )
    require(
        payloads["terminal_obstruction"]["counts"][
            "unique_order_one_obstructions"
        ]
        == 1,
        "terminal obstruction source drifted",
    )
    require(
        payloads["outer_endpoint"]["counts"][
            "hermitian_transpose_decompositions"
        ]
        == 1,
        "outer endpoint source drifted",
    )
    require(
        payloads["observation_image"]["counts"][
            "quadratic_observation_pairings"
        ]
        == 4,
        "observation image source drifted",
    )
    require(
        payloads["physical_plin"]["counts"]["physical_base_polynomials"]
        == 4,
        "physical P_lin source drifted",
    )
    require(
        payloads["full_support"]["counts"]["moving_tail_defect_identities"]
        == 1,
        "full-support source drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def audit_geometry() -> dict[str, int | str]:
    cases = 0
    for n_value in range(2, 90):
        for theta in [Fraction(0), Fraction(1, 11), Fraction(1, 2), Fraction(10, 11), Fraction(1)]:
            a_value = Fraction(n_value) + theta
            h = 1 / a_value
            require(Fraction(1, n_value) <= Fraction(3, 2) * h, "1/N bound")
            require(
                2 * h * a_value**2 / n_value**2
                <= Fraction(9, 2) * h,
                "alpha u_N/N^2 envelope",
            )
            require(Fraction(5, n_value) <= Fraction(15, 2) * h, "lower S1 envelope")
            cases += 1
    return {
        "geometry_cases": cases,
        "terminal_higher": "At mu=N, every relative-lift S_2/S_3 coefficient contains N^(-1), u_N, or alpha_Pu_N/N^2; the established boxes give N^(-1)<=(3/2)h and alpha_Pu_N/N^2<(9/2)h.",
        "lower_conditional": "At mu=1, 0<S_1^[N](alpha_P)<5/N<=(15/2)h, but the relative lifts also contain log N or log a+u_N. The lower endpoint must remain grouped with the interior.",
    }


def symbolic_certificate() -> dict:
    imaginary = sp.I
    t_n, u_n23, t_1, u_123, interior = sp.symbols(
        "T_N U_N23 T_1 U_123 I"
    )
    outer = t_n + u_n23 - t_1 - u_123 + interior
    remainder = u_n23 - t_1 - u_123 + interior
    require_zero(outer - t_n - remainder, "outer remainder split")

    lam, log_a, log_n, u_n = sp.symbols(
        "lambda log_a log_N u_N", real=True
    )
    p = sum(sp.symbols(f"p{degree}") * lam**degree for degree in range(5))
    p_dot = imaginary * (lam - log_a) * p
    p_h = sp.expand(p_dot + imaginary * u_n * p)
    p_t = sp.expand(p_dot - imaginary * u_n * p)
    relation = {u_n: log_a - log_n}
    require_zero(
        p_h.subs(relation) - imaginary * (lam - log_n) * p,
        "Hermitian relative lift",
    )
    require_zero(
        p_t.subs(relation)
        - imaginary * (lam - log_a - (log_a - log_n)) * p,
        "transpose relative lift",
    )
    require(sp.Poly(p_h, lam).degree() == 5, "Hermitian degree closure")
    require(sp.Poly(p_t, lam).degree() == 5, "transpose degree closure")

    carrier, outer_rest, moved = sp.symbols("C R M")
    dot_carrier, dot_rest, dot_moved = sp.symbols("dC dR dM")
    y_before = carrier + outer_rest
    y_after = carrier + moved + outer_rest - moved
    dy_before = dot_carrier + dot_rest
    dy_after = dot_carrier + dot_moved + dot_rest - dot_moved
    require_zero(y_after - y_before, "tie value transfer")
    require_zero(dy_after - dy_before, "tie lift transfer")
    require_zero(
        (dy_after + imaginary * u_n * y_after)
        - (dy_before + imaginary * u_n * y_before),
        "tie Hermitian transfer",
    )
    require_zero(
        (dy_after - imaginary * u_n * y_after)
        - (dy_before - imaginary * u_n * y_before),
        "tie transpose transfer",
    )

    amplitudes = {}
    rest_values = {}
    dot_values = {}
    for label in ("V", "N", "A", "Q"):
        amplitudes[label], _ = complex_symbol(f"a{label}")
        rest_values[label], _ = complex_symbol(f"y{label}")
        dot_values[label], _ = complex_symbol(f"dy{label}")
    h_cross = (
        amplitudes["V"] * sp.conjugate(dot_values["N"])
        + rest_values["V"] * sp.conjugate(-imaginary * u_n * amplitudes["N"])
        + amplitudes["N"] * sp.conjugate(dot_values["V"])
        + rest_values["N"] * sp.conjugate(-imaginary * u_n * amplitudes["V"])
        - amplitudes["A"] * sp.conjugate(dot_values["Q"])
        - rest_values["A"] * sp.conjugate(-imaginary * u_n * amplitudes["Q"])
        - amplitudes["Q"] * sp.conjugate(dot_values["A"])
        - rest_values["Q"] * sp.conjugate(-imaginary * u_n * amplitudes["A"])
    )
    t_cross = (
        amplitudes["V"] * dot_values["N"]
        + rest_values["V"] * (-imaginary * u_n * amplitudes["N"])
        + amplitudes["N"] * dot_values["V"]
        + rest_values["N"] * (-imaginary * u_n * amplitudes["V"])
        - amplitudes["A"] * dot_values["Q"]
        - rest_values["A"] * (-imaginary * u_n * amplitudes["Q"])
        - amplitudes["Q"] * dot_values["A"]
        - rest_values["Q"] * (-imaginary * u_n * amplitudes["A"])
    )
    d_h = {
        label: dot_values[label] + imaginary * u_n * rest_values[label]
        for label in rest_values
    }
    d_t = {
        label: dot_values[label] - imaginary * u_n * rest_values[label]
        for label in rest_values
    }
    h_compact = (
        amplitudes["V"] * sp.conjugate(d_h["N"])
        + amplitudes["N"] * sp.conjugate(d_h["V"])
        - amplitudes["A"] * sp.conjugate(d_h["Q"])
        - amplitudes["Q"] * sp.conjugate(d_h["A"])
    )
    t_compact = (
        amplitudes["V"] * d_t["N"]
        + amplitudes["N"] * d_t["V"]
        - amplitudes["A"] * d_t["Q"]
        - amplitudes["Q"] * d_t["A"]
    )
    require_zero(
        sp.expand_complex(sp.re(h_cross - h_compact)),
        "mixed Hermitian compression",
    )
    require_zero(sp.expand_complex(t_cross - t_compact), "mixed transpose compression")

    kappa, e_n, exp_s, s_1n, s_2n, s_3n = sp.symbols(
        "kappa e_N expS_N S_1N S_2N S_3N"
    )
    n_value, alpha = sp.symbols("N alpha", positive=True)
    p_n, p_n_prime, g_n = sp.symbols("P_N P_Nprime g_N")
    tau_0 = kappa * e_n * exp_s * s_1n
    require_zero(
        tau_0 * (-imaginary * u_n * p_n)
        - (-imaginary * u_n) * tau_0 * p_n,
        "terminal conditional lift",
    )

    def terminal_higher(value: sp.Expr, derivative_combo: sp.Expr) -> sp.Expr:
        return -kappa**2 * e_n * exp_s * (
            derivative_combo * s_2n / n_value
            + alpha * value * s_3n / n_value**2
        )

    h_value = sp.Integer(0)
    h_derivative = imaginary * p_n
    h_endpoint = terminal_higher(h_value, h_derivative)
    require_zero(
        h_endpoint
        + imaginary * kappa**2 * e_n * exp_s * p_n * s_2n / n_value,
        "terminal Hermitian higher endpoint",
    )
    t_value = -2 * imaginary * u_n * p_n
    t_derivative = imaginary * (
        p_n - 2 * u_n * (p_n_prime + g_n * p_n)
    )
    t_endpoint = terminal_higher(t_value, t_derivative)
    t_expected = -imaginary * kappa**2 * e_n * exp_s * (
        (p_n - 2 * u_n * (p_n_prime + g_n * p_n)) * s_2n / n_value
        - 2 * alpha * u_n * p_n * s_3n / n_value**2
    )
    require_zero(t_endpoint - t_expected, "terminal transpose higher endpoint")

    p_0, s_1lower = sp.symbols("P_0 S_1lower")
    lower_h = -kappa * (-imaginary * log_n * p_0) * s_1lower
    lower_t = -kappa * (-imaginary * (log_a + u_n) * p_0) * s_1lower
    require_zero(
        lower_h - imaginary * kappa * log_n * p_0 * s_1lower,
        "lower Hermitian conditional trace",
    )
    require_zero(
        lower_t - imaginary * kappa * (log_a + u_n) * p_0 * s_1lower,
        "lower transpose conditional trace",
    )

    r_n, r_l, rx_n, rx_l, c_n, c_l, q_n, q_l = sp.symbols(
        "R_N R_l Rx_N Rx_l C_N C_l Q_N Q_l"
    )
    a_n = r_n * c_n
    a_l = r_l * c_l
    n_n = r_n * q_n + rx_n * c_n
    n_l = r_l * q_l + rx_l * c_l
    b_t = c_n * n_l + n_n * c_l - a_n * q_l - q_n * a_l
    b_factored = (
        (r_l - r_n) * (c_n * q_l - q_n * c_l)
        + (rx_l + rx_n) * c_n * c_l
    )
    require_zero(b_t - b_factored, "transpose carrier bilinear factorization")
    require_zero(
        b_t.subs({r_l: r_n, rx_l: rx_n, c_l: c_n, q_l: q_n})
        - 2 * rx_n * c_n**2,
        "transpose terminal diagonal",
    )

    lam = sp.symbols("lambda", real=True)
    c_poly = sum(sp.symbols(f"c{j}") * lam**j for j in range(3))
    n_poly = sum(sp.symbols(f"n{j}") * lam**j for j in range(5))
    a_poly = sum(sp.symbols(f"a{j}") * lam**j for j in range(4))
    q_poly = sum(sp.symbols(f"q{j}") * lam**j for j in range(4))
    cv, nv, av, qv = sp.symbols("CV NV AV QV")
    b_poly = sp.expand(cv * n_poly + nv * c_poly - av * q_poly - qv * a_poly)
    require(sp.Poly(b_poly, lam).degree() == 4, "mixed bilinear degree")
    require(
        sp.Poly((lam - log_n) * b_poly, lam).degree() == 5,
        "Hermitian mixed degree",
    )
    require(
        sp.Poly((lam - log_a - u_n) * b_poly, lam).degree() == 5,
        "transpose mixed degree",
    )

    geometry = audit_geometry()
    return {
        "outer_remainder": {
            "endpoint_operators": "Write T_mu[P]=e(alpha_Plog mu)kappa A_P(mu)S_1^[N](alpha_P/mu) and U_mu[P]=-e(alpha_Plog mu)kappa^2{A_P'(mu)S_2^[N](alpha_P/mu)+(alpha_P/mu^2)A_P(mu)S_3^[N](alpha_P/mu)}.",
            "interior": "Write I_N[P]=kappa^2 sum_(r notin mathcal R_N)^sym integral_1^N e(phi_r)C_rA_Pdu; this series is absolute only after the joined twofold integration by parts.",
            "split": "The exact outer functional is mathscr O_N[P]=T_N[P]+mathscr R_N[P], where mathscr R_N[P]=U_N[P]-T_1[P]-U_1[P]+I_N[P].",
            "observations": "For X in {V,mathcal N,A,Q}, rho_X=Re{nu c_xi mathscr R_N[P_X]} and dot rho_X=Re{nu c_xi mathscr R_N[i(lambda-log a)P_X]} on a fixed roster cell.",
            "linear": "The complete linear outer functional splits without componentwise budgets as mathscr O_N[P_lin]=T_N[P_lin]+mathscr R_N[P_lin].",
        },
        "relative_lifts": {
            "terminal": "T_N[i(lambda-log a)P]=-iu_NT_N[P].",
            "hermitian": "D_Hmathscr R_N[P]:=mathscr R_N[i(lambda-log a)P]+iu_Nmathscr R_N[P]=mathscr R_N[i(lambda-log N)P].",
            "transpose": "D_Tmathscr R_N[P]:=mathscr R_N[i(lambda-log a)P]-iu_Nmathscr R_N[P]=mathscr R_N[i(lambda-log a-u_N)P].",
            "atoms": "On a retained carrier atom with lift -iu_q, D_H=-i(u_q-u_N) and D_T=-i(u_q+u_N); the q=N Hermitian relative lift vanishes exactly.",
            "ties": "At a reciprocal tie, transfer the complete mode I_P and its lift from mathscr R_N to the retained carrier. Their sum Y and both dot Y+iu_NY and dot Y-iu_NY are invariant; the complement is never differentiated alone.",
        },
        "mixed_current": {
            "complete": "Let Y_X be the retained carrier observation plus mathscr R_N[P_X], and let A_X=tau_0P_X(log N), tau_0=kappa e(alpha_Plog N)exp(S(log N))S_(1,N). Then the complete terminal/nonterminal current uses D_HY_X=dot Y_X+iu_NY_X and D_TY_X=dot Y_X-iu_NY_X.",
            "hermitian": "Re H_(a*)=Re{A_Vconjugate(D_HY_N)+A_Nconjugate(D_HY_V)-A_Aconjugate(D_HY_Q)-A_Qconjugate(D_HY_A)}.",
            "transpose": "T_(a*)=A_VD_TY_N+A_ND_TY_V-A_AD_TY_Q-A_QD_TY_A.",
            "projection": "The actual mixed real current is (1/2)Re{|nu|^2H_(a*)+nu^2c_xi^2T_(a*)}; the fixed affine row and moved-tail defect remain separate.",
        },
        "polynomial_compression": {
            "hermitian_polynomial": "B_(H,N)(lambda)=conjugate(C_N)mathcal N(lambda)+conjugate(mathcal N_N)C(lambda)-conjugate(A_N)Q(lambda)-conjugate(Q_N)A(lambda).",
            "transpose_polynomial": "B_(T,N)(lambda)=C_Nmathcal N(lambda)+mathcal N_NC(lambda)-A_NQ(lambda)-Q_NA(lambda).",
            "hermitian_functional": "The outer-remainder cross satisfies Re H_(a rho)=Re{conjugate(tau_0)mathscr R_N[i(lambda-log N)B_(H,N)(lambda)]}.",
            "transpose_functional": "The outer-remainder cross satisfies T_(a rho)=tau_0mathscr R_N[i(lambda-log a-u_N)B_(T,N)(lambda)].",
            "degree": "Both B_(H,N) and B_(T,N) have degree at most four; each relative lift has degree at most five, so the source-complete mixed cross remains inside the certified six-moment closure.",
        },
        "endpoint_audit": {
            "terminal_hermitian": "For P_H=i(lambda-log N)P, U_N[P_H]=-i kappa^2 e(alpha_Plog N)exp(S(log N))N^(-1)P(log N)S_2^[N](alpha_P/N); its terminal S_1 and S_3 values vanish.",
            "terminal_transpose": "For P_T=i(lambda-log a-u_N)P, U_N[P_T]=-i kappa^2 e(alpha_Plog N)exp(S(log N)){N^(-1)[P-2u_N(P'+gP)]S_2-2alpha_Pu_NN^(-2)P S_3} at lambda=log N.",
            "lower_hermitian": "The lower conditional piece is -T_1[P_H]=i kappa(log N)P(0)S_1^[N](alpha_P).",
            "lower_transpose": "The lower conditional piece is -T_1[P_T]=i kappa(log a+u_N)P(0)S_1^[N](alpha_P).",
        },
        "carrier_bilinear": {
            "factorization": "B_(T,N)(lambda)=[R(lambda)-R_N]{C_NQ(lambda)-Q_NC(lambda)}+[R_x(lambda)+R_(N,x)]C_NC(lambda).",
            "diagonal": "At lambda=log N, B_(T,N)=2R_(N,x)C_N^2, while the Hermitian relative factor lambda-log N is zero.",
            "meaning": "The retained q=N cross is therefore Hermitian-null and transpose-suppressed; it does not cancel the lower endpoint or interior algebraically.",
        },
        "bounds": {
            "terminal_higher": geometry["terminal_higher"],
            "lower_conditional": geometry["lower_conditional"],
        },
        "audit": {
            "geometry_cases": geometry["geometry_cases"],
            "relative_lift_families": 2,
            "mixed_polynomial_functionals": 2,
            "outer_remainder_components": 4,
        },
        "handoff": {
            "decision": "Terminal higher-endpoint and retained q=N cross terms keep explicit h mechanisms. The unresolved completed channel is the lower endpoint plus absolute interior, with both relative phase families still joined.",
            "route": "Expand the two degree-five mathscr R_N functionals into their lower U_1 and absolute-interior pieces, combine them with the finite q<N carrier sums and correction perturbation, and prove or falsify a joined Abel/Mangoldt/Vaughan or reciprocal-pair estimate before taking moduli.",
        },
    }


def build_rows(cert: dict) -> list[GateRow]:
    outer = cert["outer_remainder"]
    lifts = cert["relative_lifts"]
    mixed = cert["mixed_current"]
    poly = cert["polynomial_compression"]
    endpoint = cert["endpoint_audit"]
    carrier = cert["carrier_bilinear"]
    bounds = cert["bounds"]
    handoff = cert["handoff"]
    rows = [
        GateRow("tnr_01_endpoint", "endpoint operators", "proved", "The S_1 and S_2/S_3 endpoint rows are separated by exact linear operators.", outer["endpoint_operators"], "No endpoint is bounded componentwise."),
        GateRow("tnr_02_interior", "absolute interior", "proved", "The joined twofold interior is one exact operator.", outer["interior"], "Absolute convergence begins only after reassembly."),
        GateRow("tnr_03_split", "outer decomposition", "proved", "The terminal conditional block has an exact four-component complement.", outer["split"], "The four pieces remain joined."),
        GateRow("tnr_04_components", "remainder roster", "proved", "The nonterminal outer remainder consists of terminal higher, lower conditional, lower higher, and interior pieces.", "mathscr R_N=U_N-T_1-U_1+I_N.", "Signs are retained exactly."),
        GateRow("tnr_05_observations", "real observations", "proved", "The real rho rows are source projections of one complex remainder operator.", outer["observations"], "No eight-row norm is introduced."),
        GateRow("tnr_06_linear", "linear compression", "proved", "The linear outer correction remains one P_lin functional after the split.", outer["linear"], "The terminal and remainder pieces are not estimated separately here."),
        GateRow("tnr_07_terminal_lift", "terminal lift", "proved", "The terminal conditional block has the common -iu_N lift.", lifts["terminal"], "This is a fixed-roster identity."),
        GateRow("tnr_08_hlift", "Hermitian relative lift", "proved", "The Hermitian cross is measured relative to log N.", lifts["hermitian"], "The complement alone is not differentiated at ties."),
        GateRow("tnr_09_tlift", "transpose relative lift", "proved", "The transpose cross uses the phase-sum relative lift.", lifts["transpose"], "No phase family is discarded."),
        GateRow("tnr_10_atoms", "carrier phases", "proved", "Relative lifts reproduce the exact phase difference and phase sum on every carrier atom.", lifts["atoms"], "This does not sign either sum."),
        GateRow("tnr_11_ties", "tie transfer", "proved", "The completed relative lifts are invariant under a half-open mode transfer.", lifts["ties"], "Only the joined carrier plus remainder is differentiated."),
        GateRow("tnr_12_complete", "complete observation", "proved", "The terminal cross uses one complete nonterminal observation vector.", mixed["complete"], "The fixed endpoint edge remains external."),
        GateRow("tnr_13_hermitian", "mixed Hermitian current", "proved", "All mixed Hermitian pairings compress through D_H.", mixed["hermitian"], "Only its real projection enters."),
        GateRow("tnr_14_transpose", "mixed transpose current", "proved", "All mixed transpose pairings compress through D_T.", mixed["transpose"], "Its complex phase remains."),
        GateRow("tnr_15_projection", "source projection", "proved", "The source factors match the established H_out/T_out polarization.", mixed["projection"], "The affine and moved-tail rows are not absorbed."),
        GateRow("tnr_16_hpoly", "Hermitian polynomial", "proved", "The outer Hermitian cross has one terminal-coefficient polynomial.", poly["hermitian_polynomial"], "Its coefficients are complex."),
        GateRow("tnr_17_tpoly", "transpose polynomial", "proved", "The outer transpose cross has one terminal-coefficient polynomial.", poly["transpose_polynomial"], "It is not signed."),
        GateRow("tnr_18_hfunctional", "Hermitian functional", "proved", "Eight Hermitian cross terms reduce to one degree-five remainder functional after real projection.", poly["hermitian_functional"], "No modulus is taken."),
        GateRow("tnr_19_tfunctional", "transpose functional", "proved", "Eight transpose cross terms reduce to one degree-five remainder functional.", poly["transpose_functional"], "No modulus is taken."),
        GateRow("tnr_20_degree", "six-moment closure", "proved", "The completed mixed channel requires no seventh moment.", poly["degree"], "Degree closure is not a size bound."),
        GateRow("tnr_21_terminal_h", "terminal Hermitian endpoint", "proved", "The relative Hermitian terminal higher endpoint leaves only an S_2 term.", endpoint["terminal_hermitian"], "Its coefficient still depends on the physical row."),
        GateRow("tnr_22_terminal_t", "terminal transpose endpoint", "proved", "The relative transpose terminal higher endpoint retains explicit N^(-1) or u_N factors.", endpoint["terminal_transpose"], "No signed cancellation is inferred."),
        GateRow("tnr_23_lower_h", "lower Hermitian endpoint", "proved", "The lower conditional Hermitian trace carries log N times the small lower kernel.", endpoint["lower_hermitian"], "It must remain joined with the interior."),
        GateRow("tnr_24_lower_t", "lower transpose endpoint", "proved", "The lower conditional transpose trace carries log a+u_N times the small lower kernel.", endpoint["lower_transpose"], "It must remain joined with the interior."),
        GateRow("tnr_25_factor", "carrier bilinear", "proved", "The transpose coefficient has an exact rate/determinant factorization.", carrier["factorization"], "The factorization is not sign-definite."),
        GateRow("tnr_26_diagonal", "terminal diagonal", "proved", "The q=N Hermitian cross vanishes and the transpose row collapses to R_(N,x).", carrier["diagonal"], "This controls only the local diagonal."),
        GateRow("tnr_27_scale", "endpoint scale guard", "guard_validated", "Terminal-side higher pieces preserve h mechanisms, while lower logarithms remain visible.", bounds["terminal_higher"] + " " + bounds["lower_conditional"], "No complete h^2 estimate follows."),
        GateRow("tnr_28_route", "route decision", "open", "Prove or falsify the joined lower-endpoint/interior relative-lift estimate.", handoff["route"], "No signed completed bound is proved."),
        GateRow("tnr_29_pi", "pi provenance", "proved", "All pi factors remain Fourier-normalization factors.", "kappa=1/(2pi i) and alpha_P=xi/(2pi) are inherited; the relative lifts introduce no new pi.", "No geometric fitting is used."),
        GateRow("tnr_30_boundary", "proof boundary", "guard_validated", "This is an exact endpoint and phase compression, not a proof of RH.", "No signed lower/interior estimate, complete outer bound, quadratic residual bound, full-current theorem, contact exclusion, Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.", "This is not a proof of RH."),
    ]
    require(len(rows) == 30, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    outer = cert["outer_remainder"]
    lifts = cert["relative_lifts"]
    mixed = cert["mixed_current"]
    poly = cert["polynomial_compression"]
    endpoint = cert["endpoint_audit"]
    carrier = cert["carrier_bilinear"]
    bounds = cert["bounds"]
    handoff = cert["handoff"]
    return f"""# Full-Support Terminal/Nonterminal Relative-Lift Gate

Date: 2026-08-02

Status: exact endpoint-resolved remainder, tie-safe completed lifts, and degree-five mixed-current compression proved; signed lower/interior estimate open; not a proof of RH.

## Outer Remainder

{outer['endpoint_operators']}

{outer['interior']}

{outer['split']}

{outer['observations']}

{outer['linear']}

## Relative Lifts

{lifts['terminal']}

{lifts['hermitian']}

{lifts['transpose']}

{lifts['atoms']}

{lifts['ties']}

## Mixed Current

{mixed['complete']}

{mixed['hermitian']}

{mixed['transpose']}

{mixed['projection']}

## Polynomial Compression

{poly['hermitian_polynomial']}

{poly['transpose_polynomial']}

{poly['hermitian_functional']}

{poly['transpose_functional']}

{poly['degree']}

## Endpoint Audit

{endpoint['terminal_hermitian']}

{endpoint['terminal_transpose']}

{endpoint['lower_hermitian']}

{endpoint['lower_transpose']}

{carrier['factorization']}

{carrier['diagonal']}

{carrier['meaning']}

{bounds['terminal_higher']}

{bounds['lower_conditional']}

## Handoff

{handoff['decision']}

{handoff['route']}

## Pi Provenance

The constants `kappa=1/(2pi i)` and `alpha_P=xi/(2pi)` are inherited from the fixed Fourier character. Relative lifting introduces only logarithmic differences and sums. No circle, polygon, plotted symmetry, or fitted constant is introduced.

## Proof Boundary

This gate proves the exact four-piece outer remainder, its real observation projection, fixed-cell relative lifts, half-open tie-transfer invariance of the completed observation, exact Hermitian and transpose mixed-current compression, two degree-five polynomial functionals, terminal higher-endpoint reductions, lower conditional traces, and the transpose carrier-bilinear factorization. It proves no signed lower-endpoint/interior estimate, complete outer-complement or quadratic residual bound, signed full current, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "outer_remainder_decompositions": 1,
        "outer_remainder_components": 4,
        "real_remainder_observation_rows": 8,
        "relative_lift_families": 2,
        "tie_transfer_invariants": 4,
        "mixed_current_compressions": 2,
        "mixed_polynomial_functionals": 2,
        "maximum_mixed_polynomial_degree": 5,
        "terminal_higher_endpoint_reductions": 2,
        "lower_conditional_relative_traces": 2,
        "transpose_carrier_bilinear_factorizations": 1,
        "terminal_hermitian_diagonal_nulls": 1,
        "geometry_checks": certificate["audit"]["geometry_cases"],
        "signed_lower_interior_bounds": 0,
        "complete_outer_complement_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_full_current_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "endpoint-resolved terminal/nonterminal relative lifts proved; signed lower/interior estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact outer-remainder split, completed relative lifts and tie invariance, Hermitian/transpose mixed-current compression into two degree-five functionals, terminal endpoint reductions, lower conditional traces, and a carrier-bilinear factorization. It proves no signed lower/interior estimate, complete outer-complement or quadratic residual bound, signed full-current or Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Q209, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built full-support terminal/nonterminal relative-lift gate: "
        f"{counts['rows']} rows, {counts['outer_remainder_components']} remainder components, "
        f"{counts['relative_lift_families']} relative lifts, "
        f"{counts['mixed_polynomial_functionals']} degree-five functionals, "
        f"{counts['tie_transfer_invariants']} tie invariants, "
        f"{counts['signed_lower_interior_bounds']} signed lower/interior bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
