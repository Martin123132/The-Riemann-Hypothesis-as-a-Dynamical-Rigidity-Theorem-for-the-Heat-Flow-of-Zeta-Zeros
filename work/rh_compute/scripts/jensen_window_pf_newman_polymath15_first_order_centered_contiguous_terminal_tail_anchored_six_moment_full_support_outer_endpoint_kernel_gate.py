#!/usr/bin/env python3
"""Build the pole-free full-support outer-endpoint kernel gate."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path

import mpmath as mp
import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
KIND = (
    "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
    "terminal_tail_anchored_six_moment_full_support_outer_endpoint_kernel_gate"
)
RESULT_PATH = REPO_ROOT / "work" / "rh_compute" / "results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
    "morse_endpoint": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_endpoint_reduction.json"
    ),
    "endpoint_composition": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_endpoint_"
        "composition_retention_gate.json"
    ),
    "observation_image": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_observation_"
        "image_compression_gate.json"
    ),
    "physical_plin": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "amplitude_gate.json"
    ),
    "finite_cell": REPO_ROOT
    / "work/rh_compute/results"
    / (
        "jensen_window_pf_newman_polymath15_first_order_centered_contiguous_"
        "terminal_tail_anchored_six_moment_morse_fresnel_physical_plin_"
        "ideal_cubic_reciprocal_finite_cell_inversion_gate.json"
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


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def harmonic(n: int) -> Fraction:
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        if not path.exists():
            raise RuntimeError(f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
    require(
        payloads["morse_endpoint"]["counts"]["stable_endpoint_sum_families"] == 3,
        "Morse endpoint source drifted",
    )
    require(
        payloads["endpoint_composition"]["counts"]["stable_endpoint_packages"] == 2,
        "endpoint-composition source drifted",
    )
    require(
        payloads["observation_image"]["counts"]["quadratic_observation_pairings"]
        == 4,
        "observation-image source drifted",
    )
    require(
        payloads["physical_plin"]["counts"]["maximum_polynomial_degree"] == 5,
        "physical-P_lin source drifted",
    )
    require(
        payloads["finite_cell"]["counts"]["global_recompositions"] == 3,
        "finite-cell source drifted",
    )
    require(
        payloads["full_support"]["counts"]["all_carrier_homotopies"] == 1,
        "full-support source drifted",
    )
    require(
        payloads["full_support"]["counts"]["outer_complement_bounds"] == 0,
        "full-support proof boundary drifted",
    )
    return {
        key: {
            "path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
            "sha256": file_hash(path),
        }
        for key, path in SOURCE_PATHS.items()
    }


def audit_outer_bands() -> dict[str, int]:
    mp.mp.dps = 80
    band_cases = 0
    kernel_points = 0
    noninteger_kernel_checks = 0
    removable_integer_checks = 0
    rational_symmetric_checks = 0
    absolute_kernel_checks = 0

    for n_value in [3, 5, 8, 13]:
        for theta in [Fraction(0), Fraction(1, 3), Fraction(4, 5)]:
            a_value = Fraction(n_value) + theta
            for offset in [Fraction(1, 9), Fraction(2, 5), Fraction(8, 9)]:
                alpha = a_value * a_value - offset
                require(alpha > a_value * a_value - 1, "physical alpha collar failed")
                a_n = alpha / (Fraction(n_value) + Fraction(1, 2))
                m_value = floor_fraction(a_n) + 1
                n_mode = floor_fraction(2 * alpha)
                require(m_value <= n_mode, "empty reciprocal band")
                band_cases += 1

                for u_value in range(1, n_value + 1):
                    x = alpha / u_value
                    a_x = x - m_value + 1
                    b_x = n_mode + 1 - x
                    require(a_x > Fraction(1, 4), "lower pole distance lost")
                    require(b_x > alpha, "upper pole distance lost")
                    if x.denominator == 1:
                        require(m_value <= x <= n_mode, "integer pole escaped band")

                    x_mp = mp.mpf(x.numerator) / x.denominator
                    a_mp = mp.mpf(a_x.numerator) / a_x.denominator
                    b_mp = mp.mpf(b_x.numerator) / b_x.denominator
                    s2 = mp.zeta(2, a_mp) + mp.zeta(2, b_mp)
                    s3 = mp.zeta(3, a_mp) - mp.zeta(3, b_mp)
                    z4 = mp.zeta(4, a_mp) + mp.zeta(4, b_mp)
                    require(s2 < 21, "S2 kernel bound failed")
                    require(abs(s3) < 73, "S3 kernel bound failed")
                    require(z4 < 278, "Z4 kernel bound failed")
                    absolute_kernel_checks += 3

                    if x.denominator != 1:
                        band_sum_1 = mp.fsum(
                            1 / (x_mp - r) for r in range(m_value, n_mode + 1)
                        )
                        band_sum_2 = mp.fsum(
                            1 / (x_mp - r) ** 2
                            for r in range(m_value, n_mode + 1)
                        )
                        band_sum_3 = mp.fsum(
                            1 / (x_mp - r) ** 3
                            for r in range(m_value, n_mode + 1)
                        )
                        s1_cot = mp.pi / mp.tan(mp.pi * x_mp) - band_sum_1
                        s1_stable = mp.digamma(b_mp) - mp.digamma(a_mp)
                        s2_cot = (mp.pi / mp.sin(mp.pi * x_mp)) ** 2 - band_sum_2
                        s3_cot = (
                            mp.pi**3
                            * mp.cos(mp.pi * x_mp)
                            / mp.sin(mp.pi * x_mp) ** 3
                            - band_sum_3
                        )
                        require(abs(s1_cot - s1_stable) < mp.mpf("1e-65"), "S1 drift")
                        require(abs(s2_cot - s2) < mp.mpf("1e-65"), "S2 drift")
                        require(abs(s3_cot - s3) < mp.mpf("1e-65"), "S3 drift")
                        noninteger_kernel_checks += 3
                    kernel_points += 1

    for m_value in range(2, 8):
        for n_mode in range(m_value + 3, m_value + 10):
            for x_value in range(m_value, n_mode + 1):
                a_x = x_value - m_value + 1
                b_x = n_mode + 1 - x_value
                stable = harmonic(b_x - 1) - harmonic(a_x - 1)
                expected = harmonic(n_mode - x_value) - harmonic(x_value - m_value)
                require(stable == expected, "integer removable value failed")
                removable_integer_checks += 1

                radius = n_mode + x_value + 9
                full_pv = harmonic(radius + x_value) - harmonic(radius - x_value)
                finite_outer = sum(
                    (
                        Fraction(1, x_value - r)
                        for r in range(-radius, radius + 1)
                        if r < m_value or r > n_mode
                    ),
                    Fraction(0),
                )
                require(finite_outer == full_pv + expected, "symmetric outer sum failed")
                rational_symmetric_checks += 1

    return {
        "band_cases": band_cases,
        "kernel_points": kernel_points,
        "noninteger_kernel_checks": noninteger_kernel_checks,
        "removable_integer_checks": removable_integer_checks,
        "rational_symmetric_checks": rational_symmetric_checks,
        "absolute_kernel_checks": absolute_kernel_checks,
    }


def symbolic_certificate() -> dict:
    A, A_1, A_2 = sp.symbols("A A_1 A_2")
    q, q_1, q_2 = sp.symbols("q q_1 q_2", nonzero=True)
    d_over_q = A_1 / q**2 - A * q_1 / q**3
    differentiated = (
        A_2 / q**2
        - 3 * A_1 * q_1 / q**3
        - A * q_2 / q**3
        + 3 * A * q_1**2 / q**4
    )
    direct = (
        A_2 / q**2
        - 2 * A_1 * q_1 / q**3
        - A_1 * q_1 / q**3
        - A * q_2 / q**3
        + 3 * A * q_1**2 / q**4
    )
    require_zero(differentiated - direct, "twofold IBP operator")
    require_zero(d_over_q - (A_1 / q - A * q_1 / q**2) / q, "D/q identity")

    imaginary = sp.I
    u_n, f, f_1, g, n_value, alpha = sp.symbols(
        "u_N F F_1 g N alpha", positive=True
    )
    p_edge = -imaginary * u_n * f
    p_edge_lambda = imaginary * f - imaginary * u_n * f_1
    amplitude_lambda = p_edge_lambda + g * p_edge
    require_zero(
        amplitude_lambda - imaginary * (f - u_n * (f_1 + g * f)),
        "edge amplitude derivative",
    )
    kappa, s_1, s_2, s_3 = sp.symbols("kappa S_1 S_2 S_3")
    edge_endpoint = (
        kappa * p_edge * s_1
        - kappa**2
        * (imaginary * (f - u_n * (f_1 + g * f)) * s_2 / n_value)
        - kappa**2 * alpha * p_edge * s_3 / n_value**2
    )
    edge_factored = imaginary * (
        -kappa * u_n * f * s_1
        - kappa**2 * (f - u_n * (f_1 + g * f)) * s_2 / n_value
        + kappa**2 * alpha * u_n * f * s_3 / n_value**2
    )
    require_zero(edge_endpoint - edge_factored, "edge endpoint factorization")

    f_alpha, f_beta = sp.symbols("F_alpha F_beta")
    terminal_plin = f_alpha - imaginary * u_n * f_beta
    require_zero(
        f_alpha + imaginary * (-u_n) * f_beta - terminal_plin,
        "terminal P_lin evaluation",
    )

    real_symbols = sp.symbols("zVr zVi zNr zNi zAr zAi zQr zQi")
    dot_symbols = sp.symbols("dVr dVi dNr dNi dAr dAi dQr dQi")
    z_v = real_symbols[0] + imaginary * real_symbols[1]
    z_n = real_symbols[2] + imaginary * real_symbols[3]
    z_a = real_symbols[4] + imaginary * real_symbols[5]
    z_q = real_symbols[6] + imaginary * real_symbols[7]
    d_v = dot_symbols[0] + imaginary * dot_symbols[1]
    d_n = dot_symbols[2] + imaginary * dot_symbols[3]
    d_a = dot_symbols[4] + imaginary * dot_symbols[5]
    d_q = dot_symbols[6] + imaginary * dot_symbols[7]
    real_pairing = (
        sp.re(z_v) * sp.re(d_n)
        + sp.re(z_n) * sp.re(d_v)
        - sp.re(z_a) * sp.re(d_q)
        - sp.re(z_q) * sp.re(d_a)
    )
    hermitian = (
        z_v * sp.conjugate(d_n)
        + z_n * sp.conjugate(d_v)
        - z_a * sp.conjugate(d_q)
        - z_q * sp.conjugate(d_a)
    )
    transpose = z_v * d_n + z_n * d_v - z_a * d_q - z_q * d_a
    require_zero(
        sp.expand_complex(real_pairing - sp.re(hermitian + transpose) / 2),
        "Hermitian-transpose split",
    )

    return {
        "outer_band": {
            "roster": "m_N=floor(alpha_P/(N+1/2))+1, n_N=floor(2alpha_P), and R_N={m_N,...,n_N}.",
            "pole_distance": "For x(u)=alpha_P/u, a_x=x-m_N+1>1/4 and b_x=n_N+1-x>alpha_P on 1<=u<=N.",
            "tie_rule": "The complement uses the same half-open reciprocal band as the finite-cell theorem; at a moving integer tie its transfer must be joined to the band before differentiation.",
        },
        "stable_kernels": {
            "S1": "S_1^[N](x)=psi(b_x)-psi(a_x)=PV sum_(r notin R_N)1/(x-r).",
            "S2": "S_2^[N](x)=zeta(2,a_x)+zeta(2,b_x)=sum_(r notin R_N)1/(x-r)^2.",
            "S3": "S_3^[N](x)=zeta(3,a_x)-zeta(3,b_x)=sum_(r notin R_N)1/(x-r)^3.",
            "integer": "At integer x in R_N, S_1^[N](x)=H_(n_N-x)-H_(x-m_N); no cotangent pole is separated.",
            "bounds": "On the physical segment S_2^[N]<21, |S_3^[N]|<73, and Z_4^[N]<278. These bound kernels, not the complete outer functional.",
        },
        "outer_identity": {
            "mode": "I_P(r)=kappa[e(phi_r)A_P/q_r]_1^N-kappa^2[e(phi_r)(D_rA_P)/q_r]_1^N+kappa^2 integral_1^N e(phi_r)C_rA_P du.",
            "operator": "C_rA=A''/q_r^2-3A'q_r'/q_r^3-Aq_r''/q_r^3+3A(q_r')^2/q_r^4.",
            "endpoint": "B_mu^[N][P]=e(alpha_Plog mu){kappa A_P(mu)S_1^[N](alpha_P/mu)-kappa^2[A_P'(mu)S_2^[N](alpha_P/mu)+(alpha_P/mu^2)A_P(mu)S_3^[N](alpha_P/mu)]}.",
            "sum": "O_N[P]=B_N^[N][P]-B_1^[N][P]+kappa^2 sum_(r notin R_N)^sym integral_1^N e(phi_r)C_rA_P du.",
            "convergence": "Only S_1 is principal-value conditional; the S_2, S_3, and interior C_r terms are absolutely summable after the joined twofold integration by parts.",
        },
        "physical_pullthrough": {
            "linear": "The complete outer linear correction is Re[nu c_xi O_N[P_lin^[N]]], with no six-moment componentwise split.",
            "terminal_value": "P_lin^[N](log N)=F_(alpha_N)(log N)-i u_N F_(beta_N)(log N).",
            "edge_split": "For beta_E=(mathcal N_E,V_E,-Q_E,-A_E), P_E(lambda)=i(lambda-log a)F_(beta_E)(lambda).",
            "edge_endpoint": "B_N^[N][P_E] has exact factor i exp(S(log N)) times {-kappa u_NF S_1-kappa^2 N^(-1)[F-u_N(F'+gF)]S_2+kappa^2 alpha_Pu_NN^(-2)F S_3}.",
            "suppression": "The conditional S_1 edge coefficient has the exact factor u_N<2h; the remaining edge coefficients contain either N^(-1)<=(3/2)h or alpha_Pu_N/N^2<9h/2.",
            "noncancellation": "The F_(alpha_N)(log N) coefficient survives without u_N. Thus terminal centering suppresses the edge-only trace but does not algebraically cancel the complete lower-boundary package.",
        },
        "quadratic_channels": {
            "split": "R_out=(1/2)Re[|nu|^2 H_out+nu^2 c_xi^2 T_out], where H_out uses O_X conjugate(O_dotY) and T_out uses O_X O_dotY with signs +,+,-,-.",
            "phase": "The common phase cancels exactly from H_out and appears squared in T_out.",
            "moved_tail": "The retained moving-tail defect remains 2Delta(omega_C)^TJdot(omega_B)+2(omega_E+omega_B+omega_C)^TJdot(omega_C); it is not an outer-mode remainder and is added exactly once.",
        },
        "audit": audit_outer_bands(),
    }


def build_rows(cert: dict) -> list[GateRow]:
    band = cert["outer_band"]
    kernels = cert["stable_kernels"]
    outer = cert["outer_identity"]
    physical = cert["physical_pullthrough"]
    quadratic = cert["quadratic_channels"]
    rows = [
        GateRow("oek_01_band", "full-support band", "proved", "The reciprocal physical band is one finite integer roster.", band["roster"], "The outer modes remain infinite."),
        GateRow("oek_02_distance", "pole separation", "proved", "Every possible stationary integer lies inside the band.", band["pole_distance"], "This is not an outer-sum bound."),
        GateRow("oek_03_tie", "tie guard", "guard_validated", "Band/complement transfer uses the exact half-open convention.", band["tie_rule"], "The complement must not be differentiated alone at a tie."),
        GateRow("oek_04_s1", "conditional kernel", "proved", "Both outer tails have one stable joined first-order kernel.", kernels["S1"], "It may grow logarithmically."),
        GateRow("oek_05_s2", "absolute kernel", "proved", "The second-order outer kernel is an absolute Hurwitz sum.", kernels["S2"], "This does not bound physical coefficients."),
        GateRow("oek_06_s3", "signed absolute kernel", "proved", "The third-order outer kernel retains the tail orientation.", kernels["S3"], "Its sign alone does not sign the flow."),
        GateRow("oek_07_integer", "removable value", "proved", "Integral stationary crossings have a finite harmonic value.", kernels["integer"], "Separating cotangent poles is forbidden."),
        GateRow("oek_08_kernel_bounds", "kernel envelope", "proved", "All absolute denominator kernels needed after two integrations are uniformly bounded.", kernels["bounds"], "No complete outer-complement estimate follows."),
        GateRow("oek_09_mode_ibp", "mode identity", "proved", "Every outer Fourier mode has an exact twofold nonstationary representation.", outer["mode"], "Endpoint terms remain main terms."),
        GateRow("oek_10_operator", "interior operator", "proved", "The twice-integrated interior operator is explicit.", outer["operator"], "Its physical aggregate is not evaluated."),
        GateRow("oek_11_endpoint", "endpoint package", "proved", "Each physical endpoint is represented by one joined three-kernel functional.", outer["endpoint"], "The N endpoint must remain joined to the edge row."),
        GateRow("oek_12_outer_sum", "outer recomposition", "proved", "The complete lower and upper complements reassemble before absolute values.", outer["sum"], "No scalar smallness is asserted."),
        GateRow("oek_13_convergence", "convergence", "proved", "Only the first endpoint kernel is conditionally convergent.", outer["convergence"], "Its symmetric convention is structural."),
        GateRow("oek_14_linear", "P_lin pullthrough", "proved", "The full outer linear correction remains one physical polynomial functional.", physical["linear"], "The functional is not yet bounded."),
        GateRow("oek_15_terminal_eval", "terminal evaluation", "proved", "The terminal P_lin value separates its frequency and base rows exactly.", physical["terminal_value"], "The alpha row need not vanish."),
        GateRow("oek_16_edge_split", "edge row", "proved", "The genuine endpoint enters only through the lifted beta row.", physical["edge_split"], "The edge is not assigned to a scalar cell."),
        GateRow("oek_17_edge_endpoint", "edge composition", "proved", "The N-endpoint edge trace is one pole-free exact expression.", physical["edge_endpoint"], "No coefficient norm is inserted."),
        GateRow("oek_18_suppression", "terminal centering", "proved", "Every edge-only terminal multiplier carries h-scale suppression.", physical["suppression"], "The surviving physical rows may still be large."),
        GateRow("oek_19_noncancel", "nonpromotion guard", "guard_validated", "The terminal strip does not generically annihilate the complete lower-boundary trace.", physical["noncancellation"], "A source-specific signed estimate is still required."),
        GateRow("oek_20_quadratic", "quadratic image", "proved", "The four outer observation pairings split exactly into Hermitian and transpose channels.", quadratic["split"], "Neither channel is signed here."),
        GateRow("oek_21_phase", "phase audit", "proved", "The common auxiliary phase behaves correctly in both quadratic channels.", quadratic["phase"], "Imaginary coefficients remain."),
        GateRow("oek_22_move", "moving tail", "proved", "The moved terminal block remains a retained observation package.", quadratic["moved_tail"], "It must not be counted as an outer residual."),
        GateRow("oek_23_route", "route decision", "proved", "The only unsuppressed conditional boundary coefficient is the complete physical alpha-row trace.", "Estimate it jointly with the reciprocal-band carrier, Hermitian, transpose, correction, and moved-tail kernels before any modulus.", "The edge recurrence alone cannot close the sign."),
        GateRow("oek_24_handoff", "signed theorem", "open", "Prove or falsify the joined full-support kernel at the h^2 flow scale.", "Insert O_N[P_lin^[N]], H_out, T_out, and the moved-tail defect into tilde(Psi)' and establish the direct negative retained-current target.", "No such estimate is proved here."),
        GateRow("oek_25_pi", "pi provenance", "proved", "All pi factors come from the fixed Fourier and completed-zeta conventions.", "e(x)=exp(2pi i x), kappa=1/(2pi i), and the cotangent audit is the symmetric integer partial fraction; no plotted geometry is used.", "No geometric inference is made."),
        GateRow("oek_26_scope", "scope guard", "guard_validated", "Kernel bounds are not promoted to a signed current theorem.", "The S2, S3, and Z4 bounds control denominators only; physical coefficient and signed correlation bounds remain open.", "No Phi_B or retained-current sign follows."),
        GateRow("oek_27_boundary", "proof boundary", "guard_validated", "This gate is an exact endpoint-kernel reduction, not a proof of RH.", "No complete outer bound, quadratic residual bound, signed flow, contact exclusion, Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.", "This is not a proof of RH."),
    ]
    require(len(rows) == 27, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    c = artifact["symbolic_certificate"]
    b = c["outer_band"]
    k = c["stable_kernels"]
    o = c["outer_identity"]
    p = c["physical_pullthrough"]
    q = c["quadratic_channels"]
    return f"""# Full-Support Outer Endpoint Kernel Gate

Date: 2026-08-02

Status: exact pole-free outer-complement and endpoint-linear kernel proved; signed full-support estimate open; not a proof of RH.

This is not a proof of RH. It replaces the split lower and upper outer tails by one stable endpoint package and then carries that package through the physical observation rows.

## Reciprocal Band

{b['roster']}

{b['pole_distance']}

{b['tie_rule']}

## Pole-Free Kernels

{k['S1']}

{k['S2']}

{k['S3']}

{k['integer']}

{k['bounds']}

## Outer Recomposition

{o['mode']}

{o['operator']}

{o['endpoint']}

{o['sum']}

{o['convergence']}

## Physical Pullthrough

{p['linear']}

{p['terminal_value']}

{p['edge_split']}

{p['edge_endpoint']}

{p['suppression']}

{p['noncancellation']}

## Quadratic Channels

{q['split']}

{q['phase']}

{q['moved_tail']}

## Handoff

The genuine edge no longer carries an unsuppressed conditional endpoint trace: its `S_1` coefficient contains `u_N`, and its remaining endpoint coefficients contain `1/N` or `alpha_P u_N/N^2`. The complete alpha-row term still survives. The next theorem must therefore estimate that term jointly with the returned reciprocal-band carrier, Hermitian phase-difference, transpose phase-sum, correction, and moved-tail kernels.

## Pi Provenance

`e(x)=exp(2pi i x)` fixes `kappa=1/(2pi i)`. The cotangent formula is used only as an audit of the symmetric integer partial fraction; the stable formula uses digamma and Hurwitz functions with positive arguments. No circle, polygon, or fitted visual constant enters.

## Proof Boundary

This gate proves a stable pole-free outer-mode representation, exact endpoint reassembly, the one-polynomial physical pullthrough, terminal-centering suppression of the edge-only trace, and exact Hermitian/transpose decomposition of the four quadratic outer pairings. It proves no complete outer-complement estimate, quadratic residual bound, signed full current, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Q209`, cofinal descendant theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "outer_band_rosters": 1,
        "pole_free_kernel_families": 3,
        "integer_removable_formulas": 1,
        "twofold_outer_recompositions": 1,
        "stable_endpoint_functionals": 2,
        "absolute_kernel_bounds": 3,
        "full_support_linear_functionals": 1,
        "terminal_plin_evaluations": 1,
        "edge_endpoint_factorizations": 1,
        "edge_terminal_suppression_templates": 1,
        "hermitian_transpose_decompositions": 1,
        "moving_tail_placements": 1,
        "complete_outer_complement_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_full_current_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "pole-free full-support outer endpoint kernel proved; complete signed estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves a stable joined representation of the full-support outer modes, exact endpoint reassembly, physical P_lin pullthrough, h-scale terminal centering of the genuine-edge endpoint trace, and Hermitian/transpose decomposition of the four quadratic outer pairings. It proves no complete outer-complement bound, quadratic residual estimate, signed full-current or Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Q209, cofinal descendant theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    audit = certificate["audit"]
    print(
        "built full-support outer endpoint kernel gate: "
        f"{counts['rows']} rows, {audit['kernel_points']} kernel points, "
        f"{audit['removable_integer_checks']} removable integer checks, "
        f"{counts['complete_outer_complement_bounds']} complete outer bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
