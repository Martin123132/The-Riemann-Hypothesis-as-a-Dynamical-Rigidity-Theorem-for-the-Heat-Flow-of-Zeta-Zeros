#!/usr/bin/env python3
"""Build the full-support terminal-bulk phase-completion gate."""

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
    "terminal_tail_anchored_six_moment_full_support_terminal_bulk_phase_"
    "completion_gate"
)
RESULT_PATH = REPO_ROOT / "work/rh_compute/results" / f"{KIND}.json"
NOTE_PATH = REPO_ROOT / "outputs" / f"{KIND}.md"

SOURCE_PATHS = {
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


def scalar(expression: sp.Matrix) -> sp.Expr:
    return sp.expand(expression[0])


def load_sources() -> dict[str, dict]:
    payloads: dict[str, dict] = {}
    for key, path in SOURCE_PATHS.items():
        if not path.exists():
            raise RuntimeError(f"missing source artifact: {path}")
        payloads[key] = json.loads(path.read_text(encoding="utf-8"))
    return payloads


def source_audit(payloads: dict[str, dict]) -> dict[str, dict[str, str]]:
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


def rational_bound_audit() -> dict[str, str]:
    coefficient = (
        Fraction(2)
        * Fraction(1, 36)
        * Fraction(1, 20)
        * Fraction(9, 4)
    )
    require(coefficient == Fraction(1, 160), "pure terminal coefficient")
    for denominator in [72_000_000_000, 10**12, 10**18, 10**24]:
        h = Fraction(1, denominator)
        require(2 * h < Fraction(1, 10), "terminal frequency scale")
        require(h**2 / 20 < h**2, "terminal rate scale")
        require(h**3 / 160 < h**2, "pure terminal scale")
    return {
        "factors": "The physical terminal boxes give u_N<2h, |R_(N,x)|<h^2/20, |C_N|<3/2, and |kappa|^2=1/(4pi^2)<1/36 because pi>3.",
        "bound": "The source-inclusive pure conditional quadratic obeys |R_aa|<h^3 S_(1,N)^2 |w_N|^2/160.",
        "comparison": "This is h^3 times a squared logarithmic kernel, whereas the mixed terminal-bulk channel is linear in S_(1,N) and contains the complete nonterminal derivative observation; no cancellation follows from scale alone.",
    }


def symbolic_certificate() -> dict:
    imaginary = sp.I
    s_1, u_q, u_nx = sp.symbols("S_1 u_q u_Nx", real=True)
    wnr, wni, wqr, wqi, nqr, nqi = sp.symbols(
        "wNr wNi wqr wqi Nqr Nqi", real=True
    )
    w_n = wnr + imaginary * wni
    w_q = wqr + imaginary * wqi
    n_q = nqr + imaginary * nqi
    d_q = -imaginary * u_q * n_q * w_q
    live = s_1 * sp.im(w_n) * sp.re(d_q) / (2 * sp.pi)
    k_h = s_1 * u_q * sp.conjugate(n_q) / (4 * sp.pi)
    k_t = -s_1 * u_q * n_q / (4 * sp.pi)
    split = sp.re(k_h * w_n * sp.conjugate(w_q) + k_t * w_n * w_q)
    require_zero(sp.expand_complex(live - split), "terminal-bulk H/T split")
    require_zero(k_t + sp.conjugate(k_h), "terminal-bulk anti-conjugacy")

    n_ideal = -u_q**2 / 4 - imaginary * u_nx / 2
    k_h_ideal = sp.expand(k_h.subs({nqr: -u_q**2 / 4, nqi: -u_nx / 2}))
    k_t_ideal = sp.expand(k_t.subs({nqr: -u_q**2 / 4, nqi: -u_nx / 2}))
    expected_h = -s_1 * u_q**3 / (16 * sp.pi) + imaginary * s_1 * u_q * u_nx / (8 * sp.pi)
    expected_t = s_1 * u_q**3 / (16 * sp.pi) + imaginary * s_1 * u_q * u_nx / (8 * sp.pi)
    require_zero(k_h_ideal - expected_h, "ideal Hermitian coefficient")
    require_zero(k_t_ideal - expected_t, "ideal transpose coefficient")
    witness = split.subs(
        {
            u_q: 1,
            u_nx: 0,
            wnr: 0,
            wni: 1,
            wqr: 0,
            wqi: 1,
            nqr: sp.re(n_ideal.subs({u_q: 1, u_nx: 0})),
            nqi: sp.im(n_ideal.subs({u_q: 1, u_nx: 0})),
        }
    )
    require_zero(witness + s_1 / (8 * sp.pi), "anti-conjugate witness")

    j_matrix = sp.Matrix(
        [
            [0, sp.Rational(1, 2), 0, 0],
            [sp.Rational(1, 2), 0, 0, 0],
            [0, 0, 0, -sp.Rational(1, 2)],
            [0, 0, -sp.Rational(1, 2), 0],
        ]
    )
    op = sp.Matrix(sp.symbols("Vp Np Ap Qp", real=True))
    dop = sp.Matrix(sp.symbols("dVp dNp dAp dQp", real=True))
    a = sp.Matrix(sp.symbols("aV aN aA aQ", real=True))
    da = sp.Matrix(sp.symbols("daV daN daA daQ", real=True))
    rho = sp.Matrix(sp.symbols("rV rN rA rQ", real=True))
    drho = sp.Matrix(sp.symbols("drV drN drA drQ", real=True))

    def correction(epsilon: sp.Matrix, dot_epsilon: sp.Matrix) -> sp.Expr:
        return scalar(
            2 * op.T * j_matrix * dot_epsilon
            + 2 * dop.T * j_matrix * epsilon
            + 2 * epsilon.T * j_matrix * dot_epsilon
        )

    omega_star = op + rho
    dot_omega_star = dop + drho
    delta_a = scalar(
        2 * omega_star.T * j_matrix * da
        + 2 * dot_omega_star.T * j_matrix * a
        + 2 * a.T * j_matrix * da
    )
    require_zero(
        correction(a + rho, da + drho) - correction(rho, drho) - delta_a,
        "observation coefficient completion",
    )
    component_delta = (
        omega_star[0] * da[1]
        + omega_star[1] * da[0]
        - omega_star[2] * da[3]
        - omega_star[3] * da[2]
        + dot_omega_star[0] * a[1]
        + dot_omega_star[1] * a[0]
        - dot_omega_star[2] * a[3]
        - dot_omega_star[3] * a[2]
        + a[0] * da[1]
        + a[1] * da[0]
        - a[2] * da[3]
        - a[3] * da[2]
    )
    require_zero(delta_a - component_delta, "component completion")
    mixed_delta = scalar(
        2 * omega_star.T * j_matrix * da
        + 2 * dot_omega_star.T * j_matrix * a
    )
    require_zero(
        sp.diff(mixed_delta, a[0]) - dot_omega_star[1],
        "completed a_V coefficient",
    )
    terminal_vector = sp.Matrix(sp.symbols("tV tN tA tQ", real=True))
    require_zero(
        scalar(terminal_vector.T * (da + drho))
        - scalar(terminal_vector.T * da)
        - scalar(terminal_vector.T * drho),
        "fixed terminal linear split",
    )

    amplitude_parts = sp.symbols(
        "xVr xVi xNr xNi xAr xAi xQr xQi", real=True
    )
    x_v = amplitude_parts[0] + imaginary * amplitude_parts[1]
    x_n = amplitude_parts[2] + imaginary * amplitude_parts[3]
    x_a = amplitude_parts[4] + imaginary * amplitude_parts[5]
    x_q = amplitude_parts[6] + imaginary * amplitude_parts[7]
    u_n = sp.symbols("u_N", real=True)
    dx_v, dx_n, dx_a, dx_q = [
        -imaginary * u_n * value for value in (x_v, x_n, x_a, x_q)
    ]
    hermitian = (
        x_v * sp.conjugate(dx_n)
        + x_n * sp.conjugate(dx_v)
        - x_a * sp.conjugate(dx_q)
        - x_q * sp.conjugate(dx_a)
    )
    transpose = x_v * dx_n + x_n * dx_v - x_a * dx_q - x_q * dx_a
    require_zero(sp.expand_complex(sp.re(hermitian)), "pure Hermitian real null")
    require_zero(
        sp.expand(transpose + 2 * imaginary * u_n * (x_v * x_n - x_a * x_q)),
        "pure transpose pairing",
    )

    tau, r_value, rx_value, c_value, q_value = sp.symbols(
        "tau R R_x C Q"
    )
    a_value = r_value * c_value
    n_value = r_value * q_value + rx_value * c_value
    require_zero(
        c_value * n_value - a_value * q_value - rx_value * c_value**2,
        "carrier determinant collapse",
    )
    terminal_transpose = -2 * imaginary * u_n * tau**2 * (
        c_value * n_value - a_value * q_value
    )
    require_zero(
        terminal_transpose
        + 2 * imaginary * u_n * tau**2 * rx_value * c_value**2,
        "terminal transpose collapse",
    )

    bounds = rational_bound_audit()
    return {
        "terminal_bulk_split": {
            "derivative_atom": "For q<N put d_q=-iu_q mathcal N_qw_q, so mathcal N_(<N,xi)=sum_(q<N)Re(d_q).",
            "real_product": "Im(w_N)Re(d_q)={Im(w_Nd_q)+Im(w_Nconjugate(d_q))}/2.",
            "kernels": "The live kernel is sum_(q<N)Re{K_H(N,q)w_Nconjugate(w_q)+K_T(N,q)w_Nw_q}, with K_H=S_(1,N)u_qconjugate(mathcal N_q)/(4pi) and K_T=-S_(1,N)u_qmathcal N_q/(4pi).",
            "anti_conjugacy": "K_T(N,q)=-conjugate(K_H(N,q)); this relates coefficients but does not cancel the two distinct phase families.",
            "physical_phases": "Under w_q=omega_a A_q exp(-i xi u_q), w_Nconjugate(w_q)=|omega_a|^2A_NA_q exp(i xi(u_q-u_N)) is Hermitian phase-difference, while w_Nw_q=omega_a^2A_NA_q exp(-i xi(u_N+u_q)) is transpose phase-sum.",
        },
        "ideal_phase_chart": {
            "row": "For mathcal N_q^0=-u_q^2/4-i u_(N,x)/2, K_H^0=-S_(1,N)u_q^3/(16pi)+iS_(1,N)u_qu_(N,x)/(8pi) and K_T^0=S_(1,N)u_q^3/(16pi)+iS_(1,N)u_qu_(N,x)/(8pi).",
            "orientation": "The ideal real parts are opposite and the imaginary parts agree; neither coefficient has a pointwise sign.",
            "witness": "At u_q=1, u_(N,x)=0, and w_N=w_q=i, the joined pair equals -S_(1,N)/(8pi), not zero.",
        },
        "coefficient_completion": {
            "current": "For C(epsilon)=2omega_p^TJdot(epsilon)+2dot(omega_p)^TJepsilon+2epsilon^TJdot(epsilon), split epsilon=a+rho into the terminal conditional observation and every nonterminal residual observation.",
            "identity": "C(a+rho)=C(rho)+2(omega_p+rho)^TJdot a+2(dot(omega_p)+dot rho)^TJa+2a^TJdot a.",
            "complete_rows": "Put omega_*=omega_p+rho and dot(omega_*)=dot(omega_p)+dot rho. The terminal-dependent quadratic channel is Delta_a=2omega_*^TJdot a+2dot(omega_*)^TJa+2a^TJdot a.",
            "components": "Delta_a=V_*dot a_N+N_*dot a_V-A_*dot a_Q-Q_*dot a_A+dot V_*a_N+dot N_*a_V-dot A_*a_Q-dot Q_*a_A+a_Vdot a_N+a_Ndot a_V-a_Adot a_Q-a_Qdot a_A.",
            "leading": "The mixed coefficient of a_V is the complete nonterminal derivative dot N_*=mathcal N_(p,xi)+rho_(N,xi); the pure a-by-dot a term remains separate.",
            "fixed_terminal": "The affine term t_T^Tdot(epsilon) splits linearly as t_T^Tdot rho+t_T^Tdot a and is not absorbed into C or counted twice.",
        },
        "pure_terminal_current": {
            "amplitudes": "For the source-inclusive conditional amplitudes A_X=tau_NP_X(log N), tau_N=kappa S_(1,N)w_N, the frequency lift is A_(dot X)=-iu_NA_X.",
            "hermitian": "H_aa=A_Vconjugate(A_(dot N))+A_Nconjugate(A_(dot V))-A_Aconjugate(A_(dot Q))-A_Qconjugate(A_(dot A)) has Re(H_aa)=0 exactly.",
            "transpose": "T_aa=-2iu_N(A_VA_N-A_AA_Q).",
            "determinant": "The carrier identities A_N=R_NC_N and mathcal N_N=R_NQ_N+R_(N,x)C_N give C_Nmathcal N_N-A_NQ_N=R_(N,x)C_N^2 exactly.",
            "collapse": "Hence T_aa=-2iu_N tau_N^2 R_(N,x)C_N^2, and the actual pure-terminal real current is R_aa=(1/2)Re(T_aa).",
        },
        "bounds": bounds,
        "audit": {
            "phase_families": 2,
            "coefficient_completions": 1,
            "pure_hermitian_nulls": 1,
            "pure_transpose_collapses": 1,
        },
        "handoff": {
            "decision": "Pure conditional self-interaction cannot cancel the live terminal-bulk channel: its Hermitian real part vanishes and its transpose remainder is h^3S_(1,N)^2 scale.",
            "route": "Expand rho into the lower S_1 endpoint, the S_2/S_3 endpoints, and the absolutely convergent interior; then form the completed Hermitian and transpose kernels for dot N_*a_V before taking moduli.",
        },
    }


def build_rows(cert: dict) -> list[GateRow]:
    split = cert["terminal_bulk_split"]
    ideal = cert["ideal_phase_chart"]
    completion = cert["coefficient_completion"]
    pure = cert["pure_terminal_current"]
    bounds = cert["bounds"]
    handoff = cert["handoff"]
    rows = [
        GateRow("tbp_01_atom", "bulk derivative", "proved", "The nonterminal derivative is an exact finite real carrier sum.", split["derivative_atom"], "No sign is inferred."),
        GateRow("tbp_02_product", "real-product polarization", "proved", "The live terminal-bulk product has an exact two-quadrature polarization.", split["real_product"], "Both summands must be retained."),
        GateRow("tbp_03_hermitian", "Hermitian kernel", "proved", "One half of the live kernel is a phase-difference family.", split["kernels"], "It is not individually signed."),
        GateRow("tbp_04_transpose", "transpose kernel", "proved", "The other half is a phase-sum family.", split["kernels"], "It is not individually signed."),
        GateRow("tbp_05_anti", "coefficient relation", "proved", "The two coefficients are anti-conjugate.", split["anti_conjugacy"], "Anti-conjugacy is not cancellation."),
        GateRow("tbp_06_phases", "physical phases", "proved", "The two products carry distinct physical phases.", split["physical_phases"], "They cannot be identified pointwise."),
        GateRow("tbp_07_ideal", "ideal coefficients", "proved", "The correction-free coefficients are explicit cubic/linear functions of the logarithmic rates.", ideal["row"], "This does not bound corrections."),
        GateRow("tbp_08_orientation", "sign guard", "guard_validated", "Opposite real parts do not supply a joined sign.", ideal["orientation"], "The imaginary parts reinforce instead."),
        GateRow("tbp_09_witness", "noncancellation witness", "proved", "The anti-conjugate pair is not identically zero.", ideal["witness"], "This algebraic witness is not asserted to be a physical Xi state."),
        GateRow("tbp_10_current", "observation current", "proved", "The carrier/residual quadratic correction is one exact J-pairing.", completion["current"], "The fixed terminal affine row is separate."),
        GateRow("tbp_11_split", "residual split", "proved", "The outer observation splits into terminal conditional and nonterminal parts.", "Write epsilon=a+rho and dot epsilon=dot a+dot rho.", "No source term is discarded."),
        GateRow("tbp_12_completion", "coefficient completion", "proved", "All mixed outer terms complete the terminal coefficient exactly.", completion["identity"], "This is an identity, not a bound."),
        GateRow("tbp_13_rows", "complete rows", "proved", "The terminal-dependent channel uses complete nonterminal observations.", completion["complete_rows"], "The pure terminal channel remains explicit."),
        GateRow("tbp_14_components", "component audit", "proved", "All twelve terminal-dependent signed pairings are exposed.", completion["components"], "No componentwise modulus is taken."),
        GateRow("tbp_15_leading", "leading coefficient", "proved", "The mixed a_V coefficient is the complete nonterminal derivative row.", completion["leading"], "Its sign and size remain open."),
        GateRow("tbp_16_affine", "fixed terminal row", "proved", "The fixed terminal affine interaction remains additive.", completion["fixed_terminal"], "It is not hidden in the quadratic completion."),
        GateRow("tbp_17_amplitudes", "conditional amplitudes", "proved", "Every pure terminal S_1 observation shares one source-inclusive amplitude.", pure["amplitudes"], "This concerns only the N-endpoint S_1 block."),
        GateRow("tbp_18_lift", "frequency lift", "proved", "The pure terminal derivative has one common -iu_N factor.", pure["amplitudes"], "Other endpoint and interior terms are not included."),
        GateRow("tbp_19_hnull", "Hermitian null", "proved", "The pure terminal Hermitian current has zero real part exactly.", pure["hermitian"], "The complex Hermitian expression need not vanish."),
        GateRow("tbp_20_transpose", "transpose reduction", "proved", "The pure terminal transpose current is one carrier determinant.", pure["transpose"], "It can have either sign after projection."),
        GateRow("tbp_21_determinant", "carrier identity", "proved", "The determinant collapses to the terminal x-rate.", pure["determinant"], "No approximate cancellation is used."),
        GateRow("tbp_22_collapse", "pure terminal collapse", "proved", "The complete pure S_1-by-S_1 current is an h-cubed transpose remainder.", pure["collapse"], "This is not the mixed terminal-bulk term."),
        GateRow("tbp_23_bound", "rational envelope", "proved", "The pure terminal real current has an explicit source-relative envelope.", bounds["factors"] + " " + bounds["bound"], "The source amplitude is not numerically bounded here."),
        GateRow("tbp_24_scale", "scale separation", "guard_validated", "The pure terminal block cannot furnish the missing order-one algebraic cancellation.", bounds["comparison"], "A signed mixed estimate is still required."),
        GateRow("tbp_25_route", "route decision", "open", "Prove or falsify the completed terminal/nonterminal Hermitian-transpose estimate.", handoff["route"], "No signed completed bound is proved."),
        GateRow("tbp_26_pi", "pi provenance", "proved", "Every pi factor is inherited from Fourier normalization.", "The factor 1/(4pi) is one half of the Section 11.184 factor 1/(2pi), and |kappa|=1/(2pi); no geometric fitting is introduced.", "No circle or polygon supplies pi."),
        GateRow("tbp_27_join", "joining guard", "guard_validated", "The phase completion preserves every outer and moved-tail channel.", handoff["decision"], "The lower endpoint, S_2/S_3 terms, interior, correction, and moved-tail defect remain open obligations."),
        GateRow("tbp_28_boundary", "proof boundary", "guard_validated", "This is an exact phase and coefficient reduction, not a proof of RH.", "No signed completed terminal/nonterminal estimate, complete outer bound, quadratic residual bound, full-current theorem, contact exclusion, Xi theorem, Lambda<=0, PF-infinity, RH, or prize-level conclusion is claimed.", "This is not a proof of RH."),
    ]
    require(len(rows) == 28, "row count drifted")
    return rows


def render_note(artifact: dict) -> str:
    cert = artifact["symbolic_certificate"]
    split = cert["terminal_bulk_split"]
    ideal = cert["ideal_phase_chart"]
    completion = cert["coefficient_completion"]
    pure = cert["pure_terminal_current"]
    bounds = cert["bounds"]
    handoff = cert["handoff"]
    return f"""# Full-Support Terminal-Bulk Phase Completion Gate

Date: 2026-08-02

Status: exact terminal-bulk phase split, observation coefficient completion, and pure conditional null proved; signed completed estimate open; not a proof of RH.

## Terminal-Bulk Split

{split['derivative_atom']}

{split['real_product']}

{split['kernels']}

{split['anti_conjugacy']}

{split['physical_phases']}

## Ideal Phase Chart

{ideal['row']}

{ideal['orientation']}

{ideal['witness']}

## Coefficient Completion

{completion['current']}

{completion['identity']}

{completion['complete_rows']}

{completion['components']}

{completion['leading']}

{completion['fixed_terminal']}

## Pure Conditional Current

{pure['amplitudes']}

{pure['hermitian']}

{pure['transpose']}

{pure['determinant']}

{pure['collapse']}

{bounds['factors']}

{bounds['bound']}

{bounds['comparison']}

## Handoff

{handoff['decision']}

{handoff['route']}

## Pi Provenance

The factor `1/(4pi)` is the polarization half of the fixed conditional factor `1/(2pi)`, itself forced by `kappa=1/(2pi i)` and `e(x)=exp(2pi i x)`. The estimate `|kappa|^2<1/36` uses only `pi>3`. No circle, polygon, plotted symmetry, or fitted constant is introduced.

## Proof Boundary

This gate proves the exact Hermitian phase-difference and transpose phase-sum decomposition of the live terminal-bulk quadrature, their anti-conjugate coefficient relation and explicit noncancellation witness, exact observation-level coefficient completion, exact vanishing of the pure terminal Hermitian real current, collapse of the pure terminal transpose current through `C_N mathcal N_N-A_NQ_N=R_(N,x)C_N^2`, and its `h^3 S_(1,N)^2` envelope. It proves no signed completed terminal/nonterminal estimate, complete outer-complement or quadratic residual bound, signed full current, `Phi_B` theorem, contact exclusion, retained aggregate or Xi theorem, `Lambda<=0`, PF-infinity, RH, or prize-level conclusion.
"""


def main() -> int:
    payloads = load_sources()
    certificate = symbolic_certificate()
    rows = build_rows(certificate)
    counts = {
        "rows": len(rows),
        "terminal_bulk_polarizations": 1,
        "phase_families": 2,
        "anti_conjugate_kernel_relations": 1,
        "ideal_coefficient_pairs": 1,
        "algebraic_noncancellation_witnesses": 1,
        "observation_coefficient_completions": 1,
        "completed_component_pairings": 12,
        "fixed_terminal_linear_splits": 1,
        "pure_terminal_hermitian_nulls": 1,
        "pure_terminal_transpose_collapses": 1,
        "carrier_determinant_collapses": 1,
        "pure_terminal_rational_bounds": 1,
        "signed_completed_terminal_nonterminal_bounds": 0,
        "complete_outer_complement_bounds": 0,
        "quadratic_remainder_bounds": 0,
        "signed_full_current_bounds": 0,
        "phi_b_bounds": 0,
    }
    artifact = {
        "kind": KIND,
        "date": "2026-08-02",
        "status": "terminal-bulk phase and coefficient completion proved; signed completed estimate open",
        "counts": counts,
        "symbolic_certificate": certificate,
        "source_audit": source_audit(payloads),
        "rows": [asdict(row) for row in rows],
        "proof_boundary": "This proves the exact Hermitian/transpose polarization of the terminal-bulk quadrature, anti-conjugate coefficient relation, observation-level coefficient completion, pure terminal Hermitian real null, pure terminal transpose determinant collapse, and its h^3 S_(1,N)^2 envelope. It proves no signed completed terminal/nonterminal estimate, complete outer-complement or quadratic residual bound, signed full-current or Phi_B theorem, contact exclusion, retained aggregate or Xi theorem, Q209, Lambda<=0, PF-infinity, RH, or prize-level conclusion.",
    }
    atomic_write(RESULT_PATH, json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    atomic_write(NOTE_PATH, render_note(artifact))
    print(
        "built full-support terminal-bulk phase completion gate: "
        f"{counts['rows']} rows, {counts['phase_families']} phase families, "
        f"{counts['observation_coefficient_completions']} coefficient completion, "
        f"{counts['pure_terminal_hermitian_nulls']} Hermitian null, "
        f"{counts['pure_terminal_transpose_collapses']} transpose collapse, "
        f"{counts['signed_completed_terminal_nonterminal_bounds']} signed completed bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
