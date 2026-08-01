#!/usr/bin/env python3
"""Build the endpoint Schur-Cohn stability and first-jet guard."""

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
    "endpoint_schur_cohn_first_jet_guard"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "phase_cylinder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "phase_cylinder_jacobian_transport_guard.json"
    ),
    "joined_odd_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "joined_dyadic_odd_prefix_first_jet_reduction.json"
    ),
    "adjacent_recurrence": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_saddle_recurrence.json"
    ),
    "normalized_prefix": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "normalized_prefix_phase_flux_reduction.json"
    ),
}


@dataclass(frozen=True)
class SchurRow:
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
        "phase_cylinder": (
            "F_N(z,tau)=P_N(z,omega_*,1)+tau*r_0",
            "three-cylinder crossing budget",
        ),
        "joined_odd_prefix": (
            "P_N(z_*,omega_*,1)=H_0",
            "H_0,H_1,H_2,D_0,D_1",
            "Delta Z_0=q_n+j_0=Q_N/f_1",
        ),
        "adjacent_recurrence": (
            "Q_N=f_(N+1)+kappa_N*J_a",
            "Delta A_a=Re[-s_*'*ell*f_(N+1)",
        ),
        "normalized_prefix": (
            "Gamma_ell=mathsf_X+i*mathsf_A/ell",
            "partial_x arg(Gamma_ell)",
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


def conjugate_square(value: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.conjugate(value) * value)


def reverse_conjugate(coeffs: list[sp.Expr]) -> list[sp.Expr]:
    return [sp.conjugate(value) for value in reversed(coeffs)]


def outside_step(coeffs: list[sp.Expr]) -> list[sp.Expr]:
    top = coeffs[-1]
    constant_bar = sp.conjugate(coeffs[0])
    reversed_coeffs = reverse_conjugate(coeffs)
    transformed = [
        sp.expand(constant_bar * value - top * reverse)
        for value, reverse in zip(coeffs, reversed_coeffs)
    ]
    if sp.simplify(transformed[-1]) != 0:
        raise RuntimeError("outside Schur step did not cancel top degree")
    return transformed[:-1]


def schur_step_audit() -> dict[str, str]:
    z = sp.symbols("z")
    coeffs = list(sp.symbols("b_0:5"))
    f = sum(value * z**index for index, value in enumerate(coeffs))
    transformed = outside_step(coeffs)
    t = sum(value * z**index for index, value in enumerate(transformed))
    t_star = sum(
        value * z**index
        for index, value in enumerate(reverse_conjugate(transformed))
    )
    delta = conjugate_square(coeffs[0]) - conjugate_square(coeffs[-1])
    inverse_residual = sp.expand(
        coeffs[0] * t + coeffs[-1] * z * t_star - delta * f
    )
    for index in range(5):
        if sp.simplify(inverse_residual.coeff(z, index)) != 0:
            raise RuntimeError("outside Schur inverse identity failed")

    alpha = sp.symbols("alpha")
    middle = list(sp.symbols("a_1:4"))
    normalized = [sp.Integer(1), *middle, alpha]
    normalized_star = reverse_conjugate(normalized)
    denominator = 1 - conjugate_square(alpha)
    reduced = [
        sp.simplify(
            (value - alpha * reverse) / denominator
        )
        for value, reverse in zip(
            normalized[:-1], normalized_star[:-1]
        )
    ]
    if reduced[0] != 1:
        raise RuntimeError("normalized Schur step lost constant one")
    reduced_star = reverse_conjugate(reduced)
    reconstructed = [
        reduced[index]
        if index < len(reduced)
        else sp.Integer(0)
        for index in range(len(normalized))
    ]
    for index, value in enumerate(reduced_star, start=1):
        reconstructed[index] += alpha * value
    for actual, expected in zip(reconstructed, normalized):
        if sp.simplify(actual - expected) != 0:
            raise RuntimeError("normalized Schur inverse failed")

    return {
        "reverse": (
            "For degree d, F#(z)=z^d*conj(F(1/conj(z)))="
            "sum_(k=0)^d conj(b_(d-k))*z^k."
        ),
        "step": (
            "S_dF(z)=conj(b_0)F(z)-b_dF#(z)="
            "sum_(k=0)^(d-1)"
            "[conj(b_0)b_k-b_dconj(b_(d-k))]z^k."
        ),
        "pivot": "Delta_d=|b_0|^2-|b_d|^2.",
        "inverse": (
            "Delta_d F(z)=b_0 S_dF(z)+"
            "b_d z (S_dF)#(z), where the latter # has degree d-1."
        ),
        "normalization": (
            "For p_d(0)=1 and alpha_d=[z^d]p_d, "
            "p_(d-1)=(p_d-alpha_d p_d#)/(1-|alpha_d|^2), "
            "and p_d=p_(d-1)+alpha_d z p_(d-1)#."
        ),
    }


def frechet_audit() -> dict[str, str]:
    degree = 3
    x, nu = sp.symbols("x nu", real=True)
    b = list(sp.symbols("b_0:4"))
    u = list(sp.symbols("u_0:4"))
    v = list(sp.symbols("v_0:4"))
    w = list(sp.symbols("w_0:4"))

    def differential(base: list[sp.Expr], tangent: list[sp.Expr]) -> list[sp.Expr]:
        return [
            sp.expand(
                sp.conjugate(tangent[0]) * base[k]
                + sp.conjugate(base[0]) * tangent[k]
                - tangent[degree] * sp.conjugate(base[degree - k])
                - base[degree] * sp.conjugate(tangent[degree - k])
            )
            for k in range(degree)
        ]

    def bilinear(
        first: list[sp.Expr], second: list[sp.Expr]
    ) -> list[sp.Expr]:
        return [
            sp.expand(
                sp.conjugate(first[0]) * second[k]
                + sp.conjugate(second[0]) * first[k]
                - first[degree] * sp.conjugate(second[degree - k])
                - second[degree] * sp.conjugate(first[degree - k])
            )
            for k in range(degree)
        ]

    perturbed = [
        b[k] + x * u[k] + nu * v[k] + x * nu * w[k]
        for k in range(degree + 1)
    ]
    transformed = outside_step(perturbed)
    x_actual = [
        sp.diff(value, x).subs({x: 0, nu: 0})
        for value in transformed
    ]
    nu_actual = [
        sp.diff(value, nu).subs({x: 0, nu: 0})
        for value in transformed
    ]
    mixed_actual = [
        sp.diff(value, x, nu).subs({x: 0, nu: 0})
        for value in transformed
    ]
    x_expected = differential(b, u)
    nu_expected = differential(b, v)
    mixed_expected = [
        sp.expand(a + c)
        for a, c in zip(differential(b, w), bilinear(u, v))
    ]
    for actual, expected in zip(x_actual, x_expected):
        if sp.simplify(actual - expected) != 0:
            raise RuntimeError("Schur x tangent audit failed")
    for actual, expected in zip(nu_actual, nu_expected):
        if sp.simplify(actual - expected) != 0:
            raise RuntimeError("Schur companion tangent audit failed")
    for actual, expected in zip(mixed_actual, mixed_expected):
        if sp.simplify(actual - expected) != 0:
            raise RuntimeError("Schur mixed tangent audit failed")

    return {
        "differential": (
            "For V=sum v_kz^k, "
            "(DS_b[V])_k=conj(v_0)b_k+conj(b_0)v_k"
            "-v_dconj(b_(d-k))-b_dconj(v_(d-k))."
        ),
        "bilinear": (
            "D2S_b[U,V]_k=conj(u_0)v_k+conj(v_0)u_k"
            "-u_dconj(v_(d-k))-v_dconj(u_(d-k))."
        ),
        "mixed": (
            "For b_(x,nu)=b+xU+nuV+xnuW, "
            "partial_(x,nu)S(b_(x,nu))|_0="
            "DS_b[W]+D2S_b[U,V]."
        ),
    }


def degree_guards() -> dict[str, dict[str, str]]:
    z = sp.symbols("z")

    stable_linear = [sp.Integer(3), sp.Integer(2)]
    stable_root = sp.solve(stable_linear[0] + stable_linear[1] * z, z)[0]
    if stable_root != -sp.Rational(3, 2):
        raise RuntimeError("degree-one outside root audit failed")
    stable_pivot = conjugate_square(stable_linear[0]) - conjugate_square(
        stable_linear[1]
    )
    if stable_pivot != 5:
        raise RuntimeError("degree-one pivot audit failed")

    quadratic = [
        sp.Integer(1),
        -sp.Rational(9, 4),
        sp.Rational(1, 2),
    ]
    roots = set(sp.solve(sum(v * z**k for k, v in enumerate(quadratic)), z))
    if roots != {sp.Rational(1, 2), sp.Integer(4)}:
        raise RuntimeError("degree-two roots audit failed")
    first_reduction = outside_step(quadratic)
    first_pivot = conjugate_square(quadratic[0]) - conjugate_square(
        quadratic[-1]
    )
    second_pivot = conjugate_square(
        first_reduction[0]
    ) - conjugate_square(first_reduction[-1])
    if first_pivot != sp.Rational(3, 4):
        raise RuntimeError("degree-two first pivot audit failed")
    if second_pivot != -sp.Rational(45, 64):
        raise RuntimeError("degree-two second pivot audit failed")

    sqrt_55 = sp.sqrt(55)
    phased = [
        sp.Integer(1),
        (-23 + 3 * sp.I * sqrt_55) / 40,
        (-17 - 3 * sp.I * sqrt_55) / 40,
    ]
    phased_norms = [sp.simplify(conjugate_square(value)) for value in phased]
    if phased_norms != [
        sp.Integer(1),
        sp.Rational(16, 25),
        sp.Rational(49, 100),
    ]:
        raise RuntimeError("decreasing-phase norm audit failed")
    if sp.simplify(sum(phased)) != 0:
        raise RuntimeError("decreasing-phase boundary zero audit failed")
    phased_reduction = outside_step(phased)
    if sp.simplify(phased_reduction[0] + phased_reduction[1]) != 0:
        raise RuntimeError("decreasing-phase reduced boundary zero failed")

    endpoint = [
        -sp.Rational(2, 3),
        sp.Rational(2, 3),
    ]
    endpoint_root = sp.solve(endpoint[0] + endpoint[1] * z, z)[0]
    endpoint_pivot = conjugate_square(endpoint[0]) - conjugate_square(
        endpoint[1]
    )
    if endpoint_root != 1 or endpoint_pivot != 0:
        raise RuntimeError("endpoint pivot guard failed")

    return {
        "degree_one_orientation": {
            "polynomial": "F(z)=3+2z",
            "root": "z=-3/2",
            "pivot": "Delta_1=9-4=5>0",
            "interpretation": (
                "The constant-dominant orientation corresponds to a root "
                "outside the closed unit disk."
            ),
        },
        "first_pivot_not_sufficient": {
            "polynomial": "F(z)=1-(9/4)z+(1/2)z^2",
            "roots": "z=1/2 and z=4",
            "first_pivot": "Delta_2=3/4>0",
            "reduced_polynomial": "S_2F(z)=3/4-(9/8)z",
            "second_pivot": "Delta_1=-45/64<0",
            "interpretation": (
                "A strict first pivot is necessary but not sufficient; "
                "every recursive pivot is needed for full disk exclusion."
            ),
        },
        "decreasing_modulus_phase_guard": {
            "coefficients": (
                "b_0=1, b_1=(-23+3i*sqrt(55))/40, "
                "b_2=(-17-3i*sqrt(55))/40"
            ),
            "moduli": "|b_0|=1>|b_1|=4/5>|b_2|=7/10",
            "boundary_zero": "F(1)=b_0+b_1+b_2=0",
            "first_pivot": "Delta_2=51/100>0",
            "interpretation": (
                "Strictly decreasing coefficient moduli do not replace "
                "Schur-Cohn once the coefficient phases are free."
            ),
        },
        "endpoint_pivot_guard": {
            "prefix": "P(z)=1+(2/3)z",
            "endpoint": "r_0=-5/3",
            "augmented": "F(z)=-2/3+(2/3)z",
            "boundary_zero": "F(1)=0",
            "pivot": "Delta_1=0",
            "interpretation": (
                "The endpoint enters the first pivot and can place a zero "
                "on the unit circle even when the prefix is disk-zero-free."
            ),
        },
    }


def build_exact() -> dict:
    algebra = schur_step_audit()
    tangents = frechet_audit()
    guards = degree_guards()
    return {
        "pi_provenance": (
            "The pi in a^2=x/(4*pi)+t/16 comes from the completed-zeta "
            "normalization and the Riemann-Siegel saddle; 4*pi*(2N+1) "
            "is the cutoff-cell difference and 2*pi is the period of "
            "exp(i*theta). Schur-Cohn uses the unit circle |z|=1 and "
            "introduces no new pi."
        ),
        "endpoint_polynomial": (
            "At xi=omega_* and mu=1 put "
            "B_(k,r)=sum_(m<=M_k,m odd)(log(2^k m))^r "
            "a_(k,m)exp(i omega_*log m)c_(2^k m), "
            "E_(k,r)=sum_(m<=M_k,m odd)(log(2^k m))^r "
            "delta_(2^k m)a_(k,m)exp(i omega_*log m)c_(2^k m). "
            "Then P_r(z)=sum_k B_(k,r)z^k, D_r(z)=sum_k E_(k,r)z^k, "
            "and at z_*=exp(i omega_*log 2), P_r(z_*)=H_r and "
            "D_r(z_*)=D_r. Define F_N(z)=r_0+P_0(z)=sum_(k=0)^K b_kz^k "
            "with b_0=r_0+B_(0,0) and b_k=B_(k,0) for k>=1."
        ),
        "reverse": algebra["reverse"],
        "outside_step": algebra["step"],
        "pivot": algebra["pivot"],
        "inverse": algebra["inverse"],
        "outside_stability_theorem": (
            "After trimming exact zero top coefficients, F_d has no zero "
            "in |z|<=1 iff Delta_d>0 and S_dF_d has no zero in |z|<=1. "
            "Indeed, on |z|=1 one has |F_d#|=|F_d|. Rouche applied to "
            "S_dF_d=conj(b_0)F_d-b_dF_d# gives the forward zero count; "
            "the inverse identity and |z(S_dF_d)#|=|S_dF_d| give the "
            "reverse count. Iteration ends at a nonzero constant."
        ),
        "normalized_reflection_recursion": algebra["normalization"],
        "reflection_criterion": (
            "Normalize each effective degree-d polynomial by its nonzero "
            "constant and write alpha_d for its leading coefficient. "
            "Then F_N is zero-free in |z|<=1 iff |alpha_d|<1 at every "
            "effective stage d=K,...,1. These alpha_d are the exact "
            "outside-disk reflection coefficients."
        ),
        "boundary_margin": (
            "The inverse normalized recursion gives, for |z|=1, "
            "|p_d(z)|>=(1-|alpha_d|)|p_(d-1)(z)|. Consequently "
            "min_(|z|=1)|F_N(z)|>=|b_0|*product_(d=1)^K"
            "(1-|alpha_d|)>0 whenever every pivot is strict. This is a "
            "quantified sufficient boundary margin, not a claim that the "
            "actual Xi reflection coefficients satisfy the inequalities."
        ),
        "coefficient_x_tangent": (
            "Inside one fixed-N chart, "
            "(B_(k,r))_x=-s_*'B_(k,r+1)+E_(k,r)"
            "-i*k*log(2)*omega_*'B_(k,r). Thus the fixed-z tangent is "
            "X(z)=F_(N,x)|_z=r_(0,x)-s_*'P_1(z)+D_0(z)"
            "-i*omega_*'log(2)zP_0'(z). At z=z_*, adding "
            "z_*'=i*omega_*'log(2)z_* gives "
            "X(z_*)+z_*'F_N'(z_*)=r_(0,x)-s_*'H_1+D_0=Z_(0,x)."
        ),
        "centered_companion": (
            "Let A=log(a) and "
            "G(z)=r_A-s_*'[P_1(z)-A P_0(z)]. Then G(z_*)=Z_A. Its "
            "fixed-z x tangent is W(z)=r_(A,x)-s_*''[P_1-AP_0]"
            "-s_*'{-s_*'[P_2-AP_1]+D_1-AD_0-A_xP_0"
            "-i*omega_*'log(2)z[P_1'-AP_0']}. After adding z_*'G'(z_*), "
            "the phase terms cancel and one recovers the exact "
            "Z_(A,x) formula containing H_2,D_0,D_1,r_(A,x)."
        ),
        "frechet_recursion": tangents["differential"],
        "frechet_bilinear": tangents["bilinear"],
        "mixed_recursion": tangents["mixed"],
        "tangent_application": (
            "Apply DS at every strict Schur stage to U=X and V=G. Apply "
            "DS[W]+D2S[X,G] to the mixed x/companion tangent. This "
            "propagates Z_0,Z_A,Z_(0,x),Z_(A,x) through the same exact "
            "coefficient recursion before any modulus or triangle "
            "inequality is taken."
        ),
        "first_xi_pivot": (
            "The first physical pivot is "
            "Delta_K=|r_0+B_(0,0)|^2-|B_(K,0)|^2, equivalently "
            "|alpha_K|=|B_(K,0)/(r_0+B_(0,0))|<1. Phase-frozen dyadic "
            "contraction compares neighboring positive B_k before odd "
            "phases, d_n corrections, and endpoint addition. It does not "
            "prove this endpoint-sensitive inequality."
        ),
        "degree_one_orientation": guards["degree_one_orientation"],
        "first_pivot_guard": guards["first_pivot_not_sufficient"],
        "decreasing_phase_guard": guards["decreasing_modulus_phase_guard"],
        "endpoint_pivot_guard": guards["endpoint_pivot_guard"],
        "cutoff_update": (
            "At N->N+1 write n=2^v m, z_*=exp(i omega_*log 2), and "
            "u_n=q_n/z_*^v. The exact polynomial update is "
            "F_(N+1)(z)-F_N(z)=j_0+u_n z^v, with "
            "u_(n,x)/u_n=-s_*'log n+delta_n"
            "-i*v*log(2)*omega_*'. At z=z_* this is "
            "q_n+j_0=Q_N/f_1. Hence only the constant and valuation-v "
            "coefficients change, but every downstream reflection "
            "coefficient may change. The certified adjacent bounds are "
            "real projections; no complex pivot-continuity bound is "
            "inferred from them."
        ),
        "route_decision": (
            "Retain Schur-Cohn as a sufficient endpoint-complete "
            "zero-free route and as an exact coordinate for falsification. "
            "Reject coefficient-modulus monotonicity, the first pivot "
            "alone, and prefix zero-freeness without the endpoint. The "
            "next arithmetic task is to attack Delta_K first using the "
            "actual recurrent endpoint and terminal dyadic layer. Only if "
            "that survives should the reduced coefficients be expanded "
            "into nested odd-prefix sums and tested for subsequent strict "
            "pivots."
        ),
        "weaker_fallback": (
            "Full unit-disk stability is stronger than final linked-point "
            "nonvanishing. If an actual Xi pivot fails, the recursion still "
            "identifies the exact failure stage and its inverse formula. A "
            "weaker route must then prove nonvanishing at z_* or a signed "
            "three-cylinder crossing budget while retaining the endpoint "
            "and complete first-jet tangents; it may not promote selected "
            "finite pivots to a disk theorem."
        ),
        "live_q_ge_1_target": (
            "On q=2tL^2>=1, first prove or refute uniformly "
            "|r_0+B_(0,0)|>|B_(K,0)| in every fixed-N chart, with the "
            "actual d_n correction and endpoint recurrence. If proved, "
            "derive the next normalized coefficient array exactly and "
            "seek a structural bound |alpha_d|<=1-epsilon_d whose product "
            "margin is strong enough to feed "
            "|mathsf_X|<=delta_L => |mathcal C_N|>A_L+epsilon_term and "
            "the successor trap 0<=kappa_j<1."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 layer remains a separate "
            "multiplicity-compatible parabolic/Hermite first-jet chart. "
            "Schur-Cohn disk stability is not imposed there, and no "
            "endpoint simplicity, chart join, or contact exclusion follows "
            "from this q>=1 algebra."
        ),
        "proof_boundary": (
            "The endpoint polynomial, outside-disk Schur step and inverse, "
            "recursive reflection criterion, conditional boundary margin, "
            "fixed-z coefficient current, centered companion, Frechet and "
            "mixed tangent recursions, cutoff update, and four algebraic "
            "audits are exact. No actual Xi pivot inequality, full "
            "unit-disk stability, endpoint-complete Xi lower bound, strict "
            "successor flux upper bound, q<1 closure, finite connector or "
            "chart join, contact exclusion, Lambda<=0, PF-infinity, RH "
            "proof, or Clay-prize conclusion is asserted."
        ),
    }


def build_rows(exact: dict) -> list[SchurRow]:
    return [
        SchurRow(
            "escfjg_00_pi_provenance",
            "definition_provenance",
            "certified",
            "The endpoint Schur coordinate introduces no unexplained pi.",
            exact["pi_provenance"],
            "The unit circle is a normalized complex domain, not a new source of pi.",
        ),
        SchurRow(
            "escfjg_01_endpoint_polynomial",
            "exact_reindexing",
            "ready_to_apply",
            "The endpoint-augmented physical prefix is an explicit dyadic-layer polynomial.",
            exact["endpoint_polynomial"],
            "Its physical coefficients are complex and endpoint-sensitive.",
        ),
        SchurRow(
            "escfjg_02_outside_step",
            "exact_identity",
            "ready_to_apply",
            "The constant-dominant Schur step lowers the degree exactly.",
            exact["reverse"] + " " + exact["outside_step"] + " " + exact["pivot"],
            "Exact zero top coefficients must be trimmed before the next stage.",
        ),
        SchurRow(
            "escfjg_03_stability_theorem",
            "exact_algebraic_theorem",
            "certified",
            "Strict outside Schur pivots are equivalent to closed-unit-disk zero exclusion.",
            exact["inverse"] + " " + exact["outside_stability_theorem"],
            "The theorem is conditional on verifying every physical Xi pivot.",
        ),
        SchurRow(
            "escfjg_04_reflection_recursion",
            "exact_identity",
            "ready_to_apply",
            "The normalized recursion exposes one reflection coefficient per effective degree.",
            exact["normalized_reflection_recursion"] + " " + exact["reflection_criterion"],
            "No actual Xi reflection coefficient is yet bounded.",
        ),
        SchurRow(
            "escfjg_05_boundary_margin",
            "exact_conditional_corollary",
            "ready_to_apply",
            "Strict reflection bounds give an explicit unit-circle modulus product.",
            exact["boundary_margin"],
            "The product may be small and its Xi factors remain unproved.",
        ),
        SchurRow(
            "escfjg_06_coefficient_x_tangent",
            "exact_differential_reduction",
            "ready_to_apply",
            "The fixed-z coefficient derivative rejoins the physical H_1,D_0 current.",
            exact["coefficient_x_tangent"],
            "The moving linked point must be restored after fixed-z differentiation.",
        ),
        SchurRow(
            "escfjg_07_centered_companion",
            "exact_differential_reduction",
            "ready_to_apply",
            "A companion polynomial carries Z_A and its complete first x derivative.",
            exact["centered_companion"],
            "H_2,D_0,D_1 and endpoint derivatives are retained.",
        ),
        SchurRow(
            "escfjg_08_frechet_recursion",
            "exact_identity",
            "ready_to_apply",
            "The Schur map transports arbitrary coefficient tangents before absolute values.",
            exact["frechet_recursion"],
            "This is a tangent identity, not a tangent lower bound.",
        ),
        SchurRow(
            "escfjg_09_mixed_recursion",
            "exact_identity",
            "ready_to_apply",
            "The mixed x/centered tangent has an exact Hessian correction.",
            exact["frechet_bilinear"] + " " + exact["mixed_recursion"] + " " + exact["tangent_application"],
            "Dropping the bilinear term would corrupt Z_(A,x).",
        ),
        SchurRow(
            "escfjg_10_first_xi_pivot",
            "new_inequality_target",
            "not_ready_to_apply",
            "The first missing Xi inequality is endpoint versus terminal-layer dominance.",
            exact["first_xi_pivot"],
            "Dyadic contraction at the synchronized base does not prove it.",
        ),
        SchurRow(
            "escfjg_11_degree_one_orientation",
            "exact_convention_audit",
            "guard_validated",
            "Degree one fixes the outside-zero orientation.",
            json.dumps(exact["degree_one_orientation"], sort_keys=True),
            "This audit fixes signs and conjugations only.",
            exact["degree_one_orientation"],
        ),
        SchurRow(
            "escfjg_12_first_pivot_guard",
            "countermodel",
            "guard_validated",
            "One strict pivot does not imply disk zero-freeness.",
            json.dumps(exact["first_pivot_guard"], sort_keys=True),
            "Every reduced pivot remains an independent obligation.",
            exact["first_pivot_guard"],
        ),
        SchurRow(
            "escfjg_13_decreasing_phase_guard",
            "countermodel",
            "guard_validated",
            "Decreasing coefficient moduli do not survive arbitrary phases as a zero-free theorem.",
            json.dumps(exact["decreasing_phase_guard"], sort_keys=True),
            "Xi-specific phase structure could help only through a proved inequality.",
            exact["decreasing_phase_guard"],
        ),
        SchurRow(
            "escfjg_14_endpoint_pivot_guard",
            "countermodel",
            "guard_validated",
            "Endpoint addition can collapse the first Schur pivot.",
            json.dumps(exact["endpoint_pivot_guard"], sort_keys=True),
            "The countermodel is generic, not an Xi endpoint evaluation.",
            exact["endpoint_pivot_guard"],
        ),
        SchurRow(
            "escfjg_15_cutoff_update",
            "exact_chart_transition",
            "ready_to_apply",
            "A cutoff change updates two polynomial coefficients and then all recursive pivots.",
            exact["cutoff_update"],
            "Only certified real projected adjacent bounds are retained.",
        ),
        SchurRow(
            "escfjg_16_route_decision",
            "route_decision",
            "guard_validated",
            "Test the first physical pivot before expanding the full recursion.",
            exact["route_decision"],
            "Selected finite pivots cannot prove a uniform Xi theorem.",
        ),
        SchurRow(
            "escfjg_17_weaker_fallback",
            "conditional_route",
            "not_ready_to_apply",
            "A failed full-disk route must fall back to linked nonvanishing or signed degree.",
            exact["weaker_fallback"],
            "The endpoint and complete first jet may not be removed in the fallback.",
        ),
        SchurRow(
            "escfjg_18_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 route begins with a uniform endpoint-terminal pivot theorem.",
            exact["live_q_ge_1_target"],
            "The first pivot and all successor bounds remain open.",
        ),
        SchurRow(
            "escfjg_19_q_lt_1_nonpromotion",
            "open_theorem_target",
            "not_ready_to_apply",
            "The small-q multiplicity chart remains separate.",
            exact["q_lt_1_target"],
            "No RH or Newman conclusion follows from the q>=1 Schur algebra.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    rows = build_rows(exact)
    return {
        "kind": STEM,
        "date": "2026-07-27",
        "status": (
            "exact endpoint outside-disk Schur-Cohn recursion, conditional "
            "boundary margin, complete first-jet tangent propagation, "
            "cutoff update, and three route guards; all physical Xi pivot "
            "inequalities remain open"
        ),
        "proof_boundary": exact["proof_boundary"],
        "sources": {
            key: str(path.relative_to(REPO_ROOT))
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman Endpoint Schur-Cohn First-Jet Guard",
            "",
            "Date: 2026-07-27",
            "",
            "Status: exact endpoint stability coordinate and first-jet",
            "recursion with countermodels; not a proof of `Lambda<=0`,",
            "PF-infinity, RH, or a Clay-prize result.",
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
            "## Endpoint Polynomial",
            "",
            "```text",
            exact["endpoint_polynomial"],
            "```",
            "",
            "## Outside-Disk Schur Step",
            "",
            "```text",
            exact["reverse"],
            exact["outside_step"],
            exact["pivot"],
            exact["inverse"],
            "```",
            "",
            exact["outside_stability_theorem"],
            "",
            "## Reflection Recursion",
            "",
            "```text",
            exact["normalized_reflection_recursion"],
            exact["reflection_criterion"],
            "```",
            "",
            "## Conditional Boundary Margin",
            "",
            "```text",
            exact["boundary_margin"],
            "```",
            "",
            "## Physical Coefficient Current",
            "",
            "```text",
            exact["coefficient_x_tangent"],
            "```",
            "",
            "## Centered Companion",
            "",
            "```text",
            exact["centered_companion"],
            "```",
            "",
            "## Tangent Recursion",
            "",
            "```text",
            exact["frechet_recursion"],
            exact["frechet_bilinear"],
            exact["mixed_recursion"],
            exact["tangent_application"],
            "```",
            "",
            "## First Physical Pivot",
            "",
            "```text",
            exact["first_xi_pivot"],
            "```",
            "",
            "## Exact Guards",
            "",
            "### Degree-one orientation",
            "",
            "```json",
            json.dumps(exact["degree_one_orientation"], indent=2, sort_keys=True),
            "```",
            "",
            "### First pivot is not sufficient",
            "",
            "```json",
            json.dumps(exact["first_pivot_guard"], indent=2, sort_keys=True),
            "```",
            "",
            "### Decreasing moduli with physical-style phases",
            "",
            "```json",
            json.dumps(exact["decreasing_phase_guard"], indent=2, sort_keys=True),
            "```",
            "",
            "### Endpoint pivot collapse",
            "",
            "```json",
            json.dumps(exact["endpoint_pivot_guard"], indent=2, sort_keys=True),
            "```",
            "",
            "## Cutoff Update",
            "",
            "```text",
            exact["cutoff_update"],
            "```",
            "",
            "## Route Decision",
            "",
            "```text",
            exact["route_decision"],
            exact["weaker_fallback"],
            "```",
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
        "built Newman endpoint Schur-Cohn first-jet guard: "
        "20 rows, 1 exact outside-stability theorem, "
        "1 conditional reflection-product margin, "
        "3 exact route guards, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
