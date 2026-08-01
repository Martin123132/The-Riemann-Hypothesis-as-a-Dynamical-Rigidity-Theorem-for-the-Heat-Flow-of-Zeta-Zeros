#!/usr/bin/env python3
"""Build the physical q=1 saddle-phase and phase-variance reduction."""

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
    "contiguous_terminal_tail_anchored_five_moment_q1_saddle_"
    "phase_variance_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "critical_frame": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_residual_reduction.json"
    ),
    "reciprocal_saddle": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "reciprocal_saddle_self_duality_gate.json"
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
    "endpoint_centering": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "mangoldt_contact_centering_gate.json"
    ),
    "turan": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "contiguous_terminal_tail_anchored_five_moment_logarithmic_"
        "phase_turan_reduction.json"
    ),
}


@dataclass(frozen=True)
class SaddleVarianceRow:
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
    if sp.simplify(sp.expand(expression)) != 0:
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
    frame = payloads["critical_frame"].get("exact", {}).get(
        "critical_frame", {}
    )
    for key, marker in (
        ("definitions", "a^2=x/(4*pi)+t/16"),
        ("alpha_real", "Re(alpha)=L/2"),
        ("alpha_imag", "Im(alpha)=-(1/2)atan(x)"),
    ):
        if marker not in frame.get(key, ""):
            raise RuntimeError(f"critical-frame source drifted: {key}")

    reciprocal = payloads["reciprocal_saddle"].get("exact", {})
    if "T_0-omega" not in reciprocal.get("phase_saddle_comparison", ""):
        raise RuntimeError("reciprocal saddle source drifted")

    anchor = payloads["phase_anchor"].get("exact", {})
    if "eta=f_1/|f_1|" not in anchor.get("unit_anchor", ""):
        raise RuntimeError("absolute phase anchor drifted")

    projection = payloads["direct_projection"].get("exact", {})
    for key, marker in (
        ("saddle_amplitude", "B_a=1/2-t*delta_a/2"),
        ("saddle_bounds", "B_a=1/2-t*delta_a/2>49/100"),
        ("endpoint_projection", "B_0=beta/(|M_t(s)|*|f_1|)>0"),
    ):
        if marker not in projection.get(key, ""):
            raise RuntimeError(f"projection source drifted: {key}")

    endpoint = payloads["endpoint_centering"].get("exact", {})
    if "kappa=(-1)^N B_0 real" not in endpoint.get("endpoint_factor", ""):
        raise RuntimeError("endpoint parity source drifted")

    turan = payloads["turan"]
    leading = turan.get("leading_current_certificate", {}).get(
        "leading_turan", ""
    )
    if "P_bulk^(0)" not in leading or "u_x" not in leading:
        raise RuntimeError("leading Turan source drifted")
    if turan.get("counts", {}).get("signed_joint_bounds") != 0:
        raise RuntimeError("upstream Turan proof boundary drifted")

    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "source_kinds": {
            key: payload.get("kind", "") for key, payload in payloads.items()
        },
        "imported_signed_physical_bounds": 0,
        "imported_phi_b_bounds": 0,
    }


def physical_q1_certificate() -> dict:
    L, x = sp.symbols("L x", positive=True, real=True)
    t = 1 / (2 * L**2)
    y = x**-2
    alpha_real = L / 2 + sp.log(1 + y) / 4 - 1 / (1 + x**2)
    alpha_imag = -sp.atan(x) / 2 + 3 * x / (1 + x**2)

    sigma = sp.Rational(1, 2) + t * alpha_real / 2
    sigma_expected = (
        sp.Rational(1, 2)
        + 1 / (8 * L)
        + sp.log(1 + y) / (16 * L**2)
        - 1 / (4 * L**2 * (1 + x**2))
    )
    require_zero(sigma - sigma_expected, "q=1 sigma identity failed")

    omega = x / 2 - t * alpha_imag / 2
    omega_atan_x = (
        x / 2
        + sp.atan(x) / (8 * L**2)
        - 3 * x / (4 * L**2 * (1 + x**2))
    )
    require_zero(omega - omega_atan_x, "q=1 frequency identity failed")

    t0 = x / 2 + sp.pi / (16 * L**2)
    epsilon = (
        sp.atan(1 / x) / (8 * L**2)
        + 3 * x / (4 * L**2 * (1 + x**2))
    )
    omega_reciprocal = omega_atan_x.xreplace(
        {sp.atan(x): sp.pi / 2 - sp.atan(1 / x)}
    )
    require_zero(
        omega_reciprocal - (t0 - epsilon),
        "q=1 saddle-frequency defect failed",
    )

    z = sp.pi / (8 * L**2 * x)
    delta_0 = 1 / (1 + x**2) - sp.log(1 + y) / 4
    delta_a = delta_0 + sp.log(1 + z) / 2
    b_a = sp.Rational(1, 2) - delta_a / (4 * L**2)
    log_a = L / 2 + sp.log(1 + z) / 2
    require_zero(
        delta_a - (log_a - alpha_real),
        "q=1 saddle-amplitude defect failed",
    )

    u = sp.symbols("u", nonnegative=True, real=True)
    profile_exponent = b_a * u + u**2 / (8 * L**2)
    profile_rate = sp.diff(profile_exponent, u)
    require_zero(
        profile_rate - (b_a + u / (4 * L**2)),
        "q=1 amplitude profile derivative failed",
    )

    return {
        "domain": (
            "On q=2tL^2=1 with L>=50, t=1/(2L^2), "
            "x=4*pi*exp(L), a^2=exp(L)+1/(32L^2), and "
            "T_0=2*pi*a^2=x/2+pi/(16L^2)."
        ),
        "sigma": (
            "sigma=Re(s_*)=1/2+1/(8L)+log(1+x^(-2))/(16L^2)"
            "-1/[4L^2(1+x^2)]."
        ),
        "sigma_bounds": (
            "1/2+1/(8L)-1/(4L^2x^2)<sigma<"
            "1/2+1/(8L)-1/(16L^2x^2)."
        ),
        "frequency": (
            "Omega=-Im(s_*)=x/2+atan(x)/(8L^2)"
            "-3x/[4L^2(1+x^2)]=T_0-epsilon."
        ),
        "epsilon": (
            "epsilon=atan(1/x)/(8L^2)+3x/[4L^2(1+x^2)]."
        ),
        "epsilon_bounds": (
            "3/(8L^2x)<epsilon<7/(8L^2x)."
        ),
        "amplitude_defect": (
            "delta_a=log(a)-Re(alpha)=delta_0+(1/2)log(1+z), "
            "delta_0=1/(1+x^2)-(1/4)log(1+x^(-2)), "
            "z=pi/(8L^2x)."
        ),
        "amplitude_defect_bounds": (
            "0<delta_a<1/x^2+pi/(16L^2x)."
        ),
        "b_a": "B_a=1/2-delta_a/(4L^2).",
        "b_a_bounds": (
            "1/2-1/(4L^2x^2)-pi/(64L^4x)<B_a<1/2; "
            "in particular B_a>49/100."
        ),
        "profile": (
            "For u_n=log(a/n), A_n=A_a*exp[B_a*u_n+u_n^2/(8L^2)], "
            "A_a=|eta/S_a|*exp[t*log(a)^2/4-sigma*log(a)]>0."
        ),
        "profile_rate": (
            "partial_u log A(u)=B_a+u/(4L^2)>49/100 for u>=0. "
            "Thus n<m<a implies A_n>A_m."
        ),
        "symbolic_differences": {
            "sigma": str(sp.simplify(sigma - sigma_expected)),
            "frequency": str(sp.simplify(omega - omega_atan_x)),
            "saddle_frequency": str(
                sp.simplify(omega_reciprocal - (t0 - epsilon))
            ),
            "amplitude_defect": str(
                sp.simplify(delta_a - (log_a - alpha_real))
            ),
            "profile_rate": str(
                sp.simplify(profile_rate - (b_a + u / (4 * L**2)))
            ),
        },
        "elementary_bound_proof": (
            "Put y=x^(-2) in (0,1). The inequalities log(1+y)<y, "
            "y/(1+y)>y/2, atan(1/x)<1/x, and "
            "1/(2x)<x/(1+x^2)<1/x give all displayed bounds."
        ),
    }


def phase_anchor_certificate() -> dict:
    return {
        "real_scale": (
            "The retained endpoint scale is S_a=kappa*T_0 with "
            "kappa=(-1)^N*B_0, B_0>0, and T_0>0. Hence "
            "S_a is real and sign(S_a)=(-1)^N."
        ),
        "common_phase": (
            "eta/S_a=[(-1)^N*eta]/(B_0*T_0), so the physical common "
            "unit phase is omega_c=(-1)^N*eta="
            "(-1)^N*phi*(1+d_1)/|1+d_1|."
        ),
        "distance_phase": (
            "Since tau=Im(s_*)=-Omega, w_n=omega_c*A_n*exp(i*Omega*log n). "
            "With omega_a=omega_c*exp(i*Omega*log a), this is exactly "
            "w_n=omega_a*A_n*exp(-i*Omega*u_n)."
        ),
        "phase_boundary": (
            "The common phase is fixed by the Xi source; it is not a free "
            "rotation in a physical sign argument."
        ),
    }


def saddle_transfer_certificate() -> dict:
    return {
        "moments": (
            "Define H_k(xi)=sum_(n<=B)u_n^k*A_n*exp(-i*xi*u_n) and "
            "M_k=omega_a*H_k(Omega), 0<=k<=4."
        ),
        "transfer": (
            "For 0<=k<=4, |H_k(Omega)-H_k(T_0)|<="
            "epsilon*sum_(n<=B)u_n^(k+1)A_n."
        ),
        "proof": (
            "Termwise use |exp(-i*Omega*u)-exp(-i*T_0*u)|<="
            "|Omega-T_0|u=epsilon*u. No cancellation or denominator is used."
        ),
        "interpretation": (
            "The actual phase is exponentially close to the exact cutoff "
            "saddle phase in every unnormalized moment. Relative moments "
            "still require a lower bound on H_0 and are not inferred."
        ),
    }


def two_carrier_certificate() -> dict:
    a1, a2, u1, u2 = sp.symbols(
        "A_1 A_2 u_1 u_2", positive=True, real=True
    )
    c1, c2, s1, s2 = sp.symbols("c_1 c_2 s_1 s_2", real=True)
    f = a1 * c1 + a2 * c2
    fp = a1 * u1 * s1 + a2 * u2 * s2
    fpp = -(a1 * u1**2 * c1 + a2 * u2**2 * c2)
    four_t = sp.expand(f * fpp - fp**2)
    expanded = (
        -a1**2 * u1**2 * (c1**2 + s1**2)
        - a2**2 * u2**2 * (c2**2 + s2**2)
        - a1
        * a2
        * ((u1**2 + u2**2) * c1 * c2 + 2 * u1 * u2 * s1 * s2)
    )
    require_zero(four_t - expanded, "two-carrier Turan expansion failed")
    upper = (
        -a1**2 * u1**2
        - a2**2 * u2**2
        + a1 * a2 * (u1**2 + u2**2)
    )
    factored_upper = (a2 - a1) * (a1 * u1**2 - a2 * u2**2)
    require_zero(upper - factored_upper, "two-carrier upper factor failed")
    return {
        "expansion": (
            "For w_j=A_j(c_j+i s_j), c_j^2+s_j^2=1, the pure Turan "
            "numerator T_2=[f*f''-(f')^2]/4 satisfies "
            "4T_2=-A_1^2u_1^2-A_2^2u_2^2-A_1A_2[(u_1^2+u_2^2)c_1c_2"
            "+2u_1u_2s_1s_2]."
        ),
        "ordered_bound": (
            "If A_1>A_2>0 and u_1>u_2>0, then "
            "4T_2<=(A_2-A_1)(A_1u_1^2-A_2u_2^2)<0."
        ),
        "proof": (
            "The cross bracket is the bilinear form diag(u_1^2+u_2^2,"
            "2u_1u_2) on two unit vectors, so it is at least "
            "-(u_1^2+u_2^2)."
        ),
        "old_guard_exclusion": (
            "The Section 11.161 witness has u_1=2log(2)>u_2=log(2) but "
            "A_1=1/2<A_2=1, opposite to the physical q=1 ordering. It is "
            "therefore excluded by the exact physical amplitude law."
        ),
        "boundary": (
            "This signs the two-carrier pure Turan numerator. The full "
            "P_bulk^(0) also contains u_(N,x)f*g/2; the old real guard has "
            "g=0, but no arbitrary-phase full-current theorem is claimed."
        ),
        "symbolic_differences": {
            "expansion": str(sp.simplify(four_t - expanded)),
            "upper_factor": str(sp.simplify(upper - factored_upper)),
        },
    }


def phase_variance_certificate() -> dict:
    f, g, p, q, r, s, u_x = sp.symbols(
        "f g p q r s u_x", real=True
    )
    w = f + sp.I * g
    m1 = p + sp.I * q
    m2 = r + sp.I * s
    radius2 = sp.expand(f**2 + g**2)
    a_num = sp.expand_complex(sp.re(m1 * sp.conjugate(w)))
    c_num = sp.expand((m2 * w - m1**2) * sp.conjugate(w) ** 2)
    c_real = sp.expand_complex(sp.re(c_num))
    c_imag = sp.expand_complex(sp.im(c_num))
    discriminant = sp.expand(
        f**2 * c_real
        + a_num**2 * radius2
        - (2 * u_x * radius2**2 + c_imag) * f * g
    )
    four_current = sp.expand(-f * r - q**2 + 2 * u_x * f * g)
    require_zero(
        discriminant + four_current * radius2**2,
        "branch-free phase-variance discriminant failed",
    )
    return {
        "relative_moments": (
            "On H_0(Omega)!=0 put mu_k=H_k(Omega)/H_0(Omega) and "
            "v=mu_2-mu_1^2. For W(y)=omega_a*H_0(Omega+y)="
            "R(y)exp(i theta(y)), W'/W=-i*mu_1, "
            "(log R)''=-Re(v), theta'=-Re(mu_1), and theta''=-Im(v)."
        ),
        "phase_identity": (
            "Where f=R cos(theta)!=0, "
            "4P_bulk^(0)/f^2=-Re(v)-[Re(mu_1)]^2 sec(theta)^2"
            "+[2u_(N,x)+Im(v)]tan(theta)."
        ),
        "phase_discriminant": (
            "Define D_phase=Re(v)+[Re(mu_1)]^2 sec(theta)^2"
            "-[2u_(N,x)+Im(v)]tan(theta). Then "
            "P_bulk^(0)=-f^2 D_phase/4."
        ),
        "branch_free": (
            "For W=M_0=f+ig, M_1, M_2, R2=|W|^2, "
            "A=Re(M_1*conj(W)), and C=(M_2W-M_1^2)conj(W)^2, put "
            "D_hat=f^2Re(C)+A^2R2-[2u_(N,x)R2^2+Im(C)]fg. "
            "On W!=0, 4P_bulk^(0)R2^2=-D_hat, so P_bulk^(0)<=0 iff "
            "D_hat>=0."
        ),
        "exceptional_fibre": (
            "At W=0 the relative moments are undefined, but the primary "
            "identity gives P_bulk^(0)=-(Im M_1)^2/4<=0."
        ),
        "symbolic_difference": str(
            sp.simplify(discriminant + four_current * radius2**2)
        ),
    }


def ordered_profile_guard() -> dict:
    d = sp.log(2)
    limiting_amplitudes = (sp.Integer(4), 2 * sp.sqrt(2), sp.sqrt(2))
    distances = (4 * d, 3 * d, d)
    weights = (
        -limiting_amplitudes[0],
        limiting_amplitudes[1],
        limiting_amplitudes[2],
    )
    f = sp.simplify(sum(weights))
    s2 = sp.simplify(
        sum(distance**2 * weight for distance, weight in zip(distances, weights))
    )
    current = sp.simplify(-f * s2 / 4)
    expected = (-185 + 134 * sp.sqrt(2)) * d**2 / 2
    require_zero(current - expected, "ordered profile limiting guard failed")

    exp_argument = sp.Rational(3381, 10000)
    exp_taylor = (
        1
        + exp_argument
        + exp_argument**2 / 2
        + exp_argument**3 / 6
    )
    exp_margin = sp.factor(exp_taylor - sp.Rational(7, 5))
    f_ratio_lower = sp.factor(
        sp.Rational(7071, 10000)
        * (1 - sp.Rational(343, 2000000))
        + sp.Rational(3535, 10000)
        * (1 - sp.Rational(735, 2000000))
        - 1
    )
    curvature_ratio_lower = sp.factor(
        16 - 9 * sp.Rational(5, 7) - sp.Rational(125, 343)
    )
    if exp_margin <= 0:
        raise RuntimeError("ordered profile r lower bound failed")
    if f_ratio_lower <= sp.Rational(3, 50):
        raise RuntimeError("ordered profile value margin failed")
    if curvature_ratio_lower <= 9:
        raise RuntimeError("ordered profile curvature margin failed")

    return {
        "support": (
            "Take a=16k and n=k,2k,8k, so u=(4d,3d,d), d=log(2). "
            "Use the exact q=1 correction-free profile "
            "A(u)=A_a exp(B_a u+c u^2), 49/100<B_a<1/2, "
            "0<c=1/(8L^2)<=1/20000."
        ),
        "phase": (
            "The artificial common logarithmic half-turn tau=pi/d with "
            "omega=-exp(i*tau*log k) gives weights (-A_1,+A_2,+A_3)."
        ),
        "finite_margin": (
            "Writing r=exp(B_a d), gamma=exp(c d^2), the normalized "
            "amplitudes are (r^4 gamma^16,r^3 gamma^9,r gamma). "
            "Uniformly on the q=1 box, f/A_1>3/50 and "
            "f''/(A_1 d^2)>9, while f'=g=0. Hence "
            "P_bulk^(0)>27 A_1^2 d^2/200>0."
        ),
        "limiting_guard": (
            "At B_a=1/2,c=0 the amplitudes are (4,2sqrt(2),sqrt(2)) "
            "and P_bulk^(0)=[(-185+134sqrt(2))/2]log(2)^2>0."
        ),
        "rational_certificate": {
            "exp_taylor_minus_7_over_5": str(exp_margin),
            "f_over_a1_lower": str(f_ratio_lower),
            "f_over_a1_minus_3_over_50": str(
                sp.factor(f_ratio_lower - sp.Rational(3, 50))
            ),
            "fpp_over_a1d2_lower": str(curvature_ratio_lower),
            "limiting_current": str(current),
        },
        "scope": (
            "This guard obeys the exact finite-q=1 amplitude profile, strict "
            "amplitude ordering, dyadic integer support, and one common "
            "logarithmic phase. Its tau and common phase are deliberately "
            "not the physical Omega and omega_c. It proves that amplitude "
            "ordering and profile shape alone cannot close the current."
        ),
        "pi_provenance": (
            "The guard's pi is only the half-turn condition tau*log(2)=pi. "
            "The physical pi in T_0=2*pi*a^2 comes independently from the "
            "completed-zeta/Riemann-Siegel and Poisson normalization."
        ),
    }


def exact_payload() -> dict:
    return {
        "full_current": (
            "Retain the exact upstream decomposition "
            "h^2*Phi_B=E+P_bulk^(0)+R_corr. On W!=0 this is "
            "h^2*Phi_B=E-D_hat/(4|W|^4)+R_corr. On W=0 it is "
            "h^2*Phi_B=E-(Im M_1)^2/4+R_corr."
        ),
        "ordinary_target": (
            "On W!=0 prove the physical source inequality "
            "E+R_corr-D_hat/(4|W|^4)<=h^2/400, preferably h^2/800, "
            "using Omega=T_0-epsilon, the exact A(u), and omega_a."
        ),
        "exceptional_target": (
            "On W=0 prove E+R_corr<=h^2/400+(Im M_1)^2/4. "
            "No relative moment or phase division is allowed there."
        ),
        "next_analytic_stage": (
            "First seek a division-free estimate for D_hat and the composed "
            "E+R_corr at the exact saddle phase T_0, then transfer the five "
            "unnormalized moments to Omega with the epsilon bounds. The "
            "endpoint and transpose phase-sum family must remain composed."
        ),
        "proof_boundary": (
            "This proves exact physical q=1 parameter and common-phase laws, "
            "strict correction-free amplitude ordering, five saddle-phase "
            "transfer bounds, a two-carrier pure-Turan ordering theorem, the "
            "ordinary-fibre phase-variance identity, its branch-free "
            "discriminant, and one ordered-profile nonpromotion guard. It "
            "does not prove D_hat>=0 for the Xi source, bound E or R_corr, "
            "upper-bound Phi_B, establish a signed joint Type-I/II or Vaughan "
            "estimate, exclude contact, transfer the retained model to Xi, "
            "prove Q209, Lambda<=0, PF-infinity, RH, or a prize-level result."
        ),
    }


def build_rows(
    physical: dict,
    anchor: dict,
    transfer: dict,
    two: dict,
    variance: dict,
    guard: dict,
    exact: dict,
) -> list[SaddleVarianceRow]:
    return [
        SaddleVarianceRow("spv_01_domain", "physical q=1 domain", "proved", "The source parameters have one exact q=1 specialization.", physical["domain"], "The fixed physical terminal/bulk chart is retained."),
        SaddleVarianceRow("spv_02_sigma", "real saddle exponent", "proved", "The q=1 real exponent and exponentially small correction are explicit.", physical["sigma"] + " " + physical["sigma_bounds"], physical["elementary_bound_proof"], physical["symbolic_differences"]),
        SaddleVarianceRow("spv_03_frequency", "physical logarithmic frequency", "proved", "The positive physical logarithmic frequency is explicit.", physical["frequency"] + " " + physical["epsilon"], "Omega denotes -Im(s_*), while the previous tau is Im(s_*)=-Omega."),
        SaddleVarianceRow("spv_04_lock", "saddle-frequency lock", "proved", "The physical phase is exponentially close to the cutoff saddle phase.", physical["epsilon_bounds"], "Floors and endpoint recurrence are not identified from this smallness."),
        SaddleVarianceRow("spv_05_profile", "physical amplitude profile", "proved", "The correction-free amplitudes have one exact saddle-distance profile.", physical["amplitude_defect"] + " " + physical["b_a"] + " " + physical["profile"], "The d_n correction factors belong to R_corr and are not inserted into A_n."),
        SaddleVarianceRow("spv_06_monotone", "strict amplitude ordering", "proved", "Every correction-free bulk amplitude decreases strictly with n.", physical["amplitude_defect_bounds"] + " " + physical["b_a_bounds"] + " " + physical["profile_rate"], "This is an amplitude theorem, not a phase or full-current sign theorem."),
        SaddleVarianceRow("spv_07_anchor", "parity-fixed common phase", "proved", "The physical common phase is fixed by the branch-free Xi anchor and endpoint parity.", anchor["real_scale"] + " " + anchor["common_phase"], anchor["phase_boundary"]),
        SaddleVarianceRow("spv_08_distance", "physical distance Fourier form", "proved", "The actual carriers form one positive-amplitude distance-frequency polynomial.", anchor["distance_phase"] + " " + transfer["moments"], "The common phase cancels from relative moments but remains in the real projection."),
        SaddleVarianceRow("spv_09_transfer", "saddle-phase moment transfer", "proved", "All five unnormalized moments transfer from Omega to T_0 with an explicit positive-mass error.", transfer["transfer"] + " " + transfer["proof"], transfer["interpretation"]),
        SaddleVarianceRow("spv_10_old_guard", "old-guard exclusion", "proved", "The Section 11.161 two-atom witness violates physical amplitude ordering.", two["old_guard_exclusion"], "This removes that witness, not every possible multi-carrier obstruction."),
        SaddleVarianceRow("spv_11_two", "ordered two-carrier Turan sign", "proved", "Similarly ordered amplitudes and distances force the two-carrier pure Turan numerator negative for arbitrary phases.", two["expansion"] + " " + two["ordered_bound"] + " " + two["proof"], two["boundary"], two["symbolic_differences"]),
        SaddleVarianceRow("spv_12_variance", "complex phase variance", "proved", "The ordinary-fibre logarithmic derivatives are controlled by one complex variance.", variance["relative_moments"] + " " + variance["phase_identity"], "This coordinate requires H_0(Omega)!=0 and f!=0."),
        SaddleVarianceRow("spv_13_discriminant", "branch-free ordinary discriminant", "proved", "The phase-variance sign criterion has a branch-free moment form.", variance["phase_discriminant"] + " " + variance["branch_free"] + " " + variance["exceptional_fibre"], "The original division-free current remains primary at W=0.", {"difference": variance["symbolic_difference"]}),
        SaddleVarianceRow("spv_14_guard", "ordered-profile phase guard", "guard_validated", "Physical q=1 amplitude shape and ordering still do not sign the current without the actual phase lock.", guard["support"] + " " + guard["phase"] + " " + guard["finite_margin"], guard["scope"], guard["rational_certificate"]),
        SaddleVarianceRow("spv_15_full", "exact full-current split", "proved", "The phase discriminant rejoins the endpoint and correction kernels without loss.", exact["full_current"] + " " + exact["ordinary_target"] + " " + exact["exceptional_target"], "Neither E nor R_corr is bounded here."),
        SaddleVarianceRow("spv_16_target", "physical q=1 theorem target", "open", "The next theorem is an endpoint-composed saddle-phase estimate for the actual Xi source.", exact["next_analytic_stage"], exact["proof_boundary"]),
    ]


def render_note(payload: dict) -> str:
    physical = payload["physical_q1_certificate"]
    anchor = payload["phase_anchor_certificate"]
    transfer = payload["saddle_transfer_certificate"]
    two = payload["two_carrier_certificate"]
    variance = payload["phase_variance_certificate"]
    guard = payload["ordered_profile_guard"]
    exact = payload["exact"]
    return f"""# Physical q=1 Saddle Phase And Variance Reduction

Date: 2026-08-01

Status: exact source reduction with `0 signed physical bounds` and `0 Phi_B
bounds`. This is not a proof of RH; the endpoint-composed phase estimate is
open.

## Physical Parameters

```text
{physical['domain']}

{physical['sigma']}
{physical['sigma_bounds']}

{physical['frequency']}
{physical['epsilon']}
{physical['epsilon_bounds']}
```

Here `Omega=-Im(s_*)>0`; the `tau` of Section 11.161 is `-Omega`.

## Amplitude And Anchor

```text
{physical['amplitude_defect']}
{physical['amplitude_defect_bounds']}
{physical['b_a']}
{physical['b_a_bounds']}
{physical['profile']}
{physical['profile_rate']}

{anchor['real_scale']}
{anchor['common_phase']}
{anchor['distance_phase']}
```

Thus the physical common phase is fixed, while the relative moments depend
only on the positive amplitude profile and `Omega`.

## Saddle Transfer

```text
{transfer['moments']}
{transfer['transfer']}
```

{transfer['proof']} {transfer['interpretation']}

## Two Carriers

```text
{two['expansion']}
{two['ordered_bound']}
```

{two['proof']} {two['old_guard_exclusion']}

{two['boundary']}

## Phase Variance

```text
{variance['relative_moments']}

{variance['phase_identity']}

{variance['phase_discriminant']}
```

The branch-free ordinary-fibre form is

```text
{variance['branch_free']}
```

{variance['exceptional_fibre']}

## Stronger Nonpromotion Guard

```text
{guard['support']}
{guard['phase']}
{guard['finite_margin']}
{guard['limiting_guard']}
```

{guard['scope']}

### Pi Provenance

{guard['pi_provenance']}

## Exact Remaining Target

```text
{exact['full_current']}

{exact['ordinary_target']}

{exact['exceptional_target']}
```

{exact['next_analytic_stage']}

## Boundary

{exact['proof_boundary']}

## Reproduce

```text
python work/rh_compute/scripts/{STEM}.py
python work/rh_compute/scripts/check_{STEM}.py
```
"""


def build_payload() -> dict:
    payloads = load_sources()
    source = source_audit(payloads)
    physical = physical_q1_certificate()
    anchor = phase_anchor_certificate()
    transfer = saddle_transfer_certificate()
    two = two_carrier_certificate()
    variance = phase_variance_certificate()
    guard = ordered_profile_guard()
    exact = exact_payload()
    rows = build_rows(physical, anchor, transfer, two, variance, guard, exact)
    return {
        "kind": STEM,
        "date": "2026-08-01",
        "status": (
            "exact physical q=1 saddle parameter, common-phase, amplitude "
            "ordering, five-moment saddle transfer, two-carrier Turan, "
            "phase-variance, and ordered-profile guard reduction; physical "
            "joint bound open; no Phi_B bound, Xi-level theorem, or RH"
        ),
        "proof_boundary": exact["proof_boundary"],
        "source_audit": source,
        "physical_q1_certificate": physical,
        "phase_anchor_certificate": anchor,
        "saddle_transfer_certificate": transfer,
        "two_carrier_certificate": two,
        "phase_variance_certificate": variance,
        "ordered_profile_guard": guard,
        "exact": exact,
        "counts": {
            "rows": len(rows),
            "physical_q1_parameter_laws": 4,
            "common_phase_laws": 1,
            "amplitude_monotonicity_theorems": 1,
            "saddle_phase_transfer_bounds": 5,
            "two_carrier_order_signs": 1,
            "phase_variance_identities": 1,
            "ordinary_fibre_discriminants": 1,
            "ordered_profile_guards": 1,
            "signed_physical_targets": 1,
            "signed_physical_bounds": 0,
            "phi_b_bounds": 0,
            "xi_level_current_theorems": 0,
        },
        "rows": [asdict(row) for row in rows],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    payload = build_payload()
    atomic_write(args.out, json.dumps(payload, indent=2) + "\n")
    atomic_write(args.note, render_note(payload))
    print(
        "wrote physical q=1 saddle-phase variance reduction: "
        f"{payload['counts']['rows']} rows, 5 saddle transfers, "
        "1 ordered-profile guard, 0 signed physical bounds"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
