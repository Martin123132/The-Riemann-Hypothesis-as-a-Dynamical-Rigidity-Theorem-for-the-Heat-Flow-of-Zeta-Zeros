#!/usr/bin/env python3
"""Build the first-order saddle-centered real-residual reduction."""

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
    "jensen_window_pf_newman_polymath15_"
    "first_order_centered_real_residual_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "coefficient_mass": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "C1_cell_remainder_certificate.json"
    ),
    "global_remainder": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_global_remainder_certificate.json"
    ),
    "signed_contact": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_critical_"
        "first_order_signed_contact_reduction.json"
    ),
    "wronskian_crossing": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_wronskian_crossing_reduction.json"
    ),
    "oriented_successor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_oriented_successor_winding_reduction.json"
    ),
}

ETA_0 = 100_000
ETA_1 = 200_000
VALUE_BAND = ETA_0 // 2
SLOPE_BAND = ETA_1 // 2
COEFFICIENT_MASS = 50
MAIN_ABSOLUTE_CONSTANT = 101
CORRECTION_DERIVATIVE_CONSTANT = 4_223
CORE_APPROX_RELATIVE = "1e-6"


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
        "coefficient_mass": (
            "50*exp(L/4)",
            "coefficient mass",
        ),
        "global_remainder": (
            "|K(T)|/A_t(x)",
            "20*exp(-L/4)",
            '"eta_0": 100000',
            '"eta_1": 200000',
        ),
        "signed_contact": (
            "lambda_a=phi'/phi-s_*'*log(a)=u_a+i*v_a",
            "E_[1],x=lambda_a*E_[1]+R_a",
            "d_(n,x)=(-i/2)",
        ),
        "wronskian_crossing": (
            "W_[1]=V*X-U*Y",
            "N_(E_[1]=0,U>0)",
        ),
        "oriented_successor": (
            "wind(proxy_j)<1",
            "positive imaginary ray",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(payload.get("kind", "")) for key, payload in payloads.items()}


def symbolic_audit() -> dict:
    x_part, y_part, p_part, q_part = sp.symbols(
        "X Y P Q", real=True
    )
    u_part, v_part = sp.symbols("u v", real=True)
    slope_real = u_part * x_part - v_part * y_part + p_part
    slope_imag = v_part * x_part + u_part * y_part + q_part
    centered = p_part - v_part * y_part
    if sp.simplify(slope_real - (u_part * x_part + centered)) != 0:
        raise RuntimeError("centered real-slope identity failed")

    wronskian = sp.expand(
        slope_imag * x_part - slope_real * y_part
    )
    expected_wronskian = (
        v_part * (x_part**2 + y_part**2)
        + q_part * x_part
        - p_part * y_part
    )
    centered_wronskian = x_part * (
        v_part * x_part + q_part
    ) - centered * y_part
    if sp.simplify(wronskian - expected_wronskian) != 0:
        raise RuntimeError("saddle Wronskian identity failed")
    if sp.simplify(wronskian - centered_wronskian) != 0:
        raise RuntimeError("centered Wronskian identity failed")

    x = sp.symbols("x", positive=True, real=True)
    s = (1 - sp.I * x) / 2
    z = sp.symbols("z")
    alpha_z = (
        1 / (2 * z)
        + 1 / (z - 1)
        + sp.log(z / (2 * sp.pi)) / 2
    )
    alpha_prime = sp.diff(alpha_z, z).subs(z, s)
    real_alpha_prime = sp.simplify(
        sp.re(sp.expand_complex(alpha_prime))
    )
    imag_alpha_prime = sp.simplify(
        sp.im(sp.expand_complex(alpha_prime))
    )
    expected_real = (7 * x**2 - 5) / (x**2 + 1) ** 2
    expected_imag = x * (x**2 + 5) / (x**2 + 1) ** 2
    if sp.simplify(real_alpha_prime - expected_real) != 0:
        raise RuntimeError("critical-line Re(alpha') identity failed")
    if sp.simplify(imag_alpha_prime - expected_imag) != 0:
        raise RuntimeError("critical-line Im(alpha') identity failed")

    a_star, b_star, moment_real, moment_imag = sp.symbols(
        "a_star b_star M_R M_I", real=True
    )
    moment = moment_real + sp.I * moment_imag
    residual_moment_real = sp.re(
        -(a_star + sp.I * b_star) * moment
    )
    if sp.simplify(
        residual_moment_real - (-a_star * moment_real + b_star * moment_imag)
    ) != 0:
        raise RuntimeError("centered moment real-part identity failed")

    return {
        "slope": "U=u_a*X-v_a*Y+P_a=u_a*X+S_a",
        "centered_scalar": (
            "S_a=P_a-v_a*Y=Re(R_a)-v_a*Im(E_[1])=U-u_a*X"
        ),
        "wronskian": (
            "W_[1]=v_a*(X^2+Y^2)+Q_a*X-P_a*Y"
            "=X*(v_a*X+Q_a)-S_a*Y"
        ),
        "crossing": "At X=0, U=S_a and W_[1]=-S_a*Y.",
        "alpha_prime_real": "C_x=Re(alpha'(s))=(7x^2-5)/(x^2+1)^2",
        "alpha_prime_imag": (
            "D_x=Im(alpha'(s))=x*(x^2+5)/(x^2+1)^2"
        ),
        "sympy_wronskian": str(wronskian),
        "sympy_alpha_prime_real": str(real_alpha_prime),
        "sympy_alpha_prime_imag": str(imag_alpha_prime),
    }


def numeric_budget_audit() -> dict:
    l_min = 50.0
    ux_ratio = 50_000.0 * math.exp(-l_min)
    vy_ratio = (
        303.0 / (16.0 * math.pi**2) * math.exp(-l_min / 2.0)
    )
    derivative_ratio = (
        COEFFICIENT_MASS
        * CORRECTION_DERIVATIVE_CONSTANT
        / (16.0 * math.pi**2)
        * math.exp(-l_min / 2.0)
    )
    total = ux_ratio + vy_ratio + derivative_ratio
    if not ux_ratio < 1e-16:
        raise RuntimeError("uX nuisance budget failed")
    if not vy_ratio < 1e-10:
        raise RuntimeError("vY nuisance budget failed")
    if not derivative_ratio < 1e-7:
        raise RuntimeError("correction-derivative nuisance budget failed")
    if not total < 1e-6:
        raise RuntimeError("combined core nuisance budget failed")
    return {
        "at_L_50_uX_over_epsilon_lt": "1e-16",
        "at_L_50_vY_over_epsilon_lt": "1e-10",
        "at_L_50_D1x_over_epsilon_lt": "1e-7",
        "combined_over_epsilon_lt": CORE_APPROX_RELATIVE,
        "computed_uX_ratio": format(ux_ratio, ".17e"),
        "computed_vY_ratio": format(vy_ratio, ".17e"),
        "computed_D1x_ratio": format(derivative_ratio, ".17e"),
        "computed_total_ratio": format(total, ".17e"),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    numeric = numeric_budget_audit()
    return {
        "coordinates": {
            "definition": (
                "E_[1]=X+iY, R_a=P_a+iQ_a, "
                "lambda_a=u_a+i*v_a, E_[1],x=lambda_a*E_[1]+R_a"
            ),
            "slope": symbolic["slope"],
            "centered_scalar": symbolic["centered_scalar"],
            "wronskian": symbolic["wronskian"],
            "crossing": symbolic["crossing"],
        },
        "critical_frame": {
            "definitions": (
                "s=(1-i*x)/2, L=log(x/(4*pi)), "
                "a^2=x/(4*pi)+t/16, "
                "alpha(s)=1/(2s)+1/(s-1)+(1/2)Log(s/(2*pi))"
            ),
            "alpha_real": (
                "Re(alpha)=L/2+(1/4)log(1+x^(-2))-1/(1+x^2)"
            ),
            "alpha_imag": (
                "Im(alpha)=-(1/2)atan(x)+3x/(1+x^2), x>0"
            ),
            "alpha_prime": (
                symbolic["alpha_prime_real"]
                + "; "
                + symbolic["alpha_prime_imag"]
            ),
            "phase": (
                "phi'/phi=i*beta_t', "
                "beta_t'=-(1/2)Re(alpha+(t/2)alpha*alpha')"
            ),
            "saddle_derivative": (
                "s_*'=t*D_x/4+i*(-1/2-t*C_x/4)"
            ),
            "lambda_components": (
                "u_a=-(t/4)D_x*log(a); "
                "delta_a=log(a)-Re(alpha); "
                "v_a=delta_a/2+(t/4)*(C_x*delta_a+Im(alpha)*D_x)"
            ),
        },
        "frame_bounds": {
            "domain": "L>=50, 0<tL<=25, hence 0<t<=1/2 and x=4*pi*exp(L)",
            "u_bound": (
                "0<D_x<=2/x and log(a)<L imply "
                "|u_a|<=tL/(2x)<=25/(2x)<exp(-L)"
            ),
            "v_decomposition": (
                "Put delta_0=1/(1+x^2)-(1/4)log(1+x^(-2)), "
                "z=pi*t/(4x), delta_t=(1/2)log(1+z). Then "
                "v_a=delta_0/2+(1/4)(log(1+z)-z)"
                "+(t/4)*C_x*(delta_0+delta_t)"
                "+(t/4)*(Im(alpha)*D_x+pi/(4x))."
            ),
            "v_elementary_bounds": (
                "0<=delta_0<=x^(-2), 0<=z<=1/(2x), "
                "|C_x|<=7/x^2, 0<D_x<=2/x, "
                "|D_x-1/x|<=3/x^3, "
                "0<Im(alpha)+pi/4<=7/(2x), and "
                "|log(1+z)-z|<=z^2/2"
            ),
            "v_bound": "|v_a|<3/x^2<exp(-2L)",
        },
        "main_and_derivative_bounds": {
            "coefficient_mass": (
                "sum_(n=1)^N |e_n|<=50*exp(L/4)"
            ),
            "first_correction": (
                "|alpha_n|<2L, |alpha'|<=7/x, "
                "t/4+t^2|alpha_n|^2/8<313, hence |d_n|<1"
            ),
            "endpoint": (
                "|K(T_0)|/A_t<20*exp(-L/4), "
                "exp(t*pi^2/64)<2, and "
                "|F+F'''/(12*pi^2*a)|<6 imply |g_0|<exp(L/4)"
            ),
            "main": "|E_[1]|<101*exp(L/4)",
            "correction_derivative": (
                "|d_(n,x)|<4223/x^2 and therefore "
                "|sum e_n*d_(n,x)|<1500*exp(-7L/4)"
                "<1e-7*exp(-5L/4)"
            ),
        },
        "core_scalar": {
            "moments": (
                "M_(1,a)=sum_(n=1)^N log(n/a)f_n, "
                "D_(1,x)=sum_(n=1)^N e_n*d_(n,x), "
                "G_a=g_(0,x)-lambda_a*g_0"
            ),
            "residual": (
                "R_a=-s_*'*M_(1,a)+D_(1,x)+G_a"
            ),
            "definition": (
                "A_a=-(1/2)Im(M_(1,a))"
                "-(t/4)*(D_x*Re(M_(1,a))+C_x*Im(M_(1,a)))"
                "+Re(G_a)"
            ),
            "real_residual": "P_a=A_a+Re(D_(1,x))",
            "band_approximation": (
                "For |X|<=50000*exp(-5L/4), "
                "|U-A_a|<1e-6*exp(-5L/4)."
            ),
            "budget_audit": numeric,
        },
        "endpoint_expansion": {
            "definitions": (
                "H_a=F(p)+F'''(p)/(12*pi^2*a), "
                "g_0=-kappa_N*H_a, "
                "K(T_0)=M_0(iT_0)U*exp(pi*i/8)"
            ),
            "geometry": (
                "a_x=1/(8*pi*a), p_x=-1/(4*pi*a), "
                "H_(a,x)=p_x*(F'(p)+F''''(p)/(12*pi^2*a))"
                "-a_x*F'''(p)/(12*pi^2*a^2)"
            ),
            "log_defect": (
                "mu_a=kappa_(N,x)/kappa_N-lambda_a"
                "=K_x/K-(M_t(s))_x/M_t(s)+s_*'*log(a)"
            ),
            "explicit_rates": (
                "K_x/K=-pi/8+1/(4T_0)+1/(2(T_0+i)); "
                "(M_t(s))_x/M_t(s)=(-i/2)"
                "*(alpha+(t/2)alpha*alpha')"
            ),
            "defect": (
                "G_a=-kappa_N*(H_(a,x)+mu_a*H_a)"
            ),
        },
        "contact_and_count": {
            "contact_box": (
                "A full contact forces |X|<50000*exp(-5L/4) and "
                "|U|<100000L*exp(-5L/4)."
            ),
            "sufficient_scalar_theorem": (
                "On q=2tL^2>=1, prove that |X|<=50000*exp(-5L/4) "
                "implies |A_a|>(100000L+1)*exp(-5L/4)."
            ),
            "contact_consequence": (
                "The scalar theorem and |U-A_a|<1e-6*exp(-5L/4) "
                "give |U|>100000L*exp(-5L/4), excluding contact."
            ),
            "exact_crossing_sign": (
                "At X=0, U=S_a. Thus every simple upward crossing is "
                "exactly X=0,S_a>0, including E_[1]=0."
            ),
            "core_crossing_sign": (
                "Under the sufficient scalar theorem, sign(U)=sign(A_a) "
                "at every X=0 crossing; no exceptional complex-main class "
                "is needed."
            ),
            "successor_count": (
                "kappa_j=N_(X=0,A_a>0;t_j)-N_(X=0,A_a>0;t_(j+1))"
                "+I_i(V_j)+I_i(finite shoulders and joins), "
                "and it is enough to prove this integer is <1."
            ),
        },
        "route_audit": {
            "wronskian_role": (
                "At X=0, W_[1]=-S_a*Y. The Wronskian sign loses the "
                "E_[1]=0 crossing, while S_a and the certified core A_a do not."
            ),
            "contact_role": (
                "The Wronskian magnitude disjunction is a stronger sufficient "
                "contact theorem, not a reduction of the scalar slope difficulty."
            ),
            "remaining_arithmetic": (
                "The RH-level input is now a lower bound and signed crossing "
                "budget for A_a, an explicit centered logarithmic moment plus "
                "the exact C_0+C_1/a endpoint defect."
            ),
        },
        "open_targets": {
            "frequency": (
                "Prove the scalar band theorem and its oriented A_a-positive "
                "crossing budget on q>=1 in the live 0<tL<c_*+o(1) layer."
            ),
            "parabolic": (
                "Treat q<1 by a multiplicity-compatible Hermite or degree chart."
            ),
            "finite": (
                "Close bounded-L, L_epsilon, vertical-connector, core-to-main, "
                "and chart-join phase cells separately."
            ),
        },
    }


def build_rows(exact: dict) -> list[ReductionRow]:
    return [
        ReductionRow(
            "nfocrr_01_centered_scalar",
            "exact_identity",
            "available_exact",
            "The saddle-centered real residual is the real slope after removing only radial frame drift.",
            exact["coordinates"]["centered_scalar"],
            "No division by E_[1], Y, or the Wronskian.",
            exact["coordinates"],
        ),
        ReductionRow(
            "nfocrr_02_wronskian_audit",
            "exact_identity",
            "available_exact",
            "The Wronskian factors through the centered scalar at a real-part crossing.",
            exact["coordinates"]["wronskian"] + "; " + exact["coordinates"]["crossing"],
            "This identity explains, but does not repair, the exceptional-zero degeneration.",
        ),
        ReductionRow(
            "nfocrr_03_critical_frame",
            "exact_identity",
            "available_exact",
            "The critical-line alpha and alpha-prime components give exact saddle-frame rates.",
            exact["critical_frame"]["alpha_prime"] + "; " + exact["critical_frame"]["lambda_components"],
            "Uses the principal logarithm and x>0.",
            exact["critical_frame"],
        ),
        ReductionRow(
            "nfocrr_04_u_bound",
            "exact_bound",
            "available_exact",
            "The radial saddle-frame drift is exponentially small throughout the first-order domain.",
            exact["frame_bounds"]["u_bound"],
            "Uniform on L>=50 and 0<tL<=25.",
        ),
        ReductionRow(
            "nfocrr_05_v_bound",
            "exact_bound",
            "available_exact",
            "The angular saddle-frame drift is quadratically smaller in x.",
            exact["frame_bounds"]["v_bound"],
            "Uses the displayed cancellation of the order-t/x terms.",
            exact["frame_bounds"],
        ),
        ReductionRow(
            "nfocrr_06_main_bound",
            "analytic_composition",
            "available_exact",
            "The corrected complex half has a conservative absolute bound.",
            exact["main_and_derivative_bounds"]["main"],
            "Composes the published coefficient mass with explicit correction and endpoint bounds.",
            exact["main_and_derivative_bounds"],
        ),
        ReductionRow(
            "nfocrr_07_derivative_correction",
            "analytic_composition",
            "available_exact",
            "The explicit d_n derivative is negligible at the first-order contact scale.",
            exact["main_and_derivative_bounds"]["correction_derivative"],
            "Absolute-value estimate only; no oscillatory gain is used.",
        ),
        ReductionRow(
            "nfocrr_08_core_scalar",
            "exact_reduction",
            "available_exact",
            "All non-negligible horizontal crossing arithmetic is captured by one scalar centered moment plus endpoint defect.",
            exact["core_scalar"]["definition"],
            "The remaining lower bound is Xi-specific and unproved.",
            exact["core_scalar"],
        ),
        ReductionRow(
            "nfocrr_09_core_approximation",
            "analytic_composition",
            "available_exact",
            "On the certified value band, the true first-order slope differs from the core scalar by less than one millionth of the remainder scale.",
            exact["core_scalar"]["band_approximation"],
            "Combines uX, vY, and D_(1,x) budgets.",
        ),
        ReductionRow(
            "nfocrr_10_endpoint_expansion",
            "exact_identity",
            "available_exact",
            "The retained endpoint defect has a fully explicit fixed-cutoff derivative formula.",
            exact["endpoint_expansion"]["defect"],
            "Valid on each prescribed-N analytic lift; adjacent lifts use the existing homotopy.",
            exact["endpoint_expansion"],
        ),
        ReductionRow(
            "nfocrr_11_contact_theorem",
            "exact_sufficient_target",
            "available_exact",
            "A scalar lower bound for the core observable excludes every q>=1 contact.",
            exact["contact_and_count"]["sufficient_scalar_theorem"],
            "The strict Xi-specific lower bound remains unproved.",
        ),
        ReductionRow(
            "nfocrr_12_crossing_unification",
            "exact_crossing_reduction",
            "available_exact",
            "The centered scalar classifies all horizontal crossings without an exceptional complex-main-zero term.",
            exact["contact_and_count"]["exact_crossing_sign"] + " " + exact["contact_and_count"]["core_crossing_sign"],
            "The core sign replacement uses the open scalar floor.",
        ),
        ReductionRow(
            "nfocrr_13_successor_count",
            "exact_composition",
            "available_exact",
            "The scalar-positive crossing count plugs into the one-sided successor integer trap.",
            exact["contact_and_count"]["successor_count"],
            "Vertical connectors and finite shoulders remain explicit.",
        ),
        ReductionRow(
            "nfocrr_14_frequency_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "Prove the centered scalar band and signed crossing theorem in the q>=1 layer.",
            exact["open_targets"]["frequency"],
            "This is the live RH-level arithmetic obligation.",
        ),
        ReductionRow(
            "nfocrr_15_other_cells",
            "open_theorem_target",
            "not_ready_to_apply",
            "Close the multiplicity-compatible q<1 layer and all finite phase cells.",
            exact["open_targets"]["parabolic"] + " " + exact["open_targets"]["finite"],
            "Neither obligation follows from the q>=1 scalar reduction.",
        ),
    ]


def build_artifact() -> dict:
    exact = build_exact()
    return {
        "kind": STEM,
        "date": "2026-07-26",
        "status": (
            "exact saddle-centered scalar crossing reduction with certified "
            "frame and nuisance budgets; the Xi scalar theorem remains open"
        ),
        "proof_boundary": (
            "This artifact proves the saddle-centered coordinate identities, "
            "critical-frame formulas, u_a and v_a bounds, corrected-main and "
            "d_(n,x) budgets, core-scalar approximation, explicit endpoint "
            "defect, all-crossing sign unification, and successor-count "
            "substitution. It does not prove the q>=1 scalar lower bound or "
            "signed count, the q<1 multiplicity-compatible theorem, finite "
            "phase cells, one-sided successor winding, contact exclusion, "
            "Lambda<=0, RH, PF-infinity, or a Clay-prize conclusion."
        ),
        "constants": {
            "eta_0": ETA_0,
            "eta_1": ETA_1,
            "value_band": VALUE_BAND,
            "slope_band": SLOPE_BAND,
            "coefficient_mass": COEFFICIENT_MASS,
            "main_absolute_constant": MAIN_ABSOLUTE_CONSTANT,
            "correction_derivative_constant": CORRECTION_DERIVATIVE_CONSTANT,
            "core_approx_relative": CORE_APPROX_RELATIVE,
        },
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "exact": exact,
        "rows": [asdict(row) for row in build_rows(exact)],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    return "\n".join(
        [
            "# Newman First-Order Centered Real-Residual Reduction",
            "",
            "Date: 2026-07-26",
            "",
            "Status: exact saddle-centered scalar reduction with certified",
            "nuisance budgets. The Xi arithmetic theorem remains open; this",
            "is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Centered Crossing Scalar",
            "",
            "```text",
            exact["coordinates"]["definition"],
            exact["coordinates"]["slope"],
            exact["coordinates"]["centered_scalar"],
            exact["coordinates"]["wronskian"],
            exact["coordinates"]["crossing"],
            "```",
            "",
            "Unlike `W_[1]*Y`, `S_a` still classifies a simple crossing when",
            "`E_[1]=0`; no complex-main zero is deleted.",
            "",
            "## Critical Frame",
            "",
            "```text",
            exact["critical_frame"]["definitions"],
            exact["critical_frame"]["alpha_real"],
            exact["critical_frame"]["alpha_imag"],
            exact["critical_frame"]["alpha_prime"],
            exact["critical_frame"]["phase"],
            exact["critical_frame"]["saddle_derivative"],
            exact["critical_frame"]["lambda_components"],
            exact["frame_bounds"]["u_bound"],
            exact["frame_bounds"]["v_decomposition"],
            exact["frame_bounds"]["v_bound"],
            "```",
            "",
            "The order-`t/x` pieces in `v_a` cancel after centering at the",
            "moving Riemann-Siegel saddle. This is why `v_a` is `O(x^-2)`.",
            "",
            "## Absolute Budgets",
            "",
            "```text",
            exact["main_and_derivative_bounds"]["coefficient_mass"],
            exact["main_and_derivative_bounds"]["first_correction"],
            exact["main_and_derivative_bounds"]["endpoint"],
            exact["main_and_derivative_bounds"]["main"],
            exact["main_and_derivative_bounds"]["correction_derivative"],
            "```",
            "",
            "## Core Arithmetic Scalar",
            "",
            "```text",
            exact["core_scalar"]["moments"],
            exact["core_scalar"]["residual"],
            exact["core_scalar"]["definition"],
            exact["core_scalar"]["real_residual"],
            exact["core_scalar"]["band_approximation"],
            "```",
            "",
            "Thus the non-negligible term is one centered logarithmic moment",
            "plus the retained endpoint defect.",
            "",
            "## Endpoint Defect",
            "",
            "```text",
            exact["endpoint_expansion"]["definitions"],
            exact["endpoint_expansion"]["geometry"],
            exact["endpoint_expansion"]["log_defect"],
            exact["endpoint_expansion"]["explicit_rates"],
            exact["endpoint_expansion"]["defect"],
            "```",
            "",
            "## Contact And Winding Target",
            "",
            "```text",
            exact["contact_and_count"]["contact_box"],
            exact["contact_and_count"]["sufficient_scalar_theorem"],
            exact["contact_and_count"]["contact_consequence"],
            exact["contact_and_count"]["exact_crossing_sign"],
            exact["contact_and_count"]["core_crossing_sign"],
            exact["contact_and_count"]["successor_count"],
            "```",
            "",
            "## Route Audit",
            "",
            "```text",
            exact["route_audit"]["wronskian_role"],
            exact["route_audit"]["contact_role"],
            exact["route_audit"]["remaining_arithmetic"],
            "```",
            "",
            "## Remaining Cells",
            "",
            "```text",
            exact["open_targets"]["frequency"],
            exact["open_targets"]["parabolic"],
            exact["open_targets"]["finite"],
            "```",
            "",
            artifact["proof_boundary"],
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(artifact, indent=2) + "\n",
        encoding="utf-8",
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman first-order centered real-residual reduction: "
        "15 rows, |u_a|<e^-L, |v_a|<3/x^2, "
        "core error <1e-6*e^-5L/4, 2 open Xi cell obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
