#!/usr/bin/env python3
"""Build the direct anchored-projection regime reduction."""

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
    "first_order_centered_direct_projection_regime_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_absolute_phase_anchor_reduction.json"
    ),
    "centered_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_real_residual_reduction.json"
    ),
    "adjacent_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_adjacent_chart_stability_certificate.json"
    ),
    "bulk_pair_gate": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_bulk_pair_transfer_gate.json"
    ),
    "complex_zero_scout": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_centered_complex_zero_scout.json"
    ),
}

L_MIN = 50
D_ABSOLUTE_CONSTANT = 2_189
D_X_CONSTANT = 4_223
NORMALIZED_D_CONSTANT = 2e-7
FRAME_CONSTANT = 2
SPLIT_RATIO_DENOMINATOR = 100
SPLIT_SQRT_CONSTANT = 0.9999


@dataclass(frozen=True)
class ProjectionRow:
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
        "phase_anchor": (
            "eta=f_1/|f_1|",
            "Z_0=P_0+r_0=E_[1]/f_1",
            "Z_A=-s_*'*P_1+r_A=C_a/f_1",
        ),
        "centered_reduction": (
            "25/(2x)<exp(-L)",
            "|v_a|<3/x^2",
            "|d_(n,x)|<4223/x^2",
        ),
        "adjacent_stability": (
            "|Delta X|<1100*exp(-5L/4)",
            "|A_(a,N+1)-A_(a,N)|<5000*exp(-7L/4)",
        ),
        "bulk_pair_gate": (
            "log(rho_n)<-(2/5)h_n",
            "|d_n|<2189/x<1/2",
        ),
        "complex_zero_scout": (
            "The two-variable corrected complex zero is numerically transverse.",
            "The centered core scalar retains the nonzero crossing slope.",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(value.get("kind", "")) for key, value in payloads.items()}


def symbolic_audit() -> dict[str, str]:
    amplitude, phase = sp.symbols(
        "amplitude phase", positive=True, real=True
    )
    d_n, d_1 = sp.symbols("d_n d_1")
    mod_n, mod_1 = sp.symbols(
        "mod_n mod_1", positive=True, real=True
    )
    relative_unit = (
        sp.exp(sp.I * phase)
        * (1 + d_n)
        / mod_n
        * mod_1
        / (1 + d_1)
    )
    reconstructed = sp.cancel(
        amplitude * mod_n / mod_1 * relative_unit
    )
    expected = sp.cancel(
        amplitude
        * sp.exp(sp.I * phase)
        * (1 + d_n)
        / (1 + d_1)
    )
    if sp.simplify(reconstructed - expected) != 0:
        raise RuntimeError("amplitude/unit-phase reconstruction failed")

    c, b, ell, cosine, sine = sp.symbols(
        "c b ell cosine sine", real=True
    )
    carrier = cosine + sp.I * sine
    projected_slope = sp.expand_complex(
        sp.re(-(c + sp.I * b) * ell * carrier)
    )
    if sp.simplify(
        projected_slope - ell * (-c * cosine + b * sine)
    ) != 0:
        raise RuntimeError("direct slope projection failed")

    log_a, distance, delta, time = sp.symbols(
        "log_a distance delta time", real=True
    )
    log_n = log_a - distance
    re_alpha = log_a - delta
    re_s_star = sp.Rational(1, 2) + time * re_alpha / 2
    exponent = time * log_n**2 / 4 - re_s_star * log_n
    edge_exponent = time * log_a**2 / 4 - re_s_star * log_a
    centered_exponent = (
        edge_exponent
        + (sp.Rational(1, 2) - time * delta / 2) * distance
        + time * distance**2 / 4
    )
    if sp.expand(exponent - centered_exponent) != 0:
        raise RuntimeError("saddle-amplitude centering failed")

    z_0, z_a, z_0x, nu, lam, d_over_f = sp.symbols(
        "Z_0 Z_A Z_0x nu lambda d_over_f"
    )
    derivative_identity = sp.expand(
        nu * z_0 + z_0x - (lam * z_0 + z_a + d_over_f)
    )
    solved = z_0x + (nu - lam) * z_0 - d_over_f
    if sp.expand(
        derivative_identity.subs(z_a, solved)
    ) != 0:
        raise RuntimeError("relative-shape derivative identity failed")

    x_part, y_part, p_part, q_part = sp.symbols(
        "mathsf_X mathsf_Y p q", real=True
    )
    ratio_projection = sp.re(
        (p_part + sp.I * q_part) * (x_part + sp.I * y_part)
    )
    if sp.expand_complex(
        ratio_projection - (p_part * x_part - q_part * y_part)
    ) != 0:
        raise RuntimeError("ordinary ratio projection failed")

    x, t = sp.symbols("x t", real=True)
    heat_model = t - x**2 / 2 + sp.I * x
    if sp.diff(heat_model, t) + sp.diff(heat_model, x, 2) != 0:
        raise RuntimeError("heat countermodel PDE failed")
    jacobian = sp.Matrix(
        [
            [sp.diff(sp.re(heat_model), x), sp.diff(sp.re(heat_model), t)],
            [sp.diff(sp.im(heat_model), x), sp.diff(sp.im(heat_model), t)],
        ]
    )
    if jacobian.subs({x: 0, t: 0}).det() != -1:
        raise RuntimeError("heat countermodel Jacobian failed")
    if sp.re(sp.diff(heat_model, x)).subs({x: 0, t: 0}) != 0:
        raise RuntimeError("heat countermodel directional contact failed")

    return {
        "carrier_decomposition": (
            "sigma_*=Re(s_*), omega_*=-Im(s_*); "
            "r_n=exp[t*log(n)^2/4-sigma_*log(n)]"
            "*|1+d_n|/|1+d_1|>0; "
            "zeta_n=phi*exp(i*omega_*log(n))"
            "*(1+d_n)/|1+d_n|, |zeta_n|=1; "
            "eta*q_n=r_n*zeta_n"
        ),
        "saddle_amplitude": (
            "u_n=log(a/n)=-ell_n, delta_a=log(a)-Re(alpha), "
            "B_a=1/2-t*delta_a/2, "
            "E_a=t*log(a)^2/4-Re(s_*)*log(a); "
            "r_n=exp[E_a+B_a*u_n+t*u_n^2/4]"
            "*|1+d_n|/|1+d_1|"
        ),
        "endpoint_projection": (
            "B_0=beta/(|M_t(s)|*|f_1|)>0, "
            "J_a=H_(a,x)+mu_a*H_a; "
            "Re(eta*r_0)=(-1)^N*B_0*T_0*H_a; "
            "Re(eta*r_A)=(-1)^N*B_0"
            "*[T_0*Re(J_a)-Im(J_a)]"
        ),
        "direct_projection": (
            "C_n=Re(zeta_n), S_n=Im(zeta_n), "
            "s_*'=c+i*b; "
            "mathsf_X:=Re(eta*Z_0)="
            "(-1)^N*B_0*T_0*H_a+sum_n r_n*C_n; "
            "mathsf_A:=Re(eta*Z_A)="
            "(-1)^N*B_0*[T_0*Re(J_a)-Im(J_a)]"
            "+sum_n ell_n*r_n*(-c*C_n+b*S_n)"
        ),
        "relative_derivative": (
            "nu_1=f_(1,x)/f_1=phi'/phi+d_(1,x)/(1+d_1); "
            "Z_A=Z_(0,x)+(nu_1-lambda_a)*Z_0-D_(1,x)/f_1"
        ),
        "anchored_derivative": (
            "W_0=eta*Z_0=E_[1]/|f_1|, "
            "W_A=eta*Z_A=C_a/|f_1|, "
            "rho_1=partial_x log|f_1|="
            "Re[d_(1,x)/(1+d_1)]; "
            "W_A=W_(0,x)+(rho_1-lambda_a)*W_0"
            "-D_(1,x)/|f_1|"
        ),
        "real_derivative": (
            "For W_0=mathsf_X+i*mathsf_Y and "
            "lambda_a=u_a+i*v_a, "
            "mathsf_A=partial_x mathsf_X"
            "+(rho_1-u_a)*mathsf_X+v_a*mathsf_Y"
            "-Re(D_(1,x))/|f_1|"
        ),
        "ordinary_identity": (
            "If Z_0!=0, write h=Z_A/Z_0=p+i*q and "
            "R=|Z_0|. Then mathsf_A=p*mathsf_X-q*mathsf_Y, "
            "mathsf_X^2+mathsf_Y^2=R^2, and "
            "q=-Delta_anchor/R^2."
        ),
        "ordinary_log_current": (
            "h=W_(0,x)/W_0+rho_1-lambda_a-D_(1,x)/E_[1]. "
            "With radial current r_x=partial_x log|W_0| and "
            "branch-free phase current "
            "theta_x=Im(W_(0,x)*conj(W_0))/|W_0|^2, "
            "p=r_x+rho_1-u_a-Re(D_(1,x)/E_[1]) and "
            "q=theta_x-v_a-Im(D_(1,x)/E_[1])."
        ),
        "ordinary_sufficient": (
            "If |mathsf_X|<=delta<tau<=R, |q|>=kappa, "
            "and |p|<=K, then "
            "|mathsf_A|>=kappa*sqrt(tau^2-delta^2)-K*delta."
        ),
        "exceptional_sufficient": (
            "If |Z_0|<tau, then "
            "|mathsf_A-partial_x mathsf_X|"
            "<2*exp(-L)*tau+2e-7*exp(-5L/4). "
            "Thus |partial_x mathsf_X| exceeding the target by this "
            "amount is sufficient."
        ),
        "zero_fiber": (
            "At Z_0=0, "
            "mathsf_A=partial_x mathsf_X"
            "-Re(D_(1,x))/|f_1|."
        ),
        "heat_countermodel": (
            "W_*(x,t)=t-x^2/2+i*x satisfies "
            "partial_t W_*=-partial_x^2 W_*. At (0,0), "
            "W_*=0, |W_*,x|=1, and "
            "det D_(x,t)(Re W_*,Im W_*)=-1, but "
            "partial_x Re W_*=0. Complex-zero simplicity and a nonzero "
            "two-variable Jacobian do not imply the required directional "
            "real-slope bound."
        ),
    }


def numeric_audit() -> dict[str, str]:
    x_min = 4 * math.pi * math.exp(L_MIN)
    d_upper = D_ABSOLUTE_CONSTANT / x_min
    f_lower = 1 - d_upper
    rho_upper = D_X_CONSTANT / x_min**2 / f_lower
    lambda_upper = math.exp(-L_MIN) + 3 / x_min**2
    frame_upper = rho_upper + lambda_upper
    if not d_upper < 0.5:
        raise RuntimeError("first-coefficient disk failed")
    if not frame_upper < FRAME_CONSTANT * math.exp(-L_MIN):
        raise RuntimeError("relative derivative-frame bound failed")

    delta_upper = 1 / (4 * x_min)
    time_upper = 0.5
    saddle_linear_lower = 0.5 - time_upper * delta_upper / 2
    if not saddle_linear_lower > 0.49:
        raise RuntimeError("saddle linear-amplitude coefficient failed")

    normalized_d = 1e-7 / f_lower
    if not normalized_d < NORMALIZED_D_CONSTANT:
        raise RuntimeError("normalized D correction failed")
    split_ratio = 100_000 * math.exp(-L_MIN / 4) / L_MIN
    split_target_ratio = (
        2
        * (100_000 * L_MIN + 1)
        * math.exp(-L_MIN / 4)
        / L_MIN
    )
    if not split_ratio < 1 / SPLIT_RATIO_DENOMINATOR:
        raise RuntimeError("explicit split width failed")
    if not math.sqrt(1 - split_ratio**2) > SPLIT_SQRT_CONSTANT:
        raise RuntimeError("explicit split square-root margin failed")
    return {
        "x_min": format(x_min, ".17e"),
        "d_1_upper": format(d_upper, ".17e"),
        "f_1_lower": format(f_lower, ".17e"),
        "rho_1_upper": format(rho_upper, ".17e"),
        "lambda_upper": format(lambda_upper, ".17e"),
        "rho_minus_lambda_upper": format(frame_upper, ".17e"),
        "two_exp_minus_L": format(
            FRAME_CONSTANT * math.exp(-L_MIN), ".17e"
        ),
        "saddle_B_lower": format(saddle_linear_lower, ".17e"),
        "normalized_D_constant_raw": format(normalized_d, ".17e"),
        "split_delta_over_tau_upper": format(split_ratio, ".17e"),
        "split_target_over_tau_upper": format(
            split_target_ratio, ".17e"
        ),
        "split_sqrt_lower": format(
            math.sqrt(1 - split_ratio**2), ".17e"
        ),
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    numeric = numeric_audit()
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a)"
        ),
        "carrier_decomposition": symbolic["carrier_decomposition"],
        "saddle_amplitude": symbolic["saddle_amplitude"],
        "saddle_bounds": (
            "0<=delta_a<1/(4x), 0<=t<=1/2, hence "
            "B_a=1/2-t*delta_a/2>49/100. "
            "The amplitudes are explicit and positive; all unresolved "
            "oscillation is in the unit carriers zeta_n."
        ),
        "endpoint_projection": symbolic["endpoint_projection"],
        "direct_projection": symbolic["direct_projection"],
        "relative_derivative": symbolic["relative_derivative"],
        "anchored_derivative": symbolic["anchored_derivative"],
        "real_derivative": symbolic["real_derivative"],
        "derivative_frame_bound": (
            "|rho_1|<8446/x^2, |u_a|<exp(-L), "
            "|v_a|<3/x^2, hence "
            "|rho_1-lambda_a|<2*exp(-L). "
            "Also |D_(1,x)|/|f_1|"
            "<2e-7*exp(-5L/4)."
        ),
        "ordinary_identity": symbolic["ordinary_identity"],
        "ordinary_log_current": symbolic["ordinary_log_current"],
        "ordinary_current_error": (
            "On |Z_0|>=tau, "
            "|D_(1,x)/E_[1]|"
            "<2e-7*exp(-5L/4)/tau, while "
            "|rho_1-u_a|<2*exp(-L) and "
            "|v_a|<exp(-2L)."
        ),
        "ordinary_sufficient": symbolic["ordinary_sufficient"],
        "explicit_split": (
            "Choose tau_L=L*exp(-L). Since |f_1|>1/2, "
            "delta_L/tau_L"
            "<100000*exp(-L/4)/L<1/100 and "
            "sqrt(tau_L^2-delta_L^2)>0.9999*tau_L. "
            "Also A_L/tau_L"
            "<2*(100000L+1)*exp(-L/4)/L. "
            "Thus the ordinary branch closes if "
            "0.9999*kappa_L"
            "-K_L*100000*exp(-L/4)/L"
            ">2*(100000L+1)*exp(-L/4)/L."
        ),
        "exceptional_sufficient": symbolic["exceptional_sufficient"],
        "zero_fiber": symbolic["zero_fiber"],
        "heat_countermodel": symbolic["heat_countermodel"],
        "chart_scope": (
            "mathsf_X and mathsf_A inherit the adjacent real-projection "
            "bounds, but h=Z_A/Z_0 and its real/imaginary parts do not "
            "have a certified complex adjacent-chart bound. Use h only "
            "inside one prescribed canonical chart and transfer the final "
            "real inequalities."
        ),
        "live_target": (
            "With delta_L=50000*exp(-5L/4)/|f_1| and "
            "A_L=(100000L+1)*exp(-5L/4)/|f_1|, choose "
            "tau_L=L*exp(-L). On |Z_0|=R>=tau_L prove the direct "
            "Xi-specific inequality "
            "|Im(Z_A/Z_0)|*sqrt(R^2-delta_L^2)"
            "-|Re(Z_A/Z_0)|*delta_L>A_L. Uniform kappa_L,K_L "
            "bounds satisfying the displayed coarse condition are one "
            "stronger sufficient route, not a required formulation. "
            "On |Z_0|<tau_L prove the signed direct "
            "bound for partial_x mathsf_X, including the displayed "
            "frame and D corrections. Then count mathsf_A-positive "
            "crossings below one successor turn."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 boundary still requires a separate "
            "multiplicity-compatible parabolic/Hermite chart; neither "
            "ordinary ratio division nor a uniform zero-fiber slope may "
            "be imposed at t=0."
        ),
        "numeric_audit": numeric,
    }


def build_rows(exact: dict) -> list[ProjectionRow]:
    return [
        ProjectionRow(
            "nfocdprr_01_carriers",
            "exact_normal_form",
            "ready_to_apply",
            "Every relative coefficient splits into a positive amplitude and a branch-free unit carrier.",
            exact["carrier_decomposition"],
            "No argument branch or phase monotonicity is used.",
        ),
        ProjectionRow(
            "nfocdprr_02_saddle_amplitude",
            "exact_normal_form",
            "ready_to_apply",
            "The positive amplitudes have one exact saddle-distance form.",
            exact["saddle_amplitude"] + "; " + exact["saddle_bounds"],
            "This is an amplitude identity, not an oscillatory lower bound.",
            exact["numeric_audit"],
        ),
        ProjectionRow(
            "nfocdprr_03_endpoint_projection",
            "exact_identity",
            "ready_to_apply",
            "Both endpoint real projections are explicit in the same branch-free frame.",
            exact["endpoint_projection"],
            "The endpoint scalar J_a remains complex through mu_a.",
        ),
        ProjectionRow(
            "nfocdprr_04_direct_pair",
            "exact_normal_form",
            "ready_to_apply",
            "The live value and centered-slope observables are direct cosine/sine carrier sums.",
            exact["direct_projection"],
            "The unit-carrier cancellation remains the arithmetic problem.",
        ),
        ProjectionRow(
            "nfocdprr_05_relative_derivative",
            "exact_identity",
            "ready_to_apply",
            "The centered relative shape is the derivative of the value shape plus an explicit frame correction.",
            exact["relative_derivative"],
            "Division by f_1 is safe; division by Z_0 is not used.",
        ),
        ProjectionRow(
            "nfocdprr_06_anchored_derivative",
            "exact_identity",
            "ready_to_apply",
            "The anchored centered projection differs from the derivative of the anchored value by two certified nuisance terms.",
            exact["anchored_derivative"] + "; " + exact["real_derivative"],
            "Exact on one prescribed first-order chart.",
        ),
        ProjectionRow(
            "nfocdprr_07_frame_bound",
            "asymptotic_certificate",
            "certified",
            "The derivative-frame correction is exponentially smaller than unit scale.",
            exact["derivative_frame_bound"],
            "Valid on L>=50 and 0<=tL<=25.",
            exact["numeric_audit"],
        ),
        ProjectionRow(
            "nfocdprr_08_ordinary_identity",
            "exact_reduction",
            "ready_to_apply",
            "Away from the complex-zero fiber the direct slope is controlled by the complex shape ratio.",
            (
                exact["ordinary_identity"]
                + " "
                + exact["ordinary_log_current"]
            ),
            "The ratio exists only when Z_0 is nonzero.",
        ),
        ProjectionRow(
            "nfocdprr_09_ordinary_condition",
            "conditional_inequality",
            "ready_to_apply",
            "Explicit real and imaginary ratio bounds would close the ordinary regime.",
            (
                exact["ordinary_sufficient"]
                + " "
                + exact["ordinary_current_error"]
                + " "
                + exact["explicit_split"]
            ),
            "The required Xi-specific ratio bounds are not proved here.",
        ),
        ProjectionRow(
            "nfocdprr_10_exceptional_condition",
            "asymptotic_reduction",
            "certified",
            "Near the complex-zero fiber the theorem reduces directly to directional real transversality.",
            exact["exceptional_sufficient"],
            "This does not furnish the missing lower bound.",
        ),
        ProjectionRow(
            "nfocdprr_11_zero_fiber",
            "exact_identity",
            "ready_to_apply",
            "At a complex-main zero the centered scalar is the normalized real crossing derivative up to the tiny D correction.",
            exact["zero_fiber"],
            "No Wronskian, determinant, phase ratio, or division remains.",
        ),
        ProjectionRow(
            "nfocdprr_12_heat_guard",
            "countermodel",
            "guard_validated",
            "Complex-zero simplicity, derivative magnitude, and a nonzero two-variable Jacobian do not imply the required directional real slope.",
            exact["heat_countermodel"],
            "Exact generic backward-heat model, not an Xi counterexample.",
        ),
        ProjectionRow(
            "nfocdprr_13_chart_guard",
            "covariance_guard",
            "guard_validated",
            "The ratio theorem must be proved canonically and transferred only through the real projections.",
            exact["chart_scope"],
            "No complex chart-stability statement is promoted.",
        ),
        ProjectionRow(
            "nfocdprr_14_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 problem is now an explicit ordinary ratio theorem plus directional exceptional transversality.",
            exact["live_target"],
            "Neither branch is proved arithmetically.",
        ),
        ProjectionRow(
            "nfocdprr_15_q_lt_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1 multiplicity-compatible chart remains separate.",
            exact["q_lt_1_target"],
            "Endpoint simplicity is not assumed.",
        ),
        ProjectionRow(
            "nfocdprr_16_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "The exact regime split is kept separate from its missing Xi-specific inequalities and RH.",
            (
                "direct projection normal form and conditional split "
                "!= anchored projection lower bound"
            ),
            (
                "No q>=1 ratio/transversality theorem, q<1 theorem, "
                "finite phase closure, contact exclusion, Lambda<=0, "
                "PF-infinity, RH, or Clay-prize conclusion."
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
            "exact direct anchored-projection saddle normal form and "
            "ordinary/exceptional transversality reduction; both "
            "Xi-specific lower-bound branches remain open"
        ),
        "proof_boundary": (
            "This artifact proves the positive-amplitude/unit-carrier "
            "normal form, explicit endpoint projections, relative and "
            "anchored derivative identities, the small derivative-frame "
            "bound, the ordinary shape-ratio inequality, the exceptional "
            "directional reduction, and an exact backward-heat guard. It "
            "does not prove the q>=1 ratio or directional-transversality "
            "bounds, their signed crossing count, the q<1 "
            "multiplicity-compatible theorem, finite phase cells, "
            "one-sided successor winding, contact exclusion, Lambda<=0, "
            "RH, PF-infinity, or a Clay-prize conclusion."
        ),
        "sources": {
            key: str(path.relative_to(REPO_ROOT)).replace("\\", "/")
            for key, path in SOURCES.items()
        },
        "source_sha256": source_hashes(),
        "source_audit": source_audit(),
        "constants": {
            "L_min": L_MIN,
            "d_absolute": D_ABSOLUTE_CONSTANT,
            "d_x": D_X_CONSTANT,
            "normalized_D": NORMALIZED_D_CONSTANT,
            "frame": FRAME_CONSTANT,
            "split_ratio_denominator": SPLIT_RATIO_DENOMINATOR,
            "split_sqrt": SPLIT_SQRT_CONSTANT,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    numeric = exact["numeric_audit"]
    return "\n".join(
        [
            "# Newman First-Order Centered Direct-Projection Regime Reduction",
            "",
            "Date: 2026-07-27",
            "",
            "Status: exact direct saddle-phase coordinates and an",
            "ordinary/exceptional theorem split. The Xi-specific lower",
            "bounds remain open; this is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Direct Unit Carriers",
            "",
            "The relative coefficient chain has the branch-free split",
            "",
            "```text",
            exact["carrier_decomposition"],
            "```",
            "",
            "In the saddle-distance coordinate `u_n=log(a/n)`,",
            "",
            "```text",
            exact["saddle_amplitude"],
            exact["saddle_bounds"],
            "```",
            "",
            f"The numerical audit at `L=50` gives",
            f"`B_a>{numeric['saddle_B_lower']}`. Positivity and amplitude",
            "ordering are explicit; the unresolved arithmetic is entirely",
            "in the unit carriers.",
            "",
            "The endpoint enters the same frame without an argument branch:",
            "",
            "```text",
            exact["endpoint_projection"],
            "```",
            "",
            "Therefore the two live real observables are",
            "",
            "```text",
            exact["direct_projection"],
            "```",
            "",
            "## Derivative Identity",
            "",
            "Differentiating `E_[1]=f_1 Z_0` and comparing with the",
            "centered decomposition gives",
            "",
            "```text",
            exact["relative_derivative"],
            exact["anchored_derivative"],
            exact["real_derivative"],
            "```",
            "",
            "The complete nuisance coefficient is tiny:",
            "",
            "```text",
            exact["derivative_frame_bound"],
            "```",
            "",
            "This is stronger bookkeeping than the previous qualitative",
            "ordinary/exceptional split.",
            "",
            "## Ordinary Regime",
            "",
            "For `Z_0!=0`,",
            "",
            "```text",
            exact["ordinary_identity"],
            exact["ordinary_log_current"],
            "```",
            "",
            "Consequently,",
            "",
            "```text",
            exact["ordinary_sufficient"],
            exact["ordinary_current_error"],
            exact["explicit_split"],
            "```",
            "",
            "The missing input is now explicit: a lower bound for the",
            "phase current, with enough control of the radial current.",
            "Since the imaginary ratio is the normalized",
            "determinant, this does not resurrect the rejected generic",
            "Wronskian route; it identifies exactly what Xi-specific phase",
            "information would have to add.",
            "",
            "At `L=50`, the coarse ratios are",
            f"`delta_L/tau_L<{numeric['split_delta_over_tau_upper']}`",
            f"and `A_L/tau_L<{numeric['split_target_over_tau_upper']}`.",
            "",
            "## Exceptional Regime",
            "",
            "When `|Z_0|<tau`, no ratio is used:",
            "",
            "```text",
            exact["exceptional_sufficient"],
            exact["zero_fiber"],
            "```",
            "",
            (
                "At the zero fiber the problem is directional real "
                "transversality, not merely simplicity of the complex zero."
            ),
            "The exact guard is",
            "",
            "```text",
            exact["heat_countermodel"],
            "```",
            "",
            "Thus neither `|Z_(0,x)|>0` nor a nonzero `(x,t)` Jacobian",
            "closes the exceptional branch. The direction selected by the",
            "physical real projection must be controlled.",
            "",
            "## Live Theorem",
            "",
            "```text",
            exact["live_target"],
            "```",
            "",
            "Chart use is restricted by",
            "",
            "```text",
            exact["chart_scope"],
            "```",
            "",
            "The q<1 obligation remains",
            "",
            "```text",
            exact["q_lt_1_target"],
            "```",
            "",
            "This artifact proves no ratio lower bound, exceptional",
            "directional transversality, signed successor count, q<1",
            "chart, contact exclusion, `Lambda<=0`, PF-infinity, RH, or",
            "Clay-prize conclusion.",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built Newman centered direct-projection regime reduction: "
        "16 rows, 2 exact projection normal forms, "
        "1 heat-direction countermodel, 2 open Xi regimes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
