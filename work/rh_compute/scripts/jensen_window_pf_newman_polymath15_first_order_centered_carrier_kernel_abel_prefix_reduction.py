#!/usr/bin/env python3
"""Build the centered carrier-kernel and Abel-prefix reduction."""

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
    "carrier_kernel_abel_prefix_reduction"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCES = {
    "direct_projection": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "direct_projection_regime_reduction.json"
    ),
    "centered_reduction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_residual_reduction.json"
    ),
    "phase_anchor": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "absolute_phase_anchor_reduction.json"
    ),
    "bulk_pair": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "bulk_pair_transfer_gate.json"
    ),
    "wronskian_crossing": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "first_order_wronskian_crossing_reduction.json"
    ),
    "component_wronskian": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_component_wronskian_gate.json"
    ),
}

L_MIN = 50
D_X_CONSTANT = 4_223
TERMINAL_BAND_CONSTANT = 25_000
TERMINAL_RELATIVE_AT_L_MIN = 2e-35


@dataclass(frozen=True)
class KernelRow:
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
        "direct_projection": (
            "eta*q_n=r_n*zeta_n",
            "mathsf_A=partial_x mathsf_X",
            "Choose tau_L=L*exp(-L)",
        ),
        "centered_reduction": (
            "|d_(n,x)|<4223/x^2",
            "|v_a|<3/x^2",
        ),
        "phase_anchor": (
            "Z_0=P_0+r_0=E_[1]/f_1",
            "Z_A=-s_*'*P_1+r_A=C_a/f_1",
        ),
        "bulk_pair": (
            "F_k=sum_(n<=k)f_n",
            "log(rho_n)<-(2/5)h_n",
            "det(B_n)>=h_n^2/16",
        ),
        "wronskian_crossing": (
            "W_[1]=V*X-U*Y",
            "N_up=N_(X=0,E_[1]!=0",
        ),
        "component_wronskian": (
            "W_E=Im(E'*conj(E))",
            "inertia (1,1,m-2)",
        ),
    }
    for key, required in markers.items():
        text = json.dumps(payloads[key], sort_keys=True)
        for marker in required:
            if marker not in text:
                raise RuntimeError(f"{key} source marker missing: {marker}")
    return {key: str(value.get("kind", "")) for key, value in payloads.items()}


def symbolic_audit() -> dict[str, str]:
    x_height, heat_time, cutoff_index = sp.symbols(
        "x heat_time cutoff_index", real=True
    )
    pi_constant = sp.pi
    shifted_height = x_height / 2 + pi_constant * heat_time / 8
    saddle_squared = shifted_height / (2 * pi_constant)
    if sp.simplify(
        saddle_squared - (x_height / (4 * pi_constant) + heat_time / 16)
    ) != 0:
        raise RuntimeError("Riemann-Siegel saddle pi provenance failed")
    cell_left = 4 * pi_constant * (
        cutoff_index**2 - heat_time / 16
    )
    cell_right = 4 * pi_constant * (
        (cutoff_index + 1) ** 2 - heat_time / 16
    )
    if sp.simplify(
        cell_right - cell_left
        - 4 * pi_constant * (2 * cutoff_index + 1)
    ) != 0:
        raise RuntimeError("cutoff-cell pi provenance failed")

    c, b, u_1, u_2, r_1, r_2, theta = sp.symbols(
        "c b u_1 u_2 r_1 r_2 theta", real=True
    )
    s_prime = c + sp.I * b
    z_1 = r_1 * (sp.cos(theta) + sp.I * sp.sin(theta))
    z_2 = r_2
    bulk_sum = z_1 + z_2
    bulk_moment = u_1 * z_1 + u_2 * z_2
    kernel = sp.expand_complex(
        s_prime * bulk_moment * sp.conjugate(bulk_sum)
    )
    expected_p = (
        c * (u_1 * r_1**2 + u_2 * r_2**2)
        + r_1
        * r_2
        * (
            c * (u_1 + u_2) * sp.cos(theta)
            - b * (u_1 - u_2) * sp.sin(theta)
        )
    )
    expected_q = (
        b * (u_1 * r_1**2 + u_2 * r_2**2)
        + r_1
        * r_2
        * (
            b * (u_1 + u_2) * sp.cos(theta)
            + c * (u_1 - u_2) * sp.sin(theta)
        )
    )
    if sp.trigsimp(sp.re(kernel) - expected_p) != 0:
        raise RuntimeError("pairwise radial kernel failed")
    if sp.trigsimp(sp.im(kernel) - expected_q) != 0:
        raise RuntimeError("pairwise tangential kernel failed")

    z1, z2, z3, u3, h1, h2 = sp.symbols(
        "z1 z2 z3 u3 h1 h2"
    )
    u1 = u3 + h1 + h2
    u2 = u3 + h2
    total = z1 + z2 + z3
    prefix1 = z1
    prefix2 = z1 + z2
    moment = u1 * z1 + u2 * z2 + u3 * z3
    abel = u3 * total + h1 * prefix1 + h2 * prefix2
    if sp.expand(moment - abel) != 0:
        raise RuntimeError("finite Abel prefix identity failed")

    endpoint, endpoint_slope = sp.symbols("endpoint endpoint_slope")
    w0 = endpoint + total
    wn = (
        endpoint_slope
        + s_prime * moment
    )
    terminal = endpoint_slope - s_prime * u3 * endpoint
    wn_prefix = (
        s_prime * u3 * w0
        + terminal
        + s_prime * (h1 * prefix1 + h2 * prefix2)
    )
    if sp.expand(wn - wn_prefix) != 0:
        raise RuntimeError("endpoint-complete prefix identity failed")

    x, y, vr, vi, fr, fi, h, u = sp.symbols(
        "X Y V_R V_I F_R F_I h u", real=True
    )
    w = x + sp.I * y
    v = vr + sp.I * vi
    f = fr + sp.I * fi
    wa = s_prime * u * w + v + s_prime * h * f
    contact_scalar = vr - b * u * y + h * (c * fr - b * fi)
    if sp.simplify(
        sp.expand_complex(
            sp.re(wa) - contact_scalar - c * u * x
        )
    ) != 0:
        raise RuntimeError("all-fiber crossing scalar failed")

    fw = sp.expand_complex(f * sp.conjugate(w))
    vw = sp.expand_complex(v * sp.conjugate(w))
    radius_sq = x**2 + y**2
    prefix_p = (
        c * u * radius_sq
        + sp.re(vw)
        + h * (c * sp.re(fw) - b * sp.im(fw))
    )
    prefix_q = (
        b * u * radius_sq
        + sp.im(vw)
        + h * (b * sp.re(fw) + c * sp.im(fw))
    )
    direct_kernel = sp.expand_complex(wa * sp.conjugate(w))
    if sp.simplify(sp.re(direct_kernel) - prefix_p) != 0:
        raise RuntimeError("prefix radial-current identity failed")
    if sp.simplify(sp.im(direct_kernel) - prefix_q) != 0:
        raise RuntimeError("prefix tangential-current identity failed")

    rho_minus_u, v_a, d_real = sp.symbols(
        "rho_minus_u v_a d_real", real=True
    )
    x_derivative = (
        sp.re(wa)
        - rho_minus_u * x
        - v_a * y
        + d_real
    )
    expected_x_derivative = (
        contact_scalar
        + (c * u - rho_minus_u) * x
        - v_a * y
        + d_real
    )
    if sp.simplify(
        sp.expand_complex(x_derivative - expected_x_derivative)
    ) != 0:
        raise RuntimeError("exceptional derivative scalar failed")

    ordinary = ordinary_countermodel()
    if ordinary["kernel_zero"] != "0":
        raise RuntimeError("ordinary countermodel kernel did not vanish")
    zero_fiber = zero_fiber_countermodel()
    if zero_fiber["W_0"] != "0" or zero_fiber["Re_W_A"] != "0":
        raise RuntimeError("zero-fiber countermodel failed")

    return {
        "pi_provenance": (
            "pi is the universal constant in the completed-zeta "
            "normalization xi(s)=s*(s-1)*pi^(-s/2)*Gamma(s/2)*zeta(s)/2. "
            "Equivalently, with theta(y)=sum_(m in Z)exp(-pi*m^2*y), "
            "pi^(-s/2)*Gamma(s/2)*zeta(s)="
            "integral_0^infinity[(theta(y)-1)/2]*y^(s/2)*dy/y. "
            "For critical height T=x/2 and the heat-shifted height "
            "T_0=T+pi*t/8, the Riemann-Siegel saddle gives "
            "a^2=T_0/(2*pi)=x/(4*pi)+t/16. Hence at t=0, "
            "L=log(x/(4*pi))=log(a^2), and the fixed-t N-cell "
            "N<=a<N+1 has endpoints "
            "x_N=4*pi*(N^2-t/16), "
            "x_(N+1)=4*pi*((N+1)^2-t/16), so "
            "Delta x=4*pi*(2*N+1). The later 2*pi phase turn is the "
            "period exp(i*(theta+2*pi))=exp(i*theta). The prefix polygon "
            "is only the complex partial-sum path and defines none of "
            "these appearances of pi."
        ),
        "coordinates": (
            "z_n=r_n*zeta_n, u_n=log(a/n)>=0, "
            "S=sum_(n=1)^N z_n, U=sum_(n=1)^N u_n*z_n, "
            "e=eta*r_0=(-1)^N*B_0*(T_0+i)*H_a, "
            "g=eta*r_A=(-1)^N*B_0*(T_0+i)*J_a, "
            "W_0=e+S, W_A=g+s_*'*U, s_*'=c+i*b"
        ),
        "bulk_pairwise": (
            "For theta_nm=arg(z_n)-arg(z_m), "
            "P_B=Re(s_*'*U*conj(S))="
            "c*sum_n u_n*r_n^2"
            "+sum_(n<m)r_n*r_m*[c*(u_n+u_m)*cos(theta_nm)"
            "-b*(u_n-u_m)*sin(theta_nm)]; "
            "Q_B=Im(s_*'*U*conj(S))="
            "b*sum_n u_n*r_n^2"
            "+sum_(n<m)r_n*r_m*[b*(u_n+u_m)*cos(theta_nm)"
            "+c*(u_n-u_m)*sin(theta_nm)]"
        ),
        "endpoint_pairwise": (
            "P=Re(W_A*conj(W_0))="
            "Re(g*conj(e))+P_B"
            "+sum_n Re(g*conj(z_n)+s_*'*u_n*z_n*conj(e)); "
            "Q=Im(W_A*conj(W_0))="
            "Im(g*conj(e))+Q_B"
            "+sum_n Im(g*conj(z_n)+s_*'*u_n*z_n*conj(e)); "
            "g*conj(e)=B_0^2*(T_0^2+1)*H_a*J_a"
        ),
        "abel_prefix": (
            "F_k=sum_(n=1)^k z_n, h_k=log((k+1)/k)>0, "
            "u_n=u_N+sum_(k=n)^(N-1)h_k; "
            "U=u_N*S+sum_(k=1)^(N-1)h_k*F_k"
        ),
        "endpoint_complete_prefix": (
            "V_N=g-s_*'*u_N*e; "
            "W_A=s_*'*u_N*W_0+V_N"
            "+s_*'*sum_(k=1)^(N-1)h_k*F_k"
        ),
        "prefix_currents": (
            "For R=|W_0|, "
            "P=c*u_N*R^2+Re(V_N*conj(W_0))"
            "+sum_k h_k*[c*Re(F_k*conj(W_0))"
            "-b*Im(F_k*conj(W_0))]; "
            "Q=b*u_N*R^2+Im(V_N*conj(W_0))"
            "+sum_k h_k*[b*Re(F_k*conj(W_0))"
            "+c*Im(F_k*conj(W_0))]"
        ),
        "contact_scalar": (
            "For W_0=mathsf_X+i*mathsf_Y define "
            "mathcal_C_N=Re(V_N)-b*u_N*mathsf_Y"
            "+sum_k h_k*[c*Re(F_k)-b*Im(F_k)]. "
            "Then mathsf_A=Re(W_A)=mathcal_C_N"
            "+c*u_N*mathsf_X. At mathsf_X=0, "
            "mathsf_A=mathcal_C_N and "
            "Q=-mathsf_Y*mathcal_C_N, including mathsf_Y=0."
        ),
        "exceptional_derivative": (
            "With k_1=rho_1-u_a, "
            "partial_x mathsf_X=mathcal_C_N"
            "+(c*u_N-k_1)*mathsf_X-v_a*mathsf_Y"
            "+Re(D_(1,x))/|f_1|. At W_0=0, "
            "partial_x mathsf_X=mathcal_C_N"
            "+Re(D_(1,x))/|f_1|."
        ),
    }


def ordinary_countermodel() -> dict[str, str]:
    cosine = sp.Rational(-35, 37)
    sine = sp.Rational(-12, 37)
    unit = cosine + sp.I * sine
    s_prime = sp.Rational(13, 1056) - sp.I / 2
    z = (
        sp.Integer(1),
        sp.Rational(4, 5) * unit,
        sp.Rational(7, 10) * unit,
    )
    u = (sp.Integer(3), sp.Integer(2), sp.Integer(1))
    bulk_sum = sp.simplify(sum(z))
    bulk_moment = sp.simplify(sum(ui * zi for ui, zi in zip(u, z)))
    kernel = sp.simplify(
        sp.im(s_prime * bulk_moment * sp.conjugate(bulk_sum))
    )

    def q_at(cosine_value: sp.Expr) -> sp.Expr:
        unit_value = cosine_value
        zv = (
            sp.Integer(1),
            sp.Rational(4, 5) * unit_value,
            sp.Rational(7, 10) * unit_value,
        )
        sv = sum(zv)
        uv = sum(ui * zi for ui, zi in zip(u, zv))
        return sp.simplify(
            sp.im(s_prime * uv * sp.conjugate(sv))
        )

    return {
        "s_prime": "13/1056-i/2",
        "amplitudes": "1,4/5,7/10",
        "distances": "3,2,1",
        "shared_lower_phase_cos": "-35/37",
        "shared_lower_phase_sin": "-12/37",
        "aligned_Q": str(q_at(sp.Integer(1))),
        "opposed_Q": str(q_at(sp.Integer(-1))),
        "kernel_zero": str(kernel),
        "radius_squared_at_zero": str(
            sp.simplify(bulk_sum * sp.conjugate(bulk_sum))
        ),
        "interpretation": (
            "Strictly decreasing positive amplitudes and distances allow "
            "Q<0, Q>0, and Q=0 with |W_0|^2=61/148>0, "
            "even for c>0 and b=-1/2."
        ),
    }


def zero_fiber_countermodel() -> dict[str, str]:
    s_prime = sp.Rational(1, 16) - sp.I / 2
    omega = (-8 + sp.I) / sp.sqrt(65)
    z = (
        omega,
        -sp.Rational(3, 5) * omega,
        -sp.Rational(2, 5) * omega,
    )
    u = (sp.Integer(3), sp.Integer(2), sp.Integer(1))
    w0 = sp.simplify(sum(z))
    moment = sp.simplify(sum(ui * zi for ui, zi in zip(u, z)))
    wa = sp.simplify(s_prime * moment)
    return {
        "s_prime": "1/16-i/2",
        "global_unit_anchor": "(-8+i)/sqrt(65)",
        "amplitudes": "1,3/5,2/5",
        "distances": "3,2,1",
        "W_0": str(w0),
        "U": str(moment),
        "W_A": str(wa),
        "Re_W_A": str(sp.re(wa)),
        "abs_W_A_squared": str(sp.simplify(wa * sp.conjugate(wa))),
        "interpretation": (
            "At an exact complex-main zero, strict positive amplitude and "
            "distance ordering permit a nonzero purely imaginary W_A and "
            "zero physical real-direction slope. The Xi anchor and endpoint "
            "must enter any successful lower bound."
        ),
    }


def numeric_audit() -> dict[str, str]:
    l_value = float(L_MIN)
    x_min = 4 * math.pi * math.exp(l_value)
    carrier_rotation_lower = (
        0.5 * math.log(2)
        - 4 * D_X_CONSTANT / x_min**2
    )
    if not carrier_rotation_lower > 1 / 3:
        raise RuntimeError("relative carrier rotation lower bound failed")
    terminal_ratio = math.exp(-1.5 * l_value) / (4 * l_value)
    if not terminal_ratio < TERMINAL_RELATIVE_AT_L_MIN:
        raise RuntimeError("terminal band ratio failed")
    ordinary = ordinary_countermodel()
    if ordinary["radius_squared_at_zero"] != "61/148":
        raise RuntimeError("ordinary radius audit failed")
    if ordinary["aligned_Q"] != "-53/8":
        raise RuntimeError("ordinary aligned sign audit failed")
    if ordinary["opposed_Q"] != "7/40":
        raise RuntimeError("ordinary opposed sign audit failed")
    zero_fiber = zero_fiber_countermodel()
    if zero_fiber["abs_W_A_squared"] != "637/1280":
        raise RuntimeError("zero-fiber magnitude audit failed")
    return {
        "relative_carrier_rotation_lower": format(
            carrier_rotation_lower, ".17e"
        ),
        "terminal_relative_at_L50": format(terminal_ratio, ".17e"),
        "ordinary_radius_squared": ordinary["radius_squared_at_zero"],
        "ordinary_aligned_Q": ordinary["aligned_Q"],
        "ordinary_opposed_Q": ordinary["opposed_Q"],
        "zero_fiber_WA_squared": zero_fiber["abs_W_A_squared"],
    }


def build_exact() -> dict:
    symbolic = symbolic_audit()
    numeric = numeric_audit()
    return {
        "domain": (
            "L=log(x/(4*pi))>=50, 0<=tL<=25, "
            "a^2=x/(4*pi)+t/16, N=floor(a); "
            "the all-fiber q>=1 target additionally has q=2tL^2>=1"
        ),
        "pi_provenance": symbolic["pi_provenance"],
        "coordinates": symbolic["coordinates"],
        "bulk_pairwise": symbolic["bulk_pairwise"],
        "endpoint_pairwise": symbolic["endpoint_pairwise"],
        "abel_prefix": symbolic["abel_prefix"],
        "endpoint_complete_prefix": symbolic["endpoint_complete_prefix"],
        "prefix_currents": symbolic["prefix_currents"],
        "contact_scalar": symbolic["contact_scalar"],
        "terminal_band_bound": (
            "On |mathsf_X|<=delta_L, "
            "0<=c=t*D_x/4<exp(-L)/4 and "
            "0<=u_N=log(a/N)<2*exp(-L/2), so "
            "|c*u_N*mathsf_X|"
            "<25000*exp(-11L/4)/|f_1|"
            "<exp(-3L/2)/(4L)*A_L"
            "<2e-35*A_L for L>=50."
        ),
        "all_fiber_sufficient": (
            "Let epsilon_term=25000*exp(-11L/4)/|f_1|. "
            "On |mathsf_X|<=delta_L, "
            "|mathcal_C_N|>A_L+epsilon_term implies "
            "|mathsf_A|>A_L. At mathsf_X=0, "
            "sign(mathsf_A)=sign(mathcal_C_N) exactly, even when W_0=0."
        ),
        "ordinary_kernel_target": (
            "For R=|W_0|>=tau_L>delta_L, "
            "P+i*Q=W_A*conj(W_0) and "
            "W_A/W_0=(P+i*Q)/R^2. The previous ordinary inequality is "
            "equivalent to "
            "|Q|*sqrt(R^2-delta_L^2)-|P|*delta_L>A_L*R^2. "
            "At a real crossing Q=-mathsf_Y*mathcal_C_N; the kernel loses "
            "the zero fiber only through multiplication by mathsf_Y."
        ),
        "exceptional_derivative": symbolic["exceptional_derivative"],
        "exceptional_error": (
            "On |W_0|<tau_L, "
            "|partial_x mathsf_X-mathcal_C_N|"
            "<epsilon_term+2*exp(-L)*tau_L"
            "+2e-7*exp(-5L/4). At W_0=0 only the final D term remains."
        ),
        "ordinary_countermodel": ordinary_countermodel(),
        "zero_fiber_countermodel": zero_fiber_countermodel(),
        "prefix_sector_guard": (
            "The positive h_k make a weighted prefix-sector or aggregate "
            "oscillatory estimate a sufficient Xi input, but termwise sign "
            "is not inherited from r_(n+1)<r_n, u_(n+1)<u_n, or the common "
            "coefficient direction s_*'. The exact countermodels forbid "
            "such a generic promotion."
        ),
        "relative_carrier_rotation": (
            "For chi_n=zeta_n*conj(zeta_1), |chi_n|=1, "
            "vartheta_(n,x)=Im(chi_(n,x)*conj(chi_n))="
            "(-b)*log(n)+Im[d_(n,x)/(1+d_n)"
            "-d_(1,x)/(1+d_1)]. Since -b>=1/2, "
            "|d_n|<2189/x<1/2, and |d_(n,x)|<4223/x^2, "
            "vartheta_(n,x)>log(2)/2-16892/x^2>1/3 for n>=2. "
            "At fixed t a complete cutoff cell has "
            "Delta x=4*pi*(2N+1), so chi_2 advances by more than "
            "Delta x/3>2*pi and takes every unit-circle value."
        ),
        "sector_rotation_guard": (
            "The actual Xi relative carriers are monotone rotating, not "
            "confined to one fixed sector on a cutoff cell. Therefore a "
            "termwise cosine/sine sign or fixed common half-plane cannot "
            "be the uniform proof. Relative rotation may instead be used "
            "inside an aggregate oscillatory estimate or an explicit "
            "phase-cell decomposition."
        ),
        "live_target": (
            "With delta_L=50000*exp(-5L/4)/|f_1|, "
            "A_L=(100000L+1)*exp(-5L/4)/|f_1|, and "
            "epsilon_term=25000*exp(-11L/4)/|f_1|, prove on q>=1 that "
            "|mathsf_X|<=delta_L implies "
            "|mathcal_C_N|>A_L+epsilon_term, and prove that the number of "
            "crossings mathsf_X=0 with mathcal_C_N>0 is below one successor "
            "turn after the prescribed boundary composition. This single "
            "division-free prefix theorem includes W_0=0."
        ),
        "q_lt_1_target": (
            "The q=2tL^2<1 boundary remains a separate "
            "multiplicity-compatible parabolic/Hermite theorem. The "
            "prefix identity remains exact there, but no uniform positive "
            "slope is imposed at t=0."
        ),
        "chart_scope": (
            "The prefixes F_k and complex kernels are canonical-chart "
            "objects. Prove the prefix theorem in one prescribed chart and "
            "transfer only mathsf_X and mathsf_A through the certified "
            "adjacent real-projection bounds."
        ),
        "numeric_audit": numeric,
    }


def build_rows(exact: dict) -> list[KernelRow]:
    return [
        KernelRow(
            "nfockapr_00_pi_provenance",
            "definition_provenance",
            "certified",
            "Every occurrence of pi in the saddle, cutoff cell, and phase turn is traced to a fixed analytic normalization.",
            exact["pi_provenance"],
            "No circle fitted to the prefix polygon and no polygonal approximation defines pi here.",
        ),
        KernelRow(
            "nfockapr_01_coordinates",
            "exact_normal_form",
            "ready_to_apply",
            "The anchored value and slope are an endpoint plus a positive-amplitude unit-carrier sum and its logarithmic moment.",
            exact["coordinates"],
            "This introduces no phase branch or sign claim.",
        ),
        KernelRow(
            "nfockapr_02_bulk_pairwise",
            "exact_kernel",
            "ready_to_apply",
            "The bulk radial and tangential currents have explicit diagonal and unordered-pair kernels.",
            exact["bulk_pairwise"],
            "The cosine and sine kernels are indefinite.",
        ),
        KernelRow(
            "nfockapr_03_endpoint_pairwise",
            "exact_kernel",
            "ready_to_apply",
            "The endpoint-complete current adds one endpoint diagonal and one cross term per carrier.",
            exact["endpoint_pairwise"],
            "No endpoint interaction is discarded.",
        ),
        KernelRow(
            "nfockapr_04_abel_prefix",
            "exact_reindexing",
            "ready_to_apply",
            "Finite Abel summation collapses the logarithmic moment to positive gap weights on prefixes.",
            exact["abel_prefix"],
            "The positive weights do not make the complex prefix projections one-signed.",
        ),
        KernelRow(
            "nfockapr_05_endpoint_complete_prefix",
            "exact_normal_form",
            "ready_to_apply",
            "The complete centered slope has one terminal block and an O(N) Abel-prefix expansion.",
            exact["endpoint_complete_prefix"],
            "The terminal block retains the exact endpoint phase.",
        ),
        KernelRow(
            "nfockapr_06_prefix_currents",
            "exact_kernel",
            "ready_to_apply",
            "Both ordinary Hermitian kernels have endpoint-complete O(N) prefix formulas.",
            exact["prefix_currents"],
            "These are the same radial and tangential ratio numerators, not a new positivity theorem.",
        ),
        KernelRow(
            "nfockapr_07_contact_scalar",
            "exact_reduction",
            "ready_to_apply",
            "One division-free prefix scalar equals the physical slope at every real crossing, including the zero fiber.",
            exact["contact_scalar"],
            "Unlike the Wronskian, this scalar does not vanish merely because W_0=0.",
        ),
        KernelRow(
            "nfockapr_08_terminal_band",
            "asymptotic_certificate",
            "certified",
            "The only contact-band difference between the prefix scalar and the physical slope is negligible.",
            exact["terminal_band_bound"],
            "Valid uniformly on L>=50 and 0<=tL<=25.",
            exact["numeric_audit"],
        ),
        KernelRow(
            "nfockapr_09_all_fiber_sufficient",
            "conditional_inequality",
            "ready_to_apply",
            "A single prefix-scalar lower bound closes the direct projection theorem on every complex fiber.",
            exact["all_fiber_sufficient"],
            "The Xi-specific lower bound is not proved here.",
        ),
        KernelRow(
            "nfockapr_10_ordinary_connection",
            "exact_reduction",
            "ready_to_apply",
            "The previous ordinary ratio target is exactly the division-free Hermitian-kernel inequality.",
            exact["ordinary_kernel_target"],
            "The tangential kernel remains the centered Wronskian away from the zero fiber.",
        ),
        KernelRow(
            "nfockapr_11_exceptional_derivative",
            "exact_reduction",
            "certified",
            "The same prefix scalar controls the exceptional real derivative with certified nuisance errors.",
            exact["exceptional_derivative"] + " " + exact["exceptional_error"],
            "This identity supplies no arithmetic lower bound.",
        ),
        KernelRow(
            "nfockapr_12_ordinary_countermodel",
            "countermodel",
            "guard_validated",
            "Positive ordered amplitudes and distances do not give a phase-current sign or magnitude floor.",
            json.dumps(exact["ordinary_countermodel"], sort_keys=True),
            "Exact generic carrier model, not an Xi counterexample.",
            exact["ordinary_countermodel"],
        ),
        KernelRow(
            "nfockapr_13_zero_fiber_countermodel",
            "countermodel",
            "guard_validated",
            "Positive ordered amplitudes and distances do not give directional real transversality on the zero fiber.",
            json.dumps(exact["zero_fiber_countermodel"], sort_keys=True),
            "The actual Xi absolute anchor and endpoint are absent and must be used by a successful theorem.",
            exact["zero_fiber_countermodel"],
        ),
        KernelRow(
            "nfockapr_14_prefix_sector_guard",
            "nonpromotion_guard",
            "guard_validated",
            "Abel positivity localizes the arithmetic input but does not make prefix projections positive.",
            exact["prefix_sector_guard"],
            "No termwise cone or common-half-plane theorem is promoted.",
        ),
        KernelRow(
            "nfockapr_15_relative_carrier_rotation",
            "asymptotic_certificate",
            "certified",
            "Every actual nonfirst Xi carrier makes a full monotone relative turn inside each complete cutoff cell.",
            (
                exact["relative_carrier_rotation"]
                + " "
                + exact["sector_rotation_guard"]
            ),
            "This is phase structure, not aggregate transversality.",
            exact["numeric_audit"],
        ),
        KernelRow(
            "nfockapr_16_q_ge_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q>=1 problem is one all-fiber anchored Abel-prefix lower bound plus its oriented crossing count.",
            exact["live_target"],
            "The Xi prefix theorem and signed count remain open.",
        ),
        KernelRow(
            "nfockapr_17_q_lt_1_target",
            "open_theorem_target",
            "not_ready_to_apply",
            "The q<1 multiplicity-compatible chart remains separate.",
            exact["q_lt_1_target"],
            "Endpoint simplicity is not assumed.",
        ),
        KernelRow(
            "nfockapr_18_nonpromotion",
            "nonpromotion_guard",
            "guard_validated",
            "The exact kernel and prefix reductions are kept separate from the missing Xi theorem and RH.",
            exact["chart_scope"],
            (
                "No prefix lower bound, signed successor count, q<1 "
                "theorem, finite phase closure, contact exclusion, "
                "Lambda<=0, PF-infinity, RH, or Clay-prize conclusion."
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
            "exact endpoint-complete carrier kernels, Abel-prefix "
            "transversality scalar, and nonpromotion guards; the "
            "Xi-specific prefix lower bound and signed count remain open"
        ),
        "proof_boundary": (
            "This artifact proves the endpoint-complete pairwise radial "
            "and tangential kernels, their finite Abel-prefix forms, one "
            "continuous all-fiber crossing scalar, the negligible terminal "
            "contact-band correction, the exceptional derivative identity, "
            "and two exact generic nonpromotion countermodels. It does not "
            "prove the q>=1 Xi prefix lower bound or signed successor count, "
            "the q<1 multiplicity-compatible theorem, finite phase cells, "
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
            "d_x": D_X_CONSTANT,
            "terminal_band_constant": TERMINAL_BAND_CONSTANT,
            "terminal_relative_at_L_min": TERMINAL_RELATIVE_AT_L_MIN,
        },
        "exact": exact,
        "rows": [asdict(row) for row in rows],
    }


def render_note(artifact: dict) -> str:
    exact = artifact["exact"]
    ordinary = exact["ordinary_countermodel"]
    zero_fiber = exact["zero_fiber_countermodel"]
    numeric = exact["numeric_audit"]
    return "\n".join(
        [
            "# Newman Centered Carrier-Kernel and Abel-Prefix Reduction",
            "",
            "Date: 2026-07-27",
            "",
            "Status: exact endpoint-complete kernels and one continuous",
            "all-fiber prefix scalar. The Xi-specific lower bound and",
            "signed crossing theorem remain open; this is not a proof of RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Pi Provenance",
            "",
            "```text",
            exact["pi_provenance"],
            "```",
            "",
            "The same universal constant is the circumference-to-diameter",
            "ratio of every Euclidean circle. No particular circle is",
            "selected here, and the complex prefix polygon is not used to",
            "construct or approximate it.",
            "",
            "## Anchored Coordinates",
            "",
            "```text",
            exact["coordinates"],
            "```",
            "",
            "## Pairwise Kernels",
            "",
            "The bulk radial and tangential numerators are",
            "",
            "```text",
            exact["bulk_pairwise"],
            "```",
            "",
            "The endpoint-complete form is",
            "",
            "```text",
            exact["endpoint_pairwise"],
            "```",
            "",
            "This is the ordinary Wronskian numerator in explicit carrier",
            "coordinates. Its value is not assumed to have one sign.",
            "",
            "## Abel Prefix Form",
            "",
            "Finite Abel summation gives",
            "",
            "```text",
            exact["abel_prefix"],
            exact["endpoint_complete_prefix"],
            "```",
            "",
            "Consequently both Hermitian kernels require only O(N) prefix",
            "projections:",
            "",
            "```text",
            exact["prefix_currents"],
            "```",
            "",
            "## Continuous Crossing Scalar",
            "",
            "Define the division-free prefix scalar",
            "",
            "```text",
            exact["contact_scalar"],
            "```",
            "",
            "The terminal real correction on the contact band obeys",
            "",
            "```text",
            exact["terminal_band_bound"],
            "```",
            "",
            "At `L=50` its ratio to the target is below",
            f"`{numeric['terminal_relative_at_L50']}`.",
            "",
            "Thus one sufficient all-fiber theorem is",
            "",
            "```text",
            exact["all_fiber_sufficient"],
            "```",
            "",
            "At a real crossing the equality is exact, so the same scalar",
            "also carries the orientation. The ordinary kernel relation is",
            "",
            "```text",
            exact["ordinary_kernel_target"],
            "```",
            "",
            "The exceptional derivative is",
            "",
            "```text",
            exact["exceptional_derivative"],
            exact["exceptional_error"],
            "```",
            "",
            "This is why `mathcal_C_N`, unlike the determinant, retains",
            "the `W_0=0` class.",
            "",
            "## Exact Route Guards",
            "",
            "The ordinary three-carrier model has",
            "",
            "```text",
            json.dumps(ordinary, indent=2, sort_keys=True),
            "```",
            "",
            "It reaches both current signs and an exact current zero with",
            "`|W_0|^2=61/148`, despite strictly decreasing positive",
            "amplitudes and distances.",
            "",
            "The zero-fiber model has",
            "",
            "```text",
            json.dumps(zero_fiber, indent=2, sort_keys=True),
            "```",
            "",
            "It has `W_0=0`, `W_A!=0`, and `Re(W_A)=0`. Therefore carrier",
            "positivity, amplitude decrease, and distance ordering cannot",
            "supply the physical directional slope. The actual Xi anchor",
            "and endpoint phase must be used.",
            "",
            "The prefix-sector guard is",
            "",
            "```text",
            exact["prefix_sector_guard"],
            "```",
            "",
            "The actual Xi relative carriers satisfy",
            "",
            "```text",
            exact["relative_carrier_rotation"],
            exact["sector_rotation_guard"],
            "```",
            "",
            "At `L=50` the certified relative-current lower bound is",
            f"`{numeric['relative_carrier_rotation_lower']}`.",
            "",
            "## Live Theorem",
            "",
            "```text",
            exact["live_target"],
            "```",
            "",
            "The q<1 and chart boundaries are",
            "",
            "```text",
            exact["q_lt_1_target"],
            exact["chart_scope"],
            "```",
            "",
            "This artifact proves no Xi prefix lower bound, signed",
            "successor count, q<1 chart, contact exclusion, `Lambda<=0`,",
            "PF-infinity, RH, or Clay-prize conclusion.",
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
        "built Newman centered carrier-kernel Abel reduction: "
        "19 rows, 3 exact kernel forms, 1 all-fiber prefix scalar, "
        "2 exact countermodels, 2 open Xi obligations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
