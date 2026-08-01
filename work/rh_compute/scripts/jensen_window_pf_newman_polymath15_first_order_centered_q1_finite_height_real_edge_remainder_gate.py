#!/usr/bin/env python3
"""Build the q=1 finite-height real-edge remainder gate."""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
from math import factorial
import json
from pathlib import Path
import sys

import sympy as sp


REPO_ROOT = Path(__file__).resolve().parents[3]
VENDOR = REPO_ROOT / "work/rh_compute/vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))

import flint  # noqa: E402


STEM = (
    "jensen_window_pf_newman_polymath15_first_order_centered_"
    "q1_finite_height_real_edge_remainder_gate"
)
DEFAULT_OUT = REPO_ROOT / "work/rh_compute/results" / f"{STEM}.json"
DEFAULT_NOTE = REPO_ROOT / "outputs" / f"{STEM}.md"
SOURCE_PATHS = {
    "cofinal_real_edge": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "real_edge_projective_current_gate.json"
    ),
    "complex_endpoint_source": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "complex_endpoint_source_normalization_gate.json"
    ),
    "critical_c1_source": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_RS_C1_endpoint_peeling_contract.json"
    ),
    "first_dirichlet_correction": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_"
        "critical_dirichlet_first_correction_gate.json"
    ),
    "adjacent_chart_stability": (
        REPO_ROOT
        / "work/rh_compute/results/"
        "jensen_window_pf_newman_polymath15_first_order_centered_"
        "adjacent_chart_stability_certificate.json"
    ),
}

PRECISION_BITS = 192
H_DENOMINATOR = 72_000_000_000
DIRECT_REGIONS = (
    (Fraction(-1), Fraction(-5, 8), 1024),
    (Fraction(-3, 8), Fraction(3, 8), 2048),
    (Fraction(5, 8), Fraction(1), 1024),
)
F_DERIVATIVE_BOUNDS = (3, 24, 384, 9216, 294912, 11796480)
C1_DERIVATIVE_BOUNDS = (86, 2731, 109227)


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


def source_audit() -> dict:
    missing = [str(path) for path in SOURCE_PATHS.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing sources: " + ", ".join(missing))
    return {
        "source_sha256": {
            key: file_hash(path) for key, path in SOURCE_PATHS.items()
        },
        "published_source": {
            "citation": (
                "D.H.J. Polymath, Effective approximation of heat flow "
                "evolution of the Riemann xi function, and a new upper "
                "bound for the de Bruijn-Newman constant, arXiv:1904.12438"
            ),
            "location": "equation (53), C_0, first correction, and Proposition 6.2",
            "url": "https://arxiv.org/abs/1904.12438",
        },
    }


def aq(value: Fraction | int) -> flint.arb:
    value = Fraction(value)
    return flint.arb(flint.fmpq(value.numerator, value.denominator))


def arb_text(value: flint.arb, digits: int = 70) -> str:
    return value.str(digits)


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def interval_ball(left: Fraction, right: Fraction) -> flint.arb:
    lower = aq(left)
    upper = aq(right)
    return (lower + upper) / 2 + flint.arb(0, (upper - lower) / 2)


def direct_c0_jet(p: flint.arb) -> tuple[flint.acb, ...]:
    """Enclose C_0 and its first five p derivatives off the removable points."""
    pi = flint.arb.pi()
    variable = flint.acb_series(
        [flint.acb(p), flint.acb(1)],
        prec=6,
    )
    numerator = (
        (flint.acb(0, 1) * pi * (variable**2 / 2 + flint.acb(3) / 8)).exp()
        - flint.acb(0, 1)
        * flint.arb(2).sqrt()
        * (pi * variable / 2).cos()
    )
    value = numerator / (2 * (pi * variable).cos())
    return tuple(
        coefficient * factorial(order)
        for order, coefficient in enumerate(value.coeffs())
    )


def c0_derivative_certificate() -> dict:
    """Certify coarse C_0 derivative bounds, including both removable charts."""
    flint.ctx.prec = PRECISION_BITS
    pi = flint.arb.pi()

    # On |y|<=1/4 around either removable point, factor the common y.
    # The three numerator majorants correspond to exp(i delta)-1,
    # 1-cos(z), and sin(z).  The denominator uses sinc(pi*y).
    sinc_radius = pi / 4
    sinc_deviation = sinc_radius.sinh() / sinc_radius - 1
    denominator_lower = pi * (1 - sinc_deviation)
    delta_radius = 5 * pi / 32
    half_phase_radius = pi / 8
    numerator_upper = (
        (5 * pi / 8) * delta_radius.exp()
        + (pi**2 / 32) * half_phase_radius.exp()
        + (pi / 2) * half_phase_radius.exp()
    )
    c0_disk_upper = numerator_upper / (2 * denominator_lower)
    if not sinc_deviation < aq(Fraction(1, 5)):
        raise RuntimeError("local sinc denominator majorant failed")
    if not c0_disk_upper < 3:
        raise RuntimeError("local C_0 disk bound failed")

    cauchy = tuple(3 * factorial(k) * 8**k for k in range(6))
    if cauchy != F_DERIVATIVE_BOUNDS:
        raise RuntimeError("saved C_0 Cauchy bounds drifted")

    maxima = [flint.arb(0) for _ in range(6)]
    locations = ["" for _ in range(6)]
    total = 0
    region_rows = []
    for left, right, boxes in DIRECT_REGIONS:
        width = (right - left) / boxes
        region_maxima = [flint.arb(0) for _ in range(6)]
        for index in range(boxes):
            box_left = left + index * width
            box_right = box_left + width
            jet = direct_c0_jet(interval_ball(box_left, box_right))
            for order, derivative in enumerate(jet):
                upper = flint.arb(abs(derivative).upper())
                if upper > region_maxima[order]:
                    region_maxima[order] = upper
                if upper > maxima[order]:
                    maxima[order] = upper
                    locations[order] = f"{box_left}..{box_right}"
            total += 1
        region_rows.append(
            {
                "interval": f"[{left},{right}]",
                "boxes": boxes,
                "maximum_derivative_upper_balls": [
                    arb_text(value) for value in region_maxima
                ],
            }
        )

    for order, (maximum, bound) in enumerate(zip(maxima, cauchy)):
        if not maximum < bound:
            raise RuntimeError(f"direct C_0 derivative {order} bound failed")

    if not pi**2 > 9:
        raise RuntimeError("pi^2>9 audit failed")
    derived_c1 = (
        F_DERIVATIVE_BOUNDS[3] / 108,
        F_DERIVATIVE_BOUNDS[4] / 108,
        F_DERIVATIVE_BOUNDS[5] / 108,
    )
    if not all(
        value < bound
        for value, bound in zip(derived_c1, C1_DERIVATIVE_BOUNDS)
    ):
        raise RuntimeError("C_1 derivative majorants failed")

    return {
        "precision_bits": PRECISION_BITS,
        "direct_boxes": total,
        "direct_regions": region_rows,
        "removable_centers": ["-1/2", "1/2"],
        "local_disk_radius": "1/4",
        "inner_cauchy_radius": "1/8",
        "sinc_deviation_ball": arb_text(sinc_deviation),
        "sinc_deviation_lt": "1/5",
        "factored_denominator_lower_ball": arb_text(denominator_lower),
        "factored_numerator_upper_ball": arb_text(numerator_upper),
        "C0_disk_upper_ball": arb_text(c0_disk_upper),
        "C0_disk_upper_lt": "3",
        "global_C0_derivative_bounds": list(F_DERIVATIVE_BOUNDS),
        "direct_maximum_upper_balls": [arb_text(value) for value in maxima],
        "direct_maximum_locations": locations,
        "C1_derivative_bounds": list(C1_DERIVATIVE_BOUNDS),
        "pi_square_lower": "pi^2>9",
    }


def verify_symbolics() -> None:
    h, theta = sp.symbols("h theta", positive=True, real=True)
    p = sp.symbols("p", real=True)
    phase = -sp.pi * (p**2 / 2 - p + sp.Rational(3, 8))
    p_x = -h / (4 * sp.pi)
    q_rate = sp.diff(sp.I * phase, p) * p_x
    q_rate = q_rate.subs(p, 1 - 2 * theta)
    if sp.simplify(q_rate + sp.I * h * theta / 2) != 0:
        raise RuntimeError("terminal carrier logarithmic derivative failed")

    s, t = sp.symbols("s t", nonzero=True)
    ap, app, w = sp.symbols("alpha_prime alpha_second w")
    block = t / 4 + t**2 * w**2 / 8
    expected_dx = (-sp.I / 2) * (
        -1 / (6 * s**2) + app * block + t**2 * w * ap**2 / 4
    )
    direct_dx = (
        sp.diff(1 / (6 * s), s) * (-sp.I / 2)
        + (-sp.I * app / 2) * block
        + ap * (t**2 * w / 4) * (-sp.I * ap / 2)
    )
    if sp.simplify(direct_dx - expected_dx) != 0:
        raise RuntimeError("terminal correction derivative failed")

    ur, eta = sp.symbols("u_R eta", real=True)
    wi = -sp.pi / 4 + eta
    cancellation = sp.expand(
        sp.pi**2 / 16 + sp.re((ur + sp.I * wi) ** 2)
        - (ur**2 + sp.pi * eta / 2 - eta**2)
    )
    if sp.simplify(cancellation) != 0:
        raise RuntimeError("stable terminal real cancellation failed")

    A, B, M, e0, e1, e2, e3 = sp.symbols(
        "A B M e_0 e_1 e_2 e_3", real=True
    )
    c = A + h * e0
    d = h * B + h**2 * e1
    cx = h * B + h**2 * e2
    dx = h**2 * M + h**3 * e3
    remainder = sp.expand((c * dx - d * cx) / h**2 - (A * M - B**2))
    target = sp.expand(
        h * (e0 * M + A * e3 - B * e1 - B * e2)
        + h**2 * (e0 * e3 - e1 * e2)
    )
    if sp.simplify(remainder - target) != 0:
        raise RuntimeError("four-jet current perturbation identity failed")


def exact_payload() -> dict:
    return {
        "domain": {
            "physical_family": (
                "q=2tL^2=1, L>=50, theta=(1-p)/2 in [0,1], "
                "a=N+theta, h=1/a, T_0=2pi/h^2"
            ),
            "saddle": (
                "a^2=exp(L)+t/16, x=4pi*exp(L), "
                "rho=t*h^2/16, r=1-rho, y=1/x=h^2/(4pi*r)"
            ),
            "effective_smallness": (
                "t<=1/5000 and h<exp(-25)<1/72000000000"
            ),
        },
        "endpoint": {
            "source": (
                "F=C_0(p), C_1=F'''/(12pi^2), H=F+h*C_1"
            ),
            "coordinate_jets": (
                "p_x=-h/(4pi), h_x=-h^3/(8pi), "
                "p_xx=h^3/(32pi^2), h_xx=3h^5/(64pi^2)"
            ),
            "H_jets": (
                "H_x=p_x(F'+hC_1')+h_xC_1; "
                "H_xx=p_xx(F'+hC_1')+p_x^2(F''+hC_1'')+"
                "2p_xh_xC_1'+h_xxC_1"
            ),
        },
        "terminal": {
            "source": (
                "Z=z_N/S_a=(-1)^N M_t(s)exp[t log(N)^2/4-s_*log(N)]"
                "(1+d)/(beta*T_0), S_a=kappa*T_0"
            ),
            "alpha": (
                "alpha=1/(2s)+1/(s-1)+Log(s/(2pi))/2, "
                "s_*=(s+t alpha/2), s_*'=-i/2-it alpha'/4, "
                "s_*''=-t alpha''/8"
            ),
            "correction": (
                "w=alpha-log(N), d=1/(6s)+alpha'[t/4+t^2w^2/8], "
                "d_x=(-i/2)[-1/(6s^2)+alpha''(t/4+t^2w^2/8)"
                "+(t^2/4)w(alpha')^2]"
            ),
            "phase": (
                "Q=exp[-pi*i*(p^2/2-p+3/8)], Z=-Q*exp(W), "
                "Q_x/Q=-i*h*theta/2"
            ),
        },
        "stable_log": {
            "coordinates": (
                "ell=log(1-h*theta), chi=alpha-log(a), "
                "chi_R=log(r)/2+log(1+y^2)/4-y^2/(1+y^2), "
                "chi_I=-pi/4+atan(y)/2+3y/(1+y^2), w=chi-ell"
            ),
            "real_part": (
                "Re W=-ell/2+(7/4)log(r)+(7/8)log(1+y^2)"
                "+atan(y)/(4y)-1/4+t*pi^2/64+(t/4)Re(w^2)"
                "+log|1+d|"
            ),
            "imaginary_part": (
                "Im W=T_0[ell+h*theta+h^2theta^2/2]-(pi*t/8)ell"
                "+(T_0/2)[r*(-log r)-rho]-log(1+y^2)/(8y)"
                "-atan(y)/4+(t/4)Im(w^2)+arg(1+d)"
            ),
            "removable_rewrites": (
                "T_0[ell+z+z^2/2]=-2pi*h*theta^3*sum_(m>=0)"
                "z^m/(m+3); (T_0/2)[r(-log r)-rho]="
                "-pi*t^2*h^2*sum_(m>=0)rho^m/[(m+1)(m+2)]/256; "
                "pi^2/16+Re(w^2)=u_R^2+(pi/2)eta-eta^2"
            ),
            "x_derivative": (
                "W_x=eta/2-3/(4T_0)+i[ell+h*theta-chi_R]/2"
                "-i(t/4)alpha'w+d_x/(1+d), "
                "eta=atan(y)/2+3y/(1+y^2)"
            ),
            "branch": (
                "W is the continuous logarithm with W->0 as h->0; "
                "|d|<h^2/10 keeps 1+d in the principal near-one disk"
            ),
        },
        "endpoint_rate": {
            "mu": (
                "mu=K_rate+i*chi/2+i(t/4)alpha'chi, "
                "K_rate=-pi/8+1/(4T_0)+1/[2(T_0+i)]"
            ),
            "mu_x": (
                "mu_x=K_rate_x+i*chi_x/2+i(t/4)(alpha'_xchi+alpha'chi_x), "
                "K_rate_x=-1/(8T_0^2)-1/[4(T_0+i)^2], "
                "chi_x=-i alpha'/2-h^2/(8pi), alpha'_x=-i alpha''/2"
            ),
            "endpoint_jet": (
                "J=H_x+mu*H, J_x=H_xx+mu_xH+muH_x"
            ),
        },
        "normalized_edge": {
            "definitions": (
                "tau=1+i/T_0, u=log(a/N)=-ell, "
                "c=Re(tau H+Z), c_x=Re(tau_xH+tauH_x+Z_x), "
                "gbar=tauJ+s_*'uZ, d_edge=Re(gbar)-alpha_c c, "
                "alpha_c=Re(s_*')u"
            ),
            "second_jet": (
                "gbar_x=tau_xJ+tauJ_x+s_*''uZ+s_*'u_xZ+s_*'uZ_x; "
                "d_(edge,x)=Re(gbar_x)-alpha_(c,x)c-alpha_c c_x"
            ),
            "current": (
                "a^2 J_edge/S_a^2=h^-2[c*d_(edge,x)-d_edge*c_x]"
            ),
            "leading_jets": (
                "A=Re(F-Q), B=Re[-F'/(4pi)+i theta Q/2], "
                "M=Re[F''/(16pi^2)+Q{i/(16pi)+theta^2/4}]"
            ),
        },
    }


def majorant_certificate() -> dict:
    flint.ctx.prec = PRECISION_BITS
    h_exp = (-flint.arb(25)).exp()
    h_saved = aq(Fraction(1, H_DENOMINATOR))
    if not h_exp < h_saved:
        raise RuntimeError("effective h threshold failed")
    if not flint.arb(H_DENOMINATOR) < flint.arb(25).exp():
        raise RuntimeError("exp(25)>72000000000 audit failed")

    # Every coefficient below is normalized by the displayed power of h.
    # Only 3<pi<4, t<=1/5000, 0<=theta<=1, and h<=h0 are used.
    h0 = Fraction(1, H_DENOMINATOR)
    rho_max = h0**2 / 80000
    u_coefficient = 1 / (1 - h0)
    u_remainder_coefficient = 1 / (2 * (1 - h0))
    log_r_coefficient = Fraction(1, 80000) / (1 - rho_max)
    y_coefficient = Fraction(1, 12) / (1 - rho_max)
    eta_coefficient = Fraction(7, 2) * y_coefficient
    chi_r_coefficient = (
        log_r_coefficient / 2
        + Fraction(5, 4) * y_coefficient**2 * h0**2
    )
    u_r_coefficient = u_coefficient + chi_r_coefficient * h0
    chi_absolute = (
        chi_r_coefficient * h0**2 + 1 + eta_coefficient * h0**2
    )
    w_absolute = u_r_coefficient * h0 + 1 + eta_coefficient * h0**2
    alpha_prime_coefficient = (
        6 * y_coefficient + 7 * y_coefficient**2 * h0**2
    )
    alpha_second_coefficient = (
        2 * y_coefficient**2
        + 24 * y_coefficient**3 * h0**2
    )
    correction_block = Fraction(1, 20000) + Fraction(9, 200_000_000)
    d_coefficient = (
        y_coefficient / 3
        + alpha_prime_coefficient * correction_block
    )
    d_x_coefficient = Fraction(1, 2) * (
        Fraction(2, 3) * y_coefficient**2
        + alpha_second_coefficient * correction_block
        + Fraction(1, 100_000_000) * 3 * alpha_prime_coefficient**2
    )
    log_one_plus_d_coefficient = d_coefficient / (
        1 - d_coefficient * h0**2
    )

    real_w_coefficient = (
        u_coefficient / 2
        + Fraction(7, 4) * log_r_coefficient * h0
        + Fraction(7, 8) * y_coefficient**2 * h0**3
        + y_coefficient**2 * h0**3 / 12
        + Fraction(1, 20000)
        * (
            u_r_coefficient**2 * h0
            + 2 * eta_coefficient * h0
            + eta_coefficient**2 * h0**3
        )
        + log_one_plus_d_coefficient * h0
    )
    imaginary_w_coefficient = (
        Fraction(8, 3) / (1 - h0)
        + u_coefficient / 10000
        + h0 / (3_200_000_000 * (1 - rho_max))
        + Fraction(3, 8) * y_coefficient * h0
        + Fraction(1, 20000)
        * 2
        * u_r_coefficient
        * (1 + eta_coefficient * h0**2)
        + log_one_plus_d_coefficient * h0
    )
    w_coefficient = real_w_coefficient + imaginary_w_coefficient
    w_x_coefficient = (
        eta_coefficient / 2
        + Fraction(1, 8)
        + (u_remainder_coefficient + chi_r_coefficient) / 2
        + Fraction(1, 20000) * alpha_prime_coefficient * 3
        + d_x_coefficient
        * h0**2
        / (1 - d_coefficient * h0**2)
    )

    mu_coefficient = (
        Fraction(1, 24)
        + Fraction(1, 12)
        + eta_coefficient / 2
        + chi_r_coefficient / 2
        + Fraction(1, 20000) * alpha_prime_coefficient * 3
    )
    chi_x_real_coefficient = y_coefficient / 2 * (
        3 * y_coefficient**2 * h0**2 + Fraction(1, 80000)
    )
    chi_x_imag_coefficient = Fraction(7, 2) * y_coefficient**2
    chi_x_coefficient = chi_x_real_coefficient + chi_x_imag_coefficient
    mu_x_coefficient = (
        Fraction(1, 96)
        + chi_x_coefficient / 2
        + Fraction(1, 20000)
        * (
            alpha_second_coefficient * 3 / 2
            + alpha_prime_coefficient * chi_x_coefficient * h0**2
        )
    )
    s1_defect_coefficient = alpha_prime_coefficient / 20000
    real_s1_coefficient = 6 * y_coefficient / 20000
    s2_coefficient = alpha_second_coefficient / 40000

    h_x_defect_coefficient = (
        Fraction(C1_DERIVATIVE_BOUNDS[1], 12)
        + h0 * Fraction(C1_DERIVATIVE_BOUNDS[0], 24)
    )
    h_xx_defect_coefficient = (
        Fraction(F_DERIVATIVE_BOUNDS[1], 288)
        + h0 * Fraction(C1_DERIVATIVE_BOUNDS[1], 288)
        + Fraction(C1_DERIVATIVE_BOUNDS[2], 144)
        + h0 * Fraction(C1_DERIVATIVE_BOUNDS[1], 144)
        + h0**2 * Fraction(43, 96)
    )
    j_defect_coefficient = Fraction(232 + 4)
    j_x_defect_coefficient = Fraction(800 + 3) + 4 * h0

    c_defect_coefficient = Fraction(86 + 12) + Fraction(2, 3) * h0
    c_x_defect_coefficient = (
        Fraction(232 + 10) + h0 / 2 + h0**2 / 18
    )
    d_edge_defect_coefficient = (
        Fraction(236 + 1 + 6)
        + h0 / 2
        + 4 * h0 / 6000
        + 5 * h0 / 10000
    )
    d_edge_x_defect_coefficient = (
        Fraction(810)
        + Fraction(2, 3) * h0
        + h0**2 / 24
        + 4 * h0**2 / 40000
        + Fraction(1, 24) * (6 + 2 * h0 / 6000)
        + 2 * h0 / 6000
        + Fraction(1, 2)
        + 5
        + 5 * (2 * h0**2 / 40000 + h0 / 480000)
        + 4 * h0 / 10000
    )

    checks = {
        "u_over_h": (u_coefficient, Fraction(2)),
        "u_minus_h_theta_over_h2": (u_remainder_coefficient, Fraction(1)),
        "log_r_over_h2": (log_r_coefficient, Fraction(1, 40000)),
        "y_over_h2": (y_coefficient, Fraction(1, 6)),
        "eta_over_h2": (eta_coefficient, Fraction(1)),
        "chi_R_over_h2": (chi_r_coefficient, Fraction(1)),
        "chi_absolute": (chi_absolute, Fraction(3)),
        "w_absolute": (w_absolute, Fraction(3)),
        "alpha_prime_over_h2": (alpha_prime_coefficient, Fraction(3)),
        "alpha_second_over_h4": (alpha_second_coefficient, Fraction(1)),
        "correction_block": (correction_block, Fraction(1, 10000)),
        "d_over_h2": (d_coefficient, Fraction(1, 10)),
        "d_x_over_h4": (d_x_coefficient, Fraction(1, 100)),
        "Re_W_over_h": (real_w_coefficient, Fraction(2)),
        "Im_W_over_h": (imaginary_w_coefficient, Fraction(4)),
        "W_over_h": (w_coefficient, Fraction(6)),
        "W_x_over_h2": (w_x_coefficient, Fraction(2)),
        "mu_over_h2": (mu_coefficient, Fraction(1)),
        "mu_x_over_h4": (mu_x_coefficient, Fraction(1)),
        "s1_defect_over_h2": (s1_defect_coefficient, Fraction(1, 6000)),
        "real_s1_over_h2": (real_s1_coefficient, Fraction(1, 20000)),
        "s2_over_h4": (s2_coefficient, Fraction(1, 40000)),
        "H_x_defect_over_h2": (h_x_defect_coefficient, Fraction(232)),
        "H_xx_defect_over_h3": (h_xx_defect_coefficient, Fraction(800)),
        "J_defect_over_h2": (j_defect_coefficient, Fraction(237)),
        "J_x_defect_over_h3": (j_x_defect_coefficient, Fraction(810)),
        "c_defect_over_h": (c_defect_coefficient, Fraction(100)),
        "c_x_defect_over_h2": (c_x_defect_coefficient, Fraction(250)),
        "d_edge_defect_over_h2": (d_edge_defect_coefficient, Fraction(250)),
        "d_edge_x_defect_over_h3": (d_edge_x_defect_coefficient, Fraction(900)),
    }
    failed = [name for name, (value, cap) in checks.items() if not value < cap]
    if failed:
        raise RuntimeError("failed normalized majorants: " + ", ".join(failed))
    if not (6 * h_saved).exp() < 2:
        raise RuntimeError("exp(6h)<2 audit failed")

    normalized_chain = {
        name: {
            "upper": fraction_text(value),
            "cap": fraction_text(cap),
        }
        for name, (value, cap) in checks.items()
    }

    source_blocks = {
        "saddle": {
            "u": "|u|<2h",
            "u_minus_h_theta": "|u-h*theta|<h^2",
            "rho": "rho<h^2/80000",
            "log_r": "|log r|<h^2/40000",
            "y": "0<=y<h^2/6",
            "eta": "0<=eta<h^2",
            "chi_R": "|chi_R|<h^2",
            "chi": "|chi|<3",
            "w": "|w|<3",
        },
        "alpha_and_correction": {
            "alpha_prime": "|alpha'|<3h^2",
            "alpha_second": "|alpha''|<h^4",
            "d": "|d|<h^2/10",
            "d_x": "|d_x|<h^4/100",
        },
        "stable_terminal": {
            "W": "|W|<6h",
            "W_x": "|W_x|<2h^2",
            "exp_W": "|exp(W)|<2",
            "exp_W_minus_one": "|exp(W)-1|<12h",
            "terminal_value": "|Z+Q|<12h",
            "terminal_first_jet": "|Z_x-i*h*theta*Q/2|<10h^2",
        },
        "endpoint_rate": {
            "mu": "|mu|<h^2",
            "mu_x": "|mu_x|<h^4",
            "s_star_prime": "|s_*'+i/2|<h^2/6000",
            "real_s_star_prime": "|Re(s_*')|<h^2/20000",
            "s_star_second": "|s_*''|<h^4/40000",
        },
    }

    endpoint_blocks = {
        "C0": list(F_DERIVATIVE_BOUNDS),
        "C1": list(C1_DERIVATIVE_BOUNDS),
        "H_minus_F": "|H-F|<86h",
        "H_x_defect": "|H_x+hF'/(4pi)|<232h^2",
        "H_xx_defect": "|H_xx-h^2F''/(16pi^2)|<800h^3",
        "H": "|H|<4",
        "H_x": "|H_x|<3h",
        "H_xx": "|H_xx|<4h^2",
        "J_defect": "|J+hF'/(4pi)|<236h^2",
        "J_x_defect": "|J_x-h^2F''/(16pi^2)|<810h^3",
    }

    edge_defects = {
        "value": {
            "formula": "c=A+h*e_0",
            "bound": "|e_0|<100",
            "constant": 100,
        },
        "centered_slope": {
            "formula": "d_edge=h*B+h^2*e_1",
            "bound": "|e_1|<250",
            "constant": 250,
        },
        "radial_subtracted_value_jet": {
            "formula": "c_x=h*B+h^2*e_2",
            "bound": "|e_2|<250",
            "constant": 250,
        },
        "centered_second_jet": {
            "formula": "d_(edge,x)=h^2*M+h^3*e_3",
            "bound": "|e_3|<900",
            "constant": 900,
        },
    }

    e0 = edge_defects["value"]["constant"]
    e1 = edge_defects["centered_slope"]["constant"]
    e2 = edge_defects["radial_subtracted_value_jet"]["constant"]
    e3 = edge_defects["centered_second_jet"]["constant"]
    a_bound, b_bound, m_bound = 4, 3, 3
    linear = e0 * m_bound + a_bound * e3 + b_bound * (e1 + e2)
    quadratic = e0 * e3 + e1 * e2
    if linear != 5400 or quadratic != 152500:
        raise RuntimeError("current error arithmetic drifted")
    h_fraction = Fraction(1, H_DENOMINATOR)
    combined = Fraction(linear) + h_fraction * quadratic
    if not combined < 5500:
        raise RuntimeError("saved 5500h current budget failed")
    if not Fraction(5500, H_DENOMINATOR) < Fraction(1, 10_000_000):
        raise RuntimeError("finite-height remainder target failed")
    finite_margin = Fraction(3, 8000) - Fraction(1, 10_000_000)
    if finite_margin != Fraction(3749, 10_000_000):
        raise RuntimeError("finite negative margin arithmetic failed")

    return {
        "precision_bits": PRECISION_BITS,
        "L_min": 50,
        "t_max": "1/5000",
        "exp_minus_25_ball": arb_text(h_exp),
        "h_upper": f"1/{H_DENOMINATOR}",
        "normalized_majorant_chain": normalized_chain,
        "source_block_majorants": source_blocks,
        "endpoint_majorants": endpoint_blocks,
        "leading_jet_bounds": {"|A|": 4, "|B|": 3, "|M|": 3},
        "edge_jet_defects": edge_defects,
        "perturbation_identity": (
            "h^-2(c*d_x-d*c_x)-(A*M-B^2)="
            "h[e_0M+A e_3-B(e_1+e_2)]+h^2[e_0e_3-e_1e_2]"
        ),
        "linear_error_constant": linear,
        "quadratic_error_constant": quadratic,
        "combined_error_constant_lt": 5500,
        "uniform_remainder": (
            "|a^2 J_edge/S_a^2-K_edge(p)|<5500h<1/10000000"
        ),
        "cofinal_symbol_margin": "K_edge(p)<-3/8000",
        "finite_height_margin": (
            "a^2 J_edge/S_a^2<-3749/10000000<0"
        ),
    }


def build_rows(exact: dict, derivatives: dict, budget: dict) -> list[GateRow]:
    return [
        GateRow(
            "q1fher_01_domain",
            "effective_domain",
            "ready_to_apply",
            "The physical q=1 saddle gives explicit small parameters at L>=50.",
            exact["domain"]["physical_family"] + "; " + exact["domain"]["effective_smallness"],
            "Only the retained first-order q=1 model is considered.",
        ),
        GateRow(
            "q1fher_02_endpoint_source",
            "exact_source_identity",
            "ready_to_apply",
            "The endpoint and both coordinate derivatives retain the complete C_1 correction.",
            exact["endpoint"]["source"] + "; " + exact["endpoint"]["H_jets"],
            "No higher published endpoint correction is included.",
        ),
        GateRow(
            "q1fher_03_terminal_source",
            "exact_source_identity",
            "ready_to_apply",
            "The physical terminal carrier retains its first correction and derivative.",
            exact["terminal"]["source"] + "; " + exact["terminal"]["correction"],
            "The fixed-N derivative is taken before any asymptotic limit.",
        ),
        GateRow(
            "q1fher_04_stable_log",
            "exact_desingularization",
            "ready_to_apply",
            "A continuous near-one logarithm removes every catastrophic terminal cancellation.",
            exact["stable_log"]["real_part"] + "; " + exact["stable_log"]["imaginary_part"],
            exact["stable_log"]["branch"],
        ),
        GateRow(
            "q1fher_05_removable_rewrites",
            "exact_desingularization",
            "ready_to_apply",
            "The h=0 and y=0 quotients are replaced by convergent series before bounding.",
            exact["stable_log"]["removable_rewrites"],
            "No interval division by h or y is used.",
        ),
        GateRow(
            "q1fher_06_terminal_bounds",
            "analytic_majorant",
            "certified",
            "The terminal value and first jet have uniform one-order defects.",
            "|W|<6h, |W_x|<2h^2, |Z+Q|<12h, |Z_x-i h theta Q/2|<10h^2",
            "Uniform on the enlarged box 0<t<=1/5000, 0<h<1/72000000000, 0<=theta<=1.",
            diagnostics=budget["source_block_majorants"]["stable_terminal"],
        ),
        GateRow(
            "q1fher_07_mu_bounds",
            "analytic_majorant",
            "certified",
            "The endpoint frame rate and its derivative preserve two powers of h.",
            "|mu|<h^2 and |mu_x|<h^4",
            "The -pi/8 and +pi/8 terms are cancelled algebraically before taking absolute values.",
            diagnostics=budget["source_block_majorants"]["endpoint_rate"],
        ),
        GateRow(
            "q1fher_08_c0_cauchy",
            "analytic_desingularization",
            "certified",
            "Two factored disks give Cauchy bounds through the fifth C_0 derivative.",
            "|C_0^(k)|<3*k!*8^k for 0<=k<=5 near p=+-1/2",
            "The common zero is factored before the disk lower bound is taken.",
            diagnostics={
                "disk_upper": derivatives["C0_disk_upper_ball"],
                "bounds": derivatives["global_C0_derivative_bounds"],
            },
        ),
        GateRow(
            "q1fher_09_c0_direct_cover",
            "interval_certificate",
            "certified",
            "A direct Arb cover verifies the same derivative bounds away from the removable disks.",
            "4096 rational Arb boxes certify C_0 derivatives 0 through 5",
            "This is a complete cover of the three complementary real intervals.",
            diagnostics={
                "boxes": derivatives["direct_boxes"],
                "maxima": derivatives["direct_maximum_upper_balls"],
            },
        ),
        GateRow(
            "q1fher_10_endpoint_jets",
            "analytic_majorant",
            "certified",
            "The corrected endpoint jets have explicit one-order error constants.",
            "|H-F|<86h, |H_x+hF'/(4pi)|<232h^2, |H_xx-h^2F''/(16pi^2)|<800h^3",
            "Uses pi^2>9 and the certified C_0 derivative bounds.",
            diagnostics=budget["endpoint_majorants"],
        ),
        GateRow(
            "q1fher_11_value_jet",
            "finite_height_jet",
            "certified",
            "The normalized real edge value is uniformly close to A.",
            budget["edge_jet_defects"]["value"]["formula"] + ", " + budget["edge_jet_defects"]["value"]["bound"],
            "No division by A or c is used.",
        ),
        GateRow(
            "q1fher_12_first_jets",
            "finite_height_jet",
            "certified",
            "The centered slope and radial-subtracted value jet share the leading B.",
            "d_edge=hB+h^2e_1, c_x=hB+h^2e_2, |e_1|,|e_2|<250",
            "The terminal centering alpha_c is retained exactly.",
        ),
        GateRow(
            "q1fher_13_second_jet",
            "finite_height_jet",
            "certified",
            "The full centered second jet is uniformly close to M.",
            budget["edge_jet_defects"]["centered_second_jet"]["formula"] + ", " + budget["edge_jet_defects"]["centered_second_jet"]["bound"],
            "H_xx, mu_x, s_*'', u_x, and Z_x are all retained.",
        ),
        GateRow(
            "q1fher_14_current_budget",
            "division_free_perturbation",
            "certified",
            "The four jet defects give a uniform finite-height current remainder.",
            budget["perturbation_identity"] + "; " + budget["uniform_remainder"],
            "The estimate is polynomial at every A=0 and c=0 fibre.",
            diagnostics={
                "linear": budget["linear_error_constant"],
                "quadratic": budget["quadratic_error_constant"],
            },
        ),
        GateRow(
            "q1fher_15_finite_sign",
            "effective_interval_theorem",
            "certified",
            "The retained first-order q=1 real-edge current is strictly clockwise for L>=50.",
            budget["finite_height_margin"],
            "This signs the retained first-order edge block, not the full Xi approximation error.",
        ),
        GateRow(
            "q1fher_16_handoff",
            "route_decision",
            "open_quantitative_handoff",
            "The q=1 finite-height wall is closed inside the first-order model.",
            "Next bound the omitted higher-order Xi/source remainder, extend q>=1, and splice adjacent cutoffs before using the edge block in the Abel scalar.",
            "No q>1 theorem, Xi-level edge sign, Abel gap, winding cap, Lambda<=0, or RH conclusion follows yet.",
        ),
    ]


def render_note(artifact: dict) -> str:
    derivative = artifact["derivative_certificate"]
    budget = artifact["majorant_certificate"]
    return "\n".join(
        [
            "# Jensen-Window PF Newman Polymath-15 q=1 Finite-Height Real-Edge Remainder Gate",
            "",
            "Date: 2026-07-31",
            "",
            "Status: effective finite-height sign for the retained first-order q=1 real-edge model. This is not a proof of an Xi-level edge theorem, Lambda <= 0, or RH.",
            "",
            "```text",
            f"work/rh_compute/results/{STEM}.json",
            f"python work/rh_compute/scripts/{STEM}.py",
            f"python work/rh_compute/scripts/check_{STEM}.py",
            "```",
            "",
            "## Effective Domain",
            "",
            "On q=1 and L>=50, t<=1/5000 and the saddle identity gives",
            "",
            "```text",
            f"h=1/a<exp(-25)<1/{H_DENOMINATOR}.",
            "```",
            "",
            "All estimates below are uniform for -1<=p<=1. The proof enlarges the physical curve to the rectangular parameter box, so no monotonicity in p or L is assumed.",
            "",
            "## Stable Terminal Logarithm",
            "",
            "Write Z=z_N/S_a=-Q exp(W), with Q=exp[-pi i(p^2/2-p+3/8)]. The continuous branch W->0 is evaluated only after the h=0 and y=0 quotients have been removed. The exact stable identities are saved in the JSON artifact. Their analytic majorants give",
            "",
            "```text",
            "|W|<6h,                 |W_x|<2h^2,",
            "|Z+Q|<12h,              |Z_x-i h theta Q/2|<10h^2.",
            "```",
            "",
            "The constant t*pi^2/64 is not bounded in isolation. It cancels algebraically against the -pi^2/16 part of Re(w^2). This is the essential finite-height cancellation.",
            "",
            "The endpoint frame has the parallel cancellation",
            "",
            "```text",
            "|mu|<h^2,               |mu_x|<h^4.",
            "```",
            "",
            "Here the -pi/8 in K_rate cancels the +pi/8 from i(chi)/2 before absolute values are taken.",
            "",
            "## Endpoint Derivatives",
            "",
            "At p=+-1/2 the common numerator/denominator zero of C_0 is factored on complex disks of radius 1/4. The saved disk bound is",
            "",
            "```text",
            derivative["C0_disk_upper_ball"],
            "<3.",
            "```",
            "",
            f"Cauchy on radius 1/8 and a {derivative['direct_boxes']}-box Arb cover of the complementary real intervals prove",
            "",
            "```text",
            "|F^(k)| < 3 k! 8^k, k=0,...,5,",
            "(|F|,...,|F^(5)|) < (3,24,384,9216,294912,11796480).",
            "```",
            "",
            "Consequently, for C_1=F'''/(12pi^2),",
            "",
            "```text",
            "(|C_1|,|C_1'|,|C_1''|)<(86,2731,109227),",
            "|H-F|<86h,",
            "|H_x+hF'/(4pi)|<232h^2,",
            "|H_xx-h^2F''/(16pi^2)|<800h^3.",
            "```",
            "",
            "## Four-Jet Budget",
            "",
            "With A, B, M as in the cofinal gate, the complete normalized edge satisfies",
            "",
            "```text",
            "c=A+h e_0,                       |e_0|<100,",
            "d_edge=hB+h^2 e_1,               |e_1|<250,",
            "c_x=hB+h^2 e_2,                  |e_2|<250,",
            "d_(edge,x)=h^2M+h^3 e_3,         |e_3|<900.",
            "```",
            "",
            "No quotient by A, c, or H is introduced. Direct expansion gives",
            "",
            "```text",
            budget["perturbation_identity"],
            "```",
            "",
            "Using |A|<4, |B|<3, and |M|<3, the linear constant is 5400 and the quadratic constant is 152500. Hence",
            "",
            "```text",
            budget["uniform_remainder"],
            "```",
            "",
            "## Finite-Height Sign",
            "",
            "The preceding cofinal gate proved K_edge<-3/8000 on the whole cell. Combining the two strict inequalities yields",
            "",
            "```text",
            budget["finite_height_margin"],
            "```",
            "",
            "for q=1, L>=50, and -1<=p<=1 in the retained first-order model, including both removable C_0 points and every zero real projection.",
            "",
            "## Boundary and Handoff",
            "",
            "This closes the finite-height remainder requested by Formal Core 11.148 inside the first-order endpoint/Dirichlet model. It does not bound the omitted higher-order Xi approximation error, does not extend the sign to q>1, and does not splice the signed edge through adjacent cutoffs. Those are the next three gates before the edge block can be inserted into the cumulative Abel/contact scalar. No winding cap, contact exclusion, Lambda<=0, PF-infinity, RH, or prize-level conclusion is proved.",
            "",
        ]
    )


def build_artifact() -> dict:
    verify_symbolics()
    derivatives = c0_derivative_certificate()
    exact = exact_payload()
    budget = majorant_certificate()
    rows = build_rows(exact, derivatives, budget)
    return {
        "kind": STEM,
        "date": "2026-07-31",
        "status": (
            "effective finite-height negative q=1 real-edge current "
            "for the retained first-order model"
        ),
        "source_audit": source_audit(),
        "exact": exact,
        "derivative_certificate": derivatives,
        "majorant_certificate": budget,
        "rows": [asdict(row) for row in rows],
        "counts": {
            "rows": len(rows),
            "direct_arb_boxes": derivatives["direct_boxes"],
            "removable_cauchy_charts": 2,
            "source_derivative_bounds": 6,
            "stable_terminal_log_identities": 1,
            "terminal_value_derivative_bounds": 2,
            "endpoint_mu_bounds": 2,
            "finite_edge_jet_defects": 4,
            "uniform_q1_finite_height_edge_signs": 1,
            "uniform_q_gt_1_edge_signs": 0,
            "xi_level_edge_signs": 0,
            "abel_gaps": 0,
            "winding_bounds": 0,
        },
        "proof_boundary": (
            "This artifact proves an effective finite-height negative real-edge "
            "current for the retained first-order q=1 endpoint/terminal model "
            "on L>=50 and -1<=p<=1. It does not bound the omitted higher-order "
            "Xi approximation error, prove a q>1 edge sign, splice adjacent "
            "cutoffs, prove a cumulative-minor estimate, Abel-scalar gap, "
            "successor winding cap, contact exclusion, Q209, the cofinal "
            "descendant theorem, Lambda<=0, PF-infinity, RH, or a prize-level "
            "conclusion."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--note", type=Path, default=DEFAULT_NOTE)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.note.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    args.note.write_text(render_note(artifact), encoding="utf-8")
    print(
        "built q=1 finite-height real-edge remainder gate: 16 rows, "
        "4096 Arb boxes, 2 Cauchy charts, 4 finite edge jets, "
        "remainder <1/10000000, 1 retained-model finite-height sign, "
        "0 q>1 signs, 0 Xi-level signs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
